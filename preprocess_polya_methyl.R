library(data.table)
library(progress)

################################################################################
# Input Directory setup

# start by gleening the environment
platform <- Sys.info()['sysname']

if (platform == 'Linux'){
  # we'll assume we're on the HPC if we're in a linux environ
  root_path <- "/datasets/work/hb-rna-modifications/work/imt_sc_rnamodifications/"
  # the path to the annotations from r2d may be different
  r2d_root_path <- "/datasets/work/hb-rna-modifications/work/imt_sc_rnamodifications/r2d/"
  #r2d_root_path <- "/datasets/work/sc-kno070/work/rna-transcriptomics/"
} else {
  # otherwise, for me this mapped network drive points to the same Bowen store. 
  root_path <- "Y:/bk_imt_rnamodifications/"
  # my local mapped r2d annotations directory
  #r2d_root_path <- "Y:/r2d/"
  r2d_root_path <- "Y:/bk_imt_rnamodifications/r2d/"
}

# paths from the root(s)
polyA_path <- paste0(root_path, "8b.nanopolish_polya/polya_out/")
metadata_file <- paste0(root_path, "dRNA_seq_metadata_design.csv")

# opening the file with check.names=TRUE converts spaces to '.'
# in column names.
seq_metadata <- fread(metadata_file, check.names = TRUE)

################################################################################
# Define the samples of interest

# define the samples to read. For now, it's best to list the control first.
# It is possible to define more than two sample sets here, but it's probably best
# not to try more than about 4 or 5...
sample_list <- c("sample25", "sample26")
# use the r2d annotated methylation calls (if available) otherwise 
# use the fast5 files
# TODO: this should be determined by the files existing in the r2d_annotation dirs
use_r2d_annotated <- TRUE

# reset the paths if we're using the annotated versions
if (use_r2d_annotated){
  m5c_path <- paste0(r2d_root_path, "m5c_r2d_annotation/")
  m6a_path <- paste0(r2d_root_path, "m6a_r2d_annotation/")
} else {
  m5c_path <- paste0(root_path, "m5c_solo_output/")
  m6a_path <- paste0(root_path, "m6a_solo_output/")
}

# construct the descriptions from the metadata, adjust or override as needed.
# We're assuming the metadata file has one entry per sample, and the sample IDs
# are of the form "sampleX_fast5"
sample_desc <- list()
for (s in sample_list){
  s.id <- paste0(s, "_fast5")
  s.cell_line <- seq_metadata[Sample_ID == s.id, Cell.line]
  s.virus <- seq_metadata[Sample_ID == s.id, Virus]
  s.group <- seq_metadata[Sample_ID == s.id, Group]
  s.time <- seq_metadata[Sample_ID == s.id, Time]
  
  s.desc <- paste(s.cell_line, s.virus, s.time, s.group, sep = "-")
  
  if (length(s.desc) < 1){
    stop("No metadata for ", s.id)
  }
  
  sample_desc = c(sample_desc, list(s.desc))
}

# make an R dictionary for the descriptions so we can retrieve the descriptions
# with sample_desc[["sampleX"]]
names(sample_desc) <- sample_list

# a function to split the "contig" column into its various components
split_contig <- function(dt){
  if (is.data.table(dt) & "contig" %in% colnames(dt)){
    dt <- dt[, c("ensembl_t", 
                 "ensembl_g", 
                 "havana_g", 
                 "havana_t", 
                 "transcript_id", 
                 "gene_id", 
                 "length",
                 "transcript_type") := tstrsplit(contig, "|", fixed=TRUE) ]
    
    dt <- dt[, length := as.integer(length)]
  }
  return (dt)
}

################################################################################
# Read and filter PolyA data 

# qc_tags to keep
qc_pass = c("PASS")

dt_list = list()
for (s in sample_list){
  # filename based on sample
  polyA_file <- paste0(polyA_path, s, "_polya_result.txt")
  cat("Reading", polyA_file, "...\n")
  # fread will throw if the file isn't found
  dt <- fread(polyA_file)
  # filter only the qc_tags in out qc_pass list
  dt <- dt[qc_tag %in% qc_pass]
  # append the sample description
  dt[, sample_label := sample_desc[[s]]]
  # append to the data.table list
  dt_list = c(dt_list, list(dt))
}

# we can assume the columns are the same for the polyA files,
# rbindlist is quicker than rbind.
polyA <- rbindlist(dt_list)

# clean up
rm(dt_list, dt)

# split the contig column into new columns for gene_id etc
polyA <- split_contig(polyA)

################################################################################
# Resample PolyA

# Set up to resample an equal number of reads per sample type per contig from the polyA data.
# summary counts...

polyA_counts <- polyA[, (count = .N), by = .(contig, sample_label)]

# start from an empty list
dt_list <- list()

cat("Resampling...\n")
total_contigs <- length(unique(polyA[, contig]))
pb <- progress_bar$new(total = total_contigs)

for (t in unique(polyA[, contig])){
  # for each transcript type, get the min count over all the sample types
  # (opt) take half for the resample number
  n_t <- min(polyA_counts[contig == t][,V1])# %/% 2
  
  # if there was a min, and there are some for each sample type
  if (n_t > 0 & nrow(polyA_counts[contig == t]) == length(sample_list)){
    for (s in sample_list){
      # subset to this sample label
      dt <- polyA[sample_label == sample_desc[[s]] & contig == t]
      # take n_t rows from the dt and append to the dt list
      dt_list <- c(dt_list, 
                   list(dt[sample(.N, size = n_t, replace = FALSE)]))
    }
  }
  
  pb$tick()
}

# smoosh the resampled list back into a single dataframe
polyA_resampled <- rbindlist(dt_list)

# clean up the temps
rm(dt, dt_list)

# save the RDS. Usually we want the resampled data, but the raw data can be output
# with saveRDS(polyA, file=...)
rds_name <- paste(paste(sample_list, collapse = "_"), "PolyA_resampled.rds", sep ="_")
saveRDS(polyA_resampled, file=rds_name)

################################################################################
# Read in the methylation data

# drop low coverage sites
coverage_thres <- 5

prepare_dt <- function(dt){
  # check for the exact (new) format from r2dtool
  expected_columns <- c("chromosome","start","end","name","score","strand","transcript","start.1","end.1","name.1","score.1","strand.1","motif","coverage","stoichiometry","probability","gene_id","gene_name","transcript_biotype","tx_len","cds_start","cds_end","tx_end","transcript_metacoordinate","abs_cds_start","abs_cds_end","up_junc_dist","down_junc_dist")
  
  ret <- data.table()
  
  if (!all(expected_columns %in% colnames(dt))){
    print("incoming methyl data table does not have the expected columns")
  } else {
    # drop columns we will be extracting from the transcript==contig field
    ret <- dt[, c("gene_id", "gene_name", "transcript_biotype" ) := NULL]
    # rename transcript <=> contig and start.1 <=> position
    setnames(ret, c("transcript", "start.1"), c("contig", "position")) 
  }
  return (ret)
}

dt_list = list()

if (use_r2d_annotated){
  for (s in sample_list){
    m6a_file <- paste0(m6a_path, "r2d_", s, "/", s, "_methylationCalls_annotated_lifted.bed")
    m5c_file <- paste0(m5c_path, "r2d_", s, "/", s, "_methylationCalls_annotated_lifted.bed")
    
    cat("Reading", m6a_file, "...\n")
    m6a_dt <- fread(m6a_file, check.names = TRUE)
    
    m6a_dt <- prepare_dt(m6a_dt)
    
    if (nrow(m6a_dt) > 0){
      m6a_dt <- m6a_dt[coverage > coverage_thres]
      m6a_dt <- m6a_dt[, sample_label := sample_desc[[s]]]
      m6a_dt <- m6a_dt[, meth_type := "m6A"]
    }
    
    cat("Reading", m5c_file, "...\n")
    m5c_dt <- fread(m5c_file, check.names = TRUE)
    
    m5c_dt <- prepare_dt(m5c_dt)
    
    if (nrow(m5c_dt) > 0){
      m5c_dt <- m5c_dt[coverage > coverage_thres]
      m5c_dt <- m5c_dt[, sample_label := sample_desc[[s]]]
      m5c_dt <- m5c_dt[, meth_type := "m5C"]
    }
    dt_list <- c(dt_list, list(m6a_dt, m5c_dt))
  }
} else {
  for (s in sample_list){
    m6a_file <- paste0(m6a_path, s, "_site_level_m6A_predictions.txt")
    m5c_file <- paste0(m5c_path, s, "_site_level_m5C_predictions.txt")
    
    cat("Reading", m6a_file, "...\n")
    m6a_dt <- fread(m6a_file)
    m6a_dt <- m6a_dt[coverage > coverage_thres]
    m6a_dt <- m6a_dt[, sample_label := sample_desc[[s]]]
    m6a_dt <- m6a_dt[, meth_type := "m6A"]
    
    cat("Reading", m5c_file, "...\n")
    m5c_dt <- fread(m5c_file)
    m5c_dt <- m5c_dt[coverage > coverage_thres]
    m5c_dt <- m5c_dt[, sample_label := sample_desc[[s]]]
    m5c_dt <- m5c_dt[, meth_type := "m5C"]
    
    dt_list <- c(dt_list, list(m6a_dt, m5c_dt))
  }
}

methyl_raw <- rbindlist(dt_list)

# clean up the temps
rm(dt_list, m5c_dt, m6a_dt)

# as before, add columns for the split descriptors in the contig field
methyl_raw <- split_contig(methyl_raw)

# get the actual site read span as the max - min position, per contig per meth_type per sample_label.
# This read length may be vastly different between the samples and/or the length of the transcript given in the
# contig info.
methyl_raw <- methyl_raw[, pos_span := max(position) - min(position), by = .(contig, meth_type, sample_label)]

# just the number of sites (A or C, methylated or not) per contig per meth_type per sample_label.
# this will give us a good idea of overall coverage of the contigs.
methyl_raw <- methyl_raw[, site_count := .N, by = .(contig, meth_type, sample_label) ]

# we also want to measure the "gap" between one site and the next. for this we want to make sure the
# data table is ordered by position (per contig per meth_type per sample_label of course)
setorder(methyl_raw, contig, meth_type, sample_label, position)

# a function to measure inter-row gaps  
gap = function(x, ...){ max(x) - min(x) }
methyl_raw <- methyl_raw[, pos_gap := frollapply(position, 2, gap, fill = 0), by = .(contig, meth_type, sample_label)]

################################################################################
# Filter the methylation data and compute densities

# drop reads with less than this many sites
site_count_thres <- 50
methyl <- methyl_raw[site_count > site_count_thres]

# drop methyl reads where the read_length is less than x% of the actual length from the contig info
#methyl <- methyl[pos_span > length * 0.8]

# maximum gap from one pos to the next
methyl <- methyl[, max_gap := max(pos_gap), by = .(contig, meth_type, sample_label)]

# drop gappy contigs here
max_gap_thres <- 100
#methyl <- methyl[max_gap < max_gap_thres]

# probabilty to call methylated > prob_thres or unmethylated < prob_thres
prob_thres <- 0.8
# now use the probability threshold to just mark 1's for probably methylated, 0's for probably not 
methyl <- methyl[, methylated := ifelse( (probability >= prob_thres), 1, 0)]
#methyl <- methyl[, filter := ifelse( (probability >= prob_thres), "sig", "ns")]

# window length (per row not per position) for rolling totals/averages
rolling_window_width <- 50
# as far as i know, frollmean/sum/apply can only work on one column at a time, here we're building
# up to get a measure of the rolling densities along each contig
methyl <- methyl[, rolling_pos_gap := frollsum(pos_gap, rolling_window_width, fill = 0, align = "center"), by = .(contig, meth_type, sample_label)]

methyl <- methyl[, rolling_meth_count := frollsum(methylated, rolling_window_width, fill = 0, align = "center"), by = .(contig, meth_type, sample_label)]

methyl <- methyl[, rolling_site_count := (rolling_window_width / rolling_pos_gap)]

methyl <- methyl[, rolling_meth_density := (rolling_meth_count / rolling_pos_gap)]

methyl <- methyl[, rolling_meth_density_normed := (rolling_meth_count / rolling_site_count)]

# use site_count rather than max(position) - min(position) or contig length - there are plenty of incomplete contigs
methyl <- methyl[, total_meth_density := sum(methylated) / site_count, by = .(contig, meth_type, sample_label) ]

if (use_r2d_annotated){
  # metacoordinate binning ... like R2_plotMetaTranscript
  # subdivide the [0,3] interval
  metacoord_breaks <- seq(0, 3, 0.025)
  
  # out_ratio <- methyl[, interval := cut(transcript_metacoordinate, metacoord_breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)][, .(n =  .N), by = .(interval, meth_type, sample_label)][, dcast(.SD, interval ~ filter, value.var = "n", fill = 0)][, ratio := sig / (sig + ns)]
  
  # add a new categorical interval label based on metacoordinate
  out_ratio <- methyl[, metacoord_interval := cut(transcript_metacoordinate, metacoord_breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)]
  
  # aggregate by metacoordinate, meth_type, sample_label (ie. over all contigs)
  out_ratio <- out_ratio[, .(n = .N, sig = sum(methylated)), by = .(metacoord_interval, meth_type, sample_label)]
  
  out_ratio <- na.omit(out_ratio)
  
  #out_ratio <- out_ratio[, metacoord_sig_sites := sum(methylated), by = .(metacoord_interval, meth_type, sample_label)]
  
  out_ratio <- out_ratio[, metacoord_sig_ratio := sig / n]
}

################################################################################
# Write filtered methylation as an RDS

meth_name <- paste(paste(sample_list, collapse = "_"), "methyl.rds", sep ="_")
saveRDS(methyl, file=meth_name)