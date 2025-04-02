library(data.table)

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

# define the samples to read. For now, it's best to list a control first.
sample_list = c("sample5", "sample6")
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
sample_desc = list()
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

polyA <- split_contig(polyA)

# Set up to resample an equal number of reads per sample type per contig from the polyA data.
# summary counts...

polyA_counts <- polyA[, (count = .N), by = .(contig, sample_label)]

# start from an empty list
dt_list <- list()

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
}

# smoosh the resampled list back into a single dataframe
polyA_resampled <- rbindlist(dt_list)

# clean up the temps
rm(dt, dt_list)

# save the RDS. Usually we want the resampled
rds_name <- paste(paste(sample_list, collapse = "_"), "PolyA_resampled.rds", sep ="_")
saveRDS(polyA_resampled, file=rds_name)