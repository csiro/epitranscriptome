library(data.table)
library(progress)

# convenience function to save writing print(paste0()) all the time
pprint <- function(...){
  print(paste0(...))
}

# split the "contig" column into its various components in a data.table
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

get_metadata_dt <- function(metadata_path){
  if (!file.exists(metadata_path)){
    pprint("metadata file ", metadata_path, " not found.")
    stop()
  }
  
  # opening the file with check.names=TRUE converts spaces to '.'
  # in column names.
  metadata_dt <- fread(metadata_path, check.names = TRUE)
  expected_cols <- c("Cell.line", "Virus", "Time", "Group", "Sample_ID", "Study", "polya_path", "m6a_path", "m5c_path")
  
  if (!all(expected_cols %in% colnames(metadata_dt))){
    print("metadata file does not contain expected columns:")
    print(expected_cols)
    stop()
  }
  return (metadata_dt)
}

make_sample_desc <- function(metadata_dt, sample_id){
  
  # construct a description string from the metadata, adjust or override as needed.
  # We're assuming the metadata file has one entry per sample, and the sample IDs
  # are of the form "sampleX"

  s.cell_line <- metadata_dt[Sample_ID == sample_id, Cell.line]
  s.virus <- metadata_dt[Sample_ID == sample_id, Virus]
  s.group <- metadata_dt[Sample_ID == sample_id, Group]
  s.time <- metadata_dt[Sample_ID == sample_id, Time]
  
  s.desc <- paste(s.cell_line, s.virus, s.time, s.group, sep = "-")
  
  if (length(s.desc) < 1){
    pprint("No metadata found for ", sample_id)
    stop()
  }
  return (s.desc)
}

preprocess_polya <- function(metadata_dt, sample_list){
  
  # qc_tags to keep
  qc_pass = c("PASS")
  
  dt_list = list()
  for (s in sample_list){
    # filename based on sample
    polyA_file <- metadata_dt[Sample_ID == s, polya_path]
    
    if (!file.exists(polyA_file)){
      pprint("PolyA file ", polyA_file, " not found")
      stop()
    }
    cat("Reading", polyA_file, "...\n")
    # fread will throw if the file isn't found
    dt <- fread(polyA_file)
    # filter only the qc_tags in out qc_pass list
    dt <- dt[qc_tag %in% qc_pass]
    # append the sample description
    sample_desc <- make_sample_desc(metadata_dt, s)
    dt[, sample_label := sample_desc]
    # append to the data.table list
    dt_list = c(dt_list, list(dt))
  }
  
  # we're assuming the columns are the same for the polyA files,
  # rbindlist is quicker than rbind.
  polyA <- rbindlist(dt_list)
  
  # clean up
  rm(dt_list, dt)
  
  # split the contig column into new columns for gene_id etc
  polyA <- split_contig(polyA)
  
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
        sample_desc <- make_sample_desc(metadata_dt, s)
        dt <- polyA[sample_label == sample_desc & contig == t]
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
}

prepare_methyl_dt <- function(dt){
  # check for the exact (new) format from r2dtool
  expected_columns <- c("chromosome","start","end","name","score","strand","transcript","start.1","end.1","name.1","score.1","strand.1","motif","coverage","stoichiometry","probability","gene_id","gene_name","transcript_biotype","tx_len","cds_start","cds_end","tx_end","transcript_metacoordinate","abs_cds_start","abs_cds_end","up_junc_dist","down_junc_dist")
  
  ret <- data.table()
  
  if (!all(expected_columns %in% colnames(dt))){
    print("incoming methyl data table does not have the expected columns")
    print(expected_columns)
  } else {
    # drop columns we will be extracting from the transcript==contig field
    ret <- dt[, c("gene_id", "gene_name", "transcript_biotype" ) := NULL]
    # rename transcript <=> contig and start.1 <=> position
    setnames(ret, c("transcript", "start.1"), c("contig", "position")) 
  }
  return (ret)
}

preprocess_methyl <- function(metadata_dt, sample_list, coverage_thres = 0){
  
  dt_list = list()

  for (s in sample_list){
    sample_desc <- make_sample_desc(metadata_dt, s)
    
    m6a_file <- metadata_dt[Sample_ID == s, m6a_path]
    m5c_file <- metadata_dt[Sample_ID == s, m5c_path]
    
    cat("Reading", m6a_file, "...\n")
    m6a_dt <- fread(m6a_file, check.names = TRUE)
    
    m6a_dt <- prepare_methyl_dt(m6a_dt)
    
    if (nrow(m6a_dt) > 0){
      m6a_dt <- m6a_dt[coverage > coverage_thres, ]
      m6a_dt <- m6a_dt[, sample_label := sample_desc]
      m6a_dt <- m6a_dt[, meth_type := "m6A"]
    } else {
      pprint("Empty m6a file ", m6a_file)
      stop()
    }
    
    cat("Reading", m5c_file, "...\n")
    m5c_dt <- fread(m5c_file, check.names = TRUE)
    
    m5c_dt <- prepare_methyl_dt(m5c_dt)
    
    if (nrow(m5c_dt) > 0){
      m5c_dt <- m5c_dt[coverage > coverage_thres, ]
      m5c_dt <- m5c_dt[, sample_label := sample_desc]
      m5c_dt <- m5c_dt[, meth_type := "m5C"]
    } else {
      pprint("Empty m5c file ", m5c_file)
      stop()
    }
    dt_list <- c(dt_list, list(m6a_dt, m5c_dt))
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
  
  methyl <- methyl_raw
  
  # maximum gap from one pos to the next
  methyl <- methyl[, max_gap := max(pos_gap), by = .(contig, meth_type, sample_label)]
  
  # probabilty to call methylated > prob_thres or unmethylated < prob_thres
  prob_thres <- 0.99
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
  
  # Write methyl as an RDS
  meth_name <- paste(paste(sample_list, collapse = "_"), "methyl.rds", sep ="_")
  saveRDS(methyl, file=meth_name)
}

# why not both?
preprocess_polya_methyl <- function(metadata_path, sample_list){
  metadata_dt <- get_metadata_dt(metadata_path)
  preprocess_polya(metadata_dt, sample_list)
  preprocess_methyl(metadata_dt, sample_list)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1){
  metadata_path <- as.character(args[1])
  sample_ids <- list()
  for (idx in 2:length(args)){
    sample_ids <- append(sample_ids, list(as.character(args[idx])))
  }
  preprocess_polya_methyl(metadata_path, sample_ids)
} else {
  print("Need at least 2 arguments: metadata file path, and one or more Sample_ID's referenced in the metadata file.")
  print("Metadata file should be a .csv containing fields: Cell line, Virus, Time, Group, Sample_ID, Study, polya_path, m6a_path, m5c_path")
  print("example line: Caco, SARS-Cov-2, 24h, Control, sample5, Chang2021, path/to/polya_file, path/to/m6a_file, path/to/m5c_file")
  print("See app README for more information.")
}