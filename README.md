# Epitranscriptome Explorer

A shiny app for the exploration of RNA epitranscriptomic data (polyadenylation
and m6A/m5C methylation). For now, the input data for this app must be preprocessed 
through the accompanying `preprocess_polya_methyl.R` script. 
This script combines two (or more) sample runs into RDS files representing 
merged and filtered data.tables, one for PolyA length and one for m6A/m5C methylation.
The app acts independently on PolyA and methylation data, although in most cases 
we would be comparing data from the same sample runs.

The starting point for the app are the two upload buttons for PolyA and methylation
RDS files. These upload client side RDS files to a temp location on the server. 
The app may be set to pre-load existing files on the server side through the 
`config.json` configuration file. Note that if the shiny server is running locally,
`config.json` may point to local files, however if the server is remote, the only
way shiny can access local user files is through the upload buttons. The max upload
size is set fairly large, but it may be possible that some files are just too 
large (the app tends to start to struggle with two >500MB RDSs).

The filter section subsets the data by transcript type (eg. protein coding, rRNA,
unprocessed pseudogene), particular genes (a simple comma separated list of gene
names may be uploaded) or particular transcripts. The transcript filter only becomes
active when the data has already been filtered down to a few genes. For these
input selectors, empty or blank should be interpreted as "no filter"/everything.



