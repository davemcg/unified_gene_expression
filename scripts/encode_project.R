# encode project
# https://www.encodeproject.org/search/?type=Experiment&assay_term_name=RNA-seq&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assay_title=polyA+mRNA+RNA-seq&files.file_type=fastq&status=released&files.run_type=paired-ended
load("data/encode_metadata_paired_mRNAseq.RData")
# encode must not have understood what 'paired end' was when they setup their data
# handling structure, since the forward and reverse ends have DIFFERENT
# and UNRELATED looking accession names
#
# sigh
#
# Fortunately they label f and r (1 / 2) and what the matching accession is
# for each line
# E.G ENCFF087HWA is the reverse and ENCFF365ZYO is the forward 
# so I'm going to use the forward read as the name for the kallisto folder
# E.G ENCFF365ZYO_kallisto 
# I'm going to have to mod the kallisto script to explicitly take the f and r
# reads, since it currently assumes they have the _1 _2 format that 
# everyone else uses

#  match up to main
encode$project_accession <- "ENCODE"
#colnames(encode)[which(names(encode) == "Biosample term name")] <- "comment_ena_experiment"

forward <- encode %>% filter(`Paired end`==1)
reverse <- encode %>% filter(`Paired end`==2)
merged <- left_join(forward,reverse[,c("Paired with","File accession","File download URL")],by=c("File accession"="Paired with"))


