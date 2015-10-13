

# this should be run after running 01-check-hw-props.R since it needs some of the text files
# and variables that produces.

# first  make all the --weir-fst-pop commands
fst_comms <- paste("--weir-fst-pop", file.path("intermediates", unique(Inds$SampSite)), collapse = " ")

# then compute the Fst for all the markers
CALL <- paste("source ~/.bashrc; vcftools --vcf data/10156-586.recode.vcf --out intermediates/fst-all", fst_comms)
system(CALL)


# Now, read those in and plot them, etc
FST <- read.table("intermediates/fst-all.weir.fst", header = TRUE, stringsAsFactors = FALSE) %>%
  tbl_df


ggplot(FST, aes(x = WEIR_AND_COCKERHAM_FST)) + geom_histogram(binwidth = 0.001)
