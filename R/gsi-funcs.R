

# how about a function to turn a VCF file into a gsi_sim file

# it will just take the name of the file
# note, while developing, set VCF to boing.recode.vcf
vcf2gsi_sim <- function(VCF, outf = "gs_baseline.txt") {
  # first output to a matrix of -1, 0, 1, and 2
  CALL <- paste("source ~/.bashrc; vcftools --vcf", VCF, "--out intermediates/tmp-one-two --012")
  system(CALL)
  
  # get the indivs
  gsIndivs <- scan("intermediates/tmp-one-two.012.indv", what = "character")
  N <- length(gsIndivs)
  
  # get the markers
  tmp <- scan("intermediates/tmp-one-two.012.pos", what = "character")
  gsMarkers <- paste("X", tmp[c(T,F)], tmp[c(F,T)], sep = "_")
  M <- length(gsMarkers)
  
  # now make a hashie thing to change some allele names and make it all two-columny
  alles <- c("1 1", "1 2", "2 2", "-1 -1")
  names(alles) <- c("0", "1", "2", "-1")
  
  # then read them in and make a matrix of "two-column genotypes
  genos <- alles[scan("intermediates/tmp-one-two.012", what = "character")] %>% 
    matrix(., nrow = N, byrow = TRUE)  
  
  rownames(genos) <- gsIndivs
  genos <- genos[,-1]  # note that we tear off the column that gives the index of each individual
  
  # now we need to sort the rows so that populations are together in case they arent already grouped
  # togehter.  Just a simple order on the rownames.
  genos <- genos[order(rownames(genos)),]
  
  
  # once we have done that we should modify the row names so that the first individual of each 
  # population has "POP NameOfPop\n" in its name, so we can print those out as they need to be
  # for gsi_sim.  Note that this is specific to the naming format on the lobster data
  pops <- str_sub(rownames(genos), 1, 3)
  firsties <- c(1, 1 + which(diff(as.numeric(factor(pops))) != 0))
  rownames(genos)[firsties] <- paste("POP ", pops[firsties], "\n", rownames(genos)[firsties], sep = "")
  
  
  # then we just spooge out the preamble and the genotypes.
  cat(nrow(genos), " ", ncol(genos), "\n", sep = "", file = outf)
  cat(gsMarkers, sep = "\n", file = outf, append = TRUE)
  write.table(genos, col.names = FALSE, quote = FALSE, sep = "    ", file = outf, append = TRUE)
  
  # return the name of the file written to:
  return(outf)
}


# from gsi_sim --self-assign output, slurp out the population each individual
# got assigned to, and note whether it is correct or not.
# this relies on having pop names that are exactly three letters long.
slurp_ass <- function(infile = "gs_out.txt") {
  tmp <- readLines(infile)
  tmp2 <- tmp[str_detect(tmp, "^SELF_ASSIGN_A_LA_GC_CSV")] %>%
    str_replace_all(., "^SELF_ASSIGN_A_LA_GC_CSV:/", "")
  str_split(tmp2, ";") %>% 
    sapply(., "[", 1:3) %>% 
    t %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    tbl_df %>%
    setNames(c("ID", "To", "Score")) %>%
    mutate(From = str_sub(ID, 1, 3))
}