

source("R/gsi-funcs.R")


#### Randomly choose half (as close as possible) of the individuals as Training and Test ####
set.seed(345) # for reproducibility
Split <- Inds %>%
  group_by(SampSite) %>%
  mutate(Pool = ifelse(sample(1:length(INDV)) %% 2 == 0, "Test", "Train")) %>%
  ungroup


# now make VCF files for each of those
# first we make files with the names of the individuals
Split %>% 
  filter(Pool == "Test") %>% 
  select(INDV) %>% 
  unlist %>% 
  unname  %>%
  cat(., file = "intermediates/test-guys.txt", sep = "\n")
  

Split %>% 
  filter(Pool == "Train") %>% 
  select(INDV) %>% 
  unlist %>% 
  unname  %>%
  cat(., file = "intermediates/train-guys.txt", sep = "\n")

# then we run vcftools to grab those out:
system("source ~/.bashrc; vcftools --vcf data/10156-586.recode.vcf --keep intermediates/train-guys.txt --out intermediates/train --recode")
system("source ~/.bashrc; vcftools --vcf data/10156-586.recode.vcf --keep intermediates/test-guys.txt --out intermediates/test --recode")

# Compute Fst's from the training set and store those in a variable:
fst_comms <- paste("--weir-fst-pop", file.path("intermediates", unique(Inds$SampSite)), collapse = " ")
CALL <- paste("source ~/.bashrc; vcftools --vcf intermediates/train.recode.vcf --out intermediates/train-recode", fst_comms)
system(CALL)


# get a table of marker positions sorted by Fst high to low
train_FST <- read.table("intermediates/train-recode.weir.fst", header = TRUE, stringsAsFactors = FALSE) %>%
  tbl_df %>%
  arrange(desc(WEIR_AND_COCKERHAM_FST))


#### Now a function to do THL on a given number, nL of loci.  ####
thl <- function(nL, Split, vcffile = "data/10156-586.recode.vcf") {
  # now, choose the first nL loci, write them to a file, and get all the individuals'
  # genotypes at those loci
  fn <- file.path("intermediates", "tmp_top_loci")
  train_FST[1:nL, ] %>%
    select(CHROM, POS) %>%
    write.table(., row.names = FALSE, col.names = FALSE, sep = "\t",
                quote = FALSE, file = fn)
  
  
  CALL <- paste("source ~/.bashrc; vcftools --vcf data/10156-586.recode.vcf --positions",
                fn, 
                "--recode --out",
                fn
  )
  
  system(CALL)
  
  
  # make a gsi_sim file of that:
  vcf2gsi_sim(paste(fn, "recode.vcf", sep = "."), 
              outf = paste(fn, "gsi_sim_input", sep = "."))
  
  # then run gsi_sim on it
  system(paste("source ~/.bashrc; /Users/eriq/Documents/git-repos/gsi_sim/gsisim -b", paste(fn, "gsi_sim_input", sep = "."),
               "--self-assign > intermediates/bigdump"))
  
  
  # and slurp out the assignments
  result <- slurp_ass(file.path("intermediates", "bigdump"))
  
  # join back onto these whether or not they were Train or Test and add a "Correct" column
  # and also a number of loci column.  And return that data frame.
  result %>%
    left_join(Split, ., by = c("INDV" = "ID")) %>%
    mutate(Correct = (From == To),
           NumLoci = nL)
}


#### Now that we have that function, we can make a list of results
LocNums <- c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 10156)
ResList <- lapply(LocNums, function(x) thl(x, Split))

# then make a big tidy data frame of it:
ResDF <- bind_rows(ResList)

# compute the proportion correct from each population under each scenario of
# NumLoci and Train/Test
SumDF <- ResDF %>%
  group_by(NumLoci, From, Pool) %>%
  summarize(NumCorrect = sum(Correct),
            Ppn_Correct = NumCorrect / n())


# now plot those dudes
ggplot(SumDF, aes(x = factor(NumLoci), y = Ppn_Correct, fill = Pool)) +
  geom_boxplot()


ggsave(file = "lobster_checking_boxplots.pdf", width = 10, height = 7)
