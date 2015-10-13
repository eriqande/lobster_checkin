library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)



#### First, deal with the different populations ####
# we want a long data frame of the indiv IDs with a column that says where they
# belong
#Inds <- read.table("data/10156-indivs.txt", stringsAsFactors = FALSE) %>%
#  tbl_df %>%
#  setNames(c("ID", "GrpNum")) %>%
#  mutate(SampSite = str_sub(ID, 1, 3))

# But, actually, those seem to include individuals that are not in the file...
# so, instead, we just make Inds from the individuals we find in the file.
# Let's compute their mean depth and then use the output to get all the individuals.
system("source ~/.bashrc; vcftools --vcf data/10156-586.recode.vcf --out intermediates/get-indivs --depth")

Inds <- read.table("intermediates/get-indivs.idepth", stringsAsFactors = FALSE, header = TRUE) %>%
  tbl_df %>%
  mutate(SampSite = str_sub(INDV, 1, 3))




#### Now, we want a function that takes the SampSite name and then ####
# gets just those individuals and grabs out the genotype counts for them
# and returns it in a long data frame that has a column for SampSite.

# To do that we make a few functions

# do the HW calc with vcftools on individuals from a given SampSite.
# Note hack to get my PATH variable in there.
vcf_hardy <- function(SS, Inds) {
  # first get the file with the indivs in it
  Inds %>% 
    filter(SampSite == SS) %>% 
    select(INDV) %>% 
    write.table(., file = file.path("intermediates", SS), row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  CALL <- paste("source ~/.bashrc; vcftools --vcf data/10156-586.recode.vcf ",
                "--hardy ",
                "--keep", file.path("intermediates", SS),
                "--out", file.path("intermediates", SS)
  )

  system(CALL)
}

#### A function to read vcf hardy output and put it into usable long format ####
long_hardy <- function(file) {
  x <- read.table(file, sep ="\t", header = TRUE, stringsAsFactors = FALSE) %>%
    tbl_df %>%
    separate(OBS.HOM1.HET.HOM2., c("obs_Hom1", "obs_Het", "obs_Hom2"), sep ="/") %>%
    separate(E.HOM1.HET.HOM2., c("exp_Hom1", "exp_Het", "exp_Hom2"), sep ="/")
  
  # now put it into long format
  obs <- x %>% select(CHR, POS, starts_with("obs_")) %>%
    gather(., obs_var, obs_cnt, starts_with("obs_")) %>%
    mutate(geno = str_replace(obs_var, "^obs_", "")) %>%
    mutate(obs_cnt = as.numeric(obs_cnt)) %>%
    select(-obs_var)
  
  # now add the total number of observed individuals in there
  ntot <- obs %>% group_by(CHR, POS) %>% summarise(nindiv = sum(obs_cnt))
  
  exp <- x %>% select(CHR, POS, starts_with("exp_")) %>%
    gather(., exp_var, exp_cnt, starts_with("exp_")) %>%
    mutate(geno = str_replace(exp_var, "^exp_", "")) %>%
    mutate(exp_cnt = as.numeric(exp_cnt)) %>%
    select(-exp_var)
  
  
  cnts <- inner_join(obs, exp) %>%
    mutate(geno = factor(geno, levels = c("Hom1", "Het", "Hom2"))) %>%
    inner_join(., ntot)
  
  list(cnts = cnts, p_etc = inner_join(x, ntot))
  
}



# with those two functions we can make a single one that gets what we want:
get_geno_cnts_on_pop <- function(SS, Inds) {
  vcf_hardy(SS, Inds)
  
  long_hardy(file.path("intermediates", paste(SS, "hwe", sep = ".")))$cnts %>%
    mutate(SampSite = SS)
}


# and finally to do that over all the populations we can do this:
tmp <- lapply(unique(Inds$SampSite), function(SS) get_geno_cnts_on_pop(SS, Inds))

DF <- bind_rows(tmp)

# plot it all up
set.seed(555)
big_plot <- ggplot(DF, aes(x = exp_cnt, y = obs_cnt, colour = geno)) +
  geom_jitter(alpha = 0.15) +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(SampSite ~ geno)

#ggsave(big_plot, filename = "geno_exp_v_obs.png", width = 14, height = 30)


