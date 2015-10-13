# Lobster_checkin

This started off with me wanting to check some results that Louis' group had
on Lobster population genetics in this paper: 
[RAD genotyping reveals fine-scale genetic structuring and provides powerful population assignment in a widely distributed marine species, the American lobster (Homarus americanus)](http://onlinelibrary.wiley.com/doi/10.1111/mec.13245/abstract).

I found it strange that Benestan et al found that their power for assignment decreased as they added
more markers past 3000.  This seems worrisome to me because it is a pattern characteristics of
high grading bias.  Here I use my own program `gsi_sim` to implement the Training-Holdout-Leave-one-out
(THL) procedure from 
[Anderson 2010](http://onlinelibrary.wiley.com/doi/10.1111/j.1755-0998.2010.02846.x/full) 
to assess power for assignment of lobsters to their sampling
locations of origin.  And I don't find as much power for that as Benestan and
colleagues reported.  It will be interesting to try to track down what is going on
here.  I might have a bug, or they might have a better locus selection procedure,
or they might have suffered some high-grading bias in the procedure they used, even
though they intended to eliminate any such bias.

Note that if you want to reproduce these results you will want to be on a Mac or a Unix system.
You will need to have `vcftools` installed and in your `$PATH` variable which is given in a
`~/.bashrc` file.  Likewise you will need `gsi_sim`
([https://github.com/eriqande/gsi_sim](https://github.com/eriqande/gsi_sim) commit 280734b4 or later)
compiled and on your `$PATH`.

To reproduce my work, clone this repo, install `gsi_sim` and `vcftools` and put them on your 
`$PATH` as listed in your `~/.bashrc`, then you have to download the file `10156-586.recode.vcf` from the
Benestan et al. paper's Dryad site and put it into
`./data/10156-586.recode.vcf`and then open the [Rstudio](https://www.rstudio.com/)
project and knit the file `lobster_checkin.Rmd` to an HTML within Rstudio.
