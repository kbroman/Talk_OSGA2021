# compare geno histogram
library(here)
library(qtl2)

do <- readRDS(here("Data/attieDO_v1.rds"))
cg <-compare_geno(do, cores=0)

prop_match_fig <- function(label=FALSE)
{
    par(mar=c(6.1, 1.1, 1.1, 1.1))
    cg_upper <- cg[upper.tri(cg)]*100
    hist(cg_upper, breaks=seq(0, 100, length=201),
         xlab="Percent matching genotypes", ylab="", yaxt="n", main="")
    rug(cg_upper, col="violetred")

    if(!label) return()

    textpos <- median(cg_upper[cg_upper < 35])
    text(textpos, 4000, "DO306 / DO308", col="darkslateblue")
    arrows(textpos, 3000, textpos, 1000, length=0.15, col="violetred", lwd=2)
}

pdf(here("Figs/hist_compare_geno.pdf"), height=6, width=10, pointsize=14)
prop_match_fig(TRUE)
dev.off()
