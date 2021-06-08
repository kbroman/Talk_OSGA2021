# plot of missing data
library(qtl2)
library(broman)

attie <- readRDS("../Data/attieDO_v1.rds")
pmis <- n_missing(attie, "ind", "prop")*100

pdf("../Figs/missing_data.pdf", height=6, width=10)
par(mar=c(5.1,4.1, 1.1, 1.1))
col <-c("pink", "orange", "lightblue", "green", "purple")[as.numeric(attie$covar$wave)]
grayplot(pmis, xlab="Mouse", ylab="Percent missing genotypes",
         yaxs="i", ylim=c(0, max(pmis)*1.05), bg=col,
         hlines=seq(10, 70, by=10))
text(which(pmis > 10)-5, pmis[pmis > 10], names(pmis)[pmis > 10], adj=c(1, 0.5))
dev.off()
