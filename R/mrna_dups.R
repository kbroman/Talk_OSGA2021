mrna <- readRDS("../Data/chick_mrna.rds")
sex <- substr(rownames(mrna), 1, 1)

mrna_adj <- mrna
for(x in c("M", "F")) {
    ave <- colMeans(mrna[sex==x,], na.rm=TRUE)
    mrna_adj[sex==x,] <- t(t(mrna[sex==x,]) - ave)
}

rho <- cor(t(mrna), use="pair")
rho_adj <- cor(t(mrna_adj), use="pair")


png("../Figs/mrna_dups.png", height=1200, width=2000, pointsize=32)

par(mfrow=c(2,1), mar=c(5.1, 0.1, 2.1, 2.1))

hist(rho[lower.tri(rho)], breaks=seq(-1, 1, len=301), xlab="Correlation between samples",
     ylab="", yaxt="n", main="")
rug(rho[lower.tri(rho)])
mtext(side=3, adj=1, "mRNAs", cex=1.5, col="slateblue")

hist(rho_adj[lower.tri(rho)], breaks=seq(-1, 1, len=301), xlab="Correlation between samples",
     ylab="", yaxt="n", main="")
rug(rho_adj[lower.tri(rho)])
mtext(side=3, adj=1, "mRNAs\nsex-adjusted", cex=1.5, col="slateblue")

dev.off()
