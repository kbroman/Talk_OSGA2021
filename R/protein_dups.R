prot <- readRDS("../Data/chick_protein.rds")
sex <- substr(rownames(prot), 1, 1)

prot_adj <- prot
for(x in c("M", "F")) {
    ave <- colMeans(prot[sex==x,], na.rm=TRUE)
    prot_adj[sex==x,] <- t(t(prot[sex==x,]) - ave)
}

rho <- cor(t(prot), use="pair")
rho_adj <- cor(t(prot_adj), use="pair")


png("../Figs/protein_dups.png", height=1200, width=2000, pointsize=32)

par(mfrow=c(2,1), mar=c(5.1, 0.1, 2.1, 2.1))

hist(rho[lower.tri(rho)], breaks=seq(-1, 1, len=301), xlab="Correlation between samples",
     ylab="", yaxt="n", main="")
rug(rho[lower.tri(rho)])
mtext(side=3, adj=1, "Proteins", cex=1.5, col="slateblue")

hist(rho_adj[lower.tri(rho)], breaks=seq(-1, 1, len=301), xlab="Correlation between samples",
     ylab="", yaxt="n", main="")
rug(rho_adj[lower.tri(rho)])
mtext(side=3, adj=1, "Proteins\nsex-adjusted", cex=1.5, col="slateblue")

dev.off()
