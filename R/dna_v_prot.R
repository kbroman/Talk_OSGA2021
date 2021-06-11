# study dnva vs proteins for mix-ups

library(here)
library(lineup2)
library(broman)
library(qtl2)

prot <- readRDS(here("Data/chick_protein.rds"))
probs <- readRDS(here("Data/chick_genoprobs.rds"))
probs <- probs[,"-X"]
sex_prot <- substr(rownames(prot), 1, 1)
map <- readRDS(here("Data/chick_map.rds"))

adjust4sex <-
    function(x, sex) {
        stopifnot(nrow(x) == length(sex))

        usex <- unique(sex)
        for(u in usex) {
            x[sex==u,] <- scale(x[sex==u,], scale=FALSE)
        }

        x
    }


prot <- adjust4sex(prot, sex_prot)

cache_dir <- here("R/_cache")
if(!dir.exists(cache_dir)) dir.create(cache_dir)

# genome scans
file <- file.path(cache_dir, "maxlod_prot.rds")
if(file.exists(file)) {
    maxlod_prot <- readRDS(file)
} else {
    maxlod_prot <- scan1max(probs, prot, cores=0)
    saveRDS(maxlod_prot, file)
}

# genome scans for top 100
top_prot <- order(maxlod_prot, decreasing=TRUE)[1:100]
file <- file.path(cache_dir, "peaks_top_prot.rds")
if(file.exists(file)) {
    peaks_top_prot <- readRDS(file)
} else {
    lod_prot <- scan1(probs, prot[,top_prot], cores=0)
    peaks_top_prot <- find_peaks(lod_prot, map, threshold=min(maxlod_prot[top_prot])-0.01)
    peaks_top_prot$marker <- find_marker(map, peaks_top_prot$chr, peaks_top_prot$pos)
    saveRDS(peaks_top_prot, file)
}


prot_obs <- prot[,top_prot]

file <- file.path(cache_dir, "dist_dna_prot.rds")
if(file.exists(file)) {
    dist_dna_prot <- readRDS(file)
} else {
    # get fitted values
    prot_fitted <- matrix(ncol=100, nrow=nrow(probs[[1]]))
    dimnames(prot_fitted) <- list(rownames(probs[[1]]), colnames(prot_obs))
    for(i in 1:nrow(peaks_top_prot)) {
        pr <- pull_genoprobpos(probs, marker=peaks_top_prot[i,"marker"])
        fit <- fit1(pr, prot[,peaks_top_prot[i,"lodcolumn"]], zerosum=FALSE, se=FALSE)
        prot_fitted[,i] <- pr %*% fit$coef
    }

    # calculate distances
    dist_dna_prot <- dist_betw_matrices(prot_obs, prot_fitted,
                                        "rmsd", align_cols=FALSE, cores=0)
    saveRDS(dist_dna_prot, file)
}


# plot histogram of protein/prot correlations + scatterplots of some strong ones, colored by sex
pdf(here("Figs/hist_lod_prot.pdf"), height=6, width=10, pointsize=14)
par(mar=c(5.1, 1.1, 1.1, 1.1))
hist(maxlod_prot, breaks=seq(0, max(maxlod_prot), len=201),
     xlab="Maximum LOD score, by prot",
     ylab="", yaxt="n", main="")
rug(maxlod_prot[top_prot])
mtext(side=3, adj=1, "max LOD\nadjusting for sex",
      col="darkslateblue", line=-4, cex=1.2)
dev.off()

# heat map of correlations between samples
png(here("Figs/heatmap_dna_v_prot.png"), height=900, width=2600, pointsize=36)
par(mar=c(5.1, 4.1, 1.1, 1.1), las=1)
image(1:ncol(dist_dna_prot), 1:nrow(dist_dna_prot), t(dist_dna_prot), col=rev(revgray()),
      ylab="protein sample", xlab="DNA sample")
dev.off()

# histogram of self vs non-self
pdf(here("Figs/self_nonself_dna_prot.pdf"), height=6, width=10, pointsize=10)
hist_self_nonself(dist_dna_prot, xlab="Distance")
dev.off()

# plot of best vs self
pdf(here("Figs/self_v_best_dna_prot.pdf"), height=6, width=8, pointsize=14)
self <- get_self(dist_dna_prot)
best <- get_best(dist_dna_prot)
grayplot(self, best, xlab="self-distance", ylab="minimum distance",
         xlim=c(0.4, 1.26), ylim=c(0.4, 0.9))
dev.off()

# plot of a few of the samples
pdf(here("Figs/samples_dna_prot.pdf"), height=6, width=8, pointsize=14)
par(mfcol=c(2,3), mar=c(5.1, 4.1, 2.1, 1.1))
for(i in c("M348", "M410", "M349", "F371", "F352", "M386")) {
    v <- dist_dna_prot[i,]
    grayplot(v, ylim=c(0.4, 1.7), xlab="DNA sample", ylab="Distance", main=paste("protein sample", i))
    points(which.min(v), min(v), cex=1.2, pch=21, bg="violetred")
    text(which.min(v)-20, min(v), names(v)[which.min(v)], adj=c(1, 0.5))
    points(which(names(v)==i), v[i], cex=1.2, pch=21, bg="orange3")
}
dev.off()

# extra plot: best vs self + best vs 2nd best
pdf(here("Figs/best_v_2ndbest_dna_prot.pdf"), height=6, width=11, pointsize=14)
par(mar=c(5.1, 4.1, 1.1, 1.1), pty="s")
self <- get_self(dist_dna_prot)
best <- get_best(dist_dna_prot)
secbest <- get_2ndbest(dist_dna_prot)
par(mfrow=c(1,2))
wh <- !is.na(best) & (secbest-best < 0.15)
grayplot(self, best, xlab="self-distance", ylab="minimum distance",
         xlim=c(0.4, 1.26), ylim=c(0.4, 1.26))
abline(0,1, lwd=2, lty=2)
points(self, best, pch=21, bg="lightblue")
points(self[wh], best[wh], pch=21, bg="orange2")
grayplot(secbest, best, xlab="2nd-best distance", ylab="minimum distance",
         xlim=c(0.4, 1.26), ylim=c(0.4, 1.26))
points(secbest[wh], best[wh], pch=21, bg="orange2")
abline(0,1, lwd=2, lty=2)
dev.off()
