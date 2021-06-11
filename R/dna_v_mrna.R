# study dnva vs mrna for mix-ups

library(here)
library(lineup2)
library(broman)
library(qtl2)

mrna <- readRDS(here("Data/chick_mrna.rds"))
probs <- readRDS(here("Data/chick_genoprobs.rds"))
probs <- probs[,"-X"]
sex_mrna <- substr(rownames(mrna), 1, 1)
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


mrna <- adjust4sex(mrna, sex_mrna)

cache_dir <- here("R/_cache")
if(!dir.exists(cache_dir)) dir.create(cache_dir)

# genome scans
file <- file.path(cache_dir, "maxlod_mrna.rds")
if(file.exists(file)) {
    maxlod_mrna <- readRDS(file)
} else {
    maxlod_mrna <- scan1max(probs, mrna, cores=0)
    saveRDS(maxlod_mrna, file)
}

# genome scans for top 100
top_mrna <- order(maxlod_mrna, decreasing=TRUE)[1:100]
file <- file.path(cache_dir, "peaks_top_mrna.rds")
if(file.exists(file)) {
    peaks_top_mrna <- readRDS(file)
} else {
    lod_mrna <- scan1(probs, mrna[,top_mrna], cores=0)
    peaks_top_mrna <- find_peaks(lod_mrna, map, threshold=min(maxlod_mrna[top_mrna])-0.01)
    peaks_top_mrna$marker <- find_marker(map, peaks_top_mrna$chr, peaks_top_mrna$pos)
    saveRDS(peaks_top_mrna, file)
}


mrna_obs <- mrna[,top_mrna]

file <- file.path(cache_dir, "dist_dna_mrna.rds")
if(file.exists(file)) {
    dist_dna_mrna <- readRDS(file)
} else {
    # get fitted values
    mrna_fitted <- matrix(ncol=100, nrow=nrow(probs[[1]]))
    dimnames(mrna_fitted) <- list(rownames(probs[[1]]), colnames(mrna_obs))
    for(i in 1:nrow(peaks_top_mrna)) {
        pr <- pull_genoprobpos(probs, marker=peaks_top_mrna[i,"marker"])
        fit <- fit1(pr, mrna[,peaks_top_mrna[i,"lodcolumn"]], zerosum=FALSE, se=FALSE)
        mrna_fitted[,i] <- pr %*% fit$coef
    }

    # calculate distances
    dist_dna_mrna <- dist_betw_matrices(mrna_obs, mrna_fitted,
                                        "rmsd", align_cols=FALSE, cores=0)
    saveRDS(dist_dna_mrna, file)
}


# plot histogram of protein/mrna correlations + scatterplots of some strong ones, colored by sex
pdf(here("Figs/hist_lod_mrna.pdf"), height=6, width=10, pointsize=14)
par(mar=c(5.1, 1.1, 1.1, 1.1))
hist(maxlod_mrna, breaks=seq(0, max(maxlod_mrna), len=201),
     xlab="Maximum LOD score, by mRNA",
     ylab="", yaxt="n", main="")
rug(maxlod_mrna[top_mrna])
mtext(side=3, adj=1, "max LOD\nadjusting for sex",
      col="darkslateblue", line=-4, cex=1.2)
dev.off()

# heat map of correlations between samples
png(here("Figs/heatmap_dna_v_mrna.png"), height=1000, width=2600, pointsize=36)
par(mar=c(5.1, 4.1, 1.1, 1.1), las=1)
image(1:ncol(dist_dna_mrna), 1:nrow(dist_dna_mrna), t(dist_dna_mrna), col=rev(revgray()),
      ylab="mRNA sample", xlab="DNA sample")
dev.off()

# histogram of self vs non-self
pdf(here("Figs/self_nonself_dna_mrna.pdf"), height=6, width=10, pointsize=10)
hist_self_nonself(dist_dna_mrna, xlab="Distance")
dev.off()

# plot of best vs self
pdf(here("Figs/self_v_best_dna_mrna.pdf"), height=6, width=8, pointsize=14)
self <- get_self(dist_dna_mrna)
best <- get_best(dist_dna_mrna)
grayplot(self, best, xlab="self-distance", ylab="minimum distance",
         xlim=c(0.4, 3), ylim=c(0.4, 3))
dev.off()

# plot of a few of the samples
pdf(here("Figs/samples_dna_mrna.pdf"), height=6, width=8, pointsize=14)
par(mfcol=c(2,3), mar=c(5.1, 4.1, 2.1, 1.1))
for(i in c("M348", "M410", "M349", "F371", "F352", "M386")) {
    v <- dist_dna_mrna[i,]
    grayplot(v, ylim=c(0.4, 3), xlab="DNA sample", ylab="Distance", main=paste("mRNA sample", i))
    points(which.min(v), min(v), cex=1.2, pch=21, bg="violetred")
    text(which.min(v)-20, min(v), names(v)[which.min(v)], adj=c(1, 0.5))
    points(which(names(v)==i), v[i], cex=1.2, pch=21, bg="orange3")
}
dev.off()
