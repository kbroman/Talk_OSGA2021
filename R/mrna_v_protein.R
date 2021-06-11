# study mrna vs protein for mix-ups

library(here)
library(lineup2)
library(broman)

mrna <- readRDS(here("Data/chick_mrna.rds"))
prot <- readRDS(here("Data/chick_protein.rds"))
sex_mrna <- substr(rownames(mrna), 1, 1)
sex_prot <- substr(rownames(prot), 1, 1)

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
prot <- adjust4sex(prot, sex_prot)

cache_dir <- here("R/_cache")
if(!dir.exists(cache_dir)) dir.create(cache_dir)

# corr between prot and mrna
file <- file.path(cache_dir, "cor_prot_mrna.rds")
if(file.exists(file)) {
    cor_prot_mrna <- readRDS(file)
} else {
    cor_prot_mrna <- corr_betw_matrices(prot, mrna, "bestright", cores=0)
    saveRDS(cor_prot_mrna, file)
}

prot_select <- order(cor_prot_mrna[,1], decreasing=TRUE)[1:100]
file <- file.path(cache_dir, "sim_prot_mrna.rds")
if(file.exists(file)) {
    sim_prot_mrna <- readRDS(file)
} else {
    sim_prot_mrna <- corr_betw_matrices(t(prot[,prot_select]),
                                        t(mrna[, cor_prot_mrna[prot_select,2]]),
                                        "all", align_rows=FALSE, cores=0)
    saveRDS(sim_prot_mrna, file)
}


# plot histogram of protein/mrna correlations + scatterplots of some strong ones, colored by sex

pdf(here("Figs/hist_cor_mrna_prot.pdf"), height=6, width=10, pointsize=14)
layout(rbind(1, 2:4))
par(mar=c(5.1, 1.1, 1.1, 1.1))
hist(cor_prot_mrna[,1], breaks=seq(0, 1, len=201),
     xlab="Maximum protein/mRNA correlation, by protein",
     ylab="", yaxt="n", main="")
mtext(side=3, adj=1, "max cor(mRNA, protein)\nadjusting for sex",
      col="darkslateblue", line=-4, cex=1.2)

o <- order(cor_prot_mrna[,1], decreasing=TRUE)[1:3]
purple_green <- brocolors("web")[c("purple", "green")]
proteins <- c("Crym", "Scd1", "Pm20d1")
genes <- c("Crym", "Scd1", "Pm20d1")
par(mar=c(5.1, 4.1, 2.1, 1.1))
for(i in seq_along(o)) {
    grayplot(prot[,o[i]], mrna[rownames(prot),cor_prot_mrna[o[i],3]],
             xlab="Protein level", ylab="mRNA level", main=proteins[i],
             bg=purple_green[(sex_prot=="F")+1])
}
dev.off()


# heat map of correlations between samples
png(here("Figs/heatmap_mrna_v_protein.png"), height=1000, width=2600, pointsize=36)
par(mar=c(5.1, 4.1, 1.1, 1.1), las=1)
image(1:ncol(sim_prot_mrna), 1:nrow(sim_prot_mrna), t(sim_prot_mrna), col=revgray(),
      ylab="protein sample", xlab="mRNA sample")
dev.off()

# histogram of self vs non-self
pdf(here("Figs/self_nonself_mrna_prot.pdf"), height=6, width=10, pointsize=10)
hist_self_nonself(sim_prot_mrna, xlab="Similarity", breaks=seq(-0.6, 1, len=201))
dev.off()

# plot of best vs self
pdf(here("Figs/self_v_best_mrna_prot.pdf"), height=6, width=8, pointsize=14)
self <- get_self(sim_prot_mrna)
best <- get_best(sim_prot_mrna, get_min=FALSE)
grayplot(self, best, xlab="self-similarity", ylab="maximum similarity")
dev.off()

# plot of a few of the samples
pdf(here("Figs/samples_mrna_prot.pdf"), height=6, width=8, pointsize=14)
par(mfcol=c(2,3), mar=c(5.1, 4.1, 2.1, 1.1))
for(i in c("M348", "M410", "M349", "F371", "F352", "M386")) {
    v <- sim_prot_mrna[i,]
    grayplot(v, ylim=c(-0.6, 0.85), xlab="mRNA sample", ylab="Similarity", main=paste("protein sample", i))
    points(which.max(v), max(v), cex=1.2, pch=21, bg="violetred")
    text(which.max(v)-20, max(v), names(v)[which.max(v)], adj=c(1, 0.5))
    points(which(names(v)==i), v[i], cex=1.2, pch=21, bg="orange3")
}
dev.off()
