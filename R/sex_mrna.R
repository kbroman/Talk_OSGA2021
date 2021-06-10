# check sex in mRNA samples

library(here)
library(parallel)
library(broman)
purple_green <- brocolors("web")[c("purple", "green")]

mrna <- readRDS(here("Data/chick_mrna.rds"))
sex <- substr(rownames(mrna),1,1)

# t-test for sex
nlogp <- unlist(mclapply(1:ncol(mrna),
                         function(index) -log10(t.test(mrna[,index] ~ sex)$p.value),
                         mc.cores=detectCores()))

mrna_sub <- mrna[,order(nlogp, decreasing=TRUE)[1:100]]

genes <- c(ENSMUSG00000086503="Xist",
           ENSMUSG00000066071="Cyp4a12a",
           ENSMUSG00000078799="Sult2a5",
           ENSMUSG00000075551="Cyp3a41a",
           ENSMUSG00000038656="Cyp3a16",
           ENSMUSG00000042589="Cux2")

## plot of individual genes
pdf(here("Figs/sex_mrna_indgenes.pdf"), height=6, width=10, pointsize=14)
par(mar=c(5.1,4.1, 1.6, 1.1), mfrow=c(2,3))
sexalt <- sex
sexalt[sex=="F"] <- "female"
sexalt[sex=="M"] <- "male"
for(i in 1:6) {
    dotplot(sexalt, mrna_sub[,i], bg=purple_green[(sex=="F") + 1],
            main=genes[colnames(mrna_sub)[i]],
            xlab="", ylab="mRNA level")
}
dev.off()


# get average of males and average of females
# take RMS difference from each of these
# plot vs each other, colored by sex

mave <- colMeans(mrna_sub[sex=="M",], na.rm=TRUE)
fave <- colMeans(mrna_sub[sex=="F",], na.rm=TRUE)

mdist <- apply(t(t(mrna_sub) - mave), 1, function(a) sqrt(mean(a^2, na.rm=TRUE)))
fdist <- apply(t(t(mrna_sub) - fave), 1, function(a) sqrt(mean(a^2, na.rm=TRUE)))

pdf(here("Figs/sex_mrna.pdf"), height=6, width=7, pointsize=14)
par(mar=c(5.1,4.1, 1.6, 1.1))
grayplot(fdist, mdist, bg=purple_green[(sex=="F") + 1],
         xlab="distance from female average",
         ylab="distance from male average")
legend("topright", pch=21, pt.bg=rev(purple_green), c("female", "male"), bg="gray90")
dev.off()
