# check sex in mRNA samples

library(here)
library(parallel)
library(broman)
purple_green <- brocolors("web")[c("purple", "green")]

protein <- readRDS(here("Data/chick_protein.rds"))
sex <- substr(rownames(protein),1,1)

# t-test for sex
nlogp <- unlist(mclapply(1:ncol(protein),
                         function(index) -log10(t.test(protein[,index] ~ sex)$p.value),
                         mc.cores=detectCores()))

protein_sub <- protein[,order(nlogp, decreasing=TRUE)[1:100]]

genes <- c(ENSMUSP00000138475="Slc22a27",
           ENSMUSP00000025833="Papss2",
           ENSMUSP00000031822="Abcg2",
           ENSMUSP00000073667="Mup20",
           ENSMUSP00000096407="Sult2a3",
           ENSMUSP00000028014="Fmo4")


## plot of individual genes
pdf(here("Figs/sex_protein_indgenes.pdf"), height=6, width=10, pointsize=14)
par(mar=c(5.1,4.1, 1.6, 1.1), mfrow=c(2,3))
sexalt <- sex
sexalt[sex=="F"] <- "female"
sexalt[sex=="M"] <- "male"
for(i in 1:6) {
    dotplot(sexalt, protein_sub[,i], bg=purple_green[(sex=="F") + 1],
            main=genes[colnames(protein_sub)[i]],
            xlab="", ylab="protein level")
}
dev.off()


# get average of males and average of females
# take RMS difference from each of these
# plot vs each other, colored by sex

mave <- colMeans(protein_sub[sex=="M",], na.rm=TRUE)
fave <- colMeans(protein_sub[sex=="F",], na.rm=TRUE)

mdist <- apply(t(t(protein_sub) - mave), 1, function(a) sqrt(mean(a^2, na.rm=TRUE)))
fdist <- apply(t(t(protein_sub) - fave), 1, function(a) sqrt(mean(a^2, na.rm=TRUE)))

pdf(here("Figs/sex_protein.pdf"), height=6, width=7, pointsize=14)
par(mar=c(5.1,4.1, 1.6, 1.1))
grayplot(fdist, mdist, bg=purple_green[(sex=="F") + 1],
         xlab="distance from female average",
         ylab="distance from male average")
ontop <- (fdist < 1 & sex=="M") | (fdist > 1 & sex=="F")
points(fdist[ontop], mdist[ontop], pch=21, bg=purple_green[(sex[ontop]=="F") + 1])
legend("topright", pch=21, pt.bg=rev(purple_green), c("female", "male"), bg="gray90")
dev.off()
