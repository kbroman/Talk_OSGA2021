library(fst)
library(broman)
library(qtl2)

x <- read.fst("../Data/attieDO_chrXint.fst")
y <- read.fst("../Data/attieDO_chrYint.fst")
stopifnot(ncol(x) == ncol(y), all(colnames(x) == colnames(y)))
keep <- colnames(x)[colMeans(is.na(x)) < 0.2 & colMeans(is.na(y)) < 0.02]
# I don't understand what's going on with these
keep <- keep %wnin% paste0("DO", c(171, 141, 165))
x <- x[,keep]
y <- y[,keep]


do <- readRDS("../Data/attieDO_v1.rds")
do <- do[n_missing(do, "ind", "prop") < 0.1,]
wh <- which(colnames(x) %in% ind_ids(do))
x <- x[,c(1,wh)]
y <- y[,c(1,wh)]
is_female <- do$is_female[colnames(x)[-1]]

map <- read.csv("~/Data/MUGAarrays/UWisc/gm_uwisc_v1.csv")
mapX <- map[!is.na(map$chr) & map$chr == "X",]
mapY <- map[!is.na(map$chr) & map$chr == "Y",]
xmar <- mapX$marker %win% x[,1]
ymar <- mapY$marker %win% y[,1]

xm <- colMeans(x[x[,1] %in% xmar, -1], na.rm=TRUE)
ym <- colMeans(y[y[,1] %in% ymar, -1], na.rm=TRUE)

pdf("../Figs/xydosage.pdf", height=6, width=7, pointsize=14)
par(mar=c(6.1,5.1,1.6,1.1))
purple_green <- brocolors("web")[c("purple", "green")]
grayplot(xm, ym, bg=purple_green[is_female + 1],
         xlab="ave X chr intensity",
         ylab="ave Y chr intensity")
legend("topright", pch=21, pt.bg=purple_green, c("male", "female"), bg="gray90")
dev.off()
