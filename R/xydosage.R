library(fst)
library(broman)
library(qtl2)
library(here)

x <- read.fst(here("Data/attieDO_chrXint.fst"))
y <- read.fst(here("Data/attieDO_chrYint.fst"))
stopifnot(ncol(x) == ncol(y), all(colnames(x) == colnames(y)))

do <- readRDS(here("Data/attieDO_v0.rds"))
do <- do[n_missing(do, "ind", "prop") < 0.2,]
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



pdf(here("Figs/xydosage.pdf"), height=6, width=7, pointsize=14)
par(mar=c(6.1,5.1,1.6,1.1))
purple_green <- brocolors("web")[c("purple", "green")]
grayplot(xm, ym, bg=purple_green[is_female + 1],
         xlab="ave X chr intensity",
         ylab="ave Y chr intensity")
prob <- (is_female & ym > 0.2) | (!is_female & ym < 0.2)
points(xm[prob], ym[prob], pch=21, bg=purple_green[is_female[prob] + 1])
legend("topright", pch=21, pt.bg=purple_green, c("male", "female"), bg="gray90")
dev.off()
