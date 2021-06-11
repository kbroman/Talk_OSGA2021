######################################################################
# genotype vs. expression
######################################################################

library(here)

gve_scheme <- function(version=1) {
library(here)

# data for the illustration
load(here("Data/gve_scheme_data.RData"))
z <- t(z[,-10])
# z = similarity matrix
# og = observed genotype
# ig = inferred genotype
# e = expression data

maincolor <- "darkslateblue"
f2color <- broman::brocolors("f2")
sexcolor <- broman::brocolors("sex")
CCcolor <- qtl2::CCcolors
darkgray <- "gray60"
gray <- "gray90"
cex_main <- 1.2
cex_point <- 0.8
cex_xlab <- 0.9

set.seed(35410287)
u <- rnorm(length(og), 0, 0.08)





par(mar=c(2.1,4.1,0.1,4.1), bty="n")

plot(0,0,type="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0,100), ylim=c(107,0),
     xaxs="i", yaxs="i")

xw <- 22
xgap <- c(7, 8)
x1 <- c(xgap[1], xgap[1]+xw)
x2 <- c(x1[2]+xgap[2], x1[2]+xgap[2]+xw)
ygap <- 10
yh <- c(28, 23)
y1 <- c(ygap, ygap+yh[1])
y2 <- c(ygap, ygap+yh[2])

ytgap <- 3
xtgap <- 2

rect(x1[1], y1[1], x1[2], y2[2])
rect(x1[1], y2[2]+0.5, x1[2], y1[2]+0.5)
text(mean(x1), y1[1]-ytgap, "eQTL genotypes", cex=cex_main, col=maincolor, font=2)
text(x1[1]-xtgap, mean(y1), "DNA samples", srt=90)
text(mean(x1), y1[2]+0.5+ytgap, "eQTL")

rect(x2[1], y1[1], x2[2], y2[2])
text(mean(x2), y1[1]-ytgap, "expression traits", cex=cex_main, col=maincolor, font=2)
text(x2[1]-xtgap, mean(y2), "mRNA samples", srt=90)
text(mean(x2), y2[2]+ytgap, "genes with strong eQTL")

if(version==1) return()

recw <- 0.9
rect(x1[1]+diff(x1)*0.2, y1[1], x1[1]+diff(x1)*0.2+recw, y2[2],
     col=CCcolor[1])
rect(x2[1]+diff(x2)*0.2, y1[1], x2[1]+diff(x2)*0.2+recw, y2[2],
     col=CCcolor[1])


arrowgap <- 4.5
arroww <- 5
arrows(x2[2]+arrowgap, mean(y1), x2[2]+arrowgap+arroww, mean(y1), len=0.1, angle=15, col=maincolor, lwd=2)

x3 <- c(99-xw, 99)
rect(x3[1], y1[1], x3[2], y1[2], col=gray) # box with gray background
text(x3[1]-xtgap, mean(y1), "gene expression", srt=90)
text(mean(x3), y1[2]+2.5*ytgap, "eQTL genotype")
x <- c(-1,0,1)/3*diff(range(x3))+mean(x3)
segments(x, rep(y1[2],3), x, rep(y1[2],3)+ytgap/2)
text(x, rep(y1[2],3)+ytgap*1.5, c("AA", "AB", "BB"), cex=cex_xlab, adj=c(0.5,0))
segments(x, rep(y1[1], 3), x, rep(y1[2], 3),col=darkgray, lwd=2, lend=1, ljoin=1) # vertical lines at genotypes
rect(x3[1], y1[1], x3[2], y1[2]) # black border again

og <- og+u-2
og <- og/3*diff(range(x3)) + mean(x3)
e <- mean(y1) - ((e  - mean(range(e, na.rm=TRUE)))/diff(range(e, na.rm=TRUE))*diff(range(y1))/1.1)


points(og, e, bg="lightblue", pch=21, cex=cex_point)

if(version==2) return()

                                      # part C
y3 <- y1+y1[2]+ygap+10
y4 <- c(y3[1], y3[1]+yh[2])

rect(x1[1], y3[1], x1[2], y3[2])
text(mean(x1), y3[1]-ytgap, "predicted expression", cex=cex_main, col=maincolor, font=2)
text(x1[1]-xtgap, mean(y3), "DNA samples", srt=90)
text(mean(x1), y3[2]+ytgap, "eQTL")


rect(x2[1], y4[1], x2[2], y4[2])
text(mean(x2), y4[1]-ytgap, "expression traits", cex=cex_main, col=maincolor, font=2)
text(x2[1]-xtgap, mean(y4), "mRNA samples", srt=90)
text(mean(x2), y4[2]+ytgap, "eQTL")

if(version==3) return()

arrows(x2[2]+arrowgap, mean(y3), x2[2]+arrowgap+arroww, mean(y3), len=0.1, angle=15, col=maincolor, lwd=2)

rech <- 1.3
rect(x1[1], y3[1]+diff(y3)*0.3, x1[2], y3[1]+diff(y3)*0.3+rech, col=CCcolor[1])
rect(x2[1], y4[1]+diff(y4)*0.6, x2[2], y4[1]+diff(y4)*0.6+rech, col=sexcolor[2])


x3 <- c(99-xw, 99)
y5 <- c(y3[1], y3[1]+diff(y3)*1.2)
rect(x3[1], y5[1], x3[2], y5[2])
text(x3[1]-xtgap, mean(y5), "mRNA samples", srt=90)
text(mean(x3), y5[2]+ytgap, "DNA samples")
text(mean(x3), y5[1]-ytgap, "distance matrix", cex=cex_main, col=maincolor, font=2)



                                      # colors
library(RColorBrewer)
blues<-colorRampPalette(c("white","blue"))(256)
zcol <- matrix("", ncol=ncol(z), nrow=nrow(z))
zcol[] <- blues[z[!is.na(z)]*255+1]

xgap <- 0.3
x <- x3[1]+xgap + (0:(ncol(z)-1))*(diff(x3)-xgap*2)/ncol(z)
xn <- x + diff(x[1:2])
ygap <- 0.5
y <- y5[1]+ygap + (nrow(z):1)*(diff(y5)-2*ygap)/nrow(z)
yn <- y + diff(y[1:2])
for(i in seq(along=y))
  for(j in seq(along=x))
    rect(x[j], y[i], xn[j], yn[i], col=zcol[i,j], border="white")

rect(x3[1], y5[1], x3[2], y5[2])
}

for(i in 1:4) {
    file <- here("Figs", paste0("gve_mixup_scheme", ifelse(i==4, "", i), ".pdf"))
    pdf(file, height=6, width=10, pointsize=10)
    gve_scheme(i)
    dev.off()
}
