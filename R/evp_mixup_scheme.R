######################################################################
# genotype vs. expression
######################################################################

evp_scheme <- function(version=1) {

# data for the illustration
# z = similarity matrix
# pe_dat = protein and expression data
    library(here)
    dat <- readRDS(here("Data/evp_mixup_scheme.rds"))
    z <- dat$sim - min(dat$sim)
    z <- z/diff(range(z))
    pe_dat <- dat$pe_dat

maincolor <- "darkslateblue"
f2color <- broman::brocolors("f2")
sexcolor <- broman::brocolors("sex")
CCcolor <- qtl2::CCcolors
darkgray <- "gray60"
gray <- "gray90"

par(mar=c(2.1,4.1,0.1,4.1), bty="n", cex=1)
    cex_main <- 1.2
    cex_points <- 0.8

plot(0,0,type="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0,100), ylim=c(107,0),
     xaxs="i", yaxs="i")

xw <- 26
xgap <- c(7, 8)
x1 <- c(xgap[1], xgap[1]+xw)
x2 <- c(x1[2]+xgap[2], x1[2]+xgap[2]+xw*0.6)
ygap <- 10
yh <- c(28, 23)
y1 <- c(ygap, ygap+yh[1])
y2 <- c(ygap, ygap+yh[2])

ytgap <- 3
xtgap <- 2

rect(x1[1], y1[1], x1[2], y2[2])
rect(x1[1], y2[2]+0.5, x1[2], y1[2]+0.5)
text(mean(x1), y1[1]-ytgap, "gene expression data", cex=cex_main, col=maincolor, font=2)
text(x1[1]-xtgap, mean(y1), "mRNA samples", srt=90)
text(mean(x1), y1[2]+0.5+ytgap, "genes")

rect(x2[1], y1[1], x2[2], y2[2])
text(mean(x2), y1[1]-ytgap, "proteomics data", cex=cex_main, col=maincolor, font=2)
text(x2[1]-xtgap, mean(y2), "protein samples", srt=90)
text(mean(x2), y2[2]+ytgap, "proteins")

if(version==1) return()

recw <- 1.2
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
text(mean(x3), y1[2]+ytgap, "protein level")
x <- (pe_dat[,1]-mean(range(pe_dat[,1])))/diff(range(pe_dat[,1]))/1.05*diff(range(x3))+mean(x3)
y <- (mean(range(pe_dat[,2]))-pe_dat[,2])/diff(range(pe_dat[,2]))/1.05*diff(range(y1))+mean(y1)
rect(x3[1], y1[1], x3[2], y1[2]) # black border again
points(x,y, pch=21, bg="lightblue", cex=cex_points)

if(version==2) return()

                                      # part C
y3 <- y1+y1[2]+ygap+10
y4 <- c(y3[1], y3[1]+yh[2])
xw2 <- 0.4

rect(x1[1], y3[1], x1[2], y3[2])
text(mean(x1), y3[1]-ytgap, "gene expression data", cex=cex_main, col=maincolor, font=2)
text(x1[1]-xtgap, mean(y3), "mRNA samples", srt=90)
text(mean(x1), y3[2]+ytgap, "genes")


rect(x2[1], y4[1], x2[2], y4[2])
text(mean(x2), y4[1]-ytgap, "proteomics data", cex=cex_main, col=maincolor, font=2)
text(x2[1]-xtgap, mean(y4), "protein samples", srt=90)
text(mean(x2), y4[2]+ytgap, "proteins")

rect(x1[1], y3[1], x1[1]+diff(x2)*xw2, y3[2], col="gray90")
rect(x2[1], y4[1], x2[1]+diff(x2)*xw2, y4[2], col="gray90")


if(version==3) return()

arrows(x2[2]+arrowgap, mean(y3), x2[2]+arrowgap+arroww, mean(y3), len=0.1, angle=15, col=maincolor, lwd=2)

rech <- 1.3
rect(x1[1], y3[1]+diff(y3)*0.3, x1[1]+diff(x2)*xw2, y3[1]+diff(y3)*0.3+rech, col=CCcolor[1])
rect(x2[1], y4[1]+diff(y4)*0.6, x2[1]+diff(x2)*xw2, y4[1]+diff(y4)*0.6+rech, col=sexcolor[2])


x3 <- c(99-xw, 99)
y5 <- c(y3[1], y3[1]+diff(y3)*1.2)
rect(x3[1], y5[1], x3[2], y5[2])
text(x3[1]-xtgap, mean(y5), "protein samples", srt=90)
text(mean(x3), y5[2]+ytgap, "mRNA samples")
text(mean(x3), y5[1]-ytgap, "similarity matrix", cex=cex_main, col=maincolor, font=2)



                                      # colors
library(RColorBrewer)
blues<-colorRampPalette(c("white","blue"))(256)
zcol <- matrix("", ncol=ncol(z), nrow=nrow(z))
zcol[] <- blues[z*255+1]

xgap <- 0.3
x <- x3[1]+xgap + (0:(ncol(z)-1))*(diff(x3)-xgap*2)/ncol(z)
xn <- x + diff(x[1:2])
ygap <- 0.5
y <- y5[1]+ygap + (nrow(z):1)*(diff(y5)-2*ygap)/nrow(z)
yn <- y + diff(y[1:2])

for(i in seq(along=y)) {
    for(j in seq(along=x)) {
        rect(x[j], y[i], xn[j], yn[i], col=zcol[i,j], border="white")
  }
}

rect(x3[1], y5[1], x3[2], y5[2])
}

for(i in 1:4) {
    file <- here("Figs", paste0("evp_mixup_scheme", ifelse(i==4, "", i), ".pdf"))
    pdf(file, height=6, width=10, pointsize=10)
    evp_scheme(i)
    dev.off()
}
