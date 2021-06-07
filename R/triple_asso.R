# illustration of associations

text_color <- "darkslateblue"
text_cex <- 2
arrow_color <- "violetred"
arrow_lwd <- 3

pdf("../Figs/triple_asso.pdf", height=6, width=10, pointsize=14)

par(bty="n", mar=rep(0.1,4))
plot(0,0,type="n", xlim=c(0,100), ylim=c(0,100), xaxs="i", yaxs="i",
     xlab="", ylab="", xaxt="n", yaxt="n")


pos <- list(c(50, 75), c(25, 25), c(75, 25))
lab <- c("genotypes", "clinical\nphenotypes", "mRNA levels")

for(i in seq_along(pos)) {
    text(pos[[i]][1], pos[[i]][2], lab[i], col=text_color, cex=text_cex)
}

arrows(44, 65, 30, 40, len=0.1, lwd=arrow_lwd, col=arrow_color)
arrows(56, 65, 68, 40, len=0.1, lwd=arrow_lwd, col=arrow_color)
arrows(42, 25, 58, 25, len=0.1, lwd=arrow_lwd, col=arrow_color, code=3)


dev.off()
