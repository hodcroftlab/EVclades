args <- commandArgs(trailingOnly = TRUE)

grid <- args[1]
params_set <- args[2]
out_pdf <- args[3]
nsubtrees_pdf <- args[4]

grid <- read.csv(grid)
load("results/plots/serotypes.RData")  # should contain `serotypes` vector
lenS = length(serotypes)
grid$st.diff <- log(grid$nsubtrees / lenS)

pdf(out_pdf, width=10, height=10)
par(mar=c(5,5,2,1), mfrow=c(2,2), cex=1)

plot(grid$minlen, grid$maxlen, pch=21,
     cex=sqrt(abs(grid$nsubtrees - lenS))/5,
     bg=ifelse(grid$nsubtrees > lenS, 'white', 'black'),
     xlab="Minimum divergence (v)", ylab="Maximum patristic distance (y)")
title(main="Number of subtrees", adj=0, font.main=1)
text(x=-0.02, y=2.2, label="A", cex=2, xpd=NA)

plot(jitter(grid$maxlen), grid$mean.n.types, pch=19, cex=0.5,
     ylab="Mean number per subtree", xlab="Maximum patristic distance (y)")
title(main="Mean # of different subgenotype labels", adj=0, font.main=1)
text(x=-0.35, y=9.1, label="B", cex=2, xpd=NA)

plot(grid$minlen, grid$maxlen, cex=2.5 * sqrt(1 - grid$nlabels), 
     pch=22, bg='black', col=NA,
     xlab="Minimum divergence (v)", ylab="Maximum patristic distance (y)")
points(grid$minlen[grid$nlabels < 1], grid$maxlen[grid$nlabels < 1], 
       pch=22, cex=2.5, lwd=0.5, col='grey20')
title(main="Fraction of dropped labels", adj=0, font.main=1)
text(x=-0.02, y=2.2, label="C", cex=2, xpd=NA)

res <- read.csv(params_set)
res <- res[!grepl("UNK", res$serotype), ]  # remove unknowns if needed

# Heatmap
par(mar=c(5,5,2,2))
plot(NA, xlim=c(0.75, lenS + 0.25), ylim=c(-0.25, max(res$subtree) + 1), 
     xaxt='n', yaxt='n', xlab="Subgenotype labels", ylab="Subtree", bty='n')
title(main="Distribution of labels", adj=0, font.main=1)
text(x=-3.2, y=max(res$subtree)+2, label="D", cex=2, xpd=NA)

subtree_list <- sort(unique(res$subtree))
for (i in seq_along(subtree_list)) {
  rows <- res[res$subtree == subtree_list[i], ]
  for (sero in rows$serotype) {
    j <- which(serotypes == sero)
    count <- rows[rows$serotype == sero, "count"]
    xx <- 1 - (count / sum(rows$count))
    rect(j-1, i-1, j, i, col=rgb(xx, xx, xx), border=NA)
    if (xx >= 0.5) text(j-0.5, i-0.5, adj=0.5, label=count, cex=0.5)
  }
}
axis(side=1, at=1:lenS - 0.5, label=serotypes, cex.axis=0.8, las=2)
axis(side=2, at=0:(length(subtree_list)-1) + 0.5, label=subtree_list, las=2, cex.axis=0.8)
axis(side=4, at=0:(length(subtree_list)-1) + 0.5, 
     label=sapply(split(res$count, res$subtree), sum), cex.axis=0.6, las=2, lwd=0, line=-0.5)

for (i in 0:lenS) {
  abline(v=i, col='grey80')
}
for (i in 0:length(subtree_list)) {
  abline(h=i, col='grey80')
}

invisible(dev.off())

# Optional: ridgeplot of nsubtrees
library(ggfree)  # or ggridges
pdf(nsubtrees_pdf, width=5, height=5)
par(mar=c(5,5,1,1))
ridgeplot(split(log10(grid$nsubtrees), grid$minlen), 
          fill=gg.rainbow(10, alpha=0.5),
          ylab="Minimum divergence", xlab="Log(number of subtrees)", xaxt='n',
          cex.axis=0.8, step=0.3)
axis(side=1, at=0:4, labels=10^(0:4), cex.axis=0.8)
abline(v=log10(18), lty=2)
rug(log10(grid$nsubtrees))
invisible(dev.off())
