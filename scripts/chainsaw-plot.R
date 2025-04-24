require(xtable)

args <- commandArgs(trailingOnly = TRUE)

print("Command-line arguments:")
print(args)

nsubtrees <- args[1]
label <- args[2]
out_pdf <- args[3]
out_table <- args[4]
keep_na <- args[5]

rescale <- function(x, from, to) {
  # map vector `x` to range of `to`, given it comes from range of `from`
  diff(range(to)) * (x-min(from)) / diff(range(from)) + min(to)
}

chainsaw <- read.csv(nsubtrees)
chainsaw <- chainsaw[chainsaw$cutoff > 0.007,]  # drop cutoffs that yield too many subtrees (37 trees)

# manual runs of chainsaw.py
pdf(out_pdf, width=5, height=4)
{par(mar=c(5,5,1,5))

plot(chainsaw$cutoff, chainsaw$nsubtrees, type='n', 
     xlab="Maximum internal branch length", ylab=NA,
     bty='n', yaxt='n', cex.axis=0.8)
axis(side=2, col='royalblue3', col.axis='royalblue3', las=2, cex.axis=0.8)
mtext(side=2, text="Number of subtrees", line=3, col='royalblue3')

# rescale 
y <- rescale(chainsaw$normalized, chainsaw$normalized, chainsaw$nsubtrees)
#polygon(c(min(chainsaw$cutoff), chainsaw$cutoff, max(chainsaw$cutoff)), 
#        c(0, y, 0), col=rgb(0.9,0,0,0.2), border=NA)
points(chainsaw$cutoff, y, col='firebrick3', pch=19, cex=0.33)
segments(x0=chainsaw$cutoff[which.max(y)], x1=max(chainsaw$cutoff)*1.1,
         y0=max(y), col='firebrick3', lty=2)
text(x=mean(chainsaw$cutoff[which.max(y):nrow(chainsaw)]), 
     y=max(y)*1.02, xpd=NA,
     label=round(max(chainsaw$normalized), 4), 
     cex=0.5, col='firebrick3')

p <- pretty(chainsaw$normalized)
axis(side=4, at=rescale(p, chainsaw$normalized, chainsaw$nsubtrees), 
     label=p, cex.axis=0.7, las=2,
     col='firebrick3', col.axis='firebrick3')
u <- par('usr')
text(x=u[2]+0.06, y=mean(u[3:4]), 
     label="Normalized mutual information", col='firebrick3', 
     srt=-90, xpd=NA)

lines(chainsaw$cutoff, chainsaw$nsubtrees, type='s', col='royalblue')
#points(chainsaw$cutoff, chainsaw$nsubtrees, pch=19, cex=0.33, 
#       col='royalblue')
#abline(h=18, col=rgb(0,0,0,0.2))

# there are 12 serotypes in EV-D68
load('results/plots/serotypes.RData')

segments(x0=0, x1=0.017, y0=length(serotypes), col='royalblue', lty=2)
text(x=0.0175, y=length(serotypes), label=paste0(length(serotypes)), col="royalblue", cex=0.5)
}
print(paste("saved to", out_pdf))
dev.off()


######  NA  ######
# chainsaw <- read.csv("results/chainsaw-nsubtrees-na.csv")
# 
# pdf("~/papers/fluclades/chainsaw-na.pdf", width=5, height=4)
# 
# par(mar=c(5,5,1,5))
# plot(chainsaw$cutoff, chainsaw$nsubtrees, type='n', 
#      xlab="Maximum internal branch length", ylab=NA,
#      bty='n', yaxt='n', cex.axis=0.8)
# axis(side=2, col='royalblue3', col.axis='royalblue3', las=2, cex.axis=0.8)
# mtext(side=2, text="Number of subtrees", line=3, col='royalblue3')
# 
# # rescale 
# foo <- chainsaw$normalized[chainsaw$normalized > 0.1]
# y <- rescale(chainsaw$normalized, foo, chainsaw$nsubtrees)
# #polygon(c(min(chainsaw$cutoff), chainsaw$cutoff, max(chainsaw$cutoff)), 
# #        c(0, y, 0), col=rgb(0.9,0,0,0.2), border=NA)
# points(chainsaw$cutoff, y, pch=19, col='firebrick3', cex=0.33)
# segments(x0=chainsaw$cutoff[which.max(y)], x1=max(chainsaw$cutoff)*1.1,
#          y0=max(y), col='firebrick3', lty=2)
# text(x=mean(chainsaw$cutoff[which.max(y):nrow(chainsaw)]), 
#      y=max(y)*1.02, xpd=NA,
#      label=round(max(chainsaw$normalized), 4), 
#      cex=0.5, col='firebrick3')
# 
# p <- pretty(foo)
# axis(side=4, at=rescale(p, foo, chainsaw$nsubtrees), 
#      label=p, cex.axis=0.7, las=2,
#      col='firebrick3', col.axis='firebrick3')
# u <- par('usr')
# text(x=u[2]+0.09, y=mean(u[3:4]), 
#      label="Normalized mutual information", col='firebrick3', 
#      srt=-90, xpd=NA)
# 
# lines(chainsaw$cutoff, chainsaw$nsubtrees, type='s', col='royalblue')
# segments(x0=0, x1=0.41, y0=11, col='royalblue', lty=2)
# text(x=0.44, y=11, label="11", col="royalblue", cex=0.5)
# #points(chainsaw$cutoff, chainsaw$nsubtrees, pch=19, cex=0.5, col='royalblue')
# 
# #abline(h=11, lty=2)
# dev.off()


################# Generate table plots #################
# doFull <- TRUE  # switch between full and protein-based

# examine distribution of serotype labels among subtrees
# if (doFull) {
labels <- read.csv(label) # define a variable to get the different cutoffs
pat <- "^.*_([A-Z0-9]{1,5})_.*$" # full 
# } else {
#   labels <- read.csv("flu-results/chainsaw-NA-0.41.labels.csv")
#   pat <- ".+_H[0-9]+(N[0-9]+)_.+"  # NA  
# }

labels$serotype <- gsub(pat, "\\1", labels$tip.label)
if (keep_na=="False") {
  labels$serotype[!grepl(pat, labels$tip.label)] <- NA
  labels <- na.omit(labels) # remove NAs
  
} else {
  labels$serotype[!grepl(pat, labels$tip.label)] <- "unlabeled"
}

tab <- table(labels$subtree, labels$serotype)
#xtable(tab)  # used to embed table into LaTeX document

# tab <- matrix(sample(0:50, nrow(tab) * ncol(tab), replace = TRUE),
#               nrow = nrow(tab),
#               ncol = ncol(tab),
#               dimnames = dimnames(tab))


# generate a plot
x <- tab / apply(tab, 1, sum)
# if (doFull) {
xval <- 1:length(colnames(x))  # not sure if it would be better to keep original names
translation <- rbind(1:length(colnames(x)), colnames(x))
# } else {
#   xval <- as.integer(gsub("N([0-9]+)", "\\1", colnames(x))) 
# }
xsum <- apply(x, 1, function(row) sum(xval*row))
io <- order(xsum)
jo <- order(xval)  # column index

#hc <- hclust(dist(tab))
# if (doFull) {
  pdf(out_table, width=4.5, height=4.5)  
# } else {
#   pdf("~/papers/fluclades/chainsaw-NA-table.pdf", width=4.5, height=4.5) 
# }
{
par(mar=c(5,5,1,2))

shim <- 0.4 #0.65 or 0.4
plot(NA, xlim=c(shim, ncol(x)-shim), ylim=c(shim, nrow(x)-shim), 
     xaxt='n', yaxt='n',
     xlab="Serotype labels", ylab="Subtree", bty='n')
for (i in 1:nrow(x)) {
  for (j in 1:ncol(x)) {
    xx <- 1-x[io[i], jo[j]]
    if (xx < 1) xx <- min(0.9, xx) 
    rect(j-1, i-1, j, i, col=rgb(xx, xx, xx), border=NA)
    count <- tab[io[i], jo[j]]
    if (count>0 && xx>0.1) text(j-0.5, i-0.5, adj=0.5, label=count, cex=0.5)
  }
}
for (i in 0:ncol(x)-1) {
  abline(v=i, col='grey80')
  abline(h=i, col='grey80')
}
# if (doFull) {
axis(side=1, at=1:ncol(x)-0.5, label=sort(colnames(x)), cex.axis=0.8, las=2)  
# } else {
#   axis(side=1, at=1:ncol(x)-0.5, label=paste("N", 1:11, sep=""), 
#        cex.axis=0.8, las=2)  
# }

axis(side=2, at=1:nrow(x)-0.5, label=paste("s", io-1, sep=""), 
     las=2, cex.axis=0.8)
axis(side=4, at=1:nrow(x)-0.5, label=apply(tab[io,], 1, sum),
     cex.axis=0.6, las=2, lwd=0, line=-0.9)
}
print(out_table)
dev.off()



# ############# other segments ###############
# 
# others <- read.csv("results/chainsaw-nsubtrees-others.csv")
# others <- others[order(others$gene, others$cutoff), ]
# 
# marks <- list(
#   PB2=c(5, 0.046),
#   PB1=c(4, 0.04),
#   PA=c(5, 0.043),
#   NP=c(4, 0.055),
#   M1M2=c(6, 0.045),
#   NS1NS2=c(5, 0.15)
# )
# 
# pdf("~/papers/fluclades/chainsaw-others.pdf", width=11, height=7.3)
# pal <- gg.rainbow(n=6, l=50)
# par(mfrow=c(2,3), mar=c(5,5,1,1), cex=1)
# i <- 1
# for (g in c('PB2', 'PB1', 'PA', 'NP', 'M1M2', 'NS1NS2')) {
#   rows <- others[others$gene==g, ]
#   plot(NA, xlim=c(0.01, 0.2), ylim=c(2, 25), bty='n', 
#        xlab="Maximum internal branch length",
#        ylab="Number of subtrees", )
#   lines(rows$cutoff, rows$nsubtrees, type='s', col=pal[i])
#   title(main=g, font.main=1, adj=0.9, line=-2, cex.main=1.2)
#   points(rows$cutoff, rows$nsubtrees, pch=19, cex=0.5, col=pal[i])
#   points(marks[[g]][2], marks[[g]][1]+0.7, pch=25, bg='black', cex=0.9)
#   i <- i+1
# }
# dev.off()


