suppressPackageStartupMessages({
  require(ggfree)
  require(phangorn)
  require(xtable)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript plot-trees.R <tree_file> <output_png> <output_pdf>")
} 

tree_file <- args[1]
output_png <- args[2]
output_pdf <- args[3]

# while(FALSE){
  phy <- read.tree(tree_file) #"gb-relabeled-ha.ft2-mle.mid.nwk")
  phy <- ladderize(phy)
  eL <- tree.layout(phy)  # takes a few minutes! 
  eL$nodes$sero <- gsub("^.*_([A-Z0-9]{1,5})_.*$", "\\1", eL$nodes$label)
  idx <- which(!grepl("^.*_([A-Z0-9]{1,5})_.*$", eL$nodes$label))
  eL$nodes$sero[idx] <- NA # nodes
  eL$nodes$sero <- toupper(eL$nodes$sero)
  
  # hphy <- read.tree("gb-relabeled-ha.ft2-mle.mid.nwk")
  # hphy <- ladderize(hphy)
  # hL <- tree.layout(hphy)  # takes a few minutes! 
  # hL$nodes$sero <- gsub("^.*_([Hh][0-9]+)N*[0-9]*_.*$", "\\1", hL$nodes$label)
  # idx <- which(!grepl("^.*_([Hh][0-9]+)N*[0-9]*_.*$", hL$nodes$label))
  # hL$nodes$sero[idx] <- NA
  # 
  # nphy <- read.tree("gb-relabeled-na.ft2-mle.mid.nwk")
  # nphy <- ladderize(nphy)
  # nL <- tree.layout(nphy)
  # nL$nodes$sero <- gsub("^.*_[Hh]*[0-9]*([Nn][0-9]+)_.*$", "\\1", nL$nodes$label)
  # idx <- which(!grepl("^.*_[Hh]*[0-9]*([Nn][0-9]+)_.*$", nL$nodes$label))
  # nL$nodes$sero[idx] <- NA
  # nL$nodes$sero <- toupper(nL$nodes$sero)
  
  # save(hL, nL, file="~/workspace/fluclades/plots/treeplots.RData")
  save(eL, file="results/plots/treeplots.RData")
# }

load('results/plots/treeplots.RData')

# locate edges - this takes a while, so run once and save results!
# while(FALSE) {
  # Extract matched serotype from each label
  matches <- gsub("^.*_([A-Z0-9]{1,5})_.*$", "\\1", eL$nodes$label)
  
  # Remove NAs or things that didn't match the pattern properly
  matches <- matches[!is.na(matches) & matches != eL$nodes$label]
  
  # Count frequency of each
  serotypes <- names(table(matches))
  
  # Create a list of edges for each serotype
  edges <- lapply(1:length(serotypes), function(i) {
    idx <- which(eL$nodes$sero == serotypes[i])
    find.clade(eL, tips=eL$nodes$label[idx])
  })
  
  # h.edges <- lapply(1:18, function(i) {
  #   idx <- which(hL$nodes$sero == paste("H", i, sep=""))
  #   find.clade(hL, tips=hL$nodes$label[idx])
  # })
  # n.edges <- lapply(1:11, function(i) {
  #   idx <- which(nL$nodes$sero == paste("N", i, sep=""))
  #   find.clade(nL, tips=nL$nodes$label[idx])
  # })
  # save(h.edges, n.edges, file="results/plots/edge-index.RData")
  save(edges, file="results/plots/edge-index.RData")
  save(serotypes, file ="results/plots/serotypes.RData")
# }
load('results/plots/edge-index.RData')


# plot functions
rup <- function(x, n=2) {
  sort(x, decreasing=TRUE)[min(n, length(x))]  # runner up
}

colour.clade <- function(obj, idx, col, skip=5, offset=0.03, cex=1, lwd=2) {
  draw.clade(obj, idx, col, lwd=lwd)
  node.idx <- obj$edges$child[idx]
  node.idx <- node.idx[obj$nodes$n.tips[node.idx]==0]
  
  # label clades using tip farthest from origin
  last.tip <- node.idx[which.max(obj$nodes$x[node.idx])]
  max.x <- rup(obj$nodes$x[node.idx], skip) + (offset*max(obj$nodes$x))
  mid.y <- median(obj$nodes$y[node.idx])
  text(x=max.x, y=mid.y, label=obj$nodes$sero[last.tip], 
       col=col, xpd=NA, cex=cex)
}


# prepare plot device
res <- 600
png(output_png, width=5*res, height=5*res, res=res) #width=10*res for HA
# par(mfrow=c(1,2))
par(mfrow=c(1,1), mar=c(5,5,1,1), cex=0.8)

plot(eL, label='n', lwd=1, mar=c(0,1,0,0))
add.scalebar(eL,
             y0 = 100,   # a bit above the bottom
             dy = 100,   # distance to label
             len = 0.01, # len=0.1 for HA
             cex = 0.7) 
pal <- gg.rainbow(length(serotypes))
for (i in 1:length(serotypes)) {
  skip <- ifelse(i %in% c(1,3), 20, 2)
  colour.clade(eL, edges[[i]], col=pal[i], lwd=1, skip=skip, cex=0.8)
}
# plot(nL, label='n', lwd=1, mar=c(0,0,0,1))
# add.scalebar(nL, y0=1000, len=0.1, dy=1000, cex=0.7)
# pal <- gg.rainbow(11)
# for (i in 1:11) {
#   print(i)
#   skip <- ifelse(i %in% c(1, 4, 5, 8), 20, 2)
#   colour.clade(nL, n.edges[[i]], col=pal[i], lwd=1, skip=skip, cex=0.8)
# }
invisible(dev.off())

# get.tips in ggfree has an error - overwrite with this function
get.tips <- function(parent, obj, res = c()) {
  children <- obj$edges$child[obj$edges$parent == parent]
  for (child in children) {
    if (obj$nodes$n.tips[child] == 0) {
      res <- c(res, child)
    } else {
      res <- get.tips(child, obj, res)
    }
  }
  return(res)
}

# what are the nearest neighbours for unclassified sequences?
nearest.sero <- function(cidx, L, skip=0) {
  parent <- L$edges$parent[which(L$edges$child==cidx)]
  children <- setdiff(get.tips(parent, L), cidx)
  tab <- table(L$nodes$sero[children])
  tries <- 0
  while (length(tab) == 0 || tries <= skip) {
    # all descendants of immediate ancestor are unknown
    if (tries > 20) {
      print(paste("Failed to retrieve labelled neighbour for ", cidx))
      break  # failsafe
    }
    parent <- L$edges$parent[which(L$edges$child==parent)]
    children <- setdiff(get.tips(parent, L), cidx)
    tab <- table(L$nodes$sero[children])
    tries <- tries + 1
  }
  names(sort(tab, decreasing=TRUE))[1]
}

# extract accession numbers
eL$nodes$label <- gsub("'", "", eL$nodes$label)
eL$nodes$accn <- sapply(eL$nodes$label, function(x) strsplit(x, "_")[[1]][1])

# nL$nodes$label <- gsub("'", "", nL$nodes$label)
# nL$nodes$accn <- sapply(nL$nodes$label, function(x) strsplit(x, "_")[[1]][1])


# detect misclassified sequences by vertical placement in tree
par(mar=c(5,5,1,1))

bx <- boxplot(split(eL$nodes$y, eL$nodes$sero)) 

miss <- which(eL$nodes$y %in% bx$out)

predicted <- sapply(miss, function(cidx) {
  nearest.sero(cidx, eL, skip=2)
})

toprint <- eL$nodes[miss, c('accn', 'sero')]
toprint$predicted <- predicted

print(xtable(toprint), include.rownames=FALSE)


# boxplot(split(nL$nodes$y, nL$nodes$sero)) -> bx
# miss <- which(nL$nodes$y %in% bx$out)
# predicted <- sapply(miss, function(cidx) {
#   nearest.sero(cidx, nL, skip=2)
# })
# 
# 
# toprint <- nL$nodes[miss, c('accn', 'label', 'sero')]
# toprint$label <- gsub(".+_(A/\\w+/[^/]+/[^_]+).+", "\\1", toprint$label)
# toprint$predicted <- predicted
# print(xtable(toprint), include.rownames=FALSE)


# classify sequences with unspecified subtype
unknown <- which(eL$nodes$n.tips==0 & is.na(eL$nodes$sero))
estimate <- sapply(unknown, function(x) nearest.sero(x, eL, skip=1))
tab <- table(eL$nodes$sero)
comp <- data.frame(subtype=names(tab), known=as.integer(tab))
# tab <- table(estimate)
tab <- table(unlist(estimate)) # remove NULL values
comp$unknown <- as.integer(tab)[match(comp$subtype, names(tab))]

# n.unk <- which(nL$nodes$n.tips==0 & is.na(nL$nodes$sero))
# n.est <- sapply(n.unk, function(x) nearest.sero(x, nL, skip=1))
# tab <- table(nL$nodes$sero)
# n.comp <- data.frame(subtype=names(tab), known=as.integer(tab))
# tab <- table(unlist(n.est))
# n.comp$unknown <- as.integer(tab)[match(n.comp$subtype, names(tab))]


pdf(output_pdf, width=8, height=4)
par(mfrow=c(1,2), mar=c(5,5,1,1), cex=1)
plot(unknown~known, data=comp, las=1, cex.axis=0.8, 
     type='n', xlab="Number of annotated sequences",
     ylab="Number of inferred sequences", bty='n')
text(comp$known, comp$unknown, comp$subtype, cex=0.6)
idx <- c(2,4:8,11,13,15:17)
abline(lm(unknown~known, data=comp[idx,]), lty=2, untf=T)

# plot(unknown~known, data=n.comp[-c(2,3),], las=1, cex.axis=0.8, 
#      type='n', xlab="Number of annotated sequences",
#      ylab="Number of inferred sequences", bty='n')
# text(n.comp$known, n.comp$unknown, n.comp$subtype, cex=0.6)
# abline(lm(unknown~known, data=n.comp[-c(1,4),]), lty=2, untf=T)
invisible(dev.off())
