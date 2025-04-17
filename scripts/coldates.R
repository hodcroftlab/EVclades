setwd("~/workspace/fluclades")
library("dplyr")

# these metadata cannot be released to public domain
# gisaid <- read.csv("gisaid.csv")
genbank <- read.csv("data/metadata.tsv", sep="\t")
genbank$year <- substr(genbank$date, 1, 4)
# number of sequences per year
genbank <- genbank %>% group_by(year) %>% 
  summarise(nseq = n()) %>% 
  ungroup() %>% 
  mutate(nseq = ifelse(is.na(nseq), 0, nseq))

## WHO FluNet database 
#fnt <- read.csv("VIW_FNT.csv")
#viw <- read.csv("VIW_FID_EPI.csv")
#viw <- viw[viw$AGEGROUP_CODE=='All',]
#ncases <- sapply(split(viw$ILI_CASES, viw$MMWR_YEAR), function(x) sum(x, na.rm=T))

## number of influenza A detections (all subtypes)
#flua <- sapply(split(fnt$INF_A, fnt$MMWR_YEAR), function(x) sum(x, na.rm=T))
#nsamp <- sapply(split(fnt$SPEC_RECEIVED_NB, fnt$MMWR_YEAR), function(x) sum(x, na.rm=T))
#idx <- match(genbank$year, names(flua))
#genbank$ndetect <- flua[idx]

pdf("results/plots/genbank-nseqs.pdf", width=4.5, height=4.5)

par(mar=c(5,5,1,1))
barplot(genbank$nseq[order(genbank$year)]/100, space=0, 
        border=NA, col='grey50', cex.axis=0.8, las=2,
        ylab="Number of EV sequences (hundreds)")
p <- pretty(1:nrow(genbank))
years <- sort(unique(genbank$year))
axis(side=1, at=p+0.5, labels=years[p+1], cex.axis=0.8)
title(xlab="Year of sample collection")

dev.off()

# evidence of exponential trend
fit <- lm(log(nseq)~year, data=genbank)
summary(fit)  # displays R-squared

