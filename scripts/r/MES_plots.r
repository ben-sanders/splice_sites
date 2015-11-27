#### SETUP ####

library("ggplot2")
setwd("~/Documents/Ben_project/splice_sites/refseq_sites/")

#### DATA IMPORT ####

refseq_MES_data = read.table("out.refseq_sites.mespy", sep="\t", header=TRUE, na.strings="%s")
random_MES_data = read.table("../random_sites/out.random_sites.mespy", sep="\t", header=TRUE, na.strings="%s")

# sites where MESpy could not predict a splice-site are given scores of 999, which is 
# far higher than the perfect consensus site could get. Remove these, as they are not
# wanted.
refseq_cleaned_MES_data = subset(refseq_MES_data, WT < 900)
random_cleaned_MES_data = subset(random_MES_data, WT < 900)

refseq_MES_uncalled = length(subset(refseq_MES_data, WT > 900))
random_MES_uncalled = length(subset(random_MES_data, WT > 900))

# remove uncleaned data
rm(refseq_MES_data)
rm(random_MES_data)

# add a column for source, and combine the two data frames - can now colour by source, compare both on same plot, etc.
refseq_cleaned_MES_data$src <- "REFSEQ"
random_cleaned_MES_data$src <- "RANDOM"

## NOTE: for testing, use a subset of refseq data matching random data
merged_MES_data <- merge(refseq_cleaned_MES_data, random_cleaned_MES_data, all=TRUE)

# factorise the source column, so it can be used in plotting
merged_MES_data$src <- as.factor(merged_MES_data$src)
# normalise the results between -1 and 1 (or is it 0 and 1?)
merged_MES_data$normSCORE <- (merged_MES_data$WT - min(merged_MES_data$WT)) / (max(merged_MES_data$WT) - min(merged_MES_data$WT))

#### PLOT DATA ####

## full histogram ##

pdf("plots/refseq_full.pdf")
ggplot(refseq_cleaned_MES_data, aes(WT)) + geom_histogram(fill="blue", alpha=0.7, binwidth=1) + ylab("Count") + xlab("MES Score") + ggtitle("Distribution of MES scores for RefSeq splice sites (3' and 5' combined)")
dev.off()

pdf("plots/random_full.pdf")
ggplot(random_cleaned_MES_data, aes(WT)) + geom_histogram(fill="red", alpha=0.7, binwidth=1) + ylab("Count") + xlab("MES Score") + ggtitle("Distribution of MES scores for randomly selected non-splice AG/GT dinucleotides")
dev.off()

# show both random and refseq on the same plot
pdf("plots/combined_full.pdf")
ggplot(merged_MES_data, aes(WT)) + geom_histogram(aes(fill=src), alpha=0.5, binwidth=1, position="identity") + ylab("Count") + xlab("MES Score") + ggtitle("Distribution of MES scores for RefSeq and random sites")
dev.off()

## below zero ##

# want to see the area below zero in more detail
pdf("plots/refseq_negative.pdf")
ggplot(subset(refseq_cleaned_MES_data, WT < 0), aes(WT)) + geom_histogram(fill="blue", alpha=0.7, binwidth=1) + ylab("Count") + xlab("MES Score") + ggtitle("Distribution of scores below 0 for RefSeq splice sites")
dev.off()

## above zero ##

# and might as well do above zero as well
pdf("plots/refseq_positive.pdf")
ggplot(subset(refseq_cleaned_MES_data, WT > 0), aes(WT)) + geom_histogram(fill="blue", alpha=0.7, binwidth=1) + ylab("Count") + xlab("MES Score") + ggtitle("Distribution of scores above 0 for RefSeq splice sites")
dev.off()

#### SUMMARISE DATA ####

# want to produce a summary, and print to a file.
# once random sequences are available, might want to do the sense/spec calculations here

### sensitivity/specificity

# want the range of the scores, to iterate through an calculate sens/spec for ROC curve 
MES_ROC_results = matrix(ncol=7, nrow=length(seq(from = floor(min(merged_MES_data$normSCORE)), to = ceiling(max(merged_MES_data$normSCORE)), by = 0.0001)))
i = 1
for(cutoff in seq(from = floor(min(merged_MES_data$normSCORE)), to = ceiling(max(merged_MES_data$normSCORE)), by = 0.0001))
{
  # divide everything by a large amount, as there are so many values here it was causing integer overflow when calculating MCC
  # dividing everything keeps the proportions the same, doesn't affect results.
  tp = (length(merged_MES_data$normSCORE[merged_MES_data$normSCORE > cutoff & merged_MES_data$src=="REFSEQ"]))/1000
  # false negatives include all refseq sites not called (i.e. those with scores == 990)
  fn = (refseq_MES_uncalled + length(merged_MES_data$normSCORE[merged_MES_data$normSCORE < cutoff & merged_MES_data$src=="REFSEQ"]))/1000
  # true negatives include all random sites not called
  tn = (random_MES_uncalled + length(merged_MES_data$normSCORE[merged_MES_data$normSCORE < cutoff & merged_MES_data$src=="RANDOM"]))/1000
  fp = length(merged_MES_data$normSCORE[merged_MES_data$normSCORE > cutoff & merged_MES_data$src=="RANDOM"])/1000
  
  tpr = tp / (tp + fn)
  tnr = tn / (tn + fp)
  fpr = 1 - tnr
  ppv = tp / (tp + fp)
  npv = tn / (tn + fn)
  accuracy = (tp + tn) / (tp + tn + fp + fn)
  mcc = ((tp * tn) - (fp * fn)) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  
  MES_ROC_results[i,1] <- fpr
  MES_ROC_results[i,2] <- tpr
  MES_ROC_results[i,3] <- cutoff
  MES_ROC_results[i,4] <- accuracy
  MES_ROC_results[i,5] <- mcc
  MES_ROC_results[i,6] <- ppv
  MES_ROC_results[i,7] <- npv
  i <- i +1
}
MES_ROC_results <- data.frame(MES_ROC_results)

pdf("plots/MESpy_ROC.pdf")
ggplot(MES_ROC_results) + geom_line(aes(x=X1, y=X2), size=1) + xlab("False positive rate") +xlim(0,1) + ylab("True positive rate") +ylim(0,1) + ggtitle("MaxEntScan reciever operating characterstic curve")
dev.off()

# sensitivity, specificity, and accuracy vs. cutoff
ggplot(MES_ROC_results, aes(x=X3)) + geom_line(aes(y=X2, colour="TPR"), size=1) + geom_line(aes(y=X1, colour = "FPR"), size=1) + geom_line(aes(y=X4, colour="Accuracy"), size=1) + scale_colour_manual("", breaks=c("TPR", "FPR", "Accuracy"), values=c("TPR"="blue", "FPR"="green", "Accuracy"="red")) + xlab("cutoff") +ylab("") + ggtitle("Accuracy of MaxEntScan")
