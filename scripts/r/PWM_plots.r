#### SETUP ####

library("ggplot2")
setwd("~/Documents/Ben_project/splice_sites/refseq_sites/")

#### DATA IMPORT ####

refseq_PWM_data = read.table("out.refseq_sites.pwm", sep="\t", header=TRUE, na.strings="%s")
random_PWM_data = read.table("../random_sites/out.random_sites.pwm", sep="\t", header=TRUE, na.strings="%s")

# the PWM should give a score for everything, so this isn't strictly necessary...
refseq_cleaned_PWM_data = subset(refseq_PWM_data, WT > 0)
random_cleaned_PWM_data = subset(random_PWM_data, WT > 0)

refseq_PWM_uncalled = length(subset(refseq_PWM_data, WT == 0))
random_PWM_uncalled = length(subset(random_PWM_data, WT == 0))

# remove uncleaned data
rm(refseq_PWM_data)
rm(random_PWM_data)

# add a column for source, and combine the two data fraPWM - can now colour by source, compare both on same plot, etc.
refseq_cleaned_PWM_data$src <- "REFSEQ"
random_cleaned_PWM_data$src <- "RANDOM"

## NOTE: for testing, use a subset of refseq data matching random data
merged_PWM_data <- merge(refseq_cleaned_PWM_data, random_cleaned_PWM_data, all=TRUE)

# factorise the source column, so it can be used in plotting
merged_PWM_data$src <- as.factor(merged_PWM_data$src)
merged_PWM_data$normWT <- merged_PWM_data$WT / max(merged_PWM_data$WT)
#### PLOT DATA ####

## full histogram ##

# show both random and refseq on the same plot
pdf("plots/combined_full.pdf")
ggplot(merged_PWM_data, aes(normWT)) + geom_histogram(aes(fill=src), alpha=0.5, position="identity") + xlim(0, 0.05) + ylab("Count") + xlab("PWM Score") + ggtitle("Distribution of PWM scores for RefSeq and random sites")
dev.off()

#### SUMMARISE DATA ####

# want to produce a summary, and print to a file.
# once random sequences are available, might want to do the sense/spec calculations here

### sensitivity/specificity

# want the range of the scores, to iterate through an calculate sens/spec for ROC curve 
PWM_ROC_results = matrix(ncol=7, nrow=length(seq(from = floor(min(merged_PWM_data$normWT)), to = ceiling(max(merged_PWM_data$normWT)), by = 0.0001)))
i = 1
for(cutoff in seq(from = floor(min(merged_PWM_data$normWT)), to = ceiling(max(merged_PWM_data$normWT)), by = 0.0001))
{
  # divide everything by a large amount, as there are so many values here it was causing integer overflow when calculating MCC
  # dividing everything keeps the proportions the same, doesn't affect results.
  tp = (length(merged_PWM_data$WT[merged_PWM_data$normWT > cutoff & merged_PWM_data$src=="REFSEQ"]))/1000
  # false negatives include all refseq sites not called (i.e. those with scores == 990)
  fn = (length(merged_PWM_data$WT[merged_PWM_data$normWT < cutoff & merged_PWM_data$src=="REFSEQ"]))/1000
  # true negatives include all random sites not called
  tn = (length(merged_PWM_data$WT[merged_PWM_data$normWT < cutoff & merged_PWM_data$src=="RANDOM"]))/1000
  fp = length(merged_PWM_data$WT[merged_PWM_data$normWT > cutoff & merged_PWM_data$src=="RANDOM"])/1000
  
  tpr = tp / (tp + fn)
  tnr = tn / (tn + fp)
  fpr = 1 - tnr
  ppv = tp / (tp + fp)
  npv = tn / (tn + fn)
  accuracy = (tp + tn) / (tp + tn + fp + fn)
  mcc = ((tp * tn) - (fp * fn)) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  
  PWM_ROC_results[i,1] <- fpr
  PWM_ROC_results[i,2] <- tpr
  PWM_ROC_results[i,3] <- cutoff
  PWM_ROC_results[i,4] <- accuracy
  PWM_ROC_results[i,5] <- mcc
  PWM_ROC_results[i,6] <- ppv
  PWM_ROC_results[i,7] <- npv
  i <- i + 1
}
PWM_ROC_results <- data.frame(PWM_ROC_results)

pdf("plots/PWMpy_ROC.pdf")
ggplot(PWM_ROC_results) + geom_line(aes(x=X1, y=X2), size=1) + xlab("False positive rate") +xlim(0,1) + ylab("True positive rate") +ylim(0,1) + ggtitle("PWM model reciever operating characterstic curve")
dev.off()

# sensitivity, specificity, and accuracy vs. cutoff
ggplot(PWM_ROC_results, aes(x=X3)) + geom_line(aes(y=X2, colour="TPR"), size=1) + geom_line(aes(y=X1, colour = "FPR"), size=1) + geom_line(aes(y=X4, colour="Accuracy"), size=1) + geom_line(aes(y=X5, colour="MCC"), size=1) + xlim(0, 0.01) + scale_colour_manual("", breaks=c("TPR", "FPR", "Accuracy", "MCC"), values=c("TPR"="blue", "FPR"="green", "Accuracy"="red", "MCC"="purple")) + xlab("cutoff") +ylab("") + ggtitle("Accuracy of Position Weight Matrices")
