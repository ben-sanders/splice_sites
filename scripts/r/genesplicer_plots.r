#### SETUP ####

library("ggplot2")
setwd("~/Documents/Ben_project/splice_sites/refseq_sites/")

#### DATA IMPORT ####

refseq_GS_data = read.table("out.genesplicer.refseq", sep="\t", header=TRUE, na.strings="%s")
random_GS_data = read.table("../random_sites/out.genesplicer.random", sep="\t", header=TRUE, na.strings="%s")

# add a column for source, and combine the two data frames - can now colour by source, compare both on same plot, etc.
refseq_GS_data$src <- "REFSEQ"
random_GS_data$src <- "RANDOM"

# sites where MESpy could not predict a splice-site are given scores of 999, which is 
# far higher than the perfect consensus site could get. Remove these, as they are not
# wanted.
refseq_cleaned_GS_data = subset(refseq_GS_data, SCORE < 990)
random_cleaned_GS_data = subset(random_GS_data, SCORE < 990)

# separate out RefSeq called and uncalled
refseq_called = refseq_GS_data[refseq_GS_data$SCORE < 990,]
refseq_uncalled = refseq_GS_data[refseq_GS_data$SCORE > 990,]
# add labels
refseq_called$CALL <- "called"
refseq_uncalled$CALL <- "uncalled"
# calculate percentages
refseq_pc_called <- length(refseq_cleaned_GS_data$SCORE) / length(refseq_GS_data$SCORE) * 100
refseq_pc_uncalled <- 100 - refseq_pc_called

# separate out random called and uncalled
random_called = random_GS_data[random_GS_data$SCORE < 990,]
random_uncalled = random_GS_data[random_GS_data$SCORE > 990,]
# add labels
random_called$CALL <- "called"
random_uncalled$CALL <- "uncalled"
# calculate percentages
random_pc_called <- length(random_cleaned_GS_data$SCORE) / length(random_GS_data$SCORE) * 100
random_pc_uncalled <- 100 - random_pc_called

# merge labelled data - separate for pie charts
refseq <- merge(refseq_called, refseq_uncalled, all=TRUE)
random <- merge(random_called, random_uncalled, all=TRUE)

# merge everything for histogram
# use cleaned data, as only want the called sites
merged_GS_data <- merge(refseq_cleaned_GS_data, random_cleaned_GS_data, all=TRUE)

# this should normalise scores to between 0 and 1
merged_GS_data$normSCORE <- (merged_GS_data$SCORE - min(merged_GS_data$SCORE)) / (max(merged_GS_data$SCORE) - min(merged_GS_data$SCORE))

  # remove uncleaned data
rm(refseq_GS_data)
rm(random_GS_data)
rm(refseq_cleaned_GS_data)
rm(random_cleaned_GS_data)

#### PLOT DATA ####

# theme to remove all the stuff (labels, ticks, etc.) that's left around the pie chart but makes no sense
# after switching from bar using coord_polar().
blank_theme <- theme_grey() + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x=element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank()
  )


# called vs. uncalled pie charts
ggplot(refseq, aes(x=factor(1), fill=CALL)) + geom_bar(width=1) + coord_polar(theta="y") + blank_theme + ggtitle("Proportion of RefSeq sites called by GeneSplicer")
ggplot(random, aes(x=factor(1), fill=CALL)) + geom_bar(width=1) + coord_polar(theta="y") + blank_theme + ggtitle("Proportion of random sites called by GeneSplicer")

# merged histogram (of called sites)
ggplot(merged_GS_data, aes(normSCORE)) + geom_histogram(aes(fill=src), alpha=0.5, binwidth=1, position="identity") + ylab("Count") + xlab("Score") + ggtitle("Distribution of GeneSplicer scores for RefSeq and random sites")

#### SUMMARISE + ROC CURVE ####

# total counts
# percentages
# is ROC possible when so many aren't even called?
GS_ROC_results = matrix(ncol=7, nrow=length(seq(from = floor(min(merged_GS_data$normSCORE)), to = ceiling(max(merged_GS_data$normSCORE)), by = 0.01)))
# i is for indexing into results
i = 1
for(cutoff in seq(from = floor(min(merged_GS_data$normSCORE)), to = ceiling(max(merged_GS_data$normSCORE)), by=0.01))
{
  # called sites above the cutoff
  tp = length(merged_GS_data$normSCORE[merged_GS_data$normSCORE > cutoff & merged_GS_data$src == "REFSEQ"]) / 1000
  # any uncalled in RefSeq data are false negatives, as well as any below the cutoff
  fn = length(merged_GS_data$SCORE[merged_GS_data$normSCORE < cutoff & merged_GS_data$src == "REFSEQ"]) / 1000 # (length(refseq_uncalled$SCORE) + length(refseq_called$SCORE[refseq_called$SCORE < cutoff])) / 1000
  # true negatives are the random sites scored below cutoff + random sites not called
  tn = length(merged_GS_data$SCORE[merged_GS_data$normSCORE < cutoff & merged_GS_data$src == "RANDOM"]) / 1000 # (length(random_uncalled$SCORE) + length(random_called$SCORE[random_called$SCORE < cutoff])) / 1000
  # false positives are any called in random sites
  fp = length(merged_GS_data$SCORE[merged_GS_data$normSCORE > cutoff & merged_GS_data$src == "RANDOM"]) /1000
  
  tpr = tp / (tp + fn) # sensitivity
  tnr = tn / (tn + fp) # specificity
  fpr = fp / (fp + tn) # 1 - specificity
  ppv = tp / (tp + fp)
  npv = tn / (tn + fn)
  accuracy = (tp + tn) / (tp + tn + fp + fn)
  mcc = ((tp * tn) - (fp * fn)) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  
  GS_ROC_results[i,1] <- fpr
  GS_ROC_results[i,2] <- tpr
  GS_ROC_results[i,3] <- cutoff
  GS_ROC_results[i,4] <- accuracy
  GS_ROC_results[i,5] <- mcc
  GS_ROC_results[i,6] <- ppv
  GS_ROC_results[i,7] <- npv
  i <- i +1
}
GS_ROC_results <- data.frame(GS_ROC_results)

ggplot(GS_ROC_results) + geom_line(aes(x=X1, y=X2), size=1) + xlab("False positive rate") + xlim(0,1) + ylab("True positive rate") + ylim(0,1) + ggtitle("GeneSplicer reciever operating characterstic curve")

ggplot(GS_ROC_results, aes(x=X3)) + geom_line(aes(y=X2, colour="TPR"), size=1) + geom_line(aes(y=X1, colour = "FPR"), size=1) + geom_line(aes(y=X4, colour="Accuracy"), size=1) + scale_colour_manual("", breaks=c("TPR", "FPR", "Accuracy"), values=c("TPR"="blue", "FPR"="green", "Accuracy"="red")) + xlab("cutoff") +ylab("") + ggtitle("Accuracy of GeneSplicer")
