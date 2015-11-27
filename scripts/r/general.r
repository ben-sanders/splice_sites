# so far, this seems to be the easiest way to plot the different ROC curves on the same plot.
# but will have to manually add every set of results to both the data frame and plot.
test <- data.frame(MES_FPR=MES_ROC_results$X1, MES_TPR=MES_ROC_results$X2, PWM_FPR=PWM_ROC_results$X1, PWM_TPR=PWM_ROC_results$X2, GS_FPR=GS_ROC_results$X1, GS_TPR=GS_ROC_results$X2)
ggplot(test) + geom_step(aes(x=MES_FPR, y=MES_TPR), colour="red", size=1, alpha=0.7) + geom_step(aes(x=PWM_FPR, y=PWM_TPR), colour="green", size=1, alpha=0.7) + geom_step(aes(x=GS_FPR, y=GS_TPR), colour="blue", size=1, alpha=0.7)