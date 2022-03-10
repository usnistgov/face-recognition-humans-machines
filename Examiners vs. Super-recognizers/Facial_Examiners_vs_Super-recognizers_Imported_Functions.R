# Functions ----

auc.w<-function(posScores, negScores) {
  # Calculate AUC 
  # posScores = same-id trials
  # negScores = different-id trials
  
  numPos<-length(posScores)
  numNeg<-length(negScores)
  
  wResults<-suppressWarnings(wilcox.test(posScores, negScores))
  grabAUC<-unname(wResults$statistic/(numPos*numNeg))
  
  return(grabAUC)
}#end function auc.w

anova.supplemental<-function(df_effect,df_error,ss_effect,ss_error){
  #compute mean square error and partial eta squared for ANOVAs
  
  ms.error<-ss_error/df_error
  partial.eta.sq<-ss_effect/(ss_effect+ss_error)
  
  return(list(ms.e=ms.error,p.eta.sq=partial.eta.sq))
  
}#end function anova.supplemental


do.pairwise.kappa<-function(dataFrame=analysis.data,groupName){ 
  #dataFrame = raw data
  #groupName = '1to1Face' or 'Super'
  
  #prepare
  dataFrame<-droplevels(dataFrame[dataFrame$Group==groupName,])
  
  raters.matrix<-tapply(dataFrame$Scores,list(dataFrame$Images,dataFrame$SubID),sum,check.names = FALSE)
  raters.matrix<-as.data.frame(raters.matrix)
  
  #get all possible combinations (non repeating)
  pairs.labels<-combn(as.character(unique(dataFrame$SubID)), 2) %>% 
    t() %>%
    as_tibble(.,.name_repair="universal") %>%
    suppressMessages() %>%
    rename(SubID1=...1,SubID2=...2)
  
  #initiate data frame
  kappa.df<-data.frame(kappa.coeff=numeric(nrow(pairs.labels)))
  
  #compute
  for (i in 1:nrow(pairs.labels)) {
    hold.subj1<-raters.matrix[[as.character(pairs.labels[i,1])]]
    hold.subj2<-raters.matrix[[as.character(pairs.labels[i,2])]]
    kappa.df$kappa.coeff[i]<-mean(fleiss.linear.i(cbind(hold.subj1,hold.subj2),weights = "ordinal"))
  }
  
  return(kappa.df)
  
}#end function do.pairwise.kappa


do.pairwise.kappa.boot<-function(dataFrame1,dataFrame2,pairs.labels.1,pairs.labels.2){ 
  #dataFrame1 = rating matrix for group 1
  #dataFrame2 = rating matrix for group 2
  
  ## Examiners

  #initialize matrix
  raters.1.booted<-data.frame(kappa=rep(NA,nrow(pairs.labels.1)))
  
  #get agreement between all possible pairs of participants
    for (s in 1:nrow(pairs.labels.1)) {
      hold.subj1<-dataFrame1[[as.character(pairs.labels.1[s,1])]]
      hold.subj2<-dataFrame1[[as.character(pairs.labels.1[s,2])]]
      raters.1.booted$kappa[s]<-mean(fleiss.linear.i(cbind(hold.subj1,hold.subj2),weights = "ordinal"))
    }
    
  ## Super-recognizers

  #initialize matrix
  raters.2.booted<-data.frame(kappa=rep(NA,nrow(pairs.labels.2)))
  
  #get agreement between all possible pairs of participants
  for (s in 1:nrow(pairs.labels.2)) {
    hold.subj1<-dataFrame2[[as.character(pairs.labels.2[s,1])]]
    hold.subj2<-dataFrame2[[as.character(pairs.labels.2[s,2])]]
    raters.2.booted$kappa[s]<-mean(fleiss.linear.i(cbind(hold.subj1,hold.subj2),weights = "ordinal"))
  }
  
  kappa.diff.out<-mean(raters.1.booted$kappa)-mean(raters.2.booted$kappa)
  
  
  return(kappa.diff.out)
  
  
}#end function do.pairwise.kappa.boot


# Define variables ----

colorblind.blues<-c(RColorBrewer::brewer.pal(n = 7,name="Blues")[2],RColorBrewer::brewer.pal(n = 7,name="Blues")[7])
# Produce two colorblind safe blues