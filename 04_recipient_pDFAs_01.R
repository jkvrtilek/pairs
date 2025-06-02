# analysis for 2024 paired recordings
# question 2: do bats sound different depending on recipient?
# for each bat, do a crossed pDFA to assign calls to receiver after controlling for trial (doesn't need to be nested bc we're doing a different one for each bat), for loop through bats
# Julia Vrtilek, May 2025

# load packages
library(MASS)
library(tidyverse)

setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/pairs")

# load data
batcalls <- readRDS("vocal_data_2024-pairs_transformed.RDS") %>% 
  mutate(date = as.Date(date)) %>% 
  mutate(session = case_when(date < "2024-07-22" ~ 1,
                             date > "2024-07-21" & date < "2024-07-29" ~ 2,
                             date > "2024-07-28" ~ 3))

# make list of one df per caller
bycaller <- batcalls %>% 
  group_split(caller)

# make pDFA function ----

# Mundry's program requires the following:
# one column for test factor (groups)
# one column for control factor (individuals)
# however many variable columns
# variables should be only numbers, group and individual can be named
# NO missing values
# will take smallest number of calls

#for nested design, data do not have to be balanced
#tests for difference between groups ('testfac')
#groups and subjects ('contrfac') have to numbered consecutively and with integers beginning with 1

pDFA <- function(xdata, test_fac, contr_fac, variables, n.sel = 10, nperm = 1000) {
  if (is.factor(test_fac)==F) {test_fac=as.factor(test_fac)}
  model=paste("lda(test_fac~",variables,", prior=pr_prob, subset=sel.index, data=xdata)",sep="")
  f.table=as.data.frame(table(contr_fac))
  ncf.levels=nrow(f.table)
  ntf.levels=nrow(table(test_fac))
  pr_prob=rep(1/ntf.levels, ntf.levels)#define prior probabilities to be equal for either of two groups
  #get number of cases per subject (level of contrfac):
  #get number of calls to select per subject (subject)
  n.to.sel=min(f.table$Freq)
  #set number of random selections original classification rate should be based on:
  ur.c.val=0
  ur.sel=0
  number=(1:nrow(xdata))
  #get assignment of subjects to groups
  gr=c()
  subj.gr= table(contr_fac,test_fac)
  subject=rownames(subj.gr)
  for (i in 1:ncf.levels){
    for (k in 1:ntf.levels){
      if (subj.gr[i,k]>0) {gr=c(gr,colnames(subj.gr)[k])}
    }
  }
  
  for (k in 1:n.sel){
    #make random selection of same number of cases per subject
    sel=rep(NA,nrow(xdata))#create var. for the random selection to be indicated
    for (i in 1:ncf.levels){
      sel[contr_fac==subject[i]] =sample(c(rep(1, n.to.sel), rep(0,f.table[i,2]-n.to.sel)), f.table[i,2],replace=F)
    }
    sel.index= number[sel==1]
    #do a DFA and store results in 'res':
    res=eval(parse(text=model))
    #get predictions and store them in 'pred':
    pred=predict(res,xdata,prior=pr_prob)$class
    ur.sel=ur.sel+sum((test_fac==pred)[sel==1])
    ur.c.val= ur.c.val+ sum((test_fac==pred)[sel==0])
  }
  #save number of correctly classified calls in variable 'orig.res':
  ur.sel= ur.sel/ n.sel
  ur.c.val= ur.c.val/ n.sel
  
  #set P-value to 1 (since original data should be treated as 1 permutation):
  p.sel=1
  p.c.val=1
  all.corr=matrix(NA, nrow=nperm, ncol=2)
  all.corr[1,1]=ur.sel
  all.corr[1,2]=ur.c.val
  
  if (length(gr)==ncf.levels){
    for (k in 1:(nperm-1)){
      #randomize subjects' assignments to groups:
      r.gr=sample(gr,length(gr), replace=F)
      for (i in 1:length(subject)){
        test_fac[contr_fac==subject[i]]=r.gr[i]
      }
      
      #make random selection or same number of cases per subject
      sel=rep(NA,nrow(xdata))#create var. for the random selection to be indicated
      for (i in 1:ncf.levels){
        sel[contr_fac==subject[i]] =sample(c(rep(1, n.to.sel), rep(0,f.table[i,2]-n.to.sel)), f.table[i,2],replace=F)
      }
      sel.index= number[sel==1]
      #do a DFA and store results in 'res':
      res=eval(parse(text=model))
      #get predictions and store them in 'pred':
      pred=predict(res,xdata,prior=pr_prob)$class
      ran.sel= sum((test_fac==pred)[sel==1])
      ran.c.val= sum((test_fac==pred)[sel==0])
      if (ran.sel>=ur.sel){p.sel = p.sel + 1}
      if (ran.c.val>= ur.c.val){p.c.val= p.c.val + 1}
      all.corr[k+1,1]=ran.sel
      all.corr[k+1,2]=ran.c.val
    }
    what=c("N correctly assigned, original, selected", "P for selected", "N correctly assigned, original, cross-validated", "P for cross-validated", "N groups (levels of test factor)", "N subjects (levels of control factor)", "N cases total", "N cases selected per subject","N selected total", "N permutations","N random selections")
    value=c(ur.sel,p.sel/nperm,ur.c.val,p.c.val/nperm,ntf.levels,ncf.levels,nrow(xdata),n.to.sel,n.to.sel*ncf.levels,nperm,n.sel)
    result=data.frame(what,value)
  }else{
    result="at least one subject is member of two groups; no test done"
  }
  result
  write.table(result,file=paste("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/pairs/pDFA/",xdata$caller[1]),sep="",row.names=F,col.names=T)
  all.corr[,1]#comprises the number correctly classified selected objects for all
  #permutations (with the first value being that for the original data)
  all.corr[,2]#same for the cross-validated objects
  hist(all.corr[,1])#shows the frequency distribution
}

# for loop, one pDFA per caller ----
for (i in 1:7) {
  
  xdata <- bycaller[[i]]
  
  test_fac <- as.factor(as.character(xdata$receiver))
  contr_fac <- as.factor(as.character(xdata$session))
  variables <- paste(colnames(xdata)[5:ncol(xdata)], collapse = "+")
  
  pDFA(xdata, test_fac, contr_fac, variables)
  
}












