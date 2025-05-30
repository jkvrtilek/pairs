# analysis for 2024 paired recordings
# question 2: do bats sound different depending on recipient?
# for each bat, do a DFA to assign calls to receiver, for loop through bats
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

# make results df
results <- setNames(data.frame(matrix(ncol = 12, nrow = 7)),
         c("caller","correct.cases","all.cases","accuracy",
           "range.assign.low","range.assign.high","mean.assign","median.assign",
           "range.call.low","range.call.high","mean.call","median.call"))

for (i in 1:7) {
  results$caller[i] <- bycaller[[i]]$caller[1]
  
  dfa <- lda(receiver ~ 
               duration+    
               meanfreq+    
               sd+          
               freq.median+ 
               freq.Q25+ 
               freq.Q75+    
               freq.IQR+    
               time.median+
               time.Q25+   
               time.Q75+    
               time.IQR+    
               skew+        
               kurt+
               sp.ent+      
               time.ent+  
               entropy+     
               sfm+         
               meandom+     
               mindom+     
               maxdom+      
               dfrange+    
               modindx+     
               startdom+    
               enddom+      
               dfslope+   
               meanpeakf+   
               peakf+
               maxslope+
               minslope+
               abs_minslope+
               pos_slopes+
               neg_slopes+
               turns+
               meanslope+
               segments,
             CV=T, 
             data=bycaller[[i]])
  
  # get classification matrix
  cm <- table(bycaller[[i]]$receiver, dfa$class)
  
  # get overall correct classification rate (accuracy)
  # this is the best accuracy estimate
  results$correct.cases[i] <- sum(diag(cm))
  results$all.cases[i] <- sum(cm)
  results$accuracy[i] <- results$correct.cases[i]/results$all.cases[i]
  
  # see median, mean, and range of correct assignments/all assignments to each bat
  results$range.assign.low[i] <- range(diag(cm)/colSums(cm), na.rm=T)[1]
  results$range.assign.high[i] <- range(diag(cm)/colSums(cm), na.rm=T)[2]
  results$mean.assign[i] <- mean(diag(cm)/colSums(cm), na.rm=T)
  results$median.assign[i] <- median(diag(cm)/colSums(cm), na.rm=T)
  
  # see median, mean and range of correct assignments/all calls from each bat
  results$range.call.low[i] <- range(diag(cm)/rowSums(cm), na.rm=T)[1]
  results$range.call.high[i] <- range(diag(cm)/rowSums(cm), na.rm=T)[2]
  results$mean.call[i] <- mean(diag(cm)/rowSums(cm), na.rm=T)
  results$median.call[i] <- median(diag(cm)/rowSums(cm), na.rm=T)
}

