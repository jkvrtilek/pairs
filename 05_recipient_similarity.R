# DFA containing all calls
# get average distance of ALL of B’s calls to all of A’s calls
# compared to average distance of B’s calls to A to all of A’s calls
# need the distances for EVERY CALL, not just group centroids

# load packages
library(MASS)
library(tidyverse)

setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/pairs")

# load data
batcalls <- readRDS("vocal_data_2024-pairs_transformed.RDS") %>% 
  dplyr::select(!c(date, sound.files))

dfa <- lda(caller ~ 
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
           CV=F, 
           data=batcalls)

# get classification matrix
cm <- table(batcalls$caller, dfa$class)

# get overall correct classification rate (accuracy)
# this is the best accuracy estimate
correct.cases <- sum(diag(cm))
all.cases <- sum(cm)
accuracy <- correct.cases/all.cases


# possible functions:
# mahalanobis
# D2.dist
# mahalanobis_distance()
# maha_dist()
