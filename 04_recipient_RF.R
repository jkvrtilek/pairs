# analysis for 2024 paired recordings
# question 2: do bats sound different depending on recipient?
# for each bat, do a crossed pDFA to assign calls to receiver after controlling for trial (doesn't need to be nested bc we're doing a different one for each bat), for loop through bats
# Julia Vrtilek, May 2025

# load packages
library(ranger)
library(caret)
library(tidyverse)

setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/pairs")

# load data
batcalls <- readRDS("vocal_data_2024-pairs_transformed.RDS") %>% 
  mutate(date = as.Date(date)) %>% 
  mutate(session = case_when(date < "2024-07-22" ~ 1,
                             date > "2024-07-21" & date < "2024-07-29" ~ 2,
                             date > "2024-07-28" ~ 3))

# QUESTION: control for relationship strength?
d.RF <- batcalls %>% 
  select(caller, receiver, duration:session) %>% 
  mutate(caller = as.factor(caller)) %>% 
  mutate(receiver = as.factor(receiver))

#run RF
rf <- ranger(receiver ~., data=d.RF, keep.inbag=T,
             num.trees=8000, mtry=6, min.node.size=1,
             max.depth=0, replace=FALSE, sample.fraction=0.6)
