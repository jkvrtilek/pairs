# analysis for 2024 paired recordings
# question 1: do different bats call different amounts? does calling rate depend on recipient?
# use Bayesian social relations model, strand package
# actor, receiver, and interaction automatically included
# response variable = call count
# predictor = Haley's grooming OR food-sharing? maybe just the data from before my data?
# results: effects of actor's tendency to call, receiver's tendency to receive, correlation at dyadic level between caller calling and receiver calling (do the bats I call to more call to me more, controlling for overall calling and receiving rate), correlation at group level (if I call more to the group, does the group call more to me?) and THEN finally will tell us whether bats are calling more depending on predictor variable
# how to make use of added power from 3 replicate matrices? for now, take average, but keep thinking about this
# Julia Vrtilek
# Adapted from Haley Gmutza, May 2025

# load packages
library(STRAND)
library(tidyverse)
library(igraph)

setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/pairs")

# load and wrangle vocal data ----
batcalls <- readRDS("vocal_data_2024-pairs.RDS") %>% 
  group_by(caller, receiver) %>% 
  summarize(n.calls = n()) %>% 
  ungroup() %>% 
  add_row(caller = "quark", receiver = "yikes", n.calls = 0) %>% # quark never called. this row makes names appear alphabetically
  arrange(caller)

netcalls <- graph_from_data_frame(batcalls)
mcalls <- as_adjacency_matrix(netcalls, attr= 'n.calls', sparse=F)
nets = list(call = mcalls)

bats.used <- rownames(mcalls)

# load and wrangle social data ----

# already trimmed bat donations list
donations <- read.csv('OSU_2024_social_data.csv')
donations$Actor <- tolower(donations$Actor)
donations$Receiver <- tolower(donations$Receiver)

bat.donations <- donations %>% 
  filter(Actor %in% bats.used) %>% 
  filter(Receiver %in% bats.used)

# make foodsharing matrix
rates <-
  bat.donations %>%
  filter(Behavior == "Mouthlicking") %>% 
  mutate(edge= paste(Actor, Receiver, sep="_")) %>%
  group_by(edge) %>%
  summarize(rate= sum(rate, na.rm=T)) %>%
  filter(rate>=0) %>%
  separate(edge, into=c('Actor', 'Receiver'))

net <- graph_from_data_frame(rates)
m <- as_adjacency_matrix(net, attr= 'rate', sparse=F)

# # make grooming matrix
# ratesg <-
#   bat.donations %>%
#   filter(Behavior == "Grooming") %>% 
#   mutate(edge= paste(Actor, Receiver, sep="_")) %>%
#   group_by(edge) %>%
#   summarize(rate= sum(rate, na.rm=T)) %>%
#   filter(rate>=0) %>%
#   separate(edge, into=c('Actor', 'Receiver'))
# 
# netg <- graph_from_data_frame(ratesg)
# mg <- as_adjacency_matrix(netg, attr= 'rate', sparse=F)

# individual variable - age
chars <- read.csv('campus_bat_chars.csv', stringsAsFactors = F)
chars$Bat.name <- tolower(chars$Bat.name)

batchars <- chars %>% 
  filter(Bat.name %in% bats.used) %>% 
  select(Age)

rownames(batchars) <- bats.used

# dyadic variables - kinship, foodsharing
distmat <- read.delim('Desmodus_DistanceMatrix.txt',row.names = 1)
distmat <- data.matrix(distmat)

dyad = list(#Kinship = distmat,
            Lick = m)

# Make the STRAND data structure ----
dat = make_strand_data(
  outcome = nets,
  individual_covariates = batchars,
  dyadic_covariates = dyad,
  outcome_mode = "poisson",
  link_mode = "log",
  check_standardization = F) # I don't see why we would standardize age??? but maybe we should

# run the model ----
fit =
  fit_social_relations_model(
    data=dat,
    focal_regression = ~ Age,
    target_regression = ~ Age,
    dyad_regression = ~ Lick,
    mode="mcmc",
    stan_mcmc_parameters = list(
      chains = 3,
      iter_warmup = 1500,
      iter_sampling = 1500))

# summarize results
res = summarize_strand_results(fit)









#Wrangle model results for plotting
plotdat = vector("list",length(res$summary_list))

for(k in 1:length(res$summary_list)){
  plotdat[[k]] = data.frame(res$summary_list[[k]])
  plotdat[[k]]$SubModel = names(res$summary_list)[k]
  colnames(plotdat[[k]]) = c("Variable", "Median", "LI", "HI", "Mean","SD", "SubModel")
  for(j in 2:6)
    plotdat[[k]][,j] = as.numeric(plotdat[[k]][,j])
}

df = do.call(rbind, plotdat)

colnames(df) = c("Variable", "Median", "LI", "HI", "Mean","SD", "SubModel")

df <- df[,1:7]

#PLOT USING GGPLOT - to Gerry's specifications lol
#p1 - dyadic effects
p1 <-
  df %>% 
  filter(Variable == "dyadic effects coeffs, Kinship" | Variable == "dyadic effects sd") %>% 
  mutate(Median = as.numeric(Median), low = as.numeric(LI), high = as.numeric(HI)) %>% 
  ggplot()+
  geom_point(aes(x = Variable, y = Median), size = 2) +
  geom_errorbar(aes(x = Variable, ymin = low, ymax = high, width = 0), size = 1)+
  #ylim(-0,2.75)+
  geom_hline(aes(yintercept = 0), alpha = 0.25)+
  coord_flip()+
  theme_bw()+
  xlab("covariates")+
  ylab("standard deviation")+
  ggtitle("dyadic effects: kinship")
p1

# plot focal, target, and random effects
p2 <- 
  df %>% 
  filter(Variable == "focal effects sd" | Variable == "target effects sd" 
         | Variable == "dyadic effects sd") %>% 
  mutate(Variable = case_when(
    Variable == "focal effects sd" ~ "actor effects",
    Variable == "target effects sd" ~ "receiver effects",
    Variable == "dyadic effects sd" ~ "dyadic effects",
  )) %>%
  mutate(Variable = as_factor(fct_relevel(Variable, c("dyadic effects", "receiver effects", "actor effects")))) %>% 
  mutate(Median = as.numeric(Median), low = as.numeric(LI), high = as.numeric(HI)) %>% 
  ggplot()+
  geom_point(aes(x = Variable, y = Median), size = 2) +
  geom_errorbar(aes(x = Variable, ymin = low, ymax = high, width = 0), size = 1)+
  #ylim(-0,2.75)+
  geom_hline(aes(yintercept = 0), alpha = 0.25)+
  coord_flip()+
  theme_bw()+
  xlab("random effects")+
  ylab("standard deviation")+
  ggtitle("random effects on food sharing")
p2

#individual 
p <- 
  df %>% 
  filter(Variable == "focal effects coeffs (out-degree), Age"
         | Variable == "target effects coeffs (in-degree), Age") %>% 
  mutate(Variable = case_when(
    Variable == "focal effects coeffs (out-degree), Age" ~ "actor effect, Age",
    Variable == "target effects coeffs (in-degree), Age" ~ "receiver effects, Age",
  )) %>%
  mutate(Median = as.numeric(Median), low = as.numeric(LI), high = as.numeric(HI)) %>% 
  ggplot()+
  geom_point(aes(x = Variable, y = Median), size = 2) +
  geom_errorbar(aes(x = Variable, ymin = low, ymax = high, width = 0), size = 1)+
  #ylim(-1.75,1.75)+
  geom_hline(aes(yintercept = 0), alpha = 0.25)+
  coord_flip()+
  theme_bw()+
  #ylab("posterior estimate")+
  theme(axis.title.x = element_blank())+
  xlab("effect")
p
