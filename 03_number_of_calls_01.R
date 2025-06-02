# analysis for 2024 paired recordings
# question 1: do different bats call different amounts? does calling rate depend on recipient?
# use Bayesian social relations model, strand package
# actor, receiver, and interaction automatically included
# response variable = call count
# predictor = age, kinship, foodsharing, grooming, combined rate; then foodsharing, grooming, combined rate all w/kinship
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
  add_row(caller = "quark", # quark never called; adding this row makes names appear alphabetically
          receiver = "yikes", n.calls = 0) %>% 
  arrange(caller)

netcalls <- graph_from_data_frame(batcalls)
mcalls <- as_adjacency_matrix(netcalls, attr= 'n.calls', sparse=F)
nets = list(call = mcalls)

bats.used <- rownames(mcalls)
bats.used


# load and wrangle social data ----

# function to make scale() function return a vector not a matrix
scale2 <- function(x){as.vector(scale(x, scale = FALSE))}

# already trimmed bat donations list
raw <- read.csv('OSU_2024_social_data.csv')

raw$Date <- as.Date(raw$Date, format = "%m/%d/%Y")

donations <- raw %>% 
  filter(Date < "2024-09-01")

# how many hours of video data were scored?
hours <- length(unique(donations$Date))
seconds <- hours*3600
# matches Haley's personal communication

donations$Actor <- tolower(donations$Actor)
donations$Receiver <- tolower(donations$Receiver)

bat.donations <- donations %>% 
  filter(Actor %in% bats.used) %>% 
  filter(Receiver %in% bats.used) %>% 
  mutate(rate2 = rate/seconds)

# make foodsharing matrix
ratesf <-
  bat.donations %>%
  filter(Behavior == "Mouthlicking") %>% 
  mutate(edge= paste(Actor, Receiver, sep="_")) %>%
  group_by(edge) %>%
  summarize(rate2= sum(rate2, na.rm=T)) %>%
  filter(rate2>=0) %>%
  separate(edge, into=c('Actor', 'Receiver')) %>% 
  mutate(across(.cols=rate2, .fns = scale2))

netf <- graph_from_data_frame(ratesf)
mf <- as_adjacency_matrix(netf, attr= 'rate2', sparse=F)

# make grooming matrix
ratesg <-
  bat.donations %>%
  filter(Behavior == "Grooming") %>%
  mutate(edge= paste(Actor, Receiver, sep="_")) %>%
  group_by(edge) %>%
  summarize(rate2= sum(rate2, na.rm=T)) %>%
  filter(rate2>=0) %>%
  separate(edge, into=c('Actor', 'Receiver')) %>% 
  mutate(across(.cols=rate2, .fns = scale2))

netg <- graph_from_data_frame(ratesg)
mg <- as_adjacency_matrix(netg, attr= 'rate2', sparse=F)

# make combined "affiliation rate"
# NOTE: don't know possible grooming/foodsharing seconds, but should be the same for all dyads?

ratesa <-
  bat.donations %>% 
  filter(Behavior != "Aggression") %>% 
  mutate(edge= paste(Actor, Receiver, sep="_")) %>%
  group_by(edge) %>%
  summarize(rate2= sum(rate2, na.rm=T)) %>%
  filter(rate2>=0) %>%
  separate(edge, into=c('Actor', 'Receiver')) %>% 
  mutate(across(.cols=rate2, .fns = scale2))

neta <- graph_from_data_frame(ratesa)
ma <- as_adjacency_matrix(neta, attr= 'rate2', sparse=F)


# make the STRAND data structure ----
# individual variable - age
chars <- read.csv('campus_bat_chars.csv', stringsAsFactors = F)
chars$Bat.name <- tolower(chars$Bat.name)

batchars <- chars %>% 
  filter(Bat.name %in% bats.used) %>% 
  select(Age) %>% 
  mutate(across(.cols=Age, .fns = scale2))

rownames(batchars) <- bats.used

# dyadic variables - kinship, foodsharing
raw.kin <- read.table('kinship/KING.txt')
colnames(raw.kin) <- tolower(colnames(raw.kin))
rownames(raw.kin) <- tolower(rownames(raw.kin))

distmat <- raw.kin %>% 
  select(all_of(bats.used)) %>% 
  filter(row.names(raw.kin) %in% bats.used)

distmat <- data.matrix(distmat)

dyad = list(Kinship = distmat,
            Mouthlicking = mf,
            Grooming = mg,
            Affiliation = ma)

# combine into df
dat = make_strand_data(
  outcome = nets,
  individual_covariates = batchars,
  dyadic_covariates = dyad,
  outcome_mode = "poisson",
  link_mode = "log",
  check_standardization = F) # I don't see why we would standardize age??? but maybe we should


# model 1: kinship -----
fit =
  fit_social_relations_model(
    data=dat,
    focal_regression = ~ Age,
    target_regression = ~ Age,
    dyad_regression = ~ Kinship,
    mode="mcmc",
    stan_mcmc_parameters = list(
      chains = 3,
      iter_warmup = 1500,
      iter_sampling = 1500))

# summarize results
res = summarize_strand_results(fit)

# wrangle model results for plotting
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

# dyadic effects
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

# save plot
ggsave("STRAND_results/kinship1.pdf",
  plot = p1,
  scale = 1,
  width = 7,
  height = 5,
  units = c("in", "cm", "mm", "px"),
  dpi = 300)

# focal, target, and random effects
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
  ggtitle("random effects on number of calls - kinship model")
p2

# save plot
ggsave("STRAND_results/kinship2.pdf",
       plot = p2,
       scale = 1,
       width = 7,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)

# individual effects
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
  xlab("effect - kinship model")
p

# save plot
ggsave("STRAND_results/kinship3.pdf",
       plot = p,
       scale = 1,
       width = 7,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)

# model 2: food-sharing -----
rm(fit, res, plotdat, df, p, p1, p2)

fit =
  fit_social_relations_model(
    data=dat,
    focal_regression = ~ Age,
    target_regression = ~ Age,
    dyad_regression = ~ Mouthlicking,
    mode="mcmc",
    stan_mcmc_parameters = list(
      chains = 3,
      iter_warmup = 1500,
      iter_sampling = 1500))

# summarize results
res = summarize_strand_results(fit)

# wrangle model results for plotting
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

# dyadic effects
p1 <-
  df %>% 
  filter(Variable == "dyadic effects coeffs, Mouthlicking" | Variable == "dyadic effects sd") %>% 
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
  ggtitle("dyadic effects: mouthlicking")
p1

# save plot
ggsave("STRAND_results/foodsharing1.pdf",
       plot = p1,
       scale = 1,
       width = 7,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)

# focal, target, and random effects
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
  ggtitle("random effects on number of calls - mouthlicking model")
p2

# save plot
ggsave("STRAND_results/foodsharing2.pdf",
       plot = p2,
       scale = 1,
       width = 7,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)

# individual effects
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
  xlab("effect - mouthlicking model")
p

# save plot
ggsave("STRAND_results/foodsharing3.pdf",
       plot = p,
       scale = 1,
       width = 7,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)


# model 3: grooming -----
rm(fit, res, plotdat, df, p, p1, p2)

fit =
  fit_social_relations_model(
    data=dat,
    focal_regression = ~ Age,
    target_regression = ~ Age,
    dyad_regression = ~ Grooming,
    mode="mcmc",
    stan_mcmc_parameters = list(
      chains = 3,
      iter_warmup = 1500,
      iter_sampling = 1500))

# summarize results
res = summarize_strand_results(fit)

# wrangle model results for plotting
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

# dyadic effects
p1 <-
  df %>% 
  filter(Variable == "dyadic effects coeffs, Grooming" | Variable == "dyadic effects sd") %>% 
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
  ggtitle("dyadic effects: grooming")
p1

# save plot
ggsave("STRAND_results/grooming1.pdf",
       plot = p1,
       scale = 1,
       width = 7,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)

# focal, target, and random effects
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
  ggtitle("random effects on number of calls - grooming model")
p2

# save plot
ggsave("STRAND_results/grooming2.pdf",
       plot = p2,
       scale = 1,
       width = 7,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)

# individual effects
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
  xlab("effect - grooming model")
p

# save plot
ggsave("STRAND_results/grooming3.pdf",
       plot = p,
       scale = 1,
       width = 7,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)

# model 4: affiliation rate -----
rm(fit, res, plotdat, df, p, p1, p2)

fit =
  fit_social_relations_model(
    data=dat,
    focal_regression = ~ Age,
    target_regression = ~ Age,
    dyad_regression = ~ Affiliation,
    mode="mcmc",
    stan_mcmc_parameters = list(
      chains = 3,
      iter_warmup = 1500,
      iter_sampling = 1500))

# summarize results
res = summarize_strand_results(fit)

# wrangle model results for plotting
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

# dyadic effects
p1 <-
  df %>% 
  filter(Variable == "dyadic effects coeffs, Affiliation" | Variable == "dyadic effects sd") %>% 
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
  ggtitle("dyadic effects: affiliation")
p1

# save plot
ggsave("STRAND_results/affiliation1.pdf",
       plot = p1,
       scale = 1,
       width = 7,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)

# focal, target, and random effects
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
  ggtitle("random effects on number of calls - affiliation model")
p2

# save plot
ggsave("STRAND_results/affiliation2.pdf",
       plot = p2,
       scale = 1,
       width = 7,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)

# individual effects
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
  xlab("effect - affiliation model")
p

# save plot
ggsave("STRAND_results/affiliation3.pdf",
       plot = p,
       scale = 1,
       width = 7,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)

# model 5: food-sharing and kinship -----
rm(fit, res, plotdat, df, p, p1, p2)

fit =
  fit_social_relations_model(
    data=dat,
    focal_regression = ~ Age,
    target_regression = ~ Age,
    dyad_regression = ~ Mouthlicking * Kinship,
    mode="mcmc",
    stan_mcmc_parameters = list(
      chains = 3,
      iter_warmup = 1500,
      iter_sampling = 1500))

# summarize results
res = summarize_strand_results(fit)

# wrangle model results for plotting
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

# dyadic effects
p1 <-
  df %>% 
  filter(Variable == "dyadic effects coeffs, Mouthlicking" |
           Variable == "dyadic effects coeffs, Kinship" |
           Variable == "dyadic effects sd") %>% 
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
  ggtitle("dyadic effects: mouthlicking * kinship")
p1

# save plot
ggsave("STRAND_results/food_kin1.pdf",
       plot = p1,
       scale = 1,
       width = 7,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)

# focal, target, and random effects
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
  ggtitle("random effects on number of calls - mouthlicking * kinship model")
p2

# save plot
ggsave("STRAND_results/food_kin2.pdf",
       plot = p2,
       scale = 1,
       width = 7,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)

# individual effects
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
  xlab("effect - mouthlicking * kinship model")
p

# save plot
ggsave("STRAND_results/food_kin3.pdf",
       plot = p,
       scale = 1,
       width = 7,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)


# model 6: grooming and kinship -----
rm(fit, res, plotdat, df, p, p1, p2)

fit =
  fit_social_relations_model(
    data=dat,
    focal_regression = ~ Age,
    target_regression = ~ Age,
    dyad_regression = ~ Grooming * Kinship,
    mode="mcmc",
    stan_mcmc_parameters = list(
      chains = 3,
      iter_warmup = 1500,
      iter_sampling = 1500))

# summarize results
res = summarize_strand_results(fit)

# wrangle model results for plotting
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

# dyadic effects
p1 <-
  df %>% 
  filter(Variable == "dyadic effects coeffs, Grooming" |
           Variable == "dyadic effects coeffs, Kinship" |
           Variable == "dyadic effects sd") %>% 
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
  ggtitle("dyadic effects: grooming * kinship")
p1

# save plot
ggsave("STRAND_results/groom_kin1.pdf",
       plot = p1,
       scale = 1,
       width = 7,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)

# focal, target, and random effects
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
  ggtitle("random effects on number of calls - grooming * kinship model")
p2

# save plot
ggsave("STRAND_results/groom_kin2.pdf",
       plot = p2,
       scale = 1,
       width = 7,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)

# individual effects
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
  xlab("effect - grooming * kinship model")
p

# save plot
ggsave("STRAND_results/groom_kin3.pdf",
       plot = p,
       scale = 1,
       width = 7,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)


# model 7: affiliation rate and kinship -----
rm(fit, res, plotdat, df, p, p1, p2)

fit =
  fit_social_relations_model(
    data=dat,
    focal_regression = ~ Age,
    target_regression = ~ Age,
    dyad_regression = ~ Affiliation * Kinship,
    mode="mcmc",
    stan_mcmc_parameters = list(
      chains = 3,
      iter_warmup = 1500,
      iter_sampling = 1500))

# summarize results
res = summarize_strand_results(fit)

# wrangle model results for plotting
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

# dyadic effects
p1 <-
  df %>% 
  filter(Variable == "dyadic effects coeffs, Affiliation" |
           Variable == "dyadic effects coeffs, Kinship" |
           Variable == "dyadic effects sd") %>% 
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
  ggtitle("dyadic effects: affiliation * kinship")
p1

# save plot
ggsave("STRAND_results/affiliation_kin1.pdf",
       plot = p1,
       scale = 1,
       width = 7,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)

# focal, target, and random effects
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
  ggtitle("random effects on number of calls - affiliation * kinship model")
p2

# save plot
ggsave("STRAND_results/affiliation_kin2.pdf",
       plot = p2,
       scale = 1,
       width = 7,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)

# individual effects
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
  xlab("effect - affiliation * kinship model")
p

# save plot
ggsave("STRAND_results/affiliation_kin3.pdf",
       plot = p,
       scale = 1,
       width = 7,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)
