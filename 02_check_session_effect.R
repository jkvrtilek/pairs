# script to assess effect of session on number of calls
# 3 sessions for every dyad
library(tidyverse)

setwd("/Users/jkvrtilek/Desktop/OSU/PhD/GitHub/pairs")

# load data
batcalls <- readRDS("vocal_data_2024-pairs_transformed.RDS") %>% 
  mutate(date = as.Date(date)) %>% 
  mutate(session = case_when(date < "2024-07-22" ~ 1,
                             date > "2024-07-21" & date < "2024-07-29" ~ 2,
                             date > "2024-07-28" ~ 3))

df <- batcalls %>% 
  mutate(dyad = paste(caller, receiver, sep = "-")) %>% 
  group_by(dyad, session) %>% 
  summarize(n = n()) %>% 
  ungroup()

# plot calls per dyad per session
p <- df %>%   
  ggplot(aes(x=session, y=n, group = dyad))+
  geom_point(size=1)+
  geom_line()+
  xlab("session")+
  ylab("calls")+
  theme_bw()
p

# ANOVA: effect of session on call rate with dyad as random variable
sess_aov <- aov(n ~ session + dyad, data = df)
summary(sess_aov)




