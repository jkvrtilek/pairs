# analysis for 2024 paired recordings
# Julia Vrtilek
# December 2024

# load packages


# load data


# question 1: do different bats call different amounts? does calling rate depend on recipient?
# use Bayesian social relations model, strand package
# actor, receiver, and interaction automatically included
# response variable = call count
# predictor = Haley's grooming OR food-sharing? maybe just the data from before my data?
# results: effects of actor's tendency to call, receiver's tendency to receive, correlation at dyadic level between caller calling and receiver calling (do the bats I call to more call to me more, controlling for overall calling and receiving rate), correlation at group level (if I call more to the group, does the group call more to me?) and THEN finally will tell us whether bats are calling more depending on predictor variable
# how to make use of added power from 3 replicate matrices? for now, take average, but keep thinking about this

# question 2: do bats sound different depending on recipient?
# for each bat, do a crossed pDFA to assign calls to receiver after controlling for trial (doesn't need to be nested bc we're doing a different one for each bat), for loop through bats
# maybe figure out dynamic time warping? probably not too hard with a normal number of calls

# question 3: do bats sound more like recipient? or are they labeling?

# get all the vocal labeling papers and write down a list of their methods
