## Load packages

library(igraph)
library(ggplot2)
library(ggpubr)

## Set up and prepare the model ##

nsims <- 10000                                            # number of simulations 
n = 50                                                    # number of individuals (nodes)
kingroupraw <- rep(1:10, each=5)                          # vector with kingroups (not randomised order). Each kingroup is given by a number.
typenames <- c('ook', 'oon', 'oyk', 'oyn', 'yyk', 'yyn')  # names of dyad types. Must fit with code and typeparams. o: Old; y: Young; k: kin; n: nonkin
typeparams <- c(0.33, 0.02, 0.37, 0.05, 0.27, 0.08)       # parameters from data for each dyad type. Must be in the same order as the typenames (the format could be changed for security)
noldmin <- 0                                              # minimum number of old individuals
noldmax <- 50                                             # maximum number of old individuals
outputdata <- data.frame(n=rep(n, nsims),                 # storage for output
                         n.old.individuals = NA,
                         prop.old.individuals = NA,
                         perc.old.individuals= NA, 
                         perc.ook=NA, perc.oon=NA, 
                         perc.oyk=NA, perc.oyn=NA, 
                         perc.yyk=NA, 
                         perc.yyn=NA) 


## Run simulations ##

for(cursim in 1:nsims){ # for each simulation (cursim = current simulation)

  print(cursim) # print the current simulation number to the terminal, so that we can follow the progress of the run

  ## draw and store number of group members that are old ####
  
  nold <- sample(noldmin:noldmax,1)                           # number of old individuals for this simulation
  outputdata[cursim, 'n.old.individuals'] <- nold             # store it in the output
  outputdata[cursim, 'prop.old.individuals'] <- nold/n        # store the proportion of old individuals in the output
  outputdata[cursim, 'perc.old.individuals'] <- (nold/n)*100  # store the percentage of old individuals in the output
  
  ## create node data ####
  
  agegroup <- c(rep(1,nold), rep(0, n-nold))  # vector with agegroup for each individual (old = 1, young = 0)
  kingroup <- sample(kingroupraw)             # vector with kingroup for each individual (we must randomise one of the vectors to get independence between age and kingroup, so We randomise this one). 
  
  
  ## change the format to dyad-based ####
  
  # two vectors of individual names that, when matched, give all the unique dyads once
  inames <- rep(1:(n-1),(n-1):1)           # names of all i individuals (the first individual of each dyad)
  jnames <- sequence((n-1):1, from = 2:n)  # names of all j individuals (the second individual of each dyad)
  
  # vectors with the agegroup and kingroup for all i and j individuals, respectively
  agegroupalli <- agegroup[inames]
  agegroupallj <- agegroup[jnames]
  kingroupalli <- kingroup[inames]
  kingroupallj <- kingroup[jnames]
  
  
  ## get dyadtypes ####
  
  ndyad <- length(inames)     # number of dyads
  dyadtype <- rep(NA, ndyad)  # storage for dyadtype
  
  # find and store the type of dyad for each dyad. The type is determined by their ages and kinship.
  dyadtype[agegroupalli == 1 & agegroupallj == 1 & kingroupalli == kingroupallj] <- 'ook'  # old, old, kin
  dyadtype[agegroupalli == 1 & agegroupallj == 1 & kingroupalli != kingroupallj] <- 'oon'  # old, old, nonkin
  dyadtype[agegroupalli == 1 & agegroupallj == 0 & kingroupalli == kingroupallj] <- 'oyk'  # old, young, kin
  dyadtype[agegroupalli == 1 & agegroupallj == 0 & kingroupalli != kingroupallj] <- 'oyn'  # old, young, nonkin
  dyadtype[agegroupalli == 0 & agegroupallj == 0 & kingroupalli == kingroupallj] <- 'yyk'  # young, young, kin
  dyadtype[agegroupalli == 0 & agegroupallj == 0 & kingroupalli != kingroupallj] <- 'yyn'  # young, young, nonkin
  
  
  ## get links ####
  
  link <- rep(NA, ndyad)  # storage for links for all dyads (will contain 1 for connected dyads and 0 for unconnected dyads)
  
  for (typenum in seq_along(typenames)){                                                 # for each dyadtype
    curtypename <- typenames[typenum]                                                    # name of current dyadtype
    curtypeparam <- typeparams[typenum]                                                  # parameter for current dyadtype
    ncurtype <- sum(dyadtype == curtypename)                                             # the number of dyads of the current type
    link[dyadtype == curtypename] <- rbinom(n = ncurtype, size = 1, prob = curtypeparam) # draw and store links for the current dyadtype
  }
  
  
  ## create dyad interaction dataframe
  
  dyaddata <- data.frame(inames, jnames, agegroupalli, agegroupallj, kingroupalli, kingroupallj, dyadtype, link) # dataframe with dyad data, for all dyads
  
  ## output
  
  outputdata[cursim, 'perc.ook'] <- sum(dyaddata$dyadtype=='ook' & dyaddata$link==1)/sum(dyaddata$dyadtype=='ook')
  outputdata[cursim, 'perc.oon'] <- sum(dyaddata$dyadtype=='oon' & dyaddata$link==1)/sum(dyaddata$dyadtype=='oon')
  outputdata[cursim, 'perc.oyk'] <- sum(dyaddata$dyadtype=='oyk' & dyaddata$link==1)/sum(dyaddata$dyadtype=='oyk')
  outputdata[cursim, 'perc.oyn'] <- sum(dyaddata$dyadtype=='oyn' & dyaddata$link==1)/sum(dyaddata$dyadtype=='oyn')
  outputdata[cursim, 'perc.yyk'] <- sum(dyaddata$dyadtype=='yyk' & dyaddata$link==1)/sum(dyaddata$dyadtype=='yyk')
  outputdata[cursim, 'perc.yyn'] <- sum(dyaddata$dyadtype=='yyn' & dyaddata$link==1)/sum(dyaddata$dyadtype=='yyn')
} # end of for loop

## Figure S4 ##

# Prepare the dataset for plotting

outputdata<-na.omit(outputdata) # Omit the NAs, which can occur when there is no dyad of a certain type.
summary(outputdata) # Confirm if the means are as expected

# Prepare the individual plots
Perc.ook_fig <- 
  ggplot(data=outputdata, aes(perc.ooy))+
  geom_histogram(aes(perc.ook),size=1.5) + 
  labs(x= "Proportion of old/ old kin dyads that have a link", y="Count") +
  theme_classic(base_size = 20, base_line_size=1, base_rect_size=1)

Perc.oon_fig <- 
  ggplot(data=outputdata, aes(perc.oon))+
  geom_histogram(aes(perc.oon),size=1.5) + 
  labs(x= "Proportion of old/ old nonkin dyads that have a link", y="Count") +
  theme_classic(base_size = 20, base_line_size=1, base_rect_size=1)

Perc.oyk_fig <- 
  ggplot(data=outputdata, aes(perc.ooy))+
  geom_histogram(aes(perc.ook),size=1.5) + 
  labs(x= "Proportion of old/ young kin dyads that have a link", y="Count") +
  theme_classic(base_size = 20, base_line_size=1, base_rect_size=1)

Perc.oyn_fig <- 
  ggplot(data=outputdata, aes(perc.oon))+
  geom_histogram(aes(perc.oon),size=1.5) + 
  labs(x= "Proportion of old/ young nonkin dyads that have a link", y="Count") +
  theme_classic(base_size = 20, base_line_size=1, base_rect_size=1)

Perc.yyk_fig <- 
  ggplot(data=outputdata, aes(perc.ooy))+
  geom_histogram(aes(perc.ook),size=1.5) + 
  labs(x= "Proportion of young/ young kin dyads that have a link", y="Count") +
  theme_classic(base_size = 20, base_line_size=1, base_rect_size=1)

Perc.yyn_fig <- 
  ggplot(data=outputdata, aes(perc.oon))+
  geom_histogram(aes(perc.oon),size=1.5) + 
  labs(x= "Proportion of young/ young nonkin dyads that have a link", y="Count") +
  theme_classic(base_size = 20, base_line_size=1, base_rect_size=1)

# Arrange the plots in a single figure
ggarrange(Perc.ook_fig,Perc.oon_fig,Perc.oyk_fig,Perc.oyn_fig,Perc.yyk_fig,Perc.yyn_fig,nrow = 3, ncol=2, labels="AUTO")

