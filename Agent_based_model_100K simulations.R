## Load packages

library(igraph)

## Set up and prepare the model ##

nsims <- 100000                                           # number of simulations 
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
                         mean.degree= NA, 
                         transitivity= NA,
                         diameter= NA) 


## Run simulations ##

for(cursim in 1:nsims){  # for each simulation (cursim = current simulation)
  
  print(cursim)  # print the current simulation number to the terminal, so that we can follow the progress of the run

  ## draw and store number of group members that are old ####
  
  nold <- sample(noldmin:noldmax,1)                           # number of old individuals for this simulation
  outputdata[cursim, 'n.old.individuals'] <- nold             # store the number of old individuals in the output
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
  
  
  ## create igraph network ####
  
  dyaddata <- data.frame(inames, jnames, agegroupalli, agegroupallj, kingroupalli, kingroupallj, dyadtype, link) # dataframe with dyad data, for all dyads
  edgelistdata <- dyaddata[link==1,]                                              # dataframe with edgelist in the two first columns, link attributes in the rest (i.e. dyads without links are not included here)
  nodedata <- data.frame(nodenames = 1:n, agegroup, kingroup)                     # dataframe with node data (not sure whether this is needed but maybe igraph uses it to see nodes that have no links)
  inet <-graph_from_data_frame(d=edgelistdata, vertices=nodedata, directed=FALSE) # get network in igraph format
  
  
  ## calculate and store network metrics ####
  
  curdegrees <- degree(inet) # degrees for all individuals
  
  # global network metrics
  outputdata[cursim, 'mean.degree'] <- mean(curdegrees)
  outputdata[cursim, 'transitivity'] <- transitivity(inet, type = "global") 
  outputdata[cursim, 'diameter'] <- diameter(inet, unconnected = T) 
  
} ## end of for loop
