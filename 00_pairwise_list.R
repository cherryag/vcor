# make pairwise list of experimental conditions to compare

# CREATED : 3/17/2020 by Cherry Gao


# set up
rm(list = ls()) # clear environment
setwd('/Users/Cherry/Dropbox (MIT)/Lab/Research Files (Dropbox)/02_Vcor paper/CK-transcriptomics_data/TRANSCRIPTOME_ANALYSIS/Ranal')

cond_pair_list = list() # initialize

# list of pairwise condition to compare (first number = numerator, number on top)
cond_pair_list[[1]] = c('0muc' , '0cont')
cond_pair_list[[2]] = c('10muc' , '10cont')
cond_pair_list[[3]] = c('60muc' , '60cont')
cond_pair_list[[4]] = c('10cont' , '0cont')
cond_pair_list[[5]] = c('60cont' , '0cont')
cond_pair_list[[6]] = c('60cont' , '10cont')
cond_pair_list[[7]] = c('10muc' , '0muc')
cond_pair_list[[8]] = c('60muc' , '0muc')
cond_pair_list[[9]] = c('60muc' , '10muc')

# save
save(cond_pair_list, file='pairwise_list.Rdata')
