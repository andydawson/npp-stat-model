
# MK ??: Why do we save each step?

rm(list=ls())
setwd('~/Desktop/npp-stat-model')

library(ggplot2)
library(reshape2)

# MK ??: What is the config file? Do we need it? How do we edit it?
#source('config')

#############################################################################################################
### A. Set up working environment and directory for script ##################################################
#############################################################################################################

########################################
# Variables to be adjusted for each site
########################################

# date for save purposes
date = '13Nov2019'

# set up output 
site <- 'RH'

# does site have census data available?
census = FALSE

# MK ??: What is the difference between these two version numbers?
# data version
dvers = 'v5.0'
mvers = 'v4.0'

# burn-in range
burn = 1000

# load model output data for site
data = readRDS('~/Downloads/tree_data_ROOSTER_STAN_v0.1.RDS')

# MK : in the script, I assume that there either two models when there is census data or one model when there isn't
# when census data exists, use census model first
models = c('Model RW')
#models = c('Model RW + Census', 'Model RW')
#fnames = c('ring_model_t_date_sapl_size_pdbh', 'ring_model_t_date_sapl_size_pdbh_nc')
#fnames = c('ring_model_t_date_sapl_size_pdbh_NOCOVAR', 'ring_model_t_pdbh_nc_NOCOVAR_sigd')

#fname_data = paste0('tree_data_20_', dvers)
#load(file=paste0('data/dump/', fname_data, '.rdata'))
#fname_data = paste0('tree_data_20_no_census_', dvers)
#load(file=paste0('data/dump/', fname_data, '.rdata'))

# are dbh and ab files already created? do we want to recreate them? 
run_predict  = TRUE

########################################

# create output directory
output_dir <- file.path('output',site)
if (!file.exists(output_dir)){
  dir.create(output_dir)
  dir.create(file.path(output_dir, 'figures'))
}

# MK ??: preference on set-up of directory?
# get ring width data 
nmodels = length(fnames)
post = list()
for (i in 1:length(fnames)) {
  fname_model = fnames[i]
  load(file   = paste0('data/dump/', fname_model, '_', mvers, '.Rdata'))
  post[[i]]   = out[1:2000,]
}  
niter = dim(post[[1]])[1]

# how many iterations to keep?
keep = niter - burn

# MK ??: Do species always follow this same order? Or is it through the order given in taxa match?
# pft list 
pfts = list(ACRU = data.frame(spcd=316,acronym="ACRU"),
            ACSA = data.frame(spcd=318,acronym="ACSA3"),
            BEAL = data.frame(spcd=371,acronym="BEAL"),
            BELE = data.frame(spcd=372),
            BEPO = data.frame(spcd=379,acronym="BEPO"),
            CADE = data.frame(spcd=1000,acronym="CADE"), # id is 421, but lump into mixed hardwoods because no biomass eq
            FAGR = data.frame(spcd=531,acronym="FAGR"),
            FRAM = data.frame(spcd=541, acronym="FRAM"),
            HAVI = data.frame(spcd=1000,acronym="HAVI"), # id is 585, but lump into mixed hardwoods because no biomass eq
            PIST = data.frame(spcd=129,acronym="PIST"),
            PRSE = data.frame(spcd=762,acronym="PRSE"),
            QURU = data.frame(spcd=833,acronym="QURU"),
            QUVE = data.frame(spcd=837,acronym="QUVE"),
            QUAL = data.frame(spcd=802, acronym="QUAL"),
            TSCA = data.frame(spcd=261,acronym="TSCA")
)
pft_list  = sort(unique(names(pfts)))

#############################################################################################################
### B. Process diameter and ab predictions for ring width data ##############################################
#############################################################################################################

##### Step 1: Estimate diameter and aboveground biomass for all trees from ring width data #####

# Input: NPP Stat Model ring width results for the two different models
# Output: 4 large matrices holding the dbh and aboveground biomass predictions for the two different models (dbh_p_1_raw)

# match species acronyms to level3a available species/pfts 
acro_level3a = read.csv('data/acronym_to_level3a_v0.1.csv', stringsAsFactors = FALSE)
acro_level3a = left_join(data$taxaMatch, acro_level3a, by = c('species'='acronym'))
choj = read.csv('data/level3a_to_chojnacky_v0.4.csv')
choj = left_join(acro_level3a, choj, id = 'level3a') 

# this function takes the stat model results (dbh, rw) and returns a list of dbh arrays (tree X years X iterations)
# for each pft in the data, a list of the tree indices in the order they were used in the list, and the order of the
# pfts within the list 
build_dbh_p <- function(out, rw2year, rw2tree, Nyears, tree_list, keep, taxaMatch, taxon_list){
  
  # extract colnames and number of iterations (different value types in data, want diameter)
  col_names = sapply(strsplit(colnames(out), '\\['), function(x) x[[1]])
  niter   = dim(out)[1]
  
  # get species from taxon list 
  #pft   = as.vector(taxaMatch[match(taxon_list, taxaMatch$taxon),'species'])
  pft = as.vector(taxaMatch[match(taxon_list, taxaMatch$number),'species'])
  pfts  = sort(unique(pft))
  N_pft = length(pfts) 
  
  # create list for dbh data for each pft 
  dbh_p = list(length=N_pft)
  tree_order = vector(length=0)
  
  # loop through pfts
  for (p in 1:N_pft){
    
    # identify all trees of this pft and put in list to track 
    pft_trees = which(pft == pfts[p])
    tree_order = c(tree_order, tree_list[pft_trees])
    
    # create save array to put dbh information for pft in 
    # dbh saved for keep iterations for every year and every tree of each pft 
    N_pft_trees = length(pft_trees)
    dbh_p[[p]] = array(NA, c(N_pft_trees, Nyears, keep))
    
    # loop through trees of pft
    for (i in 1:N_pft_trees) {
      print(paste0('Tree ', i))
      tree = tree_list[pft_trees[i]]
      
      # connects individual tree to place in out and which years there is data for the tree 
      tree_idx   = which(rw2tree == tree)
      tree_years = rw2year[tree_idx]
      
      # extracts all diameter measurements for tree in iterations to keep 
      mu_dbh = out[(niter-keep+1):(niter), which(col_names=='D')[tree_idx]]
      
      # adjust negative dbh values
      if (any(mu_dbh < 0)) {
        mu_dbh[mu_dbh < 0] = 0
      }
      
      # saves extracted dbh measurements into correct years for the pft
      dbh_p[[p]][i,tree_years,] = t(mu_dbh)
    }
  }
  return(list(dbh_p=dbh_p, tree_order=tree_order, pft_order = pfts))
}

# this function takes the built dbh array list and returns a biomass array (tree X years X iterations) where
# biomass is calculated using Chojnacky's allometric coefficients in the biomass equation
build_ab_choj <- function(dbh_p, Ntrees, Nyears, keep, choj_pfts){
  
  # build save array 
  ab_p = array(NA, c(Ntrees, Nyears, keep))
  
  # loop through iterations to keep 
  for (k in 1:keep){
    print(paste0("Iteration ", k))
    
    # get dbhs for all trees for all years in iteration (number of trees X years)
    dbh_list = dbh_p[,,k]
    
    # calculate the biomass contribution of each tree for each year in iteration
    # ln(biomass-kg) = beta0 + beta1 * ln(diameter-cm)
    b0 = choj$beta0[choj_pfts]
    b1 = choj$beta1[choj_pfts]
    ab_p[,,k] = exp(b0 + b1 * log(dbh_list))
  }
  
  return(ab_p=ab_p)
}

# if ab and diameter estimates haven't already been found 
if (run_predict){
  
  # for first model with census 
  # obtain diameter estimates for all trees for
  dbh_1 = build_dbh_p(post[[1]], data$x2year, data$x2tree, data$N_years, trees, keep, data$taxaMatch, data$taxon)
  dbh_p_1 = dbh_1$dbh_p
  tree_order_1 = dbh_1$tree_order
  pft_order_1  = dbh_1$pft_order
  
  # make list of arrays into one large array across species (all trees X years X iterations)
  dbh_p_1_org = abind(dbh_p_1, along=1)
  dbh_p_1 = dbh_p_1_org[sort(tree_order_1, index.return=TRUE)$ix,,]
  
  # get choj pft index for coefficients in matrix 
  choj_pfts_1 = tree_order_1
  
  # find biomass for all trees, years, and iterationss
  ab_p_1 = build_ab_choj(dbh_p_1, data$N_trees, data$N_years, keep, choj_pfts_1)
  
  # save files 
  #saveRDS(dbh_p_1, paste0('output/',site,'/dbh_p_1_', site, '_', mvers, '.rds'))
  #saveRDS(ab_p_1, paste0('output/',site,'/ab_p_1_', site, '_', mvers, '.rds'))
  
  # organize tree info for data with no census 
  trees_p = sort(unique(datat$x2tree_p)) # get all trees ids in data
  taxon_p = data$taxon[trees_p] # gets taxons of unique trees
  x2idx_p = match(x2tree_p, unique(data$x2tree_p)) # matches all trees from data to a numerical index 
  pft_list_p = data$taxaMatch$species[taxon_p] # matches all trees to a pft
  
  # is there another model to process?
  if (census){
    # get diameter estimates from ring widths given by stat model 
    dbh_2 = build_dbh_p(post[[2]], x2year_p, x2tree_p, N_years, trees_p, keep, taxaMatch, taxon_p)
    dbh_p_2 = dbh_2$dbh_p
    tree_order_2 = dbh_2$tree_order
    pft_order_2  = dbh_2$pft_order
    
    # make list of arrays into one big array (all trees X years X iterations)
    dbh_p_2_org = abind(dbh_p_2, along=1)
    dbh_p_2 = dbh_p_2_org[sort(tree_order_2, index.return=TRUE)$ix,,]
    
    # get choj pft index for coefficients in matrix 
    choj_pfts_2 = tree_order_2
    
    # find biomasss for all trees, years, and iterations
    ab_p_2 = build_ab_choj(dbh_p_2, N_trees_p, N_years, keep, choj_pfts_2)
    
    # save files 
    #saveRDS(dbh_p_2, paste0('output/',site,'/dbh_p_2_', site, '_', mvers, '.rds'))
    #saveRDS(ab_p_2, paste0('output/',site,'/ab_p_2_', site, '_', mvers, '.rds'))
    
  }
}else {
  
  dbh_p_1 = readRDS( paste0('output/dbh_p_1_', site, '_', mvers, '.rds'))
  ab_p_1  = readRDS(paste0('output/ab_p_1_',  site, '_', mvers, '.rds'))
  
  if (census){
    dbh_p_2 = readRDS( paste0('output/dbh_p_2_', site, '_', mvers, '.rds'))
    ab_p_2  = readRDS(paste0('output/ab_p_2_', site, '_', mvers, '.rds'))
  }
}

# added separate save variable in case we need to go back to this step
dbh_p_1_raw = dbh_p_1
ab_p_1_raw = ab_p_1

if (census){
  dbh_p_2_raw = dbh_p_2
  ab_p_2_raw = ab_p_2
}

##### Step 2: Save all trees (< 5cm) to RDS files for both models #####

# Input: 4 large matrices holding dbh and ab predictions for the two different models
# Output: Reduced matrices holding dbh and ab predictions, with trees with mean DBH less than 5 cm removed 

# remove trees that have a mean diameter less than 5 centimeters across all iterations for a given year for both models
for (tree in 1:data$N_trees){
  for (year in 1:data$N_years){
    dbh_mean = mean(dbh_p_1[tree, year, ], na.rm=TRUE)
    if (is.na(dbh_mean)){next}
    if (dbh_mean < 5){
      dbh_p_1[tree, year, ] = rep(NA, keep)
      ab_p_1[tree, year, ] = rep(NA, keep)
    }
  }
}

if (census){
  for (tree in 1:N_trees_p){
    for (year in 1:N_years){
      dbh_mean = mean(dbh_p_2[tree, year, ], na.rm=TRUE)
      if (is.na(dbh_mean)){next}
      if (dbh_mean < 5){
        dbh_p_2[tree, year, ] = rep(NA, keep)
        ab_p_2[tree, year, ] = rep(NA, keep)
      }
    }
  }
  
  dbh_p_2_red = dbh_p_2
  ab_p_2_red = ab_p_2
}

# added save variable in case we need to back to this step
dbh_p_1_red = dbh_p_1
ab_p_1_red = ab_p_1

##### Step 3: Adjust for mortality using census data #####

# MK ??: What do we do for sites where census data is not available? Skip this step?

# Input: Reduced diameter and biomass values for all trees for all years and iterations, census data
# Output: Diameter and biomass values for all trees for all years and iterations smoothed for mortality seen in census data 

if(census){
  
  # get the census years
  cyr = dbh$yr
  cyr[which(cyr %in% c('01', '11'))] = paste0('20', cyr[which(cyr %in% c('01', '11'))])
  cyr[which(cyr %in% c('62', '69', '75', '91'))] = paste0('19', cyr[which(cyr %in%  c('62', '69', '75', '91'))])
  census_years = as.numeric(sort(unique(cyr)))
  idx_census   = which(years %in% census_years)
  N_census_years = length(census_years)
  dbh = data.frame(dbh, cyr=cyr)
  
  # MK ??: What is the logic behind the smoothing function?
  smooth_death <- function(dbh_p, ab_p, last_time_data, last_time, N_trees, years, keep){
    
    # last_time_data is the last time each tree was measured; only adjust prior to 2011
    last_time_data[last_time_data == 1992] = 1991
    idx_adjust = which(last_time_data < 2011)
    
    dbh_p_smooth = dbh_p
    ab_p_smooth  = ab_p
    
    for (i in 1:N_trees) {
      print(i)
      
      # index of last year measured with data 
      idx_last = which(years == last_time_data[i])
      
      if (idx_last >= 52) {
        next
        
        # fill in for death that died between last measured and first census without tree
      } else if (last_time_data[i] < last_time[i]) {
        
        print('Adjusting dbh and ab due to death!')
        
        # difference in years
        sample_int =  last_time[i] - last_time_data[i]
        
        for (k in 1:keep){
          # pick random year in time span for time of death 
          mort_year = sample(sample_int,1)
          years_na = seq((last_time_data[i]+mort_year), last_time[i])
          idx_na   = which(years %in% years_na)
          
          for (idx in idx_na){
            dbh_p_smooth[i,idx,k] = NA
            ab_p_smooth[i,idx,k] = NA
          }
        }
      }
    }
    return(list(dbh_p_smooth=dbh_p_smooth, ab_p_smooth=ab_p_smooth))
  }
  
  # smooth out dbh and aboveground biomass for both models
  smooth_1 = smooth_death(dbh_p_1, ab_p_1, last_time_data, last_time_census, N_trees, years, keep)
  dbh_p_1 = smooth_1$dbh_p_smooth
  ab_p_1 = smooth_1$ab_p_smooth
  
  smooth_2 = smooth_death(dbh_p_2, ab_p_2, last_time_data_p, last_time_census_p, N_trees_p, years, keep)
  dbh_p_2 = smooth_2$dbh_p_smooth
  ab_p_2  = smooth_2$ab_p_smooth 
}

##### Step 4: Melt data matrices into data frames for easier use #####

# melt down to data frames
dimnames(dbh_p_1)[[1]] <- seq(1, 603)
dbh_melt = melt(dbh_p_1)
colnames(dbh_melt) = c("tree", "year", "iter", "dbh")
dbh_melt$year = years[dbh_melt$year]

dimnames(ab_p_1)[[1]] <- seq(1, 603)
ab_melt = melt(ab_p_1)
colnames(ab_melt) = c("tree", "year", "iter", "ab")
ab_melt$year = years[ab_melt$year]

if (census){
  dimnames(dbh_p_2)[[1]] <- trees_p
  dbh2_melt = melt(dbh_p_2)
  colnames(dbh2_melt) = c("tree", "year", "iter", "dbh")
  dbh2_melt$year = years[dbh2_melt$year]
  
  dimnames(ab_p_2)[[1]] <- trees_p
  ab2_melt = melt(ab_p_2)
  colnames(ab2_melt) = c("tree", "year", "iter", "ab")
  ab2_melt$year = years[ab2_melt$year]
}

if (!census){
  # put models together
  dbh_melt_all = data.frame(dbh_melt, model=rep("Model RW", nrow(dbh_melt)))
  #dbh_melt_all$plot = tree_site_id[dbh_melt_all$tree]
  dbh_melt_all$plot = plot_id[dbh_melt_all$treee]
  
  ab_melt_all = data.frame(ab_melt, model=rep("Model RW", nrow(ab_melt)))
  ab_melt_all$plot = plot_id[ab_melt_all$tree]
  #ab_melt_all$plot = tree_site_id[ab_melt_all$tree]
  
  ab_mean = aggregate(ab~tree+plot+year+model, ab_melt_all, mean) 
  ab_sum = aggregate(ab~year+model+plot, ab_mean, sum)
  
  # save to new variable in case we need to go back to this step
  dbh_p_1_smth = dbh_p_1
  ab_p_1_smth = ab_p_1
}else{
  
  # put models together
  dbh_melt_all = rbind(data.frame(dbh_melt, model=rep("Model RW + Census", nrow(dbh_melt))), 
                       data.frame(dbh2_melt, model=rep("Model RW", nrow(dbh2_melt))))
  #dbh_melt_all$plot = tree_site_id[dbh_melt_all$tree]
  dbh_melt_all$plot = plot_id[dbh_melt_all$treee]
  dbh_melt_all$census_id = dbh[match(dbh_melt_all$tree, dbh$stat_id),'census_id']
  dbh_melt_all$dist_census = census[match(dbh_melt_all$census_id, census$census_id),'dist_census']/3.28084
  
  ab_melt_all = rbind(data.frame(ab_melt, model=rep("Model RW + Census", nrow(ab_melt))), 
                      data.frame(ab2_melt, model=rep("Model RW", nrow(ab2_melt))))
  #ab_melt_all$plot = tree_site_id[ab_melt_all$tree]
  ab_melt_all$plot = plot_id[ab_melt_all$tree]
  
  ab_mean = aggregate(ab~tree+plot+year+model, ab_melt_all, mean) 
  ab_sum = aggregate(ab~year+model+plot, ab_mean, sum)
  
  # save to new variable in case we need to go back to this step
  dbh_p_1_smth = dbh_p_1
  ab_p_1_smth = ab_p_1
  dbh_p_2_smth = dbh_p_2
  ab_p_2_smth = ab_p_2
}


##### Step 5: Calculate biomass and biomass increment #####

# Input: Smoothed, reduced diameter and biomass estimates for all trees for all years and estimates in four matrices (ab in Kg/plot)
# Output: Melted data frames for aboveground biomass (Mg/ha) and aboveground biomass increment with smoothed and reduced data for
#         all iterations, years, and trees

# in Kg/plot, rescale so it is Mg/ha
ab_pr1 = ab_p_1/(20^2*pi) * 10

# melt predicted biomass for both models
ab_p1_melt = melt(ab_pr1)
colnames(ab_p1_melt) = c('tree_id', 'year_id', 'iter', 'ab')
#ab_p1_melt$site_id = tree_site_id[ab_p1_melt$tree_id]
ab_p1_melt$site_id = data$plot_id[ab_p1_melt$tree_id]

if (census){
  
  # in Kg/plot, rescale so it is Mg/ha
  ab_pr2 = ab_p_2/(20^2*pi) * 10
  
  # melt predicted biomass matrix
  ab_p2_melt = melt(ab_pr2)
  colnames(ab_p2_melt) = c('tree_id', 'year_id', 'iter', 'ab')
  #ab_p2_melt$site_id = tree_site_id[ab_p2_melt$tree_id]
  ab_p2_melt$site_id = data$plot_id[ab_p2_melt$tree_id]
  
  # put them together
  ab_p_melt = rbind(data.frame(ab_p1_melt, model=rep('Model RW + Census')), 
                    data.frame(ab_p2_melt, model=rep('Model RW')))
  ab_p_melt$taxon   = taxon[ab_p_melt$tree_id]
  ab_p_melt$taxon   = data$taxaMatch[match(ab_p_melt$taxon, taxaMatch$number), 'species']
  #ab_p_melt$taxon   = data$taxaMatch[match(ab_p_melt$taxon, taxaMatch$taxon), 'species']
  ab_p_melt$year = data$years[ab_p_melt$year_id]
  
}else{
  
  # create total 
  ab_p_melt = data.frame(ab_p1_melt, model=rep('Model RW'))
  ab_p_melt$taxon   = data$taxon[ab_p_melt$tree_id]
  ab_p_melt$taxon   = data$taxaMatch[match(ab_p_melt$taxon, data$taxaMatch$number), 'species']
  #ab_p_melt$taxon   = data$taxaMatch[match(ab_p_melt$taxon, data$taxaMatch$taxon), 'species']
  ab_p_melt$year = data$years[ab_p_melt$year_id]
  
}

# save aboveground biomass predictions
#saveRDS(ab_p_melt, file=paste0('output/',site,'/ab_p_', site, '_', mvers, '.RDS'))

# determine biomass increment for model 1 
abi_p1 = apply(ab_pr1, c(1,3), function(x) diff(x))
abi_p1 = aperm(abi_p1, c(2, 1, 3))
abi_p1_melt = melt(abi_p1)
colnames(abi_p1_melt) = c('tree_id', 'year_id', 'iter', 'abi')

if (census){
  
  # determine biomass increment for model 2
  abi_p2 = apply(ab_pr2, c(1,3), function(x) diff(x))
  abi_p2 = aperm(abi_p2, c(2, 1, 3))
  abi_p2_melt = melt(abi_p2)
  colnames(abi_p2_melt) = c('tree_id', 'year_id', 'iter', 'abi')
  
  # put them together
  abi_p_melt = rbind(data.frame(abi_p1_melt, model=rep('Model RW + Census')), 
                     data.frame(abi_p2_melt, model=rep('Model RW')))
  abi_p_melt$taxon   = data$taxon[abi_p_melt$tree_id]
  abi_p_melt$taxon   = data$taxaMatch[match(abi_p_melt$taxon, data$taxaMatch$number), 'species']
  #abi_p_melt$taxon   = data$taxaMatch[match(abi_p_melt$taxon, data$taxaMatch$taxon), 'species']
  #abi_p_melt$site_id = tree_site_id[abi_p_melt$tree_id]
  abi_p_melt$site_id = data$plot_id[abi_p_melt$tree_id]
  abi_p_melt$year    = data$years[abi_p_melt$year_id]
  
}else{
  
  # put them together
  abi_p_melt = data.frame(abi_p1_melt, model=rep('Model RW + Census'))
  abi_p_melt$taxon   = data$taxon[abi_p_melt$tree_id]
  abi_p_melt$taxon   = data$taxaMatch[match(abi_p_melt$taxon, data$taxaMatch$number), 'species']
  #abi_p_melt$taxon   = data$taxaMatch[match(abi_p_melt$taxon, data$taxaMatch$taxon), 'species']
  #abi_p_melt$site_id = tree_site_id[abi_p_melt$tree_id]
  abi_p_melt$site_id = data$plot_id[abi_p_melt$tree_id]
  abi_p_melt$year    = data$years[abi_p_melt$year_id]
  
}

# save aboveground biomass increment predictions
#saveRDS(abi_p_melt, file=paste0('output/',site,'/abi_p_', site, '_', mvers, '.RDS'))

##### Step 6: Prep AGB and ABI for wiki #####

# Saves single data frame containing all aboveground biomass values and aboveground biomass increment values
# for all trees, years, and iterations in one single large data frame 

abi_p_melt = abi_p_melt[,c('tree_id', 'year', 'site_id', 'taxon', 'model', 'iter', 'abi')]
ab_p_melt = ab_p_melt[,c('tree_id', 'year', 'site_id', 'taxon', 'model', 'iter', 'ab')]

## MK ??: Why 2012? Is this HF specific?
ab_p_melt = ab_p_melt[ab_p_melt$year<2013,]
abi_p_melt = abi_p_melt[abi_p_melt$year<2013,]

colnames(ab_p_melt)[which(colnames(ab_p_melt) == "ab")] = 'value'
colnames(abi_p_melt)[which(colnames(abi_p_melt) == "abi")] = 'value'

#############################################################################################################
### C. Process dbh and ab from increments and most recent dbh measurement ###################################
#############################################################################################################

##### Step 1: Track diameter back from most recent dbh measurement using increment data #####

# Input: Raw data from models on dbh increment (mm) and most recent dbh measurement (cm)
# Output: Matrix containing diameter measurements (cm) for all trees for years with increment data given 

N_samples = 200
dbh_m = array(NA, c(data$N_trees, data$N_years))

# loop through all trees 
for (i in 1:data$N_trees) {
  print(i)
  tree = data$trees[i]
  pft  = as.vector(data$taxaMatch[which(data$taxaMatch$taxon == data$taxon[tree]),'species'])
  
  # if tree doesn't have increment data, make NA for all years 
  if (length(which(data$m2tree_a == tree)) == 0) {
    print(paste0('No increments for tree ', tree, ' !'))
    
    dbh_idx   = which(dbh$stat_id == tree)
    dbh_years = dbh[dbh_idx, 'year']
    
    dbh_dat_tree = dbh[dbh_idx,]
    
    dbh_m[tree, ] = rep(NA, N_years)
    
    # if tree does have increment data, track growth 
  } else {
    
    # for now use average increments
    incr_tree  = exp(logXobs_a[which(m2tree_a == tree)])
    incr_years = years[m2ti_a[which(m2tree_a == tree)]]
    
    # get diameter measurements for tree
    dbh_idx   = which(dbh$stat_id == tree)
    dbh_years = dbh[dbh_idx, 'year']
    pdbh_idx = which(pdbh$stat_id == tree)
    pdbh_year = pdbh[pdbh_idx, 'dbh_year']
    
    # diameter data available?
    if (sum(dbh_years %in% c(incr_years, pdbh_year)) == 0){
      next
    }
    
    # get most recent diameter measurement for tree 
    if (length(pdbh_idx) > 0){
      dbh_year = pdbh_year
      dbh_tree = pdbh[pdbh_idx, 'value']
    } else {
      dbh_idx  = dbh_idx[dbh_years %in% incr_years]
      dbh_dat_tree = dbh[dbh_idx,]
      dbh_year = dbh_dat_tree[which.max(dbh_dat_tree$year),'year']
      dbh_tree = dbh_dat_tree[which.max(dbh_dat_tree$year),'value']
    }
    
    # MK ??: What's the deal with these trees?
    # tree 139
    if (tree==139){
      dbh_year = 2009
      dbh_tree = 46.3
    } 
    if (tree == 145) {
      dbh_year = 2007
      dbh_tree = 13.3
    }
    
    # MK ??: Shouldn't we base this partially on what year our diameter value is given?
    # remove last index of increment data 
    incr_years = incr_years[-(length(incr_years))]
    incr_tree = incr_tree[-(length(incr_tree))]
    tree_years = seq(min(incr_years), max(incr_years+1))
    
    # determine diameter each year based on increment 
    incr_cumsum = rev(cumsum(rev(incr_tree[incr_years %in% tree_years])))
    dbh_calc = c(dbh_tree - 2*incr_cumsum/10, dbh_tree)
    
    # line up years for data
    year_idx = match(tree_years, years)
    
    if (length(dbh_calc) > length(year_idx)){
      dbh_calc = dbh_calc[-1]
    }
    dbh_m[tree,year_idx] = dbh_calc
  }
  
  if (last_ti[i] < length(years)) {
    dbh_m[tree, (last_ti[i]+1):N_years] = rep(NA, N_years - last_ti[i])
  }
}

##### Step 2: Get biomass using Chojnacky coefficients for all trees #####

# Input: Matrix containing diameter information for all trees for all years 
# Output: 

# get choj pft index for coefficients in matrix 
choj_pfts_m = as.vector(sapply(pft_list[taxon], function(x){which(acro_level3a$acronym == x)}))

# ln(biomass-kg) = beta0 + beta1 * ln(diameter-cm)
b0 = choj$beta0[choj_pfts_m]
b1 = choj$beta1[choj_pfts_m]
ab_m = exp(b0 + b1*log(dbh_m))

# gets rid of data for trees with diameter less than 5 cm
for (n in 1:N_trees){
  for (t in 1:N_years){
    if (is.na(dbh_m[n,t])){
      next
    } else if (dbh_m[n,t] < 5) {
      dbh_m[n,t] = NA
      ab_m[n,t]  = NA
    }  
  }
}

# in Kg/plot, rescale so it is Mg/ha
ab_mr = ab_m/(20^2*pi) * 10

# melt measured
ab_m_melt = melt(ab_mr)
colnames(ab_m_melt) = c('tree_id', 'year_id', 'ab')
ab_m_melt$taxon   = data$taxon[ab_m_melt$tree_id]
#ab_m_melt$taxon   = taxaMatch[match(ab_m_melt$taxon, taxaMatch$taxon), 'species']
ab_m_melt$taxon   = data$taxaMatch[match(ab_m_melt$taxon, taxaMatch$number), 'species']
#ab_m_melt$site_id = tree_site_id[ab_m_melt$tree_id]
ab_m_melt$site_id = data$plot_id[ab_m_melt$tree_id]
ab_m_melt$year    = data$years[ab_m_melt$year_id]

## MK ??: Do we need to remove these dates for abi, median, etc?
ab_m_melt = ab_m_melt[ab_m_melt$year<2013,]

# measured
#ab_mr_median = apply(ab_mr, c(1,2), median, na.rm=TRUE)
abi_m = t(apply(ab_mr, 1, function(x) diff(x)))

# melt measured
abi_m_melt = melt(abi_m)
colnames(abi_m_melt) = c('tree_id', 'year_id', 'abi')
abi_m_melt$taxon   = data$taxon[abi_m_melt$tree_id]
#abi_m_melt$taxon   = taxaMatch[match(abi_m_melt$taxon, taxaMatch$taxon), 'species']
abi_m_melt$taxon   = data$taxaMatch[match(abi_m_melt$taxon, data$taxaMatch$number), 'species']
#abi_m_melt$site_id = tree_site_id[abi_m_melt$tree_id]
abi_m_melt$site_id = data$plot_id[abi_m_melt$tree_id]
abi_m_melt$year    = data$years[abi_m_melt$year_id]

abi_m_melt = abi_m_melt[abi_m_melt$year<2013,]

#############################################################################################################
### D. Process dbh and ab from census data ##################################################################
#############################################################################################################

if (census){

  ##### Step 1: Organize census data #####
  
  ## MK??: why are 'yr' and 'year' different?
  # get census years
  # census_years = sort(unique(dbh$year))
  cyr = dbh$yr
  cyr[which(cyr %in% c('01', '11'))] = paste0('20', cyr[which(cyr %in% c('01', '11'))])
  cyr[which(cyr %in% c('62', '69', '75', '91'))] = paste0('19', cyr[which(cyr %in%  c('62', '69', '75', '91'))])
  census_years = sort(unique(cyr))
  
  dbh = data.frame(dbh, cyr=cyr)
  N_census_years = length(census_years)
  N_samples = 200
  
  # organize data using dcast 
  dbh_c_cast = dcast(dbh, stat_id+site~cyr, value.var=c('value'))
  
  ##### Step 2: Get biomass estimates based on Chojnacky coefficients #####
  
  dbh_c = dbh_c_cast[,3:ncol(dbh_c_cast)]
  
  # get Choj pft indices for getting correct coefficients
  choj_pfts_c = as.vector(sapply(pft_list[taxon], function(x){which(acro_level3a$acronym == x)}))
  
  # ln(biomass-kg) = beta0 + beta1 * ln(diameter-cm)
  b0 = choj$beta0[choj_pfts_c]
  b1 = choj$beta1[choj_pfts_c]
  ab_c = exp(b0 + b1*log(dbh_c))
  
  # in Kg/plot, rescale so it is Mg/ha
  ab_cr = as.matrix(ab_c/(20^2*pi) * 10)
  
  # melt census
  ab_c_melt = melt(ab_cr)
  colnames(ab_c_melt) = c('tree_id', 'year_id', 'ab')
  ab_c_melt$taxon   = taxon[ab_c_melt$tree_id]
  #ab_c_melt$taxon   = taxaMatch[match(ab_c_melt$taxon, taxaMatch$taxon), 'species']
  ab_c_melt$taxon   = taxaMatch[match(ab_c_melt$taxon, taxaMatch$number), 'species']
  ab_c_melt$site_id = tree_site_id[ab_c_melt$tree_id]
  ab_c_melt$year    = as.numeric(census_years[ab_c_melt$year_id])
  
  # aboveground biomass increment computation
  census_years = as.numeric(census_years)
  abi_c = t(apply(ab_cr, 1, function(x) diff(x)))
  for (i in 1:ncol(abi_c)){
    abi_c[,i] = abi_c[,i]/diff(census_years)[i]
  }
  
  # melt 
  abi_c_melt = melt(abi_c)
  colnames(abi_c_melt) = c('tree_id', 'year', 'abi')
  abi_c_melt$taxon   = taxon[abi_c_melt$tree_id]
  abi_c_melt$taxon   = taxaMatch[match(abi_c_melt$taxon, taxaMatch$taxon), 'species']
  abi_c_melt$site_id = tree_site_id[abi_c_melt$tree_id]
    
}

#############################################################################################################
### E. Save WIKI RDS and CSV file + clean up working directory ##############################################
#############################################################################################################

# clean up working directory because the next file is huge
variables = ls()
not_rm = variables %in% c('ab_p_melt','abi_p_melt','output_dir','site',
           'ab_m_melt','abi_m_melt','ab_c_melt','abi_c_melt',
           'years', 'ab_sum')
rm = variables[!not_rm]

rm(list=rm)

# MK ??: Do we need to save all of the NA observations in this data frame?
#preds = rbind(data.frame(ab_p_melt, type="AB"), data.frame(abi_p_melt, type="ABI"))

#write.csv(preds, file=paste0('output/', output_dir, '/NPP_STAT_MODEL_', site, '.csv'), row.names=FALSE)
#saveRDS(preds, paste0('output/', output_dir, '/NPP_STAT_MODEL_', site, '.RDS'))

#############################################################################################################
### F. Data Visualizations ##################################################################################
#############################################################################################################

###### Data manipulations for plotting #####

# for aboveground biomass iterations for stat model data
ab_p_sum_by_iter = aggregate(value ~ year+site_id+iter+model, ab_p_melt, function(x) sum(x, na.rm=TRUE))
ab_p_quants = aggregate(value~year+site_id+model, data=ab_p_sum_by_iter, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
ab_p_quants = data.frame(ab_p_quants)
ab_p_quants = cbind(ab_p_quants[,1:3], ab_p_quants[,4])
colnames(ab_p_quants)[4:6] = c('ab25', 'ab50', 'ab975') 

# for aboveground biomass measurements 
ab_m_sum_by_iter = aggregate(ab ~ year+site_id, ab_m_melt, function(x) sum(x, na.rm=TRUE))
ab_m_quants = aggregate(ab~year+site_id, data=ab_m_sum_by_iter, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
ab_m_quants = data.frame(ab_m_quants)
ab_m_quants = cbind(ab_m_quants[,1:2], ab_m_quants[,3])
colnames(ab_m_quants)[3:5] = c('ab25', 'ab50', 'ab975') 

# for aboveground biomass estimates from census data
colnames(ab_c_melt)[which(colnames(ab_c_melt)=="year_id")] = "year"
ab_c_sum_by_iter = aggregate(ab ~ year+site_id, ab_c_melt, function(x) sum(x, na.rm=TRUE))
ab_c_quants = aggregate(ab~year+site_id, data=ab_c_sum_by_iter, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
ab_c_quants = data.frame(ab_c_quants)
ab_c_quants = cbind(ab_c_quants[,1:2], ab_c_quants[,3])
colnames(ab_c_quants)[3:5] = c('ab25', 'ab50', 'ab975')  

# for aboveground biomass increment estimates for stat model results
abi_p_sum_by_iter = aggregate(value ~ year+site_id+iter+model, abi_p_melt, function(x) sum(x, na.rm=TRUE))
abi_p_quants = aggregate(value~year+site_id+model, data=abi_p_sum_by_iter, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
abi_p_quants = data.frame(abi_p_quants)
abi_p_quants = cbind(abi_p_quants[,1:3], abi_p_quants[,4])
colnames(abi_p_quants)[4:6] = c('ab25', 'ab50', 'ab975')  

# sums for census and measured
abi_c_sum = aggregate(abi ~ year+site_id, abi_c_melt, function(x) sum(x, na.rm=TRUE))
abi_m_sum = aggregate(abi ~ year+site_id, abi_m_melt, function(x) sum(x, na.rm=TRUE))

## MK ??: unsure why this is necessary 
# super hack. apologies...
abi_m_sum$year = abi_m_sum$year-1

##### Plot 1: AGB plots for four different methods ######

cols = c('Model RW + Census'='#084594', 'Model RW'='#8c2d04', 'Empirical RW'='#238b45', 'Empirical Census'='black')
cols_fill = c('Model RW + Census'="#4292c6", 'Model RW'="#fdd0a2", 'Empirical RW'="white", 'Empirical Census'="white")
ggplot() +  
  geom_ribbon(data=ab_p_quants, aes(x=year, ymin=ab25, ymax=ab975, fill=model, colour=model), alpha=0.4) +
  geom_line(data=ab_p_quants, aes(x=year, y=ab50, colour=model), size=1) + 
  geom_ribbon(data=ab_m_quants, aes(x=year, ymin=ab25, ymax=ab975, colour='Empirical RW', fill='Empirical RW'), alpha=0.4) +
  geom_line(data=ab_m_quants, aes(x=year, y=ab50, colour='Empirical RW'), size=1) + 
  geom_point(data=ab_c_quants, aes(x=year, y=ab50, colour='Empirical Census', fill = 'Empirical Census'), size=2) + 
  geom_errorbar(data=ab_c_quants, aes(x=year, ymin=ab25, ymax=ab975, colour='Empirical RW'), width=0.5) +
  facet_grid(site_id~.) + 
  scale_color_manual(values=cols, name='Method') +
  scale_fill_manual(values=cols_fill, name='Method') +
  scale_x_continuous(breaks=seq(min(years), max(years), by=5)) +
  theme_bw() + 
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14)) +
  ylab("Biomass (Mg/ha)") + 
  xlab('Year')

#ggsave(file=paste0('figures/AGB_by_site_', mvers, '.pdf'))
#ggsave(file=paste0('figures/AGB_by_site_', mvers, '.png'))

ggplot() +  
  geom_ribbon(data=ab_p_quants, aes(x=year, ymin=ab25, ymax=ab975, color = model, fill = model), alpha=0.4) +
  geom_line(data=ab_p_quants, aes(x=year, y=ab50, colour=model), size=1) + 
  geom_ribbon(data=ab_m_quants, aes(x=year, ymin=ab25, ymax=ab975, colour='Empirical RW', fill='Empirical RW'), alpha=0.4) +
  geom_line(data=ab_m_quants, aes(x=year, y=ab50, colour='Empirical RW'), size=1) + 
  geom_point(data=ab_c_quants, aes(x=year, y=ab50, colour='Empirical Census', fill='Empirical Census'), size=2) + 
  geom_errorbar(data=ab_c_quants, aes(x=year, ymin=ab25, ymax=ab975), width=0.5) +
  geom_line(data=ab_p_quants, aes(x=year, y=ab25, colour=model), linetype=1, size=0.5) +
  geom_line(data=ab_p_quants, aes(x=year, y=ab975, colour=model), linetype=1, size=0.5) +
  facet_grid(site_id~., scales="free_y") + scale_color_manual(values=cols, name='Method') +
  scale_fill_manual(values=cols_fill, name='Method') +
  theme_bw() + theme(axis.title=element_text(size=14), axis.text=element_text(size=14)) +
  ylab("Biomass (Mg/ha)") + xlab('Year') +
  scale_x_continuous(breaks=seq(min(years), max(years), by=5))

#ggsave(file=paste0('figures/', figures_dir, '/AGB_by_site_', site, '.pdf'))
#ggsave(file=paste0('figures/', figures_dir, '/AGB_by_site_', site, '.png'))


##### Plot 2: Mean Biomass #####

ggplot(data=ab_sum) + 
  geom_line(data=ab_sum, aes(x=year, y=ab, colour=model)) + 
  xlim(c(1960,2012)) + 
  facet_grid(plot~.)
#ggsave('figures/NOCOVAR/sum_ab_by_plot_HF_v4.0.pdf')

##### Plot 3: AGB increment plots for four different methods #####

ggplot() +  
  geom_ribbon(data=abi_p_quants, aes(x=year, ymin=ab25, ymax=ab975, fill=model, colour=model), alpha=0.4) +
  geom_line(data=abi_p_quants, aes(x=year, y=ab50, colour=model), size=1) + 
  geom_line(data=abi_m_sum, aes(x=year, y=abi, colour='Empirical RW', fill='Empirical RW'), size=1) + 
  geom_point(data=abi_c_sum, aes(x=year, y=abi, colour='Empirical Census', fill='Empirical Census'),size=2) + 
  facet_grid(site_id~.) +
  scale_color_manual(values=cols, name='Method') +
  scale_fill_manual(values=cols_fill, name='Method') + 
  theme_bw() +
  #ylim(-5,5) +
  ylab("Biomass Increment (Mg/ha/year)") + xlab('Year') +
  scale_x_continuous(breaks=seq(min(years), max(years), by=5), limits=c(min(years), max(years)))

#ggsave(file=paste0('figures/', figures_dir, '/AGBI_by_site_', site, '.pdf'))
#ggsave(file=paste0('figures/', figures_dir, '/AGBI_by_site_', site, '.png'))

#############################################################################################################
#############################################################################################################

# things that are never used: 

#ab_m_sum = aggregate(ab ~ year+site_id, ab_m_melt, function(x) sum(x, na.rm=TRUE))

#colnames(dbh_melt_all)[which(colnames(dbh_melt_all)=="plot")] = "site_id"
#colnames(dbh_melt_all)[which(colnames(dbh_melt_all)=="tree")] = "tree_id"
#colnames(dbh_melt_all)[which(colnames(dbh_melt_all)=="year")] = "year_id"

#ab_cr_median = apply(ab_cr, c(1,2), median, na.rm=TRUE)

#colnames(dbh_mean_by_iter)[which(colnames(dbh_mean_by_iter)=="plot")] = "site_id"
#fade_fix = merge(ab_p_sum_by_iter, dbh_mean_by_iter, by=c('year', 'site_id', 'model', 'iter'))
#saveRDS(fade_fix , file=paste0('allom/fading_record_correct_', site, '_', mvers, '.RDS'))

#fade_fix_by_tree = merge(ab_p_melt, dbh_melt_all, by=c('year_id', 'tree_id','site_id', 'model', 'iter'))
#saveRDS(fade_fix_by_tree, file=paste0('data/HF_DBH_AB_tree_iter.RDS'))
