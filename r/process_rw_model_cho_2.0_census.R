# THIS SCRIPT SHOULD BE USED FOR SITES WITH AVAILABLE CENSUS DATA # 

rm(list=ls())
setwd('~/Desktop/npp-stat-model')

library(ggplot2)
library(reshape2)
library(abind)
library(dplyr)

#############################################################################################################
### A. Set up working environment and directory for script ##################################################
#############################################################################################################

########################################
# Variables to be adjusted for each site
########################################

# date for save purposes
date = '27Dec2019' 

# set up output 
site <- 'HARVARD'

# data version
dvers = 'v5.0'
mvers = 'v4.0'

# burn-in range
burn = 1000

# load model output data for site
data_dir = file.path('data','sites',site)

# get ring width data
models = c('Model RW + Census', 'Model RW')
fnames = c('ring_model_t_date_sapl_size_pdbh_NOCOVAR_v4.0.Rdata', 'ring_model_t_pdbh_nc_NOCOVAR_sigd_v4.0.Rdata')

# load tree data 
load(file.path(data_dir,'tree_data_20_no_census_v5.0.rdata'))
load(file.path(data_dir,'tree_data_20_v5.0.rdata'))

# are dbh and ab files already created? do we want to recreate them? 
run_ringwidth = TRUE

# the year ring widths were collected (this year and next years are removed)
rw_year = 2013

########################################

# create output directory
output_dir <- file.path('sites',site)
if (!file.exists(output_dir)){
  dir.create(output_dir)
  dir.create(file.path(output_dir, 'figures'))
}

# get ring width data
nmodels = length(fnames)
post = list()
for (i in 1:length(fnames)) {
  fname_model = fnames[i]
  load(file   = paste0(data_dir,'/', fname_model))
  post[[i]]   = out[1:2000,]
}  
niter = dim(post[[1]])[1]

# how many iterations to keep?
keep = niter - burn
keep = 5

# match species acronyms to level3a available species/pfts 
acro_level3a = read.csv('data/acronym_to_level3a_v0.1.csv', stringsAsFactors = FALSE)
acro_level3a = left_join(taxaMatch, acro_level3a, by = c('species'='acronym'))
choj = read.csv('data/level3a_to_chojnacky_v0.4.csv', stringsAsFactors = F)
choj = left_join(acro_level3a, choj, id = 'level3a') 

#############################################################################################################
### B. Process diameter and ab predictions for ring width data ##############################################
#############################################################################################################

##### Step 1: Estimate diameter and aboveground biomass for all trees from ring width data #####

# this function takes the stat model results (dbh, rw) and returns a list of dbh arrays (tree X years X iterations)
# for each pft in the data, a list of the tree indices in the order they were used in the list, and the order of the
# pfts within the list 
build_dbh_p <- function(out, rw2year, rw2tree, Nyears, tree_list, keep, taxaMatch, taxon_list){
  
  # extract colnames and number of iterations (different value types in data, want diameter)
  col_names = sapply(strsplit(colnames(out), '\\['), function(x) x[[1]])
  niter   = dim(out)[1]
  
  # get species from taxon list 
  pft = as.vector(taxaMatch[match(taxon_list, taxaMatch$taxon),'species'])
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
if (run_ringwidth){
  
  # for first model with census 
  # obtain diameter estimates for all trees for
  dbh_1 = build_dbh_p(post[[1]], x2year, x2tree, N_years, trees, keep, taxaMatch, taxon)
  dbh_p_1 = dbh_1$dbh_p
  tree_order_1 = dbh_1$tree_order
  pft_order_1  = dbh_1$pft_order
  
  # make list of arrays into one large array across species (all trees X years X iterations)
  dbh_p_1_org = abind(dbh_p_1, along=1)
  dbh_p_1 = dbh_p_1_org[sort(tree_order_1, index.return=TRUE)$ix,,]
  
  # get choj pft index for coefficients in matrix 
  choj_pfts_1 = taxon
  
  # find biomass for all trees, years, and iterationss
  ab_p_1 = build_ab_choj(dbh_p_1, N_trees, N_years, keep, choj_pfts_1)

  # 2: for model without census
  # organize tree info for data with no census 
  trees_p = sort(unique(x2tree_p)) # get all trees ids in data
  taxon_p = taxon[trees_p] # gets taxons of unique trees
  x2idx_p = match(x2tree_p, unique(x2tree_p)) # matches all trees from data to a numerical index 
  pft_list_p = taxaMatch$species[taxon_p] # matches all trees to a pft
  
  # get diameter estimates from ring widths given by stat model 
  dbh_2 = build_dbh_p(post[[2]], x2year_p, x2tree_p, N_years, trees_p, keep, taxaMatch, taxon_p)
  dbh_p_2 = dbh_2$dbh_p
  tree_order_2 = dbh_2$tree_order
  pft_order_2  = dbh_2$pft_order
  
  # make list of arrays into one big array (all trees X years X iterations)
  dbh_p_2_org = abind(dbh_p_2, along=1)
  dbh_p_2 = dbh_p_2_org[sort(tree_order_2, index.return=TRUE)$ix,,]
  
  # get choj pft index for coefficients in matrix 
  choj_pfts_2 = taxon_p
  
  # find biomasss for all trees, years, and iterations
  ab_p_2 = build_ab_choj(dbh_p_2, N_trees_p, N_years, keep, choj_pfts_2)
}


# added separate save variable in case we need to go back to this step
dbh_p_1_raw = dbh_p_1
ab_p_1_raw = ab_p_1
dbh_p_2_raw = dbh_p_2
ab_p_2_raw = ab_p_2

##### Step 2: Remove the data points for trees with diameters less than 5 cm #####

# remove trees that have a mean diameter less than 5 centimeters across all iterations for a given year for both models
for (tree in 1:N_trees){
  for (year in 1:N_years){
    dbh_mean = mean(dbh_p_1[tree, year, ], na.rm=TRUE)
    if (is.na(dbh_mean)){next}
    if (dbh_mean < 5){
      dbh_p_1[tree, year, ] = rep(NA, keep)
      ab_p_1[tree, year, ] = rep(NA, keep)
    }
  }
}

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
  
# added save variable in case we need to back to this step
dbh_p_1_red = dbh_p_1
ab_p_1_red = ab_p_1
dbh_p_2_red = dbh_p_2
ab_p_2_red = ab_p_2

##### Step 3: Adjust for mortality using census data #####

# get the census years
cyr = dbh$yr
cyr[which(cyr %in% c('01', '11'))] = paste0('20', cyr[which(cyr %in% c('01', '11'))])
cyr[which(cyr %in% c('62', '69', '75', '91'))] = paste0('19', cyr[which(cyr %in%  c('62', '69', '75', '91'))])
census_years = as.numeric(sort(unique(cyr)))
idx_census   = which(years %in% census_years)
N_census_years = length(census_years)
dbh = data.frame(dbh, cyr=cyr)

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

##### Step 4: Melt data matrices into data frames for easier use #####

# melt down to data frames
rownames(dbh_p_1) <- seq(1, N_trees)
dbh_melt = melt(dbh_p_1)
colnames(dbh_melt) = c("tree", "year", "iter", "dbh")
dbh_melt$year = years[dbh_melt$year]
dbh_melt$plot = tree_site_id[dbh_melt$tree]

rownames(ab_p_1) <- seq(1, N_trees)
ab_melt = melt(ab_p_1)
colnames(ab_melt) = c("tree", "year", "iter", "ab")
ab_melt$year = years[ab_melt$year]
ab_melt$plot = tree_site_id[ab_melt$tree]

rownames(dbh_p_2) <- trees_p
dbh2_melt = melt(dbh_p_2)
colnames(dbh2_melt) = c("tree", "year", "iter", "dbh")
dbh2_melt$year = years[dbh2_melt$year]
dbh2_melt$plot = tree_site_id[dbh2_melt$tree]

rownames(ab_p_2) <- trees_p
ab2_melt = melt(ab_p_2)
colnames(ab2_melt) = c("tree", "year", "iter", "ab")
ab2_melt$year = years[ab2_melt$year]
ab2_melt$plot = tree_site_id[ab2_melt$tree]

# put models together
dbh_melt_all = rbind(data.frame(dbh_melt, model=rep("Model RW + Census", nrow(dbh_melt))), 
                     data.frame(dbh2_melt, model=rep("Model RW", nrow(dbh2_melt))))
dbh_melt_all$plot = tree_site_id[dbh_melt_all$tree]
dbh_melt_all$census_id = dbh[match(dbh_melt_all$tree, dbh$stat_id),'census_id']
dbh_melt_all$dist_census = census[match(dbh_melt_all$census_id, census$census_id),'dist_census']/3.28084

ab_melt_all = rbind(data.frame(ab_melt, model=rep("Model RW + Census", nrow(ab_melt))), 
                    data.frame(ab2_melt, model=rep("Model RW", nrow(ab2_melt))))
ab_melt_all$plot = tree_site_id[ab_melt_all$tree]

# aggregate for some plots later
dbh_mean = aggregate(dbh~tree+year+plot, dbh_melt, mean, na.rm=TRUE) 
dbh_sum = aggregate(dbh~year+plot, dbh_mean, sum, na.rm=TRUE)
ab_mean = aggregate(ab~tree+plot+year+model, ab_melt_all, mean) 
ab_sum = aggregate(ab~year+model+plot, ab_mean, sum)

# save to new variable in case we need to go back to this step
dbh_p_1_smth = dbh_p_1
ab_p_1_smth = ab_p_1
dbh_p_2_smth = dbh_p_2
ab_p_2_smth = ab_p_2


##### Step 5: Calculate biomass and biomass increment #####

# MODEL 1 
# in Kg/plot, rescale so it is Mg/ha
ab_pr1 = ab_p_1/(20^2*pi) * 10

# melt predicted biomass for both models
ab_p1_melt = melt(ab_pr1)
colnames(ab_p1_melt) = c('tree_id', 'year_id', 'iter', 'ab')
ab_p1_melt$plot = tree_site_id[ab_p1_melt$tree_id]

# MODEL 2 
# in Kg/plot, rescale so it is Mg/ha
ab_pr2 = ab_p_2/(20^2*pi) * 10

# melt predicted biomass matrix
ab_p2_melt = melt(ab_pr2)
colnames(ab_p2_melt) = c('tree_id', 'year_id', 'iter', 'ab')
ab_p2_melt$plot = tree_site_id[ab_p2_melt$tree_id]

# PUT THEM TOGETHER
ab_p_melt = rbind(data.frame(ab_p1_melt, model=rep('Model RW + Census')), 
                  data.frame(ab_p2_melt, model=rep('Model RW')))
ab_p_melt$taxon   = taxon[ab_p_melt$tree_id]
ab_p_melt$taxon   = taxaMatch[match(ab_p_melt$taxon, taxaMatch$taxon), 'species']
ab_p_melt$year = years[ab_p_melt$year_id]

# determine biomass increment for model 1 
abi_p1 = apply(ab_pr1, c(1,3), function(x) diff(x))
abi_p1 = aperm(abi_p1, c(2, 1, 3))
abi_p1_melt = melt(abi_p1)
colnames(abi_p1_melt) = c('tree_id', 'year_id', 'iter', 'abi')
  
# determine biomass increment for model 2
abi_p2 = apply(ab_pr2, c(1,3), function(x) diff(x))
abi_p2 = aperm(abi_p2, c(2, 1, 3))
abi_p2_melt = melt(abi_p2)
colnames(abi_p2_melt) = c('tree_id', 'year_id', 'iter', 'abi')

# put them together
abi_p_melt = rbind(data.frame(abi_p1_melt, model=rep('Model RW + Census')), 
                   data.frame(abi_p2_melt, model=rep('Model RW')))
abi_p_melt$taxon   = taxon[abi_p_melt$tree_id]
abi_p_melt$taxon   = taxaMatch[match(abi_p_melt$taxon, taxaMatch$taxon), 'species']
abi_p_melt$plot = tree_site_id[abi_p_melt$tree_id]
abi_p_melt$year    = years[abi_p_melt$year_id]

##### Step 6: Prep AGB and ABI for wiki #####

abi_p_melt = abi_p_melt[,c('tree_id', 'year', 'plot', 'taxon', 'model', 'iter', 'abi')]
ab_p_melt = ab_p_melt[,c('tree_id', 'year', 'plot', 'taxon', 'model', 'iter', 'ab')]

ab_p_melt = ab_p_melt[ab_p_melt$year<rw_year,]
abi_p_melt = abi_p_melt[abi_p_melt$year<rw_year,]

colnames(ab_p_melt)[which(colnames(ab_p_melt) == "ab")] = 'value'
colnames(abi_p_melt)[which(colnames(abi_p_melt) == "abi")] = 'value'

#############################################################################################################
### C. Process dbh and ab from most recent dbh measurement ###################################
#############################################################################################################

##### Step 1: Track diameter back from most recent dbh measurement using increment data #####

dbh_m = array(NA, c(N_trees, N_years))

# loop through all trees 
for (i in 1:N_trees) {
  print(i)
  tree = trees[i]
  pft  = as.vector(taxaMatch[which(taxaMatch$taxon == taxon[tree]),'species'])
  
  # if tree doesn't have increment data, make NA for all years 
  if (length(which(m2tree_a == tree)) == 0) {
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
    
    # census allows for discrepancies in data products 
    # tree was dead and then re-measured as alive
    # 10 cm in census, 5 cm in next coring
    # tree 139
    if (tree==139){
      dbh_year = 2009
      dbh_tree = 46.3
    } 
    # tree 145
    if (tree == 145) {
      dbh_year = 2007
      dbh_tree = 13.3
    }
    
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

# ln(biomass-kg) = beta0 + beta1 * ln(diameter-cm)
b0 = choj$beta0[taxon]
b1 = choj$beta1[taxon]
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

##### Step 3: Melt down and plot #####

# melt to data frames
rownames(dbh_m) <- seq(1, N_trees)
dbh_m_melt = melt(dbh_m)
colnames(dbh_m_melt) = c("tree", "year", "dbh")
dbh_m_melt$year = years[dbh_m_melt$year]
dbh_m_melt$plot = tree_site_id[dbh_m_melt$tree]

# in Kg/plot, rescale so it is Mg/ha
ab_mr = ab_m/(20^2*pi) * 10

# melt measured
ab_m_melt = melt(ab_mr)
colnames(ab_m_melt) = c('tree_id', 'year_id', 'ab')
ab_m_melt$taxon   = taxon[ab_m_melt$tree_id]
ab_m_melt$taxon   = taxaMatch[match(ab_m_melt$taxon, taxaMatch$taxon), 'species']
ab_m_melt$plot = tree_site_id[ab_m_melt$tree_id]
ab_m_melt$year    = years[ab_m_melt$year_id]
ab_m_melt = ab_m_melt[ab_m_melt$year<rw_year,]

# measured biomass increment
abi_m = t(apply(ab_mr, 1, function(x) diff(x)))

# melt measured
abi_m_melt = melt(abi_m)
colnames(abi_m_melt) = c('tree_id', 'year_id', 'abi')
abi_m_melt$taxon   = taxon[abi_m_melt$tree_id]
abi_m_melt$taxon   = taxaMatch[match(abi_m_melt$taxon, taxaMatch$taxon), 'species']
abi_m_melt$plot = tree_site_id[abi_m_melt$tree_id]
abi_m_melt$year    = years[abi_m_melt$year_id]
abi_m_melt = abi_m_melt[abi_m_melt$year<rw_year,]

dbh_m_sum = aggregate(dbh~year+plot, dbh_m_melt, sum)

ggplot() + 
  geom_line(data=dbh_sum, aes(x=year, y=dbh, colour=factor(plot), group=factor(plot)), linetype = "dashed") +
  geom_line(data=dbh_m_sum, aes(x=year, y=dbh, colour=factor(plot), group=factor(plot))) + 
  labs(title ='Sum of DBH Values for DBH and RW Obs', colour = 'Plot ID')
ggsave(paste0(output_dir, '/figures/sum_dbh_by_plot_', site,'.pdf'))

#############################################################################################################
### D. Process dbh and ab from census data ##################################################################
#############################################################################################################

##### Step 1: Organize census data #####
# get census years
cyr = dbh$yr
cyr[which(cyr %in% c('01', '11'))] = paste0('20', cyr[which(cyr %in% c('01', '11'))])
cyr[which(cyr %in% c('62', '69', '75', '91'))] = paste0('19', cyr[which(cyr %in%  c('62', '69', '75', '91'))])
census_years = sort(unique(cyr))

dbh = data.frame(dbh, cyr=cyr)
N_census_years = length(census_years)

# organize data using dcast 
dbh_c_cast = dcast(dbh, stat_id+site~cyr, value.var=c('value'))

##### Step 2: Get biomass estimates based on Chojnacky coefficients #####

dbh_c = dbh_c_cast[,3:ncol(dbh_c_cast)]

# ln(biomass-kg) = beta0 + beta1 * ln(diameter-cm)
b0 = choj$beta0[taxon]
b1 = choj$beta1[taxon]
ab_c = exp(b0 + b1*log(dbh_c))

# in Kg/plot, rescale so it is Mg/ha
ab_cr = as.matrix(ab_c/(20^2*pi) * 10)

# melt census
ab_c_melt = melt(ab_cr)
colnames(ab_c_melt) = c('tree_id', 'year', 'ab')
ab_c_melt$taxon   = taxon[ab_c_melt$tree_id]
ab_c_melt$taxon   = taxaMatch[match(ab_c_melt$taxon, taxaMatch$taxon), 'species']
ab_c_melt$plot = tree_site_id[ab_c_melt$tree_id]

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
abi_c_melt$plot = tree_site_id[abi_c_melt$tree_id]

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

# get rid of NA values and save to one large data frame
ab_p_melt = ab_p_melt %>% filter(!is.na(value))
abi_p_melt = abi_p_melt %>% filter(!is.na(value))
preds = rbind(data.frame(ab_p_melt, type="AB"), data.frame(abi_p_melt, type="ABI"))

write.csv(preds, file=paste0(output_dir, '/output/','NPP_STAT_MODEL_', site, '.csv'), row.names=FALSE)
saveRDS(preds, paste0(output_dir, '/output/','NPP_STAT_MODEL_', site, '.RDS'))

#############################################################################################################
### F. Data Visualizations ##################################################################################
#############################################################################################################

###### Data manipulations for plotting #####

# for aboveground biomass iterations for stat model data
ab_p_sum_by_iter = aggregate(value ~ year+iter+plot+model, ab_p_melt, function(x) sum(x, na.rm=TRUE))
ab_p_quants = aggregate(value~year+plot+model, data=ab_p_sum_by_iter, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
ab_p_quants = data.frame(ab_p_quants)
ab_p_quants = cbind(ab_p_quants[,1:3], ab_p_quants[,4])
colnames(ab_p_quants) = c('year','plot','model','ab25', 'ab50', 'ab975') 

# for aboveground biomass measurements 
ab_m_sum_by_iter = aggregate(ab ~ year+plot, ab_m_melt, function(x) sum(x, na.rm=TRUE))
ab_m_quants = aggregate(ab~year+plot, data=ab_m_sum_by_iter, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
ab_m_quants = data.frame(ab_m_quants)
ab_m_quants = cbind(ab_m_quants[,1:2], ab_m_quants[,3])
colnames(ab_m_quants) = c('year','plot','ab25', 'ab50', 'ab975') 

# for aboveground biomass estimates from census data
ab_c_sum_by_iter = aggregate(ab ~ year+plot, ab_c_melt, function(x) sum(x, na.rm=TRUE))
ab_c_quants = aggregate(ab~year+plot, data=ab_c_sum_by_iter, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
ab_c_quants = data.frame(ab_c_quants)
ab_c_quants = cbind(ab_c_quants[,1:2], ab_c_quants[,3])
colnames(ab_c_quants) = c('year','plot','ab25', 'ab50', 'ab975')  

# for aboveground biomass increment estimates for stat model results
abi_p_sum_by_iter = aggregate(value ~ year+plot+iter+model, abi_p_melt, function(x) sum(x, na.rm=TRUE))
abi_p_quants = aggregate(value~year+plot+model, data=abi_p_sum_by_iter, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
abi_p_quants = data.frame(abi_p_quants)
abi_p_quants = cbind(abi_p_quants[,1:3], abi_p_quants[,4])
colnames(abi_p_quants) = c('year','plot','model','ab25', 'ab50', 'ab975')  

# sums for census and measured
abi_c_sum = aggregate(abi ~ year+plot, abi_c_melt, function(x) sum(x, na.rm=TRUE))
abi_m_sum = aggregate(abi ~ year+plot, abi_m_melt, function(x) sum(x, na.rm=TRUE))

# MK :: super hack. apologies...
#abi_m_sum$year = abi_m_sum$year-1

##### Plot 1: AGB plots for four different methods ######

cols = c('Model RW + Census'='#084594', 'Model RW'='#8c2d04', 'Empirical RW'='#238b45', 'Empirical Census'='black')
cols_fill = c('Model RW + Census'="#4292c6", 'Model RW'="#fdd0a2", 'Empirical RW'="white", 'Empirical Census'="white")

# same y scales 
ggplot() +  
  geom_ribbon(data=ab_p_quants, aes(x=year, ymin=ab25, ymax=ab975, fill=model, colour=model), alpha=0.4) +
  geom_line(data=ab_p_quants, aes(x=year, y=ab50, colour=model), size=1) + 
  geom_ribbon(data=ab_m_quants, aes(x=year, ymin=ab25, ymax=ab975, colour='Empirical RW', fill='Empirical RW'), alpha=0.4) +
  geom_line(data=ab_m_quants, aes(x=year, y=ab50, colour='Empirical RW'), size=1) + 
  geom_point(data=ab_c_quants, aes(x=year, y=ab50, colour='Empirical Census', fill = 'Empirical Census'), size=2) + 
  geom_errorbar(data=ab_c_quants, aes(x=year, ymin=ab25, ymax=ab975, colour='Empirical RW'), width=0.5) +
  facet_grid(plot~.) + 
  scale_color_manual(values=cols, name='Method') +
  scale_fill_manual(values=cols_fill, name='Method') +
  scale_x_continuous(breaks=seq(min(years), max(years), by=5)) +
  theme_bw() + 
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14)) +
  ylab("Biomass (Mg/ha)") + 
  xlab('Year')
ggsave(file=paste0(output_dir,'/figures/AGB_by_plot_',site,'.png'))

# free y scales 
ggplot() +  
  geom_ribbon(data=ab_p_quants, aes(x=year, ymin=ab25, ymax=ab975, color = model, fill = model), alpha=0.4) +
  geom_line(data=ab_p_quants, aes(x=year, y=ab50, colour=model), size=1) + 
  geom_ribbon(data=ab_m_quants, aes(x=year, ymin=ab25, ymax=ab975, colour='Empirical RW', fill='Empirical RW'), alpha=0.4) +
  geom_line(data=ab_m_quants, aes(x=year, y=ab50, colour='Empirical RW'), size=1) + 
  geom_point(data=ab_c_quants, aes(x=year, y=ab50, colour='Empirical Census', fill='Empirical Census'), size=2) + 
  geom_errorbar(data=ab_c_quants, aes(x=year, ymin=ab25, ymax=ab975), width=0.5) +
  geom_line(data=ab_p_quants, aes(x=year, y=ab25, colour=model), linetype=1, size=0.5) +
  geom_line(data=ab_p_quants, aes(x=year, y=ab975, colour=model), linetype=1, size=0.5) +
  facet_grid(plot~., scales="free_y") + scale_color_manual(values=cols, name='Method') +
  scale_fill_manual(values=cols_fill, name='Method') +
  theme_bw() + theme(axis.title=element_text(size=14), axis.text=element_text(size=14)) +
  ylab("Biomass (Mg/ha)") + xlab('Year') +
  scale_x_continuous(breaks=seq(min(years), max(years), by=5))
ggsave(file=paste0(output_dir,'/figures/AGB_by_plot_',site,'_freey.png'))

##### Plot 2: AGB increment plots for four different methods #####

ggplot() +  
  geom_ribbon(data=abi_p_quants, aes(x=year, ymin=ab25, ymax=ab975, fill=model, colour=model), alpha=0.4) +
  geom_line(data=abi_p_quants, aes(x=year, y=ab50, colour=model), size=1) + 
  geom_line(data=abi_m_sum, aes(x=year, y=abi, colour='Empirical RW', fill='Empirical RW'), size=1) + 
  geom_point(data=abi_c_sum, aes(x=year, y=abi, colour='Empirical Census', fill='Empirical Census'),size=2) + 
  facet_grid(plot~.) +
  scale_color_manual(values=cols, name='Method') +
  scale_fill_manual(values=cols_fill, name='Method') + 
  theme_bw() +
  ylab("Biomass Increment (Mg/ha/year)") + xlab('Year') +
  scale_x_continuous(breaks=seq(min(years), max(years), by=5), limits=c(min(years), max(years)))
ggsave(file=paste0(output_dir,'/figures/AGBI_by_plot_',site,'.png'))


#############################################################################################################
#############################################################################################################

# things that are never used: 

#ab_m_sum = aggregate(ab ~ year+plot, ab_m_melt, function(x) sum(x, na.rm=TRUE))

#colnames(dbh_melt_all)[which(colnames(dbh_melt_all)=="plot")] = "plot"
#colnames(dbh_melt_all)[which(colnames(dbh_melt_all)=="tree")] = "tree_id"
#colnames(dbh_melt_all)[which(colnames(dbh_melt_all)=="year")] = "year_id"

#ab_cr_median = apply(ab_cr, c(1,2), median, na.rm=TRUE)

#colnames(dbh_mean_by_iter)[which(colnames(dbh_mean_by_iter)=="plot")] = "plot"
#fade_fix = merge(ab_p_sum_by_iter, dbh_mean_by_iter, by=c('year', 'plot', 'model', 'iter'))
#saveRDS(fade_fix , file=paste0('allom/fading_record_correct_', site, '_', mvers, '.RDS'))

#fade_fix_by_tree = merge(ab_p_melt, dbh_melt_all, by=c('year_id', 'tree_id','plot', 'model', 'iter'))
#saveRDS(fade_fix_by_tree, file=paste0('data/HF_DBH_AB_tree_iter.RDS'))

# 
# # MK ??: Do species always follow this same order? Or is it through the order given in taxa match?
# # pft list 
# pfts = list(ACRU = data.frame(spcd=316,acronym="ACRU"),
#             ACSA = data.frame(spcd=318,acronym="ACSA3"),
#             BEAL = data.frame(spcd=371,acronym="BEAL"),
#             BELE = data.frame(spcd=372),
#             BEPO = data.frame(spcd=379,acronym="BEPO"),
#             CADE = data.frame(spcd=1000,acronym="CADE"), # id is 421, but lump into mixed hardwoods because no biomass eq
#             FAGR = data.frame(spcd=531,acronym="FAGR"),
#             FRAM = data.frame(spcd=541, acronym="FRAM"),
#             HAVI = data.frame(spcd=1000,acronym="HAVI"), # id is 585, but lump into mixed hardwoods because no biomass eq
#             PIST = data.frame(spcd=129,acronym="PIST"),
#             PRSE = data.frame(spcd=762,acronym="PRSE"),
#             QURU = data.frame(spcd=833,acronym="QURU"),
#             QUVE = data.frame(spcd=837,acronym="QUVE"),
#             QUAL = data.frame(spcd=802, acronym="QUAL"),
#             TSCA = data.frame(spcd=261,acronym="TSCA")
# )
# pft_list  = sort(unique(names(pfts)))

