#!/usr/bin/Rscript
library(plotrix)
library(dplR)
library(fields)
library(reshape2)
library(plyr)

# source("config_HMC")
dataDir = '~/Documents/projects/npp-stat-model/data'
it_wd = 'ROOSTER'
dvers = "v0.1"
mvers = "v0.1"

nPlots <- 2
ftPerMeter <- 3.2808399

lastYear  <- 2014
firstYear <- 1940
years <- firstYear:lastYear
nT <- length(years)

rwFiles <- list.files(paste0("data/", it_wd, '/', "rwl"))
rwFiles <- rwFiles[grep(".rwl$", rwFiles)]
rwData <- list()
for(fn in rwFiles) {
    id <- gsub(".rw", "", fn)
    rwData[[id]] <- t(read.tucson(file.path("data", it_wd, "rwl", fn)))  # rows are tree, cols are times
}

treeMeta = read.csv("data/ROOSTER/RoosterHillAllPlots.csv", skip=3)

incr=ldply(rwData, rbind)
rownames(incr) = as.vector(unlist(lapply(rwData, rownames)))
incr[,1] = rownames(incr)

######################################################################################################################################
## make nimble data
######################################################################################################################################
if (!file.exists('data/dump')){
  dir.create('data/dump')
} 
  
incr_data = melt(incr)
colnames(incr_data) = c('id', 'year', 'incr')
incr_data$plot   = as.numeric(substr(incr_data$id, 3, 3))
incr_data$TreeID = incr_data$id 
incr_data$id     = as.numeric(substr(incr_data$id, 4, 6))
incr_data$year = as.vector(incr_data$year)

# 
# # treeMetaNo   = treeMeta[match(incr_data$TreeID, treeMeta$TreeID), 'No.']
# # treeMetaPlot = treeMeta[match(incr_data$TreeID, treeMeta$TreeID), 'Plot']
# 
# year_start = firstYear
# year_end   = lastYear
# years      = seq(year_start, year_end)
# incr_data  = incr_data[which(incr_data$year %in% years),]
# trees_inc = sort(unique(incr_data$id))
# 
# # order by tree and year
# incr_data = incr_data[order(incr_data$id, incr_data$year),]
# 
# N_inc   = nrow(incr_data) # number of measurement 
# m2t     = incr_data$year
# m2tree  = incr_data$id
# m2treecode = incr_data$id
# m2ti    = match(m2t, years)
# Xobs    = incr_data$incr
# Xobs[Xobs==0] = 0.0001
# logXobs = log(Xobs)
# 
# # make pdbh
# pdbh = aggregate(year~id, incr_data, max)
# pdbh = pdbh[order(pdbh$id, pdbh$year),]
# pdbh$dbh = treeMeta$DBH_cm[match(pdbh$id, treeMeta$Tree)]
# 
# N_pdbh = nrow(pdbh)
# logPDobs = log(pdbh$dbh)
# pdbh_tree_id = pdbh$id
# # pdbh_day_id  = rep(100, nrow(pdbh)) / growDays #as.numeric(dbh$day) / growDays # FIXME
# pdbh_year_id = as.numeric(pdbh$year) - year_start + 1#getTimeIndex(dbh$year)
# pdbh_tree_code = pdbh$id
# 
# # ids_table = dbh[,c('stat_id', 'tree_id' )]
# # ids_table = ids_table[!duplicated(ids_table$stat_id),]
# 
# trees   = sort(unique(pdbh$id))
# N_trees = length(unique(pdbh$id))
# N_years = length(years)
# 
# last_time_data = vector(length=N_trees)
# last_time = vector(length=N_trees)
# last_time_pdbh = vector(length=N_trees)
# for (i in 1:N_trees){
#   tree = trees[i]
#   print(tree)
#   last_time[i] = max(c(pdbh$year[which(pdbh$id == tree)], incr_data$year[which(incr_data$id == tree)]), na.rm=TRUE)
#   last_time_data[i] = last_time[i]
# 
#   last_time_pdbh[i] = max(pdbh[pdbh$id == tree,'year'], last_time[i], na.rm=TRUE) 
# }
# 
# last_time = as.numeric(last_time)
# last_time_data = as.numeric(last_time_data)
# last_time_pdbh = as.numeric(last_time_pdbh)
# 
# 
# X_ord = data.frame(meas=numeric(0), tree_id=numeric(0), year=numeric(0))
# n = 1
# for (i in 1:N_trees){
#   print(i)
#   print(last_time[i])
#   #year = seq(year_start, last_time[i])
#   year = seq(year_start, last_time_pdbh[i])
#   meas = seq(n, n+length(year)-1)
#   n = n + length(year)
#   
#   X_ord = rbind(X_ord, data.frame(meas=meas, tree_id=rep(trees[i], length(year)), year=year))
# }
# 
# x2tree  = X_ord$tree_id
# x2year  = match(X_ord$year, years) 
# # last_ti = last_time-year_start +1
# last_ti = last_time_pdbh-year_start +1
# 
# N_X   = nrow(X_ord)
# N_D   = N_X
# 
# meas2x = vector(length=N_inc)
# for (i in 1:N_inc) {
#   id = incr_data$id[i]
#   year    = incr_data$year[i]
#   
#   meas2x[i] = which((X_ord$tree_id == id) & (X_ord$year == year))
# }
# 
# pdbh2d = vector(length=N_pdbh)
# for (i in 1:N_pdbh){
#   id = pdbh$id[i]
#   year = pdbh$year[i]
#   
#   print(i)
#   which((X_ord$tree_id == id) & (X_ord$year == year))
#   
#   pdbh2d[i] = which((X_ord$tree_id == id) & (X_ord$year == year))
# }
# 
# ones = rep(1, 4)
# N_taxa=length(unique(treeMeta$SpecCode))
# 
# cs_last_ti = cumsum(last_ti)
# first_ti   = c(1,cs_last_ti[1:(length(cs_last_ti)-1)]+1)
# 
# X = rep(0.1, N_X)
# D = rep(0.1, N_X)
# logX = rep(log(0.1), N_X)
# D0 = rep(3, N_trees)
# 
# beta = rep(0.1, N_trees)
# beta_t = rep(0.1, N_years)
# beta0 = 0.1
# sig_x_obs = 0.7
# sig_d_obs = 0.01
# sig_d = 0.1
# sig_d_sap = 0.1
# sig_x = 0.1
# beta_sd = 0.1
# beta_t_sd = 0.1
# beta_spp_sd = 0.1
# beta_spp = rep(0.1, N_taxa)
# beta_slope = 0.1
# 
# b0 = 0.1
# b1 = 10
# 
# nu = 0.1
# rho = 0.1
# 
# taxon = treeMeta$SpecCode
# 
# tbl <- table(taxon)
# taxaMatch <- data.frame(tbl); names(taxaMatch) <- c('species', 'count')
# taxaMatch$taxon <- 1:nrow(taxaMatch)
# 
# dump(c('N_inc', 'N_pdbh', 'N_X', 'N_D', 'N_trees', 'N_years',
#        'logXobs', 'm2t', 'm2ti', 'm2tree', 
#        'logPDobs', 'pdbh_tree_id', 'pdbh_year_id', 
#        'meas2x', 'x2year', 'x2tree', 
#        'pdbh2d',
#        'last_ti',
#        'ones',
#        'year_start', 'year_end',
#        'taxon', 'N_taxa', 
#        # 'N_saplings', 'sap2x', 'not_sap2x','max_size', 'sapling_tree_id', 'sapling_year_id',
#        'first_ti', 'cs_last_ti'),
#      file=paste0('data/dump/tree_data_IT_', dvers, '.dump'))
# 
# save(N_inc, N_pdbh, N_X, N_D, N_trees, N_years, 
#      logXobs, m2t, m2ti, m2tree,
#      logPDobs, pdbh_tree_id, pdbh_year_id, 
#      meas2x, x2year, x2tree,
#      pdbh2d,
#      last_ti, last_time, last_time_data,
#      ones,
#      year_start, year_end,
#      trees, years, 
#      taxon, N_taxa,
#      first_ti, cs_last_ti, pdbh, taxaMatch,
#      file=paste0('data/dump/tree_data_IT_', dvers, '.rdata'))


##########################################################################################################################
## STAN DATA
##########################################################################################################################
incr = incr[,-1]

year_end = max(as.numeric(incr_data$year), na.rm=TRUE)
year_start = min(as.numeric(incr_data$year), na.rm=TRUE)
N_years = year_end - year_start + 1
years = seq(year_start, year_end)

incr_data = incr_data[which(!is.na(incr_data$incr)),]
year_idx = data.frame(as.numeric(aggregate(year~id, incr_data, min)[,2]), as.numeric(aggregate(year~id, incr_data, max)[,2]))
year_idx[,2] = rep(2016, nrow(year_idx))
year_idx = year_idx - year_start + 1

incr[is.na(incr)] = -999
logXobs=incr[order(rownames(incr)),]


# make pdbh
pdbh = aggregate(year~id+plot, incr_data, max)
pdbh = pdbh[order(pdbh$id, pdbh$year),]

for (n in 1:nrow(pdbh)){
  pdbh$dbh[n] = treeMeta$DBH[which((as.numeric(substr(treeMeta$Site, 3, 3))==pdbh$plot[n])&(treeMeta$Tree.Number == pdbh$id[n]))]
}
N_pdbh = nrow(pdbh)
logPDobs = log(pdbh$dbh)
pdbh_tree_id = pdbh$id
logPDobs[is.na(logPDobs)] = -999



pdbh_year_id = rep(N_years, N_trees)

idx_stack = data.frame(meas=numeric(0), tree_id=numeric(0), year=numeric(0))
n = 1
for (tree in 1:N_trees){
  year = seq(year_idx[tree,1], year_idx[tree,2])
  meas = seq(n, n+length(year)-1)
  n = n + length(year)
  idx_stack = rbind(idx_stack, data.frame(meas=meas, tree_id=rep(tree, length(year)), year=year))
}

idx_tree = which(!duplicated(idx_stack$tree_id))
idx_tree = data.frame(idx_tree, c(idx_tree[-1]-1, nrow(idx_stack)))

x2tree  = idx_stack$tree_id
x2year  = idx_stack$year 

N_vals   = nrow(idx_stack)

# meas2x = vector(length=N_vals)
# for (i in 1:N_inc) {
#   id = incr_data$id[i]
#   year    = incr_data$year[i]
#   
#   meas2x[i] = which((X_ord$tree_id == id) & (X_ord$year == year))
# }

pdbh$year = rep(N_years, nrow(pdbh))

pdbh2val = vector(length=N_pdbh)
for (i in 1:N_pdbh){
  id = pdbh$id[i]
  year = pdbh$year[i]
  
  print(i)
  which((idx_stack$tree_id == id) & (idx_stack$year == year))
  
  pdbh2val[i] = which((idx_stack$tree_id == id) & (idx_stack$year == year))
}


saveRDS(list(N_trees=N_trees, 
           N_years=N_years,
           N_vals=N_vals,
           logXobs=logXobs, 
           logPDobs=logPDobs,
           year_idx=year_idx, 
           taxon=taxon, 
           N_taxa=N_taxa,
           pdbh_year=pdbh_year_id,
           idx_tree =idx_tree, 
           pdbh2val=pdbh2val,
           x2tree=x2tree,
           x2year=x2year,
           taxaMatch=taxaMatch,
           years = years),
           file=paste0('data/dump/tree_data_IT_STAN_', dvers, '.RDS'))
