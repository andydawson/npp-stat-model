#!/usr/bin/Rscript
library(plotrix)
library(dplR)
library(fields)
library(reshape2)

source("config_HMC")
dataDir = '~/Documents/projects/npp-stat-model/data'

plot_ids = data.frame(plot = seq(1,4), census_plot = c(7094, 7092, 7095, 7093))

nPlots <- 2

ftPerMeter <- 3.2808399

lastYear  <- 2015
firstYear <- 1960#1905
years <- firstYear:lastYear
nT <- length(years)

growDays <- 30+31+30+31+31+30+31

huron_wd <- paste0(dataDir, '/huron')


# encapsulate firstYear in getTimeIndex closure
getTimeIndex_gen <- function(firstYear) {
    function(year) {
        year - firstYear + 1
    }
}
getTimeIndex <- getTimeIndex_gen(firstYear)


## read data in
censusFile = 'HMC_census.csv'
treeMetaFile = 'HMC.field.data.csv'

censusFull <- read.csv(paste0(huron_wd, '/', censusFile), stringsAsFactors = FALSE)
treeMeta <- read.csv(paste0(huron_wd, '/', treeMetaFile), skip = 0, stringsAsFactors = FALSE)

# ringMeta <- read.csv(paste0(lyford_wd, '/', ringMetaFile), stringsAsFactors = FALSE)
# censusDates <- read.csv(paste0(lyford_wd, '/', censusSampleDatesFile), stringsAsFactors = FALSE)
# rwSampleDates <- read.csv(paste0(lyford_wd, '/', rwPlotSampleDatesFile), stringsAsFactors = FALSE)

rwFiles <- list.files(paste0(huron_wd, '/', "RingWidths"))
rwFiles <- rwFiles[grep(".rwl$", rwFiles)]
rwData <- list()
for(fn in rwFiles) {
    id <- gsub(".rw", "", fn)
    rwData[[id]] <- t(read.tucson(file.path(huron_wd, "RingWidths", fn)))  # rows are tree, cols are times
}

# ## initial processing of census sample dates
# threshold <- function(vals, minDate, lower = TRUE) {
#     if(!is(minDate, "Date")) minDate <- as.Date(minDate)
#     if(length(minDate) == 1) minDate <- rep(minDate, length(vals))
#     if(lower) {
#         vals[julian(vals) < julian(minDate)] <- minDate[julian(vals) < julian(minDate)]
#     } else vals[julian(vals) > julian(minDate)] <- minDate[julian(vals) > julian(minDate)]
#     return(vals)
# }
# 
# inds <- which(censusDates$X1987.1992_start == "")
# censusDates$X1987.1992_start[inds] <-  censusDates$X1987.1992_end[inds]
# censusDates$yr1991 <- as.numeric(gsub(".*(19[89][0-9])$", "\\1", censusDates$X1987.1992_start))
# 
# start <- as.Date("1969-3-31")
# end <- as.Date("1969-11-1")
# tmp1 <- as.Date(censusDates$X1969_start, format = "%m/%d/%Y")
# tmp1 <- threshold(tmp1, start)
# tmp1 <- threshold(tmp1, end, lower = FALSE)
# tmp2 <- as.Date(censusDates$X1969_end, format = "%m/%d/%Y")
# tmp2 <- threshold(tmp2, start)
# tmp2 <- threshold(tmp2, end, lower = FALSE)
# censusDates$day69 <- (julian(tmp1) + julian(tmp2)) / 2 - julian(start)
# 
# start <- as.Date("1975-3-31")
# end <- as.Date("1975-11-1")
# censusDates$X1975[is.na(censusDates$X1975)] <- "5/15/1975" # rough average of other blocks
# tmp1 <- as.Date(censusDates$X1975, format = "%m/%d/%Y")
# tmp1 <- threshold(tmp1, start)
# tmp1 <- threshold(tmp1, end, lower = FALSE)
# censusDates$day75 <- julian(tmp1) - julian(start)
# 
# start <- as.Date(paste0(censusDates$yr1991, "-3-31"))
# end <- as.Date(paste0(censusDates$yr1991, "-11-1"))
# tmp1 <- as.Date(censusDates$X1987.1992_start, format = "%m/%d/%Y")
# tmp1 <- threshold(tmp1, start)
# tmp1 <- threshold(tmp1, end, lower = FALSE)
# tmp2 <- as.Date(censusDates$X1987.1992_end, format = "%m/%d/%Y")
# tmp2 <- threshold(tmp2, start)
# tmp2 <- threshold(tmp2, end, lower = FALSE)
# censusDates$day91 <- (julian(tmp1) + julian(tmp2)) / 2 - julian(start)
# 
# start <- as.Date("2001-3-31")
# end <- as.Date("2001-11-1")
# tmp1 <- as.Date(censusDates$X2001_start, format = "%m/%d/%Y")
# tmp1 <- threshold(tmp1, start)
# tmp1 <- threshold(tmp1, end, lower = FALSE)
# tmp2 <- as.Date(censusDates$X2001_end, format = "%m/%d/%Y")
# tmp2 <- threshold(tmp2, start)
# tmp2 <- threshold(tmp2, end, lower = FALSE)
# censusDates$day01 <- (julian(tmp1) + julian(tmp2)) / 2 - julian(start)
# 
# start <- as.Date("2011-3-31")
# end <- as.Date("2011-11-1")
# tmp1 <- as.Date(censusDates$X2011_start, format = "%m/%d/%Y")
# tmp1 <- threshold(tmp1, start)
# tmp1 <- threshold(tmp1, end, lower = FALSE)
# tmp2 <- as.Date(censusDates$X2011_end, format = "%m/%d/%Y")
# tmp2 <- threshold(tmp2, start)
# tmp2 <- threshold(tmp2, end, lower = FALSE)
# censusDates$day11 <- (julian(tmp1) + julian(tmp2)) / 2 - julian(start)
# 
# start <- as.Date("1962-3-31")
# end <- as.Date("1962-11-1")
# tmp1 <- as.Date(censusDates$X1962, format = "%m/%d/%Y")
# censusDates$day62 <- julian(tmp1) - julian(start)

## get census trees only in plots
# census <- censusFull[which(censusFull$No. %in% treeMeta$No.),]
census <- censusFull[which(((censusFull$Plot %in% treeMeta$Plot) & (censusFull$Dist. <= 16))|
                             ((censusFull$Plot %in% treeMeta$Plot) & (censusFull$No. %in% treeMeta$No.))),]


# these trees we alive in census 14, but no mathching RW values
# stemindex 345 amd 375 below cutoff for coring?
census[which((census$Plot %in% treeMeta$Plot) & !(census$No. %in% treeMeta$No.) & (!is.na(census$D14))),]

# stemindex 381 was probably out of radius at 15.9 from center; remove
census = census[which(census$stemindex != 381),]

# no trees with ring-widths missing from census 
treeMeta[which((treeMeta$Plot %in% census$Plot) & !(treeMeta$No. %in% census$No.)),]
census[which((census$Plot %in% treeMeta$Plot) & !(census$No. %in% treeMeta$No.)),]



# reorganize census data
census = census[,c(1:14)]
census$stat_id = seq(1, nrow(census))
census_long = melt(census, id=c('No.', 'stat_id','Sp', 'Plot', 'Azi', 'Dist.', 'stemindex'))
colnames(census_long) = c('tree_id', 'stat_id', 'taxon', 'plot', 'azi', 'dist', 'stemindex', 'year', 'dbh')

census_long = census_long[!is.na(census_long$dbh),]
census_long$year = substr(census_long$year, 2,3)
census_long$year = as.vector(sapply(census_long$year, function(x) if (x %in% c('04', '09', '14')){paste0('20', x)} else {paste0('19', x)}))


# sum over multiple cores to get tree-level increment
# insert 0 for NA when match non-NA in another core

combineCores <- function(rw) {
    zeroGrowthFlag <- -1
    id <- as.numeric(substring(dimnames(rw)[[1]], 3, 5))
    rw <- rw[!is.na(id), ]
    id <- id[!is.na(id)]
    orient <- substring(dimnames(rw)[[1]], 6, 7)
    trees <- unique(id)
    plot <- unique(as.numeric(substring(dimnames(rw)[[1]], 3, 3)))
    if(length(plot) > 1) stop("multiple plots in object")
    vals <- matrix(zeroGrowthFlag, length(trees), ncol(rw)) # -1 is code for no ring 
    dimnames(vals)[[2]] <- dimnames(rw)[[2]]
    for(i in seq_along(trees)) {
        tree_rw <- rw[id == trees[i], , drop = FALSE]
        cat("processing ", plot, trees[i], substring(dimnames(tree_rw)[[1]], 6, 6), "\n")
        anyNonNA <- apply(tree_rw, 2, function(x) sum(!is.na(x)) > 0)
        # vals[i, anyNonNA] <- (colSums(tree_rw, na.rm = TRUE) / (0.5*nrow(tree_rw)))[anyNonNA]
        vals[i, anyNonNA] <- (colSums(tree_rw, na.rm = TRUE) / (nrow(tree_rw)))[anyNonNA]
        if(!anyNonNA[1]) {
            lastZero <- max(which(cumsum(vals[i,]) == -(1:ncol(rw))))
            vals[i, 1:lastZero] <- NA  # NA is for before tree existed
        }
    }
    dimnames(vals)[[1]] <- trees
    return(vals)
}

# # organize cores in the same way as combineCores, without cobining cores for individuals trees
# orgCores <- function(rw) {
#   zeroGrowthFlag <- -1
# #   id <- as.numeric(substring(dimnames(rw)[[1]], 3, 5))
#   id <- substring(dimnames(rw)[[1]], 5, 7)
#   rw <- rw[!is.na(id), ]
#   id <- id[!is.na(id)]
#   orient <- substring(dimnames(rw)[[1]], 8, 8)
#   trees <- unique(id)
#   plot <- unique(as.numeric(substring(dimnames(rw)[[1]], 4, 4)))
#   if(length(plot) > 1) stop("multiple plots in object")
#   vals <- matrix(zeroGrowthFlag, length(id), ncol(rw)) # -1 is code for no ring 
#   dimnames(vals)[[2]] <- dimnames(rw)[[2]]
#   for(i in seq_along(trees)) {
#     tree_rw <- rw[id == trees[i], , drop = FALSE]
#     cat("processing ", plot, trees[i], substring(dimnames(tree_rw)[[1]], 8, 8), "\n")
#     anyNonNA <- apply(tree_rw, 2, function(x) sum(!is.na(x)) > 0)
#     # vals[i, anyNonNA] <- (colSums(tree_rw, na.rm = TRUE) / (0.5*nrow(tree_rw)))[anyNonNA]
#     vals[i, anyNonNA] <- (colSums(tree_rw, na.rm = TRUE) / (nrow(tree_rw)))[anyNonNA]
#     if(!anyNonNA[1]) {
#       lastZero <- max(which(cumsum(vals[i,]) == -(1:ncol(rw))))
#       vals[i, 1:lastZero] <- NA  # NA is for before tree existed
#     }
#   }
#   dimnames(vals)[[1]] <- trees
# #   return(vals)
#   # vals = rw
#   dimnames(vals)[[2]] <- dimnames(rw)[[2]]
#   dimnames(vals)[[1]] <- id
# 
#   return(vals)
# }

# organize cores in the same way as combineCores, without cobining cores for individuals trees
orgCores <- function(rw) {
  zeroGrowthFlag <- -1
  #   id <- as.numeric(substring(dimnames(rw)[[1]], 3, 5))
  tree_ids <- substring(rownames(rw), 5, 7)
  core_ids <-  substring(rownames(rw), 4, 8)
  id <- substring(rownames(rw), 5, 7)
  rw <- rw[!is.na(id), ]
  id <- id[!is.na(id)]
  orient <- substring(rownames(rw), 8, 8)
  trees <- unique(tree_ids)
  plot <- unique(as.numeric(substring(rownames(rw), 4, 4)))
  # if(length(plot) > 1) stop("multiple plots in object")
  vals <- matrix(zeroGrowthFlag, length(id), ncol(rw)) # -1 is code for no ring 
  dimnames(vals)[[2]] <- dimnames(rw)[[2]]
  for(i in seq_along(core_ids)) {
    tree_rw <- rw[i, , drop = FALSE]
    cat("processing ", plot, core_ids[i], "\n")
    anyNonNA <- apply(tree_rw, 2, function(x) sum(!is.na(x)) > 0)
    # vals[i, anyNonNA] <- (colSums(tree_rw, na.rm = TRUE) / (0.5*nrow(tree_rw)))[anyNonNA]
    vals[i, anyNonNA] <- (colSums(tree_rw, na.rm = TRUE) / (nrow(tree_rw)))[anyNonNA]
    # if(!anyNonNA[1]) {
    #   lastZero <- max(which(cumsum(vals[i,]) == -(1:ncol(rw))))
    #   vals[i, 1:lastZero] <- NA  # NA is for before tree existed
    # }
  }
  # dimnames(vals)[[1]] <- core_ids
  #   return(vals)
  # vals = rw
  dimnames(vals)[[2]] <- dimnames(rw)[[2]]
  dimnames(vals)[[1]] <- core_ids
  
  return(vals)
}


# # single measurement per tree (average over cores)
# treeRW <- lapply(rwData, combineCores)

# multiple measurements per tree (no averaging)
treeRW <- lapply(rwData, orgCores)

# combine rw data to get increment dataset
ringYears <- sapply(treeRW, function(x) as.numeric(dimnames(x)[[2]]))
nTreesWithCores <- sum(sapply(treeRW, function(x) nrow(x)))
incr <- matrix(NA, nTreesWithCores, nT)
dimnames(incr)[[2]] <- years

start <- 1
for(i in seq_along(treeRW)) {
    end <-  start -1 + nrow(treeRW[[i]])
    include <- which(ringYears[[i]] <= lastYear & ringYears[[i]] >= firstYear)
    incr[start:end, as.character(ringYears[[i]][include])] <- treeRW[[i]][ , include]
 
    start <- end + 1
}
dimnames(incr)[[1]] <- unlist(sapply(treeRW, function(x) dimnames(x)[[1]]))

tbl <- table(census_long$taxon)
taxaMatch <- data.frame(tbl); names(taxaMatch) <- c('species', 'count')
taxaMatch$taxon <- 1:nrow(taxaMatch)

######################################################################################################################################
## make nimble data
######################################################################################################################################
if (!file.exists('data/dump')){
  dir.create('data/dump')
} 
  
incr_data = melt(incr)
colnames(incr_data) = c('id', 'year', 'incr')
incr_data$plot   = tolower(substr(incr_data$id, 1, 1))
incr_data$orient = tolower(substr(incr_data$id, 5, 5))
incr_data$id     = substr(incr_data$id, 2, 4)
incr_data$TreeID = paste0('HMC', incr_data$plot, incr_data$id)


treeMetaNo   = treeMeta[match(incr_data$TreeID, treeMeta$TreeID), 'No.']
treeMetaPlot = treeMeta[match(incr_data$TreeID, treeMeta$TreeID), 'Plot']


incr_data$stat_id = rep(NA, nrow(incr_data))
for (i in 1:nrow(incr_data)){
  print(i)
  if (any(which((census$No. == treeMetaNo[i]) & (census$Plot == treeMetaPlot[i])))){
    print(census[which((census$No. == treeMetaNo[i]) & (census$Plot == treeMetaPlot[i])), 'stat_id' ])
    # incr_data$stat_id[i] = censusFull[which((censusFull$No. == treeMetaNo[i]) & (censusFull$Plot == treeMetaPlot[i])), 'stat_id' ]
    incr_data$stat_id[i] = census[which((census$No. == treeMetaNo[i]) & (census$Plot == treeMetaPlot[i])), 'stat_id' ]
  }
}

incr_data = incr_data[which(!is.na(incr_data$stat_id)),]
incr_data = incr_data[which(incr_data$incr>=0),]

incr_data_orig = incr_data

year_start = firstYear
year_end   = lastYear
years      = seq(year_start, year_end)
incr_data  = incr_data[which(incr_data$year %in% years),]
trees_inc = sort(unique(incr_data$stat_id))

# order by tree and year
incr_data = incr_data[order(incr_data$stat_id, incr_data$year),]

# how many measurements for each tree for each year?
ncores = vector(length=nrow(incr_data))
for (i in 1:length(trees_inc)){
  for (t in 1:length(years)){
    idx = which((incr_data$stat_id == trees_inc[i]) & (incr_data$year == years[t]))
    ncores[idx] = rep(length(idx), length(idx))
  }
}

incr_data = cbind(incr_data, ncores)
incr_data = data.frame(incr_data, measno = seq(1, nrow(incr_data)))

head(incr_data)

N_inc   = nrow(incr_data) # number of measurement 
m2t     = incr_data$year
m2tree  = incr_data$stat_id
m2treecode = incr_data$id
m2nc = incr_data$ncores
m2ti    = match(m2t, years)
m2orient = incr_data$orient
Xobs    = incr_data$incr
Xobs[Xobs==0] = 0.0001
logXobs = log(Xobs)

dbh = census_long
dbh = dbh[which(dbh$year %in% years),]
dbh = dbh[order(dbh$stat_id, dbh$year),]

dbh$plot_paleon = plot_ids[match(dbh$plot, plot_ids$census_plot), 'plot']

N_dbh   = nrow(dbh)
logDobs = log(dbh$dbh)
dbh_tree_id = dbh$stat_id
dbh_day_id  = rep(100, nrow(dbh)) / growDays #as.numeric(dbh$day) / growDays # FIXME
dbh_year_id = as.numeric(dbh$year) - year_start + 1#getTimeIndex(dbh$year)
dbh_tree_code = dbh$stat_id

# make pdbh
pdbh = aggregate(year~TreeID + stat_id, incr_data, max)
pdbh = pdbh[order(pdbh$stat_id, pdbh$year),]
pdbh$dbh = treeMeta$DBH[match(pdbh$TreeID, treeMeta$TreeID)]

#### FIXING
pdbh[which(pdbh$year<2015),]

N_pdbh = nrow(pdbh)
logPDobs = log(pdbh$dbh)
pdbh_tree_id = pdbh$stat_id
pdbh_day_id  = rep(100, nrow(pdbh)) / growDays #as.numeric(dbh$day) / growDays # FIXME
pdbh_year_id = as.numeric(pdbh$year) - year_start + 1#getTimeIndex(dbh$year)
pdbh_tree_code = pdbh$stat_id

dbh_recent = dbh[which(dbh$year == 2014),]
compare_dbh = data.frame(pdbh=pdbh$dbh, dbh=dbh_recent[match(pdbh_tree_id, dbh_recent$stat_id), 'dbh'], stat_id = pdbh_tree_id)
plot(compare_dbh$pdbh, compare_dbh$dbh)
abline(a=0, b=1, col='red')

foo = compare_dbh[,1] - compare_dbh[,2]
hist(foo)


ids_table = dbh[,c('stat_id', 'tree_id' )]
ids_table = ids_table[!duplicated(ids_table$stat_id),]

trees   = sort(unique(dbh$stat_id))
N_trees = length(unique(dbh$stat_id))
N_years = length(years)

# get the site number
site_data    <- data.frame(stat_id=dbh$stat_id, site=dbh$plot_paleon, census_id=dbh$tree_id, stemindex=dbh$stemindex, plot=dbh$plot)
site_data    <- site_data[!duplicated(dbh$stat_id),]
tree_site_id <- site_data[order(site_data$stat_id), 'site']
# tree_census_id <- site_data[order(site_data$stat_id), 'census_id']
N_sites      <- length(unique(tree_site_id))


census_years = unique(census_long$year)

last_time = vector(length=N_trees)
last_time_pdbh = vector(length=N_trees)
for (i in 1:N_trees){
  tree = trees[i]
  print(tree)
  last_time[i] = max(c(dbh$year[which(dbh$stat_id == tree)], incr_data$year[which(incr_data$stat_id == tree)]), na.rm=TRUE)
  last_time_pdbh[i] = max(pdbh[pdbh$stat_id == tree,'dbh_year'], last_time[i], na.rm=TRUE)

  #print(last_time[i])
}
last_time = as.numeric(last_time)
last_time_pdbh = as.numeric(last_time_pdbh)


X_ord = data.frame(meas=numeric(0), tree_id=numeric(0), year=numeric(0))
n = 1
for (i in 1:N_trees){
  print(i)
  print(last_time[i])
  #year = seq(year_start, last_time[i])
  year = seq(year_start, last_time_pdbh[i])
  meas = seq(n, n+length(year)-1)
  n = n + length(year)
  
  X_ord = rbind(X_ord, data.frame(meas=meas, tree_id=rep(trees[i], length(year)), year=year))
}

x2tree  = X_ord$tree_id
x2year  = match(X_ord$year, years) 
# last_ti = last_time-year_start +1
last_ti = last_time_pdbh-year_start +1

# X_ord = aggregate(orient~year + stat_id, incr_data, function(x) length(unique(x)))
# X_ord = X_ord[order(X_ord$stat_id, X_ord$year),]
N_X   = nrow(X_ord)
N_D   = N_X

meas2x = vector(length=N_inc)
for (i in 1:N_inc) {
  stat_id = incr_data$stat_id[i]
  year    = incr_data$year[i]
  
  meas2x[i] = which((X_ord$tree_id == stat_id) & (X_ord$year == year))
}

meas2d = vector(length=N_dbh)
for (i in 1:N_dbh) {
  stat_id = dbh$stat_id[i]
  year    = dbh$year[i]
  
  meas2d[i] = which((X_ord$tree_id == stat_id) & (X_ord$year == year))
}


# N_pdbh = nrow(pdbh)
# logPDobs = log(pdbh$value)
# pdbh_year_id = pdbh$dbh_year - year_start + 1
# pdbh_tree_id = pdbh$stat_id
# pdbh_day_id  = as.numeric(pdbh$day)

pdbh2d = vector(length=N_pdbh)
for (i in 1:N_pdbh){
  stat_id = pdbh$stat_id[i]
  year    = pdbh$year[i]
  
  print(i)
  which((X_ord$tree_id == stat_id) & (X_ord$year == year))
  
  pdbh2d[i] = which((X_ord$tree_id == stat_id) & (X_ord$year == year))
}


# ncores_years = ncores$year
# ncores_tree  = ncores$stat_id
# ncores = ncores$orient
N_ncores = length(ncores)

i1core2m = vector(mode="integer")
i2core2m = vector(mode="integer")
i3core2m = vector(mode="integer")
i4core2m = vector(mode="integer")
n = 1
while (n <= length(m2nc)) {
  if (m2nc[n] == 1) {
    i1core2m = c(i1core2m, n)
    n = n + 1
  } else if (m2nc[n] == 2) {
    i2core2m = c(i2core2m, n)
    n = n + 2
  } else if (m2nc[n] == 3) {
    i3core2m = c(i3core2m, n)
    n = n + 3
  } else if (m2nc[n] == 4) {
    i4core2m = c(i4core2m, n)
    n = n + 4
  } else {
    print("WTF")
  }
}

n1cores = length(i1core2m)
n2cores = length(i2core2m)
n3cores = length(i3core2m)
n4cores = length(i4core2m)

ones = rep(1, 4)
open_dbh =25
N_taxa=length(unique(dbh$taxon))

cs_last_ti = cumsum(last_ti)
first_ti   = c(1,cs_last_ti[1:(length(cs_last_ti)-1)]+1)

X = rep(0.1, N_X)
D = rep(0.1, N_X)
logX = rep(log(0.1), N_X)
D0 = rep(3, N_trees)

beta = rep(0.1, N_trees)
beta_t = rep(0.1, N_years)
beta0 = 0.1
sig_x_obs = 0.7
sig_d_obs = 0.01
sig_d = 0.1
sig_d_sap = 0.1
sig_x = 0.1
beta_sd = 0.1
beta_t_sd = 0.1
beta_spp_sd = 0.1
beta_spp = rep(0.1, N_taxa)
beta_slope = 0.1

b0 = 0.1
b1 = 10

tau2 = 0
tau3 = 0
tau4 = 0

nu = 0.1
rho = 0.1

# D_saplings = rep(3, N_saplings)
# D_pre = rep(3, N_X-N_saplings)

# dump(c('X', 'logX', 'D0', 'beta', 'beta_t', 'beta0', 'sig_x_obs', 'sig_d_obs', 'sig_d', 'sig_d_sap',
#        'sig_x', 'beta_sd', 'beta_t_sd', 'beta_spp', 'beta_spp_sd', 'beta_slope', 
#        'tau2', 'tau3', 'tau4', 'b0', 'b1', 'D_saplings', 'D_pre', 'nu', 'rho'), 
#      file=paste0('data/dump/tree_data_20_', dvers, '_inits.dump'))


# average rw values
incr_data$incr[incr_data$incr == 0] = 0.0001
if (any(incr_data$incr==0)){
  print('Zero ring-widths!')
}

logOFmean=FALSE

incr_data = data.frame(incr_data, incr_log = log(incr_data$incr))
if (logOFmean){
  incrAvg = aggregate(incr ~ stat_id + year, incr_data, mean)
} else {
  incrAvg = aggregate(incr_log ~ stat_id + year, incr_data, mean)
}
# incrAvg = aggregate(incr_log ~ stat_id + year, incr_data, mean)
incrAvg = incrAvg[order(incrAvg$stat_id, incrAvg$year),]
incrAvg = data.frame(incrAvg, measno = seq(1, nrow(incrAvg)))

if (any(is.infinite(incr_data$incr_log))){
  print('Infinitely negative log ring-widths!')
}

N_inc_a   = nrow(incrAvg) # number of measurement 
m2t_a     = incrAvg$year
m2tree_a  = incrAvg$stat_id
m2treecode_a = incrAvg$id
m2ti_a    = match(m2t_a, years)
# Xobs_a    = incrAvg$incr
# Xobs_a[Xobs_a==0] = 0.0001
# logXobs_a = log(Xobs_a)
if (logOFmean){
  logXobs_a    = log(incrAvg$incr)
} else {
  logXobs_a    = incrAvg$incr_log
}
  
# X_ord the same; still estimate same times and trees!
meas2x_a = vector(length=N_inc_a)
for (i in 1:N_inc_a) {
  stat_id = incrAvg$stat_id[i]
  year    = incrAvg$year[i]
  
  meas2x_a[i] = which((X_ord$tree_id == stat_id) & (X_ord$year == year))
}

meas2d = vector(length=N_dbh)
for (i in 1:N_dbh) {
  stat_id = dbh$stat_id[i]
  year    = dbh$year[i]
  
  meas2d[i] = which((X_ord$tree_id == stat_id) & (X_ord$year == year))
}

dump(c('N_inc', 'N_dbh', 'N_pdbh', 'N_X', 'N_D', 'N_trees', 'N_years', 'N_sites',
       'logXobs', 'm2t', 'm2ti', 'm2nc', 'm2tree', 'ncores',
       'N_inc_a', 'logXobs_a', 'm2t_a', 'm2ti_a', 'm2tree_a', 
       'logDobs', 'dbh_tree_id', 'dbh_day_id', 'dbh_year_id',
       'logPDobs', 'pdbh_tree_id', 'pdbh_year_id', 'pdbh_day_id',
       'tree_site_id', 'tree_census_id',
       'n1cores', 'n2cores', 'n3cores', 'n4cores', 
       'i1core2m', 'i2core2m', 'i3core2m', 'i4core2m',
       'meas2x',  'meas2x_a','x2year', 'x2tree', 
       'meas2d',
       'pdbh2d',
       'last_ti',
       'ones',
       'year_start', 'year_end',
       'taxon', 'N_taxa', 'open_dbh',
       # 'N_saplings', 'sap2x', 'not_sap2x','max_size', 'sapling_tree_id', 'sapling_year_id',
       'first_ti', 'cs_last_ti'),
     file=paste0('data/dump/tree_data_HMC_', dvers, '.dump'))

save(N_inc, N_dbh, N_pdbh, N_X, N_D, N_trees, N_years, N_sites,
     logXobs, m2t, m2ti, m2nc, m2treecode, m2tree, ncores,
     N_inc_a, logXobs_a, m2t_a, m2ti_a, m2tree_a, 
     logDobs, dbh_tree_id, dbh_day_id, dbh_year_id,
     logPDobs, pdbh_tree_id, pdbh_year_id, pdbh_day_id,
     tree_site_id, tree_census_id,
     n1cores, n2cores, n3cores, n4cores, 
     i1core2m, i2core2m, i3core2m, i4core2m,
     meas2x, meas2x_a, x2year, x2tree, 
     meas2d,
     pdbh2d,
     last_ti, last_time, last_time_data,
     ones,
     year_start, year_end,
     trees, years, m2orient,
     taxon, N_taxa, taxaMatch, open_dbh,
     # N_saplings, sap2x, not_sap2x, max_size, sapling_tree_id, sapling_year_id,
     first_ti, cs_last_ti, dbh, pdbh,
     file=paste0('data/dump/tree_data_HMC_', dvers, '.rdata'))

######################################################################################################################################
## make nimble data: no census
######################################################################################################################################

trees_p   =  sort(unique(c(pdbh$stat_id, incr_data$stat_id)))
N_trees_p = length(trees_p)

# get the site number
tree_site_id_p <- site_data[match(trees_p, site_data$stat_id), 'site']
tree_census_id_p <- site_data[match(trees_p, site_data$stat_id), 'census_id']
N_sites_p      <- length(unique(tree_site_id_p))

last_ti_p = last_ti[match(trees_p, trees)]
last_time_p = last_time[match(trees_p, trees)]
last_time_data_p = last_time_data[match(trees_p, trees)]
last_time_pdbh_p = last_time_pdbh[match(trees_p, trees)]


X_ord_p = data.frame(meas=numeric(0), tree_id=numeric(0), year=numeric(0))
n = 1
for (i in 1:N_trees_p){
  print(i)
  print(last_time_p[i])
  #year = seq(year_start, last_time[i])
  year = seq(year_start, last_time_pdbh_p[i])
  meas = seq(n, n+length(year)-1)
  n = n + length(year)
  
  X_ord_p = rbind(X_ord_p, data.frame(meas=meas, tree_id=rep(trees_p[i], length(year)), year=year))
}

x2tree_p  = X_ord_p$tree_id
x2year_p  = match(X_ord_p$year, years) 
taxon_p   = taxon[X_ord_p$tree_id]
# last_ti = last_time-year_start +1

# X_ord = aggregate(orient~year + stat_id, incr_data, function(x) length(unique(x)))
# X_ord = X_ord[order(X_ord$stat_id, X_ord$year),]
N_X_p   = nrow(X_ord_p)
N_D_p   = N_X_p

meas2x_p = vector(length=N_inc)
for (i in 1:N_inc) {
  print(i)
  stat_id = incr_data$stat_id[i]
  year    = incr_data$year[i]
  
  meas2x_p[i] = which((X_ord_p$tree_id == stat_id) & (X_ord_p$year == year))
}

pdbh2d_p = vector(length=N_pdbh)
for (i in 1:N_pdbh){
  stat_id = pdbh$stat_id[i]
  year    = pdbh$year[i]
  
  print(i)
  
  pdbh2d_p[i] = which((X_ord_p$tree_id == stat_id) & (X_ord_p$year == year))
}



X_p = rep(0.1, N_X_p)
D_p = rep(0.1, N_X_p)
logX_p = rep(log(0.1), N_X_p)
D0 = rep(3, N_trees_p)

beta_p = rep(0.1, N_trees_p)

cs_last_ti_p = cumsum(last_ti_p)
first_ti_p   = c(1,cs_last_ti_p[1:(length(cs_last_ti_p)-1)]+1)

dump(c('X_p', 'logX_p', 'D0', 'beta', 'beta_t', 'beta0', 'sig_x_obs', 'sig_d_obs', 'sig_d', 'sig_d_sap',
       'sig_x', 'beta_sd', 'beta_t_sd', 'beta_spp', 'beta_spp_sd', 'beta_slope', 
       'tau2', 'b0', 'b1', 'nu', 'rho'), 
     file=paste0('data/dump/tree_data_HMC_', dvers, '_no_census_inits.dump'))

dump(c('N_inc', 'N_pdbh', 'N_X_p', 'N_D_p', 'N_trees_p', 'N_years', 'N_sites',
       'logXobs', 'm2t', 'm2ti', 'm2nc', 'm2tree', 'ncores',
       'logPDobs', 'pdbh_tree_id', 'pdbh_year_id', 'pdbh_day_id',
       'tree_site_id_p', 'tree_census_id_p',
       'n1cores', 'n2cores', 
       'i1core2m', 'i2core2m', 
       'meas2x_p',  'x2year_p', 'x2tree_p', 
       'pdbh2d_p',
       'last_ti_p',
       'ones',
       'year_start', 'year_end',
       'taxon_p', 'N_taxa', 'open_dbh',
       'first_ti_p', 'cs_last_ti_p'),
     file=paste0('data/dump/tree_data_HMC_no_census_', dvers, '.dump'))

save(N_inc, N_pdbh, N_X_p, N_D_p, N_trees_p, N_years, N_sites,
     logXobs, m2t, m2ti, m2nc, m2treecode, m2tree, ncores,
     logPDobs, pdbh_tree_id, pdbh_year_id, pdbh_day_id,
     tree_site_id_p, tree_census_id_p,
     n1cores, n2cores, 
     i1core2m, i2core2m, 
     meas2x_p, x2year_p, x2tree_p, 
     pdbh2d_p,
     last_ti_p, last_time_p, last_time_data_p,
     ones,
     year_start, year_end,
     trees, years, m2orient,
     taxon, N_taxa, taxaMatch, open_dbh,
     first_ti_p, cs_last_ti_p, dbh, pdbh,
     file=paste0('data/dump/tree_data_HMC_no_census_', dvers, '.rdata'))

# ######################################################################################################################################
# ## make nimble data: single core dataset
# ######################################################################################################################################
# incr_dataOrig = incr_data
# 
# #year_start = 1960
# # year_end   = 2014#max(incr_data$year)
# 
# year_start = firstYear
# year_end   = lastYear
# 
# # reset years for now
# # years = seq(year_start, year_end)
# years = seq(firstYear, lastYear)
# incr_data  = incr_data[which(incr_data$year %in% years),]
# trees_inc = sort(unique(incr_data$stat_id))
# 
# # order by tree and year
# incr_data = incr_data[order(incr_data$stat_id, incr_data$year),]
# 
# # how many measurements for each tree for each year?
# ncores = vector(length=nrow(incr_data))
# for (i in 1:length(trees_inc)){
#   for (t in 1:length(years)){
#     idx = which((incr_data$stat_id == trees_inc[i]) & (incr_data$year == years[t]))
#     ncores[idx] = rep(length(idx), length(idx))
#   }
# }
# 
# incr_data = cbind(incr_data, ncores)
# incr_data = data.frame(incr_data, measno = seq(1, nrow(incr_data)))
# 
# head(incr_data)
# 
# N_inc   = nrow(incr_data) # number of measurement 
# m2t     = incr_data$year
# m2tree  = incr_data$stat_id
# m2treecode = incr_data$id
# m2nc = incr_data$ncores
# m2ti    = match(m2t, years)
# m2orient = incr_data$orient
# Xobs    = incr_data$incr
# Xobs[Xobs==0] = 0.0001
# logXobs = log(Xobs)
# 
# dbh = dbh[which(dbh$year %in% years),]
# dbh = dbh[order(dbh$stat_id, dbh$year),]
# N_dbh   = nrow(dbh)
# logDobs = log(dbh$value)
# dbh_tree_id = dbh$stat_id
# dbh_day_id  = as.numeric(dbh$day) / growDays
# dbh_year_id = dbh$year - year_start + 1#getTimeIndex(dbh$year)
# dbh_tree_code = dbh$id
# 
# ids_table = dbh[,c(5,6,7)]
# ids_table = ids_table[!duplicated(ids_table$stat_id),]
# 
# trees   = sort(unique(dbh$stat_id))
# N_trees = length(unique(dbh$stat_id))
# N_years = length(years)
# 
# # get the site number
# site_data    <- data.frame(stat_id=dbh$stat_id, site=dbh$site, census_id=dbh$census_id)
# site_data    <- site_data[!duplicated(dbh$stat_id),]
# tree_site_id <- site_data[order(site_data$stat_id), 'site']
# tree_census_id <- site_data[order(site_data$stat_id), 'census_id']
# N_sites      <- length(unique(tree_site_id))
# 
# 
# census_years = c(1969, 1975, 1991, 2001, 2011)
# 
# last_time_data = vector(length=N_trees)
# last_time = vector(length=N_trees)
# last_time_pdbh = vector(length=N_trees)
# for (i in 1:N_trees){
#   tree = trees[i]
#   print(tree)
#   last_time[i] = max(c(dbh$year[which(dbh$stat_id == tree)], incr_data$year[which(incr_data$stat_id == tree)]), na.rm=TRUE)
#   last_time_data[i] = last_time[i]
#   if (last_time[i] %in% census_years) {
#     if (which(census_years == last_time[i]) == length(census_years)) {
#       last_time[i] = max(years)
#     } else if (last_time_data[i] == 1991) {
#       last_time[i] = 2001
#     } else {
#       last_time[i] = census_years[which(census_years == last_time[i]) + 1]
#     }
#   } else if (last_time_data[i] == 1992) {
#     last_time[i] = 2001
#   }
#   
#   #print(pdbh[pdbh$stat_id == tree,'dbh_year'])
#   last_time_pdbh[i] = max(pdbh[pdbh$stat_id == tree,'dbh_year'], last_time[i], na.rm=TRUE) 
#   
#   #print(last_time[i])
# }
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
# # X_ord = aggregate(orient~year + stat_id, incr_data, function(x) length(unique(x)))
# # X_ord = X_ord[order(X_ord$stat_id, X_ord$year),]
# N_X   = nrow(X_ord)
# N_D   = N_X
# 
# meas2x = vector(length=N_inc)
# for (i in 1:N_inc) {
#   stat_id = incr_data$stat_id[i]
#   year    = incr_data$year[i]
#   
#   meas2x[i] = which((X_ord$tree_id == stat_id) & (X_ord$year == year))
# }
# 
# meas2d = vector(length=N_dbh)
# for (i in 1:N_dbh) {
#   stat_id = dbh$stat_id[i]
#   year    = dbh$year[i]
#   
#   meas2d[i] = which((X_ord$tree_id == stat_id) & (X_ord$year == year))
# }
# 
# 
# N_pdbh = nrow(pdbh)
# logPDobs = log(pdbh$value)
# pdbh_year_id = pdbh$dbh_year - year_start + 1
# pdbh_tree_id = pdbh$stat_id
# pdbh_day_id  = as.numeric(pdbh$day)
# 
# pdbh2d = vector(length=N_pdbh)
# for (i in 1:N_pdbh){
#   stat_id = pdbh$stat_id[i]
#   year    = pdbh$dbh_year[i]
#   
#   print(i)
#   which((X_ord$tree_id == stat_id) & (X_ord$year == year))
#   
#   pdbh2d[i] = which((X_ord$tree_id == stat_id) & (X_ord$year == year))
# }
# 
# 
# # ncores_years = ncores$year
# # ncores_tree  = ncores$stat_id
# # ncores = ncores$orient
# N_ncores = length(ncores)
# 
# i1core2m = vector(mode="integer")
# i2core2m = vector(mode="integer")
# i3core2m = vector(mode="integer")
# i4core2m = vector(mode="integer")
# n = 1
# while (n <= length(m2nc)) {
#   if (m2nc[n] == 1) {
#     i1core2m = c(i1core2m, n)
#     n = n + 1
#   } else if (m2nc[n] == 2) {
#     i2core2m = c(i2core2m, n)
#     n = n + 2
#   } else if (m2nc[n] == 3) {
#     i3core2m = c(i3core2m, n)
#     n = n + 3
#   } else if (m2nc[n] == 4) {
#     i4core2m = c(i4core2m, n)
#     n = n + 4
#   } else {
#     print("WTF")
#   }
# }
# 
# n1cores = length(i1core2m)
# n2cores = length(i2core2m)
# n3cores = length(i3core2m)
# n4cores = length(i4core2m)
# 
# ones = rep(1, 4)
# open_dbh =25
# N_taxa=nTaxa
# 
# cs_last_ti = cumsum(last_ti)
# first_ti   = c(1,cs_last_ti[1:(length(cs_last_ti)-1)]+1)
# 
# N_saplings = nSaplings 
# sap2x = vector(length=N_saplings)
# for (i in 1:N_saplings) {
#   sap2x[i] = which((X_ord$tree_id == sapling_tree_id[i]) & (X_ord$year == years[sapling_year_id[i]]))
# }
# 
# not_sap2x = which(!(seq(1,N_X) %in% sap2x))
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
# tau2 = 0
# tau3 = 0
# tau4 = 0
# 
# nu = 0.1
# rho = 0.1
# 
# D_saplings = rep(3, N_saplings)
# D_pre = rep(3, N_X-N_saplings)
# 
# dump(c('X', 'logX', 'D0', 'beta', 'beta_t', 'beta0', 'sig_x_obs', 'sig_d_obs', 'sig_d', 'sig_d_sap',
#        'sig_x', 'beta_sd', 'beta_t_sd', 'beta_spp', 'beta_spp_sd', 'beta_slope', 
#        'tau2', 'tau3', 'tau4', 'b0', 'b1', 'D_saplings', 'D_pre', 'nu', 'rho'), 
#      file=paste0('data/dump/tree_data_20_', dvers, '_inits.dump'))
# 
# 
# # average rw values
# incr_data$incr[incr_data$incr == 0] = 0.0001
# if (any(incr_data$incr==0)){
#   print('Zero ring-widths!')
# }
# 
# logOFmean=FALSE
# 
# incr_data = data.frame(incr_data, incr_log = log(incr_data$incr))
# if (logOFmean){
#   incrAvg = aggregate(incr ~ stat_id + year, incr_data, mean)
# } else {
#   incrAvg = aggregate(incr_log ~ stat_id + year, incr_data, mean)
# }
# # incrAvg = aggregate(incr_log ~ stat_id + year, incr_data, mean)
# incrAvg = incrAvg[order(incrAvg$stat_id, incrAvg$year),]
# incrAvg = data.frame(incrAvg, measno = seq(1, nrow(incrAvg)))
# 
# if (any(is.infinite(incr_data$incr_log))){
#   print('Infinitely negative log ring-widths!')
# }
# 
# N_inc_a   = nrow(incrAvg) # number of measurement 
# m2t_a     = incrAvg$year
# m2tree_a  = incrAvg$stat_id
# m2treecode_a = incrAvg$id
# m2ti_a    = match(m2t_a, years)
# # Xobs_a    = incrAvg$incr
# # Xobs_a[Xobs_a==0] = 0.0001
# # logXobs_a = log(Xobs_a)
# if (logOFmean){
#   logXobs_a    = log(incrAvg$incr)
# } else {
#   logXobs_a    = incrAvg$incr_log
# }
# 
# # X_ord the same; still estimate same times and trees!
# meas2x_a = vector(length=N_inc_a)
# for (i in 1:N_inc_a) {
#   stat_id = incrAvg$stat_id[i]
#   year    = incrAvg$year[i]
#   
#   meas2x_a[i] = which((X_ord$tree_id == stat_id) & (X_ord$year == year))
# }
# 
# meas2d = vector(length=N_dbh)
# for (i in 1:N_dbh) {
#   stat_id = dbh$stat_id[i]
#   year    = dbh$year[i]
#   
#   meas2d[i] = which((X_ord$tree_id == stat_id) & (X_ord$year == year))
# }
# 
# # census dates
# aggregate(year~site, dbh, function(x) sort(unique(x)))
# 
# 
# dump(c('N_inc', 'N_dbh', 'N_pdbh', 'N_X', 'N_D', 'N_trees', 'N_years', 'N_sites',
#        'logXobs', 'm2t', 'm2ti', 'm2nc', 'm2tree', 'ncores',
#        'N_inc_a', 'logXobs_a', 'm2t_a', 'm2ti_a', 'm2tree_a', 
#        'logDobs', 'dbh_tree_id', 'dbh_day_id', 'dbh_year_id',
#        'logPDobs', 'pdbh_tree_id', 'pdbh_year_id', 'pdbh_day_id',
#        'tree_site_id', 'tree_census_id',
#        'n1cores', 'n2cores', 'n3cores', 'n4cores', 
#        'i1core2m', 'i2core2m', 'i3core2m', 'i4core2m',
#        'meas2x',  'meas2x_a','x2year', 'x2tree', 
#        'meas2d',
#        'pdbh2d',
#        'last_ti',
#        'ones',
#        'year_start', 'year_end',
#        'taxon', 'N_taxa', 'open_dbh',
#        'N_saplings', 'sap2x', 'not_sap2x','max_size', 'sapling_tree_id', 'sapling_year_id',
#        'first_ti', 'cs_last_ti'),
#      file=paste0('data/dump/tree_data_HMC_average_', dvers, '.rdata'))
# 
# save(N_inc, N_dbh, N_pdbh, N_X, N_D, N_trees, N_years, N_sites,
#      logXobs, m2t, m2ti, m2nc, m2treecode, m2tree, ncores,
#      N_inc_a, logXobs_a, m2t_a, m2ti_a, m2tree_a, 
#      logDobs, dbh_tree_id, dbh_day_id, dbh_year_id,
#      logPDobs, pdbh_tree_id, pdbh_year_id, pdbh_day_id,
#      tree_site_id, tree_census_id,
#      n1cores, n2cores, n3cores, n4cores, 
#      i1core2m, i2core2m, i3core2m, i4core2m,
#      meas2x, meas2x_a, x2year, x2tree, 
#      meas2d,
#      pdbh2d,
#      last_ti, last_time, last_time_data,
#      ones,
#      year_start, year_end,
#      trees, years, m2orient,
#      taxon, N_taxa, taxaMatch, open_dbh,
#      N_saplings, sap2x, not_sap2x, max_size, sapling_tree_id, sapling_year_id,
#      first_ti, cs_last_ti, dbh, pdbh,
#      file=paste0('data/dump/tree_data_HMC_average_', dvers, '.rdata'))