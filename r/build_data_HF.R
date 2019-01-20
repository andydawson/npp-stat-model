#!/usr/bin/Rscript
library(plotrix)
library(dplR)
library(fields)
library(reshape2)

source("config")

centroids <- read.csv('data/centroids.csv')
nPlots <- nrow(centroids)

ftPerMeter <- 3.2808399

lastYear  <- 2015
firstYear <- 1960#1905
years <- firstYear:lastYear
nT <- length(years)

growDays <- 30+31+30+31+31+30+31

# wd <- paste0(dataDir, '/lyford')
lyford_wd <- paste0(dataDir, '/lyford20')


# encapsulate firstYear in getTimeIndex closure
getTimeIndex_gen <- function(firstYear) {
    function(year) {
        year - firstYear + 1
    }
}
getTimeIndex <- getTimeIndex_gen(firstYear)


## read data in

censusFull <- read.csv(paste0(lyford_wd, '/', censusFile), stringsAsFactors = FALSE)
ringMeta <- read.csv(paste0(lyford_wd, '/', ringMetaFile), stringsAsFactors = FALSE)
treeMeta <- read.csv(paste0(lyford_wd, '/', treeMetaFile), skip = 2, stringsAsFactors = FALSE)
censusDates <- read.csv(paste0(lyford_wd, '/', censusSampleDatesFile), stringsAsFactors = FALSE)
rwSampleDates <- read.csv(paste0(lyford_wd, '/', rwPlotSampleDatesFile), stringsAsFactors = FALSE)

# get the sample year for pdbh measurements
# treeMeta$sample_date = rep(rwSampleDates$revised[2])
# treeMeta$sample_date[which((treeMeta$site==1)&(treeMeta$tree.number<22))] = rep(rwSampleDates$revised[1])
# treeMeta$sample_date[which((treeMeta$site==2)&(treeMeta$tree.number<32))] = rep(rwSampleDates$revised[1])
# treeMeta$sample_date[which((treeMeta$site==3)&(treeMeta$tree.number<28))] = rep(rwSampleDates$revised[1])

# treeMeta$sample_year = rep(2014)
# treeMeta$sample_year[which((treeMeta$site==1)&(treeMeta$tree.number<22))] = rep(2013)
# treeMeta$sample_year[which((treeMeta$site==2)&(treeMeta$tree.number<32))] = rep(2013)
# treeMeta$sample_year[which((treeMeta$site==3)&(treeMeta$tree.number<28))] = rep(2013)

# we are going to say the dbh_year and sample_year are the same
treeMeta$dbh_year = rep(2014)
treeMeta$dbh_year[which((treeMeta$site==1)&(treeMeta$tree.number<22))] = rep(2013)
treeMeta$dbh_year[which((treeMeta$site==2)&(treeMeta$tree.number<32))] = rep(2013)
treeMeta$dbh_year[which((treeMeta$site==3)&(treeMeta$tree.number<28))] = rep(2013)


# # read in plot data to get live dead status
# rwStatus <- list(length=nPlots)
# for (j in 1:nPlots){
#   rwStatus[[j]] <- read.csv(paste0(lyford_wd, '/', 'Lyford_Data_Full20m_Final/LyFordPlot', j, '.csv'), 
#                             skip=7,header=TRUE,stringsAsFactors = FALSE)
# }
    
# remove second tree that has Audrey id of 1863; this one is not consistent with Audrey's data
treeMeta <- treeMeta[!(treeMeta$Site == "LF2" & treeMeta$Tree.Number == 27), ]
ringMeta <- ringMeta[!(ringMeta$SITE == "LF2" & ringMeta$TREE == 27), ]

treeMeta <- treeMeta[!(treeMeta$Site == "LF2" & treeMeta$Tree.Number == 45), ]

rwFiles <- list.files(paste0(lyford_wd, '/', "RW"))
rwFiles <- rwFiles[grep(".rwl$", rwFiles)]
rwData <- list()
for(fn in rwFiles) {
    id <- gsub(".rw", "", fn)
    rwData[[id]] <- t(read.tucson(file.path(lyford_wd, "RW", fn)))  # rows are tree, cols are times
}

## initial processing of input datasets
names(censusFull)[names(censusFull) == "treeid"] <- "census_id"

ringMeta$SITE <- as.numeric(gsub("LF", "", ringMeta$SITE))
ringMeta$X <-  NULL
names(ringMeta) <- tolower(names(ringMeta))

names(treeMeta) <- tolower(names(treeMeta))
treeMeta$site <- as.numeric(gsub("LF", "", treeMeta$site))
wh <- which(names(treeMeta) == "tag")
names(treeMeta)[wh] <- "census_id" # to match census file
treeMeta$id <- treeMeta$site*100 + treeMeta$tree.number


idx_prob = c()
treeMeta$dbh_year = NA
lastTreeIdFirstSampling <- c(21, 31, 37)
for (i in 1:nrow(treeMeta)){
  if (treeMeta$status[i] == 'Li'){
  if (treeMeta$tree.number[i]<lastTreeIdFirstSampling[treeMeta$site[i]]){
    treeMeta$dbh_year[i] = 2013
  } else {
    treeMeta$dbh_year[i] = 2014
  }
  } else {
    if (any((ringMeta$site==treeMeta$site[i]) & 
                           (ringMeta$tree==treeMeta$tree.number[i]))){
    treeMeta$dbh_year[i] = max(ringMeta[which((ringMeta$site==treeMeta$site[i]) & 
                           (ringMeta$tree==treeMeta$tree.number[i])),'meas_outer'])
    } else {
      print(i)
      idx_prob = c(idx_prob, i)
    }
  }
}


## initial processing of census sample dates
threshold <- function(vals, minDate, lower = TRUE) {
    if(!is(minDate, "Date")) minDate <- as.Date(minDate)
    if(length(minDate) == 1) minDate <- rep(minDate, length(vals))
    if(lower) {
        vals[julian(vals) < julian(minDate)] <- minDate[julian(vals) < julian(minDate)]
    } else vals[julian(vals) > julian(minDate)] <- minDate[julian(vals) > julian(minDate)]
    return(vals)
}

inds <- which(censusDates$X1987.1992_start == "")
censusDates$X1987.1992_start[inds] <-  censusDates$X1987.1992_end[inds]
censusDates$yr1991 <- as.numeric(gsub(".*(19[89][0-9])$", "\\1", censusDates$X1987.1992_start))

start <- as.Date("1969-3-31")
end <- as.Date("1969-11-1")
tmp1 <- as.Date(censusDates$X1969_start, format = "%m/%d/%Y")
tmp1 <- threshold(tmp1, start)
tmp1 <- threshold(tmp1, end, lower = FALSE)
tmp2 <- as.Date(censusDates$X1969_end, format = "%m/%d/%Y")
tmp2 <- threshold(tmp2, start)
tmp2 <- threshold(tmp2, end, lower = FALSE)
censusDates$day69 <- (julian(tmp1) + julian(tmp2)) / 2 - julian(start)

start <- as.Date("1975-3-31")
end <- as.Date("1975-11-1")
censusDates$X1975[is.na(censusDates$X1975)] <- "5/15/1975" # rough average of other blocks
tmp1 <- as.Date(censusDates$X1975, format = "%m/%d/%Y")
tmp1 <- threshold(tmp1, start)
tmp1 <- threshold(tmp1, end, lower = FALSE)
censusDates$day75 <- julian(tmp1) - julian(start)

start <- as.Date(paste0(censusDates$yr1991, "-3-31"))
end <- as.Date(paste0(censusDates$yr1991, "-11-1"))
tmp1 <- as.Date(censusDates$X1987.1992_start, format = "%m/%d/%Y")
tmp1 <- threshold(tmp1, start)
tmp1 <- threshold(tmp1, end, lower = FALSE)
tmp2 <- as.Date(censusDates$X1987.1992_end, format = "%m/%d/%Y")
tmp2 <- threshold(tmp2, start)
tmp2 <- threshold(tmp2, end, lower = FALSE)
censusDates$day91 <- (julian(tmp1) + julian(tmp2)) / 2 - julian(start)

start <- as.Date("2001-3-31")
end <- as.Date("2001-11-1")
tmp1 <- as.Date(censusDates$X2001_start, format = "%m/%d/%Y")
tmp1 <- threshold(tmp1, start)
tmp1 <- threshold(tmp1, end, lower = FALSE)
tmp2 <- as.Date(censusDates$X2001_end, format = "%m/%d/%Y")
tmp2 <- threshold(tmp2, start)
tmp2 <- threshold(tmp2, end, lower = FALSE)
censusDates$day01 <- (julian(tmp1) + julian(tmp2)) / 2 - julian(start)

start <- as.Date("2011-3-31")
end <- as.Date("2011-11-1")
tmp1 <- as.Date(censusDates$X2011_start, format = "%m/%d/%Y")
tmp1 <- threshold(tmp1, start)
tmp1 <- threshold(tmp1, end, lower = FALSE)
tmp2 <- as.Date(censusDates$X2011_end, format = "%m/%d/%Y")
tmp2 <- threshold(tmp2, start)
tmp2 <- threshold(tmp2, end, lower = FALSE)
censusDates$day11 <- (julian(tmp1) + julian(tmp2)) / 2 - julian(start)

start <- as.Date("1962-3-31")
end <- as.Date("1962-11-1")
tmp1 <- as.Date(censusDates$X1962, format = "%m/%d/%Y")
censusDates$day62 <- julian(tmp1) - julian(start)




## get census trees only in plots
census <- list(); length(census) <- nPlots
for(i in seq_len(nPlots)) {
    dist <- rdist(centroids[i, 2:3], censusFull[ , c('xsite', 'ysite')])
    census[[i]] <- censusFull[dist < plotRadius*ftPerMeter, ]
    census[[i]]$site <- i
    census[[i]]$dist_census <- dist[dist < plotRadius*ftPerMeter]
}

# check that we are getting the right trees
par(mfrow=c(1,1))
plot(censusFull$xsite, censusFull$ysite, asp=1, col="lightgrey")#, xlim=c(-500,500))
points(centroids$x, centroids$y, col="red", pch='X')

for (i in 1:nPlots){
  points(census[[i]]$xsite, census[[i]]$ysite, pch=19, col='blue')
  draw.circle(centroids$x[i], centroids$y[i], plotRadius*ftPerMeter)
}
# looks right to me!


census <- do.call(rbind, census)
dbhCols <- grep("dbh", names(census))
condCols <- grep("cond", names(census))

existFun <- function(x) {
    dbhs <- x[ , dbhCols]
    conds <- x[ , condCols]
    numAlive <- sum(conds == "L", na.rm = TRUE)
    numAlive[is.na(numAlive)] <- 0
    numAlive[]
    return(numAlive > 0)
}

# remove trees that we deem as out
out_trees <- c(3829, 1866, 1890, 2111, 995, 1089)
census = census[which(!(census$census_id %in% out_trees)),]

numAlive <- apply(census, 1, function(x) sum(x[condCols] == "L", na.rm = TRUE))
noDbh <- apply(census[ , dbhCols], 1, function(x) sum(is.na(x)) == length(dbhCols))
numAlive[noDbh] <- 0
census <- census[numAlive > 0, ]
# for now throw out trees never alive in a census - assume these were never big so limited impact on biomass increment

censusYears <- strsplit(censusYears, ",")[[1]]
censusYears2digit <- substring(censusYears, 3, 4)
nCensus <- length(censusYears)

# this is where pdbh gets added
census <- merge(census, treeMeta, by=c('census_id', 'site'), all.x = TRUE, all.y = FALSE)

census <- merge(census, censusDates[ , c('Block', 'yr1991', paste0('day', censusYears2digit))], by.x = 'block', by.y = 'Block', all.x = TRUE, all.y = FALSE)

census <- census[order(census$site, census$id), ]

census$stat_id <- 1:nrow(census)

# points(census[which(census$census_id==1890),'xsite'], census[which(census$census_id==1890),'ysite'], pch=19, col='red')
# points(census[which(census$census_id==2111),'xsite'], census[which(census$census_id==2111),'ysite'], pch=19, col='red')

## extract dbh values for living trees

dbhCols  <- paste0('dbh', censusYears2digit)
condCols <- paste0('cond', censusYears2digit)
dayCols  <- paste0('day', censusYears2digit)

colnames(census)[colnames(census) == "species.x"] = "species"

# get the dbh values from the RW plots
cols <- c("census_id", "site", "id", "stat_id", "species", "dbh", "dbh_year")
tmp <- census[ , cols]
pdbh <- melt(tmp, id.vars = c("species","census_id", "site", "stat_id", "id","dbh_year"))
pdbh <- pdbh[!is.na(pdbh$value),]
pdbh$yr = substring(pdbh$dbh_year,3,4)
pdbh$day  = rwSampleDates[match(pdbh$dbh_year, rwSampleDates$year),'date']

# FIXME
for (i in 1:nrow(pdbh)){
  if (!is.na(pdbh$dbh_year[i])){
    print(i)
    start <- as.Date(paste0(pdbh$dbh_year[i], "-3-31"))
    end   <- as.Date(paste0(pdbh$dbh_year[i], "-11-1"))
    if (pdbh$dbh_year[i] %in% c(2013,2014)) {
      tmp1 <- as.Date(pdbh$day[i])
      tmp1 <- threshold(tmp1, start)
      tmp1 <- threshold(tmp1, end, lower = FALSE)
    } else{
      tmp1   <- as.Date(paste0(pdbh$dbh_year[i], "-11-1"))
    }
      pdbh$day[i] <- julian(tmp1) - julian(start)
  } else {
    next
  }
}

pdbh = pdbh[!is.na(pdbh$day),]

# get the dbh values from the census
cols <- c("census_id", "site", "id", "stat_id", "species", dbhCols, condCols)
tmp <- census[ , cols]
dbh <- melt(tmp, id.vars = c("species","census_id", "site", "stat_id", "id", condCols))
dbh$yr=substring(dbh$variable,4,5)

cols <- c("census_id", "site", "id", "stat_id", "species", dayCols, condCols)
tmp <- census[ , cols]
day <- melt(tmp, id.vars = c("species","census_id", "site", "stat_id", "id", condCols))
day$yr=substring(day$variable,4,5)

names(day)[names(day) == 'value'] <- 'day'

dbh <- merge(dbh, day[ , c('stat_id', 'day', 'yr')], all.x = TRUE, all.y = FALSE)
# get in terms of 4-digit year, not 2-digit
match <- data.frame(yr = censusYears2digit, year = as.numeric(censusYears)) 
dbh <- merge(dbh, match, all.x = TRUE, all.y = FALSE)

dbh_good = dbh

colMatch = data.frame(censusYears2digit, which(names(dbh) %in% condCols))
names(colMatch) = c('yr','col')
dbh2 = merge(dbh, colMatch, by.x= 'yr', by.y = 'yr')
dbh2$status <- dbh2[cbind(1:nrow(dbh2), dbh2$col)]
dbh2 <- merge(dbh2, census[ , c('stat_id', 'yr1991')], all.x = TRUE, all.y = FALSE)
dbh2$year[dbh2$year == 1991] <- dbh2$yr1991[dbh2$year == 1991]

dbh <- subset(dbh2, status == "L", c('yr', 'year', 'species', 'day', 'census_id', 'id', 'stat_id', 'site', 'value'))
dbh <- subset(dbh, !is.na(value))


#   # comment this out unless doing data cleaning!
#   dbh = rbind(data.frame(dbh, type=rep('census')),
#               data.frame(yr=pdbh$yr,
#                          year=pdbh$dbh_year,
#                          species=pdbh$species,
#                          day=pdbh$day,
#                          census_id=pdbh$census_id,
#                          id=pdbh$id,
#                          stat_id=pdbh$stat_id,
#                          site=pdbh$site,
#                          value=pdbh$value,
#                          type=rep('rw')))#[,('yr', 'year','species', 'day', 'census_id', 'id', 'stat_id','site','value')])
# 
# 
# dbh_sub = dbh[which(dbh$year>=2011),]
# 
# dbh_sub = dcast(dbh_sub, stat_id~type)
# dbh_sub = transform(dbh_sub, census=as.numeric(census), rw=as.numeric(rw))
# 
# 
# dbh_sub[which(abs(dbh_sub$census - dbh_sub$rw) > 8),]
# dbh[which((dbh$stat_id %in% dbh_sub[which(abs(dbh_sub$census - dbh_sub$rw) > 8),'stat_id']) & (dbh$year >= 2011)), ]
# 
# 
# ggplot(dbh_sub) + geom_point(aes(x=census, y=rw)) + geom_abline(aes(slope=1, intercept=0))
# ggsave(file='figures/dbh_census_vs_rw.pdf')


# dbh_sub$diff = dbh_sub$rw - dbh_sub$census
# ggplot(dbh_sub) + geom_point(aes(x=census, y=diff)) + geom_hline(aes(yintercept=0)) + xlab('Census DBH') + ylab('RW - Census DBH')
# ggsave(file='figures/dbh_diff_vs_census.pdf')

####

## dbh is a core data object going into the stat model

saplings <- subset(dbh2, (is.na(status) & (yr != 62)))
saplings <- rbind(saplings, subset(dbh2, (is.na(status) & (yr == 62) & (site == 1))))
## saplings is a core data object going into the stat model

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
# 
# # organize cores in the same way as combineCores, without cobining cores for individuals trees
# orgCores <- function(rw) {
#   zeroGrowthFlag <- -1
#   id <- as.numeric(substring(dimnames(rw)[[1]], 3, 5))
#   rw <- rw[!is.na(id), ]
#   id <- id[!is.na(id)]
#   orient <- substring(dimnames(rw)[[1]], 6, 7)
#   trees <- unique(id)
#   plot <- unique(as.numeric(substring(dimnames(rw)[[1]], 3, 3)))
#   if(length(plot) > 1) stop("multiple plots in object")
#   vals <- matrix(zeroGrowthFlag, length(trees), ncol(rw)) # -1 is code for no ring 
#   dimnames(vals)[[2]] <- dimnames(rw)[[2]]
#   for(i in seq_along(trees)) {
#     tree_rw <- rw[id == trees[i], , drop = FALSE]
#     cat("processing ", plot, trees[i], substring(dimnames(tree_rw)[[1]], 6, 6), "\n")
#     anyNonNA <- apply(tree_rw, 2, function(x) sum(!is.na(x)) > 0)
#     vals[i, anyNonNA] <- (colSums(tree_rw, na.rm = TRUE) / (0.5*nrow(tree_rw)))[anyNonNA]
#     if(!anyNonNA[1]) {
#       lastZero <- max(which(cumsum(vals[i,]) == -(1:ncol(rw))))
#       vals[i, 1:lastZero] <- NA  # NA is for before tree existed
#     }
#   }
#   dimnames(vals)[[1]] <- trees
#   return(vals)
# }


# organize cores in the same way as combineCores, without cobining cores for individuals trees
orgCores <- function(rw) {
  zeroGrowthFlag <- -1
#   id <- as.numeric(substring(dimnames(rw)[[1]], 3, 5))
  id <- substring(dimnames(rw)[[1]], 3, 6)
  rw <- rw[!is.na(id), ]
  id <- id[!is.na(id)]
  orient <- substring(dimnames(rw)[[1]], 6, 7)
  trees <- unique(id)
  plot <- unique(as.numeric(substring(dimnames(rw)[[1]], 3, 3)))
  if(length(plot) > 1) stop("multiple plots in object")
  vals <- matrix(zeroGrowthFlag, length(id), ncol(rw)) # -1 is code for no ring 
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
#   return(vals)
  # vals = rw
  dimnames(vals)[[2]] <- dimnames(rw)[[2]]
  dimnames(vals)[[1]] <- id

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

# need to figure out lastDate and put in census
census$lastYear <- 2012
census$lastYear[census$cond11 != "L" & !is.na(census$cond11)] <- 2010
census$lastYear[census$cond01 != "L" & !is.na(census$cond01)] <- 2000
census$lastYear[census$cond91 != "L" & !is.na(census$cond91)] <- census$yr1991[census$cond91 != "L" & !is.na(census$cond91)] - 1
census$lastYear[census$cond75 != "L" & !is.na(census$cond75)] <- 1974
census$lastYear[census$cond69 != "L" & !is.na(census$cond69)] <- 1968
census$lastYear[census$cond62 != "L" & !is.na(census$cond62)] <- 1962

# now walk through trees with rings
firstNoGrowth <- function(x) {
    tmp <- which(x == -1)
    if(length(tmp)) return(as.numeric(names(x)[min(tmp)])) else return(NA)
}
zeroRing <- apply(incr, 1, firstNoGrowth)
zeroRing[is.na(zeroRing)] <- lastYear + 1 # kludge so that 2012 is given as last year in code below for trees with ring in 2012

censusMatches <-  which(census$id %in% substr(names(zeroRing),1,3))

census$lastYear[censusMatches] <- zeroRing[as.character(census$id[censusMatches])] - 1


# temporary restrict to one site to quicken calcs
if(F){
census <- subset(census, site == 1)
dbh <- dbh[substring(dbh$id, 1,1)=="1", ]
wh <- substring(dimnames(incr)[[1]], 1 ,1) == "1"
incr <- incr[wh,]
}

tbl <- table(census$species)
taxaMatch <- data.frame(tbl); names(taxaMatch) <- c('species', 'count')
taxaMatch$taxon <- 1:nrow(taxaMatch)

# get the site number
site_data    <- data.frame(stat_id=dbh$stat_id, site=dbh$site)
site_data    <- site_data[!duplicated(dbh$stat_id),]
tree_site_id <- site_data[order(site_data$stat_id), 'site']
nSites       <- length(unique(tree_site_id))

# figure out plot number business!
plot_data <- dbh
plot_data <- plot_data[!duplicated(plot_data$stat_id),]


dbh$value <- as.numeric(dbh$value)

nDBH <- nrow(dbh)
logDobs <- log(dbh$value)
dbh_tree_id <- dbh$stat_id
dbh_day_id <- as.numeric(dbh$day) / growDays
dbh_year_id <- getTimeIndex(dbh$year)

nSaplings <- nrow(saplings)
sapling_tree_id <- saplings$stat_id
sapling_year_id <- getTimeIndex(saplings$year)
max_size <- ifelse(saplings$year == 1969, 4.5, 5)

incrMelted <- melt(incr)
names(incrMelted) <- c('id', 'year', 'incr')

incrMelted$orient = tolower(substr(incrMelted$id, 4, 4))

incrMelted$id     = substr(incrMelted$id, 1, 3)

incrData <- merge(incrMelted, census[ , c('id', 'stat_id')])
incrData <- subset(incrData, !is.na(incr))
incrData <- subset(incrData, incr != -1)

nWidths <- nrow(incrData)
incr_tree_id <- incrData$stat_id
incr_year_id <- getTimeIndex(incrData$year)
logXobs <- log(incrData$incr)
# fudge 0 incr for now
logXobs[logXobs < log(.02)] <- log(.02)

tmp <- merge(census[ , c('stat_id', 'species')], taxaMatch[ , c('species', 'taxon')], all.x = TRUE, all.y = FALSE)
tmp <- tmp[order(tmp$stat_id), ]
taxon <- tmp$taxon
nTaxa <- nrow(taxaMatch)

n <- nrow(census)
dead <- census$lastYear < lastYear

last_time <- getTimeIndex(census$lastYear)

last_time = vector(length=n)
for (i in 1:n){
  last_time[i] = max(c(dbh$year[which(dbh$stat_id == i)], incrData$year[which(incrData$stat_id == i)]))
  print(last_time[i])
}

last_ti = match(years, last_time)

save(incrData, incr, treeMeta, ringMeta,
     file = paste0('/home/adawson/Documents/projects/npp/data/meas/incr_data_', dvers,'.Rda'))

var_wt = apply(incr, 1, var, na.rm=TRUE)
pdf('figures/variance_within_tree.pdf')
hist(var_wt, breaks=20)
dev.off()

######################################################################################################################################
## make nimble data
######################################################################################################################################
if (!file.exists('data/dump')){
  dir.create('data/dump')
} 
  


incrDataOrig = incrData

#year_start = 1960
# year_end   = 2014#max(incrData$year)

year_start = firstYear
year_end   = lastYear

# reset years for now
# years = seq(year_start, year_end)
years = seq(firstYear, lastYear)
incrData  = incrData[which(incrData$year %in% years),]
trees_inc = sort(unique(incrData$stat_id))

# order by tree and year
incrData = incrData[order(incrData$stat_id, incrData$year),]

# how many measurements for each tree for each year?
ncores = vector(length=nrow(incrData))
for (i in 1:length(trees_inc)){
  for (t in 1:length(years)){
    idx = which((incrData$stat_id == trees_inc[i]) & (incrData$year == years[t]))
    ncores[idx] = rep(length(idx), length(idx))
  }
}

incrData = cbind(incrData, ncores)
incrData = data.frame(incrData, measno = seq(1, nrow(incrData)))

head(incrData)

N_inc   = nrow(incrData) # number of measurement 
m2t     = incrData$year
m2tree  = incrData$stat_id
m2treecode = incrData$id
m2nc = incrData$ncores
m2ti    = match(m2t, years)
m2orient = incrData$orient
Xobs    = incrData$incr
Xobs[Xobs==0] = 0.0001
logXobs = log(Xobs)

dbh = dbh[which(dbh$year %in% years),]
dbh = dbh[order(dbh$stat_id, dbh$year),]
N_dbh   = nrow(dbh)
logDobs = log(dbh$value)
dbh_tree_id = dbh$stat_id
dbh_day_id  = as.numeric(dbh$day) / growDays
dbh_year_id = dbh$year - year_start + 1#getTimeIndex(dbh$year)
dbh_tree_code = dbh$id

ids_table = dbh[,c(5,6,7)]
ids_table = ids_table[!duplicated(ids_table$stat_id),]

trees   = sort(unique(dbh$stat_id))
N_trees = length(unique(dbh$stat_id))
N_years = length(years)

# get the site number
site_data    <- data.frame(stat_id=dbh$stat_id, site=dbh$site, census_id=dbh$census_id)
site_data    <- site_data[!duplicated(dbh$stat_id),]
tree_site_id <- site_data[order(site_data$stat_id), 'site']
tree_census_id <- site_data[order(site_data$stat_id), 'census_id']
N_sites      <- length(unique(tree_site_id))


census_years = c(1969, 1975, 1991, 2001, 2011)

last_time_data = vector(length=N_trees)
last_time = vector(length=N_trees)
last_time_pdbh = vector(length=N_trees)
for (i in 1:N_trees){
  tree = trees[i]
  print(tree)
  last_time[i] = max(c(dbh$year[which(dbh$stat_id == tree)], incrData$year[which(incrData$stat_id == tree)]), na.rm=TRUE)
  last_time_data[i] = last_time[i]
  if (last_time[i] %in% census_years) {
    if (which(census_years == last_time[i]) == length(census_years)) {
      last_time[i] = max(years)
    } else if (last_time_data[i] == 1991) {
      last_time[i] = 2001
    } else {
      last_time[i] = census_years[which(census_years == last_time[i]) + 1]
    }
  } else if (last_time_data[i] == 1992) {
    last_time[i] = 2001
  }
  
  #print(pdbh[pdbh$stat_id == tree,'dbh_year'])
  last_time_pdbh[i] = max(pdbh[pdbh$stat_id == tree,'dbh_year'], last_time[i], na.rm=TRUE) 
  
  #print(last_time[i])
}
last_time = as.numeric(last_time)
last_time_data = as.numeric(last_time_data)
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

# X_ord = aggregate(orient~year + stat_id, incrData, function(x) length(unique(x)))
# X_ord = X_ord[order(X_ord$stat_id, X_ord$year),]
N_X   = nrow(X_ord)
N_D   = N_X

meas2x = vector(length=N_inc)
for (i in 1:N_inc) {
  stat_id = incrData$stat_id[i]
  year    = incrData$year[i]
  
  meas2x[i] = which((X_ord$tree_id == stat_id) & (X_ord$year == year))
}

meas2d = vector(length=N_dbh)
for (i in 1:N_dbh) {
  stat_id = dbh$stat_id[i]
  year    = dbh$year[i]
  
  meas2d[i] = which((X_ord$tree_id == stat_id) & (X_ord$year == year))
}


N_pdbh = nrow(pdbh)
logPDobs = log(pdbh$value)
pdbh_year_id = pdbh$dbh_year - year_start + 1
pdbh_tree_id = pdbh$stat_id
pdbh_day_id  = as.numeric(pdbh$day)

pdbh2d = vector(length=N_pdbh)
for (i in 1:N_pdbh){
  stat_id = pdbh$stat_id[i]
  year    = pdbh$dbh_year[i]
  
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
N_taxa=nTaxa

cs_last_ti = cumsum(last_ti)
first_ti   = c(1,cs_last_ti[1:(length(cs_last_ti)-1)]+1)

N_saplings = nSaplings 
sap2x = vector(length=N_saplings)
for (i in 1:N_saplings) {
  sap2x[i] = which((X_ord$tree_id == sapling_tree_id[i]) & (X_ord$year == years[sapling_year_id[i]]))
}

not_sap2x = which(!(seq(1,N_X) %in% sap2x))

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

D_saplings = rep(3, N_saplings)
D_pre = rep(3, N_X-N_saplings)

dump(c('X', 'logX', 'D0', 'beta', 'beta_t', 'beta0', 'sig_x_obs', 'sig_d_obs', 'sig_d', 'sig_d_sap',
       'sig_x', 'beta_sd', 'beta_t_sd', 'beta_spp', 'beta_spp_sd', 'beta_slope', 
       'tau2', 'tau3', 'tau4', 'b0', 'b1', 'D_saplings', 'D_pre', 'nu', 'rho'), 
     file=paste0('data/dump/tree_data_20_', dvers, '_inits.dump'))


# average rw values
incrData$incr[incrData$incr == 0] = 0.0001
if (any(incrData$incr==0)){
  print('Zero ring-widths!')
}

logOFmean=FALSE

incrData = data.frame(incrData, incr_log = log(incrData$incr))
if (logOFmean){
  incrAvg = aggregate(incr ~ stat_id + year, incrData, mean)
} else {
  incrAvg = aggregate(incr_log ~ stat_id + year, incrData, mean)
}
# incrAvg = aggregate(incr_log ~ stat_id + year, incrData, mean)
incrAvg = incrAvg[order(incrAvg$stat_id, incrAvg$year),]
incrAvg = data.frame(incrAvg, measno = seq(1, nrow(incrAvg)))

if (any(is.infinite(incrData$incr_log))){
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

# census dates
aggregate(year~site, dbh, function(x) sort(unique(x)))


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
       'N_saplings', 'sap2x', 'not_sap2x','max_size', 'sapling_tree_id', 'sapling_year_id',
       'first_ti', 'cs_last_ti'),
     file=paste0('data/dump/tree_data_20_', dvers, '.dump'))

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
     N_saplings, sap2x, not_sap2x, max_size, sapling_tree_id, sapling_year_id,
     first_ti, cs_last_ti, dbh, pdbh,
     file=paste0('data/dump/tree_data_20_', dvers, '.rdata'))

######################################################################################################################################
## make nimble data: no census
######################################################################################################################################

trees_p   =  sort(unique(c(pdbh$stat_id, incrData$stat_id)))
N_trees_p = length(trees_p)

# get the site number
tree_site_id_p <- site_data[trees_p, 'site']
tree_census_id_p <- site_data[trees_p, 'census_id']
N_sites_p      <- length(unique(tree_site_id_p))

last_ti_p = last_ti[trees_p]
last_time_p = last_time[trees_p]
last_time_data_p = last_time_data[trees_p]
last_time_pdbh_p = last_time_pdbh[trees_p]


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

# X_ord = aggregate(orient~year + stat_id, incrData, function(x) length(unique(x)))
# X_ord = X_ord[order(X_ord$stat_id, X_ord$year),]
N_X_p   = nrow(X_ord_p)
N_D_p   = N_X_p

meas2x_p = vector(length=N_inc)
for (i in 1:N_inc) {
  print(i)
  stat_id = incrData$stat_id[i]
  year    = incrData$year[i]
  
  meas2x_p[i] = which((X_ord_p$tree_id == stat_id) & (X_ord_p$year == year))
}

pdbh2d_p = vector(length=N_pdbh)
for (i in 1:N_pdbh){
  stat_id = pdbh$stat_id[i]
  year    = pdbh$dbh_year[i]
  
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
       'tau2', 'tau3', 'tau4', 'b0', 'b1', 'nu', 'rho'), 
     file=paste0('data/dump/tree_data_20_', dvers, '_no_census_inits.dump'))

dump(c('N_inc', 'N_pdbh', 'N_X_p', 'N_D_p', 'N_trees_p', 'N_years', 'N_sites',
       'logXobs', 'm2t', 'm2ti', 'm2nc', 'm2tree', 'ncores',
       'logPDobs', 'pdbh_tree_id', 'pdbh_year_id', 'pdbh_day_id',
       'tree_site_id_p', 'tree_census_id_p',
       'n1cores', 'n2cores', 'n3cores', 'n4cores', 
       'i1core2m', 'i2core2m', 'i3core2m', 'i4core2m',
       'meas2x_p',  'x2year_p', 'x2tree_p', 
       'pdbh2d_p',
       'last_ti_p',
       'ones',
       'year_start', 'year_end',
       'taxon_p', 'N_taxa', 'open_dbh',
       'first_ti_p', 'cs_last_ti_p'),
     file=paste0('data/dump/tree_data_20_no_census_', dvers, '.dump'))

save(N_inc, N_pdbh, N_X_p, N_D_p, N_trees_p, N_years, N_sites,
     logXobs, m2t, m2ti, m2nc, m2treecode, m2tree, ncores,
     logPDobs, pdbh_tree_id, pdbh_year_id, pdbh_day_id,
     tree_site_id_p, tree_census_id_p,
     n1cores, n2cores, n3cores, n4cores, 
     i1core2m, i2core2m, i3core2m, i4core2m,
     meas2x_p, x2year_p, x2tree_p, 
     pdbh2d_p,
     last_ti_p, last_time_p, last_time_data_p,
     ones,
     year_start, year_end,
     trees, years, m2orient,
     taxon, N_taxa, taxaMatch, open_dbh,
     first_ti_p, cs_last_ti_p, dbh, pdbh,
     file=paste0('data/dump/tree_data_20_no_census_', dvers, '.rdata'))

# ######################################################################################################################################
# ## make nimble data: single core dataset
# ######################################################################################################################################
# incrDataOrig = incrData
# 
# #year_start = 1960
# # year_end   = 2014#max(incrData$year)
# 
# year_start = firstYear
# year_end   = lastYear
# 
# # reset years for now
# # years = seq(year_start, year_end)
# years = seq(firstYear, lastYear)
# incrData  = incrData[which(incrData$year %in% years),]
# trees_inc = sort(unique(incrData$stat_id))
# 
# # order by tree and year
# incrData = incrData[order(incrData$stat_id, incrData$year),]
# 
# # how many measurements for each tree for each year?
# ncores = vector(length=nrow(incrData))
# for (i in 1:length(trees_inc)){
#   for (t in 1:length(years)){
#     idx = which((incrData$stat_id == trees_inc[i]) & (incrData$year == years[t]))
#     ncores[idx] = rep(length(idx), length(idx))
#   }
# }
# 
# incrData = cbind(incrData, ncores)
# incrData = data.frame(incrData, measno = seq(1, nrow(incrData)))
# 
# head(incrData)
# 
# N_inc   = nrow(incrData) # number of measurement 
# m2t     = incrData$year
# m2tree  = incrData$stat_id
# m2treecode = incrData$id
# m2nc = incrData$ncores
# m2ti    = match(m2t, years)
# m2orient = incrData$orient
# Xobs    = incrData$incr
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
#   last_time[i] = max(c(dbh$year[which(dbh$stat_id == tree)], incrData$year[which(incrData$stat_id == tree)]), na.rm=TRUE)
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
# # X_ord = aggregate(orient~year + stat_id, incrData, function(x) length(unique(x)))
# # X_ord = X_ord[order(X_ord$stat_id, X_ord$year),]
# N_X   = nrow(X_ord)
# N_D   = N_X
# 
# meas2x = vector(length=N_inc)
# for (i in 1:N_inc) {
#   stat_id = incrData$stat_id[i]
#   year    = incrData$year[i]
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
# incrData$incr[incrData$incr == 0] = 0.0001
# if (any(incrData$incr==0)){
#   print('Zero ring-widths!')
# }
# 
# logOFmean=FALSE
# 
# incrData = data.frame(incrData, incr_log = log(incrData$incr))
# if (logOFmean){
#   incrAvg = aggregate(incr ~ stat_id + year, incrData, mean)
# } else {
#   incrAvg = aggregate(incr_log ~ stat_id + year, incrData, mean)
# }
# # incrAvg = aggregate(incr_log ~ stat_id + year, incrData, mean)
# incrAvg = incrAvg[order(incrAvg$stat_id, incrAvg$year),]
# incrAvg = data.frame(incrAvg, measno = seq(1, nrow(incrAvg)))
# 
# if (any(is.infinite(incrData$incr_log))){
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
#      file=paste0('data/dump/tree_data_20_', dvers, '.dump'))
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
#      file=paste0('data/dump/tree_data_20_', dvers, '.rdata'))