library(PEcAn.allometry)
library(ggplot2)
library(reshape2)
library(abind)

run_predict = FALSE
sample_allom = TRUE
site = "ROOSTER"

output_dir = 'ROOSTER'
if (!file.exists(paste0('output/', output_dir))){
  dir.create(paste0('output/', output_dir))
  dir.create(paste0('output/', output_dir, '/figures'))
}

dvers = "v0.1"
mvers = "v0.1"

fname_data = paste0('tree_data_ROOSTER_STAN_', dvers)
dat = readRDS(paste0('data/dump/', fname_data, '.RDS'))

N_trees = dat$N_trees
N_years = dat$N_years
N_vals = dat$N_vals
logXobs = dat$logXobs
logPDobs = dat$logPDobs
year_idx = dat$year_idx
m2taxon = dat$m2taxon
N_taxa = dat$N_taxa
pdbh_year = dat$pdbh_year
idx_tree = dat$idx_tree
pdbh2val = dat$pdbh2val
x2tree = dat$x2tree
x2year = dat$x2year
taxaMatch = dat$taxaMatch
years = dat$years
plot_id = dat$plot_id
taxon = dat$taxon
m2tree = dat$m2tree
m2t = dat$m2t
pdbh_year = dat$pdbh_year

trees = seq(1, N_trees)

# fnames = c('ring_model_t_pdbh_IT')
fnames = "ring_model_t_pdbh"
models = c('Model RW')

post= list()
for (i in 1:length(fnames)) {
  fname_model = fnames[i]
  load(file   = paste0('output/', fname_model, '_', site, '_', mvers, '.Rdata'))
  post[[i]]   = post
}  

data(allom.components)
allom.components

pfts = list(ACRU = data.frame(spcd=316,acronym="ACRU"),
            BEPA = data.frame(spcd=375,acronym="BEPA"),
            FAGR = data.frame(spcd=531,acronym="FAGR"),
            PCRU = data.frame(spcd=97,acronym="PCRU"),
            PIST = data.frame(spcd=129, acronym="PIST"),
            PRSE = data.frame(spcd=762,acronym="PRSE"),
            QURU = data.frame(spcd=833,acronym="QURU"))
pft_list  = sort(unique(names(pfts)))


# breaks if there is no dir called allom
if (!file.exists("allom")){
  dir.create("dir")
}

allom.fit   = load.allom('allom/')
if ((!(all(names(pfts) %in% names(allom.fit))))|sample_allom){
  allom.stats = AllomAve(pfts, ngibbs=500, components=6, outdir='allom')
  allom.fit   = load.allom('allom/')
}


# tbl <- table(taxon)
# taxaMatch <- data.frame(tbl); names(taxaMatch) <- c('species', 'count')
# taxaMatch$taxon <- 1:nrow(taxaMatch)

taxaMatch$group = c(1, 2, 3, 4, 5, 6, 7)
group_name      = c('Acer', 'Betula', 'Fagus', 'Picea', 'Pinus', 'Prunus', 'Quercus')
taxaMatch

niter = length(post$lp__)
keep = 500

array_or_vector_p <- function(object) {
  ## Test if object is what we consider an array for the purposes of
  ## the code below.
  ##
  ## Lists are not arrays, but vectors are, even if they don't have a
  ## dim attribute. The following give an error:
  ## - objects which are neither lists, (non-list) vectors, or arrays
  ## - vectors and arrays with zero length
  ## - arrays with rank lower than 1
  if (is.list(object))
    FALSE
  else {
    stopifnot(is.array(object) || is.vector(object))
    stopifnot(length(object) > 0)
    if (is.array(object))
      stopifnot(length(dim(object)) > 0)
    TRUE
  }
}
subarray <- function(array, index) {
  ## equivalent to array[index, , ...]
  stopifnot(array_or_vector_p(array))
  if (is.vector(array))
    array[index]
  else
    do.call("[",c(list(array,index),rep(TRUE,length(dim(array))-1)))
}
thin_subarrays <- function(subarrays, n=30) {
  ## thin subarrays. for a single n, draw that many, otherwise use that as indexes.
  indexes <- if (length(n)==1) 1:n else n
  lapply(subarrays, function(array) subarray(array, indexes))
}

post[[1]] = thin_subarrays(post[[1]], n=keep)

nmodels = length(fnames)

#dbh_1 = build_dbh_p(post[[1]], x2year, x2tree, N_years, trees, keep, taxaMatch, taxon)
build_dbh_p <- function(out, x2year, x2tree, N_years, trees, keep, taxaMatch, taxon){
  
  # col_names = sapply(strsplit(colnames(out), '\\['), function(x) x[[1]])
  col_names = names(out)
  niter     = length(out[[1]])
  
  pft   = as.vector(taxaMatch[match(taxon, taxaMatch$number),'species'])
  pfts  = sort(unique(pft))
  N_pft = length(pfts) 
  
  dbh_p = list(length=N_pft)
  tree_order = vector(length=0)
  
  for (p in 1:N_pft){
    pft_trees   = which(pft == pfts[p])
    tree_order = c(tree_order, trees[pft_trees])
    N_pft_trees = length(pft_trees)
    dbh_p[[p]] = array(NA, c(N_pft_trees, N_years, keep))
    
    for (i in 1:N_pft_trees) {
      print(paste0('Tree ', i))
      tree = trees[pft_trees[i]]
      
      tree_idx   = which(x2tree == tree)
      tree_years = x2year[tree_idx]
      
      mu_dbh = out$D[(niter-keep+1):(niter),tree_idx]
      
      # negative dbh values
      if (any(mu_dbh < 0)) {
        mu_dbh[mu_dbh < 0] = 0
      }
      
      dbh_p[[p]][i,tree_years,] = t(mu_dbh)
    }
  }
  return(list(dbh_p=dbh_p, tree_order=tree_order, pft_order = pfts))
}

build_ab_p <- function(allom.fit, dbh_p, Ntrees, Nyears, keep, pfts){
  
  ab_p = array(NA, c(Ntrees, Nyears, keep))
  
  for (k in 1:keep){
    print(paste0("Iteration ", k))
    dbh_list = sapply(dbh_p, function(x) x[,,k])
    
    pred = allom.predict(allom.fit[pfts],
                         dbh = dbh_list,
                         pft = pfts,
                         component = 6,
                         use = "Bg",
                         interval = "prediction", 
                         single.tree = FALSE)
    
    ngibbs = dim(pred[[1]])[1]
    pred_draw = sapply(pred, function(x) x[sample(seq(1,ngibbs), 1),,])
    
    pred_draw_array = array(NA, c(0, Nyears))
    for (i in 1:length(pred_draw)){
      pred_draw_array = rbind(pred_draw_array, pred_draw[[i]]) 
    }
    ab_p[,,k] = pred_draw_array
  }
  
  return(ab_p=ab_p)
}

if (run_predict){
  
  dbh_1 = build_dbh_p(post[[1]], x2year, x2tree, N_years, trees, keep, taxaMatch, taxon)
  dbh_p_1 = dbh_1$dbh_p
  tree_order_1 = dbh_1$tree_order
  pft_order_1  = dbh_1$pft_order
  
  ab_p_1 = build_ab_p(allom.fit, dbh_p_1, N_trees, N_years, keep, pft_list)
  ab_p_1 = ab_p_1[sort(tree_order_1, index.return=TRUE)$ix,,]
  
  dbh_p_1_org = abind(dbh_p_1, along=1)
  dbh_p_1 = dbh_p_1_org[sort(tree_order_1, index.return=TRUE)$ix,,]
  
  saveRDS(dbh_p_1, file=paste0('output/', site, '/dbh_p_',site, '_', mvers, '.rds'))
  saveRDS(ab_p_1, paste0('output/', site, '/ab_p_', site, '_', mvers, '.rds'))
  
} else {
  dbh_p_1 = readRDS(paste0('output/', output_dir, '/dbh_p_', site, '_', mvers, '.rds'))
  ab_p_1  = readRDS(paste0('output/', output_dir, '/ab_p_', site, '_', mvers, '.rds'))
}

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

rownames(dbh_p_1) <- seq(1, N_trees)
dbh_melt = melt(dbh_p_1)
colnames(dbh_melt) = c("tree", "year", "iter", "dbh")
dbh_melt$year = years[dbh_melt$year]
dbh_melt$plot = plot_id[dbh_melt$tree]

# saveRDS(dbh_melt, 'data/IT_DBH_iterations.RDS') 

dbh_mean = aggregate(dbh~tree+year+plot, dbh_melt, mean, na.rm=TRUE) 
dbh_sum = aggregate(dbh~year+plot, dbh_mean, sum, na.rm=TRUE)

ggplot(data=dbh_sum) + 
  geom_line(data=dbh_sum, aes(x=year, y=dbh, colour=plot, group=plot)) #+ xlim(c(1960,2012))
ggsave(paste0('output/', output_dir, '/figures/sum_dbh_by_plot_', site, '_', mvers, '.pdf'))

######################################################################################################################################
## get dbh and ab from increments and most recent dbh measurement 
######################################################################################################################################

N_samples = keep
dbh_m = array(NA, c(N_trees, N_years))
ab_m = array(NA, c(N_trees, N_years, N_samples))

for (i in 1:N_trees) {
  print(i)
  
  tree = trees[i]
  tree_idx = which(trees == tree)
  pft  = as.vector(taxaMatch[which(taxaMatch$number == taxon[tree_idx]),'species'] )
  
  # for now use average increments
  incr_tree  = exp(logXobs[which(m2tree == tree)])
  incr_years = years[m2t[which(m2tree == tree)]]
  
  incr_mean = aggregate(incr_tree~incr_years, FUN=median)
  
  incr_tree = incr_mean[,2]
  incr_years = incr_mean[,1]
  
  pdbh_year = pdbh_year[i]
  
  dbh_year = pdbh_year
  dbh_tree = exp(logPDobs[tree])
  
  incr_years = incr_years[-(length(incr_years))]
  incr_tree = incr_tree[-(length(incr_tree))]
  
  tree_years = seq(min(incr_years), max(incr_years+1))
  
  incr_cumsum = rev(cumsum(rev(incr_tree[incr_years %in% tree_years])))
  dbh_calc = c(dbh_tree - 2*incr_cumsum/10, dbh_tree)
  
  year_idx = match(tree_years, years)
  
  if (length(dbh_calc) > length(year_idx)){
    dbh_calc = dbh_calc[-1]
  }
  
  dbh_m[tree_idx,year_idx] = dbh_calc
  
  # if (last_ti[i] < length(years)) {
  #   # dbh_m[tree, (last_ti[i]+1):N_years] = rep(0, N_years - last_ti[i])
  #   dbh_m[tree_idx, (last_ti[i]+1):N_years] = rep(NA, N_years - last_ti[i])
  # }
}

# 
# pft_numbers = unique(taxon)
# pfts = toupper(as.vector(taxaMatch[match(pft_numbers, taxaMatch$species),'species'] ))
dbh_m_list = as.list(rep(NA, length(pfts)))
names(dbh_m_list) = names(pfts)
idx_m = c()

pft_list = names(pfts)
taxon_code = taxaMatch$species[match(taxon, taxaMatch$number)]

for (p in 1:length(pft_list)){
  idx_m = c(idx_m, which(taxon_code==pft_list[p]))
  dbh_m_list[[p]] = dbh_m[which(taxon_code==pft_list[p]),]
}

for (k in 1:N_samples){
  pred = allom.predict(allom.fit[pft_list],
                       dbh = dbh_m_list,
                       pft = pft_list,
                       component = 6,
                       use = "Bg",
                       interval = "prediction",
                       single.tree = FALSE)
  
  ngibbs = dim(pred[[1]])[1]
  pred_draw = sapply(pred, function(x) x[sample(seq(1,ngibbs), 1),,])
  
  pred_draw_array = array(NA, c(0, N_years))
  for (i in 1:length(pred_draw)){
    pred_draw_array = rbind(pred_draw_array, pred_draw[[i]]) 
  }
  ab_m[,,k] = pred_draw_array
}

ab_m = ab_m[sort(idx_m, index.return=TRUE)$ix,,]


for (n in 1:N_trees){
  for (t in 1:N_years){
    if (is.na(dbh_m[n,t])){
      next
    } else if (dbh_m[n,t] < 5) {
      dbh_m[n,t] = NA
      ab_m[n,t,]  = rep(NA, N_samples)
    }  
  }
}

rownames(dbh_m) <- seq(1, N_trees)
dbh_m_melt = melt(dbh_m)
colnames(dbh_m_melt) = c("tree", "year", "dbh")
dbh_m_melt$year = years[dbh_m_melt$year]
dbh_m_melt$plot = plot_id[dbh_m_melt$tree]

# dbh_m_mean = aggregate(dbh~tree+year+plot, dbh_m_melt, mean) 
dbh_m_sum = aggregate(dbh~year+plot, dbh_m_melt, sum)

ggplot() + 
  geom_line(data=dbh_sum, aes(x=year, y=dbh, colour=plot, group=plot)) + #+ xlim(c(1960,2012))
  geom_line(data=dbh_m_sum, aes(x=year, y=dbh, colour=plot, group=plot))
ggsave(paste0('output/', output_dir, '/figures/sum_dbh_by_plot_', site, '_', mvers, '.pdf'))

#########################################################################################################################################
## make the AGB plot
#########################################################################################################################################

# in Kg/plot, rescale so it is Mg/ha
ab_mr = ab_m/(30^2*pi) * 10
ab_pr1 = ab_p_1/(30^2*pi) * 10

# melt measured
ab_m_melt = melt(ab_mr)
colnames(ab_m_melt) = c('tree_id', 'year_id', 'iter', 'ab')
ab_m_melt$taxon   = taxon[ab_m_melt$tree_id]
ab_m_melt$taxon   = taxaMatch[match(ab_m_melt$taxon, taxaMatch$taxon), 'species']
ab_m_melt$year    = years[ab_m_melt$year_id]
ab_m_melt$plot    = plot_id[ab_m_melt$tree_id] 

saveRDS(ab_m_melt, file=paste0('output/', output_dir, '/ab_m_melt_', site, '.RDS'))

# melt predicted
ab_p1_melt = melt(ab_pr1)
colnames(ab_p1_melt) = c('tree_id', 'year_id', 'iter', 'ab')

# put them together
ab_p_melt = data.frame(ab_p1_melt, model=rep('Model RW'))

ab_p_melt$taxon   = taxon[ab_p_melt$tree_id]
ab_p_melt$taxon   = taxaMatch[match(ab_p_melt$taxon, taxaMatch$species), 'species']
# ab_p_melt$site_id = tree_site_id[ab_p_melt$tree_id]
ab_p_melt$year    = years[ab_p_melt$year_id]
ab_p_melt$plot    = plot_id[ab_p_melt$tree_id] 
saveRDS(ab_p_melt, paste0('output/', output_dir, '/ab_p_', site, '_', mvers, '.RDS'))

ab_p_sum_by_iter = aggregate(ab ~ year+iter+plot, ab_p_melt, function(x) sum(x, na.rm=TRUE))
# agb_sum2 = aggregate(ab~year_id+iter, data=posts_ab, function(x) sum(x, na.rm=TRUE))
ab_p_quants = aggregate(ab~year+plot, data=ab_p_sum_by_iter, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
ab_p_quants = data.frame(ab_p_quants)
ab_p_quants = cbind(ab_p_quants[,c(1,2)], ab_p_quants[,3])
colnames(ab_p_quants) = c('year', 'plot', 'ab25', 'ab50', 'ab975')  

ab_m_sum_by_iter = aggregate(ab ~ year+iter+plot, ab_m_melt, function(x) sum(x, na.rm=TRUE))
ab_m_quants = aggregate(ab~year+plot, data=ab_m_sum_by_iter, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
ab_m_quants = data.frame(ab_m_quants)
ab_m_quants = cbind(ab_m_quants[,c(1,2)], ab_m_quants[,3])
colnames(ab_m_quants) = c('year', 'plot', 'ab25', 'ab50', 'ab975')  

ab_m_sum = aggregate(ab ~ year+plot, ab_m_melt, function(x) sum(x, na.rm=TRUE))

ab_p_quants = as.data.frame(ab_p_quants)
ab_m_quants = as.data.frame(ab_m_quants)

ggplot() +  geom_ribbon(data=ab_p_quants, aes(x=year, ymin=ab25, ymax=ab975), alpha=0.4) +
  geom_line(data=ab_p_quants, aes(x=year, y=ab50), size=1) + 
  geom_line(data=ab_m_quants, aes(x=year, y=ab50, colour='Empirical RW', fill='Empirical RW'), size=1) + 
  geom_line(data=ab_p_quants, aes(x=year, y=ab25), linetype=1, size=0.5) +
  geom_line(data=ab_p_quants, aes(x=year, y=ab975), linetype=1, size=0.5) +
  # scale_fill_manual(values=cols_fill, name='Method')+
  theme_bw() + theme(axis.title=element_text(size=14), axis.text=element_text(size=14)) +
  ylab("Biomass (Mg/ha)") + xlab('Year') + facet_grid(plot~., scales="free_y") #+
  # scale_x_continuous(breaks=seq(min(years), max(years), by=5))
ggsave(file=paste0('output/', output_dir, '/figures/AGB_by_site_', site, '.pdf'))
ggsave(file=paste0('output/', output_dir, '/figures/AGB_by_site_', site, '.png'))

#########################################################################################################################################
## make some plots
#########################################################################################################################################

abi_p1 = apply(ab_pr1, c(1,3), function(x) diff(x))
abi_p1 = aperm(abi_p1, c(2, 1, 3))
abi_p1_melt = melt(abi_p1)
colnames(abi_p1_melt) = c('tree_id', 'year_id', 'iter', 'abi')

# put them together
abi_p_melt = data.frame(abi_p1_melt, model=rep('Model RW'))

# 
# abi_p = apply(ab_pr, c(1,3,4), function(x) diff(x))
# abi_p = aperm(abi_p, c(2, 1, 3, 4))
# abi_p_melt = melt(abi_p)
# colnames(abi_p_melt) = c('tree_id', 'year_id', 'iter', 'model', 'abi')

abi_p_melt$taxon   = taxon[abi_p_melt$tree_id]
abi_p_melt$taxon   = taxaMatch[match(abi_p_melt$taxon, taxaMatch$taxon), 'species']
abi_p_melt$year    = years[abi_p_melt$year_id]
abi_p_melt$plot    = plot_id[abi_p_melt$tree_id]

saveRDS(abi_p_melt, file=paste0('data/abi_p_', site, '_', mvers, '.RDS'))

# measured
ab_mr_median = apply(ab_mr, c(1,2), median, na.rm=TRUE)
abi_m = t(apply(ab_mr_median, c(1), function(x) diff(x)))

# melt measured
abi_m_melt = melt(abi_m)
colnames(abi_m_melt) = c('tree_id', 'year_id', 'abi')
abi_m_melt$taxon   = taxon[abi_m_melt$tree_id]
abi_m_melt$taxon   = taxaMatch[match(abi_m_melt$taxon, taxaMatch$taxon), 'species']
abi_m_melt$site_id = tree_site_id[abi_m_melt$tree_id]
abi_m_melt$year    = years[abi_m_melt$year_id]
abi_m_melt$plot    = plot_id[abi_m_melt$tree_id]
# ab_m_melt$iter    = rep(NA)
# ab_m_melt$model   = rep(NA)

# abi_p_melt = abi_p_melt[abi_p_melt$year<2013,]
# 
# ab_p_m1 = ab_p_melt[which(ab_p_melt$model == 1),]
# saveRDS(ab_p_m1, file=paste0('r/allom/ab_p_m1_', mvers, '.RDS'))

abi_p_sum_by_iter = aggregate(abi ~ year+plot+iter+model, abi_p_melt, function(x) sum(x, na.rm=TRUE))
# agb_sum2 = aggregate(ab~year_id+iter, data=posts_ab, function(x) sum(x, na.rm=TRUE))
abi_p_quants = aggregate(abi~year+plot+model, data=abi_p_sum_by_iter, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
abi_p_quants = data.frame(abi_p_quants)
abi_p_quants = cbind(abi_p_quants[,1:3], abi_p_quants[,4])
colnames(abi_p_quants)[4:6] = c('ab25', 'ab50', 'ab975')  

abi_m_sum = aggregate(abi ~ year+plot, abi_m_melt, function(x) sum(x, na.rm=TRUE))

cols = c('Model RW'='#8c2d04', 'Empirical RW'='#238b45')
cols_fill = c('Model RW'="#fdd0a2", 'Empirical RW'="white")

# super hack. apologies...
abi_m_sum$year = abi_m_sum$year-1

ggplot() +  geom_ribbon(data=abi_p_quants, aes(x=year, ymin=ab25, ymax=ab975, fill=model, colour=model), alpha=0.4) +
  geom_line(data=abi_p_quants, aes(x=year, y=ab50, colour=model), size=1) + 
  geom_line(data=abi_m_sum, aes(x=year, y=abi, colour='Empirical RW', fill='Empirical RW'), size=1) + 
  facet_grid(plot~.) + scale_color_manual(values=cols, name='Method')+#, labels=c('RW + Census', 'RW')) + 
  scale_fill_manual(values=cols_fill, name='Method')+#, labels=c('RW + Census', 'RW')) + 
  theme_bw()+
  ylab("Biomass Increment (Mg/ha/year)") + xlab('Year') + 
  scale_x_continuous(breaks=seq(min(years), 2015, by=10), limits=c(min(years)+4,2015)) + 
  ylim(c(0, 3))
ggsave(file=paste0('output/', output_dir, '/figures/AGBI_by_site_', site, '.pdf'))
ggsave(file=paste0('output/', output_dir, '/figures/AGBI_by_site_', site, '.png'))


#########################################################################################################################################
## combine ab and abi predictions into single data frame
#########################################################################################################################################

abi_p_melt$taxon = taxaMatch$species[match(taxon[abi_p_melt$tree_id], taxaMatch$number)]
ab_p_melt$taxon = taxaMatch$species[match(taxon[ab_p_melt$tree_id], taxaMatch$number)]

colnames(ab_p_melt)[which(colnames(ab_p_melt) == "ab")] = 'value'
colnames(abi_p_melt)[which(colnames(abi_p_melt) == "abi")] = 'value'

preds = rbind(data.frame(ab_p_melt, type="AB"),
              data.frame(abi_p_melt, type="ABI"))
write.csv(preds, file=paste0('output/', output_dir, '/NPP_STAT_MODEL_', site, '.csv'), row.names=FALSE)
saveRDS(preds, paste0('output/', output_dir, '/NPP_STAT_MODEL_', site, '.RDS'))
