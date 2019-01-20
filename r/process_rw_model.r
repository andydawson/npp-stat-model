library(PEcAn.allometry)
library(ggplot2)
library(reshape2)
library(abind)

source('config')
run_predict  = TRUE
sample_allom = FALSE

figures_dir = 'NOCOVAR'
if (!file.exists(figures_dir)){
  dir.create(paste0("figures/", figures_dir))
}

fname_data = paste0('tree_data_20_', dvers)
load(file=paste0('data/dump/', fname_data, '.rdata'))
fname_data = paste0('tree_data_20_no_census_', dvers)
load(file=paste0('data/dump/', fname_data, '.rdata'))

# fname_model = 'ring_model_t_date_sapl_size_pdbh'
# load(file = paste0('output/', fname_model, '_', mvers, '.Rdata'))
# post = out
# col_names = sapply(strsplit(colnames(post), '\\['), function(x) x[[1]])

# fnames = c('ring_model_t_date_sapl_size_pdbh', 'ring_model_t_date_sapl_size_pdbh_nc')
fnames = c('ring_model_t_date_sapl_size_pdbh_NOCOVAR', 'ring_model_t_pdbh_nc_NOCOVAR_sigd')

# models = c('RW + Census', 'RW')
models = c('Model RW + Census', 'Model RW')

post= list()

for (i in 1:length(fnames)) {
  fname_model = fnames[i]
  load(file   = paste0('output/', fname_model, '_', mvers, '.Rdata'))
  post[[i]]   = out[1:2000,]
}  



data(allom.components)
allom.components

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

# pfts = list(ACSA = data.frame(spcd=318,acronym="ACSA3"))

# breaks if there is no dir called allom
if (!file.exists("allom")){
  dir.create("dir")
}

allom.fit   = load.allom('allom/')
if ((!(all(names(pfts) %in% names(allom.fit))))|sample_allom){
  allom.stats = AllomAve(pfts, ngibbs=200, components=6, outdir='allom')
  allom.fit   = load.allom('allom/')
}

taxaMatch$group = c(1, 1, 2, 2, 2, 3, 4, 5, 6, 7, 8, 9, 9, 9, 10)
group_name      = c('Acer', 'Betula', 'Castanea', 'Fagus', 'Fraxinus', 'Havi', 'Pinus', 'Prunus', 'Quercus', 'Tsuga')
taxaMatch

niter = dim(post[[1]])[1]
burn = 1750
keep = niter - burn

nmodels = length(fnames)

trees_p = sort(unique(x2tree_p))
taxon_p = taxon[trees_p]
x2idx_p = match(x2tree_p, unique(x2tree_p))
pft_list_p = sort(unique(pft_list[taxon_p]))

# dbh_2 = build_dbh_p(post[[2]], x2year_p, x2tree_p, trees_p, N_trees_p, N_years, keep, taxaMatch, taxon_p)
# trees_p = sort(unique(x2tree_p))
# taxon_p = taxon[trees_p]
# dbh_2 = build_dbh_p(post[[2]], x2year_p, x2tree_p, trees_p, keep, taxaMatch, taxon_p)
# dbh_p_2 = dbh_2$dbh_p
# tree_order_2 = dbh_2$tree_order
# pft_order_2  = dbh_2$pft_order

# out= post[[2]]
# rw2year=x2year_p
# rw2tree=x2tree_p
# tree_list = trees_p
# taxon_list= taxon_p

build_dbh_p <- function(out, rw2year, rw2tree, Nyears, tree_list, keep, taxaMatch, taxon_list){

  col_names = sapply(strsplit(colnames(out), '\\['), function(x) x[[1]])
  niter     = dim(out)[1]
  
  pft   = as.vector(taxaMatch[match(taxon_list, taxaMatch$taxon),'species'])
  pfts  = sort(unique(pft))
  N_pft = length(pfts) 
  
  dbh_p = list(length=N_pft)
  tree_order = vector(length=0)
  
  for (p in 1:N_pft){
    pft_trees   = which(pft == pfts[p])
    tree_order = c(tree_order, tree_list[pft_trees])
    N_pft_trees = length(pft_trees)
    dbh_p[[p]] = array(NA, c(N_pft_trees, Nyears, keep))
    
    for (i in 1:N_pft_trees) {
      print(paste0('Tree ', i))
      tree = tree_list[pft_trees[i]]

      tree_idx   = which(rw2tree == tree)
      tree_years = rw2year[tree_idx]
      
      mu_dbh = out[(niter-keep+1):(niter), which(col_names=='D')[tree_idx]]
      
      # negative dbh values
      if (any(mu_dbh < 0)) {
        mu_dbh[mu_dbh < 0] = 0
      }

      dbh_p[[p]][i,tree_years,] = t(mu_dbh)
    }
  }
  return(list(dbh_p=dbh_p, tree_order=tree_order, pft_order = pfts))
}

# 
# dbh_1 = build_dbh_p(post[[1]], x2year, x2tree, trees, N_trees, N_years, keep, taxaMatch, taxon)
# dbh_p_1 = dbh_1$dbh_p
# tree_order_1 = dbh_1$tree_order
# 
# trees_p = sort(unique(x2tree_p))
# taxon_p = taxon[trees_p]
# dbh_2 = build_dbh_p(post[[2]], x2year_p, x2tree_p, trees_p, N_trees_p, N_years, keep, taxaMatch, taxon_p)
# dbh_p_2 = dbh_2$dbh_p
# tree_order_2 = dbh_2$tree_order

# ab_p_2 = build_ab_p(allom.fit, dbh_p_2, N_trees_p, N_years, keep, taxaMatch, taxon_p, pft_list_p)

# ab_p_2 = build_ab_p(allom.fit, dbh_p_2, N_trees_p, N_years, keep, pft_list_p)

# ab_p_1 = build_ab_p(allom.fit, dbh_p_1, N_trees, N_years, keep, pft_list)
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

trees_p = sort(unique(x2tree_p))
if (run_predict){
  
  dbh_1 = build_dbh_p(post[[1]], x2year, x2tree, N_years, trees, keep, taxaMatch, taxon)
  dbh_p_1 = dbh_1$dbh_p
  tree_order_1 = dbh_1$tree_order
  pft_order_1  = dbh_1$pft_order
  
  ab_p_1 = build_ab_p(allom.fit, dbh_p_1, N_trees, N_years, keep, pft_list)
  ab_p_1 = ab_p_1[sort(tree_order_1, index.return=TRUE)$ix,,]
  
  dbh_p_1_org = abind(dbh_p_1, along=1)
  dbh_p_1 = dbh_p_1_org[sort(tree_order_1, index.return=TRUE)$ix,,]
  
  saveRDS(dbh_p_1, paste0('allom/dbh_p_1_', location, '_', mvers, '.rds'))
  saveRDS(ab_p_1, paste0('allom/ab_p_1_', location, '_', mvers, '.rds'))
  
  trees_p = sort(unique(x2tree_p))
  taxon_p = taxon[trees_p]
  dbh_2 = build_dbh_p(post[[2]], x2year_p, x2tree_p, N_years, trees_p, keep, taxaMatch, taxon_p)
  dbh_p_2 = dbh_2$dbh_p
  tree_order_2 = dbh_2$tree_order
  pft_order_2  = dbh_2$pft_order

  ab_p_2 = build_ab_p(allom.fit, dbh_p_2, N_trees_p, N_years, keep, pft_list_p)
  ab_p_2 = ab_p_2[sort(tree_order_2, index.return=TRUE)$ix,,]
  
  dbh_p_2_org = abind(dbh_p_2, along=1)
  dbh_p_2 = dbh_p_2_org[sort(tree_order_2, index.return=TRUE)$ix,,]
  
  saveRDS(dbh_p_2, paste0('allom/dbh_p_2_', location, '_', mvers, '.rds'))
  saveRDS(ab_p_2, paste0('allom/ab_p_2_', location, '_', mvers, '.rds'))
  
} else {
  dbh_p_1 = readRDS( paste0('allom/dbh_p_1_', location, '_', mvers, '.rds'))
  ab_p_1  = readRDS(paste0('allom/ab_p_1_',  location, '_', mvers, '.rds'))
  
  dbh_p_2 = readRDS( paste0('allom/dbh_p_2_', location, '_', mvers, '.rds'))
  ab_p_2  = readRDS(paste0('allom/ab_p_2_', location, '_', mvers, '.rds'))
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

dimnames(dbh_p_1)[[1]] <- seq(1, 603)
dbh_melt = melt(dbh_p_1)
colnames(dbh_melt) = c("tree", "year", "iter", "dbh")
dbh_melt$year = years[dbh_melt$year]

saveRDS(dbh_melt, 'data/HF_DBH_iterations.RDS') 

dimnames(dbh_p_2)[[1]] <- trees_p
dbh2_melt = melt(dbh_p_2)
colnames(dbh2_melt) = c("tree", "year", "iter", "dbh")
dbh2_melt$year = years[dbh2_melt$year]

dbh_melt_all = rbind(data.frame(dbh_melt, model=rep("Model RW + Census", nrow(dbh_melt))), 
                     data.frame(dbh2_melt, model=rep("Model RW", nrow(dbh2_melt))))
dbh_melt_all$plot = tree_site_id[dbh_melt_all$tree]
dbh_melt_all = dbh_melt_all[which(dbh_melt_all$dbh>=5),]

dbh_mean = aggregate(dbh~tree+plot+year+model, dbh_melt_all, mean) 
dbh_sum = aggregate(dbh~year+model+plot, dbh_mean, sum)

ggplot(data=dbh_sum) + geom_line(data=dbh_sum, aes(x=year, y=dbh, colour=model)) + xlim(c(1960,2012)) + facet_grid(plot~.)
ggsave('figures/NOCOVAR/sum_dbh_by_plot_HF_v4.0.pdf')


foo = dbh_melt_all[which((dbh_melt_all$model=="Model RW")&(dbh_melt_all$plot==1)),]

bar = dbh_melt_all[which((dbh_melt_all$model=="Model RW + Census")&(dbh_melt_all$plot==1)),]


hist(foo$dbh)

hist(bar$dbh)

## Make mean DBH by year for Istem
dbh_mean_by_iter = aggregate(dbh~plot+year+model+iter, dbh_melt_all, median) 
ggplot(data=dbh_mean_by_iter) + 
  geom_line(data=dbh_mean_by_iter, aes(x=year, y=dbh, colour=model, group=interaction(model, iter))) + 
  xlim(c(1960,2012)) + facet_grid(plot~.)
ggsave('figures/NOCOVAR/dbh_median_by_iter_HF_v4.0.pdf')

######################################################################################################################################
## predict from dbh mean
######################################################################################################################################
# first check the trace plots

plot(dbh_p_2[1,1,], type='l')

######################################################################################################################################
## smooth out death
######################################################################################################################################

# get the census years
cyr = dbh$yr
cyr[which(cyr %in% c('01', '11'))] = paste0('20', cyr[which(cyr %in% c('01', '11'))])
cyr[which(cyr %in% c('62', '69', '75', '91'))] = paste0('19', cyr[which(cyr %in%  c('62', '69', '75', '91'))])
census_years = sort(unique(cyr))
census_years = as.numeric(census_years)
idx_census   = which(years %in% census_years)
N_census_years = length(census_years)
dbh = data.frame(dbh, cyr=cyr)

# # for now don't smooth the without census data...
# smooth_1 = smooth_death(dbh_p_1, ab_p_1, last_time_data, last_time, N_trees, years, keep)

smooth_death <- function(dbh_p, ab_p, last_time_data, last_time, N_trees, years, keep){

  last_time_data[last_time_data == 1992] = 1991
  idx_adjust = which(last_time_data < 2011)

  dbh_p_smooth = dbh_p
  ab_p_smooth  = ab_p

  for (i in 1:N_trees) {
    print(i)

    idx_last = which(years == last_time_data[i])

    if (idx_last >= 52) {
      next
    } else if (last_time_data[i] < last_time[i]) {
      print('Adjusting dbh and ab due to death!')

      sample_int =  last_time[i] - last_time_data[i]

      for (k in 1:keep){
        # print(k)
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


smooth_1 = smooth_death(dbh_p_1, ab_p_1, last_time_data, last_time, N_trees, years, keep)
dbh_p_1 = smooth_1$dbh_p_smooth
ab_p_1 = smooth_1$ab_p_smooth

smooth_2 = smooth_death(dbh_p_2, ab_p_2, last_time_data_p, last_time_p, N_trees_p, years, keep)
dbh_p_2 = smooth_2$dbh_p_smooth
ab_p_2  = smooth_2$ab_p_smooth

# for now take mean values
# dbh_p_smooth = apply(dbh_p_smooth, c(1,2,4), mean, na.rm=TRUE)
# ab_p_smooth = apply(ab_p_smooth, c(1,2,4), mean, na.rm=TRUE)
# dbh_p_smooth = apply(dbh_p_smooth, c(1,2,4), function(x) quantile(x, probs=c(0.5), na.rm=TRUE))
# ab_p_smooth  = apply(ab_p_smooth, c(1,2,4),  function(x) quantile(x, probs=c(0.5), na.rm=TRUE))

dimnames(dbh_p_1)[[1]] <- seq(1, 603)
dbh_melt = melt(dbh_p_1)
colnames(dbh_melt) = c("tree", "year", "iter", "dbh")
dbh_melt$year = years[dbh_melt$year]

# saveRDS(dbh_melt, 'data/HF_DBH_iterations.RDS') 

dimnames(dbh_p_2)[[1]] <- trees_p
dbh2_melt = melt(dbh_p_2)
colnames(dbh2_melt) = c("tree", "year", "iter", "dbh")
dbh2_melt$year = years[dbh2_melt$year]

dbh_melt_all = rbind(data.frame(dbh_melt, model=rep("Model RW + Census", nrow(dbh_melt))), 
                     data.frame(dbh2_melt, model=rep("Model RW", nrow(dbh2_melt))))
dbh_melt_all$plot = tree_site_id[dbh_melt_all$tree]
dbh_melt_all = dbh_melt_all[which(dbh_melt_all$dbh>=5),]

dbh_sum = aggregate(dbh~year+model+plot+iter, dbh_melt_all, sum, na.rm=TRUE)
dbh_mean = aggregate(dbh~plot+year+model, dbh_sum, mean) 

ggplot(data=dbh_mean) + geom_line(data=dbh_mean, aes(x=year, y=dbh, colour=model)) + xlim(c(1960,2012)) + facet_grid(plot~.)
ggsave('figures/NOCOVAR/mean_sum_dbh_by_plot_HF_v4.0.pdf')

## Make mean DBH by year for Istem
dbh_mean_by_iter = aggregate(dbh~plot+year+model+iter, dbh_melt_all, median) 
ggplot(data=dbh_mean_by_iter) + 
  geom_line(data=dbh_mean_by_iter, aes(x=year, y=dbh, colour=model, group=interaction(model, iter))) + 
  xlim(c(1960,2012)) + facet_grid(plot~.)
ggsave('figures/NOCOVAR/dbh_mean_by_iter_HF_v4.0.pdf')

saveRDS(dbh_mean_by_iter, 'data/HF_DBH_median_by_plot_iter.RDS') 

######################################################################################################################################
## get dbh and ab from increments and most recent dbh measurement 
######################################################################################################################################

N_samples = 200
dbh_m = array(NA, c(N_trees, N_years))
ab_m = array(NA, c(N_trees, N_years, N_samples))

for (i in 1:N_trees) {
  print(i)
  
  tree = trees[i]
  pft  = as.vector(taxaMatch[which(taxaMatch$taxon == taxon[tree]),'species'] )
  
  if (length(which(m2tree_a == tree)) == 0) {
    print(paste0('No increments for tree ', tree, ' !'))
    
    dbh_idx   = which(dbh$stat_id == tree)
    dbh_years = dbh[dbh_idx, 'year']
    
    dbh_dat_tree = dbh[dbh_idx,]
    
    if (max(dbh_dat_tree$value) < 10) {
      print(paste0('Tree dbh is ', max(dbh_dat_tree$value), '; set ab to zero!'))
    }
    
    # dbh_m[tree, ] = rep(0, N_years)
    dbh_m[tree, ] = rep(NA, N_years)
    
  } else {
    
    # for now use average increments
    incr_tree  = exp(logXobs_a[which(m2tree_a == tree)])
    incr_years = years[m2ti_a[which(m2tree_a == tree)]]
    
    dbh_idx   = which(dbh$stat_id == tree)
    dbh_years = dbh[dbh_idx, 'year']
    
    pdbh_idx = which(pdbh$stat_id == tree)
    pdbh_year = pdbh[pdbh_idx, 'dbh_year']
    
    if (sum(dbh_years %in% c(incr_years, pdbh_year)) == 0){
      next
    }
    
    if (length(pdbh_idx) > 0){
      dbh_year = pdbh_year
      dbh_tree = pdbh[pdbh_idx, 'value']
    } else {
      dbh_idx  = dbh_idx[dbh_years %in% incr_years]
      dbh_dat_tree = dbh[dbh_idx,]
      dbh_year = dbh_dat_tree[which.max(dbh_dat_tree$year),'year']
      dbh_tree = dbh_dat_tree[which.max(dbh_dat_tree$year),'value']
    }
    
    # tree 139
    if (tree==139){
      dbh_year = 2009
      dbh_tree = 46.3
    } 
    
    if (tree == 145) {
      dbh_year = 2007
      dbh_tree = 13.3
    }
    
    incr_years = incr_years[-(length(incr_years))]
    incr_tree = incr_tree[-(length(incr_tree))]
    
    tree_years = seq(min(incr_years), max(incr_years+1))
    
    # if (dbh_year == max(incr_years)){
    #   incr_years = incr_years[-length(incr_years)]
    #   
    #   tree_years = seq(min(incr_years), dbh_year)
    #   
    # }
    
    # years
    # 
    
    # 
    incr_cumsum = rev(cumsum(rev(incr_tree[incr_years %in% tree_years])))
    dbh_calc = c(dbh_tree - 2*incr_cumsum/10, dbh_tree)
    
    year_idx = match(tree_years, years)
    
    if (length(dbh_calc) > length(year_idx)){
      dbh_calc = dbh_calc[-1]
    }
    
    # year_new = seq(min(incr_years), dbh_year)
    
    # dbh_calc = c(dbh_tree - incr_cumsum/10, dbh_tree)
    # dbh_calc[which(dbh_calc < 0)] = 0
    
    # if (min(tree_years)==year_start) {
    #   dbh_calc = dbh_calc[-1]
    #   # year_idx=year_idx[-1]
    # } else {
    #   year_idx = c(min(year_idx), year_idx)
    # }
    
    dbh_m[tree,year_idx] = dbh_calc
  }
  
  if (last_ti[i] < length(years)) {
    # dbh_m[tree, (last_ti[i]+1):N_years] = rep(0, N_years - last_ti[i])
    dbh_m[tree, (last_ti[i]+1):N_years] = rep(NA, N_years - last_ti[i])
  }
  
  # dbh_in = which((!is.na(dbh_m[tree,])) & (dbh_m[tree,] >= 0))  
  # if (length(dbh_in)==0) {
  #   next()
  # } else {
  #   pred = allom.predict(allom.fit[pft],dbh = dbh_m[tree, dbh_in],pft = pft,component = 6,use = "Bg",interval = "prediction", single.tree=TRUE)
  #   # pred = allom.predict(allom.fit[pft],dbh = dbh_m[tree, dbh_in],pft = pft,component = 6,use = "Bg",interval = "confidence")
  #   # pred_mean = colMeans(pred)
  #   pred_draw=pred[sample(seq(1,nrow(pred)), 1),]
  #   # pred_mean = apply(pred, 2, function(x) quantile(x, probs=0.5))
  #   # PI = apply(pred,2,quantile,c(0.025,0.5,0.975),na.rm=TRUE)
  #   ab_m[tree,dbh_in] = pred_draw
  # }
}


pft_numbers = unique(taxon)
pfts = toupper(as.vector(taxaMatch[match(pft_numbers, taxaMatch$taxon),'species'] ))
dbh_m_list = as.list(rep(NA, length(pft_numbers)))
names(dbh_m_list) = toupper(pfts)
idx_m = c()

for (p in 1:length(pft_numbers)){
  idx_m = c(idx_m, which(taxon==pft_numbers[p]))
  dbh_m_list[[p]] = dbh_m[which(taxon==pft_numbers[p]),]
}

for (k in 1:N_samples){
  print(k)
  pred = allom.predict(allom.fit[pfts],
                       dbh = dbh_m_list,
                       pft = pfts,
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

#####################################################################################################################################
## ab from census years alone
#####################################################################################################################################
# census_years = sort(unique(dbh$year))
census_years = sort(unique(as.vector(dbh$year)))

cyr = dbh$yr
cyr[which(cyr %in% c('01', '11'))] = paste0('20', cyr[which(cyr %in% c('01', '11'))])
cyr[which(cyr %in% c('62', '69', '75', '91'))] = paste0('19', cyr[which(cyr %in%  c('62', '69', '75', '91'))])
census_years = sort(unique(cyr))

dbh = data.frame(dbh, cyr=cyr)

N_census_years = length(census_years)

N_samples = 200
dbh_c = array(NA, c(N_trees, N_census_years))
ab_c = array(NA, c(N_trees, N_census_years, N_samples))

for (i in 1:N_trees) {
  print(i)
  
  tree = trees[i]
  tree_idx = which(trees == tree)
  
  pft  = as.vector(taxaMatch[which(taxaMatch$taxon == taxon[tree_idx]),'species'] )
  
  dbh_idx   = which(dbh$stat_id == tree)
  dbh_years = dbh[dbh_idx, 'cyr']
  dbh_dat_tree = dbh[dbh_idx,]
  
  ab_idx = match(dbh_years, census_years)
  
  dbh_c[tree_idx, ab_idx] = dbh_dat_tree$value
  
}

dbh_c_cast = dcast(dbh, stat_id+site~cyr, value.var=c('value'))

pft_numbers = unique(taxon)
pfts = toupper(as.vector(taxaMatch[match(pft_numbers, taxaMatch$taxon),'species'] ))
dbh_c_list = as.list(rep(NA, length(pft_numbers)))
names(dbh_c_list) = toupper(pfts)
idx_c = c()

for (p in 1:length(pft_numbers)){
  idx_c = c(idx_c, dbh_c_cast[which(taxon==pft_numbers[p]),1])
  dbh_c_list[[p]] = dbh_c_cast[which(taxon==pft_numbers[p]),3:ncol(dbh_c_cast)]
}

for (k in 1:N_samples){
  print(k)
  pred = allom.predict(allom.fit[pfts],
                       dbh = dbh_c_list,
                       pft = pfts,
                       component = 6,
                       use = "Bg",
                       interval = "prediction",
                       single.tree = FALSE)
  
  ngibbs = dim(pred[[1]])[1]
  pred_draw = sapply(pred, function(x) x[sample(seq(1,ngibbs), 1),,])
  
  pred_draw_array = array(NA, c(0, N_census_years))
  for (i in 1:length(pred_draw)){
    pred_draw_array = rbind(pred_draw_array, pred_draw[[i]]) 
  }
  ab_c[,,k] = pred_draw_array
}

ab_c = ab_c[sort(idx_c, index.return=TRUE)$ix,,]
#########################################################################################################################################
## make the AGB plot
#########################################################################################################################################

# in Kg/plot, rescale so it is Mg/ha
ab_mr = ab_m/(20^2*pi) * 10
ab_cr = ab_c/(20^2*pi) * 10

ab_pr1 = ab_p_1/(20^2*pi) * 10
ab_pr2 = ab_p_2/(20^2*pi) * 10

# # in Kg/plot, rescale so it is Kg/m2
# ab_mr = ab_m/(20^2*pi) 
# ab_pr = ab_p_smooth/(20^2*pi) 
# ab_cr = ab_c/(20^2*pi) 

# melt measured
ab_m_melt = melt(ab_mr)
colnames(ab_m_melt) = c('tree_id', 'year_id', 'iter', 'ab')
ab_m_melt$taxon   = taxon[ab_m_melt$tree_id]
ab_m_melt$taxon   = taxaMatch[match(ab_m_melt$taxon, taxaMatch$taxon), 'species']
ab_m_melt$site_id = tree_site_id[ab_m_melt$tree_id]
ab_m_melt$year    = years[ab_m_melt$year_id]

# melt census
ab_c_melt = melt(ab_cr)
colnames(ab_c_melt) = c('tree_id', 'year_id', 'iter', 'ab')
ab_c_melt$taxon   = taxon[ab_c_melt$tree_id]
ab_c_melt$taxon   = taxaMatch[match(ab_c_melt$taxon, taxaMatch$taxon), 'species']
ab_c_melt$site_id = tree_site_id[ab_c_melt$tree_id]
ab_c_melt$year    = as.numeric(census_years[ab_c_melt$year_id])

saveRDS(ab_c_melt, file=paste0('data/ab_c_melt_', location, '.RDS'))

# melt predicted
ab_p1_melt = melt(ab_pr1)
colnames(ab_p1_melt) = c('tree_id', 'year_id', 'iter', 'ab')
ab_p1_melt$site_id = tree_site_id[ab_p1_melt$tree_id]

ab_p2_melt = melt(ab_pr2)
colnames(ab_p2_melt) = c('tree_id', 'year_id', 'iter', 'ab')
ab_p2_melt$tree_id = trees_p[ab_p2_melt$tree_id]
ab_p2_melt$site_id = tree_site_id[ab_p2_melt$tree_id]

# put them together
ab_p_melt = rbind(data.frame(ab_p1_melt, model=rep('Model RW + Census')), 
                  data.frame(ab_p2_melt, model=rep('Model RW')))

ab_p_melt$taxon   = taxon[ab_p_melt$tree_id]
ab_p_melt$taxon   = taxaMatch[match(ab_p_melt$taxon, taxaMatch$taxon), 'species']
# ab_p_melt$site_id = tree_site_id[ab_p_melt$tree_id]
ab_p_melt$year    = years[ab_p_melt$year_id]
saveRDS(ab_p_melt, file=paste0('allom/ab_p_', location, '_', mvers, '.RDS'))

ab_p_melt = ab_p_melt[ab_p_melt$year<2013,]
ab_m_melt = ab_m_melt[ab_m_melt$year<2013,]


# merge 


ab_p_sum_by_iter = aggregate(ab ~ year+site_id+iter+model, ab_p_melt, function(x) sum(x, na.rm=TRUE))
# agb_sum2 = aggregate(ab~year_id+iter, data=posts_ab, function(x) sum(x, na.rm=TRUE))
ab_p_quants = aggregate(ab~year+site_id+model, data=ab_p_sum_by_iter, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
ab_p_quants = data.frame(ab_p_quants)
ab_p_quants = cbind(ab_p_quants[,1:3], ab_p_quants[,4])
colnames(ab_p_quants)[4:6] = c('ab25', 'ab50', 'ab975')  

colnames(dbh_mean_by_iter)[which(colnames(dbh_mean_by_iter)=="plot")] = "site_id"
fade_fix = merge(ab_p_sum_by_iter, dbh_mean_by_iter, by=c('year', 'site_id', 'model', 'iter'))
saveRDS(fade_fix , file=paste0('allom/fading_record_correct_', location, '_', mvers, '.RDS'))


ab_m_sum_by_iter = aggregate(ab ~ year+site_id+iter, ab_m_melt, function(x) sum(x, na.rm=TRUE))
ab_m_quants = aggregate(ab~year+site_id, data=ab_m_sum_by_iter, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
ab_m_quants = data.frame(ab_m_quants)
ab_m_quants = cbind(ab_m_quants[,1:2], ab_m_quants[,3])
colnames(ab_m_quants)[3:5] = c('ab25', 'ab50', 'ab975')  

ab_c_sum_by_iter = aggregate(ab ~ year+site_id+iter, ab_c_melt, function(x) sum(x, na.rm=TRUE))
ab_c_quants = aggregate(ab~year+site_id, data=ab_c_sum_by_iter, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
ab_c_quants = data.frame(ab_c_quants)
ab_c_quants = cbind(ab_c_quants[,1:2], ab_c_quants[,3])
colnames(ab_c_quants)[3:5] = c('ab25', 'ab50', 'ab975')  

ab_m_sum = aggregate(ab ~ year+site_id, ab_m_melt, function(x) sum(x, na.rm=TRUE))

saveRDS(ab_c_sum, 'data/ab_c_sum.RDS')

cols = c('#084594', '#8c2d04', '#238b45', 'black')
cols_fill = c('#4292c6', '#fdd0a2', 'white', 'white')

cols = c('Model RW + Census'='#084594', 'Model RW'='#8c2d04', 'Empirical RW'='#238b45', 'Empirical Census'='black')
cols_fill = c('Model RW + Census'="#4292c6", 'Model RW'="#fdd0a2", 'Empirical RW'="white", 'Empirical Census'="white")

ggplot() +  geom_ribbon(data=ab_p_quants, aes(x=year, ymin=ab25, ymax=ab975, fill=model, colour=model), alpha=0.4) +
  geom_line(data=ab_p_quants, aes(x=year, y=ab50, colour=model), size=1) + 
   geom_ribbon(data=ab_m_quants, aes(x=year, ymin=ab25, ymax=ab975, colour='Empirical RW', fill='Empirical RW'), alpha=0.4)+
  geom_line(data=ab_m_quants, aes(x=year, y=ab50, colour='Empirical RW', fill='Empirical RW'), size=1) + 
  geom_point(data=ab_c_quants, aes(x=year, y=ab50, colour='Empirical Census', fill='Empirical Census'), size=2) + 
  geom_errorbar(data=ab_c_quants, aes(x=year, ymin=ab25, ymax=ab975), width=0.5)+
  facet_grid(site_id~.) + scale_color_manual(values=cols, name='Method')+
  scale_fill_manual(values=cols_fill, name='Method')+
  theme_bw() + theme(axis.title=element_text(size=14), axis.text=element_text(size=14)) +
  ylab("Biomass (Mg/ha)") + xlab('Year') + 
  scale_x_continuous(breaks=seq(min(years), max(years), by=5))
ggsave(file=paste0('figures/AGB_by_site_', mvers, '.pdf'))
ggsave(file=paste0('figures/AGB_by_site_', mvers, '.png'))

ggplot() +  geom_ribbon(data=ab_p_quants, aes(x=year, ymin=ab25, ymax=ab975, fill=model), alpha=0.4) +
  geom_line(data=ab_p_quants, aes(x=year, y=ab50, colour=model), size=1) + 
  # geom_ribbon(data=ab_m_quants, aes(x=year, ymin=ab25, ymax=ab975, colour='Empirical RW', fill='Empirical RW'), alpha=0.4)+
  geom_line(data=ab_m_quants, aes(x=year, y=ab50, colour='Empirical RW', fill='Empirical RW'), size=1) + 
  geom_point(data=ab_c_quants, aes(x=year, y=ab50, colour='Empirical Census', fill='Empirical Census'), size=2) + 
  # geom_errorbar(data=ab_c_quants, aes(x=year, ymin=ab25, ymax=ab975), width=0.5)+
  geom_line(data=ab_p_quants, aes(x=year, y=ab25, colour=model), linetype=1, size=0.5) +
  geom_line(data=ab_p_quants, aes(x=year, y=ab975, colour=model), linetype=1, size=0.5) +
  facet_grid(site_id~., scales="free_y") + scale_color_manual(values=cols, name='Method')+
  scale_fill_manual(values=cols_fill, name='Method')+
  theme_bw() + theme(axis.title=element_text(size=14), axis.text=element_text(size=14)) +
  ylab("Biomass (Mg/ha)") + xlab('Year') +
  scale_x_continuous(breaks=seq(min(years), max(years), by=5))
ggsave(file=paste0('figures/', figures_dir, '/AGB_by_site_', location, '.pdf'))
ggsave(file=paste0('figures/', figures_dir, '/AGB_by_site_', location, '.png'))

#########################################################################################################################################
## plot mean biomass
#########################################################################################################################################
dimnames(ab_p_1)[[1]] <- seq(1, 603)
ab_melt = melt(ab_p_1)
colnames(ab_melt) = c("tree", "year", "iter", "ab")
ab_melt$year = years[ab_melt$year]

# saveRDS(dbh_melt, 'data/HF_DBH_iterations.RDS') 

dimnames(ab_p_2)[[1]] <- trees_p
ab2_melt = melt(ab_p_2)
colnames(ab2_melt) = c("tree", "year", "iter", "ab")
ab2_melt$year = years[ab2_melt$year]

ab_melt_all = rbind(data.frame(ab_melt, model=rep("Model RW + Census", nrow(ab_melt))), 
                     data.frame(ab2_melt, model=rep("Model RW", nrow(ab2_melt))))
ab_melt_all$plot = tree_site_id[ab_melt_all$tree]
# ab_melt_all = ab_melt_all[which(dbh_melt_all$dbh>=5),]

ab_mean = aggregate(ab~tree+plot+year+model, ab_melt_all, mean) 
ab_sum = aggregate(ab~year+model+plot, ab_mean, sum)

ggplot(data=ab_sum) + geom_line(data=ab_sum, aes(x=year, y=ab, colour=model)) + xlim(c(1960,2012)) + facet_grid(plot~.)
ggsave('figures/NOCOVAR/sum_ab_by_plot_HF_v4.0.pdf')

#########################################################################################################################################
## make some plots
#########################################################################################################################################

abi_p1 = apply(ab_pr1, c(1,3), function(x) diff(x))
abi_p1 = aperm(abi_p1, c(2, 1, 3))
abi_p1_melt = melt(abi_p1)
colnames(abi_p1_melt) = c('tree_id', 'year_id', 'iter', 'abi')

saveRDS(abi_p1, paste0("allom/abi_p1_",mvers, ".rds" ))

abi_p2 = apply(ab_pr2, c(1,3), function(x) diff(x))
abi_p2 = aperm(abi_p2, c(2, 1, 3))
abi_p2_melt = melt(abi_p2)
colnames(abi_p2_melt) = c('tree_id', 'year_id', 'iter', 'abi')
abi_p2_melt$tree_id = trees_p[abi_p2_melt$tree_id]

# put them together
abi_p_melt = rbind(data.frame(abi_p1_melt, model=rep('Model RW + Census')), data.frame(abi_p2_melt, model=rep('Model RW')))

# 
# abi_p = apply(ab_pr, c(1,3,4), function(x) diff(x))
# abi_p = aperm(abi_p, c(2, 1, 3, 4))
# abi_p_melt = melt(abi_p)
# colnames(abi_p_melt) = c('tree_id', 'year_id', 'iter', 'model', 'abi')

abi_p_melt$taxon   = taxon[abi_p_melt$tree_id]
abi_p_melt$taxon   = taxaMatch[match(abi_p_melt$taxon, taxaMatch$taxon), 'species']
abi_p_melt$site_id = tree_site_id[abi_p_melt$tree_id]
abi_p_melt$year    = years[abi_p_melt$year_id]

saveRDS(abi_p_melt, file=paste0('allom/abi_p_', location, '_', mvers, '.RDS'))

# measured
ab_mr_median = apply(ab_mr, c(1,2), median, na.rm=TRUE)
abi_m = t(apply(ab_mr, 1, function(x) diff(x)))

# census
ab_cr_median = apply(ab_cr, c(1,2), median, na.rm=TRUE)
census_years = as.numeric(census_years)
abi_c = t(apply(ab_cr, 1, function(x) diff(x)))
for (i in 1:ncol(abi_c)){
  abi_c[,i] = abi_c[,i]/diff(census_years)[i]
}

# melt measured
abi_m_melt = melt(abi_m)
colnames(abi_m_melt) = c('tree_id', 'year_id', 'abi')
abi_m_melt$taxon   = taxon[abi_m_melt$tree_id]
abi_m_melt$taxon   = taxaMatch[match(abi_m_melt$taxon, taxaMatch$taxon), 'species']
abi_m_melt$site_id = tree_site_id[abi_m_melt$tree_id]
abi_m_melt$year    = years[abi_m_melt$year_id]
# ab_m_melt$iter    = rep(NA)
# ab_m_melt$model   = rep(NA)

# melt census
abi_c_melt = melt(abi_c)
colnames(abi_c_melt) = c('tree_id', 'year_id', 'abi')
abi_c_melt$taxon   = taxon[abi_c_melt$tree_id]
abi_c_melt$taxon   = taxaMatch[match(abi_c_melt$taxon, taxaMatch$taxon), 'species']
abi_c_melt$site_id = tree_site_id[abi_c_melt$tree_id]
abi_c_melt$year    = as.numeric(census_years[abi_c_melt$year_id])
# ab_c_melt$iter    = rep(NA)
# ab_c_melt$model   = rep(NA)


abi_p_melt = abi_p_melt[abi_p_melt$year<2013,]
abi_m_melt = abi_m_melt[abi_m_melt$year<2013,]

# 
# ab_p_m1 = ab_p_melt[which(ab_p_melt$model == 1),]
# saveRDS(ab_p_m1, file=paste0('r/allom/ab_p_m1_', mvers, '.RDS'))

abi_p_sum_by_iter = aggregate(abi ~ year+site_id+iter+model, abi_p_melt, function(x) sum(x, na.rm=TRUE))
# agb_sum2 = aggregate(ab~year_id+iter, data=posts_ab, function(x) sum(x, na.rm=TRUE))
abi_p_quants = aggregate(abi~year+site_id+model, data=abi_p_sum_by_iter, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
abi_p_quants = data.frame(abi_p_quants)
abi_p_quants = cbind(abi_p_quants[,1:3], abi_p_quants[,4])
colnames(abi_p_quants)[4:6] = c('ab25', 'ab50', 'ab975')  

abi_c_sum = aggregate(abi ~ year+site_id, abi_c_melt, function(x) sum(x, na.rm=TRUE))
abi_m_sum = aggregate(abi ~ year+site_id, abi_m_melt, function(x) sum(x, na.rm=TRUE))

# cols = c('#084594', '#8c2d04')
# cols_fill = c('#4292c6', '#feedde')

cols = c('Model RW + Census'='#084594', 'Model RW'='#8c2d04', 'Empirical RW'='#238b45', 'Empirical Census'='black')
cols_fill = c('Model RW + Census'="#4292c6", 'Model RW'="#fdd0a2", 'Empirical RW'="white", 'Empirical Census'="white")

# super hack. apologies...
abi_m_sum$year = abi_m_sum$year-1

ggplot() +  geom_ribbon(data=abi_p_quants, aes(x=year, ymin=ab25, ymax=ab975, fill=model, colour=model), alpha=0.4) +
  geom_line(data=abi_p_quants, aes(x=year, y=ab50, colour=model), size=1) + 
  geom_line(data=abi_m_sum, aes(x=year, y=abi, colour='Empirical RW', fill='Empirical RW'), size=1) + 
  # geom_point(data=abi_c_sum, aes(x=year, y=abi, colour='Empirical Census', fill='Empirical Census'),size=2) + 
  facet_grid(site_id~.) + scale_color_manual(values=cols, name='Method')+#, labels=c('RW + Census', 'RW')) + 
  scale_fill_manual(values=cols_fill, name='Method')+#, labels=c('RW + Census', 'RW')) + 
  theme_bw()+
  ylab("Biomass Increment (Mg/ha/year)") + xlab('Year') +
  scale_x_continuous(breaks=seq(min(years), max(years), by=5))
ggsave(file=paste0('figures/', figures_dir, '/AGBI_by_site_', location, '.pdf'))
ggsave(file=paste0('figures/', figures_dir, '/AGBI_by_site_', location, '.png'))


# #########################################################################################################################################
# ## make some plots
# #########################################################################################################################################
# 
# 
# ########################################################################################
# 
# dbh_p_melt = melt(dbh_p_smooth)
# colnames(dbh_p_melt) = c('tree_id', 'year_id', 'model', 'dbh')
# dbh_p_melt = dbh_p_melt[,c(3,1,2,4)]
# 
# ab_p_melt = melt(ab_pr)
# colnames(ab_p_melt) = c('tree_id', 'year_id', 'model', 'ab')
# ab_p_melt = ab_p_melt[,c(3,1,2,4)]
# 
# ab_p_melt$taxon = taxon[ab_p_melt$tree_id]
# ab_p_melt$taxon = taxaMatch[match(ab_p_melt$taxon, taxaMatch$taxon), 'species']
# 
# saveRDS(ab_p_melt, file=paste0('r/allom/ab_p', mvers, '.RDS'))
# saveRDS(taxon, file='r/allom/tree_to_taxon.RDS')
# saveRDS(taxaMatch, file='r/allom/taxon_table.RDS')
# 
# dbh_m_melt = melt(dbh_m)
# colnames(dbh_m_melt) = c('tree_id', 'year_id', 'dbh')
# 
# ab_m_melt = melt(ab_mr)
# colnames(ab_m_melt) = c('tree_id', 'year_id', 'ab')
# 
# dbh_c_melt = melt(dbh_c)
# colnames(dbh_c_melt) = c('tree_id', 'year_id', 'dbh')
# 
# ab_c_melt = melt(ab_cr)
# colnames(ab_c_melt) = c('tree_id', 'year_id', 'ab')
# 
# 
# dbh_all = rbind(cbind(dbh_p_melt, type=rep('Predicted', nrow(dbh_p_melt))),
#                 cbind(model=rep('mean', nrow(dbh_m_melt)), dbh_m_melt, type=rep('Measured', nrow(dbh_m_melt))))
# dbh_all$site_id=tree_site_id[dbh_all$tree_id]
# dbh_all$year=years[dbh_all$year_id]
# 
# ab_all = rbind(cbind(ab_p_melt, type=rep('Predicted', nrow(ab_p_melt))),
#                cbind(model=rep('mean', nrow(ab_m_melt)), ab_m_melt, type=rep('Measured', nrow(ab_m_melt))))
# ab_all$site_id=tree_site_id[ab_all$tree_id]
# ab_all$year=years[ab_all$year_id]
# 
# ab_c_melt = data.frame(ab_c_melt, site_id=tree_site_id[ab_c_melt$tree_id], year=census_years[ab_c_melt$year_id])
# # ab_c_melt = ab_c_melt[ab_c_melt$year_id < 53,]
# ab_c_melt = data.frame(ab_c_melt, taxon=taxon[ab_c_melt$tree_id])
# 
# #########################################################################################################################################
# ## increment
# #########################################################################################################################################
# 
# abi_m = t(apply(ab_mr, 1, function(x) diff(x)))
# 
# 
# abi_p = apply(ab_pr, c(1,3,4), function(x) diff(x))
# abi_p_melt = melt(abi_p)
# colnames(abi_p_melt) = c('year', 'tree_id', 'iter', 'model', 'abi')
# 
# abi_p_melt$site_id = tree_site_id[abi_p_melt$tree_id]
# abi_all_site  = aggregate(abi  ~ site_id + year + iter + model, abi_p_melt, sum, na.rm=TRUE)
# 
# abi_all_site_quants = aggregate(abi ~ site_id + year + model, abi_all_site, mean, na.rm=TRUE)
# 
# ggplot() + geom_line(data=abi_all_site_quants, aes(x=year, y=abi, group=factor(model), colour=factor(model))) + facet_grid(site_id~.)
# 
# ### to here
# 
# # # abi_p = t(apply(ab_pr, 1, function(x) diff(x)))
# # abi_p = t(apply(ab_pr[3,,], 1, function(x) diff(x)))
# # abi_pl = t(apply(ab_pr[2,,], 1, function(x) diff(x)))
# # abi_pu = t(apply(ab_pr[4,,], 1, function(x) diff(x)))
# 
# # abi_p_m1 = t(apply(ab_pr[,,1], c(1), function(x) diff(x)))
# # abi_p_m2 = t(apply(ab_pr[,,2], c(1), function(x) diff(x)))
# # abi_p_m3 = t(apply(ab_pr[,,3], c(1), function(x) diff(x)))
# 
# abi_c = t(apply(ab_cr, 1, function(x) diff(x)))
# for (i in 1:ncol(abi_c)){
#   abi_c[,i] = abi_c[,i]/diff(census_years)[i]
# }
# 
# 
# # plot abi by plot
# ggplot() + geom_line(data=ab_all_site, aes(x=year, y=ab, colour=factor(type), lintype=factor(type)), size=1) + 
#   # geom_ribbon(data=ab_site_r, aes(x=year, ymin=ab25, ymax=ab975), alpha = 0.4, fill="lightgreen") + 
#   scale_shape_manual(values=c(19,1,15)) + 
#   # xlab("Year") + ylab('Stem biomass (Kg / m^2)') + 
#   xlab("Year") + ylab('Stem biomass (Mg / ha)') + 
#   facet_grid(site_id~.) +
#   scale_x_continuous(breaks=seq(min(years), max(years), by=5)) + 
#   labs(colour="Method",shape="Method") +
#   theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14)) + 
#   theme(legend.title=element_text(size=14), 
#         legend.text=element_text(size=14))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# abi_m_melt = data.frame(rep('rw'), melt(abi_m))
# colnames(abi_m_melt) = c('model', 'tree_id', 'year_id', 'abi')
# 
# abi_p_melt1 = data.frame(rep('all data'), melt(abi_p_m1))
# colnames(abi_p_melt1) = c('model', 'tree_id', 'year_id', 'abi')
# abi_p_melt2 = data.frame(rep('all data'), melt(abi_p_m2))
# colnames(abi_p_melt2) = c('model', 'tree_id', 'year_id', 'abi')
# abi_p_melt3 = data.frame(rep('one core'), melt(abi_p_m3))
# colnames(abi_p_melt3) = c('model', 'tree_id', 'year_id', 'abi')
# 
# abi_c_melt = data.frame(rep('census'), melt(abi_c))
# colnames(abi_m_melt) = c('tree_id', 'year_id', 'abi')
# 
# abi_all =rbind(data.frame(abi_m_melt, type=rep('Measured', nrow(abi_m_melt))), 
#                data.frame(abi_p_melt1, type=rep('Predicted', nrow(abi_p_melt))),
#                data.frame(abi_p_melt2, type=rep('Predicted', nrow(abi_p_melt))),
#                data.frame(abi_p_melt3, type=rep('Predicted', nrow(abi_p_melt))))
# abi_all = data.frame(abi_all, 
#                      site_id=tree_site_id[abi_all$tree_id], 
#                      year=years[abi_all$year_id])
# 
# abi_c_melt = data.frame(abi_c_melt, type=rep('C', nrow(abi_c_melt)))
# colnames(abi_c_melt) = c('tree_id', 'year_id', 'abi', 'type')
# abi_c_melt = data.frame(abi_c_melt, site_id=tree_site_id[abi_c_melt$tree_id], year=census_years[abi_c_melt$year_id])
# 
# 
# breaks = c(0, 10, 20, 100)
# 
# dbh_all$size = cut(dbh_all$dbh, breaks=breaks, labels=FALSE)
# dbh_all$size = factor(dbh_all$size)
# levels(dbh_all$size) = c('small', 'medium', 'large')
# 
# ab_all$size  = dbh_all$size
# 
# dbh_all$taxon = taxon[dbh_all$tree_id]
# ab_all$taxon  = dbh_all$taxon
# 
# dbh_all$group = taxaMatch[match(dbh_all$taxon, taxaMatch$taxon), 'group']
# ab_all$group  = taxaMatch[match(ab_all$taxon, taxaMatch$taxon), 'group']
# 
# abi_all = merge(abi_all, ab_all, by=c('tree_id', 'year_id', 'type', 'site_id'))
# abi_all$year=years[abi_all$year_id]
# 
# #####################################################################################################################################
# ## make some plots
# #####################################################################################################################################
# 
# dbh_all_site = aggregate(dbh  ~ site_id + year + size + type + model, dbh_all, sum, na.rm=TRUE)
# ab_all_site  = aggregate(ab  ~ site_id + year + size + type + model, ab_all, sum, na.rm=TRUE)
# 
# # plot total biomass by site by type
# ab_all_site  = aggregate(ab  ~ site_id + year + type + model, ab_all, sum, na.rm=TRUE)
# 
# # ab_site_r = reshape(ab_all_site[ab_all_site$type=='Predicted',c(1,2,4,5)], timevar="quant", idvar=c('year', 'site_id'), direction="wide")
# # colnames(ab_site_r)[3:7] = c("abmean", "absd", "ab25", "ab50", "ab975") 
# 
# ggplot() + geom_line(data=ab_all_site, aes(x=year, y=ab, colour=factor(type), lintype=factor(type)), size=1) + 
#   # geom_ribbon(data=ab_site_r, aes(x=year, ymin=ab25, ymax=ab975), alpha = 0.4, fill="lightgreen") + 
#   scale_shape_manual(values=c(19,1,15)) + 
#   # xlab("Year") + ylab('Stem biomass (Kg / m^2)') + 
#   xlab("Year") + ylab('Stem biomass (Mg / ha)') + 
#   facet_grid(site_id~.) +
#   scale_x_continuous(breaks=seq(min(years), max(years), by=5)) + 
#   labs(colour="Method",shape="Method") +
#   theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14)) + 
#   theme(legend.title=element_text(size=14), 
#         legend.text=element_text(size=14))
# 
# ggsave(paste0('figures/ab_all_mod_', suff, '.pdf'))
# 
# # plot abi by site
# abi_all_site  = aggregate(abi  ~ site_id + year + type + quant, abi_all, sum, na.rm=TRUE)
# 
# abi_site_r = reshape(abi_all_site[abi_all_site$type=='Predicted',c(1,2,4,5)], timevar="quant", idvar=c('year', 'site_id'), direction="wide")
# colnames(abi_site_r)[3:7] = c("abmean", "absd", "ab25", "ab50", "ab975") 
# 
# p <- ggplot() + geom_line(data=abi_all_site, aes(x=year, y=abi, group=type, colour=type), size=1, alpha=0.5) + 
#   geom_point(data=abi_all_site, aes(x=year, y=abi, colour=type, shape=type)) 
# p <- p + geom_ribbon(data=abi_site_r, aes(x=year, ymin=ab25, ymax=ab975), alpha=0.4, fill="lightpink")
# # p <- p + geom_segment(data=abi_c_site, aes(x=x, y=abi, xend=xend,  yend=abi), 
# #                       colour=cols[2], size=1)
# p <- p + facet_grid(site_id~.)
# p <- p + scale_shape_manual(values=c(19,1)) + 
#   xlab("Year") +  ylab('ABI (Mg / ha)') + 
#   # xlab("Year") + ylab('ABI (Kg / m^2)') + 
#   labs(colour="Method",shape="Method") + 
#   theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14)) #+ 
# print(p)
# 
# ggsave(paste0('figures/abi_by_site_', suff, '.pdf'))
# 
# 
# # plot abi by size
# abi_all_size  = aggregate(abi  ~ site_id + year + type + quant + size, abi_all, sum, na.rm=TRUE)
# 
# ggplot() + geom_line(data=abi_all_size[abi_all_size$quant=='mean',], aes(x=year, y=abi, colour=size, linetype=type), size=0.9) + 
#   # ylab('Aboveground biomass increment (Kg / m^2)') +
#   ylab('Aboveground biomass increment (Mg / ha)') +
#   facet_grid(site_id~.)
# ggsave(paste0('figures/abi_by_site_size_', suff, '.pdf'))
# 
# 
# # ab by taxon 
# ab_taxon = aggregate(ab ~ taxon + year + site_id + quant + type, ab_all, function(x) sum(x, na.rm=TRUE))
# 
# ggplot() + geom_line(data=ab_taxon[ab_taxon$quant=='mean',], aes(x=year, y=ab, colour=factor(taxon), linetype=type), size=0.9) + 
#   # ylab('Aboveground biomass increment (Kg / m^2)') +
#   ylab('Aboveground biomass (Mg / ha)') +
#   facet_grid(site_id~.)
# ggsave(paste0('figures/ab_by_taxon_', suff, '.pdf'))
# 
# # ab by pft 
# ab_group = aggregate(ab ~ group + year + site_id + quant + type, ab_all, function(x) sum(x, na.rm=TRUE))
# 
# ggplot() + geom_line(data=ab_group[ab_group$quant=='mean',], aes(x=year, y=ab, colour=factor(group_name[group]), linetype=type), size=0.9) + 
#   # ylab('Aboveground biomass increment (Kg / m^2)') +
#   ylab('Aboveground biomass (Mg / ha)') +
#   facet_grid(site_id~.)
# ggsave(paste0('figures/ab_by_species_group_', suff, '.pdf'))
# 
# 
# ggplot(data=ab_pft[ab_group$quant=='mean',]) + geom_line(data=ab_group[ab_group$quant=='mean',], aes(x=year, y=ab, colour=factor(site_id), linetype=type), size=0.9) + 
#   # ylab('Aboveground biomass increment (Kg / m^2)') +
#   ylab('Aboveground biomass (Mg / ha)') +
#   facet_grid(pft~., scales='free')
# 
# 
# # dbh_binned = aggregate(dbh ~ site_id + year + size + type + quant, dat, sum, na.rm=TRUE)
# # ab_binned  = aggregate(ab  ~ site_id + year + size + type + quant, dat, sum, na.rm=TRUE)
# 
# ###########################################################################################################################
# # for ann and mike
# ###########################################################################################################################
# ab_group$name = group_name[ab_group$group]
# 
# ab_group_sub    = ab_group[ab_group$type == 'Predicted', ]
# ab_group_sub$ab = ab_group_sub$ab/10 # change units to Kg/m2
# 
# saveRDS(ab_group_sub, file="data/lyford_ab_group_v1.RDS")
# 
# ###########################################################################################################################
# # garbage
# ###########################################################################################################################
# 
# # figure out which don't have measurements, they are so different
# dbh_p_smooth[1,,52]
# dbh_m[,52]
# 
# sum(dbh_p_smooth[1,,52],na.rm=TRUE)
# sum(dbh_m[,52],na.rm=TRUE)
# 
# sum_dbh_p = colSums(dbh_p_smooth[1,,], na.rm=TRUE)
# sum_dbh_m = colSums(dbh_m, na.rm=TRUE)
# 
# plot(years, sum_dbh_p, ylim=c(0,6000))
# lines(years, sum_dbh_m)
# 
# 
# miss = apply(dbh_m, 1, function(x) all(is.na(x)))
# sum(miss)
# 
# saplings = which(miss & (dbh_p_smooth[1,,52]<10))
# med_miss = which(miss & (dbh_p_smooth[1,,52]>10) & (dbh_p_smooth[1,,52]<20))
# big_miss = which(miss & (dbh_p_smooth[1,,52]>=20))
# 
# dbh[dbh$stat_id %in% big_miss, ]
# dbh[which((dbh$stat_id %in% big_miss) & (dbh$year == 2011)), ]
# 
# pdf(file='figures/missing_3829.pdf')
# plot(census$xsite, census$ysite, asp=1, main='Missing tree 3829')
# points(census[census$census_id == 3829,'xsite'], census[census$census_id == 3829,'ysite'], col='blue', pch=19)
# dev.off()
# 
# dbh[dbh$stat_id %in% med_miss, ]
# dbh[which((dbh$stat_id %in% med_miss) & (dbh$year == 2011)), ]
# 
# 
# plot(census[census$census_id %in% dbh$census_id[which((dbh$stat_id %in% big_miss) & (dbh$year == 2011))],'xsite'], 
#      census[census$census_id %in% dbh$census_id[which((dbh$stat_id %in% big_miss) & (dbh$year == 2011))],'ysite'], asp=1)
# 
# plot(census[census$census_id %in% dbh$census_id[which((dbh$stat_id %in% med_miss) & (dbh$year == 2011))],'xsite'], 
#      census[census$census_id %in% dbh$census_id[which((dbh$stat_id %in% med_miss) & (dbh$year == 2011))],'ysite'], asp=1)
# 
# 
# ## get census trees only in plots
# med_miss_cid = unique(dbh[which(dbh$stat_id %in% med_miss), 'census_id'])
# 
# plotRadius = 20
# ftPerMeter = 3.28084
# check <- vector()
# for(i in 1:3) {
#   dist <- rdist(centroids[i, 2:3], census[, c('xsite', 'ysite')])
#   check <- c(check, census[which((census$census_id %in% med_miss_cid) & (dist > 10*ftPerMeter) & 
#                                    (dist < plotRadius*ftPerMeter)), 'census_id'])
#   #   outer[[i]]$site <- i
#   #   outer[[i]]$dist_census <- dist[which((census$census_id %in% med_miss_cid) & (dist > 10*ftPerMeter) & (dist < plotRadius*ftPerMeter))]
# }
# 
# med_miss_cid[which(!(med_miss_cid %in% check))]
# 
# # from sites 2,2,3
# dist <- rdist(centroids[, 2:3], census[which(census$census_id %in%  med_miss_cid[which(!(med_miss_cid %in% check))]), c('xsite', 'ysite')])
# 
# # idx_missing = which(!is.na(dbh_p[which(is.na(dbh_m[,52])),52]))
# 
# # i=2
# # 
# # count = 0
# # for (i in 1:length(idx_missing)){
# #   if (all(is.na(dbh_m[idx_missing[i],]))){
# #     print(idx_missing[i])
# #     print(dbh_p[idx_missing[i],52])
# #   } else {
# #     count = count + 1
# #   }
# # }
# 
# # count
# idx_missing = which(is.na(dbh_m[,52]))
# dbh_p_smooth[idx_missing, 52]
# 
# for (k in 1:N_sites){
#   meas = dbh_m[idx_missing[which(tree_site_id[idx_missing] == k)],]
#   miss = apply(meas, 1, function(x) all(is.na(x)))
#   pred = dbh_p[idx_missing[which(tree_site_id[idx_missing] == k)],]
#   missp = apply(pred, 1, function(x) all(is.na(x)))
#   pred[miss,52]
# }
# 
# idx_m = which(!is.na(dbh_m[which(tree_site_id==1),52]))
# length(idx_m)
# 
# idx_p = which(!is.na(dbh_p[which(tree_site_id==1),52]))
# length(idx_p)
# 
# sum(dbh_p[which(tree_site_id==1),52]>10, na.rm=TRUE)
# 
# meas = apply(dbh_m, 1, function(x) any(!is.na(x)))
# # 133
# # only 37 trees measured from site 1? correct?
# 
# 
# 
# sum(dbh_p[,52][idx_m],na.rm=TRUE)
# sum(dbh_m[,52][idx_m],na.rm=TRUE)
# 
# ########################################################################################################################################
# ## ab contributions by taxon
# ########################################################################################################################################
# 
# ab_all_raw =rbind(data.frame(ab_c_melt, type=rep('C', nrow(ab_c_melt))), 
#                   data.frame(ab_m_melt, type=rep('M', nrow(ab_m_melt))), 
#                   data.frame(ab_melt, type=rep('P', nrow(ab_melt))))
# 
# # ab_taxon = aggregate(ab ~ taxon + year + site_id + type, ab_all_raw, function(x) sum(x, na.rm=TRUE))
# ab_taxon = aggregate(ab ~ taxon + year + type + site_id, ab_all_raw, function(x) sum(x, na.rm=TRUE))
# 
# ab_taxon_per = ab_taxon
# for (i in 1:nrow(ab_taxon)) {
#   total = sum(ab_taxon[which( (ab_taxon_per$year[i] == ab_taxon$year) & 
#                                 (ab_taxon_per$type[i] == ab_taxon$type) & 
#                                 (ab_taxon_per$site_id[i] == ab_taxon$site_id) ), 'ab'])
#   ab_taxon_per$ab[i] = ab_taxon_per$ab[i] / total
# }
# 
# # check it!
# foo=aggregate(ab ~ year + type + site_id, ab_taxon_per, sum)
# foo
# 
# ab_taxon_per = ab_taxon_per[ab_taxon_per$year %in% census_years,]
# ab_taxon_per$taxon = taxaMatch$species[ab_taxon_per$taxon]
# 
# any(is.na(ab_taxon_per$ab))
# 
# ab_taxon_per = ab_taxon_per[which(ab_taxon_per$type == 'C'), ]
# 
# p <- ggplot(data=ab_taxon_per) + geom_bar(data=ab_taxon_per, aes(x=factor(1), y=ab, fill=factor(taxon)), stat='identity') 
# p <- p + facet_grid(year~site_id)
# p <- p + coord_polar(theta="y") 
# p <- p + scale_fill_brewer(palette='Spectral')
# print(p)
# 
# p <- ggplot(data=ab_taxon_per) + geom_bar(data=ab_taxon_per, aes(x=year, y=ab, fill=taxon), stat='identity', position='dodge') 
# p <- p + facet_grid(site_id~.)
# # p <- p + coord_polar(theta="y") 
# p <- p + scale_fill_brewer(palette='Spectral')
# print(p)
# 
# 
# taxaMatch = data.frame(taxaMatch, group=c(2, 5, 4, 4, 5, 3, 5, 5, 1, 5))
# ab_taxon_per$taxon = taxaMatch[match(ab_taxon_per$taxon, taxaMatch$species), 4]
# ab_taxon_agg = aggregate(ab ~ taxon + year + type + site_id, ab_taxon_per, sum)
# 
# ab_taxon_agg$taxon <- factor(ab_taxon_agg$taxon, levels=unique(ab_taxon_agg$taxon))
# levels(ab_taxon_agg$taxon) <- c('Q. rubra', 'A. rubrum', 'F. gradifolia', 'Betula spp.', 'Other')
# 
# p <- ggplot(data=ab_taxon_agg) + geom_bar(data=ab_taxon_agg, aes(x=factor(1), y=ab, fill=factor(taxon)), stat='identity') 
# p <- p + facet_grid(year~site_id)
# p <- p + coord_polar(theta="y") 
# p <- p + scale_fill_brewer(palette='Set1')
# print(p)
# 
# p <- ggplot(data=ab_taxon_agg) + geom_bar(data=ab_taxon_agg, aes(x=year, y=ab, fill=factor(taxon)), stat='identity', position='dodge') 
# p <- p + facet_grid(site_id~.)
# # p <- p + coord_polar(theta="y") 
# p <- p + scale_fill_brewer(palette='Spectral')
# print(p)
# 
# ab_taxon_agg = ab_taxon_agg[ab_taxon_agg$year %in% c(1969, 2011),]
# p <- ggplot(data=ab_taxon_agg) + geom_bar(data=ab_taxon_agg, aes(x=taxon, y=ab, fill=factor(year)), stat='identity', position='dodge') 
# p <- p + facet_grid(site_id~.)
# # p <- p + coord_polar(theta="y") 
# p <- p + scale_fill_brewer(palette='Set1')
# print(p)
# 
# 
# ab_taxon$taxon = taxaMatch[match(ab_taxon$taxon, taxaMatch$taxon), 4]
# ab_taxon_agg = aggregate(ab ~ taxon + year + type + site_id, ab_taxon, sum)
# # ab_taxon_count = aggregate(ab ~ taxon + year + type + site_id, ab_taxon, NROW)
# 
# ab_all_raw$taxon = taxaMatch[match(ab_all_raw$taxon, taxaMatch$taxon), 4]
# ab_taxon_count = aggregate(ab ~ taxon + year + type + site_id, ab_all_raw, NROW)
# 
# ab_taxon_agg$taxon <- factor(ab_taxon_agg$taxon, levels=unique(ab_taxon_agg$taxon))
# levels(ab_taxon_agg$taxon) <- c('Q. rubra', 'A. rubrum', 'F. grandifolia', 'Betula spp.', 'Other')
# 
# ab_taxon_count$taxon <- factor(ab_taxon_count$taxon, levels=unique(ab_taxon_count$taxon))
# levels(ab_taxon_count$taxon) <- c('Q. rubra', 'A. rubrum', 'F. grandifolia', 'Betula spp.', 'Other')
# 
# # p <- ggplot(data=ab_taxon_agg) + geom_bar(data=ab_taxon_agg, aes(x=factor(1), y=ab, fill=factor(taxon)), stat='identity') 
# # p <- p + facet_grid(year~site_id)
# # p <- p + coord_polar(theta="y") 
# # p <- p + scale_fill_brewer(palette='Set1')
# # print(p)
# 
# p <- ggplot(data=ab_taxon_agg) + geom_bar(data=ab_taxon_agg, aes(x=year, y=ab, fill=factor(taxon)), stat='identity', position='dodge') 
# p <- p + facet_grid(site_id~.)
# # p <- p + coord_polar(theta="y") 
# p <- p + scale_fill_brewer(palette='Spectral')
# print(p)
# 
# ab_taxon_agg = ab_taxon_agg[ab_taxon_agg$year %in% c(1969, 2011),]
# ab_taxon_agg$ab = ab_taxon_agg$ab * 10 / 13^2 / pi 
# ab_taxon_count = ab_taxon_count[ab_taxon_count$year %in% c(1969, 2011),]
# 
# ab_taxon_agg_c = data.frame(ab_taxon_agg[which(ab_taxon_agg$type == 'C'),], count=ab_taxon_count[which(ab_taxon_count$type == 'C'),]$ab)
# # ab_taxon_count_c = ab_taxon_count[which(ab_taxon_count$type == 'C'),]
# ab_taxon_agg_c = rbind(ab_taxon_agg_c, 
#                        data.frame(taxon='Other', year=1969, type='C', site_id=1, ab=0, count=0),
#                        data.frame(taxon='F. grandifolia', year=1969, type='C', site_id=2, ab=0, count=0))
# 
# dodgewidth <- position_dodge(width=0.9)
# p <- ggplot(data=ab_taxon_agg_c, aes(x=taxon, y=ab, fill=factor(year))) + 
#   geom_bar(data=ab_taxon_agg_c, stat='identity', position='dodge') 
# p <- p + facet_grid(site_id~.)
# # p <- p + coord_polar(theta="y") 
# # p <- p + scale_fill_brewer(palette='Set1')
# p <- p + geom_text(aes(label=count), position=position_dodge(width=0.9), vjust=-0.25)
# p <- p + ylim(0, 275)
# # xlab("Year") + ylab('Stem biomass (Mg / ha)') + 
# p <- p + xlab("") + ylab('') + 
#   # theme(panel.grid = element_blank()) + 
#   labs(fill="Year") +
#   theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14)) + 
#   theme(strip.background = element_blank(),
#         strip.text.y = element_blank(), 
#         legend.title=element_text(size=14), 
#         legend.text=element_text(size=14))
# # p <- p + stat_bin(geom="text", position= dodgewidth, aes(label=count), vjust=-1)
# print(p)
# ggsave('ab_taxon_1969_2011.pdf')
# ggsave('ab_taxon_1969_2011.png')
# ggsave('ab_taxon_1969_2011.eps')
# 
# ########################################################################################################################################
# ## fading record; number dying, etc.
# ########################################################################################################################################
# 
# # number of trees contributing over time
# dbh_c
# dbh_m
# 
# counts_p = apply(dbh_p[,,keep], 2, function(x) sum(!is.na(x)))
# counts_m = apply(dbh_m, 2, function(x) sum(!is.na(x)))
# 
# counts_c = rep(NA, N_years)
# counts_c[match(census_years, years)] = apply(dbh_c, 2, function(x) sum(!is.na(x)))
# 
# counts_dat = melt(data.frame(years, counts_p, counts_m, counts_c), id.vars='years')
# 
# ggplot(data=counts_dat) + geom_point(data=counts_dat, aes(x=years, y=value, colour=variable), size=3) + 
#   ylab('Number of trees') + xlab('Years')
# 
# # by site
# 
# dbh_p_melt = melt(dbh_p)
# colnames(dbh_p_melt) <- c('tree_id', 'year', 'dbh')
# dbh_p_melt$year = years[dbh_p_melt$year]
# dbh_p_melt = data.frame(dbh_p_melt, site_id = tree_site_id[dbh_p_melt$tree_id], type=rep('predicted', nrow(dbh_p_melt)))
# 
# dbh_m_melt = melt(dbh_m)
# colnames(dbh_m_melt) <- c('tree_id', 'year', 'dbh')
# dbh_m_melt$year = years[dbh_m_melt$year]
# dbh_m_melt = data.frame(dbh_m_melt, site_id = tree_site_id[dbh_m_melt$tree_id], type=rep('measured', nrow(dbh_m_melt)))
# 
# dbh_c_melt = melt(dbh_c)
# colnames(dbh_c_melt) <- c('tree_id', 'year', 'dbh')
# dbh_c_melt$year = census_years[dbh_c_melt$year]
# dbh_c_melt = data.frame(dbh_c_melt, site_id = tree_site_id[dbh_c_melt$tree_id], type=rep('census', nrow(dbh_c_melt)))
# 
# # don't forget the saplings!
# dbh_sap_melt = data.frame(tree_id = sapling_tree_id, year=years[sapling_year_id], dbh = rep(5, N_saplings), 
#                           site_id = tree_site_id[sapling_tree_id], type=rep('census', N_saplings))
# dbh_sap_melt$year[dbh_sap_melt$year == 1992] = 1991
# dbh_sap_melt$year[dbh_sap_melt$year == 1987] = 1991
# 
# dbh_all = rbind(dbh_p_melt, dbh_m_melt, dbh_c_melt, dbh_sap_melt)
# counts_dat = aggregate(dbh ~ year + site_id + type, dbh_all, function(x) sum(!is.na(x)))
# 
# ggplot(data=counts_dat) + geom_point(data=counts_dat, aes(x=year, y=dbh, colour=type), size=3) + 
#   ylab('Number of trees') + xlab('Years') + facet_grid(site_id~.)
# ggsave(paste0(file='figures/number_trees', suff, '.pdf'))
# 
# ########################################################################################################################################
# ## dbh distributions
# ########################################################################################################################################
# 
# dbh_sub = dbh_all[dbh_all$year %in% census_years,]
# counts_sub = counts_dat[counts_dat$year %in% census_years,]
# # dbh_sub = dbh_all[dbh_all$year %in% census_years,]
# 
# ggplot(data=dbh_sub) + geom_histogram(data=dbh_sub, aes(x=dbh), binwidth=2) +facet_grid(year~type)
# ggsave('figures/dbh_hists_census_years.pdf')
# 
# ggplot(data=dbh_sub) + geom_density(data=dbh_sub, aes(x=dbh, fill=type), alpha=0.2) +facet_grid(year~type)
# ggsave('figures/dbh_density_census_years.pdf')
# 
# # what about for the live vs dead trees
# 
# # note last observed
# last_ti
# 
# # record year of death and dbh
# mort = array(NA, c(0, 5))
# for (i in 1:N_trees) {
#   tree = i
#   dbh_tree = dbh_sub[which( (dbh_sub$tree_id == tree) & (dbh_sub$type == 'predicted') ),]
#   dbh_tree = dbh_tree[!is.na(dbh_tree$dbh), ]
#   if (max(dbh_tree$year) < 2011) { 
#     idx_dead = which.max(dbh_tree$year)
#     mort = rbind(mort, dbh_tree[idx_dead,])
#   }
# }
# 
# live = dbh_sub[!(dbh_sub$tree_id %in% mort$tree_id),]
# live$type = rep('live', nrow(live))
# mort$type = rep('dead', nrow(mort))
# 
# status = rbind(live, mort)
# 
# ggplot(data=status) + geom_histogram(data=status, aes(x=dbh, fill=type), binwidth=2) #+facet_grid(year~type)
# 
# ggplot(data=status) + geom_histogram(data=status, aes(x=dbh, fill=type), binwidth=2) +facet_grid(year~type)
# ggsave('figures/dbh_hists_status.pdf')
# 
# ggplot(data=status) + geom_density(data=status, aes(x=dbh, fill=type), alpha=0.2) +facet_grid(year~.)
# ggsave('figures/dbh_density_status.pdf')
# 
# # death by taxon
# status = data.frame(status, taxon=taxaMatch[taxon[status$tree_id], 'species'])
# 
# ggplot(data=status) + geom_histogram(data=status, aes(x=dbh, fill=type), binwidth=2) +facet_grid(taxon~year)
# ggsave('figures/dbh_hists_by_status_taxon_year.pdf')
# 
# # ggplot(data=status) + geom_histogram(data=status, aes(x=dbh, fill=type), binwidth=2) +facet_grid(~year)
# 
# ggplot(data=status) + geom_histogram(data=status, aes(x=dbh, fill=type), binwidth=2) + facet_wrap(~taxon, nrow=2)
# ggsave('figures/dbh_hists_by_status_taxon.pdf')
# 
# # for the poster
# levels(status$taxon) <- list(Acer = c('ACSA', 'ACRU'), Betula=c('BEAL', 'BELE'), Castanea='CADE', Fagus='FAGR', 
#                              Pine='PIST', Quercus='QURU', Tsuga='TSCA', Other='HAVI')
# status$type <- factor(status$type)
# levels(status$type) <- list(Dead='dead', Live='live')
# colnames(status[5]) = 'Status'
# ggplot(data=status) + geom_histogram(data=status, aes(x=dbh, fill=type), binwidth=4, alpha=0.4) + facet_grid(taxon~.) + 
#   xlab('Dbh (cm)') + ylab('Number of trees') #+ scale_fill_manual(values=c('red', 'blue'))
# ggsave('figures/dbh_hists_by_status_taxon.pdf')
# 
# ggplot(data=status) + geom_histogram(data=status, aes(x=dbh, fill=type), binwidth=4, alpha=0.9) + facet_wrap(~taxon, nrow=2) + 
#   xlab('Dbh (cm)') + ylab('Number of trees') + scale_fill_manual(values=c('#fc8d59', '#67a9cf')) + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14)) + 
#   theme(strip.text.x = element_text(size=14), 
#         legend.title=element_text(size=14), 
#         legend.text=element_text(size=14)) +  labs(fill="Status")  + xlab("") + ylab('') 
# ggsave('figures/dbh_hists_by_status_taxon.pdf')
# 
# ggplot(data=status) + geom_histogram(data=status, aes(x=dbh, y=..density.., fill=type), binwidth=2, alpha=0.6) + facet_wrap(~taxon, nrow=2) + 
#   xlab('Dbh (cm)') + ylab('Density of trees') + scale_fill_manual(values=c('#fc8d59', '#67a9cf'))
# ggsave('figures/dbh_density_by_status_taxon.pdf', width=14, height=6)
# 
# # ggplot(data=status) + geom_histogram(data=status, aes(x=dbh), binwidth=4, alpha=0.4) + facet_grid(type~taxon) + 
# #   xlab('Dbh (cm)') + ylab('Number of trees') #+ scale_fill_manual(values=c('red', 'blue'))
# # 
# # ggplot(data=status) + geom_boxplot(data=status, aes(x=taxon, y=dbh, fill=type), alpha=0.4)# + facet_grid(taxon~.) + 
# #   xlab('Dbh (cm)') + ylab('Number of trees') #+ scale_fill_manual(values=c('red', 'blue'))
# 
# # ggplot(data=status) + geom_density(data=status, aes(x=dbh, fill=type), alpha=0.2) + facet_grid(taxon~.) + 
# # xlab('Dbh (cm)') + ylab('Number of trees')
# 
# # ggplot(data=status) + geom_density(data=status, aes(x=dbh, fill=type), alpha=0.2) +facet_grid(taxon~., nrow=2)
# 
# # death by location
# xloc = census$xsite[match(status$tree_id, census$stat_id)]
# yloc = census$ysite[match(status$tree_id, census$stat_id)]
# status = data.frame(status, xloc, yloc)
# # check that we are getting the right trees
# # par(mfrow=c(1,1))
# # plot(status$xloc[which(status$type == 'live')], status$yloc[which(status$type == 'live')], asp=1, col="lightgrey")#, xlim=c(-500,500))
# # points(status$xloc[which(status$type == 'dead')], status$yloc[which(status$type == 'dead')], asp=1, col="blue")
# # points(centroids$x, centroids$y, col="red", pch='X')
# # 
# # for (i in 1:nPlots){
# #   points(census[[i]]$xsite, census[[i]]$ysite, pch=19, col='blue')
# #   draw.circle(centroids$x[i], centroids$y[i], plotRadius*ftPerMeter)
# # }
# 
# ggplot(status) + geom_point(data=status, aes(x=xloc, y=yloc, colour=type)) + coord_fixed() + facet_grid(year~site_id)#, nrow=5)#, scales="free")
# ggsave(file="figures/spatial_patterns_status.pdf")
# 
# # how many died at each site; should match count data
# nmort = aggregate(tree_id ~ year + site_id, mort, function(x) length(unique(x)))
# nmort = rbind(nmort, c(1969, 1, 0))
# 
# ggplot(data=nmort) + geom_bar(data=nmort, 
#                               aes(x=year, y=tree_id, fill=factor(site_id)), stat='identity', width=0.5, position='dodge') + 
#   xlab('Year') + ylab('Number of tree deaths') + scale_x_discrete(labels=c("1969 - 1975", "1975 - 1991", "1991 - 2001")) + 
#   scale_y_continuous(breaks=0:12) + 
#   guides(fill=guide_legend(title='Site'))#+ facet_grid(site_id~.)
# ggsave('figures/mort_counts_by_site.pdf')
# 
# # ggplot(data=nmort) + geom_point(data=nmort, aes(y=year, x=tree_id, colour=factor(site_id)), size=4)# + facet_grid(site_id~.)
# 
# ########################################################################################################################################
# ## any ring widths from dead trees?
# ########################################################################################################################################
# 
# 
# ########################################################################################################################################
# ## why are measured biomass estimates so much higher than predicted?? expect opposite....
# ########################################################################################################################################
# 
# # dbh differences? is this a model problem?
# dbh_diff = dbh_p - dbh_m
# 
# par(mfrow=c(1,1))
# plot(c(0,0), xlim=c(year_start,year_end), ylim=c(min(dbh_diff, na.rm=TRUE), max(dbh_diff, na.rm=TRUE)), type='n')
# for (i in 1:N_trees){
#   points(years, dbh_diff[i,])
# }
# 
# pdf('figures/dbh_diffs.pdf', width=10,height=8)
# par(mfrow=c(1,1))
# for (i in 1:N_trees){
#   plot(c(0,0), xlim=c(year_start,year_end), ylim=c(min(dbh_p, na.rm=TRUE), max(dbh_p, na.rm=TRUE)), type='n')
#   points(years, dbh_p[i,])
#   points(years, dbh_m[i,], col='blue', pch=19)
# }
# dev.off()
# 
# ########################################################################################################################################
# ## plot census data and computed dbh
# ########################################################################################################################################
# 
# # plot dbh_m
# 
# pdf('figures/dbh_reconstruct.pdf', width=10,height=8)
# for (i in 1:N_trees){
#   tree = trees[i]
#   tree_idx = which(x2tree == tree)
#   tree_years = x2year[tree_idx]
#   tree_years = years[tree_years]
#   
#   idx_dbh = which(dbh_tree_id == tree)
#   yrs = dbh_year_id[idx_dbh]
#   
#   D_cen = exp(logDobs[idx_dbh])
#   D_calc = dbh_m[tree,]
#   
#   plot(c(0,0), xlim=c(min(tree_years),max(tree_years)), ylim=c(0,max(c(D_cen,D_calc), na.rm=TRUE)),#c(min(c(D_dat,D_mu)), max(c(D_dat,D_mu))), 
#        xlab='Year', ylab='dbh', type='n', main=paste0('Tree ', tree))
#   
#   lines(years, D_calc, col='blue')#, lty=2)
#   points(years[yrs], D_cen, pch=19, col='black')  
# }
# dev.off()
# 
# 
# ########################################################################################################################################
# ## abi take 2
# ########################################################################################################################################
# 
# ab_both_wide <- reshape(ab_all, 
#                         timevar = "year",
#                         idvar = c("site_id", "type"),
#                         direction = "wide")
# 
# abi_both = t(apply(ab_both_wide[,3:ncol(ab_both_wide)], 1, function(x) diff(x) ))
# abi_both = cbind(ab_both_wide[,1:2], abi_both)
# 
# abi_melt = melt(abi_both, id.vars=c('site_id', 'type'))
# colnames(abi_melt)[3] = 'year'
# colnames(abi_melt)[4] = 'abi'
# abi_melt$year = substr(abi_melt$year, 4, 8)
# 
# abi_melt = abi_melt[abi_melt$type != 'Census',]
# 
# abi_melt$type <- factor(abi_melt$type)
# levels(abi_melt$type) <- c('Measured', 'Predicted')
# 
# ggplot(data=abi_melt) + geom_point(data=abi_melt, aes(x=year, y=abi, colour=factor(site_id), shape=type)) + 
#   scale_shape_manual(values=c(19,3)) + xlab("Year") + ylab('Stem biomass increment (Mg / ha)') + facet_grid(site_id~.)
# 
# ggplot(data=abi_melt) + geom_line(aes(x=year, y=abi, group=type, colour=type)) + 
#   geom_point(data=abi_melt, aes(x=year, y=abi, colour=type, shape=type)) +
#   scale_shape_manual(values=c(19,3)) + xlab("Year") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
#   ylab('Stem biomass increment (Mg / ha)') + facet_grid(site_id~.)
# ggsave('abi_both.pdf')
# 
# ##############################################################################################################
# # Plot centers: 42.53065 N, -72.18346 W (LF1); 42.53128 N, -72.18271 W (LF2); 42.53008 N, -72.18246 W (LF3)
# 
# centroids
# 
# for (i in 1:nrow(treeMeta)){
#   
#   plot_zero = centroids[treeMeta[i, 'site'],c('x', 'y')]
#   
#   az = (90 - treeMeta[i, 'azimuth']) * pi / 180
#   r  = ftPerMeter * treeMeta[i, 'distance']
#   
#   treeMeta[i, 'x'] = r * cos(az) + plot_zero[1]
#   treeMeta[i, 'y'] = r * sin(az) + plot_zero[2]
#   
# }
# 
# pdf(file='figures/missing_trees.pdf')
# par(mfrow=c(1,1))
# plot(census$xsite, census$ysite, col='grey', asp=1)
# points(treeMeta$x, treeMeta$y, col='lightblue', pch=19)
# points(census[census$census_id %in% dbh$census_id[which((dbh$stat_id %in% big_miss) & (dbh$year == 2011))],'xsite'], 
#        census[census$census_id %in% dbh$census_id[which((dbh$stat_id %in% big_miss) & (dbh$year == 2011))],'ysite'], 
#        col='darksalmon', pch=19)
# text(census[census$census_id %in% dbh$census_id[which((dbh$stat_id %in% big_miss) & (dbh$year == 2011))],'xsite'], 
#      census[census$census_id %in% dbh$census_id[which((dbh$stat_id %in% big_miss) & (dbh$year == 2011))],'ysite'], 
#      labels=census[census$census_id %in% dbh$census_id[which((dbh$stat_id %in% big_miss) & (dbh$year == 2011))],'census_id'],
#      pos=c(3, 1, 3, 1, 2, 1, 3, 4), cex=1)
# dev.off()
# 
# census[census$census_id %in% dbh$census_id[which((dbh$stat_id %in% big_miss) & (dbh$year == 2011))], c('census_id', 'dbh11')]
