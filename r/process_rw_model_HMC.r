library(PEcAn.allometry)
library(ggplot2)
library(reshape2)
library(abind)

run_predict = TRUE
sample_allom = FALSE

figures_dir = 'NOCOVAR'
if (!file.exists(figures_dir)){
  dir.create(paste0("figures/", figures_dir))
}

source('config_HMC')

fname_data = paste0('tree_data_HMC_', dvers)
load(file=paste0('data/dump/', fname_data, '.rdata'))
fname_data = paste0('tree_data_HMC_no_census_', dvers)
load(file=paste0('data/dump/', fname_data, '.rdata'))

fnames = c('ring_model_t_pdbh_HMC_NOCOVAR', 'ring_model_t_pdbh_nc_HMC_NOCOVAR_sigd')
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
            BELU = data.frame(spcd=372), # FIXME, subbing in for BELU
            BEPA = data.frame(spcd=375,acronym="BEPA"),
            OSVI = data.frame(spcd=701,acronym="OSVI"),
            THOC = data.frame(spcd=241,acronym="THOC"),
            POGR = data.frame(spcd=743,acronym="POGR"),
            TIAM = data.frame(spcd=951,acronym="TIAM"),
            TSCA = data.frame(spcd=261,acronym="TSCA")
)
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

taxaMatch$group = c(1, 1, 4, 4, 6, 2, 5, 3, 7)
group_name      = c('Acer', 'Aspen',  'Basswood', 'Betula', 'Cedar', 'Ostrya', 'Tsuga')
taxaMatch

niter = dim(post[[1]])[1]
burn = 1800
keep = niter - burn

nmodels = length(fnames)

x2idx_p = match(x2tree_p, unique(x2tree_p))

taxon = taxaMatch[match(dbh[match(trees, dbh$stat_id), 'taxon'], taxaMatch$species),'taxon']
trees_p = sort(unique(x2tree_p))
taxon_p = taxon[match(trees_p, trees)]
pft_list_p = sort(unique(pft_list[taxon_p]))


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
  
  saveRDS(dbh_p_1, file=paste0('allom/dbh_p_1_',location, '_', mvers, '.rds'))
  saveRDS(ab_p_1, paste0('allom/ab_p_1_', location, '_', mvers, '.rds'))
  
  trees_p = sort(unique(x2tree_p))
  # taxon_p = taxon[trees_p]
  dbh_2 = build_dbh_p(post[[2]], x2year_p, x2tree_p, N_years, trees_p, keep, taxaMatch, taxon_p)
  dbh_p_2 = dbh_2$dbh_p
  tree_order_2 = dbh_2$tree_order
  pft_order_2  = dbh_2$pft_order
  
  # #FIXME: is this the problem?
  # pft_list_p = toupper(sort(unique(as.vector(taxaMatch[match(taxon_p, taxaMatch$taxon),'species']))))
  
  ab_p_2 = build_ab_p(allom.fit, dbh_p_2, N_trees_p, N_years, keep, pft_list_p)
  ab_p_2 = ab_p_2[sort(tree_order_2, index.return=TRUE)$ix,,]
  
  dbh_p_2_org = abind(dbh_p_2, along=1)
  dbh_p_2 = dbh_p_2_org[sort(tree_order_2, index.return=TRUE)$ix,,]
  
  saveRDS(dbh_p_2, paste0('allom/dbh_p_2_', location, '_', mvers, '.rds'))
  saveRDS(ab_p_2, paste0('allom/ab_p_2_', location, '_', mvers, '.rds'))
  
} else {
  dbh_p_1 = readRDS( paste0('allom/dbh_p_1_', location, '_', mvers, '.rds'))
  ab_p_1  = readRDS(paste0('allom/ab_p_1_', mvers, '.rds'))
  
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

dimnames(dbh_p_1)[[1]] <- seq(1, 279)
dbh_melt = melt(dbh_p_1)
colnames(dbh_melt) = c("tree", "year", "iter", "dbh")
dbh_melt$year = years[dbh_melt$year]

saveRDS(dbh_melt, 'data/HMC_DBH_iterations.RDS') 

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
ggsave('figures/NOCOVAR/sum_dbh_by_plot_HMC_v4.0.pdf')


######################################################################################################################################
## smooth out death
######################################################################################################################################

# get the census years
census_years = as.numeric(sort(unique(dbh$year)))
idx_census   = which(years %in% census_years)
N_census_years = length(census_years)

# # for now don't smooth the without census data...
# smooth_death <- function(dbh_p, ab_p, last_time_data, last_time, N_trees, years, keep){
#   
#   last_time_data[last_time_data == 1992] = 1991
#   idx_adjust = which(last_time_data < 2011)
#   
#   dbh_p_smooth = dbh_p
#   ab_p_smooth  = ab_p
#   
#   for (i in 1:N_trees) {
#     print(i)
#     
#     idx_last = which(years == last_time_data[i])
#     
#     if (idx_last >= 52) {
#       next
#     } else if (last_time_data[i] < last_time[i]) {
#       print('Adjusting dbh and ab due to death!')
#       
#       sample_int =  last_time[i] - last_time_data[i]
#       
#       for (k in 1:keep){
#         mort_year = sample(sample_int,1)
#         years_na = seq((last_time_data[i]+mort_year), last_time[i])
#         idx_na   = which(years %in% years_na)  
#         
#         
#         dbh_p_smooth[i,idx_na,k] = rep(NA, length(idx_na))
#         ab_p_smooth[i,idx_na,k] = rep(NA, length(idx_na))
#       }
#     }
#   }
#   return(list(dbh_p_smooth=dbh_p_smooth, ab_p_smooth=ab_p_smooth))
# }
# 
#   
#   
# smooth_1 = smooth_death(dbh_p_1, ab_p_1, last_time_data, last_time, N_trees, years, keep)
# dbh_p_1 = smooth_1$dbh_p_smooth
# ab_p_1 = smooth_1$ab_p_smooth

# for now take mean values
# dbh_p_smooth = apply(dbh_p_smooth, c(1,2,4), mean, na.rm=TRUE)
# ab_p_smooth = apply(ab_p_smooth, c(1,2,4), mean, na.rm=TRUE)
# dbh_p_smooth = apply(dbh_p_smooth, c(1,2,4), function(x) quantile(x, probs=c(0.5), na.rm=TRUE))
# ab_p_smooth  = apply(ab_p_smooth, c(1,2,4),  function(x) quantile(x, probs=c(0.5), na.rm=TRUE))

######################################################################################################################################
## get dbh and ab from increments and most recent dbh measurement 
######################################################################################################################################

N_samples = 200#00
dbh_m = array(NA, c(N_trees, N_years))
ab_m = array(NA, c(N_trees, N_years, N_samples))

for (i in 1:N_trees) {
  print(i)
  
  tree = trees[i]
  tree_idx = which(trees == tree)
  pft  = as.vector(taxaMatch[which(taxaMatch$taxon == taxon[tree_idx]),'species'] )
  
  if (length(which(m2tree_a == tree)) == 0) {
    print(paste0('No increments for tree ', tree, ' !'))
    
    dbh_idx   = which(dbh$stat_id == tree)
    dbh_years = dbh[dbh_idx, 'year']
    
    dbh_dat_tree = dbh[dbh_idx,]
    
    if (max(dbh_dat_tree$dbh) < 10) {
      print(paste0('Tree dbh is ', max(dbh_dat_tree$dbh), '; set ab to zero!'))
    }
    
    # dbh_m[tree, ] = rep(0, N_years)
    dbh_m[tree_idx, ] = rep(NA, N_years)
    
  } else {
    
    # for now use average increments
    incr_tree  = exp(logXobs_a[which(m2tree_a == tree)])
    incr_years = years[m2ti_a[which(m2tree_a == tree)]]
    
    dbh_idx   = which(dbh$stat_id == tree)
    dbh_years = dbh[dbh_idx, 'year']
    
    pdbh_idx = which(pdbh$stat_id == tree)
    pdbh_year = pdbh[pdbh_idx, 'year']
    
    if (sum(dbh_years %in% c(incr_years, pdbh_year)) == 0){
      next
    }
    
    if (length(pdbh_idx) > 0){
      dbh_year = pdbh_year
      dbh_tree = pdbh[pdbh_idx, 'dbh']
    } else {
      dbh_idx  = dbh_idx[dbh_years %in% incr_years]
      dbh_dat_tree = dbh[dbh_idx,]
      dbh_year = dbh_dat_tree[which.max(dbh_dat_tree$year),'year']
      dbh_tree = dbh_dat_tree[which.max(dbh_dat_tree$year),'dbh']
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
    
    dbh_m[tree_idx,year_idx] = dbh_calc
  }
  
  if (last_ti[i] < length(years)) {
    # dbh_m[tree, (last_ti[i]+1):N_years] = rep(0, N_years - last_ti[i])
    dbh_m[tree_idx, (last_ti[i]+1):N_years] = rep(NA, N_years - last_ti[i])
  }
  
  # dbh_in = which((!is.na(dbh_m[tree_idx,])) & (dbh_m[tree_idx,] >= 0))  
  # if (length(dbh_in)==0) {
  #   next()
  # } else {
  #   pred = allom.predict(allom.fit[toupper(pft)],
  #                        dbh = dbh_m[tree_idx, dbh_in],
  #                        pft = toupper(pft),
  #                        component = 6,
  #                        use = "Bg",
  #                        interval = "prediction", 
  #                        single.tree=TRUE)
  #   # pred = allom.predict(allom.fit[pft],dbh = dbh_m[tree, dbh_in],pft = pft,component = 6,use = "Bg",interval = "confidence")
  #   # pred_mean = colMeans(pred)
  #   pred_draw=pred[sample(seq(1,nrow(pred)), 1),]
  #   # pred_mean = apply(pred, 2, function(x) quantile(x, probs=0.5))
  #   # PI = apply(pred,2,quantile,c(0.025,0.5,0.975),na.rm=TRUE)
  #   ab_m[tree_idx,dbh_in] = pred_draw
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
N_samples = 200
dbh_c = array(NA, c(N_trees, N_census_years))
ab_c = array(NA, c(N_trees, N_census_years, N_samples))

for (i in 1:N_trees) {
  print(i)
  
  tree = trees[i]
  tree_idx = which(trees == tree)
  
  pft  = as.vector(taxaMatch[which(taxaMatch$taxon == taxon[tree_idx]),'species'] )
  
  dbh_idx   = which(dbh$stat_id == tree)
  dbh_years = dbh[dbh_idx, 'year']
  dbh_dat_tree = dbh[dbh_idx,]
  
  ab_idx = match(dbh_years, census_years)
  
  dbh_c[tree_idx, ab_idx] = dbh_dat_tree$dbh

}

dbh_c_cast = dcast(dbh, stat_id+plot_paleon~year, value.var=c('dbh'))

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

# melt measured
ab_m_melt = melt(ab_mr)
colnames(ab_m_melt) = c('tree_id', 'year_id', 'iter', 'ab')
ab_m_melt$taxon   = taxon[ab_m_melt$tree_id]
ab_m_melt$taxon   = taxaMatch[match(ab_m_melt$taxon, taxaMatch$taxon), 'species']
ab_m_melt$site_id = tree_site_id[ab_m_melt$tree_id]
ab_m_melt$year    = years[ab_m_melt$year_id]

saveRDS(ab_m_melt, file=paste0('data/ab_m_melt_', location, '.RDS'))

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
ab_p_melt = rbind(data.frame(ab_p1_melt, model=rep('Model RW + Census')), data.frame(ab_p2_melt, model=rep('Model RW')))

ab_p_melt$taxon   = taxon[ab_p_melt$tree_id]
ab_p_melt$taxon   = taxaMatch[match(ab_p_melt$taxon, taxaMatch$taxon), 'species']
# ab_p_melt$site_id = tree_site_id[ab_p_melt$tree_id]
ab_p_melt$year    = years[ab_p_melt$year_id]
saveRDS(ab_p_melt, paste0('allom/ab_p_', location, '_', mvers, '.RDS'))

# ab_p_melt = ab_p_melt[ab_p_melt$year<2013,]
# ab_m_melt = ab_m_melt[ab_m_melt$year<2013,]

ab_p_sum_by_iter = aggregate(ab ~ year+site_id+iter+model, ab_p_melt, function(x) sum(x, na.rm=TRUE))
# agb_sum2 = aggregate(ab~year_id+iter, data=posts_ab, function(x) sum(x, na.rm=TRUE))
ab_p_quants = aggregate(ab~year+site_id+model, data=ab_p_sum_by_iter, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
ab_p_quants = data.frame(ab_p_quants)
ab_p_quants = cbind(ab_p_quants[,1:3], ab_p_quants[,4])
colnames(ab_p_quants)[4:6] = c('ab25', 'ab50', 'ab975')  

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

cols = c('#084594', '#8c2d04', '#238b45', 'black')
cols_fill = c('#4292c6', '#fdd0a2', 'white', 'white')

cols = c('Model RW + Census'='#084594', 'Model RW'='#8c2d04', 'Empirical RW'='#238b45', 'Empirical Census'='black')
cols_fill = c('Model RW + Census'="#4292c6", 'Model RW'="#fdd0a2", 'Empirical RW'="lightgreen", 'Empirical Census'="white")

# cols = c('Model RW + Census'='grey18', 
#          'Model RW'='grey63', 
#          'Empirical RW'='black', 
#          'Empirical Census'='black')
# cols_fill = c('Model RW + Census'="grey24", 
#               'Model RW'="grey59", 
#               'Empirical RW'="white", 
#               'Empirical Census'="white")

ggplot() +  geom_ribbon(data=ab_p_quants, aes(x=year, ymin=ab25, ymax=ab975, fill=model), alpha=0.4) +
  geom_line(data=ab_p_quants, aes(x=year, y=ab50, colour=model), size=1) + 
   geom_ribbon(data=ab_m_quants, aes(x=year, ymin=ab25, ymax=ab975, colour='Empirical RW', fill='Empirical RW'), alpha=0.4)+
  geom_line(data=ab_m_quants, aes(x=year, y=ab50, colour='Empirical RW', fill='Empirical RW'), size=1) + 
  geom_point(data=ab_c_quants, aes(x=year, y=ab50, colour='Empirical Census', fill='Empirical Census'), size=2) + 
  geom_errorbar(data=ab_c_quants, aes(x=year, ymin=ab25, ymax=ab975), width=0.5)+
  geom_line(data=ab_p_quants, aes(x=year, y=ab25, colour=model), linetype=1, size=0.5) +
  geom_line(data=ab_p_quants, aes(x=year, y=ab975, colour=model), linetype=1, size=0.5) +
  facet_grid(site_id~., scales="free_y") + scale_color_manual(values=cols, name='Method')+
  scale_fill_manual(values=cols_fill, name='Method')+
  theme_bw() + theme(axis.title=element_text(size=14), axis.text=element_text(size=14)) +
  ylab("Biomass (Mg/ha)") + xlab('Year') +
  scale_x_continuous(breaks=seq(min(years), max(years), by=5))
ggsave(file=paste0('figures/AGB_by_site_', location, '.pdf'))
ggsave(file=paste0('figures/AGB_by_site_', location, '.png'))



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
## make some plots
#########################################################################################################################################

abi_p1 = apply(ab_pr1, c(1,3), function(x) diff(x))
abi_p1 = aperm(abi_p1, c(2, 1, 3))
abi_p1_melt = melt(abi_p1)
colnames(abi_p1_melt) = c('tree_id', 'year_id', 'iter', 'abi')

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
abi_m = t(apply(ab_mr_median, c(1), function(x) diff(x)))

# census
ab_cr_median = apply(ab_cr, c(1,2), median, na.rm=TRUE)
census_years = as.numeric(census_years)
abi_c = t(apply(ab_cr_median, 1, function(x) diff(x)))
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
  geom_point(data=abi_c_sum, aes(x=year, y=abi, colour='Empirical Census', fill='Empirical Census'),size=2) + 
  facet_grid(site_id~.) + scale_color_manual(values=cols, name='Method')+#, labels=c('RW + Census', 'RW')) + 
  scale_fill_manual(values=cols_fill, name='Method')+#, labels=c('RW + Census', 'RW')) + 
  theme_bw()+
  ylab("Biomass Increment (Mg/ha/year)") + xlab('Year') + 
  scale_x_continuous(breaks=seq(min(years), 2012, by=5), limits=c(min(years)+4,2012)) + 
  ylim(c(0.3, 8.5))
ggsave(file=paste0('figures/', figures_dir, '/AGBI_by_site_', location, '.pdf'))
ggsave(file=paste0('figures/', figures_dir, '/AGBI_by_site_', location, '.png'))
