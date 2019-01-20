library(PEcAn.allometry)
library(ggplot2)
library(reshape2)
library(abind)

run_predict = TRUE
sample_allom = FALSE
site = "IT"

figures_dir = 'NOCOVAR'
if (!file.exists(figures_dir)){
  dir.create(paste0("figures/", figures_dir))
}

# source('config_HMC')

dvers = "v0.1"
mvers = "v0.1"

# fname_data = paste0('tree_data_IT_', dvers)
# load(file=paste0('data/dump/', fname_data, '.rdata'))
# 
fname_data = paste0('tree_data_IT_STAN_', dvers)
dat = readRDS(paste0('data/dump/', fname_data, '.RDS'))

N_trees = dat$N_trees
N_years = dat$N_years
N_vals = dat$N_vals
logXobs = dat$logXobs
logPDobs = dat$logPDobs
year_idx = dat$year_idx
taxon = dat$taxon
N_taxa = dat$N_taxa
pdbh_year = dat$pdbh_year
idx_tree = dat$idx_tree
pdbh2val = dat$pdbh2val
x2tree = dat$x2tree
x2year = dat$x2year
taxaMatch = dat$taxaMatch
years= dat$years

trees = seq(1, N_trees)

# fnames = c('ring_model_t_pdbh_IT')
fnames = "ring_model_t_pdbh_IT_stan_v2"
models = c('Model RW')

post= list()
for (i in 1:length(fnames)) {
  fname_model = fnames[i]
  load(file   = paste0('output/', fname_model, '_', site, '_', mvers, '.Rdata'))
  post[[i]]   = post
}  

data(allom.components)
allom.components

pfts = list(ACRA = data.frame(spcd=316,acronym="ACRA"),
            ACSA = data.frame(spcd=318,acronym="ACSA3"),
            BEPA = data.frame(spcd=375,acronym="BEPA"),
            PIRE = data.frame(spcd=125, acronym="PIRE"),
            PIST = data.frame(spcd=129, acronym="PIST"))
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

taxaMatch$group = c(1, 1, 2, 3, 3)
group_name      = c('Acer', 'Acer', 'Betula', 'Pinus', 'Pinus')
taxaMatch

niter = length(post$lp__)
burn = 0
keep = niter - burn

nmodels = length(fnames)

taxon = as.vector(taxon)

# dbh_1 = build_dbh_p(post[[1]], x2year, x2tree, N_years, trees, keep, taxaMatch, taxon)
build_dbh_p <- function(out, x2year, x2tree, N_years, trees, keep, taxaMatch, taxon){
  
  # col_names = sapply(strsplit(colnames(out), '\\['), function(x) x[[1]])
  col_names = names(out)
  niter     = length(out[[1]])
  
  pft   = as.vector(taxaMatch[match(taxon, taxaMatch$species),'species'])
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
  
  saveRDS(dbh_p_1, file=paste0('allom/dbh_p_',site, '_', mvers, '.rds'))
  saveRDS(ab_p_1, paste0('allom/ab_p_', site, '_', mvers, '.rds'))
  
} else {
  dbh_p_1 = readRDS(paste0('allom/dbh_p_', site, '_', mvers, '.rds'))
  ab_p_1  = readRDS(paste0('allom/ab_p_', site, '_', mvers, '.rds'))
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

saveRDS(dbh_melt, 'data/IT_DBH_iterations.RDS') 

dbh_mean = aggregate(dbh~tree+year, dbh_melt, mean, na.rm=TRUE) 
dbh_sum = aggregate(dbh~year, dbh_mean, sum, na.rm=TRUE)

ggplot(data=dbh_sum) + 
  geom_line(data=dbh_sum, aes(x=year, y=dbh)) #+ xlim(c(1960,2012))
ggsave('figures/NOCOVAR/sum_dbh_by_plot_IT_v1.0.pdf')

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
  pft  = as.vector(taxaMatch[which(taxaMatch$species == taxon[tree_idx]),'species'] )
  
  # for now use average increments
  incr_tree  = exp(logXobs[which(m2tree == tree)])
  incr_years = years[m2ti[which(m2tree == tree)]]
  
  pdbh_idx = which(pdbh$id == tree)
  pdbh_year = pdbh[pdbh_idx, 'year']
  
  if (length(pdbh_idx) > 0){
    dbh_year = pdbh_year
    dbh_tree = pdbh[pdbh_idx, 'dbh']
  } 
  
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

for (p in 1:length(pft_list)){
  idx_m = c(idx_m, which(taxon==pft_list[p]))
  dbh_m_list[[p]] = dbh_m[which(taxon==pft_list[p]),]
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
dbh_melt$year = years[dbh_melt$year]

dbh_m_mean = aggregate(dbh~tree+year, dbh_m_melt, mean) 
dbh_m_sum = aggregate(dbh~year, dbh_m_mean, sum)
dbh_m_sum$year = years[dbh_m_sum$year]

ggplot() + 
  geom_line(data=dbh_sum, aes(x=year, y=dbh)) + #+ xlim(c(1960,2012))
  geom_line(data=dbh_m_sum, aes(x=year, y=dbh))
ggsave('figures/NOCOVAR/sum_dbh_by_plot_IT_v1.0.pdf')

#########################################################################################################################################
## make the AGB plot
#########################################################################################################################################

# in Kg/plot, rescale so it is Mg/ha
ab_mr = ab_m/(15^2*pi) * 10
ab_pr1 = ab_p_1/(15^2*pi) * 10

# melt measured
ab_m_melt = melt(ab_mr)
colnames(ab_m_melt) = c('tree_id', 'year_id', 'iter', 'ab')
ab_m_melt$taxon   = taxon[ab_m_melt$tree_id]
ab_m_melt$taxon   = taxaMatch[match(ab_m_melt$taxon, taxaMatch$taxon), 'species']
ab_m_melt$year    = years[ab_m_melt$year_id]

saveRDS(ab_m_melt, file=paste0('data/ab_m_melt_', site, '.RDS'))

# melt predicted
ab_p1_melt = melt(ab_pr1)
colnames(ab_p1_melt) = c('tree_id', 'year_id', 'iter', 'ab')

# put them together
ab_p_melt = ab_p1_melt

ab_p_melt$taxon   = taxon[ab_p_melt$tree_id]
ab_p_melt$taxon   = taxaMatch[match(ab_p_melt$taxon, taxaMatch$species), 'species']
# ab_p_melt$site_id = tree_site_id[ab_p_melt$tree_id]
ab_p_melt$year    = years[ab_p_melt$year_id]
saveRDS(ab_p_melt, paste0('data/ab_p_', site, '_', mvers, '.RDS'))

ab_p_sum_by_iter = aggregate(ab ~ year+iter, ab_p_melt, function(x) sum(x, na.rm=TRUE))
# agb_sum2 = aggregate(ab~year_id+iter, data=posts_ab, function(x) sum(x, na.rm=TRUE))
ab_p_quants = aggregate(ab~year, data=ab_p_sum_by_iter, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
ab_p_quants = data.frame(ab_p_quants)
ab_p_quants = cbind(ab_p_quants[,1], ab_p_quants[,2])
colnames(ab_p_quants) = c('year', 'ab25', 'ab50', 'ab975')  

ab_m_sum_by_iter = aggregate(ab ~ year+iter, ab_m_melt, function(x) sum(x, na.rm=TRUE))
ab_m_quants = aggregate(ab~year, data=ab_m_sum_by_iter, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
ab_m_quants = data.frame(ab_m_quants)
ab_m_quants = cbind(ab_m_quants[,1], ab_m_quants[,2])
colnames(ab_m_quants) = c('year', 'ab25', 'ab50', 'ab975')  

ab_m_sum = aggregate(ab ~ year, ab_m_melt, function(x) sum(x, na.rm=TRUE))

# cols = c('#084594', '#8c2d04', '#238b45', 'black')
# cols_fill = c('#4292c6', '#fdd0a2', 'white', 'white')
# 
# cols = c('Model RW + Census'='#084594', 'Model RW'='#8c2d04', 'Empirical RW'='#238b45', 'Empirical Census'='black')
# cols_fill = c('Model RW + Census'="#4292c6", 'Model RW'="#fdd0a2", 'Empirical RW'="lightgreen", 'Empirical Census'="white")
# 
# ggplot() +  geom_ribbon(data=ab_p_quants, aes(x=year, ymin=ab25, ymax=ab975, fill=model), alpha=0.4) +
#   geom_line(data=ab_p_quants, aes(x=year, y=ab50, colour=model), size=1) + 
#   geom_ribbon(data=ab_m_quants, aes(x=year, ymin=ab25, ymax=ab975, colour='Empirical RW', fill='Empirical RW'), alpha=0.4)+
#   geom_line(data=ab_m_quants, aes(x=year, y=ab50, colour='Empirical RW', fill='Empirical RW'), size=1) + 
#   geom_point(data=ab_c_quants, aes(x=year, y=ab50, colour='Empirical Census', fill='Empirical Census'), size=2) + 
#   geom_errorbar(data=ab_c_quants, aes(x=year, ymin=ab25, ymax=ab975), width=0.5)+
#   geom_line(data=ab_p_quants, aes(x=year, y=ab25, colour=model), linetype=1, size=0.5) +
#   geom_line(data=ab_p_quants, aes(x=year, y=ab975, colour=model), linetype=1, size=0.5) +
#   facet_grid(site_id~., scales="free_y") + scale_color_manual(values=cols, name='Method')+
#   scale_fill_manual(values=cols_fill, name='Method')+
#   theme_bw() + theme(axis.title=element_text(size=14), axis.text=element_text(size=14)) +
#   ylab("Biomass (Mg/ha)") + xlab('Year') +
#   scale_x_continuous(breaks=seq(min(years), max(years), by=5))
# ggsave(file=paste0('figures/AGB_by_site_', site, '.pdf'))
# ggsave(file=paste0('figures/AGB_by_site_', site, '.png'))
# 

ab_p_quants = as.data.frame(ab_p_quants)
ab_m_quants = as.data.frame(ab_m_quants)

ggplot() +  geom_ribbon(data=ab_p_quants, aes(x=year, ymin=ab25, ymax=ab975), alpha=0.4) +
  geom_line(data=ab_p_quants, aes(x=year, y=ab50), size=1) + 
  geom_line(data=ab_m_quants, aes(x=year, y=ab50, colour='Empirical RW', fill='Empirical RW'), size=1) + 
  geom_line(data=ab_p_quants, aes(x=year, y=ab25), linetype=1, size=0.5) +
  geom_line(data=ab_p_quants, aes(x=year, y=ab975), linetype=1, size=0.5) +
  # scale_fill_manual(values=cols_fill, name='Method')+
  theme_bw() + theme(axis.title=element_text(size=14), axis.text=element_text(size=14)) +
  ylab("Biomass (Mg/ha)") + xlab('Year') #+
  # scale_x_continuous(breaks=seq(min(years), max(years), by=5))
ggsave(file=paste0('figures/', figures_dir, '/AGB_by_site_', site, '.pdf'))
ggsave(file=paste0('figures/', figures_dir, '/AGB_by_site_', site, '.png'))

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

saveRDS(abi_p_melt, file=paste0('data/abi_p_', site, '_', mvers, '.RDS'))

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
ggsave(file=paste0('figures/', figures_dir, '/AGBI_by_site_', site, '.pdf'))
ggsave(file=paste0('figures/', figures_dir, '/AGBI_by_site_', site, '.png'))
