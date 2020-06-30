library(nimble)

source('config_HF')

dvers="v5.0"
mvers="v5.0"
site = "HARVARD"

# fname_data = paste0('tree_data_HMC_', dvers)

# load(paste0('sites/', site, '/data/tree_data_20_', site ,'_NIMBLE_', dvers, '.rdata'))
load(paste0('sites/', site, '/data/tree_data_', site, '_NIMBLE_', dvers, '.rdata'))

# check initial values
D_init = rep(0, N_X)
for (i in 1:N_trees){
  D_init[first_ti[i]:cs_last_ti[i]] = 0.0 + 2 * 0.5/ 10 * seq(1, last_ti[i])
}

D_init[sap2x] < max_size
all(D_init[sap2x] < max_size)
#######################################################################################################################################
# MODEL 1: date, sapling, size, pdbh, zero covariance 
#######################################################################################################################################
fname_model = 'ring_model_t_date_sapl_size_pdbh_NOCOVAR'
# load(paste0('data/dump/', fname_data, '.rdata'))
load(paste0('sites/', site, '/data/tree_data_', site ,'_NIMBLE_', dvers, '.rdata'))
source(paste0('sites/', site, '/models/', fname_model, '.R'))

m <- nimbleModel(body(fname_model),
                 constants = list(N_inc = N_inc,
                                  N_dbh = N_dbh,
                                  N_pdbh = N_pdbh,
                                  N_trees = N_trees,
                                  N_years = N_years,
                                  N_X = N_X,
                                  N_taxa = N_taxa,
                                  dbh_tree_id = dbh_tree_id,
                                  dbh_year_id = dbh_year_id,
                                  dbh_day_id = dbh_day_id,
                                  pdbh_day_id = pdbh_day_id,
                                  open_dbh = open_dbh,
                                  m2tree = m2tree,
                                  m2ti = m2ti,
                                  m2nc = m2nc,
                                  n1cores = n1cores,
                                  n2cores=n2cores,
                                  n3cores=n3cores,
                                  n4cores=n4cores,
                                  i1core2m = i1core2m,
                                  i2core2m = i2core2m,
                                  i3core2m = i3core2m,
                                  i4core2m = i4core2m,
                                  x2tree = x2tree, 
                                  x2year = x2year,
                                  ones = ones,
                                  meas2d = meas2d,
                                  meas2x = meas2x, 
                                  pdbh2d = pdbh2d,
                                  last_ti = last_ti,
                                  first_ti = first_ti,
                                  cs_last_ti = cs_last_ti,
                                  taxon = taxon,
                                  sap2x = sap2x,
                                  max_size = max_size,
                                  N_saplings = N_saplings,
                                  tree_site_id = tree_site_id
                                  
                 ), check=FALSE)


#debug(m$setData)

m$setData(list(logDobs = logDobs, logXobs = logXobs, logPDobs = logPDobs, sapling_constraint = rep(1, N_saplings)))

inits = list(X = rep(0.5, N_X), 
             D0 = rep(0, N_trees),
             # D = rep(.1, N_X), 
             beta0 = .1, 
             beta   = rep(.1, N_trees),
             beta_t = rep(.1, N_years), 
             sig_x_obs = .5,
             sig_d_obs = .1,
             beta_sd   = .1, 
             beta_t_sd = .1, 
             sig_x     = .3)
inits = c(inits, b0 = 0.5, b1 = 10)
inits = c(inits, list(beta_slope=0.1))

m$setInits(inits)

spec <- configureMCMC(m, thin = 5) 
spec$addMonitors(c('D', 'X', 'beta', 'beta_t', 'beta_t_sd', 'beta_sd', 'sig_x_obs', 'sig_d_obs'))

Rmcmc <- buildMCMC(spec)
cm    <- compileNimble(m)
Cmcmc <- compileNimble(Rmcmc, project = m)
Cmcmc$run(10000) #25000
out <- as.matrix(Cmcmc$mvSamples)#[101:500, ]


saveRDS(out, file = paste0('sites/', site, '/output/', fname_model, '_', mvers, '.RDS'))
# #######################################################################################################################################
# # MODEL 2: size, pdbh, zero covariance 
# #######################################################################################################################################
# fname_model = 'ring_model_t_size_pdbh_NOCOVAR'
# load(paste0('sites/', site, '/data/tree_data_', site ,'_NIMBLE_', dvers, '.rdata'))
# # load(paste0('data/dump/', fname_data, '.rdata'))
# source(paste0('sites/', site, '/models/', fname_model, '.R'))
# 
# m <- nimbleModel(body(fname_model),
#                  constants = list(N_inc = N_inc,
#                                   N_dbh = N_dbh,
#                                   N_pdbh = N_pdbh,
#                                   N_trees = N_trees,
#                                   N_years = N_years,
#                                   N_X = N_X,
#                                   N_taxa = N_taxa,
#                                   dbh_tree_id = dbh_tree_id,
#                                   dbh_year_id = dbh_year_id,
#                                   dbh_day_id = dbh_day_id,
#                                   pdbh_day_id = pdbh_day_id,
#                                   open_dbh = open_dbh,
#                                   m2tree = m2tree,
#                                   m2ti = m2ti,
#                                   m2nc = m2nc,
#                                   n1cores = n1cores,
#                                   n2cores=n2cores,
#                                   n3cores=n3cores,
#                                   n4cores=n4cores,
#                                   i1core2m = i1core2m,
#                                   i2core2m = i2core2m,
#                                   i3core2m = i3core2m,
#                                   i4core2m = i4core2m,
#                                   x2tree = x2tree, 
#                                   x2year = x2year,
#                                   ones = ones,
#                                   meas2d = meas2d,
#                                   meas2x = meas2x, 
#                                   pdbh2d = pdbh2d,
#                                   last_ti = last_ti,
#                                   first_ti = first_ti,
#                                   cs_last_ti = cs_last_ti,
#                                   taxon = taxon,
#                                   tree_site_id = tree_site_id
#                                   
#                  ), check=FALSE)
# 
# 
# #debug(m$setData)
# 
# m$setData(list(logDobs = logDobs, logXobs = logXobs, logPDobs = logPDobs))
# 
# inits = list(X = rep(0.5, N_X), 
#              D0 = rep(0, N_trees),
#              # D = rep(.1, N_X), 
#              beta0 = .1, 
#              beta   = rep(.1, N_trees),
#              beta_t = rep(.1, N_years), 
#              sig_x_obs = .5,
#              sig_d_obs = .1,
#              beta_sd   = .1, 
#              beta_t_sd = .1, 
#              sig_x     = .3)
# inits = c(inits, list(beta_slope=0.1))
# 
# m$setInits(inits)
# 
# spec <- configureMCMC(m, thin = 5) 
# spec$addMonitors(c('D', 'X', 'beta', 'beta_t', 'beta_t_sd', 'beta_sd', 'sig_x_obs', 'sig_d_obs'))
# 
# Rmcmc <- buildMCMC(spec)
# cm    <- compileNimble(m)
# Cmcmc <- compileNimble(Rmcmc, project = m)
# Cmcmc$run(10000) #25000
# out <- as.matrix(Cmcmc$mvSamples)#[101:500, ]
# 
# saveRDS(out, file =paste0('sites/', site, '/output/', fname_model, '_', mvers, '.RDS'))
# 
# #######################################################################################################################################
# # MODEL 3: pdbh, zero covariance 
# #######################################################################################################################################
# fname_model = 'ring_model_t_pdbh_NOCOVAR'
# load(paste0('sites/', site, '/data/tree_data_', site ,'_NIMBLE_', dvers, '.rdata'))
# # load(paste0('data/dump/', fname_data, '.rdata'))
# source(paste0('sites/', site, '/models/', fname_model, '.R'))
# 
# m <- nimbleModel(body(fname_model),
#                  constants = list(N_inc = N_inc,
#                                   N_dbh = N_dbh,
#                                   N_pdbh = N_pdbh,
#                                   N_trees = N_trees,
#                                   N_years = N_years,
#                                   N_X = N_X,
#                                   N_taxa = N_taxa,
#                                   dbh_tree_id = dbh_tree_id,
#                                   dbh_year_id = dbh_year_id,
#                                   dbh_day_id = dbh_day_id,
#                                   pdbh_day_id = pdbh_day_id,
#                                   open_dbh = open_dbh,
#                                   m2tree = m2tree,
#                                   m2ti = m2ti,
#                                   m2nc = m2nc,
#                                   n1cores = n1cores,
#                                   n2cores=n2cores,
#                                   n3cores=n3cores,
#                                   n4cores=n4cores,
#                                   i1core2m = i1core2m,
#                                   i2core2m = i2core2m,
#                                   i3core2m = i3core2m,
#                                   i4core2m = i4core2m,
#                                   x2tree = x2tree, 
#                                   x2year = x2year,
#                                   ones = ones,
#                                   meas2d = meas2d,
#                                   meas2x = meas2x, 
#                                   pdbh2d = pdbh2d,
#                                   last_ti = last_ti,
#                                   first_ti = first_ti,
#                                   cs_last_ti = cs_last_ti,
#                                   taxon = taxon,
#                                   tree_site_id = tree_site_id
#                                   
#                  ), check=FALSE)
# 
# 
# #debug(m$setData)
# 
# m$setData(list(logDobs = logDobs, logXobs = logXobs, logPDobs = logPDobs))
# 
# inits = list(X = rep(0.5, N_X), 
#              D0 = rep(0, N_trees),
#              # D = rep(.1, N_X), 
#              beta0 = .1, 
#              beta   = rep(.1, N_trees),
#              beta_t = rep(.1, N_years), 
#              sig_x_obs = .5,
#              sig_d_obs = .1,
#              beta_sd   = .1, 
#              beta_t_sd = .1, 
#              sig_x     = .3)
# inits = c(inits, b0 = 0.5, b1 = 10)
# inits = c(inits, list(beta_slope=0.1))
# 
# m$setInits(inits)
# 
# spec <- configureMCMC(m, thin = 5) 
# spec$addMonitors(c('D', 'X', 'beta', 'beta_t', 'beta_t_sd', 'beta_sd', 'sig_x_obs', 'sig_d_obs'))
# 
# Rmcmc <- buildMCMC(spec)
# cm    <- compileNimble(m)
# Cmcmc <- compileNimble(Rmcmc, project = m)
# Cmcmc$run(10000) #25000
# out <- as.matrix(Cmcmc$mvSamples)#[101:500, ]
# 
# saveRDS(out, file =  paste0('sites/', site, '/output/', fname_model, '_', mvers, '.RDS'))
#########################################################################################################################################
## MODEL 4: RW ONLY, pdbh, zero covariance; input sig_d_obs from CENSUS + RW MODEL 
#########################################################################################################################################
fname_model = fname_model = 'ring_model_t_date_sapl_size_pdbh_NOCOVAR'
# 'ring_model_t_date_sapl_size_pdbh_NOCOVAR'
out = readRDS(paste0('sites/', site, '/output/', fname_model, '_', mvers, '.RDS'))

col_names = sapply(strsplit(colnames(out), '\\['), function(x) x[[1]])
hist(out[,which(col_names=="sig_d_obs")])
sig_d_obs = mean(out[,which(col_names=="sig_d_obs")])

fname_model = 'ring_model_t_pdbh_nc_NOCOVAR_sigd'

load(paste0('sites/', site, '/data/tree_data_no_census_', site ,'_NIMBLE_', dvers, '.rdata'))
source(paste0('sites/', site, '/models/', fname_model, '.R'))

# check initial values
D_init = rep(0, N_X_p)
for (i in 1:N_trees_p){
  # for (j in 1:last_ti[i]){
  D_init[first_ti_p[i]:cs_last_ti_p[i]] = 0.0 + 2 * 0.5/ 10 * seq(1, last_ti_p[i])
  # }
}
all(D_init>0)

x2idx_p = match(x2tree_p, unique(x2tree_p))

# full model
m <- nimbleModel(body(fname_model),
                 constants = list(N_inc = N_inc,
                                  N_pdbh = N_pdbh,
                                  N_trees = N_trees_p,
                                  N_years = N_years,
                                  N_X = N_X_p,
                                  N_taxa = N_taxa,
                                  pdbh_day_id = pdbh_day_id,
                                  open_dbh = open_dbh,
                                  m2tree = m2tree,
                                  m2ti = m2ti,
                                  m2nc = m2nc,
                                  n1cores = n1cores,
                                  n2cores=n2cores,
                                  n3cores=n3cores,
                                  n4cores=n4cores,
                                  i1core2m = i1core2m,
                                  i2core2m = i2core2m,
                                  i3core2m = i3core2m,
                                  i4core2m = i4core2m,
                                  x2tree = x2idx_p, 
                                  x2year = x2year_p,
                                  ones = ones,
                                  meas2x = meas2x_p, 
                                  pdbh2d = pdbh2d_p,
                                  last_ti = last_ti_p,
                                  first_ti = first_ti_p,
                                  cs_last_ti = cs_last_ti_p,
                                  sig_d_obs = sig_d_obs,
                                  tree_site_id = tree_site_id
                 ), check=FALSE)


#debug(m$setData)

m$setData(list(logXobs = logXobs, logPDobs = logPDobs))#, tau2_constraint=c(1), tau3_constraint=c(1)))#not_cens = rep(1, N_saplings)))#s

inits = list(X = rep(0.5, N_X_p), 
             D0 = rep(0, N_trees_p),
             # D = rep(.1, N_X_p), 
             beta0 = .1, 
             beta   = rep(.1, N_trees_p),
             beta_t = rep(.1, N_years), 
             sig_x_obs = .5,
             #sig_d_obs = .1,
             beta_sd   = .1, 
             beta_t_sd = .1, 
             sig_x     = .3)

m$setInits(inits)

spec <- configureMCMC(m, thin = 5) 
spec$addMonitors(c('D', 'X', 'beta', 'beta_t', 'beta_t_sd', 'beta_sd', 'sig_x_obs'))

Rmcmc <- buildMCMC(spec)
cm    <- compileNimble(m)
Cmcmc <- compileNimble(Rmcmc, project = m)
Cmcmc$run(10000) #25000
out <- as.matrix(Cmcmc$mvSamples)#[101:500, ]

saveRDS(out, file = paste0('sites/', site, '/output/', fname_model, '_', mvers, '.RDS'))

# #########################################################################################################################################
# ## MODEL 5: RW ONLY, pdbh, zero covariance; try to estimate sig_d_obs 
# #########################################################################################################################################
# fname_model = 'ring_model_t_pdbh_nc_NOCOVAR'
# 
# load(paste0('sites/', site, '/data/tree_data_', site ,'_no_census_NIMBLE_', dvers, '.rdata'))
# source(paste0('sites/', site, '/models/', fname_model, '.R'))
# 
# # check initial values
# D_init = rep(0, N_X_p)
# for (i in 1:N_trees_p){
#   # for (j in 1:last_ti[i]){
#   D_init[first_ti_p[i]:cs_last_ti_p[i]] = 0.0 + 2 * 0.5/ 10 * seq(1, last_ti_p[i])
#   # }
# }
# all(D_init>0)
# 
# x2idx_p = match(x2tree_p, unique(x2tree_p))
# 
# # full model
# m <- nimbleModel(body(fname_model),
#                  constants = list(N_inc = N_inc,
#                                   N_pdbh = N_pdbh,
#                                   N_trees = N_trees_p,
#                                   N_years = N_years,
#                                   N_X = N_X_p,
#                                   N_taxa = N_taxa,
#                                   pdbh_day_id = pdbh_day_id,
#                                   open_dbh = open_dbh,
#                                   m2tree = m2tree,
#                                   m2ti = m2ti,
#                                   m2nc = m2nc,
#                                   n1cores = n1cores,
#                                   n2cores=n2cores,
#                                   n3cores=n3cores,
#                                   n4cores=n4cores,
#                                   i1core2m = i1core2m,
#                                   i2core2m = i2core2m,
#                                   i3core2m = i3core2m,
#                                   i4core2m = i4core2m,
#                                   x2tree = x2idx_p, 
#                                   x2year = x2year_p,
#                                   ones = ones,
#                                   meas2x = meas2x_p, 
#                                   pdbh2d = pdbh2d_p,
#                                   last_ti = last_ti_p,
#                                   first_ti = first_ti_p,
#                                   cs_last_ti = cs_last_ti_p#,
#                                   #sig_d_obs = sig_d_obs#,
#                                   # tree_site_id = tree_site_id
#                  ), check=FALSE)
# 
# 
# #debug(m$setData)
# 
# m$setData(list(logXobs = logXobs, logPDobs = logPDobs))#, tau2_constraint=c(1), tau3_constraint=c(1)))#not_cens = rep(1, N_saplings)))#s
# 
# inits = list(X = rep(0.5, N_X_p), 
#              D0 = rep(0, N_trees_p),
#              # D = rep(.1, N_X_p), 
#              beta0 = .1, 
#              beta   = rep(.1, N_trees_p),
#              beta_t = rep(.1, N_years), 
#              sig_x_obs = .5,
#              sig_d_obs = .1,
#              beta_sd   = .1, 
#              beta_t_sd = .1, 
#              sig_x     = .3)
# 
# m$setInits(inits)
# 
# spec <- configureMCMC(m, thin = 5) 
# spec$addMonitors(c('D', 'X', 'beta', 'beta_t', 'beta_t_sd', 'beta_sd', 'sig_x_obs', 'sig_d_obs'))
# 
# Rmcmc <- buildMCMC(spec)
# cm    <- compileNimble(m)
# Cmcmc <- compileNimble(Rmcmc, project = m)
# Cmcmc$run(10000) #25000
# out <- as.matrix(Cmcmc$mvSamples)#[101:500, ]
# 
# saveRDS(out, file = paste0('sites/', site, '/output/', fname_model, '_', mvers, '.RDS'))
# 
# #######################################################################################################################################
# # MODEL 6: average ring width
# #######################################################################################################################################
# fname_data = paste0('tree_data_HMC_', dvers)
# load(paste0('sites/', site, '/data/tree_data_', site ,'_NIMBLE_', dvers, '.rdata'))
# # load(paste0('data/dump/', fname_data, '.rdata'))
# 
# fname_model = 'ring_model_t_date_size_pdbh_avg'
# source(paste0('sites/', site, '/models/', fname_model, '.R'))
# 
# # m <- nimbleModel(body(fname_model),
# #                  constants = list(N_inc = N_inc,
# #                                   N_dbh = N_dbh,
# #                                   N_pdbh = N_pdbh,
# #                                   N_trees = N_trees,
# #                                   N_years = N_years,
# #                                   N_X = N_X,
# #                                   N_taxa = N_taxa,
# #                                   dbh_tree_id = dbh_tree_id,
# #                                   dbh_year_id = dbh_year_id,
# #                                   dbh_day_id = dbh_day_id,
# #                                   pdbh_day_id = pdbh_day_id,
# #                                   open_dbh = open_dbh,
# #                                   m2tree = m2tree,
# #                                   m2ti = m2ti,
# #                                   m2nc = m2nc,
# #                                   n1cores = n1cores,
# #                                   n2cores=n2cores,
# #                                   n3cores=n3cores,
# #                                   n4cores=n4cores,
# #                                   i1core2m = i1core2m,
# #                                   i2core2m = i2core2m,
# #                                   i3core2m = i3core2m,
# #                                   i4core2m = i4core2m,
# #                                   x2tree = x2tree, 
# #                                   x2year = x2year,
# #                                   ones = ones,
# #                                   meas2d = meas2d,
# #                                   meas2x = meas2x, 
# #                                   pdbh2d = pdbh2d,
# #                                   last_ti = last_ti,
# #                                   first_ti = first_ti,
# #                                   cs_last_ti = cs_last_ti,
# #                                   taxon = taxon,
# #                                   tree_site_id = tree_site_id
# #                                   
# #                  ), check=FALSE)
# 
# 
# m <- nimbleModel(body(fname_model),
#                  constants = list(N_inc = N_inc_a,
#                                   N_dbh = N_dbh,
#                                   N_pdbh = N_pdbh,
#                                   N_trees = N_trees,
#                                   N_years = N_years,
#                                   N_X = N_X,
#                                   N_taxa = N_taxa,
#                                   dbh_tree_id = dbh_tree_id,
#                                   dbh_year_id = dbh_year_id,
#                                   dbh_day_id = dbh_day_id,
#                                   pdbh_day_id = pdbh_day_id,
#                                   open_dbh = open_dbh,
#                                   m2tree = m2tree_a,
#                                   m2ti = m2ti_a,
#                                   x2tree = x2tree,
#                                   x2year = x2year,
#                                   meas2x = meas2x_a,
#                                   meas2d = meas2d,
#                                   pdbh2d = pdbh2d,
#                                   last_ti = last_ti,
#                                   first_ti = first_ti,
#                                   cs_last_ti = cs_last_ti,
#                                   taxon = taxon,
#                                   tree_site_id = tree_site_id
#                  ), check=FALSE)
# 
# 
# #debug(m$setData)
# 
# # m$setData(list(logDobs = logDobs, logXobs = logXobs_a))#, tau2_constraint=c(1), tau3_constraint=c(1)))#not_cens = rep(1, N_saplings)))#s
# m$setData(list(logDobs = logDobs, logXobs = logXobs_a, logPDobs = logPDobs))
# 
# inits = list(X = rep(0.5, N_X),
#              D0 = rep(0, N_trees),
#              D = rep(.1, N_X),
#              beta0 = .1,
#              beta   = rep(.1, N_trees),
#              beta_t = rep(.1, N_years),
#              sig_x_obs = .5,
#              sig_d_obs = .1,
#              beta_sd   = .1,
#              beta_t_sd = .1,
#              sig_x     = .3)
# inits = c(inits, b0 = 0.5, b1 = 10)
# inits = c(inits, list(beta_slope=0.1))
# 
# m$setInits(inits)
# 
# spec <- configureMCMC(m, thin = 5)
# spec$addMonitors(c('D', 'X', 'beta', 'beta_t', 'beta_t_sd', 'beta_sd', 'sig_x_obs', 'sig_d_obs'))
# 
# Rmcmc <- buildMCMC(spec)
# cm    <- compileNimble(m)
# Cmcmc <- compileNimble(Rmcmc, project = m)
# Cmcmc$run(10000) #25000
# out <- as.matrix(Cmcmc$mvSamples)#[101:500, ]
# 
# saveRDS(out, file =  paste0('sites/', site, '/output/', fname_model, '_', mvers, '.RDS'))

###


# #######################################################################################################################################
# # single core
# #######################################################################################################################################
# 
# # load('data/dump/tree_full.rdata')
# suff = '1c'
# load(paste0('data/dump/tree_data_20_', dvers, '_', suff, '.rdata'))
# 
# # check initial values
# D_init = rep(0, N_X)
# for (i in 1:N_trees){
#   # for (j in 1:last_ti[i]){
#   D_init[first_ti[i]:cs_last_ti[i]] = 0.0 + 2 * 0.5/ 10 * seq(1, last_ti[i])
#   # }
# }
# 
# D_init[sap2x] < max_size
# all(D_init[sap2x] < max_size)
# 
# fname_model = 'ring_model_t_date_sapl_size_avg'
# # fname_model = 'ring_model_t_date_sapl_size_constraint'
# # fname_model = 'ring_model_t_date_sapl_size_negtau'
# 
# source(paste0('r/models/', fname_model, '.R'))
# # full model
# m <- nimbleModel(body(fname_model),
#                  constants = list(N_inc = N_inc,
#                                   N_dbh = N_dbh,
#                                   N_trees = N_trees,
#                                   N_years = N_years,
#                                   N_X = N_X,
#                                   N_taxa = N_taxa,
#                                   dbh_tree_id = dbh_tree_id,
#                                   dbh_year_id = dbh_year_id,
#                                   dbh_day_id = dbh_day_id,
#                                   open_dbh = open_dbh,
#                                   m2tree = m2tree,
#                                   m2ti = m2ti,
#                                   x2tree = x2tree, 
#                                   x2year = x2year,
#                                   meas2d = meas2d,
#                                   meas2x = meas2x, 
#                                   last_ti = last_ti,
#                                   first_ti = first_ti,
#                                   cs_last_ti = cs_last_ti,
#                                   taxon = taxon,
#                                   sap2x = sap2x,
#                                   max_size = max_size,
#                                   N_saplings = N_saplings,
#                                   tree_site_id = tree_site_id
#                                   
#                  ), check=FALSE)
# 
# 
# #debug(m$setData)
# 
# # m$setData(list(logDobs = logDobs, logXobs = logXobs, sapling_constraint = rep(1, N_saplings), tau2_constraint=c(1), tau3_constraint=c(1)))#not_cens = rep(1, N_saplings)))#s
# 
# m$setData(list(logDobs = logDobs, logXobs = logXobs, sapling_constraint = rep(1, N_saplings)))#, tau2_constraint=c(1), tau3_constraint=c(1)))#not_cens = rep(1, N_saplings)))#s
# 
# # m$setData(list(logDobs = logDobs, logXobs = logXobs, sapling_constraint = rep(1, N_saplings), tau2_constraint=c(1), tau3_constraint=c(1)))#not_cens = rep(1, N_saplings)))#s
# 
# 
# 
# inits = list(X = rep(0.5, N_X), 
#              D0 = rep(0, N_trees),
#              # D = rep(.1, N_X), 
#              beta0 = .1, 
#              beta   = rep(.1, N_trees),
#              beta_t = rep(.1, N_years), 
#              sig_x_obs = .5,
#              sig_d_obs = .1,
#              beta_sd   = .1, 
#              beta_t_sd = .1, 
#              sig_x     = .3)
# inits = c(inits, b0 = 0.5, b1 = 10)
# 
# 
# if (fname_model == 'ring_model_t_date_taxon') {
#   inits = c(inits, list(b0 = 0.5, b1 = 10, beta_spp_sd = 0.1, beta_slope=0.1, beta_spp=rep(0.1, N_taxa)))
# }
# 
# if ((fname_model == 'ring_model_t_date_sapl_size_avg')|(fname_model == 'ring_model_t_date_sapl_size_constraint')) {
#   inits = c(inits, list(beta_slope=0.1))
# }
# 
# 
# m$setInits(inits)
# 
# spec <- configureMCMC(m, thin = 5) 
# spec$addMonitors(c('D', 'X', 'beta', 'beta_t', 'beta_t_sd', 'beta_sd', 'sig_x_obs', 'sig_d_obs'))
# 
# Rmcmc <- buildMCMC(spec)
# cm    <- compileNimble(m)
# Cmcmc <- compileNimble(Rmcmc, project = m)
# Cmcmc$run(10000) #25000
# out <- as.matrix(Cmcmc$mvSamples)#[101:500, ]
# 
# save(out, file = paste0('output/', fname_model, '_',  suff, '_nimble_', machine, '_', mvers, '.Rdata'))

# # full model
# fname_model = 'ring_model_t_date_size_pdbh_HMC'
# source(paste0('models/', fname_model, '.R'))
# m <- nimbleModel(body(fname_model),
#                  constants = list(N_inc = N_inc,
#                                   N_dbh = N_dbh,
#                                   N_pdbh = N_pdbh,
#                                   N_trees = N_trees,
#                                   N_years = N_years,
#                                   N_X = N_X,
#                                   N_taxa = N_taxa,
#                                   dbh_tree_id = dbh_tree_id,
#                                   dbh_year_id = dbh_year_id,
#                                   dbh_day_id = dbh_day_id,
#                                   pdbh_day_id = pdbh_day_id,
#                                   open_dbh = open_dbh,
#                                   m2tree = m2tree,
#                                   m2ti = m2ti,
#                                   m2nc = m2nc,
#                                   n1cores = n1cores,
#                                   n2cores=n2cores,
#                                   # n3cores=n3cores,
#                                   # n4cores=n4cores,
#                                   i1core2m = i1core2m,
#                                   i2core2m = i2core2m,
#                                   # i3core2m = i3core2m,
#                                   # i4core2m = i4core2m,
#                                   x2tree = x2tree, 
#                                   x2year = x2year,
#                                   ones = ones,
#                                   meas2d = meas2d,
#                                   meas2x = meas2x, 
#                                   pdbh2d = pdbh2d,
#                                   last_ti = last_ti,
#                                   first_ti = first_ti,
#                                   cs_last_ti = cs_last_ti,
#                                   # taxon = taxon,
#                                   # sap2x = sap2x,
#                                   # max_size = max_size,
#                                   # N_saplings = N_saplings,
#                                   tree_site_id = tree_site_id
#                                   
#                  ), check=FALSE)
# 
# 
# #debug(m$setData)
# 
# m$setData(list(logDobs = logDobs, logXobs = logXobs, logPDobs = logPDobs))#, sapling_constraint = rep(1, N_saplings)))
# 
# inits = list(X = rep(0.5, N_X), 
#              D0 = rep(0, N_trees),
#              # D = rep(.1, N_X), 
#              beta0 = .1, 
#              beta   = rep(.1, N_trees),
#              beta_t = rep(.1, N_years), 
#              sig_x_obs = .5,
#              sig_d_obs = .1,
#              beta_sd   = .1, 
#              beta_t_sd = .1, 
#              sig_x     = .3,
#              tau2      = -.1)#,
#              # tau3      = -.1,
#              # tau4      = -.1)
# inits = c(inits, b0 = 0.5, b1 = 10)
# inits = c(inits, list(beta_slope=0.1))
# 
# m$setInits(inits)
# 
# spec <- configureMCMC(m, thin = 5) 
# spec$addMonitors(c('D', 'X', 'beta', 'beta_t', 'beta_t_sd', 'beta_sd', 'sig_x_obs', 'tau2', 'sig_d_obs'))
# 
# Rmcmc <- buildMCMC(spec)
# cm    <- compileNimble(m)
# Cmcmc <- compileNimble(Rmcmc, project = m)
# Cmcmc$run(10000) #25000
# out <- as.matrix(Cmcmc$mvSamples)#[101:500, ]
# 
# save(out, file = paste0('output/', fname_model, '_v', mvers, '.Rdata'))
# 
# #########################################################################################################################################
# ## ring model with size and pdbh; no date
# #########################################################################################################################################
# 
# load(paste0('data/dump/', fname_data, '.rdata'))
# 
# # # check initial values
# # D_init = rep(0, N_X)
# # for (i in 1:N_trees){
# #   # for (j in 1:last_ti[i]){
# #   D_init[first_ti[i]:cs_last_ti[i]] = 0.0 + 2 * 0.5/ 10 * seq(1, last_ti[i])
# #   # }
# # }
# # 
# # D_init[sap2x] < max_size
# # all(D_init[sap2x] < max_size)
# 
# # full model
# fname_model = 'ring_model_t_size_pdbh_HMC'
# source(paste0('models/', fname_model, '.R'))
# m <- nimbleModel(body(fname_model),
#                  constants = list(N_inc = N_inc,
#                                   N_dbh = N_dbh,
#                                   N_pdbh = N_pdbh,
#                                   N_trees = N_trees,
#                                   N_years = N_years,
#                                   N_X = N_X,
#                                   N_taxa = N_taxa,
#                                   dbh_tree_id = dbh_tree_id,
#                                   dbh_year_id = dbh_year_id,
#                                   dbh_day_id = dbh_day_id,
#                                   pdbh_day_id = pdbh_day_id,
#                                   open_dbh = open_dbh,
#                                   m2tree = m2tree,
#                                   m2ti = m2ti,
#                                   m2nc = m2nc,
#                                   n1cores = n1cores,
#                                   n2cores=n2cores,
#                                   # n3cores=n3cores,
#                                   # n4cores=n4cores,
#                                   i1core2m = i1core2m,
#                                   i2core2m = i2core2m,
#                                   # i3core2m = i3core2m,
#                                   # i4core2m = i4core2m,
#                                   x2tree = x2tree, 
#                                   x2year = x2year,
#                                   ones = ones,
#                                   meas2d = meas2d,
#                                   meas2x = meas2x, 
#                                   pdbh2d = pdbh2d,
#                                   last_ti = last_ti,
#                                   first_ti = first_ti,
#                                   cs_last_ti = cs_last_ti,
#                                   # taxon = taxon,
#                                   # sap2x = sap2x,
#                                   # max_size = max_size,
#                                   # N_saplings = N_saplings,
#                                   tree_site_id = tree_site_id
#                                   
#                  ), check=FALSE)
# 
# 
# #debug(m$setData)
# 
# m$setData(list(logDobs = logDobs, logXobs = logXobs, logPDobs = logPDobs))#, sapling_constraint = rep(1, N_saplings)))
# 
# inits = list(X = rep(0.5, N_X), 
#              D0 = rep(0, N_trees),
#              # D = rep(.1, N_X), 
#              beta0 = .1, 
#              beta   = rep(.1, N_trees),
#              beta_t = rep(.1, N_years), 
#              sig_x_obs = .5,
#              sig_d_obs = .1,
#              beta_sd   = .1, 
#              beta_t_sd = .1, 
#              sig_x     = .3,
#              tau2      = -.1)#,
# # tau3      = -.1,
# # tau4      = -.1)
# inits = c(inits, b0 = 0.5, b1 = 10)
# inits = c(inits, list(beta_slope=0.1))
# 
# m$setInits(inits)
# 
# spec <- configureMCMC(m, thin = 5) 
# spec$addMonitors(c('D', 'X', 'beta', 'beta_t', 'beta_t_sd', 'beta_sd', 'sig_x_obs', 'tau2', 'sig_d_obs'))
# 
# Rmcmc <- buildMCMC(spec)
# cm    <- compileNimble(m)
# Cmcmc <- compileNimble(Rmcmc, project = m)
# Cmcmc$run(10000) #25000
# out <- as.matrix(Cmcmc$mvSamples)#[101:500, ]
# 
# save(out, file = paste0('output/', fname_model, '_v', mvers, '.Rdata'))
# 
# #########################################################################################################################################
# ## run without census data (same but no dbh, only pdbh)
# #########################################################################################################################################
# load(paste0('data/dump/tree_data_HMC_no_census_', dvers, '.rdata'))
# 
# # check initial values
# D_init = rep(0, N_X_p)
# for (i in 1:N_trees_p){
#   # for (j in 1:last_ti[i]){
#   D_init[first_ti_p[i]:cs_last_ti_p[i]] = 0.0 + 2 * 0.5/ 10 * seq(1, last_ti_p[i])
#   # }
# }
# 
# # D_init[sap2x] < max_size
# # all(D_init[sap2x] < max_size)
# 
# fname_model = 'ring_model_t_pdbh_HMC_nc'
# 
# x2idx_p = match(x2tree_p, unique(x2tree_p))
# 
# source(paste0('models/', fname_model, '.R'))
# # full model
# m <- nimbleModel(body(fname_model),
#                  constants = list(N_inc = N_inc,
#                                   N_pdbh = N_pdbh,
#                                   N_trees = N_trees_p,
#                                   N_years = N_years,
#                                   N_X = N_X_p,
#                                   N_taxa = N_taxa,
#                                   pdbh_day_id = pdbh_day_id,
#                                   open_dbh = open_dbh,
#                                   m2tree = m2tree,
#                                   m2ti = m2ti,
#                                   m2nc = m2nc,
#                                   n1cores = n1cores,
#                                   n2cores=n2cores,
#                                   i1core2m = i1core2m,
#                                   i2core2m = i2core2m,
#                                   x2tree = x2idx_p, 
#                                   x2year = x2year_p,
#                                   ones = ones,
#                                   meas2x = meas2x_p, 
#                                   pdbh2d = pdbh2d_p,
#                                   last_ti = last_ti_p,
#                                   first_ti = first_ti_p,
#                                   cs_last_ti = cs_last_ti_p,
#                                   #sap2x = sap2x,
#                                   #max_size = max_size,
#                                   tree_site_id = tree_site_id
#                                   #taxon = taxon,
#                                   #tree_site_id = tree_site_id
#                                   
#                  ), check=FALSE)
# 
# 
# #debug(m$setData)
# 
# m$setData(list(logXobs = logXobs, logPDobs = logPDobs))#, tau2_constraint=c(1), tau3_constraint=c(1)))#not_cens = rep(1, N_saplings)))#s
# 
# inits = list(X = rep(0.5, N_X_p), 
#              D0 = rep(0, N_trees_p),
#              # D = rep(.1, N_X_p), 
#              beta0 = .1, 
#              beta   = rep(.1, N_trees_p),
#              beta_t = rep(.1, N_years), 
#              sig_x_obs = .5,
#              sig_d_obs = .1,
#              beta_sd   = .1, 
#              beta_t_sd = .1, 
#              sig_x     = .3,
#              tau2      = -.1)
# # inits = c(inits, b0 = 0.5, b1 = 10)
# # inits = c(inits, list(beta_slope=0.1))
# 
# m$setInits(inits)
# 
# spec <- configureMCMC(m, thin = 5) 
# spec$addMonitors(c('D', 'X', 'beta', 'beta_t', 'beta_t_sd', 'beta_sd', 'sig_x_obs', 'tau2', 'sig_d_obs'))
# 
# Rmcmc <- buildMCMC(spec)
# cm    <- compileNimble(m)
# Cmcmc <- compileNimble(Rmcmc, project = m)
# Cmcmc$run(10000) #25000
# out <- as.matrix(Cmcmc$mvSamples)#[101:500, ]
# 
# save(out, file = paste0('output/', fname_model, '_', mvers, '.Rdata'))
# 
# # #########################################################################################################################################
# # ## run without census data (same but no dbh, no correction for date)
# # #########################################################################################################################################
# # # load('data/dump/tree_full.rdata')
# # load(paste0('data/dump/tree_data_20_no_census_', dvers, '.rdata'))
# # # load(paste0('data/dump/tree_full_20_no_census_v5.rdata'))
# # 
# # # check initial values
# # D_init = rep(0, N_X_p)
# # for (i in 1:N_trees_p){
# #   # for (j in 1:last_ti[i]){
# #   D_init[first_ti_p[i]:cs_last_ti_p[i]] = 0.0 + 2 * 0.5/ 10 * seq(1, last_ti_p[i])
# #   # }
# # }
# # 
# # # D_init[sap2x] < max_size
# # # all(D_init[sap2x] < max_size)
# # 
# # fname_model = 'ring_model_t_size_pdbh_nc'
# # 
# # x2idx_p = match(x2tree_p, unique(x2tree_p))
# # 
# # source(paste0('models/', fname_model, '.R'))
# # # full model
# # m <- nimbleModel(body(fname_model),
# #                  constants = list(N_inc = N_inc,
# #                                   N_pdbh = N_pdbh,
# #                                   N_trees = N_trees_p,
# #                                   N_years = N_years,
# #                                   N_X = N_X_p,
# #                                   N_taxa = N_taxa,
# #                                   pdbh_day_id = pdbh_day_id,
# #                                   open_dbh = open_dbh,
# #                                   m2tree = m2tree,
# #                                   m2ti = m2ti,
# #                                   m2nc = m2nc,
# #                                   n1cores = n1cores,
# #                                   n2cores=n2cores,
# #                                   n3cores=n3cores,
# #                                   n4cores=n4cores,
# #                                   i1core2m = i1core2m,
# #                                   i2core2m = i2core2m,
# #                                   i3core2m = i3core2m,
# #                                   i4core2m = i4core2m,
# #                                   x2tree = x2idx_p, 
# #                                   x2year = x2year_p,
# #                                   ones = ones,
# #                                   meas2x = meas2x_p, 
# #                                   pdbh2d = pdbh2d_p,
# #                                   last_ti = last_ti_p,
# #                                   first_ti = first_ti_p,
# #                                   cs_last_ti = cs_last_ti_p,
# #                                   #sap2x = sap2x,
# #                                   max_size = max_size,
# #                                   N_saplings = N_saplings,
# #                                   tree_site_id = tree_site_id
# #                                   #taxon = taxon,
# #                                   #tree_site_id = tree_site_id
# #                                   
# #                  ), check=FALSE)
# # 
# # 
# # #debug(m$setData)
# # 
# # m$setData(list(logXobs = logXobs, logPDobs = logPDobs))#, tau2_constraint=c(1), tau3_constraint=c(1)))#not_cens = rep(1, N_saplings)))#s
# # 
# # inits = list(X = rep(0.5, N_X_p), 
# #              D0 = rep(0, N_trees_p),
# #              # D = rep(.1, N_X_p), 
# #              beta0 = .1, 
# #              beta   = rep(.1, N_trees_p),
# #              beta_t = rep(.1, N_years), 
# #              sig_x_obs = .5,
# #              sig_d_obs = .1,
# #              beta_sd   = .1, 
# #              beta_t_sd = .1, 
# #              sig_x     = .3,
# #              tau2      = -.1,
# #              tau3      = -.1,
# #              tau4      = -.1)
# # inits = c(inits, list(beta_slope=0.1))
# # 
# # m$setInits(inits)
# # 
# # spec <- configureMCMC(m, thin = 5) 
# # spec$addMonitors(c('D', 'X', 'beta', 'beta_t', 'beta_t_sd', 'beta_sd', 'sig_x_obs', 'tau2', 'tau3', 'tau4', 'sig_d_obs'))
# # 
# # Rmcmc <- buildMCMC(spec)
# # cm    <- compileNimble(m)
# # Cmcmc <- compileNimble(Rmcmc, project = m)
# # Cmcmc$run(5000) #25000
# # out <- as.matrix(Cmcmc$mvSamples)#[101:500, ]
# # 
# # save(out, file = paste0('output/', fname_model, '_', mvers, '.Rdata'))
# # 
# #########################################################################################################################################
# ## run without census data (same but no dbh, no correction for date, no size)
# #########################################################################################################################################
# # load('data/dump/tree_full.rdata')
# load(paste0('data/dump/tree_data_20_no_census_', dvers, '.rdata'))
# # load(paste0('data/dump/tree_full_20_no_census_v5.rdata'))
# 
# # check initial values
# D_init = rep(0, N_X_p)
# for (i in 1:N_trees_p){
#   # for (j in 1:last_ti[i]){
#   D_init[first_ti_p[i]:cs_last_ti_p[i]] = 0.0 + 2 * 0.5/ 10 * seq(1, last_ti_p[i])
#   # }
# }
# 
# # D_init[sap2x] < max_size
# # all(D_init[sap2x] < max_size)
# 
# fname_model = 'ring_model_t_pdbh_nc'
# 
# x2idx_p = match(x2tree_p, unique(x2tree_p))
# 
# source(paste0('models/', fname_model, '.R'))
# # full model
# m <- nimbleModel(body(fname_model),
#                  constants = list(N_inc = N_inc,
#                                   N_pdbh = N_pdbh,
#                                   N_trees = N_trees_p,
#                                   N_years = N_years,
#                                   N_X = N_X_p,
#                                   N_taxa = N_taxa,
#                                   pdbh_day_id = pdbh_day_id,
#                                   open_dbh = open_dbh,
#                                   m2tree = m2tree,
#                                   m2ti = m2ti,
#                                   m2nc = m2nc,
#                                   n1cores = n1cores,
#                                   n2cores=n2cores,
#                                   n3cores=n3cores,
#                                   n4cores=n4cores,
#                                   i1core2m = i1core2m,
#                                   i2core2m = i2core2m,
#                                   i3core2m = i3core2m,
#                                   i4core2m = i4core2m,
#                                   x2tree = x2idx_p,
#                                   x2year = x2year_p,
#                                   ones = ones,
#                                   meas2x = meas2x_p,
#                                   pdbh2d = pdbh2d_p,
#                                   last_ti = last_ti_p,
#                                   first_ti = first_ti_p,
#                                   cs_last_ti = cs_last_ti_p,
#                                   max_size = max_size,
#                                   tree_site_id = tree_site_id
#                  ), check=FALSE)
# 
# 
# #debug(m$setData)
# 
# m$setData(list(logXobs = logXobs, logPDobs = logPDobs))#, tau2_constraint=c(1), tau3_constraint=c(1)))#not_cens = rep(1, N_saplings)))#s
# 
# inits = list(X = rep(0.5, N_X_p),
#              D0 = rep(0, N_trees_p),
#              # D = rep(.1, N_X_p),
#              beta0 = .1,
#              beta   = rep(.1, N_trees_p),
#              beta_t = rep(.1, N_years),
#              sig_x_obs = .5,
#              sig_d_obs = .1,
#              beta_sd   = .1,
#              beta_t_sd = .1,
#              sig_x     = .3,
#              tau2      = -.1)
# m$setInits(inits)
# 
# spec <- configureMCMC(m, thin = 5)
# spec$addMonitors(c('D', 'X', 'beta', 'beta_t', 'beta_t_sd', 'beta_sd', 'sig_x_obs', 'tau2', 'sig_d_obs'))
# 
# Rmcmc <- buildMCMC(spec)
# cm    <- compileNimble(m)
# Cmcmc <- compileNimble(Rmcmc, project = m)
# Cmcmc$run(10000) #25000
# out <- as.matrix(Cmcmc$mvSamples)#[101:500, ]
# 
# save(out, file = paste0('output/', fname_model, '_', mvers, '.Rdata'))
