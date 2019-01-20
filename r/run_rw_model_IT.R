library(nimble)

# source('config_HMC')
dvers="v0.1"
mvers="v0.1"

fname_data = paste0('tree_data_IT_', dvers)

load(paste0('data/dump/', fname_data, '.rdata'))


#######################################################################################################################################
# full model but with zero covariance; not efficient
#######################################################################################################################################

fname_model = 'ring_model_t_date_sapl_size_pdbh_NOCOVAR'
load(paste0('output/', fname_model, '_v4.0.Rdata'))

col_names = sapply(strsplit(colnames(out), '\\['), function(x) x[[1]])
hist(out[,which(col_names=="sig_d_obs")])
sig_d_obs = mean(out[,which(col_names=="sig_d_obs")])

#
fname_model = 'ring_model_t_pdbh_IT'
load(paste0('data/dump/', fname_data, '.rdata'))
source(paste0('models/', fname_model, '.R'))

m <- nimbleModel(body(fname_model),
                 constants = list(N_inc = N_inc,
                                  N_pdbh = N_pdbh,
                                  N_trees = N_trees,
                                  N_years = N_years,
                                  N_X = N_X,
                                  N_taxa = N_taxa,
                                  m2tree = m2tree,
                                  m2ti = m2ti,
                                  x2tree = x2tree, 
                                  x2year = x2year,
                                  ones = ones,
                                  meas2x = meas2x, 
                                  pdbh2d = pdbh2d,
                                  last_ti = last_ti,
                                  first_ti = first_ti,
                                  cs_last_ti = cs_last_ti,
                                  taxon = taxon,
                                  sig_d_obs = sig_d_obs
                 ), check=FALSE)


#debug(m$setData)

m$setData(list(logXobs = logXobs, logPDobs = logPDobs))

inits = list(X = rep(0.5, N_X), 
             D0 = rep(0, N_trees),
             # D = rep(.1, N_X), 
             beta0 = .1, 
             beta   = rep(.1, N_trees),
             beta_t = rep(.1, N_years), 
             sig_x_obs = .5,
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

save(out, file = paste0('output/', fname_model, '_', mvers, '.Rdata'))


#######################################################################################################################################
# full model but with zero covariance; not efficient
#######################################################################################################################################

#
fname_model = 'ring_model_t_pdbh_IT'
load(paste0('data/dump/', fname_data, '.rdata'))
source(paste0('models/', fname_model, '.R'))

m <- nimbleModel(body(fname_model),
                 constants = list(N_inc = N_inc,
                                  N_pdbh = N_pdbh,
                                  N_trees = N_trees,
                                  N_years = N_years,
                                  N_X = N_X,
                                  N_taxa = N_taxa,
                                  m2tree = m2tree,
                                  m2ti = m2ti,
                                  x2tree = x2tree, 
                                  x2year = x2year,
                                  ones = ones,
                                  meas2x = meas2x, 
                                  pdbh2d = pdbh2d,
                                  last_ti = last_ti,
                                  first_ti = first_ti,
                                  cs_last_ti = cs_last_ti,
                                  taxon = taxon#,
                                  #sig_d_obs = sig_d_obs
                 ), check=FALSE)


#debug(m$setData)

m$setData(list(logXobs = logXobs, logPDobs = logPDobs))

inits = list(X = rep(0.5, N_X), 
             D0 = rep(0, N_trees),
             # D = rep(.1, N_X), 
             beta0 = .1, 
             beta   = rep(.1, N_trees),
             beta_t = rep(.1, N_years), 
             sig_x_obs = .5,
             sig_x_obs = .5,
             beta_sd   = .1, 
             beta_t_sd = .1, 
             sig_x     = .3)
m$setInits(inits)

spec <- configureMCMC(m, thin = 5) 
spec$addMonitors(c('D', 'X', 'beta', 'beta_t', 'beta_t_sd', 'beta_sd', 'sig_x_obs', 'sig_d_obs'))

Rmcmc <- buildMCMC(spec)
cm    <- compileNimble(m)
Cmcmc <- compileNimble(Rmcmc, project = m)
Cmcmc$run(10000) #25000
out <- as.matrix(Cmcmc$mvSamples)#[101:500, ]

save(out, file = paste0('output/', fname_model, '_', mvers, '.Rdata'))
