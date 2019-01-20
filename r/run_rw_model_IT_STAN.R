library(rstan)

# source('config_HMC')
dvers="v0.1"
mvers="v0.1"
site = "IT"

fname_data = paste0('tree_data_IT_STAN_', dvers)
fname_model = "ring_model_t_pdbh_IT_stan_v2"

dat = readRDS(paste0('data/dump/', fname_data, '.RDS'))


#######################################################################################################################################
# full model but with zero covariance; not efficient
#######################################################################################################################################

compiled <- stan_model(file = 'models/stan/ring_model_t_pdbh_IT_v2.stan')

fit <- sampling(compiled, 
                data = dat, 
                iter = 200, 
                chains = 1)
rm(compiled)

post=rstan::extract(fit)
rm(fit)


save(post, file = paste0('output/', fname_model, '_', site, '_', mvers, '.Rdata'))
