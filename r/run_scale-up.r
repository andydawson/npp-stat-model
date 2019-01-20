library(nimble)
library(reshape)
library(ggplot2)

source('config')
# source('config_HMC')

figures_dir = 'NOCOVAR'

get_cols <- function(x, varname){
  col_names = sapply(strsplit(colnames(x), '\\['), function(x) x[[1]])
  return(which(col_names == varname))
}

compile_model <- function(thin=thin, b=b){
  # first run for log(ab)
  source(paste0('models/scale_up.R'))
  # full model
  m <- nimbleModel(body(scale_up),
                   constants = list(N_sites = 3,
                                    T = T,
                                    sig_mu = 100), check=FALSE)
  
  #debug(m$setData)
  m$setData(list(logb = log(b)))
  
  inits = list(mu = rep(0.5, T), 
               sigma = 0.5)
  
  m$setInits(inits)
  spec <- configureMCMC(m, thin = thin) 
  spec$addMonitors(c('mu', 'sigma'))
  Rmcmc <- buildMCMC(spec)
  cm    <- compileNimble(m)
  Cmcmc <- compileNimble(Rmcmc, project = m)
  
  return(list(compiled_model = cm, compiled_mcmc = Cmcmc))
}


run_scale_model <- function(compiled, b, thin=thin, niter_scale=niter_scale){
  
  compiled_model = compiled$compiled_model
  compiled_mcmc = compiled$compiled_mcmc
  
  # compiled_mcmc$reset()
  compiled_model$setData(list(logb = log(b)))
  
  compiled_mcmc$run(niter_scale, reset=TRUE) 
  post <- as.matrix(compiled_mcmc$mvSamples)
  
  return(post)
}

if (location=='HF'){
  fname_data = paste0('tree_data_20_', dvers)
} else if (location=='HMC'){
  fname_data = paste0('tree_data_HMC_', dvers)
}
load(file=paste0('data/dump/', fname_data, '.rdata'))

if (location=='HF'){
  fname_data = paste0('tree_data_20_no_census_', dvers)
} else if (location=='HMC'){
  fname_data = paste0('tree_data_HMC_no_census_', dvers)
}
load(file=paste0('data/dump/', fname_data, '.rdata'))

ab  = readRDS(file=paste0('allom/ab_p_', location, '_', mvers, '.RDS'))
colnames(ab)[which(colnames(ab) == 'ab')] = 'value'

abi = readRDS(file=paste0('allom/abi_p_', location, '_', mvers, '.RDS'))
colnames(abi)[which(colnames(abi) == 'abi')] = 'value'

dat = list(ab=ab, abi=abi)

for (var in c('ab', 'abi')){
  
  bio=dat[[var]]
  
  bio=bio[which(bio$iter<20),]
  
  niter=length(unique(bio$iter))
  niter_scale = 1000
  thin = 5
  nsamps= niter_scale/thin
  burn = 0#250
  nkeep = 10
  
  T = 53#length(years)-2
  models  = as.vector(unique(bio$model))
  nmodels = length(models)
  mu_post = array(NA, c(nmodels, niter*nkeep, T))
  
  for (n in 1:nmodels){
    bio_model = bio[which(bio$model==models[n]),]
    for (i in 1:niter){
      # organize the data
      bio_agg = aggregate(value~year_id+site_id+iter, bio_model[which(bio_model$iter==i),], function(x) sum(x, na.rm=TRUE))
      bio_agg = bio_agg[,c('year_id', 'site_id', 'value')]
      bio_wide = reshape(bio_agg, timevar="site_id", idvar = c("year_id"), direction="wide")
      
      b = t(bio_wide[,2:4])
      
      # remove the last year
      b = b[,1:T]
      
      if (i==1){
        compiled <- compile_model(thin=thin, b=b)
      }
      
      post <- run_scale_model(compiled, b, thin=thin, niter_scale=niter_scale)
      
      post <- post[burn:nsamps,]
      mu_post[n,((i-1)*nkeep+1):(i*nkeep),] = post[sample.int(nrow(post),nkeep),get_cols(post,'mu')]
    }
  }
  
  mu_post_melt = melt(mu_post)
  colnames(mu_post_melt) = c('model', 'iter', 'year', 'value')
  mu_post_melt$model = models[mu_post_melt$model]
  mu_post_melt$year = years[mu_post_melt$year]
  mu_post_melt$value = exp(mu_post_melt$value)
  
  saveRDS(mu_post_melt, file=paste0('allom/', var, '_scale_up_samples_', location, '.RDS'))
}


mu_ab = readRDS(paste0('allom/ab_scale_up_samples_', location, '.RDS'))
mu_quants = aggregate(value~year+model, mu_ab, function(x) quantile(x, c(0.025, 0.5, 0.975)))
mu_quants = data.frame(mu_quants[,1:2], mu_quants[,3])
colnames(mu_quants) = c('year', 'model', 'lb', 'median', 'ub')

# ab_census = readRDS(file=paste0('data/ab_c_melt_', location, '.RDS'))
# 
# ab_census_w = reshape(ab_census, timevar='quant', idvar=c('year'), direction='wide')
# colnames(ab_census_w) = c('year', 'lb', 'median', 'ub')

ggplot() + geom_ribbon(data=mu_quants, aes(x=year, ymin = lb, ymax = ub, fill=model), alpha=0.4) + 
  geom_line(data=mu_quants, aes(x=year, y=median, group=model, colour=model)) +  
  # geom_point(data=subset(ab_census,quant %in% c('median')), aes(x=year, y=value)) +
  # geom_errorbar(data=ab_census_w, aes(x=year, ymin=lb, ymax=ub),width=1.5) +
  theme_bw()+
  ylab("Biomass (Mg/ha)") + xlab('Year') +
  scale_x_continuous(breaks=seq(min(years[seq(1,T)]), max(years[seq(1,T)]), by=5)) + 
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14))
ggsave(paste0('figures/', figures_dir, '/AGB_scaled_up_', location, '.pdf'))

mu_abi = readRDS(paste0('allom/abi_scale_up_samples_', location, '.RDS'))
mu_quants = aggregate(value~year+model, mu_abi, function(x) quantile(x, c(0.025, 0.5, 0.975)))
mu_quants = data.frame(mu_quants[,1:2], mu_quants[,3])
colnames(mu_quants) = c('year', 'model', 'lb', 'median', 'ub')

ggplot() + geom_ribbon(data=mu_quants, aes(x=year, ymin = lb, ymax = ub, fill=model), alpha=0.4) + 
  geom_line(data=mu_quants, aes(x=year, y=median, group=model, colour=model)) +  
  theme_bw()+
  ylab("Biomass Increment (Mg/ha)") + xlab('Year') +
  scale_x_continuous(breaks=seq(min(years[seq(1,T)]), max(years[seq(1,T)]), by=5)) + 
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14))
ggsave(paste0('figures/', figures_dir, '/ABI_scaled_up_', location, '.pdf'))


# plot(mu_post[,8], type='l')
# 
# mu_post = post[,get_cols(post,'mu')]
# mu = apply(mu_post, 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
# colnames(mu) = seq(1, T)
# 
# plot(mu_post[,10], type='l')
# 
# mu_melt = melt(mu)
# colnames(mu_melt) = c('quant', 'year_id', 'mu')
# 
# mu_rib = data.frame(mu_melt[mu_melt$quant == '2.5%', c('year_id', 'mu')], mu_melt[mu_melt$quant == '97.5%', c('mu')])
# colnames(mu_rib) = c('year_id', 'mu_lb', 'mu_ub')
# 
# ggplot(mu_melt) + geom_ribbon(data=mu_rib, aes(x=year_id, ymin = mu_lb, ymax = mu_ub), fill = "grey90") + 
#   geom_line(aes(x=year_id, y=mu, colour=quant)) + scale_colour_manual(values=c('dodgerblue', 'black', 'dodgerblue'))
# 
# sigma_post = post[,get_cols(post,'sigma')]
# plot(sigma_post, type='l')
