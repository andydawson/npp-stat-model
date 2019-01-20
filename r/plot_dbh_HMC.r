library(reshape2)
library(ggplot2)

source('r/utils/plot_funs.r')
source('config_HMC')

figures_dir = 'NOCOVAR'
if (!file.exists(figures_dir)){
  dir.create(paste0("figures/", figures_dir))
}

# fnames = c('ring_model_t_date_size_pdbh_HMC', 'ring_model_t_pdbh_HMC_nc') 
fnames = c('ring_model_t_pdbh_HMC_NOCOVAR', 'ring_model_t_pdbh_nc_HMC_NOCOVAR_sigd')
models = c('Model RW + Census', 'Model RW')

fname_data = paste0('tree_data_HMC_', dvers)
load(file=paste0('data/dump/', fname_data, '.rdata'))
fname_data = paste0('tree_data_HMC_no_census_', dvers)
load(file=paste0('data/dump/', fname_data, '.rdata'))

post= list()
for (i in 1:length(fnames)) {
  fname_model = fnames[i]
  load(file   = paste0('output/', fname_model, '_', mvers, '.Rdata'))
  post[[i]]   = out#[1:500,]
}  

####

burn=25
niter=dim(out)[1]

########################################################################################################################################
library(RColorBrewer)
darkcols <- brewer.pal(4, "Set1")

nc=FALSE
x2idx_p = match(x2tree_p, unique(x2tree_p))
trees = sort(unique(x2tree))
trees_p = sort(unique(x2tree_p))

pdf(paste0('figures/', figures_dir, '/growth_results_both_HMC_', mvers, '.pdf'), width=12, height=10)
for (i in 1:N_trees){
  
  print(paste0('Tree ', i))
  
  tree = trees[i]
  tree_p = trees_p[trees_p==tree]
  
  tree_idx = which(x2tree == tree)
  tree_idx_p =  which(x2tree_p == tree_p)
  
  # estimated increment
  X_qs = apply(post[[1]][,get_cols(post[[1]],'X')[tree_idx]], 2, 
               function(x) quantile(x, probs=c(0.025,0.5, 0.975), na.rm=TRUE))
  X_qs_p = apply(post[[2]][,get_cols(post[[2]],'X')[tree_idx_p]], 2, 
                 function(x) quantile(x, probs=c(0.025,0.5, 0.975), na.rm=TRUE))
  
  
  tree_years = x2year[tree_idx]
  tree_years_p = x2year_p[tree_idx_p]
  
  # raw data
  dat_idx   = which(m2tree == tree)
  core_nums = unique(m2orient[dat_idx])
  rws = exp(logXobs[dat_idx])
  
  if (length(rws)>0){
    ymin = min(c(rws,X_qs[1,]))
    ymax = max(c(rws,X_qs[3,]))
  } else {
    ymin = min(X_qs[1,])
    ymax = max(X_qs[3,])
  }
  
  
  tree_dat = data.frame(value=numeric(0), year=numeric(0), var=character(0), type=character(0), subtype=character(0))
  tree_dat = rbind(tree_dat, data.frame(value=X_qs[2,], year=years[tree_years], var=rep('RW'), type=rep('model'), subtype=rep('Both')))
  
  ribbon_dat = data.frame(L=numeric(0), U=numeric(0), year=numeric(0), var=character(0), subtype=character(0))
  ribbon_dat = rbind(ribbon_dat, data.frame(L=X_qs[1,], U=X_qs[3,], year=years[tree_years], var=rep('RW'), subtype=rep('Both')))
  
  if (length(tree_idx_p)>0){
    tree_dat = rbind(tree_dat, data.frame(value=X_qs_p[2,], year=years[tree_years_p], var=rep('RW'), type=rep('model'), subtype=rep('RW')))
    ribbon_dat = rbind(ribbon_dat, data.frame(L=X_qs_p[1,], U=X_qs_p[3,], year=years[tree_years_p], var=rep('RW'), subtype=rep('RW')))
  }
  
  if (length(core_nums) > 1) {
    idx_a = which(m2tree_a == tree)
    dat   = exp(logXobs_a[idx_a])
    yrs   = m2t_a[idx_a] 
    
    tree_dat = rbind(tree_dat, data.frame(value=dat, year=yrs, var=rep('RW'), type=rep('data'), subtype=rep('raw avg')))
  }
  
  for (core in core_nums){
    idx = which((m2tree == tree) & (m2orient == core))
    yrs = m2t[idx] 
    tree_dat = rbind(tree_dat, data.frame(value=exp(logXobs[idx]), year=yrs, var=rep('RW'), type=rep('data'), subtype=rep(core)))
  }
  
  # now plot D!
  D_qs = apply(post[[1]][,get_cols(post[[1]],'D')[tree_idx]], 2, function(x) quantile(x, probs=c(0.025,0.5, 0.975), na.rm=TRUE))
  D_qs_p = apply(post[[2]][,get_cols(post[[2]],'D')[tree_idx_p]], 2, function(x) quantile(x, probs=c(0.025,0.5, 0.975), na.rm=TRUE))
  
  tree_dat = rbind(tree_dat, data.frame(value=D_qs[2,], year=years[tree_years], var=rep('DBH'), type=rep('model'), subtype=rep('Both')))
  
  ribbon_dat = rbind(ribbon_dat, data.frame(L=D_qs[1,], U=D_qs[3,], year=years[tree_years], var=rep('DBH'), subtype=rep('Both')))
  
  if (length(tree_idx_p)>0){
    tree_dat = rbind(tree_dat, data.frame(value=D_qs_p[2,], year=years[tree_years_p], var=rep('DBH'), type=rep('model'), subtype=rep('RW')))
    ribbon_dat = rbind(ribbon_dat, data.frame(L=D_qs_p[1,], U=D_qs_p[3,], year=years[tree_years_p], var=rep('DBH'), subtype=rep('RW')))
  }
  
  D_dat = NA
  if (!nc){
    idx_dbh = which(dbh_tree_id == tree)
    yrs = dbh_year_id[idx_dbh]
    
    D_dat = exp(logDobs[idx_dbh])
  }
  
  ymax = max(c(D_dat,D_qs[2,]), na.rm=TRUE)
  
  tree_dat = rbind(tree_dat, data.frame(value=D_dat, year=years[yrs], var=rep('DBH'), type=rep('data'), subtype=rep('census')))
  
  
  if (any(pdbh$stat_id == tree)){
    
    idx_pdbh = which(pdbh$stat_id == tree)
    yrs = pdbh_year_id[idx_pdbh]
    
    PD_dat = exp(logPDobs[idx_pdbh])
    
    ymax = max(c(PD_dat,ymax))
    
    tree_dat = rbind(tree_dat, data.frame(value=PD_dat, year=years[yrs], var=rep('DBH'), type=rep('data'), subtype=rep('paleon')))
    
  } else {
    tree_dat = rbind(tree_dat, data.frame(value=NA, year=NA, var=rep('DBH'), type=rep('data'), subtype=rep('paleon')))
    
  }
  
  cols = c('#084594', '#8c2d04')
  cols_fill = c('#4292c6', 'coral2')
  
  census_id = dbh[which(dbh$stat_id == tree), 'stemindex']
  hf_id = dbh[which(dbh$stat_id == tree), 'tree_id']
  
  p <- ggplot(tree_dat) + geom_ribbon(data=ribbon_dat, aes(x=year, ymin=L, ymax=U, fill=subtype), alpha=0.4) + 
    geom_line(data=subset(tree_dat, type %in% c('model')), aes(x=year, y=value, colour=subtype), size=1) + 
    geom_point(data=subset(tree_dat, (type %in% c('data')) & (var %in% c('RW')) & (subtype %in% c('raw avg'))), 
               aes(x=year, y=value, group=subtype), colour='brown', size=4, shape=20, alpha=0.5) +
    geom_line(data=subset(tree_dat, (type %in% c('data')) & (var %in% c('RW')) & (!(subtype %in% c('raw avg')))), 
              aes(x=year, y=value, group=subtype), alpha=0.7,  colour='black', linetype=2, size=0.8, show.legend=FALSE) +
    geom_point(data=subset(tree_dat, (type %in% c('data')) & (var %in% c('DBH'))) , 
               aes(x=year, y=value, shape=subtype), size=3) +
    # xlim(c(1960, 2015)) + 
    scale_color_manual(values=cols, name='Data', labels=c('RW + Census', 'RW')) + 
    scale_fill_manual(values=cols_fill, name='Data', labels=c('RW + Census', 'RW')) + 
    scale_shape_manual(values=c(19, 8, 10), guide='none') +
    theme_bw()+
    theme(axis.title.y=element_blank()) + 
    scale_x_continuous(breaks=seq(1960, 2015, by=5)) + 
    facet_grid(var~., scales="free_y") + 
    ggtitle(paste0('Stat id: ', i , '; Census id: ', census_id, '; PalEON id :', hf_id)) #+ 
  print(p)
}
dev.off()

####################################################################################################################################
## trace and density plots of parameters
####################################################################################################################################
burn=200
plot_sig(post, burn, figure_dir, location, mvers)
#########################################################################################################################################

b0 = sapply(post, function(x) x[,get_cols(x,'b0')][burn:nrow(x)])
colMeans(b0)
b0 = melt(b0)
colnames(b0) =c('iter', 'model', 'value')

b1 = sapply(post, function(x) x[,get_cols(x,'b1')][burn:nrow(x)])
colMeans(b1)
b1 = melt(b1)
colnames(b1) =c('iter', 'model', 'value')

b = rbind(data.frame(b0, par=rep('b0')),
          data.frame(b1, par=rep('b1')))

b$model = models[b$model]

ggplot(data=b) + geom_line(aes(x=iter, y=value, colour=factor(model))) + facet_grid(par~., scales="free_y")+ 
  labs(colour='Data')
ggsave(file=paste0('figures/b_trace_', suff ,'.pdf'))


#########################################################################################################################################

# tau2 = sapply(post, function(x) x[,get_cols(x,'tau2')][burn:nrow(x)])
# colMeans(tau2)
# tau2 = melt(tau2)
# colnames(tau2) = c('iter', 'model', 'value')
# 
# tau = data.frame(tau2, par=rep('tau2'))
# 
# tau$model = models[tau$model]
# 
# ggplot(data=tau) + geom_line(aes(x=iter, y=value, colour=factor(model))) + facet_grid(par~., scales="free_y")+ 
#   labs(colour='Data')
# ggsave(file=paste0('figures/tau_trace_', suff, '.pdf'))
# 
# 
# ggplot(data=tau, aes(x=value, y=..scaled..)) + #geom_histogram(aes(value, colour=model, fill=model), binwidth=0.02, stat='bin') + 
#   geom_density(aes(colour=model, fill=model),alpha=.2)+ facet_grid(par~., scales="free_y") + 
#   labs(colour='Method') + 
#   labs(fill='Method') + ylab('Density') + xlab('Value') + 
#   theme(strip.text = element_text(size = 12), legend.text = element_text(size=12),legend.title = element_text(size=12), axis.title = element_text(size=12), axis.text = element_text(size=12))
# ggsave(file=paste0('figures/tau_density.pdf'))
# ggsave(file=paste0('figures/tau_density.png'))

