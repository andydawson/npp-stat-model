
source("config")

load(paste0('data/lyford20/incr_data_', dvers ,'.Rda'))

# work with incr
incrData$incr[incrData$incr == 0] = 0.000001
incr_log = log(incrData$incr)

incrData = data.frame(incrData, incr_log = incr_log)
angles = data.frame(dir=c('e', 'n', 'w', 's'), angle=c(0, 90, 180, 360))
ids = unique(incrData$id)

incrMean = data.frame(id   = numeric(0), 
                      year = numeric(0), 
                      angle = numeric(0),
                      incr1 = numeric(0), 
                      incr2 = numeric(0), 
                      incr_mean = numeric(0), 
                      incr_var = numeric(0), 
                      stat_id  = numeric(0), 
                      nmeas = numeric(0),
                      ncores   = numeric(0),
                      core_one = numeric(0),
                      core_two = numeric(0),
                      core_three = numeric(0),
                      core_four  = numeric(0))

for (i in 1:length(ids)){
  print(i)
  id = ids[i]
  
  years = unique(incrData$year[which(incrData$id == id)])
  
  for (year in years) {
    
    dat = incrData[(incrData$id == id) & (incrData$year == year),]
    dat = dat[dat$orient %in% c('n', 's', 'e', 'w'), ] 
    
    if (nrow(dat) < 2) {
      next
    }
    
    stat_id = dat$stat_id[1]
    nmeas   = nrow(dat)
    combs   = combn(seq(1,nmeas), 2)
    
    comb_sep = rep(NA, ncol(combs))
    for (i in 1:ncol(combs)){
      comb_sep[i] = abs(diff(angles[angles$dir %in% dat$orient[combs[,i]], 2])) %% 180
    }
    
    for (i in 1:ncol(combs)) {
      
      comb = combs[,i]
      
      incrMean = rbind(incrMean, 
                       data.frame(id       = id, 
                                  year     = year, 
                                  angle    = comb_sep[i],
                                  incr1    = dat$incr_log[comb[1]],
                                  incr2    = dat$incr_log[comb[2]],
                                  incr_mean = mean(dat$incr_log[comb]),
                                  incr_var  = var(dat$incr[comb]), 
                                  stat_id   = stat_id, 
                                  nmeas     = nmeas,
                                  ncores    = 2,
                                  core_one  = 1 %in% comb,
                                  core_two  = 2 %in% comb,
                                  core_three = 3 %in% comb,
                                  core_four  = 4 %in% comb)
      )
      
    }
  }
}

# take differences
years = sort(unique(incrMean$year))
ids   = as.numeric(levels(unique(incrMean$id))) # convert from "factor" to numeric list

incr_diff = data.frame(id    = numeric(0), 
                       year  = numeric(0), 
                       angle = numeric(0),
                       diff  = numeric(0), 
                       stat_id  = numeric(0), 
                       nmeas    = numeric(0),
                       ncores   = numeric(0))#,
# common   = numeric(0))


for (i in 1:nrow(incrMean)){
  print(i)
  
  incr_dat1 = data.frame(incrMean[i,c(1,2,3)], diff=incrMean[i,'incr_mean'] - incrMean[i,'incr1'], incrMean[i,c(8,9,10)])
  incr_dat2 = data.frame(incrMean[i,c(1,2,3)], diff=incrMean[i,'incr_mean'] - incrMean[i,'incr2'], incrMean[i,c(8,9,10)])
  
  incr_diff = rbind(incr_diff, incr_dat1, incr_dat2)
  
}

# using all increments
diffall = incr_diff[which(incr_diff$nmeas==2),]
eall = var(diffall$diff)
eall
median((diffall$diff - mean(diffall$diff))^2)

pdf('figures/abs_diff_means_2core.pdf')
hist(abs(diffall$diff), breaks=100)
dev.off()

#mad
eall = 1.4826*median(abs(diffall$diff - median(diffall$diff)))
eall

# using all increments
diffall = incr_diff[which(incr_diff$nmeas==3),]
eall = var(diffall$diff)
eall
median((diffall$diff - mean(diffall$diff))^2)

pdf('figures/abs_diff_means_3core.pdf')
hist(abs(diffall$diff), breaks=100)
dev.off()

#mad
eall = 1.4826*median(abs(diffall$diff - median(diffall$diff)))
eall

sig_x_obs = 0.6

tau2 = sig_x_obs^2 - 2*eall 
tau2

##
## think about angle of separation
##

diff90 = incr_diff[incr_diff$angle == 90,]
e90 = var(diff90$diff)
e90
median((diff90$diff - mean(diff90$diff))^2)

#mad
e90 = 1.4826*median(abs(diff90$diff - median(diff90$diff)))

diff180 = incr_diff[incr_diff$angle == 0,]
e180 = var(diff180$diff)
e180

median((diff180$diff - mean(diff180$diff))^2)

#mad
e180 = 1.4826*median(abs(diff180$diff - median(diff180$diff)))

sig_x_obs = 0.5

tau2 = sig_x_obs^2 - 2*eall 
tau2

tau2 = sig_x_obs^2 - 2*e180 
tau2

tau2 = sig_x_obs^2 - 2*e90 
tau2

# interpret tau's as the correlation
tau2 = 1 - 2 * eall / (sig_x_obs^2) 
tau2

# interpret tau's as the correlation
tau2 = 1 - 2 * e180 / (sig_x_obs^2) 
tau2

# interpret tau's as the correlation
tau2 = 1 - 2 * e90 / (sig_x_obs^2) 
tau2

pdf('figures/diffs_by_angle.pdf')
par(mfrow=c(2,1))
hist(diff90$diff, breaks=60, xlim=c(-6,6), freq=FALSE)
hist(diff180$diff, breaks=60, xlim=c(-6,6), freq=FALSE)
dev.off()

pdf('figures/tau2_by_angle.pdf')
par(mfrow=c(1,1))
sig = seq(0.32,2,by=0.01)
plot(sig, 1 - 2 * e90 / (sig^2), type='l', ylim=c(-1,1), xlab='sigma', ylab='tau2') 
lines(sig, 1 - 2 * e180 / (sig^2), col='blue')
# lines(sig, 1 - 2 * e21b / (sig^2), col='red')
# lines(sig, 1 - 6 * e32b/ (sig^2) , col='green')
legend('bottomright', c('90', '180', 'All'), lty=c(1,1), col=c('black', 'blue', 'red'))
dev.off()

diff90 = incr_diff[which((incr_diff$angle == 90) & (incr_diff$nmeas==2)),]
e90 = var(diff90$diff)
e90
median((diff90$diff - mean(diff90$diff))^2)

#mad
e90 = 1.4826*median(abs(diff90$diff - median(diff90$diff)))
e90

diff180 = incr_diff[which((incr_diff$angle == 0) & (incr_diff$nmeas==2)),]
e180 = var(diff180$diff)
e180

median((diff180$diff - mean(diff180$diff))^2)

#mad
e180 = 1.4826*median(abs(diff180$diff - median(diff180$diff)))
e180

sig_x_obs = 0.5

tau2 = sig_x_obs^2 - 2*e180 
tau2

tau2 = sig_x_obs^2 - 2*e90 
tau2

# interpret tau's as the correlation
tau2 = 1 - 2 * e180 / (sig_x_obs^2) 
tau2

# interpret tau's as the correlation
tau2 = 1 - 2 * e90 / (sig_x_obs^2) 
tau2

pdf('figures/diffs_by_angle_2core.pdf')
par(mfrow=c(2,1))
hist(diff90$diff, breaks=60, xlim=c(-6,6), freq=FALSE)
hist(diff180$diff, breaks=60, xlim=c(-6,6), freq=FALSE)
dev.off()

pdf('figures/tau2_by_angle_2core.pdf')
sig = seq(0.32,2,by=0.01)
plot(sig, 1 - 2 * e90 / (sig^2), type='l', ylim=c(-1,1), xlab='sigma', ylab='tau2') 
lines(sig, 1 - 2 * e180 / (sig^2), col='blue')
lines(sig, 1 - 2 * e21b / (sig^2), col='red')
# lines(sig, 1 - 6 * e32b/ (sig^2) , col='green')
legend('bottomright', c('90', '180', 'All'), lty=c(1,1), col=c('black', 'blue', 'red'))
dev.off()
