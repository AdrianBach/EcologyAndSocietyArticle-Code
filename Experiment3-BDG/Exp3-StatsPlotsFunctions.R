# Adrian Bach
# PhD : Using AI to improve decision-making in conservation conflicts
# University of Stirling

#### libraries ####
library(plotly)

#### functions needed ####

boot_sd_ci <- function(x, confidence = 95, itr = 1000) {
  
  # init iterations and sample
  i <- itr
  bt_avg <- NULL
  
  # loop over iterations
  while (i > 0) {
    # sample randomly from x
    spl <- sample(x, length(x), replace = TRUE)
    
    # store mean
    bt_avg <- c(bt_avg, sum(spl)/length(spl))
    
    # decrement i
    i <- i-1
  }
  
  # mean over the bootstrapped samples
  bt_est <- sum(bt_avg)/itr
  
  # compute standard deviation
  bt_sd <- sqrt((1/(length(x)-1)) * sum((bt_avg-bt_est)^2))
  
  # compute confidence interval
  # sorting bt_avg numerically
  st_avg <- sort(bt_avg)
  
  # get the first value after the 2.5 first centiles
  bt_95ci_inf <- st_avg[floor(0.5*(1-0.01*confidence)*itr)+1]
  
  # get the last value before the 2.5 last centiles
  bt_95ci_sup <- st_avg[floor((0.01*confidence+0.5*(1-0.01*confidence))*itr)-1]
  
  res <- c(bt_sd, bt_95ci_inf, bt_95ci_sup)
  return(res) 
  
}

OTI_stats <- function(df, ts, omit.extinction = FALSE) {
  
  df <- as.data.frame(df)
  
  # levels of ratio, UT and BB parameters 
  bgt_rat <- levels(as.factor(df$ratio))
  upd_thr <- levels(as.factor(df$at))
  bud_bon <- levels(as.factor(df$bb))
  
  # number of samples for bootstrap 
  nbs <- 1000
  
  if (omit.extinction == "TRUE") { 
    
    # a loop to calculate extinction freq
    sub <- subset(df, at == 0 & bb == 0 & ratio == bgt_rat[1])
    
    # NULL tab
    ext_freq <- NULL
    sd_ci <- boot_sd_ci(sub$extinct, itr = nbs)
    res <- c(sum(sub$extinct)/dim(sub)[1], sd_ci[1], sd_ci[2],sd_ci[3])
    ext_freq <- rbind(ext_freq, res)
    
    # loop over BRs just for control strat
    for (i in 2:length(bgt_rat)) {
      sub <- subset(df, at == upd_thr[1] & bb == bud_bon[1] & ratio == bgt_rat[i])
      sd_ci <- boot_sd_ci(sub$extinct, itr = nbs)
      res <- c(sum(sub$extinct)/dim(sub)[1], sd_ci[1], sd_ci[2],sd_ci[3])
      ext_freq <- rbind(ext_freq, res)
    }
    
    # loop over the rest
    for (i in 1:length(bgt_rat)) {
      for (j in 2:length(upd_thr)) {
        for (k in 1:length(bud_bon)) {
          sub <- subset(df, at == upd_thr[j] & bb == bud_bon[k] & ratio == bgt_rat[i])
          sd_ci <- boot_sd_ci(sub$extinct, itr = nbs)
          res <- c(sum(sub$extinct)/dim(sub)[1], sd_ci[1], sd_ci[2],sd_ci[3])
          ext_freq <- rbind(ext_freq, res)
        }
      }
    }
    
    print("Ommiting replicates that resulted in Resource extinction")
    df <- subset(df, extinct == 0)
    
    # levels after ommission
    bgt_rat <- levels(as.factor(df$ratio))
    upd_thr <- levels(as.factor(df$at))
    bud_bon <- levels(as.factor(df$bb))
    
    # initiate a count for later
    zz <- 1
    
  } # end if loop on omit.extinction 
  
  # create empty tab
  res_tab <- NULL
  
  ## Special subset for UT = 0, for which there was only BB = 0
  # initiate a string
  res_str <- NULL
  
  # subset
  sub <- subset(df, at == upd_thr[1] & bb == bud_bon[1] & ratio == bgt_rat[1])
  
  # number of replicates
  nbrep <- dim(sub)[1]
  
  # extinction frequency
  if (omit.extinction == TRUE) {
    ext <- ext_freq[zz,]
  } else {
    sd_ci <- boot_sd_ci(sub$extinct, itr = nbs)
    ext <- c(sum(sub$extinct)/dim(sub)[1], sd_ci[1], sd_ci[2], sd_ci[3])
  }
  
  # start filling the string in
  res_str <- c(nbrep, mean(sub$memory), mean(sub$budget), bgt_rat[1], upd_thr[1], bud_bon[1], ext)
  
  ## mean, sd and 95ci for each proxy
  
  # Actual Resource population deviation from Manager's target
  sd_ci <- boot_sd_ci(sub$act_dev, itr = nbs)
  res_str <- c(res_str, mean(sub$act_dev), sd_ci[1], sd_ci[2], sd_ci[3])
  
  # Absolute value of the actual Resource population deviation from Manager's target
  sd_ci <- boot_sd_ci(sub$abs_act_dev, itr = nbs)
  res_str <- c(res_str, mean(sub$abs_act_dev), sd_ci[1], sd_ci[2], sd_ci[3])
  
  # Users' total final yield
  sd_ci <- boot_sd_ci(sub$fin_yield/40000, itr = nbs)
  res_str <- c(res_str, mean(sub$fin_yield/40000), sd_ci[1], sd_ci[2], sd_ci[3])
  
  # Difference between the highest and the lowest yield
  sd_ci <- boot_sd_ci(sub$max_diff_yield, itr = nbs)
  res_str <- c(res_str, mean(sub$max_diff_yield), sd_ci[1], sd_ci[2], sd_ci[3])
  
  # Percentage of time steps of non-updating
  sd_ci <- boot_sd_ci(sub$inac_ts, itr = nbs)
  res_str <- c(res_str, mean(sub$inac_ts), sd_ci[1], sd_ci[2], sd_ci[3])
  
  # Sum of absolute deviation from target
  sd_ci <- boot_sd_ci(sub$SumAbsDev/sub$final_ts, itr = nbs)
  res_str <- c(res_str, mean(sub$SumAbsDev/sub$final_ts), sd_ci[1], sd_ci[2], sd_ci[3])
  
  # binding the string to the tab
  res_tab <- rbind(res_tab, as.numeric(res_str))
  
  # loop over BRs just for control strat
  for (i in 2:length(bgt_rat)) {
    
    # increment tracker
    if (omit.extinction == TRUE) {
      zz <- zz + 1
    }
    
    # initiate a string
    res_str <- NULL
    
    # subset
    sub <- subset(df, ratio == bgt_rat[i] & at == upd_thr[1] & bb == bud_bon[1])
    
    # number of replicates
    nbrep <- dim(sub)[1]
    
    # extinction frequency
    if (omit.extinction == TRUE) {
      ext <- ext_freq[zz,]
    } else {
      sd_ci <- boot_sd_ci(sub$extinct, itr = nbs)
      ext <- c(sum(sub$extinct)/dim(sub)[1], sd_ci[1], sd_ci[2], sd_ci[3])
    }
    
    # start filling the string in
    res_str <- c(nbrep,mean(sub$memory), mean(sub$budget), bgt_rat[i], upd_thr[1], bud_bon[1], ext)
    
    # avoid problems if there is only one replicate
    if (nbrep >= 2) {
      
      ## mean, sd and 95ci for each proxy
      
      # Actual Resource population deviation from Manager's target
      sd_ci <- boot_sd_ci(sub$act_dev, itr = nbs)
      res_str <- c(res_str, mean(sub$act_dev), sd_ci[1], sd_ci[2], sd_ci[3])
      
      # Absolute value of the Actual Resource population deviation from Manager's target
      sd_ci <- boot_sd_ci(sub$abs_act_dev, itr = nbs)
      res_str <- c(res_str, mean(sub$abs_act_dev), sd_ci[1], sd_ci[2], sd_ci[3])
      
      # Users' total final yield
      sd_ci <- boot_sd_ci(sub$fin_yield/40000, itr = nbs)
      res_str <- c(res_str, mean(sub$fin_yield/40000), sd_ci[1], sd_ci[2], sd_ci[3])
      
      # Difference between the highest and the lowest yield
      sd_ci <- boot_sd_ci(sub$max_diff_yield, itr = nbs)
      res_str <- c(res_str, mean(sub$max_diff_yield), sd_ci[1], sd_ci[2], sd_ci[3])
      
      # Percentage of time steps of non-updating
      sd_ci <- boot_sd_ci(sub$inac_ts, itr = nbs)
      res_str <- c(res_str, mean(sub$inac_ts), sd_ci[1], sd_ci[2], sd_ci[3])
      
      # Sum of absolute deviation from target
      sd_ci <- boot_sd_ci(sub$SumAbsDev/sub$final_ts, itr = nbs)
      res_str <- c(res_str, mean(sub$SumAbsDev/sub$final_ts), sd_ci[1], sd_ci[2], sd_ci[3])
      
    } else {
      
      print(paste("parameter set with budget ratio = ", as.numeric(bgt_rat)[i]*100, "%, UT = ", as.numeric(upd_thr[j])*100, "% and BB = ", as.numeric(bud_bon[k])*100, "% has less than 2 replicates"))
      
      # Actual Resource population deviation from Manager's target
      res_str <- c(res_str, sub$act_dev, 0, sub$act_dev, sub$act_dev)
      
      # Absolute value of the Actual Resource population deviation from Manager's target
      res_str <- c(res_str, sub$abs_act_dev, 0, sub$abs_act_dev, sub$abs_act_dev)
      
      # Users' total final yield
      res_str <- c(res_str, sub$fin_yield/40000, 0, sub$fin_yield/40000, sub$fin_yield/40000)
      
      # Difference between the highest and the lowest yield
      res_str <- c(res_str, sub$max_diff_yield, 0, sub$max_diff_yield, sub$max_diff_yield)
      
      # Percentage of time steps of non-updating
      res_str <- c(res_str, sub$inac_ts, 0, sub$inac_ts, sub$inac_ts)
      
      # Sum of squared deviation from target
      res_str <- c(res_str, sub$SumAbsDev/sub$final_ts, 0, sub$SumAbsDev, sub$SumAbsDev)
    } # end else loop on nbrep
    
    # binding the string to the tab
    res_tab <- rbind(res_tab, as.numeric(res_str))
    
  }
  
  ## loop over the other OTI parameters
  # for each ratio value
  for (i in 1:length(bgt_rat)) {
    
    # for each positive at value
    for (j in 2:length(upd_thr)) {
      
      # for each bb value
      for (k in 1:length(bud_bon)) {
        
        # increment tracker
        if (omit.extinction == TRUE) {
          zz <- zz + 1
        }
        
        # initiate a string
        res_str <- NULL
        
        # subset
        sub <- subset(df, ratio == bgt_rat[i] & at == upd_thr[j] & bb == bud_bon[k])
        
        # number of replicates
        nbrep <- dim(sub)[1]
        
        # extinction frequency
        if (omit.extinction == TRUE) {
          ext <- ext_freq[zz,]
        } else {
          sd_ci <- boot_sd_ci(sub$extinct, itr = nbs)
          ext <- c(sum(sub$extinct)/dim(sub)[1], sd_ci[1], sd_ci[2], sd_ci[3])
        }
        
        # start filling the string in
        res_str <- c(nbrep,mean(sub$memory), mean(sub$budget), bgt_rat[i], upd_thr[j], bud_bon[k], ext)
        
        # avoid problems if there is only one replicate
        if (nbrep >= 2) {
          
          ## mean, sd and 95ci for each proxy
          
          # Actual Resource population deviation from Manager's target
          sd_ci <- boot_sd_ci(sub$act_dev, itr = nbs)
          res_str <- c(res_str, mean(sub$act_dev), sd_ci[1], sd_ci[2], sd_ci[3])
          
          # Absolute value of the Actual Resource population deviation from Manager's target
          sd_ci <- boot_sd_ci(sub$abs_act_dev, itr = nbs)
          res_str <- c(res_str, mean(sub$abs_act_dev), sd_ci[1], sd_ci[2], sd_ci[3])
          
          # Users' total final yield
          sd_ci <- boot_sd_ci(sub$fin_yield/40000, itr = nbs)
          res_str <- c(res_str, mean(sub$fin_yield/40000), sd_ci[1], sd_ci[2], sd_ci[3])
          
          # Difference between the highest and the lowest yield
          sd_ci <- boot_sd_ci(sub$max_diff_yield, itr = nbs)
          res_str <- c(res_str, mean(sub$max_diff_yield), sd_ci[1], sd_ci[2], sd_ci[3])
          
          # Percentage of time steps of non-updating
          sd_ci <- boot_sd_ci(sub$inac_ts, itr = nbs)
          res_str <- c(res_str, mean(sub$inac_ts), sd_ci[1], sd_ci[2], sd_ci[3])
          
          # Sum of absolute deviation from target
          sd_ci <- boot_sd_ci(sub$SumAbsDev/sub$final_ts, itr = nbs)
          res_str <- c(res_str, mean(sub$SumAbsDev/sub$final_ts), sd_ci[1], sd_ci[2], sd_ci[3])
          
        } else {
          
          print(paste("parameter set with budget ratio = ", as.numeric(bgt_rat)[i]*100, "%, UT = ", as.numeric(upd_thr[j])*100, "% and BB = ", as.numeric(bud_bon[k])*100, "% has less than 2 replicates"))
          
          # Actual Resource population deviation from Manager's target
          res_str <- c(res_str, sub$act_dev, 0, sub$act_dev, sub$act_dev)
          
          # Absolute value of the Actual Resource population deviation from Manager's target
          res_str <- c(res_str, sub$abs_act_dev, 0, sub$abs_act_dev, sub$abs_act_dev)
          
          # Users' total final yield
          res_str <- c(res_str, sub$fin_yield/40000, 0, sub$fin_yield/40000, sub$fin_yield/40000)
          
          # Difference between the highest and the lowest yield
          res_str <- c(res_str, sub$max_diff_yield, 0, sub$max_diff_yield, sub$max_diff_yield)
          
          # Percentage of time steps of non-updating
          res_str <- c(res_str, sub$inac_ts, 0, sub$inac_ts, sub$inac_ts)
          
          # Sum of squared deviation from target
          res_str <- c(res_str, sub$SumAbsDev/sub$final_ts, 0, sub$SumAbsDev, sub$SumAbsDev)
        } # end else loop on nbrep
        
        # binding the string to the tab
        res_tab <- rbind(res_tab, as.numeric(res_str))

      } # end for loop on budget bonus
    } # end for loop on update threshold
  } # end for loop on budget ratio
  
  # Array of column names
  colnames(res_tab) <- c("rep", "memory", "budget", "ratio", "at", "bb", "ext_prob", "ext_prob_sd", "ext_prob_95ci_inf", "ext_prob_95ci_sup", "act_dev", "act_dev_sd", "act_dev_95ci_inf", "act_dev_95ci_sup", "abs_act_dev", "abs_act_dev_sd", "abs_act_dev_95ci_inf", "abs_act_dev_95ci_sup", "fin_yield", "fin_yield_sd", "fin_yield_95ci_inf", "fin_yield_95ci_sup", "max_diff_yield", "max_diff_yield_sd", "max_diff_yield_95ci_inf", "max_diff_yield_95ci_sup", "inac_ts", "inac_ts_sd", "inac_ts_95ci_inf", "inac_ts_95ci_sup", "SumAbsDev", "SumAbsDev_sd", "SumAbsDev_95ci_inf", "SumAbsDev_95ci_sup")
  
  res_tab <- as.data.frame(res_tab)
  
  return(res_tab)
  
} # end function

#### import data ####

path <- "~/Desktop/PhD/GitKraken/gmse_fork_RQ1/Budget-ratio-large-batch/"

setwd(path)

dir.name <- "Mem-res"

brut <- as.data.frame(read.csv(paste(dir.name, "/merged-res/mem-budget-ratio-merged.csv", sep = "")))

stat <- OTI_stats(df = brut, ts = 20, omit.extinction = F) 
woe_stat <- OTI_stats(df = brut, ts = 20, omit.extinction = T)

# in case of problem
stat2 <- stat[-which(is.na(stat$at)),]
stat2 <- stat2[-c(1:6),]

write.csv(stat, file = paste(dir.name, "/merged-res/bgt-ratio-mem-stats.csv", sep=""))
write.csv(woe_stat, file = paste(dir.name,"/bgt-ratio-noMem-woExt-stats.csv", sep=""))

#### change dev in case of extinctions ####

brut <- read.csv("~/Desktop/PhD/GitKraken/gmse_fork_RQ1/Budget-ratio-large-batch/Mem-res/merged-res/mem-budget-ratio-merged.csv")
brut.mod <- brut
brut.mod[which(brut.mod$extinct == 1), 8] <- -1

stat.mod <- OTI_stats(df = brut.mod, ts = 20, omit.extinction = F)

######## Contour figures ########

stat <- read.csv("~/Desktop/PhD/GitKraken/gmse_fork_RQ1/Budget-ratio-large-batch/Mem-res/merged-res/bgt-ratio-mem-stats.csv")

#### Extinction frequency ####

# build results matrix

d <- subset(stat, at == 0.5)[,-1]

d$bb <- d$bb*100
d$ratio <- d$ratio*1000

bubo <- levels(as.factor(d$bb))
bura <- levels(as.factor(d$ratio))

resmat <- matrix(data = d$ext_prob, ncol = length(bura), nrow = length(bubo))

fig <- plot_ly(
  x = bubo, 
  y = bura, 
  # z = matrix(data = d$ext_prob, ncol = length(bubo), nrow = length(bura)), 
  z = t(resmat),
  type = "contour",
  colorscale = list(c(0, 0.5, 1), c('green', 'orange', 'red')),
  autocontour = F,
  # contours = list(showlabels = TRUE),
  contours = list(
    start = 0,
    end = 1,
    size = 0.1,
    showlabels = T
  )
)

xlab <- list(
  title = "Budget bonus (%)"#,
  # titlefont = f
)
ylab <- list(
  title = "Initial budget"#,
  # titlefont = f
)

fig <- fig %>% colorbar(title = "Extinction \n frequency")
fig <- fig %>% layout(xaxis = xlab, yaxis = ylab)
fig

#### Deviation from target ####

# build results matrix

d <- subset(stat, at == 0.5)

d$bb <- d$bb*100

bubo <- levels(as.factor(d$bb))
bura <- levels(as.factor(d$ratio))

resmat <- matrix(data = 100*d$act_dev, nrow = length(bubo), ncol = length(bura))

fig <- plot_ly(
  x = bubo, 
  y = bura, 
  z = t(resmat), 
  type = "contour",
  colorscale = list(c(0, 0.35, 0.5, 0.65, 1), c('red', 'orange', 'green', 'orange', 'red')),
  autocontour = F,
  # contours = list(showlabels = TRUE)#,
  contours = list(
    start = -100,
    end = 100,
    size = 10,
    showlabels = T
  )
)

xlab <- list(
  title = "Budget bonus (%)"#,
  # titlefont = f
)
ylab <- list(
  title = "Manager's budget"#,
  # titlefont = f
)

fig <- fig %>% colorbar(title = "Dev. from \n target (%)")
fig <- fig %>% layout(xaxis = xlab, yaxis = ylab)
fig

#### Deviation from target at last ts ####

# build results matrix

d <- subset(stat.mod, at == 0.5)

d$bb <- d$bb*100

bubo <- levels(as.factor(d$bb))
bura <- levels(as.factor(d$ratio))

resmat <- matrix(data = 100*d$act_dev, nrow = length(bubo), ncol = length(bura))
# resmat <- t(resmat)

fig <- plot_ly(
  x = bubo, 
  y = bura, 
  z = t(resmat), 
  type = "contour",
  colorscale = list(c(0, 0.35, 0.5, 0.65, 1), c('red', 'orange', 'green', 'orange', 'red')),
  autocontour = F,
  # contours = list(showlabels = TRUE)#,
  contours = list(
    start = -100,
    end = 100,
    size = 10,
    showlabels = T
  )
)

xlab <- list(
  title = "Budget bonus (%)"#,
  # titlefont = f
)
ylab <- list(
  title = "Manager's budget"#,
  # titlefont = f
)

fig <- fig %>% colorbar(title = "Dev. from \n target (%)")
fig <- fig %>% layout(xaxis = xlab, yaxis = ylab)
fig

#### users yield ####
# build results matrix

d <- subset(stat, at == 0.5)[,-1]

d$bb <- d$bb*100
d$ratio <- d$ratio*100

bubo <- levels(as.factor(d$bb))
bura <- levels(as.factor(d$ratio))

resmat <- matrix(data = 100*d$fin_yield, ncol = length(bura), nrow = length(bubo))

fig <- plot_ly(
  x = bubo, 
  y = bura, 
  z = t(resmat), 
  type = "contour",
  colorscale = list(c(0, 0.5, 1), c('red', 'orange', 'green')),
  autocontour = F,
  # contours = list(showlabels = TRUE)#,
  contours = list(
    start = 50,
    end = 100,
    size = 5,
    showlabels = T
  )
)

xlab <- list(
  title = "Budget bonus (%)"#,
  # titlefont = f
)
ylab <- list(
  title = "Manager's budget"#,
  # titlefont = f
)

fig <- fig %>% colorbar(title = "final\nyield (%)")
fig <- fig %>% layout(xaxis = xlab, yaxis = ylab)
fig

#### Yield inequity ####

# build results matrix

d <- subset(stat, at == 0.5)

d$bb <- d$bb*100

bubo <- levels(as.factor(d$bb))
bura <- levels(as.factor(d$ratio))

resmat <- matrix(data = 100*d$max_diff_yield, ncol = length(bura), nrow = length(bubo))

fig <- plot_ly(
  x = bubo, 
  y = bura, 
  # z = matrix(data = d$ext_prob, ncol = length(bubo), nrow = length(bura)), 
  z = t(resmat),
  type = "contour",
  colorscale = list(c(0, 0.5, 1), c('green', 'orange', 'red')),
  autocontour = F,
  # contours = list(showlabels = TRUE),
  contours = list(
    start = 0,
    end = 20,
    size = 2,
    showlabels = T
  )
)

xlab <- list(
  title = "Budget bonus (%)"#,
  # titlefont = f
)
ylab <- list(
  title = "Manager's budget (%)"#,
  # titlefont = f
)

fig <- fig %>% colorbar(title = "Equity\n(%)")
fig <- fig %>% layout(xaxis = xlab, yaxis = ylab)
fig

#### time waiting ####

# build results matrix

d <- subset(stat, at == 0.5)

d$bb <- d$bb*100

bubo <- levels(as.factor(d$bb))
bura <- levels(as.factor(d$ratio))

resmat <- matrix(data = 100*d$inac_ts, ncol = length(bura), nrow = length(bubo))

fig <- plot_ly(
  x = bubo, 
  y = bura, 
  # z = matrix(data = d$ext_prob, ncol = length(bubo), nrow = length(bura)), 
  z = t(resmat),
  type = "contour",
  #colorscale = list(c(0, 0.5, 1), c('green', 'orange', 'red')),
  autocontour = F,
  # contours = list(showlabels = TRUE),
  contours = list(
    start = 0,
    end = 50,
    size = 5,
    showlabels = T
  )
)

xlab <- list(
  title = "Budget bonus (%)"#,
  # titlefont = f
)
ylab <- list(
  title = "Manager's budget (%)"#,
  # titlefont = f
)

fig <- fig %>% colorbar(title = "Waiting\n(%)")
fig <- fig %>% layout(xaxis = xlab, yaxis = ylab)
fig

#### effect of BB:ratio on extinction frequency according to UT ####

{d <- subset(stat, at == 0.5)
  
  d$at <- d$at*100
  d$bb <- d$bb*100
  
  # get max, min and average of each UT
  upth <- levels(as.factor(d$at))
  bubo <- levels(as.factor(d$bb))
  bura <- levels(as.factor(d$ratio))
  
  sub <- as.data.frame(subset(d, ratio == bura[1]))
  ext <- sub$ext_prob
  sd <- sub$ext_prob_sd
  ci_inf <- sub$ext_prob_95ci_inf
  ci_sup <- sub$ext_prob_95ci_sup
  
  for (i in 2:length(bura)) {
    sub <- as.data.frame(subset(d, ratio == bura[i]))
    ext <- rbind(ext, sub$ext_prob)
    sd <- rbind(sd, sub$ext_prob_sd)
    ci_inf <- rbind(ci_inf, sub$ext_prob_95ci_inf)
    ci_sup <- rbind(ci_sup, sub$ext_prob_95ci_sup)
  }
  
  # Without management?
  no.mgmt <- read.csv("~/Desktop/PhD/GitKraken/gmse_fork_RQ1/batch6-perfObs/manager_budget_is_1.csv", sep = "\t", header = FALSE)
  no.mgmt.var <- sum(no.mgmt[,5])/dim(no.mgmt)[1]
  
  # Without humans?
  no.hum <- read.csv("~/Desktop/PhD/GitKraken/gmse_fork_RQ1/batch6-perfObs/user_and_manager_budget_is_1.csv", sep = "\t", header = FALSE)
  no.hum.var <- sum(no.hum[,5])/dim(no.hum)[1]
}

# plot and export in pdf
{pdf(file = "mem-UT50-bgtRatio-extfreq.pdf", width = par('din')[1], height = par('din')[2])

  {# # enlarge margins
  # par(mar = c(5, 5, 1, 1))
  
  # points cex
  pts <- 0.5
  # pts <- 1
    
  # trunk ext if necessary
  ext <- ext[-(1:2),]
  ci_inf <- ci_inf[-(1:2),]
  ci_sup <- ci_sup[-(1:2),]
  
  xadj <- seq(-2,2,length.out=dim(ext)[1])
  colo <- rainbow(dim(ext)[1])
  
  # plot base
  plot(1, type = "n",
       ylim = c(0, 1),
       xlim = c(0, 105),
       ylab = "Extinction frequency", #
       xlab = "Budget Bonus (%)", #cex.lab = 1.5, cex.axis = 1.5, cex = 1.5,
       cex.lab = pts + 0.2, cex.axis = pts + 0.2)
  
  # # Control band
  # polygon(c(xtendrange,rev(xtendrange)),c(rep(ci_sup[1,1], length(xtendrange)),rev(rep(ci_inf[1,1], length(xtendrange)))),col="lightgrey", border = "grey") 
  
  k <- 1
  for (i in 1:dim(ext)[1]) {
    arrows(x0 = as.numeric(bubo) + xadj[k], y0 = ci_inf[i,], x1 = as.numeric(bubo) + xadj[k], y1 = ci_sup[i,], length=0.02, angle=90, code=3, col = colo[i])
    points(x = as.numeric(bubo) + xadj[k], y = ext[i,], type = "b", cex = pts, lwd = pts-0.2, col = colo[i], pch = 20);
    k <- k+1
  }
  
  # legend
  legend( # 110, 0.5,             # Location of legend 
         "topright",
         # xpd = TRUE,                          # Allow drawing outside plot area
         ncol = 2,
         # xjust = 0,                           # Left justify legend box on x
         # yjust = 0.5,                          # Center legend box on y
         legend = bura[3:6],
           # c(paste("BR = ", bura[1]),
           #          paste("BR = ", bura[2]),
           #          paste("BR = ", bura[3]),
           #          paste("BR = ", bura[4]),
           #          paste("BR = ", bura[5]),
           #          paste("BR = ", bura[6]),
           #          paste("BR = ", bura[7]),
           #          paste("BR = ", bura[8]),
           #          paste("BR = ", bura[9])),
         col = colo,
           # c(colo[1],
           #       colo[2],
           #       colo[3],
           #       colo[4],
           #       colo[5],
           #       colo[6],
           #       colo[7],
           #       colo[8],
           #       colo[9]),          
         pch = rep(16, dim(ext)[1]),
           # c(23,
           #       17,
           #       15,
           #       19,
           #       21,
           #       20),                      # Legend Element Styles          
         lty = rep(1, dim(ext)[1]),
           # c(3,
           #       2,
           #       1,
           #       1,
           #       2,
           #       1),       
         cex = pts-0.2,
         # cex = 0.6,
         title = "Ratios") #,                  # Legend Title
         # title.col = gray(.2) ,                # Legend title color
         # box.lty = 1,                         # Legend box line type
         # box.lwd = 1)                         # Legend box line width
  }
dev.off()
}

#### effect of BB:ratio on deviation from target according to UT ####

{d <- subset(stat, at == 0.3)

d$at <- d$at*100
d$bb <- d$bb*100
d$act_dev <- d$act_dev*100
d$act_dev_sd <- d$act_dev_sd*100
d$act_dev_95ci_inf <- d$act_dev_95ci_inf*100
d$act_dev_95ci_sup <- d$act_dev_95ci_sup*100

# get max, min and average of each UT
upth <- levels(as.factor(d$at))
bubo <- levels(as.factor(d$bb))
bura <- levels(as.factor(d$ratio))

sub <- as.data.frame(subset(d, ratio == bura[1]))
ext <- sub$act_dev
sd <- sub$act_dev_sd
ci_inf <- sub$act_dev_95ci_inf
ci_sup <- sub$act_dev_95ci_sup

for (i in 2:length(bura)) {
  sub <- as.data.frame(subset(d, ratio == bura[i]))
  ext <- rbind(ext, sub$act_dev)
  sd <- rbind(sd, sub$act_dev_sd)
  ci_inf <- rbind(ci_inf, sub$act_dev_95ci_inf)
  ci_sup <- rbind(ci_sup, sub$act_dev_95ci_sup)
}

# # Without management?
# no.mgmt <- read.csv("~/Desktop/PhD/GitKraken/gmse_fork_RQ1/batch6-perfObs/manager_budget_is_1.csv", sep = "\t", header = FALSE)
# no.mgmt.var <- sum(no.mgmt[,5])/dim(no.mgmt)[1]
# 
# # Without humans?
# no.hum <- read.csv("~/Desktop/PhD/GitKraken/gmse_fork_RQ1/batch6-perfObs/user_and_manager_budget_is_1.csv", sep = "\t", header = FALSE)
# no.hum.var <- sum(no.hum[,5])/dim(no.hum)[1]

# # plotting convenience
# xadj1 <- as.numeric(upth[-1]) - 1;
# xadj2 <- as.numeric(upth[-1]) + 1;
# xtendrange <- seq(-10,110,1)
}

#### Details of BB effect and waiting strategy on dev according to ratio ####

{d <- subset(stat, at == 0 & ratio == 0.8 | at == 0.3 & ratio == 0.8 | at == 0.5 & ratio == 0.8)

d$at <- d$at*100
d$bb <- d$bb*100

bubo <- levels(as.factor(d$bb))

sub <- subset(d, at != 50)
ext <- sub$act_dev*100
sd <- sub$act_dev_sd*100
ci_inf <- sub$act_dev_95ci_inf*100
ci_sup <- sub$act_dev_95ci_sup*100

sub <- subset(d, at != 30)
ext <- cbind(ext, sub$act_dev*100)
sd <- cbind(sd, sub$act_dev_sd*100)
ci_inf <- cbind(ci_inf, sub$act_dev_95ci_inf*100)
ci_sup <- cbind(ci_sup, sub$act_dev_95ci_sup*100)
}

# plot and export in pdf
{pdf(file = "mem-BR08-bgtRatio-dev.pdf", width = par('din')[1], height = par('din')[2])
  
  { # points cex
    pts <- 0.5
    
    # plot base
    plot(1, type = "n",
         ylim = c(-100, 20),
         xlim = c(0, 100),
         ylab = "Deviation from target (%)", #
         xlab = "Budget Bonus (%)", #cex.lab = 1.5, cex.axis = 1.5, cex = 1.5,
         cex.lab = pts + 0.2, cex.axis = pts + 0.2)
    
    # Control band
    xtendrange <- seq(-10,110,1)
    
    polygon(c(xtendrange,rev(xtendrange)),c(rep(ci_sup[1,1], length(xtendrange)),rev(rep(ci_inf[1,1], length(xtendrange)))),col="lightgrey", border = "grey")
    abline(h = ext[1,1], lwd = pts)
    
    # best possible
    abline(h = 0, lty = 2, lwd = pts, col = "darkgreen")
    
    # UT30
    arrows(x0 = as.numeric(bubo)-1, y0 = ci_inf[-1,1], x1 = as.numeric(bubo)-1, y1 = ci_sup[-1,1], length=0.02, angle=90, code=3, col = 'blue')
    points(x = as.numeric(bubo)-1, y = ext[-1,1], type = "b", cex = pts, lwd = pts, col = 'blue', pch = 20);
    
    # UT50
    arrows(x0 = as.numeric(bubo)+1, y0 = ci_inf[-1,2], x1 = as.numeric(bubo)+1, y1 = ci_sup[-1,2], length=0.02, angle=90, code=3, col = 'violet')
    points(x = as.numeric(bubo)+1, y = ext[-1,2], type = "b", cex = pts, lwd = pts, col = 'violet', pch = 4);
    
    # legend
    legend( # 110, 0.5,             # Location of legend 
      "topright",
      # xpd = TRUE,                          # Allow drawing outside plot area
      # ncol = 2,
      # xjust = 0,                           # Left justify legend box on x
      # yjust = 0.5,                          # Center legend box on y
      legend = c("Control",
                 "UT 30%",
                 "UT 50%"),
      col = c("black",
              "blue",
              "violet"),        
      pch = c(NA_integer_,
              20,
              4),                    # Legend Element Styles          
      lty = c(1, 
              1,
              1),     
      cex = pts-0.2,
      # cex = 0.6,
      title = "Strategies - 0.8 ratio") #,                  # Legend Title
    # title.col = gray(.2) ,                # Legend title color
    # box.lty = 1,                         # Legend box line type
    # box.lwd = 1)                         # Legend box line width
}
  dev.off()
}



#### Details of BB effect and waiting strategy on ext freq according to ratio ####

stat <- read.csv(file = "~/Desktop/PhD/GitKraken/gmse_fork_RQ1/Budget-ratio-large-batch/Mem-res/merged-res/bgt-ratio-mem-stats.csv", header = T, stringsAsFactors = T)[,-1]

{d <- subset(stat, at == 0 & ratio == 0.8 | at == 0.5 & ratio == 0.8) # | at == 0.5 & ratio == 0.7

d$at <- d$at*100
d$bb <- d$bb*100

bubo <- levels(as.factor(d$bb))

sub <- subset(d, at != 50)
ext <- sub$ext_prob
sd <- sub$ext_prob_sd
ci_inf <- sub$ext_prob_95ci_inf
ci_sup <- sub$ext_prob_95ci_sup

# sub <- subset(d, at != 30)
# ext <- cbind(ext, sub$ext_prob)
# sd <- cbind(sd, sub$ext_prob_sd)
# ci_inf <- cbind(ci_inf, sub$ext_prob_95ci_inf)
# ci_sup <- cbind(ci_sup, sub$ext_prob_95ci_sup)

sub <- subset(d, at != 0)
ext <- c(ext, sub$ext_prob)
sd <- c(sd, sub$ext_prob_sd)
ci_inf <- c(ci_inf, sub$ext_prob_95ci_inf)
ci_sup <- c(ci_sup, sub$ext_prob_95ci_sup)
}

# plot and export in pdf
{pdf(file = "mem-BR08-UT50-bgtRatio-extfreq.pdf", width = par('din')[1], height = par('din')[2])
  
  { # points cex
    pts <- 0.5
    
    # plot base
    plot(1, type = "n",
         ylim = c(0, 1),
         xlim = c(0, 100),
         ylab = "Extinction frequency", #
         xlab = "Budget Bonus (%)", #cex.lab = 1.5, cex.axis = 1.5, cex = 1.5,
         cex.lab = pts + 0.2, cex.axis = pts + 0.2)
    
    # Control band
    xtendrange <- seq(-10,110,1)
    
    # polygon(c(xtendrange,rev(xtendrange)),c(rep(ci_sup[1,1], length(xtendrange)),rev(rep(ci_inf[1,1], length(xtendrange)))),col="lightgrey", border = "grey")
    # abline(h = ext[1,1], lwd = pts)
    
    polygon(c(xtendrange,rev(xtendrange)),c(rep(ci_sup[1], length(xtendrange)),rev(rep(ci_inf[1], length(xtendrange)))),col="lightgrey", border = "grey")
    abline(h = ext[1], lwd = pts)
    
    # best possible
    abline(h = 0, lty = 2, lwd = pts, col = "darkgreen")
    
    # # UT30
    # arrows(x0 = as.numeric(bubo)-1, y0 = ci_inf[-1,1], x1 = as.numeric(bubo)-1, y1 = ci_sup[-1,1], length=0.02, angle=90, code=3, col = 'blue')
    # points(x = as.numeric(bubo)-1, y = ext[-1,1], type = "b", cex = pts, lwd = pts, col = 'blue', pch = 20);
    
    # UT50
    # arrows(x0 = as.numeric(bubo)+1, y0 = ci_inf[-1,2], x1 = as.numeric(bubo)+1, y1 = ci_sup[-1,2], length=0.02, angle=90, code=3, col = 'violet')
    # points(x = as.numeric(bubo)+1, y = ext[-1,2], type = "b", cex = pts, lwd = pts, col = 'violet', pch = 4);
    arrows(x0 = as.numeric(bubo), y0 = ci_inf[-1], x1 = as.numeric(bubo), y1 = ci_sup[-1], length=0.02, angle=90, code=3, col = 'violet')
    points(x = as.numeric(bubo), y = ext[-1], type = "b", cex = pts, lwd = pts, col = 'violet', pch = 4);
    
    # legend
    legend( # 110, 0.5,             # Location of legend 
      "topright",
      # xpd = TRUE,                          # Allow drawing outside plot area
      # ncol = 2,
      # xjust = 0,                           # Left justify legend box on x
      # yjust = 0.5,                          # Center legend box on y
      legend = c("Control",
                 # "UT 30%",
                 "UT 50%"),
      col = c("black",
              # "blue",
              "violet"),        
      pch = c(NA_integer_,
              # 20,
              4),                    # Legend Element Styles          
      lty = c(1, 
              # 1,
              1),     
      cex = pts-0.2,
      # cex = 0.6,
      title = "Strategies") #,          - 0.7 ratio         # Legend Title
    # title.col = gray(.2) ,                # Legend title color
    # box.lty = 1,                         # Legend box line type
    # box.lwd = 1)                         # Legend box line width
}
  dev.off()
}

# plot and export in png
setwd("~/Desktop/PhD/GitKraken/gmse_fork_RQ1/Budget-ratio-large-batch/")
{png(file = "extfreq-TRJ-PT50-BR800-bigLabs.png", width = 1000, height = 1000)
  
  { # enlarge margins
    par(mar = c(5, 6, 1, 1)+0.1)
    
    # sizes
    pts <- 3
    lb <- 3
    ax <- 2
    
    # plot base
    plot(1, type = "n",
         ylim = c(0, 1),
         xlim = c(0, 100),
         ylab = "Extinction frequency", #
         xlab = "Budget Bonus (%)", #cex.lab = 1.5, cex.axis = 1.5, cex = 1.5,
         cex.lab = lb, cex.axis = ax)
    
    # Control band
    xtendrange <- seq(-10,110,1)
    
    # polygon(c(xtendrange,rev(xtendrange)),c(rep(ci_sup[1,1], length(xtendrange)),rev(rep(ci_inf[1,1], length(xtendrange)))),col="lightgrey", border = "grey")
    # abline(h = ext[1,1], lwd = pts)
    
    polygon(c(xtendrange,rev(xtendrange)),c(rep(ci_sup[1], length(xtendrange)),rev(rep(ci_inf[1], length(xtendrange)))),col="lightgrey", border = "grey")
    abline(h = ext[1], lwd = pts)
    
    # best possible
    abline(h = 0, lty = 2, lwd = pts, col = "darkgreen")
    
    # # UT30
    # arrows(x0 = as.numeric(bubo)-1, y0 = ci_inf[-1,1], x1 = as.numeric(bubo)-1, y1 = ci_sup[-1,1], length=0.02, angle=90, code=3, col = 'blue')
    # points(x = as.numeric(bubo)-1, y = ext[-1,1], type = "b", cex = pts, lwd = pts, col = 'blue', pch = 20);
    
    # UT50
    # arrows(x0 = as.numeric(bubo)+1, y0 = ci_inf[-1,2], x1 = as.numeric(bubo)+1, y1 = ci_sup[-1,2], length=0.02, angle=90, code=3, col = 'violet')
    # points(x = as.numeric(bubo)+1, y = ext[-1,2], type = "b", cex = pts, lwd = pts, col = 'violet', pch = 4);
    arrows(x0 = as.numeric(bubo), y0 = ci_inf[-1], x1 = as.numeric(bubo), y1 = ci_sup[-1], length=0.1, angle=90, code=3, col = 'darkmagenta', lwd = pts-1)
    points(x = as.numeric(bubo), y = ext[-1], type = "b", cex = pts, lwd = pts-1, col = 'darkmagenta', pch = 17, lty = 3);
    
    # legend
    legend( # 110, 0.5,             # Location of legend 
      "topright",
      # xpd = TRUE,                          # Allow drawing outside plot area
      # ncol = 2,
      # xjust = 0,                           # Left justify legend box on x
      # yjust = 0.5,                          # Center legend box on y
      legend = c("Control",
                 # "UT 30%",
                 "TRJ - PT 50%"),
      col = c("black",
              # "blue",
              "darkmagenta"),        
      pch = c(NA_integer_,
              # 20,
              17),                    # Legend Element Styles          
      lty = c(1, 
              # 1,
              3),     
      cex = pts-1,
      # cex = 0.6,
      title = "Strategies") #,          - 0.7 ratio         # Legend Title
    # title.col = gray(.2) ,                # Legend title color
    # box.lty = 1,                         # Legend box line type
    # box.lwd = 1)                         # Legend box line width
}
  dev.off()
}

#### Details of BB effect and waiting strategy on yield according to ratio ####

{d <- subset(stat, at == 0 & ratio == 0.8 | at == 0.5 & ratio == 0.8) # | at == 0.5 & ratio == 0.7

d$at <- d$at*100
d$bb <- d$bb*100
d$fin_yield <-                   d$fin_yield*100
d$fin_yield_sd <-             d$fin_yield_sd*100
d$fin_yield_95ci_inf <- d$fin_yield_95ci_inf*100
d$fin_yield_95ci_sup <- d$fin_yield_95ci_sup*100

bubo <- levels(as.factor(d$bb))

sub <- subset(d, at != 50)
ext <-    sub$fin_yield
sd <-     sub$fin_yield_sd
ci_inf <- sub$fin_yield_95ci_inf
ci_sup <- sub$fin_yield_95ci_sup

# sub <- subset(d, at != 30)
# ext <- cbind(ext, sub$ext_prob)
# sd <- cbind(sd, sub$ext_prob_sd)
# ci_inf <- cbind(ci_inf, sub$ext_prob_95ci_inf)
# ci_sup <- cbind(ci_sup, sub$ext_prob_95ci_sup)

sub <- subset(d, at != 0)
ext <- c(ext, sub$fin_yield)
sd <- c(sd, sub$fin_yield_sd)
ci_inf <- c(ci_inf, sub$fin_yield_95ci_inf)
ci_sup <- c(ci_sup, sub$fin_yield_95ci_sup)
}

# plot and export in pdf
{pdf(file = "BDG-PT50-BM800-finYie.pdf", width = par('din')[1], height = par('din')[2])
  
  { # points cex
    # pts <- 0.5
    pts <- 1
    
    # plot base
    plot(1, type = "n",
         ylim = c(80, 100),
         xlim = c(0, 100),
         ylab = "Farmers final yield (%)", #
         xlab = "Budget Bonus (%)", #cex.lab = 1.5, cex.axis = 1.5, cex = 1.5,
         cex.lab = pts + 0.2, cex.axis = pts + 0.2)
    
    # Control band
    xtendrange <- seq(-10,110,1)
    
    # polygon(c(xtendrange,rev(xtendrange)),c(rep(ci_sup[1,1], length(xtendrange)),rev(rep(ci_inf[1,1], length(xtendrange)))),col="lightgrey", border = "grey")
    # abline(h = ext[1,1], lwd = pts)
    
    polygon(c(xtendrange,rev(xtendrange)),c(rep(ci_sup[1], length(xtendrange)),rev(rep(ci_inf[1], length(xtendrange)))),col="lightgrey", border = "grey")
    abline(h = ext[1], lwd = pts)
    
    # best possible
    abline(h = 40, lty = 2, lwd = pts, col = "darkgreen")
    
    # # UT30
    # arrows(x0 = as.numeric(bubo)-1, y0 = ci_inf[-1,1], x1 = as.numeric(bubo)-1, y1 = ci_sup[-1,1], length=0.02, angle=90, code=3, col = 'blue')
    # points(x = as.numeric(bubo)-1, y = ext[-1,1], type = "b", cex = pts, lwd = pts, col = 'blue', pch = 20);
    
    # UT50
    # arrows(x0 = as.numeric(bubo)+1, y0 = ci_inf[-1,2], x1 = as.numeric(bubo)+1, y1 = ci_sup[-1,2], length=0.02, angle=90, code=3, col = 'violet')
    # points(x = as.numeric(bubo)+1, y = ext[-1,2], type = "b", cex = pts, lwd = pts, col = 'violet', pch = 4);
    arrows(x0 = as.numeric(bubo), y0 = ci_inf[-1], x1 = as.numeric(bubo), y1 = ci_sup[-1], length=0.02, angle=90, code=3, col = 'violet')
    points(x = as.numeric(bubo), y = ext[-1], type = "b", cex = pts, lwd = pts, col = 'violet', pch = 20);
    
    # legend
    legend( # 110, 0.5,             # Location of legend 
      "bottomright",
      # xpd = TRUE,                          # Allow drawing outside plot area
      # ncol = 2,
      # xjust = 0,                           # Left justify legend box on x
      # yjust = 0.5,                          # Center legend box on y
      legend = c("Control",
                 # "UT 30%",
                 "Trajectory"),
      col = c("black",
              # "blue",
              "violet"),        
      pch = c(NA_integer_,
              # 20,
              20),                    # Legend Element Styles          
      lty = c(1, 
              # 1,
              1),     
      # cex = pts-0.2,
      cex = 0.6,
      title = "Strategies") #,          - 0.7 ratio         # Legend Title
    # title.col = gray(.2) ,                # Legend title color
    # box.lty = 1,                         # Legend box line type
    # box.lwd = 1)                         # Legend box line width
}
  dev.off()
}

# plot and export in png
setwd("~/Desktop/PhD/GitKraken/gmse_fork_RQ1/Budget-ratio-large-batch/")
{png(file = "finYie-TRJ-PT50-BR800-bigLabs.png", width = 1000, height = 1000)
  
  { # enlarge margins
    par(mar = c(5, 6, 1, 1)+0.1)
    
    # sizes
    pts <- 3
    lb <- 3
    ax <- 2
    
    # plot base
    plot(1, type = "n",
         ylim = c(85, 100),
         xlim = c(0, 100),
         ylab = "Farmers final yield (%)", #
         xlab = "Budget Bonus (%)", #cex.lab = 1.5, cex.axis = 1.5, cex = 1.5,
         cex.lab = lb, cex.axis = ax)

    # Control band
    xtendrange <- seq(-10,110,1)
    
    # polygon(c(xtendrange,rev(xtendrange)),c(rep(ci_sup[1,1], length(xtendrange)),rev(rep(ci_inf[1,1], length(xtendrange)))),col="lightgrey", border = "grey")
    # abline(h = ext[1,1], lwd = pts)
    
    polygon(c(xtendrange,rev(xtendrange)),c(rep(ci_sup[1], length(xtendrange)),rev(rep(ci_inf[1], length(xtendrange)))),col="lightgrey", border = "grey")
    abline(h = ext[1], lwd = pts)
    
    # best possible
    abline(h = 100, lty = pts, lwd = pts, col = "darkgreen")
    
    # # UT30
    # arrows(x0 = as.numeric(bubo)-1, y0 = ci_inf[-1,1], x1 = as.numeric(bubo)-1, y1 = ci_sup[-1,1], length=0.02, angle=90, code=3, col = 'blue')
    # points(x = as.numeric(bubo)-1, y = ext[-1,1], type = "b", cex = pts, lwd = pts, col = 'blue', pch = 20);
    
    # UT50
    # arrows(x0 = as.numeric(bubo)+1, y0 = ci_inf[-1,2], x1 = as.numeric(bubo)+1, y1 = ci_sup[-1,2], length=0.02, angle=90, code=3, col = 'violet')
    # points(x = as.numeric(bubo)+1, y = ext[-1,2], type = "b", cex = pts, lwd = pts, col = 'violet', pch = 4);
    arrows(x0 = as.numeric(bubo), y0 = ci_inf[-1], x1 = as.numeric(bubo), y1 = ci_sup[-1], length=0.1, angle=90, code=3, col = 'darkmagenta', lwd = pts-1)
    points(x = as.numeric(bubo), y = ext[-1], type = "b", cex = pts, lwd = pts-1, col = 'darkmagenta', pch = 17, lty = 3);
    
    # legend
    legend( # 110, 0.5,             # Location of legend 
      "bottomright",
      # xpd = TRUE,                          # Allow drawing outside plot area
      # ncol = 2,
      # xjust = 0,                           # Left justify legend box on x
      # yjust = 0.5,                          # Center legend box on y
      legend = c("Control",
                 # "UT 30%",
                 "TRJ - PT 50%"),
      col = c("black",
              # "blue",
              "darkmagenta"),        
      pch = c(NA_integer_,
              # 20,
              17),                    # Legend Element Styles          
      lty = c(1, 
              # 1,
              3),     
      cex = pts-1,
      # cex = 0.6,
      title = "Strategies") #,          - 0.7 ratio         # Legend Title
    # title.col = gray(.2) ,                # Legend title color
    # box.lty = 1,                         # Legend box line type
    # box.lwd = 1)                         # Legend box line width
    }
  dev.off()
}
    
#### Details of BB effect and waiting strategy on time waiting according to ratio ####

{d <- subset(stat, at == 0 & ratio == 0.8 | at == 0.5 & ratio == 0.8) # | at == 0.5 & ratio == 0.7

d$at <- d$at*100
d$bb <- d$bb*100
d$inac_ts          <- d$inac_ts*100
d$inac_ts_sd       <- d$inac_ts_sd*100
d$inac_ts_95ci_inf <- d$inac_ts_95ci_inf*100
d$inac_ts_95ci_sup <- d$inac_ts_95ci_sup*100

bubo <- levels(as.factor(d$bb))

sub <- subset(d, at != 50)
ext <-    sub$inac_ts
sd <-     sub$inac_ts_sd
ci_inf <- sub$inac_ts_95ci_inf
ci_sup <- sub$inac_ts_95ci_sup

# sub <- subset(d, at != 30)
# ext <- cbind(ext, sub$ext_prob)
# sd <- cbind(sd, sub$ext_prob_sd)
# ci_inf <- cbind(ci_inf, sub$ext_prob_95ci_inf)
# ci_sup <- cbind(ci_sup, sub$ext_prob_95ci_sup)

sub <- subset(d, at != 0)
ext    <- c(ext,    sub$inac_ts)
sd     <- c(sd,     sub$inac_ts_sd)
ci_inf <- c(ci_inf, sub$inac_ts_95ci_inf)
ci_sup <- c(ci_sup, sub$inac_ts_95ci_sup)
}

# plot and export in pdf
{pdf(file = "mem-BR08-UT50-bgtRatio-tw-large.pdf", width = par('din')[1], height = par('din')[2])
  
  { # points cex
    # pts <- 0.5
    pts <- 1
    
    # plot base
    plot(1, type = "n",
         ylim = c(0, 40),
         xlim = c(0, 100),
         ylab = "Average waiting time (%)", #
         xlab = "Budget Bonus (%)", #cex.lab = 1.5, cex.axis = 1.5, cex = 1.5,
         cex.lab = pts + 0.2, cex.axis = pts + 0.2)
    
    # # Control band
    # xtendrange <- seq(-10,110,1)
    
    # polygon(c(xtendrange,rev(xtendrange)),c(rep(ci_sup[1,1], length(xtendrange)),rev(rep(ci_inf[1,1], length(xtendrange)))),col="lightgrey", border = "grey")
    # abline(h = ext[1,1], lwd = pts)
    
    # polygon(c(xtendrange,rev(xtendrange)),c(rep(ci_sup[1], length(xtendrange)),rev(rep(ci_inf[1], length(xtendrange)))),col="lightgrey", border = "grey")
    # abline(h = ext[1], lwd = pts)
    
    # # best possible
    # abline(h = 40, lty = 2, lwd = pts, col = "darkgreen")
    
    # # UT30
    # arrows(x0 = as.numeric(bubo)-1, y0 = ci_inf[-1,1], x1 = as.numeric(bubo)-1, y1 = ci_sup[-1,1], length=0.02, angle=90, code=3, col = 'blue')
    # points(x = as.numeric(bubo)-1, y = ext[-1,1], type = "b", cex = pts, lwd = pts, col = 'blue', pch = 20);
    
    # UT50
    # arrows(x0 = as.numeric(bubo)+1, y0 = ci_inf[-1,2], x1 = as.numeric(bubo)+1, y1 = ci_sup[-1,2], length=0.02, angle=90, code=3, col = 'violet')
    # points(x = as.numeric(bubo)+1, y = ext[-1,2], type = "b", cex = pts, lwd = pts, col = 'violet', pch = 4);
    arrows(x0 = as.numeric(bubo), y0 = ci_inf[-1], x1 = as.numeric(bubo), y1 = ci_sup[-1], length=0.02, angle=90, code=3, col = 'violet')
    points(x = as.numeric(bubo), y = ext[-1], type = "b", cex = pts, lwd = pts, col = 'violet', pch = 20);
    
#     # legend
#     legend( # 110, 0.5,             # Location of legend 
#       "bottomright",
#       # xpd = TRUE,                          # Allow drawing outside plot area
#       # ncol = 2,
#       # xjust = 0,                           # Left justify legend box on x
#       # yjust = 0.5,                          # Center legend box on y
#       legend = c("Control",
#                  # "UT 30%",
#                  "Trajectory"),
#       col = c("black",
#               # "blue",
#               "violet"),        
#       pch = c(NA_integer_,
#               # 20,
#               20),                    # Legend Element Styles          
#       lty = c(1, 
#               # 1,
#               1),     
#       # cex = pts-0.2,
#       cex = 0.6,
#       title = "Strategies") #,          - 0.7 ratio         # Legend Title
#     # title.col = gray(.2) ,                # Legend title color
#     # box.lty = 1,                         # Legend box line type
#     # box.lwd = 1)                         # Legend box line width
  }
  dev.off()
}

#### effect of UT:BB on poplation deviation from target ####

{d <- stat

d$at <- d$at*100
d$bb <- d$bb*100
d$act_dev <- d$act_dev*100
d$act_dev_95ci_inf <- d$act_dev_95ci_inf*100
d$act_dev_95ci_sup <- d$act_dev_95ci_sup*100

# get max, min and average of each UT
upth <- levels(as.factor(d$at))
bubo <- levels(as.factor(d$bb))

sub <- as.data.frame(subset(d, at == upth[1]))
var <- c(sub$act_dev, rep(NA,length(upth)-1))
sd <- c(sub$act_dev_sd, rep(NA,length(upth)-1))
ci_inf <- c(sub$act_dev_95ci_inf, rep(NA,length(upth)-1))
ci_sup <- c(sub$act_dev_95ci_sup, rep(NA,length(upth)-1))

for (i in 2:length(upth)) {
  sub <- as.data.frame(subset(d, at == upth[i]))
  var <- rbind(var, sub$act_dev)
  sd <- rbind(sd, sub$act_dev_sd)
  ci_inf <- rbind(ci_inf, sub$act_dev_95ci_inf)
  ci_sup <- rbind(ci_sup, sub$act_dev_95ci_sup)
}

# Without management?
no.mgmt <- read.csv("~/Desktop/PhD/GitKraken/gmse_fork_RQ1/batch6-perfObs/manager_budget_is_1.csv", sep = "\t", header = FALSE)
no.mgmt.var <- sum(no.mgmt[,6])/dim(no.mgmt)[1]

# Without humans?
no.hum <- read.csv("~/Desktop/PhD/GitKraken/gmse_fork_RQ1/batch6-perfObs/user_and_manager_budget_is_1.csv", sep = "\t", header = FALSE)
no.hum.var <- sum(no.hum[,6])/dim(no.hum)[1]

# plotting convenience
xadj1 <- as.numeric(upth[-1]) - 1;
xadj2 <- as.numeric(upth[-1]) + 1;
xtendrange <- seq(-10,110,1)
}

# plot and export in pdf
{pdf(file = "T1Q1-dev-small.pdf", width = par('din')[1], height = par('din')[2])
  
  {# # enlarge margins
    # par(mar = c(5, 5, 1, 1))
    
    # points cex
    pts <- 0.5
    # pts <- 1
    
    # plot base
    plot(1, type = "n",
         ylim = c(-100, 100),
         xlim = c(0, 100),
         ylab = "Average deviation from target", #
         xlab = "Update threshold (%)", #cex.lab = 1.5, cex.axis = 1.5, cex = 1.5,
         cex.lab = pts + 0.2, cex.axis = pts + 0.2)
    
    # Control band
    polygon(c(xtendrange,rev(xtendrange)),c(rep(ci_sup[1,1], length(xtendrange)),rev(rep(ci_inf[1,1], length(xtendrange)))),col="lightgrey", border = "grey") 
    
    # Put the other budget bonuses in a faint grey to show but de-emphasise
    for (i in 3:(dim(var)[2]-1)) {
      points(x = as.numeric(upth[-1]), y = var[-1,i], type = "b", cex = pts-0.2, lwd = pts-0.2, lty = "solid",
             col = "grey", pch = 20);
    }
    
    # Plot the results
    arrows(0, ci_inf[1,1], 0, ci_sup[1,1], length=0.02, angle=90, code=3, col = "black")
    arrows(xadj1, ci_inf[-1,2], xadj1, ci_sup[-1,2], length=0.02, angle=90, code=3, col = "darkred")
    arrows(xadj2, ci_inf[-1,11], xadj2, ci_sup[-1,11], length=0.02, angle=90, code=3, col = "darkred")
    points(x = xadj1, y = var[-1,2], type = "b", pch = 16, col = "darkred", cex = pts, lwd = pts)
    points(x = xadj2, y = var[-1,11], type = "b", pch = 21, lty = "dashed", col = "darkred", cex = pts, lwd = pts)
    points(x = 0, y = var[1,1], pch = 15, cex = pts, lwd = pts)
    abline(h=var[1,1], lty = 1, col = "black", lwd = pts)
    
    # best possible
    abline(h = 0, col = "darkgreen", lty = 2)
    
    # No management
    points(y = no.mgmt.var, x = 0, pch = 17, col = "black", cex = pts)
    abline(h = no.mgmt.var, lty = 2, col = "black", lwd = pts)
    
    # No humans
    points(y = no.hum.var, x = 0, pch = 23, col = "black", cex = pts)
    abline(h = no.hum.var, lty = 3, col = "black", lwd = pts)
    
    # legend
    legend( # 110, 0.5,             # Location of legend 
      "bottomright",
      # xpd = TRUE,                          # Allow drawing outside plot area
      # ncol = 2,
      # xjust = 0,                           # Left justify legend box on x
      # yjust = 0.5,                          # Center legend box on y
      legend = c("No humans",
                 "No management",
                 "Null deviation",
                 "Control", 
                 "ATI - 0% BB",
                 "ATI - 100% BB", 
                 "other BB values"),
      col = c("black",                 # Legend Element Colors
              "black",
              "darkgreen",
              "black",
              "darkred",
              "darkred",
              "lightgrey"),          
      pch = c(23,
              17,
              NA_integer_,
              15,
              19,
              21,
              20),                      # Legend Element Styles          
      lty = c(3,
              2,
              2,
              1,
              1,
              2,
              1),       
      cex = pts-0.2,
      # cex = 0.6,
      title = "Strategies") #,                  # Legend Title
    # title.col = gray(.2) ,                # Legend title color
    # box.lty = 1,                         # Legend box line type
    # box.lwd = 1)                         # Legend box line width
}
  dev.off()
}

#### effect of UT:BB on the users final yield ####

{d <- stat

d$at <- d$at*100
d$bb <- d$bb*100
d$fin_yield <- d$fin_yield/1000
d$fin_yield_95ci_inf <- d$fin_yield_95ci_inf/1000
d$fin_yield_95ci_sup <- d$fin_yield_95ci_sup/1000

# get max, min and average of each UT
upth <- levels(as.factor(d$at))
bubo <- levels(as.factor(d$bb))

sub <- as.data.frame(subset(d, at == upth[1]))
dev <- c(sub$fin_yield, rep(NA,length(upth)-1))
sd <- c(sub$fin_yield_sd, rep(NA,length(upth)-1))
ci_inf <- c(sub$fin_yield_95ci_inf, rep(NA,length(upth)-1))
ci_sup <- c(sub$fin_yield_95ci_sup, rep(NA,length(upth)-1))

for (i in 2:length(upth)) {
  sub <- as.data.frame(subset(d, at == upth[i]))
  dev <- rbind(dev, sub$fin_yield)
  sd <- rbind(sd, sub$fin_yield_sd)
  ci_inf <- rbind(ci_inf, sub$fin_yield_95ci_inf)
  ci_sup <- rbind(ci_sup, sub$fin_yield_95ci_sup)
}

xadj1 <- as.numeric(upth[-1]) - 3;
xadj2 <- as.numeric(upth[-1]) - 1;
xtendrange <- seq(-10,110,1)

par(mar = c(5, 5, 1, 1))
plot(1, type = "n",
     ylab = "Sum of users' final yield (10^3 a.b.u +/- 95%CI)", cex.lab = 1.5, cex.axis = 1.5, cex = 1.5,
     ylim = c(20, 40), xlim = c(0, 100), lwd = 2)
polygon(c(xtendrange,rev(xtendrange)),c(rep(ci_sup[1,1], length(xtendrange)),rev(rep(ci_inf[1,1], length(xtendrange)))),col="lightgrey", border = "grey") 

# Put the other budget bonuses in a faint grey to show but de-emphasise
for (i in 3:(dim(dev)[2]-1)) {
  points(x = as.numeric(upth[-1])-2, y = dev[-1,i], type = "b", cex = 0.6, lwd = 0.8, lty = "solid",
         col = "grey", pch = 20);
}
abline(h = 40, col = "darkgreen", lty = 2)

arrows(0, ci_inf[1,1], 0, ci_sup[1,1], length=0.05, angle=90, code=3, col = "black")
arrows(xadj1, ci_inf[-1,2], xadj1, ci_sup[-1,2], length=0.05, angle=90, code=3, col = "darkblue")
arrows(xadj2, ci_inf[-1,11], xadj2, ci_sup[-1,11], length=0.05, angle=90, code=3, col = "darkblue")
points(x = xadj1, y = dev[-1,2], xlab = "Update threshold", type = "b", pch = 20,
       cex = 1.5,
       lwd = 2, col = "blue");
points(x = xadj2, y = dev[-1,11], type = "b", cex = 1, lwd = 2, lty = "dashed", col = "blue")
points(x = 0, y = dev[1,1], pch = 15)
abline(h=dev[1,1], lty = 1, col = "black")

# Without management?
no.mgmt <- read.csv("~/Desktop/PhD/GitKraken/gmse_fork_RQ1/batch6-perfObs/manager_budget_is_1.csv", sep = "\t", header = FALSE)

# no.mgmt.extfreq <- length(which(no.mgmt$time_step < 20))/100
no.mgmt.finyie <- (sum(no.mgmt[,8])/dim(no.mgmt)[1])/1000
boot <- boot_sd_ci(no.mgmt[,8])/1000
arrows(0, boot[2], 0, boot[3], length=0.05, angle=90, code=3, col = "black")
points(y = no.mgmt.finyie, x = 0, pch = 17, col = "black")
# abline(h = no.mgmt.finyie, lty = 2, lwd = 1, col = "black")

# Without humans?
no.hum <- read.csv("~/Desktop/PhD/GitKraken/gmse_fork_RQ1/batch6-perfObs/user_and_manager_budget_is_1.csv", sep = "\t", header = FALSE)

no.hum.finyie <- sum(no.hum[,8])/dim(no.hum)[1] * 1/1000
boot <- boot_sd_ci(no.hum[,8]) * 1/1000
arrows(0, boot[2], 0, boot[3], length=0.05, angle=90, code=3, col = "black")
points(y = no.hum.finyie, x = 0, pch = 23, col = "black")
abline(h = no.hum.finyie, lty = 3, lwd = 1, col = "black")

# legend
legend("bottomright",             # Location of legend 
       # ncol = 2,
       # xpd = TRUE,                          # Allow drawing outside plot area
       # xjust = 0,                           # Left justify legend box on x
       # yjust = .5,                          # Center legend box on y
       legend = c("No humans",
                  "No management",
                  "Maximum yield",
                  "Control", 
                  "ATI - 0% BB",
                  "ATI - 100% BB", 
                  "other BB values"),
       col = c("black",                 # Legend Element Colors
               "black",
               "darkgreen",
               "black",
               "blue",
               "blue",
               "lightgrey"),          
       pch = c(23,
               17,
               NA_integer_,
               15,
               20,
               21,
               20),                      # Legend Element Styles          
       lty = c(3,
               NA_integer_,
               2,
               1,
               2,
               1),       
       cex = 0.7,
       title = "Strategies",                  # Legend Title
       title.col = gray(.2),                # Legend title color
       box.lty = 1,                         # Legend box line type
       box.lwd = 1)                         # Legend box line width  
}

#### Plot the evolution of population, costs and actions in FLI and OTI strategy ####

multi.mean.trj <- function(df, upth, bubo, bura, tmax, yaxis) {
  
  # sizes
  pts <- 7
  lb <- 3
  ax <- 2
  
  # subsetting
  dd <- subset(df, UT == upth & BB == bubo & ratio == bura)
  
  # zz <- subset(df, UT == 0 & BB == 0 & ratio == bura)
  
  max.dd <- max(dd[10:dim(dd)[2]])+0.2*max(dd[10:dim(dd)[2]])
  # max.zz <- max(zz[10:dim(dd)[2]])+0.2*max(zz[10:dim(dd)[2]])
  
  # plot
  xrange <- seq(0, tmax, 1)
  plot.new()
  plot(1, type = "n", xlab = "Time", ylab = yaxis, ylim = c(0,max.dd), xlim = c(0,20), xaxt = 'n',
       cex.lab = lb, cex.axis=ax) # , main = paste("UT = ", upth, " BB = ",bubo)
  # plot(1, type = "n", xlab = "time", ylab = yaxis, ylim = c(0,max(max.zz,max.dd)), xlim = c(0,20)) # , main = paste("UT = ", upth, " BB = ",bubo)
  axis(side=1, at=seq(0,length(xrange)-1,1), labels=xrange, cex.axis=ax)
  
  # max with FLI strategy
  if (str_detect(yaxis, "Cost")){
    abline(h = (dd[1,3]*1000/10)+10, lty = 2, pch = 2, col = "black")
    abline(h = 10, lty = 2, pch = 2, col = "black")
  } else if (str_detect(yaxis, "budget")){
    abline(h = dd[1,3]*1000, lty = 2, pch = 2, col = "black")
  } else if (str_detect(yaxis, "action")){
    abline(h = dd[1,2]/10, lty = 2, pch = 2, col = "black")
  } else if (str_detect(yaxis, "Population size")){
    xtendrange <- seq(-1,tmax+1,1)
    trans <- adjustcolor("lightgreen",alpha.f=0.5)
    upperUT <- rep((1+upth)*dd[1,8], length(xtendrange))
    lowerUT <- rep((1-upth)*dd[1,8], length(xtendrange))
    
    polygon(c(xtendrange,rev(xtendrange)),c(upperUT,rev(lowerUT)),col=trans, border = "green")
    abline(h = dd[1,8], col = "darkgreen", lty = 2, lwd = pts-2)
  }
  
  # loop over replicates
  for (i in which(dd$rep %% 10 == 0)) {
    points(x = xrange, y = dd[i,9:dim(dd)[2]], type = 'l', lwd = pts-2, col = 'grey')
  }
  
  moy <- mean(dd[,10])
  infci <- NULL
  supci <- NULL
  for (i in 11:(dim(dd)[2])) {
    moy <- c(moy, mean(dd[,i]))
    infci <- c(infci, boot_sd_ci(dd[,i])[2])
    supci <- c(supci, boot_sd_ci(dd[,i])[3])
  }
  
  # confidence interval
  arrows(xrange[-c(1,2)], infci, xrange[-c(1,2)], supci, length=0.1, angle=90, code=3, col = "black", lwd = pts-4)
  
  # trajectory
  points(x = xrange[-1], y = moy, type = 'l', pch = 1, col = "black", lwd = pts)
}

popul <- read.csv(file = "~/Desktop/PhD/GitKraken/gmse_fork_RQ1/Budget-ratio-large-batch/Mem-res/merged-res/pop-mem-budget-ratio-merged.csv")
library(stringr)
par(mar = c(5,6,1,1)+0.1)
multi.mean.trj(df = popul, upth = 0.5, bubo = 1.0, bura = 0.8, tmax = 20, yaxis = "Population size")
setwd(dir = "~/Desktop/PhD/GitKraken/gmse_fork_RQ1/Budget-ratio-large-batch/")
dev.print(device = png, file = "pop-TRJ-PT50-BB100-BM800-bigLabs.png", width = 1000)

#### plot only one replicate ####
uniq <- function(df, upth, bubo, tmax, yaxis, replicate) {
  # subsetting
  d <- subset(df, UT == upth)
  dd <- subset(df, UT == upth & BB == bubo & rep == replicate) 
  
  # plot
  xrange <- seq(1, tmax, 1)
  plot(1, type = "n", xlab = "time", ylab = yaxis, ylim = c(0,max(d[9:dim(d)[2]])), xlim = c(0,20)) # , main = paste("UT = ", upth, " BB = ",bubo)
  
  # max with FLI strategy
  if (str_detect(yaxis, "Cost")){
    abline(h = (dd[1,1]/10)+10, lty = 2, pch = 2, col = "black")
  } else if (str_detect(yaxis, "budget")){
    abline(h = dd[1,1], lty = 2, pch = 2, col = "black")
  } else if (str_detect(yaxis, "action")){
    abline(h = dd[1,1]/10, lty = 2, pch = 2, col = "black")
  } else if (str_detect(yaxis, "population")){
    xtendrange <- seq(-1,tmax+1,1)
    trans <- adjustcolor("lightgreen",alpha.f=0.5)
    upperUT <- rep((1+upth)*dd[1,6], length(xtendrange))
    lowerUT <- rep((1-upth)*dd[1,6], length(xtendrange))
    
    polygon(c(xtendrange,rev(xtendrange)),c(upperUT,rev(lowerUT)),col=trans, border = "green")
    abline(h = dd[1,6], col = "darkgreen", lty = 2)
  }
  
  points(x = xrange, y = dd[,8:dim(dd)[2]], type = "b", pch = 1, col = "black")
}

#### trajectories panel ####

# find max in 

multi.mean.pop.mod <- function(df, upth, bubo, bura, tmax) {
  
  # subsetting
  dd <- subset(df, UT == upth & BB == bubo & ratio == bura)
  
  # sizes
  pts <- 7
  lb <- 3
  ax <- 2  
  ymax <- 4200
  
  # plot
  xrange <- seq(0, tmax, 1)
  plot(1, type = "n", xlab = "", ylab = "", ylim = c(0,ymax), xlim = c(0,20), xaxt = 'n',
       cex.lab = lb, cex.axis=ax) # , main = paste("UT = ", upth, " BB = ",bubo)
  # plot(1, type = "n", xlab = "time", ylab = yaxis, ylim = c(0,max(max.zz,max.dd)), xlim = c(0,20)) # , main = paste("UT = ", upth, " BB = ",bubo)
  axis(side=1, at=seq(0,length(xrange)-1,1), labels=xrange, cex.axis=ax)
  
  xtendrange <- seq(-1,tmax+1,1)
  trans <- adjustcolor("lightgreen",alpha.f=0.5)
  upperUT <- rep((1+upth)*dd[1,8], length(xtendrange))
  lowerUT <- rep((1-upth)*dd[1,8], length(xtendrange))
  
  polygon(c(xtendrange,rev(xtendrange)),c(upperUT,rev(lowerUT)),col=trans, border = "green")
  abline(h = dd[1,8], col = "darkgreen", lty = 2, lwd = pts-3)
  
  # loop over replicates
  for (i in which(dd$rep %% 10 == 0)) {
    points(x = xrange, y = dd[i,9:dim(dd)[2]], type = 'l', lwd = pts-3, col = 'grey')
  }
  
  moy <- mean(dd[,10])
  infci <- NULL
  supci <- NULL
  for (i in 11:(dim(dd)[2])) {
    moy <- c(moy, mean(dd[,i]))
    infci <- c(infci, boot_sd_ci(dd[,i])[2])
    supci <- c(supci, boot_sd_ci(dd[,i])[3])
  }
  
  # confidence interval
  arrows(xrange[-c(1,2)], infci, xrange[-c(1,2)], supci, length=0.05, angle=90, code=3, col = "black", lwd = pts-5)
  
  # trajectory
  points(x = xrange[-1], y = moy, type = 'l', pch = 1, col = "black", lwd = pts-1)
}

# plot.new()
{
  # divide into four boxes
  layout(matrix(c(1,2), nrow = 1, byrow = F))
  par(mar = c(3, 3, 1, 1))
  
  # set space for x and y titles
  par(oma = c(3, 2.5, 0, 0))
  
  # upper left: CTL
  multi.mean.pop.mod(df = popul, upth = 0, bubo = 0, bura = 0.8, tmax = 20)
  mtext("CTL", side = 3, line = -2, adj = 0.98, cex = 2, col = "grey40")
  # upper right: ATI
  multi.mean.pop.mod(df = popul, upth = 0.5, bubo = 0.3, bura = 0.8, tmax = 20)
  mtext("TRJ", side = 3, line = -2, adj = 0.98, cex = 2, col = "grey40")
  # # lower right: TRJ
  # multi.mean.pop.mod(popul.trj, upth = 0.3, bubo = 0, tmax = 20)
  # mtext("TRJ", side = 3, line = -2, adj = 0.98, cex = 2, col = "grey40")
  # lower left: legend
  # legend( 1, 1,             # Location of legend
  #   # "right",
  #   xpd = TRUE,                          # Allow drawing outside plot area
  #   # ncol = 2,
  #   # xjust = 0,                           # Left justify legend box on x
  #   # yjust = 0.5,                          # Center legend box on y
  #   legend = c("Target",
  #              "PT range",
  #              "Replicate", 
  #              "Mean"),
  #           col = c("darkgreen",                 # Legend Element Colors
  #                   "green",
  #                   "grey",
  #                   "black"),          
  #           pch = c(NA_integer_,
  #                   NA_integer_,
  #                   NA_integer_,
  #                   NA_integer_ ),                      # Legend Element Styles          
  #           lty = c(2,
  #                   NA_integer_,
  #                   1,
  #                   1),       
  #           cex = 2,
  #           fill = "green"
  #           # cex = 0.6,
  #           # title = "Strategies") #,                  # Legend Title
  #           # title.col = gray(.2) ,                # Legend title color
  #           # box.lty = 1,                         # Legend box line type
  #           # box.lwd = 1)                         # Legend box line width
  # )
  
  # x title
  mtext("Time", outer = TRUE, cex = 2.5, side = 1, line = 0.5, adj = 0.52)
  # y title
  mtext("Population size", outer = TRUE, cex = 2.5, side = 2, line = 0, adj = 0.5)
  
}
# dev.off()

setwd(dir = "~/Desktop/PhD/GitKraken/gmse_fork_RQ1/")
dev.print(device = png, file = "pop-BDG-bigLabs.png", width = 1000, height = 600)
