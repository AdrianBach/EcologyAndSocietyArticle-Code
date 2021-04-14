# Adrian Bach
# PhD : Using AI to improve decision-making in conservation conflicts
# University of Stirling

#### Bootstrapping function ####

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

#### Stats function ####

OTI_stats <- function(df, ts, omit.extinction = FALSE) {
  
  df <- as.data.frame(df)
  
  # levels of OTI parameters
  upd_thr <- levels(as.factor(df$at))
  bud_bon <- levels(as.factor(df$bb))
  
  # number of samples for bootstrap 
  nbs <- 1000
  
  if (omit.extinction == "TRUE") { 
    
    # a loop to calculate extinction freq
    sub <- subset(df, at == 0 & bb == 0)
    
    # NULL tab
    ext_freq <- NULL
    sd_ci <- boot_sd_ci(sub$extinct, itr = nbs)
    res <- c(sum(sub$extinct)/dim(sub)[1], sd_ci[1], sd_ci[2],sd_ci[3])
    ext_freq <- rbind(ext_freq, res)
    
    for (i in 2:length(upd_thr)) {
      for (j in 1:length(bud_bon)) {
        sub <- subset(df, at == upd_thr[i] & bb == bud_bon[j])
        sd_ci <- boot_sd_ci(sub$extinct, itr = nbs)
        res <- c(sum(sub$extinct)/dim(sub)[1], sd_ci[1], sd_ci[2],sd_ci[3])
        ext_freq <- rbind(ext_freq, res)
      }
    }
    
    print("Ommiting replicates that resulted in Resource extinction")
    df <- subset(df, extinct == 0)
    
    # levels of OTI parameters
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
  sub <- subset(df, at == upd_thr[1] & bb == bud_bon[1])
  
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
  res_str <- c(nbrep,mean(sub$budget),upd_thr[1],bud_bon[1], ext)
  
  ## mean, sd and 95ci for each proxy
  
  # Actual Resource population deviation from Manager's target
  sd_ci <- boot_sd_ci(sub$act_dev, itr = nbs)
  res_str <- c(res_str, mean(sub$act_dev), sd_ci[1], sd_ci[2], sd_ci[3])
  
  # Absolute value of the Actual Resource population deviation from Manager's target
  sd_ci <- boot_sd_ci(sub$abs_act_dev, itr = nbs)
  res_str <- c(res_str, mean(sub$abs_act_dev), sd_ci[1], sd_ci[2], sd_ci[3])
  
  # Users' total final yield
  # sd_ci <- boot_sd_ci(sub$fin_yield, itr = nbs)
  # res_str <- c(res_str, mean(sub$fin_yield), sd_ci[1], sd_ci[2], sd_ci[3])
  sd_ci <- boot_sd_ci(sub$fin_yield/40000, itr = nbs)
  res_str <- c(res_str, mean(sub$fin_yield/40000), sd_ci[1], sd_ci[2], sd_ci[3])
  
  # Difference between the highest and the lowest yield
  sd_ci <- boot_sd_ci(sub$max_diff_yield, itr = nbs)
  res_str <- c(res_str, mean(sub$max_diff_yield), sd_ci[1], sd_ci[2], sd_ci[3])
  
  # Percentage of time steps of non-updating
  sd_ci <- boot_sd_ci(sub$inac_ts, itr = nbs)
  res_str <- c(res_str, mean(sub$inac_ts), sd_ci[1], sd_ci[2], sd_ci[3])
  
  # Sum of squared deviation from target
  sd_ci <- boot_sd_ci(sub$SAD/sub$final_ts, itr = nbs)
  res_str <- c(res_str, mean(sub$SAD/sub$final_ts), sd_ci[1], sd_ci[2], sd_ci[3])
  
  # binding the string to the tab
  res_tab <- rbind(res_tab, as.numeric(res_str))
  
  ## loop over the other OTI parameters
  # for each at value
  for (i in 2:length(upd_thr)) {
    
    # for each bb value
    for (j in 1:length(bud_bon)) {
      
      # increment tracker
      if (omit.extinction == TRUE) {
        zz <- zz + 1
      }
      
      # initiate a string
      res_str <- NULL
      
      # subset
      sub <- subset(df, at == upd_thr[i] & bb == bud_bon[j])
      
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
      res_str <- c(nbrep,mean(sub$budget),upd_thr[i],bud_bon[j], ext)
      
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
        # sd_ci <- boot_sd_ci(sub$fin_yield, itr = nbs)
        # res_str <- c(res_str, mean(sub$fin_yield), sd_ci[1], sd_ci[2], sd_ci[3])
        sd_ci <- boot_sd_ci(sub$fin_yield/40000, itr = nbs)
        res_str <- c(res_str, mean(sub$fin_yield/40000), sd_ci[1], sd_ci[2], sd_ci[3])
        
        # Difference between the highest and the lowest yield
        sd_ci <- boot_sd_ci(sub$max_diff_yield, itr = nbs)
        res_str <- c(res_str, mean(sub$max_diff_yield), sd_ci[1], sd_ci[2], sd_ci[3])
        
        # Percentage of time steps of non-updating
        sd_ci <- boot_sd_ci(sub$inac_ts, itr = nbs)
        res_str <- c(res_str, mean(sub$inac_ts), sd_ci[1], sd_ci[2], sd_ci[3])
        
        # # Percentage of time steps of K overshooting
        # sd_ci <- boot_sd_ci(sub$overK, itr = nbs)
        # res_str <- c(res_str, mean(sub$overK), sd_ci[1], sd_ci[2], sd_ci[3])
        
        # Sum of squared deviation from target
        sd_ci <- boot_sd_ci(sub$SAD/sub$final_ts, itr = nbs)
        res_str <- c(res_str, mean(sub$SAD/sub$final_ts), sd_ci[1], sd_ci[2], sd_ci[3])
        
      } else {
        
        print(paste("parameter set with UT = ", as.numeric(upd_thr[i])*100, "% and BB = ", as.numeric(bud_bon[j])*100, "% has less than 2 replicates"))
        
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
        
        # # Percentage of time steps of K overshooting
        # res_str <- c(res_str, sub$overK, 0, sub$overK, sub$overK)
        
        # Sum of squared deviation from target
        res_str <- c(res_str, sub$SAD/sub$final_ts, 0, sub$SAD, sub$SAD)
      } # end else loop on nbrep
      
      # binding the string to the tab
      res_tab <- rbind(res_tab, as.numeric(res_str))
      
    } # end for loop on budget bonus
  } # end for loop on update threshold
  
  # Array of column names
  colnames(res_tab) <- c("rep", "budget", "at", "bb", "ext_prob", "ext_prob_sd", "ext_prob_95ci_inf", "ext_prob_95ci_sup", "act_dev", "act_dev_sd", "act_dev_95ci_inf", "act_dev_95ci_sup", "abs_act_dev", "abs_act_dev_sd", "abs_act_dev_95ci_inf", "abs_act_dev_95ci_sup", "fin_yield", "fin_yield_sd", "fin_yield_95ci_inf", "fin_yield_95ci_sup", "max_diff_yield", "max_diff_yield_sd", "max_diff_yield_95ci_inf", "max_diff_yield_95ci_sup", "inac_ts", "inac_ts_sd", "inac_ts_95ci_inf", "inac_ts_95ci_sup", "SumAbsDev", "SumAbsDev_sd", "SumAbsDev_95ci_inf", "SumAbsDev_95ci_sup")
  
  res_tab <- as.data.frame(res_tab)
  
  return(res_tab)
  
} # end function

#### change dev, yield in case of extinction ####

brut <- read.csv("~/Desktop/PhD/GitKraken/gmse_fork_RQ1/mem-noreset-results/mem-noreset-merged-results/ATI-mem-noreset-merged-results.csv")[,-1]
brut.mod <- brut
brut.mod[which(brut.mod$extinct == 1), 6] <- -1
brut.mod[which(brut.mod$extinct == 1), 8] <- 40000

stat.mod <- OTI_stats(df = brut.mod, ts = 20, omit.extinction = F)