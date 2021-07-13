# Adrian Bach
# PhD : Using AI to improve decision-making in conservation conflicts
# University of Stirling

#### libraries ####
library("colorspace")
library("stringr")

### functions needed ####

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

#### Figure to compare alternative strategies against CTL strategy for each management outcome ####

OTI_vs_FLI_plot <- function(df, upth, goal = c(0:5), variance = c("sd","ci"), nb_replicates, nemar = FALSE) {
  
  # text sizes
  labsize = 3
  ticksize = 2
  linewd = 2
  pointsize = 2
  barsize = 0.1
  
  # subsetting
  fli <- subset(df, at == 0)
  oti <- subset(df, at == upth)
  
  trans <- adjustcolor("lightgrey",alpha.f=0.7)
  
  # extinction frequency 
  if (goal == 0) {
    
    # get values
    fli_avg <- fli$ext_prob
    oti_avg <- oti$ext_prob
    
    ## macnemar test for comparison of the frequencies for paired samples, for each BB value to the FLI
    if (nemar == TRUE) {
      
      # initiate vector for pvalues
      pv <- NULL
      for (i in 1:length(oti_avg)) {
        # table of contingency
        tc <- matrix(c(fli$ext_prob+oti_avg[i], (1-fli$ext_prob)+oti_avg[i], fli$ext_prob+(1-oti_avg[i]), (1-fli$ext_prob)+(1-oti_avg[i]))*nb_replicates, nrow = 2, ncol = 2)
        if (tc[1,1]>10 & tc[1,2]>10 & tc[2,1]>10 & tc[2,2]>10) {
          pv <- c(pv, mcnemar.test(tc, correct = F)$p.value)
        } else {
          pv <- c(pv, mcnemar.test(tc, correct = T)$p.value)
        }      
      }
      
      # convert p-values in stars
      # initiate a vector for the stars
      st <- rep(NA, length(pv))
      for (i in 1:length(pv)) {
        st[i] <- ""
        if(pv[i] < 0.05 & pv[i] >= 0.01) {
          st[i] <- "*"
        }
        if (pv[i] < 0.01 & pv[i] >= 0.001) {
          st[i] <- "**"
        }
        if (pv[i] < 0.001) {
          st[i] <- "***"
        }
      }
    }
    
    if (variance == "sd") {
      
      # FLI strat
      fli_sd <- fli$ext_prob_sd
      
      # OTI strat
      oti_sd <- oti$ext_prob_sd
      
      # plotting
      xrange <- seq(0,100,10)
      xtendrange <- seq(-10,110,1)
      flimax <- rep(fli_avg+fli_sd, length(xtendrange))
      flimin <- rep(fli_avg-fli_sd, length(xtendrange))
      xoti <- oti$bb*100
      
      # diag <- barplot(oti_avg, col = "lightblue", space = 1, width = 4, names.arg = xoti,
      #                 xlab = "Budget bonus\n(in % of initial budget)", ylab = paste("Extinction frequency (N =", nb_replicates,")"), ylim = c(0, max(fli$ext_prob, max(oti_avg))+0.1)) 
      #                 # ylim = c(0,1), xlim = c(0,100)) # ,
      #                 # main = paste("UT = ", upth*100,"%"))
      # abline(h = fli$ext_prob, lty = 1, lwd = 2, col = "black")
      
      # plot base
      plot(1, type = "n", xlab = "Budget bonus\n(in % of initial budget)", ylab = "Extinction frequency (+/- SD)",
           ylim = c(0,1), xlim = c(0,100)) # , max(max(fli_avg+fli_sd),max(oti_avg+oti_sd))
      # main = paste("UT = ", upth*100,"%"))
      polygon(c(xtendrange,rev(xtendrange)),c(flimax,rev(flimin)),col=trans, border = "grey")
      abline(h=fli_avg, lwd=linewd)
      
      points(y = oti_avg, x = xoti, pch = 16, col = "blue", cex = pointsize)
      # arrows(xoti, oti_avg-oti_sd_neg, xoti, oti_avg+oti_sd, length=0.05, angle=90, code=3, col = "blue")
      arrows(xoti, oti_avg-oti_sd, xoti, oti_avg+oti_sd, length=barsize, angle=90, code=3, col = "darkblue")
      
      if (nemar == TRUE) {
        # add the stars above the bars
        text(diag,oti_avg+oti_sd+0.05,as.character(st),cex=1)
      }
    }  
    
    if (variance == "ci") {
      
      # FLI strat
      fli_ci_inf <- fli$ext_prob_95ci_inf
      fli_ci_sup <- fli$ext_prob_95ci_sup
      
      # OTI strat
      oti_ci_inf <- oti$ext_prob_95ci_inf
      oti_ci_sup <- oti$ext_prob_95ci_sup
      
      # plotting
      xrange <- seq(0,100,10)
      xtendrange <- seq(-10,110)
      flimax <- rep(fli_ci_sup, length(xtendrange))
      flimin <- rep(fli_ci_inf, length(xtendrange))
      trans <- adjustcolor("lightgrey",alpha.f=0.7)
      xoti <- oti$bb*100
      
      # diag <- barplot(oti_avg, col = "lightblue", space = 1, width = 4, names.arg = xoti,
      #                 xlab = "Budget bonus\n(in % of initial budget)", ylab = paste("Extinction frequency (N =", nb_replicates,")"), ylim = c(0, max(fli$ext_prob, max(oti_avg))+0.1)) 
      #                 # ylim = c(0,1), xlim = c(0,100)) # ,
      #                 # main = paste("UT = ", upth*100,"%"))
      # abline(h = fli$ext_prob, lty = 1, lwd = 2, col = "black")
      
      # plot base
      plot(1, type = "n", xlab = "n", ylab = "Extinction frequency", cex.lab = labsize, cex.axis = ticksize,
           ylim = c(0, max(max(oti$ext_prob_95ci_sup)+0.1, max(flimax)+0.1)), xlim = c(0,100)) # ,max(fli_ci_sup,max(oti_ci_sup)
      # main = paste("UT = ", upth*100,"%"))
      polygon(c(xtendrange,rev(xtendrange)),c(flimax,rev(flimin)),col=trans, border = "grey")
      abline(h=fli_avg, lwd=linewd)
      
      points(y = oti_avg, x = xoti, pch = 16, col = "blue", cex = pointsize)
      # arrows(xoti, oti_avg-oti_sd_neg, xoti, oti_avg+oti_sd, length=0.05, angle=90, code=3, col = "blue")
      arrows(xoti, oti_ci_inf, xoti, oti_ci_sup, length=barsize, angle=90, code=3, col = "darkblue")
      
      if (nemar == TRUE) {
        # add the stars above the bars
        text(diag,oti_avg+oti_sd+0.05,as.character(st),cex=1)
      }
    }  
  }
  
  # Deviation from target
  if (goal == 1) {
    
    # get means, sd, and 95ci and plot
    fli_avg <- fli$act_dev*100
    oti_avg <- oti$act_dev*100
    
    if (variance == "sd") {
      
      # FLI strat
      fli_sd <- fli$act_dev_sd*100
      # # prevent the sd range to go over the borders
      # fli_sd_neg <- ifelse(test = fli_avg-fli_sd < -100, fli_sd+(fli_avg-fli_sd+100), fli_sd)
      
      # OTI strat
      oti_sd <- oti$act_dev_sd*100
      # # prevent the sd range to go over the borders
      # oti_sd_neg <- ifelse(test = oti_avg-oti_sd < -100, oti_sd+(oti_avg-oti_sd+100), oti_sd)
      
      # plotting
      xrange <- seq(0,100,10)
      xtendrange <- seq(-10,110)
      flimax <- rep(fli_avg+fli_sd, length(xtendrange))
      # flimin <- rep(fli_avg-fli_sd_neg, length(xtendrange))
      flimin <- rep(fli_avg-fli_sd, length(xtendrange))
      trans <- adjustcolor("lightgrey",alpha.f=0.7)
      xoti <- oti$bb*100
      
      # plot base
      plot(1, type = "n", xlab = "Budget bonus\n(in % of initial budget)", ylab = "Resource population deviation from MT\n(in %, mean +/- SD)",
           ylim = c(-100,ifelse(max(max(fli_avg+fli_sd)+10,max(oti_avg+oti_sd)+10)<0, 0, max(max(fli_avg+fli_sd)+10,max(oti_avg+oti_sd)+10))), xlim = c(0,100)) # ,
      # main = paste("UT = ", upth*100,"%"))
      polygon(c(xtendrange,rev(xtendrange)),c(flimax,rev(flimin)),col=trans, border = "grey")
      abline(h = fli_avg, lwd = 2)
      abline(h = 0, lty = 2, lwd = 1.5, col = "darkgreen")
      # abline(h = -100, lty = 2, lwd = 1.2, col = "red")
      
      points(y = oti_avg, x = xoti, pch = 16, col = "blue", cex = pointsize)
      # arrows(xoti, oti_avg-oti_sd_neg, xoti, oti_avg+oti_sd, length=0.05, angle=90, code=3, col = "blue")
      arrows(xoti, oti_avg-oti_sd, xoti, oti_avg+oti_sd, length=barsize, angle=90, code=3, col = "darkblue")
      
    }
    
    if (variance == "ci") {
      
      # FLI strat
      # fli_95ci <- fli$act_dev_95ci*100
      # # prevent the 95ci range to go over the borders
      # fli_95ci_neg <- ifelse(test = fli_avg-fli_95ci < -100, fli_95ci+(fli_avg-fli_95ci+100), fli_95ci)
      fli_95ci_inf <- fli$act_dev_95ci_inf*100
      fli_95ci_sup <- fli$act_dev_95ci_sup*100
      
      # OTI strat
      # oti_95ci <- oti$act_dev_95ci*100
      # # prevent the 95ci range to go over the borders
      # oti_95ci_neg <- ifelse(test = oti_avg-oti_95ci < -100, oti_95ci+(oti_avg-oti_95ci+100), oti_95ci)
      oti_95ci_inf <- oti$act_dev_95ci_inf*100
      oti_95ci_sup <- oti$act_dev_95ci_sup*100
      
      # plotting
      xrange <- seq(0,100,10)
      xtendrange <- seq(-10,110)
      flimax <- rep(fli_95ci_sup, length(xtendrange))
      flimin <- rep(fli_95ci_inf, length(xtendrange))
      trans <- adjustcolor("lightgrey",alpha.f=0.7)
      xoti <- oti$bb*100
      
      # plot base
      plot(1, type = "n", xlab = "Budget bonus (%)", ylab = "Deviation from target (%)", cex.lab = labsize, cex.axis = ticksize,
           ylim = c(min(min(flimin)-5, min(oti_95ci_inf)-5),ifelse(max(max(fli_95ci_sup)+5,max(oti_95ci_sup)+5)<0, 0, max(max(fli_95ci_sup)+5,max(oti_95ci_sup)+5))), xlim = c(0,100)) # ,
      # main = paste("UT = ", upth*100,"%"))
      polygon(c(xtendrange,rev(xtendrange)),c(flimax,rev(flimin)),col=trans, border = "grey")
      abline(h = fli_avg, lwd = 2)
      abline(h = 0, lty = 2, lwd = 1.5, col = "darkgreen")
      # abline(h = -100, lty = 2, lwd = 1.2, col = "red")
      
      points(y = oti_avg, x = xoti, pch = 16, col = "blue", cex = pointsize)
      arrows(xoti, oti_95ci_inf, xoti, oti_95ci_sup, length=barsize, angle=90, code=3, col = "darkblue")
    }
  }
  
  # Final yield
  if (goal == 2) {
    
    # get means, sd, and 95ci and plot
    fli_avg <- fli$fin_yield*100
    oti_avg <- oti$fin_yield*100
    
    if (variance == "sd") {
      
      # FLI strat
      fli_sd <- fli$fin_yield_sd*100
      # # prevent the sd range to go over the borders
      # fli_sd_neg <- ifelse(test = fli_avg-fli_sd < 0, fli_sd+(fli_avg-fli_sd), fli_sd)
      # fli_sd_pos <- ifelse(test = fli_avg+fli_sd > 100, fli_sd-(fli_avg+fli_sd-100), fli_sd)
      
      # OTI strat
      oti_sd <- oti$fin_yield_sd*100
      # # prevent the sd range to go over the borders
      # oti_sd_neg <- ifelse(test = oti_avg-oti_sd < 0, oti_sd+(oti_avg-oti_sd+100), oti_sd)
      # oti_sd_pos <- ifelse(test = oti_avg+oti_sd > 100, oti_sd-(oti_avg+oti_sd-100), oti_sd)
      
      # plotting
      xrange <- seq(0,100,10)
      xtendrange <- seq(-10,110)
      # flimax <- rep(fli_avg+fli_sd_pos, length(xtendrange))
      # flimin <- rep(fli_avg-fli_sd_neg, length(xtendrange))
      flimax <- rep(fli_avg+fli_sd, length(xtendrange))
      flimin <- rep(fli_avg-fli_sd, length(xtendrange))
      trans <- adjustcolor("lightgrey",alpha.f=0.7)
      xoti <- oti$bb*100
      
      # plot base
      plot(1, type = "n", xlab = "Budget bonus\n(in % of initial budget)", ylab = "Sum of Users final budgets\n(in k-a.b.u, mean +/- SD)",
           ylim = c(min(min(fli_avg+fli_sd)-5,min(oti_avg+oti_sd))-5,100), xlim = c(0,100)) # ,
      # main = paste("UT = ", upth*100,"%"))
      polygon(c(xtendrange,rev(xtendrange)),c(flimax,rev(flimin)),col=trans, border = "grey")
      abline(h = fli_avg, lwd = 2)
      abline(h = 100, lty = 2, lwd = 1.5, col = "darkgreen")
      
      points(y = oti_avg, x = xoti, pch = 16, col = "blue", cex = pointsize)
      # arrows(xoti, oti_avg-oti_sd_neg, xoti, oti_avg+oti_sd_pos, length=0.05, angle=90, code=3, col = "blue")
      arrows(xoti, oti_avg-oti_sd, xoti, oti_avg+oti_sd, length=barsize, angle=90, code=3, col = "darkblue")
    }
    
    if (variance == "ci") {
      
      # FLI strat
      # fli_95ci <- fli$fin_yield_95ci/100
      # # prevent the 95ci range to go over the borders
      # fli_95ci_neg <- ifelse(test = fli_avg-fli_95ci < 0, fli_95ci+(fli_avg-fli_95ci), fli_95ci)
      # fli_95ci_pos <- ifelse(test = fli_avg+fli_95ci > 100, fli_95ci-(fli_avg+fli_95ci-100), fli_95ci)
      fli_95ci_inf <- fli$fin_yield_95ci_inf*100
      fli_95ci_sup <- fli$fin_yield_95ci_sup*100
      
      # OTI strat
      # oti_95ci <- oti$fin_yield_95ci/100
      # # prevent the 95ci range to go over the borders
      # oti_95ci_neg <- ifelse(test = oti_avg-oti_95ci < 0, oti_95ci+(oti_avg-oti_95ci), oti_95ci)
      # oti_95ci_pos <- ifelse(test = oti_avg+oti_95ci > 100, oti_95ci-(oti_avg+oti_95ci-100), oti_95ci)
      oti_95ci_inf <- oti$fin_yield_95ci_inf*100
      oti_95ci_sup <- oti$fin_yield_95ci_sup*100
      
      # plotting
      xrange <- seq(0,100,10)
      xtendrange <- seq(-10,110)
      flimax <- rep(fli_95ci_sup, length(xtendrange))
      flimin <- rep(fli_95ci_inf, length(xtendrange))
      trans <- adjustcolor("lightgrey",alpha.f=0.7)
      xoti <- oti$bb*100
      
      # plot base
      plot(1, type = "n", xlab = "Budget bonus (%)", ylab = "Farmers' final yield (%)", cex.lab = labsize, cex.axis = ticksize,
           ylim = c(min(min(fli_95ci_inf)-5,min(oti_95ci_inf))-5,100), xlim = c(0,100)) # , 
      # main = paste("UT = ", upth*100,"%"))
      polygon(c(xtendrange,rev(xtendrange)),c(flimax,rev(flimin)),col=trans, border = "grey")
      abline(h = fli_avg, lwd = 2)
      abline(h = 100, lty = 2, lwd = 1.5, col = "darkgreen")
      
      points(y = oti_avg, x = xoti, pch = 16, col = "blue", cex = pointsize)
      arrows(xoti, oti_95ci_inf, xoti, oti_95ci_sup, length=barsize, angle=90, code=3, col = "darkblue")
    }
  }
  
  if (goal == 3) {
    
    # get means, sd, and 95ci and plot
    fli_avg <- fli$max_diff_yield*100
    oti_avg <- oti$max_diff_yield*100
    
    if (variance == "sd") {
      
      # FLI strat
      fli_sd <- fli$max_diff_yield_sd*100
      # # prevent the sd range to go over the borders
      # fli_sd_neg <- ifelse(test = fli_avg-fli_sd < 0, fli_sd+(fli_avg-fli_sd), fli_sd)
      # fli_sd_pos <- ifelse(test = fli_avg+fli_sd > 100, fli_sd-(fli_avg+fli_sd-100), fli_sd)
      
      # OTI strat
      oti_sd <- oti$max_diff_yield_sd*100
      # # prevent the sd range to go over the borders
      # oti_sd_neg <- ifelse(test = oti_avg-oti_sd < 0, oti_sd+(oti_avg-oti_sd), oti_sd)
      # oti_sd_pos <- ifelse(test = oti_avg+oti_sd > 100, oti_sd-(oti_avg+oti_sd-100), oti_sd)
      
      # plotting
      xrange <- seq(0,100,10)
      # xtendrange <- seq(-10,110)
      # flimax <- rep(fli_avg+fli_sd_pos, length(xtendrange))
      # flimin <- rep(fli_avg-fli_sd_neg, length(xtendrange))
      flimax <- rep(fli_avg+fli_sd, length(xtendrange))
      flimin <- rep(fli_avg-fli_sd, length(xtendrange))
      trans <- adjustcolor("lightgrey",alpha.f=0.7)
      xoti <- oti$bb*100
      
      # plot base
      plot(1, type = "n", xlab = "Budget bonus (%)", ylab = "Maximum difference between Users yields\n(in % of the highest yield, mean +/- SD)",
           ylim = c(0,max(max(fli_avg+fli_sd),max(oti_avg+oti_sd))+5), xlim = c(0,100)) #,
      # main = paste("UT = ", upth*100,"%"))
      polygon(c(xtendrange,rev(xtendrange)),c(flimax,rev(flimin)),col=trans, border = "grey")
      abline(h = fli_avg, lwd = 2)
      abline(h = 0, lty = 2, lwd = 1.5, col = "darkgreen")
      
      # points(y = oti_avg, x = xoti, pch = 16, col = "blue")
      # arrows(xoti, oti_avg-oti_sd_neg, xoti, oti_avg+oti_sd_pos, length=0.05, angle=90, code=3, col = "blue")
      arrows(xoti, oti_avg-oti_sd, xoti, oti_avg+oti_sd, length=barsize, angle=90, code=3, col = "darkblue")
    }
    
    if (variance == "ci") {
      
      # FLI strat
      # fli_95ci <- fli$max_diff_yield_95ci*100
      # # prevent the 95ci range to go over the borders
      # fli_95ci_neg <- ifelse(test = fli_avg-fli_95ci < 0, fli_95ci+(fli_avg-fli_95ci), fli_95ci)
      # fli_95ci_pos <- ifelse(test = fli_avg+fli_95ci > 100, fli_95ci-(fli_avg+fli_95ci-100), fli_95ci)
      fli_95ci_inf <- fli$max_diff_yield_95ci_inf*100
      fli_95ci_sup <- fli$max_diff_yield_95ci_sup*100
      
      # OTI strat
      # oti_95ci <- oti$max_diff_yield_95ci*100
      # # prevent the 95ci range to go over the borders
      # oti_95ci_neg <- ifelse(test = oti_avg-oti_95ci < 0, oti_95ci+(oti_avg-oti_95ci), oti_95ci)
      # oti_95ci_pos <- ifelse(test = oti_avg+oti_95ci > 100, oti_95ci-(oti_avg+oti_95ci-100), oti_95ci)
      oti_95ci_inf <- oti$max_diff_yield_95ci_inf*100
      oti_95ci_sup <- oti$max_diff_yield_95ci_sup*100
      
      # plotting
      xrange <- seq(0,100,10)
      xtendrange <- seq(-10,110)
      # flimax <- rep(fli_avg+fli_95ci_pos, length(xtendrange))
      # flimin <- rep(fli_avg-fli_95ci_neg, length(xtendrange))
      flimax <- rep(fli_95ci_sup, length(xtendrange))
      flimin <- rep(fli_95ci_inf, length(xtendrange))
      trans <- adjustcolor("lightgrey",alpha.f=0.7)
      xoti <- oti$bb*100
      
      # plot base
      plot(1, type = "n", xlab = "Budget bonus\n(in % of initial budget)", ylab = "Maximum difference between Users yields\n(in % of the highest yield, mean +/- 95CI)",
           # ylim = c(0,max(max(fli_avg+fli_sd),max(oti_avg+oti_sd))+5), xlim = c(0,100)) # ,
           ylim = c(0,max(max(fli_95ci_sup)+1,max(oti_95ci_sup))+1), xlim = c(0,100)) # ,
      # main = paste("UT = ", upth*100,"%"))
      polygon(c(xtendrange,rev(xtendrange)),c(flimax,rev(flimin)),col=trans, border = "grey")
      abline(h = fli_avg, lwd = 2)
      abline(h = 0, lty = 2, lwd = 1.5, col = "darkgreen")
      
      points(y = oti_avg, x = xoti, pch = 16, col = "blue", cex = pointsize)
      # arrows(xoti, oti_avg-oti_95ci_neg, xoti, oti_avg+oti_95ci_pos, length=0.05, angle=90, code=3, col = "blue")
      arrows(xoti, oti_95ci_inf, xoti, oti_95ci_sup, length=barsize, angle=90, code=3, col = "darkblue")
    }
  }
  
  # time spent waiting
  if (goal == 4) {
    
    # get means, sd, 95ci and plot
    oti_avg <- oti$inac_ts*100
    
    if (variance == "sd") {
      
      # OTI strat
      oti_sd <- oti$inac_ts_sd*100
      # # prevent the sd range to go over the borders
      # oti_sd_neg <- ifelse(test = oti_avg-oti_sd < 0, oti_sd+(oti_avg-oti_sd), oti_sd)
      # oti_sd_pos <- ifelse(test = oti_avg+oti_sd > 100, oti_sd-(oti_avg+oti_sd-100), oti_sd)
      
      # plotting
      xoti <- oti$bb*100
      
      # barplot
      # diag = barplot(oti_avg, col = "lightblue", space = 1, width = 4, names.arg = xoti,
      #                xlab = "Budget bonus\n(in % of initial budget)", ylab = "Time steps without updating\n(in %, mean +/- SD)",
      #                ylim = c(0,max(oti_avg+oti_sd))) # ,
      # ylim = c(0,max(oti_avg+oti_sd_pos))) # ,
      # ylim = c(0,100), xlim = c(0,100)) # ,
      # main = paste("UT = ", upth*100,"%"))
      # arrows(diag, oti_avg-oti_sd_neg, diag, oti_avg+oti_sd_pos, length=0.03, angle=90, code=3, col = "black")
      
      # plot base
      plot(1, type = "n", xlab = "Budget bonus\n(in % of initial budget)", ylab = "Time steps without updating (%)", cex.lab = labsize, cex.axis = ticksize,
           # ylim = c(0,max(max(fli_avg+fli_sd),max(oti_avg+oti_sd))+5), xlim = c(0,100)) # ,
           ylim = c(0,max(max(fli_95ci_sup)+1,max(oti_95ci_sup))+1), xlim = c(0,100)) # ,
      # main = paste("UT = ", upth*100,"%"))
      polygon(c(xtendrange,rev(xtendrange)),c(flimax,rev(flimin)),col=trans, border = "grey")
      abline(h = fli_avg, lwd = 2)
      abline(h = 0, lty = 2, lwd = 1.5, col = "darkgreen")
      
      points(y = oti_avg, x = xoti, pch = 16, col = "blue", cex = pointsize)
      
      arrows(diag, oti_avg-oti_sd, diag, oti_avg+oti_sd, length=0.03, angle=90, code=3, col = "black")
    }
    
    if (variance == "ci") {
      
      # OTI strat
      # oti_95ci <- oti$inac_ts_95ci*100
      # # prevent the 95ci range to go over the borders
      # oti_95ci_neg <- ifelse(test = oti_avg-oti_95ci < 0, oti_95ci+(oti_avg-oti_95ci), oti_95ci)
      # oti_95ci_pos <- ifelse(test = oti_avg+oti_95ci > 100, oti_95ci-(oti_avg+oti_95ci-100), oti_95ci)
      oti_95ci_inf <- oti$inac_ts_95ci_inf*100
      oti_95ci_sup <- oti$inac_ts_95ci_sup*100
      
      # plotting
      xoti <- oti$bb*100
      
      # barplot
      # diag = barplot(oti_avg, col = "lightblue", space = 1, width = 4, names.arg = xoti,
      #                xlab = "Budget bonus\n(in % of initial budget)", ylab = "Time steps without updating\n(in %, mean +/- 95CI)",
      #                ylim = c(0,max(oti_95ci_sup)))
      # ylim = c(0,max(oti_avg+oti_95ci_pos))) 
      # , xlim = c(0,100)) # ,
      # main = paste("UT = ", upth*100,"%"))
      
      # plot base
      plot(1, type = "n", xlab = "Budget bonus (%)", ylab = "Time steps without updating (%)", cex.lab = labsize, cex.axis = ticksize,
           # ylim = c(0,max(max(fli_avg+fli_sd),max(oti_avg+oti_sd))+5), xlim = c(0,100)) # ,
           ylim = c(0,max(oti_95ci_sup)+5), xlim = c(0,100)) # ,
      
      points(y = oti_avg, x = xoti, pch = 16, col = "blue", cex = pointsize)
      arrows(xoti, oti_95ci_inf, xoti, oti_95ci_sup, length=barsize, angle=90, code=3, col = "darkblue")
      # arrows(diag, oti_avg-oti_95ci_neg, diag, oti_avg+oti_95ci_pos, length=0.03, angle=90, code=3, col = "black")
    }
  }
  
  if (goal == 5) {
    
    # get means, sd, 95ci and plot
    oti_avg <- oti$SumAbsDev/10^6
    
    if (variance == "sd") {
      
      # OTI strat
      oti_sd <- oti$SumAbsDev_sd/10^6
      # # prevent the sd range to go over the borders
      # oti_sd_neg <- ifelse(test = oti_avg-oti_sd < 0, oti_sd+(oti_avg-oti_sd), oti_sd)
      # oti_sd_pos <- ifelse(test = oti_avg+oti_sd > 100, oti_sd-(oti_avg+oti_sd-100), oti_sd)
      
      # plotting
      xoti <- oti$bb*100
      
      # barplot
      diag = barplot(oti_avg, col = "lightblue", space = 1, width = 4, names.arg = xoti,
                     xlab = "Budget bonus\n(in % of initial budget)", ylab = "Sum of sq. deviation from target\n(ind.10^6 +/- SD)",
                     ylim = c(0,max(oti_avg+oti_sd))) # ,
      # ylim = c(0,max(oti_avg+oti_sd_pos))) # ,
      # ylim = c(0,100), xlim = c(0,100)) # ,
      # main = paste("UT = ", upth*100,"%"))
      # arrows(diag, oti_avg-oti_sd_neg, diag, oti_avg+oti_sd_pos, length=0.03, angle=90, code=3, col = "black")
      arrows(diag, oti_avg-oti_sd, diag, oti_avg+oti_sd, length=0.03, angle=90, code=3, col = "black")
    }
    
    if (variance == "ci") {
      
      # OTI strat
      # oti_95ci <- oti$inac_ts_95ci*100
      # # prevent the 95ci range to go over the borders
      # oti_95ci_neg <- ifelse(test = oti_avg-oti_95ci < 0, oti_95ci+(oti_avg-oti_95ci), oti_95ci)
      # oti_95ci_pos <- ifelse(test = oti_avg+oti_95ci > 100, oti_95ci-(oti_avg+oti_95ci-100), oti_95ci)
      oti_95ci_inf <- oti$SumAbsDev_95ci_inf/10^6
      oti_95ci_sup <- oti$SumAbsDev_95ci_sup/10^6
      
      # plotting
      xoti <- oti$bb*100
      
      # barplot
      diag = barplot(oti_avg, col = "lightblue", space = 1, width = 4, names.arg = xoti,
                     xlab = "Budget bonus\n(in % of initial budget)", ylab = "Sum of sq. deviation from target\n (ind.10^6 +/- 95%CI)", cex = 1.5,
                     ylim = c(0,max(oti_95ci_sup)))
      # ylim = c(0,max(oti_avg+oti_95ci_pos))) 
      # , xlim = c(0,100)) # ,
      # main = paste("UT = ", upth*100,"%"))
      arrows(diag, oti_95ci_inf, diag, oti_95ci_sup, length=0.03, angle=90, code=3, col = "black")
      # arrows(diag, oti_avg-oti_95ci_neg, diag, oti_avg+oti_95ci_pos, length=0.03, angle=90, code=3, col = "black")
    }
  }
}

#### Compile results of the comparison between the FLI and ATI strategies ####

OTI_diagnostic <- function(df, upth, variance = c("sd", "ci"), nb_replicates, omit.extinction = FALSE) {
  
  # divide into four boxes
  layout(matrix(c(1,2,3,4), nrow = 2, byrow = T))
  par(mar = c(2, 6, 1, 1))
  # par(mfrow=c(2,2))
  
  # set space for a x title
  par(oma = c(3, 0, 0, 0))
  # layout.show(n = 4)
  
  # upper left: extinction frequency (goal 0)
  OTI_vs_FLI_plot(df, upth, goal = 0, variance, nb_replicates)
  mtext(letters[1], side = 3, line = -2, adj = 0.98, cex = 2, col = "grey40")
  # upper right: Resource population deviation from MT (goal 1)
  OTI_vs_FLI_plot(df, upth, goal = 1, variance, nb_replicates)
  mtext(letters[2], side = 3, line = -2, adj = 0.98, cex = 2, col = "grey40")
  # lower right: sum of Users final yield
  OTI_vs_FLI_plot(df, upth, goal = 2, variance, nb_replicates)
  mtext(letters[3], side = 3, line = -3, adj = 0.98, cex = 2, col = "grey40")
  # lower left: time steps without intervention
  OTI_vs_FLI_plot(df, upth, goal = 4, variance, nb_replicates)
  mtext(letters[4], side = 3, line = -2, adj = 0.98, cex = 2, col = "grey40")
  
  # Global title
  # mtext(paste("UT =", upth*100, "%"), outer = TRUE, cex = 1, line = 1)
  # mtext(paste("UT =", upth*100, "% -", ifelse(omit.extinction == FALSE,"including", "excluding"), "extinctions"), outer = TRUE, cex = 1, line = 1.5)
  
  # x title
  mtext("Budget bonus (%)", outer = TRUE, cex = 2.5, side = 1, line = 1.5, adj = 0.53)
  
  # par(mfrow=c(1,1))
}

# Examples
OTI_diagnostic(df = stat.mod, upth = 0.3, variance = "ci", nb_replicates = 100, omit.extinction = F)
setwd(dir = "~/Desktop/PhD/GitKraken/gmse_fork_RQ1/mem-noreset-results/mem-noreset-merged-results/")
dev.print(device = png, file = "TRJ-PT30-bigLabs.png", width = 1000)

OTI_diagnostic(df = woe_stat, upth = 0.3, variance = "ci", nb_replicates = 100, omit.extinction = T)

######## Contour figures ########

library(plotly)
 
#### Extinction frequency ####

# build results matrix

stat <- read.csv("~/Desktop/PhD/GitKraken/gmse_fork_RQ1/mem-noreset-results/mem-noreset-merged-results/ATI-mem-noreset-stats.csv")

d <- stat[,-1]

if (is.na(subset(d, at == 0.3 & bb == 0.3)$ext_prob) == TRUE) {
  print("UT30BB30 is missing, extrapolating")
  d$ext_prob[27] <-  (d$ext_prob[26]+d$ext_prob[28])/2
  d$act_dev[27] <-  (d$act_dev[26]+d$act_dev[28])/2
  d$fin_yield[27] <-  (d$fin_yield[26]+d$fin_yield[28])/2
  d$max_diff_yield[27] <-  (d$max_diff_yield[26]+d$max_diff_yield[28])/2
  d$inac_ts[27] <-  (d$inac_ts[26]+d$inac_ts[28])/2
}

d$at <- d$at*100
d$bb <- d$bb*100

bubo <- levels(as.factor(d$bb))
upth <- levels(as.factor(d$at))

control <- d[1,]
for (i in 2:length(bubo)) {
  zz <- d[1,]
  zz[4] <- bubo[i]
  control <- rbind(control, zz)
}

d <- rbind(control, d[-1,])

# resmat <- matrix(data = seq(1:(length(bubo)*length(bura))), ncol = length(bura), nrow = length(bubo))
resmat <- matrix(data = d$ext_prob, ncol = length(upth), nrow = length(bubo))

fig <- plot_ly(
  x = bubo, 
  y = upth, 
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
  title = "Permissiveness (%)"#,
  # titlefont = f
)

fig <- fig %>% colorbar(title = "Extinction \n frequency")
fig <- fig %>% layout(xaxis = xlab, yaxis = ylab)
fig

#### Users' final yield ####

# resmat <- matrix(data = seq(1:(length(bubo)*length(bura))), ncol = length(bura), nrow = length(bubo))
resmat <- matrix(data = 100*d$fin_yield, ncol = length(upth), nrow = length(bubo))

fig <- plot_ly(
  x = bubo, 
  y = upth, 
  # z = matrix(data = d$ext_prob, ncol = length(bubo), nrow = length(bura)), 
  z = t(resmat),
  type = "contour",
  colorscale = list(c(0, 0.5, 1), c('red', 'orange', 'green')),
  autocontour = F,
  # contours = list(showlabels = TRUE),
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
  title = "Permissiveness (%)"#,
  # titlefont = f
)

fig <- fig %>% colorbar(title = "Final \n yield (%)")
fig <- fig %>% layout(xaxis = xlab, yaxis = ylab)
fig

#### Users' equity ####

resmat <- matrix(data = 100*d$max_diff_yield, ncol = length(upth), nrow = length(bubo))

fig <- plot_ly(
  x = bubo, 
  y = upth, 
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
  title = "Permissiveness (%)"#,
  # titlefont = f
)

fig <- fig %>% colorbar(title = "Equity\n(%)")
fig <- fig %>% layout(xaxis = xlab, yaxis = ylab)
fig

#### Distance to target ####

resmat <- matrix(data = 100*d$act_dev, ncol = length(upth), nrow = length(bubo))

fig <- plot_ly(
  x = bubo, 
  y = upth, 
  # z = matrix(data = d$ext_prob, ncol = length(bubo), nrow = length(bura)), 
  z = t(resmat),
  type = "contour",
  colorscale = list(c(0, 0.35, 0.5, 0.65, 1), c('red', 'orange', 'green', 'orange', 'red')),
  autocontour = F,
  # contours = list(showlabels = TRUE),
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
  title = "Permissiveness (%)"#,
  # titlefont = f
)

fig <- fig %>% colorbar(title = "Dev. from\ntarget(%)")
fig <- fig %>% layout(xaxis = xlab, yaxis = ylab)
fig

#### Distance to target at the LAST ts ####

d <- stat.mod

if (is.na(subset(d, at == 0.3 & bb == 0.3)$ext_prob) == TRUE) {
  print("UT30BB30 is missing, extrapolating")
  d$ext_prob[27] <-  (d$ext_prob[26]+d$ext_prob[28])/2
  d$act_dev[27] <-  (d$act_dev[26]+d$act_dev[28])/2
  d$fin_yield[27] <-  (d$fin_yield[26]+d$fin_yield[28])/2
  d$max_diff_yield[27] <-  (d$max_diff_yield[26]+d$max_diff_yield[28])/2
  d$inac_ts[27] <-  (d$inac_ts[26]+d$inac_ts[28])/2
}

d$at <- d$at*100
d$bb <- d$bb*100

bubo <- levels(as.factor(d$bb))
upth <- levels(as.factor(d$at))

control <- d[1,]
for (i in 2:length(bubo)) {
  zz <- d[1,]
  zz[4] <- bubo[i]
  control <- rbind(control, zz)
}

d <- rbind(control, d[-1,])

# resmat <- matrix(data = seq(1:(length(bubo)*length(bura))), ncol = length(bura), nrow = length(bubo))
resmat <- matrix(data = 100*d$act_dev, ncol = length(upth), nrow = length(bubo))

fig <- plot_ly(
  x = bubo, 
  y = upth, 
  # z = matrix(data = d$ext_prob, ncol = length(bubo), nrow = length(bura)), 
  z = t(resmat),
  type = "contour",
  colorscale = list(c(0, 0.35, 0.5, 0.65, 1), c('red', 'orange', 'green', 'orange', 'red')),
  autocontour = F,
  # contours = list(showlabels = TRUE),
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
  title = "Permissiveness (%)"#,
  # titlefont = f
)

fig <- fig %>% colorbar(title = "Dev. from\ntarget(%)")
fig <- fig %>% layout(xaxis = xlab, yaxis = ylab)
fig

#### sum yield at the LAST ts ####

d <- stat.mod

if (is.na(subset(d, at == 0.3 & bb == 0.3)$ext_prob) == TRUE) {
  print("UT30BB30 is missing, extrapolating")
  d$ext_prob[27] <-  (d$ext_prob[26]+d$ext_prob[28])/2
  d$act_dev[27] <-  (d$act_dev[26]+d$act_dev[28])/2
  d$fin_yield[27] <-  (d$fin_yield[26]+d$fin_yield[28])/2
  d$max_diff_yield[27] <-  (d$max_diff_yield[26]+d$max_diff_yield[28])/2
  d$inac_ts[27] <-  (d$inac_ts[26]+d$inac_ts[28])/2
}

d$at <- d$at*100
d$bb <- d$bb*100

bubo <- levels(as.factor(d$bb))
upth <- levels(as.factor(d$at))

control <- d[1,]
for (i in 2:length(bubo)) {
  zz <- d[1,]
  zz[4] <- bubo[i]
  control <- rbind(control, zz)
}

d <- rbind(control, d[-1,])

# resmat <- matrix(data = seq(1:(length(bubo)*length(bura))), ncol = length(bura), nrow = length(bubo))
resmat <- matrix(data = 100*d$fin_yield, ncol = length(upth), nrow = length(bubo))

fig <- plot_ly(
  x = bubo, 
  y = upth, 
  # z = matrix(data = d$ext_prob, ncol = length(bubo), nrow = length(bura)), 
  z = t(resmat),
  type = "contour",
  colorscale = list(c(0, 0.5, 1), c('red', 'orange', 'green')),
  autocontour = F,
  # contours = list(showlabels = TRUE),
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
  title = "Permissiveness (%)"#,
  # titlefont = f
)

fig <- fig %>% colorbar(title = "Final\nyield(%)")
fig <- fig %>% layout(xaxis = xlab, yaxis = ylab)
fig

#### time spend waiting ####

resmat <- matrix(data = 100*d$inac_ts, ncol = length(upth), nrow = length(bubo))

fig <- plot_ly(
  x = bubo, 
  y = upth, 
  # z = matrix(data = d$ext_prob, ncol = length(bubo), nrow = length(bura)), 
  z = t(resmat),
  type = "contour",
  #colorscale = list(c(0, 1), c('orange', 'green')),
  autocontour = F,
  # contours = list(showlabels = TRUE),
  contours = list(
    start = 0,
    end = 50,
    size = 10,
    showlabels = T
  )
)

xlab <- list(
  title = "Budget bonus (%)"#,
  # titlefont = f
)
ylab <- list(
  title = "Permissiveness (%)"#,
  # titlefont = f
)

fig <- fig %>% colorbar(title = "Waiting\n(%)")
fig <- fig %>% layout(xaxis = xlab, yaxis = ylab)
fig

######## plot trajectories in one parameter combo #########

multi.mean.trj <- function(df, upth, bubo, tmax, yaxis) {
  # subsetting
  dd <- subset(df, UT == upth & BB == bubo)
  
  # zz <- subset(df, UT == 0 & BB == 0 & ratio == bura)
  
  max.dd <- max(dd[10:dim(dd)[2]])+0.2*max(dd[10:dim(dd)[2]])
  # max.zz <- max(zz[10:dim(dd)[2]])+0.2*max(zz[10:dim(dd)[2]])
  
  # plot
  xrange <- seq(0, tmax, 1)
  plot.new()
  plot(1, type = "n", xlab = "Time", ylab = yaxis, ylim = c(0,max.dd), xlim = c(0,20)
       , xaxt = 'n', cex.axis = 2, cex.lab = 3) # , main = paste("UT = ", upth, " BB = ",bubo)
  # plot(1, type = "n", xlab = "time", ylab = yaxis, ylim = c(0,max(max.zz,max.dd)), xlim = c(0,20)) # , main = paste("UT = ", upth, " BB = ",bubo)
  axis(side=1, at=seq(0,length(xrange)-1,1), labels=xrange, cex.axis=2)
  
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
    upperUT <- rep((1+upth)*dd[1,6], length(xtendrange))
    lowerUT <- rep((1-upth)*dd[1,6], length(xtendrange))
    
    polygon(c(xtendrange,rev(xtendrange)),c(upperUT,rev(lowerUT)),col=trans, border = "green")
    abline(h = dd[1,6], col = "darkgreen", lty = 2, lwd = 5)
  }
  
  # loop over replicates
  for (i in which(dd$rep %% 10 == 0)) {
    points(x = xrange, y = dd[i,7:dim(dd)[2]], type = 'l', lwd = 5, col = 'grey')
  }
  
  moy <- mean(dd[,8])
  # infci <- NULL
  # supci <- NULL
  for (i in 9:(dim(dd)[2])) {
    moy <- c(moy, mean(dd[,i]))
  #   infci <- c(infci, boot_sd_ci(dd[,i])[2])
  #   supci <- c(supci, boot_sd_ci(dd[,i])[3])
  }
  
  # # confidence interval
  # arrows(xrange[-c(1,2)], infci, xrange[-c(1,2)], supci, length=0.03, angle=90, code=3, col = "black")
  
  # trajectory
  points(x = xrange[-1], y = moy, type = 'l', lwd = 5 + 2, col = "black")
}
 
# Example
popul <- read.csv(file = "~/Desktop/PhD/GitKraken/gmse_fork_RQ1/mem-noreset-results/mem-noreset-merged-results/pop-ATI-mem-noreset-merged.csv")[,-1]
library(stringr)
par(mar = c(5,6,2,2)+0.1)
multi.mean.trj(df = popul, upth = 0, bubo = 0, tmax = 20, yaxis = "Population size")
setwd(dir = "~/Desktop/PhD/GitKraken/gmse_fork_RQ1/mem-noreset-results/mem-noreset-merged-results/")
dev.print(device = png, file = "pop-CTL-bigLabs.png", width = 1000)

#### trajectories panel ####

# find max in 

multi.mean.pop.mod <- function(df, upth, bubo, tmax) {
  # subsetting
  dd <- subset(df, UT == upth & BB == bubo)
  
  ymax <- 4100
  lnwdt <- 4
  
  # plot
  xrange <- seq(0, tmax, 1)
  plot(1, type = "n", xlab = "", ylab = "", ylim = c(0,ymax), xlim = c(0,20)
       , xaxt = 'n', cex.axis = 2) # , cex.lab = 3, main = paste("UT = ", upth, " BB = ",bubo)
  # plot(1, type = "n", xlab = "time", ylab = yaxis, ylim = c(0,max(max.zz,max.dd)), xlim = c(0,20)) # , main = paste("UT = ", upth, " BB = ",bubo)
  axis(side=1, at=seq(0,length(xrange)-1,1), labels=xrange, cex.axis=2)

  xtendrange <- seq(-1,tmax+1,1)
  trans <- adjustcolor("lightgreen",alpha.f=0.5)
  upperUT <- rep((1+upth)*dd[1,6], length(xtendrange))
  lowerUT <- rep((1-upth)*dd[1,6], length(xtendrange))
    
  polygon(c(xtendrange,rev(xtendrange)),c(upperUT,rev(lowerUT)),col=trans, border = "green")
  abline(h = dd[1,6], col = "darkgreen", lty = 2, lwd = lnwdt)
  
  # loop over replicates
  for (i in which(dd$rep %% 10 == 0)) {
    points(x = xrange, y = dd[i,7:dim(dd)[2]], type = 'l', lwd = lnwdt, col = 'grey')
  }
  
  moy <- mean(dd[,8])
  # infci <- NULL
  # supci <- NULL
  for (i in 9:(dim(dd)[2])) {
    moy <- c(moy, mean(dd[,i]))
    #   infci <- c(infci, boot_sd_ci(dd[,i])[2])
    #   supci <- c(supci, boot_sd_ci(dd[,i])[3])
  }
  
  # # confidence interval
  # arrows(xrange[-c(1,2)], infci, xrange[-c(1,2)], supci, length=0.03, angle=90, code=3, col = "black")
  
  # trajectory
  points(x = xrange[-1], y = moy, type = 'l', lwd = lnwdt + 2, col = "black")
}

popul.ati <- read.csv("~/Desktop/PhD/GitKraken/gmse_fork_RQ1/noreset-results-save/noreset-merged-results/pop-ATI-noreset-merged.csv")[,-1]
popul.trj <- read.csv(file = "~/Desktop/PhD/GitKraken/gmse_fork_RQ1/mem-noreset-results/mem-noreset-merged-results/pop-ATI-mem-noreset-merged.csv")[,-1]

{
  plot.new()
  
  # divide into four boxes
  layout(matrix(c(1,2,3,0), nrow = 2, byrow = T))
  par(mar = c(3, 5, 1, 1))
  
  # set space for x and y titles
  par(oma = c(3, 2.5, 0, 0))
  
  # upper left: CTL
  multi.mean.pop.mod(df = popul.trj, upth = 0, bubo = 0, tmax = 20)
  mtext("CTL", side = 3, line = -2, adj = 0.98, cex = 2, col = "grey40")
  # upper right: ATI
  multi.mean.pop.mod(popul.ati, upth = 0.3, bubo = 0.1, tmax = 20)
  mtext("ATI", side = 3, line = -2, adj = 0.98, cex = 2, col = "grey40")
  # lower right: TRJ
  multi.mean.pop.mod(popul.trj, upth = 0.3, bubo = 0, tmax = 20)
  mtext("TRJ", side = 3, line = -2, adj = 0.98, cex = 2, col = "grey40")
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
dev.print(device = png, file = "pop-Strat-bigLabs.png", width = 1000)
