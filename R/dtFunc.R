# quick check whether data are loadable (i.e., radio buttons and persistence distribution are all fine)
DTparmLoadable<-function(parms){
  loadable<-T
  if (! (parms$persistence_distn %in% pdnames)) loadable<-F
  if ((parms$phifix!=1 & parms$phifix!=0) | (parms$ffix!=1 & parms$ffix!=0) | (parms$Ifix!=1 & parms$Ifix!=0)) loadable<-F
  if (!basicMode){
    if (! (parms$arrfun %in% c("Uniform", "Compound"))) loadable<-F
    if (sum(parms$arrcomponents %in% 0:1) != 3) loadable<-F
    if (loadable){# the radios and listboxes can be filled. Are the arrival parameters OK?
      val<-suppressWarnings(as.numeric(parms$arrstart))
      if (length(val) != 1 || nchar(val)==0 || is.na(val) || val<0 || val>=365 || abs(val-round(val))>0.001) loadable<-F
      lwr.u<-suppressWarnings(as.numeric(parms$lwr.u))
      if (length(lwr.u) != 1 || nchar(lwr.u)==0 || is.na(lwr.u) || lwr.u<0 || lwr.u>=364 || abs(lwr.u-round(lwr.u))>0.001) loadable<-F
      lwr.p1<-suppressWarnings(as.numeric(parms$lwr.p1))
      if (length(lwr.p1) != 1 || nchar(lwr.p1)==0 || is.na(lwr.p1) || lwr.p1<0 || lwr.p1>=364 || abs(lwr.p1-round(lwr.p1))>0.001) loadable<-F
      lwr.p2<-suppressWarnings(as.numeric(parms$lwr.p2))
      if (length(lwr.p2) != 1 || nchar(lwr.p2)==0 || is.na(lwr.p2) || lwr.p2<0 || lwr.p2>=364 || abs(lwr.p2-round(lwr.p2))>0.001) loadable<-F
    }
    if (loadable){
      upr.u<-suppressWarnings(as.numeric(parms$upr.u))
      if (length(upr.u) != 1 || nchar(upr.u)==0 || is.na(upr.u) || upr.u<=0 || upr.u>=364 || abs(upr.u-round(upr.u))>0.001) loadable<-F
      if (lwr.u >= upr.u) loadable<-F
      upr.p1<-suppressWarnings(as.numeric(parms$upr.p1))
      if (length(upr.p1) != 1 || nchar(upr.p1)==0 || is.na(upr.p1) || upr.p1<=0 || upr.p1>=364 || abs(upr.p1-round(upr.p1))>0.001) loadable<-F
      if (lwr.p1 >= upr.p1) loadable<-F
      upr.p2<-suppressWarnings(as.numeric(parms$upr.p2))
      if (length(upr.p2) != 1 || nchar(upr.p2)==0 || is.na(upr.p2) || upr.p2<=0 || upr.p2>=364 || abs(upr.p2-round(upr.p2))>0.001) loadable<-F
      if (lwr.p2 >= upr.p2) loadable<-F
    }
    if (loadable){
      wt.u<-suppressWarnings(as.numeric(parms$wt.u))
      if (length(wt.u) != 1 || nchar(wt.u)==0 || is.na(wt.u) || wt.u<0 || wt.u>1) loadable<-F
      wt.p1<-suppressWarnings(as.numeric(parms$wt.p1))
      if (length(wt.p1) != 1 || nchar(wt.p1)==0 || is.na(wt.p1) || wt.p1<0 || wt.p1>1) loadable<-F
      wt.p2<-suppressWarnings(as.numeric(parms$wt.p2))
      if (length(wt.p2) != 1 || nchar(wt.p2)==0 || is.na(wt.p2) || wt.p2<0 || wt.p2>1) loadable<-F
    }
    if (loadable){
      totwt<-wt.u+wt.p1+wt.p2
      parms$wt.u<-parms$wt.u/totwt
      parms$wt.p1<-parms$wt.p1/totwt
      parms$wt.p2<-parms$wt.p2/totwt
      a.p1<-suppressWarnings(as.numeric(parms$a.p1))
      if (length(a.p1) != 1 || nchar(a.p1)==0 || is.na(a.p1) || a.p1<=0) loadable<-F
      b.p1<-suppressWarnings(as.numeric(parms$b.p1))
      if (length(b.p1) != 1 || nchar(b.p1)==0 || is.na(b.p1) || b.p1<=0) loadable<-F
      a.p2<-suppressWarnings(as.numeric(parms$a.p2))
      if (length(a.p2) != 1 || nchar(a.p2)==0 || is.na(a.p2) || a.p2<=0) loadable<-F
      b.p2<-suppressWarnings(as.numeric(parms$b.p2))
      if (length(b.p2) != 1 || nchar(b.p2)==0 || is.na(b.p2) || b.p2<=0) loadable<-F
    }
  }
  if (!loadable){
    tkmessageBox(message="Error in data. Cannot load.",  icon="error", type="ok")
    return(F)
  }
  return(T)
}

drawdt<-function(dtdat){
  with(dtdat, {
  # parameters have been error-checked before the Calculate button is enabled, so just write tk's to dtPrevious
    alpha<-1-crlev
    grinc <- 120   # resolution of the graphs
    if (ffix==1){
      nse<-1
    } else {
      f<-seq(fmin, fmax, length = grinc)
      nse <- grinc
    }
    if (phifix==1){
      nphi <-1
    } else {
      phi <- seq(phimin, phimax, length = grinc)
      nphi <- grinc
    }
    if (Ifix == 0) {
        # then, there are a number of search schedules to consider. How many?
      if (ffix + phifix + Ifix == 0){ # then draw graphs for five representatives of search intervals (shortest, longest, middle, and two others)
        if (Imax-Imin>4){
         si<-unique(ceiling(((Imax/Imin)^.25)^(0:4)*Imin))    # this gives integer search intervals
         si[length(si)]<-Imax
         nsi<-length(si)
        } else {
          si<-round(seq(Imin,Imax,len=5),2)
          nsi<-5
        }
      } else {
        ni <- (span/Imin):(span/Imax)
        nsi <- length(ni)
        si <- span/ni
        sitest <- floor(si/Imin) * Imin
        sitest <- c(0, sitest)
        si <- si[sitest[2:nsi] - sitest[1:(nsi - 1)] > 0]
        nsi <- length(si)
        if (nsi > grinc) {
            si <- si[c(1, 1 + which(diff(round((1:nsi) * (grinc - 2)/nsi)) > 0), nsi)]
            nsi <- length(si)
        }
      }
    } else {
        nsi <- 1
        si<-Isam
    }
    garray <- array(dim = c(nse, nphi, nsi))
  ########## graphs and calculations
  # upon initiation of the design frame, determine the size of the default graphics window
    desw<-7*.Rvar$charSizeAdjust; desh<-8.5*.Rvar$charSizeAdjust # width and height of design graph windows
    des3sz<-c(desw,desh)
#    dev.new(width=desw,height=desh, noRStudioGD = T)
    dessz<-par('din'); dev.off()
    des3w<-11*.Rvar$charSizeAdjust; des3h<-7*.Rvar$charSizeAdjust;
    des3sz<-c(des3w,des3h)
#    dev.new(width=des3w,height=des3h, noRStudioGD = T); des3sz<-par('din'); dev.off()
    ### step 0: preliminaries
    ### step 1: define space for array of g values to plot
    if (!basicMode){
      if (arrfun=='Compound'){
        if (span > 365){
          tkmessageBox(message="
            Analysis is with compound arrival function is limited to one year.
              (1) Use a monitoring period of one year or less, or
              (2) Use uniform arrival function.
            ",icon='error')
          return(F)
        }
      syr<-as.numeric(format(as.Date(firstsearch),"%Y"))
      arryr<-ifelse(as.Date(firstsearch)-as.Date(paste0(syr,"-01-01")) < arrstart,syr-1,syr)
      arryr<-arryr
      s0<-as.numeric(as.Date(firstsearch)-as.Date(arrstart,origin=paste0(arryr,"-01-01")))
      FullYear<-paste0(format(as.Date(arrstart,origin=paste0(arryr,"-01-01")),"%d-%b-%Y"), " through ", format(as.Date(arrstart-1,origin=paste0(arryr+1,"-01-01")),"%d-%b-%Y"))
      FullYear<-FullYear
      # integral of arrival function from arrstart to monitoring start
      # monitoring ends at:
      sf<-s0+span
      if (sf > 365){ # then the monitoring season extends over one year past the start of arrivals...define new arrival date
        tkmessageBox(message="Monitoring season extends more than one year after arrivals begin.\nReset arrival function start date or monitoring schedule.",icon='error')
        return(F)
      }
      } else {
        s0<-0
        sf<-span
      }
      MonitoredPeriod<-paste0(format(as.Date(firstsearch),"%d-%b-%Y")," through ",format(as.Date(as.numeric(as.Date(firstsearch))+span,origin="1970-01-01"),"%d-%b-%Y"))
      MonitoredPeriod<-MonitoredPeriod
      if (arrfun == "Uniform"){
        arrmiss0<-0
        arrmissf<-0
      } else {
    ### fraction of carcasses arriving before the monitoring begins:
        # integral starts at 0
        # monitoring start is at:
        # integral of arrival function from arrstart to monitoring start: parts for each component
        arrmiss0<- # fraction of carcass that arrive before the monitoring season and are not accounted for in estimate for total arriving in monitored period
          ifelse(s0>lwr.u,arrcomponents[1]*wt.u*(s0-lwr.u)/(upr.u-lwr.u),0) +
          ifelse(s0>lwr.p1,arrcomponents[2]*wt.p1*pbeta((s0-lwr.p1)/(upr.p1-lwr.p1),a.p1,b.p1),0) +
          ifelse(s0>lwr.p2,arrcomponents[3]*wt.p2*pbeta((s0-lwr.p2)/(upr.p2-lwr.p2),a.p2,b.p2),0)
        arrmissf<- # fraction of carcasses arriving after monitoring ends:
          ifelse(sf<upr.u,arrcomponents[1]*wt.u*(upr.u-sf)/(upr.u-lwr.u),0) +
          ifelse(sf<upr.p1,arrcomponents[2]*wt.p1*(1-pbeta((sf-lwr.p1)/(upr.p1-lwr.p1),a.p1,b.p1)),0) +
          ifelse(sf<upr.p2,arrcomponents[3]*wt.p2*(1-pbeta((sf-lwr.p2)/(upr.p2-lwr.p2),a.p2,b.p2)),0)
      }
    }
    for (sii in 1:nsi) {
        #### calculate g for the search ith search schedule (assuming phi = 1)
        days <- seq(0, span, by = si[sii])
#        if (span/si[sii] - floor(span/si[sii]) > 0.8) days <- c(days, span)
        nsearch <- length(days)-1
        nt<-nsearch
        ind1 <- rep(1:nt, times = nt:1)
        ind2 <- ind1 + 1
        ind3 <- unlist(lapply(1:nt, function(x) x:nt)) + 1
        schedule <- cbind(days[ind1], days[ind2], days[ind3])
        schedule.index <- cbind(ind1, ind2, ind3)  #columns for arrival interval and search number
        nmiss <- schedule.index[, 3] - schedule.index[, 2]
        maxmiss <- max(nmiss)
        if (!basicMode){
          if (arrfun == "Uniform"){
            arrvec <- (schedule[, 2] - schedule[, 1])/span
          } else if (arrfun == "Compound"){
            # fraction of carcasses arriving in each search interval
            # integrate the arrival function by components
            arrvec0.u<-numeric(nt) # fraction of carcasses falling in each search interval
            # several different possibilities:
            # 1. search interval does not overlap uniform arrivals -> arrvec0.u == 0
            # 2. search interval entirely contained within uniform arrivals:
            ind<-which(s0+days[1:nsearch]>=lwr.u & s0+days[2:length(days)]<upr.u) # indices for which search intervals are entirely within uniform arrival period
            arrvec0.u[ind]<-(days[ind+1]-days[ind])/(upr.u-lwr.u)
            # 3. search interval straddles start of uniform arrivals (at most one interval)
            ind<-which(s0+days[1:nsearch]<=lwr.u & s0+days[2:length(days)]>lwr.u & s0+days[2:length(days)]<=upr.u)
            arrvec0.u[ind]<-(days[ind+1]+s0-lwr.u)/(upr.u-lwr.u)
            # 4. search interval straddles end of uniform arrivals
            ind<-which(s0+days[1:nsearch]<=upr.u & s0+days[2:length(days)]>upr.u)
            arrvec0.u[ind]<-(upr.u-(days[ind]+s0))/(upr.u-lwr.u)
            # 5. one search interval entirely covers uniform arrivals
            ind<-which(s0+days[1:nsearch]<=lwr.u & s0+days[2:length(days)]>upr.u)
            arrvec0.u[ind]<-1

            arrvec0.p1<-numeric(nsearch) # fraction of carcasses falling in each search interval
            # several different possibilities:
            # 1. search interval does not overlap p1 arrivals -> arrvec0.p1 == 0
            # 2. search interval entirely contained within p1 arrivals:
            ind<-which(s0+days[1:nsearch]>=lwr.p1 & s0+days[2:length(days)]<=upr.p1) # indices for which search intervals are entirely within p1 arrival period
            arrvec0.p1[ind]<-(pbeta((s0+days[ind+1]-lwr.p1)/(upr.p1-lwr.p1),a.p1,b.p1)-pbeta((s0+days[ind]-lwr.p1)/(upr.p1-lwr.p1),a.p1,b.p1))
            # 3. search interval straddles start of p1 arrivals (at most one interval)
            ind<-which(s0+days[1:nsearch]<=lwr.p1 & s0+days[2:length(days)]>lwr.p1 & s0+days[2:length(days)]<=upr.p1)
            arrvec0.p1[ind]<-pbeta((days[ind+1]+s0-lwr.p1)/(upr.p1-lwr.p1),a.p1,b.p1)
            # 4. search interval straddles end of p1 arrivals (at most one interval)
            ind<-which(s0+days[1:nsearch]<=upr.p1 & s0+days[2:length(days)]>upr.p1)
            arrvec0.p1[ind]<-1-pbeta(((days[ind]+s0)-lwr.p1)/(upr.p1-lwr.p1),a.p1,b.p1)
            # 5. one search interval entirely covers p1 arrivals
            ind<-which(s0+days[1:nsearch]<=lwr.p1 & s0+days[2:length(days)] > upr.p1)
            arrvec0.p1[ind]<-1

            arrvec0.p2<-numeric(nsearch) # fraction of carcasses falling in each search interval
            # several different possibilities:
            # 1. search interval does not overlap p2 arrivals -> arrvec0.p2 == 0
            # 2. search interval entirely contained within p2 arrivals:
            ind<-which(s0+days[1:nsearch]>=lwr.p2 & s0+days[2:length(days)]<=upr.p2) # indices for which search intervals are entirely within p2 arrival period
            arrvec0.p2[ind]<-(pbeta((s0+days[ind+1]-lwr.p2)/(upr.p2-lwr.p2),a.p2,b.p2)-pbeta((s0+days[ind]-lwr.p2)/(upr.p2-lwr.p2),a.p2,b.p2))
            # 3. search interval straddles start of p2 arrivals (at most one interval)
            ind<-which(s0+days[1:nsearch]<=lwr.p2 & s0+days[2:length(days)]>lwr.p2 & s0+days[2:length(days)]<=upr.p2)
            arrvec0.p2[ind]<-pbeta((days[ind+1]+s0-lwr.p2)/(upr.p2-lwr.p2),a.p2,b.p2)
            # 4. search interval straddles end of p2 arrivals (at most one interval)
            ind<-which(s0+days[1:nsearch]<=upr.p2 & s0+days[2:length(days)]>upr.p2)
            arrvec0.p2[ind]<-1-pbeta(((days[ind]+s0)-lwr.p2)/(upr.p2-lwr.p2),a.p2,b.p2)
            # 5. one search interval entirely covers p2 arrivals
            ind<-which(s0+days[1:nsearch]<=lwr.p2 & s0+days[2:length(days)] > upr.p2)
            arrvec0.p2[ind]<-1

            arrvec0<-arrvec0.u*wt.u+arrvec0.p1*wt.p1+arrvec0.p2*wt.p2
            arrvec0<-arrvec0/sum(arrvec0) # fraction of total carcasses arriving in each interval...scaled to just the monitored period
            # to extrapolate to the whole year, multiply the final prob_obs by arrmiss0 + arrmissf to incorporate the temporal coverage
            arrvec<-rep(arrvec0,nsearch:1)
          }
        } else {
          arrvec <- (schedule[, 2] - schedule[, 1])/span #basic mode parallels uniform (but adjusts by v upon final calculation of g)
          arrmiss0<-0; arrmissf<-1-dtdat$tarr
        }
        if (ffix == 1) {
            f0 <- f[1]
            powk <- cumprod(c(1, rep(k, maxmiss)))  # vector of k^i's
            notfind <- cumprod(1 - f0 * powk[-length(powk)])
            nvec <- c(1, notfind) * f0
            pfind.si <- nvec * powk  # conditional probability of finding a carcass on the ith search (row) after arrival for given (simulated) searcher efficiency (column)
            # persistences
            intxsearch <- unique(round(cbind(schedule[, 2] - schedule[, 1], schedule[, 3] - schedule[, 2]), 3), MAR = 1)
            ppersu <- ppersist(persistence_distn, t_arrive0 = 0, t_arrive1 = intxsearch[, 1], t_search = intxsearch[, 1] + intxsearch[, 2], pda = pda0, pdb = pdb0)
            prob_obs <- numeric(dim(schedule)[1])
            for (poi in 1:length(prob_obs)) {
                prob_obs[poi] <- pfind.si[nmiss[poi] + 1] * ppersu[which(abs(intxsearch[, 1] - (schedule[poi, 2] - schedule[poi, 1])) < 0.01 & abs(intxsearch[, 2] - (schedule[poi, 3] - schedule[poi, 2])) < 0.01)] *
                    arrvec[poi]
            }
            prob_obs <- sum(prob_obs) * phi # probability of observing a carcass that arrives in the monitored period
            garray[1, , sii] <- prob_obs * (1-(arrmiss0+arrmissf)) # probability of observing a carcass that arrives in the year
        } else {
            # searcher efficiencies
            if (maxmiss > 0) {
                powk <- cumprod(c(1, rep(k, maxmiss)))  # vector of k^i's
                notfind <- apply(1 - f %o% powk[-length(powk)], FUN = "cumprod", MARGIN = 1)
                if (maxmiss == 1)
                    nvec <- cbind(rep.int(1, grinc), notfind) * f else nvec <- cbind(rep.int(1, grinc), t(notfind)) * f
                pfind.si <- t(t(nvec) * powk)  # conditional probability of finding a carcass on the ith search (row) after arrival for given (simulated) searcher efficiency (column)
            } else {
                pfind.si <- f
            }
            intxsearch <- unique(round(cbind(schedule[, 2] - schedule[, 1], schedule[, 3] - schedule[, 2]), 3), MAR = 1)
            ppersu <- ppersist(persistence_distn, t_arrive0 = 0, t_arrive1 = intxsearch[, 1], t_search = intxsearch[, 1] + intxsearch[, 2], pda = pda0, pdb = pdb0)
            # arrivals if uniform arrivals
  #          if (arrfun == "Uniform")
  #              arrvec <- (schedule[, 2] - schedule[, 1])/span else if (arrfun == "Beta")
  #              arrvec <- pbeta(schedule[, 2]/span, shape1 = afa, shape2 = afb) - pbeta(schedule[, 1]/span, shape1 = afa, shape2 = afb)
            # add the probabilities
            prob_obs <- numeric(grinc)
            if (maxmiss > 0) {
                for (i in 1:dim(schedule)[1]) {
                    prob_obs <- prob_obs + pfind.si[, nmiss[i] + 1] * ppersu[which(abs(intxsearch[, 1] - (schedule[i, 2] - schedule[i, 1])) < 0.01 & abs(intxsearch[, 2] - (schedule[i, 3] - schedule[i, 2])) < 0.01),
                      ] * arrvec[i]
                }
            } else {
                for (i in 1:dim(schedule)[1]) {
                    prob_obs <- prob_obs + pfind.si[nmiss[i] + 1] * ppersu[which(abs(intxsearch[, 1] - (schedule[i, 2] - schedule[i, 1])) < 0.01 & abs(intxsearch[, 2] - (schedule[i, 3] - schedule[i, 2])) < 0.01),
                      ] * arrvec[i]
                }
            }
            garray[, , sii] <- outer(prob_obs, phi) * (1-(arrmiss0+arrmissf))
        }
    }
    newwin<-T
    new3win<-T
    if (ffix+phifix+Ifix == 0){
      newwin<-F
      if (length(dev.list()) > 0){
        for (i in dev.list()){
          dev.set(i)
          if (sum((par('din')-des3sz)^2==0)) {
            new3win<-F
            break
          }
        }
      }
    } else if (ffix+phifix+Ifix == 1 | ffix+phifix+Ifix == 2){
      new3win<-F
      if (length(dev.list()) > 0){
        for (i in dev.list()){
          dev.set(i)
          if (sum((par('din')-dessz)^2==0)) {
            newwin<-F
            break
          }
        }
      }
    } else {
      newwin<-F
      new3win<-F
    }

    if (newwin){
      if (.Rvar$platform == 'windows') windows.options(width = desw, height = desh)
      if (.Rvar$platform == 'mac') quartz.options(width = desw, height = desh)
      if (.Rvar$platform == 'linux') X11.options(width = desw, height = desh)
    }
    if (new3win){
      if (.Rvar$platform == 'windows') windows.options(width = des3w, height = des3h)
      if (.Rvar$platform == 'mac') quartz.options(width = des3w, height = des3h)
      if (.Rvar$platform == 'linux') X11.options(width = des3w, height = des3h)
    }
    dev.new(noRStudioGD = T)
#   image(f,phi,pMgtA,col=colorRampPalette(colr)(200),breaks=breaks,xlab="Searcher efficiency (p)\nP(observing carcass | carcass present)",ylab="Sampling Coverage (a)")
    if (dfoption == 'a'){
      coloroption<-'new'
      if (coloroption=='old'){
        colr<-c(colors()[c(490,552,498,652)],"#FFFFFFBF")
        breaks<-c(seq(0,0.05,len=200/5),seq(0.05,0.1,len=200/5),seq(0.1,.3,len=200/5),seq(0.3,.5,len=200/5),seq(0.5,1,len=200/5+1))  # colors for alpha
      } else {
        colr<-c("#FFFFFFBF", colors()[c(652, 498, 552, 490)])
        breaks<-c(seq(0, 0.01, len=200/5),seq(0.01,0.1,len=200/5),seq(0.1,.3,len=200/5),seq(0.3,.55,len=200/5),seq(0.55,1,len=200/5+1))
      }
    } else if (dfoption == 'g') {
      colr <- colors()[c(490, 125, 577, 653, 53, 552)]
      seqlen <- 200/length(colr)
      breaks <- c(seq(0, 0.1, len = seqlen), seq(0.1, 0.2, len = seqlen), seq(0.2, 0.35, len = seqlen), seq(0.35, 0.6, len = seqlen), seq(0.6, 0.8, len = seqlen), seq(0.8, 1, len = seqlen))
    }
    if (ffix + phifix + Ifix == 0){ # all parameters vary; draw graphs for five different search intervals
      par(mar=c(4,4,2,1),mgp=c(2,.7,0),mfrow=c(2,3))
      par(oma=c(0,0,3,0))
      xvals<-f; yvals<-phi
      xlbl<-"Searcher Efficiency (p)"; ylbl<-"Search Coverage (a)"
      for (sii in 1:nsi){
        gplo<-garray[,,sii]
        if (dfoption == 'a'){
          gplo <- pmax(round(gplo,4), 0.0001)
          gplo <- pmin(gplo, 0.9999)
          guni<-unique(as.vector(gplo))
          for (gi in 1:length(guni)){
            mmax<-fmmax(X,guni[gi])
            pM<-diff(sqrt(0:(mmax+1))); pM<-pM/sum(pM)
            pMgX<-dbinom(X, size = 0:mmax, prob = guni[gi])*pM; pMgX<-pMgX/sum(pMgX)
            gplo[gplo==guni[gi]]<-ifelse(mmax>=tau, 1-sum(pMgX[1:(round(tau)+1)]), 0)
          }
          lvl<-c(0.001,0.01,0.025,0.05,0.1,0.25,0.5,.75)
          lvl<-lvl[abs(lvl/alpha-1)>.01]
        } else {
          rng<-diff(range(gplo))
          if (rng<0.05) {
            lvl<-pretty(gplo,nlevels=10)
            lincr<-lvl[2]-lvl[1]
          } else if (rng<0.15) {
            lincr<-0.01
            lvl<-seq(round(min(gplo),2),round(max(gplo),2),by=lincr)
          } else if (rng<0.3) {
            lincr<-0.02
            lvl<-seq(round(min(gplo),2),round(max(gplo),2),by=lincr)
          } else if (rng<0.4) {
            lincr<-0.03
            bot<-seq(0.02,1,by=lincr)
            bot<-bot[min(which(bot>min(gplo)))]
            lvl<-seq(bot,round(max(gplo),2),by=lincr)
          } else if (rng<0.6) {
            lincr<-0.04
            bot<-seq(0.02,1,by=lincr)
            bot<-bot[min(which(bot>min(gplo)))]
            lvl<-seq(bot,round(max(gplo),2),by=lincr)
          } else if (rng<0.7) {
            lincr<-0.05
            bot<-seq(0.05,1,by=lincr)
            bot<-bot[min(which(bot>min(gplo)))]
            lvl<-seq(bot,round(max(gplo),2),by=lincr)
          } else {
            lincr<-0.1
            bot<-seq(0.1,1,by=lincr)
            bot<-bot[min(which(bot>min(gplo)))]
            lvl<-seq(bot,round(max(gplo),2),by=lincr)
          }
          lvl<-lvl[abs(lvl-target_g)>lincr/4]
        }
        image(xvals, yvals, gplo, col = colorRampPalette(colr)(length(breaks) - 1), breaks = breaks, xlab = xlbl, ylab = ylbl)
        if (dfoption =='a'){
          if (coloroption == 'old'){
            contour(xvals, yvals, gplo,levels=lvl,add=T,labcex=.Rvar$charSizeAdjust,col=ifelse(lvl>0.05, colors()[125], colors()[124]), lab = 1-lvl)
            contour(xvals, yvals, gplo,add=T,col=ifelse(alpha<.15, colors()[1], colors()[68]),levels=alpha,lwd=2,labcex=.9*.Rvar$charSizeAdjust, lab = 1-alpha)
          } else {
            contour(xvals, yvals, gplo,levels=lvl,add=T,labcex=.Rvar$charSizeAdjust,
#              col=ifelse(lvl<0.101,colors()[123],ifelse(lvl<.51, colors()[43], colors()[121])), lab = 1-lvl)
              col=ifelse(lvl<0.101,colors()[123], ifelse(lvl<0.51, colors()[43], colors()[124])), lab = 1-lvl)
#            contour(xvals, yvals, gplo,add=T,col=ifelse(alpha<.19, colors()[91], colors()[13]),levels=alpha,lwd=2,labcex=.9*.Rvar$charSizeAdjust, lab = 1-alpha)
            contour(xvals, yvals, gplo,add=T,col=ifelse(alpha < 1 - .899, 'black', colors()[13]),levels=alpha,lwd=2,labcex=.9*.Rvar$charSizeAdjust, lab = 1-alpha)

          }
        } else {
          contour(xvals, yvals, gplo, add = T, col = ifelse(lvl >= 0.19, colors()[125], colors()[123]), labcex = .8*.Rvar$charSizeAdjust, levels = lvl,method='flattest')
          contour(xvals,yvals,gplo,add=T,col=colors()[1],levels=target_g,lwd=2,labcex=.9*.Rvar$charSizeAdjust)
        }
        box()
        title(paste0("Search interval (I) = ",si[sii]))
        if(.Rvar$platform == "windows") bringToTop()
      }
      if (dfoption == 'a') tit0 <- 'Credibility level (1 - \u03b1)' else tit0 <- "Detection probability"
      mtext(paste0(tit0, " vs. searcher efficiency, coverage, and search interval"),side=3, line=1, outer=T, cex=1.5*.Rvar$charSizeAdjust)
#tit0 <- "Detection probability"
#mtext(paste0(tit0, " vs. searcher efficiency, coverage, and search interval"),side=3, line=1, outer=T, cex=1.25*.Rvar$charSizeAdjust,col=4)

      par(mar=c(0,0,0,0)) # space for summary of parameters
      plot(0,0,xlim=c(0,1),ylim=c(0,1),type='n',axes=F)
      lines(par('usr')[1:2],rep(.99,2),lty=3)
      par(family='serif')
      labat=0.01
      labcex<-.Rvar$charSizeAdjust
      if (dfoption == 'g') mtext(side=3, line=-4,  adj=0, at=labat, text=paste0("Target g = ", target_g), cex=labcex)
      if (dfoption == 'a') mtext(side=3, line=-4,  adj=0, at=labat, text=paste0("1 - \u03b1 = P(M <= ", tau, " | X = ", X, ")"),cex=labcex)
      mtext(side=3, line=-5.5, adj=0, at=labat, text=paste0("span = ", span),cex=labcex)
      mtext(side=3, line=-7, adj=0, at=labat, text=paste0("k = ", round(k,3)),cex=labcex)
      perslbl<-paste0(persistence_distn, " persistence")
      if (persistence_distn=="Exponential"){
        perslbl1<-paste0("   mean CP = ", round(pdb0,2), ", \u03bb = ",round(1/pdb0,4))
        perslbl2<-''
      } else {
        perslbl1<-paste0("   shape = ", round(pda0,4), ", scale = ", round(pdb0,4))
        CPr<-rCPgab(persistence_distn,pda0,pdb0,si[1])
        perslbl2<-paste0('   mean CP = ',round(CPr[1],2), " days")
        if (Ifix==1) perslbl2<-paste0(perslbl2,", r = ", round(CPr[2],3))
      }
      mtext(side=3, line=-8.5, adj=0, at=labat,text=perslbl,cex=labcex)
      mtext(side=3, line=-10, adj=0, at=labat,text=perslbl1,cex=labcex)
      mtext(side=3, line=-11.5, adj=0, at=labat,text=perslbl2,cex=labcex)
      arrlbl<-paste0(arrfun, " arrivals")
      startline=-13
      if (arrfun == "Compound"){
        arrlbl<-paste0(arrlbl, " with ", sum(arrcomponents), " components")
        arrsublbls<-character(3)
        if (arrcomponents[1]) arrsublbls[1]<-paste0('   uniform (', round(wt.u/sum(c(wt.u,wt.p1,wt.p2)*arrcomponents),3),'), ', format(as.Date(arrstart+lwr.u,origin="1970-01-01"),"%b %d"),
          ' to ', format(as.Date(arrstart+upr.u,origin="1970-01-01"),"%b %d"))
        if (arrcomponents[2]) arrsublbls[2]<-paste0('   pulse1 (',round(wt.p1/sum(c(wt.u,wt.p1,wt.p2)*arrcomponents),3),'), ',
          format(as.Date(arrstart+lwr.p1,origin="1970-01-01"),"%b %d"), ' to ',
          format(as.Date(arrstart+upr.p1,origin="1970-01-01"),"%b %d"), ', a = ',signif(a.p1,4),', b = ',signif(b.p1,4))
        if (arrcomponents[3]) arrsublbls[3]<-paste0('   pulse2 (',round(wt.p2/sum(c(wt.u,wt.p1,wt.p2)*arrcomponents),3),'), ',
          format(as.Date(arrstart+lwr.p2,origin="1970-01-01"),"%b %d"), ' to ',
          format(as.Date(arrstart+upr.p2,origin="1970-01-01"),"%b %d"), ', a = ',signif(a.p2,4),', b = ',signif(b.p2,4))
        mtext(side=3, line=startline, adj=0, at=labat,text=arrlbl,cex=labcex)
        for (i in 1:3) if (arrcomponents[i]) mtext(side=3, line=startline-1.5*sum(arrcomponents[1:i]), adj=0, at=labat,text=arrsublbls[i],cex=labcex)
        mtext(side=3,line=startline-1.5*sum(arrcomponents)-1.5,adj=0,at=labat, text=paste0("First search on ",firstsearch),cex=labcex)
        mtext(side=3,line=startline-1.5*sum(arrcomponents)-3,adj=0,at=labat, text=paste0("Temporal coverage = ",round(1-(arrmiss0+arrmissf),3)),cex=labcex)
      } else {
        mtext(side=3, line=startline, adj=0, at=labat,text=arrlbl,cex=labcex)
      }
      par(family='sans')
      if(.Rvar$platform == "windows") bringToTop()
    }
    if (ffix + phifix + Ifix == 1){
      par(mar=c(4,4,2,1),mgp=c(2,.7,0))
      if (Ifix == 1){
        xvals<-f; yvals<-phi; gplo<-garray[,,1]
        xlbl<-"Searcher Efficiency (p)"; ylbl<-"Search Coverage (a)"
        fixlbl<-paste0("Search interval: I = ",round(si,2),", span = ", round(span,2))
      } else if (phifix==1){
        xvals<-f; yvals<-si; gplo<-garray[,1,]
        xlbl<-"Searcher Efficiency (p)"; ylbl<-"Search Interval (I)\n(days)"
        fixlbl<-paste0("Search coverage: a = ",round(phi,3))
      } else if (ffix==1){
        xvals<-phi; yvals<-si; gplo<-garray[1,,]
        xlbl<-"Search Coverage (a)"; ylbl<-"Search Interval (I)\n(days)"
        fixlbl<-paste0("Searcher efficiency: p = ",round(f,3))
      }
      if (dfoption == 'a'){
        gplo <- pmax(round(gplo,4), 0.0001)
        gplo <- pmin(gplo, 0.9999)
        guni<-unique(as.vector(gplo))
        for (gi in 1:length(guni)){
          mmax<-fmmax(X,guni[gi])
          pM<-diff(sqrt(0:(mmax+1))); pM<-pM/sum(pM)
          pMgX<-dbinom(X, size = 0:mmax, prob = guni[gi])*pM; pMgX<-pMgX/sum(pMgX)
          gplo[gplo==guni[gi]]<-ifelse(mmax>=tau, 1-sum(pMgX[1:(round(tau)+1)]), 0)
        }
        lvl<-c(0.001,0.01,0.025,0.05,0.1,0.25,0.5,.75)
        lvl<-lvl[abs(lvl/alpha-1)>.01]
      } else {
        rng<-diff(range(gplo))
        if (rng<0.05) {
          lvl<-pretty(gplo,nlevels=10)
          lincr<-lvl[2]-lvl[1]
        } else if (rng<0.15) {
          lincr<-0.01
          lvl<-seq(round(min(gplo),2),round(max(gplo),2),by=lincr)
        } else if (rng<0.3) {
          lincr<-0.02
          lvl<-seq(round(min(gplo),2),round(max(gplo),2),by=lincr)
        } else if (rng<0.4) {
          lincr<-0.03
          bot<-seq(0.02,1,by=lincr)
          bot<-bot[min(which(bot>min(gplo)))]
          lvl<-seq(bot,round(max(gplo),2),by=lincr)
        } else if (rng<0.6) {
          lincr<-0.04
          bot<-seq(0.02,1,by=lincr)
          bot<-bot[min(which(bot>min(gplo)))]
          lvl<-seq(bot,round(max(gplo),2),by=lincr)
        } else if (rng<0.7) {
          lincr<-0.05
          bot<-seq(0.05,1,by=lincr)
          bot<-bot[min(which(bot>min(gplo)))]
          lvl<-seq(bot,round(max(gplo),2),by=lincr)
        } else {
          lincr<-0.1
          bot<-seq(0.1,1,by=lincr)
          bot<-bot[min(which(bot>min(gplo)))]
          lvl<-seq(bot,round(max(gplo),2),by=lincr)
        }
        lvl<-lvl[abs(lvl-target_g)>lincr/4]
      }
      par(fig=c(0,1,.2,1)) #space for figure
      image(xvals, yvals, gplo, col = colorRampPalette(colr)(length(breaks) - 1), breaks = breaks, xlab = xlbl, ylab = ylbl)
#      contour(xvals, yvals, gplo, add = T, col = ifelse(lvl >= 0.19, colors()[125], colors()[123]), labcex = .8, levels = lvl,method='flattest')
#      contour(xvals,yvals,gplo,add=T,col=colors()[1],levels=target_g,lwd=2,labcex=.9)
      if (dfoption =='a'){
        if (coloroption == 'old'){
          contour(xvals,yvals,gplo,levels=lvl,add=T,labcex=1*.Rvar$charSizeAdjust,col=ifelse(lvl>0.05, colors()[125], colors()[123]), lab = 1-lvl)
          contour(xvals,yvals,gplo,add=T,col=ifelse(alpha<.15, colors()[1], colors()[68]),levels=alpha,lwd=2,labcex=.9*.Rvar$charSizeAdjust, lab=1-alpha)
        } else {
          contour(xvals, yvals, gplo,levels=lvl,add=T,labcex=.Rvar$charSizeAdjust,
#              col=ifelse(lvl<0.101,colors()[123],ifelse(lvl<.51, colors()[43], colors()[121])), lab = 1-lvl)
#            col=ifelse(lvl<0.101,colors()[123], colors()[43]), lab = 1-lvl)
            col=ifelse(lvl<0.101,colors()[123], ifelse(lvl<0.51, colors()[43], colors()[124])), lab = 1-lvl)
#          contour(xvals, yvals, gplo,add=T,col=ifelse(alpha<.19, colors()[91], colors()[13]),levels=alpha,lwd=2,labcex=.9*.Rvar$charSizeAdjust, lab = 1-alpha)
          contour(xvals, yvals, gplo,add=T,col=ifelse(alpha < 1 - .899, 'black', colors()[13]),levels=alpha,lwd=2,labcex=.9*.Rvar$charSizeAdjust, lab = 1-alpha)
        }
      } else {
        contour(xvals, yvals, gplo, add = T, col = ifelse(lvl >= 0.19, colors()[125], colors()[124]), labcex = .8*.Rvar$charSizeAdjust, levels = lvl,method='flattest')
        contour(xvals,yvals,gplo,add=T,col=colors()[1],levels=target_g,lwd=2,labcex=.9*.Rvar$charSizeAdjust)
      }
      box()

      if (dfoption == 'g') title("Detection probability (g)") else title("Credibility level (1 - \u03b1)")
      par(mar=c(0,0,0,0),fig=c(0,1,0,.2),new=T) # space for summary of parameters
      plot(0,0,xlim=c(0,1),ylim=c(0,1),type='n',axes=F)
      lines(par('usr')[1:2],rep(.99,2),lty=3)
      # parameter values to write depend on which options are active
      par(family='serif')
      labat<-0.01
      labcex<-96/(par('cra')/par('cin'))[2]#.Rvar$charSizeAdjust
      startline=-1*labcex
      if (dfoption == 'g') mtext(side=3, line=labcex*-2,  adj=0, at=labat, text=paste0("Target g = ",target_g),cex=labcex)
      if (dfoption == 'a') mtext(side=3, line=labcex*-2,  adj=0, at=labat, text=paste0("1 - \u03b1 = P(M <= ",tau, " | X = ",X, ")"),cex=labcex)
      mtext(side=3, line=labcex*-3,  adj=0, at=labat, text=fixlbl,cex=labcex)
      mtext(side=3, line=labcex*-4, adj=0, at=labat, text=paste0("k = ", round(k,3)),cex=labcex)
      perslbl<-paste0(persistence_distn, " persistence")
      if (persistence_distn=="Exponential"){
        perslbl1<-paste0("   mean CP = ", round(pdb0,2), ", \u03bb = ",round(1/pdb0,4))
        perslbl2<-''
      } else {
        perslbl1<-paste0("   shape (\u03b1) = ", round(pda0,4), ", scale (\u03b2) = ", round(pdb0,4))
        CPr<-rCPgab(persistence_distn,pda0,pdb0,si[1])
        perslbl2<-paste0('   mean CP = ',round(CPr[1],2), " days")
        if (Ifix==1) perslbl2<-paste0(perslbl2,", r = ", round(CPr[2],3))
      }
      mtext(side=3, line=labcex*-5, adj=0, at=labat,text=perslbl,cex=labcex)
      mtext(side=3, line=labcex*-6, adj=0, at=labat,text=perslbl1,cex=labcex)
      mtext(side=3, line=labcex*-7, adj=0, at=labat,text=perslbl2,cex=labcex)
      arrlbl<-paste0(arrfun, " arrivals")
      if (arrfun == "Uniform") {
        mtext(side=3, line=labcex*-8, adj=0, at=labat+(!Ifix)*.2, text=arrlbl,cex=labcex)
      } else if (arrfun == "Compound") {
        arrlbl<-paste0(arrlbl,' with ', sum(arrcomponents),' components:')
        arrsublbls<-character(3)
        if (arrcomponents[1]) arrsublbls[1]<-paste0('   uniform (', round(wt.u/sum(c(wt.u,wt.p1,wt.p2)*arrcomponents),3),'), ', format(as.Date(arrstart+lwr.u,origin="1970-01-01"),"%b %d"),
          ' to ', format(as.Date(arrstart+upr.u,origin="1970-01-01"),"%b %d"))
        if (arrcomponents[2]) arrsublbls[2]<-paste0('   pulse1 (',round(wt.p1/sum(c(wt.u,wt.p1,wt.p2)*arrcomponents),3),'), ',
          format(as.Date(arrstart+lwr.p1,origin="1970-01-01"),"%b %d"), ' to ',
          format(as.Date(arrstart+upr.p1,origin="1970-01-01"),"%b %d"), ', a = ',signif(a.p1,4),', b = ',signif(b.p1,4))
        if (arrcomponents[3]) arrsublbls[3]<-paste0('   pulse2 (',round(wt.p2/sum(c(wt.u,wt.p1,wt.p2)*arrcomponents),3),'), ',
          format(as.Date(arrstart+lwr.p2,origin="1970-01-01"),"%b %d"), ' to ',
          format(as.Date(arrstart+upr.p2,origin="1970-01-01"),"%b %d"), ', a = ',signif(a.p2,4),', b = ',signif(b.p2,4))
        mtext(side=3, line=labcex*(startline-1), adj=0, at=0.5,text=arrlbl,cex=labcex)
        for (i in 1:3) if (arrcomponents[i]) mtext(side=3, line=labcex*(startline-1-sum(arrcomponents[1:i])), adj=0, at=0.5,text=arrsublbls[i],cex=labcex)
        mtext(side=3,line=labcex*(startline-1-sum(arrcomponents)-1),adj=0,at=0.5, text=paste0("First search on ",firstsearch, ", span = ", round(span,2)),cex=labcex)
        mtext(side=3,line=labcex*(startline-1-sum(arrcomponents)-2),adj=0,at=0.5, text=paste0("Temporal coverage = ",round(1-(arrmiss0+arrmissf),3)),cex=labcex)
      }
      if (Ifix==0 && arrfun!="Compound") mtext(side=3,line=-8*labcex, adj=0, at=labat,text=paste0("span = ",round(span,2)),cex=labcex)
      par(family='sans')
      if(.Rvar$platform == "windows") bringToTop()
    }
    if (ffix + phifix + Ifix == 2){
      par(mar=c(4,4,2,1),mgp=c(2,.7,0))
      labcex<-96/(par('cra')/par('cin'))[2]
      ylbl<-ifelse(dfoption == 'g', "Overall Detection Probability (g)", paste0("1 - \u03b1 = P(M <= ", tau, " | X = ", X, ")"))
      tit0<-ifelse(dfoption == 'g', "Overall detection probability", "Credibility level (1 - \u03b1)")
      if (Ifix == 0){ # I varies
        xvals<-si; gplo<-garray[1,1,]
        if (dfoption == 'a') {
          gplo <- pmax(round(gplo,4), 0.0001)
          gplo <- pmin(gplo, 0.9999)
          guni<-unique(as.vector(gplo))
          for (gi in 1:length(guni)){
            mmax<-fmmax(X,guni[gi])
            pM<-diff(sqrt(0:(mmax+1))); pM<-pM/sum(pM)
            pMgX<-dbinom(X, size = 0:mmax, prob = guni[gi])*pM; pMgX<-pMgX/sum(pMgX)
            gplo[gplo==guni[gi]]<-ifelse(mmax>=tau, 1-sum(pMgX[1:(round(tau)+1)]), 0)
          }
          gplo[is.na(gplo)]
        }
        xlbl<-"Search Interval (I)\n(days)"
        tit<-paste0(tit0, " vs. search interval")
        fixlbl1<- paste0("Coverage (a) = ", phi, ", searcher efficiency (p) = ",f)#coverage + searcher efficiency + span
        fixlbl2<-paste0("span = ", span)#coverage + searcher efficiency + span
      } else if (phifix==0){ #coverage varies
        xvals<-phi; gplo<-garray[1,,1]
        if (dfoption == 'a') {
          gplo <- pmax(round(gplo,4), 0.0001)
          gplo <- pmin(gplo, 0.9999)
          guni<-unique(as.vector(gplo))
          for (gi in 1:length(guni)){
            mmax<-fmmax(X,guni[gi])
            pM<-diff(sqrt(0:(mmax+1))); pM<-pM/sum(pM)
            pMgX<-dbinom(X, size = 0:mmax, prob = guni[gi])*pM; pMgX<-pMgX/sum(pMgX)
            gplo[gplo==guni[gi]]<-ifelse(mmax>=tau, 1-sum(pMgX[1:(round(tau)+1)]), 0)
          }
        }
        xlbl<-"Search Coverage (a)"
        tit<-paste0(tit0, " vs. search coverage")
        fixlbl1<- paste0("Searcher efficiency (p) = ", f)#searcher efficiency + search interval + span
        fixlbl2<-paste0("search interval = ", si, ", span = ", span)#coverage + searcher efficiency + span
      } else if (ffix==0){ #searcher efficiency varies
        xvals<-f; yvals<-si; gplo<-garray[,1,1]
        if (dfoption == 'a') {
          gplo <- pmax(round(gplo,4), 0.0001)
          gplo <- pmin(gplo, 0.9999)
          guni<-unique(as.vector(gplo))
          for (gi in 1:length(guni)){
            mmax<-fmmax(X,guni[gi])
            pM<-diff(sqrt(0:(mmax+1))); pM<-pM/sum(pM)
            pMgX<-dbinom(X, size = 0:mmax, prob = guni[gi])*pM; pMgX<-pMgX/sum(pMgX)
            gplo[gplo==guni[gi]]<-ifelse(mmax>=tau, 1-sum(pMgX[1:(round(tau)+1)]), 0)
          }
        }
        xlbl<-"Searcher Efficiency (p)"
        fixlbl1<- paste0("Search coverage (a) = ",phi)#coverage + search interval + span
        fixlbl2<-paste0("search interval = ", si, ", span = ", span)#coverage + searcher efficiency + span
        tit<-paste0(tit0, " vs. searcher efficiency")
      }
      par(fig=c(0,1,.2,1)) #space for figure
      typ<-ifelse(Ifix!=0,'l','n')
      if (dfoption == 'g'){
        plot(xvals,gplo, xlab='',ylab='',type=typ)
        title(xlab=xlbl,line=3); title(ylab=ylbl,line=2.5)
        if (max(gplo)>=target_g){
          if (Ifix == 0) {
            ind<-min(which(gplo>=target_g))
            frac<-(gplo[ind]-target_g)/(gplo[ind]-gplo[ind+1])
            xbd<-xvals[ind]-frac*(xvals[ind]-xvals[ind+1])
            if (ind>1) lines(xvals[1:ind],gplo[1:ind],col=2)
            lines(c(xvals[ind],xbd),c(gplo[ind],target_g),col=2)
            lines(c(xbd,xvals[ind+1]),c(target_g,gplo[ind+1]),col=1)
            lines(xvals[(ind+1):nsi],gplo[(ind+1):nsi])
            points(xbd,target_g,col=2,pch=16)
          } else {
            ind<-min(which(gplo>=target_g))
            lines(xvals[ind:length(gplo)],gplo[ind:length(gplo)],col=2)
            points(xvals[ind],gplo[ind],pch=16,col=2)
          }
        } else {
          if (Ifix == 0) lines(xvals,gplo)
        }
      } else {
        gplo<-1-gplo
        plot(xvals,gplo, xlab='',ylab='',type=typ)
        title(xlab=xlbl,line=3); title(ylab=ylbl,line=2.5)
        if (max(gplo) >= 1-alpha){
          if (Ifix == 0) {
            ind<-max(which(gplo >= 1-alpha))
            frac<-(gplo[ind]-(1-alpha))/(gplo[ind]-gplo[ind+1])
            xbd<-xvals[ind]-frac*(xvals[ind]-xvals[ind+1])
            if (ind>1) lines(xvals[1:ind],gplo[1:ind],col=2)
            lines(c(xvals[ind],xbd),c(gplo[ind],1-alpha),col=2)
            lines(c(xbd,xvals[ind+1]),c(1-alpha,gplo[ind+1]),col=1)
            lines(xvals[(ind+1):nsi],gplo[(ind+1):nsi])
            points(xbd,1-alpha,col=2,pch=16)
          } else {
            ind<-min(which(gplo >= 1-alpha))
            lines(xvals[ind:length(gplo)],gplo[ind:length(gplo)],col=2)
            points(xvals[ind],gplo[ind],pch=16,col=2)
          }
        } else {
          if (Ifix == 0) lines(xvals,gplo)
        }
      }
      title(tit)
      par(fig=c(0,1,0,.2),mar=c(0,0,0,0),new=T) #space for parameter list
      plot(0,0,type='n',xlim=c(0,1),ylim=c(0,1),axes=F)
      lines(par('usr')[1:2],rep(.99,2),lty=3)
      par(family='serif')
      startline<- -2

      if (dfoption == 'g'){
        glab<-paste0("Target g = ",target_g)
        if (max(gplo)<target_g){ # then other parameters need to be improved to attain target g
          if (Ifix==0){
            glab<-paste0(glab," cannot be attained with search interval >= ",si[1]," and the given parameters:")
          } else {
            glab<-paste0(glab," cannot be attained with ",ifelse(phifix==0, "coverage", "searcher efficiency")," <= ",xvals[1]," and the given parameters:")
          }
        } else { # g can be attained
          if (Ifix==1){
            if (ind > 1){# then there is a minimum value of xval that gives proper target g
              #linear interpolation of minimum value required
                frac<-(target_g-gplo[ind-1])/(gplo[ind]-gplo[ind-1])
                xbd<-xvals[ind-1]+frac*(xvals[ind]-xvals[ind-1])
                glab<-paste0(glab," attained when ",ifelse(phifix==0,"coverage","searcher efficiency")," exceeds ",round(xbd,3))
            } else {
              glab<-paste0(glab," attained with minimum ",ifelse(phifix==0,"coverage","searcher efficiency")," checked (",round(xvals[1],3),")")
            }
          } else { # minimum value in the list is sufficient for target g (and maybe a good deal smaller would be OK too)
            if (ind < length(si)){# then there is a minimum value of xval that gives proper target g
              #linear interpolation of minimum value required
                glab<-paste0(glab," attained when search interval <= ", round(xbd,1))
            } else {
              glab<-paste0(glab," attained with longest search interval checked (",round(xvals[1],1),")")
            }
          }
        }
      } else {
        glab<-paste0("Credibility level (1 - \u03b1) = ", 1-alpha)
        if (max(gplo)<1-alpha){ # then other parameters need to be improved to attain target g
          if (Ifix==0){
            glab<-paste0(glab," cannot be attained with search interval >= ",si[1]," and the given parameters:")
          } else {
            glab<-paste0(glab," cannot be attained with ",ifelse(phifix==0, "coverage", "searcher efficiency")," <= ",xvals[1]," and the given parameters:")
          }
        } else { # g can be attained
          if (Ifix==1){
            if (ind > 1){# then there is a minimum value of xval that gives proper target g
              #linear interpolation of minimum value required
                frac<-(1-alpha-gplo[ind-1])/(gplo[ind]-gplo[ind-1])
                xbd<-xvals[ind-1]+frac*(xvals[ind]-xvals[ind-1])
                glab<-paste0(glab," attained when ",ifelse(phifix==0,"coverage","searcher efficiency")," exceeds ",round(xbd,3))
            } else {
              glab<-paste0(glab," attained with minimum ",ifelse(phifix==0,"coverage","searcher efficiency")," checked (",round(xvals[1],3),")")
            }
          } else { # minimum value in the list is sufficient for target g (and maybe a good deal smaller would be OK too)
            if (ind < length(si)){# then there is a minimum value of xval that gives proper target g
              #linear interpolation of minimum value required
                glab<-paste0(glab," attained when search interval <= ", round(xbd,1))
            } else {
              glab<-paste0(glab," attained with longest search interval checked (",round(xvals[1],1),")")
            }
          }
        }
      }
      labat<-0.0
      mtext(side=3, line=startline*labcex,  adj=0, at=labat, text=glab,cex=labcex)
      mtext(side=3, line=labcex*(startline-1),  adj=0, at=labat, text=fixlbl1,cex=labcex)
      mtext(side=3, line=labcex*(startline-2),  adj=0, at=labat, text=fixlbl2,cex=labcex)
      mtext(side=3, line=labcex*(startline-3), adj=0, at=labat, text=paste0("k = ", round(k,3)),cex=labcex)
      perslbl<-paste0(persistence_distn, " persistence: ")
      if (persistence_distn=="Exponential"){
        perslbl1<-paste0("  mean CP = ", round(pdb0,2))
        if (Ifix==1){
          r<-rCPgab(persistence_distn,pda0,pdb0,si[1])[2]
        }
        perslbl1<-paste0(perslbl1, ", r = ",round(r, 3))
        perslbl2<-''
      } else {
        perslbl1<-paste0("  shape (\u03b1) = ", round(pda0,4), ", scale (\u03b2) = ", round(pdb0,4))
        CPr<-rCPgab(persistence_distn,pda0,pdb0,si[1])
        perslbl2<-paste0('  mean CP = ',round(CPr[1],2), " days")
        if (Ifix==1) perslbl2<-paste0(perslbl2,", r = ", round(CPr[2],3))
      }
      mtext(side=3, line=labcex*(startline-4), adj=0, at=labat,text=perslbl,cex=labcex)
      mtext(side=3, line=labcex*(startline-5), adj=0, at=labat,text=perslbl1,cex=labcex)
      mtext(side=3, line=labcex*(startline-6), adj=0, at=labat,text=perslbl2,cex=labcex)
      arrlbl<-paste0(arrfun, " arrivals")
      if (arrfun=="Compound") {
        arrlbl<-paste0(arrlbl,' with ', sum(arrcomponents),' components:')  # this needs to be filled in...
        arrsublbls<-character(3)
        if (arrcomponents[1]) arrsublbls[1]<-paste0('   uniform (', round(wt.u/sum(c(wt.u,wt.p1,wt.p2)*arrcomponents),3),'), ', format(as.Date(arrstart+lwr.u,origin="1970-01-01"),"%b %d"),
          ' to ', format(as.Date(arrstart+upr.u,origin="1970-01-01"),"%b %d"))
        if (arrcomponents[2]) arrsublbls[2]<-paste0('   pulse1 (',round(wt.p1/sum(c(wt.u,wt.p1,wt.p2)*arrcomponents),3),'), ',
          format(as.Date(arrstart+lwr.p1,origin="1970-01-01"),"%b %d"), ' to ',
          format(as.Date(arrstart+upr.p1,origin="1970-01-01"),"%b %d"), ', a = ',signif(a.p1,4),', b = ',signif(b.p1,4))
        if (arrcomponents[3]) arrsublbls[3]<-paste0('   pulse2 (',round(wt.p2/sum(c(wt.u,wt.p1,wt.p2)*arrcomponents),3),'), ',
          format(as.Date(arrstart+lwr.p2,origin="1970-01-01"),"%b %d"), ' to ',
          format(as.Date(arrstart+upr.p2,origin="1970-01-01"),"%b %d"), ', a = ',signif(a.p2,4),', b = ',signif(b.p2,4))
        mtext(side=3, line=labcex*(startline-1), adj=0, at=0.5,text=arrlbl,cex=labcex)
        for (i in 1:3) if (arrcomponents[i]) mtext(side=3, line=labcex*(startline-1-sum(arrcomponents[1:i])), adj=0, at=0.5,text=arrsublbls[i],cex=labcex)
        mtext(side=3,line=labcex*(startline-1-sum(arrcomponents)-1),adj=0,at=0.5, text=paste0("First search on ",firstsearch),cex=labcex)
        mtext(side=3,line=labcex*(startline-1-sum(arrcomponents)-2),adj=0,at=0.5, text=paste0("Temporal coverage = ",round(1-(arrmiss0+arrmissf),3)),cex=labcex)
      } else {
        mtext(side=3, line=labcex*(startline-1), adj=0, at=0.5,text=arrlbl,cex=labcex)
      }
      par(family='sans')
      if(.Rvar$platform == "windows") bringToTop()
    }
    # if all three parameters are fixed, then just repeat the parameters with a fixed g
    if (ffix + phifix + Ifix == 3){
      graphics.off()
      glab<-paste0("Overall Detection Probability: g =  ", round(garray[1,1,1],3))
      if (dfoption == 'a') {
        mmax<-fmmax(X,garray[1,1,1])
        pM<-diff(sqrt(0:(mmax+1))); pM<-pM/sum(pM)
        pMgX<-dbinom(X, size = 0:mmax, prob = garray[1,1,1])*pM; pMgX<-pMgX/sum(pMgX)
        crlev <- ifelse(mmax > tau, round(sum(pMgX[1:(round(tau)+1)]), 3), 1)
      }
      while(1){ if (sink.number()==0) break else sink() }
      sink(paste0(.Rvar$datadir,"/output"))
      cat(paste0(glab,'\n'))
      if (dfoption == 'a') cat(paste0("Credibility level: 1 - \u03b1 = P(M <= ", tau, " | X = ", X, ") = ", crlev,'\n'))
      cat(paste0("Coverage: a = ", phi, "\n"))
      cat(paste0("Searcher efficiency: p = ",f,"\n"))
      cat(paste0("Search schedule: First search = ", firstsearch,", I = ",si,", span = ",span,"\n"))
      if (arrfun=="Compound")    cat(paste0("   Temporal coverage: ", round(1-(arrmiss0+arrmissf),3),'\n'))
      cat(paste0("k = ", round(k,3),'\n'))
      perslbl1<-paste0(persistence_distn, " persistence")
      if (persistence_distn=="Exponential"){
        perslbl1<-paste0(perslbl1,"mean CP = ", round(pdb0,2))
        r<-rCPgab(persistence_distn,pda0,pdb0,si[1])[2]
        perslbl1<-paste0(perslbl1, ", r = ",round(r,3))
      } else {
        perslbl1<-paste0(perslbl1,"shape = ", round(pda0,4), ", scale = ", round(pdb0,4))
        CPr<-rCPgab(persistence_distn,pda0,pdb0,si[1])
        perslbl2<-paste0('    mean persistence time = ',round(CPr[1],2), " days")
        if (Ifix==1) perslbl2<-paste0(perslbl2,", r = ", round(CPr[2],3), " for Ir = ",si)
      }
      cat(paste0(perslbl1,'\n'))
      if (persistence_distn!="Exponential") cat(paste0(perslbl2,'\n'))
      arrlbl<-paste0(arrfun, " arrivals")
      if (arrfun=="Compound") {
        arrlbl<-paste0(arrlbl,' with ', sum(arrcomponents),' components:\n')  # this needs to be filled in...
        if (arrcomponents[1]) arrlbl<-paste0(arrlbl,'     uniform arrivals (', round(wt.u/sum(c(wt.u,wt.p1,wt.p2)*arrcomponents),3),') from ', format(as.Date(arrstart+lwr.u,origin="1970-01-01"),"%b %d"),
          ' to ', format(as.Date(arrstart+upr.u,origin="1970-01-01"),"%b %d"),'\n')
        if (arrcomponents[2]) arrlbl<-paste0(arrlbl, '     beta pulse (',round(wt.p1/sum(c(wt.u,wt.p1,wt.p2)*arrcomponents),3),') from ',
          format(as.Date(arrstart+lwr.p1,origin="1970-01-01"),"%b %d"), ' to ',
          format(as.Date(arrstart+upr.p1,origin="1970-01-01"),"%b %d"), ' with a = ',signif(a.p1,4),' and b = ',signif(b.p1,4),'\n')
        if (arrcomponents[3]) arrlbl<-paste0(arrlbl, '     beta pulse (',round(wt.p2/sum(c(wt.u,wt.p1,wt.p2)*arrcomponents),3),') from ',
          format(as.Date(arrstart+lwr.p2,origin="1970-01-01"),"%b %d"), ' to ',
          format(as.Date(arrstart+upr.p2,origin="1970-01-01"),"%b %d"), ' with a = ',signif(a.p2,4),' and b = ',signif(b.p2,4),'\n')
      }
      cat(arrlbl)
      sink()
      file.show(paste0(.Rvar$datadir,"/output"),delete.file=T,title="Estimated detection probability (g)")
    }
  })
}
dtsave<-function(dt){
  save(dt, file = paste0(.Rvar$datadir,"/dtPrevious.Rdata"))
}

dtmvg<-function(X, crlev){
  .Rvar$viewMvg <- tktoplevel()
  tktitle(.Rvar$viewMvg) <- "M* vs. g"
  tkgrab.set(.Rvar$viewMvg);  tkfocus(.Rvar$viewMvg)
  tkconfigure(.Rvar$viewMvg,width=1000)
  tkwm.resizable(.Rvar$viewMvg,0,0)
  tkwm.deiconify(.Rvar$viewMvg)
  gBounds <- tkframe(.Rvar$viewMvg)
  barwidth<-5
  wbord<-1
  scwid<-900
  relief<-'flat'
  mvgdatFrame<-tkframe(.Rvar$viewMvg)
  tkX<-tclVar(0)
  tkcrlev<-tclVar(.8)
  X.edit<-tkentry(mvgdatFrame, textvariable = tkX, width = 4, bg = 'white', justify = 'right')
  X.lbl<-tklabel(mvgdatFrame, text='X')
  crlev.edit<-tkentry(mvgdatFrame, textvariable= tkcrlev, width =4, bg='white', justify = 'right')
  crlev.lbl<-tklabel(mvgdatFrame, text='Credibility level (1 - \u03b1)')
  tkgrid(X.lbl, X.edit, crlev.lbl, crlev.edit, padx = 5)
  tkgrid(mvgdatFrame, sticky = 'w', pady = c(15,0), padx=30)
  tkXok <- tclVar("T")
  tkaok <- tclVar("T")
  .Rvar$dogo<-T
  tkgmin<-tclVar(0.01)
  tkgmax<-tclVar(0.95)
  plotmvg <- function(X, crlev) { # not currently scaled to endpoints
    if (.Rvar$dogo){
      ng<-200
      g<-seq(toR(tkgmin), toR(tkgmax), length = ng)
      mmax<-Vectorize(fmmax, vectorize.args = "g")(x = X, g)
      Mstar<-numeric(ng)
      for (gi in 1:ng){
        pM<-diff(sqrt(0:(mmax[gi]+1)))
        pXgM<-dbinom(X, 0:mmax[gi], g[gi])
        pMgX<-pXgM*pM/sum(pXgM*pM)
        Mstar[gi]<-min(which(cumsum(pMgX) >= crlev)) - 1
      }
      par(mar=c(4,4,1,1), mgp = c(2.5,.8, 0))
      plot(g,Mstar, type='s', ylab = 'M*', lwd = 2, col = colors()[73],xaxs = 'i', yaxs = 'i', ylim = c(0,max(Mstar)*1.05))
    } else {
      par(mar=c(4,4,1,1), mgp = c(2.5,.8, 0))
      plot(0,0,type='n', ylab = 'M*', lwd = 2, col = colors()[73],xaxs = 'i', yaxs = 'i', ylim = c(0,100), xlim=c(0,1))
    }
  }
  mvgfig <- tkrplot::tkrplot(.Rvar$viewMvg, fun = function() plotmvg(X = toR(tkX), crlev = toR(tkcrlev)), hscale=2.5, vscale=1.4)
  onChange <- function(...) {
    tkrplot::tkrreplot(mvgfig)
  }
  gminSlider <- tkscale(gBounds, from = 0.005, to = 1, variable = tkgmin, orient = "horizontal",
    length = scwid,
    width = barwidth,
    command = function(...) {
      if (toR(tkgmin) >= toR(tkgmax)){
        tclvalue(tkgmin)<-toR(tkgmax)-0.001
      }
      onChange()
    },
    resolution = .003,
    sliderrelief=relief,
    sliderlength=8,
    troughcolor=colors()[125],
    showvalue=F,
    borderwidth=wbord
  )
  gmaxSlider <- tkscale(gBounds, from = 0.005, to = 1, variable = tkgmax, orient = "horizontal",
    length = scwid,
    width = barwidth,
    command = function(...) {
      if (toR(tkgmin) >= toR(tkgmax)){
        tclvalue(tkgmax)<-toR(tkgmin)+0.001
      }
      onChange()
    },
    resolution = .003,
    sliderrelief=relief,
    sliderlength=8,
    troughcolor=colors()[624],
    showvalue=F,
    borderwidth=wbord
  )
  tkgrid(mvgfig)
  tkgrid(gminSlider)
  tkgrid(gmaxSlider)
  tkgrid(gBounds,pady=15)
  tkbind(X.edit, "<KeyRelease>", function(){
    X<-suppressWarnings(as.numeric(toR(tkX)))
    if (is.na(X) || X < 0 || X != round(X)) {
      tclvalue(tkXok) <- "F"
      tkconfigure(X.edit, bg = 'yellow')
    } else {
      tclvalue(tkXok) <- "T"
      tkconfigure(X.edit, bg = 'white')
    }
    if (tclvalue(tkXok) == "T" & tclvalue(tkaok) == "T") {
      .Rvar$dogo <- T
    } else {
      .Rvar$dogo <- F
    }
    tkrplot::tkrreplot(mvgfig)
  })
  tkbind(crlev.edit, "<KeyRelease>", function(){
    crlev<-suppressWarnings(as.numeric(toR(tkcrlev)))
    if (is.na(crlev) || crlev <= 0 || crlev >= 1) {
      tclvalue(tkaok) <- "F"
      tkconfigure(crlev.edit, bg = 'yellow')
    } else {
      tclvalue(tkaok) <- "T"
      tkconfigure(crlev.edit, bg = 'white')
    }
    if (tclvalue(tkXok) == "T" & tclvalue(tkaok) == "T") {
      .Rvar$dogo <- T
    } else {
      .Rvar$dogo <- F
    }
    tkrplot::tkrreplot(mvgfig)
  })
}