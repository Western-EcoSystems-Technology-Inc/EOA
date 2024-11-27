trig<-function(x, Eg, Vg, pBa, pBb, Tau, lta){
  mmax<-ifelse (Vg < 0.00001, fmmax(x,Eg),fmmax.ab(x,pBa,pBb))
  M<-x:mmax
 # posterior for M | x
  if (Vg > 0.00001){
    pXgM<-VGAM::dbetabinom.ab(x,size=M,shape1=pBa,shape2=pBb) # the probabilities of X for M = 0:mmax
  } else {
    pXgM<-dbinom(x,size=M,prob=Eg) # the probabilities of X for M = x:mmax
  }
  pM<-diff((sqrt(x:(mmax+1))))
  pMgX<-pXgM*pM; pMgX<-pMgX/sum(pMgX) # posterior distribution for M (ignoring M < x, which has probability = zero)
  mstar<-(x:mmax)[min(which(cumsum(pMgX)>1-lta))]
  return(list(mstar = mstar, trigFired = ifelse( mstar > Tau, T, F)))
}

scexCalc<-function(scexdat, viewone = 0, scex1frm = NA){
  if (viewone) {
    scexoutTable <- scex1frm$scexoutTable
    scexoutData <- scex1frm$scexoutData
  }
  Mstarx<-NA
  dofast<-NA
  with(scexdat, {
    gcum<-Xmov<-Mmov<-gmov<-numeric(nyr)+NA
    res<-array(dim=c(nyr,13,nsim))
    ## columns in the results array correspond to the following variables:
    #1. M = number of fatalities
    #2. gi = indicator of monitoring protocol with 1 => intensive and 2 => non-intensive (indexing directly into g array)
    #3. x
    #4. Mstar
    #5. lams
    #6. lamr
    #7. AMAi
    #8. stT
    #9. ltT
    #10. rtT
    #11. rho
    #12. pBa (cumulative, not short-term)
    #13. pBb
    dofast<<-!(iAMA & stT) & !rT & !viewone & (gpost==0 | g3i==2)
    if (dofast){ # then do exact calculation (fast)
      res[,7,]<-0     # baseline AMA
      res[,7:10,]<-0 # triggers are presumed not to have fired until they fire
      res[,11,]<-1
      g<-  cbind(c(g1, g2, 0),c(g1lwr, g2lwr, 0),c(g1upr, g2upr, 0))
      if (gpost){
        g[3,]<-g[g3i,]
      }
      res[,2,]<-c(rep(1,yrs),rep(2,nyr-yrs)) # index into g array
      mu<-maxp<-minp<-Vg<-numeric(nyr)
      mu[1:yrs]<-g[1,1] # tmp[1:yi,2] is index into rows of g
      minp[1:yrs]<-g[1,2]
      maxp[1:yrs]<-g[1,3]
      if (nyr>yrs){
        mu[(yrs+1):nyr]<-g[2,1] # tmp[1:yi,2] is index into rows of g
        minp[(yrs+1):nyr]<-g[2,2]
        maxp[(yrs+1):nyr]<-g[2,3]
      }
      sig2<-((maxp-minp)/4)^2 # estimated means and uncertainties for g
      Eg<-cumsum(mu)/(1:nyr) # average, expected G is the weighted average of the observed g's
      gcum<-Eg
      Vg<-cumsum(sig2)/(1:nyr)^2
      pBa<-Eg^2/Vg*(1-Eg)-Eg; pBb<-pBa*(1/Eg-1) # these are the two shape parameters for the beta distribution underlying the beta-binomial for X | M
      res[,12,]<-rep(pBa,nsim)
      res[,13,]<-rep(pBb,nsim)
      minxTrig<-numeric(nyr)
      yi<-1
      # find minimum x that triggers for yi = 1
      if (ltT){
        for (x in 0:(Tau+1)){
          if (trig(x, Eg[yi], Vg[yi], pBa[yi], pBb[yi], Tau, lta)$trigFired){
            minxTrig[yi]<-x
            break
          }
        }
        dev.new(width = 7, height = 0.75)
        par(mar=c(0,0,2,0))
        plot(0,0,type='n',xlim=c(0,1),xaxs='i',ylim=c(0,1),yaxs='i',axes=F,xlab='years')
        title('Calculating...',adj=0)
        points(seq(1/nsim,1,length=100),rep(.5,100),pch='.')
        colr<-c(4,3,7,2)
        if(.Rvar$platform == "windows") bringToTop()
        polygon(1/nyr+c(-1,-1,0,0)/nyr,c(0,1,1,0),col=colorRampPalette(colr)( nyr )[1],border=NA) # showing progress of calculations
        for (yi in 2:nyr){
  #       if (minxTrig[yi-1]) triggers, then try a smaller x
          polygon(yi/nyr+c(-1,-1,0,0)/nyr,c(0,1,1,0),col=colorRampPalette(colr)( nyr )[yi],border=NA) # showing progress of calculations
          x <- minxTrig[yi-1]
          # check whether x[yi-1] triggers for yi
          if (trig(x, Eg[yi], Vg[yi], pBa[yi], pBb[yi], Tau, lta)$trigFired){
            # need to find a smaller x
            for (x in (minxTrig[yi-1]-1):0){
              if (x==0){
                minxTrig[yi]<-0
              } else {
                if (!trig(x, Eg[yi], Vg[yi], pBa[yi], pBb[yi], Tau, lta)$trigFired){
                  minxTrig[yi]<-x+1
                  break
                }
              }
            }
          } else {
            # need to find a larger x
            for (x in minxTrig[yi-1]:Tau+1){
              if (trig(x, Eg[yi], Vg[yi], pBa[yi], pBb[yi], Tau, lta)$trigFired){
                minxTrig[yi]<-x
                break
              }
            }
          }
        }
      }
      # simulate:
      M<-rpois(nsim*nyr, lambda)
      res[,1,]<-M
      res[,3,]<-rbinom(nsim*nyr, size = M, prob = c(rep(g1, yrs), rep(g2, nyr-yrs)))
      Xcum<-matrixStats::colCumsums(res[,3,])
      Mstarx<<-list()
      if (ltT){
        polygon(c(0,0,nyr,nyr),c(0,1,1,0),col='white',border=NA) # showing progress of calculations
        for (yi in 1:nyr){
          Mstarx[[yi]]<-numeric(max(Xcum[yi,])+1)
          for (x in 1:length(Mstarx[[yi]])-1){
            Mstarx[[yi]][x+1]<-trig(x, Eg[yi], Vg[yi], pBa[yi], pBb[yi], Tau, lta)$mstar
          }
          polygon(yi/nyr+c(-1,-1,0,0)/nyr,c(0,1,1,0),col=3,border=NA) # showing progress of calculations
        }
        dev.off()
      } else {
        for (yi in 1:nyr){
          Mstarx[[yi]]<-numeric(max(Xcum[yi,])+1)
          for (x in 1:length(Mstarx[[yi]])-1){
            Mstarx[[yi]][x+1]<-trig(x, Eg[yi], Vg[yi], pBa[yi], pBb[yi], Tau, lta)$mstar
          }
        }
      }
      Mstar<-array(numeric(nsim*2),dim=c(nsim,2)) # Mstar at triggering and at end of project
      if (ltT){
        avoid <- (Xcum>=minxTrig)*(1:nyr)
        avoid[avoid==0] <- nyr + 1
        avoid<-matrixStats::colMins(avoid)
        res[,9,][((1:nsim)-1)*nyr+pmin(avoid,nyr)]<-1
        res[nyr,9,avoid==nyr+1]<-0
        curt<-sapply(avoid, function(x) {x<-min(x,nyr); c(numeric(x), 1+numeric(nyr-x))})
        res[,1,]<-res[,1,]-curt*res[,1,]+rpois(length(res[,1,]), lambda = lambda * rhoinf * curt)

        if (!gpost){ # no monitoring after trigger, so x's are zero after avoidance
          res[,3,]<-res[,3,]*(1-curt)
          Xcum<-matrixStats::colCumsums(res[,3,])
        } else if (g3i==2){ # monitoring is at non-intensive level following avoidance
          res[,3,]<-res[,3,]*(1-curt)+rbinom(nsim*nyr, size=res[,1,]*curt, prob=g[2,1])
        }
        yind<-pmin(avoid,nyr)
        Mstar[,1]<-sapply(1:nsim, function(x) Mstarx[[yind[x]]][Xcum[yind[x],x]+1]) # Mstar for the cumulative X at the triggering year
      } else {
        res[,9,][((1:nsim)-1)*nyr+nyr]<-1
        Mstar[,1]<-sapply(1:nsim, function(x) Mstarx[[nyr]][Xcum[nyr,x]+1]) # Mstar for the cumulative X at the triggering year
      }
      Mstar[,2]<-Mstarx[[nyr]][Xcum[nyr,]+1]
    } else { # rT, iAMA, or g3i so do brute force calculation
      if (stT & iAMA) steps <- dim(scexdat$iAMAschedule)[1]
    # data are presumed to have already been error-checked and fed into R
        g<-  cbind(c(g1, g2, 0),c(g1lwr, g2lwr, 0),c(g1upr, g2upr, 0))
        if (gpost){
          g[3,]<-g[g3i,]
        }
        if (viewone==0){
          tmp<-array(dim=c(nyr,13))
          graphics.off()
          dev.new(width =7,height=.75) # figure for watching progress of calculations
          par(mar=c(0,0,2,0))
          plot(0,0,type='n',xlim=c(0,1),xaxs='i',ylim=c(0,1),yaxs='i',axes=F,xlab='years')
          title('Calculating...',adj=0)
          points(seq(1/nsim,1,length=100),rep(.5,100),pch='.')
          colr<-c(4,3,7,2)
          if(.Rvar$platform == "windows") bringToTop()
        } else {
          nsim<-1
        }
      #  numeric(nyr)
        for (simi in 1:nsim){
          tmp<-array(dim=c(nyr,13))    # tmp is an array to store results for one simulation run; attached to full res array at end of simulation draw
          tmp[1:yrs,2]<-1 # intensive monitoring in the first yrs of the project
          if (nyr < yrs) tmp[(yrs+1),2]<-2 # non-intensive monitoring in the last nyr-yrs of the project
          tmp[1,7]<-0     # baseline AMA
          reverted<-F   # indicator whether reversion has occurred
          noncompliance<-F # indicator whether ITP has been exceeded
          tmp[,7:10]<-0 # triggers are presumed not to have fired until they fire
          tmp[,11]<-1
          for (yi in 1:nyr){
          ### determine rate
            lam<-lambda*tmp[yi,11] # rate was determined at the end of last year's monitoring (with default = 1 but changes with triggering)
          ### determine detection probability
            if (!is.na(tmp[yi,2])){ # then the monitoring rate has been determined by earlier short-term triggering or because project is still in initial yrs
              gi<-tmp[yi,2]
            } else {
              gi<-2 # after the initial years of intensive monitoring, the default is non-intensive monitoring (unless non-compliance)
              tmp[yi,2]<-gi
            }
          ### generate carcasses
            tmp[yi,1]<-rpois(1,lam)
          ### count carcasses
            tmp[yi,3]<-ifelse(tmp[yi,1]==0, 0, rbinom(1,tmp[yi,1],g[tmp[yi,2],1]))
          ### run long-term test
          # calculate cumulative g
            # create an array to derive pBa and pBb from
            mu<-g[tmp[1:yi,2],1] # tmp[1:yi,2] is index into rows of g
            minp<-g[tmp[1:yi,2],2]; maxp<-g[tmp[1:yi,2],3]; a<-tmp[1:yi,11]/sum(tmp[1:yi,11])
            x<-sum(tmp[1:yi,3])
            sig2<-((maxp-minp)/4)^2 # estimated means and uncertainties for g
            Eg<-sum(mu*a)/sum(a) # average, expected G is the weighted average of the observed g's
            gcum[yi]<-Eg
          # Lest<-max(x, 0.5)/Eg
            Vg<-sum(sig2*a^2)#+(sum(a*(mu^2+sig2))-Eg^2)/Lest
            pBa<-Eg^2/Vg*(1-Eg)-Eg; pBb<-pBa*(1/Eg-1) # these are the two shape parameters for the beta distribution underlying the beta-binomial for X | M
            tmp[yi,12]<-pBa; tmp[yi,13]<-pBb
          # M needs to be calculated up to a reasonable maximum, but the final distribution is clipped at P(M > m) < 0.0001
            mmax<-NA # stotmp the variable for use outside the loop
            mmax<-ifelse (Vg < 0.00001, fmmax(x,Eg),fmmax.ab(x,pBa,pBb))
            M<-x:mmax
           # posterior for M | x
            if (Vg > 0.00001){
              pBa<-Eg^2/Vg*(1-Eg)-Eg; pBb<-pBa*(1/Eg-1) # these are the two shape parameters for the beta distribution underlying the beta-binomial for X | M
              pXgM<-VGAM::dbetabinom.ab(x,size=M,shape1=pBa,shape2=pBb) # the probabilities of X for M = 0:mmax
            } else {
              pXgM<-dbinom(x,size=M,prob=Eg) # the probabilities of X for M = x:mmax
            }
            pM<-diff((sqrt(c(x-1,x:mmax)+1)))
            pMgX<-pXgM*pM; pMgX<-pMgX/sum(pMgX) # posterior distribution for M (ignoring M < x, which has probability = zero)

            tmp[yi,4]<-M[min(which(1-cumsum(pMgX)<lta))]
            if (ltT){
              if (tmp[yi,4]>Tau){ # then long-term trigger is fired
                noncompliance<-T
                tmp[yi,9]<-1
                if (yi < nyr) tmp[yi+1,7]<-ifelse(iAMA & stT, steps, 1)
              }
            }
            # test short-term trigger (and, if noncompliance already, then break afterwards)
            if (stT){
              minyi<-max(1,yi-sty+1)
              x3<-sum(tmp[minyi:yi,3])
              Xmov[yi]<-x3; Mmov[yi]<-sum(tmp[minyi:yi,1])
              mu<-g[tmp[minyi:yi,2],1] # where gi
              minp<-g[tmp[minyi:yi,2],2]; maxp<-g[tmp[minyi:yi,2],3];
              if (sum(tmp[minyi:yi,11]) > 0) a<-tmp[minyi:yi,11]/sum(tmp[minyi:yi,11]) else a<-1
              sig2<-((maxp-minp)/4)^2 # estimated means and uncertainties for g
              Eg3<-sum(mu*a)/sum(a) # average, expected G is the weighted average of the observed g's
              gmov[yi]<-Eg3
              Lest<-max(x3/3, 0.5)/Eg3
              Vg3<-sum(sig2*a^2)+(sum(a*(mu^2+sig2))-Eg3^2)/Lest
              if (Vg3<0.000001){ #then g is essentially constant, so use binomial rather than betabinomial
                tmp[yi,5]<-pLgX(ifelse(stTon=='same',Tau/nyr, tau)*sty,x3,Eg3,fmmax(x3,Eg3))
              } else { #use beta distribution (after calculating betaa and betab--don't use pBa and pBb here because those names are reserved for cumulative rather than short-term
                betaa<-Eg3^2/Vg3*(1-Eg3)-Eg3; betab<-betaa*(1/Eg3-1) # these are the two shape parameters for the beta distribution underlying the beta-binomial for X | M
                tmp[yi,5]<-pLgX.ab(ifelse(stTon=='same',Tau/nyr, tau)*sty,x3,betaa,betab,fmmax.ab(x3,betaa,betab))
              }
              if (tmp[yi,5]<sta){ # then the short-term trigger is fired
                tmp[yi,8]<-1   # marker for stT firing
                if (iAMA){
                  if (yi<nyr){ # then make adjustments for future years
                    if (reverted) { # reversion trigger is active, so reset rho to 1, steps to 0, and g indicator to first step
                      if (!noncompliance | gpost){ #monitoring to continue in future
                        tmp[(yi+1):nyr,7]<-0  # step is set to 0
                        tmp[(yi+1):nyr,11]<-1 # rho is set to 1
                        tmp[(yi+1):nyr, 2]<-2 # sampling is set to non-intensive
                      }
                      reverted<-F
                    } else { # operations are not in 'reverted' mode, so increment iAMA
                      if (!noncompliance | gpost) {
                        tmp[(yi+1):nyr,7]<-min(steps,tmp[yi,7]+1) # increment in the iAMA level (unless it is already maxed out)
                      }
                      if (tmp[yi+1,7] >= steps){ # then short-term trigger has fired maximum number of times to signal noncompliance
                        noncompliance<-T
                        res[yi+1,2,simi]<-0
                      }
                      if (!noncompliance){
                        tmp[(yi+1):nyr,11]<-iAMAschedule[tmp[yi+1,7],1] # the rho for the following years is adjusted according to the AMA schedule
                        if (iAMAschedule[tmp[yi,7]+1,2] > 0){
                        # i.e. if the AMA schedule requires additional intensive monitoring after triggering
                        #      then schedule that monitoring in the db
      #                    yrsam<-min(nyr,yi+toR(iAMAdata[[tmp[yi+1,7],2]]))
                          yrsam<-min(nyr,yi+iAMAschedule[tmp[yi+1,7],2])
                          tmp[(yi+1):yrsam,2]<-1 # ...intensive sampling for the next (few) years (as dictated by the short-term trigger)
                          #tmp[yi+1] gives the AMA level; 2 gives the number of years of intensive sampling required at that level
                          if (yrsam < nyr) tmp[(yrsam+1):nyr, 2]<-2
                        }
                      }
                      if (noncompliance & gpost){  # triggers have maxed out but monitoring and estimation continue
                        tmp[(yi+1):nyr, 2] <- 3
                        tmp[(yi+1):nyr, 7] <- steps
                        tmp[(yi+1):nyr, 11] <- rhoinf
                      }
                    }
                  }
                }
              }
              if (tmp[yi,8]==1){
                 if (viewone==0){
                    polygon(simi/nsim+c(-1,-1,0,0)/nsim,c((yi-1)/nyr,yi/nyr,yi/nyr,(yi-1)/nyr),col=colorRampPalette(colr)( nsim )[simi],border=NA) # showing progress of calculations
                 }
                 if (!noncompliance | gpost) next
              }
            }
            if (noncompliance){ #then pull the plug if no monitoring later, but continue estimating if monitoring continues
              if (yi < nyr){
                tmp[(yi+1):nyr,11]<-rhoinf
                tmp[(yi+1):nyr,2]<-3
                if (!gpost){ # no monitoring or estimation after final AMA
                  if (rhoinf == 0){
                    tmp[(yi+1):nyr,1]<-0
                    tmp[(yi+1):nyr,3]<-0
                  } else {
                    tmp[(yi+1):nyr,1]<-rpois(nyr-yi, lambda*rhoinf)
                    tmp[(yi+1):nyr,3]<-rbinom(nyr-yi,tmp[(yi+1):nyr,1],g[3,1])
                  }
                } else { # yes monitoring after final AMA
      #            tmp[(yi+1):nyr,3]<-
      #            lam<-lambda*tmp[yi,11] # rate was determined at the end of last year's monitoring (with default = 1 but changes with triggering)
                }
              }
              if (viewone==0) polygon(simi/nsim+c(-1,-1,0,0)/nsim,c((yi-1)/nyr,1,1,(yi-1)/nyr),col=colorRampPalette(colr)( nsim )[simi],border=NA) # showing progress of calculations
              if (!gpost) break # no need to test for reversion or continue monitoring and testing in future years
            }
            # check reversion trigger
            if (rT && !noncompliance){#test for reversion; if there is non-compliance, then no reversion is possible
              if (!ltT){ #then we need to calculate pBa and pBb
                mu<-g[tmp[1:yi,2],1] # where gi
                minp<-g[tmp[1:yi,2],2]; maxp<-g[tmp[1:yi,2],3]; a<-tmp[1:yi,11]/sum(tmp[1:yi,11])
                x<-sum(tmp[1:yi,3])
                sig2<-((maxp-minp)/4)^2 # estimated means and uncertainties for g
                Eg<-sum(mu*a)/sum(a) # average, expected G is the weighted average of the observed g's
    #            Lest<-max(x, 0.5)/Eg
                Vg<-sum(sig2*a^2)#+(sum(a*(mu^2+sig2))-Eg^2)/Lest
                pBa<-Eg^2/Vg*(1-Eg)-Eg; pBb<-pBa*(1/Eg-1) # these are the two shape parameters for the beta distribution underlying the beta-binomial for X | M
                tmp[yi,12]<-pBa; tmp[yi,13]<-pBb
              }
              if (pBa*pBb/((pBa+pBa)^2*(pBa+pBa+1))<0.000001){
                tmp[yi,6]<-pLgX(rhorev*Tau/nyr*yi,x,Eg,fmmax(x,Eg))
              } else {
                tmp[yi,6]<-pLgX.ab(rhorev*Tau/nyr*yi,x,pBa,pBb,fmmax.ab(x,pBa,pBb))
              }
              if (tmp[yi,6]>1-rta){# then reversion
                reverted<-T
                tmp[yi:nyr, 7]<-0 # iAMA is reset
                tmp[yi,10]<-1
                if (yi<nyr){
                  tmp[(yi+1):nyr,11]<-1/rhorev
                }
              }
            }
          #  polygon(yi+c(-1,-1,0,0),c(0,1,1,0),col=7) # showing progress of calculations
               if (viewone==0) polygon(simi/nsim+c(-1,-1,0,0)/nsim,c((yi-1)/nyr,yi/nyr,yi/nyr,(yi-1)/nyr),col=colorRampPalette(colr)( nsim )[simi],border=NA) # showing progress of calculations
      #         print(paste(c(yi,tmp[yi,],'\n'),collapse=' ')); flush.console() # not printed after stT fires?
            }
          res[,,simi]<-tmp
        }
      }
      # M at conclusion of project
      if (viewone==0) {
        M<-colSums(res[,1,]) #apply(res[,1,],F=sum,M=2)
      } else {
        M<-sum(res[,1,])
      }

      # years before avoidance
      if (stT & iAMA){ #then avoidance when AMAi = steps
        junk<-(res[,7,]==steps)*(1:nyr)-1
      } else if (rT | (ltT & gpost==1 & g3i == 1)) { # avoidance when ltT = 1 [NOTE: when !rT, this quantity ('junk') has already been calculated as 'avoid']
        junk<-(res[,9,]==1)*(1:nyr)
      }
      if (viewone==0){
        if (rT | (stT & iAMA)| (ltT & gpost==1 & g3i == 1)){
          junk[junk==0 | junk==-1]<-nyr+1
  #        avoid<-apply(junk,F=min,M=2)
          avoid<-matrixStats::colMins(junk)
        }
        if (rT | (stT & iAMA) | (ltT & gpost==1 & g3i == 1)){
          Mstar<-array(numeric(nsim*2),dim=c(nsim,2)) # Mstar at triggering and at end of project
          for (simi in 1:nsim) {
            Mstar[simi,1]<-res[min(avoid[simi],nyr),4,simi]
            Mstar[simi,2]<-max(na.omit(res[,4,simi]))
          }
        }
        # distribution of short-term trigger firings
        if (stT){
          if (rT | (stT & iAMA) | (ltT & gpost==1 & g3i == 1)){
            stTf<-apply(res[,8,],F=sum,M=2)
          } else {
            ### calculate the number of times the short-term trigger fires for each simulation run
            ## step 1: calculate the running totals for x3 and g3
            x3<-array(dim=c(nyr, nsim))
            g3<-array(dim=c(nyr, 2)) # pBa and pBb parameters for running averages
            Eg3<-g[1,1] # average, expected G is the weighted average of the observed g's
            Vg3<-((g[1,3]-g[1,2])/4)^2
            Vg3<-max(0.000001, Vg3)
            g3[1,1]<-Eg3^2/Vg3*(1-Eg3)-Eg3; g3[1,2]<-g3[1,1]*(1/Eg3-1)
            tmpg<-c(rep(1, yrs), rep(2, nyr-yrs))
            X<-res[,3,]
            x3[1,]<-X[1,]
            for (yi in 2:nyr){
              minyi<-max(1,yi-sty+1)
              if(minyi<yi) x3[yi,]<-colSums(X[minyi:yi,]) else x3[yi,]<-X[minyi:yi,]
              mu<-g[tmpg[minyi:yi],1]
              minp<-g[tmpg[minyi:yi],2]; maxp<-g[tmpg[minyi:yi],3];
              sig2<-((maxp-minp)/4)^2 # estimated means and uncertainties for g
              a<-rep(1/(yi-minyi+1),yi-minyi+1)
              Eg3<-sum(mu*a)/sum(a) # average, expected G is the weighted average of the observed g's
              Vg3<-sum(sig2*a^2)
              # NOTE: this measure assumes that relative weights remain the same for all years in the window
              # If ltT reduces rate by factor of rho, this assumption is wrong and needs to be corrected for the few years following ltT
              Vg3<-max(0.000001, Vg3)
              g3[yi,1]<-Eg3^2/Vg3*(1-Eg3)-Eg3; g3[yi,2]<-g3[yi,1]*(1/Eg3-1)
              # these are the two shape parameters for the beta distribution underlying the beta-binomial for X | M
              # rarely, there may be an error (when the weights are wrong for yrs-1 years after a ltT
              # the error only comes into play when the stT fires in the yrs-1 years after ltT is fired
              # error is corrected on a case-by-case basis below...
            }
            ## step 2: detemine which x3's fire stT for which years
            ## for initial years (g1)...
            Tix<-0 # minimum x that triggers short-term for the initial monitoring intensity
            if (yrs>0){
              Eg<-g[1,1] # average, expected G is the weighted average of the observed g's
              Vg<-max(0.000001,((g[1,3]-g[1,2])/4)^2)
              # NOTE: this measure assumes that relative weights remain the same for all years in the window
              # If ltT reduces rate by factor of rho, this assumption is wrong and needs to be corrected for the few years following ltT
              betaa<-Eg^2/Vg*(1-Eg)-Eg; betab<-betaa*(1/Eg-1)
              while (1) {
                if (pLgX.ab(ifelse(stTon=='same',Tau/nyr, tau)*sty, Tix, betaa, betab, fmmax.ab(Tix,betaa,betab)) < sta) break
                Tix<-Tix+1
              }
            }
            Tfx<-0 # minimum x that triggers short-term for the less-intensive monitoring intensity
            if (yrs>0){
              Eg<-g[2,1] # average, expected G is the weighted average of the observed g's
              Vg<-max(0.000001,((g[2,3]-g[2,2])/4)^2)
              # NOTE: this measure assumes that relative weights remain the same for all years in the window
              # If ltT reduces rate by factor of rho, this assumption is wrong and needs to be corrected for the few years following ltT
              betaa<-Eg^2/Vg*(1-Eg)-Eg; betab<-betaa*(1/Eg-1)
              while (1) {
                if (pLgX.ab(ifelse(stTon=='same',Tau/nyr, tau)*sty, Tfx, betaa, betab, fmmax.ab(Tfx,betaa,betab)) < sta) break
                Tfx<-Tfx+1
              }
            }
            Tiix<-NULL
            if (yrs > 0 & sty > 1 ){ #then there is transition years between intensive and non-intensive g3's
              Tiix<-numeric(sty-1) # array of x's for transition between intensive and non-intensive
              for (i in 1:length(Tiix)){
                ind<-c(rep(1, sty-i), rep(2, i))
                mu<-g[ind,1]
                minp<-g[ind,2]; maxp<-g[ind,3];
                sig2<-((maxp-minp)/4)^2 # estimated means and uncertainties for g
                a<-rep(1/sty,sty)
                Eg3<-sum(mu*a)/sum(a) # average, expected G is the weighted average of the observed g's
                Vg3<-sum(sig2*a^2)
                # NOTE: this measure assumes that relative weights remain the same for all years in the window
                # If ltT reduces rate by factor of rho, this assumption is wrong and needs to be corrected for the few years following ltT
                Vg3<-max(0.000001, Vg3)
                betaa<-Eg3^2/Vg3*(1-Eg3)-Eg3; betab<-betaa*(1/Eg3-1)
                while(1){
                  if (pLgX.ab(ifelse(stTon=='same',Tau/nyr, tau)*sty, Tiix[i], betaa, betab, fmmax.ab(Tiix[i],betaa,betab)) < sta) break
                  Tiix[i]<-Tiix[i]+1
                }
              }
            }
            if (nyr-length(Tiix)-yrs>0) {
              Tx<-c(rep(Tix, yrs), Tiix, rep(Tfx, nyr-length(Tiix)-yrs))
            } else {
              tadd<-dim(x3)[1]-length(Tix)*yrs
              if (tadd>0) Tx<-c(rep(Tix, yrs), Tiix[1:tadd]) else Tx<-rep(Tix, yrs)
            }
            stTf<-colSums(x3>=Tx)
          }
        }
        # distribution of reversion trigger firings
        if (rT){
          rTftmp<-rbind(res[,10,],31)
          rTf<-numeric(nsim)
          for (simi in 1:nsim) rTf[simi]<-min(which(rTftmp[,simi]>0))
        }

        .Rvar$res <- res
        wd<-7.5
        graphics.off()
        dev.new(height=wd,width=wd)
        layout(matrix(c(1,2,3),1,3),wid=c(1,1,1/3))
        par(mar=c(3,3,1,.5),mgp=c(2,.7,0),family='sans',oma=c(18,0,0,1))#,mfrow=c(1,3))
        if (sum(is.na(Mstar[,2]))>0) ploM <- quantile(Mstar[,1],c(0.05,0.9), type=3) else ploM <- quantile(Mstar[,2],c(0.05,0.9), type = 2)
        plot(0,0,type='n',
          xlim=c(.5,2.5), xaxs='i',
          ylim=range(c(Tau,quantile(M,c(.05,.95), type = 2),ploM)),
          xlab='',axes=F,
          ylab='Number of Fatalities')
        axis(2)
        axis(1,at=1:2,lab=c("Actual\n(M)",paste0("Estimated\n(M* for \u03b1 = ", lta, ")")),padj=.5)
        polygon(par('usr')[c(1,2,2,1)],c(par('usr')[c(3,3)],c(Tau,Tau)),col=colors()[600],border=NA)
        polygon(par('usr')[c(1,2,2,1)],c(par('usr')[c(4,4)],c(Tau,Tau)),col=colors()[394],border=NA)
        qtls<-quantile(M, probs=c(.05,.1,.25,.5,.75,.9,.95), type = 2)
        polygon(1+c(-.25,.25,.25,-.25),qtls[c(3,3,5,5)])
        lines(1+c(-.25,.25),rep(qtls[4],2))
        lines(rep(1,2),qtls[1:2],lty=3)
        lines(rep(1,2),qtls[2:3],lty=1)
        lines(rep(1,2),qtls[5:6],lty=1)
        lines(rep(1,2),qtls[6:7],lty=3)
        points(1,mean(M),pch=16)
        if(.Rvar$platform == "windows") bringToTop()
#        if (F & gpost){# && noncompliance){
        if (gpost){# && noncompliance){
        # gpost = T if monitoring continues after trigger
        # noncompliance = T if
        # distribution of Mstar at time of final AMA
          qtls<-quantile(Mstar[,1], probs=c(.05,.1,.25,.5,.75,.9,.95), type = 2)
          polygon(2.4+0.2*c(-.25,.25,.25,-.25),qtls[c(3,3,5,5)], border = colors()[623])
          lines(2.4+0.2*c(-.25,.25),rep(qtls[4],2), col = colors()[623])
          lines(0.4+rep(2,2),qtls[1:2],lty=3, col = colors()[623])
          lines(0.4+rep(2,2),qtls[2:3],lty=1, col = colors()[623])
          lines(0.4+rep(2,2),qtls[5:6],lty=1, col = colors()[623])
          lines(0.4+rep(2,2),qtls[6:7],lty=3, col = colors()[623])
          par(family='serif',xpd=T)
          points(0.4+2,mean(Mstar[,1]),pch=16, col = colors()[623])
          text(2.3, mean(Mstar[,1]), lab = "at trigger", srt=90, col = colors()[623], cex=.Rvar$charSizeAdjust)
          qtls<-quantile(Mstar[,2],probs=c(.05,.1,.25,.5,.75,.9,.95), type = 2)
          polygon(2+c(-.25,.25,.25,-.25),qtls[c(3,3,5,5)])
          lines(2+c(-.25,.25),rep(qtls[4],2))
          lines(rep(2,2),qtls[1:2],lty=3)
          lines(rep(2,2),qtls[2:3],lty=1)
          lines(rep(2,2),qtls[5:6],lty=1)
          lines(rep(2,2),qtls[6:7],lty=3)
          par(family='serif',xpd=T)
          text(.6,Tau,'\u03a4',adj=c(.5,.4), cex=.Rvar$charSizeAdjust)
          text(2,par('usr')[3]-0.15*(diff(par('usr')[3:4])),paste('\u03b1 = ',lta,sep=''), cex=.Rvar$charSizeAdjust)
          par(family='sans',xpd=F)
          points(2,mean(Mstar[,2]),pch=16)
        } else {
          qtls<-quantile(Mstar[,1],probs=c(.05,.1,.25,.5,.75,.9,.95), type = 2)
          polygon(2+c(-.25,.25,.25,-.25),qtls[c(3,3,5,5)])
          lines(2+c(-.25,.25),rep(qtls[4],2))
          lines(rep(2,2),qtls[1:2],lty=3)
          lines(rep(2,2),qtls[2:3],lty=1)
          lines(rep(2,2),qtls[5:6],lty=1)
          lines(rep(2,2),qtls[6:7],lty=3)
          par(family='serif',xpd=T)
          text(.6,Tau,'\u03a4',adj=c(.5,.4), cex=.Rvar$charSizeAdjust)
          text(2,par('usr')[3]-0.15*(diff(par('usr')[3:4])),paste('\u03b1 = ',lta,sep=''), cex=.Rvar$charSizeAdjust)
          par(family='sans',xpd=F)
          points(2,mean(Mstar[,1]),pch=16)
        }
        par(family='serif',xpd=T)
        text(1,par('usr')[3]-0.15*(diff(par('usr')[3:4])),paste('\u03a4 = ',Tau,sep=''), cex=.Rvar$charSizeAdjust)
        par(family='sans',xpd=F)
        box()

        # years of operation before avoidance
        plot(0,0,type='n',xlim=c(.5,2.5),xaxs='i',ylim=c(0,nyr+1),yaxs='i',xlab='',axes=F,ylab='Years of operation before...')
        clrs<-c(colors()[254],colors()[254],colors()[11],colors()[114],colors()[250])
        if (ltT | iAMA){
          qtls<-c(0,quantile(pmin(avoid,nyr+1),seq(.01,1,by=.01), type = 2))
          for(i in 1:(length(qtls)-1)){
            polygon(1+.5*c(-1,1,1,-1),c(rep(qtls[i],2),rep(qtls[i+1],2)),col=colorRampPalette(clrs)( 100 )[i],border=NA)
          }
          qtls<-quantile(pmin(avoid,nyr+1),probs=c(.05,.1,.25,.5,.75,.9,.95), type = 2)
          polygon(1+c(-.25,.25,.25,-.25),qtls[c(3,3,5,5)])
          lines(1+c(-.25,.25),rep(qtls[4],2))
          lines(rep(1,2),qtls[1:2],lty=3)
          lines(rep(1,2),qtls[2:3],lty=1)
          lines(rep(1,2),qtls[5:6],lty=1)
          lines(rep(1,2),qtls[6:7],lty=3)
#          points(1,mean(pmin(avoid,nyr)),pch=16)
          par(family='serif',xpd=T)
          text(1,par('usr')[3]-0.15*(diff(par('usr')[3:4])),bquote(rho[infinity]==.(rhoinf)) , cex=.Rvar$charSizeAdjust)
          par(family='sans',xpd=F)
        } else {
          polygon(1+.5*c(-1,-1,1,1),par('usr')[c(3,4,4,3)],col=colors()[250])
          text(1,mean(par('usr')[c(3,4)]),'NA',cex=2*.Rvar$charSizeAdjust)
        }
        # years of operation before reversion
        #plot(0,0,type='n',xlim=c(.5,2.5),xaxs='i',ylim=c(0,nyr),yaxs='i',xlab='',axes=F,ylab='Years of operation before...')
        if (rT){
          rTf[rTf==0]<-nyr+1
          qtls<-c(0,quantile(rTf,seq(.01,1,by=.01), type = 2))
          clrs<-clrs[length(clrs):1]
          for(i in 1:(length(qtls)-1)){
            polygon(2+.5*c(-1,1,1,-1),c(rep(qtls[i],2),rep(qtls[i+1],2)),col=colorRampPalette(clrs)( 100 )[i],border=NA)
          }
          qtls<-quantile(rTf,probs=c(.05,.1,.25,.5,.75,.9,.95), type = 2)
          polygon(2+c(-.25,.25,.25,-.25),qtls[c(3,3,5,5)])
          lines(2+c(-.25,.25),rep(qtls[4],2))
          lines(rep(2,2),qtls[1:2],lty=3)
          lines(rep(2,2),qtls[2:3],lty=1)
          lines(rep(2,2),qtls[5:6],lty=1)
          lines(rep(2,2),qtls[6:7],lty=3)
#          points(2,mean(rTf),pch=16)
          par(family='serif',xpd=T)
          text(2,par('usr')[3]-0.15*(diff(par('usr')[3:4])),paste('\u03b1 = ',rta,'\n\u03c10 = ',rhorev,sep='') , cex=.Rvar$charSizeAdjust)
          par(family='sans',xpd=F)
        } else {
          polygon(2+.5*c(-1,-1,1,1),par('usr')[c(3,4,4,3)],col=colors()[250])
          text(2,mean(par('usr')[c(3,4)]),'NA',cex=2*.Rvar$charSizeAdjust)
        }
        axis(2)
        axis(1,at=1:2,lab=c("Avoidance","Reversion"),padj=.5)
        lines(rep(1.5,2),par('usr')[3:4])
        box()
        if(.Rvar$platform == "windows") bringToTop()
        # number of times short-term trigger fires
        plot(0,0,type='n',xlim=c(0,1),xaxs='i',ylim=c(0,1),yaxs='i',ylab="Fraction of projects",axes=F,xlab='')
        clrs<-c(3,7,2,colors()[440],8)
        polygon(c(0,1,1,0),c(0,0,1,1),col=clrs[5])
        axis(2)
        axis(2,at=seq(.1,.9,by=.2),lab=F)
        axis(1,at=.5,lab="Frequency\nof short-term\n trigger firing",padj=.7)
        if(stT){
          ff<-c(0,cumsum(tabulate(stTf+1))/nsim)
          for (i in 1:length(ff)){
            polygon(c(0,1,1,0),ff[c(i,i,i+1,i+1)],col=clrs[min(i,5)])
          }
        #  for (i in 1:5){
        #    polygon(c(.9,1,1,.9),c(i,i,i-1,i-1)*.2,col=clrs[i])
        #  }
        } else {
          text(.5,.5,"NA",cex=2*.Rvar$charSizeAdjust)
        }
        par(family='serif',xpd=T)
        ln<-1
        mtext(side=1,line=ln,'--------------------------------------------------------------------------------------------------------------------------',outer=T,adj=0, cex=.Rvar$charSizeAdjust); ln<-ln+1
        mtext(side=1,line=ln,paste('  Take limit (\u03a4) = ',Tau,' in ',nyr, ' yrs,  Baseline fatality rate (\u03bb) = ',lambda,', number of reps = ',nsim,sep=''),outer=T,adj=0, cex=.Rvar$charSizeAdjust); ln<-ln+1.2
        mtext(side=1,line=ln,paste('  Monitoring: ',yrs,' years of intensive monitoring at start of permit',sep=''),outer=T,adj=0, cex=.Rvar$charSizeAdjust);ln<-ln+1.2
        mtext(side=1,line=ln,paste('  Detection probability, intensive monitoring: g = ',g[1,1],' with 95% CI = [',g[1,2],', ',g[1,3], ']',sep=''),outer=T,adj=0, cex=.Rvar$charSizeAdjust);ln<-ln+1.2
        mtext(side=1,line=ln,paste('  Detection probability, non-intensive monitoring: g = ',g[2,1],' with 95% CI = [',g[2,2],', ',g[2,3], ']',sep=''),outer=T,adj=0, cex=.Rvar$charSizeAdjust);ln<-ln+1.2
        lab1<-'Long-term trigger: '
        if (ltT){
          lab2<-paste('\u03b1 = ',lta,sep='')
        } else {
          lab2<-'NA'
        }
        mtext(side=1,line=ln,paste('  ',lab1,lab2,sep=''),adj=0,outer=T, cex=.Rvar$charSizeAdjust); ln<-ln+1.2
        lab1<-'Reversion trigger: '
        if (rT){
          lab2<-paste('\u03c10 = ',rhorev,', \u03b1 = ',rta,sep='')
        } else {
          lab2<-'NA'
        }
        mtext(side=1,line=ln,paste('  ',lab1,lab2,sep=''),adj=0,outer=T, cex=.Rvar$charSizeAdjust); ln<-ln+1.2

        lab1<-'Short-term trigger: '
        if (stT){
          lab2<-paste('\u03b1 = ',sta,' for rolling average rate over ',sty,' years, ', "\u03c4 = ", ifelse(stTon=='same',round(Tau/nyr,2), tau),sep='')
        } else {
          lab2<-"NA"
        }
        mtext(side=1,line=ln,paste('  ',lab1,lab2,sep=''),adj=0,outer=T, cex=.Rvar$charSizeAdjust); ln<-ln+1.2
        if(stT){
          if (iAMA){
            ln<-ln+1.4
  #          lab2<-paste0('Incremental AMA: \u03c1 = ',toR(iAMAdata[[1,1]]))
            lab2<-paste0('Incremental AMA: \u03c1 = ',iAMAschedule[1,1])
            for (i in 2:steps){
  #            lab2<-paste(lab2,toR(iAMAdata[[i,1]]),sep=', ')
              lab2<-paste(lab2,iAMAschedule[i,1],sep=', ')
            }
            lab2<-paste(lab2,'\n   and',sep='')
            for (i in 1:steps){
  #            lab2<-paste(lab2,' ',toR(iAMAdata[[i,2]]),sep='')
              lab2<-paste(lab2,' ',iAMAschedule[i,2],sep='')
            }
            lab2<-paste(lab2,' years of intensive monitoring following each trigger firing',sep='')
          } else {
            lab2<-' with no intensive monitoring after triggering'
          }
        }
        mtext(side=1,line=ln,paste('  ',lab2),outer=T,adj=0, cex=.Rvar$charSizeAdjust);ln<-ln+1.2
        mtext(side=1,line=ln,'--------------------------------------------------------------------------------------------------------------------------',outer=T,adj=0, cex=.Rvar$charSizeAdjust); ln<-ln+1.2
        mtext(side=1,line=ln,'  Boxplots show IQR with median, and whiskers to show 5th, 10th, 90th, and 95th percentiles; \u25cf = mean',outer=T,adj=0, cex=.Rvar$charSizeAdjust); ln<-ln+1.2
        if (stT){
          mtext(side=1,line=ln,'  Number of short-term trigger firings in project: ',adj=0,outer=T, cex=.Rvar$charSizeAdjust)
          xst<--7; xincr<-1;clrs<-c(3,7,2,colors()[440],8)
          layout(1)
          par(new=T,mar=c(0,0,0,0),oma=c(0,0,0,0),fig=c(.42,.97,.005,.035))
          plot(1,1,axes=F,xlab='',ylab='',pch=16,col=5,xlim=c(0,1),xaxs='i',ylim=c(0,1),xaxs='i',type='n')
          for (xi in 1:5){
            polygon(xincr*c(xi,xi,xi-1,xi-1)/5,c(0,1,1,0),col=clrs[xi])
            text((xi-.5)/5,.5,ifelse(xi<5,xi-1,'4+'), cex=.Rvar$charSizeAdjust)
          }
        }
        par(family='sans',xpd=F)
      } else {# viewone==1
        endo<-min(min(which(res[,11,1]==rhoinf)-1,nyr))
        empg<-sum(res[1:endo,3,1])/sum(res[1:endo,1,1])
      # create new form to show results (two tables: inputs and outputs)
      # output table (not editable)
        for (i in 1:nyr){
          scexoutData[[i,0]]<-i
          scexoutData[[i,1]]<-res[i,1,1]
          scexoutData[[i,2]]<-res[i,3,1]
          scexoutData[[i,3]]<-g[res[i,2,1],1]
          scexoutData[[i,4]]<-round(res[i,11,1], 2)
          scexoutData[[i,5]]<-round(res[i,11,1]*lambda,2)
          scexoutData[[i,6]]<-sum(res[1:i,1,1])
          scexoutData[[i,7]]<-sum(res[1:i,3,1])
          if (!is.na(gcum[i])) scexoutData[[i,8]]<-round(gcum[i],3) else  scexoutData[[i,8]]<-"NA"
          if (ltT && !is.na(res[i,4,1])) scexoutData[[i,9]]<-res[i,4,1] else scexoutData[[i,9]]<-"NA"
          if (rT && !is.na(res[i,6,1]))  scexoutData[[i,10]]<-round(res[i,6,1],3) else  scexoutData[[i,10]]<-"NA"
          if (stT && !is.na(Mmov[i])){
            scexoutData[[i,11]]<-Mmov[i]
            scexoutData[[i,12]]<-Xmov[i]
            scexoutData[[i,13]]<-round(gmov[i],3)
            scexoutData[[i,14]]<-round(1-res[i,5,1],3)
          } else {
            scexoutData[[i,11]]<-"NA"
            scexoutData[[i,12]]<-"NA"
            scexoutData[[i,13]]<-"NA"
            scexoutData[[i,14]]<-"NA"
          }
        }
        for (i in 1:nyr){
          val<-suppressWarnings(as.numeric(toR(scexoutData[[i,6]])))
          if (!is.na(val) && val > Tau) {
            tcl(scexoutTable, "tag", "celltag", "exceeds", cellind(i,6))
          } else {
            tcl(scexoutTable, "tag", "celltag", "delete", cellind(i,6)) # removes all tags for the given cell
          }
          val<-suppressWarnings(as.numeric(toR(scexoutData[[i,9]])))
          if (!is.na(val) && val > Tau) {
            tcl(scexoutTable, "tag", "celltag", "triggered", cellind(i,9))
          } else {
            tcl(scexoutTable, "tag", "celltag", "delete", cellind(i,9)) # removes all tags for the given cell
          }
          val<-suppressWarnings(as.numeric(toR(scexoutData[[i,14]])))
          if (!is.na(val) && val > 1-sta) {
            tcl(scexoutTable, "tag", "celltag", "triggered", cellind(i,14))
          } else {
            tcl(scexoutTable, "tag", "celltag", "delete", cellind(i,14)) # removes all tags for the given cell
          }
          val<-suppressWarnings(as.numeric(toR(scexoutData[[i,10]])))
          if (!is.na(val) && val >= 1-rta){
            tcl(scexoutTable, "tag", "celltag", "reverted", cellind(i,10))
          } else {
            tcl(scexoutTable, "tag", "celltag", "delete", cellind(i,10)) # removes all tags for the given cell
          }
        }
      }
    })
}




