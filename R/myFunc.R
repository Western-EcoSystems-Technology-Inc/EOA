myCalc<-function(mydat, sides){ # different graphs, depending on which options are selected
  binog<-F
# first stage is to write the values to
  aCI<-mydat$aCI[1]
  if (mydat$option=="M" & (mydat$Mtype %in% c("T","P"))){
    # all M options share initial data setup
    yrs<-length(mydat$years) # number of years of monitoring
    nyr<-mydat$nyr
    gcum<-array(dim=c(yrs,3)) # cumulative detection probability (with beta parameters)
    Mest.yi<-numeric(yrs) # estimated M for the given aM for each year
    #### initial years
    X<-mydat$X;  x<-sum(X)
    pBa0<-mydat$Ba; pBb0<-mydat$Bb
    a<-mydat$rel_wt/sum(mydat$rel_wt)
    aM<-1-mydat$crlev
    Tau<-mydat$Tau
    mu<-pBa0/(pBa0+pBb0)
    sig2<-pBa0*pBb0/((pBa0+pBb0)^2*(pBa0+pBb0+1))
    nq<-200 # number of quantiles to display from the posterior of M (going from 50th to 95th)
    Mstar<-array(dim=c(yrs,nq)) # quantiles of posterior from 50th to 95th
    aseq<-seq(.05,.95,length=nq)
    Eg<-cumsum(mu*a)/cumsum(a) # average, expected G is the weighted average of the observed g's
    .Rvar$Eg<-Eg
    if (yrs > 1)
      gcum[1:yrs,1]<-Eg #else gcum[1]<-Eg
#    lambda<-pmax(cumsum(X), 0.5)/Eg
    Vg<-numeric(yrs)
    Vg[1]<-sig2[1]
    if (yrs>1){
      for (i in 2:yrs){
        aa<-a[1:i]/sum(a[1:i])
        Vg[i]<-sum(sig2[1:i]*aa[1:i]^2)#+(sum(aa[1:i]*(mu[1:i]^2+sig2[1:i]))-Eg[i]^2)/lambda[i]
      }
    }
    pBa<-Eg^2/Vg*(1-Eg)-Eg; pBb<-pBa*(1/Eg-1) # these are the two shape parameters for the beta distribution underlying the beta-binomial for X | M
    .Rvar$pBa<-pBa; .Rvar$pBb<-pBb
#    if (yrs > 1){
      gcum[1:yrs,2]<-pBa
      gcum[1:yrs,3]<-pBb
#    } else {
#      gcum[2]<-pBa
#      gcum[3]<-pBb
#    }
    # M needs to be calculated up to a reasonable maximum, but the final distribution is clipped at P(M > m) < 0.0001
    mmax<-NA # declaration of the variable
    .Rvar$myTrack<-array(dim=c(yrs, 11))
    # columns:
    # 1    => M*
    # 2-10 => quantiles of posterior: 0.025, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.975
    resalpha<-c(0.025, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.975)
    # 11   => mean(M)
    for (yi in 1:yrs){
      mmax<-ifelse (Vg[yi] < 0.00001, fmmax(sum(X[1:yi]),.Rvar$Eg[yi]),fmmax.ab(sum(X[1:yi]),.Rvar$pBa[yi],.Rvar$pBb[yi]))
      x<-sum(X[1:yi])
      M<-x:mmax
    ############
        if (Vg[yi]>0.00001){
    #      pBa<-Eg^2/Vg*(1-Eg)-Eg; pBb<-pBa*(1/Eg-1) # these are the two shape parameters for the beta distribution underlying the beta-binomial for X | M
          pXgM<-VGAM::dbetabinom.ab(x,size=M,shape1=.Rvar$pBa[yi],shape2=.Rvar$pBb[yi]) # the probabilities of X for M = 0:mmax
        } else {
          pXgM<-dbinom(x,size=M,prob=.Rvar$Eg[yi]) # the probabilities of X for M = x:mmax
        }
        pM<-diff((sqrt(c(x-1,x:mmax)+1)))
        pMgX<-pXgM*pM; pMgX<-pMgX/sum(pMgX) # posterior distribution for M (ignoring M < x, which has probability = zero)
    #####################
      for (ai in 1:nq) Mstar[yi,ai]<-M[min(which(1-cumsum(pMgX)<1-aseq[ai]))]
      Mest.yi[yi]<-M[min(which(1-cumsum(pMgX)<aM))]
      .Rvar$myTrack[yi,1]<-Mest.yi[yi]
      for (ai in 1:length(resalpha)) .Rvar$myTrack[yi, 1 + ai]<-M[min(which(1-cumsum(pMgX)<1-resalpha[ai]))]
      .Rvar$myTrack[yi,length(resalpha)+2]<-sum(M*pMgX)
    }
    x<-sum(X)
    M0<-M
    pMgX0<-pMgX
    if (mydat$Mtype=="T"){ # tracking
      # Tau    -- total permitted take over the course of the permit
      # nyr    -- length of the permit (or the number of years to show on the graph)
      # aM   -- alpha level for testing long-term trigger
      ############ tracking graph
      ### estimated
      graphics.off() # closes old graphics window and draws new on of the proper size
      ht<-7*.Rvar$charSizeAdjust
      dev.new(height=ht,width=2*ht, noRStudioGD =T)

      par(mar=c(3,3,6,.5),mgp=c(2,.7,0),family='sans')
      gtop<-max(c(Tau,Mstar))
      plot(0,0,type='n',xlim=c(0,nyr)+.5,ylim=c(0,gtop),xlab='Year',ylab='Fatalities',xaxs='i',axes=F)
      axis(1)
      axis(1,at=1:nyr,lab=F,tck=-0.01)
      if (par('usr')[4]<110){
        axis(2,at=seq(0,par('usr')[4],by=10))
        axis(2,at=seq(5,par('usr')[4],by=10),lab=F)
        axis(2,at=seq(1,par('usr')[4],by=1),lab=F, tck=-0.01)
      } else if (par('usr')[4]<150){
        axis(2,at=seq(0,par('usr')[4],by=20))
        axis(2,at=seq(10,par('usr')[4],by=10),lab=F)
        axis(2,at=seq(5,par('usr')[4],by=5),lab=F, tck=-0.01)
      } else {
        axis(2)
        if (par('usr')[4]<1100 & (round(par('yaxp')[2]/par('yaxp')[3])%%10==0))   axis(2,at=seq(0,par('yaxp')[2],by=10),lab=F, tck=-0.01)
      }
      title('Cumulative mortality (M)', cex.main = .Rvar$charSizeAdjust*14/12, adj=0)

      # estimation...
      par(family='serif',xpd=T)
      # the data
      colr<-c(colors()[c(400, 128, 31, 552, 552, 498, 652)],"#FFFFE7")
      medind<-round(mean(which(abs(aseq-.5)==min(abs(aseq-.5)))))
.Rvar$Mstar<-Mstar
.Rvar$aseq<-aseq
      for (yi in 1:yrs){
        xi<-yi+c(-.5,-.5,.5,.5)
        for (ai in 2:nq){
          polygon(xi,c(Mstar[yi,ai-1],Mstar[c(yi,yi),ai],Mstar[yi,ai-1]),col=colorRampPalette(colr)( nq )[ai],border=colorRampPalette(colr)( nq )[ai])
        }
        #if (abs(aM-.5)>=0.01)
        points(yi, Mest.yi[yi], pch=8)#, cex=.Rvar$charSizeAdjust*.7)
        lines(yi+c(-.5,.5),rep(Mstar[yi,1],2),lty=3)
        lines(yi+c(-.5,.5),rep(Mstar[yi,nq],2),lty=3)
        lines(yi+c(-.5,.5),rep(Mstar[yi,medind],2),col=colors()[259],lwd=1.5)
      }
      lines(par('usr')[1:2],rep(Tau,2),col=4)
      # legend
      par(family='serif')
      mtext(side = 3, at=sum(par('usr')[1:2]*c(.8,.2))-.2*nyr/30, line = 5, 'Posterior distribution of M:',adj=0, cex=.Rvar$charSizeAdjust, family = 'serif')
      tmp <- par('usr')[4]+seq((1-par('plt')[4])*0.05, 0.7*(1-par('plt')[4]), length = nq)*diff(par('usr')[3:4])/diff(par('plt')[3:4])
      xi<-sum(par('usr')[1:2]*c(.8,.2))+.3*nyr/30*c(1,1,-1,-1)
      for (i in nq:1)  polygon(xi,tmp[c(1,i,i,1)],col=colorRampPalette(colr)( nq )[i],border=colorRampPalette(colr)( nq )[i])
      lines(xi[c(1,3)],rep(tmp[round(.4*nq)],2),col=colors()[259], lwd=1.5)
      if (abs(aM-.5)>=0.01){
        text(xi[2]+.5*nyr/30, tmp[round(.7*nq)],adj=0,lab=paste0("M* for 1 - \u03b1 = ",1-aM), cex=.Rvar$charSizeAdjust)
        text(xi[2]+.5*nyr/30, tmp[round(.4*nq)],adj=0,lab='Median', cex=.Rvar$charSizeAdjust)
        points(mean(xi[2:3]),tmp[round(.7*nq)],pch=8, cex= .Rvar$charSizeAdjust*.7)
      } else {
        text(xi[2]+.5*nyr/30, tmp[round(.4*nq)],adj=0,lab='Median = M*', cex=.Rvar$charSizeAdjust)
        points(mean(xi[2:3]),tmp[round(.4*nq)],pch=8, cex= .Rvar$charSizeAdjust*.7)
      }
#      text(xi[2]+.5*nyr/30, tmp[40],adj=0,lab=paste('M* for \u03b1 = ',aM,sep=''), cex=.Rvar$charSizeAdjust)
      lines(xi[c(1,3)], rep(max(tmp), 2),lty=3)
      text(xi[2]+.5*nyr/30, max(tmp),adj=0,lab='95th percentile',cex=.Rvar$charSizeAdjust)
      lines(xi[c(1,3)], rep(min(tmp), 2),lty=3)
      text(xi[2]+.5*nyr/30, tmp[max(round(0.03*nq),1)],adj=0,lab='5th percentile',cex=.Rvar$charSizeAdjust)
      par(family='serif')
      text(1,Tau+diff(par('usr')[3:4])*.005,'\u03a4',adj=c(0.5,0),cex=.Rvar$charSizeAdjust)
      mtext(text = "Overall detection probability (cumulative):   ",
        side = 3, line = 4, at = par('usr')[2], adj = 1, cex=.Rvar$charSizeAdjust
      )
      mtext(text = paste0(" g = ", round(.Rvar$pBa[yrs]/(.Rvar$pBa[yrs]+.Rvar$pBb[yrs]),4), " 95% CI = [",round(qbeta(.025,.Rvar$pBa[yrs],.Rvar$pBb[yrs]),4),", ",round(qbeta(.975,.Rvar$pBa[yrs],.Rvar$pBb[yrs]),4),"]"),
        side = 3, line = 3, at = par('usr')[2], adj = 1, cex=.Rvar$charSizeAdjust
      )
      par(family='sans')
      box()
      if (min(M)>0){
        minm<-min(M)
        M<-c(0:(minm-1),M)
        pMgX<-c(rep(0,minm),pMgX)
      }
      ### estimate lambda
      ctprob<-0.0001 # A Jeffreys prior is used for the distribution of lambda. What is a reasonable maximum (and minimum) for lambda?
      # Pick a mmax = minimum M value such that P(X <= x | M) < ctprob. Then, define the Lmax to be minimum such that the P(lambda > Lmax) < ctprob.
      # p is assumed to be distributed as a beta RV with mean = pmean and variance = ((maxp-minp)/4)^2.
      # This determines beta parameters as follows:
#      Lest<-max(x, 0.5)/Eg
      gmin<-gmax<-0
      Vg<-sum(sig2*a^2)#+(sum(a*(mu^2+sig2))-Eg^2)/Lest
      if (Vg > 0.000001){
        mmax<-fmmax.ab(x,.Rvar$pBa[yrs],.Rvar$pBb[yrs])
        if (ppois(mmax,100)<ctprob){
          lvec<-seq(0,100,by=0.01)
          Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))] # this is the biggest that lambda could credibly be and still have the biggest M conceivable
        } else {
          i<-1
          while (1){
            i<-i+1
            if (ppois(mmax,100*i)<ctprob){
              lvec<-seq((i-1)*100,i*100,by=0.01)
              Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))]
              break
            }
          }
        }
      } else {
        mmax<-fmmax(x,.Rvar$Eg[yrs])
        if (ppois(mmax,100)<ctprob){
          lvec<-seq(0,100,by=0.01)
          Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))] # this is the biggest that lambda could credibly be and still have the biggest M conceivable
        } else {
          i<-1
          while (1){
            i<-i+1
            if (ppois(mmax,100*i)<ctprob){
              lvec<-seq((i-1)*100,i*100,by=0.01)
              Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))]
              break
            }
          }
        }
      }
      if (Vg > 0.000001){
        binog <- F
        cpLgX<-function(L,x,pBa,pBb,mmax) 1-pLgX.ab(L,x,pBa,pBb,mmax)
        meanL<-integrate(Vectorize(cpLgX,"L"),lower=0.0001,upper=Lmax,x=x,.Rvar$pBa[yrs],.Rvar$pBb[yrs],mmax=mmax)$val/sum(mydat$rel_wt)#styr
        CI<-c(optimize(f=optqL.ab,interval=c(0.00001,Lmax),x=x,.Rvar$pBa[yrs],.Rvar$pBb[yrs],mmax=mmax,p=.025)$minimum, optimize(f=optqL.ab,interval=c(0.00001,Lmax),x=x,.Rvar$pBa[yrs],.Rvar$pBb[yrs],mmax=mmax,p=.975)$minimum)/sum(mydat$rel_wt)#styr
      } else {
        binog <- T
        cpLgX<-function(L,x,g,mmax) 1-pLgX(L,x,g,mmax)
        meanL<-integrate(Vectorize(cpLgX,"L"),lower=0.0001,upper=Lmax,x=x,g=.Rvar$Eg[yrs],mmax=mmax)$val/sum(mydat$rel_wt)#styr
        CI<-c(optimize(f=optqL,interval=c(0.00001,Lmax),x=x,g=.Rvar$Eg[yrs],mmax=mmax,p=.025)$minimum, optimize(f=optqL,interval=c(0.00001,Lmax),x=x,g=.Rvar$Eg[yrs],mmax=mmax,p=.975)$minimum)/sum(mydat$rel_wt)#styr
      }
      par(family='serif')
      mtext(text = "Estimated baseline annual fatality rate:   ",
        side = 3, line = 1.5, adj = 1,cex = .Rvar$charSizeAdjust
      )
      mtext(text = paste0("\u03bb\u0302 = ", signif(meanL,3), " 95% CI = [",signif(CI[1],3),", ",signif(CI[2],3),"]"),
        side = 3, line = 0.5, adj = 1,cex = .Rvar$charSizeAdjust
      )
      text(1,Tau+diff(par('usr')[3:4])*.005,'\u03a4',adj=c(0.5,0),cex=.Rvar$charSizeAdjust)
      par(family='sans')
      while(1){ if (sink.number()==0) break else sink() }
      sink(paste0(.Rvar$datadir, "/output"))
      cat(paste("Summary statistics for mortality estimates through ", length(mydat$years), " years\n",sep=''))
      cat("------------------------------------\n")
      cat("Results\n")
#$$$$$$$$$$$$$$$
      yrs<-length(mydat$years) # number of years of monitoring
      nyr<-mydat$nyr
      gcum<-array(dim=c(yrs,3)) # cumulative detection probability (with beta parameters)
      Mest.yi<-numeric(yrs) # estimated M for the given aM for each year
      #### initial years
      ### The function calculates an overall average g along with estimated variance.
      ### Posterior is betabinomial on g and is calculated based on a uniform prior.
      X<-mydat$X; nyr<-length(X)
      x<-sum(X)
      pBa0<-mydat$Ba; pBb0<-mydat$Bb
      a<-mydat$rel_wt/sum(mydat$rel_wt)
      rel_wt<-mydat$rel_wt
      # p is assumed to be distributed as a beta RV with mean = pmean and variance = ((maxp-minp)/4)^2.
      # This determines beta parameters as follows:
      aM<-1-mydat$crlev
      Tau<-mydat$Tau
      mu<-pBa0/(pBa0+pBb0)
      sig2<-pBa0*pBb0/((pBa0+pBb0)^2*(pBa0+pBb0+1))
      Eg<-sum(mu*a) # average, expected G is the weighted average of the observed g's
#      lambda<-max(sum(X), 0.5)/Eg
      Vg<-sum(sig2*a^2)#+(sum(a*(mu^2+sig2))-Eg^2)/lambda
      pBa<-Eg^2/Vg*(1-Eg)-Eg; pBb<-pBa*(1/Eg-1) # these are the two shape parameters for the beta distribution underlying the beta-binomial for X | M
      pMgX<-postM.ab(x, pBa, pBb)
      cims<-MCI(pMgX, .95)
      cat('\n')
      if (sides == 2){
        cat(paste0(100*(1-aM),"% CI for M = [",cims[1],", ",  cims[2],"]\n"))
      } else {
        tmpmst<-calcMstar(pMgX, aM)
        cat(paste0("M* = ",tmpmst, " for 1 - \u03b1 = ",1-aM,", i.e., P(M <= ", tmpmst ,") >= ", (1-aM)*100, "%", '\n'))
      }
      cat(paste0("Estimated overall detection probability: g = ", signif(Eg, 3),
        ", 95% CI = [",signif(qbeta(.025, pBa, pBb),3),", ",signif(qbeta(.975, pBa, pBb),3),"]\n",
        "   Ba = ", signif(pBa, 5), ", Bb = ", signif(pBb, 5), "\n"))
#      CI<-c(L[min(which(pL>(1-aCI)/2))]/sum(mydat$rel_wt),L[min(which(pL>aCI+(1-aCI)/2))]/sum(mydat$rel_wt))
#      CI<-c(L[min(which(pL>.025))]/sum(mydat$rel_wt),L[min(which(pL>.975))]/sum(mydat$rel_wt))
#      par(family='serif')
      cat(paste0("Estimated baseline fatality rate (for rho = 1): lambda = ",signif(meanL,4),", 95% CI = [",signif(CI[1],3),", ",signif(CI[2],3),"]\n"))
      cat("\n")
      cat(paste0("Cumulative Mortality Estimates\n"))
      cat("                                              mean\n")
      cat("Year          X    g     M*  median  95% CI  lambda       95% CI\n")
      for (yi in 1:length(mydat$years)){
        postmi<-postM.ab(sum(mydat$X[1:yi]), .Rvar$pBa[yi], .Rvar$pBb[yi])
        mcii<-MCI(postmi, .95)
        lstat<-postL.abCI(sum(mydat$X[1:yi]), .Rvar$pBa[yi], .Rvar$pBb[yi])
        cat(sprintf("%-12s%3.0f %5.3f %4.0f  %4.0f     [%.0f, %.0f] %6.4g   [%.4g, %.4g]\n",
          mydat$years[yi],
          sum(mydat$X[1:yi]),
          .Rvar$Eg[yi],
          .Rvar$myTrack[yi,1],
          .Rvar$myTrack[yi,6],
          mcii[1],
          mcii[2],
          signif(lstat$meanL,4),
          signif(lstat$CI[1],4),
          signif(lstat$CI[2],4)))
      }

      cat(paste0("\n\nAnnual Mortality Estimates\n"))
      cat("                                              mean\n")
      cat("Year          X    g     M*  median  95% CI  lambda       95% CI\n")
      ht<-7*.Rvar$charSizeAdjust
      lstat<-list(meanL=numeric(length(mydat$years)), CI=array(dim=c(length(mydat$years), 2)))
      for (yi in 1:length(mydat$years)){
        postmi<-postM.ab(mydat$X[yi], mydat$Ba[yi], mydat$Bb[yi])
        mcii<-MCI(postmi, .95)
        junk<-postL.abCI(mydat$X[yi], mydat$Ba[yi], mydat$Bb[yi])
        lstat$meanL[yi]<-junk$meanL
        lstat$CI[yi,]<-junk$CI
        cat(sprintf("%-12s%3.0f %5.3f %4.0f  %4.0f     [%.0f, %.0f] %6.4f   [%.4f, %.4f]\n",
          mydat$years[yi],
          mydat$X[yi],
          round(mydat$Ba[yi]/(mydat$Ba[yi]+mydat$Bb[yi]),3),
          calcMstar(postmi, 1-mydat$crlev),
          calcMstar(postmi, 0.5),
          mcii[1],
          mcii[2],
          signif(lstat$meanL[yi],4),
          signif(lstat$CI[yi,1],4),
          signif(lstat$CI[yi,2],4)))
      }
      cat(paste0("\n\nTest of assumed relative weights (rho) and potential bias\n"))
      rhosumry<-rhotest(mydat)
      cat("              Fitted rho\n")
      cat("Assumed rho     95% CI \n")
      for (yi in 1:length(rhosumry$rho.ass)) {
        cat(sprintf(
          " %6.3g     [%5.3f, %5.3f]\n", rhosumry$rho.ass[yi], rhosumry$rhoqtls[1, yi], rhosumry$rhoqtls[5, yi]))
      }
#      for (yi in 1:length(rhosumry$rho.ass)) cat(sprintf(" %5.3f", rhosumry$rho.ass[yi]))
#      cat("\nFitted rho:  ")
#      for (yi in 1:length(rhosumry$rho.obs)) cat(sprintf(" %5.3f", rhosumry$rho.obs[yi]*sum(rhosumry$rho.ass)))
      cat(paste0("\np = ", round(rhosumry$pval, 5), ' for likelihood ratio test of H0: assumed rho = true rho\n'))
      cat(paste0("Quick test of relative bias: ", round(rhosumry$quicktest, 3),'\n\n'))
      dev.new(noRStudioGD = T)
      plot(0,0, type='n', xlim=c(1,length(mydat$years))+.25*c(-1,1), ylim=range(lstat), axes=F, xlab='Year', ylab = '\u03bb')
      axis(1, at=1:length(mydat$years), lab=mydat$years)
      axis(2)
      box()
      for (i in 1:length(mydat$years)){
#        lines(i+.15*c(-1,1),rep(mean(postL.ab(mydat$X[i],mydat$Ba[i], mydat$Bb[i], 0.001)$CI),2), lwd=2)
#        points(i,mean(postL.ab(mydat$X[i],mydat$Ba[i], mydat$Bb[i], 0.001)$CI), pch='+', cex=1.5)
#        lines(rep(i,2), lstat$CI[i,1:2], lwd=2)
        qtls<-postL.sumry(mydat$X[i], mydat$Ba[i], mydat$Bb[i])
        polygon(i+.15*c(-1, 1, 1, -1), qtls[c(2,2,4,4)], col=NA)
        lines(i+0.15*c(-1,1), qtls[c(3,3)])
        lines(rep(i,2), qtls[1:2])
        lines(rep(i,2), qtls[4:5])
      }
      title("Annual posterior median \u03bb with IQR and 95% CI")
      cat(paste0(unlist(rep("=",80)),collapse=''))
      cat("\nInput\n")
      cat("Year (or period)  rho    X    Ba     Bb   ghat    95% CI\n")
      for (i in 1:nyr) cat(sprintf("%-12s     %5.3f %3.0f %6.4g %6.4g %5.3f [%5.3f, %5.3f]\n",mydat$years[i], rel_wt[i], X[i], pBa0[i], pBb0[i], pBa0[i]/(pBa0[i]+pBb0[i]), qbeta(.025, pBa0[i], pBb0[i]), qbeta(.975, pBa0[i], pBb0[i])))
      sink()
      file.show(paste0(.Rvar$datadir, "/output"),delete.file=T,title=paste0("Mortality over ", yrs, " years"))
    } else if (mydat$Mtype=="P"){
    # Projections:
      yrs<-length(mydat$X) # number of years of previous data
      nyr<-mydat$nyr
      # Tau    -- total permitted take over the course of the permit
      # nyr    -- length of the permit (or the number of years to show on the graph)
      # rho    -- future rate relative to this year (actual)
      #wt<-rho # future rate relative to this year (assumed)...EoA version 1.1 assumes wt = rho
      # alev   -- alpha level for testing long-term trigger
      alev<-1-mydat$crlev
      Tau<-mydat$Tau
      nsim<-10000
      ############### -------------> HERE:
      if (mydat$Ptype=="I") {
        g<-mydat$Ba[yrs]/(mydat$Ba[yrs]+mydat$Bb[yrs])
        s2<-mydat$Ba[yrs]*mydat$Bb[yrs]/((mydat$Ba[yrs]+mydat$Bb[yrs])^2*((mydat$Ba[yrs]+mydat$Bb[yrs])+1))
        gmin<-g-2*sqrt(s2); gmax<-g+2*sqrt(s2)
        gfut<-cbind(rep(g,nyr-yrs),rep(gmin,nyr-yrs),rep(gmax,nyr-yrs))
        rho<-rep(mydat$rel_wt[yrs],nyr-yrs)
      }
      if (mydat$Ptype=="C"){
        gfut<-cbind(rep(mydat$g, nyr-yrs),rep(mydat$glwr, nyr-yrs),rep(mydat$gupr, nyr-yrs))
        rho<-numeric(nyr-yrs)+mydat$prho
      }
      if (mydat$Ptype=="V"){
        gfut<-cbind(mydat$projg,mydat$projglwr,mydat$projgupr)
        rho<-mydat$projrho
      }
      gcum<-array(dim=c(nyr,3)) # cumulative detection probability (with beta parameters)
      mu<-mydat$Ba/(mydat$Ba+mydat$Bb)
      s2<-mydat$Ba*mydat$Bb/((mydat$Ba+mydat$Bb)^2*(mydat$Ba+mydat$Bb+1))
      minp<-mu-2*sqrt(s2); maxp<-mu+2*sqrt(s2)
      gproj<-rbind(cbind(mu,minp,maxp),gfut) ##################
      #### initial years...

      X<-mydat$X;  x<-sum(X)# NOTE:
      a<-mydat$rel_wt/sum(mydat$rel_wt)
      Mest.yi<-numeric(yrs) # estimated M for the given alev for each year
      nq<-200 # number of quantiles to display from the posterior of M (going from 5th to 95th)
      Mstar<-array(dim=c(yrs,nq)) # quantiles of posterior from 5th to 95th
      aseq<-seq(.05,.95,length=nq)
      sig2<-((maxp-minp)/4)^2 # estimated means and uncertainties for g
      Eg<-cumsum(mu*a)/cumsum(a) # average, expected G is the weighted average of the observed g's
      .Rvar$Eg<-Eg
      gcum[1:yrs,1]<-Eg
#      lambda<-pmax(cumsum(X), 0.5)/Eg
      Vg<-numeric(yrs)
      for (i in 1:yrs){
        aa<-a[1:i]/sum(a[1:i])
        Vg[i]<-sum(sig2[1:i]*aa[1:i]^2)#+(sum(aa[1:i]*(mu[1:i]^2+sig2[1:i]))-Eg[i]^2)/lambda[i]
      }
      pBa<-Eg^2/Vg*(1-Eg)-Eg; pBb<-pBa*(1/Eg-1) # these are the two shape parameters for the beta distribution underlying the beta-binomial for X | M
      .Rvar$pBa<-pBa; .Rvar$pBb<-pBb
      gcum[1:yrs,2]<-pBa; gcum[1:yrs,3]<-pBb
      # M needs to be calculated up to a reasonable maximum, but the final distribution is clipped at P(M > m) < 0.0001
      mmax<-NA # stores the variable for use outside the loop
      dev.new(width = 7*.Rvar$charSizeAdjust, height = .75) # figure for watching progress of calculations
      par(mar=c(0,0,2,0))
      plot(0,0,type='n',xlim=c(0,nyr),xaxs='i',ylim=c(0,1),yaxs='i',axes=F,xlab='years')
      title('Calculating...',adj=0)
      for (yi in 1:yrs){
        mmax<-ifelse (Vg[yi] < 0.00001, fmmax(sum(X[1:yi]),.Rvar$Eg[yi]),fmmax.ab(sum(X[1:yi]),.Rvar$pBa[yi],.Rvar$pBb[yi]))
        x<-sum(X[1:yi])
        M<-x:mmax
      #########
         # posterior for M | x
          if (Vg[yi]>0.00001){
      #      pBa<-Eg^2/Vg*(1-Eg)-Eg; pBb<-pBa*(1/Eg-1) # these are the two shape parameters for the beta distribution underlying the beta-binomial for X | M
            pXgM<-VGAM::dbetabinom.ab(x,size=M,shape1=.Rvar$pBa[yi],shape2=.Rvar$pBb[yi]) # the probabilities of X for M = 0:mmax
          } else {
            pXgM<-dbinom(x,size=M,prob=.Rvar$Eg[yi]) # the probabilities of X for M = x:mmax
          }

          pM<-diff((sqrt(c(x-1,x:mmax)+1)))
          pMgX<-pXgM*pM; pMgX<-pMgX/sum(pMgX) # posterior distribution for M (ignoring M < x, which has probability = zero)
      #########
        for (ai in 1:nq) Mstar[yi,ai]<-M[min(which(1-cumsum(pMgX)<1-aseq[ai]))]
        Mest.yi[yi]<-M[min(which(1-cumsum(pMgX)<alev))]
        polygon(yi+c(-1,-1,0,0),c(0,1,1,0),col=7) # showing progress of calculations
      }
      x<-sum(X)
      M0<-M
      pMgX0<-pMgX
      ##### projection: future fatality is projected from the posterior on lambda (scaled to annual rate for most recent year)
      #### M values to be based on simulation, drawing from posterior of M for first yrs + Poisson draws from random draws from lambda posterior in future years
      ###1. generate random lambdas from projected distribution
      ##a. initial years (use posterior for lambda thus far)
################# MUST use Jeffreys prior, NOT uniform
      # calculate posterior for lambda
      ctprob<-0.0001 # A Jeffreys prior is used for the distribution of lambda. What is a reasonable maximum (and minimum) for lambda?
      # Pick a mmax = minimum M value such that P(X <= x | M) < ctprob. Then, define the Lmax to be minimum such that the P(lambda > Lmax) < ctprob.
      # p is assumed to be distributed as a beta RV with mean = pmean and variance = ((maxp-minp)/4)^2.
      # This determines beta parameters as follows:
#      Lest<-max(x, 0.5)/Eg
#      gmin<-gmax<-0
      Vg<-sum(sig2*a^2)#+(sum(a*(mu^2+sig2))-Eg^2)/Lest
      if (Vg > 0.000001){
        mmax<-fmmax.ab(x,.Rvar$pBa[yrs],.Rvar$pBb[yrs])
        if (ppois(mmax,100)<ctprob){
          lvec<-seq(0,100,by=0.01)
          Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))] # this is the biggest that lambda could credibly be and still have the biggest M conceivable
        } else {
          i<-1
          while (1){
            i<-i+1
            if (ppois(mmax,100*i)<ctprob){
              lvec<-seq((i-1)*100,i*100,by=0.01)
              Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))]
              break
            }
          }
        }
      } else {
        mmax<-fmmax(x,.Rvar$Eg[yrs])
        if (ppois(mmax,100)<ctprob){
          lvec<-seq(0,100,by=0.01)
          Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))] # this is the biggest that lambda could credibly be and still have the biggest M conceivable
        } else {
          i<-1
          while (1){
            i<-i+1
            if (ppois(mmax,100*i)<ctprob){
              lvec<-seq((i-1)*100,i*100,by=0.01)
              Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))]
              break
            }
          }
        }
      }
      L<-seq(.0001,Lmax, length=100)
      if (Vg > 0.000001){
        cpLgX<-function(L,x,pBa,pBb,mmax) 1-pLgX.ab(L,x,pBa,pBb,mmax)
        meanL<-integrate(Vectorize(cpLgX,"L"),lower=0.0001,upper=Lmax,x=x, pBa = .Rvar$pBa[yrs],pBb = .Rvar$pBb[yrs],mmax=mmax)$val/sum(mydat$rel_wt)#styr
        CI<-c(optimize(f=optqL.ab,interval=c(0.00001,Lmax),x=x,pBa = .Rvar$pBa[yrs],pBb = .Rvar$pBb[yrs],mmax=mmax,p=.025)$minimum, optimize(f=optqL.ab,interval=c(0.00001,Lmax),x=x,.Rvar$pBa[yrs],.Rvar$pBb[yrs],mmax=mmax,p=.975)$minimum)/sum(mydat$rel_wt)#styr
        pL<-Vectorize(pLgX.ab, "L")(L, x, .Rvar$pBa[yrs], .Rvar$pBb[yrs], mmax)
      } else {
        cpLgX<-function(L,x,g,mmax) 1-pLgX(L,x,g,mmax)
        meanL<-integrate(Vectorize(cpLgX,"L"),lower=0.0001,upper=Lmax,x=x,g=.Rvar$Eg[yrs],mmax=mmax)$val/sum(mydat$rel_wt)#styr
        CI<-c(optimize(f=optqL,interval=c(0.00001,Lmax),x=x,g=.Rvar$Eg[yrs],mmax=mmax,p=.025)$minimum, optimize(f=optqL,interval=c(0.00001,Lmax),x=x,g=.Rvar$Eg[yrs],mmax=mmax,p=.975)$minimum)/sum(mydat$rel_wt)#styr
        pL<-Vectorize(pLgX, "L")(L, x, .Rvar$Eg[yrs], mmax)
      }
      # generate carcasses for initial years
      Msim<-M[findInterval(runif(nsim),c(0,cumsum(pMgX)))] #simulated number of fatalities through the first yrs years
      # generate carcasses in future years
        u<-runif(nsim) # generate a new sequence of random variables because future years are a separate process from the past years (reusing the RVs --> variance too high)
      # sample lambda for projection
      # finding exact values of sampled L is very slow: optimize(f=optqL.ab,interval=c(0.00001,Lmax),x=x,pBa = pBa[yrs],pBb = pBb[yrs],mmax=mmax,p=runif(1))$minimum
      # quicker would be to do linear interpolation OR to sample M as a negative binomial (details?)
      Li<-seq(0, Lmax, by = Lmax/1000)
      pLi<-Vectorize(pLgX.ab,"L")(Li, x, .Rvar$pBa[length(pBa)], .Rvar$pBb[length(pBb)], mmax); pLi<-pLi/max(pLi)
      u<-runif(nsim)
      ind<-findInterval(u, pLi)
      Lsim<-(Li[ind]+(Li[ind+1]-Li[ind])*(pLi[ind+1]-u)/(pLi[ind+1]-pLi[ind]))/sum(mydat$rel_wt) # random draw from L (assuming rho = 1)
      if (yrs < nyr){
        Mproj<-array(rpois(nsim*(nyr-yrs),as.vector(rho%o%Lsim)),dim=c(nyr-yrs,nsim)) # rho is the change in fatality rate in a future year compared with a rho = 1 case
        # note: the Mproj formula works for both scalar rho (gtype = 1 or 2) and vector rho (gtype = 3)
      }
      # what about the future estimated fatalities?
      # generate annual x's based on M's and a schedule of search
      # assume remaining years are like the last year of data
      mug<-gfut[,1]
      sig2g<-((gfut[,3]-gfut[,2])/4)^2
      betaa<-mug^2/sig2g*(1-mug)-mug; betab<-betaa*(1/mug-1)
      Xproj<-apply(array(rbinom(length(Mproj),size=as.vector(Mproj),prob=rbeta(nsim, betaa, betab)),dim=c(nyr-yrs,nsim)), F = cumsum, M = 2) + sum(X) # total cumulative X through year i
      if (nyr-yrs == 1) Xproj<-array(Xproj, dim=c(1, length(Xproj)))
      # predicting estimated take: year-by-year calculate g and X [each year and simulation run --> posterior --> 50th, 80th, 90th percentiles]
      Mhatproj<-array(dim=c(nyr-yrs,nsim)) # for each year and simulation run, the projected Mstar is given (for user-defined alpha = alev)
      a<-c(mydat$rel_wt,rho)
      for (yi in 1:(nyr-yrs)){
        aa<-a[1:(yrs+yi)]; aa<-aa/sum(aa)
        sig2<-((gproj[1:(yrs+yi),3]-gproj[1:(yrs+yi),2])/4)^2 # estimated means and uncertainties for g
        Eg<-sum(gproj[1:(yrs+yi),1]*aa) # average, expected G is the weighted average of the observed g's
        for (xx in unique(Xproj[yi,])){
          # calculate g (for the cumulative)
#          lambda<-max(x, 0.5)/Eg
          Vg<-sum(sig2*aa^2)#+(sum(aa*(gproj[1:(yrs+yi),1]^2+sig2))-Eg^2)/lambda
          pBa<-Eg^2/Vg*(1-Eg)-Eg; pBb<-pBa*(1/Eg-1) # these are the two shape parameters for the beta distribution underlying the beta-binomial for X | M
          # calculate posterior of M
          mmax<-ifelse (Vg < 0.00001, fmmax(xx,Eg),fmmax.ab(xx,pBa,pBb))
          M<-xx:mmax
      #    pM<-rep(1/length(M),length(M))
      #    if (Vg >= 0.00001){
      #      pXgM<-VGAM::dbetabinom.ab(x,size=M,shape1=pBa,shape2=pBb) # the probabilities of X for M = 0:mmax
      #    } else {
      #      pXgM<-dbinom(x,size=M,prob=Eg) # the probabilities of X for M = 0:mmax
      #    }
      #    pMgX<-numeric(length(M))
      #    pMgX<-pXgM*pM/sum(pXgM*pM)
          if (Vg>0.00001){
      #      pBa<-Eg^2/Vg*(1-Eg)-Eg; pBb<-pBa*(1/Eg-1) # these are the two shape parameters for the beta distribution underlying the beta-binomial for X | M
            pXgM<-VGAM::dbetabinom.ab(xx,size=M,shape1=pBa,shape2=pBb) # the probabilities of X for M = 0:mmax
          } else {
            pXgM<-dbinom(xx,size=M,prob=Eg) # the probabilities of X for M = xx:mmax
          }
          pM<-diff((sqrt(c(xx-1,xx:mmax)+1)))
          pMgX<-pXgM*pM; pMgX<-pMgX/sum(pMgX) # posterior distribution for M (ignoring M < xx, which has probability = zero)
          ###########
          Mhatproj[yi,Xproj[yi,]==xx]<-M[min(which(cumsum(pMgX)>1-alev))]
        }
        polygon(yrs+yi+c(-1,-1,0,0),c(0,1,1,0),col=8)
      }
      if(.Rvar$platform == "windows") bringToTop()
      dev.off() # progress bar
      ############ projection graph
      ### estimated
      ht<-7*getData('charSizeAdjust')
      dev.new(height=ht,width=2*ht, noRStudioGD = T)
      charSizeAdjust<-par('din')[1]/12
      par(mar=c(9,3,3,.5),mgp=c(2,.7,0),family='sans')
      gtop<-max(c(Tau,Mstar,quantile(Msim+apply(Mproj,F=sum,M=2),.95),quantile(Mhatproj[nyr-yrs,],.95)))
      plot(0,0,type='n',xlim=c(0,nyr)+.5,ylim=c(0,gtop),xlab='Year',ylab='Cumulative Fatalities',xaxs='i',axes=F)
      axis(1)
      axis(1,at=1:nyr,lab=F,tck=-0.01)
      if (par('usr')[4]<110){
        axis(2,at=seq(0,par('usr')[4],by=10))
        axis(2,at=seq(5,par('usr')[4],by=10),lab=F)
        axis(2,at=seq(1,par('usr')[4],by=1),lab=F, tck=-0.01)
      } else if (par('usr')[4]<150){
        axis(2,at=seq(0,par('usr')[4],by=20))
        axis(2,at=seq(10,par('usr')[4],by=10),lab=F)
        axis(2,at=seq(5,par('usr')[4],by=5),lab=F, tck=-0.01)
      } else {
        axis(2)
        if (par('usr')[4]<1100&((par('yaxp')[2]/par('yaxp')[3])%%10==0))   axis(2,at=seq(0,par('usr')[4],by=10),lab=F, tck=-0.01)
      }
      title('Cumulative mortality',adj=0,line=1.5, cex.main = charSizeAdjust * 14/12)

      ## projection...
      qcollist<-colors()[c(245,235,224)]
      tmpM<-tmpMstar<-array(dim=c(nyr-yrs, 8))
      for (yi in 1:(nyr-yrs)){
        xi<-yrs+yi+c(-.5,-.5,.5,.5)
        mMstar<-mean(Mhatproj[yi,])
      #  mMstar<-median(Mhatproj[yi,])
        qtls<-quantile(Mhatproj[yi,],probs=c(0.05,.1,.25,.5,.75,.9,0.95), type = 3)
        # predicted estimates (note: if alpha = 0.5, the estimated fatality is slightly lower than actual on average because the median underestimates...mean is closer to unbiased when M is large)
        polygon(xi,qtls[c(1,7,7,1)],col=qcollist[1],border=NA)
        polygon(xi,qtls[c(2,6,6,2)],col=qcollist[2],border=NA)
        polygon(xi,qtls[c(3,5,5,3)],col=qcollist[3],border=NA)
        lines(yrs+yi+c(-.5,.5), rep(qtls[4],2), col = qcollist[2])
      #  points(yrs+yi,qtls[4],pch=3)
        points(yrs+yi,mMstar,col=2,pch=18)
        tmpMstar[yi,]<-c(mMstar, qtls)
        # projected fatalities
        if (yi == 1) projM.yi<-Msim+Mproj[1,]
        if (yi > 1)  projM.yi<-Msim+apply(Mproj[1:yi,],F=sum, M=2)
        meanM<-mean(projM.yi)
        qtls<-quantile(projM.yi,probs=c(0.05,.1,.25,.5,.75,.9,0.95), type = 3)
        polygon(yrs+yi+0.07*c(1,1,-1,-1),qtls[c(3,5,5,3)])
        lines(yrs+yi+0.07*c(-1,1),rep(qtls[4],2))
        points(yrs+yi,meanM,pch=20)
        lines(rep(yrs+yi,2),qtls[2:3]);lines(rep(yrs+yi,2),qtls[5:6])
        lines(rep(yrs+yi,2),qtls[1:2],lty=3);lines(rep(yrs+yi,2),qtls[6:7],lty=3)
        tmpM[yi,]<-c(meanM, qtls)
      }
      .Rvar$Mhatproj.sumry<-tmpMstar
      .Rvar$Mproj.sumry<-tmpM
      # estimation...
      par(family='serif',xpd=T)
      # the data
#      colr<-c(colors()[c(552,498,652)],"#FFFFE7")
      colr<-c(colors()[c(400, 128, 31, 552, 552, 498, 652)],"#FFFFE7")
      medind<-round(mean(which(abs(aseq-.5)==min(abs(aseq-.5)))))
      for (yi in 1:yrs){
        xi<-yi+c(-.5,-.5,.5,.5)
        for (ai in 2:nq){
          polygon(xi,c(Mstar[yi,ai-1],Mstar[c(yi,yi),ai],Mstar[yi,ai-1]),col=colorRampPalette(colr)( nq )[ai],border=colorRampPalette(colr)( nq )[ai])
        }
        #if (abs(alev-0.5)> 0.01)
        points(yi,Mest.yi[yi],pch=8, cex= charSizeAdjust*.7)
        lines(yi+c(-.5,.5),rep(Mstar[yi,medind],2),col=colors()[69],lwd=1.5)
        lines(yi+c(-.5,.5),rep(Mstar[yi,nq],2),lty=3)
        lines(yi+c(-.5,.5),rep(Mstar[yi,1],2),lty=3)
      }
      lines(par('usr')[1:2],rep(Tau,2),col=4)
      lines(c(.8, 1.2), rep(Tau,2),col=colors()[1])
      par(family='serif',xpd=T)
      mtext(text = 'Posterior distribution of M:',
        side = 1, line = 2.5, at = par('usr')[1]-diff(par('usr')[1:2])*.01, adj = 0, cex = charSizeAdjust
      )
      mtext(text = 'Estimated baseline annual fatality rate (\u03bb for \u03c1 = 1):',
        side = 1, line = 6, at = par('usr')[1]+diff(par('usr')[1:2])*.12, adj = 0, cex = charSizeAdjust
      )
      mtext(text = paste0('Mean = ',signif(meanL,3),', 95% CI = [',signif(CI[1],3),', ',signif(CI[2],3),']'),
        side = 1, line = 7, adj=0,cex = charSizeAdjust, at = par('usr')[1]+diff(par('usr')[1:2])*.12,
      )

      hTot <- diff(par('usr')[3:4])/diff(par('plt')[3:4]) # total height of plot window (in usr coords)
      bottom <- par('usr')[3] - hTot*par('plt')[3]
      nq<-100
      tmp <- bottom + seq(0.15, 0.55, length = nq) * (par('usr')[3] - bottom)
      xi<-par('usr')[1]+.3*nyr/30*c(1,1,-1,-1)
      for (i in nq:1)  polygon(xi,tmp[c(1,i,i,1)],col=colorRampPalette(colr)( nq )[i],border=colorRampPalette(colr)( nq )[i])
      lines(xi[c(1,3)],rep(tmp[40],2),col=colors()[43], lwd=1.5)
#      text(xi[2]+.5*nyr/30, tmp[40],adj=0,lab=paste('M* for \u03b1 = ',alev,sep=''), cex=charSizeAdjust)
      if (abs(alev-.5)>=0.01){
        points(mean(xi[2:3]),tmp[70],pch=8, cex= charSizeAdjust*.7)
        text(xi[2]+.5*nyr/30, tmp[70],adj=0,lab=paste0("M* for 1 - \u03b1 = ",1-alev), cex=charSizeAdjust)
        text(xi[2]+.5*nyr/30, tmp[40],adj=0,lab='Median', cex=charSizeAdjust)
      } else {
        points(mean(xi[2:3]),tmp[40],pch=8, cex= charSizeAdjust*.7)
        text(xi[2]+.5*nyr/30, tmp[40],adj=0,lab='Median = M*', cex=charSizeAdjust)
      }
      lines(xi[c(1,3)], rep(max(tmp), 2),lty=3)
      text(xi[2]+.5*nyr/30, max(tmp),adj=0,lab='95th percentile',cex=charSizeAdjust)
      lines(xi[c(1,3)], rep(min(tmp), 2),lty=3)
      par(family='serif')
      text(1,Tau,'\u03a4',adj=c(0.5,.5),cex=charSizeAdjust)
      par(family='sans')
      box()
      par(xpd=T,family='serif')
      polygon(c(par('usr')[1],yrs+.5,yrs+.5,par('usr')[1]),par('usr')[4]+c(0,0,rep(0.05*diff(par('usr')[3:4]),2)),col=colors()[652])
      text((yrs+.5)/2,par('usr')[4],'estimated',adj=c(.5,-.3),cex=charSizeAdjust)
      polygon(c(par('usr')[2],yrs+.5,yrs+.5,par('usr')[2]),par('usr')[4]+c(0,0,rep(0.05*diff(par('usr')[3:4]),2)),col=qcollist[1])
      text((yrs+.5+par('usr')[2])/2,par('usr')[4],'projected',adj=c(.5,-.3),cex= charSizeAdjust)

      xsz<-nyr/30;
      ysz<-0.75*(diff(par('usr')[3:4]))/30
      top<-bottom + .6 * (par('usr')[3]-bottom)
      polygon(x=par('usr')[1]+xsz*c(25,25,27,27),y=top+ysz*c(-1,0,0,-1))
      points(x=par('usr')[1]+xsz*26.3,y=top-ysz*.5,pch=16)
      text(x=par('usr')[1]+xsz*26.3,y=top+ysz*.4,'mean',adj=c(.5,0),cex=charSizeAdjust)
      lines(x=par('usr')[1]+xsz*c(25,24),y=top+ysz*rep(-.5,2))
      lines(x=par('usr')[1]+xsz*c(23,24),y=top+ysz*rep(-.5,2),lty=3)
      lines(par('usr')[1]+xsz*rep(26,2),top+ysz*c(-1,0))
      lines(x=par('usr')[1]+xsz*c(27,28),y=top+ysz*rep(-.5,2))
      lines(x=par('usr')[1]+xsz*c(28,29),y=top+ysz*rep(-.5,2),lty=3)
      par(family='serif')
      text(par('usr')[1]+xsz*23,top-ysz,'5th',adj=c(0.5,1.4),cex=charSizeAdjust)
      text(par('usr')[1]+xsz*24,top-ysz,'10th',adj=c(0.5,1.4),cex=charSizeAdjust)
      text(par('usr')[1]+xsz*25,top-ysz,'25th',adj=c(0.5,1.4),cex=charSizeAdjust)
      text(par('usr')[1]+xsz*26,top-ysz,'50th',adj=c(0.5,1.4),cex=charSizeAdjust)
      text(par('usr')[1]+xsz*27,top-ysz,'75th',adj=c(0.5,1.4),cex=charSizeAdjust)
      text(par('usr')[1]+xsz*28,top-ysz,'90th',adj=c(0.5,1.4),cex=charSizeAdjust)
      text(par('usr')[1]+xsz*29,top-ysz,'95th',adj=c(0.5,1.4),cex=charSizeAdjust)
      text(par('usr')[1]+xsz*22.5,top-ysz,'Projected mortality (M)',adj=c(1,0),family='serif',cex=charSizeAdjust)
      text(par('usr')[1]+xsz*22.5,top-ysz,'percentiles:',adj=c(1,1.4),family='serif',cex=charSizeAdjust)
      top<-bottom + .2 * (par('usr')[3]-bottom)
      polygon(x=par('usr')[1]+xsz*c(23,23,29,29),y=top-ysz*c(1,-1,-1,1),col=qcollist[1],border=NA)
      polygon(x=par('usr')[1]+xsz*c(24,24,28,28),y=top-ysz*c(1,-1,-1,1),col=qcollist[2],border=NA)
      polygon(x=par('usr')[1]+xsz*c(25,25,27,27),y=top-ysz*c(1,-1,-1,1),col=qcollist[3],border=NA)
      lines(x=par('usr')[1]+xsz*c(26,26),y=top-ysz*c(1,-1), col=qcollist[2], lwd=1.5)
      text(par('usr')[1]+xsz*22.5,top-ysz,paste('Projected mortality estimates (M* for 1 - \u03b1 = ',1-alev,')\nwith median (line) and 50%, 80%, and 90% PIs',sep=''),adj=c(1,0),family='serif',cex=charSizeAdjust)
      points(par('usr')[1]+xsz*26,top,pch=16,col=2)
      text(par('usr')[1]+xsz*26,top + ysz, 'mean', cex=charSizeAdjust)
      par(xpd=F)
      par(family='sans')
      if(.Rvar$platform == "windows") bringToTop()
      if (max(Mest.yi)<=Tau){ # no triggering in initial years
      # CDF graph of fraction of projects triggering and fraction of projects exceeding...
      # fraction of projects triggering within yi years
        dev.new(height=4, width=6, noRStudioGD = T)
        par (mar=c(3,3,1,1), mgp=c(2,.7,0), tck=-.01)
        tmp<-(yrs+(1:(nyr-yrs)))*(Mhatproj>Tau)
        tmp[tmp==0]<-nyr+1
        junk<-apply(tmp,F=min,M=2)
#        pTrig<-cumsum(tabulate(junk))[1:nyr]/nsim
        pTrig<-1+numeric(nyr)
        pTrig[1:min(max(junk), nyr)]<-cumsum(tabulate(junk))[1:min(max(junk), nyr)]/nsim
        pExceed<-numeric(nyr)
        # generate carcasses according to lambda for the initial years
        Mi<-array(apply(array(rpois(nsim*yrs,as.vector(mydat$rel_wt%o%Lsim)),dim=c(yrs,nsim)),F=cumsum,M=2), dim=c(yrs, nsim))
        for (yi in 1:min(yrs, nyr)){
          pExceed[yi]<-sum(Mi[yi,]>Tau)/nsim
        }
        Mtot.yi<-apply(Mproj,F=cumsum,M=2)
        if (nyr-yrs == 1) Mtot.yi<-array(dim=c(1, length(Mtot.yi)))
        if (nyr>yrs){
          for (yi in 1:(nyr-yrs)){
            pExceed[yi+yrs]<-sum(Mtot.yi[yi,]+Msim>Tau)/nsim
          }
        }
        plot(0,0,type='n',xlim=c(1,nyr),ylim=c(0,1),ylab='Posterior predictive probability',xlab='Year')
        lines(1:nyr,pExceed,type='s',lwd=2)
        lines(1:nyr,pTrig,type='s',col=2)
        if (pTrig[round(nyr/2)]>.8) loc<-"bottomright" else loc<-"topleft"
        legend(x=loc,legend=c("p(M > Tau) [exceedance]","p(M* > Tau) [triggering]"),lty=c(1,1),col=1:2,lwd=2:1)
      }
    #### write parameter set and summary results to a data file
      while(1){ if (sink.number()==0) break else sink() }
      sink(paste0(.Rvar$datadir, "/output"))
      cat("\n================================================================================\n")
      cat(paste("Summary statistics from posterior predictive distributions for ",nsim, " simulated projects\n",sep=''))
      cat("------------------------------------\n")
      cat(paste("Estimated annual baseline fatality rate (lambda for rho = 1): mean = ", signif(meanL,3),", 95% CI = [",signif(CI[1],3),", ",signif(CI[2],3),"]\n",sep=''))
      if (max(Mest.yi)<=Tau){ # Why cut this routine out if max(Mest.yi) > Tau
        cat("\nProjected fatalities and fatality estimates...\n")
        cat(paste("p(M > Tau within ",nyr," years) = ",sum(apply(Mproj,M=2,F=sum)+Msim>Tau)/nsim,"  [exceedance]\n",sep=""))
        trigYr<-(Mhatproj>Tau)*(1:(nyr-yrs)); trigYr[trigYr==0]<-nyr-yrs+1; trigYr<-apply(trigYr,F=min,M=2)# years into projection years (i.e., years after past data were collected)
        cat(paste("p(M* > Tau within ",nyr," years) = ",sum(trigYr<=nyr-yrs)/nsim,"  [triggering]\n",sep=''))
        cat(paste0("M* based on credibility level 1 - alpha = ", mydat$crlev, "\n"))
        Mtot<-numeric(nsim)
        for (simi in 1:nsim) Mtot[simi]<-Msim[simi]+sum(Mproj[1:min(trigYr[simi],nyr-yrs),simi]) # total fatality at time of triggering or end of project
        qtls<-quantile(Mtot[trigYr<=nyr-yrs],c(0.25,0.5,0.75), type = 3) # quantiles of totals in projects that trigger
        cat(sprintf("\nAmong projects with triggering (%1.2f%%), mean(M) = %1.2f at time of triggering, with median = %g and IQR = [%g, %g]\n",sum(Mhatproj[nyr-yrs,]>Tau)/nsim*100,mean(Mtot[trigYr<=nyr-yrs]),qtls[2],qtls[1],qtls[3]))
        qtls<-quantile(Mtot[trigYr>nyr-yrs],c(0.25,0.5,0.75), type = 3) # quantiles of totals in projects that trigger
        cat(sprintf("Among projects with no triggering (%1.2f%%), mean(M) = %1.2f at end of %i years, with median = %g and IQR = [%g, %g]\n",sum(Mhatproj[nyr-yrs,]<=Tau)/nsim*100,mean(Mtot[trigYr>nyr-yrs]),nyr,qtls[2],qtls[1],qtls[3]))
        cat('\n')
        cat('Years of operations without triggering:\n')
        tsum<-trigYr; tsum[tsum>nyr-yrs]<-nyr-yrs
        qtls<-quantile(tsum,c(0.25,.5,0.75), type = 3)+yrs
        cat(sprintf(' Mean = %0.2f, with median = % g and IQR = [%g, %g]\n',yrs+mean(tsum),qtls[2],qtls[1],qtls[3]))
      } else {
        cat(paste0("\nTriggering in year ", mydat$years[min(which(Mest.yi>Tau))], '\n'))
      }
      cat('\n----------------------------------------')
      cat("\nSummary statistics for projection years\n")
      cat('----------------------------------------------------------------------------------------------------\n')
      cat(sprintf("Yr   Mean             %-39s| %-40s\n", "quantiles of M", "quantiles of M*"))
      cat("        M       M*  0.05  0.10  0.25  0.50  0.75  0.90  0.95 |  0.05  0.10  0.25  0.50  0.75  0.90  0.95\n")
      cat('----------------------------------------------------------------------------------------------------\n')
      m1<-.Rvar$Mproj.sumry; m2<-.Rvar$Mhatproj.sumry
      for (yi in 1:(nyr-yrs)){
        cat(sprintf("%-4.0f%7.1f%7.1f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f |%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f%6.0f\n",
        yi, signif(m1[yi,1],4), signif(m2[yi,1],4), m1[yi,2], m1[yi,3], m1[yi,4], m1[yi,5], m1[yi,6], m1[yi,7], m1[yi,8], m2[yi,2], m2[yi,3], m2[yi,4], m2[yi,5], m2[yi,6], m2[yi,7], m2[yi,8] ))
      }
      cat("\n================================================================================\n")
      cat(paste("\nGoverning parameters: Tau = ",Tau,", alpha = ",alev,"\n\n",sep=''))
      cat(paste("Data for ",yrs," years of monitoring:\n",sep=''))
      cat("            yr    x    g    glwr   gupr  rho  M*\n")
      for (yi in 1:yrs){
        cat(sprintf("%15s %3.0f %.4f %.4f %.4f  %.3g %3.0f\n",mydat$years[yi],mydat$X[yi],mu[yi],minp[yi],maxp[yi],mydat$rel_wt[yi],Mest.yi[yi]))
      }
      cat('\n')
      cat("Parameters for future monitoring and operations:\n")
      cat(paste("  Number of years: ",nyr-yrs,'\n'))
      if (mydat$Ptype=="I" | mydat$Ptype=="C"){
        cat(sprintf("  g = %.4g, 95%% CI  [%.4g, %.4g]\n",gfut[1,1],gfut[1,2],gfut[1,3]))
        cat(sprintf("  Relative weight (rho): %.3g",rho[1])) # rho is denormalized to match format of user input
      } else if (mydat$Ptype=='V'){
        cat(" yr    g      glwr   gupr   rho\n")
        for (yi in 1:(nyr-yrs)){
          cat(sprintf("%3.0f  %.4f  %.4f %.4f  %.3g\n",yi+yrs,gfut[yi,1],gfut[yi,2],gfut[yi,3],rho[yi]))#*mydat$rel_wt[yrs])) # rho is denormalized to match format of user input
        }
      } else {
        tkmessageBox(message="Bad Ptype = ", mydat$Ptype, " in myCalc. Aborting calculation")
        return(F)
      }
      cat("\n***************************************************************************************\n")
#        ############## past results
        ctprob<-0.0001 # A Jeffreys prior is used for the distribution of lambda. What is a reasonable maximum (and minimum) for lambda?
        # Pick a mmax = minimum M value such that P(X <= x | M) < ctprob. Then, define the Lmax to be minimum such that the P(lambda > Lmax) < ctprob.
        # p is assumed to be distributed as a beta RV with mean = pmean and variance = ((maxp-minp)/4)^2.
        # This determines beta parameters as follows:
  #      Lest<-max(x, 0.5)/Eg
        gmin<-gmax<-0
        Vg<-sum(sig2*a^2)#+(sum(a*(mu^2+sig2))-Eg^2)/Lest
        if (Vg > 0.000001){
          mmax<-fmmax.ab(x,.Rvar$pBa[yrs],.Rvar$pBb[yrs])
          if (ppois(mmax,100)<ctprob){
            lvec<-seq(0,100,by=0.01)
            Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))] # this is the biggest that lambda could credibly be and still have the biggest M conceivable
          } else {
            i<-1
            while (1){
              i<-i+1
              if (ppois(mmax,100*i)<ctprob){
                lvec<-seq((i-1)*100,i*100,by=0.01)
                Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))]
                break
              }
            }
          }
        } else {
          mmax<-fmmax(x,.Rvar$Eg[yrs])
          if (ppois(mmax,100)<ctprob){
            lvec<-seq(0,100,by=0.01)
            Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))] # this is the biggest that lambda could credibly be and still have the biggest M conceivable
          } else {
            i<-1
            while (1){
              i<-i+1
              if (ppois(mmax,100*i)<ctprob){
                lvec<-seq((i-1)*100,i*100,by=0.01)
                Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))]
                break
              }
            }
          }
        }
        if (Vg > 0.000001){
          cpLgX<-function(L,x,pBa,pBb,mmax) 1-pLgX.ab(L,x,pBa,pBb,mmax)
          meanL<-integrate(Vectorize(cpLgX,"L"),lower=0.0001,upper=Lmax,x=x,.Rvar$pBa[yrs],.Rvar$pBb[yrs],mmax=mmax)$val/sum(mydat$rel_wt)#styr
          CI<-c(optimize(f=optqL.ab,interval=c(0.00001,Lmax),x=x,.Rvar$pBa[yrs],.Rvar$pBb[yrs],mmax=mmax,p=.025)$minimum, optimize(f=optqL.ab,interval=c(0.00001,Lmax),x=x,.Rvar$pBa[yrs],.Rvar$pBb[yrs],mmax=mmax,p=.975)$minimum)/sum(mydat$rel_wt)#styr
        } else {
          cpLgX<-function(L,x,g,mmax) 1-pLgX(L,x,g,mmax)
          meanL<-integrate(Vectorize(cpLgX,"L"),lower=0.0001,upper=Lmax,x=x,g=.Rvar$Eg[yrs],mmax=mmax)$val/sum(mydat$rel_wt)#styr
          CI<-c(optimize(f=optqL,interval=c(0.00001,Lmax),x=x,g=.Rvar$Eg[yrs],mmax=mmax,p=.025)$minimum, optimize(f=optqL,interval=c(0.00001,Lmax),x=x,g=.Rvar$Eg[yrs],mmax=mmax,p=.975)$minimum)/sum(mydat$rel_wt)#styr
        }
        cat(paste("Summary statistics for mortality estimates through ", length(mydat$years), " years\n",sep=''))
        cat("------------------------------------\n")
        cat("Results\n")
        cat(paste0("Totals through ", length(mydat$years), " years\n"))
  #$$$$$$$$$$$$$$$
        yrs<-length(mydat$years) # number of years of monitoring
        nyr<-mydat$nyr
        gcum<-array(dim=c(yrs,3)) # cumulative detection probability (with beta parameters)
        Mest.yi<-numeric(yrs) # estimated M for the given aM for each year
        #### initial years
        ### The function calculates an overall average g along with estimated variance.
        ### Posterior is betabinomial on g and is calculated based on a uniform prior.
        X<-mydat$X; nyr<-length(X)
        x<-sum(X)
        pBa0<-mydat$Ba; pBb0<-mydat$Bb
        a<-mydat$rel_wt/sum(mydat$rel_wt)
        rel_wt<-mydat$rel_wt
        # p is assumed to be distributed as a beta RV with mean = pmean and variance = ((maxp-minp)/4)^2.
        # This determines beta parameters as follows:
        aM<-1-mydat$crlev
        Tau<-mydat$Tau
        mu<-pBa0/(pBa0+pBb0)
        sig2<-pBa0*pBb0/((pBa0+pBb0)^2*(pBa0+pBb0+1))
        Eg<-sum(mu*a) # average, expected G is the weighted average of the observed g's
  #      lambda<-max(sum(X), 0.5)/Eg
        Vg<-sum(sig2*a^2)#+(sum(a*(mu^2+sig2))-Eg^2)/lambda
        pBa<-Eg^2/Vg*(1-Eg)-Eg; pBb<-pBa*(1/Eg-1) # these are the two shape parameters for the beta distribution underlying the beta-binomial for X | M
        pMgX<-postM.ab(x, pBa, pBb)
        cims<-MCI(pMgX, .95)
        cat('\n')
        if (sides == 2){
          cat(paste0(100*(1-aM),"% CI for M = [",cims[1],", ",  cims[2],"]\n"))
        } else {
          tmpmst<-calcMstar(pMgX, aM)
          cat(paste0("M* = ",tmpmst, " for 1 - alpha = ",1 - aM,", i.e., P(M <= ", tmpmst ,") >= ", (1-aM)*100, "%", '\n'))
        }
        cat(paste0("Estimated overall detection probability: ",
           "g = ", signif(Eg, 3),", 95% CI = [",signif(qbeta(.025, pBa, pBb),3),", ",signif(qbeta(.975, pBa, pBb),3),"]\n ",
           "   Ba = ", signif(pBa, 5), ", Bb = ", signif(pBb, 5), "\n"))
  #      CI<-c(L[min(which(pL>(1-aCI)/2))]/sum(mydat$rel_wt),L[min(which(pL>aCI+(1-aCI)/2))]/sum(mydat$rel_wt))
  #      CI<-c(L[min(which(pL>.025))]/sum(mydat$rel_wt),L[min(which(pL>.975))]/sum(mydat$rel_wt))
  #      par(family='serif')
        cat(paste0("Estimated baseline fatality rate (for rho = 1): lambda = ",signif(meanL,4),", 95% CI = [",signif(CI[1],3),", ",signif(CI[2],3),"]\n"))
        cat("\n")

        cat(paste0("Cumulative Mortality Estimates\n"))
        cat("Year         M*   median   95% CI   mean(lambda) 95% CI\n")
        for (yi in 1:length(mydat$years)){
          postmi<-postM.ab(sum(mydat$X[1:yi]), .Rvar$pBa[yi], .Rvar$pBb[yi])
          mcii<-MCI(postmi, .95)
          lstat<-postL.abCI(sum(mydat$X[1:yi]), .Rvar$pBa[yi], .Rvar$pBb[yi])
          cat(sprintf("%-12s  %-5.0f %-7.0f[%.0f, %.0f]  %6.4f   [%6.4g, %6.4g]\n",
            mydat$years[yi],
          .Rvar$myTrack[yi,1],
          .Rvar$myTrack[yi,6],
            mcii[1],
            mcii[2],
            signif(lstat$meanL,4),
            signif(lstat$CI[1],4),
            signif(lstat$CI[2],4)))
        }

        cat(paste0("\nAnnual Mortality Estimates\n"))
        cat("Year         M*   median   95% CI   mean(lambda) 95% CI\n")
        for (i in 1:length(mydat$years)){
          postmi<-postM.ab(mydat$X[i], mydat$Ba[i], mydat$Bb[i])
          mcii<-MCI(postmi, .95)
          lstat<-postL.abCI(mydat$X[i], mydat$Ba[i], mydat$Bb[i])
          cat(sprintf("%-12s  %-5.0f %-7.0f[%.0f, %.0f]  %6.4f   [%6.4g, %6.4g]\n",
            mydat$years[i],
            calcMstar(postmi, 1 - mydat$crlev),
            calcMstar(postmi, 0.5),
            mcii[1],
            mcii[2],
            signif(lstat$meanL,4),
            signif(lstat$CI[1],4),
            signif(lstat$CI[2],4)))
        }
        cat(paste0("\n\nTest of assumed relative weights (rho) and potential bias"))
        rhosumry<-rhotest(mydat)
        cat("              Fitted rho\n")
        cat("Assumed rho     95% CI \n")
        for (yi in 1:length(rhosumry$rho.ass)) {
          cat(sprintf(
            " %6.3g     [%5.3f, %5.3f]\n", rhosumry$rho.ass[yi], rhosumry$rhoqtls[1, yi], rhosumry$rhoqtls[5, yi]))
        }
#        cat("\nAssumed rho: ")
#        for (yi in 1:length(rhosumry$rho.ass)) cat(sprintf(" %5.3f", rhosumry$rho.ass[yi]))
#        cat("\nFitted rho:  ")
#        for (yi in 1:length(rhosumry$rho.obs)) cat(sprintf(" %5.3f", rhosumry$rho.obs[yi]*sum(rhosumry$rho.ass)))
        cat(paste0("\np = ", round(rhosumry$pval, 5), ' for likelihood ratio test of H0: assumed rho = true rho\n'))
        cat(paste0("Quick test of relative bias: ", round(rhosumry$quicktest, 3),'\n\n'))
#        cat(paste0(", significance level = ", round(rhosumry$qtpval, 5)),'\n\n')
        cat(paste0(unlist(rep("=",80)),collapse=''))
        cat("\nInput\n")
        cat("Year (or period) rel_wt  X    Ba     Bb   ghat    95% CI\n")
        for (i in 1:nyr) cat(sprintf("%-12s     %5.3f %3.0f %6.4g %6.4g %5.3f [%5.3f, %5.3f]\n",mydat$years[i], rel_wt[i], X[i], pBa0[i], pBb0[i], pBa0[i]/(pBa0[i]+pBb0[i]), qbeta(.025, pBa0[i], pBb0[i]), qbeta(.975, pBa0[i], pBb0[i])))
#        ##############

        sink()
        file.show(paste0(.Rvar$datadir, "/output"),delete.file=T,title="Summary of projected fatality and triggering")
        ##########
      }
    } else if (mydat$option=="M" & mydat$Mtype=="C"){# estimate current cumulative carcass count-based total fatalities
      yrs<-length(mydat$years) # number of years of monitoring
      nyr<-mydat$nyr
      gcum<-array(dim=c(yrs,3)) # cumulative detection probability (with beta parameters)
      Mest.yi<-numeric(yrs) # estimated M for the given aM for each year
      #### initial years
      ### The function calculates an overall average g along with estimated variance.
      ### Posterior is betabinomial on g and is calculated based on a uniform prior.
      X<-mydat$X; nyr<-length(X)
      x<-sum(X)
      pBa0<-mydat$Ba; pBb0<-mydat$Bb
      a<-mydat$rel_wt/sum(mydat$rel_wt)
      rel_wt<-mydat$rel_wt
      # p is assumed to be distributed as a beta RV with mean = pmean and variance = ((maxp-minp)/4)^2.
      # This determines beta parameters as follows:
      aM<-1-mydat$crlev
      Tau<-mydat$Tau
      mu<-pBa0/(pBa0+pBb0)
      sig2<-pBa0*pBb0/((pBa0+pBb0)^2*(pBa0+pBb0+1))
      Eg<-sum(mu*a) # average, expected G is the weighted average of the observed g's
#      lambda<-max(sum(X), 0.5)/Eg
      Vg<-sum(sig2*a^2)#+(sum(a*(mu^2+sig2))-Eg^2)/lambda
      pBa<-Eg^2/Vg*(1-Eg)-Eg; pBb<-pBa*(1/Eg-1) # these are the two shape parameters for the beta distribution underlying the beta-binomial for X | M
      gmin<-qbeta(.025,shape1=pBa,shape2=pBb); gmax<-qbeta(.975,shape1=pBa,shape2=pBb)
      ### M needs to be calculated up to a reasonable maximum, but the final distribution is clipped at P(M > m) < 0.0001
      mmax<-ifelse (Vg < 0.00001, fmmax(x,Eg),fmmax.ab(x,pBa,pBb))
      M<-x:mmax
      if (Vg>0.00001){
      #      pBa<-Eg^2/Vg*(1-Eg)-Eg; pBb<-pBa*(1/Eg-1) # these are the two shape parameters for the beta distribution underlying the beta-binomial for X | M
        pXgM<-VGAM::dbetabinom.ab(x,size=M,shape1=pBa,shape2=pBb) # the probabilities of X for M = 0:mmax
      } else {
        pXgM<-dbinom(x,size=M,prob=Eg) # the probabilities of X for M = x:mmax
      }
      pM<-diff((sqrt(c(x-1,x:mmax)+1)))
      pMgX<-pXgM*pM; pMgX<-pMgX/sum(pMgX) # posterior distribution for M (ignoring M < x, which has probability = zero)
      M<-c(rep(0,x), M)
      pMgX<-c(rep(0,x),pMgX)
      cs<-cumsum(pMgX)
      pMem<-pMgX[1:min(which(cs>.999))]
      pMem<-pMem/sum(pMem)
      cs<-cumsum(pMem)
      if (sides=='1') py<-c(1,1-cumsum(pMem[-length(pMem)])) else py <- pMem
      m<-1:length(py)-1
      cut.off<-min(which(py<=aM))-1
      if (cut.off > 100) {
         wid<-1
      } else if (cut.off < 50){
        wid<-3
      } else {
         wid<-2
      }
      plot(m,py,
        xlab='m',
        ylab=ifelse(sides == '1', paste0("Prob( M \u2265 m | X = ", x,")"), paste0("Prob( M = m | X = ", x,")")),
        cex.axis=.9*.Rvar$charSizeAdjust,
        cex.lab=.9*.Rvar$charSizeAdjust,
        xlim = c(ifelse(sides=='1', 0, m[min(which(cumsum(py)>0.001))]),max(m))
      )
      col.use=rep("black", length(py))
      if (sides=='1'){
         cut.off.cr<-min(which(cs>=1-aM))
         col.use[m[1]:(cut.off.cr)]='red'
         for(i in 1:length(py)){
           lines(x=rep(m[i],2), y= c(0,py[i]), col=col.use[i], lwd=wid)
         }
      } else {
         lwrbnd<-min(which(cs > aM/2))-1
         lwrArea<-ifelse(lwrbnd == 0, 0, cs[lwrbnd])
         uprbnd<-min(which(cs > 1-aM + lwrArea))-1
         reds<-lwrbnd:uprbnd+1
         col.use[reds]=colors()[90]
         for(i in 1:length(py)){
           lines(x=rep(m[i],2), y= c(0,py[i]), col=col.use[i], lwd=wid)
         }
      }
      top<-par('usr')[4]/1.04
      if(.Rvar$platform == "windows") bringToTop()
      ## estimate annual rate
      x<-sum(X)
      mmax<-ifelse (Vg< 0.00001, fmmax(x,Eg),fmmax.ab(sum(X),pBa,pBb))
      ctprob<-0.0001 # A Jeffreys prior is used for the distribution of lambda. What is a reasonable maximum (and minimum) for lambda?
      # Pick a mmax = minimum M value such that P(X <= x | M) < ctprob. Then, define the Lmax to be minimum such that the P(lambda > Lmax) < ctprob.
      # p is assumed to be distributed as a beta RV with mean = pmean and variance = ((maxp-minp)/4)^2.
      # This determines beta parameters as follows:
#      Lest<-max(x, 0.5)/Eg
      #gmin<-gmax<-0
      Vg<-sum(sig2*a^2)#+(sum(a*(mu^2+sig2))-Eg^2)/Lest
      if (Vg > 0.000001){
        mmax<-fmmax.ab(x,pBa,pBb)
        if (ppois(mmax,100)<ctprob){
          lvec<-seq(0,100,by=0.01)
          Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))] # this is the biggest that lambda could credibly be and still have the biggest M conceivable
        } else {
          i<-1
          while (1){
            i<-i+1
            if (ppois(mmax,100*i)<ctprob){
              lvec<-seq((i-1)*100,i*100,by=0.01)
              Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))]
              break
            }
          }
        }
      } else {
        mmax<-fmmax(x,Eg)
        if (ppois(mmax,100)<ctprob){
          lvec<-seq(0,100,by=0.01)
          Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))] # this is the biggest that lambda could credibly be and still have the biggest M conceivable
        } else {
          i<-1
          while (1){
            i<-i+1
            if (ppois(mmax,100*i)<ctprob){
              lvec<-seq((i-1)*100,i*100,by=0.01)
              Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))]
              break
            }
          }
        }
      }
      if (Vg > 0.000001){
        cpLgX<-function(L,x,pBa,pBb,mmax) 1-pLgX.ab(L,x,pBa,pBb,mmax)
        meanL<-integrate(Vectorize(cpLgX,"L"),lower=0.0001,upper=Lmax,x=x,pBa,pBb,mmax=mmax)$val/sum(mydat$rel_wt)#styr
        CI<-c(optimize(f=optqL.ab,interval=c(0.00001,Lmax),x=x,pBa,pBb,mmax=mmax,p=.025)$minimum, optimize(f=optqL.ab,interval=c(0.00001,Lmax),x=x,pBa,pBb,mmax=mmax,p=.975)$minimum)/sum(mydat$rel_wt)#styr
      } else {
        cpLgX<-function(L,x,g,mmax) 1-pLgX(L,x,g,mmax)
        meanL<-integrate(Vectorize(cpLgX,"L"),lower=0.0001,upper=Lmax,x=x,g=Eg,mmax=mmax)$val/sum(mydat$rel_wt)#styr
        CI<-c(optimize(f=optqL,interval=c(0.00001,Lmax),x=x,g=Eg,mmax=mmax,p=.025)$minimum, optimize(f=optqL,interval=c(0.00001,Lmax),x=x,g=Eg,mmax=mmax,p=.975)$minimum)/sum(mydat$rel_wt)#styr
      }
      text(max(m),top,paste("Estimated baseline annual fatality rate: \u03bb\u0303 = ",signif(meanL,3), "\n 95% CI = [ ",signif(CI[1],3),", ",signif(CI[2],4),"]",sep=""),adj=1,cex=.85*.Rvar$charSizeAdjust)
      text(max(m),.9*top,paste("Overall detection probability: g = ",signif(Eg,4), "\n 95% CI = [ ",signif(gmin,4),", ",signif(gmax,4),"]",sep=""),adj=1,cex=.85*.Rvar$charSizeAdjust)
      if (sides=='1'){
        text(max(m),.8*top,paste("M* = ",m[cut.off], ", i.e., P(M \u2264 ",  m[cut.off],") \u2265 ", (1-aM)*100, "%",sep=''),adj=1,cex=.85*.Rvar$charSizeAdjust)
      } else {
        text(max(m),.8*top,paste0(100*(1-aM),"% CI for total M = [",m[min(reds)],", ",  m[max(reds)],"]  "),adj=1,cex=.85*.Rvar$charSizeAdjust)
      }
      mtext(side=3,paste('Credibility level (1 - \u03b1) = ',1-aM),adj=0,family='serif', cex=.Rvar$charSizeAdjust)
      title(paste("Posterior Distribution of Total Fatalities over ",yrs," years",sep=""))
      M<-1:length(pMem)-1
      CpMem<-1-cumsum(pMem)
### write a table of results:
      while(1){ if (sink.number()==0) break else sink() }
      sink(paste0(.Rvar$datadir, "/output"))
      cat(paste("Summary statistics for total mortality through ", length(mydat$years), " years\n",sep=''))
      cat("------------------------------------\n")
      cat("Results\n")
      if (sides == 2){
        cat(paste0(100*(1-aM),"% CI for M = [",m[min(reds)],", ",  m[max(reds)],"]\n\n"))
      } else {
        cat(paste0("M* = ",m[cut.off], " for 1 - \u03b1 = ",1-aM,", i.e., P(M <= ",  m[cut.off],") >= ", (1-aM)*100, "%", '\n\n'))
      }
      cat(paste0("Estimated overall detection probability: g = ", signif(Eg, 3),", 95% CI = [",signif(gmin,3),", ",signif(gmax,3),"]\n",
        "   Ba = ", signif(pBa,5), ", Bb = ", signif(pBb,5), "\n\n"))
#      CI<-c(L[min(which(pL>(1-aCI)/2))]/sum(mydat$rel_wt),L[min(which(pL>aCI+(1-aCI)/2))]/sum(mydat$rel_wt))
#      CI<-c(L[min(which(pL>.025))]/sum(mydat$rel_wt),L[min(which(pL>.975))]/sum(mydat$rel_wt))
#      par(family='serif')
      cat(paste0("Estimated baseline fatality rate: lambda = ",signif(meanL,4),", 95% CI = [",signif(CI[1],3),", ",signif(CI[2],3),"]\n"))
      cat(paste0("\n\nTest of assumed relative weights (rho) and potential bias"))
      rhosumry<-rhotest(mydat)
      cat("              Fitted rho\n")
      cat("Assumed rho     95% CI \n")
      for (yi in 1:length(rhosumry$rho.ass)) {
        cat(sprintf(
          " %6.3g     [%5.3f, %5.3f]\n", rhosumry$rho.ass[yi], rhosumry$rhoqtls[1, yi], rhosumry$rhoqtls[5, yi]))
      }
#      cat("\nAssumed rho: ")
#      for (yi in 1:length(rhosumry$rho.ass)) cat(sprintf(" %5.3f", rhosumry$rho.ass[yi]))
#      cat("\nFitted rho:  ")
#      for (yi in 1:length(rhosumry$rho.obs)) cat(sprintf(" %5.3f", rhosumry$rho.obs[yi]*sum(rhosumry$rho.ass)))
      cat(paste0("\np = ", round(rhosumry$pval, 5), ' for likelihood ratio test of H0: assumed rho = true rho\n'))
      cat(paste0("Quick test of relative bias: ", round(rhosumry$quicktest, 3),'\n\n'))
      .Rvar$MposteriorPDF <- cbind(1:length(pMgX)-1, pMgX)
      cat('Posterior distribution of M\n')
      cat('m     p(M = m) p(M > m)\n')
      Mcdf<-cumsum(.Rvar$MposteriorPDF[,2])
      for (i in 1:length(.Rvar$MposteriorPDF[,2])) {
        cat(sprintf("%-5.0f  %6.4f   %6.4f\n",i-1, round(.Rvar$MposteriorPDF[i,2],4), round(1-Mcdf[i],4)))
      }
      cat(paste0(unlist(rep("=",80)),collapse=''))
      cat("\nInput\n")
      cat("Year (or period) rel_wt  X    Ba     Bb   ghat    95% CI\n")
      for (i in 1:nyr) cat(sprintf("%-12s     %5.3f %3.0f %6.4g %6.4g %5.3f [%5.3f, %5.3f]\n",mydat$years[i], rel_wt[i], X[i], pBa0[i], pBb0[i], pBa0[i]/(pBa0[i]+pBb0[i]), qbeta(.025, pBa0[i], pBb0[i]), qbeta(.975, pBa0[i], pBb0[i])))
      sink()
      file.show(paste0(.Rvar$datadir, "/output"),delete.file=T,title=paste0("Mortality over ", yrs, " years"))
    } else if (mydat$option=="L" & mydat$Ltype=="stT"){
    # short-term trigger
      styr<-mydat$styr
      X<-mydat$X; nyr<-length(X)
      yind<-(nyr-styr+1):nyr
      yrs<-mydat$years[yind]
      X<-X[yind]
      x<-sum(X)
      pBa0<-mydat$Ba[yind]; pBb0<-mydat$Bb[yind]
      rel_wt<-mydat$rel_wt[yind]
      a<-mydat$rel_wt[yind]/sum(mydat$rel_wt[yind])
      aL<-mydat$aL
      # p is assumed to be distributed as a beta RV with mean = pmean and variance = ((maxp-minp)/4)^2.
      # This determines beta parameters as follows:
      tau<-mydat$tau
      mu<-pBa0/(pBa0+pBb0)
      sig2<-pBa0*pBb0/((pBa0+pBb0)^2*(pBa0+pBb0+1))
      Eg<-sum(mu*a) # average, expected G is the weighted average of the observed g's
#      Lest<-max(sum(X), 0.5)/Eg
      Vg<-sum(sig2*a^2)#+(sum(a*(mu^2+sig2))-Eg^2)/Lest
      pBa<-Eg^2/Vg*(1-Eg)-Eg; pBb<-pBa*(1/Eg-1) # these are the two shape parameters for the beta distribution underlying the beta-binomial for X | M
      ## governing parameters for the calculation
      ctprob<-0.0001 # A Jeffreys prior is used for the distribution of lambda. What is a reasonable maximum (and minimum) for lambda?
      # Pick a mmax = minimum M value such that P(X <= x | M) < ctprob. Then, define the Lmax to be minimum such that the P(lambda > Lmax) < ctprob.
      # p is assumed to be distributed as a beta RV with mean = pmean and variance = ((maxp-minp)/4)^2.

      # This determines beta parameters as follows:
#      Lest<-max(x, 0.5)/Eg
      gmin<-gmax<-0
      Vg<-sum(sig2*a^2)#+(sum(a*(mu^2+sig2))-Eg^2)/Lest
      if (Vg > 0.000001){
        mmax<-fmmax.ab(x,pBa,pBb)
        if (ppois(mmax,100)<ctprob){
          lvec<-seq(0,100,by=0.01)
          Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))] # this is the biggest that lambda could credibly be and still have the biggest M conceivable
        } else {
          i<-1
          while (1){
            i<-i+1
            if (ppois(mmax,100*i)<ctprob){
              lvec<-seq((i-1)*100,i*100,by=0.01)
              Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))]
              break
            }
          }
        }
      #  Lalpha<-optimize(f=optqL.ab,interval=c(0.001,Lmax),x=x,pBa=pBa,pBb=pBb,mmax=mmax,p=alpha)$minimum
        Lalpha<-optimize(optqL.ab,interval=c(0.00001,Lmax),x=x,pBa=pBa,pBb=pBb,mmax=mmax,p=aL)$min/styr
        gmin<-qbeta(0.025,pBa,pBb); gmax<-qbeta(0.975,pBa,pBb)
      } else {
        mmax<-fmmax(x,Eg)
        if (ppois(mmax,100)<ctprob){
          lvec<-seq(0,100,by=0.01)
          Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))] # this is the biggest that lambda could credibly be and still have the biggest M conceivable
        } else {
          i<-1
          while (1){
            i<-i+1
            if (ppois(mmax,100*i)<ctprob){
              lvec<-seq((i-1)*100,i*100,by=0.01)
              Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))]
              break
            }
          }
        }
        Lalpha<-optimize(f=optqL,interval=c(0.00001,Lmax),x=x,g=Eg,mmax=mmax,p=aL)$minimum/styr
        gmin<-Eg; gmax<-Eg
      }
      if (Vg > 0.000001){
        cpLgX<-function(L,x,pBa,pBb,mmax) 1-pLgX.ab(L,x,pBa,pBb,mmax)
        meanL<-integrate(Vectorize(cpLgX,"L"),lower=0.0001,upper=Lmax,x=x,pBa,pBb,mmax=mmax)$val/styr
        CI<-c(optimize(f=optqL.ab,interval=c(0.00001,Lmax),x=x,pBa,pBb,mmax=mmax,p=.025)$minimum, optimize(f=optqL.ab,interval=c(0.00001,Lmax),x=x,pBa,pBb,mmax=mmax,p=.975)$minimum)/styr
      } else {
        cpLgX<-function(L,x,g,mmax) 1-pLgX(L,x,g,mmax)
        meanL<-integrate(Vectorize(cpLgX,"L"),lower=0.0001,upper=Lmax,x=x,g=Eg,mmax=mmax)$val/styr
        CI<-c(optimize(f=optqL,interval=c(0.00001,Lmax),x=x,g=Eg,mmax=mmax,p=.025)$minimum, optimize(f=optqL,interval=c(0.00001,Lmax),x=x,g=Eg,mmax=mmax,p=.975)$minimum)/styr
      }

      ### write results:
      while(1){ if (sink.number()==0) break else sink() }
      sink(paste0(.Rvar$datadir, "/output"))
      cat(paste0("Short-term trigger: Test of average fatality rate (lambda) over ", styr, " years\n"))
      cat(paste0("Years: ",yrs[1]," - ",yrs[styr],"\n"))
      cat(paste0(unlist(rep("=",80)),collapse=''))
      cat("\nResults\n")
      cat(paste0("Estimated overall detection probability: g = ",
       signif(Eg, 3),", 95% CI = [",signif(gmin,3),", ",signif(gmax,3),"]\n   Ba = ", signif(pBa, 5), ", Bb = ", signif(pBb, 5), "\n\n"))
#      CI<-c(L[min(which(pL>(1-aCI)/2))]/sum(mydat$rel_wt),L[min(which(pL>aCI+(1-aCI)/2))]/sum(mydat$rel_wt))
#      CI<-c(L[min(which(pL>.025))]/sum(mydat$rel_wt),L[min(which(pL>.975))]/sum(mydat$rel_wt))
#      par(family='serif')
      cat(paste0("Estimated annual fatality rate over the past ", styr," years: lambda = ",signif(meanL,4),", 95% CI = [",signif(CI[1],3),", ",signif(CI[2],3),"]\n"))
      cat(paste0("   P(lambda > ", tau, ") = ", 1-round(pLgX.ab(tau*styr, x, pBa, pBb, mmax),4),sep=''))
      cat(paste0("\n   ",ifelse(Lalpha > tau,paste0("Exceedance: lambda > ",tau, " with ",(1-aL)*100,"% credibility"),paste0("Compliance: Cannot infer lambda > ", tau, " with ",(1-aL)*100, "% credibility"))))
      cat("\n")
      cat(paste0(unlist(rep("_",80)),collapse=''))
      cat("\n")
      cat("Input\n")
      cat("Threshold for short-term rate (tau) = ", tau, " per year\n\n")
      cat("Period       rel_wt  X    Ba     Bb   ghat    95% CI\n")
      for (i in 1:styr) cat(sprintf("%-12s %5.3f %3.0f %6.4g %6.4g %5.3f [%5.3f, %5.3f]\n",yrs[i], rel_wt[i], X[i], pBa0[i], pBb0[i], pBa0[i]/(pBa0[i]+pBb0[i]), qbeta(.025, pBa0[i], pBb0[i]), qbeta(.975, pBa0[i], pBb0[i])))
      sink()
      file.show(paste0(.Rvar$datadir, "/output"),delete.file=T,title="Short-term Trigger")
    } else if (mydat$option=="L" & mydat$Ltype=="rT"){ # reversion trigger
    ### reversion trigger
      X<-mydat$X; nyr<-length(X)
      yrs<-mydat$years
      X<-X
      x<-sum(X)
      pBa0<-mydat$Ba; pBb0<-mydat$Bb
      rel_wt<-mydat$rel_wt
      a<-mydat$rel_wt/sum(mydat$rel_wt)
      aR<-mydat$aR
      rho<-mydat$rho
      tau<-mydat$tau
      mu<-pBa0/(pBa0+pBb0)
      sig2<-pBa0*pBb0/((pBa0+pBb0)^2*(pBa0+pBb0+1))
      Eg<-sum(mu*a) # average, expected G is the weighted average of the observed g's
#      Lest<-max(sum(X), 0.5)/Eg
      Vg<-sum(sig2*a^2)#+(sum(a*(mu^2+sig2))-Eg^2)/Lest
      pBa<-Eg^2/Vg*(1-Eg)-Eg; pBb<-pBa*(1/Eg-1) # these are the two shape parameters for the beta distribution underlying the beta-binomial for X | M
      ## governing parameters for the calculation
      ctprob<-0.0001 # A Jeffreys prior is used for the distribution of lambda. What is a reasonable maximum (and minimum) for lambda?
#      Lest<-max(x, 0.5)/Eg
      gmin<-gmax<-0
      Vg<-sum(sig2*a^2)#+(sum(a*(mu^2+sig2))-Eg^2)/Lest
      if (Vg > 0.000001){
        mmax<-fmmax.ab(x,pBa,pBb)
        if (ppois(mmax,100)<ctprob){
          lvec<-seq(0,100,by=0.01)
          Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))] # this is the biggest that lambda could credibly be and still have the biggest M conceivable
        } else {
          i<-1
          while (1){
            i<-i+1
            if (ppois(mmax,100*i)<ctprob){
              lvec<-seq((i-1)*100,i*100,by=0.01)
              Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))]
              break
            }
          }
        }
        Lalpha<-optimize(optqL.ab,interval=c(0.00001,Lmax),x=x,pBa=pBa,pBb=pBb,mmax=mmax,p=1-aR)$min/nyr
        gmin<-qbeta(0.025,pBa,pBb); gmax<-qbeta(0.975,pBa,pBb)
      } else {
        mmax<-fmmax(x,Eg)
        if (ppois(mmax,100)<ctprob){
          lvec<-seq(0,100,by=0.01)
          Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))] # this is the biggest that lambda could credibly be and still have the biggest M conceivable
        } else {
          i<-1
          while (1){
            i<-i+1
            if (ppois(mmax,100*i)<ctprob){
              lvec<-seq((i-1)*100,i*100,by=0.01)
              Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))]
              break
            }
          }
        }
        Lalpha<-optimize(f=optqL,interval=c(0.00001,Lmax),x=x,g=Eg,mmax=mmax,p=1-aR)$minimum/nyr
        gmin<-Eg; gmax<-Eg
      }
      if (Vg > 0.000001){
        cpLgX<-function(L,x,pBa,pBb,mmax) 1-pLgX.ab(L,x,pBa,pBb,mmax)
        meanL<-integrate(Vectorize(cpLgX,"L"),lower=0.0001,upper=Lmax,x=x,pBa,pBb,mmax=mmax)$val/nyr
        CI<-c(optimize(f=optqL.ab,interval=c(0.00001,Lmax),x=x,pBa,pBb,mmax=mmax,p=.025)$minimum, optimize(f=optqL.ab,interval=c(0.00001,Lmax),x=x,pBa,pBb,mmax=mmax,p=.975)$minimum)/nyr
        pLgtrhotau<-pLgX.ab(rho*tau*nyr, x, pBa, pBb, mmax)
      } else {
        cpLgX<-function(L,x,g,mmax) 1-pLgX(L,x,g,mmax)
        meanL<-integrate(Vectorize(cpLgX,"L"),lower=0.0001,upper=Lmax,x=x,g=Eg,mmax=mmax)$val/nyr
        CI<-c(optimize(f=optqL,interval=c(0.00001,Lmax),x=x,g=Eg,mmax=mmax,p=.025)$minimum, optimize(f=optqL,interval=c(0.00001,Lmax),x=x,g=Eg,mmax=mmax,p=.975)$minimum)/nyr
      }
      while(1){ if (sink.number()==0) break else sink() }
      sink(paste0(.Rvar$datadir, "/output"))
      cat(paste0("Test whether mortality rate (stochastic) over ", nyr, " years is lower than rho * tau = ",rho," * ",tau, " = ", rho*tau,"\n"))
      cat(paste0("Years: ",yrs[1]," - ",yrs[nyr],"\n"))
      cat(paste0(unlist(rep("=",80)),collapse=''))
      cat("\nResults\n")
      cat(paste0("Total number of carcasses recovered: ", sum(X),"\n"))
      cat(paste0("Estimated overall detection probability: g = ",
        signif(Eg, 3),", 95% CI = [",signif(gmin,3),", ",signif(gmax,3),"]\n   Ba = ", signif(pBa, 5), ", Bb = ", signif(pBb, 5), "\n\n"))
      cat(paste0("Estimated average annual fatality rate: lambda = ",signif(meanL,4),", 95% CI = [",signif(CI[1],3),", ",signif(CI[2],3),"]\n"))
      cat(paste0("   p(lambda <= ", tau*rho, ") = ", round(pLgX.ab(tau*rho*nyr, sum(X), pBa, pBb, mmax),4)))
      cat(paste0("\n   ",ifelse(pLgX.ab(tau*rho*nyr, sum(X), pBa, pBb, mmax) >= 1-aR,paste0("Reversion: lambda < rho * tau = ",tau*rho, " with ",(1-aR)*100,"% credibility"),paste0("No reversion: Cannot infer lambda < rho * tau = ", tau*rho, " with ",(1-aR)*100, "% credibility"))))
      cat("\n")
      cat(paste0(unlist(rep("_",80)),collapse=''))
      cat("\n")
      cat("Input\n")
      cat(paste0("Threshold for short-term rate (tau)   = ", tau, " per year\n\n"))
      cat(paste0("Effect of AMA to reverse (rho) = ", rho*100, "% expected reduction in fatality rate\n\n"))
      cat("Year (or period) rel_wt  X    Ba     Bb   ghat    95% CI\n")
      for (i in 1:nyr) cat(sprintf("%-12s     %5.3f %3.0f %6.4g %6.4g %5.3f [%5.3f, %5.3f]\n",yrs[i], rel_wt[i], X[i], pBa0[i], pBb0[i], pBa0[i]/(pBa0[i]+pBb0[i]), qbeta(.025, pBa0[i], pBb0[i]), qbeta(.975, pBa0[i], pBb0[i])))
      sink()
      file.show(paste0(.Rvar$datadir, "/output"),delete.file=T,title="Reversion Trigger")
    } else if (mydat$option=="L" & mydat$Ltype=="ciT"){ # confidence interval for short-term rate
      X<-mydat$X; nyr<-length(X)
      yrs<-mydat$years
      X<-X
      x<-sum(X)
      pBa0<-mydat$Ba; pBb0<-mydat$Bb
      rel_wt<-mydat$rel_wt
      a<-mydat$rel_wt/sum(mydat$rel_wt)
      aL<-mydat$aL
      # p is assumed to be distributed as a beta RV with mean = pmean and variance = ((maxp-minp)/4)^2.
      # This determines beta parameters as follows:
      tau<-mydat$tau
      mu<-pBa0/(pBa0+pBb0)
      sig2<-pBa0*pBb0/((pBa0+pBb0)^2*(pBa0+pBb0+1))
      Eg<-sum(mu*a) # average, expected G is the weighted average of the observed g's
#      Lest<-max(sum(X), 0.5)/Eg
      Vg<-sum(sig2*a^2)#+(sum(a*(mu^2+sig2))-Eg^2)/Lest
      pBa<-Eg^2/Vg*(1-Eg)-Eg; pBb<-pBa*(1/Eg-1) # these are the two shape parameters for the beta distribution underlying the beta-binomial for X | M
      ## governing parameters for the calculation
      ctprob<-0.0001 # A Jeffreys prior is used for the distribution of lambda. What is a reasonable maximum (and minimum) for lambda?
      # Pick a mmax = minimum M value such that P(X <= x | M) < ctprob. Then, define the Lmax to be minimum such that the P(lambda > Lmax) < ctprob.
      # p is assumed to be distributed as a beta RV with mean = pmean and variance = ((maxp-minp)/4)^2.
      # This determines beta parameters as follows:
#      Lest<-max(x, 0.5)/Eg
      gmin<-gmax<-0
      Vg<-sum(sig2*a^2)#+(sum(a*(mu^2+sig2))-Eg^2)/Lest
      if (Vg > 0.000001){
        mmax<-fmmax.ab(x,pBa,pBb)
        if (ppois(mmax,100)<ctprob){
          lvec<-seq(0,100,by=0.01)
          Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))] # this is the biggest that lambda could credibly be and still have the biggest M conceivable
        } else {
          i<-1
          while (1){
            i<-i+1
            if (ppois(mmax,100*i)<ctprob){
              lvec<-seq((i-1)*100,i*100,by=0.01)
              Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))]
              break
            }
          }
        }
      #  Lalpha<-optimize(f=optqL.ab,interval=c(0.001,Lmax),x=x,pBa=pBa,pBb=pBb,mmax=mmax,p=alpha)$minimum
        Lalpha<-optimize(optqL.ab,interval=c(0.00001,Lmax),x=x,pBa=pBa,pBb=pBb,mmax=mmax,p=aL)$min/nyr
        gmin<-qbeta(0.025,pBa,pBb); gmax<-qbeta(0.975,pBa,pBb)
      } else {
        mmax<-fmmax(x,Eg)
        if (ppois(mmax,100)<ctprob){
          lvec<-seq(0,100,by=0.01)
          Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))] # this is the biggest that lambda could credibly be and still have the biggest M conceivable
        } else {
          i<-1
          while (1){
            i<-i+1
            if (ppois(mmax,100*i)<ctprob){
              lvec<-seq((i-1)*100,i*100,by=0.01)
              Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))]
              break
            }
          }
        }
        Lalpha<-optimize(f=optqL,interval=c(0.00001,Lmax),x=x,g=Eg,mmax=mmax,p=aL)$minimum/nyr
        gmin<-Eg; gmax<-Eg
      }
      if (Vg > 0.000001){
        binog<-F
        cpLgX<-function(L,x,pBa,pBb,mmax) 1-pLgX.ab(L,x,pBa,pBb,mmax)
        meanL<-integrate(Vectorize(cpLgX,"L"),lower=0.0001,upper=Lmax,x=x,pBa,pBb,mmax=mmax)$val/nyr
        CI<-c(optimize(f=optqL.ab,interval=c(0.00001,Lmax),x=x,pBa,pBb,mmax=mmax,p=(1-aCI)/2)$minimum, optimize(f=optqL.ab,interval=c(0.00001,Lmax),x=x,pBa,pBb,mmax=mmax,p=1-(1-aCI)/2)$minimum)/nyr
      } else {
        binog<-T
        cpLgX<-function(L,x,g,mmax) 1-pLgX(L,x,g,mmax)
        meanL<-integrate(Vectorize(cpLgX,"L"),lower=0.0001,upper=Lmax,x=x,g=Eg,mmax=mmax)$val/nyr
        CI<-c(optimize(f=optqL,interval=c(0.00001,Lmax),x=x,g=Eg,mmax=mmax,p=(1-aCI)/2)$minimum, optimize(f=optqL,interval=c(0.00001,Lmax),x=x,g=Eg,mmax=mmax,p=1-(1-aCI)/2)$minimum)/nyr
      }
      ### write results:
      while(1){ if (sink.number()==0) break else sink() }
      rhosumry<-rhotest(mydat)
      sink(paste0(.Rvar$datadir, "/output"))
      cat(paste0("Estimation of mortality rate (stochastic) over ", nyr, " years\n"))
      cat(paste0("Years: ",yrs[1]," - ",yrs[nyr],"\n"))
      cat(paste0(unlist(rep("=",80)),collapse=''))
      cat("\nResults\n")
      cat(paste0("Total number of carcasses recovered: ", sum(X),"\n"))
      if (binog){
        cat(paste0("Overall detection probability, g = ", signif(Eg, 3), "\n\n"))
      } else {
      cat(paste0("Estimated overall detection probability, g = ",
        signif(Eg, 3),", 95% CI = [",signif(gmin,3),", ",signif(gmax,3),"]\n   Ba = ", signif(pBa, 5), ", Bb = ",signif(pBb, 5), "\n\n"))
      }
      cat(paste0(
        "Estimated annual fatality rate: \n lambda = ", signif(meanL,3),
        " with ", 100*aCI, "% CI = [",signif(CI[1],3),", ",signif(CI[2],3),"]",sep=''))
      cat("\n")
      cat(paste0(unlist(rep("_",80)),collapse=''))
      cat("\n")
      cat("Input\n")
      cat("Threshold for short-term rate (tau) = ", tau, " per year\n\n")
      cat("Year (or period) rel_wt  X    Ba     Bb   ghat    95% CI\n")
      for (i in 1:nyr) cat(sprintf("%-12s     %5.3f %3.0f %6.4g %6.4g %5.3f [%5.3f, %5.3f]\n",yrs[i], rel_wt[i], X[i], pBa0[i], pBb0[i], pBa0[i]/(pBa0[i]+pBb0[i]), qbeta(.025, pBa0[i], pBb0[i]), qbeta(.975, pBa0[i], pBb0[i])))
      sink()
      file.show(paste0(.Rvar$datadir, "/output"),delete.file=T,title="Estimation of Mortality Rate (stochastic)")
    } else {
      return(F)
    }
    if (mydat$option == "L"){
      if (Vg > 0.000001){
        .Rvar$LposteriorCDF <- Vectorize(function(L) pLgX.ab(L, x=x, pBa=pBa, pBb = pBb, mmax = mmax), vectorize.args='L')
      } else {
        .Rvar$LposteriorCDF <- Vectorize(function(L) pLgX(L, x=x, g = g, mmax = mmax), vectorize.args='L')
      }
    }
}
pchk<-function(a){ #probability? i.e., numerical value in [0, 1]
  if (length(a) != 1 || is.na(a) || nchar(a) == 0 || a<=0 || a>=1) return(F)
  return(T)
}
pintchk<-function(n){ # positive integer (non-zero)?
  if (length(n) != 1 || is.na(n) || nchar(n) == 0 || n<=0 || n != round(n)) return(F)
  return(T)
}
gt0chk<-function(x){
  if (length(x) != 1 || is.na(x) || nchar(x) == 0 || x<=0) return(F)
  return(T)
}
