#compcolor<-colors()[626]
checkok1<-function(dat){
  pd<-dat$persistence_distn
  ck1<-(Xok & aok & SEnok & kok & Iok & nsearchok & ((pd=="Exponential" & lamok) | (pd!="Exponential" & pdaok & pdbok)))
  ck2<- (syAok & SExok & bminok & bmaxok & prok)
  ck1 & ck2 & startok
}
checkgok1<-function(dat){ # dat is the database that would be used to calculate g if all data are valid
  pd<-dat$persistence_distn
  ck1<-(aok & SEnok & kok & Iok & nsearchok & ((pd=="Exponential" & lamok) | (pd!="Exponential" & pdaok & pdbok)))
  ck2<- (SExok & bminok & bmaxok)
  ck1 & ck2 & startok
}
MpriorOK<-function(prior){
  if (is.numeric(prior)){ # numeric => must be two-dimensional array with probabilities starting at m = 0
    if (length(dim(prior)) != 2){
      warning("error in prior")
      return(F)
    }
    if (abs(sum(diff(prior[,1])-rep(1,length(prior[,1])-1)))>0.00000001){ # differences in M values not uniformly = 1
      warning("error in prior")
      return(F)
    }
    if (prior[1,1] != 0){
      warning("error in prior: m[1] must be zero")
      return(F)
    }
    if (abs(sum(prior[,2])-1) > 0.00001){
      warning("prior: probabilities must sum to 1")
      return(F)
    }
  } else {
    return(F)
  }
  return(T)
}

feedR<-function(dat){ #dat is a list whose data will be made available without having to reference the long name of the list
  for (i in 1:length(dat)){
    assign(names(dat)[i],dat[[i]],pos=1)
  }
  if (dat$samtype=="Custom") {
    Isam<<-round(max(dat$days)/(length(dat$days)-1),1)
    span<<-max(dat$days)
    nsearch<<-length(days)-1
    nt<<-nsearch
    days<<-dat$days
  } else {
    Isam<<-dat$Isam
    nt<<-nsearch    # not including search at t = 0
    days<<-(0:nsearch)*Isam
  }
}
fmmax.ab<-function(x,pBa,pBb){ # find the maximum m to sum over in calculating posterior (for beta-binomial)
  if (VGAM::pbetabinom.ab(x, 1e5, pBa, pBb) > 0.0001) { # if too huge of M's are required, return a large one and give warning
    g <- pBa/(pBa + pBb)
    warning(paste0("P(X <= ", x," | g = ", g, ", m = 10000) = ", signif(pbinom(x,10000,g),6), ". Taking mmax = 10,000..."
  ))
    return(1e5)
  }
  mmax<-x
  while (1){
    m<-mmax:(mmax+100)
    mmax<-mmax+100
    if (VGAM::pbetabinom.ab(x,size=mmax,shape1=pBa,shape2=pBb)<0.0001){
      mmax<-m[min(which(VGAM::pbetabinom.ab(x,size=m,shape1=pBa,shape2=pBb)<0.0001))]
      break
    }
  }
  mmax
}

fmmax<-function(x,g){ # find the maximum m to sum over in calculating posterior
  if (pbinom(x,1e5,g) > 0.0001) { # if too huge of M's are required, return a large one and give warning
    warning(paste0("P(X <= ", x," | g = ", g, ", m = 100000) = ", signif(pbinom(x,100000,g),6), ". Using mmax = 100,000..."    ))
    return(1e5)
  }
  mmax<-x
  while (1){
    m<-mmax:(mmax+100)
    mmax<-mmax+100
    if (pbinom(x,mmax,g)<0.0001){ # pick mmax large enough so that it is very unlikely to have a higher m
      mmax<-m[min(which(pbinom(x,m,g)<0.0001))]
      break
    }
  }
  mmax
}
postM.ab<-function(x, Ba, Bb, prior="IbinRef", mmax=NULL){
# choices for prior are:
# IbinRef = integrated reference prior for binomial
# binRef = reference prior for binomial (gives Inf for P(M = 0 | X = 0)
# IbetabinRef = integrated reference prior for betabinomial
# betabinRef = reference prior for betabinomial (gives Inf for P(M = 0 | X = 0)
# uniform
# custom prior = numeric array with columns m and p(M = m)
  # error-checking
  suppressWarnings({
    if (length(x) * length(Ba) * length(Bb) != 1){
      print("error in data")
      return(F)
    }
    if (is.na(as.numeric(x)) * is.na(as.numeric(Ba)) * is.na(as.numeric(Bb)) != 0){
      print("error in data: all values must be numeric")
      return(F)
    }
    if (x < 0 || abs(round(x)-x) > 0.000001){
      print("error in data: x must be a non-negative integer")
      return(F)
    }
    if (Ba < 0.000001){
      print("error in data: Ba must be positive")
      return(F)
    }
    if (Bb < 0.000001){
      print("error in data: Bb must be positive")
      return(F)
    }
  })
  x<-round(x)
  # check whether the mmax provided is valid
  if (!missing(mmax)){ # then mmax must be a non-negative integer
    if (!is.numeric(mmax)){
      warning('mmax must be a non-negative integer')
      return(F)
    }
    if (abs(round(mmax)-mmax)>0.000001){
      warning('mmax must be an integer')
      return(F)
    }
    if (mmax< -0.000001){
      warning('mmax must be non-negative')
      return(F)
    }
  }
  # check prior
  # if a custom prior is entered, check the format; if string is entered, calculated prior
  if (!is.numeric(prior) && !is.character(prior)){ # format is neither custom (numeric) nor named (string)
    warning("error in prior")
    return(F)
  } else {
    if (is.numeric(prior)){ # numeric => custom prior; must be two-dimensional array with probabilities starting at m = 0
      if (MpriorOK(prior)){
        if (missing(mmax)){
          mmax<-max(prior[,1])
        } else {
          if (round(mmax) != round(max(prior[,1]))) {
            warning("mmax != max of custom prior. Aborting calculation.")
            return(F)
          } else {
            mmax<-round(max(prior[,1]))
            M<-x:mmax
            prior_M<-prior[,2]
          }
        }
      } else { # custom prior is not correctly formatted (numeric but not two columns etc.)
        return(F)
      }
    } else { # named priors (following different philosophies of uninformed)
      if (prior %in% c('IbinRef','binRef', 'IbetabinRef', 'betabinRef', 'uniform')) {
        mmax<-round(ifelse(missing(mmax), fmmax.ab(x, Ba, Bb), mmax))
        M<-x:mmax
        prior_M<-switch(prior,
          IbinRef = diff(sqrt(0:(mmax+1)))/sum(diff(sqrt(0:(1+mmax)))),
          binRef = 1/sqrt(0:mmax),
          IbetabinRef = diff(log(sqrt(0:(mmax+1)+Ba+Bb)+sqrt(0:(mmax+1)))),
          betabinRef = 1/sqrt(0:mmax*(0:mmax+Ba+Bb)),
          uniform = rep(1, mmax + 1)
        )
      } else {
          warning("error in prior")
          return(F)
      }
    }
  }
#  prior_M <- prior_M/sum(prior_M)
  if (prior_M[1] == Inf & x == 0) return (1)
  pXgM<-VGAM::dbetabinom.ab(x, M, Ba, Bb)
  pM<-prior_M[x:mmax+1]
  pMgX<-pXgM*pM; pMgX<-pMgX/sum(pMgX) # posterior distribution for M (ignoring M < X, which has probability = zero)
  pMgX<-c(rep(0,x),pMgX)
  pMgX
}

postM<-function(x, g, prior='IbinRef', mmax = NA){
 # choices for prior are:
 # IbinRef = integrated reference prior for binomial
 # binRef = reference prior for binomial (gives Inf for P(M = 0 | X = 0)
 # uniform
 # custom prior = numeric array with columns m and p(M = m)
  if (length(x) * length(g) != 1){
    warning("error in data: x and g must be scalars")
    return(F)
  }
  if (is.na(as.numeric(x)) * is.na(as.numeric(g)) != 0){
    warning("error in data: x and g must be numeric")
    return(F)
  }
  if (x < -0.000001 || abs(round(x)-x) > 0.000001){
    warning("error in data: x must be a non-negative integer")
    return(F)
  }
  if (g < 0.00001){
    warning("error in data: g must be strictly greater than 0")
    return(F)
  }
  if (g > 0.99999) return (x)
  if (!is.na(mmax) && (!is.numeric(mmax) || abs(round(mmax)-x) < 0.000001 || round(mmax) < 1)){
    print("problem is here")
    warning("error in mmax")
    return(F)
  }
  x<-round(x)
  if (!is.numeric(prior) && !is.character(prior)){ # custom prior
    warning("error in prior")
    return(F)
  } else {
    if (is.numeric(prior)){ # numeric => must be two-dimensional array with probabilities starting at m = 0
      if (MpriorOK(prior)){
        if (is.na(mmax)){
          mmax<-max(prior[,1])
        } else {
          if (round(mmax) != round(max(prior[,1]))) {
            warning("mmax != max of custom prior. Aborting calculation.")
            return(F)
          } else {
            mmax<-round(max(prior[,1]))
            M<-x:mmax
            prior_M<-prior[,2]/sum(prior[,2])
          }
        }
      } else { # custom prior is not correctly formatted (numeric but not two columns etc.)
        return(F)
      }
    } else { # named priors (following different philosophies of uninformed)
      if (prior %in% c('IbinRef','binRef', 'uniform')) {
        mmax<-round(ifelse(is.na(mmax), fmmax(x,g), mmax))
        M<-x:mmax
        prior_M<-switch(prior,
          IbinRef = diff(sqrt(0:(mmax+1)))/sum(diff(sqrt(0:(1+mmax)))),
          binRef = 1/sqrt(0:mmax),
          uniform = rep(1, mmax + 1)
        )
      } else {
          warning("error in prior")
          return(F)
      }
    }
  }
  if (prior_M[1] == Inf & x == 0) return (1)
  pXgM<-dbinom(x, M, g)
  pM<-prior_M[x:mmax+1]
  pMgX<-pXgM*pM; pMgX<-pMgX/sum(pMgX) # posterior distribution for M (ignoring M < X, which has probability = zero)
  c(rep(0,x),pMgX) # full posterior, counting from m = 0 to mmax
}

calcMstar<-function(pMgX, alpha){
  min(which(cumsum(pMgX)>=1-alpha))-1
}

MCI<-function(pMgX, crlev=0.95){
  cs<-cumsum(pMgX)
  aM<-1-crlev
  lwrbnd<-min(which(cs > aM/2))-1
  lwrArea<-ifelse(lwrbnd == 0, 0, cs[lwrbnd])
  uprbnd<-min(which(cs > 1-aM + lwrArea))-1
  c(lwrbnd,uprbnd)
}

postMstar.ab<-function(x, Ba, Bb, conf.level = 0.95){
  calcMstar(postM.ab(x=x, Ba=Ba, Bb=Bb), alpha=1-conf.level)
}
postMstar<-function(x, g, conf.level = 0.95){
  calcMstar(postM(x=x, g=g), alpha=1-conf.level)
}
postMCI.ab<-function(x, Ba, Bb, conf.level = 0.95){
  MCI(postM.ab(x=x, Ba=Ba, Bb=Bb), crlev=conf.level)
}
postMCI<-function(x, g, conf.level = 0.95){
  MCI(postM(x=x, g=g), crlev=conf.level)
}

postL.abCI<-function(x, Ba, Bb, conf.lev = 0.95){
  suppressWarnings({
    if (length(x) * length(Ba) * length(Bb) * length(conf.lev) != 1){
      print("error in data")
      return(F)
    }
    if (is.na(as.numeric(x)) * is.na(as.numeric(Ba)) * is.na(as.numeric(Bb)) * is.na(as.numeric(conf.lev)) != 0){
      print("error in data: all values must be numeric")
      return(F)
    }
    if (x < 0 || abs(round(x)-x) > 0.000001){
      print("error in data: x must be a non-negative integer")
      return(F)
    }
    if (Ba < 0.000001){
      print("error in data: Ba must be positive")
      return(F)
    }
    if (Bb < 0.000001){
      print("error in data: Bb must be positive")
      return(F)
    }
    if (conf.lev >= 0.999999 || conf.lev <= 0.000001){
      print("error in data: conf.lev must be strictly between 0 and 1")
      return(F)
    }
  })
  ctprob<-0.0001
  mmax<-fmmax.ab(x,Ba,Bb)
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
  cpLgX<-function(L,x,Ba,Bb,mmax) 1-pLgX.ab(L,x,Ba,Bb,mmax)
  meanL<-integrate(Vectorize(cpLgX,"L"),lower=0.0001,upper=Lmax,x=x,Ba,Bb,mmax=mmax)$val
  CI<-c(optimize(f=optqL.ab,interval=c(0.00001,Lmax),x=x,Ba,Bb,mmax=mmax,p=(1-conf.lev)/2)$minimum, optimize(f=optqL.ab,interval=c(0.00001,Lmax),x=x,Ba,Bb,mmax=mmax,p=1-(1-conf.lev)/2)$minimum)#styr
  list(meanL=meanL, CI=CI)
}
postL.CI<-function(x, g, conf.lev = 0.95){
  suppressWarnings({
    if (length(x) * length(g) * length(conf.lev) != 1){
      print("error in data")
      return(F)
    }
    if (is.na(as.numeric(x)) * is.na(as.numeric(g)) * is.na(as.numeric(conf.lev)) != 0){
      print("error in data: all values must be numeric")
      return(F)
    }
    if (x < 0 || abs(round(x)-x) > 0.000001){
      print("error in data: x must be a non-negative integer")
      return(F)
    }
    if (g < 0.00001 || g > 0.99999){
      print("error in data: g must be in interval (0, 1)")
      return(F)
    }
    if (conf.lev >= 0.999999 || conf.lev <= 0.000001){
      print("error in data: conf.lev must be strictly between 0 and 1")
      return(F)
    }
  })
  ctprob<-0.0001
  mmax<-fmmax(x,g)
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
  cpLgX<-function(L,x,g,mmax) 1-pLgX(L,x,g,mmax)
  meanL<-integrate(Vectorize(cpLgX,"L"),lower=0.0001,upper=Lmax,x=x,g,mmax=mmax)$val
  CI<-c(optimize(f=optqL,interval=c(0.00001,Lmax),x=x,g,mmax=mmax,p=(1-conf.lev)/2)$minimum, optimize(f=optqL,interval=c(0.00001,Lmax),x=x,g,mmax=mmax,p=1-(1-conf.lev)/2)$minimum)#styr
  list(meanL=meanL, CI=CI)
}
posteriorL.ab<-function(x, Ba, Bb, pL='jeffreys'){
  suppressWarnings({
    if (length(x) * length(Ba) * length(Bb)!= 1){
      print("error in data")
      return(F)
    }
    if (is.na(as.numeric(x)) * is.na(as.numeric(Ba)) * is.na(as.numeric(Bb)) != 0){
      print("error in data: all values must be numeric")
      return(F)
    }
    if (x < 0 || abs(round(x)-x) > 0.000001){
      print("error in data: x must be a non-negative integer")
      return(F)
    }
    if (Ba < 0.000001){
      print("error in data: Ba must be positive")
      return(F)
    }
    if (Bb < 0.000001){
      print("error in data: Bb must be positive")
      return(F)
    }
  })
  if (!is.character(pL) & !is.function(pL)){
    warning("error in prior")
    return(F)
  }
  ctprob<-0.0001
  mmax<-fmmax.ab(x, Ba, Bb)
  force(x); force(Ba); force(Bb)
  if (is.function(pL)){    # custom prior provided by user
    # minimal error-checking on the prior
    # improper (integrates to infinity) and non-normalized (does not integrate to 1) priors are allowed
    # prior is truncated above lambda where observing x or smaller would be extremely unlikely
    #   NOTE: this does not provide protection against pathological priors where all the weight
    #   is well above anything remotely compatible with the data
    ctprob<-1e-7
    if (ppois(mmax,100)<ctprob){
      lvec<-seq(0,100,by=0.01)
      Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))] # this is the biggest that lambda could credibly be and still have the biggest M conceivable
    } else {
      i<-0
      while (1){
        i<-i+1
        if (ppois(mmax,100*i)<ctprob){
          lvec<-seq((i-1)*100,i*100,by=0.01)
          Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))]
          break
        }
      }
    }
    Lmin<-1e-5
    deno<-integrate(pL, lower=Lmin, Lmax)$val
    pL0<-function(L) pL(L)/deno   # normalized version of pL, valid in interval [Lmin, Lmax]. Posterior will be set = 0 outside [Lmin, Lmax]
    ans0<-function(L){ # calculate the posterior
    # the numerator of the posterior
      numerator<-colSums(outer(x:mmax, L, F=dpois)*VGAM::dbetabinom.ab(x, x:mmax, Ba, Bb))*pL0(L)
    # the denominator of the posterior
      denominator<-0
      for (m in x:mmax){
        denominator<-denominator+integrate(function(L) colSums(outer(m, L, F=dpois)*VGAM::dbetabinom.ab(x, m, Ba, Bb))*pL0(L), lower=Lmin, upper=Lmax)$val
      }
      retval<-numeric(length(L))
      retval[L>=Lmin & L<=Lmax]<-numerator[L>=Lmin & L<=Lmax]/denominator
      retval
    }
    ans<-Vectorize(function(L) integrate(ans0, lower=Lmin, upper=max(L))$val, "L")
  } else { # jeffreys or uniform prior
    if (pL == 'uniform'){
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
      # for uniform distribution, denominator simplifies:
      # integral(sum(dbetabinom * dpois * pL)) = sum(dbetabinom * integral(dpois)) because pL cancels from numerator and denominator and dbetabinom is not a function of L
      # = sum(dbetabinom) because integral(dpois) = 1 regardless of what M is
      ans<-function(L) 1-colSums(outer(x:mmax, L, ppois) * VGAM::dbetabinom.ab(x, x:mmax, Ba, Bb))/sum(VGAM::dbetabinom.ab(x, x:mmax, Ba, Bb))
    } else if (pL == 'jeffreys'){
      lfac<-beta(x:mmax+0.5,0.5)/sqrt(pi)
      deno<-sum(VGAM::dbetabinom.ab(x,size=x:mmax,shape1=Ba,shape2=Bb)*lfac)
      first.part<-VGAM::dbetabinom.ab(x,x:mmax,shape1=Ba,shape2=Bb)*lfac
      ans<-function(L) {
        colSums(first.part*t(outer(L,x:mmax+0.5, FUN="pgamma")))/deno
      }
    } else {
      warning("error in prior")
      return (F)
    }
  }
  ans
}

calcLstar.ab<-function(x, Ba, Bb, alpha){
  ctprob<-0.0001
  mmax<-fmmax.ab(x, Ba, Bb)
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
  cpLgX<-function(L,x,Ba,Bb,mmax) 1-pLgX.ab(L,x,Ba,Bb,mmax)
  optimize(f=optqL.ab,interval=c(0.00001, Lmax), x=x ,Ba, Bb, mmax=mmax, p=1-alpha)$minimum
}
posteriorLpdf.ab<-function(x, Ba, Bb, pL='jeffreys'){
  suppressWarnings({
    if (length(x) * length(Ba) * length(Bb)!= 1){
      print("error in data")
      return(F)
    }
    if (is.na(as.numeric(x)) * is.na(as.numeric(Ba)) * is.na(as.numeric(Bb)) != 0){
      print("error in data: all values must be numeric")
      return(F)
    }
    if (x < 0 || abs(round(x)-x) > 0.000001){
      print("error in data: x must be a non-negative integer")
      return(F)
    }
    if (Ba < 0.000001){
      print("error in data: Ba must be positive")
      return(F)
    }
    if (Bb < 0.000001){
      print("error in data: Bb must be positive")
      return(F)
    }
  })
  if (!is.character(pL) & !is.function(pL)){
    warning("error in prior")
    return(F)
  }
  ctprob<-0.0001
  mmax<-fmmax.ab(x, Ba, Bb)
  force(x); force(Ba); force(Bb)
  if (is.function(pL)){    # custom prior provided by user
    # minimal error-checking on the prior
    # improper (integrates to infinity) and non-normalized (does not integrate to 1) priors are allowed
    # prior is truncated above lambda where observing x or smaller would be extremely unlikely
    #   NOTE: this does not provide protection against pathological priors where all the weight
    #   is well above anything remotely compatible with the data
    ctprob<-1e-7
    if (ppois(mmax,100)<ctprob){
      lvec<-seq(0,100,by=0.01)
      Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))] # this is the biggest that lambda could credibly be and still have the biggest M conceivable
    } else {
      i<-0
      while (1){
        i<-i+1
        if (ppois(mmax,100*i)<ctprob){
          lvec<-seq((i-1)*100,i*100,by=0.01)
          Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))]
          break
        }
      }
    }
    Lmin<-1e-5
    deno<-integrate(pL, lower=Lmin, Lmax)$val
    pL0<-function(L) pL(L)/deno   # normalized version of pL, valid in interval [Lmin, Lmax]. Posterior will be set = 0 outside [Lmin, Lmax]
    denominator<-0
    for (m in x:mmax){
      denominator<-denominator+integrate(function(L) colSums(outer(m, L, F=dpois)*VGAM::dbetabinom.ab(x, m, Ba, Bb))*pL0(L), lower=Lmin, upper=Lmax)$val
    }
    ans<-function(L){ # calculate the posterior
    # the numerator of the posterior
    # the denominator of the posterior
      numerator<-colSums(outer(x:mmax, L, F=dpois)*VGAM::dbetabinom.ab(x, x:mmax, Ba, Bb))*pL0(L)
      retval<-numeric(length(L))
      retval[L>=Lmin & L<=Lmax]<-numerator[L>=Lmin & L<=Lmax]/denominator
      retval
    }
  } else { # jeffreys or uniform prior
    if (pL == 'uniform'){
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
      # for uniform distribution, denominator simplifies:
      # integral(sum(dbetabinom * dpois * pL)) = sum(dbetabinom * integral(dpois)) because pL cancels from numerator and denominator and dbetabinom is not a function of L
      # = sum(dbetabinom) because integral(dpois) = 1 regardless of what M is
      ans<-function(L) colSums(outer(x:mmax, L, dpois) * VGAM::dbetabinom.ab(x, x:mmax, Ba, Bb))/sum(VGAM::dbetabinom.ab(x, x:mmax, Ba, Bb))
    } else if (pL == 'jeffreys'){
      lfac<-beta(x:mmax+0.5,0.5)/sqrt(pi)
      deno<-sum(VGAM::dbetabinom.ab(x,size=x:mmax,shape1=Ba,shape2=Bb)*lfac)
      first.part<-VGAM::dbetabinom.ab(x,x:mmax,shape1=Ba,shape2=Bb)*lfac
      ans<-function(L) {
        colSums(first.part*t(outer(L,x:mmax+0.5, FUN="dgamma")))/deno
      }
    } else {
      warning("error in prior")
      return (F)
    }
  }
  ans
}
pkfit<-function(pkdat, burn = 2000, n.iter = 2000){
# function to fit pk model to data from searcher efficiency field trials
  ## error-checking:
  if (!is.list(pkdat)){
    warning("pkdat must be list with elements n, M, and X")
    return
  }
  if (sum(c('n','M','X') %in% names(pkdat)) != 3){
    warning("pkdat must be list with elements n, M, and X")
    return
  }
  if (length(pkdat$M) != pkdat$n || length(pkdat$X) != pkdat$n){
    warning("pkdat$M and pkdat$X must be vectors of length pkdat$n")
    return
  }
  if (!is.numeric(pkdat$M) || !is.numeric(pkdat$X)){
    warning("pkdat$M and pkdat$X must be numeric")
    return
  }
  with (pkdat, {
    if (sum(M<1) > 0 || sum(round(M)-M)>0){
      warning("pkdat$M must be a vector of positive integers")
      return
    }
    if (sum(X<0) > 0 || sum(round(X)-X)>0){
      warning("pkdat$X must be a vector of non-negative integers")
      return
    }
  })
  pkjags <- rjags::jags.model(
    textConnection(eoa::pkmod),
    data = pkdat,
    inits = with(pkdat,list(p = X[1]/M[1],k = max(min((X[2]/M[2])/(X[1]/M[1]),.99),.01)))
  )
  update(pkjags, burn)
  as.matrix(rjags::coda.samples(pkjags, variable.names=c('p','k'), n.iter=n.iter)[[1]][,2:1])
}
################################################
cpfit<-function(cpdat, persistence_distn){
  if (!(persistence_distn %in% eoa::pdnames)){
    msg<-"error. persistence_distn must be one of the following: "
    for (nm in eoa::pdnames) msg<-paste0(msg, "\n   '",nm, "'")
    warning(paste0(msg, "\n\nAborting calculation."))
  }
  if (dim(cpdat)[2] < 2){
    warning(paste0("Error in data (",fileName,")\nRequired: two columns with data for CPmin and CPmax.\nCheck file."))
    return(F)
  } else {
    cpdat<-cpdat[,1:2]
    if (!is.data.frame(cpdat)) cpdat<-data.frame(cpdat)
    names(cpdat)<-c("CPmin", "CPmax")
  }
  if (!is.numeric(cpdat[,1]) || !is.numeric(cpdat[,2]) || sum(cpdat[,1] < 0) > 0 || sum(cpdat[,1] > cpdat[,2]) > 0) {
    tkmessageBox(message = "Error in cpdat. Cannot calculate.")
    return(F)
  }
  xind<-which(cpdat$CPmin == 0 & cpdat$CPmax == Inf)
  if (length(xind)>0){
    cpdat$CPmin<-cpdat$CPmin[-xind]
    cpdat$CPmax<-cpdat$CPmax[-xind]
  }
  cpdat$CPmin <- pmax(cpdat$CPmin, 0.001)
  event<-ifelse(cpdat$CPmin == cpdat$CPmax, 1, ifelse(cpdat$CPmax == Inf, 0, 3))
  left<-cpdat$CPmin
  right<-cpdat$CPmax
  right[event==0]<-cpdat$CPmin[event==0]

  if (sum(cpdat$CPmax==Inf)==length(cpdat$CPmax)){
    warning("No carcasses removed in persistence trials. Cannot fit model using 'survival' package. Use CP analysis submodule of Single Class module to fit model.")
    return(F)
  }
  surv<- survival::Surv(time=left, time2=right, event=event, type=c('interval'))
  # fit survival models to persistence data and plot
  return(switch(persistence_distn,
    "Exponential"   = survival::survreg(surv~1, dist="exponential"),
    "Weibull"       = survival::survreg(surv~1, dist="weibull"),
    "Lognormal"     = survival::survreg(surv~1, dist="lognormal"),
    "Log-Logistic"  = survival::survreg(surv~1, dist="loglogistic")
  ))
}
cpsim<-function(cpmod, nsim, option = "parms"){
  CPab<-array(dim=c(nsim,2))
  if (cpmod$dist == "exponential"){
    CPab[,2]<-exp(rnorm(nsim,mean=cpmod$coef, sd=sqrt(cpmod$var[1])))
    CPab[,1]<-1/CPab[,2]
  } else if (cpmod$dist == "weibull"){
    CPparms<-MASS::mvrnorm(nsim,c(cpmod$coef[1],cpmod$scale),cpmod$var)
    CPab[,1]<-1/CPparms[,2] #shape
    CPab[,2]<-exp(CPparms[,1]) # scale
    i0<-which(CPab[,1]<=0 | (3/CPab[,2])^CPab[,1]<1e-320)
    while(length(i0)>0){
      CPparms[i0,]<-MASS::mvrnorm(length(i0),c(cpmod$coef[1],cpmod$scale), cpmod$var)
      CPab[i0,1]<-1/CPparms[i0,2]
      CPab[i0,2]<-exp(CPparms[i0,1])
      i0<-which(CPab[,1]<=0 | (3/CPab[,2])^CPab[,1]<1e-320)
    }
  } else if (cpmod$dist == "loglogistic"){
    # NOTE: log-logistic a = shape = 1/mod.ll$scale, b = scale = exp(mod.ll$coef[1])
    CPparms<-MASS::mvrnorm(nsim,c(cpmod$coef[1],cpmod$scale),cpmod$var)
    CPab[,1]<-1/CPparms[,2] #shape
    CPab[,2]<-exp(CPparms[,1]) # scale
    i0<-which(CPab[,1]<=0)
    while(length(i0)>0){
      CPparms[i0,]<-MASS::mvrnorm(length(i0),c(cpmod$coef[1],cpmod$scale), cpmod$var)
      CPab[i0,1]<-1/CPparms[i0,2]
      CPab[i0,2]<-exp(CPparms[i0,1])
      i0<-which(CPab[,1]<=0)
    }
  } else if (cpmod$dist == "lognormal"){
    # NOTE: lognormal a = sdlog^2 = mod.ln$scale^2, b = meanlog = mod.ln$coef[1]
    CPparms<-MASS::mvrnorm(nsim,c(cpmod$coef[1],cpmod$scale),cpmod$var)
    CPab[,1]<- CPparms[,2]^2 #shape
    CPab[,2]<- CPparms[,1] # scale
  }
  if (option == "parms") return(CPab)
  return(switch(cpmod$dist,
    "exponential" = rexp(nsim, CPab[,1]),
    "weibull"     = rweibull(nsim, shape = CPab[,1], scale = CPab[,2]),
    "loglogistic" = actuar::rllogis(nsim, shape = CPab[,1], scale = CPab[,2]),
    "lognormal"   = rlnorm(nsim, sdlog = sqrt(CPab[,1]), meanlog = CPab[,2])
  ))
}

posteriorL<-function(x, g, pL='jeffreys'){#CDF
# possibilities for prior:
# 'jeffreys', 'uniform', or a custom function(L)
  if (length(x) * length(g) != 1){
    print("error in data")
    return(F)
  }
  if (is.na(as.numeric(x)) * is.na(as.numeric(g)) != 0){
    print("error in data: all values must be numeric")
    return(F)
  }
  if (x < 0 || abs(round(x)-x) > 0.000001){
    print("error in data: x must be a non-negative integer")
    return(F)
  }
  if (g < 0.00001 || g > 0.99999){
    print("error in data: g must be in interval (0, 1)")
    return(F)
  }
  if (!is.character(pL) && !is.function(pL)){
    warning("error in prior")
    return(F)
  }
  ctprob<-0.0001
  mmax<-fmmax(x, g)
  force(x); force(g)
  if (is.function(pL)){    # custom prior provided by user
    # minimal error-checking on the prior
    # improper (integrates to infinity) and non-normalized (does not integrate to 1) priors are allowed
    # prior is truncated above lambda where observing x or smaller would be extremely unlikely
    #   NOTE: this does not provide protection against pathological priors where all the weight
    #   is well above anything remotely compatible with the data
    ctprob<-1e-7
    if (ppois(mmax,100)<ctprob){
      lvec<-seq(0,100,by=0.01)
      Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))] # this is the biggest that lambda could credibly be and still have the biggest M conceivable
    } else {
      i<-0
      while (1){
        i<-i+1
        if (ppois(mmax,100*i)<ctprob){
          lvec<-seq((i-1)*100,i*100,by=0.01)
          Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))]
          break
        }
      }
    }
    Lmin<-1e-5
    deno<-integrate(pL, lower=Lmin, Lmax)$val
    pL0<-function(L) pL(L)/deno   # normalized version of pL, valid in interval [Lmin, Lmax]. Posterior will be set = 0 outside [Lmin, Lmax]
    ans0<-function(L){
    # the numerator of the posterior
      numerator<-(g/(1-g))^x*exp(-L)/factorial(x)*colSums(exp(outer(x:mmax, log((1-g)*L))-lfactorial(x:mmax-x)))*pL0(L)
    # the denominator of the posterior
      denominator<-0
      for (m in x:mmax){
        denominator<-denominator+integrate(function(L) exp(m*log((1-g)*L)-L-lfactorial(m-x)+log(pL0(L))), lower=Lmin, upper=Lmax)$val
      }
      denominator<-denominator*(g/(1-g))^x/factorial(x)
      retval<-numeric(length(L))
      retval[L>=Lmin & L<=Lmax]<-numerator[L>=Lmin & L<=Lmax]/denominator
      retval
    }
    ans<-Vectorize(function(L) integrate(ans0, lower=Lmin, upper=max(L))$val, "L")
  } else if (is.character(pL) && pL %in% c('jeffreys', 'uniform')){
    ans<-switch(pL,
      jeffreys = function(L){
          lfac<-beta(x:mmax + 0.5, 0.5)/sqrt(pi)
          deno<-sum(dbinom(x, size=x:mmax, prob=g)*lfac)
          colSums(dbinom(x, x:mmax, g)*t(outer(L,x:mmax+0.5, FUN="pgamma"))*lfac)/deno },
      uniform = function(L) as.vector(g*(outer(L, 1+x:mmax, FUN="pgamma")%*%dbinom(x,x:mmax,g)))
    )
  } else {
    warning("error in prior")
    return(F)
  }
  ans
}

posteriorLpdf<-function(x, g, pL='jeffreys'){ #CDF
# possibilities for prior:
# 'jeffreys', 'uniform', or a custom function(L)
  if (length(x) * length(g) != 1){
    print("error in data")
    return(F)
  }
  if (is.na(as.numeric(x)) * is.na(as.numeric(g)) != 0){
    print("error in data: all values must be numeric")
    return(F)
  }
  if (x < 0 || abs(round(x)-x) > 0.000001){
    print("error in data: x must be a non-negative integer")
    return(F)
  }
  if (g < 0.00001 || g > 0.99999){
    print("error in data: g must be in interval (0, 1)")
    return(F)
  }
  if (!is.character(pL) && !is.function(pL)){
    warning("error in prior")
    return(F)
  }
  ctprob<-0.0001
  mmax<-fmmax(x, g)
  force(x); force(g)
  if (is.function(pL)){    # custom prior provided by user
    # minimal error-checking on the prior
    # improper (integrates to infinity) and non-normalized (does not integrate to 1) priors are allowed
    # prior is truncated above lambda where observing x or smaller would be extremely unlikely
    #   NOTE: this does not provide protection against pathological priors where all the weight
    #   is well above anything remotely compatible with the data
    ctprob<-1e-7
    if (ppois(mmax,100)<ctprob){
      lvec<-seq(0,100,by=0.01)
      Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))] # this is the biggest that lambda could credibly be and still have the biggest M conceivable
    } else {
      i<-0
      while (1){
        i<-i+1
        if (ppois(mmax,100*i)<ctprob){
          lvec<-seq((i-1)*100,i*100,by=0.01)
          Lmax<-lvec[min(which(ppois(mmax,lvec)<ctprob))]
          break
        }
      }
    }
    Lmin<-1e-5
    deno<-integrate(pL, lower=Lmin, Lmax)$val
    pL0<-function(L) pL(L)/deno   # normalized version of pL, valid in interval [Lmin, Lmax]. Posterior will be set = 0 outside [Lmin, Lmax]
    ans<-function(L){
    # the numerator of the posterior
      numerator<-colSums(outer(x:mmax, L, dpois) * dbinom(x, x:mmax, g)) * pL0(L)
    # the denominator of the posterior
      denominator<-0
      for (m in x:mmax){
        denominator<-denominator+integrate(function(L) dpois(m, L)*dbinom(x, m, g)*pL0(L), lower=Lmin, upper=Lmax)$val
      }
      retval<-numeric(length(L))
      retval[L>=Lmin & L<=Lmax]<-numerator[L>=Lmin & L<=Lmax]/denominator
      retval
    }
  } else if (is.character(pL) && pL %in% c('jeffreys', 'uniform')){
    ans<-switch(pL,
      jeffreys = function(L){
          lfac<-beta(x:mmax + 0.5, 0.5)/sqrt(pi)
          deno<-sum(dbinom(x, size=x:mmax, prob=g)*lfac)
          colSums(dbinom(x, x:mmax, g)*t(outer(L,x:mmax+0.5, FUN="dgamma"))*lfac)/deno },
      uniform = function(L) as.vector(g*(outer(L, 1+x:mmax, FUN="dgamma")%*%dbinom(x,x:mmax,g)))
    )
  } else {
    warning("error in prior")
    return(F)
  }
  ans
}

plotPrior<-function(prior_M, prtype){
#  graphics.off()
  if (.Rvar$platform == 'windows') windows.options(reset = T)
  if (.Rvar$platform == 'mac') quartz.options(reset = T)
  if (.Rvar$platform == 'linux') X11.options(reset = T)
  dev.new(noRStudioGD = T)
  par(fig=c(0,1,0,1))
  par(mar=c(5.1,4.1,4.1,2.1))
  x<-0:(length(prior_M)-1)
  y<-prior_M
  y<-y/sum(y)
  cs<-cumsum(prior_M)
  if (prtype=="Objective"){
    plot(x,y,type='h',xlab = "Fatalities (m)",ylab = "P(M = m)",lwd=2, axes=F, xlim=range(x)+c(0,diff(range(x))*0.04))#,col=cols)
    axis(1)
    axis(2, lab=F)
    points(max(x)+(1:3)*(par('usr')[2]-max(x))/4, (rep(min(y),3)+par('usr')[3])/2, pch=20)
    box()
  } else {
    plot(x,y,type='h',xlab = "Fatalities (m)",ylab = "P(M = m)",lwd=2)#,col=cols)
  }
  assign('devToClose', dev.cur(), env = .Rvar)
  assign('sizeToClose', par('usr'), env = .Rvar)
  if(.Rvar$platform == "windows") bringToTop()
  points(x,y)
  title(paste(prtype,"Prior"))
}
optqL<-function(L,x,g,mmax,p){ #objective function for optimizing L to find the 1-alpha quantile of the posterior
  (pLgX(L,x,g,mmax)-p)^2
}
optqL.ab<-function(L,x,pBa,pBb,mmax,p){ #objective function for optimizing L to find the 1-alpha quantile of the posterior
  (pLgX.ab(L,x,pBa,pBb,mmax)-p)^2
}
pLgX<-function(L,x,g,mmax){ #CDF for posterior (using Jeffrey's prior)
  lfac<-beta(x:mmax+0.5,0.5)/sqrt(pi)
  deno<-sum(dbinom(x,x:mmax,g)*lfac)
  sum(dbinom(x,x:mmax,g)*pgamma(L,x:mmax+0.5)*lfac)/deno
}
pLgX.ab<-function(L,x,pBa,pBb,mmax){ #CDF for posterior (using Jeffrey's prior), where pBa and pBb are parameters for beta distribution estimate of g
  lfac<-beta(x:mmax+0.5,0.5)/sqrt(pi)
  deno<-sum(VGAM::dbetabinom.ab(x,size=x:mmax,shape1=pBa,shape2=pBb)*lfac)
  sum(VGAM::dbetabinom.ab(x,x:mmax,shape1=pBa,shape2=pBb)*pgamma(L,x:mmax+0.5)*lfac)/deno
}
plotPersDist<-function(persistence_distn,Ir,pda,pdb,bmax=0, bmin=0){
# persistence distribution is plotted: mean (red) and uncertainty (dotted)
   if (persistence_distn=="Exponential") {
      x<-seq(0,min(qexp(.99,pda),10*Ir),length=1000)
      y<-1-pexp(x,pda)
      rhat<-signif(1-integrate(pexp,lower=0,upper=Ir,rate=pda)$val/Ir,3)
      subt<-substitute(list(.dist,lambda==.pda,"mean CP"==.mn, .rlab,"I "==.I),list(.dist="Exponential",.pda=signif(pda,3),.mn=signif(1/pda,1),.rlab=bquote(r[I]==.(rhat)),.I=Ir))
      if (bmax!=0){
        y.low<-1-pexp(x,1/bmin)
        y.high<-1-pexp(x,1/bmax)
      }
   } else if (persistence_distn=="Weibull") {
      x<-seq(0,min(qweibull(.995,pda,pdb),10*Ir),length=1000)
      y<-1-pweibull(x,pda,pdb)
      rhat<-signif(1-integrate(pweibull,lower=0,upper=Ir,shape=pda,scale=pdb)$val/Ir,3)
      subt<-substitute(list(.dist,alpha==.pda,beta==.pdb,"mean CP"==.mn, .rlab,"I "==.I),list(.dist="Weibull",.pda=signif(pda,3),.pdb=signif(pdb,3),.mn=round(pdb*gamma(1+1/pda),1),.rlab=bquote(r[I]==.(rhat)),.I=Ir))
      if (bmax!=0){
        y.low<-1-pweibull(x,shape=pda,scale=bmin)
        y.high<-1-pweibull(x,shape=pda,scale=bmax)
      }
   } else if (persistence_distn=="Log-Logistic") {
      x<-seq(0,min(actuar::qllogis(.995,shape=pda,scale=pdb),10*Ir),length=1000)
      y<-1-actuar::pllogis(x,shape=pda,scale=pdb)
      rhat<-signif(1-integrate(actuar::pllogis,lower=0,upper=Ir,shape=pda,scale=pdb)$val/Ir,3)
      if (pda>1){
        subt<-substitute(list(.dist,alpha==.pda,beta==.pdb,"mean CP"==.mn, .rlab,"I "==.I),list(.dist="Log-Logistic",.pda=signif(pda,3),.pdb=signif(pdb,3),.mn=round(pdb*pi/pda/sin(pi/pda),1),.rlab=bquote(r[I]==.(rhat)),.I=Ir))
      } else {
        subt<-substitute(list(.dist,alpha==.pda,beta==.pdb,"mean CP = Inf", .rlab,"I "==.I),list(.dist="Log-Logistic",.pda=signif(pda,3),.pdb=signif(pdb,3),.rlab=bquote(r[I]==.(rhat)),.I=Ir))
      }
      if (bmax!=0){
        y.low<-1-actuar::pllogis(x,shape=pda,scale=bmin)
        y.high<-1-actuar::pllogis(x,shape=pda,scale=bmax)
      }

   } else if (persistence_distn=="Lognormal") {
      x<-seq(0,min(qlnorm(.995,meanlog=pdb, sdlog=sqrt(pda)),10*Ir),length=1000)
      y<-1-plnorm(x,meanlog=pdb, sdlog=sqrt(pda))
      rhat<-signif(1-integrate(plnorm,lower=0,upper=Ir,meanlog=pdb,sdlog=sqrt(pda))$val/Ir,3)
      subt<-substitute(list(.dist,alpha==.pda,beta==.pdb,"mean CP"==.mn, .rlab,"I "==.I),list(.dist="Lognormal",.pda=pda,.pdb=pdb,.mn=round(exp(pdb+pda/2),1),.rlab=bquote(r[I]==.(rhat)),.I=Ir))
      if (bmax!=0){
        y.low<-1-plnorm(x,sdlog=sqrt(pda),meanlog=bmin)
        y.high<-1-plnorm(x,sdlog=sqrt(pda),meanlog=bmax)
      }

   }
   dev.new(noRStudioGD = T)
   plot(x,y,type='l',xlab="Days",ylab="Fraction of carcasses remaining",lwd=2,col=2,ylim=c(0,1),xaxs='i',yaxs='i', xlim=c(0,max(x)))
   if (bmax!=0) {
     lines(x,y.low,lty=3)
     lines(x,y.high,lty=3)
   }
   title("Persistence Distribution")
   if(.Rvar$platform == "windows") bringToTop()
   mtext(subt,cex=.8*.Rvar$charSizeAdjust)
}

feedR.sysc<-function(){
# translate the appropriate TCL data to R variables for analysis of single year (single class) data
# error-checking is assumed to have occurred before calling feedR.sy
# function is used for calculation either from sy or symc module
  for (i in 1:length(syscPrevious)){
    if (!(names(syscPrevious)[i] %in% c("prior_f", "prior_M", "persistence_distn", "arrcomponents")))
      assign(names(syscPrevious[i]),toR(get(paste0("tk.",names(syscPrevious[i])))),pos=1)
  }
  samtype<<-tclvalue(tk.samtype)
  if (samtype=='Formula'){
    Isam<<-as.numeric(tclvalue(tk.Isam))
    nsearch<<-as.numeric(tclvalue(tk.nsearch))
    days<<-(0:nsearch)*Isam
  } else {
    days<<-numeric(length(tk.days))
    for (i in 1:length(tk.days)) days[i]<<-as.numeric(tclvalue(tk.days[[i-1]]))
    Isam<<-round(max(days)/(length(days)-1),1)
    nsearch<<-length(days)-1
  }
  firstsearch<<-tclvalue(tk.firstsearch)
#  persistence_distn<<-persistence_distn
}
chkArr<-function(){
  if (tclvalue(tk.arrfun)=="Uniform"){
    arrok<<-toR(tk.firstsearch)
    return(arrok)
  }
  if (length(tk.arrcomponents) != 3)return(F)
  arrcomponents<-toR(tk.arrcomponents)
  if (is.na(sum(suppressWarnings(as.numeric(arrcomponents))))>0) return(F)
  if (!(arrcomponents[1] %in% 0:1)) return(F)
  if (!(arrcomponents[2] %in% 0:1)) return(F)
  if (!(arrcomponents[3] %in% 0:1)) return(F)
  if (sum(arrcomponents) == 0) return(F)
  if (arrcomponents[1]){
    x<-toR(tk.lwr.u)
    if (length(x)!=1) return(F)
    if (!is.numeric(x))  return(F)
    if (x<0 | x>364) return(F)
    y<-toR(tk.upr.u)
    if (length(y)!=1) return(F)
    if (!is.numeric(y)) return(F)
    if (y<0 | y>364) return(F)
    if (x>=y) return(F)
    w0<-toR(tk.wt.u)
    if (length(w0)!=1) return(F)
    if (!is.numeric(w0)) return(F)
    if (w0<0 | w0>1) return(F)
  }
  if (arrcomponents[2]){
    x<-toR(tk.lwr.p1)
    if (length(x)!=1) return(F)
    if (!is.numeric(x)) return(F)
    if (x<0 | x>364) return(F)
    y<-toR(tk.upr.p1)
    if (length(y)!=1) return(F)
    if (!is.numeric(y)) return(F)
    if (y<0 | y>364) return(F)
    if (x>=y) return(F)
    w1<-toR(tk.wt.p1)
    if (length(w1)!=1)  return(F)
    if (!is.numeric(w1)) return(F)
    if (w1<0 | w1>1) return(F)
    x<-toR(tk.a.p1)
    if (length(x)!=1) return(F)
    if (!is.numeric(x)) return(F)
    if (x<=0) return(F)
    x<-toR(tk.b.p1)
    if (length(x)!=1) return(F)
    if (!is.numeric(x)) return(F)
    if (x<=0) return(F)
  }
  if (arrcomponents[3]){
    x<-toR(tk.lwr.p2)
    if (length(x)!=1) return(F)
    if (!is.numeric(x)) return(F)
    if (x<0 | x>364) return(F)
    y<-toR(tk.upr.p2)
    if (length(y)!=1) return(F)
    if (!is.numeric(y)) return(F)
    if (y<0 | y>364) return(F)
    if (x>=y) return(F)
    w2<-toR(tk.wt.p2)
    if (length(w2)!=1) return(F)
    if (!is.numeric(w2)) return(F)
    if (w2<0 | w2>1) return(F)
    x<-toR(tk.a.p2)
    if (length(x)!=1) return(F)
    if (!is.numeric(x)) return(F)
    if (x<=0) return(F)
    x<-toR(tk.b.p2)
    if (length(x)!=1) return(F)
    if (!is.numeric(x)) return(F)
    if (x<=0) return(F)
  }
  tclvalue(tk.wt.u)<-w0*arrcomponents[1]/(c(w0,w1,w2)%*%arrcomponents)
  tclvalue(tk.wt.p1)<-w1*arrcomponents[2]/(c(w0,w1,w2)%*%arrcomponents)
  tclvalue(tk.wt.p2)<-w2*arrcomponents[3]/(c(w0,w1,w2)%*%arrcomponents)
  return(startchk(toR(tk.firstsearch)))
}


getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}
cellind<-function(i,j) as.tclObj(paste0(i, ',', j))
# returns a table cell index as a tclvalue
# useful for tagging cells via, e.g., tcl(classTable,"tag","celltag","error", cellind(2,3))
abgrCP<-function(pd,r,Ir,CP){
   amax<-0
   if (pd=="Weibull"){
      CP<-max(CP,r*Ir/.99)
      for (a in 1:40){ #find an appropriate upper bound on alpha
         if(1-integrate(pweibull,lower=0,upper=Ir,shape=a,scale=CP/gamma(1+1/a))$val/Ir>r){
            amax<-a
            break
         }
      }
      a<-optimize(f=f.wei.a,si=Ir,mu=CP,rI=r,interval=c(0.01,amax))$min
      b<-CP/gamma(1+1/a)
   }  else if (pd=="Lognormal"){
      CP<-max(CP,r*Ir/.999)
      a<-seq(0.001,20,length=1000)
      imax<-2
      for (i in 1:length(a)){ #find an appropriate bounds on alpha
         if(1-integrate(plnorm,lower=0,upper=Ir,meanlog=log(CP)-a[i]/2,sdlog=sqrt(a[i]))$val/Ir<r){
            imax<-i
            break
         }
      }
      a<-optimize(f=f.lnorm.a,si=Ir,mu=CP,rI=r,interval=c(a[imax-1],a[imax]))$min
      b<-log(CP)-a/2
   } else if (pd=="Log-Logistic"){
      CP<-max(CP,r*Ir/.999)
      for (a in 1:40){ #find an appropriate upper bound on alpha
         if(1-integrate(actuar::pllogis,lower=0,upper=Ir,shape=a,scale=CP*sin(pi/a)*a/pi)$val/Ir>r){
            amax<-a
            break
         }
      }
      a<-optimize(f=f.ll.a,si=Ir,mu=CP,rI=r,interval=c(max(0.01,amax-1),amax))$min
      b<-CP*sin(pi/a)*(a/pi)
   } else if (pd=="Exponential"){
      a<-optimize(f=f.e.a,si=Ir,mu=CP,rI=r,interval=c(0.01,20))$min
      b<-1/a
   }
   return(c(a,b))
}
f.wei.a<-function(a,si,mu,rI){# objective function to minimize in finding Weibull alpha parameter given r and meanCP
   abs(1-integrate(pweibull,lower=0,upper=si,shape=a,scale=mu/gamma(1+1/a))$val/si-rI)
}
f.lnorm.a<-function(a,si,mu,rI){# objective function to minimize in finding Weibull alpha parameter given r and meanCP
   abs(1-integrate(plnorm,lower=0,upper=si,meanlog=log(mu)-a/2,sdlog=sqrt(a))$val/si-rI)
}
f.ll.a<-function(a,si,mu,rI){# objective function to minimize in finding Weibull alpha parameter given r and meanCP
   abs(1-integrate(actuar::pllogis,lower=0,upper=si,shape=a,scale=mu*a/pi*sin(pi/a))$val/si-rI)
}
f.e.a<-function(a,si,mu,rI){# objective function to minimize in finding Weibull alpha parameter given r and meanCP
   abs(1-integrate(pexp,lower=0,upper=si,rate=a)$val/si-rI)
}
setT<-function(tabdat,rowm,colm) tcl(tabdat,"tag", "celltag", "cellok",as.tclObj(paste0(rowm,',', colm)))  #
setF<-function(tabdat,rowm,colm) tcl(tabdat,"tag", "celltag", "error",as.tclObj(paste0(rowm,',', colm)))
##########################
toR<-function(dat){ # convert a numeric tcl array, vector, or scalar to an R variable of the same dimension
  n<-length(dat)
  if (n==0) return(NA)
  if (is.numeric(dat)) return(dat)
  if ("tclArray" %in% class(dat)){
    if (length(grep(',',names(dat))) == 0){ # vector
      ans<-numeric(n)
      for (i in 1:n)
        ans[i]<-ifelse(is.na(suppressWarnings(as.numeric(tclvalue(dat[[i-1]])))),
          ifelse(class(try(tclvalue(dat[[i-1]]),silent=T))=="try-error",
            NA, tclvalue(dat[[i-1]])),
          as.numeric(dat[[i-1]]))
      return(ans)
    }
    if (length(grep(',', names(dat))) == n){ # 2-d array
      dims<-apply(array(unlist(lapply(strsplit(names(dat),','),"as.numeric")),dim=c(2,length(dat))),F=max,M=1)+1
      ans<-array(dim=dims)
      for (i in 1:dims[1])
        for (j in 1:dims[2])
          ans[i,j]<-ifelse(is.na(suppressWarnings(as.numeric(dat[[i-1,j-1]]))),
            tclvalue(dat[[i-1,j-1]]), as.numeric(dat[[i-1,j-1]]))
      return(ans)
    }
    return(NA) # neither vector or array
  }
  # single number (or word)
  ans0<-try(tclvalue(dat),silent=T)
  if (class(ans0)=="try-error") return(NA)
  if (nchar(ans0)==0) return(NA)
  ans<-ifelse(is.na(suppressWarnings(as.numeric(ans0))), ans0, as.numeric(ans0))
  return(ans)
}
saveparm<-function(dat){
# write the given data set to an .rds file for later reading (R list)
# no error-checking...this enables users to save bad parameter sets to revisit later
# the parameters for the different persistence distributions are not saved
# saved data are read from R list (which not been error-checked)
  if (!exists('csvpath', env = .Rvar)) assign('csvpath', getwd(), env = .Rvar)
  filename <- tclvalue(tkgetSaveFile(filetypes = "{{R data files} {.rds}}", defaultextension = ".rds", initialfile = '.rds', title = "Save", initialdir = .Rvar$csvpath))
  tmp<-unlist(strsplit(filename,'/')); pathname<-paste(tmp[-length(tmp)],collapse='/')
  if (nchar(pathname)>0) .Rvar$csvpath <- pathname
  if (filename == "") return(FALSE)
  saveRDS(dat, file=filename)
  .Rvar$dataFileTitle <- tmp[length(tmp)]
  return(TRUE)
}
plot_arrivals <- function(arriv.dist, arrstart, arrcomponents,
  lwr.u, upr.u, wt.u,
  lwr.p1,upr.p1,wt.p1, a.p1, b.p1,
  lwr.p2,upr.p2,wt.p2, a.p2, b.p2){
  if (arriv.dist =='Uniform') {
    return(F)
  } else if (sum(arrcomponents)==0){
    tkmessageBox(message="No model componenents have been defined.\n Click 'Edit' to build a model or select 'Uniform' for simple model.")
  } else {         #  otherwise... compound arrival function
    xlen <- 1000
    xx<-seq(0,364,length=xlen)
    yy.u<-wt.u/(upr.u-lwr.u)*(xx>=lwr.u & xx<=upr.u)
    yy.p1<-numeric(xlen); ind<-which(xx>=lwr.p1 & xx<=upr.p1)
    yy.p1[ind]<-wt.p1*dbeta((xx[ind]-lwr.p1)/(upr.p1-lwr.p1), shape1=a.p1, shape2=b.p1)/(upr.p1-lwr.p1)
    yy.p2<-numeric(xlen); ind<-which(xx>=lwr.p2 & xx<=upr.p2)
    yy.p2[ind]<-wt.p2*dbeta((xx[ind]-lwr.p2)/(upr.p2-lwr.p2), shape1=a.p2, shape2=b.p2)/(upr.p2-lwr.p2)
    yy<-yy.u*arrcomponents[1]+yy.p1*arrcomponents[2]+yy.p2*arrcomponents[3]
    ymax<-max(yy)
    par(mar=c(4,2.5,.5,1.5),mgp=c(2,.7,0))
    plot(xx,yy,type='n',xlab = "Date",ylab = "",axes=F,xaxs='i',yaxs='i',ylim=c(0,ymax*1.06))
    mtext(side=2,line=.7,"Relative Arrival Rate", cex=.Rvar$charSizeAdjust)
    if (arrcomponents[1]) lines(xx,yy.u,col=uclr)
    if (arrcomponents[2]) lines(xx,yy.p1,col=p1clr)
    if (arrcomponents[3]) lines(xx,yy.p2,col=p2clr)
    lines(xx,yy,lwd=2)
    ats <- (ats0-arrstart)%%365
    axis(1,at=ats,lab=figlab)
    if(.Rvar$platform == "windows") bringToTop()
    box();
  }
}
rCPgab<-function(pd,a,b,Ir){
   if (pd=="Weibull"){
      CP<-b*gamma(1+1/a)
      r<-1-integrate(pweibull,lower=0,upper=Ir,shape=a,scale=b)$val/Ir
   } else if (pd=="Lognormal") {
      CP<-exp(b+a/2)
      r<-1-integrate(plnorm,lower=0,upper=Ir,meanlog=b,sdlog=sqrt(a))$val/Ir
   } else if (pd=="Log-Logistic") {
      if (a>1) {
        CP<-b/a*pi/sin(pi/a)
      } else {
        CP<- Inf
      }
      r<-1-integrate(actuar::pllogis,lower=0,upper=Ir,shape=a,scale=b)$val/Ir
   }
   if (pd=="Exponential"){
    CP<-b;
    r<-1-integrate(pexp,lower=0,upper=Ir,rate=1/b)$val/Ir
   }
   return(c(CP,r))
}

checkmyr<-function(myr){
# myr matches the format of the tcltable in the edit priors frame (i.e., first row has names, first column has years)
# the years are not included in the error check data array
  if (is.null(dim(myr))){
    tkmessageBox(message="Error: data must be 2-d array with columns for year, rel_wt, carcass count, Ba, and Bb")
    return(FALSE)
  }
  rowm<-dim(myr)[1]; colm<-dim(myr)[2]
  if (colm != 5) {
    tkmessageBox(message="Error: must be columns for year, rel_wt, carcass count, Ba, and Bb")
    return(FALSE)
  }
  NAtot<-sum(is.na(suppressWarnings(as.numeric(myr))))
  if (NAtot==colm & rowm==1){
    tkmessageBox(message="Error: no data")
    return(FALSE)
  }
  if (NAtot > 0){ # parse the missing data: full rows or additional data too?
    if(NAtot%%colm==0){ # data must be numeric (but empty rows at the end are OK)
      if (sum(is.na(as.vector(t(myr))[(length(myr)-NAtot+1):length(myr)]))==NAtot){ # empty rows at bottom of table but no other missing values
        myr<-myr[1:(rowm-NAtot/colm),] #DEBUG <<   :production <
      }
      if (sum(is.na(myr))>0){
        tkmessageBox(message="Error: Missing or non-numeric data")
        return(FALSE)
      }
    } else {
      tkmessageBox(message="Error: Missing or non-numeric data")
#      rowind<-min(which(is.na(myr),arr.ind=T)[,1])
#      if (length(rowind)==1) colind<-which(is.na(myr),arr.ind=T)[2] else colind<-min(which(is.na(myr),arr.ind=T)[rowind,2]) # indices for the first NA
      return(FALSE)
    }
  }
  if (is.null(dim(myr))) myr<-array(myr,dim=c(1,length(myr))) # turns vector into 1-d array (not sure what the point is...)
  if (sum(myr[,2] < 0) > 0){
    tkmessageBox(message="Error: Relative weights (rel_wt) must be non-negative")
    return(FALSE)
  }
  if (sum(myr[,3]< 0) > 0 | sum(abs(myr[,3]-round(myr[,3]))>0.0001)){ # carcass counts must be non-negative integers
    tkmessageBox(message="Error: Carcass counts must be whole numbers")
    return(FALSE)
  } else {
    myr[,3]<-round(myr[,3])
  }
  if (sum(myr[,4]<=0) > 0){ # Ba must be positive
    tkmessageBox(message="Error: Ba parameters must be positive")
    return(FALSE)
  }
  if (sum(myr[,5]<=0)){ # Bb parameters must be positive
    tkmessageBox(message="Error: Bb parameters must be positive")
    return(FALSE)
  }
  return(myr)
}
conversionsCalculator<-function(){
  w<-8
  chkg<-function(v){
    g<-suppressWarnings(as.numeric(tclvalue(v)))
    if (is.null(g) || is.na(g) || length(g) != 1 || nchar(g)==0 || g <= 0 || g >= 1){
      return(F)
    } else {
      return(T)
    }
  }
  chkab<-function(v){
    ab<-suppressWarnings(as.numeric(tclvalue(v)))
    if (is.null(ab) || is.na(ab) || length(ab) != 1 || nchar(ab)==0 || ab <= 0){
      return(F)
    } else {
      return(T)
    }
  }
  calculator<-tktoplevel()
  tkgrab.set(calculator);  tkfocus(calculator)
  tkwm.title(calculator,paste0("EoA, v", .Rvar$VER," - Parameter Conversion"))
  tkwm.resizable(calculator,0,0)
  tkwm.deiconify(calculator)

  gcoFrame<-ttklabelframe(calculator, text="   Detection Probability Parameters")
  glbl<-tklabel(gcoFrame, text = "g\u0302")
  cilbl<-tklabel(gcoFrame, text = "95% CI")
  Balbl<-tklabel(gcoFrame, text = "Ba")
  Bblbl<-tklabel(gcoFrame, text = "Bb")
  tk.gco<-tclVar(0.3)
  gco.edit<-tkentry(gcoFrame, textvariable = tk.gco, bg='white', justify = 'right',width=w)
  tk.glco<-tclVar(0.15)
  glco.edit<-tkentry(gcoFrame, textvariable = tk.glco, bg='white', justify = 'right',width=w)
  tk.guco<-tclVar(0.45)
  mu<-toR(tk.gco)
  sig2<-((toR(tk.guco)-toR(tk.glco))/4)^2
  betaa<-mu^2/sig2*(1-mu)-mu; betab<-betaa*(1/mu-1)
  guco.edit<-tkentry(gcoFrame, textvariable = tk.guco, bg='white', justify = 'right',width=w)
  tk.Ba<-tclVar(signif(betaa,5))
  Ba.edit<-tkentry(gcoFrame, textvariable = tk.Ba, bg='white', justify = 'right',width=w)
  tk.Bb<-tclVar(signif(betab,5))
  Bb.edit<-tkentry(gcoFrame, textvariable = tk.Bb, bg='white', justify = 'right',width=w)
  ######
  tkgrid(glbl, column = 1, row = 0)
  tkgrid(gco.edit, column = 2, row = 0, padx=c(3,10))
  tkgrid(cilbl, column = 3, row = 0)
  tkgrid(glco.edit, column = 4, row = 0)
  tkgrid(guco.edit, column = 5, row = 0, padx=c(0,10))
  tkgrid(Balbl, column = 1, row = 1)
  tkgrid(Ba.edit, column = 2, row = 1, padx=c(3,10))
  tkgrid(Bblbl, column = 3, row = 1, sticky='e',padx=c(10,3))
  tkgrid(Bb.edit, column = 4, row = 1)
  tkgrid(gcoFrame,pady=10,sticky='w')

  tkbind(gco.edit,"<KeyRelease>", function(){
    if (!chkg(tk.gco)){
      tkconfigure(gco.edit,bg=colors()[652])
      tclvalue(tk.Ba)<-"NA"
      tclvalue(tk.Bb)<-"NA"
      tkconfigure(Ba.edit, bg = 'yellow')
      tkconfigure(Bb.edit, bg = 'yellow')
      return(F)
    }
    tkconfigure(gco.edit,bg='white')
    g<-toR(tk.gco)
    if (chkg(tk.glco) && chkg(tk.guco)){ # are glco and guco reals in (0, 1)?
      glco<-toR(tk.glco)
      guco<-toR(tk.guco)
      if (glco < g && guco > g){
        tkconfigure(gco.edit,bg='white')
        tkconfigure(glco.edit,bg='white')
        tkconfigure(guco.edit,bg='white')
        mu<-toR(tk.gco)
        sig2<-((toR(tk.guco)-toR(tk.glco))/4)^2
        betaa<-mu^2/sig2*(1-mu)-mu; betab<-betaa*(1/mu-1)
        tclvalue(tk.Ba)<-signif(betaa,5)
        tclvalue(tk.Bb)<-signif(betab,5)
        tkconfigure(Ba.edit, bg = 'white')
        tkconfigure(Bb.edit, bg = 'white')
        return(T)
      } else {
        tclvalue(tk.Ba)<-"NA"
        tclvalue(tk.Bb)<-"NA"
        tkconfigure(Ba.edit, bg = 'yellow')
        tkconfigure(Bb.edit, bg = 'yellow')
        tkconfigure(gco.edit,bg=colors()[400])
        tkconfigure(glco.edit,bg=colors()[400])
        tkconfigure(guco.edit,bg=colors()[400])
      }
    }
  })

  tkbind(glco.edit,"<KeyRelease>", function(){
    if (!chkg(tk.glco)){
      tkconfigure(glco.edit,bg=colors()[652])
      tclvalue(tk.Ba)<-"NA"
      tclvalue(tk.Bb)<-"NA"
      tkconfigure(Ba.edit, bg = 'yellow')
      tkconfigure(Bb.edit, bg = 'yellow')
      return(F)
    }
    tkconfigure(glco.edit,bg='white')
    glco<-toR(tk.glco)
    if (chkg(tk.gco) && chkg(tk.guco)){
      g<-toR(tk.gco)
      guco<-toR(tk.guco)
      if (glco < g && guco > g){
        tkconfigure(gco.edit,bg='white')
        tkconfigure(glco.edit,bg='white')
        tkconfigure(guco.edit,bg='white')
        mu<-toR(tk.gco)
        sig2<-((toR(tk.guco)-toR(tk.glco))/4)^2
        betaa<-mu^2/sig2*(1-mu)-mu; betab<-betaa*(1/mu-1)
        tclvalue(tk.Ba)<-signif(betaa,5)
        tclvalue(tk.Bb)<-signif(betab,5)
        tkconfigure(Ba.edit, bg = 'white')
        tkconfigure(Bb.edit, bg = 'white')
        return(T)
      } else {
        tclvalue(tk.Ba)<-"NA"
        tclvalue(tk.Bb)<-"NA"
        tkconfigure(Ba.edit, bg = 'yellow')
        tkconfigure(Bb.edit, bg = 'yellow')
        tkconfigure(gco.edit,bg=colors()[400])
        tkconfigure(glco.edit,bg=colors()[400])
        tkconfigure(guco.edit,bg=colors()[400])
      }
    }
  })
  tkbind(guco.edit,"<KeyRelease>", function(){
    if (!chkg(tk.guco)){
      tkconfigure(guco.edit,bg=colors()[652])
      tclvalue(tk.Ba)<-"NA"
      tclvalue(tk.Bb)<-"NA"
      tkconfigure(Ba.edit, bg = 'yellow')
      tkconfigure(Bb.edit, bg = 'yellow')
      return(F)
    }
    tkconfigure(guco.edit,bg='white')
    guco<-toR(tk.guco)
    if (chkg(tk.gco) && chkg(tk.glco)){
      glco<-toR(tk.glco)
      g<-toR(tk.gco)
      if (glco < g && guco > g){
        tkconfigure(gco.edit,bg='white')
        tkconfigure(glco.edit,bg='white')
        tkconfigure(guco.edit,bg='white')
        mu<-toR(tk.gco)
        sig2<-((toR(tk.guco)-toR(tk.glco))/4)^2
        betaa<-mu^2/sig2*(1-mu)-mu; betab<-betaa*(1/mu-1)
        tclvalue(tk.Ba)<-signif(betaa,5)
        tclvalue(tk.Bb)<-signif(betab,5)
        tkconfigure(Ba.edit, bg = 'white')
        tkconfigure(Bb.edit, bg = 'white')
        return(T)
      } else {
        tclvalue(tk.Ba)<-"NA"
        tclvalue(tk.Bb)<-"NA"
        tkconfigure(Ba.edit, bg = 'yellow')
        tkconfigure(Bb.edit, bg = 'yellow')
        tkconfigure(gco.edit,bg=colors()[400])
        tkconfigure(glco.edit,bg=colors()[400])
        tkconfigure(guco.edit,bg=colors()[400])
      }
    }
  })

  tkbind(Ba.edit,"<KeyRelease>", function(){
    if (!chkab(tk.Ba)){
      tkconfigure(Ba.edit,bg=colors()[652])
      tclvalue(tk.gco)<-'NA'
      tclvalue(tk.glco)<-'NA'
      tclvalue(tk.guco)<-'NA'
      tkconfigure(gco.edit,bg='yellow')
      tkconfigure(glco.edit,bg='yellow')
      tkconfigure(guco.edit,bg='yellow')
      return(F)
    }
    Ba<-toR(tk.Ba)
    tkconfigure(Ba.edit,bg='white')
    if (chkab(tk.Bb)){
      Bb<-toR(tk.Bb)
      tclvalue(tk.gco) <-signif(Ba/(Ba+Bb),5)
      tclvalue(tk.glco)<-signif(qbeta(.025,Ba,Bb),5)
      tclvalue(tk.guco)<-signif(qbeta(.975,Ba,Bb),5)
      tkconfigure(gco.edit,bg='white')
      tkconfigure(glco.edit,bg='white')
      tkconfigure(guco.edit,bg='white')
    }
  })

  tkbind(Bb.edit,"<KeyRelease>", function(){
    if (!chkab(tk.Bb)){
      tkconfigure(Bb.edit,bg=colors()[652])
      tclvalue(tk.gco)<-'NA'
      tclvalue(tk.glco)<-'NA'
      tclvalue(tk.guco)<-'NA'
      tkconfigure(gco.edit,bg='yellow')
      tkconfigure(glco.edit,bg='yellow')
      tkconfigure(guco.edit,bg='yellow')
      return(F)
    }
    tkconfigure(Bb.edit,bg='white')
    Bb<-toR(tk.Bb)
    if (chkab(tk.Ba)){
      Ba<-toR(tk.Ba)
      tclvalue(tk.gco) <-signif(Ba/(Ba+Bb),5)
      tclvalue(tk.glco)<-signif(qbeta(.025,Ba,Bb),5)
      tclvalue(tk.guco)<-signif(qbeta(.975,Ba,Bb),5)
      tkconfigure(gco.edit,bg='white')
      tkconfigure(glco.edit,bg='white')
      tkconfigure(guco.edit,bg='white')
    }
  })
  OKbutton<-tkbutton(calculator, text="Done", command=function()tkdestroy(calculator))
  tkgrid(OKbutton, row=1, column=1)
}

val.numeric<-function(v) if (length(v) != 1 || is.na(v) || nchar(v) == 0) return(F) else return(T)
val.integer<-function(v) if (abs(round(v)-v) > 0.0000001) return(F) else return(T)
val.gt0<-function(v)  if (v <= 0.0000001) return(F) else return(T)
val.gte0<-function(v) if (v < 0) return(F) else return(T)
val.lte1<-function(v) if (v > 1) return(F) else return(T)
val.lt1<-function(v)  if (v >= .9999999) return(F) else return(T)

postL.sumry<-function (x, Ba, Bb){
    ctprob <- 1e-04
    mmax <- fmmax.ab(x, Ba, Bb)
    if (ppois(mmax, 100) < ctprob) {
        lvec <- seq(0, 100, by = 0.01)
        Lmax <- lvec[min(which(ppois(mmax, lvec) < ctprob))]
    }
    else {
        i <- 1
        while (1) {
            i <- i + 1
            if (ppois(mmax, 100 * i) < ctprob) {
                lvec <- seq((i - 1) * 100, i * 100, by = 0.01)
                Lmax <- lvec[min(which(ppois(mmax, lvec) < ctprob))]
                break
            }
        }
    }
    cpLgX <- function(L, x, Ba, Bb, mmax) 1 - pLgX.ab(L, x, Ba,
        Bb, mmax)
    c(
      optimize(f = optqL.ab, interval = c(1e-05, Lmax), x = x, Ba, Bb, mmax = mmax, p = .025)$min,
      optimize(f = optqL.ab, interval = c(1e-05, Lmax), x = x, Ba, Bb, mmax = mmax, p = 0.25)$min,
      optimize(f = optqL.ab, interval = c(1e-05, Lmax), x = x, Ba, Bb, mmax = mmax, p = 0.5)$min,
      optimize(f = optqL.ab, interval = c(1e-05, Lmax), x = x, Ba, Bb, mmax = mmax, p = 0.75)$min,
      optimize(f = optqL.ab, interval = c(1e-05, Lmax), x = x, Ba, Bb, mmax = mmax, p = 0.975)$min
    )
}
combg<-function(Ba, Bb){
  n<-length(Ba)
  mu<-sum(Ba/(Ba+Bb))/n
  sig2<-sum(Ba*Bb/((Ba+Bb)^2*(Ba+Bb+1)))/n^2
  Bab<-numeric(2)
  Bab[1]<-mu^2/sig2*(1-mu)-mu; Bab[2]<-Bab[1]*(1/mu-1)
  Bab
}
#combg.w<-function(Ba, Bb, rho){
#  rho<-rho/sum(rho)
#  mu<-Ba/(Ba+Bb)
#  sig2<-Ba*Bb/((Ba+Bb)^2*(Ba+Bb+1))
#  Eg<-sum(mu*rho)
#  Vg<-sum(sig2*rho^2)#+(sum(aa[1:i]*(mu[1:i]^2+sig2[1:i]))-Eg[i]^2)/lambda[i]
#  Bab<-numeric(2)
#  Bab[1]<-Eg^2/Vg*(1-Eg)-Eg; Bab[2]<-Bab[1]*(1/Eg-1) # these are the two shape parameters for the beta distribution underlying the beta-binomial for X | M
#  Bab
#}
combg.w<-function(Ba, Bb, rho){
# input: parallel vectors for detection probability beta parameters (Ba and Bb) and relative weights (rho)
# alternatively: rho can be an array of vectors with dim = c(nsim, length(Ba))
  if (!is.array(rho)) {
    scalar<-T
    rho<-array(rho, dim=c(1, length(rho)))
  } else {
    scalar<-F
  }
  rho<-rho/rowSums(rho)
  mu<-Ba/(Ba+Bb)
  sig2<-Ba*Bb/((Ba+Bb)^2*(Ba+Bb+1))
  Eg<-rho%*%mu
  Vg<-rho^2%*%sig2#+(sum(aa[1:i]*(mu[1:i]^2+sig2[1:i]))-Eg[i]^2)/lambda[i]
  Bab<-array(dim=c(dim(rho)[1],2))
  Bab[,1]<-Eg^2/Vg*(1-Eg)-Eg; Bab[,2]<-Bab[,1]*(1/Eg-1) # these are the two shape parameters for the beta distribution underlying the beta-binomial for X | M
  if (scalar) as.vector(Bab) else Bab
}
breakoutProbs<-function(Tau, alpha, pBa, pBb){
  xcrit<-0 # smallest x that triggers
  while(1){
    mmax<-fmmax.ab(xcrit, pBa, pBb)
    M<-xcrit:mmax
    pXgM<-VGAM::dbetabinom.ab(xcrit,size=M,shape1=pBa,shape2=pBb) # the probabilities of X for M = 0:mmax
    pM<-diff((sqrt(c(xcrit-1,xcrit:mmax)+1)))
    pMgX<-pXgM*pM; pMgX<-pMgX/sum(pMgX) # posterior distribution for M (ignoring M < x, which has probability = zero)
    if (M[min(which(1-cumsum(pMgX)<1-alpha))] > Tau) break
    xcrit<-xcrit+1
  }
  p.notrigger <- VGAM::pbetabinom.ab(xcrit-1, Tau:mmax, pBa, pBb)
  maxind<-min(which(p.notrigger<0.01))
  list(M = Tau:M[maxind], p.notriger = p.notrigger[1:maxind])
}

rhotest<-function(mydat){
  # rho test...
  # H0: rho is defined correctly
  # likelihood ratio test
  x<-mydat$X; Ba<-mydat$Ba; Bb<-mydat$Bb; rho<-mydat$rel_wt; alpha<-1-mydat$crlev
  nx<-length(x)
  rho.norm<-rho/sum(rho)
  pXgL0<-numeric(nx)
  Bab<-combg.w(Ba, Bb, rho)
  maxLhat<-5*sum(x+.5)/qbeta(0.1, Bab[1], Bab[2])
  lik0<-ifelse(sum(x) == 0, 1,
    optimize(function(L){
      for (yi in 1:nx){
        mmax<-fmmax.ab(x[yi], Ba[yi], Bb[yi])
        pXgL0[yi]<<-sum(VGAM::dbetabinom.ab(x[yi], x[yi]:mmax, Ba[yi], Bb[yi])*dpois(x[yi]:mmax, L*rho.norm[yi]))
      }
      prod(pXgL0)
    }, lower=0.001, upper=maxLhat, maximum=T)$obj
  )
  # likelihood of fitted model (maximum likelihoods year-by-year fitting lambdas separately)
  pXgL.fit<-numeric(nx)
  for (yi in 1:nx){
    if (x[yi]==0) pXgL.fit[yi]<-1 else { # if x[yi] = 0, then lambda[yi] = 0 is MLE
      mmax<-fmmax.ab(x[yi], Ba[yi], Bb[yi])
      maxLhat<-5*x[yi]/qbeta(0.1, Ba[yi], Bb[yi])
      pXgL.fit[yi]<-optimize(function(L) sum(VGAM::dbetabinom.ab(x[yi], x[yi]:mmax, Ba[yi], Bb[yi])*dpois(x[yi]:mmax, L)), lower=0, upper=maxLhat, maximum=T)$obj
    }
  }
  lik.fit<-prod(pXgL.fit)
  pval<-1-pchisq(2*(log(lik.fit)-log(lik0)), nx-1)

  if (length(mydat$X)==1)
    list(pval = 1, rho.ass = mydat$rel_wt, rho.obs = mydat$rel_wt, quicktest = 0, qtpval = 1) # quicktest is a point estimate of the bias due to misspecification of rho
  nsim<-10000
  Lposts<-list()
  Lmax<-numeric(nx)
  L.obs<-numeric(nx)
  rho.obs<-numeric(nx)
  for (yi in 1:nx){ # calculate posterior distributions
    Lposts[[yi]]<-posteriorL.ab(x[yi], Ba[yi], Bb[yi])
    Lmax[yi]<-uniroot(function(L) Lposts[[yi]](L) - 0.999999, interval=c(0, 1e7))$root
    L.obs[yi]<-uniroot(function(L) Lposts[[yi]](L) - (1-alpha), interval=c(0, Lmax[yi]))$root
  }
  rho.obs<-L.obs/sum(L.obs)
#  rho.teststat<-sum((rho.obs-rho/sum(rho))^2)
#  # if H0, then overall posterior lambda can be calculated by combining classes according to assumed rho
  Bab<-combg.w(Ba, Bb, rho)
  Lpost0<-posteriorL.ab(sum(x), Bab[1], Bab[2]) # posteriorL under H0

  # simulated rho's
  nsim<-1000
  Lsim<-array(dim=c(nsim, length(mydat$X)))
  for (i in 1:length(mydat$X)) Lsim[,i]<-rL(nsim, Lposts[[i]])
  # scale to sum to 1
  Lsim<-Lsim/rowSums(Lsim)
  rhosim<-Lsim*sum(mydat$rel_wt)

#  # sample from the H0 posterior and disperse the rate among the nx years according to H0
#  Lsim0<-numeric(nsim)
#  Lsim<-array(dim=c(nsim,nx))
#  Lmax0<-uniroot(function(L) Lpost0(L) - 0.999999, interval=c(0, 1e10))$root
#  Lvec<-seq(0, Lmax0, length=1000)
#  ends<-Lpost0(Lvec)
#  ends[length(ends)]<-1 # truncate the distribution
#  r<-runif(nsim)
#  rint<-findInterval(r, ends)
#  Lsim0<-Lvec[rint]+(r-ends[rint])*(Lvec[rint+1]-Lvec[rint])/(ends[rint+1]-ends[rint])

  # bias quick test
  Bab.f<-combg.w(Ba, Bb, rho.obs)
  Lpostf<-posteriorL.ab(sum(x), Bab.f[1], Bab.f[2])
  lf<-uniroot(function(L) Lpostf(L)-(1-alpha), interval=c(0,1e6))$root
  l0<-uniroot(function(L) Lpost0(L)-(1-alpha), interval=c(0,1e6))$root
#.Rvar$Lsimf<-Lsimf
#.Rvar$Lsim0<-Lsim0
  list(pval = pval, rho.ass = rho, rho.obs = rho.obs, rhoqtls = apply(rhosim, F=quantile, M=2, prob=c(0.025, 0.25, 0.5, 0.75, 0.975), type=3), quicktest = l0/lf) # quicktest is a "point estimate" of the bias due to misspecification of rho
}

rL<-function(nsim, Lcdf, exact = F){
  r<-runif(nsim)
  Lmax<-uniroot(function(L) Lcdf(L) - 0.999999, interval=c(0, 1e10))$root
  if (!exact){
    Lvec<-seq(0, uniroot(function(L) Lcdf(L) - 0.999999, interval=c(0, 1e10))$root, length=1000) # sample from a discrete representation of L
    ends<-Lcdf(Lvec)
    ends[length(ends)]<-1 # truncate the distribution
    rint<-findInterval(r, ends)
    return(Lvec[rint]+(r-ends[rint])*(Lvec[rint+1]-Lvec[rint])/(ends[rint+1]-ends[rint]))
  }
  ans<-numeric(nsim)
  for (i in 1:nsim){
    ans[i]<-ifelse(r[i] < 0.9999989, uniroot(function(L) Lcdf(L) - r[i], interval=c(0, 1e10))$root, Lmax)
  }
  ans
}
g2ab<-function(g, glwr, gupr){
  if (glwr < g && gupr > g){
    mu<-g
    sig2<-((gupr-glwr)/4)^2
    betaa<-mu^2/sig2*(1-mu)-mu; betab<-betaa*(1/mu-1)
    return(c(betaa, betab))
  } else {
    print("error in data")
    return(F)
  }
}
getData<-function(datname) return(get(datname,env = .Rvar))
arrivalModel<-R6::R6Class("arrivalModel",
  public = list(
    arrcomponents = c(T, T, T), # which components of the compound arrival function are included?
    lwr.u=0,    upr.u =365, wt.u =1/4, # parameters for uniform component
    lwr.p1=36,  upr.p1=146, wt.p1=1/4, a.p1=2.5, b.p1=2.5, # parameters for beta1 component
    lwr.p2=243, upr.p2=320, wt.p2=1/2, a.p2=3.5, b.p2=4.5, # parameters for beta2 component
    monLwr = 0, monUpr = 365, Mdo=F,# monitoring period (on arrivals' scale)
#    searchDays = NULL, # search schedule: vector of increasing times beginning at zero (or NULL)
#    s0 = NULL, # beginning of monitoring season, i.e., monitoring begins s0 days after arrivals begin
    duration = NULL, # How long the arrival + monitoring lasts; may be specified expressly or calculated implicitly by internal calculation
    # NOTE: the beginning of the arrival season is assumed to be zero; season extends to 'duration' (arrivals may be zero after s0 + max(searchDays))
    arrfun = NULL,
    initialize = function(parms = NULL, duration = NULL){
    # 'parms' contains a list of parameters that match the format of compound arrival function parameters
    # Optionally, may include searchDays and/or s0
      # If parms are provided, write the values into the data
#######   NOTE: data-checking for internal consistency of parameter values is required (but not yet implemented)
# arrcomponents must be logical
# parms must be numeric, non-negative, scalars
# search days must be a strictly increasing, numeric vector that starts at zero (or NULL)
# s0 must be non-negative, numeric (or NULL)
# duration: must be at least max(the component upper bounds, searchDays + s0) [assuming that calculation makes sense b/c data are present)
      if (!missing(parms)){
        if (!is.list(parms) && !("R6" %in% inherits(parms))) {warning("Bogus parms. Aborting."); return}
        self$arrcomponents <- parms$arrcomponents
        self$lwr.u <- parms$lwr.u; self$upr.u <- parms$upr.u; self$wt.u <- parms$wt.u # parameters for uniform component
        self$lwr.p1 <- parms$lwr.p1;  self$upr.p1 <- parms$upr.p1; self$wt.p1 <- parms$wt.p1; self$a.p1=parms$a.p1; self$b.p1 <- parms$b.p1
        self$lwr.p2 <- parms$lwr.p2;  self$upr.p2 <- parms$upr.p2; self$wt.p2 <- parms$wt.p2; self$a.p2=parms$a.p2; self$b.p2 <- parms$b.p2
        self$monLwr <- ifelse(!is.null(parms$monLwr), parms$monLwr, 0)
        self$monUpr <- ifelse(!is.null(parms$monUpr), parms$monUpr, ifelse(!is.null(duration), duration, max(parms$upr.u, parms$upr.p1, parms$upr.p2)))
        self$Mdo<-ifelse(!is.null(parms$monUpr) && is.logical(parms$monUpr), parms$monUpr, F)
#        if (!("searchDays" %in% names(parms))) self$searchDays <- NULL else self$searchDays<-parms$searchDays
#        if (!("s0" %in% names(parms) || is.null(parms$s0))) self$s0<- 0
        self$duration<-max(self$upr.u, self$upr.p1, self$upr.p2, self$monUpr)
      } else { # parms are null: assign dummy values for building off of in model builder app
        if (missing(duration)) {
          self$duration <- 365
        } else {
          self$lwr.u<-0; self$upr.u<-duration # default parameters for uniform component
          self$lwr.p1<-round(duration*0.0986); self$upr.p1<-round(duration/2.5) # parameters for beta1 component
          self$lwr.p2<-round(duration*2/3); self$upr.p2<-round(duration*0.88) # parameters for beta2 component
          self$duration<-duration
          self$monLwr <- 0
          self$monUpr <- duration
        }
      }
      self$buildModel()
    },
    buildModel = function() {
      # NOTE: Updates the class variables and Returns a
      # construct a form (with sliders) for building an arrival model
      # build arrivals form to set parameters for compound arrivals
      # minimal error-checking on parms
      for (nm in names(self)) {  # check whether they are the wrong parameters to extract: don't extract environment
        if (is.character(self[[nm]]) || is.numeric(self[[nm]]) || is.logical(self[[nm]]) || is.null(self[[nm]])) {
          assign(nm, self[[nm]]) # extract parameters into function environment
        }
      }
      Udo <-arrcomponents[1]
      P1do<-arrcomponents[2]
      P2do<-arrcomponents[3]
#      duration<-max(Udo*upr.u, P1do*upr.p1, P2do*upr.p2)
      mu.p1<-tclVar(a.p1/(a.p1+b.p1)) # these are values used by TCL in the sliders; must be translated back to R upon saving
      s2.p1<-tclVar(log(12*a.p1*b.p1/((a.p1+b.p1)^2*(a.p1+b.p1+1))))
      mu.p2<-tclVar(a.p2/(a.p2+b.p2))
      s2.p2<-tclVar(log(12*a.p2*b.p2/((a.p2+b.p2)^2*(a.p2+b.p2+1))))
      LU<-tclVar(lwr.u); UU<-tclVar(upr.u) # tcl versions of the parameters
      LP1<-tclVar(lwr.p1); UP1<-tclVar(upr.p1)
      LP2<-tclVar(lwr.p2); UP2<-tclVar(upr.p2)
      WU<-tclVar(wt.u)
      WP1<-tclVar(wt.p1)
      WP2<-tclVar(wt.p2)
      LM<-tclVar(monLwr)
      UM<-tclVar(monUpr)
      .Rvar$arrProcess <- tktoplevel()
      tktitle(.Rvar$arrProcess) <- "Arrival Function Builder"
      tkgrab.set(.Rvar$arrProcess);  tkfocus(.Rvar$arrProcess)
      tkconfigure(.Rvar$arrProcess,width=1000)
      tkwm.resizable(.Rvar$arrProcess,0,0)
      tkwm.deiconify(.Rvar$arrProcess)
      barwidth<-7
      disclr<-colors()[353]
      parmclr<-'white'
      wbord<-1
      scwid<-745
      relief<-"flat"
      xx<-seq(0.005,.995,length=1000)
      # initial values of component parameters

      arrparms<-tclArray()
      columnNames<-c("","Population", "wt", "start", "end", "Ba", "Bb")
      for (i in 1:length(columnNames)) arrparms[[0,i-1]]<-strsplit(columnNames[i]," ",fixed=T)[[1]]
      # unsearched area
      arrparms[[1,0]]<-as.tclObj('',drop=T)
      arrparms[[1,1]]<-"uniform"; arrparms[[1,2]]<-round(wt.u,3) ; arrparms[[1,3]]<-round(lwr.u,3) ; arrparms[[1,4]]<-round(upr.u,3); arrparms[[1,5]]<-1   ; arrparms[[1,6]]<-1
      arrparms[[2,1]]<-"beta1";   arrparms[[2,2]]<-round(wt.p1,3); arrparms[[2,3]]<-round(lwr.p1,3); arrparms[[2,4]]<-round(upr.p1,3);
      arrparms[[2,5]]<-signif(a.p1,4); arrparms[[2,6]]<-signif(b.p1,4)
      arrparms[[3,1]]<-"beta2";   arrparms[[3,2]]<-round(wt.p2,3); arrparms[[3,3]]<-round(lwr.p2,3); arrparms[[3,4]]<-round(upr.p2,3)
      arrparms[[3,5]]<-signif(a.p2,4); arrparms[[3,6]]<-signif(b.p2,4)

      parmTable<-tcltk2::tk2table(.Rvar$arrProcess,
        rows=4,
        cols=length(columnNames),
        selectmode="extended",
        variable=arrparms,
        titlerows="1",
        titlecols="1",
        resizeborders="none",
        multiline=F,
        rowseparator='\n',
        colseparator='\t',
        state='disabled',
        fg='white'
      )
      colwidths<-c(1,10,6,6,8,8,8,7)
      for(i in 1:8) {
        tcl(parmTable, "width", i - 1, colwidths[i])
      }
      tcl(parmTable,"tag","configure", "U",bg=uclr,fg='white')
      tcl(parmTable,"tag","configure", "P1",bg=p1clr,fg='white')
      tcl(parmTable,"tag","configure", "P2",bg=p2clr,fg='white')
      tcl(parmTable,"tag","configure", "hide",bg=disclr, fg=disclr)
      if (Udo) tcl(parmTable,"tag", "rowtag", "U", "1") else tcl(parmTable,"tag", "rowtag", "hide", "1")
      if (P1do) tcl(parmTable,"tag", "rowtag", "P1","2") else tcl(parmTable,"tag", "rowtag", "hide","2")
      if (P2do) tcl(parmTable,"tag", "rowtag", "P2","3") else tcl(parmTable,"tag", "rowtag", "hide","3")
      tcl(parmTable,"tag","configure", "error",bg=colors()[652]) # cells with error have yellow background

    # start of monitoring period is assumed to be t = 0 (although in the database it is s0 in calculations)
    # start of arrival season is entered by user (with default of 0)
    # end of arrival season is entered by user (with default = max(days))
      plotarr <- function() { # not currently scaled to endpoints
        if (!Udo & !P1do & !P2do & !Mdo){
          par(mar=c(4,1.5,.5,.5),mgp=c(2,.7,0))
          plot(0,0,axes=F,xlab='Arrival Time',ylab='',type='n',xlim=c(0,self$duration), yaxs='i')
          mtext(side=1, at = 0, text=0)
          mtext(side=1, at = self$duration, text=self$duration)
          mtext(side=2,line=.7,"Relative Arrival Rate")
          box()
#          if (!is.null(searchDays) && is.numeric(searchDays)){
#            axis(1, at=s0 + searchDays, lab=F, tck=-0.02)
#            axis(1, at=s0 + searchDays, lab=F, tck=0.02)
#          }
          return(F)
        }
        lwr.u<-eoa::toR(LU);   upr.u<-eoa::toR(UU)
        lwr.p1<-eoa::toR(LP1); upr.p1<-eoa::toR(UP1)
        lwr.p2<-eoa::toR(LP2); upr.p2<-eoa::toR(UP2)
        mu.p1 <- as.numeric(tclvalue(mu.p1)); s2.p1<- exp(as.numeric(tclvalue(s2.p1))+log(1/12))
        a.p1 <<- mu.p1^2/s2.p1*(1-mu.p1)-mu.p1; b.p1 <<- a.p1*(1/mu.p1-1)
        # double assignment operator is necessary because changes in values need to be reflected outside the function
        arrparms[[2,5]]<<-signif(a.p1,4); arrparms[[2,6]]<<-signif(b.p1,4)
        mu.p2 <- as.numeric(tclvalue(mu.p2)); s2.p2<- exp(as.numeric(tclvalue(s2.p2))+log(1/12))
        a.p2 <<- mu.p2^2/s2.p2*(1-mu.p2)-mu.p2; b.p2<<-a.p2*(1/mu.p2-1)
        arrparms[[3,5]]<<-signif(a.p2,4); arrparms[[3,6]]<<-signif(b.p2,4)
        wt.u<-eoa::toR(WU)*Udo
        wt.p1<-eoa::toR(WP1)*P1do
        wt.p2<-eoa::toR(WP2)*P2do
        totw <- wt.u + wt.p1 + wt.p2
        wt.u<-wt.u/totw; wt.p1<-wt.p1/totw; wt.p2<-wt.p2/totw
        monLwr<-eoa::toR(LM)
        monUpr<-eoa::toR(UM)
        xx<-seq(0,self$duration,length=1000)
        part1<-numeric(length(xx))
        if (P1do){
          ind<-which(xx>=lwr.p1&xx<=upr.p1)
          part1[ind]<-suppressWarnings(dbeta((xx[ind]-lwr.p1)/(upr.p1-lwr.p1),shape1=a.p1,shape2=b.p1)/(upr.p1-lwr.p1)*wt.p1)
        }
        part2<-numeric(length(xx))
        if (P2do){
          ind<-which(xx>=lwr.p2&xx<=upr.p2)
          part2[ind]<-suppressWarnings(dbeta((xx[ind]-lwr.p2)/(upr.p2-lwr.p2),shape1=a.p2,shape2=b.p2)/(upr.p2-lwr.p2)*wt.p2)
        }
        part3<-numeric(length(xx))
        if (Udo){
          ind<-which(xx>=eoa::toR(lwr.u)&xx<=upr.u)
          part3[ind]<-wt.u/(upr.u-lwr.u)
        }
        arf<-part1 + part2 + part3
        ymax<-max(arf)
        if (is.finite(ymax)){
          par(mar=c(4,1.5,.5,.5),mgp=c(2,.7,0))
          plot(0,0,axes=F,xlab='Arrival Time',ylab='',type='n',xlim=c(1,365),ylim=c(0,1.05)*ymax,yaxs='i')
          mtext(side=1, at = 0, text=0)
          mtext(side=1, at = self$duration, text=self$duration)
          mtext(side=2,line=.7,"Relative Arrival Rate")
#          if (!is.null(searchDays) && is.numeric(searchDays)){
#            axis(1, at=s0 + searchDays, lab=F, tck=-0.02)
#            axis(1, at=s0 + searchDays, lab=F, tck=0.02)
#          }
          if (Mdo) {
            polygon(c(monLwr, monUpr, monUpr, monLwr),par("usr")[c(2,2,3,3)], border=NA, col=colors()[394])
            mtext(side = 1, at = c(monLwr, monUpr), text=c(monLwr, monUpr))
            arrdat<-list(arrcomponents = c(Udo, P1do, P2do),
              lwr.u=lwr.u,    upr.u =upr.u, wt.u =wt.u, # parameters for uniform component
              lwr.p1=lwr.p1,  upr.p1=upr.p1, wt.p1=wt.p1, a.p1=a.p1, b.p1=b.p1, # parameters for beta1 component
              lwr.p2=lwr.p2, upr.p2=upr.p2, wt.p2=wt.p2, a.p2=a.p2, b.p2=b.p2 # parameters for beta2 component
            )
          }
          box()
          if (Udo & lwr.u<upr.u){
            lines(c(lwr.u,upr.u),rep(wt.u/(upr.u-lwr.u),2),col=uclr,lwd=1)
            lines(rep(lwr.u,2),c(0,wt.u/(upr.u-lwr.u)),col=uclr)
            lines(rep(upr.u,2),c(0,wt.u/(upr.u-lwr.u)),col=uclr)
          }
          if (P1do & lwr.p1<upr.p1){
           xx<-seq(lwr.p1,upr.p1,len=250)
           lines(xx,(dbeta((xx-lwr.p1)/(upr.p1-lwr.p1),shape1=a.p1,shape2=b.p1)/(upr.p1-lwr.p1))*wt.p1,col=p1clr,lwd=1)
          }
          if (P2do &lwr.p2<upr.p2){
           xx<-seq(lwr.p2,upr.p2,len=250)
           mu<-(lwr.p2+upr.p2)/2; sd<-(-lwr.p2+upr.p2)/5
           lines(xx,(dbeta((xx-lwr.p2)/(upr.p2-lwr.p2),shape1=a.p2,shape2=b.p2)/(upr.p2-lwr.p2))*wt.p2,col=p2clr,lwd=1)
          }
          if (Mdo){
            # text showing integral of parts
            loc1<-round(integrate(afunc, arrdat=arrdat, lower=0, upper=monLwr)$value, 3)
            loc2<-round(integrate(afunc, arrdat=arrdat, lower=monUpr, upper=self$duration)$value, 3)
            text(par('usr')[1]+diff(par('usr')[1:2])*.03, par('usr')[4]-diff(par('usr')[3:4])*.02, loc1)
            text(par('usr')[2]-diff(par('usr')[1:2])*.03, par('usr')[4]-diff(par('usr')[3:4])*.02, loc2)
            text(mean(c(monLwr, monUpr)), par('usr')[4]-diff(par('usr')[3:4])*.02, lab = paste0("v = ", round(1-loc1-loc2,3)))
          }
          xx<-seq(0.1,364.9,length=1000)
          lines(xx,arf, lwd=2, col=1)
        } else {
          par(mar=c(4,1.5,.5,1.5),mgp=c(2,.7,0))
          plot(0,0,axes=F,xlab='time',ylab='',type='n',xlim=c(0,365),xaxs='i',yaxs='i')
          mtext(side=2,line=.7,"Relative Arrival Rate")
          mtext(side=1, at = 0, text=0)
          mtext(side=1, at = self$duration, text=self$duration)
          box()
#          if (!is.null(searchDays) && is.numeric(searchDays)) {
#            axis(1, at=searchDays, lab=F, tck=-0.02)
#            axis(1,at=searchDays, lab=F, tck=0.02)
#          }
          polygon(c(monLwr, monUpr, monUpr, monLwr),par('usr')[2,2,3,3], colors()[394], border=NA)
          return(F)
        }
      }
      arrfig <- tkrplot::tkrplot(.Rvar$arrProcess, fun = plotarr, hscale=2.2, vscale=1.3)
      onChange <- function(...) {
        tkrplot::tkrreplot(arrfig)
      }
      # radio buttons for defining model form
      modelDefinition<-tkframe(.Rvar$arrProcess) # frame for holding the model definition options
      Mlab<-tklabel(modelDefinition,text="monitoring")
      Mradio<-tclVar(ifelse(Mdo,"yes","no"))
      MyesRadio<-tkradiobutton(modelDefinition,variable=Mradio,value="yes")
      MnoRadio<-tkradiobutton(modelDefinition,variable=Mradio,value="no")
      Ulab<-tklabel(modelDefinition,text="uniform")
      Uradio<-tclVar(ifelse(arrcomponents[1],"yes","no"))
      UyesRadio<-tkradiobutton(modelDefinition,variable=Uradio,value="yes")
      UnoRadio<-tkradiobutton(modelDefinition,variable=Uradio,value="no")
      P1lab<-tklabel(modelDefinition,text="beta1")
      P1radio<-tclVar(ifelse(arrcomponents[2],"yes","no"))
      P1yesRadio<-tkradiobutton(modelDefinition,variable=P1radio,value="yes")
      P1noRadio<-tkradiobutton(modelDefinition,variable=P1radio,value="no")
      P2lab<-tklabel(modelDefinition,text="beta2")
      P2radio<-tclVar(ifelse(arrcomponents[3],"yes","no"))
      P2yesRadio<-tkradiobutton(modelDefinition,variable=P2radio,value="yes")
      P2noRadio<-tkradiobutton(modelDefinition,variable=P2radio,value="no")
      tkgrid(Mlab,MyesRadio,MnoRadio,rowspan=2)
      tkgrid(Ulab,UyesRadio,UnoRadio,rowspan=2)
      tkgrid(P1lab,P1yesRadio,P1noRadio,rowspan=2)
      tkgrid(P2lab,P2yesRadio,P2noRadio,rowspan=2)

      processBounds<-tkframe(.Rvar$arrProcess)
      bndlab<-tklabel(processBounds, text = "Bounds")
      tkgrid(bndlab)
      Ulwr <- tkscale(processBounds, from = 0, to = self$duration, variable = LU, orient = "horizontal",
        length = scwid,
        width=barwidth,
        command = function(...) {
          if (eoa::toR(LU)>=eoa::toR(UU)){
            tclvalue(LU)<-eoa::toR(UU)-1
            arrparms[[1,3]]<<-round(eoa::toR(LU))
            return(F)
          }
          arrparms[[1,3]]<<-round(eoa::toR(LU))
          onChange()
        },
        resolution = 1,
        sliderrelief=relief,
        sliderlength=8,
        troughcolor=uclr,
        showvalue=F,
        borderwidth=wbord
      )
      Uupr <- tkscale(processBounds, from = 0, to = self$duration, variable = UU, orient = "horizontal",
        length = scwid,
        width=barwidth,
        command = function(...) {
          if (eoa::toR(LU)>=eoa::toR(UU)){
            tclvalue(UU)<-eoa::toR(LU)+1
            arrparms[[1,4]]<<-round(eoa::toR(UU))
            return(F)
          }
          arrparms[[1,4]]<<-round(eoa::toR(UU))
          onChange()
        },
        resolution = 1,
        sliderrelief=relief,
        sliderlength=8,
        troughcolor=uclr,
        showvalue=F,
        borderwidth=wbord
      )
      P1lwr <- tkscale(processBounds, from = 0, to = self$duration, variable = LP1, orient = "horizontal",
        length = scwid,
        width=barwidth,
        command = function(...) {
          if (eoa::toR(LP1)>=eoa::toR(UP1)){
            tclvalue(LP1)<-eoa::toR(UP1)-1
            arrparms[[2,3]]<<-round(eoa::toR(LP1))
            return(F)
          }
          arrparms[[2,3]]<<-round(eoa::toR(LP1))
          onChange()
        },
        resolution = 1,
        sliderrelief=relief,
        sliderlength=8,
        troughcolor=p1clr,
        showvalue=F,
        borderwidth=wbord
      )
      P1upr <- tkscale(processBounds, from = 0, to = self$duration, variable = UP1, orient = "horizontal",
        length = scwid,
        width=barwidth,
        command = function(...) {
          if (eoa::toR(LP1)>=eoa::toR(UP1)){
            tclvalue(UP1)<-eoa::toR(LP1)+1
            arrparms[[2,4]]<<-round(eoa::toR(UP1))
            return(F)
          }
          arrparms[[2,4]]<<-round(eoa::toR(UP1))
          onChange()
        },
        resolution = 1,
        sliderrelief=relief,
        sliderlength=8,
        troughcolor=p1clr,
        showvalue=F,
        borderwidth=wbord
      )
      P2lwr <- tkscale(processBounds, from = 0, to = self$duration, variable = LP2, orient = "horizontal",
        length = scwid,
        width=barwidth,
        command = function(...) {
          if (eoa::toR(LP2)>=eoa::toR(UP2)){
            tclvalue(LP2)<-eoa::toR(UP2)-1
            arrparms[[3,3]]<<-round(eoa::toR(LP2))
            return(F)
          }
          arrparms[[3,3]]<<-round(eoa::toR(LP2))
          onChange()
        },
        resolution = 1,
        sliderrelief=relief,
        sliderlength=8,
        troughcolor=p2clr,
        showvalue=F,
        borderwidth=wbord
      )
      P2upr <- tkscale(processBounds, from = 0, to = self$duration, variable = UP2, orient = "horizontal",
        length = scwid,
        width=barwidth,
        command = function(...) {
          if (eoa::toR(LP2)>=eoa::toR(UP2)){
            tclvalue(UP2)<-eoa::toR(LP2)+1
            arrparms[[3,4]]<<-round(eoa::toR(UP2))
            return(F)
          }
          arrparms[[3,4]]<<-round(eoa::toR(UP2))
          onChange()
        },
        resolution = 1,
        sliderrelief=relief,
        sliderlength=8,
        troughcolor=p2clr,
        showvalue=F,
        borderwidth=wbord
      )
      Mlwr<-tkscale(processBounds, from = 0, to = self$duration, variable = LM, orient = "horizontal",
        length = scwid,
        width=barwidth,
        command = function(...) {
          if (eoa::toR(LM)>=eoa::toR(UM)){
            tclvalue(LM)<-eoa::toR(UM)-1
            return(F)
          }
          onChange()
        },
        resolution = 1,
        sliderrelief=relief,
        sliderlength=8,
        troughcolor=monclr,
        showvalue=F,
        borderwidth=wbord
      )
      Mupr<-tkscale(processBounds, from = 0, to = self$duration, variable = UM, orient = "horizontal",
        length = scwid,
        width=barwidth,
        command = function(...) {
          if (eoa::toR(LM)>=eoa::toR(UM)){
            tclvalue(UM)<-eoa::toR(LM)+1
            return(F)
          }
          onChange()
        },
        resolution = 1,
        sliderrelief=relief,
        sliderlength=8,
        troughcolor=monclr,
        showvalue=F,
        borderwidth=wbord
      )
      tkgrid(Mlwr)
      tkgrid(Mupr)
      tkgrid(Ulwr)
      tkgrid(Uupr)
      tkgrid(P1lwr)
      tkgrid(P1upr)
      tkgrid(P2lwr)
      tkgrid(P2upr)

      ## parameter summary for bounds
      processWts<-tkframe(.Rvar$arrProcess)
      wtslids<-tkframe(processWts)
      scht<-400
      setWts<-function(){
        ww<-eoa::toR(WU)
        wo<-eoa::toR(WP1)*P1do+eoa::toR(WP2)*P2do
        arrparms[[1,2]]<<-round(ifelse(abs(wo)<0.000001,1,ww/(wo+ww)),2)
        ww<-eoa::toR(WP1)
        wo<-eoa::toR(WU)*Udo+eoa::toR(WP2)*P2do
        arrparms[[2,2]]<<-round(ifelse(abs(wo)<0.000001,1,ww/(wo+ww)),2)
        ww<-eoa::toR(WP2)
        wo<-eoa::toR(WU)*Udo+eoa::toR(WP1)*P1do
        arrparms[[3,2]]<<-round(ifelse(abs(wo)<0.000001,1,ww/(wo+ww)),2)
      }
      Uwt<-tkscale(wtslids, from = 1, to = 0, variable = WU, orient = "vertical",
        length = scht,
        width=barwidth,
        command = function(...) {
          setWts()
          onChange()
        },
        resolution = 0.01,
        sliderrelief=relief,
        sliderlength=8,
        troughcolor=uclr,
        showvalue=F,
        borderwidth=wbord
      )
      P1wt<-tkscale(wtslids, from = 1, to = 0, variable = WP1, orient = "vertical",
        length = scht,
        width=barwidth,
        command = function(...){
          setWts()
          onChange()
        },
        resolution = 0.01,
        sliderrelief=relief,
        sliderlength=8,
        troughcolor=p1clr,
        showvalue=F,
        borderwidth=wbord
      )
      P2wt<-tkscale(wtslids, from = 1, to = 0, variable = WP2, orient = "vertical",
        length = scht,
        width=barwidth,
        command = function(...){
          setWts()
          onChange()
        },
        resolution = 0.01,
        sliderrelief=relief,
        sliderlength=8,
        troughcolor=p2clr,
        showvalue=F,
        borderwidth=wbord
      )
      tkgrid(Uwt,P1wt,P2wt,pady=10, sticky='e')
      wtlab<-tklabel(processWts, text = "Weights")
      tkgrid(wtslids)
      tkgrid(wtlab)

      processMeans<-tkframe(.Rvar$arrProcess)
      mudiff<-0.015 # how close does the mean have to be in order to trigger a "sticky" uniform distribution for beta (==> a=1, b=1)
      s2diff<-1.05 # how close does variance have to be to initiate sticky uniform beta?
      P1mu<-tkscale(processMeans, from = 0, to = 1, variable = mu.p1, orient = "horizontal",
        length = scwid,
        width=barwidth,
        command = function(...){
          mu<-as.numeric(tclvalue(mu.p1))
          s2<-exp(as.numeric(tclvalue(s2.p1))+log(1/12))
          if (mu < 0.5 + mudiff & mu > 0.5-mudiff) {
            mu<-0.5
          }
          mu<-max(mu, (1-sqrt(1-4*s2))/2+0.005)
          mu<-min(mu, (1+sqrt(1-4*s2))/2-0.005)
          tclvalue(mu.p1)<-mu
          onChange()
        },
        resolution=0.005,
        sliderrelief=relief,
        sliderlength=8,
        troughcolor=p1clr,
        showvalue=F,
        borderwidth=wbord
      )
      P2mu<-tkscale(processMeans, from = 0, to = 1, variable = mu.p2, orient = "horizontal",
        length = scwid,
        width=barwidth,
        command = function(...){
          mu<-as.numeric(tclvalue(mu.p2))
          s2<-exp(as.numeric(tclvalue(s2.p2))+log(1/12))
          if (mu < 0.5 + mudiff & mu > 0.5-mudiff) {
            mu<-0.5
          }
          mu<-max(mu, (1-sqrt(1-4*s2))/2+0.005)
          mu<-min(mu, (1+sqrt(1-4*s2))/2-0.005)
          tclvalue(mu.p2)<-mu
          onChange()
        },
        resolution = 0.005,
        sliderrelief=relief,
        sliderlength=8,
        troughcolor=p2clr,
        showvalue=F,
        borderwidth=wbord
      )
      tkgrid(P1mu)
      tkgrid(P2mu)
      mulab<-tklabel(processMeans, text = "Symmetry")
      tkgrid(mulab)

      processVariances<-tkframe(.Rvar$arrProcess)
      vslids<-tkframe(processVariances)
      P1s2<-tkscale(vslids, from = -6, to = 1.05, variable = s2.p1, orient = "vertical",
        length = scht,
        width=barwidth,
        command = function(...){
          if (abs(as.numeric(tclvalue(s2.p1)))<.09) tclvalue(s2.p1)<-0
          s2<-exp(as.numeric(tclvalue(s2.p1))+log(1/12))
          mu<-as.numeric(tclvalue(mu.p1))
          s2<-min(s2, mu*(1-mu)*.99)
          tclvalue(s2.p1)<-log(12*s2)
          onChange()
        },
        resolution = 0.05,
        sliderrelief=relief,
        sliderlength=8,
        troughcolor=p1clr,
        showvalue=F,
        borderwidth=wbord
      )
      P2s2<-tkscale(vslids, from = -6, to = 1.05, variable = s2.p2, orient = "vertical",
        length = scht,
        width=barwidth,
        command = function(...){
          if (abs(as.numeric(tclvalue(s2.p2)))<.09) tclvalue(s2.p2)<-0
          s2<-exp(as.numeric(tclvalue(s2.p2))+log(1/12))
          mu<-as.numeric(tclvalue(mu.p2))
          s2<-min(s2, mu*(1-mu)*.99)
          tclvalue(s2.p2)<-log(12*s2)
          onChange()
        },
        resolution = 0.1,
        sliderrelief=relief,
        sliderlength=8,
        troughcolor=p2clr,
        showvalue=F,
        borderwidth=wbord
      )
      vlab<-tklabel(processVariances, text = "Aspect")
      tkgrid(P1s2,P2s2,pady=10)
      tkgrid(vslids)
      tkgrid(vlab)
      tkgrid(processVariances)
      arrButtFrame<-tkframe(.Rvar$arrProcess)
      arrSave<-tkbutton(arrButtFrame,text="Save", width=15,command=function(){
        # compile compound parameters
        self$arrcomponents<-as.logical(c(Udo,P1do, P2do))
        self$Mdo <- as.logical(Mdo)
        self$wt.u  <- eoa::toR(arrparms[[1,2]])
        self$lwr.u <- eoa::toR(arrparms[[1,3]])
        self$upr.u <- eoa::toR(arrparms[[1,4]])
        self$wt.p1 <- eoa::toR(arrparms[[2,2]])
        self$lwr.p1<- eoa::toR(arrparms[[2,3]])
        self$upr.p1<- eoa::toR(arrparms[[2,4]])
        self$wt.p2 <- eoa::toR(arrparms[[3,2]])
        self$lwr.p2<- eoa::toR(arrparms[[3,3]])
        self$upr.p2<- eoa::toR(arrparms[[3,4]])
        self$monLwr<- eoa::toR(LM)
        self$monUpr<- eoa::toR(UM)
        self$arrfun<-function(tt){
          part1<-numeric(length(tt))
          Udo<-force(self$arrcomponents[1])
          P1do<-force(self$arrcomponents[2])
          P2do<-force(self$arrcomponents[3])
          lwr.u<-force(self$lwr.u); upr.u<-force(self$upr.u); wt.u<-force(self$wt.u)
          lwr.p1<-force(self$lwr.p1); upr.p1<-force(self$upr.p1); wt.p1<-force(self$wt.p1)
          lwr.p2<-force(self$lwr.p2); upr.p2<-force(self$upr.p2); wt.p2<-force(self$wt.p2)
          a.p1<-force(self$a.p1); b.p1<-force(self$b.p1)
          a.p2<-force(self$a.p2); b.p2<-force(self$b.p2)
          if (P1do){
            ind<-which(tt>=lwr.p1&tt<=upr.p1)
            part1[ind]<-suppressWarnings(dbeta((tt[ind]-lwr.p1)/(upr.p1-lwr.p1),shape1=a.p1,shape2=b.p1)/(upr.p1-lwr.p1)*wt.p1)
          }
          part2<-numeric(length(tt))
          if (P2do){
            ind<-which(tt>=lwr.p2&tt<=upr.p2)
            part2[ind]<-suppressWarnings(dbeta((tt[ind]-lwr.p2)/(upr.p2-lwr.p2),shape1=a.p2,shape2=b.p2)/(upr.p2-lwr.p2)*wt.p2)
          }
          part3<-numeric(length(tt))
          if (Udo){
            ind<-which(tt>=lwr.u&tt<=upr.u)
            part3[ind]<-wt.u/(upr.u-lwr.u)
          }
          part1 + part2 + part3
        }
        tkdestroy(.Rvar$arrProcess)
      })
      arrCancel<-tkbutton(arrButtFrame,text='Cancel', width=15,command=function() tkdestroy(.Rvar$arrProcess))
      tkgrid(arrSave)
      tkgrid(arrCancel)

      tkbind(MyesRadio,"<Button-1>", function() {
        arrcomponents[1]<-1
        tkconfigure(Mlwr,state='normal',troughcolor=monclr)
        tkconfigure(Mupr,state='normal',troughcolor=monclr)
        # redraw the figure with the uniform added back in
        Mdo<<-T
        tkrplot::tkrreplot(arrfig)
      })
      tkbind(MnoRadio,"<Button-1>", function() {
    #    if (!dateset) return(F)
        tkconfigure(Mlwr,state='disabled',troughcolor=disclr)
        tkconfigure(Mupr,state='disabled',troughcolor=disclr)
        Mdo<<-F
        tkrplot::tkrreplot(arrfig)
      })

      tkbind(UyesRadio,"<Button-1>", function() {
        arrcomponents[1]<-1
        tkconfigure(Uwt,state='normal',troughcolor=uclr)
        tkconfigure(Ulwr,state='normal',troughcolor=uclr)
        tkconfigure(Uupr,state='normal',troughcolor=uclr)
        tcl(parmTable,"tag", "rowtag", "U","1")
        # redraw the figure with the uniform added back in
        Udo<<-T
        setWts()
        tkrplot::tkrreplot(arrfig)
      })
      tkbind(UnoRadio,"<Button-1>", function() {
    #    if (!dateset) return(F)
        arrcomponents[1]<<-0
        tkconfigure(Uwt,state='disabled',troughcolor=disclr)
        tkconfigure(Ulwr,state='disabled',troughcolor=disclr)
        tkconfigure(Uupr,state='disabled',troughcolor=disclr)
        tcl(parmTable,"tag", "rowtag", "hide","1")
        # redraw the figure with the uniform taken out
        Udo<<-F
        setWts()
        tkrplot::tkrreplot(arrfig)
      })
      tkbind(P1yesRadio,"<Button-1>", function() {
    #    if (!dateset) return(F)
        arrcomponents[2]<<-1
        tkconfigure(P1wt,state='normal',troughcolor=p1clr)
        tkconfigure(P1lwr,state='normal',troughcolor=p1clr)
        tkconfigure(P1upr,state='normal',troughcolor=p1clr)
        tkconfigure(P1mu,state='normal',troughcolor=p1clr)
        tkconfigure(P1s2,state='normal',troughcolor=p1clr)
        tcl(parmTable,"tag", "rowtag", "P1","2")
        # redraw the figure with the uniform added back in
        P1do<<-T
        setWts()
        tkrplot::tkrreplot(arrfig)
      })
      tkbind(P1noRadio,"<Button-1>", function() {
    #    if (!dateset) return(F)
        arrcomponents[2]<<-0
        tkconfigure(P1wt,state='disabled',troughcolor=disclr)
        tkconfigure(P1lwr,state='disabled',troughcolor=disclr)
        tkconfigure(P1upr,state='disabled',troughcolor=disclr)
        tkconfigure(P1mu,state='disabled',troughcolor=disclr)
        tkconfigure(P1s2,state='disabled',troughcolor=disclr)
        tcl(parmTable,"tag", "rowtag", "hide","2")
        # redraw the figure with the uniform taken out
        P1do<<-F
        setWts()
        tkrplot::tkrreplot(arrfig)
      })
      tkbind(P2yesRadio,"<Button-1>", function() {
    #    if (!dateset) return(F)
        arrcomponents[3]<<-1
        tkconfigure(P2wt,state='normal',troughcolor=p2clr)
        tkconfigure(P2lwr,state='normal',troughcolor=p2clr)
        tkconfigure(P2upr,state='normal',troughcolor=p2clr)
        tkconfigure(P2mu,state='normal',troughcolor=p2clr)
        tkconfigure(P2s2,state='normal',troughcolor=p2clr)
        tcl(parmTable,"tag", "rowtag", "P2","3")
        # redraw the figure with the uniform added back in
        P2do<<-T
        setWts()
        tkrplot::tkrreplot(arrfig)
      })
      tkbind(P2noRadio,"<Button-1>", function() {
    #    if (!dateset) return(F)
        arrcomponents[3]<<-0
        tkconfigure(P2wt,state='disabled',troughcolor=disclr)
        tkconfigure(P2lwr,state='disabled',troughcolor=disclr)
        tkconfigure(P2upr,state='disabled',troughcolor=disclr)
        tkconfigure(P2mu,state='disabled',troughcolor=disclr)
        tkconfigure(P2s2,state='disabled',troughcolor=disclr)
        tcl(parmTable,"tag", "rowtag", "hide","3")
        # redraw the figure with the uniform taken out
        P2do<<-F
        setWts()
        tkrplot::tkrreplot(arrfig)
      })
      ##############
      if (Mdo){
        Mstate<-'normal'
        clr<-monclr
      } else {
        Mstate<-'disabled'
        clr<-disclr
      }
      tkconfigure(Mlwr,state=Mstate,troughcolor=clr)
      tkconfigure(Mupr,state=Mstate,troughcolor=clr)
      if (Udo){
        Ustate<-'normal'
        clr<-uclr
      } else {
        Ustate<-'disabled'
        clr<-disclr
      }
      tkconfigure(Uwt,state=Ustate,troughcolor=clr)
      tkconfigure(Ulwr,state=Ustate,troughcolor=clr)
      tkconfigure(Uupr,state=Ustate,troughcolor=clr)
      if (P1do){
        P1state<-'normal'
        clr<-p1clr
      } else {
        P1state<-'disabled'
        clr<-disclr
      }
      tkconfigure(P1wt,state=P1state,troughcolor=clr)
      tkconfigure(P1lwr,state=P1state,troughcolor=clr)
      tkconfigure(P1upr,state=P1state,troughcolor=clr)
      tkconfigure(P1mu,state=P1state,troughcolor=clr)
      tkconfigure(P1s2,state=P1state,troughcolor=clr)
      if (P2do){
        P2state<-'normal'
        clr<-p2clr
      } else {
        P2state<-'disabled'
        clr<-disclr
      }
      tkconfigure(P2wt,state=P2state,troughcolor=clr)
      tkconfigure(P2lwr,state=P2state,troughcolor=clr)
      tkconfigure(P2upr,state=P2state,troughcolor=clr)
      tkconfigure(P2mu,state=P2state,troughcolor=clr)
      tkconfigure(P2s2,state=P2state,troughcolor=clr)

      # labels for sliders
      tkgrid(parmTable, row = 0, column = 0, columnspan=2,sticky='sw',pady=10,padx=10)
      tkgrid(arrButtFrame, column=1, row=0, sticky='e')
      tkgrid(modelDefinition, row = 1, column = 0, sticky='s')
      tkgrid(processBounds, column = 1, row = 1, sticky='n')
#      arrlab<-tkrplot::tkrplot(.Rvar$arrProcess, function() plot(0,0), hscale=2.2, vscale=1.3)
#      tkgrid(arrlab, row=2, column =1 ,sticky='ne')
      tkgrid(processWts, arrfig, row = 2, sticky='ne',padx=10)
      tkgrid(processVariances, column=2, row=2, sticky = 'nw', padx=c(0,5))
      tkgrid(processMeans, row = 3, column = 1, sticky='n')
      blankframe<-tklabel(.Rvar$arrProcess,text=as.tclObj('        ',drop=T))
      tkgrid(blankframe)

    }
  )
)
afunc<-function(tt, arrdat){
  # arrdat is a list with model parameters
  part1<-part2<-part3<-numeric(length(tt))
  with(arrdat, {
    if (arrcomponents[1]){
      ind<-which(tt>=lwr.u&tt<=upr.u)
      part1[ind]<<-wt.u/(upr.u-lwr.u)
    }
    if (arrcomponents[2]){
      ind<-which(tt>=lwr.p1&tt<=upr.p1)
      part2[ind]<<-suppressWarnings(dbeta((tt[ind]-lwr.p1)/(upr.p1-lwr.p1),shape1=a.p1,shape2=b.p1)/(upr.p1-lwr.p1)*wt.p1)
    }
    if (arrcomponents[3]){
      ind<-which(tt>=lwr.p2&tt<=upr.p2)
      part3[ind]<<-suppressWarnings(dbeta((tt[ind]-lwr.p2)/(upr.p2-lwr.p2),shape1=a.p2,shape2=b.p2)/(upr.p2-lwr.p2)*wt.p2)
    }
  })
  part1 + part2 + part3
}
# input format for calling the estg function...list
calcg.fixed <- function(gdat, arrdat = NULL){
###############
### arrdat specifies the arrival parameters to use
## arrdat = NULL   => uniform arrivals over the monitored period
## arrdat = vector => fraction of carcasses arriving in each search interval
## arrdat = list() => arrival season longer than monitored period with:
#   arrfun = vector, 'Uniform', or function defined (at least) on [0, duration]
#   additional arguments needed if arrfun is a vector (ignored otherwise):
#     arrmiss0 = fraction of carcasses arriving before monitoring begins
#     arrmissf = fraction of carcasses arriving after monitoring ends
#   additional arguments needed if arrfun is a function or 'Uniform' (ignored otherwise):
#     s0 = beginning of monitoring (relative to beginning of arrivals)
#     duration = length of monitoring season
#     arrSimplify = T or F to indicate whether or not to assume uniform arrivals within search intervals
###############
  
  clipProb<-0.001
  for (nm in names(gdat)) assign(nm, gdat[[nm]])

  ### preliminaries
  # search schedule
  if (samtype == "Formula") {
    days<-c(0, 1:nsearch)*Isam
  } else if (samtype == "Custom"){
    Isam <- round(max(days)/(length(days)-1),1)
    nsearch <- length(days)-1
  } else {
    warning("error in search schedule. aborting calculation.")
    return(F)
  }
  ind1<-rep(1:nsearch,times=nsearch:1)
  ind2<-ind1+1
  ind3<-unlist(lapply(1:nsearch, function(x) x:nsearch))+1
  schedule.index<-cbind(ind1,ind2,ind3)
  schedule<-cbind(days[ind1],days[ind2],days[ind3])
  nmiss<-schedule.index[,3]-schedule.index[,2]
  maxmiss<-max(nmiss)
  pobs<-numeric(dim(schedule)[1])
  # searcher efficiency
  powk<-cumprod(c(1,rep(k,maxmiss))) # vector of k^i's
  notfind<-cumprod(1-p*powk[-length(powk)])
  nvec<-c(1,notfind)*p
  pfind.si<-nvec*powk # conditional probability of finding a carcass on the ith search (row) after arrival for given (simulated) searcher efficiency (column)
  # persistences
  intxsearch<-unique(cbind(schedule[,2] - schedule[,1], schedule[,3] - schedule[,2]), MAR=1)
  if (persistence_distn == "Exponential") pda<-1/pdb
  ppersu<-eoa::ppersist(persistence_distn, t_arrive0 = 0, t_arrive1=intxsearch[,1], t_search = intxsearch[,1]+intxsearch[,2], pda = pda, pdb = pdb)
  # arrivals
  if (is.null(arrdat)) { # inference is for monitored period only, with uniform arrivals
    arrvec<-(schedule[,2]-schedule[,1])/max(days)
    arrmiss0<-0
    arrmissf<-0
    v<-1
    arrSimplify <-T
  } else if (is.numeric(arrdat)) {
    if (length(arrfun) != nsearch){
      warning(paste0('Vector of arrival probabilities is ', ifelse(length(arrfun)> nsearch, 'longer ', 'shorter '), 'than number of searches. Aborting calculation.'))
      return(F)
    }
    if (sum(is.na(arrfun)) > 0){
      warning('Missing values not allowed in vector of arrival probabilities (arrfun).')
      return(F)
    }
    arrvec<-rep(arrfun, nsearch:1)
    arrmiss0<-0
    arrmissf<-0
    v<-1
    arrSimplify<-T
  } else if (is.list(arrdat) || inherits(arrdat, "R6")){
    if (is.function(arrdat$arrfun)){
      if (arrdat$duration < arrdat$s0 + max(days)){
        warning("Error in arrdat: duration must be at least s0 + max(days). Period of inference shorter than monitored period. Aborting calculation.")
        return(F)
      }
      arrSimplify<-ifelse(is.null(arrdat$arrSimplify), F, arrdat$arrSimplify)
      deno<-integrate(arrfun, lower = 0, upper = arrdat$duration)$val
      arrmiss0<-integrate(arrfun, lower = 0, upper = arrdat$s0)$val/deno
      arrmissf<-integrate(arrfun, lower = arrdat$s0 + max(days), upper = arrdat$duration)$val/deno
      v<-1-arrmiss0-arrmissf
      if (arrSimplify){
        arrvec0<-numeric(nsearch)
        for (i in 1:nsearch) {
          arrvec0[i] <- integrate(arrfun, lower = days[i], upper = days[i+1])$val
        }
        arrvec0<-arrvec0/sum(arrvec0) # normalize the arrival vector
        arrvec<-rep(arrvec0,nsearch:1)
      }
    } else if (is.numeric(arrdat$arrfun) || arrdat$arrfun == "Uniform"){
      arrmiss0<-arrdat$arrmiss0
      arrmissf<-arrdat$arrmissf
      v<-1-arrmiss0-arrmissf
      arrSimplify <- T
      if (arrdat$arrfun=="Uniform"){
        arrvec<-(schedule[,2]-schedule[,1])/max(days)
      } else {
        arrvec<-rep(arrfun, nsearch:1)
      }
    }
  } else {
    warning("error in arrdat specification. Aborting calculation")
    return(F)
  }
  # detection probabilities
  if (arrSimplify){
    for (i in 1:length(pobs)){
      pobs[i]<-pfind.si[nmiss[i]+1] *
        ppersu[which(
          abs(intxsearch[,1]-(schedule[i,2]-schedule[i,1]))<clipProb &
          abs(intxsearch[,2]-(schedule[i,3]-schedule[i,2]))<clipProb),] *
        arrvec[i]
    }
    pobs<-sum(pobs)
  } else { # must integrate for each interval
     # calculate probability arrives when it does and persists until given search date
      arrdeno<-integrate(arrfun, lower = s0, upper = sf)$val
      ipart<-numeric(dim(schedule)[1])
      for (i in 1:dim(schedule)[1]){
        ipart[i]<-switch(persistence_distn,
          "Exponential" = integrate(function(tt)
            (1-pexp(schedule[i,3]-tt, rate=1/pdb)) *
            arrfun(tt + s0)/arrdeno, lower=schedule[i,1], upper=schedule[i,2])$val,
          "Weibull" =  integrate(function(tt)
            (1-pweibull(schedule[i,3]-tt, shape=pda, scale=pdb)) *
            arrfun(tt + s0)/arrdeno, lower=schedule[i,1], upper=schedule[i,2])$val,
          "Lognormal" = integrate(function(tt)
            (1-plnorm(schedule[i,3]-tt, meanlog=pdb, sdlog=sqrt(pda))) *
            arrfun(tt + s0)/arrdeno, lower=schedule[i,1], upper=schedule[i,2])$val,
          "Log-Logistic" = integrate(function(tt)
            (1-actuar::pllogis(schedule[i,3]-tt, shape=pda, scale=pdb)) *
            arrfun(tt + s0)/arrdeno, lower=schedule[i,1], upper=schedule[i,2])$val,
          NA)
      }
      pfind<-pfind.si[nmiss+1]
      pobs<-sum(pfind*ipart)
   }
  return(list('Full site, full year' = pobs*a*v, 'Full site, monitored period' = pobs*a, 'Searched area, monitored period' = pobs)) # detection probability for full year and site
}



### calculate a CI for a posterior lambda
LCI<-function(fun, crlev = 0.95){
  # fun is assumed to be a CDF (posterior distribution) with support x >= 0...no error-checking
  # find the lower bound
  LL <- 0.001
  if (fun(LL) > (1 - crlev)/2){
    while (fun(LL) > (1 - crlev)/2) LL <- LL/2
    low<-uniroot(function(L) fun(L)-(1-crlev)/2, lower = LL, upper = 2 * LL)$root
  } else {
    while (fun(LL) < (1 - crlev)/2) LL <- 2*LL
    low<-uniroot(function(L) fun(L)-(1-crlev)/2, lower = LL/2, upper = LL)$root
  }
  LU <- LL
  while (fun(LU) < 1-(1-crlev)/2) LU <- 2 * LU
  upp<-uniroot(function(L) fun(L)-(1-(1-crlev)/2), lower = LU/2, upper = LU)$root
  c(low, upp)
}
