compcolor<-colors()[405]
errcolor<-colors()[652]
uclr<-colors()[448]
p1clr<-colors()[123]
p2clr<-colors()[96]
lwrclr<-colors()[519]
uprclr<-colors()[652]
monclr<-colors()[94]
ats0 <- c(0,31,59,90,120,151,181,212,243,273,304,334)   # days since jan01
figlab<-format(as.Date(ats0,origin="1970-1-1"),"%b-%d")
listboxbgclr<-'white'
expoColor<-colors()[370]
weibullColor<-colors()[148]
llogisColor<-colors()[69]
lnormColor<-colors()[125]
basicMode<-T
prnames<-c("Objective", "Custom")
pdnames<-c("Exponential","Weibull","Log-Logistic","Lognormal")
singleYearDefault<-list(X = 2, a = .4, v = 0.75, SEn = 158, SEx = 98, k = .674, SEopt = 'h', # 'h' -> enter k by hand; 'f' -> estimate k via field trials
  samtype="Formula", # "Formula" or "Custom"
  Isam = 7, nsearch = 26, days = seq(0, 182, by = 7), firstsearch = "2016-03-15",
  perstype = 'hand', # persistence parameters entered by hand ('hand') or derived from field trial data ('field')
  persistence_distn = "Lognormal" , pda=4.0827, pdb = 1.1707, blwr = 0.4871, bupr= 1.854, #Ir = 7,
  pda.w =0.57792, pdb.w  = 7.1683, blwr.w = 3.919, bupr.w = 13.11,
  pda.ll=0.83052, pdb.ll = 3.3816, blwr.ll = 1.705, bupr.ll= 6.708,
  pda.ln=4.0827, pdb.ln = 1.1707, blwr.ln = 0.4871, bupr.ln= 1.854,
  pda.e=0.11594, pdb.e = 8.6251, blwr.e = 6.123, bupr.e= 12.15,
  prior_f='Objective',
  prior_M=diff(sqrt(0:101))/sum(diff(sqrt(0:101))),
  custom.prior = NA,#c(0.000906, 0.002425, 0.004328, 0.006443, 0.008638, 0.010816, 0.012909, 0.014867, 0.016658, 0.018263, 0.019672, 0.020883, 0.0219, 0.022729, 0.02338, 0.023865, 0.024198, 0.02439, 0.024455, 0.024408, 0.024259, 0.024023, 0.023709, 0.023329, 0.022892, 0.022407, 0.021884, 0.021328, 0.020748, 0.020149, 0.019536, 0.018914, 0.018287, 0.01766, 0.017035, 0.016414, 0.015802, 0.015198, 0.014606, 0.014025, 0.013459, 0.012907, 0.01237, 0.011849, 0.011344, 0.010856, 0.010384, 0.009929, 0.00949, 0.009067, 0.00866, 0.00827, 0.007894, 0.007534, 0.007189, 0.006859, 0.006542, 0.006239, 0.005949, 0.005672, 0.005407, 0.005154, 0.004912, 0.004681, 0.004461, 0.004251, 0.004051, 0.003859, 0.003677, 0.003503, 0.003338, 0.003179, 0.003029, 0.002885, 0.002748, 0.002616, 0.002491, 0.002371, 0.002256, 0.002146, 0.002041, 0.001939, 0.001841, 0.001746, 0.001654, 0.001565, 0.001478, 0.001394, 0.001312, 0.001232, 0.001153, 0.001077, 0.001003, 0.00093, 0.00086, 0.000791, 0.000725, 0.000662, 0.000601, 0.000544, 0.000489, 0.000437, 0.000389, 0.000344, 0.000302, 0.000264, 0.000229, 0.000198, 0.00017, 0.000145, 0.000122, 0.000103, 8.59E-05, 7.13E-05, 5.87E-05, 4.81E-05, 3.91E-05, 3.15E-05, 2.53E-05, 2.01E-05, 1.59E-05, 1.25E-05, 9.73E-06, 7.53E-06, 5.79E-06, 4.42E-06, 3.35E-06, 2.52E-06, 1.89E-06, 1.4E-06),
  objective.prior = diff(sqrt(0:99))/sum(diff(sqrt(0:99))),
  arrfun='Uniform', arrcomponents = c(T, T, T), arrstart = 0,
  lwr.u=5, upr.u=360, wt.u=1/4,
  lwr.p1=30, upr.p1=150, wt.p1=1/4, a.p1=2.5, b.p1=2.5,
  lwr.p2=200, upr.p2=300, wt.p2=1/2, a.p2=3.5, b.p2=4.5,
#  atemporal = 0.722, arrmiss0 = 0.206, arrmissf = 0.072, # variables that are not user-entered should not be kept in parameter set
#  Bab = c(8.8269,19.4668), BabAnn = c(10.2332,35.3502),
#  syA=0.1
  crlev = 0.9
)
dtDefault<-list(Ifix = 1, ffix = 0, phifix = 0,
  target_g = 0.25, grinc = 120, tarr = 1,
  Isam = 3, Imin = 3, Imax = 28, span = 168, firstsearch = "2016-03-15",
  phi = .85, phimin = 0.1, phimax = 1,
  f = .25, fmin = 0.05, fmax = 1,
  k=.67,
  persistence_distn = "Weibull" , pda0 = 0.7209, pdb0 = 8.118,
  pda.w=0.7209, pdb.w = 8.118,
  pda.ll=2.37246, pdb.ll = 7.32333,
  pda.ln=1.09861, pdb.ln = 1.75328,
  pda.e=1/10, pdb.e = 10,
  arrfun='Uniform', arrcomponents = c(T, T, T), arrstart = 0,
  lwr.u=5, upr.u=360, wt.u=1/4,
  lwr.p1=30, upr.p1=150, wt.p1=1/4, a.p1=2.5, b.p1=2.5,
  lwr.p2=200, upr.p2=300, wt.p2=1/2, a.p2=3.5, b.p2=4.5,
  tau = 5, crlev = 0.9, X = 0, dfoption = 'g'
)
symcDefault<-list(
  class=c("Unsearched", "Easy", "Moderate", "Difficult"),
  Xmc=c(0, 1, 0, 0),
  Ba=c(0.00001, 3.05, 4.5, 3.5),
  Bb=c(1000, 1.64, 5.5, 31.5),
  rel_wt=c(0.3, 0.15, 0.35, 0.2),
  option="M", # calculate overall detection probability only (g) or estimate total mortality (M) [the latter requires a prior and an alpha]
  ICEoption="h", # enter data for classes via serial estimation from monitoring data (e) or by hand (h)
  crlev=0.8,
  firstsearch = "2016-03-15", span = 180,
  prior_f='Objective',
  prior_M=diff(sqrt(0:63))/sum(diff(sqrt(0:63))),
  custom.prior = NA,#c(0.000906, 0.002425, 0.004328, 0.006443, 0.008638, 0.010816, 0.012909, 0.014867, 0.016658, 0.018263, 0.019672, 0.020883, 0.0219, 0.022729, 0.02338, 0.023865, 0.024198, 0.02439, 0.024455, 0.024408, 0.024259, 0.024023, 0.023709, 0.023329, 0.022892, 0.022407, 0.021884, 0.021328, 0.020748, 0.020149, 0.019536, 0.018914, 0.018287, 0.01766, 0.017035, 0.016414, 0.015802, 0.015198, 0.014606, 0.014025, 0.013459, 0.012907, 0.01237, 0.011849, 0.011344, 0.010856, 0.010384, 0.009929, 0.00949, 0.009067, 0.00866, 0.00827, 0.007894, 0.007534, 0.007189, 0.006859, 0.006542, 0.006239, 0.005949, 0.005672, 0.005407, 0.005154, 0.004912, 0.004681, 0.004461, 0.004251, 0.004051, 0.003859, 0.003677, 0.003503, 0.003338, 0.003179, 0.003029, 0.002885, 0.002748, 0.002616, 0.002491, 0.002371, 0.002256, 0.002146, 0.002041, 0.001939, 0.001841, 0.001746, 0.001654, 0.001565, 0.001478, 0.001394, 0.001312, 0.001232, 0.001153, 0.001077, 0.001003, 0.00093, 0.00086, 0.000791, 0.000725, 0.000662, 0.000601, 0.000544, 0.000489, 0.000437, 0.000389, 0.000344, 0.000302, 0.000264, 0.000229, 0.000198, 0.00017, 0.000145, 0.000122, 0.000103, 8.59E-05, 7.13E-05, 5.87E-05, 4.81E-05, 3.91E-05, 3.15E-05, 2.53E-05, 2.01E-05, 1.59E-05, 1.25E-05, 9.73E-06, 7.53E-06, 5.79E-06, 4.42E-06, 3.35E-06, 2.52E-06, 1.89E-06, 1.4E-06),
  objective.prior = diff(sqrt(0:63))/sum(diff(sqrt(0:63))),
  arrstart=0, arr0u=5, arr1u=360,
  arrcomponents=c(T,T,T),
  lwr.u = 5, upr.u = 360, wt.u = 0.25,
  lwr.p1 = 30, upr.p1 = 150, wt.p1 = 0.25, a.p1 = 2.5, b.p1 = 2.5,
  lwr.p2 = 200, upr.p2 = 300, wt.p2 = 0.5, a.p2 = 3.5, b.p2 = 4.5,
  arrfun='Uniform'
)
syscDefault<-singleYearDefault # syscPrevious is the data from the previous sy class when doing multiple classes
myDefault<-list(
  years=1942:1946,
  X=c(2,0,1,0,0),
  Ba=c(100.50,100.50,100.50,235.44, 235.44),
  Bb=c(234.50, 234.50, 234.50, 2707.56, 2707.56),
  rel_wt=rep(1,5),
  option ="M", # M, L: what type of test to run? mortality (if "M" then further options for tracking or projection), short-term rate, reversion
  Mtype="C", # C, T, P: what type of total to calculate? Cumultative total (C), tracking (T), projection (P)
  Ltype="stT", # stT, rT, ciT: what type of lambda test to run? short-term test, reversion test, confidence interval
  Ptype="I", # I, C, V: what type of monitoring and operations in future year? I = same as most recent, C = constant but not same as last, V = varies by year
  crlev = 0.5, nyr = 30, Tau = 60, # alpha value for estimating mortality, length of project, long-term take limit
  g=0.08, glwr=0.07, gupr=0.09, prho=1, # detection probabilities for future years when management options are constant
  projyears=1947:(1947+30-5-1),
  projrho=rep(1,25),
  projg=rep(0.08,25),
  projglwr=rep(0.07,25),
  projgupr=rep(0.09,25),
  aL = 0.01, styr=3, tau = 2, # alpha value for short-term test, years in moving average window, threshold for fatality rate (annual)
  aR = 0.10, rho = 0.6, # alpha value for reversion, AMA
  aCI = 0.90, # confidence level for lambda
  prior_f='Objective',
  prior_M=diff(sqrt(0:999))/sum(diff(sqrt(0:999))) # this may not run high enough if count is not small, but it is adjusted for the calculations
)

scexDefault<-list(
  ltT = 1, # is long-term trigger part of the analysis
  stT = 1, # is short-term trigger part of the analysis
  rT = 0, # is reversion trigger part of the analysis
  nyr = 30, # total years in project
  Tau = 60, # total permitted take
  lambda = 2, # baseline fatality rate (at rho = 1)
  rhoinf = 0, # effectiveness of final AMA (rhoinf = 0 => fatalities cease upon implementation of the AMA) --> NOTE: last iAMAschedule rho should be == rhoinf too
  yrs  = 3, # initial years of intensive monitoring at the start of a project
  g1 = 0.3, # detection probability under intensive monitoring (e.g., cleared plots)
  g1lwr = 0.25, # lower bound of 95% CI for detection probability under intensive monitoring
  g1upr = 0.35, # uppper bound of 95% CI for detection probability under intensive monitoring
  g2 = 0.08, # detection probability under non-intensive monitoring (e.g., roads and pads)
  g2lwr = 0.07, # lower bound of 95% CI for detection probability under non-intensive monitoring
  g2upr = 0.09, # uppper bound of 95% CI for detection probability under nonintensive monitoring
  gpost = 0,  #continue monitoring after long-term triggering?
  g3i = 2,    # monitoring continues at intensive level (1 for intensive, 2 for non-intensive, 3 for no monitoring)
  nsim = 10, # number of simulation draws (greater number gives better results but takes more time)
  lta = 0.5, # alpha value for determining LT compliance
  sta = 0.01, # alpha value for determining ST compliance
  sty = 3, # number of years in moving window for estimating short-term lambda
  stTon = 'same', # short-term trigger = Tau/n ('same') or short-term trigger != Tau/n ('different')
  tau = 2,
  iAMA = F, # is a schedule of incremental AMAs implemented?
  iAMAschedule = cbind(c(0.8, 0.7, 0),c(2,2,0)), # rho and number of years of intensive monitoring associated with each increment in AMA
  rta = 0.1, # alpha value for reversion trigger
  rhorev = 0.6 # assumed effectiveness of AMA under consideration for reversion
)

CP <- data.frame(cbind( # field trial data from actual PA site (moderate visibility class)
c(0.95,4.93,9.22,3.12,0.82,0.08,9.21,5.02,1.88,5.92,0.13,5.68,12.06,0.1,0.06,0.08,0.09,5.13,2.19,12.05,12.28,9,0.15,6.21,12.03,0,0,0,0,0,0,0,0,0,0,22.1,23.17,20.82,21.11,19.97,18.95,19.07),
c(5.99,6.85,12.03,5.14,3.12,2.03,10.34,6.19,3.08,7.07,1.15,8.71,13.11,1.14,1.19,2.02,0.96,8.13,5.19,15.13,15.29,12.05,2.98,9.19,14.78,3.11,5.94,1.16,0.9,0.11,3.28,2.92,3.06,3.3,2.96,Inf,Inf,Inf,Inf,Inf,Inf,Inf)
))
names(CP)<-c("CPmin", "CPmax")
CP$CPmin <- pmax(CP$CPmin, 0.001)
event<-ifelse(CP$CPmin == CP$CPmax, 1, ifelse(CP$CPmax == Inf, 0, 3))
left<-CP$CPmin
right<-CP$CPmax
right[event==0] <- CP$CPmin[event==0]
surv<-survival::Surv(time=left,time2=right, event=event,type=c('interval'))
mod.e<-survival::survreg(surv~1, dist="exponential")
mod.w<-survival::survreg(surv~1, dist="weibull")
mod.ll<-survival::survreg(surv~1, dist="loglogistic")
mod.ln<-survival::survreg(surv~1, dist="lognormal")
#CPr<-rCPgab("Weibull", 1/mod.w$scale, exp(mod.w$coef[1]), ifelse(singleYearDefault$samtype=="Formula", singleYearDefault$Isam, getmode(diff(singleYearDefault$days))))
CPdataDefault<-list(
  CP = CP,
  surv = surv,
  persistence_distn = "Lognormal",
  mod.e = mod.e,
  mod.w = mod.w,
  mod.ll = mod.ll,
  mod.ln = mod.ln,
  pda = 1/mod.w$scale,
  pdb = exp(mod.w$coef[1]),
  blwr = exp(mod.w$coef+qt(0.025,mod.w$df.res)*sqrt(mod.w$var[1])),
  bupr = exp(mod.w$coef+qt(0.975,mod.w$df.res)*sqrt(mod.w$var[1]))
# , meanCP = CPr[1],
#  rhat = CPr[2]
)
pkdatDefault<-list(
  n = 12,
  X = c(98,  23,  6,  7,  1,  1, 2, 1, 0, 0, 1, 0),
  M = c(158, 55, 31, 21, 13, 10, 9, 7, 5, 4, 4, 3)
)

pkmod <- "
model {
  for (i in 1:n) {
    X[i] ~ dbin(p*k^(i-1), M[i])
  }
  p ~ dbeta(1/2, 1/2)
  k ~ dbeta(1, 1)
}
"

## create names of tcl variables that need to be stored for all data that need to be read into forms
syVar<-c(
  # arrivals
  'a', 'v', 'arrfun', 'arrstart',
  'lwr.u', 'upr.u', 'wt.u',
  'lwr.p1', 'upr.p1', 'wt.p1', 'a.p1', 'b.p1',
  'lwr.p2', 'upr.p2', 'wt.p2', 'a.p2', 'b.p2',
  # search
  'samtype', 'firstsearch', 'Isam', 'SEopt', 'nsearch', 'SEn', 'SEx', 'k',
  # persistence
  'persistence_distn', 'perstype', 'pda', 'pdb', 'blwr', 'bupr',
  'pda.w', 'pdb.w', 'blwr.w', 'bupr.w',
  'pda.ll', 'pdb.ll', 'blwr.ll', 'bupr.ll',
  'pda.ln', 'pdb.ln', 'blwr.ln', 'bupr.ln',
  'pda.e', 'pdb.e', 'blwr.e', 'bupr.e',
  # estimation of M
  'X', 'crlev', 'prior_f'
)
syArray<-c(
  # arrival and search
  'arrcomponents', 'days',
  # priors
  'prior_M', 'custom.prior', 'objective.prior'
)
#############
symcVar<-c(
  'ICEoption', 'option', 'prior_f', 'crlev',
  'firstsearch', 'span',
  'arrfun', 'arrstart',
  'lwr.u', 'upr.u', 'wt.u',
  'lwr.p1', 'upr.p1', 'wt.p1', 'a.p1', 'b.p1',
  'lwr.p2', 'upr.p2', 'wt.p2', 'a.p2', 'b.p2'
)
symcArray<-c(
  'arrcomponents',
  'class', 'rel_wt', 'Xmc', 'Ba', 'Bb',
  'prior_M', 'custom.prior', 'objective.prior'
)
#############
scexVar<-c(
  'g1', 'g1lwr', 'g1upr', 'g2', 'g2lwr', 'g2upr', 'g3i', 'gpost',
  'lambda', 'lta', 'ltT', 'yrs', 'iAMA',
  'nsim', 'nyr',
  'rhoinf', 'rhorev', 'rT', 'rta', 'sta', 'stT', 'sty', 'stTon', 'tau', 'Tau'
)
scexArray<-c(
  'iAMAschedule'
)
#############
dtVar<-c(
  'grinc',
  'target_g', 'tarr',
  'ffix', 'f', 'fmin', 'fmax',
  'phifix', 'phi', 'phimin', 'phimax',
  'Ifix', 'Isam', 'Imin', 'Imax', 'span', 'firstsearch', 'k',
  'persistence_distn',
  'pda0',   'pdb0',
  'pda.w',  'pdb.w',
  'pda.ll', 'pdb.ll',
  'pda.ln', 'pdb.ln',
  'pda.e',  'pdb.e',
  'arrfun', 'arrstart',
  'lwr.u', 'upr.u', 'wt.u',
  'lwr.p1', 'upr.p1', 'wt.p1', 'a.p1', 'b.p1',
  'lwr.p2', 'upr.p2', 'wt.p2', 'a.p2', 'b.p2',
  'tau', 'crlev', 'X', 'dfoption'
)
dtArray<-c(
  'arrcomponents'
)
#############
myVar<-c(
  'aCI', 'aL', 'crlev', 'aR',
  'g', 'glwr', 'gupr', 'prho',
  'Ltype', 'Mtype', 'nyr', 'Ptype', 'option',
  'prior_f',
  'rho',
  'styr',
  'tau',
  'Tau',
  'years'
)
myArray<-c(
  'X', 'Ba', 'Bb','prior_M',
  'projyears', 'projg', 'projglwr', 'projgupr', 'projrho',
  'rel_wt'
)
gdatDefault<-list(
  # NOTE: the format for exact = F differs from format for exact = T b/c uncertainty parameters are not needed in the former case
  a = .4,  # spatial coverage = fraction of carcasses arriving in the searched area
  v = 0.75, # temporal coverage =  fraction of carcasses arriving in the monitored period
  # NOTE: v parameter is ignored if compound or custom arrival function is specified (temporal coverage is dictated by arrival function)
  p = 0.3, # searcher efficiency = conditional probability of finding a (fresh) carcass given that it is present at time of search
  k = 0.7, # factor by which searcher efficiency changes with each successive search; k ?[0, 1]
  samtype="Formula", # are search dates regular ("Formula") or irregular ("Custom")
  # NOTE: calculations are faster if dates are entered by formula (search interval and number of searches)
  Isam = 7, nsearch = 26, # format for entering searches by "Formula"
  days = c(0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 70),  # format for entering "Custom" searches
  persistence_distn = "Lognormal" , pda=4.0827, pdb = 1.1707 # format for entering persistence distribution
)
# data
parmsDefault<-list(arrcomponents = c(T, T, T), # which components of the compound arrival function are included?
  lwr.u=0, upr.u=365, wt.u=1/4, # parameters for uniform component
  lwr.p1=30, upr.p1=150, wt.p1=1/4, a.p1=2.5, b.p1=2.5, # parameters for beta1 component
  lwr.p2=200, upr.p2=300, wt.p2=1/2, a.p2=3.5, b.p2=4.5 # parameters for beta2 component
)
#############
### constants (included in production copy) -->
devtools::use_data(
  # primary R data lists
  myDefault, dtDefault, singleYearDefault, symcDefault, scexDefault,
  # secondary R data lists
  syscDefault, CPdataDefault, pkdatDefault, gdatDefault, parmsDefault,
  # vectors of names of variables corresponding to tcl variables (or arrays)
  syVar, syArray, symcVar, symcArray, scexVar, scexArray, dtVar, dtArray, myVar, myArray, #syscVar, syscArray,
  # format constants
  VER, about_text,
  compcolor, errcolor, uclr, p1clr, p2clr, lwrclr, uprclr, monclr, ats0, figlab, basicMode,
  expoColor, weibullColor, llogisColor, lnormColor,
  listboxbgclr, prnames, pdnames, pkmod,
  internal = T, overwrite = T
)
