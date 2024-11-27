SYparmLoadable<-function(parms){
# a function to check whether the basic skeleton of SY parms is valid for loading to the SY form
# check whether radio buttons and listboxes can be properly assigned
# if not, an error message is given and the function exits (F)
# if so, function returns (T)
# no error-checking is performed on other specific parameter values
  if (is.null(parms$arrfun)){
    tkmessageBox(icon='error',message="Error in data:\nArrival function ($arrfun) is missing")
    printParms(parms,filename)
    return(F)
  } else if (parms$arrfun!="Uniform" & parms$arrfun!="Compound"){
    tkmessageBox(icon='error',message="Error in data:\nArrival function ($arrfun) must be Uniform or Compound")
    printParms(parms,filename)
    return(F)
  } else {
    ### need to write in a routine for handling arrival function data
  }
  if (is.null(parms$samtype)){
    tkmessageBox(icon='error',message="Error in data:\nSearch type ($samtype) is missing")
    printParms(parms,filename)
    return(F)
  } else {
    if (parms$samtype != "Formula" & parms$samtype != "Custom") {
      tkmessageBox(icon='error',message='Error in data:\nSearch type ($samtype) must be "Formula" or "Custom"')
      printParms(parms,filename)
      return(F)
    }
  }
  if (is.null(parms$persistence_distn)){
    tkmessageBox(icon='error',message="Error in data:\nPersistence distribution ($persistence_distn) is missing")
    printParms(parms,filename)
    return(F)
  } else if (parms$persistence_distn!="Exponential" & parms$persistence_distn!="Weibull" & parms$persistence_distn!="Lognormal" & parms$persistence_distn!="Log-Logistic"){
    tkmessageBox(icon='error',message="Error in data:\nPersistence distribution ($persistence_distn) must be Exponential, Weibull, Lognormal, or Log-Logistic")
    printParms(parms,filename)
    return(F)
  }
  if (is.null(parms$prior_f)){
    tkmessageBox(icon='error',message="Error in data:\nPrior distribution ($prior_f) is missing")
    printParms(parms,filename)
    return(F)
  } else if (parms$prior_f!="Objective" & parms$prior_f!="Custom"){
    tkmessageBox(icon='error',message="Error in data:\nPrior distribution ($prior_f) must be Objective or Custom")
    printParms(parms,filename)
    return(F)
  }
  return(T)
}
arrivalFun<-function(sydat){
  .Rvar$arrProcess <- tktoplevel()
  tktitle(.Rvar$arrProcess) <- "Arrival process"
  tkgrab.set(.Rvar$arrProcess);  tkfocus(.Rvar$arrProcess)
  tkconfigure(.Rvar$arrProcess,width=1000)
  tkwm.resizable(.Rvar$arrProcess,0,0)
  tkwm.deiconify(.Rvar$arrProcess)
  arrcomponents<-toR(sydat$tkvars$arrcomponents)
  Udo <-arrcomponents[1]
  P1do<-arrcomponents[2]
  P2do<-arrcomponents[3]
  dateset<-T

  barwidth<-5
  disclr<-colors()[353]
  parmclr<-'white'
  wbord<-1
  scwid<-900
  relief<-"flat"
  xx<-seq(0.005,.995,length=1000)
  a.p1<-toR(sydat$tkvars$a.p1)
  b.p1<-toR(sydat$tkvars$b.p1)
  mu.p1<-tclVar(a.p1/(a.p1+b.p1))
  s2.p1<-tclVar(log(12*a.p1*b.p1/((a.p1+b.p1)^2*(a.p1+b.p1+1))))
  a.p2<-toR(sydat$tkvars$a.p2); b.p2<-toR(sydat$tkvars$b.p2)
  mu.p2<-tclVar(a.p2/(a.p2+b.p2))
  s2.p2<-tclVar(log(12*a.p2*b.p2/((a.p2+b.p2)^2*(a.p2+b.p2+1))))

  arrparms<-tclArray()
  columnNames<-c("","Population", "wt", "start", "end", "Ba", "Bb")
  for (i in 1:length(columnNames)) arrparms[[0,i-1]]<-strsplit(columnNames[i]," ",fixed=T)[[1]]
  # unsearched area
  arrparms[[1,0]]<-as.tclObj('',drop=T)
  arrparms[[1,1]]<-"resident"; arrparms[[1,2]]<-tclvalue(sydat$tkvars$wt.u) ; arrparms[[1,3]]<-tclvalue(sydat$tkvars$lwr.u) ; arrparms[[1,4]]<-tclvalue(sydat$tkvars$upr.u); arrparms[[1,5]]<-1   ; arrparms[[1,6]]<-1
  arrparms[[2,1]]<-"pulse1";   arrparms[[2,2]]<-round(as.numeric(tclvalue(sydat$tkvars$wt.p1)),3); arrparms[[2,3]]<-round(as.numeric(tclvalue(sydat$tkvars$lwr.p1)),3); arrparms[[2,4]]<-round(as.numeric(tclvalue(sydat$tkvars$upr.p1)),3);
  arrparms[[2,5]]<-signif(a.p1,4); arrparms[[2,6]]<-signif(b.p1,4)
  arrparms[[3,1]]<-"pulse2";   arrparms[[3,2]]<-round(as.numeric(tclvalue(sydat$tkvars$wt.p2)),3); arrparms[[3,3]]<-round(as.numeric(tclvalue(sydat$tkvars$lwr.p2)),3); arrparms[[3,4]]<-round(as.numeric(tclvalue(sydat$tkvars$upr.p2)),3)
  arrparms[[3,5]]<-signif(a.p2,4); arrparms[[3,6]]<-signif(b.p2,4)

  parmTable<-tk2table(.Rvar$arrProcess,
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
  if (Udo) tcl(parmTable,"tag", "rowtag", "U","1") else tcl(parmTable,"tag", "rowtag", "hide","1")
  if (P1do) tcl(parmTable,"tag", "rowtag", "P1","2") else tcl(parmTable,"tag", "rowtag", "hide","2")
  if (P2do) tcl(parmTable,"tag", "rowtag", "P2","3") else tcl(parmTable,"tag", "rowtag", "hide","3")
  tcl(parmTable,"tag","configure", "error",bg=colors()[652]) # cells with error have yellow background

#  ats0 <- c(0,31,59,90,120,151,181,212,243,273,304,334)
  figlab<-format(as.Date(ats0,origin="1970-01-01"),"%b-%d")
#  startDate <- tclVar("0")
  ats <- (ats0-as.numeric(tclvalue(sydat$tkvars$arrstart)))%%365
  startDateFrame<-tkframe(.Rvar$arrProcess)
  startlbl <- tklabel(.Rvar$arrProcess,
    text = paste0("Starting date: ",format(as.Date(as.numeric(tclvalue(sydat$tkvars$arrstart)),
    origin="1970-01-01"),"%b-%d")),
    width=20,
    anchor="w")
  # A function that changes the label
  startDateOnSlide <- function(...) {
    tkconfigure(startlbl, text = paste0("Starting date: ",format(as.Date(as.numeric(tclvalue(sydat$tkvars$arrstart)),origin="1970-01-01"),"%b-%d")))
    ats <<- (ats0-as.numeric(tclvalue(sydat$tkvars$arrstart)))%%365
    assign('ats',(ats0-as.numeric(tclvalue(sydat$tkvars$arrstart)))%%365, .Rvar)
    tkrplot::tkrreplot(arrfig)
  }
  # Add the slider
  startDateSlider <- tkscale(startDateFrame, from = 0, to = 364, variable = sydat$tkvars$arrstart, orient = "horizontal", length = scwid,
    command = startDateOnSlide,
    resolution = 1,
    showvalue=F,
    troughcolor=colors()[652],
    relief='flat',
    background=colors()[641],
    active=colors()[641],
    bd=0
  )
  plotarr <- function() { # not currently scaled to endpoints
    if(!dateset){
      par(mar=c(4,1.5,.5,1.5),mgp=c(2,.7,0))
      plot(0,0,axes=F,xlab='Date',ylab='',type='n',xlim=c(0,365))
      mtext(side=2,line=.7,"                                  ")
      lines(-rep(as.numeric(tclvalue(sydat$tkvars$arrstart)),2)%%365,par('usr')[3:4],lty=3, cex=.Rvar$charSizeAdjust)
      box()
      ats <<- (ats0-as.numeric(tclvalue(sydat$tkvars$arrstart)))%%365
      axis(1, at=ats, lab=figlab)
      return(F)
    }
    if (!Udo & !P1do & !P2do){
      par(mar=c(4,1.5,.5,1.5),mgp=c(2,.7,0))
      plot(0,0,axes=F,xlab='Date',ylab='',type='n',xlim=c(0,365),xaxs='i',yaxs='i')
      mtext(side=2,line=.7,"Relative Arrival Rate", cex=.Rvar$charSizeAdjust)
      box()
      axis(1, at=ats, lab=figlab)
      return(F)
    }
    lwr.u<-toR(sydat$tkvars$lwr.u);   upr.u<-toR(sydat$tkvars$upr.u)
    lwr.p1<-toR(sydat$tkvars$lwr.p1); upr.p1<-toR(sydat$tkvars$upr.p1)
    lwr.p2<-toR(sydat$tkvars$lwr.p2); upr.p2<-toR(sydat$tkvars$upr.p2)
    mu.p1 <- as.numeric(tclvalue(mu.p1)); s2.p1<- exp(as.numeric(tclvalue(s2.p1))+log(1/12))
    a.p1 <<- mu.p1^2/s2.p1*(1-mu.p1)-mu.p1; b.p1<<-a.p1*(1/mu.p1-1)
    arrparms[[2,5]]<<-signif(a.p1,4); arrparms[[2,6]]<<-signif(b.p1,4)
    mu.p2 <- as.numeric(tclvalue(mu.p2)); s2.p2<- exp(as.numeric(tclvalue(s2.p2))+log(1/12))
    a.p2 <<- mu.p2^2/s2.p2*(1-mu.p2)-mu.p2; b.p2<<-a.p2*(1/mu.p2-1)
    arrparms[[3,5]]<<-signif(a.p2,4); arrparms[[3,6]]<<-signif(b.p2,4)
    wt.u<-toR(sydat$tkvars$wt.u)*Udo
    wt.p1<-toR(sydat$tkvars$wt.p1)*P1do
    wt.p2<-toR(sydat$tkvars$wt.p2)*P2do
    totw <- wt.u + wt.p1 + wt.p2
    wt.u<-wt.u/totw; wt.p1<-wt.p1/totw; wt.p2<-wt.p2/totw
    xx<-seq(0.1,364.9,length=1000)
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
      ind<-which(xx>=lwr.u&xx<=upr.u)
      part3[ind]<-wt.u/(upr.u-lwr.u)
    }
    arf<-part1 + part2 + part3
    ymax<-max(arf)
    if (is.finite(ymax)){
      par(mar=c(4,1.5,.5,1.5),mgp=c(2,.7,0))
      plot(0,0,axes=F,xlab='Date',ylab='',type='n',xlim=c(1,365),xaxs='i',ylim=c(0,1.05)*ymax,yaxs='i')
      mtext(side=2,line=.7,"Relative Arrival Rate", cex=.Rvar$charSizeAdjust)
      box()
      axis(1, at=ats, lab=figlab)

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
      xx<-seq(0.1,364.9,length=1000)
      lines(xx,arf, lwd=2, col=1)
    } else {
      par(mar=c(4,1.5,.5,1.5),mgp=c(2,.7,0))
      plot(0,0,axes=F,xlab='time',ylab='',type='n',xlim=c(0,365),xaxs='i',yaxs='i')
      mtext(side=2,line=.7,"Relative Arrival Rate", cex=.Rvar$charSizeAdjust)
      box()
      axis(1, at=ats, lab=figlab)
      return(F)
    }
  }
  arrfig <- tkrplot::tkrplot(.Rvar$arrProcess, fun = plotarr, hscale=2.2*.Rvar$charSizeAdjust, vscale=1.3*.Rvar$charSizeAdjust)
  onChange <- function(...) {
    tkrplot::tkrreplot(arrfig)
  }
  # radio buttons for defining model form
  modelDefinition<-tkframe(.Rvar$arrProcess) # frame for holding the model definition options
  Ulab<-tklabel(modelDefinition,text="uniform")
  Uradio<-tclVar(ifelse(arrcomponents[1],"yes","no"))
  UyesRadio<-tkradiobutton(modelDefinition,variable=Uradio,value="yes")
  UnoRadio<-tkradiobutton(modelDefinition,variable=Uradio,value="no")
  P1lab<-tklabel(modelDefinition,text="pulse1")
  P1radio<-tclVar(ifelse(arrcomponents[2],"yes","no"))
  P1yesRadio<-tkradiobutton(modelDefinition,variable=P1radio,value="yes")
  P1noRadio<-tkradiobutton(modelDefinition,variable=P1radio,value="no")
  P2lab<-tklabel(modelDefinition,text="pulse2")
  P2radio<-tclVar(ifelse(arrcomponents[3],"yes","no"))
  P2yesRadio<-tkradiobutton(modelDefinition,variable=P2radio,value="yes")
  P2noRadio<-tkradiobutton(modelDefinition,variable=P2radio,value="no")
  tkgrid(Ulab,UyesRadio,UnoRadio,rowspan=2)
  tkgrid(P1lab,P1yesRadio,P1noRadio,rowspan=2)
  tkgrid(P2lab,P2yesRadio,P2noRadio,rowspan=2)

  processBounds<-tkframe(.Rvar$arrProcess)

  Ulwr <- tkscale(processBounds, from = 0, to = 365, variable = sydat$tkvars$lwr.u, orient = "horizontal",
    length = scwid,
    width=barwidth,
    command = function(...) {
      if (toR(sydat$tkvars$lwr.u)>=toR(sydat$tkvars$upr.u)){
        tclvalue(sydat$tkvars$lwr.u)<-toR(sydat$tkvars$upr.u)-1
        arrparms[[1,3]]<<-round(toR(sydat$tkvars$lwr.u))
        return(F)
      }
      arrparms[[1,3]]<<-round(toR(sydat$tkvars$lwr.u))
      onChange()
    },
    resolution = 1,
    sliderrelief=relief,
    sliderlength=8,
    troughcolor=uclr,
    showvalue=F,
    borderwidth=wbord
  )
  Uupr <- tkscale(processBounds, from = 0, to = 365, variable = sydat$tkvars$upr.u, orient = "horizontal",
    length = scwid,
    width=barwidth,
    command = function(...) {
      if (toR(sydat$tkvars$lwr.u)>=toR(sydat$tkvars$upr.u)){
        tclvalue(sydat$tkvars$upr.u)<-toR(sydat$tkvars$lwr.u)+1
        arrparms[[1,4]]<<-round(toR(sydat$tkvars$lwr.u))
        return(F)
      }
      arrparms[[1,4]]<<-round(toR(sydat$tkvars$upr.u))
      onChange()
    },
    resolution = 1,
    sliderrelief=relief,
    sliderlength=8,
    troughcolor=uclr,
    showvalue=F,
    borderwidth=wbord
  )
  P1lwr <- tkscale(processBounds, from = 0, to = 365, variable = sydat$tkvars$lwr.p1, orient = "horizontal",
    length = scwid,
    width=barwidth,
    command = function(...) {
      if (toR(sydat$tkvars$lwr.p1)>=toR(sydat$tkvars$upr.p1)){
        tclvalue(sydat$tkvars$lwr.p1)<-toR(sydat$tkvars$upr.p1)-1
        arrparms[[2,3]]<<-round(toR(sydat$tkvars$lwr.p1))
        return(F)
      }
      arrparms[[2,3]]<<-round(toR(sydat$tkvars$lwr.p1))
      onChange()
    },
    resolution = 1,
    sliderrelief=relief,
    sliderlength=8,
    troughcolor=p1clr,
    showvalue=F,
    borderwidth=wbord
  )
  P1upr <- tkscale(processBounds, from = 0, to = 365, variable = sydat$tkvars$upr.p1, orient = "horizontal",
    length = scwid,
    width=barwidth,
    command = function(...) {
      if (toR(sydat$tkvars$lwr.p1)>=toR(sydat$tkvars$upr.p1)){
        tclvalue(sydat$tkvars$upr.p1)<-toR(sydat$tkvars$lwr.p1)+1
        arrparms[[2,4]]<<-round(toR(sydat$tkvars$upr.p1))
        return(F)
      }
      arrparms[[2,4]]<<-round(toR(sydat$tkvars$upr.p1))
      onChange()
    },
    resolution = 1,
    sliderrelief=relief,
    sliderlength=8,
    troughcolor=p1clr,
    showvalue=F,
    borderwidth=wbord
  )
  P2lwr <- tkscale(processBounds, from = 0, to = 365, variable = sydat$tkvars$lwr.p2, orient = "horizontal",
    length = scwid,
    width=barwidth,
    command = function(...) {
      if (toR(sydat$tkvars$lwr.p2)>=toR(sydat$tkvars$upr.p2)){
        tclvalue(sydat$tkvars$lwr.p2)<-toR(sydat$tkvars$upr.p2)-1
        arrparms[[3,3]]<<-round(toR(sydat$tkvars$lwr.p2))
        return(F)
      }
      arrparms[[3,3]]<<-round(toR(sydat$tkvars$lwr.p2))
      onChange()
    },
    resolution = 1,
    sliderrelief=relief,
    sliderlength=8,
    troughcolor=p2clr,
    showvalue=F,
    borderwidth=wbord
  )
  P2upr <- tkscale(processBounds, from = 0, to = 365, variable = sydat$tkvars$upr.p2, orient = "horizontal",
    length = scwid,
    width=barwidth,
    command = function(...) {
      if (toR(sydat$tkvars$lwr.p2)>=toR(sydat$tkvars$upr.p2)){
        tclvalue(sydat$tkvars$upr.p2)<-toR(sydat$tkvars$lwr.p2)+1
        arrparms[[3,4]]<<-round(toR(sydat$tkvars$upr.p2))
        return(F)
      }
      arrparms[[3,4]]<<-round(toR(sydat$tkvars$upr.p2))
      onChange()
    },
    resolution = 1,
    sliderrelief=relief,
    sliderlength=8,
    troughcolor=p2clr,
    showvalue=F,
    borderwidth=wbord
  )
  tkgrid(Ulwr)
  tkgrid(Uupr)
  tkgrid(P1lwr)
  tkgrid(P1upr)
  tkgrid(P2lwr)
  tkgrid(P2upr)

  ## parameter summary for bounds
  processWts<-tkframe(.Rvar$arrProcess)
  scht<-420
  setWts<-function(){
    ww<-toR(sydat$tkvars$wt.u)
    wo<-toR(sydat$tkvars$wt.p1)*P1do+toR(sydat$tkvars$wt.p2)*P2do
    arrparms[[1,2]]<<-round(ifelse(wo==0,1,ww/(wo+ww)),2)
    ww<-toR(sydat$tkvars$wt.p1)
    wo<-toR(sydat$tkvars$wt.u)*Udo+toR(sydat$tkvars$wt.p2)*P2do
    arrparms[[2,2]]<<-round(ifelse(wo==0,1,ww/(wo+ww)),2)
    ww<-toR(sydat$tkvars$wt.p2)
    wo<-toR(sydat$tkvars$wt.u)*Udo+toR(sydat$tkvars$wt.p1)*P1do
    arrparms[[3,2]]<<-round(ifelse(wo==0,1,ww/(wo+ww)),2)
  }
  Uwt<-tkscale(processWts, from = 1, to = 0, variable = sydat$tkvars$wt.u, orient = "vertical",
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
  P1wt<-tkscale(processWts, from = 1, to = 0, variable = sydat$tkvars$wt.p1, orient = "vertical",
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
  P2wt<-tkscale(processWts, from = 1, to = 0, variable = sydat$tkvars$wt.p2, orient = "vertical",
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
  tkgrid(Uwt,P1wt,P2wt,pady=10)

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
  processVariances<-tkframe(.Rvar$arrProcess)
  P1s2<-tkscale(processVariances, from = -6, to = 1.05, variable = s2.p1, orient = "vertical",
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
  P2s2<-tkscale(processVariances, from = -6, to = 1.05, variable = s2.p2, orient = "vertical",
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
  tkgrid(P1s2,P2s2,pady=10)

  arrButtFrame<-tkframe(.Rvar$arrProcess)
  arrReset<-tkbutton(.Rvar$arrProcess,text="Reset Start Date",command=function(){
    tkconfigure(arrSave,state='disabled')
    tkconfigure(arrRestore,state='disabled')
  # data
    dateset<<-F
    tclvalue(sydat$tkvars$arrstart)<-0
    tkconfigure(UyesRadio,state='disabled')
    tkconfigure(UnoRadio,state='disabled')
    tkconfigure(P1yesRadio,state='disabled')
    tkconfigure(P1noRadio,state='disabled')
    tkconfigure(P2yesRadio,state='disabled')
    tkconfigure(P2noRadio,state='disabled')
    tkgrid(startOK, startDateSlider, pady=5,padx=10) # this bad boy gets gridded when reset button is hit
    tkgrid(startDateFrame,sticky='nw',columnspan=4,column=1,row=0)
    tkgrid.remove(arrReset)
  })
  arrRestore<-tkbutton(arrButtFrame,text="Restore Defaults", width=15,command=function(){
    adat<-singleYearDefault
    dateset<<-T
    Udo<<-T; P1do<<-T; P2do<<-T
    tclvalue(sydat$tkvars$arrstart)<-adat$arrstart
    tclvalue(sydat$tkvars$lwr.u)<-adat$lwr.u
    tclvalue(sydat$tkvars$upr.u)<-adat$upr.u
    tclvalue(sydat$tkvars$lwr.p1)<-adat$lwr.p1
    tclvalue(sydat$tkvars$upr.p1)<-adat$upr.p1
    tclvalue(sydat$tkvars$lwr.p2)<-adat$lwr.p2
    tclvalue(sydat$tkvars$upr.p2)<-adat$upr.p2
    a.p1<-adat$a.p1; b.p1<-adat$b.p1
    tclvalue(mu.p1)<-a.p1/(a.p1+b.p1)
    tclvalue(s2.p1)<-log(12*a.p1*b.p1/((a.p1+b.p1)^2*(a.p1+b.p1+1)))
    a.p2<-adat$a.p2; b.p2<-adat$b.p2
    tclvalue(mu.p2)<-a.p2/(a.p2+b.p2)
    tclvalue(s2.p2)<-log(12*a.p2*b.p2/((a.p2+b.p2)^2*(a.p2+b.p2+1)))
    tclvalue(sydat$tkvars$wt.u)<-adat$wt.u
    tclvalue(sydat$tkvars$wt.p1)<-adat$wt.p1
    tclvalue(sydat$tkvars$wt.p2)<-adat$wt.p2
    arrparms[[1,0]]<-as.tclObj('',drop=T)
    arrparms[[1,1]]<-"resident"; arrparms[[1,2]]<-tclvalue(sydat$tkvars$wt.u) ; arrparms[[1,3]]<-tclvalue(sydat$tkvars$lwr.u) ; arrparms[[1,4]]<-tclvalue(sydat$tkvars$upr.u); arrparms[[1,5]]<-1   ; arrparms[[1,6]]<-1
    arrparms[[2,1]]<-"pulse1";   arrparms[[2,2]]<-round(toR(sydat$tkvars$wt.p1),3); arrparms[[2,3]]<-round(toR(sydat$tkvars$lwr.p1),3); arrparms[[2,4]]<-round(toR(sydat$tkvars$upr.p1),3);
    arrparms[[2,5]]<-signif(a.p1,4); arrparms[[2,6]]<-signif(b.p1,4)
    arrparms[[3,1]]<-"pulse2";   arrparms[[3,2]]<-round(toR(sydat$tkvars$wt.p2),3); arrparms[[3,3]]<-round(toR(sydat$tkvars$lwr.p2),3); arrparms[[3,4]]<-round(toR(sydat$tkvars$upr.p2),3)
    arrparms[[3,5]]<-signif(a.p2,4); arrparms[[3,6]]<-signif(b.p2,4)
  ### form
  # set scales colors to gray
    tkconfigure(Uwt,state='normal',troughcolor=uclr)
    tkconfigure(Ulwr,state='normal',troughcolor=uclr)
    tkconfigure(Uupr,state='normal',troughcolor=uclr)
    tkconfigure(P1wt,state='normal',troughcolor=p1clr)
    tkconfigure(P1lwr,state='normal',troughcolor=p1clr)
    tkconfigure(P1upr,state='normal',troughcolor=p1clr)
    tkconfigure(P1mu,state='normal',troughcolor=p1clr)
    tkconfigure(P1s2,state='normal',troughcolor=p1clr)
    tkconfigure(P2wt,state='normal',troughcolor=p2clr)
    tkconfigure(P2lwr,state='normal',troughcolor=p2clr)
    tkconfigure(P2upr,state='normal',troughcolor=p2clr)
    tkconfigure(P2mu,state='normal',troughcolor=p2clr)
    tkconfigure(P2s2,state='normal',troughcolor=p2clr)
  # hide table values
    tcl(parmTable,"tag", "rowtag", "U","1")
    tcl(parmTable,"tag", "rowtag", "P1","2")
    tcl(parmTable,"tag", "rowtag", "P2","3")
  # disable radio buttons (enable again upon OK for starting date)
    tclvalue(Uradio)<-"yes"
    tclvalue(P1radio)<-"yes"
    tclvalue(P2radio)<-"yes"
    tkrplot::tkrreplot(arrfig)
  })
  arrSave<-tkbutton(arrButtFrame,text="Save", width=15,command=function(){
    tclvalue(sydat$tkvars$arrfun)<<-"Compound"
    sydat$tkvars$arrcomponents[[0]]<<-Udo; sydat$tkvars$arrcomponents[[1]]<<-P1do; sydat$tkvars$arrcomponents[[2]]<<-P2do
    sydat$arrcomponents<<-c(Udo, P1do, P2do)
    tclvalue(sydat$tkvars$arrstart) <<- tclvalue(sydat$tkvars$arrstart)
    tclvalue(sydat$tkvars$wt.u)  <<- tclvalue(arrparms[[1,2]])
    tclvalue(sydat$tkvars$lwr.u) <<- tclvalue(arrparms[[1,3]])
    tclvalue(sydat$tkvars$upr.u) <<- tclvalue(arrparms[[1,4]])
    tclvalue(sydat$tkvars$wt.p1) <<- tclvalue(arrparms[[2,2]])
    tclvalue(sydat$tkvars$lwr.p1)<<- tclvalue(arrparms[[2,3]])
    tclvalue(sydat$tkvars$upr.p1)<<- tclvalue(arrparms[[2,4]])
    tclvalue(sydat$tkvars$wt.p2) <<- tclvalue(arrparms[[3,2]])
    tclvalue(sydat$tkvars$lwr.p2)<<- tclvalue(arrparms[[3,3]])
    tclvalue(sydat$tkvars$upr.p2)<<- tclvalue(arrparms[[3,4]])
    mu<-as.numeric(tclvalue(mu.p1))
    s2<-exp(as.numeric(tclvalue(s2.p1))+log(1/12))
    tclvalue(sydat$tkvars$a.p1)<<-mu^2/s2*(1-mu)-mu
    tclvalue(sydat$tkvars$b.p1)<<-toR(sydat$tkvars$a.p1)*(1/mu-1)
    mu<-as.numeric(tclvalue(mu.p2))
    s2<-exp(as.numeric(tclvalue(s2.p2))+log(1/12))
    tclvalue(sydat$tkvars$a.p2)<<-mu^2/s2*(1-mu)-mu
    tclvalue(sydat$tkvars$b.p2)<<-toR(sydat$tkvars$a.p2)*(1/mu-1)

    tkrplot::tkrreplot(sydat$arrfig.mini)
    tkdestroy(.Rvar$arrProcess)
  })
  arrCancel<-tkbutton(arrButtFrame,text='Cancel', width=15,command=function()tkdestroy(.Rvar$arrProcess))
  tkgrid(arrRestore)
  tkgrid(arrSave)
  tkgrid(arrCancel)

  tkbind(UyesRadio,"<Button-1>", function() {
    if (!dateset) return(F)
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
    if (!dateset) return(F)
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
    if (!dateset) return(F)
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
    if (!dateset) return(F)
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
    if (!dateset) return(F)
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
    if (!dateset) return(F)
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

  OnStartOK<-function(){
    dateset<<-T
    tkgrid.remove(startDateFrame)
    tkconfigure(UyesRadio,state='normal')
    tkconfigure(UnoRadio,state='normal')
    tkconfigure(P1yesRadio,state='normal')
    tkconfigure(P1noRadio,state='normal')
    tkconfigure(P2yesRadio,state='normal')
    tkconfigure(P2noRadio,state='normal')
    tkconfigure(arrSave,state='normal')
    tkconfigure(arrRestore,state='normal')
    tkgrid(arrReset)
    tkrplot::tkrreplot(arrfig)
  }
  startOK<-tkbutton(startDateFrame,text="OK",command=OnStartOK)
  if (sydat$name != "design") tkgrid(startlbl,arrReset,sticky='w') else tkgrid(startlbl, sticky='w')

  tkgrid(parmTable,columnspan=2,sticky='sw',pady=10,padx=10)
  tkgrid(arrButtFrame, column=1,row=1)
  tkgrid(modelDefinition,processBounds)
  tkgrid(processWts, arrfig, processVariances, sticky='ne',padx=10)
  tkgrid(processMeans, column=1,sticky='n')
  blankframe<-tklabel(.Rvar$arrProcess,text=as.tclObj('        ',drop=T))
  tkgrid(blankframe)
}
checkprior<-function(v){
  if (sum(is.na(v))>0){ # missing data
    tkmessageBox(icon='error',message=paste0("\n\nError: Missing data"))
    tkconfigure(cprOK,state='disabled')
    return(FALSE)
  }
  if (dim(v)[2] != 2){ # two columns?
    tkmessageBox(icon='error',message="\n\nError...Required: two columns [m and p(M = m)]")
    tkconfigure(cprOK,state='disabled')
    return(FALSE)
  }
  if (!is.numeric(v[1,1])|!is.numeric(v[1,2])){ # header? ...ignore header
    v<-v[-1,] #data (possibly) contains header, so delete the first row
  }
  prtmp<-v[,2]
  if (!is.numeric(prtmp)){ # are the probabilities numeric?
    tkmessageBox(icon='error',message=paste0("\n\nError in data\nRequired: numeric probabilities"))
    tkconfigure(cprOK,state='disabled')
    return(FALSE)
  }
  if (!is.numeric(v[,1])){ # are the m values numeric?
    tkmessageBox(icon='error',message=paste0("\n\nError in data\nRequired: numeric m values"))
    tkconfigure(cprOK,state='disabled')
    return(FALSE)
  }
  if (sum(prtmp<0)>0){  # probabilities must all be non-negative
    tkmessageBox(icon='error',message=paste0("\n\nError in data\nRequired: probabilities must be non-negative"))
    tkconfigure(cprOK,state='disabled')
    return(FALSE)
  }
  if (abs(sum(prtmp)-1)>0.0001){  # probabilities must sum to 1
    tkmessageBox(icon='error',message=paste0("\n\nError in data\nRequired: probabilities must sum to 1"))
    tkconfigure(cprOK,state='disabled')
    return(FALSE)
  } else {
    prtmp<-prtmp/sum(prtmp)
  }
  if(sum(v[,1]<0)>0){ # m values must all be non-negative
    tkmessageBox(icon='error',message=paste0("\n\nError in data\nRequired: m values must be non-negative"))
    tkconfigure(cprOK,state='disabled')
    return(FALSE)
  }
  if (sum(abs(round(v[,1])-v[,1]))>0.0001){ # m values must be integers
    tkmessageBox(icon='error',message=paste0("\n\nError in data\nRequired: m values must be integers"))
    tkconfigure(cprOK,state='disabled')
    return(FALSE)
  } else {
    v[,1]<-round(v[,1])
  }
  uval <-unique(diff(v[,1])) # m values must be sequential
  if (length(uval)>1){
    tkmessageBox(icon='error',message=paste0("\n\nError in data\nRequired: m values must be sequential"))
    tkconfigure(cprOK,state='disabled')
    return(FALSE)
  } else if (uval != 1){
    tkmessageBox(icon='error',message=paste0("\n\nError in data\nRequired: m values must be sequential"))
    tkconfigure(cprOK,state='disabled')
    return(FALSE)
  }
  if (length(v[,1])!=max(v[,1])+1){ # if the m values do not start with zero, fill in the start of the array with zeros
    v2<-array(numeric(2*(max(v[,1]+1))),dim=c(max(v[,1])+1,2))
    v2[(min(v[,1])+1):length(v2[,1]),2]<-v[,2]
    v2[,1]<-0:(length(v2[,1])-1)
    v<-v2
  }
#  prtmp<<-prtmp
  return(TRUE)
}
