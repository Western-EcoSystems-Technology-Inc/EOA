# initial options:
tcltk::tcl("option","add","*tearOff",0)
syForm<-R6::R6Class("syForm",
  portable = FALSE,
  public = list(
    # tcl variables:
    tkvars = list(), dynLbl=list(),  sides = NA,
    # R variables:
    Isam = 7, Ir = 7, days = NA, fba = NA, fbb = NA,
    persistence_distn = 'Weibull', # ...for singleYear data w/hand CP parameters (for field CP, tkvars$persistence_distn and .Rvar$CPdat$tkvars$persistence_distn)
    rhat = NA, meanCP = NA,
    prior_f = "Objective",
    prior_M = c(0.8, 0.1, 0.05, 0.03, 0.02),
    myr = NA,
    custom.prior = NA,
    objective.prior = NA,
    prdat = NA,
    arrcomponents = c(1, 0, 0),
    syModule = NA, cdFrame = NA, cprFrame = NA, fcpModule = NA, pkModule = NA,
    rCPhtext1 = NA, rCPhtext2 = NA,
    name = 'sy',
    arrfig.mini = NA,
    Xok = T, aok = T, vok = T, SEnok = T, SExok = T, SExleg = T,
    kok = T, Iok = T, nsearchok = T, Irok = T,
    startok = T, pdaok = T, pdbok = T, bminok = T, bminleg = T, bmaxok = T, bmaxleg = T,
    syAok = T, lamok = T, gook = T,  prok = T, blwrok = T, buprok = T,
    pMgX = NA, pMgX.raw = NA, pMgX.ann = NA, pLpost = NA, M = NA, partial = F,
    # form variables that should be available in other routines?
    prior.lbl = NA, priorViewButton = NA, pkci.lbl = NA, pkr.lbl = NA, cuslab = NA,
    priorEditButton = NA, prior.lbox = NA, SEci.lbl = NA, distrcpr.txt = NA, distr.txt = NA,
    pdtypeHand = NA, pd.lbox = NA, sclbnd = NA, pda.edit = NA, rCP.lbl = NA, distrcpr.lbl = NA,
    # initial values
    initialize = function(sydat, partial = F) { # 'partial' excludes arrival and prior...for sysc case
      # create the required tcl variables from the data input
      partial<<-partial
      if (partial) name <<- 'symc' else name <<- 'sy'
      for (nm in syVar){
        tkvars[[nm]] <<- tclVar(sydat[[nm]])
      }
      for (nm in syArray){
        tmp<-tclArray()
        if (length(sydat[[nm]]) == 0){ # missing data
          tmp[[0]] <- as.tclObj('', drop=T)
        } else if (length(sydat[[nm]]) == 1) { # scalar
          tmp[[0]] <- as.tclObj(ifelse(!is.na(sydat[[nm]]), sydat[[nm]], ''), drop=T)
        } else if (is.null(dim(sydat[[nm]]))) { # vector
          for (j in 1:length(sydat[[nm]])) tmp[[j-1]] <- as.tclObj(sydat[[nm]][j],drop=T)
        } else { # matrix
          for (rowi in 1:dim(sydat[[nm]])[1])
            for (coli in 1:dim(sydat[[nm]])[2])
              tmp[[rowi-1, coli-1]] <- as.tclObj(sydat[[nm]][rowi, coli], drop = T)
        }
       tkvars[[nm]] <<- tmp
      }
      dynLbl[['pkci']]<<-tclVar('pkci')
      dynLbl[['pkr']]<<-tclVar('pkr')
      dynLbl[['distr']]<<-tclVar('distr')
      dynLbl[['distrcpr']]<<-tclVar('distrcpr')
      dynLbl[['seci']]<<-tclVar('seci')
      chkForm.sy<-function(){
        colr<-ifelse(Xchk(toR(tkvars$X)),'white',colors()[652]);       tkconfigure(X.edit, bg=colr)
        colr<-ifelse(SEnchk(toR(tkvars$SEn)),'white',colors()[652]);   tkconfigure(SEn.edit, bg=colr)
        colr<-ifelse(SExchk(toR(tkvars$SEx)),'white',colors()[652]);   tkconfigure(SEx.edit, bg=colr)
        colr<-ifelse(kchk(toR(tkvars$k)),'white',colors()[652]);       tkconfigure(k.edit, bg=colr)
        colr<-ifelse(pdachk(toR(tkvars$pda)),'white',colors()[652]);   tkconfigure(pda.edit, bg=colr)
        colr<-ifelse(pdbchk(toR(tkvars$pdb)),'white',colors()[652]);   tkconfigure(pdb.edit, bg=colr)
        colr<-ifelse(blwrchk(toR(tkvars$blwr)),'white',colors()[652]); tkconfigure(blwr.edit, bg=colr)
        colr<-ifelse(buprchk(toR(tkvars$bupr)),'white',colors()[652]); tkconfigure(bupr.edit, bg=colr)
        colr<-ifelse(achk(toR(tkvars$a)),'white',colors()[652]);       tkconfigure(a.edit, bg=colr)
        colr<-ifelse(vchk(toR(tkvars$v)),'white', colors()[652]);      tkconfigure(v.edit, bg=colr)
        colr<-ifelse(startchk(toR(tkvars$firstsearch)),'white',colors()[652]); tkconfigure(start.edit, bg=colr)
        if (tclvalue(tkvars$samtype)=='Formula'){
          colr<-ifelse(Ichk(toR(tkvars$Isam)),'white',colors()[652]); tkconfigure(I.edit, bg=colr)
          colr<-ifelse(nsearchchk(toR(tkvars$nsearch)),'white',colors()[652]); tkconfigure(nsearch.edit, bg=colr)
        } else {
          days<-toR(tkvars$days)
          if (days[1] != 0) {
            tkmessageBox(message=paste0("Error in search dates"))
            tclvalue(tkvars$samtype)<-"Formula"
            tclvalue(tkvars$Isam)<-singleYearPrevious$Isam
            tclvalue(tkvars$nsearch)<-singleYearPrevious$nsearch
          } else if (min(diff(days)) <= 0){
            tkmessageBox(message=paste0("Error in search dates"))
            tclvalue(tkvars$samtype)<-"Formula"
            tclvalue(tkvars$Isam)<-singleYearPrevious$Isam
            tclvalue(tkvars$nsearch)<-singleYearPrevious$nsearch
          }
        }
        colr<-ifelse(syAchk(toR(tkvars$crlev)),'white',colors()[652]); tkconfigure(crlev.edit, bg=colr)
      #  Ichk(toR(tkvars$Ir))
        if (!basicMode){
          strt<-suppressWarnings(as.numeric(toR(tkvars$arrstart)))
          if (tclvalue(tkvars$arrfun)=="Compound"){
            for (i in 1){
              if (is.na(strt) || strt < 0 || strt >= 365) {
                changeArrFun<-T
                break
              } else {
                arrcomponents<-toR(tkvars$arrcomponents)
                if (length(arrcomponents)!=3)     {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                if (!(arrcomponents[1] %in% 0:1)) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                if (!(arrcomponents[2] %in% 0:1)) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                if (!(arrcomponents[3] %in% 0:1)) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                if (sum(arrcomponents) == 0) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                if (arrcomponents[1]){
                  x<-toR(tkvars$lwr.u)
                  if (length(x)!=1) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (!is.numeric(x)) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (x<0 | x>364) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  y<-toR(tkvars$upr.u)
                  if (length(y)!=1) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (!is.numeric(y)) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (y<0 | y>364) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (x>=y) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  w0<-toR(tkvars$wt.u)
                  if (length(w0)!=1) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (!is.numeric(w0)) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (w0<0 | w0>1) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                }
                if (arrcomponents[2]){
                  x<-toR(tkvars$lwr.p1)
                  if (length(x)!=1) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (!is.numeric(x)) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (x<0 | x>364) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  y<-toR(tkvars$upr.p1)
                  if (length(y)!=1) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (!is.numeric(y)) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (y<0 | y>364) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (x>=y) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  w1<-toR(tkvars$wt.p1)
                  if (length(w1)!=1) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (!is.numeric(w1)) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (w1<0 | w1>1) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  x<-toR(tkvars$a.p1)
                  if (length(x)!=1) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (!is.numeric(x)) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (x<=0) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  x<-toR(tkvars$b.p1)
                  if (length(x)!=1) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (!is.numeric(x)) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (x<=0) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                }
                if (arrcomponents[3]){
                  x<-toR(tkvars$lwr.p2)
                  if (length(x)!=1) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (!is.numeric(x)) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (x<0 | x>364) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  y<-toR(tkvars$upr.p2)
                  if (length(y)!=1) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (!is.numeric(y)) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (y<0 | y>364) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (x>=y) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  w2<-toR(tkvars$wt.p2)
                  if (length(w2)!=1) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (!is.numeric(w2)) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (w2<0 | w2>1) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  x<-toR(tkvars$a.p2)
                  if (length(x)!=1) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (!is.numeric(x)) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (x<=0) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  x<-toR(tkvars$b.p2)
                  if (length(x)!=1) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (!is.numeric(x)) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                  if (x<=0) {tclvalue(tkvars$arrfun)<-"Uniform"; break}
                }
                tclvalue(tkvars$wt.u)<-w0*arrcomponents[1]/(c(w0,w1,w2)%*%arrcomponents)
                tclvalue(tkvars$wt.p1)<-w1*arrcomponents[1]/(c(w0,w1,w2)%*%arrcomponents)
                tclvalue(tkvars$wt.p2)<-w2*arrcomponents[1]/(c(w0,w1,w2)%*%arrcomponents)
              }
            }
          }
          if (tclvalue(tkvars$arrfun)=="Uniform"){
            tclvalue(tkvars$arrstart)<-singleYearPrevious$arrstart
            tclvalue(tkvars$lwr.u)<-singleYearPrevious$lwr.u
            tclvalue(tkvars$upr.u)<-singleYearPrevious$upr.u
            tclvalue(tkvars$wt.u)<-singleYearPrevious$wt.u
            tclvalue(tkvars$lwr.p1)<-singleYearPrevious$lwr.p1
            tclvalue(tkvars$upr.p1)<-singleYearPrevious$upr.p1
            tclvalue(tkvars$wt.p1)<-singleYearPrevious$wt.p1
            tclvalue(tkvars$a.p1)<-singleYearPrevious$a.p1
            tclvalue(tkvars$b.p1)<-singleYearPrevious$b.p1
            tclvalue(tkvars$lwr.p2)<-singleYearPrevious$lwr.p2
            tclvalue(tkvars$upr.p2)<-singleYearPrevious$upr.p2
            tclvalue(tkvars$wt.p2)<-singleYearPrevious$wt.p2
            tclvalue(tkvars$a.p2)<-singleYearPrevious$a.p2
            tclvalue(tkvars$b.p2)<-singleYearPrevious$b.p2
          }
          if (!is.numeric(prior_M)){
            prior_f<<-"Objective"
          } else if (min(prior_M)<0 | sum(prior_M)>1.0001 | sum(prior_M)< 0.9999){
            prior_f<<-"Objective"
          }
          if (prior_f=="Custom") custom.prior<<-prior_M
          if (prior_f=="Objective") {
            prior_M<<-singleYearPrevious$objective.prior
            objective.prior<<-prior_M
          }
        }
      }
      ### R variables stored in the syForm list
      persistence_distn <<- sydat$persistence_distn
      tclvalue(tkvars$persistence_distn) <<- .Rvar$CPdataPrevious$persistence_distn
      if (!basicMode){
        prior_f <<- sydat$prior_f
        prior_M <<- sydat$prior_M
        custom.prior <<- sydat$custom.prior
        objective.prior <<- sydat$objective.prior
        arrcomponents <<- sydat$arrcomponents
      }
      assign('pkdat', .Rvar$pkdatPrevious, env = .Rvar)
      assign('CPdata', .Rvar$CPdataPrevious, env = .Rvar)
      ### build the components of the form
      ## toplevel and menus
      syModule <<- tktoplevel()
      tkwm.protocol(syModule, "WM_DELETE_WINDOW", function(){
        tkdestroy(syModule)
        if (!partial) {
          sysave(.Rvar$singleYearPrevious, .Rvar$CPdataPrevious, .Rvar$pkdatPrevious, .Rvar$csvpath)
          tkwm.deiconify(.Rvar$EoA)
        } else {
          sysave(.Rvar$syscPrevious, .Rvar$CPdataPrevious, .Rvar$pkdatPrevious, .Rvar$csvpath, T)
        }
      })
      if (!partial) tkwm.withdraw(syModule)
      tcl("option","add","*tearOff",0)
      if (partial) tkwm.title(syModule,paste0("EoA, v", .Rvar$VER, " - Search Class")) else tkwm.title(syModule,paste0("EoA, v", .Rvar$VER, " - Single Class Module"))
      tkwm.resizable(syModule,0,0)
      buttwid<-6
      SYtopMenu  <- tkmenu(syModule)
      SYeditMenu <- tkmenu(SYtopMenu,activebackground=colors()[125],activeforeground=colors()[109]); tkconfigure(syModule,menu=SYtopMenu)
      tkadd(SYeditMenu,"command",label="Restore defaults",command=function(){
        if (partial) {
          sydat<-syscDefault
          for (nm in names(sydat)){
            if (nm %in% names(.Rvar$symcdat)){
              sydat[[nm]]<-.Rvar$symcdat[[nm]]
            }
          }
          assign("syscPrevious", sydat, env = .Rvar)
          assign("CPdataPrevious", CPdataDefault, env = .Rvar)
          assign("pkdatPrevious", pkdatDefault, env = .Rvar)
          suppressWarnings(file.remove(paste0(.Rvar$datadir,"/syscPrevious.Rdata")))
          suppressWarnings(file.remove(paste0(.Rvar$datadir,"/CPdataPrevious.Rdata")))
          suppressWarnings(file.remove(paste0(.Rvar$datadir,"/pkdatPrevious.Rdata")))
        } else {
          sydat<-singleYearDefault
          assign("singleYearPrevious", sydat, env = .Rvar)
          assign("CPdataPrevious", CPdataDefault, env = .Rvar)
          assign("pkdatPrevious", pkdatDefault, env = .Rvar)
          suppressWarnings(file.remove(paste0(.Rvar$datadir,"/singleYearPrevious.Rdata")))
          suppressWarnings(file.remove(paste0(.Rvar$datadir,"/CPdataPrevious.Rdata")))
          suppressWarnings(file.remove(paste0(.Rvar$datadir,"/pkdatPrevious.Rdata")))
        }
        tkdestroy(syModule)
        initialize(sydat, partial = partial)
      })
      tkadd(SYeditMenu,"command",label="Restore previous",command=function(){ # restores the data set that was loaded at the start of the session
        tkdestroy(syModule)
        initialize(.Rvar$singleYearPrevious)
      })
      tkadd(SYeditMenu,"command",label="Save to file (.Rdata)",command = function() {
        if (!exists('csvpath', env = .Rvar)) assign('csvpath', getwd(), env = .Rvar)
        filename <- tclvalue(tkgetSaveFile(filetypes = "{{R images} {.Rdata}}", defaultextension = ".Rdata", initialfile = '.Rdata', title = "Save", initialdir = .Rvar$csvpath))
        tmp<-unlist(strsplit(filename,'/')); pathname<-paste(tmp[-length(tmp)],collapse='/')
        if (nchar(pathname)>0) .Rvar$csvpath <- pathname
        if (filename == "") return(FALSE)
        syProvisional<-sy2rds()
        pkdat<-.Rvar$pkdat
        CPdata<-.Rvar$CPdata
        save(syProvisional, CPdata, pkdat, file = filename)
        save(csvpath, file = paste0(.Rvar$datadir,"/csvpath.Rdata"))
        if (partial) tkwm.title(syModule,paste0(tmp[length(tmp)], " - EoA, v", .Rvar$VER, " - Search Class")) else tkwm.title(syModule,paste0(tmp[length(tmp)], " - EoA, v", .Rvar$VER, " - Single Class Module"))
      })
      tkadd(SYeditMenu,"command",label="Read from file (.Rdata)",command= SYreadparm)
      tkadd(SYeditMenu,"command",label="Save", command = function() { # <--- test
      # write tk data to singleYearPrevious and save to .Rdata
        if (syChkAll()) {
          updatePrevious()
          if (!partial) {
            sysave(.Rvar$singleYearPrevious, .Rvar$CPdataPrevious, .Rvar$pkdatPrevious, .Rvar$csvpath)
          } else {
            sysave(.Rvar$syscPrevious, .Rvar$CPdataPrevious, .Rvar$pkdatPrevious, .Rvar$csvpath, T)
          }
        } else {
          tkmessageBox(message="Sorry. Bad data. No save")
        }
      })
      tkadd(SYeditMenu,"command",label="Close",command= function(){
        graphics.off()
        if (!partial) {
          sysave(.Rvar$singleYearPrevious, .Rvar$CPdataPrevious, .Rvar$pkdatPrevious, .Rvar$csvpath)
          tkwm.deiconify(.Rvar$EoA)
        } else {
          sysave(.Rvar$syscPrevious, .Rvar$CPdataPrevious, .Rvar$pkdatPrevious, .Rvar$csvpath, T)
        }
        tkdestroy(syModule)
      })
      tkadd(SYtopMenu,"cascade", label="Edit", menu=SYeditMenu)
#      allok1()  <-- WARNING! Production version needs this (or at least have it taken care of)
      SYhelpMenu <- tkmenu(SYtopMenu,activebackground=colors()[125],activeforeground=colors()[109])
      tkadd(SYhelpMenu, "command", label="Parameter Conversion", command = conversionsCalculator)
      tkadd(SYhelpMenu,"command",label="About",command=function()tkmessageBox(title='Evidence of Absence (EoA)',message=about_text))
      tkadd(SYtopMenu,"cascade",label="Help",menu=SYhelpMenu)
      ##### single year module: master frame skeleton
      #### single year module, first column
      sy1<-tkframe(syModule)
      ### sy1: gFrame = detection probability frame
      w=7
      gFrame<-ttklabelframe(sy1,text="     Detection Probability (g)",padding=10)
      datesFrame<-ttklabelframe(gFrame,text="Search Schedule")
      ## sy1: gFrame : coverage
      aFrame<-tkframe(datesFrame)
      a.lbl<-tklabel(aFrame,text=ifelse(partial,"Fraction of total carcasses\narriving in class (dwp)", "Spatial coverage (a) "))
      a.edit<-tkentry(aFrame,width=6,textvariable=tkvars$a,justify='right',bg='white')
      if (!basicMode) tkgrid(a.lbl,a.edit,sticky='e')
      if (basicMode){
        v.lbl<-tklabel(aFrame,text="Temporal coverage (v) ")
        v.edit<-tkentry(aFrame,width=6,textvariable=tkvars$v,justify='right',bg='white')
        if (!partial) tkgrid(a.lbl,a.edit,sticky='e')
        tkgrid(v.lbl, v.edit, sticky='e')
      } else {
      ###
        arrfFrame<-ttklabelframe(gFrame,text="Arrival Function")
        arrUniformRadio<-tkradiobutton(arrfFrame, variable=tkvars$arrfun, text="Uniform",  value="Uniform", command = function(){
          tclvalue(tkvars$arrfun)<-"Uniform"
          tkrplot::tkrreplot(arrfig.mini)
          tkconfigure(editArrivalsbutton,state='disabled')
          if (name=='symc') {
            tkconfigure(firstsearch.lbl,state='disabled')
            tkconfigure(firstsearch.edit,state='normal')
            tkconfigure(span.lbl,state='normal')
            tkconfigure(span.edit,state='disabled')
          }
        })
        arrCompoundRadio<-tkradiobutton(arrfFrame,variable=tkvars$arrfun, text="Compound", value="Compound", command = function(){
          if (tclvalue(tkcget(arrUniformRadio,"-state"))=="disabled") return(F)
          tclvalue(tkvars$arrfun)<-"Compound"
          tkrplot::tkrreplot(arrfig.mini)
          tkconfigure(editArrivalsbutton,state='normal')
          if (name=='symc') {
            tkconfigure(span.lbl,state='disabled')
            tkconfigure(span.edit,state='disabled')
          }
        })
        buttwid<-6
        viewArrivalsbutton<-tkbutton(arrfFrame,text="View",width=buttwid,command=function() {
        #  arrcomponents<-toR(tkvars$arrcomponents)
          plot_arrivals(tclvalue(tkvars$arrfun),
            toR(tkvars$arrstart), arrcomponents,
            toR(tkvars$lwr.u), toR(tkvars$upr.u), toR(tkvars$wt.u),
            toR(tkvars$lwr.p1),toR(tkvars$upr.p1),toR(tkvars$wt.p1), toR(tkvars$a.p1), toR(tkvars$b.p1),
            toR(tkvars$lwr.p2),toR(tkvars$upr.p2),toR(tkvars$wt.p2), toR(tkvars$a.p2), toR(tkvars$b.p2))
        })
        editArrivalsbutton<-tkbutton(arrfFrame, text="Edit")
        if (tclvalue(tkvars$arrfun)=="Uniform") tkconfigure(editArrivalsbutton, state="disabled")
        tkbind(editArrivalsbutton,"<ButtonPress-1>", function(X,Y){  # edit arrivals button just brings up a window to edit the arrival variables
          if (tclvalue(tcl(editArrivalsbutton,"cget","-state"))=="disabled") return(F)
          assign('geo1', as.numeric(c(X,Y)), env = .Rvar)
          arrivalFun(self)
          tkrplot::tkrreplot(arrfig.mini)
        })
        arrfig.mini <<- tkrplot::tkrplot(arrfFrame, fun = plotarr.mini, hscale=.25, vscale=.08)
        # construct arrival frame
        tkgrid(arrUniformRadio,sticky='w')
        tkgrid(arrCompoundRadio,sticky='w')
        tkgrid(arrfig.mini, column=1,row=0,rowspan=2,padx=10)
        tkgrid(editArrivalsbutton,column=2,row=0,rowspan=2, padx=10)
        if (partial) {
          tkconfigure(arrUniformRadio, state='disabled')
          tkconfigure(arrCompoundRadio, state='disabled')
          tkconfigure(arrfig.mini, state='disabled')
          tkconfigure(editArrivalsbutton, state='disabled')
        }
      }
    ##########################
    ## sy1: gFrame : searcher efficiency
      SEFrame<-ttklabelframe(gFrame,text="     Searcher Efficiency", padding=10, width = 600)
      hSE<-tkradiobutton(SEFrame, text = "Carcasses removed after one search", variable = tkvars$SEopt, value = 'h', command = function(){
        for (obj in unlist(strsplit(tclvalue(tkwinfo("children",fSEframe)), ' '))) tkconfigure(obj,state='disabled')
        for (obj in unlist(strsplit(tclvalue(tkwinfo("children",hSEframe)), ' '))) tkconfigure(obj,state='normal')
      })
      fSE<-tkradiobutton(SEFrame, text = "Carcasses available for several searches", variable = tkvars$SEopt, value = 'f', command = function(){
        for (obj in unlist(strsplit(tclvalue(tkwinfo("children",hSEframe)), ' '))) tkconfigure(obj,state='disabled')
        for (obj in unlist(strsplit(tclvalue(tkwinfo("children",fSEframe)), ' '))) tkconfigure(obj,state='normal')
      })

      hSEframe<-tkframe(SEFrame)
      SEn.lbl<-tklabel(hSEframe,text="Carcasses available ");   SEn.edit<-tkentry(hSEframe,width=4,textvariable=tkvars$SEn,justify='right',bg='white')
      SEx.lbl<-tklabel(hSEframe,text="Carcasses found ");  SEx.edit<-tkentry(hSEframe,width=4,textvariable=tkvars$SEx,justify='right',bg='white')
      SEci.lbl<<-tklabel(hSEframe,width=35)
      if (SEnchk(toR(tkvars$SEn)) && SExchk(toR(tkvars$SEx))){
        SEx<-as.numeric(tclvalue(tkvars$SEx)); SEn<-as.numeric(tclvalue(tkvars$SEn))
        fba<-SEx+0.5; fbb<-SEn-SEx+.5
        tclvalue(dynLbl$seci)<-paste0("p\u0302 = ", round(SEx/SEn,3)," with 95% CI = [" ,round(qbeta(0.025, fba, fbb),3),", ",round(qbeta(0.975, fba, fbb),3),"]")
      } else {
        tclvalue(dynLbl$seci)<-paste0("p\u0302 = NA with 95% CI = [NA, NA]")
      }
      tkconfigure(SEci.lbl,textvariable=dynLbl$seci)
      k.lbl<-tklabel(hSEframe,text="Factor by which searcher\nefficiency changes with\neach search (k) "); k.edit<-tkentry(hSEframe,width=6,textvariable=tkvars$k,justify='right',bg='white')
      tkgrid(SEn.lbl,SEn.edit,sticky='e')
      tkgrid(SEx.lbl,SEx.edit,sticky='e')
      tkgrid(SEci.lbl,columnspan=2,sticky='e')
      tkgrid(k.lbl,k.edit,sticky='e')
      fSEframe<-tkframe(SEFrame)
      .Rvar$pkdat<-.Rvar$pkdatPrevious
      capture.output({
        .Rvar$pkjags <- rjags::jags.model(
          textConnection(pkmod),
          data = .Rvar$pkdat,
          inits = with(.Rvar$pkdat,list(p = X[1]/M[1],k = max(min((X[2]/M[2])/(X[1]/M[1]),.99),.01)))
        )
        update(.Rvar$pkjags, 2000)
        .Rvar$tmppk<-rjags::coda.samples(.Rvar$pkjags, variable.names=c('p','k'), n.iter=2000)[[1]][,2:1]
        },
        file = paste0(.Rvar$datadir, '/NULL')
      )
      pkstat <- list(phat = mean(.Rvar$tmppk[,1]), CIp = quantile(.Rvar$tmppk[,1],c(0.025, 0.975)), khat = mean(.Rvar$tmppk[,2]), CIk = quantile(.Rvar$tmppk[,2],c(0.025, 0.975)), r = cor(.Rvar$tmppk[,1:2])[2])
      .Rvar$pkres<-array(dim=c(dim(.Rvar$tmppk)[1],.Rvar$pkdat$n+2))
      .Rvar$pkstat<-pkstat
      .Rvar$pkres[,1:2]<-.Rvar$tmppk
      .Rvar$pkres[,3]<-.Rvar$tmppk[,1]
      for (i in 4:(.Rvar$pkdat$n+2)){
        .Rvar$pkres[,i]<-.Rvar$pkres[,i-1]*.Rvar$pkres[,2]
      }
      tclvalue(dynLbl$pkci) <- paste0("95% CIs: p \u2208 [", round(pkstat$CIp[1],3), ", ", round(pkstat$CIp[2],3),"], ",
        "k \u2208 [", round(pkstat$CIk[1],3), ", ", round(pkstat$CIk[2],3),"]")
      pkci.lbl<<-tklabel(fSEframe, textvariable = dynLbl$pkci)
      tclvalue(dynLbl$pkr)<<-paste0("p\u0302 = ",round(.Rvar$pkdat$X[1]/.Rvar$pkdat$M[1],3), ", k\u0302 = ", round(.Rvar$pkstat$khat,3))
      pkr.lbl<<-tklabel(fSEframe, textvariable = dynLbl$pkr)
      pkEditsy<-tkbutton(fSEframe, text = "Edit", command = function() {
        pkdatEnter(.Rvar$pkdat)
      })
      pkViewsy<-tkbutton(fSEframe, text = "View", command = function(){
        junk<-.Rvar$pkdat
        npk<-dim(.Rvar$pkres)[2]-2
        if (length(dev.list()) > 0 && !prod(par()$fig==c(0,1,0,1))) {
           wd<-6.5*.Rvar$charSizeAdjust; ht<-7*.Rvar$charSizeAdjust
          if (.Rvar$platform == 'windows') windows.options(width = wd, height = ht)
          if (.Rvar$platform == 'mac') quartz.options(width = wd, height = ht)
          if (.Rvar$platform == 'linux') X11.options(width = wd, height = ht)
           dev.new(noRStudioGD = T)
        }
        par(mar=c(4, 4, 2.5, 1.5),mgp=c(2,.7,0))
        plot(1:npk,rep(.5,npk), type='n',
          xlim = c(.75, npk+.25),
          ylim=range(apply(.Rvar$pkres[,-(1:2)], F = quantile, M = 2, probs=c(.975, 0.025))),
          xlab = "Search Occasion (i)",
          ylab = expression(Searcher ~Efficiency~ (p[~i]))
        )
        if(.Rvar$platform == "windows") bringToTop()
        lines(1:npk, apply(.Rvar$pkres[,-(1:2)], F = mean, M = 2),lwd=2)
        lines(1:npk, apply(.Rvar$pkres[,-(1:2)], F = quantile, M = 2, probs=.975),lty=3, col=2)
        lines(1:npk, apply(.Rvar$pkres[,-(1:2)], F = quantile, M = 2, probs=.025),lty=3, col=2)
        points(1:npk, junk$X/junk$M, pch = 18, cex = 3*.Rvar$charSizeAdjust*sqrt(junk$M)/sqrt(max(junk$M)), col = colors()[125])
        pklabs<-list()
        for (i in 1:npk) pklabs[[i]]<-paste0(junk$X[i],"/",junk$M[i])
        text(1:npk, junk$X/junk$M, lab = pklabs, col=colors()[125], adj = -.25*c(1,1), cex = .75*.Rvar$charSizeAdjust)
        legend(x='topright', legend = c('Empirical X[i]/M[i]', 'Fitted values', '95% credibility bands'),
          lty=c(NA, 1, 3), lwd=c(NA, 2,1), pch = c(18, NA, NA), pt.cex = 1.5, col = c(colors()[125], 1, 2)
        )
        lhs<-par('usr')[1]+diff(par('usr')[1:2])*0.015
        mtext(text = paste0(
          "95% CIs: p \u2208 [", round(quantile(.Rvar$pkres[,1],0.025),3),", ", round(quantile(.Rvar$pkres[,1],0.975),3),"], ",
          "k \u2208 [", round(quantile(.Rvar$pkres[,2],0.025),3), ", ", round(quantile(.Rvar$pkres[,2],0.975),3),"]"),
          side = 1, line = -2.5, at=lhs, adj=0, cex = 0.8*.Rvar$charSizeAdjust)
#        mtext(text = bquote(hat(r)["p,k"]== .(round(cor(.Rvar$pkres[,1],.Rvar$pkres[,2]),4))), side = 1, line = -1.25, at=lhs, adj=0, cex = 0.8*.Rvar$charSizeAdjust)
        mtext(text = paste0("p\u0302 = ",round(.Rvar$pkdat$X[1]/.Rvar$pkdat$M[1],3), ", k\u0302 = ", round(.Rvar$pkstat$khat,3)), side = 1, line = -1.25, at=lhs, adj=0, cex = 0.8*.Rvar$charSizeAdjust)
        title(expression(Estimation~of~p~and~k:~~p[i]==p*k^"i - 1"))
      })
      tkgrid(pkci.lbl, columnspan=3)
      tkgrid(pkr.lbl, pkViewsy, pkEditsy)

      if (tclvalue(tkvars$SEopt)=='h'){
        for (obj in unlist(strsplit(tclvalue(tkwinfo("children",fSEframe)), ' '))) tkconfigure(obj,state='disabled')
        for (obj in unlist(strsplit(tclvalue(tkwinfo("children",hSEframe)), ' '))) tkconfigure(obj,state='normal')
      } else {
        for (obj in unlist(strsplit(tclvalue(tkwinfo("children",hSEframe)), ' '))) tkconfigure(obj,state='disabled')
        for (obj in unlist(strsplit(tclvalue(tkwinfo("children",fSEframe)), ' '))) tkconfigure(obj,state='normal')
      }

      tkgrid(fSE, sticky='w')
      tkgrid(fSEframe, pady=c(0,10))
      tkgrid(hSE, sticky='w', pady=c(10,0))
      tkgrid(hSEframe)
      ## sy1: gFrame : search dates
      start.lbl<-tklabel(datesFrame,text="Start of monitoring\n   (yyyy-mm-dd)",anchor='w')
      start.edit<-tkentry(datesFrame,width=10,textvariable=tkvars$firstsearch,justify='left',bg='white')
#      tkconfigure(start.edit, state = ifelse(partial & sydat$arrfun == "Uniform", "disabled", "normal"))
      formulaRadio<-tkradiobutton(datesFrame,variable=tkvars$samtype,text="Formula",value="Formula", command = function(){
        tkconfigure(I.lbl,state="normal");    tkconfigure(I.edit,state="normal")
        tkconfigure(nsearch.lbl,state="normal"); tkconfigure(nsearch.edit,state="normal")
        tkconfigure(datesEditView.button,state='disabled')
        tkconfigure(cuslab, state='disabled')
        Ir<<-round(toR(tkvars$Isam), 1)
        updateCPhandLabel()
        updateCPfieldLabel(.Rvar$CPdata)
      })
      customdatesRadio<-tkradiobutton(datesFrame,variable=tkvars$samtype,text="Custom",value="Custom", command = function(){
        tkconfigure(I.lbl,state="disabled");    tkconfigure(I.edit,state="disabled")
        tkconfigure(nsearch.lbl,state="disabled"); tkconfigure(nsearch.edit,state="disabled")
        tkconfigure(datesEditView.button,state='normal')
        tkconfigure(cuslab, state='normal')
        dztmp<-toR(tkvars$days)
        Ir<<-round(max(dztmp)/(length(dztmp)-1), 1)
        updateCPhandLabel()
        updateCPfieldLabel(.Rvar$CPdata)
      })
      datesEditView.button<-tkbutton(datesFrame,text="Edit/View",command=function() cd_form())
      I.lbl<-tklabel(datesFrame,text="    Search interval (I)",anchor='w');  I.edit<-tkentry(datesFrame,width=4,textvariable=tkvars$Isam,justify='right',bg='white')
      nsearch.lbl<-tklabel(datesFrame,text="    Number of searches",anchor='w');  nsearch.edit<-tkentry(datesFrame,width=4,textvariable=tkvars$nsearch,justify='right',bg='white')
      cuslab<<-tklabel(datesFrame, text = paste0("    span = ", max(toR(tkvars$days)), ", I (mean) = ", round(max(toR(tkvars$days))/(length(toR(tkvars$days))-1),1)))
      tkgrid(start.lbl, start.edit, columnspan=2,sticky='w',pady=c(0,8))
      tkgrid(formulaRadio,column=0,rowspan=2,sticky='nw',columnspan=2)
      tkgrid(I.lbl,I.edit,sticky='w',columnspan=2)
      tkgrid(nsearch.lbl,nsearch.edit,sticky='w',columnspan=2);
      tkgrid(customdatesRadio,datesEditView.button,pady=c(10, 4), sticky='w')
      tkgrid(cuslab, row = 6, sticky='new', columnspan = 2, pady=c(0, 10))
      if (tclvalue(tkvars$samtype)=="Formula"){
        tkconfigure(cuslab, state='disabled')
        tkconfigure(datesEditView.button,state="disabled")
 #        Ir<-.Rvar$singleYearPrevious$Isam
         Ir<<-sydat$Isam
      }
      if (tclvalue(tkvars$samtype)=="Custom"){
        tkconfigure(cuslab, state='normal')
#        tkconfigure(start.lbl,state="disabled"); tkconfigure(start.edit,state="disabled")
        tkconfigure(I.lbl,state="disabled"); tkconfigure(I.edit,state="disabled")
        tkconfigure(nsearch.lbl,state="disabled"); tkconfigure(nsearch.edit,state="disabled")
#        Ir<-getmode(diff(.Rvar$singleYearPrevious$days))
        Ir<<-round(max(sydat$days)/(length(sydat$days)-1), 1)
      }
      ## sy1: gFrame : persistenceFrame
      persistenceFrame<-ttklabelframe(gFrame,text="   Persistence Distribution",padding=10)
      # sy1: gFrame : persistenceFrame: pdpargFrame
      # prepare summary data for display in field persistence labels...
      nsim<-1000
      CPab<-array(dim=c(nsim, 2))
      CPr<-rCPgab(.Rvar$CPdata$persistence_distn, a = .Rvar$CPdata$pda, b = .Rvar$CPdata$pdb, Ir = Ir)
      if (.Rvar$CPdataPrevious$persistence_distn == "Exponential"){
        hCP.ci<-exp(qnorm(c(0.025, 0.975), mean = .Rvar$CPdata$mod.e$coef, sd = sqrt(.Rvar$CPdata$mod.e$var[1])))
        rlwr<-rCPgab("Exponential", a=1/exp(.Rvar$CPdata$mod.e$coef[1]),b=exp(.Rvar$CPdata$mod.e$coef[1]+qt(.025,.Rvar$CPdata$mod.e$df.res)*sqrt(.Rvar$CPdata$mod.e$var[1])),Ir=Ir)[2]
        rupr<-rCPgab("Exponential", a=1/exp(.Rvar$CPdata$mod.e$coef[1]),b=exp(.Rvar$CPdata$mod.e$coef[1]+qt(.975,.Rvar$CPdata$mod.e$df.res)*sqrt(.Rvar$CPdata$mod.e$var[1])),Ir=Ir)[2]
        hr.ci<-c(rlwr, rupr)
        bci <- exp(.Rvar$CPdata$mod.e$coef+qt(c(0.025, 0.975),.Rvar$CPdata$mod.e$df.res)*sqrt(.Rvar$CPdata$mod.e$var[1]))
        hshape<-exp(1/.Rvar$CPdata$mod.e$coef)
        hscale<-exp(.Rvar$CPdata$mod.e$coef)
      } else if (.Rvar$CPdata$persistence_distn == "Weibull"){
        CPparms<-MASS::mvrnorm(nsim,c(.Rvar$CPdata$mod.w$coef[1],log(.Rvar$CPdata$mod.w$scale)),.Rvar$CPdata$mod.w$var)
        CPab[,1]<-1/exp(CPparms[,2]) #shape
        CPab[,2]<-exp(CPparms[,1]) # scale
        i0<-which(CPab[,1]<=0 | (Isam/CPab[,2])^CPab[,1]<1e-320)
        while(length(i0)>0){
          CPparms[i0,]<-MASS::mvrnorm(length(i0),c(.Rvar$CPdata$mod.w$coef[1],log(.Rvar$CPdata$mod.w$scale)), .Rvar$CPdata$mod.w$var)
          CPab[i0,1]<-1/exp(CPparms[i0,2])
          CPab[i0,2]<-exp(CPparms[i0,1])
          i0<-which(CPab[,1]<=0 | (Isam/CPab[,2])^CPab[,1]<1e-320)
        }
        hCP.ci<-quantile(CPab[,2]*gamma(1+1/CPab[,1]), c(0.025, 0.975))
        bci <- exp(.Rvar$CPdata$mod.w$coef[1]+qt(c(0.025, 0.975),.Rvar$CPdata$mod.w$df.res)*sqrt(.Rvar$CPdata$mod.w$var[1]))
        hshape<-1/.Rvar$CPdata$mod.w$scale
        hscale<-exp(.Rvar$CPdata$mod.w$coef[1])
      } else if (.Rvar$CPdata$persistence_distn == "Log-Logistic"){
        CPparms<-MASS::mvrnorm(nsim,c(.Rvar$CPdata$mod.ll$coef[1],log(.Rvar$CPdata$mod.ll$scale)),.Rvar$CPdata$mod.ll$var)
        CPab[,1]<-1/exp(CPparms[,2]) #shape
        CPab[,2]<-exp(CPparms[,1]) # scale
        i0<-which(CPab[,1]<=0)
        while(length(i0)>0){
          CPparms[i0,]<-MASS::mvrnorm(length(i0),c(.Rvar$CPdata$mod.ll$coef[1],log(.Rvar$CPdata$mod.ll$scale)), .Rvar$CPdata$mod.ll$var)
          CPab[i0,1]<-1/exp(CPparms[i0,2])
          CPab[i0,2]<-exp(CPparms[i0,1])
          i0<-which(CPab[,1]<=0)
        }
        llmcp<-numeric(nsim)
        llmcp[CPab[,1]<=1]<-Inf
        ind<-which(CPab[,1]>1)
        llmcp[ind]<-(CPab[ind,2]*pi/CPab[ind,1])/sin(pi/CPab[ind,1])
        hCP.ci<-quantile(llmcp, c(0.025, 0.975))
        bci <- exp(.Rvar$CPdata$mod.ll$coef[1]+qt(c(0.025, 0.975),.Rvar$CPdata$mod.ll$df.res)*sqrt(.Rvar$CPdata$mod.ll$var[1]))
        hshape<-1/.Rvar$CPdata$mod.ll$scale
        hscale<-exp(.Rvar$CPdata$mod.ll$coef[1])
      } else if (.Rvar$CPdata$persistence_distn == "Lognormal"){
        CPparms<-MASS::mvrnorm(nsim,c(.Rvar$CPdata$mod.ln$coef[1], log(.Rvar$CPdata$mod.ln$scale)),.Rvar$CPdata$mod.ln$var)
        CPab[,1]<- exp(CPparms[,2])^2 #shape
        CPab[,2]<- CPparms[,1] # scale
        hCP.ci<-quantile(exp(CPab[,2]+CPab[,1]/2),c(0.025, 0.975))
        bci <- .Rvar$CPdata$mod.ln$coef[1] + qt(c(0.025, 0.975),.Rvar$CPdata$mod.ln$df.res)* sqrt(.Rvar$CPdata$mod.ln$var[1])
        hshape<-.Rvar$CPdata$mod.ln$scale^2
        hscale<-.Rvar$CPdata$mod.ln$coef[1]
      }
      if (.Rvar$CPdata$persistence_distn != "Exponential"){
        rr <- as.vector(ppersist(.Rvar$CPdata$persistence_distn,t_arrive0 = 0, t_arrive1=Ir, t_search = Ir, pda = CPab[,1], pdb = CPab[,2]))
        rr[is.na(rr)] <- max(rr[!is.na(rr)])
        hr.ci<-quantile(rr, probs = c(0.025, 0.975))
      }

      pdfieldFrame<-tkframe(persistenceFrame)
      distr.txt <<- paste0("Distribution: ", .Rvar$CPdata$persistence_distn,
        " with shape (\u03b1) = ", signif(hshape,4)," and scale (\u03b2) = ", signif(hscale, 4))
      tclvalue(dynLbl$distr) <-distr.txt
      distr.lbl<-tklabel(pdfieldFrame, textvariable = dynLbl$distr)
      distrcpr.lbl<<-tklabel(pdfieldFrame, textvariable = dynLbl$distrcpr)
      updateCPfieldLabel(.Rvar$CPdata)
      pdfieldEdit<-tkbutton(pdfieldFrame, text = "View/Edit", command = fcp_form)
      pdtypeField<-tkradiobutton(pdfieldFrame,
        text = "Use field trials to estimate parameters",
        variable=tkvars$perstype,
        value = "field",
        command=function(){
          # normalize the 'field' labels
          tkconfigure(distr.lbl, state='normal')
          tkconfigure(distrcpr.lbl, state='normal')
          tkconfigure(pdfieldEdit, state='normal')
          # disable the 'hand' labels (and edit boxes)
          tkconfigure(pda.lbl,state='disabled'); tkconfigure(pda.edit, state='disabled')
          tkconfigure(pdb.lbl,state='disabled'); tkconfigure(pdb.edit, state='disabled')
          tkconfigure(blwr.lbl,state='disabled'); tkconfigure(blwr.edit,state='disabled')
          tkconfigure(bupr.lbl,state='disabled'); tkconfigure(bupr.edit,state='disabled')
          tkconfigure(rCP.lbl,state='disabled')
          tkconfigure(pd.lbox, state='disabled', bg='gray93')
          tkconfigure(pdView, state = 'disabled')
          .Rvar$CPdat$persistence_distn <- tclvalue(tkvars$persistence_distn) #tkvars$persistence_distn changes in fcp_form
          # open the fcp form
          # source("fcp_form.Rea")
      })
      tkgrid(pdtypeField, sticky='w')
      tkgrid(pdfieldEdit, row = 0, column = 1, sticky='w')
      tkgrid(distr.lbl, sticky='w', padx=15, columnspan = 2)
      tkgrid(distrcpr.lbl, sticky='w', padx=15, columnspan = 2)
      pdparmFrame<-ttklabelframe(persistenceFrame,text='Parameters')
      pdtypeHand<<-tkradiobutton(persistenceFrame,
        text = "Enter parameter estimates manually",
        variable = tkvars$perstype,
        value = "hand",
        command=function(){
          # disable the 'field' labels
          tkconfigure(distr.lbl, state='disabled')
          tkconfigure(distrcpr.lbl, state='disabled')
          tkconfigure(pdfieldEdit, state='disabled')
          # normalize the 'hand' labels (and edit boxes)
          tkconfigure(pda.lbl,state='normal');
          ast<-ifelse(tclvalue(tkget(pd.lbox,tkcurselection(pd.lbox)))=="Exponential", "disabled", "normal")
          tkconfigure(pda.edit, state=ast)
          tkconfigure(pdb.lbl,state='normal'); tkconfigure(pdb.edit, state='normal')
          tkconfigure(blwr.lbl,state='normal'); tkconfigure(blwr.edit,state='normal')
          tkconfigure(bupr.lbl,state='normal'); tkconfigure(bupr.edit,state='normal')
          tkconfigure(rCP.lbl,state='normal')
          tkconfigure(pd.lbox, state='normal', bg='white')
          tkconfigure(pdView, state = 'normal')
          persistence_distn <<- tclvalue(tkget(pd.lbox,tkcurselection(pd.lbox)))
      })
      shapescaleFrame<-tkframe(pdparmFrame,width=5)
      pda.lbl<-tklabel(shapescaleFrame,text="shape (\u03b1) "); pda.edit<<-tkentry(pdparmFrame,width=w,textvariable=tkvars$pda,justify='right',bg='white')
      pdb.lbl<-tklabel(shapescaleFrame,text="scale (\u03b2) "); pdb.edit<-tkentry(pdparmFrame,width=w,textvariable=tkvars$pdb,justify='right',bg='white')
      blwr.lbl<-tklabel(shapescaleFrame,text="lwr ");           blwr.edit<-tkentry(pdparmFrame,width=w,textvariable=tkvars$blwr,justify='right',bg='white')
      bupr.lbl<-tklabel(shapescaleFrame,text="upr ");           bupr.edit<-tkentry(pdparmFrame,width=w,textvariable=tkvars$bupr,justify='right',bg='white')
      tkconfigure(pda.lbl,width=10)
      tkconfigure(pdb.lbl,width=10)
      rCP.lbl<<-tklabel(pdparmFrame, width=40)
      if (pdaok & pdbok & ifelse(toR(tkvars$samtype)=="Formula", Ichk(toR(tkvars$Isam)),Ichk(getmode(diff(toR(tkvars$days)))))){
        pda<-as.numeric(tclvalue(tkvars$pda))
        pdb<-as.numeric(tclvalue(tkvars$pdb))
        rCP<-rCPgab(persistence_distn, pda, pdb, Ir)
        rCPhtext1<<-paste0("r = ",signif(rCP[2],3), " for Ir = ", Ir)

        tkconfigure(rCP.lbl,text=paste0(rCPhtext1, rCPhtext2))

        meanCP<<-rCP[1]; rhat<<-rCP[2]
        if (bminok & bmaxok){
          blwr<-as.numeric(tclvalue(tkvars$blwr))
          bupr<-as.numeric(tclvalue(tkvars$bupr))
          rCPlwr<-rCPgab(persistence_distn, pda, blwr, Ir)
          rCPupr<-rCPgab(persistence_distn, pda, bupr, Ir)
          rCPhtext2<<-paste0(", with 95% CI: r \u2208 [",signif(rCPlwr[2],3),", ",signif(rCPupr[2],3),"]")
          tkconfigure(rCP.lbl,text=paste0(rCPhtext1, rCPhtext2))
        }
      } else {
        rCPhtext1<<-paste0("r = ",NA, " for Ir = ",NA)
        meanCP<<-NA; rhat<<-NA
      }
      if (! (pdaok & pdbok & bminok & bmaxok)) {
        rCPhtext2<<-paste0(", with 95% CI: r \u2208 [",NA,", ",NA,"]")
        tkconfigure(rCP.lbl,text=paste0(rCPhtext1, rCPhtext2))
      }
      if (persistence_distn=="Exponential"){
        tkconfigure(pda.edit,state='disabled')
      }
      # sy1: gFrame : persistenceFrame: pdlboxFrame
      pdlboxFrame<-tkframe(persistenceFrame)
      pdView<-tkbutton(persistenceFrame,text="View",width=buttwid,command=function() {
        if (pdaok & pdbok & bminok & bmaxok & Irok){
          if (length(dev.list()) > 0 && !prod(par()$fig==c(0,1,0,1))) {
            if (.Rvar$platform == 'windows') windows.options(width = 7*.Rvar$charSizeAdjust, height = 7*.Rvar$charSizeAdjust)
            if (.Rvar$platform == 'mac') quartz.options(width = 7*.Rvar$charSizeAdjust, height = 7*.Rvar$charSizeAdjust)
            if (.Rvar$platform == 'linux') X11.options(width = 7*.Rvar$charSizeAdjust, height = 7*.Rvar$charSizeAdjust)
            dev.new(noRStudioGD = FALSE)
          }
          plotPersDist(persistence_distn,
            Ir,#as.numeric(tclvalue(tkvars$Ir)),
            as.numeric(tclvalue(tkvars$pda)),
            as.numeric(tclvalue(tkvars$pdb)),
            as.numeric(tclvalue(tkvars$blwr)),
            as.numeric(tclvalue(tkvars$bupr)))
        } else {
          tkmessageBox(message="Error in data")
          return(F)
        }
      })
      pd.lbox<<-tklistbox(pdlboxFrame,height=4,width=12,selectmode="single",exportselection="0",bg=listboxbgclr)
      for (i in 1:length(pdnames)) tkinsert(pd.lbox,"end", pdnames[i])
      tkselection.set(pd.lbox,which(pdnames==persistence_distn)-1)
      tkactivate(pd.lbox,which(pdnames==persistence_distn)-1)
      if (tclvalue(tkvars$perstype)=='field'){
          tkconfigure(distr.lbl, state='normal')
          tkconfigure(distrcpr.lbl, state='normal')
          tkconfigure(pdfieldEdit, state='normal')
          # disable the 'hand' labels (and edit boxes)
          tkconfigure(pda.lbl,state='disabled'); tkconfigure(pda.edit, state='disabled')
          tkconfigure(pdb.lbl,state='disabled'); tkconfigure(pdb.edit, state='disabled')
          tkconfigure(blwr.lbl,state='disabled'); tkconfigure(blwr.edit,state='disabled')
          tkconfigure(bupr.lbl,state='disabled'); tkconfigure(bupr.edit,state='disabled')
          tkconfigure(rCP.lbl,state='disabled')
          tkconfigure(pd.lbox, state='disabled', bg='gray93')
      } else {
          tkconfigure(distr.lbl, state='disabled')
          tkconfigure(distrcpr.lbl, state='disabled')
          tkconfigure(pdfieldEdit, state='disabled')
          # disable the 'hand' labels (and edit boxes)
          tkconfigure(pda.lbl,state='normal')
          ast<-ifelse(tclvalue(tkget(pd.lbox,tkcurselection(pd.lbox)))=="Exponential", "disabled", "normal")
          tkconfigure(pda.edit, state=ast)
          tkconfigure(pdb.lbl,state='normal'); tkconfigure(pdb.edit, state='normal')
          tkconfigure(blwr.lbl,state='normal'); tkconfigure(blwr.edit,state='normal')
          tkconfigure(bupr.lbl,state='normal'); tkconfigure(bupr.edit,state='normal')
          tkconfigure(rCP.lbl, state='normal')
          tkconfigure(pd.lbox, state='normal', bg='white')
      }

      tkgrid(pda.lbl,pda.edit)
      tkgrid(pdb.lbl,pdb.edit,blwr.lbl,blwr.edit,bupr.lbl,bupr.edit)
      tkgrid(shapescaleFrame)
      tkgrid(rCP.lbl,sticky='w',columnspan=6)
      tkgrid(pdfieldFrame, sticky='w', columnspan=2, pady=c(0, 10))
      tkgrid(pdtypeHand, sticky='w', columnspan=2, pady = c(10,0))
      tkgrid(pdView, column=1, row = 1)
      tkgrid(pd.lbox, pady=7)
      tkgrid(pdlboxFrame, padx = 10, column = 0, sticky = 'n', pady = 15)
      tkgrid(pdparmFrame, padx = 10, column = 1, row = 2, sticky = 'n', pady = c(15, 5))
      if(!basicMode) tkgrid(arrfFrame,row=0,column=0, sticky = 'nw', pady = c(10, 5))
      tkgrid(datesFrame,row=0,column=0,sticky='nw',pady=10, ipadx = 10)
      syCalcg<-tkbutton(gFrame,text='Estimate g',width=12, command=function(){SYcalcg(writeg=T)})
      tkgrid(syCalcg)
      tkgrid(aFrame,row = 7,column=0,sticky='w',pady=c(5, 10),columnspan=3)
#      tkgrid.propagate(SEFrame, 0)
      tkconfigure(SEFrame, width = 300, height = 280)
      tkgrid(SEFrame,row=0,column=1, rowspan = 3, sticky = 'n', pady = 10)
      tkgrid(persistenceFrame, row = 0, column = 2, rowspan=3, sticky = 'nw', pady=10)
      tkgrid(gFrame, columnspan = 2, pady=15)
        ###  sy2 = second column on single year frame
      #sy2<-tkframe(syModule)
      ## sy2: fatality frame
      Mframe<-ttklabelframe(sy1,text="   Fatality estimation (M, \u03bb)",padding=10);
      ## sy2: observations (X) frame
      Xlbl<-tklabel(Mframe,text="Carcass Count (X) ")
      X.edit<-tkentry(Mframe,width=5,textvariable=tkvars$X,justify='right',bg='white') # define input box for X
      crlev.lbl<-tklabel(Mframe,text="Credibility level (1 - \u3b1) ");  crlev.edit<-tkentry(Mframe,width=5,textvariable=tkvars$crlev,justify='right',bg='white')
      MoptFrame<-tkframe(Mframe)
      syCalcM<-tkbutton(MoptFrame,text='Estimate M',width=12, command=function(){SYcalcPost(writepost=T)})
      syCalcL<-tkbutton(Mframe, text="Estimate \u03bb",width=12, command=function(){SYcalcLambda()})
      sides <<- tclVar(1)
      oneside <- tkradiobutton(MoptFrame, text = "One-sided CI (M*)", variable = sides, value = 1)
      twoside <- tkradiobutton(MoptFrame, text = "Two-sided CI", variable = sides, value = 2)
      tkgrid(Xlbl,X.edit)
      tkgrid(crlev.lbl,crlev.edit)
      tkgrid(syCalcM)
      tkgrid(oneside, row=0, column=1)
      tkgrid(twoside, row=0, column=2)
      tkgrid(MoptFrame, row=0, column=2, padx=10)
      tkgrid(syCalcL, row=1, column=2, padx=10, sticky='w')

      if (partial == F){
        tkgrid(oneside, sticky='w', padx=c(30,0))
        tkgrid(twoside, sticky='w', padx=c(30,0))
      }
      if (!basicMode){
        ## sy2: Mframe: priorFrame
        if (prior_f=="Custom")   prior_M<<-toR(tkvars$custom.prior)
        if (prior_f=="Objective") prior_M<<-toR(tkvars$objective.prior)
        priorFrame<-ttklabelframe(Mframe,text="Prior Distribution")
        prior.lbox<<-tklistbox(priorFrame,height=3,width=12,selectmode="single",exportselection="0",bg=listboxbgclr);
        ttPriortmp1<-tkframe(priorFrame)
        priorEditButton<<-tkbutton(ttPriortmp1,text="Edit",width=buttwid,command=cpr_form)

        priorViewButton<<-tkbutton(ttPriortmp1,text="View",width=buttwid,command=function(){
          # this button is disabled unless prior_M, prior_f
        #  prior_M<-toR(tkvars$prior_M) prior_M is directly handled without tk mediator
            plotPrior(self$prior_M, self$prior_f)
        })
        tkgrid(priorEditButton); tkgrid(priorViewButton)
        tkgrid(prior.lbox,ttPriortmp1,sticky='w')
        for (i in 1:length(prnames)) tkinsert(prior.lbox,i-1,prnames[i])
        if (prior_f=="Objective"){
          tmplbl<-"          ...           "
        } else if (prior_f=="Custom") {
          tmplbl<-paste0("95th percentile = ", min(which(cumsum(toR(tkvars$custom.prior))>=0.95))-1)
        }
        tkselection.set(prior.lbox,which(prnames==prior_f)-1)
        tkactivate(prior.lbox,which(prnames==prior_f)-1)
        prior.lbl<<-tklabel(priorFrame,text=tmplbl)
        tkgrid(prior.lbl,columnspan=2)
        if (prior_f!="Objective") {
          tkconfigure(priorEditButton, state='normal')
        } else {
          tkconfigure(priorEditButton, state='disabled')
        }
        tkgrid(priorFrame, rowspan=4, column=2, row = 0, sticky='nw', padx=15)
      }
      #tkgrid(prandarFrame,columnspan=2,sticky='n')
      ## sy2: syButtonsFrame
      #syButtonsFrame<-tkframe(gFrame)
      syButtonsFrame<-tkframe(sy1)
      # buttons
      syClose<-tkbutton(syButtonsFrame, text = "Close", width =12, command = function(){
        graphics.off()
        if (!partial) {
          sysave(.Rvar$singleYearPrevious, .Rvar$CPdataPrevious, .Rvar$pkdatPrevious, .Rvar$csvpath)
          tkwm.deiconify(.Rvar$EoA)
        } else {
          sysave(.Rvar$syscPrevious, .Rvar$CPdataPrevious, .Rvar$pkdatPrevious, .Rvar$csvpath, T)
        }
        tkdestroy(syModule)
      })
      tkgrid(syClose)
      tkgrid(Mframe, sticky='nw')
      tkgrid(syButtonsFrame, column = 1, row = 1)
      tkgrid(sy1)

      # singleYearForm actions
      #tkbind(arrUniformRadio,"<Button-1>", function() {updatearrrad("Uniform")})
      #tkbind(arrCompoundRadio,"<Button-1>", function()  {updatearrrad("Compound")})
      if (!partial){
        tkbind(prior.lbox,"<<ListboxSelect>>", function() {
           v<-tclvalue(tkget(prior.lbox,tkcurselection(prior.lbox)))
            prior_f<<-v
            if (v=='Objective') {
              tkconfigure(prior.lbl,text="          ...           ")
              tkconfigure(priorEditButton, state='disabled')
              tkconfigure(priorViewButton, state='normal')
              prior_M<<-toR(tkvars$objective.prior)
              prok<<-T
            }
            if (v=='Custom') {
              custom.prior<-toR(tkvars$custom.prior)
              if (length(custom.prior)<=1){
                pr95<-NA
              } else {
                pr95<-min(which(cumsum(custom.prior)>=0.95))-1
              }
              tkconfigure(prior.lbl,text=paste0("95th percentile = ", pr95))
              tkconfigure(priorEditButton, state='normal')
              buttstate<-ifelse(sum(is.na(suppressWarnings(as.numeric(custom.prior)))), "disabled", "normal")
              tkconfigure(priorViewButton,state=buttstate)
              prior_M<<-custom.prior
              if (buttstate=='disabled') prok<<-F else prok<<-T
            }
          })
      }
      tkbind(pd.lbox,"<<ListboxSelect>>", function(){
        persistence_distn<<-tclvalue(tkget(pd.lbox,tkcurselection(pd.lbox)))
        v<-persistence_distn
        if (v=='Exponential') {
          tkconfigure(pda.lbl,text='rate');
          pdb.e<-toR(tkvars$pdb.e)
          tclvalue(tkvars$pda)<-signif(1/pdb.e,5); tkconfigure(pda.edit,state='disabled')
          tclvalue(tkvars$pdb)<-tclvalue(tkvars$pdb.e)
          tclvalue(tkvars$blwr)<-tclvalue(tkvars$blwr.e)
          tclvalue(tkvars$bupr)<-tclvalue(tkvars$bupr.e)
        } else if (v=='Weibull' |  v=='Log-Logistic' | v=='Lognormal') {
          tkconfigure(pda.lbl,text='shape (\u03b1) '); tkconfigure(pda.edit,state='normal')
          if (v=='Weibull'){
            tclvalue(tkvars$pda)<-tclvalue(tkvars$pda.w)
            tclvalue(tkvars$pdb)<-tclvalue(tkvars$pdb.w)
            tclvalue(tkvars$blwr)<-tclvalue(tkvars$blwr.w)
            tclvalue(tkvars$bupr)<-tclvalue(tkvars$bupr.w)
          } else if (v=='Log-Logistic') {
            tclvalue(tkvars$pda)<-tclvalue(tkvars$pda.ll)
            tclvalue(tkvars$pdb)<-tclvalue(tkvars$pdb.ll)
            tclvalue(tkvars$blwr)<-tclvalue(tkvars$blwr.ll)
            tclvalue(tkvars$bupr)<-tclvalue(tkvars$bupr.ll)
          } else if (v=='Lognormal') {
            tclvalue(tkvars$pda)<-tclvalue(tkvars$pda.ln)
            tclvalue(tkvars$pdb)<-tclvalue(tkvars$pdb.ln)
            tclvalue(tkvars$blwr)<-tclvalue(tkvars$blwr.ln)
            tclvalue(tkvars$bupr)<-tclvalue(tkvars$bupr.ln)
          }
        }
        pda<-toR(tkvars$pda)
        pdb<-toR(tkvars$pdb)
      #    Ir<-as.numeric(tclvalue(tkvars$Ir))
        blwr<-toR(tkvars$blwr)
        bupr<-toR(tkvars$bupr)
        rCP<-rCPgab(v, pda, pdb, Ir)
        meanCP<<-rCP[1]; rhat<<-rCP[2]
        rCPlwr<-rCPgab(v, pda, blwr, Ir) # not quite...need to find values for new distribution that looks like the old
        rCPupr<-rCPgab(v, pda, bupr, Ir)
        rCPhtext1<<-paste0("r = ",signif(rCP[2],3), " for Ir = ", Ir)
        rCPhtext2<<-paste0(", with 95% CI: r \u2208 [",signif(rCPlwr[2],3),", ",signif(rCPupr[2],3),"]")
        tkconfigure(rCP.lbl,text=paste0(rCPhtext1, rCPhtext2))
        # when listbox selection is changed, then the parameter set that is loaded is proper, so...
        #  -- all the backgrounds are whited
        #  -- the OK and leg flags are set to true
        #  -- the "view" button is set to true
        #  -- the full suite of flags is checked for possible activation of the 'calculate' buttons
        tkconfigure(pda.edit, bg='white'); pdaok<<-T
        tkconfigure(pdb.edit, bg='white'); pdbok<<-T
        tkconfigure(blwr.edit,bg='white'); bminok<<-T; bminleg<<-T
        tkconfigure(bupr.edit,bg='white'); bmaxok<<-T; bmaxleg<<-T
      })
      # live error-checking
      tkbind(X.edit,"<KeyRelease>", function() {
        X<-try(suppressWarnings(as.numeric(tclvalue(tkvars$X))),silent=T)
        if (class(X) == "try-error") Xok<<-F else Xok<<-Xchk(X)
        if (Xok) tkconfigure(X.edit, bg='white') else tkconfigure(X.edit, bg=colors()[652])
      })
      tkbind(a.edit,"<KeyRelease>", function() {
        a<-try(suppressWarnings(as.numeric(tclvalue(tkvars$a))),silent=T)
        if (class(a) == "try-error") aok<<-F else aok<<-achk(a)
        if (aok) tkconfigure(a.edit, bg='white') else tkconfigure(a.edit, bg=colors()[652])

      })
      tkbind(v.edit,"<KeyRelease>", function() {
        v<-try(suppressWarnings(as.numeric(tclvalue(tkvars$v))),silent=T)
        if (class(v) == "try-error") vok<<-F else vok<<-vchk(v)
        if (vok) tkconfigure(v.edit, bg='white') else tkconfigure(v.edit, bg=colors()[652])
      })
      tkbind(SEn.edit,"<KeyRelease>", function() {
        SEn<-try(suppressWarnings(as.numeric(tclvalue(tkvars$SEn))),silent=T)
        if (class(SEn) == "try-error") SEnok<<-F else SEnok<<-SEnchk(SEn)
        if (SEnok){
          tkconfigure(SEn.edit, bg='white')
        } else {
          tkconfigure(SEn.edit, bg=colors()[652])
          tclvalue(dynLbl$seci)<-"p\u0302 = NA, with 95% CI = [NA, NA] "
        }
        if (SExok & SEnok){
          SEx<-toR(tkvars$SEx)
          fba<<-SEx+.5; fbb<<-SEn-SEx+.5
          tclvalue(dynLbl$seci)<-paste0("p\u0302 = ", round(SEx/SEn,3),", with 95% CI = [" ,round(qbeta(0.025, fba, fbb),3),", ",round(qbeta(0.975, fba, fbb),3),"] ")
        } else {
          tclvalue(dynLbl$seci)<-"p\u0302 = NA, with 95% CI = [NA, NA] "
        }
        if (SExok) tkconfigure(SEx.edit, bg='white') else tkconfigure(SEx.edit, bg=colors()[652])
      })
      tkbind(SEx.edit,"<KeyRelease>", function() {
        SEx<-try(suppressWarnings(as.numeric(tclvalue(tkvars$SEx))),silent=T)
        if (class(SEx) == "try-error") SExok<<-F else SExok<<-SExchk(SEx)
        if (SExok) {
          tkconfigure(SEx.edit, bg='white')
        } else {
          tkconfigure(SEx.edit, bg=colors()[652])
          tclvalue(dynLbl$seci) <- "p\u0302 = NA, with 95% CI = [NA, NA] "
        }
        if (SExok & SEnok){
          SEn<-toR(tkvars$SEn)
          fba<<-SEx+.5; fbb<<-SEn-SEx+.5
          tclvalue(dynLbl$seci)<-paste0("p\u0302 = ", round(SEx/SEn,3),", with 95% CI = [" ,round(qbeta(0.025, fba, fbb),3),", ",round(qbeta(0.975, fba, fbb),3),"] ")
        } else {
          tclvalue(dynLbl$seci)<-"p\u0302 = NA, 95% CI = [NA, NA] "
          return(T)
        }
      })
      tkbind(k.edit,"<KeyRelease>", function() {
        k<-try(suppressWarnings(as.numeric(tclvalue(tkvars$k))),silent=T)
        if (class(k) == "try-error") kok<<-F else kok<<-kchk(k)
        if (kok) tkconfigure(k.edit, bg='white') else tkconfigure(k.edit, bg=colors()[652])
      })
      tkbind(I.edit,"<KeyRelease>", function() {
        Ir<-try(suppressWarnings(as.numeric(tclvalue(tkvars$Isam))),silent=T)
        if (class(Ir) == "try-error")  Iok<<-F else Iok<<-Ichk(Ir)
        if (Iok) {
          Ir<<-Ir
          tkconfigure(I.edit, bg='white')
          if (!basicMode && nsearchok & startok & (tclvalue(tkvars$arrfun)=="Uniform")) tkrplot::tkrreplot(arrfig.mini)
              # if persistence parameters are OK, then label for r/CP should be changed to Iparm, r/CP
          if (pdaok & pdbok){
            pda<-as.numeric(tclvalue(tkvars$pda)); pdb<-as.numeric(tclvalue(tkvars$pdb))
            rCP<-rCPgab(persistence_distn,pda,pdb,Ir)
            rCPhtext1<<-paste0("r = ",signif(rCP[2],3), " for Ir = ",Ir)
            meanCP<<-rCP[1]; rhat<<-rCP[2]
            if (bminok & bmaxok){
              blwr<-as.numeric(tclvalue(tkvars$blwr)); bupr<-as.numeric(tclvalue(tkvars$bupr))
              rCPlwr<-rCPgab(persistence_distn, pda, blwr, Ir)
              rCPupr<-rCPgab(persistence_distn, pda, bupr, Ir)
              rCPhtext2<<-paste0(", with 95% CI: r \u2208 [",signif(rCPlwr[2],3),", ",signif(rCPupr[2],3),"]")
            } else {
              rCPhtext2<<-paste0(", with 95% CI: r \u2208 [",NA,", ",NA,"]")
            }
          } else {
            rCPhtext1<<-paste0("r = ",NA, " for Ir = ",Ir)
            rCPhtext2<<-paste0(", with 95% CI: r \u2208 [",NA,", ",NA,"]")
            meanCP<<-NA; rhat<<-NA

          }
        } else {
          tkconfigure(I.edit, bg=colors()[652])
          rCPhtext1<<-paste0("r = ",NA, " for Ir = ",NA)
          rCPhtext2<<-paste0(", with 95% CI: r \u2208 [",NA,", ",NA,"]")
        }
        tkconfigure(rCP.lbl, text = paste0(rCPhtext1, rCPhtext2))
        updateCPfieldLabel(.Rvar$CPdata)
      })
      tkbind(nsearch.edit,"<KeyRelease>", function() {
        nsearch<-try(suppressWarnings(as.numeric(tclvalue(tkvars$nsearch))),silent=T)
        if (class(nsearch) == "try-error") nsearchok<<-F else nsearchok<<-nsearchchk(nsearch)
        if (nsearchok){
          tkconfigure(nsearch.edit, bg='white')
          if (!basicMode && Iok & startok & (tclvalue(tkvars$arrfun)=="Uniform")) tkrplot::tkrreplot(arrfig.mini)
        } else {
          tkconfigure(nsearch.edit, bg=colors()[652])
        }
      })
      tkbind(start.edit,"<KeyRelease>",function(){
        firstsearch<-tclvalue(tkvars$firstsearch)
        if (startchk(firstsearch)){
          tkconfigure(start.edit, bg='white')
          if (!basicMode && Iok & nsearchok & (tclvalue(tkvars$arrfun)=="Uniform")) tkrplot::tkrreplot(arrfig.mini)
        } else {
          tkconfigure(start.edit, bg=colors()[652])
        }
      })
      tkbind(pda.edit,"<KeyRelease>", function() {
        pda<-try(suppressWarnings(as.numeric(tclvalue(tkvars$pda))), silent=T)
        if (class(pda) == "try-error") pdaok <<- F else pdaok <<- pdachk(pda)
        if (pdaok){
          tkconfigure(pda.edit, bg='white')
          if (pdaok & pdbok & ifelse(toR(tkvars$samtype)=="Formula", Ichk(toR(tkvars$Isam)), Ichk(getmode(diff(days))))){
            pda<-as.numeric(tclvalue(tkvars$pda)); pdb<-as.numeric(tclvalue(tkvars$pdb))
            rCP<-rCPgab(persistence_distn, pda, pdb, Ir)
            meanCP<<-rCP[1]; rhat<<-rCP[2]
            rCPhtext1<<-paste0("r = ",signif(rCP[2],3), " for Ir = ",Ir)
            if (bminok & bmaxok){
              blwr<-as.numeric(tclvalue(tkvars$blwr))
              bupr<-as.numeric(tclvalue(tkvars$bupr))
              rCPlwr<-rCPgab(persistence_distn, pda, blwr, Ir)
              rCPupr<-rCPgab(persistence_distn, pda, bupr, Ir)
              rCPhtext2<<-paste0(", with 95% CI: r \u2208 [",signif(rCPlwr[2],3),", ",signif(rCPupr[2],3),"]")
            }
          }
        } else {
          tkconfigure(pda.edit, bg=colors()[652])
          rCPhtext1<<-paste0("r = ", NA, " for Ir = ", Ir)
          rCPhtext2<<-paste0(", with 95% CI: r \u2208 [",NA,", ",NA,"]")
        }
        tkconfigure(rCP.lbl, text = paste0(rCPhtext1, rCPhtext2))
      })
      tkbind(pdb.edit,"<KeyRelease>", function() {
        pdb<-try(suppressWarnings(as.numeric(tclvalue(tkvars$pdb))),silent=T)
        if (class(pdb) == "try-error") pdbok<<-F else pdbok<<-pdbchk(pdb)
        if (pdbok) {
          tkconfigure(pdb.edit, bg='white')
          if (bminleg) { #pmin is a legitimate value (i.e. in (0, 1]) and entered into provisional parameter list; it's O.K. if pmin < pdb
            blwr<-as.numeric(tclvalue(tkvars$blwr))
            if (blwr < pdb){
              bminok<<-T
              tkconfigure(blwr.edit,bg='white')
            } else {
              bminok<<-F
              tkconfigure(blwr.edit,bg=colors()[652])
            }
          }
          if (bmaxleg) { #pmin is a legitimate value (i.e. in (0, 1]) and entered into provisional parameter list; it's O.K. if pmin < pdb
            bupr<-as.numeric(tclvalue(tkvars$bupr))
            if (bupr > pdb){
              bmaxok<<-T
              tkconfigure(bupr.edit,bg='white')
            } else {
              bmaxok<<-F
              tkconfigure(bupr.edit,bg=colors()[652])
            }
          }
          if (bmaxok & bminok & pdaok & pdbok){
            rCPlwr<-rCPgab(persistence_distn, pda, blwr, Ir)
            rCPupr<-rCPgab(persistence_distn, pda, bupr, Ir)
            rCPhtext2<<-paste0(", with 95% CI: r \u2208 [",signif(rCPlwr[2],3),", ",signif(rCPupr[2],3),"]")
          } else {
            rCPhtext2<<-paste0(", with 95% CI: r \u2208 [",NA,", ",NA,"]")
          }
        } else {
          tkconfigure(pdb.edit, bg=colors()[652])
        }
        if (bminok) tkconfigure(blwr.edit, bg='white') else tkconfigure(blwr.edit, bg=colors()[652])
        if (bmaxok) tkconfigure(bupr.edit, bg='white') else tkconfigure(bupr.edit, bg=colors()[652])
        if (pdaok & pdbok){
      #    Ir<<-as.numeric(tclvalue(tkvars$Ir))
          rCP<-rCPgab(persistence_distn, toR(tkvars$pda), toR(tkvars$pdb), Ir)
          meanCP<<-rCP[1]; rhat<<-rCP[2]
          rCPhtext1<<-paste0("r = ",signif(rCP[2],3), " for Ir = ", Ir)
        } else {
          rCPhtext1<<-paste0("r = NA for Ir = ", Ir)
        }
        tkconfigure(rCP.lbl, text = paste0(rCPhtext1, rCPhtext2))
      })

      tkbind(blwr.edit,"<KeyRelease>", function() {
        blwr<-try(suppressWarnings(as.numeric(tclvalue(tkvars$blwr))),silent=T)
        if (class(blwr) == "try-error") blwrok<<-F else blwrok<<-blwrchk(blwr)
        if (bminok) tkconfigure(blwr.edit, bg='white') else tkconfigure(blwr.edit, bg=colors()[652])
        if (bmaxok) tkconfigure(bupr.edit, bg='white') else tkconfigure(bupr.edit, bg=colors()[652])
        if (bmaxok & bminok & pdaok & pdbok){
          pda<-as.numeric(tclvalue(tkvars$pda))
          bupr<-as.numeric(tclvalue(tkvars$bupr))
      #    Ir<-as.numeric(tclvalue(tkvars$Ir))
          rCPlwr<-rCPgab(persistence_distn, pda, blwr, Ir)
          rCPupr<-rCPgab(persistence_distn, pda, bupr, Ir)
          rCPhtext2<<-paste0(", with 95% CI: r \u2208 [",signif(rCPlwr[2],3),", ",signif(rCPupr[2],3),"]")
        } else {
          rCPhtext2<<-paste0(", with 95% CI: r \u2208 [",NA,", ",NA,"]")
        }
        tkconfigure(rCP.lbl, text = paste0(rCPhtext1, rCPhtext2))
      })
      tkbind(bupr.edit,"<KeyRelease>", function() {
        bupr<-try(suppressWarnings(as.numeric(tclvalue(tkvars$bupr))),silent=T)
        if (class(bupr) == "try-error") buprok<<-F else buprok<<-buprchk(bupr)
        if (bminok) tkconfigure(blwr.edit, bg='white') else tkconfigure(blwr.edit, bg=colors()[652])
        if (bmaxok) tkconfigure(bupr.edit, bg='white') else tkconfigure(bupr.edit, bg=colors()[652])
        if (bmaxok & bminok & pdaok & pdbok){
          pda<-as.numeric(tclvalue(tkvars$pda))
          blwr<-as.numeric(tclvalue(tkvars$blwr))
      #    Ir<-as.numeric(tclvalue(tkvars$Ir))
          rCPlwr<-rCPgab(persistence_distn, pda, blwr, Ir)
          rCPupr<-rCPgab(persistence_distn, pda, bupr, Ir)
          rCPhtext2<<-paste0(", with 95% CI: r \u2208 [",signif(rCPlwr[2],3),", ",signif(rCPupr[2],3),"]")
        } else {
          rCPhtext2<<-paste0(", with 95% CI: r \u2208 [",NA,", ",NA,"]")
        }
        tkconfigure(rCP.lbl, text = paste0(rCPhtext1, rCPhtext2))
      })
      tkbind(crlev.edit,"<KeyRelease>", function() {
        crlev<-try(suppressWarnings(as.numeric(tclvalue(tkvars$crlev))),silent=T)
        if (class(crlev) == "try-error") syAok<<-F else syAok<<-syAchk(crlev)
        if (syAok) tkconfigure(crlev.edit, bg='white') else tkconfigure(crlev.edit, bg=colors()[652])
      })
      tkbind(syModule,"<Destroy>",function() {
        try(save(singleYearPrevious,paste0(.Rvar$datadir,'singleYearPrevious.Rdata')),silent=T)
        tkdestroy(syModule)
      })
      tkwm.deiconify(syModule)
      tkgrab.set(syModule)
      tkfocus(syModule)

      if (partial){
        # disable the arrival and prior frames
        if (!basicMode){
          for (obj in unlist(strsplit(tclvalue(tkwinfo("children", arrfFrame)), ' '))) tkconfigure(obj,state='disabled')
          tkconfigure(prior.lbox, state='disabled')
          tkconfigure(priorEditButton, state='disabled')
          tkconfigure(priorViewButton, state='disabled')
          tkconfigure(prior.lbox, bg = 'gray93')
        }
        tkconfigure(crlev.lbl, state='disabled')
        tkconfigure(crlev.edit, state='disabled', bg='gray93')
        tkconfigure(syCalcL, state = 'disabled')
        tkconfigure(syClose, text = "Cancel", command = function(){
          tkdestroy(syModule)
        })
        tkconfigure(syCalcg, command=function(){SYcalcg(writeg=F)})
      }
    },
    cd_form = function(){
#      datestmp<-NA # why store as R variable instead of tclday as tmp and then tkvars$days as quasi-permanent?
      datestmp <- toR(tkvars$days)
      tclday<-tclArray()
      for (i in 1:length(tkvars$days)) tclday[[i-1, 0]]<-tclvalue(tkvars$days[[i-1]])
      pasteFromClipboard<-function(){
        suppressWarnings(datestmp<-as.numeric(readClipboard(1)))
        if (length(datestmp)<2){
          tkmessageBox(icon='error',message="Error in data\n\nRequired: numeric column",type='ok')
          tkconfigure(dateOK,state='disabled')
          return (FALSE)
        } else { # do the rest of the error check
          if (sum(is.na(datestmp)>1)){
            tkmessageBox(icon='error',message="Error in data\n\nRequired: numeric column",type='ok')
            tkconfigure(dateOK,state='disabled')
            return (FALSE)
          } else if (is.na(datestmp[1])==1){
            if (!is.na(datestmp[1])){
              tkconfigure(dateOK,state='disabled')
              tkmessageBox(icon='error',message="Error in data\n\nRequired: numeric column",type='ok')
              return(FALSE)
            } else {
              datestmp<-datestmp[-1]
            }
          }
          if (sum(is.na(datestmp))>0){
            tkmessageBox(icon='error',message="Error in data\n\nRequired: numeric column",type='ok')
            tkconfigure(dateOK,state='disabled')
            return (FALSE)
          }
          if (datestmp[1]!=0) {
            tkmessageBox(icon="error",message=paste0("First search date must be 0 but it's = ",datestmp[1]),type='ok')
            tkconfigure(dateOK,state='disabled')
            return (FALSE)
          }
          if (sum(datestmp<0)>0) {
            tkmessageBox(icon="error",message="Dates must be positive",type='ok')
            tkconfigure(dateOK,state='disabled')
            return (FALSE)
          }
          if (length(unique(datestmp))<length(datestmp)){
            tkmessageBox(icon="error",message="Replicate dates")
            tkconfigure(dateOK,state='disabled')
            return (FALSE)
          }
          if (sum(diff(datestmp)<0)>0) {
            tkmessageBox(icon="error",message="Dates must be in increasing order")
            tkconfigure(dateOK,state='disabled')
            return (FALSE)
          }
          # if the error check goes through, then paste the dates onto a tktable and copy to singleYearProvisional
        }
        nold<-as.numeric(tclvalue(tcl(dateTable,"index","end","row")));
        nnew<-length(datestmp)-1
        tkconfigure(dateTable,state='normal',flashmode=T,flashtime=1)
        if (nnew < nold) tkdelete(dateTable,"rows",1, nold-nnew+1)
        for (i in 1:length(datestmp)) tclday[[i-1, 0]]<<-datestmp[i]
        tkconfigure(dateTable,rows=length(datestmp),state='disabled')
        tkconfigure(dateOK,state='normal')
        return (TRUE)
      }
      readCSV<-function(){
        fileName<-tclvalue(tkgetOpenFile(initialdir=.Rvar$csvpath))
        if (!nchar(fileName)) {
          return (FALSE)
        }
        tmp<-unlist(strsplit(fileName,'/'))
        .Rvar$csvpath<-paste(tmp[-length(tmp)],collapse='/')
        v<-try(read.csv(fileName,sep=",",as.is=T),silent=T)
        if (class(v)=="try-error" | dim(v)[2] != 1){
          tkmessageBox(icon='error',message=paste0(fileName,"\n\nError in file\nRequired: 1 column of dates"))
          tkconfigure(dateOK,state='disabled')
          return(FALSE)
        }
        datestmp<-v[,1]

        if (!is.numeric(datestmp)){
          tkmessageBox(icon='error',message=paste0(fileName,"\n\nError in data\nRequired: numeric dates"))
          tkconfigure(dateOK,state='disabled')
          return(FALSE)
        }
        if (datestmp[1]!=0){
          tkmessageBox(icon='error',message=paste0(fileName,"\n\nError in data\nRequired: header + search days with first date = 0"))
          tkconfigure(dateOK,state='disabled')
          return(FALSE)
        }
        if (sum(diff(datestmp)==0)>0){
          tkmessageBox(icon='error',message=paste0(fileName,"\n\nError in data\nDuplicate search dates"))
          tkconfigure(dateOK,state='disabled')
          return(FALSE)
        }
        if (sum(datestmp<0) > 0){
          tkmessageBox(icon='error',message=paste0(fileName,"\n\nError in data\nDates must be non-negative"))
          tkconfigure(dateOK,state='disabled')
          return(FALSE)
        }
        nold<-as.numeric(tclvalue(tcl(dateTable,"index","end","row")));
        nnew<-length(datestmp)
        tkconfigure(dateTable,state='normal',flashmode=T,flashtime=1)
        if (nnew < nold) {
          tkdelete(dateTable,"rows",1,nold - nnew + 1)
        }
        for (i in 1:length(datestmp)){
          tclday[[i-1,0]]<<-datestmp[i]
        }
        tkconfigure(dateTable,rows=length(datestmp),state='disabled')
        tkconfigure(dateOK,state='normal')
      }

      cdFrame <<- tktoplevel()
      tkgrab.set(cdFrame);  tkfocus(cdFrame)
      topMenu <- tkmenu(cdFrame); tkconfigure(cdFrame,menu=topMenu)
      tkwm.title(cdFrame,paste0("EoA, v", .Rvar$VER, " - Custom Dates"))
      tkwm.resizable(cdFrame,0,0)
      tkwm.deiconify(cdFrame)
      #try(tkwm.iconify(singleYearModule),silent=T)
      # menus
      editMenu <- tkmenu(topMenu,tearoff=FALSE,activebackground=colors()[125],activeforeground=colors()[109])
        tkadd(editMenu,"command",label="Read from file",command=readCSV)
        tkadd(editMenu,"command",label="Paste",command=pasteFromClipboard)
        tkadd(editMenu,"command",label="Copy", command=function() writeClipboard(paste(c("dates",datestmp),sep='\n')))
        tkadd(topMenu,"cascade",label="Edit",menu=editMenu)
      helpMenu <- tkmenu(topMenu,tearoff=FALSE,activebackground=colors()[125],activeforeground=colors()[109])
        tkadd(helpMenu,"command",label="About",command=function()tkmessageBox(title='Evidence of Absence (EoA)',message=about_text))
        tkadd(topMenu,"cascade",label="Help",menu=helpMenu)
      tk.cstart<-tclVar(tclvalue(tkvars$firstsearch))
      cstart.lbl<-tklabel(cdFrame,text='Start of monitoring')
      cstart.edit<-tkentry(cdFrame,width=10,textvariable=tk.cstart,bg='white')
      tkbind(cstart.edit,"<KeyRelease>", function() {
        cst<-suppressWarnings(as.numeric(tclvalue(tk.cstart)))
        test<-try(as.Date(tclvalue(tk.cstart)),silent=T)
        if (class(test)=="try-error"){
            tkconfigure(cstart.edit,bg=colors()[652])
            tkconfigure(dateOK,state='disabled')
            return(F)
        }
        tkconfigure(dateOK, state='normal')
        tkconfigure(cstart.edit,bg='white')
      })# create table for dates and enter dates from singleYearPrevious$days
      schedule.lbl<-tklabel(cdFrame,text='Search schedule\n(days after first search)',justify='left')
      tblF<-tkframe(cdFrame)
      dateTable<-tcltk2::tk2table(tblF,
        rows=length(tclday),
        cols=1,
        selectmode="extended",
        colwidth="5",
        variable=tclday,
        resizeborders="none",
        height="5",
        selecttitle = 1#,
#        sparsearray = 1
      )
      yscr <- ttkscrollbar(tblF, orient='vertical', command=function(...) tkyview(dateTable,...))
      tkconfigure(dateTable,yscrollcommand=function(...) tkset(yscr,...))
      tkconfigure(dateTable,state='disabled')# control whether rows and/or columns can be edited
      tkconfigure(dateTable,multiline="0")# prevent line-wrapping within cells
      tkconfigure(dateTable,rowseparator="\"\n\"",colseparator="\"\t\"") # separation characters for copying into clipboard
      tkgrid(dateTable,row=0,column=0)
      tkgrid(yscr,row=0,column=1,sticky='ns')#row=1,column=1,sticky='news')
      dateButtonsFrame<-tkframe(cdFrame)
      dateOK<-tkbutton(dateButtonsFrame,text='OK',width=7, command=function() {
        nold<-length(tkvars$days)
        tclday$active<-NULL
        if (nold > length(tclday)){
          for (i in (length(tclday)-1):(nold-1)) tkvars$days[[i]]<-NULL
        }
        for (i in 1:length(tclday)){
          tkvars$days[[i-1]]<<-tclvalue(tclday[[i-1,0]])
        }
        tclvalue(tkvars$firstsearch)<-tclvalue(tk.cstart)
        Ir<<-round(max(toR(tclday))/(length(tclday)-1),1)
        updateCPfieldLabel(.Rvar$CPdata)
        updateCPhandLabel()
        tkconfigure(cuslab, text = paste0("    span = ", max(toR(tkvars$days)), ", I (mean) = ", round(max(toR(tkvars$days))/(length(toR(tkvars$days))-1),1)))
        tkdestroy(cdFrame)
      })
      dateCancel<-tkbutton(dateButtonsFrame,text="Cancel",width=7,command=function(){tkdestroy(cdFrame)})
      tkgrid(dateOK)
      tkgrid(dateCancel)
      tkgrid(cstart.lbl,cstart.edit)
      tkgrid(schedule.lbl,sticky='n')
      tkgrid(tblF, column=1,row=1,pady=12,rowspan=2)
      tkgrid(dateButtonsFrame,row=0,column=2,sticky='n',padx=10)
    },
    cpr_form = function(){
      tclcpr<-tclArray()
      tclcpr[[0,0]]<-"m"
      tclcpr[[0,1]]<- as.tclObj("P(M = m)",drop=T)
      custom.prior<-toR(tkvars$custom.prior)
      v<-cbind(1:length(custom.prior)-1,custom.prior)
      if (sum(is.na(as.numeric(custom.prior))) == 0 && custom.prior != ''){
        # create table for prior and enter initial values from singleYearProvisional$custom.prior (NOTE: in the singleYearProvisional db, the p(M) values start for m = 0
        for (i in 1:length(custom.prior)){
          tclcpr[[i,0]]<-i
          tclcpr[[i,1]]<-round(custom.prior[i],6)
          Buttstate<-'normal'
        }
      } else {
        tclcpr[[1,0]]<-"NA"
        tclcpr[[1,1]]<-"NA"
        Buttstate<- 'disabled'
      }
      pasteFromClipboardcpr<-function(){
        junk<-suppressWarnings(try(read.table(file='clipboard',as.is=T)))
        if (class(junk) != "try-error") v<-na.omit(junk) else {tkmessageBox(icon='error',message="Error in data.\nRequired: two columns [m and p(M = m)].\nCheck clipboard."); return(F) }
        if (length(dim(junk)) != 2 || dim(junk)[2]!=2) {tkmessageBox(icon='error',message="Error in data.\nRequired: two columns [m and p(M = m)].\nCheck clipboard."); return(F) }
        v<-suppressWarnings(array(as.numeric(cbind(junk$V1,junk$V2)),dim=dim(junk)))
        if (sum(is.na(v[1,]))>0) v<-v[-1,]
        if (checkprior(v)) {
          entercpr(v)
        } else {
          tkconfigure(cprOK,state='disabled')
          tkconfigure(cprView,state='disabled')
        }
      }
      readCSVcpr<-function(){ # need to read in and parse two columns
        fileName<-tclvalue(tkgetOpenFile(initialdir=.Rvar$csvpath))
        if (!nchar(fileName)) {
          tkmessageBox(message="No file was selected!")
          return (FALSE)
        }
        tmp<-unlist(strsplit(fileName,'/'))
        .Rvar$csvpath<-paste(tmp[-length(tmp)],collapse='/')
        v<-na.omit(try(read.csv(fileName,sep=",",as.is=T),silent=T)[,1:2]) # debug: <<     production: <
      #  v<<-try(read.csv(fileName,sep=",",as.is=T))
        if (class(v) != "try-error") v<-na.omit(v) else {tkmessageBox(icon='error',message=paste0("Error in data (",fileName,"\nRequired: two columns [m and p(M = m)].\nCheck file.")); return(F) }
        if (checkprior(v)) entercpr(v)
      }
      entercpr<-function(v){
        custom.prior<<-v[,2]
        tclcpr<<-tclArray()
        tclcpr[[0,0]]<-'m'; tclcpr[[0,1]]<-as.tclObj("p(M = m)",drop=T)
        for (i in 1:length(custom.prior)) {tclcpr[[i,0]]<- i-1; tclcpr[[i,1]] <- custom.prior[i]}
        tkconfigure(cprTable,rows=length(custom.prior)+1,cols=2,variable=tclcpr)
        tkconfigure(cprOK,state='normal')
        tkconfigure(cprView,state='normal')
      }
      cprFrame <<- tktoplevel()
      tkgrab.set(cprFrame);  tkfocus(cprFrame)
      topMenu <- tkmenu(cprFrame); tkconfigure(cprFrame,menu=topMenu)
      tkwm.title(cprFrame,paste0("EoA, v", .Rvar$VER, " - Custom Prior"))
      tkwm.resizable(cprFrame,0,0)
      tkwm.deiconify(cprFrame)
      # menus
      editMenu <- tkmenu(topMenu,activebackground=colors()[125],activeforeground=colors()[109])
        tkadd(editMenu,"command",label="Read CSV",command=readCSVcpr)
        tkadd(editMenu,"command",label="Copy", command=function() write.table(v,file='clipboard',row.names=F,col.names=c('m','p(M = m)'),sep='\t'))
        tkadd(editMenu,"command",label="Paste",command=pasteFromClipboardcpr)
        tkadd(topMenu,"cascade",label="Edit",menu=editMenu)
      helpMenu <- tkmenu(topMenu,activebackground=colors()[125],activeforeground=colors()[109])
        tkadd(helpMenu,"command",label="About",command=function()tkmessageBox(title='Evidence of Absence (EoA)',message=about_text))
        tkadd(topMenu,"cascade",label="Help",menu=helpMenu)
      cprTable<-tkwidget(cprFrame,"table",rows=length(custom.prior)+1,cols=2, selectmode="extended",colwidth="10",
        variable=tclcpr, resizeborders="none",titlerows=1,height=10, selecttitle = 1)
      yscr <- tkscrollbar(cprFrame,orient='vertical', command=function(...)tkyview(cprTable,...))
      tkconfigure(cprTable,yscrollcommand=function(...) tkset(yscr,...))
      tkgrid(cprTable,sticky='n')
      tkgrid(yscr,sticky='nse',column=2,row=0)
      tkconfigure(cprTable,state='disabled')# control whether rows and/or columns can be edited
      tkconfigure(cprTable,multiline="0")# prevent line-wrapping within cells
      tkconfigure(cprTable,rowseparator="\"\n\"",colseparator="\"\t\"") # separation characters for copying into clipboard
      tkconfigure(cprTable,variable=tclcpr)
      cprButtFrame<-tkframe(cprFrame)
      cprOK<-tkbutton(cprButtFrame,text='OK',command=function() {
        tkvars$custom.prior<<-tclArray()
        for (i in 1:length(custom.prior)){
          tkvars$custom.prior[[i-1]]<<-custom.prior[i]
        }
        prior_f<<-"Custom"
        prior_M<<-custom.prior
        custom.prior<<-custom.prior
        prok<<-T
        tkdestroy(cprFrame)
        updateprlab
      })
      cprCancel<-tkbutton(cprButtFrame,text="Cancel",command=function(){ # if the old custom prior is not workable, revert to objective prior
        prok<<-T
        custom.prior<-toR(tkvars$custom.prior)
        if (sum(is.na(suppressWarnings(as.numeric(custom.prior))))) {
          tkconfigure(priorViewButton,state='normal')
          tkconfigure(priorEditButton,state='disabled')
          tkselection.clear(prior.lbox,which(prnames=='Custom')-1)
          tkselection.set(prior.lbox,which(prnames=="Objective")-1)
          tkactivate(prior.lbox,which(prnames=="Objective")-1)
          prior_f<<-"Objective"
          prior_M<<-toR(tkvars$objective.prior)
          tkconfigure(prior.lbl,text="           ...           ")
        } else if (sum(custom.prior)==1) {
          tkconfigure(priorViewButton,state='normal')
          tkconfigure(priorEditButton,state='normal')
        } else {
          tkconfigure(priorViewButton,state='normal')
          tkconfigure(priorEditButton,state='disabled')
          tkselection.set(prior.lbox,which(prnames=="Objective")-1)
          tkactivate(prior.lbox,which(prnames=="Objective")-1)
          prior_f<<-"Objective"
          prior_M<<-toR(tkvars$objective.prior)
          tkconfigure(prior.lbl,text="           ...           ")
        }
        tkdestroy(cprFrame)
      })
      cprView<-tkbutton(cprButtFrame,text="View",command=function(){
        # cprView button is disabled unless data in table are OK
        plotPrior(custom.prior,"Custom")
      })
      ButtWidth<-12
      tkconfigure(cprOK,width=ButtWidth, state = Buttstate)
      tkconfigure(cprView,width=ButtWidth, state = Buttstate)
      tkconfigure(cprCancel,width=ButtWidth)
      tkgrid(cprOK)
      tkgrid(cprView)
      tkgrid(cprCancel)
      tkgrid(cprButtFrame,row=0,column=3,sticky='n')
      if (sum(is.na(toR(tkvars$custom.prior))>0)){
        tkconfigure(cprOK,state='disabled')
        tkconfigure(cprView,state='disabled')
      }
    },
    SYcalcPost = function(writepost){ # first need to do error-checking on the data set
      # estimate g
      if (!SYcalcg(writeg=F)) return(F) # this gives Bab as beta parameters for distribution of ghat
#      prior_M <- NA
      with(.Rvar$singleYearPrevious,{
        # calculate posterior
        if (prior_f == "Objective"){
          if (.Rvar$syresult$Bab[1]>0){
            mmax <-fmmax.ab(X, .Rvar$syresult$Bab[1], .Rvar$syresult$Bab[2]) # find the maximum m to sum over in calculating posterior (for beta-binomial)
          } else {
            mmax<-fmmax(X, .Rvar$syresult$ghat)
          }
          prior_M <- diff(sqrt(0:mmax))/sum(diff(sqrt(0:mmax)))
        } else {
          if (sum(is.na(as.numeric(prior_M))) > 0) {
            tkmessageBox(message = "Error in prior distribution. Missing values or non-numeric data")
            return(F)
          }
          if (sum(prior_M < 0) > 0){
            tkmessageBox(message = "Error in prior distribution. Negative probabilities")
            return(F)
          }
          if (abs(sum(prior_M)-1)>0.00001){
            tkmessageBox(message = "Error in prior distribution. Probabilities do not sum to 1")
            return(F)
          }
        }
        mmax<-length(prior_M)-1
      #  Bab<-singleYearProvisional$Bab # this has been calculated and stored in Bab on main page
        M<-X:mmax
       # posterior for M | X
        if (.Rvar$syresult$BabRaw[1]>0){
          pXgM<-VGAM::dbetabinom.ab(X,size=M,shape1=.Rvar$syresult$BabRaw[1],shape2=.Rvar$syresult$BabRaw[2]) # the probabilities of X for M = 0:mmax
        } else {
          pXgM<-dbinom(X, size=M, prob=.Rvar$syresult$ghat/a) # the probabilities of X for M = X:mmax
        }
        pM<-prior_M[X:mmax+1]
        pMgX<-pXgM*pM; pMgX<-pMgX/sum(pMgX) # posterior distribution for M (ignoring M < X, which has probability = zero)
        pMgX<-c(rep(0,X),pMgX)
        if (pMgX[length(pMgX)]<0.001) pMgX<-pMgX[cumsum(pMgX)<0.999]
        pMgX.raw<<-pMgX/sum(pMgX)
        .Rvar$pMgX.raw<-pMgX.raw

        M<-X:mmax
       # posterior for M | X
        if (.Rvar$syresult$Bab[1]>0){
          pXgM<-VGAM::dbetabinom.ab(X,size=M,shape1=.Rvar$syresult$Bab[1],shape2=.Rvar$syresult$Bab[2]) # the probabilities of X for M = 0:mmax
        } else {
          pXgM<-dbinom(X,size=M,prob=.Rvar$syresult$ghat) # the probabilities of X for M = X:mmax
        }
        pM<-prior_M[X:mmax+1]
        pMgX<-pXgM*pM; pMgX<-pMgX/sum(pMgX) # posterior distribution for M (ignoring M < X, which has probability = zero)
        pMgX<-c(rep(0,X),pMgX)
        if (pMgX[length(pMgX)]<0.001) pMgX<-pMgX[cumsum(pMgX)<0.999]
        pMgX<<-pMgX/sum(pMgX)
        .Rvar$pMgX<-pMgX
        # extrapolate to full year
#        if (tclvalue(tkvars$arrfun)=="Compound"){
          Bab<-.Rvar$syresult$BabAnn
          if (Bab[1]>0){
            mmax <-fmmax.ab(X, Bab[1], Bab[2]) # find the maximum m to sum over in calculating posterior (for beta-binomial)
          } else {
            mmax<-fmmax(X, .Rvar$syresult$ghat)
          }
          prior_M <- diff(sqrt(0:mmax))/sum(diff(sqrt(0:mmax)))
          mmax<-length(prior_M)-1
          M<-X:mmax
      #    Bab<-singleYearProvisional$BabAnn
          if (Bab[1]>0){
            pXgM<-VGAM::dbetabinom.ab(X,size=M,shape1=Bab[1],shape2=Bab[2]) # the probabilities of X for M = 0:mmax
          } else {
            pXgM<-dbinom(X,size=M,prob=.Rvar$syresult$ghat*v) # the probabilities of X for M = X:mmax
          }
          pM<-prior_M[X:mmax+1]
          pMgX<-pXgM*pM; pMgX<-pMgX/sum(pMgX) # posterior distribution for M (ignoring M < X, which has probability = zero)
          pMgX<-c(rep(0,X),pMgX)
          if (pMgX[length(pMgX)]<0.001) pMgX<-pMgX[cumsum(pMgX)<0.999]
          pMgX.ann<<-pMgX/sum(pMgX)
          .Rvar$pMgX.ann<-pMgX.ann
#        }
        .Rvar$singleYearPrevious$prior_M<-prior_M
      })
#      .Rvar$singleYearPrevious$prior_M<-prior_M
#      prior_M<<-prior_M
      # plot posterior
      plotPost()
      # table of results
      writePostResults()
    },
    fixArrData = function(parms){
      for (i in 1:3) tkvars$arrcomponents[[i-1]] <- parms$arrcomponents[i]
      tclvalue(tkvars$arrstart) <- parms$arrstart
      tclvalue(tkvars$lwr.u) <- parms$lwr.u
      tclvalue(tkvars$upr.u) <- parms$upr.u
      tclvalue(tkvars$wt.u) <- parms$wt.u
      tclvalue(tkvars$lwr.p1) <- parms$lwr.p1
      tclvalue(tkvars$upr.p1) <-parms$upr.p1
      tclvalue(tkvars$wt.p1)<-parms$wt.p1
      tclvalue(tkvars$a.p1)<-parms$a.p1
      tclvalue(tkvars$b.p1)<-parms$b.p1
      tclvalue(tkvars$lwr.p2)<-parms$lwr.p2
      tclvalue(tkvars$upr.p2)<-parms$upr.p2
      tclvalue(tkvars$wt.p2)<-parms$wt.p2
      tclvalue(tkvars$a.p2)<-parms$a.p2
      tclvalue(tkvars$b.p2)<-parms$b.p2
    },
    updatearrrad = function(arf){
      if (arf=="Uniform"){
        tkconfigure(editArrivalsbutton,state='disabled')
        tkconfigure(viewArrivalsbutton,state='disabled')
      } else if (arf=="Compound") {
        tkconfigure(editArrivalsbutton,state='normal')
        tkconfigure(viewArrivalsbutton,state='normal')
      }
    },
    updatePrevious = function(){
      tmp<-.Rvar$singleYearPrevious
      for (nm in syVar) tmp[[nm]]<-toR(tkvars[[nm]])
      for (nm in syArray) tmp[[nm]]<-toR(tkvars[[nm]])
      tmp$persistence_distn <- persistence_distn
      tmp$prior_f<-prior_f
      tmp$prior_M<-prior_M
      .Rvar$CPdata$persistence_distn<-tclvalue(tkvars$persistence_distn)
      if (tclvalue(tkvars$perstype) == "field") .Rvar$CPdataPrevious<-.Rvar$CPdata
      if (tclvalue(tkvars$SEopt) == "f") .Rvar$pkdatPrevious<-.Rvar$pkdat
      if (!partial) .Rvar$singleYearPrevious<-tmp else .Rvar$syscPrevious<-tmp
    },
    syChkAll = function(){
    # no check on arrival function or prior distribution (which are assumed to be correct b/c they are edited via separate forms)
    # this function gives a quick error check on data that have been entered by hand or edited from a preloaded data set
    # ...not sufficient for checking data loaded from a file
      if (! (Xok & aok & SEnok & SExok & kok & pdaok & pdbok & bminok & bmaxok & syAok)) return(F)
      samtype<-tclvalue(tkvars$samtype) # the radio buttons should be correct by default...no need to check
      if (samtype=="Formula"){ # check number of searches and search interval
        Isam<-try(as.numeric(toR(tkvars$Isam)),silent=T)
        if (class(Isam)=="try-error") return(F) else if (!Ichk(Isam)) return(F)
        nsearch<-try(as.numeric(toR(tkvars$nsearch)),silent=T)
        if (class(nsearch)=="try-error") return(F) else if (!nsearchchk(nsearch)) return(F)
      } else { # check specific search days
        days<-toR(tkvars$days)
        if (days[1] != 0) {
          tkmessageBox(message=paste0("Error. Initial search date (", days[1], ") should 0."))
          return(F)
        }
        if (min(diff(days)) <= 0){
            tkmessageBox(message=paste0("Error in search dates"))
            return(F)
        }
      }
      firstsearch<-try(tclvalue(tkvars$firstsearch),silent=T)
      if (class(firstsearch)=="try-error") return(F) else if (!startchk(firstsearch)) return(F)
      return(T)
    },
    SYcalcg = function(writeg){
    # check for errors before calculating
    # if no errors, translate the appropriate tcl variables into R variables for calculation
      if (!syChkAll()){ # then error-message and abort
        tkmessageBox(message="Error in data. Aborting calculation...", icon='error')
        return(F)
      }
    # no errors, so write data to R (assumes prior is OK)
      updatePrevious()
      if (tclvalue(tkvars$perstype) == 'hand'){
        if (persistence_distn == "Exponential"){
          tclvalue(tkvars$pda.e) <<- tclvalue(tkvars$pda)
          tclvalue(tkvars$pdb.e) <<- tclvalue(tkvars$pdb)
          tclvalue(tkvars$blwr.e) <<- tclvalue(tkvars$blwr)
          tclvalue(tkvars$bupr.e) <<- tclvalue(tkvars$bupr)
        } else if (persistence_distn=="Weibull") {
          tclvalue(tkvars$pda.w) <<-  tclvalue(tkvars$pda)
          tclvalue(tkvars$pdb.w) <<-  tclvalue(tkvars$pdb)
          tclvalue(tkvars$blwr.w) <<- tclvalue(tkvars$blwr)
          tclvalue(tkvars$bupr.w) <<- tclvalue(tkvars$bupr)
        } else if (persistence_distn=="Log-Logistic") {
          tclvalue(tkvars$pda.ll) <<- tclvalue(tkvars$pda)
          tclvalue(tkvars$pdb.ll) <<- tclvalue(tkvars$pdb)
          tclvalue(tkvars$blwr.ll) <<-tclvalue(tkvars$blwr)
          tclvalue(tkvars$bupr.ll) <<-tclvalue(tkvars$bupr)
        } else if (persistence_distn=="Lognormal") {
          tclvalue(tkvars$pda.ln) <<-tclvalue(tkvars$pda)
          tclvalue(tkvars$pdb.ln) <<-tclvalue(tkvars$pdb)
          tclvalue(tkvars$blwr.ln) <<-tclvalue(tkvars$blwr)
          tclvalue(tkvars$bupr.ln) <<-tclvalue(tkvars$bupr)
        }
      }
      for (nm in names(tkvars)){
        if (!(nm %in% c("prior_f", "prior_M", "persistence_distn", "arrcomponents"))) assign(nm,toR(tkvars[[nm]]))
      }
      if (perstype=='field'){
        blwr <- .Rvar$CPdataPrevious$blwr
        bupr <- .Rvar$CPdataPrevious$bupr
        pda  <- .Rvar$CPdataPrevious$pda
        pdb  <- .Rvar$CPdataPrevious$pdb
        persistence_distn <- .Rvar$CPdataPrevious$persistence_distn
      }
      samtype <- tclvalue(tkvars$samtype)
      if (samtype=='Formula'){
        Isam <- as.numeric(tclvalue(tkvars$Isam))
        nsearch <- as.numeric(tclvalue(tkvars$nsearch))
        days <- (0:nsearch)*Isam
      } else {
        days <- toR(tkvars$days)
        Isam <- round(max(days)/(length(days)-1),1)
        nsearch <- length(days)-1
      }
      firstsearch <- tclvalue(tkvars$firstsearch)
      prior_f <<- prior_f
      if (tclvalue(tkvars$SEopt)=='f'){
        SEx<-.Rvar$pkdat$X[1]
        SEn<-.Rvar$pkdat$M[1]
      }
      fba<<-SEx+0.5
      fbb<<-SEn-SEx+.5
      if (tclvalue(tkvars$arrfun)=='Compound' & !basicMode){
        if (max(days) > 365){
          tkmessageBox(message="
            Analysis with compound arrival function is limited to one year.
              (1) Use a monitoring period of one year or less, or
              (2) Use uniform arrival function.
            ",icon='error')
          return(F)
        }
        syr<-as.numeric(format(as.Date(firstsearch),"%Y"))
        arryr <- ifelse(as.Date(firstsearch)-as.Date(paste0(syr,"-01-01")) < arrstart,syr-1,syr)
        s0<-as.numeric(as.Date(firstsearch)-as.Date(arrstart,origin=paste0(arryr,"-01-01")))
        FullYear <- paste0(format(as.Date(arrstart,origin=paste0(arryr,"-01-01")),"%d-%b-%Y"), " through ", format(as.Date(arrstart-1,origin=paste0(arryr+1,"-01-01")),"%d-%b-%Y"))

        # integral of arrival function from arrstart to monitoring start
        # monitoring ends at:
        sf<-s0+max(days)
        if (sf > 365){ # then the monitoring season extends over one year past the start of arrivals...define new arrival date
          tkmessageBox(message="Monitoring season extends more than one year after arrivals begin.\nReset arrival function start date or monitoring schedule.",icon='error')
          return(F)
        }
      } else {
        FullYear<-NA
      }
      MonitoredPeriod <- paste0(format(as.Date(firstsearch),"%d-%b-%Y")," through ",format(as.Date(as.numeric(as.Date(firstsearch))+max(days),origin="1970-01-01"),"%d-%b-%Y"))

      # check whether arrival start date and monitoring schedule are all completed within a year
      # need s0 +
      pda0<-NA; pdb0<-NA
      if (perstype == "hand"){
        pdb0<-pdb
        pda0<-pda
      } else {
        with(.Rvar$CPdataPrevious, {
          if (persistence_distn == "Exponential"){
            # NOTE: exponential b = meanCP = exp(.Rvar$CPdataPrevious$mod.e$coef), variance of estimated mod.e$coef = mod.e$var
            pdb0 <<- exp(mod.e$coef[1])
            pda0 <<- 1/pdb0
          } else if (persistence_distn == "Weibull"){
            # NOTE: weibull a = shape = 1/mod.w$scale, b = scale = exp(mod.w$coef[1])
            pdb0 <<- exp(mod.w$coef[1])
            pda0 <<- 1/mod.w$scale
          } else if (persistence_distn == "Log-Logistic"){
            # NOTE: log-logistic a = shape = 1/mod.ll$scale, b = scale = exp(mod.ll$coef[1])
            pdb0 <<- exp(mod.ll$coef[1])
            pda0 <<- 1/mod.ll$scale
          } else if (persistence_distn == "Lognormal"){
            # NOTE: lognormal a = sdlog^2 = mod.ln$scale^2, b = meanlog = mod.ln$coef[1]
            pdb0 <<- mod.ln$coef[1]
            pda0 <<- mod.ln$scale^2
          }
        })
      }
      ###1. setting estimation control parameters
      ##  a. search limit: required number of searches after arrival to include in estimate of searcher efficiency
      ##   [when the number of searches is high, including them all in the estimation is calculation intensive but does not contribute signficantly to the result]
      f0<-SEx/SEn
      if (tclvalue(tkvars$SEopt)=='f') k <- .Rvar$pkdat$M[1]*.Rvar$pkdat$X[2]/(.Rvar$pkdat$M[2]*.Rvar$pkdat$X[1])
      ind1<-rep(1:nsearch,times=nsearch:1)
      ind2<-ind1+1
      ind3<-unlist(lapply(1:nsearch, function(x) x:nsearch))+1
      schedule.index<-cbind(ind1,ind2,ind3)
      schedule<-cbind(days[ind1],days[ind2],days[ind3])
      nmiss<-schedule.index[,3]-schedule.index[,2]
      maxmiss<-max(nmiss)
      powk<-cumprod(c(1,rep(k,maxmiss))) # vector of k^i's
      notfind<-cumprod(1-f0*powk[-length(powk)])
      nvec<-c(1,notfind)*f0
      pfind.si<-nvec*powk # conditional probability of finding a carcass on the ith search (row) after arrival for given (simulated) searcher efficiency (column)
      # persistences
      intxsearch<-unique(cbind(schedule[,2] - schedule[,1], schedule[,3] - schedule[,2]), MAR=1)
      ppersu<-ppersist(persistence_distn, t_arrive0=0, t_arrive1=intxsearch[,1], t_search = intxsearch[,1]+intxsearch[,2], pda = pda0, pdb = pdb0)
      ppersu[is.na(ppersu)] <- max(ppersu[!is.na(ppersu)])
      if (tclvalue(tkvars$arrfun) =="Uniform" | basicMode) { # fraction of total carcasses that arrive in each search interval...inference is limited to the monitoring period
        arrvec<-(schedule[,2]-schedule[,1])/max(days) #
      } else if (tclvalue(tkvars$arrfun)=="Compound" & !basicMode){
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
        # fraction of carcasses arriving in each search interval
        # integrate the arrival function by components
        arrvec0.u<-numeric(nsearch) # fraction of carcasses falling in each search interval
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
      prob_obs<-numeric(dim(schedule)[1])
      for (i in 1:length(prob_obs)){
        prob_obs[i]<-pfind.si[nmiss[i]+1] *
          ppersu[which(
            abs(intxsearch[,1]-(schedule[i,2]-schedule[i,1]))<0.001 &
            abs(intxsearch[,2]-(schedule[i,3]-schedule[i,2]))<0.001),] *
          arrvec[i]
      }
      prob_obs<-prob_obs*a
      ggnm<-numeric(maxmiss+1)
      for (missi in 0:maxmiss){
        ggnm[missi+1]<-sum(prob_obs[nmiss==missi])
      }
      if (nsearch>10){
        iskip<-min(which(cumsum(ggnm)/sum(ggnm)>.99))+1 # cutting off the end introduces a bias. Correct by multiplying the final g's by sum(ggnm)/ggnm[iskip]
        gadj<-sum(ggnm)/sum(ggnm[1:iskip])
      } else {
        iskip<-maxmiss
        gadj<-1
      }
      ##  b. determine how many simulation draws are needed: choose nsim large enough so that the estimated ghat is known within 1% [or SE(ghat) < 0.005*ghat]
      r0<-ppersist(persistence_distn,t_arrive0 = 0, t_arrive1=round(Isam), t_search = round(Isam), pda = pda0, pdb = pdb0) # need an "if"
      f0<-SEx/SEn
      g0<-f0*r0
      # estimated variance
      nsim<-1000
      if (perstype == "hand"){
        if (persistence_distn=="Weibull" | persistence_distn=="Exponential"){
          pdb<-rlnorm(nsim,sdlog=(log(bupr)-log(blwr))/(2*1.96),meanlog=log(pdb0)) # from a lognormal distribution with mean = log(pdb) and sd = (log(maxb)-log(minb))/4
        } else if (persistence_distn=="Lognormal"){
#          pdb<-pmax(rnorm(nsim,mean=pdb0,sd=(bupr-blwr)/(2*1.96)),.1)
          pdb<-rnorm(nsim,mean=pdb0,sd=(bupr-blwr)/(2*1.96))
        } else if (persistence_distn=="Log-Logistic"){
          pdb<-rlnorm(nsim,meanlog=log(pdb0),sdlog=(log(bupr)-log(blwr))/(2*1.96)) # not entirely satisfactory...doesn't look great for highly skewed distribution of b
        }
        rr <- as.vector(ppersist(persistence_distn,t_arrive0 = 0, t_arrive1=round(Isam), t_search = round(Isam), pda = pda0, pdb = pdb))
        rr[is.na(rr)] <- max(rr[!is.na(rr)])
        Vhatr<-var(rr)
      } else {
        CPab<-array(dim=c(nsim, 2))
        if (persistence_distn == "Exponential"){
          CPab[,2]<-exp(rnorm(nsim,mean=.Rvar$CPdataPrevious$mod.e$coef, sd=sqrt(.Rvar$CPdataPrevious$mod.e$var[1])))
          CPab[,1]<-pda0
        } else if (persistence_distn == "Weibull"){
          CPparms<-MASS::mvrnorm(nsim,c(.Rvar$CPdata$mod.w$coef[1],log(.Rvar$CPdata$mod.w$scale)),.Rvar$CPdata$mod.w$var)
          CPab[,1]<-1/exp(CPparms[,2]) #shape
          CPab[,2]<-exp(CPparms[,1]) # scale
          i0<-which(CPab[,1]<=0 | (Isam/CPab[,2])^CPab[,1]<1e-320)
          while(length(i0)>0){
            CPparms[i0,]<-MASS::mvrnorm(length(i0),c(.Rvar$CPdata$mod.w$coef[1],log(.Rvar$CPdata$mod.w$scale)), .Rvar$CPdata$mod.w$var)
            CPab[i0,1]<-1/exp(CPparms[i0,2])
            CPab[i0,2]<-exp(CPparms[i0,1])
            i0<-which(CPab[,1]<=0 | (Isam/CPab[,2])^CPab[,1]<1e-320)
          }
        } else if (persistence_distn == "Log-Logistic"){
          # NOTE: log-logistic a = shape = 1/mod.ll$scale, b = scale = exp(mod.ll$coef[1])
          CPparms<-MASS::mvrnorm(nsim,c(.Rvar$CPdata$mod.ll$coef[1],log(.Rvar$CPdata$mod.ll$scale)),.Rvar$CPdata$mod.ll$var)
          CPab[,1]<-1/exp(CPparms[,2]) #shape
          CPab[,2]<-exp(CPparms[,1]) # scale
          i0<-which(CPab[,1]<=0)
          while(length(i0)>0){
            CPparms[i0,]<-MASS::mvrnorm(length(i0),c(.Rvar$CPdata$mod.ll$coef[1],log(.Rvar$CPdata$mod.ll$scale)), .Rvar$CPdata$mod.ll$var)
            CPab[i0,1]<-1/log(CPparms[i0,2])
            CPab[i0,2]<-exp(CPparms[i0,1])
            i0<-which(CPab[,1]<=0)
          }
        } else if (persistence_distn == "Lognormal"){
          # NOTE: lognormal a = sdlog^2 = mod.ln$scale^2, b = meanlog = mod.ln$coef[1]
          CPparms<-MASS::mvrnorm(nsim,c(.Rvar$CPdata$mod.ln$coef[1], log(.Rvar$CPdata$mod.ln$scale)),.Rvar$CPdata$mod.ln$var)
          CPab[,1]<- exp(CPparms[,2])^2 #shape
          CPab[,2]<- CPparms[,1] # scale

        }
        rr<-as.vector(ppersist(persistence_distn,t_arrive0 = 0, t_arrive1=round(Isam), t_search = round(Isam), pda = CPab[,1], pdb = CPab[,2]))
        rr[is.na(rr)]<-max(rr[!is.na(rr)])
        Vhatr<-var(rr)
      }
      Vhatf<-fba*fbb/((fba+fbb)^2*(fba+fbb+1))
      shatg<-sqrt(f0^2*Vhatr+r0^2*Vhatf+Vhatf*Vhatr)
      nsim<-min(round((shatg/(g0*0.003))^2),20000)
      ###2. estimation of g
      #a. subset the search schedule (to ignore probabilities of detection carcasses after they have been missed several times):
      if (nsim<=1 | shatg < 0.001){
        ghat<-sum(ggnm)
        prob_obs<-sum(ggnm)
        Bab<-c(-1, -1)    # this is an indicator of fixed g
        BabAnn<-c(-1, -1)
        BabRaw<-c(-1, -1)
      } else {
        schedule<-cbind(days[ind1],days[ind2],days[ind3])[ind2>=ind3-iskip+1,] # these are defined in section 1
        schedule.index<-cbind(ind1,ind2,ind3)[ind2>=ind3-iskip+1,] #columns for arrival interval and search number
        nmiss<-schedule.index[,3]-schedule.index[,2]
        maxmiss<-max(nmiss)
        # searcher efficiencies
        if (tclvalue(tkvars$SEopt)=='h'){
          f<-rbeta(nsim,shape1=fba,shape2=fbb)
          if (maxmiss>0){
            powk<-cumprod(c(1,rep(k,maxmiss))) # vector of k^i's
            notfind<-apply(1-f%o%powk[-length(powk)],FUN="cumprod",MARGIN=1)
            if (maxmiss==1) nvec<-cbind(rep.int(1,nsim),notfind)*f else nvec<-cbind(rep.int(1,nsim),t(notfind))*f
            pfind.si<-t(t(nvec)*powk) # conditional probability of finding a carcass on the ith search (row) after arrival for given (simulated) searcher efficiency (column)
          } else {
            pfind.si<-f
          }
        } else if (tclvalue(tkvars$SEopt)=='f'){
          #capture.output(pk<-as.array(rstan::sampling(pkmod, .Rvar$pkdat, chains = 1, iter = nsim+100, warmup = 100, verbose=F, show_messages = F))[,1,1:2], file=paste0(.Rvar$datadir,'/NULL'))
          # if pk are fit by field trials, then the pk model is fit when the form opens and is updated when parameters are updated.
          # all that remains is resampling the distribution
          capture.output({
            .Rvar$pkjags <- rjags::jags.model(
              textConnection(pkmod),
              data = .Rvar$pkdat,
              inits = with(.Rvar$pkdat,list(p = X[1]/M[1],k = max(min((X[2]/M[2])/(X[1]/M[1]),.99),.01)))
            )
            update(.Rvar$pkjags, 1000)
            .Rvar$tmppk<-rjags::coda.samples(.Rvar$pkjags, variable.names=c('p','k'), n.iter=nsim)[[1]][,2:1]
            },
            file = paste0(.Rvar$datadir, '/NULL')
          )
          .Rvar$pkstat<-list(phat = mean(.Rvar$tmppk[,1]), CIp = quantile(.Rvar$tmppk[,1],c(0.025, 0.975)), khat = mean(.Rvar$tmppk[,2]), CIk = quantile(.Rvar$tmppk[,2],c(0.025, 0.975)), r = cor(.Rvar$tmppk[,1:2])[2])
          .Rvar$pkres<-array(dim=c(dim(.Rvar$tmppk)[1],.Rvar$pkdat$n+2))
          .Rvar$pkres[,1:2]<-.Rvar$tmppk
          .Rvar$pkres[,3]<-.Rvar$tmppk[,1]
          for (i in 4:(.Rvar$pkdat$n+2)){
            .Rvar$pkres[,i]<-.Rvar$pkres[,i-1]*.Rvar$pkres[,2]
          }

          if (maxmiss == 0) pfind.si<-.Rvar$tmppk[,1]
          if (maxmiss == 1) pfind.si<-cbind(.Rvar$tmppk[,1], (1-.Rvar$tmppk[,1])*.Rvar$tmppk[,2]*.Rvar$tmppk[,1])
          if (maxmiss > 1){
            powk<-array(dim=c(nsim,maxmiss+1))
            powk[,1]<-1
            for (i in 1:maxmiss+1) powk[,i] <- powk[,i-1]*.Rvar$tmppk[,2]
            pfind.si<-.Rvar$tmppk[,1]*powk*cbind(rep(1,nsim),t(apply(1-(.Rvar$tmppk[,1]*powk)[,1:maxmiss], F = cumprod, M = 1)))
          }
        } else {
          tkmessageBox(message='Holy cow! tclvalue(tkvars$SEopt) = ', tclvalue(tkvars$SEopt), '\nAborting SYcalcg() in sy_function.R...')
          return(F)
        }
        intxsearch<-unique(cbind(schedule[,2] - schedule[,1], schedule[,3] - schedule[,2]), MAR=1)
        CPab<-array(dim=c(nsim, 2))
        if (perstype == "hand"){
          if (persistence_distn=="Weibull" | persistence_distn=="Exponential"){
            CPab[,2]<-rlnorm(nsim,sdlog=(log(bupr)-log(blwr))/(2*1.96),meanlog=log(pdb0)) # from a lognormal distribution with mean = log(pdb) and sd = (log(maxb)-log(minb))/4
          } else if (persistence_distn=="Lognormal"){
#            CPab[,2]<-pmax(rnorm(nsim,mean=pdb0,sd=(bupr-blwr)/(2*1.96)),.1)
            CPab[,2]<-rnorm(nsim,mean=pdb0,sd=(bupr-blwr)/(2*1.96))
          } else if (persistence_distn=="Log-Logistic"){
            CPab[,2]<-rlnorm(nsim,meanlog=log(pdb0),sdlog=(log(bupr)-log(blwr))/(2*1.96)) # not entirely satisfactory...doesn't look great for highly skewed distribution of b
          }
          CPab[,1]<-pda0
        } else {
          if (persistence_distn == "Exponential"){
            CPab[,2]<-exp(rnorm(nsim,mean=.Rvar$CPdataPrevious$mod.e$coef, sd=sqrt(.Rvar$CPdataPrevious$mod.e$var[1])))
            CPab[,1]<-pda0
          } else if (persistence_distn == "Weibull"){
            CPparms<-MASS::mvrnorm(nsim,c(.Rvar$CPdata$mod.w$coef[1],log(.Rvar$CPdata$mod.w$scale)), .Rvar$CPdata$mod.w$var)
            CPab[,1]<-1/exp(CPparms[,2]) #shape
            CPab[,2]<-exp(CPparms[,1]) # scale
            i0<-which(CPab[,1]<=0 | (Isam/CPab[,2])^CPab[,1]<1e-320)
            len0<-length(i0)
            while(length(i0)>0){
              CPparms[i0,]<-MASS::mvrnorm(length(i0),c(.Rvar$CPdata$mod.w$coef[1],log(.Rvar$CPdata$mod.w$scale)), .Rvar$CPdata$mod.w$var)
              CPab[i0,1]<-1/exp(CPparms[i0,2])
              CPab[i0,2]<-exp(CPparms[i0,1])
              i0<-which(CPab[,1]<=0 | (Isam/CPab[,2])^CPab[,1]<1e-320)
            }
          } else if (persistence_distn == "Log-Logistic"){
            # NOTE: log-logistic a = shape = 1/mod.ll$scale, b = scale = exp(mod.ll$coef[1])
            CPparms<-MASS::mvrnorm(nsim,c(.Rvar$CPdata$mod.ll$coef[1],log(.Rvar$CPdata$mod.ll$scale)), .Rvar$CPdata$mod.ll$var)
            CPab[,1]<-1/exp(CPparms[,2]) #shape
            CPab[,2]<-exp(CPparms[,1]) # scale
            i0<-which(CPab[,1]<=0)
            len0<-length(i0)
            while(length(i0)>0){
              CPparms[i0,]<-MASS::mvrnorm(length(i0),c(.Rvar$CPdata$mod.ll$coef[1],log(.Rvar$CPdata$mod.ll$scale)), .Rvar$CPdata$mod.ll$var)
              CPab[i0,1]<-1/exp(CPparms[i0,2])
              CPab[i0,2]<-exp(CPparms[i0,1])
              i0<-which(CPab[,1]<=0)
            }
          } else if (persistence_distn == "Lognormal"){
            # NOTE: lognormal a = sdlog^2 = mod.ln$scale^2, b = meanlog = mod.ln$coef[1]
            CPparms<-MASS::mvrnorm(nsim,c(.Rvar$CPdata$mod.ln$coef[1], log(.Rvar$CPdata$mod.ln$scale)), .Rvar$CPdata$mod.ln$var)
            CPab[,1]<- exp(CPparms[,2])^2 #shape
            CPab[,2]<- CPparms[,1] # scale
          }
        }
        ppersu<-ppersist(persistence_distn,t_arrive0=0, t_arrive1=intxsearch[,1], t_search = intxsearch[,1]+intxsearch[,2], pda = CPab[,1], pdb = CPab[,2])
        ppersu[is.na(ppersu)] <- max(ppersu[!is.na(ppersu)])
        # arrivals
        # if uniform arrivals
        if (tclvalue(tkvars$arrfun)=="Uniform" | basicMode){
          arrvec<-(schedule[,2]-schedule[,1])/(nsearch*Isam)
         } else if (tclvalue(tkvars$arrfun)=="Compound" & !basicMode) {
          arrvec<-c(rep(arrvec0[1:(length(arrvec0)-maxmiss-1)],maxmiss+1),rep(arrvec0[(length(arrvec0)-maxmiss):length(arrvec0)],(maxmiss+1):1))
         }
        # add the probabilities
        prob_obs<-numeric(nsim)
        if (maxmiss>0){
          for (i in 1:dim(schedule)[1]){
            prob_obs<-prob_obs + pfind.si[,nmiss[i]+1] * ppersu[which(abs(intxsearch[,1]-(schedule[i,2]-schedule[i,1]))<0.001 & abs(intxsearch[,2]-(schedule[i,3]-schedule[i,2]))<0.001),] * arrvec[i]
          }
        } else {
          for (i in 1:dim(schedule)[1]){
            prob_obs<-prob_obs + pfind.si[nmiss[i]+1] * ppersu[which(abs(intxsearch[,1]-(schedule[i,2]-schedule[i,1]))<0.001 & abs(intxsearch[,2]-(schedule[i,3]-schedule[i,2]))<0.001),] * arrvec[i]
          }
        }
        # g for monitored period
        if (max(prob_obs) * gadj < 1) prob_obs <- prob_obs * gadj
        muB<-mean(prob_obs); sig2B<-var(prob_obs)
        Ba<-muB^2/sig2B*(1-muB)-muB; Bb<-Ba*(1/muB-1)
        BabRaw<-suppressWarnings(MASS::fitdistr(prob_obs,'beta',start=list(shape1=Ba,shape2=Bb),control=list(parscale=c(Ba,Bb)))$estimate)
        BabAnn<-NA
        prob_obs <- prob_obs * ifelse(partial, 1, a)
        muB<-mean(prob_obs); sig2B<-var(prob_obs)
        Ba<-muB^2/sig2B*(1-muB)-muB; Bb<-Ba*(1/muB-1)
        Bab<-suppressWarnings(MASS::fitdistr(prob_obs,'beta',start=list(shape1=Ba,shape2=Bb),control=list(parscale=c(Ba,Bb)))$estimate)
        prob_obs<-prob_obs
        if (!basicMode){
          if (tclvalue(tkvars$arrfun)=="Compound"){
            # g for the entire year
            atemporal<-1-(arrmiss0+arrmissf)
            arrmiss0 <- arrmiss0
            arrmissf <- arrmissf
            atemporal<- atemporal
            prob_obs<-prob_obs * atemporal
            muB<-mean(prob_obs); sig2B<-var(prob_obs)
            Ba<-muB^2/sig2B*(1-muB)-muB; Bb<-Ba*(1/muB-1)
            BabAnn<-suppressWarnings(MASS::fitdistr(prob_obs,'beta',start=list(shape1=Ba,shape2=Bb),control=list(parscale=c(Ba,Bb)))$estimate)
          } else {
            BabAnn <- NA
            atemporal <- NA
            arrmiss0 <- NA
            arrmissf <- NA
          }
        } else {
            atemporal<-toR(tkvars$v)
            prob_obs<-prob_obs*atemporal
            muB<-mean(prob_obs); sig2B<-var(prob_obs)
            Ba<-muB^2/sig2B*(1-muB)-muB; Bb<-Ba*(1/muB-1)
            if (sig2B > 0.00001){
              BabAnn<-suppressWarnings(MASS::fitdistr(prob_obs,'beta',start=list(shape1=Ba,shape2=Bb),control=list(parscale=c(Ba,Bb)))$estimate)
            } else {
              BabAnn<-c(Ba, Bb)
            }
        }
      }
.Rvar$BabAnn<-BabAnn
      if (perstype=='hand'){
        rCP<-rCPgab(persistence_distn, toR(tkvars$pda), toR(tkvars$pdb), Ir)
      } else {
        pdi<-(persistence_distn==c("Exponential", "Weibull", "Log-Logistic", "Lognormal"))
        with(.Rvar$CPdataPrevious,{
          pda <- sum(pdi * c(exp(1/mod.e$coef), 1/mod.w$scale, 1/mod.ll$scale, mod.ln$scale^2))
          pdb <- sum(pdi * c(exp(mod.e$coef), exp(mod.w$coef), exp(mod.ll$coef), mod.ln$coef))
          if (persistence_distn=='Exponential'){
            blwr <- exp(mod.e$coef+qt(0.025,mod.e$df.res)*sqrt(mod.e$var[1]))
            bupr <- exp(mod.e$coef+qt(0.975,mod.e$df.res)*sqrt(mod.e$var[1]))
          } else if (persistence_distn=='Weibull'){
            blwr <- exp(mod.w$coef+qt(0.025,mod.w$df.res)*sqrt(mod.w$var[1]))
            bupr <- exp(mod.w$coef+qt(0.975,mod.w$df.res)*sqrt(mod.w$var[1]))
          } else if (persistence_distn=='Log-Logistic'){
            blwr <- exp(mod.ll$coef+qt(0.025,mod.ll$df.res)*sqrt(mod.ll$var[1]))
            bupr <- exp(mod.ll$coef+qt(0.975,mod.ll$df.res)*sqrt(mod.ll$var[1]))
          } else if (persistence_distn=='Lognormal'){
            blwr <- mod.ln$coef+qt(0.025,mod.ln$df.res)*sqrt(mod.ln$var[1])
            bupr <- mod.ln$coef+qt(0.975,mod.ln$df.res)*sqrt(mod.ln$var[1])
          }
        })
        rCP <- rCPgab(persistence_distn, pda, pdb, Ir)
        pda<-pda
        pdb<-pdb
        blwr <-blwr
        bupr<-bupr
      }
      meanCP<-rCP[1]; rhat<-rCP[2]
      if (basicMode) {arrmiss0<-NA; arrmissf<-NA}
      .Rvar$syresult<-list(
        prob_obs=prob_obs, Bab = Bab, ghat = ifelse(nsim <= 1,ghat, 0), FullYear = FullYear, BabAnn = BabAnn, BabRaw = BabRaw,
        MonitoredPeriod = MonitoredPeriod,  atemporal = toR(tkvars$v), arrmiss0 = arrmiss0, arrmissf = arrmissf, X = toR(tkvars$X), a = toR(tkvars$a))
      if (writeg) writegResults(.Rvar$syresult)
      if (partial == T) {.Rvar$symcWindow$addsymc(); tkdestroy(syModule)}
      return(T)
    },
    writegResults = function(res){
      #### write parameter set and summary results to a data file
      with(.Rvar$singleYearPrevious, {
        with(res, {
          if (tclvalue(tkvars$perstype)=="field"){
            pda<-.Rvar$CPdataPrevious$pda
            pdb<-.Rvar$CPdataPrevious$pdb
            bupr<-.Rvar$CPdataPrevious$bupr
            blwr<-.Rvar$CPdataPrevious$blwr
            persistence_distn<-.Rvar$CPdataPrevious$persistence_distn
            rhat<-rCPgab(persistence_distn, pda, pdb, Ir)[2]
          } else {
            pda<-.Rvar$singleYearPrevious$pda
            pdb<-.Rvar$singleYearPrevious$pdb
            bupr<-.Rvar$singleYearPrevious$bupr
            blwr<-.Rvar$singleYearPrevious$blwr
            persistence_distn<-persistence_distn
            rhat<-rCPgab(persistence_distn, pda, pdb, Ir)[2]
          }
          if (sum(is.na(prob_obs))==0 & !is.null(prob_obs)){
            while(1){ if (sink.number()==0) break else sink() }
            sink(paste0(.Rvar$datadir,"/output"))
            cat("Summary statistics for estimation of detection probability (g)\n")
            cat(paste0(unlist(rep("=",80)),collapse=''))
            cat('\nResults:\n')
            if (tclvalue(tkvars$arrfun) == "Uniform" & !basicMode){
              cat("\nFull year:\n  Full year arrival function not provided. Cannot extrapolate to full year.\n")
            } else if (tclvalue(tkvars$arrfun) != "Uniform" | basicMode){
              cat(paste0("\nFull site for full year ",ifelse(basicMode, '', paste0(', ',FullYear)),'\n'))
              if (BabAnn[1]>0){
                cat(paste0('   Estimated g = ', signif(BabAnn[1]/sum(BabAnn),3),", 95% CI = [", signif(qbeta(0.025,BabAnn[1],BabAnn[2]),3),", ",signif(qbeta(0.975,BabAnn[1],BabAnn[2]),3),"]\n",sep=''))
                cat(paste0("   Fitted beta distribution parameters for estimated g: Ba = ", round(BabAnn[1],4),", Bb = ",round(BabAnn[2],4),"\n"))
              } else {
                cat(paste0('   Estimated g = ', round(ghat,3),"\n",sep=''))
                cat("   NOTE: Variance of ghat extremely small. Assuming g is fixed and known.")
              }
            }
            cat(paste0("\n\nFull site for monitored period, ", MonitoredPeriod,"\n"))
        #    Bab<-singleYearProvisional$Bab
            if (Bab[1]>0){
              cat(paste0("   Estimated g = ", signif(Bab[1]/sum(Bab),3),", 95% CI = [", signif(qbeta(0.025,Bab[1],Bab[2]),3),", ",signif(qbeta(0.975,Bab[1],Bab[2]),3),"]\n"))
              cat(paste0("   Fitted beta distribution parameters for estimated g: Ba = ", round(Bab[1],4),", Bb = ",round(Bab[2],4),"\n"))
              cat(paste0('   Temporal coverage (within year) = ',round(atemporal,3),'\n'))

            } else {
              cat(paste("   Estimated g = ", signif(ghat,3),"\n",sep=''))
              cat("   NOTE: Variance of ghat extremely small. Assuming g is fixed and known.")
            }
            cat('\n\n')
            cat(paste0('Searched area for monitored period, ', MonitoredPeriod,'\n'))
            if (BabRaw[1]>0){
              cat(paste0('   Estimated g = ', round(BabRaw[1]/sum(BabRaw),3),", 95% CI = [", signif(qbeta(0.025,BabRaw[1],BabRaw[2]),3),", ",signif(qbeta(0.975,BabRaw[1],BabRaw[2]),3),"]\n"))
              cat(paste0("   Fitted beta distribution parameters for estimated g: Ba = ", round(BabRaw[1],4),", Bb = ",round(BabRaw[2],4),"\n"))
              cat('\n')
            } else {
              cat(paste0('   Estimated g = ', round(ghat/(atemporal*toR(tkvars$a)),3),"\n"))
              cat('\n')
            }
            cat(paste0(unlist(rep("=",80)),collapse=''))
            cat("\nInput:\n")
            cat(paste0("Search parameters\n"))
            if (SEopt == 'h'){
              cat(paste0("   trial carcasses placed = ",SEn, ", carcasses found = ", SEx,"\n"))
              cat(paste0("   estimated searcher efficiency: p = ",signif(SEx/SEn,3), ", 95% CI = [",signif(qbeta(.025,fba,fbb),3),", ",signif(qbeta(.975,fba,fbb),3),"]\n"))
              cat(paste0("   k = ",k,"\n"))
            } else {
              cat(paste0("   Searcher efficiency trials\n"))
              cat(paste0("     carcasses available: "))
              for (i in 1:.Rvar$pkdat$n) cat(paste0(" ",.Rvar$pkdat$M[i]))
              cat('\n')
              cat(paste0("     carcasses discovered:"))
              for (i in 1:.Rvar$pkdat$n) cat(paste0(" ",.Rvar$pkdat$X[i]))
              cat('\n')
              cat("     searcher efficiency: ")
              cat(paste0("estimated p = ",round(.Rvar$pkstat$phat,3), ", "))
              cat(paste0("estimated k = ", round(.Rvar$pkstat$khat,3)))
              cat(paste0("\n       95% CIs:",
                " p in [", round(.Rvar$pkstat$CIp[1],3), ", ", round(.Rvar$pkstat$CIp[2],3), "]",
                " k in [", round(.Rvar$pkstat$CIk[1],3), ", ", round(.Rvar$pkstat$CIk[2],3), "]"
              ))
              cat('\n')
            }
            cat(paste0("   Search schedule: "))
            if (samtype=="Formula"){
              cat(paste0("Search interval (I) = ", Isam, ", number of searches = ",nsearch,", span = ",Isam*nsearch,"\n"))
            } else {
              cat(days)
              cat('\n')
            }
            cat(paste0("     spatial coverage: ", tclvalue(tkvars$a)))
            cat(paste0("     temporal coverage: ", tclvalue(tkvars$v)))
            cat('\n')
            cat(paste0(unlist(rep("_",80)),collapse=''))
            cat("\nCarcass persistence:\n")
            cat(paste0("   ",persistence_distn," persistence distribution\n"))
            rlwr <- rCPgab(persistence_distn, pda, blwr, Ir)[2]
            rupr <- rCPgab(persistence_distn, pda, bupr, Ir)[2]
            if (persistence_distn=="Exponential"){
              persparm.lbl<-paste0("     scale (\u03b2) = ",round(pdb,3))
              persparm.lbl<-paste0(persparm.lbl,"\n     95% CI \u03b2 = [", round(blwr,3) ,", ", round(bupr,3),"]")
              persparm.lbl<-paste0(persparm.lbl," and r = ",round(rhat,3)," for Ir = ",Ir, " with 95% CI = [", round(rlwr, 3), ", ", round(rupr, 3),"]")
            } else {
              persparm.lbl<-paste0("     shape (\u03b1) = ",signif(pda,4), " and scale (\u03b2) = ",round(pdb,3))
              persparm.lbl<-paste0(persparm.lbl,"\n     95% CI \u03b2 = [",round(blwr,3) ,", ",round(bupr,3),"]\n")
              persparm.lbl<-paste0(persparm.lbl,"     r = ",round(rhat,3), " for Ir = ",Ir, " with 95% CI = [", round(rlwr, 3), ", ", round(rupr, 3),"]")
            }
            persparm.lbl<-paste0(persparm.lbl,"\n     ",ifelse(tclvalue(tkvars$perstype)=="field", paste0("n = ", dim(.Rvar$CPdata$CP)[1]), "Parameters entered manually"))
            msg<-persparm.lbl
            cat(msg); cat('\n')
            arrivals.lbl<-paste0("   ",tclvalue(tkvars$arrfun)," arrivals")
            if (tclvalue(tkvars$arrfun)=="Uniform"){
              arrivals.lbl<-paste0(arrivals.lbl,'\n')
            } else {
              arrivals.lbl<-paste0(arrivals.lbl,' with ', sum(arrcomponents),' components:\n')  # this needs to be filled in...
              if (arrcomponents[1]) arrivals.lbl<-paste0(arrivals.lbl,'     uniform arrivals (', round(wt.u/sum(c(wt.u,wt.p1,wt.p2)*arrcomponents),3),') from ', format(as.Date(arrstart+lwr.u,origin="1970-01-01"),"%b %d"),
                ' to ', format(as.Date(arrstart+upr.u,origin="1970-01-01"),"%b %d"),'\n')
              if (arrcomponents[2]) arrivals.lbl<-paste0(arrivals.lbl, '     beta pulse (',round(wt.p1/sum(c(wt.u,wt.p1,wt.p2)*arrcomponents),3),')  ',
                format(as.Date(arrstart+lwr.p1,origin="1970-01-01"),"%b %d"), ' to ',
                format(as.Date(arrstart+upr.p1,origin="1970-01-01"),"%b %d"), ' with a = ',signif(a.p1,4),' and b = ',signif(b.p1,4),'\n')
              if (arrcomponents[3]) arrivals.lbl<-paste0(arrivals.lbl, '     beta pulse (',round(wt.p2/sum(c(wt.u,wt.p1,wt.p2)*arrcomponents),3),') from ',
                format(as.Date(arrstart+lwr.p2,origin="1970-01-01"),"%b %d"), ' to ',
                format(as.Date(arrstart+upr.p2,origin="1970-01-01"),"%b %d"), ' with a = ',signif(a.p2,4),' and b = ',signif(b.p2,4),'\n')
            }
            cat(arrivals.lbl)
            cat(paste0(unlist(rep("_",80)),collapse=''))
          } else {
            while(1){ if (sink.number()==0) break else sink() }
            sink(paste0(.Rvar$datadir,"/output"))
            cat(paste0(unlist(rep("=",80)),collapse=''))
            sink()
            file.show(paste0(.Rvar$datadir,"/output"),delete.file=T,title="Estimated detection probability (g)")
          }
            sink()
            file.show(paste0(.Rvar$datadir,"/output"),delete.file=T,title="Estimated detection probability (g)")
        })
      })
    },
    writePostResults = function(){
      #### write parameter set and summary results to a data file
      with(.Rvar$singleYearPrevious,{
        with(.Rvar$syresult,{
          if (tclvalue(tkvars$perstype)=="field"){
            pda<-.Rvar$CPdataPrevious$pda
            pdb<-.Rvar$CPdataPrevious$pdb
            bupr<-.Rvar$CPdataPrevious$bupr
            blwr<-.Rvar$CPdataPrevious$blwr
            persistence_distn<-.Rvar$CPdataPrevious$persistence_distn
            rhat<-rCPgab(persistence_distn, pda, pdb, Ir)[2]
          } else {
            pda<-.Rvar$singleYearPrevious$pda
            pdb<-.Rvar$singleYearPrevious$pdb
            bupr<-.Rvar$singleYearPrevious$bupr
            blwr<-.Rvar$singleYearPrevious$blwr
            persistence_distn<-persistence_distn
            rhat<-rCPgab(persistence_distn, pda, pdb, Ir)[2]
          }
          if (sum(is.na(prob_obs))==0 & !is.null(prob_obs)){
            while(1){ if (sink.number()==0) break else sink() }
            sink(paste0(.Rvar$datadir,"/output"))
            cat("Summary statistics for fatality estimation (M)\n")
            cat(paste0(unlist(rep("=",80)),collapse=''))
            cat("\nResults:\n")
            cat(paste0("\nCarcasses discovered: X = ",X,'\n'))
            if (tclvalue(tkvars$arrfun)=="Uniform" & !basicMode){
              cat("\nResults for full year:\n  Full year arrival function not provided. Cannot extrapolate to full year.\n")
            } else {
            cat("\nFull site for full year", ifelse (basicMode, "", paste0(", ",FullYear)),"\n")
        #        Bab<-singleYearProvisional$BabAnn
              if (tclvalue(sides) == '1'){
                cat(paste0("   M* = ",min(which(cumsum(pMgX.ann)>crlev))-1," for 1 - alpha = ",crlev,"\n"))
              } else {
                lwrbnd<-min(which(cumsum(pMgX.ann)> (1-crlev)/2))-1
                lwrArea<-ifelse(lwrbnd == 0, 0, sum(pMgX.ann[1:lwrbnd]))
                cat(paste0("   ",100*crlev,"% CI for M = [", lwrbnd,", ", min(which(cumsum(pMgX.ann) > crlev + lwrArea))-1,"], median = ", min(which(cumsum(pMgX.ann) > 0.5))-1, "\n"))
              }
              if (BabAnn[1]>0){
                cat(paste0('   Estimated g: ', signif(BabAnn[1]/sum(BabAnn),3),", 95% CI = [", signif(qbeta(0.025,BabAnn[1],BabAnn[2]),3),", ",signif(qbeta(0.975,BabAnn[1],BabAnn[2]),3),"]\n",sep=''))
                cat(paste0("   Fitted beta distribution parameters for estimated g: Ba = ", round(BabAnn[1],4),", Bb = ",round(BabAnn[2],4),"\n"))
              } else {
                cat(paste0('   Estimated g: ', signif(.Rvar$syresult$ghat,3),"]\n",sep=''))
              }
              cat(paste0('   Temporal coverage (within year) = ',round(atemporal,3),'\n'))
            }
            cat(paste0("\nFull site for monitored period, ",MonitoredPeriod,"\n"))
            if (tclvalue(sides) == '1'){
              cat(paste0("   M* = ",min(which(cumsum(pMgX)>crlev))-1," for 1 - alpha = ",crlev,"\n"))
            } else {
              lwrbnd<-min(which(cumsum(pMgX)> (1-crlev)/2))-1
              lwrArea<-ifelse(lwrbnd == 0, 0, sum(pMgX[1:lwrbnd]))
              cat(paste0("   ",100*crlev,"% CI for M = [", lwrbnd,", ",min(which(cumsum(pMgX) > crlev + lwrArea))-1,"], median = ", min(which(cumsum(pMgX) > 0.5))-1,"\n"))
            }
        #    Bab<-singleYearProvisional$Bab
            if (Bab[1]>0){
              cat(paste("   Estimated g: ", signif(Bab[1]/sum(Bab),3),", 95% CI = [", signif(qbeta(0.025,Bab[1],Bab[2]),3),", ",signif(qbeta(0.975,Bab[1],Bab[2]),3),"]\n",sep=''))
              cat(paste0("   Fitted beta distribution parameters for estimated g: Ba = ", round(Bab[1],4),", Bb = ",round(Bab[2],4)))
              cat("\n")
            } else {
              cat(paste("   Estimated g: ", signif(ghat,3),", 95% CI = [NA, NA]\n",sep=''))
              cat("   Fitted beta distribution parameters for estimated g: Ba = NA, Bb = NA")
              cat("\n   NOTE: Variance of ghat extremely small. Assuming g is fixed and known.")
              cat("\n")
            }
            cat('\n')
            cat(paste0('Searched area for monitored period, ', MonitoredPeriod,'\n'))
            if (tclvalue(sides) == '1'){
              cat(paste0("   M* = ",min(which(cumsum(pMgX.raw)>crlev))-1," for 1 - alpha = ",crlev,"\n"))
            } else {
              lwrbnd<-min(which(cumsum(pMgX.raw)> (1-crlev)/2))-1
              lwrArea<-ifelse(lwrbnd == 0, 0, sum(pMgX.raw[1:lwrbnd]))
              cat(paste0("   ",100*crlev,"% CI for M = [", lwrbnd,", ",min(which(cumsum(pMgX.raw) > crlev + lwrArea))-1,"], median = ", min(which(cumsum(pMgX.raw) > 0.5))-1,"\n"))
            }
            if (Bab[1]>0){
              cat(paste0('   Estimated g = ', round(BabRaw[1]/sum(BabRaw),3),", 95% CI = [", signif(qbeta(0.025,BabRaw[1],BabRaw[2]),3),", ",signif(qbeta(0.975,BabRaw[1],BabRaw[2]),3),"]\n"))
              cat(paste0("   Fitted beta distribution parameters for estimated g: Ba = ", round(BabRaw[1],4),", Bb = ",round(BabRaw[2],4),"\n"))
            } else {
              cat(paste0('   Estimated g = ', round(.Rvar$syresult$ghat,3),"\n"))
              cat(paste0("   Fitted beta distribution parameters for estimated g: Ba = NA, Bb = NA\n"))
            }
            cat('\n')
            cat(paste0(unlist(rep("_",80)),collapse='')); cat("\n")
            cat("   Posterior distribution of M\n")
            if (tclvalue(tkvars$arrfun)=="Uniform" & !basicMode){
              cat("    m  p(M = m) p(M > m)\n")
              for (m in 1:length(pMgX)){
                cat(sprintf("%5i %8.4f %8.4f\n",m-1,round(pMgX[m],4),round(1-sum(pMgX[1:m]),4)))
              }
            } else {
              cat("         Full site for         Full site for       Searched area for\n")
              cat("           full year         monitored period       monitored period\n")
              cat("    m  p(M = m) p(M > m)     p(M = m) p(M > m)      p(M = m) p(M > m)\n")
              for (m in 1:length(pMgX.ann)){
                cat(sprintf("%5i %8.4f %8.4f     %8.4f %8.4f     %8.4f %8.4f\n",
                  m-1,
                  round(pMgX.ann[m],4),
                  round(1-sum(pMgX.ann[1:m]),4),
                  ifelse(m<=length(pMgX),round(pMgX[m],4),0),
                  ifelse(m<=length(pMgX),round(1-sum(pMgX[1:m]),4),0),
                  ifelse(m<=length(pMgX.raw),round(pMgX.raw[m],4),0),
                  ifelse(m<=length(pMgX.raw),round(1-sum(pMgX.raw[1:m]),4),0))
                )
              }
            }
            cat(paste0(unlist(rep("=",80)),collapse=''))
            cat("\n")
            cat("\nInput:\n")

            cat(paste0("Search parameters\n"))
            if (SEopt == 'h'){
              cat(paste0("   trial carcasses placed = ",SEn, ", carcasses found = ", SEx,"\n"))
              cat(paste0("   estimated searcher efficiency: p = ",signif(SEx/SEn,3), ", 95% CI = [",signif(qbeta(.025,fba,fbb),3),", ",signif(qbeta(.975,fba,fbb),3),"]\n"))
              cat(paste0("   k = ",k,", spatial coverage: a = ",a,"\n"))
            } else {
              cat(paste0("   Searcher efficiency trials\n"))
              cat(paste0("     carcasses available: "))
              for (i in 1:.Rvar$pkdat$n) cat(paste0(" ",.Rvar$pkdat$M[i]))
              cat('\n')
              cat(paste0("     carcasses discovered:"))
              for (i in 1:.Rvar$pkdat$n) cat(paste0(" ",.Rvar$pkdat$X[i]))
              cat('\n')
              cat("     searcher efficiency: ")
              cat(paste0("estimated p = ",round(.Rvar$pkstat$phat,3), ", "))
              cat(paste0("estimated k = ", round(.Rvar$pkstat$khat,3)))
              cat(paste0("\n       95% CIs:",
                " p in [", round(.Rvar$pkstat$CIp[1],3), ", ", round(.Rvar$pkstat$CIp[2],3), "]",
                " k in [", round(.Rvar$pkstat$CIk[1],3), ", ", round(.Rvar$pkstat$CIk[2],3), "]"
              ))
              cat('\n')
            }
            cat(paste0("   Search schedule: "))
            if (samtype=="Formula"){
              cat(paste0("Search interval (I) = ", Isam, ", number of searches = ",nsearch,", span = ",Isam*nsearch,"\n"))
            } else {
              cat(days)
              cat('\n')
            }
            cat(paste0(unlist(rep("_",80)),collapse=''))
            cat("\nCarcass persistence:\n")
            cat(paste0("   ",persistence_distn," persistence distribution"))
            rlwr <- rCPgab(persistence_distn, pda, blwr, Ir)[2]
            rupr <- rCPgab(persistence_distn, pda, bupr, Ir)[2]
            if (persistence_distn=="Exponential"){
              persparm.lbl<-paste0(" with scale (\u03b2) = ",round(pdb,3))
              persparm.lbl<-paste0(persparm.lbl,"\n     95% CI \u03b2 = [", round(blwr,3) ,", ", round(bupr,3),"]")
              persparm.lbl<-paste0(persparm.lbl," and r = ",round(rhat,3)," for Ir = ",Ir, " with 95% CI = [", round(rlwr, 3), ", ", round(rupr, 3),"]")
            } else {
              persparm.lbl<-paste0(" with shape (\u03b1) = ",signif(pda,4), " and scale (\u03b2) = ",round(pdb,3))
              persparm.lbl<-paste0(persparm.lbl,"\n     95% CI \u03b2 = [",round(blwr,3) ,", ",round(bupr,3),"]\n")
              persparm.lbl<-paste0(persparm.lbl,"     r = ",round(rhat,3), " for Ir = ",Ir, " with 95% CI = [", round(rlwr, 3), ", ", round(rupr, 3),"]")
            }
            cat(paste0(persparm.lbl,"\n"))
            arrivals.lbl<-paste0("   ",tclvalue(tkvars$arrfun)," arrivals")
            if (tclvalue(tkvars$arrfun)=="Uniform"){
              arrivals.lbl<-paste0(arrivals.lbl,'\n')
            } else {
              arrivals.lbl<-paste0(arrivals.lbl,'\n')
        #      arrivals.lbl<-paste0(arrivals.lbl,' with a = ',afa,' and b = ',afb,"\n")
            }
            cat(arrivals.lbl)
            cat(paste0(unlist(rep("_",80)),collapse=''))
            cat("\nOther\n")
            cat(paste0("  Integrated reference prior for binomial detection probability\n"))
            cat(paste0("    p(M = m) proportional to sqrt(m+1)-sqrt(m)\n"))
            cat(paste0("    Prior distribution truncated at m = ", length(prior_M)-1,'\n'))
            cat("\n")
          } else {
            while(1){ if (sink.number()==0) break else sink() }
            sink(paste0(.Rvar$datadir,"/output"))
            cat(paste0(unlist(rep("=",80)),collapse=''))
            sink()
            file.show(paste0(.Rvar$datadir,"/output"),delete.file=T,title="Fatality estimation (M)")
            return(0)
          }
        })
      })
      sink()
      file.show(paste0(.Rvar$datadir,"/output"),delete.file=T,title="Fatality estimation (M)")
    },
    printParms = function(parms=singleYearPrevious,filename=''){
        while(1){ if (sink.number()==0) break else sink() }
        sink(paste0(.Rvar$datadir,"/output"))
        cat(paste0(filename,'\n'))
        print(parms)
        sink()
        file.show(paste0(.Rvar$datadir,"/output"),delete.file=T,title="Parameter set")
    },
    syClose = function(){
      if (syChkAll()) updatePrevious()
      tkdestroy(syModule)
      if (!partial) tkwm.deiconify(.Rvar$EoA)
    },
    plotPost = function(){
      fullYearBarcolor1<-colors()[123]
      fullYearBarcolor2<-  colors()[257]
      monitoredBarcolor1<-"red"
      monitoredBarcolor2<- colors()[641]
        # This function plots the pmf of M, with red bars covering at least crlev and black bars covering at most (1-crlev).
        # The distribution of M is passed as a vector pMemgX, with p.m.f. given as P(M = 0), P(M = 1), etc.
        # Included in the graph labels are a description of the active parameter set.
        # This function is called from Excel after loading the active parameter set and posterior distribution to R.
        # Open graphing window
      cexlab<-96/(par('cra')/par('cin'))[2]
      with(.Rvar$singleYearPrevious, {
        with(.Rvar$syresult,{
          graphics.off()
          compound <- (tclvalue(tkvars$arrfun)=="Compound")
          if (!compound & !basicMode){ # graph is just like EoA 1's (with slight changes in labels)
             wd<-6.5*.Rvar$charSizeAdjust; ht<-7*.Rvar$charSizeAdjust
              if (.Rvar$platform == 'windows') windows.options(width = wd, height = ht)
              if (.Rvar$platform == 'mac') quartz.options(width = wd, height = ht)
              if (.Rvar$platform == 'linux') X11.options(width = wd, height = ht)
             dev.new(noRStudioGD = T)
             par(mar=c(4,4,4,.5),mgp=c(2.5,.7,0)) #@
          } else { # need to have a bigger graphics window to accommodate three figs
             wd<-8*.Rvar$charSizeAdjust; ht<-10*.Rvar$charSizeAdjust
              if (.Rvar$platform == 'windows') windows.options(width = wd, height = ht)
              if (.Rvar$platform == 'mac') quartz.options(width = wd, height = ht)
              if (.Rvar$platform == 'linux') X11.options(width = wd, height = ht)
             dev.new(noRStudioGD = T)
             par(mar=c(4,4,4,.5),mgp=c(2.5,.7,0)) #@
         }
        # Plot posterior distribution for monitored period
          if (!compound & !basicMode){
             m<-1:length(pMgX)-1
             maxplot<-min(length(m),min(which(cumsum(pMgX)>.999)))
             m<-m[1:maxplot]; pMgX<-pMgX[1:maxplot];
             if (tclvalue(sides)=='1') py<-c(1,1-cumsum(pMgX[-length(pMgX)])) else py <- pMgX
             if (length(dev.list()) > 0 && !prod(par()$fig==c(0,1,0,1))){
              wd<-Rvar$charSizeAdjust * 7; ht = .Rvar$charSizeAdjust * 7
                if (.Rvar$platform == 'windows') windows.options(width = wd, height = ht)
                if (.Rvar$platform == 'mac') quartz.options(width = wd, height = ht)
                if (.Rvar$platform == 'linux') X11.options(width = wd, height = ht)
              dev.new(noRStudioGD = T)
             }
             par(fig=c(0,1,.18,1))
             plot(m,py,
              xlab='m',
              ylab=ifelse (tclvalue(sides) == '1', paste0("Prob( M \u2265 m | X = ", X ,")"), paste0("Prob( M = m | X = ", X ,")")),
              xlim = c(ifelse(tclvalue(sides) == '1', 0, m[min(which(cumsum(pMgX) > 0.001))]),max(m))
             )
             if(.Rvar$platform == "windows") bringToTop()

             cs<-cumsum(pMgX)
             if (tclvalue(sides) == '2'){
               col.use=rep('black', length(m))
               lwrbnd<-min(which(cs > (1-crlev)/2))-1
               lwrArea<-ifelse(lwrbnd == 0, 0, cs[lwrbnd])
               uprbnd<-min(which(cs > crlev + lwrArea))-1
               reds<-lwrbnd:uprbnd+1
               col.use[reds]=monitoredBarcolor2
               for(i in 1:length(py)){
                 lines(x=rep(m[i],2), y= c(0,py[i]), col=col.use[i], lwd=3)
               }
             } else {
               col.use=rep("black", length(m))
               reds<-c(0,which(cs < crlev))+1
               col.use[reds]=monitoredBarcolor1
               for(i in 1:length(py)){
                 lines(x=rep(m[i],2), y= c(0,py[i]), col=col.use[i], lwd=3)
               }
             }
             title(paste0("Posterior Distribution of M (full site, monitored period)"), cex.main = .Rvar$charSizeAdjust*1.2)
             mtext(side=3,MonitoredPeriod, cex=.Rvar$charSizeAdjust)
             if (Bab[1]>0){
               mtext(side=3,line=-2,adj=1,cex=.85*.Rvar$charSizeAdjust,paste("Overall detection probability: g\u0302 = ",signif(Bab[1]/sum(Bab),3), "  \n 95% CI = [ ",round(qbeta(0.025,Bab[1],Bab[2]),3),", ",round(qbeta(0.975,Bab[1],Bab[2]),3),"]  ",sep=""))
             } else {
               mtext(side=3,line=-2,adj=1,cex=.85*.Rvar$charSizeAdjust,paste("Overall detection probability: g\u0302 = ",signif(ghat,3),"  ",sep=""))
             }
        #     mtext(side=3,line=-3,adj=1,cex=.85,paste("P(M \u2264 ",  m[cut.off.cr],") \u2265 ", (1-syA)*100, "%  ",sep=''))
             if (tclvalue(sides) == '1'){
               mtext(side=3,line=-3,adj=1,cex=.85*.Rvar$charSizeAdjust,paste0("M* = ",m[max(reds)],", i.e., P(M \u2264 ",  m[max(reds)],") \u2265 ", crlev*100, "%  ",sep=''))
             } else {
               mtext(side=3,line=-3,adj=1,cex=.85*.Rvar$charSizeAdjust,paste0(100*crlev,"% CI for M = [",m[min(reds)], ", ", m[max(reds)],"]  "))
             }
          } else {
            # arrivals and monitoring
            if (!basicMode){
              xlen <- 1000
              xx<-seq(0,364,length=xlen)
              yy.u<-wt.u/(upr.u-lwr.u)*(xx>=lwr.u & xx<=upr.u)
              yy.p1<-numeric(xlen); ind<-which(xx>=lwr.p1 & xx<=upr.p1)
              yy.p1[ind]<-wt.p1*dbeta((xx[ind]-lwr.p1)/(upr.p1-lwr.p1), shape1=a.p1, shape2=b.p1)/(upr.p1-lwr.p1)
              yy.p2<-numeric(xlen); ind<-which(xx>=lwr.p2 & xx<=upr.p2)
              yy.p2[ind]<-wt.p2*dbeta((xx[ind]-lwr.p2)/(upr.p2-lwr.p2), shape1=a.p2, shape2=b.p2)/(upr.p2-lwr.p2)
              yy<-yy.u*arrcomponents[1]+yy.p1*arrcomponents[2]+yy.p2*arrcomponents[3]
              ymax<-max(yy)
              par(mar=c(3,3,2,.5), mgp=c(2,.7,0), cex = .Rvar$charSizeAdjust, tck= -.04) #smaller margins
              par(fig=c(0,1,0.75,1))
              plot(xx,yy,type='n',
                xlab = "Date",ylab = "",
                axes=F, xaxs='i', yaxs='i',ylim=c(0,ymax*1.06), cex.axis=.9*.Rvar$charSizeAdjust,
                cex.lab=.85*.Rvar$charSizeAdjust)
              s0<-as.numeric(as.Date(format(as.Date(firstsearch),"1970-%m-%d"))-as.Date(arrstart,origin="1970-01-01"))
              sf<-s0+max(days)
              polygon(x=c(par('usr')[1],s0,s0,par('usr')[1]),y=c(par('usr')[3],par('usr')[3],par('usr')[4],par('usr')[4]),col=colors()[350],border=NA)
              polygon(x=c(par('usr')[2],sf,sf,par('usr')[2]),y=c(par('usr')[4],par('usr')[4],par('usr')[3],par('usr')[3]),col=colors()[350],border=NA)
              mtext(side=2,line=.7,"Relative\nArrival Rate", cex=.Rvar$charSizeAdjust*0.8)
              if (arrcomponents[1]) lines(xx,yy.u,col=uclr)
              if (arrcomponents[2]) lines(xx,yy.p1,col=p1clr)
              if (arrcomponents[3]) lines(xx,yy.p2,col=p2clr)
              lines(xx,yy,lwd=2)
              ats <- (ats0-arrstart)%%365
              axis(1, at=ats, lab=figlab, cex.axis = 0.9*.Rvar$charSizeAdjust)
              axis(1, at=s0+days, lab=F, tck=.04)
              mtext(side=3,adj=0, "Unmonitored period shown in gray. Search dates indicated by tick marks on interior of x-axis",cex=96/(par('cra')/par('cin'))[2])
              if (arrmiss0 > 0.01) mtext(text=paste0(round(arrmiss0*100,1),"%"),side=3,line=-2,at=(arrstart+s0)/2,cex=.8*.Rvar$charSizeAdjust)
              if (arrmiss0 > 0.01 | arrmissf > 0.01) mtext(text=paste0(round(atemporal*100,1),"%"),side=3,line=-2,at=(s0+sf)/2,cex=.8*.Rvar$charSizeAdjust)
              if (arrmissf > 0.01) mtext(text=paste0(round(arrmissf*100,1),"%"),side=3,line=-2,at=(sf+par('usr')[2])/2,cex=.8*.Rvar$charSizeAdjust)
              box();
            }
            # mortality during monitored period
            pMgX<-c(pMgX,rep(0,length(pMgX.ann)-length(pMgX)))
            m<-1:length(pMgX)-1
            maxplot<-min(length(m),min(which(cumsum(pMgX.ann)>.999)))
            m<-m[1:maxplot]
            pMgX<-pMgX[1:maxplot]
            pMgX.ann<-pMgX.ann[1:maxplot]
            if (tclvalue(sides)=='1') py<-c(1,1-cumsum(pMgX[-length(pMgX)])) else py <- pMgX
            par(mar=c(3,3,1.5,.5), mgp=c(2,.7,0), tck=-.03) #smaller margins
            if(!basicMode) {
              suppressWarnings(par(fig=c(0,1,.45,.75), new=T, cex = .Rvar$charSizeAdjust))
            } else {
              suppressWarnings(par(fig=c(0,1,.15,.57), new=T, cex = .Rvar$charSizeAdjust))
            }
            plot(m,py,
              xlab='m',
              ylab=ifelse (tclvalue(sides) == '1', paste0("Prob( M \u2265 m | X = ", X ,")"), paste0("Prob( M = m | X = ", X ,")")),
              cex.axis = 0.9,
              cex.lab = 0.9,
              xlim = c(ifelse(tclvalue(sides)=='1', 0, m[min(which(cumsum(pMgX)>0.001))]),max(m))
            )
            if(.Rvar$platform == "windows") bringToTop()

            cs<-cumsum(pMgX)
            if (tclvalue(sides) == '2'){
              col.use=rep('black', length(m))
               lwrbnd<-min(which(cs > (1-crlev)/2))-1
               lwrArea<-ifelse(lwrbnd == 0, 0, cs[lwrbnd])
               uprbnd<-min(which(cs > crlev + lwrArea))-1
               reds<-lwrbnd:uprbnd+1
              col.use[reds]=monitoredBarcolor2
              for(i in 1:length(py)){
                lines(x=rep(m[i],2), y= c(0,py[i]), col=col.use[i], lwd=3)
              }
            } else {
              col.use=rep("black", length(m))
              reds<-c(1, which(cs<crlev) + 1)
              col.use[reds]=monitoredBarcolor1
              for(i in 1:length(py)){
                lines(x=rep(m[i],2), y= c(0,py[i]), col=col.use[i], lwd=3)
              }
            }
            mtext(side=3,adj=0, paste0("Posterior Distribution of M (full site, monitored period), ",MonitoredPeriod), cex=cexlab)
            if (Bab[1]>0){
              mtext(side=3,line=-2,adj=1,cex=.85*cexlab,paste("Overall detection probability: g\u0302 = ",signif(Bab[1]/sum(Bab),3), "  \n 95% CI = [ ",round(qbeta(0.025,Bab[1],Bab[2]),3),", ",round(qbeta(0.975,Bab[1],Bab[2]),3),"]  ",sep=""))
            } else {
              mtext(side=3,line=-2,adj=1,cex=.85*cexlab,paste("Overall detection probability: g\u0302 = ",signif(ghat,3),"  ",sep=""))
            }
            if (tclvalue(sides)=='1'){
               mtext(side=3,line=-3,adj=1,cex=.85*cexlab,paste0("M* = ",m[max(reds)],", i.e., P(M \u2264 ",  m[max(reds)],") \u2265 ", crlev*100, "%  ",sep=''))
            } else {
               mtext(side=3,line=-3,adj=1,cex=.85*cexlab,paste0(100*crlev,"% CI for M = [",m[min(reds)],", ",  m[max(reds)],"]  "))
            }
            # full-year mortality
             cs<-cumsum(pMgX.ann)
            if (tclvalue(sides)=='1') py<-c(1,1-cumsum(pMgX.ann[-length(pMgX.ann)])) else py<-pMgX.ann
             par(mar=c(3,3,1.5,.5), mgp=c(2,.7,0), cex = .Rvar$charSizeAdjust, tck=-.03) #smaller margins
             suppressWarnings(par(fig=c(0,1,.57,1),new=T))
            plot(m,py,
              xlab='m',
              ylab=ifelse(tclvalue(sides) == '1', paste0("Prob( M \u2265 m | X = ", X ,")"), paste0("Prob( M = m | X = ", X ,")")),
              cex.axis=.9,
              cex.lab=.9,
              xlim = c(ifelse(tclvalue(sides)=='1', 0, m[min(which(cumsum(pMgX.ann)>0.001))]),max(m))
            )
             col.use=rep("black", length(m))
            if (tclvalue(sides)=='1'){
               cut.off.cr<-min(which(cs>=crlev))
               col.use[m[1]:(cut.off.cr)]=fullYearBarcolor1
               for(i in 1:length(py)){
                 lines(x=rep(m[i],2), y= c(0,py[i]), col=col.use[i], lwd=3)
               }
            } else {
#              reds<-which(cs<=1-syA/2 & cs > syA/2) + 1
#              reds<-c(min(reds)-1, reds)
               lwrbnd<-min(which(cs > (1-crlev)/2))-1
               lwrArea<-ifelse(lwrbnd == 0, 0, cs[lwrbnd])
               uprbnd<-min(which(cs > crlev + lwrArea))-1
               reds<-lwrbnd:uprbnd+1
               col.use[reds]=fullYearBarcolor2
              for(i in 1:length(py)){
                lines(x=rep(m[i],2), y= c(0,py[i]), col=col.use[i], lwd=3)
              }
            }
            mtext(side=3,adj=0, paste0("Posterior Distribution of M (full site, full year)", ifelse(basicMode, '', FullYear)), cex=96/(par('cra')/par('cin'))[2])
            if (BabAnn[1]>0){
               mtext(side=3,line=-2,adj=1,cex=.85*cexlab,paste("Overall detection probability: g\u0302 = ",signif(BabAnn[1]/sum(BabAnn),3), "  \n 95% CI = [ ",round(qbeta(0.025,BabAnn[1],BabAnn[2]),3),", ",round(qbeta(0.975,BabAnn[1],BabAnn[2]),3),"]  ",sep=""))
             } else {
               mtext(side=3,line=-2,adj=1,cex=.85*cexlab,paste("Overall detection probability: g\u0302 = ",signif(ghat*atemporal,3),"  ",sep=""))
             }
             #text(max(m),.95*top,paste("Black bars cover less than ",(1-crlev)*100, "% of distribution",sep=''),adj=1,cex=.75)
            if (tclvalue(sides) == 1){
              mtext(side=3,line=-3,adj=1,cex=.85*cexlab,paste("M* = ",m[cut.off.cr],", i.e., P(M \u2264 ",  m[cut.off.cr],") \u2265 ", crlev*100, "%  ",sep=''))
            } else {
              mtext(side=3,line=-3,adj=1,cex=.85*cexlab,paste0(100*crlev,"% CI for M = [",m[min(reds)],", ",  m[max(reds)],"]  "))
            }
          }
           #################
           # graph labels
          if (tclvalue(tkvars$perstype)=="field"){
            pda<-.Rvar$CPdataPrevious$pda
            pdb<-.Rvar$CPdataPrevious$pdb
            bupr<-.Rvar$CPdataPrevious$bupr
            blwr<-.Rvar$CPdataPrevious$blwr
            persistence_distn<-.Rvar$CPdataPrevious$persistence_distn
            rhat<-rCPgab(persistence_distn, pda, pdb, Ir)[2]
          } else {
            pda<-.Rvar$singleYearPrevious$pda
            pdb<-.Rvar$singleYearPrevious$pdb
            bupr<-.Rvar$singleYearPrevious$bupr
            blwr<-.Rvar$singleYearPrevious$blwr
            persistence_distn<-persistence_distn
            rhat<-rCPgab(persistence_distn, pda, pdb, Ir)[2]
          }
           pdlab1<-paste("Carcass count: X = ",X)
           if (tclvalue(tkvars$samtype)=='Formula') {
             pdlab4<-paste("Monitoring: I = ", Isam, ", ",nsearch, " searches, span = ",nsearch*Isam,sep="")
           } else {
             pdlab4<-paste("Monitoring: mean(I) = ", mean(days), ", ", length(days), " searches, span = ",max(days),sep="")
           }
           if (tclvalue(tkvars$SEopt)=='f'){
              pdlab5<-substitute(
                list("Search efficiency": hat(p) == x, hat(k)== y),
                list(x = round(.Rvar$pkstat$phat,3), y=round(.Rvar$pkstat$khat,3))
              )
              pdlab6<-paste0("   95% CIs:",
                " p \u2208 [", round(.Rvar$pkstat$CIp[1],3), ", ", round(.Rvar$pkstat$CIp[2],3), "]",
                " k \u2208 [", round(.Rvar$pkstat$CIk[1],3), ", ", round(.Rvar$pkstat$CIk[2],3), "]"
              )
              pdlab12<-paste("Coverage: a = ",a ,", Prior: ",prior_f,sep='')
           } else {
             pdlab5<-paste("Search efficiency: p = ", ifelse(SEopt == 'h', signif(SEx/SEn,3), ),
            ", 95% CI = [", signif(qbeta(0.025,shape1=fba,shape2=fbb),3), ", ",signif(qbeta(0.975,shape1=fba,shape2=fbb),3), "]",sep='')
             pdlab6<-paste0("k = ",k, ", Coverage: a = ",a)
             pdlab12<-paste("Prior: ",prior_f,sep='')
           }
           laba<-paste("Persistence: ", persistence_distn,sep='')
           if (persistence_distn=="Exponential"){
              pdlab8<-substitute(list(.laba,alpha==.pda, r==.r),list(.laba=laba,.pda=signif(pda,3),.r=signif(rhat,3)))
           } else {
              pdlab8<-substitute(list(.laba,alpha==.pda,beta==.pdb, r==.r),list(.laba=laba,.pda=signif(pda,3),.pdb=signif(pdb,3),.r=signif(rhat,3)))
           }
           laba<-paste("Arrivals: ",tclvalue(tkvars$arrfun),sep='')
           if (tclvalue(tkvars$arrfun)=="Uniform"){
              pdlab11<-laba
           } else {
              pdlab11<-paste0(laba,", with ",sum(arrcomponents), " components")
              complab<-''
              if (arrcomponents[1]) complab<-paste0(complab,"  uniform (",wt.u,"), ", format(as.Date(lwr.u+arrstart,origin="1970-01-01"),"%b-%d")," to ",format(as.Date(upr.u+arrstart,origin="1970-01-01"),"%b-%d"))
              if (arrcomponents[2]) complab<-paste0(complab,"\n  beta (",wt.p1,"), ",
                format(as.Date(lwr.p1+arrstart,origin="1970-01-01"),"%b-%d")," to ",
                format(as.Date(upr.p1+arrstart,origin="1970-01-01"),"%b-%d"),
                ", \u03b1 = ",signif(a.p1,3),", \u03b2 = ",signif(b.p1,3))
              if (arrcomponents[3]) complab<-paste0(complab,"\n  beta (",wt.p2,"), ",
                format(as.Date(lwr.p2+arrstart,origin="1970-01-01"),"%b-%d")," to ",
                format(as.Date(upr.p2+arrstart,origin="1970-01-01"),"%b-%d"),
                ", \u03b1 = ",signif(a.p2,3),", \u03b2 = ",signif(b.p2,3))
           }
           par(mar=c(0,0,0,0)) #smaller margins
           suppressWarnings(par(fig=c(0,1,0,ifelse(compound | basicMode, .15, .18)),new=T))
           plot(0,axes=F,xlab='',ylab='',type='n')
           if(.Rvar$platform == "windows") bringToTop()
           lines(par('usr')[1:2],par('usr')[c(4,4)])
           lhs1<-par('usr')[1]+diff(par('usr')[1:2])*.02
           mtext(pdlab1,side=3,line=-1,adj=0,at=lhs1,cex=.75*cexlab)
           mtext(pdlab4,side=3,line=-2,adj=0,at=lhs1,cex=.75*cexlab)
           mtext(pdlab5,side=3,line=-3,adj=0,at=lhs1,cex=.75*cexlab)
           mtext(pdlab6,side=3,line=-4,adj=0,at=lhs1,cex=.75*cexlab)
           mtext(pdlab12,side=3,line=-5,adj=0,at=lhs1,cex=.75*cexlab)
           lhs2<-sum(par('usr')[1:2])/2+diff(par('usr')[1:2])*.04
           mtext(pdlab8,side=3,line=-1,adj=0,at=lhs2,cex=.75*cexlab)
           mtext(pdlab11,side=3,line=-2,adj=0,at=lhs2,cex=.75*cexlab)
          if (tclvalue(tkvars$arrfun) != "Uniform"){
            mtext(paste0("  Temporal coverage = ", round(atemporal,3)),side=3,line=-2.7,adj=0, at=lhs2,cex=.75*cexlab)
            mtext(complab,side=3,line=-3.1,adj=0,at=lhs2,cex=.75*cexlab,padj=1)
          }
        })
      })
    },
    SYcalcLambda = function(){
      if (!SYcalcg(F)) return(F) # this gives Bab as beta parameters for distribution of ghat
      if (tclvalue(tkvars$arrfun)=="Uniform" & !basicMode){
        pBa<-.Rvar$syresult$Bab[1]; pBb<-.Rvar$syresult$Bab[2]
      } else {
        pBa<-.Rvar$syresult$BabAnn[1]; pBb<-.Rvar$syresult$BabAnn[2]
      }
      x<-toR(tkvars$X)
      mmax <- fmmax.ab(x,pBa,pBb)
      ctprob<-0.0001
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
      if (.Rvar$platform == 'windows') windows.options(reset = T)
      if (.Rvar$platform == 'mac') quartz.options(reset = T)
      if (.Rvar$platform == 'linux') X11.options(reset = T)
      graphics.off()
      par(fig=c(0,1,0,1))
      L<-seq(
          optimize(optqL.ab,interval=c(0,mmax),x=x,pBa=pBa,pBb=pBb,mmax=mmax,p=.005)$min,
          optimize(optqL.ab,interval=c(0,mmax+5*sqrt(mmax)),x=x,pBa=pBa,pBb=pBb,mmax=mmax,p=.995)$min,
          length=1000)
      plot(L,Vectorize(pLgX.ab,"L")(L,x=x, pBa=pBa, pBb=pBb, mmax=mmax), type='l',xlab='\u03bb',ylab='Posterior CDF of \u03bb')
      if(.Rvar$platform == "windows") bringToTop()
      cpLgX<-function(L,x,pBa,pBb,mmax) 1-pLgX.ab(L,x,pBa,pBb,mmax)
      q05<-optimize(f=optqL.ab,interval=c(0.00001,Lmax),x=x,pBa,pBb,mmax=mmax,p=.05)$minimum
      q50<-optimize(f=optqL.ab,interval=c(0.00001,Lmax),x=x,pBa,pBb,mmax=mmax,p=.50)$minimum
      q95<-optimize(f=optqL.ab,interval=c(0.00001,Lmax),x=x,pBa,pBb,mmax=mmax,p=.95)$minimum
      lines(rep(q05,2),c(0,.05),lty=3); lines(c(q05,par('usr')[1]),rep(.05,2),lty=3)
      lines(rep(q50,2),c(0,.5),lty=3); lines(c(q50,par('usr')[1]),rep(.5,2),lty=3)
      lines(rep(q95,2),c(0,.95),lty=3); lines(c(q95,par('usr')[1]),rep(.95,2),lty=3)
      text(par('usr')[1]+diff(par('usr')[1:2])*.03,0.05,paste0("5th percentile = ", round(q05,3)),adj=c(0,-.5), cex=.Rvar$charSizeAdjust)
      text(par('usr')[1]+diff(par('usr')[1:2])*.03,0.5,paste0("median = ",round(q50,3)),adj=c(0,-.5), cex=.Rvar$charSizeAdjust)
      text(par('usr')[1]+diff(par('usr')[1:2])*.03,0.95,paste0("95th percentile = ",round(q95,3)),adj=c(0,-.5), cex=.Rvar$charSizeAdjust)
      tmp<-postL.abCI(x, pBa, pBb, toR(tkvars$crlev))
      title("Posterior of \u03bb", family='serif')
      mtext(side=3,text=paste0("mean = ",signif(tmp$mean,4), ", ", signif(100*toR(tkvars$crlev),3),"% CI = [", signif(tmp$CI[1],3), ", ", signif(tmp$CI[2],3), "]"), cex=.Rvar$charSizeAdjust)
      text(par('usr')[1]+diff(par('usr')[1:2])*.95, .1, adj=1, lab=paste0("x = ", x, ", Ba = ", signif(pBa,3), ", Bb = ", signif(pBb, 3)))
      text(par('usr')[1]+diff(par('usr')[1:2])*.95, 0.05, adj=1, lab=paste0("g\u0302 = ", signif(pBa/(pBa+pBb),3), ", 95% CI = [", signif(qbeta(.025, pBa, pBb),3), ", ", signif(qbeta(.975, pBa, pBb),3), ']'))
    },
    updateCPfieldLabel = function(CPdat){
      nsim<-1000
      CPab<-array(dim=c(nsim, 2))
      if (CPdat$persistence_distn == "Exponential"){
        CPr<-rCPgab(CPdat$persistence_distn,a=1/exp(CPdat$mod.e$coef[1]),b=exp(CPdat$mod.e$coef[1]),Ir=Ir)
        hCP.ci<-exp(qnorm(c(0.025, 0.975), mean = CPdat$mod.e$coef, sd = sqrt(CPdat$mod.e$var[1])))
        rlwr<-rCPgab("Exponential", a=1/exp(CPdat$mod.e$coef[1]),b=exp(CPdat$mod.e$coef[1]+qt(.025,CPdat$mod.e$df.res)*sqrt(CPdat$mod.e$var[1])),Ir=Ir)[2]
        rupr<-rCPgab("Exponential", a=1/exp(CPdat$mod.e$coef[1]),b=exp(CPdat$mod.e$coef[1]+qt(.975,CPdat$mod.e$df.res)*sqrt(CPdat$mod.e$var[1])),Ir=Ir)[2]
        hr.ci<-c(rlwr, rupr)
        bci <- exp(CPdat$mod.e$coef+qt(c(0.025, 0.975),CPdat$mod.e$df.res)*sqrt(CPdat$mod.e$var[1]))
        hshape<-1/CPr[1]
        hscale<-CPr[1]
      } else if (CPdat$persistence_distn == "Weibull"){
        CPr<-rCPgab(CPdat$persistence_distn,a=1/CPdat$mod.w$scale,b=exp(CPdat$mod.w$coef[1]),Ir=Ir)
        CPparms<-MASS::mvrnorm(nsim,c(CPdat$mod.w$coef[1],log(CPdat$mod.w$scale)),CPdat$mod.w$var)
        CPab[,1]<-1/exp(CPparms[,2]) #shape
        CPab[,2]<-exp(CPparms[,1]) # scale
        i0<-which(CPab[,1]<=0 | (Isam/CPab[,2])^CPab[,1]<1e-320)
        while(length(i0)>0){
          CPparms[i0,]<-MASS::mvrnorm(length(i0),c(.Rvar$CPdata$mod.w$coef[1],log(.Rvar$CPdata$mod.w$scale)), .Rvar$CPdata$mod.w$var)
          CPab[i0,1]<-1/exp(CPparms[i0,2])
          CPab[i0,2]<-exp(CPparms[i0,1])
          i0<-which(CPab[,1]<=0 | (Isam/CPab[,2])^CPab[,1]<1e-320)
        }
        hCP.ci<-quantile(CPab[,2]*gamma(1+1/CPab[,1]), c(0.025, 0.975))
        bci <- exp(CPdat$mod.w$coef[1]+qt(c(0.025, 0.975),CPdat$mod.w$df.res)*sqrt(CPdat$mod.w$var[1]))
        hshape<-1/CPdat$mod.w$scale
        hscale<-exp(CPdat$mod.w$coef[1])
      } else if (CPdat$persistence_distn == "Log-Logistic"){
        CPr<-rCPgab(CPdat$persistence_distn, a = 1/CPdat$mod.ll$scale,b=exp(CPdat$mod.ll$coef[1]),Ir=Ir)
        CPparms<-MASS::mvrnorm(nsim,c(CPdat$mod.ll$coef[1],log(CPdat$mod.ll$scale)),CPdat$mod.ll$var)
        CPab[,1]<-1/exp(CPparms[,2]) #shape
        CPab[,2]<-exp(CPparms[,1]) # scale
        i0<-which(CPab[,1]<=0)
        while(length(i0)>0){
          CPparms[i0,]<-MASS::mvrnorm(length(i0),c(.Rvar$CPdata$mod.ll$coef[1],log(.Rvar$CPdata$mod.ll$scale)), .Rvar$CPdata$mod.ll$var)
          CPab[i0,1]<-1/exp(CPparms[i0,2])
          CPab[i0,2]<-exp(CPparms[i0,1])
          i0<-which(CPab[,1]<=0)
        }
        llmcp<-numeric(nsim)
        llmcp[CPab[,1]<=1]<-Inf
        ind<-which(CPab[,1]>1)
        llmcp[ind]<-(CPab[ind,2]*pi/CPab[ind,1])/sin(pi/CPab[ind,1])
        hCP.ci<-quantile(llmcp, c(0.025, 0.975))
        bci <- exp(CPdat$mod.ll$coef[1]+qt(c(0.025, 0.975),CPdat$mod.ll$df.res)*sqrt(CPdat$mod.ll$var[1]))
        hshape<-1/CPdat$mod.ll$scale
        hscale<-exp(CPdat$mod.ll$coef[1])
      } else if (CPdat$persistence_distn == "Lognormal"){
        CPr<-rCPgab(CPdat$persistence_distn, a = CPdat$mod.ln$scale^2,b=CPdat$mod.ln$coef[1],Ir=Ir)
        CPparms<-MASS::mvrnorm(nsim,c(CPdat$mod.ln$coef[1],log(CPdat$mod.ln$scale)),CPdat$mod.ln$var)
        CPab[,1]<- exp(CPparms[,2])^2 #shape
        CPab[,2]<- CPparms[,1] # scale
        hCP.ci<-quantile(exp(CPab[,2]+CPab[,1]/2),c(0.025, 0.975))
        bci <- CPdat$mod.ln$coef[1] + qt(c(0.025, 0.975),CPdat$mod.ln$df.res)* sqrt(CPdat$mod.ln$var[1])
        hshape<-CPdat$mod.ln$scale^2
        hscale<-CPdat$mod.ln$coef[1]
      }
      if (CPdat$persistence_distn != "Exponential"){
        rr <- as.vector(ppersist(CPdat$persistence_distn,t_arrive0 = 0, t_arrive1=Ir, t_search = Ir, pda = CPab[,1], pdb = CPab[,2]))
        rr[is.na(rr)] <- max(rr[!is.na(rr)])
        hr.ci<-quantile(rr, probs = c(0.025, 0.975), na.rm=T)
      }
      tclvalue(dynLbl$distr)<-paste0("Distribution: ", CPdat$persistence_distn,
        " with shape (\u03b1) = ", signif(hshape,4)," and scale (\u03b2) = ", signif(hscale, 4))
      tclvalue(dynLbl$distrcpr)<-paste0("r = ", signif(CPr[2],3), " for Ir = ", Ir,
        ", with 95% CIs: r = [",round(hr.ci[1],3),", ", round(hr.ci[2],3), "], \u03b2 = [", round(bci[1],4), ", ",round(bci[2],4),"]")
    },
#################
    updateCPhandLabel = function(){
      if (pdaok & pdbok){
        pda<-toR(tkvars$pda); pdb<-toR(tkvars$pdb)
        rCP<-rCPgab(persistence_distn, pda, pdb, Ir)
        rCPhtext1<<-paste0("r = ",signif(rCP[2],3), " for Ir = ", Ir)
        tkconfigure(rCP.lbl,text=paste0(rCPhtext1, rCPhtext2))
        meanCP<<-rCP[1]; rhat<<-rCP[2]
        if (bminok & bmaxok){
          blwr<-as.numeric(tclvalue(tkvars$blwr))
          bupr<-as.numeric(tclvalue(tkvars$bupr))
          rCPlwr<-rCPgab(persistence_distn, pda, blwr, Ir)
          rCPupr<-rCPgab(persistence_distn, pda, bupr, Ir)
          rCPhtext2<<-paste0(", with 95% CI: r \u2208 [",signif(rCPlwr[2],3),", ",signif(rCPupr[2],3),"]")
          tkconfigure(rCP.lbl,text=paste0(rCPhtext1, rCPhtext2))
        }
      }
    },
#      if (pdaok & pdbok & ifelse(toR(tkvars$samtype)=="Formula", Ichk(toR(tkvars$Isam)),Ichk(getmode(diff(toR(tkvars$days)))))){
#        pda<-as.numeric(tclvalue(tkvars$pda))
#        pdb<-as.numeric(tclvalue(tkvars$pdb))
#        rCP<-rCPgab(persistence_distn, pda, pdb, Ir)
#        rCPhtext1<<-paste0("r = ",signif(rCP[2],3), " for Ir = ", Ir)
#
#        tkconfigure(rCP.lbl,text=paste0(rCPhtext1, rCPhtext2))
#
#        meanCP<<-rCP[1]; rhat<<-rCP[2]
#        if (bminok & bmaxok){
#          blwr<-as.numeric(tclvalue(tkvars$blwr))
#          bupr<-as.numeric(tclvalue(tkvars$bupr))
#          rCPlwr<-rCPgab(persistence_distn, pda, blwr, Ir)
#          rCPupr<-rCPgab(persistence_distn, pda, bupr, Ir)
#          rCPhtext2<<-paste0(", with 95% CI: r \u2208 [",signif(rCPlwr[2],3),", ",signif(rCPupr[2],3),"]")
#          tkconfigure(rCP.lbl,text=paste0(rCPhtext1, rCPhtext2))
#        }
#      } else {
#        rCPhtext1<<-paste0("r = ",NA, " for Ir = ",NA)
#        meanCP<<-NA; rhat<<-NA
#      }
#      if (! (pdaok & pdbok & bminok & bmaxok)) {
#        rCPhtext2<<-paste0(", with 95% CI: r \u2208 [",NA,", ",NA,"]")
#        tkconfigure(rCP.lbl,text=paste0(rCPhtext1, rCPhtext2))
#      }
#      if (persistence_distn=="Exponential"){
#        tkconfigure(pda.edit,state='disabled')
#      }
#    }
##################
    fcp_form = function(){
      fcpModule <<- tktoplevel()
      tkgrab.set(fcpModule)
      tkfocus(fcpModule)
      tkwm.title(fcpModule,paste0("EoA - Evidence of Absence, ",.Rvar$VER," - Carcass Persistence"))
      tkwm.resizable(fcpModule,0,0)
      tkwm.deiconify(fcpModule)

      # assign .Rvar$CPdata variables locally for interim calculations; recollate upon OK-ing
      for (nm in names(.Rvar$CPdata)) assign(nm, .Rvar$CPdata[[nm]])
#      for (nm in names(.Rvar$CPdataPrevious)) assign(nm, CPdataDefault[[nm]])
      coln<-6
      rown<-5
      Ir<-ifelse(tclvalue(tkvars$samtype) == 'Formula', toR(tkvars$Isam), getmode(diff(toR(tkvars$days))))
      colname<-c("Distribution", "\u0394AIC", paste0("     r\n(Ir = ", Ir,")"), "\u03b1", "\u03b2", "95% CI for \u03b2")
      fcpResults<-tclArray()
      for (i in 1:coln) fcpResults[[0,i-1]]<-as.tclObj(colname[i],drop=T)
      fcpResults[[1,0]]<-"Exponential"
      fcpResults[[2,0]]<-"Weibull"
      fcpResults[[3,0]]<-"Log-Logistic"
      fcpResults[[4,0]]<-"Lognormal"
      fcpTopFrame<-tkframe(fcpModule)
      fcpTable<-tcltk2::tk2table(fcpTopFrame,
        rows = rown,
        cols = coln,
        selectmode = "extended",
        variable = fcpResults,
        titlerows = 1,
        resizeborders = "none",
        multiline = 1,
        state='disabled',
        selecttitle = 1,
        flashmode = T,
        flashtime = 2,
        rowseparator='\n',
        colseparator='\t'
      )
      tkconfigure(fcpTopFrame, padx=15,pady=15)
      columnwidths<-c(12, 7, 9, 9, 9, 15)
      for (i in 1:length(columnwidths)) tcl(fcpTable, "width", i-1, columnwidths[i])
      tcl(fcpTable,"height",0, 2)
      buttonFrame<-tkframe(fcpTopFrame)
      #dummyButton<-tkradiobutton(buttonFrame, bg='gray93', fg='gray93', selectcolor='gray93', relief='flat', borderwidth=0, highlightcolor='gray93' )
      dummyButton<-tklabel(buttonFrame)
      tkconfigure(dummyButton, borderwidth=2, height=2)
      tclvalue(tkvars$persistence_distn)<-persistence_distn
      expoButton<-tkradiobutton(buttonFrame, variable=tkvars$persistence_distn, value= 'Exponential')
      weibullButton<-tkradiobutton(buttonFrame, variable=tkvars$persistence_distn, value= 'Weibull')
      llogisButton<-tkradiobutton(buttonFrame, variable=tkvars$persistence_distn, value= 'Log-Logistic')
      lnormButton<-tkradiobutton(buttonFrame, variable=tkvars$persistence_distn, value= 'Lognormal')
      tkconfigure(expoButton, borderwidth=0, width=2, height=1, highlightthickness=0, anchor='w')
      tkconfigure(weibullButton, borderwidth=0, width=2, height=1, highlightthickness=0, anchor='w')
      tkconfigure(llogisButton,borderwidth=0, width=2, height=1, highlightthickness=0, anchor='w')
      tkconfigure(lnormButton, borderwidth=0, width=2, height=1, highlightthickness=0, anchor='w')
      tkgrid(dummyButton,sticky='e')
      tkgrid(expoButton,sticky='e')
      tkgrid(weibullButton,sticky='e')
      tkgrid(llogisButton,sticky='e')
      tkgrid(lnormButton,sticky='e')
      warningColor<-colors()[552]
      cautionColor<-colors()[652]
      OKColor<-colors()[257]
      tcl(fcpTable,"tag","configure", "warning", bg=warningColor, fg = 'white')
      tcl(fcpTable,"tag","configure", "caution", bg=cautionColor, fg = lnormColor)
      tcl(fcpTable,"tag","configure", "recommended", bg = OKColor, fg = 'white')
      tcl(fcpTable,"tag","configure", "expo",bg=expoColor,fg='white')
      tcl(fcpTable,"tag","configure", "weibull",bg=weibullColor,fg='black')
      tcl(fcpTable,"tag","configure", "llogis",bg=llogisColor,fg=lnormColor)
      tcl(fcpTable,"tag","configure", "lnorm",bg=lnormColor, fg='white')
      tcl(fcpTable,"tag","configure", "flash", bg = 'green')#colors()[124])
      tcl(fcpTable,"tag", "rowtag", "expo","1")
      tcl(fcpTable,"tag", "rowtag", "weibull","2")
      tcl(fcpTable,"tag", "rowtag", "llogis","3")
      tcl(fcpTable,"tag", "rowtag", "lnorm","4")

      ## first time the window opens, use .Rvar$CPdataPrevious to fill the form and graph

      fcpFillTable<-function(){
        aic<-c(extractAIC(mod.e)[2], extractAIC(mod.w)[2], extractAIC(mod.ll)[2], extractAIC(mod.ln)[2])
        aic<-aic-min(aic)
        pda<-c(1/exp(mod.e$coef[1]), 1/mod.w$scale, 1/mod.ll$scale, mod.ln$scale^2)
        pdb<-c(exp(mod.e$coef[1]), exp(mod.w$coef[1]), exp(mod.ll$coef[1]), mod.ln$coef[1])
      #  Ir<-max(round(pdb[1],1))
        CPr<-array(dim=c(4,2))
        CPr[1,]<-rCPgab("Exponential", pda[1], pdb[1], Ir)
        CPr[2,]<-rCPgab("Weibull", pda[2], pdb[2], Ir)
        CPr[3,]<-rCPgab("Log-Logistic", pda[3], pdb[3], Ir)
        CPr[4,]<-rCPgab("Lognormal", pda[4], pdb[4], Ir)
        b95<-array(dim=c(4,2))
        b95[1,] <- exp(mod.e$coef+qt(c(0.025, 0.975),mod.e$df.res)*sqrt(mod.e$var[1]))
        b95[2,] <- exp(mod.w$coef[1]+qt(c(0.025, 0.975),mod.w$df.res)*sqrt(mod.w$var[1]))
        b95[3,] <- exp(mod.ll$coef[1]+qt(c(0.025, 0.975),mod.ll$df.res)*sqrt(mod.ll$var[1]))
        b95[4,] <- mod.ln$coef[1] + qt(c(0.025, 0.975),mod.ln$df.res)* sqrt(mod.ln$var[1])
        for (i in 1:4){
          fcpResults[[i,1]]<-round(aic[i], 2)
      #    fcpResults[[i,2]]<-round(CPr[i,1],2)
          fcpResults[[i,2]]<-round(CPr[i,2],3)
          fcpResults[[i,3]]<-signif(pda[i],5)
          fcpResults[[i,4]]<-signif(pdb[i],5)
          fcpResults[[i,5]]<-as.tclObj(paste0("[", signif(b95[i,1],4), ", ", signif(b95[i,2],4), "]"),drop = T)
        }
        if (aic[1] >= 2) {
          tkconfigure(expoButton, bg = warningColor, text="\u2718")
        } else if (aic[1] >= 1) {
          tkconfigure(expoButton, bg = cautionColor, text="\u2757")
        } else {
          tkconfigure(expoButton, bg = OKColor, text="\u2714")
          tclvalue(tkvars$persistence_distn)<-"Exponential"
        } ##########
        if (aic[2] >= 2){
          tkconfigure(weibullButton, bg = warningColor, text="\u2718")
        } else if (aic[2] > 0){
          tkconfigure(weibullButton, bg = cautionColor, text="\u2757")
        } else {
          tkconfigure(weibullButton, bg = OKColor, text="\u2714")
          tclvalue(tkvars$persistence_distn)<-"Weibull"
        } ##########
        if (aic[3] >= 1){
          tkconfigure(llogisButton, bg = warningColor, text="\u2718")
        } else if (aic[3] == 0 && min(aic[-3]) >= 2){
          tclvalue(tkvars$persistence_distn)<-"Log-Logistic"
          tkconfigure(llogisButton, bg = OKColor, text="\u2714")
        } else {
          tkconfigure(llogisButton, bg = cautionColor, text="\u2757")
        } ##########
        if (aic[4] >= 2){
          tkconfigure(lnormButton, bg = warningColor, text="\u2718")
        } else if (aic[4] > 0){
          tkconfigure(lnormButton, bg = cautionColor, text="\u2757")
        } else {
          tkconfigure(lnormButton, bg = OKColor, text="\u2714")
          tclvalue(tkvars$persistence_distn)<-"Lognormal"
        }
        if (aic[3] == 0 & min(aic[-3])<2){
          nextbest<-which(aic[-3]==min(aic[-3]))
          if (aic[-3][nextbest] == 1){
            tkconfigure(expoButton, bg = OKColor, text="\u2714")
            tclvalue(tkvars$persistence_distn)<-"Exponential"
          } else if (aic[-3][nextbest] == 2){
            tkconfigure(weibullButton, bg = OKColor, text="\u2714")
            tclvalue(tkvars$persistence_distn)<-"Weibull"
          } else {
            tkconfigure(lnormButton, bg = OKColor, text="\u2714")
            tclvalue(tkvars$persistence_distn)<-"Lognormal"
          }
        }
      }
      fcpFillTable()
      tclvalue(tkvars$persistence_distn)<-.Rvar$CPdata$persistence_distn
      actionsFrame<-tkframe(fcpModule)
      fcpSave<-tkbutton(actionsFrame, text = "Save results\nto file", command=function(){
        # ask for file name, path, and format
        filename <- tclvalue(tkgetSaveFile(filetypes = "
          {{Windows bitmap} {.bmp}}
          {{Encapsulated Postscript} {.eps}}
          {{Joint Photographic Experts Group} {.jpg .jpeg}}
          {{Portable Document Format} {.pdf}}
          {{Portable Network Graphics} {.png}}
          {{Postscript} {.ps}}
          {{Tagged Image File Format} {.tiff .tif}}
          {{Enhanced metafiles} {.emf .wmf}} ",
          defaultextension = T,  initialdir=.Rvar$csvpath,
          title = "Save"
        ))
        if (filename == "") return(FALSE)
        tmp<-unlist(strsplit(filename,'/'))
        pathname<-paste(tmp[-length(tmp)],collapse='/')
        if (nchar(pathname)>0) .Rvar$csvpath<-pathname
        ext<-tolower(tail(unlist(strsplit(tail(tmp,1),"[.]")),1))
        # draw graph
        if (ext == 'bmp') bmp(filename=filename)
        if (ext == 'eps') {setEPS(); postscript(file = filename)}
        if (ext %in% c('jpeg', 'jpg')) jpeg(filename=filename)
        if (ext == 'pdf') pdf(file=filename)
        if (ext == 'png') png(filename=filename)
        if (ext == 'ps') {setPS(); postscript(file=filename)}
        if (ext %in% c('tif', 'tiff')) tiff(filename=filename)
        if (ext %in% c('emf', 'wmf')) win.metafile(filename=filename)

        par(mar=c(10, 4, 2, 1), mgp=c(2,.7,0), tck=-.015, family = 'sans')
        if (length(dev.list()) > 0 && !prod(par()$fig==c(0,1,0,1))) {
          wd = 7*.Rvar$charSizeAdjust; ht = 7*.Rvar$charSizeAdjust
          if (.Rvar$platform == 'windows') windows.options(width = wd, height = ht)
          if (.Rvar$platform == 'mac') quartz.options(width = wd, height = ht)
          if (.Rvar$platform == 'linux') X11.options(width = wd, height = ht)
          dev.new(noRStudioGD = T)
        }
        plot(survival::survfit(surv~1,data=CP),xlab='Days', ylab = 'Fraction of carcasses persisting', xaxs='i', yaxs='i')
        xx<-seq(par('usr')[1],par('usr')[2],length=1000)
        lines(xx,1-pexp(xx, rate = 1/exp(mod.e$coef)),col=expoColor, lwd=2)
        lines(xx,1-actuar::pllogis(xx,shape=1/mod.ll$scale, scale = exp(mod.ll$coef[1])),col=llogisColor, lwd=2)
        lines(xx,1-plnorm(xx,meanlog=mod.ln$coef[1], sdlog = mod.ln$scale),col=lnormColor, lwd = 2)
        lines(xx,1-pweibull(xx, 1/mod.w$scale, exp(mod.w$coef[1])), col=weibullColor, lwd = 2)
        box()
        legend(x='topright', legend=c("Exponential", "Weibull", "Log-logistic", "Lognormal"),
          lty = 1, lwd = 2, col = c(expoColor, weibullColor, llogisColor, lnormColor))
        title("Persistence Distribution")
        # print label
        par(family='mono')
        just<-par('usr')[1]-0.07*diff(par('usr')[1:2])
        mtext(side=1, line = 4, at = just, adj=0, text = sprintf("%-12s  %6s  %4s     %4s      %4s      %s","Distribution", "\u0394AIC","r", "\u03b1", "\u03b2", "95% CI for \u03b2" ), cex=.Rvar$charSizeAdjust)
        for (i in 1:4) mtext(side=1, line = 4+i, at = just, adj=0, text = sprintf("%-12s %6.2g  %6.3g  %6.4f   %6.4f  %15s",
            toR(fcpResults[[i,0]]),
            toR(fcpResults[[i,1]]),
            toR(fcpResults[[i,2]]),
            toR(fcpResults[[i,3]]),
            toR(fcpResults[[i,4]]),
            toR(fcpResults[[i,5]])), cex=.Rvar$charSizeAdjust
        )
        mtext(side=1, line = 9, at = just, adj=0, text = sprintf("%-15s  %5s   %10s)","", "",paste0("(for Ir = ", Ir)), cex=.8*.Rvar$charSizeAdjust)
        dev.off()
      })
      fcpOK<-tkbutton(actionsFrame, text = "OK", command=function(){
        .Rvar$CPdata$CP<-CP
        .Rvar$CPdata$surv <- surv
        .Rvar$CPdata$persistence_distn<-tclvalue(tkvars$persistence_distn)
        .Rvar$CPdata$mod.e<-mod.e
        .Rvar$CPdata$mod.w<-mod.w
        .Rvar$CPdata$mod.ll<-mod.ll
        .Rvar$CPdata$mod.ln<-mod.ln
        if (.Rvar$CPdata$persistence_distn=='Exponential'){
          .Rvar$CPdata$pda <- 1/exp(mod.e$coef)
          .Rvar$CPdata$pdb <- exp(mod.e$coef)
          .Rvar$CPdata$blwr <- exp(mod.e$coef+qt(0.025,mod.e$df.res)*sqrt(mod.e$var[1]))
          .Rvar$CPdata$bupr <- exp(mod.e$coef+qt(0.975,mod.e$df.res)*sqrt(mod.e$var[1]))
          CPr <- rCPgab("Exponential",.Rvar$CPdata$pda, .Rvar$CPdata$pdb, Ir)
        } else if (.Rvar$CPdata$persistence_distn=='Weibull'){
          .Rvar$CPdata$pda <- 1/mod.w$scale
          .Rvar$CPdata$pdb <- exp(mod.w$coef[1])
          .Rvar$CPdata$blwr <- exp(mod.w$coef+qt(0.025,mod.w$df.res)*sqrt(mod.w$var[1]))
          .Rvar$CPdata$bupr <- exp(mod.w$coef+qt(0.975,mod.w$df.res)*sqrt(mod.w$var[1]))
          CPr <- rCPgab("Weibull", .Rvar$CPdata$pda, .Rvar$CPdata$pdb, Ir)
        } else if (.Rvar$CPdata$persistence_distn=='Log-Logistic'){
          .Rvar$CPdata$pda <- 1/mod.ll$scale
          .Rvar$CPdata$pdb <- exp(mod.ll$coef[1])
          .Rvar$CPdata$blwr <- exp(mod.ll$coef+qt(0.025,mod.ll$df.res)*sqrt(mod.ll$var[1]))
          .Rvar$CPdata$bupr <- exp(mod.ll$coef+qt(0.975,mod.ll$df.res)*sqrt(mod.ll$var[1]))
          CPr <- rCPgab("Log-Logistic", .Rvar$CPdata$pda, .Rvar$CPdata$pdb, Ir)
        } else if (.Rvar$CPdata$persistence_distn=='Lognormal'){
          .Rvar$CPdata$pda <- mod.ln$scale^2
          .Rvar$CPdata$pdb <- mod.ln$coef[1]
          .Rvar$CPdata$blwr <- mod.ln$coef+qt(0.025,mod.ln$df.res)*sqrt(mod.ln$var[1])
          .Rvar$CPdata$bupr <- mod.ln$coef+qt(0.975,mod.ln$df.res)*sqrt(mod.ln$var[1])
          CPr <-rCPgab("Lognormal", .Rvar$CPdata$pda, .Rvar$CPdata$pdb, Ir)
        }
        .Rvar$CPdata$meanCP<-CPr[1]
        .Rvar$CPdata$rhat<-CPr[2]
        .Rvar$CPdata$persistence_distn <- tclvalue(tkvars$persistence_distn)
        updateCPfieldLabel(.Rvar$CPdata)
        tkdestroy(fcpModule)
      })
      fcpCancel<-tkbutton(actionsFrame, text = "Cancel", command=function(){
#        tclvalue(tkvars$perstype)<-'hand'
#        tcl(pdtypeHand,'invoke')
        tkdestroy(fcpModule)
      })
      getCPdataButton<-tkbutton(actionsFrame, text = 'Get data from\nfield trials', command = function() {
        if (getCPdata()) {
          tkrplot::tkrreplot(persFig)
          fcpFillTable()
        } else {
          tkdestroy(fcpModule)
          tclvalue(tkvars$persistence_distn)<-"Exponential"
          tkinvoke(pdtypeHand)
          tkselection.clear(pd.lbox, 1, 3)
          tkselection.set(pd.lbox, which(pdnames=='Exponential')-1)
          tkactivate(pd.lbox, which(pdnames=='Exponential')-1)
          self$persistence_distn <- "Exponential"
          pdb.e<-sclbnd[1]
          tclvalue(tkvars$pdb.e)<-pdb.e
          tclvalue(tkvars$pda)<-signif(1/pdb.e,5); tkconfigure(pda.edit,state='disabled')
          tclvalue(tkvars$pdb)<-tclvalue(tkvars$pdb.e)
          tclvalue(tkvars$blwr)<-tclvalue(tkvars$blwr.e)<-sclbnd[2]
          tclvalue(tkvars$bupr)<-tclvalue(tkvars$bupr.e)<-sclbnd[3]
        }
      })
      tkconfigure(fcpOK,width=12)
      tkconfigure(fcpCancel,width=12)
      tkconfigure(getCPdataButton,width=12)
      tkconfigure(fcpSave, width=12)
      tkgrid(getCPdataButton, padx=5)
      tkgrid(fcpSave, padx=5)
      tkgrid(fcpOK, padx=5)
      tkgrid(fcpCancel, padx=5)
      persFig <- tkrplot::tkrplot(fcpModule, hscale=1.5, vscale = 1.2,  fun = function(){
        par(mar=c(4, 4, 0.5, 0.5), mgp=c(2,.7,0), tck=-.015)
        plot(survival::survfit(surv~1,data=CP),xlab='Days', ylab = 'Fraction of carcasses persisting')
        xx<-seq(par('usr')[1],par('usr')[2],length=1000)
        lines(xx,1-pexp(xx, rate = 1/exp(mod.e$coef)),col=expoColor, lwd=2)
        lines(xx,1-actuar::pllogis(xx,shape=1/mod.ll$scale, scale = exp(mod.ll$coef[1])),col=llogisColor, lwd=2)
        lines(xx,1-plnorm(xx,meanlog=mod.ln$coef[1], sdlog = mod.ln$scale),col=lnormColor, lwd = 2)
        lines(xx,1-pweibull(xx, 1/mod.w$scale, exp(mod.w$coef[1])), col=weibullColor, lwd = 2)
        box()
        mtext(side = 1, line = -1, at = 0, adj = 0, text=paste0("  n = ", summary(mod.e)$n))
      #  legend(x='topright', title = 'Distribution:', title.adj=0.05,
      #  legend = c('Weibull', 'Lognormal', 'Loglogistic', 'Exponential', 'Kaplan-Meyer'),
      #  col = c(colors()[33], 4, 5, 3, 1),
      #  lwd=c(2,2,2,2,1))
      })
      getCPdata<-function(){
        fileName<-tclvalue(tkgetOpenFile(initialdir=.Rvar$csvpath))
        if (!nchar(fileName)) {
          return (FALSE)
        }
        tmp<-unlist(strsplit(fileName,'/'))
        .Rvar$csvpath<-paste(tmp[-length(tmp)],collapse='/')
        #  v<-na.omit(try(read.csv(fileName,sep=",",as.is=T),silent=T)[,1:2]) # debug: <<     production: <
        CP<-try(read.csv(fileName,sep=",",as.is=T))
        # CP data file has two columns: CPmin and CPmax, which bracket the time when the carcass disappeared
        #   if carcass disappears before first check, CPmin = 0
        #   if carcass scavenging event was observed, CPmin = CPmax = time since carcass placement
        #   if carcass does not disappear before end of study, CPmin = time of last check, CPmax = Inf (not "Inf" or inf or "inf" or "Infinity" or 100000 or...)
        if (class(CP) != "try-error") CP<-na.omit(CP) else {tkmessageBox(icon='error',message=paste0("Error in data (",fileName,"\nRequired: two columns with data for CPmin and CPmax.\nCheck file.")); return(F) }
        if (dim(CP)[2] < 2){
          tkmessageBox(icon='error',message=paste0("Error in data (",fileName,"\nRequired: two columns with data for CPmin and CPmax.\nCheck file.")); return(F)
        } else {
          CP<-CP[,1:2]
          names(CP)<-c("CPmin", "CPmax")
        }
        if (!is.numeric(CP[,1]) || !is.numeric(CP[,2]) || sum(CP[,1] < 0) > 0 || sum(CP[,1] > CP[,2]) > 0) {
          tkmessageBox(message = "Error in CP data. Cannot calculate.")
          return(F)
        }

        xind<-which(CP$CPmin == 0 & CP$CPmax == Inf)
        if (length(xind)>0){
          CP$CPmin<-CP$CPmin[-xind]
          CP$CPmax<-CP$CPmax[-xind]
        }

        CP$CPmin <- pmax(CP$CPmin, 0.001)
        event<-ifelse(CP$CPmin == CP$CPmax, 1, ifelse(CP$CPmax == Inf, 0, 3))
        left<-CP$CPmin
        right<-CP$CPmax
        right[event==0]<-CP$CPmin[event==0]

        if (sum(CP$CPmax==Inf)==length(CP$CPmax)){
          n<-length(CP$CPmax)
          trial.period<-mean(CP$CPmin)
          deno<-gsl::hyperg_2F1(1/2, 1/2 - n, 3/2, 1)
          rbnd<-suppressWarnings(c(
            uniroot(function(p) gsl::hyperg_2F1(1/2, 1/2 - n, 3/2, 1-p)*sqrt(1-p)/deno - (1-pexp(1,1)), interval=c(0,1))$root,
            uniroot(function(p) gsl::hyperg_2F1(1/2, 1/2 - n, 3/2, 1-p)*sqrt(1-p)/deno - 0.975, interval=c(0,1))$root,
            uniroot(function(p) gsl::hyperg_2F1(1/2, 1/2 - n, 3/2, 1-p)*sqrt(1-p)/deno - 0.025, interval=c(0,1))$root))
          sclbnd<<-signif(trial.period/-log(rbnd), 3)
          tkmessageBox(message=paste0("No carcasses removed in persistence trials. Use exponential persistence distribution and enter parameters manually.\n",
            "scale = ", signif(sclbnd[1], 3),
            "\nlwr = ", signif(sclbnd[2], 3),
            "\nupr = ", signif(sclbnd[3], 3)))
          return (F)
        }
        surv<- survival::Surv(time=left, time2=right, event=event, type=c('interval'))
        # fit survival models to persistence data and plot
        mod.e<<-survival::survreg(surv~1, dist="exponential")
        # NOTE: exponential b = meanCP = exp(mod.e$coef), variance of estimated mod.e$coef = mod.e$var

        mod.w<<-survival::survreg(surv~1, dist="weibull")
        # NOTE: weibull a = shape = 1/mod.w$scale, b = scale = exp(mod.w$coef[1])

        mod.ll<<-survival::survreg(surv~1, dist="loglogistic")
        # NOTE: log-logistic a = shape = 1/mod.ll$scale, b = scale = exp(mod.ll$coef[1])

        mod.ln<<-survival::survreg(surv~1, dist="lognormal")
        # NOTE: lognormal a = sdlog^2 = mod.ln$scale^2, b = meanlog = mod.ln$coef[1]
        surv<<-surv
        CP<<-CP
        persistence_distn<<-'Weibull'
        pda <<- 1/mod.w$scale
        pdb <<- exp(mod.w$coef[1])
        blwr <<- exp(mod.w$coef+qt(0.025,mod.w$df.res)*sqrt(mod.w$var[1]))
        bupr <<- exp(mod.w$coef+qt(0.975,mod.w$df.res)*sqrt(mod.w$var[1]))
        return(T)
      }
      tkgrid(buttonFrame, fcpTable, actionsFrame)
      tkgrid(fcpTopFrame)
      #tkgrid(IrFrame,getCPdataButton)
      tkgrid(persFig, columnspan = 3)
      tkbind(fcpTable,"<Enter>", function(){
        tkconfigure(fcpTable, cursor = "left_ptr")
      })
    },
    plotarr.mini = function() { # not currently scaled to endpoints

      par(mar=c(0,0,0,0))
      if (tclvalue(tkvars$arrfun)=="Uniform"){
        plot(0,0,type='n',xlim=c(0,364),ylim=c(0,1),yaxs='i',xaxs='i')
        x0<-as.numeric(as.Date(tclvalue(tkvars$firstsearch))-as.Date(paste0(format(as.Date(tclvalue(tkvars$firstsearch)),"%Y"),"-01-01")))
        if (tclvalue(tkvars$samtype)=="Formula"){
          span <- toR(tkvars$nsearch) * toR(tkvars$Isam)
        } else if (tclvalue(tkvars$samtype)=="Custom"){
          span <- max(toR(tkvars$days))
        } else {
          span<-toR(tkvars$span)
        }
        x1<-x0 + span
        lht<-.3
        polygon(c(x0,x0,x1,x1),c(0,lht, lht, 0))
      } else {
        if (sum(arrcomponents)==0){
          tkmessageBox(message="No model components have been defined.\n Click 'Edit' to build a model or select 'Uniform' for simple model.")
        } else {         #  otherwise... compound arrival function
          arrstart<-toR(tkvars$arrstart);# arrcomponents<-toR(tkvars$arrcomponents)
          lwr.u<-toR(tkvars$lwr.u); upr.u<-toR(tkvars$upr.u); wt.u<-toR(tkvars$wt.u)
          lwr.p1<-toR(tkvars$lwr.p1); upr.p1<-toR(tkvars$upr.p1); wt.p1<-toR(tkvars$wt.p1); a.p1<-toR(tkvars$a.p1); b.p1<-toR(tkvars$b.p1)
          lwr.p2<-toR(tkvars$lwr.p2);upr.p2<-toR(tkvars$upr.p2);wt.p2<-toR(tkvars$wt.p2); a.p2<-toR(tkvars$a.p2); b.p2<-toR(tkvars$b.p2)
          xlen <- 1000
          xx<-seq(0,364,length=xlen)
          yy.u<-wt.u/(upr.u-lwr.u)*(xx>=lwr.u & xx<=upr.u)
          yy.p1<-numeric(xlen); ind<-which(xx>=lwr.p1 & xx<=upr.p1)
          yy.p1[ind]<-wt.p1*dbeta((xx[ind]-lwr.p1)/(upr.p1-lwr.p1), shape1=a.p1, shape2=b.p1)/(upr.p1-lwr.p1)
          yy.p2<-numeric(xlen); ind<-which(xx>=lwr.p2 & xx<=upr.p2)
          yy.p2[ind]<-wt.p2*dbeta((xx[ind]-lwr.p2)/(upr.p2-lwr.p2), shape1=a.p2, shape2=b.p2)/(upr.p2-lwr.p2)
          yy<-yy.u*arrcomponents[1]+yy.p1*arrcomponents[2]+yy.p2*arrcomponents[3]
          ymax<-max(yy)
          plot(0,0,type='n',xlim=c(0,364),ylim=c(0,ymax*1.06),yaxs='i',xaxs='i')
          lines(xx,yy)
        }
      }
    },
    sy2rds = function(){
      junk <- list()
      for (nm in names(tkvars)) {
        junk[[nm]]<-toR(tkvars[[nm]])
      }
      if (junk$samtype=='Formula'){
        junk$Isam<-toR(tkvars$Isam)
        junk$nsearch<-toR(tkvars$nsearch)
        junk$days<-.Rvar$singleYearPrevious$days
      } else {
        junk$Isam<-.Rvar$singleYearPrevious$Isam
        junk$nsearch<-.Rvar$singleYearPrevious$nsearch
        junk$days<-toR(tkvars$days)
      }
#      if (exists("arrcomponents", env = .Rvar) && is.na(.Rvar$arrcomponents)) junk$arrcomponents<-c(T,T,T) else junk$arrcomponents<-.Rvar$arrcomponents
      junk$prior_f <- prior_f
      junk$prior_M <- prior_M
      junk$custom.prior<-ifelse(prior_f=="Custom", prior_M, .Rvar$singleYearPrevious$custom.prior)
      junk$objective.prior <- ifelse(prior_f=="Objective", prior_M, .Rvar$singleYearPrevious$objective.prior)
      junk$persistence_distn <- persistence_distn
      return(junk)
    },
    pkdatEnter = function(pkdat){
      pkModule <<- tktoplevel()
      pkData<-tclArray()
      tkwm.title(pkModule,paste0("EoA, v",.Rvar$VER," - Enter p-k Data"))
      tkgrab.set(pkModule)
      tkfocus(pkModule)
      tkwm.resizable(pkModule,0,0)
      tkwm.deiconify(pkModule)

      pkData[[0,0]]<-as.tclObj("\nSearch",drop=T)
      pkData[[0,1]]<-as.tclObj("Carcasses\nAvailable",drop=T)
      pkData[[0,2]]<-as.tclObj("Carcasses\nDiscovered",drop=T)
      for (i in 1:pkdat$n){
        for (j in 0:2){
          pkData[[i,0]]<-i
          pkData[[i,1]]<-pkdat$M[i]
          pkData[[i,2]]<-pkdat$X[i]
        }
      }
      setT<-function(tabdat,rowm,colm) tcl(tabdat,"tag", "celltag", "cellok",as.tclObj(paste0(rowm,',', colm)))  #
      setF<-function(tabdat,rowm,colm) tcl(tabdat,"tag", "celltag", "error",as.tclObj(paste0(rowm,',', colm)))
      pkCellChk<-function(val, rowm, colm){
        val<-suppressWarnings(as.numeric(val))
        if (length(val)==0 || is.na(val)){
          setF(pkTable, rowm, colm)
          return(F)
        }
        if (val < 0 || val != round(val)){
          setF(pkTable, rowm, colm)
          return(F)
        } else {
          setT(pkTable, rowm, colm)  # tag the cell
          return(T)
        }
      }
      pkTableChk<-function(){
        goodtab<-T
        for (i in 1:as.numeric(tclvalue(tcl(pkTable,"index","end","row")))){
          if (!pkCellChk(tclvalue(pkData[[i,1]]), i, 1)){
            tcl(pkTable, "tag", "celltag", "error", as.tclObj(paste0(i, ',', 1)))
            goodtab<-F
          }
          if (!pkCellChk(tclvalue(pkData[[i,2]]), i, 2)) {
            tcl(pkTable, "tag", "celltag", "error", as.tclObj(paste0(i, ',', 2)))
            goodtab<-F
          }
          if (toR(pkData[[i,1]]) < toR(pkData[[i,2]])){
            tcl(pkTable, "tag", "celltag", "adderr", as.tclObj(paste0(i, ',', 1)))
            tcl(pkTable, "tag", "celltag", "adderr", as.tclObj(paste0(i, ',', 2)))
            goodtab<-F
          }
        }
        return(goodtab)
      }
      valCharpk<-function(S){
        actcol<-as.numeric(tclvalue(tkindex(pkTable,"active","col")))
        actrow<-as.numeric(tclvalue(tkindex(pkTable,"active","row")))
        if (length(grep("\n",S))>0){ #
          return(tcl("expr", FALSE))
        } else {
          pkCellChk(S, actrow, actcol)
          return(tcl("expr", TRUE))# but other kinds of space are not--> error-checking if good value is added
        }
      }

      pkTable<-tcltk2::tk2table(pkModule,
        rows=length(pkdat$M)+1,
        cols=3,
        selectmode="extended",
        variable=pkData,
        resizeborders="none",
        titlerows = 1,
        titlecols = 1,
        rowseparator='\n',
        colseparator='\t',
        background='white',
        validate = 1,
        vcmd=function(S) valCharpk(S),
        selecttitle = 1
      )
      colwidths<-c(7, 10, 10)
      for(i in 1:3) {
        tcl(pkTable, "width", i - 1, colwidths[i])
      }
      tcl(pkTable, "height", 0, 2)
      ## state configurations
      tcl(pkTable, "tag", "configure", "error",bg=colors()[652]) # cells with error have yellow background
      tcl(pkTable, "tag", "configure", "adderr", bg=colors()[400]) # cells that have valid format but incompatible with other values
      tcl(pkTable, "tag", "configure", "cellok",bg='white')
      tcl(pkTable, "tag", "configure", "active",fg='black',relief='groove')

      #### table features
      # disable standard paste command for the table
      # tkevent.delete("<<Paste>>","<Control-v>") # reactivate when window is destroyed
      # add row
      tkbind(pkTable,"<Control-KeyPress-a>", function(){
      # add a row
        tkinsert(pkTable, "rows", tclvalue(tkindex(pkTable, "active", "row")), 1)
        tkevent.generate(pkTable,"<KeyPress-Down>")
        rowind<-as.numeric(tclvalue(tcl(pkTable,"index","active","row")))
        for (i in 1:2){
          tcl(pkTable,"tag","celltag","error", as.tclObj(paste0(rowind,',', i)))
        }
      # update the row index labels
        for (i in 1:as.numeric(tclvalue(tcl(pkTable,"index","end","row")))) pkData[[i,0]]<-i
      })
      # add row if 'return' on last row
      tkbind(pkTable,"<Return>",function(){
        if (tclvalue(tkindex(pkTable,"active","row")) == tclvalue(tkindex(pkTable,"end","row"))){
          tkevent.generate(pkTable,"<Control-KeyPress-a>")
          return(T)
        } else {
          tkevent.generate(pkTable,"<KeyPress-Down>")
        }
      })
      # delete row
      tkbind(pkTable,"<Control-KeyPress-d>", function(){
      # delete current row (if there is only one row, just erase the data)
        actrow<-tclvalue(tkindex(pkTable, "active", "row"))
        if (as.numeric(actrow)>1) tkdelete(pkTable, "rows", actrow, 1)
      # update the row index labels
        for (i in 1:as.numeric(tclvalue(tcl(pkTable,"index","end","row")))) pkData[[i,0]]<-i
      })
      tkbind(pkTable,"<Control-Key-v>", function(){
      # after pasting data from the clipboard, read the prior table data into a buffer
      # if there are errors in the new data, give error message and replace
        junk<-readClipboard()
        ncols<-max(nchar(gsub('[^\t]','',junk)))
        if (ncols+as.numeric(tclvalue(tcl(pkTable,"index","active","col"))) > 2) {
          tkmessageBox(message="Error: would-be paste region has too many columns. Cannot paste",icon='error')
          return(F)
        }
        rown <- as.numeric(tclvalue(tkindex(pkTable,"end","row")))
        tkconfigure(pkTable,rows=max(rown+1,as.numeric(tclvalue(tcl(pkTable,"index","active","row")))+length(junk)))
        tkevent.generate(pkTable,"<<Paste>>")
        # check table for errors and missing values...
        pkTableChk()
      })
      # check cell values
      # check table
      pkButFrame<-ttkframe(pkModule)
      pkView<-tkbutton(pkButFrame, text = "View", width=6, command = function(){
        if (!pkTableChk()){
          tkmessageBox(message="Error in data. Cannot calculate...")
          return(F)
        }
        # fix dimensions of pkdat array
        junk<-list()
        # add new data to pkdat
        junk$n<-as.numeric(tclvalue(tcl(pkTable,"index","end","row")))
        junk$X<-numeric(junk$n)
        junk$M<-numeric(junk$n)
        for (i in 1:junk$n){
          junk$M[i]<-toR(pkData[[i,1]])
          junk$X[i]<-toR(pkData[[i,2]])
        }
        pkdat<-junk
        capture.output({
          .Rvar$pkjags <- rjags::jags.model(
            textConnection(pkmod),
            data = pkdat,
            inits = with(pkdat,list(p = X[1]/M[1],k = max(min((X[2]/M[2])/(X[1]/M[1]),.99),.01)))
          )
          update(.Rvar$pkjags, 1000)
          .Rvar$tmppk<-rjags::coda.samples(.Rvar$pkjags, variable.names=c('p','k'), n.iter=2000)[[1]][,2:1]
          },
         file = paste0(.Rvar$datadir, '/NULL')
        )
        pkstat <- list(phat = mean(.Rvar$tmppk[,1]), CIp = quantile(.Rvar$tmppk[,1],c(0.025, 0.975)), khat = mean(.Rvar$tmppk[,2]), CIk = quantile(.Rvar$tmppk[,2],c(0.025, 0.975)), r = cor(.Rvar$tmppk[,1:2])[2])
        .Rvar$pkres<-array(dim=c(dim(.Rvar$tmppk)[1],pkdat$n+2))
        pkstat<-pkstat
        .Rvar$pkres[,1:2]<-.Rvar$tmppk
        .Rvar$pkres[,3]<-.Rvar$tmppk[,1]
        for (i in 4:(pkdat$n+2)){
          .Rvar$pkres[,i]<-.Rvar$pkres[,i-1]*.Rvar$pkres[,2]
        }
        npk<-dim(.Rvar$pkres)[2]-2
        par(mar=c(4, 4, 2.5, 1.5),mgp=c(2,.7,0))
        if (length(dev.list()) > 0 && !prod(par()$fig==c(0,1,0,1))) {
          wd = 7*.Rvar$charSizeAdjust; ht = 7*.Rvar$charSizeAdjust
          if (.Rvar$platform == 'windows') windows.options(width = wd, height = ht)
          if (.Rvar$platform == 'mac') quartz.options(width = wd, height = ht)
          if (.Rvar$platform == 'linux') X11.options(width = wd, height = ht)
          dev.new(noRStudioGD =T)
        }
        plot(1:npk,rep(.5,npk), type='n',
          xlim = c(.75, npk + .25),
          ylim=range(apply(.Rvar$pkres[,-(1:2)], F = quantile, M = 2, probs=c(.975, 0.025))),
          xlab = "Search Occasion (i)",
          ylab = expression(Searcher ~Efficiency~ (p[~i]))
        )
        if(.Rvar$platform == "windows") bringToTop()
      #  for (i in 1:npk){
      #    lines(rep(i,2), qbeta(c(0.05, 0.95), junk$X[i]+.5, junk$M[i]-junk$X[i]+.5), col=colors()[125])
      #  }
        lines(1:npk, apply(.Rvar$pkres[,-(1:2)], F = mean, M = 2),lwd=2)
        lines(1:npk, apply(.Rvar$pkres[,-(1:2)], F = quantile, M = 2, probs=.975),lty=3, col=2)
        lines(1:npk, apply(.Rvar$pkres[,-(1:2)], F = quantile, M = 2, probs=.025),lty=3, col=2)
        points(1:npk, junk$X/junk$M, pch = 18, cex = 3*.Rvar$charSizeAdjust*sqrt(junk$M)/sqrt(max(junk$M)), col = colors()[125])
        pklabs<-list()
        for (i in 1:npk) pklabs[[i]]<-paste0(junk$X[i],"/",junk$M[i])
        text(1:npk, junk$X/junk$M, lab = pklabs, col=colors()[125], adj = -.25*c(1,1), cex = .75*.Rvar$charSizeAdjust)
        legend(x='topright', legend = c('Empirical X[i]/M[i]', 'Fitted values', '95% credibility bands'),
          lty=c(NA, 1, 3), lwd=c(NA, 2,1), pch = c(18, NA, NA), pt.cex = 1.5, col = c(colors()[125], 1, 2)
        )
        lhs<-par('usr')[1]+diff(par('usr')[1:2])*0.015
        mtext(text = paste0(
          "95% CIs: p \u2208 [", round(quantile(.Rvar$pkres[,1],0.025),3),", ", round(quantile(.Rvar$pkres[,1],0.975),3),"], ",
          "k \u2208 [", round(quantile(.Rvar$pkres[,2],0.025),3), ", ", round(quantile(.Rvar$pkres[,2],0.975),3),"]"),
          side = 1, line = -2.5, at=lhs, adj=0, cex = 0.8*.Rvar$charSizeAdjust)
#        mtext(text = bquote(hat(r)["p,k"]== .(round(cor(.Rvar$pkres[,1],.Rvar$pkres[,2]),4))), side = 1, line = -1.25, at=lhs, adj=0, cex = 0.8*.Rvar$charSizeAdjust)
        mtext(text = paste0("p\u0302 = ",round(.Rvar$pkdat$X[1]/.Rvar$pkdat$M[1],3), ", k\u0302 = ", round(mean(.Rvar$pkres[,2]),3)), side = 1, line = -1.25, at=lhs, adj=0, cex = 0.8*.Rvar$charSizeAdjust)
        title(expression(Estimation~of~p~and~k:~~p[i]==p*k^"i - 1"))
      })
      pkSave<-tkbutton(pkButFrame, text = "OK", width = 6, command = function(){
        # error check
        if (!pkTableChk()){
          tkmessageBox(message="Error in data. Cannot calculate...")
          return(F)
        }
        # fix dimensions of pkdat array
        junk<-list()
        # add new data to pkdat
        junk$n<-as.numeric(tclvalue(tcl(pkTable,"index","end","row")))
        junk$X<-numeric(junk$n)
        junk$M<-numeric(junk$n)
        for (i in 1:junk$n){
          junk$M[i]<-toR(pkData[[i,1]])
          junk$X[i]<-toR(pkData[[i,2]])
        }

        pkdat<-junk
        capture.output({
          .Rvar$pkjags <- rjags::jags.model(
            textConnection(pkmod),
            data = pkdat,
            inits = with(pkdat,list(p = X[1]/M[1],k = max(min((X[2]/M[2])/(X[1]/M[1]),.99),.01)))
          )
          update(.Rvar$pkjags, 1000)
          .Rvar$tmppk<-rjags::coda.samples(.Rvar$pkjags, variable.names=c('p','k'), n.iter=2000)[[1]][,2:1]
          },
          file = paste0(.Rvar$datadir, '/NULL')
        )
        pkstat <- list(phat = mean(.Rvar$tmppk[,1]), CIp = quantile(.Rvar$tmppk[,1],c(0.025, 0.975)), khat = mean(.Rvar$tmppk[,2]), CIk = quantile(.Rvar$tmppk[,2],c(0.025, 0.975)), r = cor(.Rvar$tmppk[,1:2])[2])
        .Rvar$pkres<-array(dim=c(dim(.Rvar$tmppk)[1],pkdat$n+2))
        pkstat<-pkstat
        .Rvar$pkres[,1:2]<-.Rvar$tmppk
        .Rvar$pkres[,3]<-.Rvar$tmppk[,1]
        for (i in 4:(pkdat$n+2)){
          .Rvar$pkres[,i]<-.Rvar$pkres[,i-1]*.Rvar$pkres[,2]
        }
        .Rvar$pkdat<-junk
        pkstat<-list(phat = mean(.Rvar$pkres[,1]),
          CIp = quantile(.Rvar$pkres[,1],c(0.025, 0.975)),
          khat = mean(.Rvar$pkres[,2]), CIk = quantile(.Rvar$pkres[,2],c(0.025, 0.975)),
          r = cor(.Rvar$pkres[,1], .Rvar$pkres[,2])
        )
        tclvalue(dynLbl$pkci)<<-paste0("95% CIs: p \u2208 [", round(pkstat$CIp[1],3), ", ", round(pkstat$CIp[2],3),"], ",
          "k \u2208 [", round(pkstat$CIk[1],3), ", ", round(pkstat$CIk[2],3),"]")
#        tclvalue(dynLbl$pkr) <- paste0("correlation r\u0302(p,k) = ",round(pkstat$r,3))
        tclvalue(dynLbl$pkr)<<-paste0("p\u0302 = ",round(.Rvar$pkdat$X[1]/.Rvar$pkdat$M[1],3), ", k\u0302 = ", round(pkstat$khat,3))
        .Rvar$pkstat <- pkstat
        tkdestroy(pkModule)
        return(T)
      })
      pkCancel<-tkbutton(pkButFrame, text = "Cancel", width=6, command = function(){
      # set radio in symodule back to 'estimate p only' mode
        tkdestroy(pkModule)
      })

      tkgrid(pkView)
      tkgrid(pkSave)
      tkgrid(pkCancel)
      tkgrid(pkTable)
      tkgrid(pkButFrame, row = 0, column = 1, padx = 15, pady = 15, sticky = 'nw')
      tkconfigure(pkModule, width=500)
    },
    SYreadparm = function(){
      filename <- tclvalue(tkgetOpenFile(filetypes = "{{R images} {.Rdata}}",defaultextension = ".Rdata", initialfile = '.Rdata', title = "Read", initialdir = .Rvar$csvpath))
      tmp<-unlist(strsplit(filename,'/')); pathname<-paste(tmp[-length(tmp)],collapse='/')
      if (filename == "") return(FALSE)
      tryCatch(load(filename, env = .Rvar), error=function(){tkmessageBox(icon='error',message=paste0("Error: Unable to read file:\n\n",filename)); return(F)})
      assign("csvpath", pathname, env = .Rvar)
      # check whether search schedule type, persistence distribution, prior, and arrival function are defined properly
      if (SYparmLoadable(.Rvar$syProvisional)){ # if loadable, then redraw form
        .Rvar$CPdataPrevious<-.Rvar$CPdata
        .Rvar$pkdatPrevious<-.Rvar$pkdat
        tkdestroy(syModule)
        if (partial) {#then read the arrival and prior info from symc rather than from loaded file
          for (nm in names(singleYearDefault)){
            if (nm %in% names(symcDefault)){
              .Rvar$syProvisional[[nm]]<-.Rvar$symcdat[[nm]]
            }
          }
        }
        initialize(.Rvar$syProvisional, ifelse(name == 'symc', T, F)) # loads the values to the form
        if (partial) tkwm.title(syModule,paste0(tmp[length(tmp)], " - EoA, v", .Rvar$VER, " - Search Class")) else tkwm.title(syModule,paste0(tmp[length(tmp)], " - EoA, v", .Rvar$VER, " - Single Class Module"))
#        syChkAll()
#        chkForm.sy()
      }
    },
    Xchk = function(X) {
      X<-suppressWarnings(as.numeric(X))
      if (is.na(X) || length(X) != 1 || nchar(X) == 0){
        Xok<<-F
        return(F)
      }
      if (is.na(X) || X < 0 || round(X)!=X) {
        Xok<<-F
        return(F)
      } else {
        Xok<<-T
      }
      return(T)
    },
    achk = function(acoverage){
      acoverage<-suppressWarnings(as.numeric(acoverage))
      if (is.na(acoverage) || length(acoverage) != 1 || nchar(acoverage)==0){
        aok<<-F
        return(F)
      } else if ( is.na(acoverage) || acoverage<=0 || acoverage>1){
        aok<<-F
        return(F)
      } else {
        aok<<-T
      }
      return(T)
    },
    vchk = function(tcoverage){
      tcoverage<-suppressWarnings(as.numeric(tcoverage))
      if (is.na(tcoverage) || length(tcoverage) != 1 || nchar(tcoverage)==0){
        vok<<-F
        return(F)
      } else if ( is.na(tcoverage) || tcoverage<=0 || tcoverage>1){
        vok<<-F
        return(F)
      } else {
        vok<<-T
      }
      return(T)
    },
    SEnchk = function(SEn){
      SEn<-suppressWarnings(as.numeric(SEn))
      if (is.na(SEn) || length(SEn) != 1 ||nchar(SEn)==0){
        SEnok<<-F
        return(F)
      } else if (is.na(SEn) || SEn <= 0 || round(SEn) != SEn) {
        SEnok<<-F
        return(F)
      } else {
        SEnok<<-T
        if (SExleg) { #SEx is a legitimate value (i.e. positive integer)
          SEx<-as.numeric(tclvalue(tkvars$SEx))
          if (SEx <= SEn){
            SExok<<-T
          } else {
            SExok<<-F
          }
        }
      }
      return(T)
    },
    SExchk = function(SEx) {
      SEx<-suppressWarnings(as.numeric(SEx))
      if (length(SEx) != 1 || nchar(SEx)==0 || is.na(SEx)){
        SExok<<-F
        SExleg<<-F
        return(F)
      } else if (SEx<=0 || round(SEx)!=SEx) {
        SExok<<-F
        SExleg<<-F
        return(F)
      } else {
        SExleg<<-T
        if (SEnok) {
          SEn<-as.numeric(tclvalue(tkvars$SEn))
          if (SEx<=SEn){
            SExok<<-T
          } else {
            SExok<<-F
          }
        }
      }
      if (SExok)  return(T) else return(F)
    },
    kchk = function(kparm) {          #
      kparm<-suppressWarnings(as.numeric(kparm))
      if (is.na(kparm) || length(kparm)!=1 || nchar(kparm)==0){
        kok<<-F
        return(F)
      } else if (is.na(kparm) || kparm < 0 || kparm > 1) {
        kok<<-F
        return(F)
      } else {
        kok<<-T
      }
      return(T)
    },
    Ichk = function(Iparm){ # no real restrictions...needs to be positive
      Iparm<-suppressWarnings(as.numeric(Iparm))
      if (is.na(Iparm) || length(Iparm) != 1 || nchar(Iparm)==0){
        Iok<<-F
        return(F)
      } else if (is.na(Iparm) || Iparm<=0) {
        rhat<-NA
        Iok<<-F
        return(F)
      } else {
        Iok<<-T
        Isam<<-Iparm
        Ir<<-Iparm
    #    tclvalue(tkvars$Ir)<<-Ir
      }
      return(T)
    },
    nsearchchk = function(nsearchparm){
      nsearchparm<-suppressWarnings(as.numeric(nsearchparm))
      if (is.na(nsearchparm) || length(nsearchparm) != 1 || nchar(nsearchparm)==0){
        nsearchok<<-F
        return(F)
      } else if (is.na(nsearchparm) || nsearchparm <= 0 || nsearchparm!=round(nsearchparm)) {
        nsearchok<<-F
        nsearchleg<<-F
        return(F)
      } else {
        nsearchok<<-T
      }
      return(T)
    },
    startchk = function(startparm){
      if (is.na(startparm) || length(startparm)!=1 || nchar(startparm)==0 || class(try(as.Date(startparm),silent=T))=="try-error"){ # text box is empty OR not a date
        startok<<-F
        return(F)
      } else {
        startok<<-T
    #    singleYearProvisional$firstsearch<<-startparm
      }
      return(T)
    },
    pdachk = function(pda) { # need to update labels for r/CP if necessary
      pda<-suppressWarnings(as.numeric(pda))
      if (is.na(pda) || length(pda) != 1 || nchar(pda)==0){ # all the pda's need to be greater than zero
        pdaok<<-F
        return(F)
      } else if (pda<=0) {
        pdaok<<-F
        return(F)
      } else {
        pdaok<<-T
      }
      return(T)
    },
    pdbchk = function(pdb){
    # pdb's need to be greater than zero except for lognormal's may be negative
      pdb<-suppressWarnings(as.numeric(pdb))
      if (pdaok) pda<-as.numeric(tclvalue(tkvars$pda))
      if (is.na(pdb) || (pdb<=0 & persistence_distn != "Lognormal")) {
        pdbok<<-F
        return(F)
      } else {
        pdbok<<-T
        if (persistence_distn=="Exponential"){
          pdaok<<-T
          pda<-1/pdb
          tclvalue(tkvars$pda)<-signif(pda,3)
        }
      }
      return(T)
    },
    blwrchk = function(blwr){
      blwr<-suppressWarnings(as.numeric(blwr))
      if (is.na(blwr) || length(blwr)!=1 || nchar(blwr)==0){
        bminok<<-F
        bminleg<<-F
        return(F)
      } else if (blwr<=0 & persistence_distn!="Lognormal" ) {
        bminok<<-F
        bminleg<<-F
        return(F)
      } else {
        bminleg<<-T
        if (pdbok) {
          pdb<-as.numeric(tclvalue(tkvars$pdb))
          if (blwr < pdb){
            bminok<<-T
          } else {
            bminok<<-F
          }
        }
      }
      return(T)
    },
    buprchk = function(bupr) {
      bupr<-suppressWarnings(as.numeric(bupr))
      if (is.na(bupr) || length(bupr)!=1 || nchar(bupr)==0){
        bmaxok<<-F
        bmaxleg<<-F
        return(F)
      } else if (bupr<=0 & persistence_distn!="Lognormal") {
        bmaxok<<-F
        bmaxleg<<-F
        return(F)
      } else {
        bmaxleg<<-T
        if (pdbok) {
          pdb<-as.numeric(tclvalue(tkvars$pdb))
          if (bupr > pdb){
            bmaxok<<-T
          } else {
            bmaxok<<-F
          }
        }
      }
      return(T)
    },
    syAchk = function(syA) {          #
      syA<-suppressWarnings(as.numeric(syA))
      if (length(syA) != 1){
        syAok<<-F
        return(F)
      }
      if (is.na(syA) || nchar(syA)==0 || syA<=0 || syA>=1) {
        syAok<<-F
        return(F)
      } else {
        syAok<<-T
      }
      return(T)
    }
  ),
  active = list(
    allok1 = function(){
      Xok<<-T
      aok<<-T
      SEnok<<-T;  SExok<<-T;  SExleg<<-T
      kok<<-T
      Iok<<-T;  nsearchok<<-T; startok<<-T
      lamok<<-T
      pdaok<<-T
      pdbok<<-T; bminok<<-T; bmaxok<<-T; bminleg<<-T; bmaxleg<<-T
      Irok<<-T
      syAok<<-T
      gook<<-T
      prok<<-T
    },
    updateprlab = function(){
      tkconfigure(prior.lbl,text=paste0("95th percentile = ",min(which(cumsum(toR(tkvars$custom.prior))>=0.95))-1))
      tkconfigure(priorViewButton,state="normal")
    }
  )
)

sysave<-function(sy, CPdataPrevious, pkdatPrevious, csvpath, sysc = F){
  if (sysc) {
    filename<-paste0(.Rvar$datadir, "/syscPrevious.Rdata")
    syscPrevious<-sy
    save(syscPrevious, CPdataPrevious, pkdatPrevious, file = filename)
  } else {
    filename<-paste0(.Rvar$datadir,"/singleYearPrevious.Rdata")
    singleYearPrevious<-sy
    save(singleYearPrevious, CPdataPrevious, pkdatPrevious, file = filename)
  }
  save(csvpath, file = paste0(.Rvar$datadir,"/csvpath.Rdata"))
}