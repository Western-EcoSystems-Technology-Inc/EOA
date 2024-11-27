scexForm <- R6::R6Class("scexForm",
  portable = FALSE,
  public = list(
    # tcl variables:
    scexModule = NA, tkvars = list(),

    # R variables:
    wid = 5, # width of (many of) the data entry boxes
    iAMAschedule = NA, iAMAtable = NA, iAMAdata = NA,
    scextmp = list(),
    # data checks...
    tableok = T, nyrok = T, Tauok = T, nsimok = T, lambdaok = T, rhoinfok = T,
    yrsok = T, g1ok = T, g1lwrok = T, g1uprok = T, g2ok = T, g2lwrok = T, g2uprok = T,
    rtaok = T, rhorevok = T, ltaok = T, styok = T, staok = T,
    viewone = 0, doneonce = F,

    # initial values
    initialize = function(scexdat) {
      # create the required tcl variables from the data input
#      allscexok() # assuming a valid data set (should be true if loading scexDefault or scexPrevious)
                  # error-checking happens later
      for (i in 1:length(scexVar)){
        ind <- which(names(scexdat)==scexVar[i])
        tkvars[[scexVar[i]]] <<- tclVar(as.tclObj(ifelse(length(ind) > 0, scexdat[[ind]],""),drop=T))
      }
      narrvar<-length(scexArray)
      if (narrvar > 0){
        for (li in 1:narrvar){
          ind <- which(names(scexdat)==scexArray[li])
          tmp<-tclArray()
          if (length(ind) == 0){ # missing data
            tmp[[0]] <- as.tclObj('', drop=T)
          } else if (length(scexdat[[ind]]) == 1) { # scalar
            tmp[[0]] <- as.tclObj(ifelse(!is.na(scexdat[[ind]]), scexdat[[ind]], ''), drop=T)
          } else if (is.null(dim(scexdat[[ind]]))) { # vector
            for (j in 1:length(scexdat[[ind]])) tmp[[j-1]] <- as.tclObj(scexdat[[ind]][j],drop=T)
          } else { # matrix
            for (rowi in 1:dim(scexdat[[ind]])[1])
              for (coli in 1:dim(scexdat[[ind]])[2])
                tmp[[rowi-1, coli-1]] <- as.tclObj(scexdat[[ind]][rowi, coli], drop = T)
          }
          tkvars[[scexArray[li]]] <<- tmp
        }
      }
      iAMAschedule<<-scexdat$iAMAschedule

      scexok<-function(){
      # must be true for all options:
        if(! (nyrok & Tauok & lambdaok & rhoinfok & yrsok)){
          tkmessageBox(message="Error in data. Cannot calculate.")
          return(F)
        }
        if ((toR(tkvars$ltT) & !ltaok) || (toR(tkvars$stT) & (!staok | !styok)) || (toR(tkvars$rT) & (!rtaok | !rhorevok))){
          tkmessageBox(message="Error in data. Cannot calculate.")
          return(F)
        }
        if (toR(tkvars$stT) & toR(tkvars$iAMA) & !iAMAtablechk()){
          tkmessageBox(message="Error in data. Cannot calculate.")
          return(F)
        }
        if (toR(tkvars$yrs) > toR(tkvars$nyr)){
          tkconfigure(yrs.edit,bg=colors()[400])
          tkconfigure(nyr.edit,bg=colors()[400])
          tkmessageBox(message="Years of intensive monitoring cannot exceed years in permit\n\nCannot calculate.")
          return(F)
        }
        if (toR(tkvars$yrs) > 0){
          if (g1ok & g1lwrok & g1uprok){
            if (toR(tkvars$g1)<=toR(tkvars$g1lwr)){
              tkconfigure(g1.edit,bg=colors()[400])
              tkconfigure(g1lwr.edit,bg=colors()[400])
              tkmessageBox(message="mean g must be within 95% CI for g")
              return(F)
            }
            if (toR(tkvars$g1)>=toR(tkvars$g1upr)){
              tkconfigure(g1.edit,bg=colors()[400])
              tkconfigure(g1upr.edit,bg=colors()[400])
              tkmessageBox(message="mean g must be within 95% CI for g")
              return(F)
            }
          } else {
            tkmessageBox(message="Error in data. Cannot calculate.")
            return(F)
          }
        }
        if (toR(tkvars$yrs) < toR(tkvars$nyr)){
          if (g2ok & g2lwrok & g2uprok){
            if (toR(tkvars$g2)<=toR(tkvars$g2lwr)){
              tkconfigure(g2.edit,bg=colors()[400])
              tkconfigure(g2lwr.edit,bg=colors()[400])
              tkmessageBox(message="mean g must be within 95% CI for g")
              return(F)
            }
            if (toR(tkvars$g2)>=toR(tkvars$g2upr)){
              tkconfigure(g2.edit,bg=colors()[400])
              tkconfigure(g2upr.edit,bg=colors()[400])
              tkmessageBox(message="mean g must be within 95% CI for g")
              return(F)
            }
          } else {
            tkmessageBox(message="Error in data. Cannot calculate.")
            return(F)
          }
        }
        return(T)
      }
      tcltk::tcl("option","add","*tearOff",0)
      scexModule<<-tktoplevel(); tkgrab.set(scexModule);  tkfocus(scexModule)
      tkwm.title(scexModule,paste0("EoA, v", .Rvar$VER, " - Scenario Explorer"))
      tkwm.resizable(scexModule,0,0)
      SCEXtopMenu <- tkmenu(scexModule); tkconfigure(scexModule,menu=SCEXtopMenu)
      SCEXhelpMenu <- tkmenu(SCEXtopMenu,activebackground=colors()[125],activeforeground=colors()[109])
      SCEXeditMenu <- tkmenu(SCEXtopMenu,activebackground=colors()[125],activeforeground=colors()[109])
      tkadd(SCEXeditMenu,"command",label="Restore defaults",command=function(){
        .Rvar$scexPrevious <- scexDefault
        suppressWarnings(file.remove(paste0(.Rvar$datadir,"/scexPrevious.Rdata")))
        tkwm.state(scexModule, 'withdraw')
        initialize(.Rvar$scexPrevious)
      })
      tkadd(SCEXeditMenu,"command",label="Restore previous",command=function(){
        tkdestroy(scexModule)
        initialize(.Rvar$scexPrevious)
      })
      tkadd(SCEXeditMenu,"command",label="Save to file (.rds)",command=function() {
        junk<-list()
        for (nm in names(tkvars)) junk[[nm]]<-toR(tkvars[[nm]])
        saveparm(junk)
        tkwm.title(scexModule,paste0(.Rvar$dataFileTitle, " - EoA, v", .Rvar$VER, " - Scenario Explorer"))
      })
      tkadd(SCEXeditMenu,"command",label="Read from file (.rds)",command=function() {
        filename <- tclvalue(tkgetOpenFile(filetypes = "{{R data files} {.rds}}",defaultextension = ".rds", initialfile = '.rds', title = "Read"))
        tmp<-unlist(strsplit(filename,'/')); pathname<-paste(tmp[-length(tmp)],collapse='/')
        if (filename == "") return(FALSE)
        parms <- tryCatch(readRDS(filename), error=function(){tkmessageBox(icon='error',message=paste0("Error: Unable to read file:\n\n",filename)); return(F)})
        .Rvar$csvpath<-pathname
        tkdestroy(scexModule)
        initialize(parms)
        tkwm.title(scexModule,paste0(tmp[length(tmp)], " - EoA, v", .Rvar$VER, " - Scenario Explorer"))
      })
      tkadd(SCEXtopMenu, "cascade", label="Edit", menu=SCEXeditMenu)
      tkadd(SCEXhelpMenu, "command", label="Parameter Conversion", command = conversionsCalculator)
      tkadd(SCEXhelpMenu, "command", label="About", command=function() tkmessageBox(title='Evidence of Absence (EoA)',message=about_text))
      tkadd(SCEXtopMenu, "cascade", label="Help", menu=SCEXhelpMenu)
      ###
      scexGenFrame<-ttklabelframe(scexModule,text="General Framework",padding=10)
      scexFieldFrame<-ttklabelframe(scexModule,text="Field Parameters")
        scexg1Frame<-ttklabelframe(scexFieldFrame,text="Detection probability (g), intensive monitoring")
        scexg2Frame<-ttklabelframe(scexFieldFrame, text="Detection probability (g), non-intensive monitoring")
        scexg3Frame<-ttklabelframe(scexFieldFrame, text="Post-trigger Monitoring")
      scexGovFrame<-ttklabelframe(scexModule,text="Governing Parameters")
        scexLTFrame<-ttklabelframe(scexGovFrame,text="Long-Term Trigger")
        scexSTFrame<-ttklabelframe(scexGovFrame,text="Short-Term Trigger")
          scexAMAiFrame<-ttklabelframe(scexSTFrame,text="Schedule of incremental AMAs")
        scexRTFrame<-ttklabelframe(scexGovFrame,text="Reversion Trigger")

      # labels and entry widgets for variables
      wid<-5
      ltT.box<-tkcheckbutton(scexLTFrame, variable = tkvars$ltT) # include long-term trigger in analysis?
      stT.box<-tkcheckbutton(scexSTFrame,text="", variable = tkvars$stT) # include short-term trigger in analysis?
      rT.box<-tkcheckbutton(scexRTFrame,text="", variable = tkvars$rT) # include reversion trigger in analysis?
      nyr.lbl<-tklabel(scexGenFrame, text="Total years in permit (n)")
      nyr.edit<-tkentry(scexGenFrame, textvariable = tkvars$nyr, justify = 'right', width=wid, bg = 'white')
      Tau.lbl<-tklabel(scexGenFrame, text="Total permitted take (\u03a4)")
      Tau.edit<-tkentry(scexGenFrame, textvariable = tkvars$Tau, justify = 'right', width=wid, bg = 'white')
      nsim.lbl<-tklabel(scexGenFrame, text="Number of simulation draws (nsim)")
      nsim.edit<-tkentry(scexGenFrame, textvariable = tkvars$nsim, justify='right', width=wid, bg = 'white')
      lambda.lbl<-tklabel(scexFieldFrame, text="Baseline fatality rate (\u03bb)")
      lambda.edit<-tkentry(scexFieldFrame, textvariable = tkvars$lambda, justify='right', width=wid, bg = 'white')
      rhoinf.lbl<-tklabel(scexFieldFrame, text="Effect of avoidance AMA (\u03c1\u221E)")
      rhoinf.edit<-tkentry(scexFieldFrame, textvariable = tkvars$rhoinf, justify='right', width=wid, bg = 'white')
      yrs.lbl<-tklabel(scexFieldFrame, text="Initial years of intensive monitoring")
      yrs.edit<-tkentry(scexFieldFrame, textvariable = tkvars$yrs, justify='right', width=wid, bg = 'white')
      g1.lbl<-tklabel(scexg1Frame, text="mean g")
      g1.edit<-tkentry(scexg1Frame, textvariable = tkvars$g1, justify='right', width=wid, bg = 'white')
      g1ci.lbl<-tklabel(scexg1Frame, text = "95% CI")
      g1lwr.edit<-tkentry(scexg1Frame, textvariable = tkvars$g1lwr, justify='right', width=wid, bg = 'white')
      g1upr.edit<-tkentry(scexg1Frame, textvariable = tkvars$g1upr, justify='right', width=wid, bg = 'white')
      g2.lbl<-tklabel(scexg2Frame, text="mean g")
      g2.edit<-tkentry(scexg2Frame, textvariable = tkvars$g2, justify='right', width=wid, bg = 'white')
      g2ci.lbl<-tklabel(scexg2Frame, text = "95% CI")
      g2lwr.edit<-tkentry(scexg2Frame, textvariable = tkvars$g2lwr, justify='right', width=wid, bg = 'white')
      g2upr.edit<-tkentry(scexg2Frame, textvariable = tkvars$g2upr, justify='right', width=wid, bg = 'white')
#      tkvars$gpost<-tclVar(.Rvar$scexPrevious$gpost)
#      tkvars$g3i<-tclVar(.Rvar$scexPrevious$g3i)
      gpost.box<-tkcheckbutton(scexg3Frame, variable = tkvars$gpost, text = "Monitoring continues\n after long-term trigger")
      minigframe<-tkframe(scexg3Frame)
      g3i.intenseRadio<-tkradiobutton(minigframe, variable = tkvars$g3i, text = "Intensive", value = 1)
      g3i.nonintenseRadio<-tkradiobutton(minigframe, variable = tkvars$g3i, text = "Non-intensive", value = 2)
      tkgrid(g3i.intenseRadio, sticky='w')
      tkgrid(g3i.nonintenseRadio, sticky='w')
      rta.lbl<-tklabel(scexRTFrame, text="\u03b1")
      rta.edit<-tkentry(scexRTFrame, textvariable = tkvars$rta, justify='right', width=wid, bg = 'white')
      rhorev.lbl<-tklabel(scexRTFrame, text="\u03c1\u2080")
      rhorev.edit<-tkentry(scexRTFrame, textvariable = tkvars$rhorev, justify='right', width=wid, bg = 'white')
      lta.lbl<-tklabel(scexLTFrame, text="\u03b1")
      lta.edit<-tkentry(scexLTFrame, textvariable = tkvars$lta, justify='right', width=wid, bg = 'white')
      sty.lbl<-tklabel(scexSTFrame, text="Term (y)")
      sty.edit<-tkentry(scexSTFrame, textvariable = tkvars$sty, justify='right', width=wid, bg = 'white')
      sta.lbl<-tklabel(scexSTFrame, text="\u03b1")
      sta.edit<-tkentry(scexSTFrame, textvariable = tkvars$sta, justify='right', width=wid, bg = 'white')
      stTonyesRadio<-tkradiobutton(scexSTFrame, variable = tkvars$stTon, value = 'same', text = "\u03c4 = \u03a4/n", command=function() tkgrid.remove(tau.edit))
      stTonnoRadio<-tkradiobutton(scexSTFrame, variable = tkvars$stTon, value = 'different', text = "\u03c4 \u2260 \u03a4/n", command=function() tkgrid(tau.edit, row = 1, column = 4))
      tau.edit<-tkentry(scexSTFrame, textvariable = tkvars$tau, justify = 'right', width=wid, bg='white')
      iAMA.box<-tkcheckbutton(scexAMAiFrame, text="Include incremental AMA", variable = tkvars$iAMA)
      buttFrame<-tkframe(scexModule)
      bw<-15
      calcButton<-tkbutton(buttFrame, text="Calculate", width=bw, command=function(){
        if (scexok()){
          v<-suppressWarnings(as.numeric(tclvalue(tkvars$nsim)))
          if (!val.numeric(v) || !val.integer(v) || !val.gt0(v)){
            tkmessageBox(message="Number of simulation draws must be a positive integer. Cannot calculate.")
            return(F)
          } else if (v > 10000 && (toR(tkvars$rT) | (toR(tkvars$stT) & toR(tkvars$iAMA)))){
            tkconfigure(nsim.edit,bg=colors()[542])
            if (!tclvalue(tkmessageBox(message=paste0("Simulation with ", v, " draws will run until the cows come home.\n\nAre you sure you want to go through with this?"),        type="yesno"))=="yes"){
              return(F)
            }
          }
        } else {
          return(F)
        }
        .Rvar$scexPrevious<-feedR.scex(self$tkvars)
        scexCalc(.Rvar$scexPrevious)
      })
      viewOneButton<-tkbutton(buttFrame, text="View one example", width=bw, command=function(){
        if (!scexok()) return(F)
        tkwm.state(scexModule, "withdraw")
        .Rvar$scexPrevious<-feedR.scex(self$tkvars)
        .Rvar$scex1mod<-scex1Form$new(.Rvar$scexPrevious)
        scexCalc(.Rvar$scexPrevious, viewone = 1, .Rvar$scex1mod)
      })
      closeButton<-tkbutton(buttFrame,text="Close", width=bw, command=function(){
        # destroy the window
        graphics.off()
        tkdestroy(scexModule)
        scexsave(.Rvar$scexPrevious)
        tkwm.deiconify(.Rvar$EoA)
      })
      tkgrid(calcButton, viewOneButton, closeButton)
      #### step 1: create table for holding, displaying, editing data
      steps<-dim(iAMAschedule)[1]
      iAMAdata<<-tclArray()
      columnNames<-c("","\n\u03c1", "Years of intensive monitoring")
      for (i in 1:length(columnNames)) iAMAdata[[0,i-1]]<-as.tclObj(strsplit(columnNames[i]," ",fixed=T)[[1]],drop=T)
      for (i in 1:steps)
      for (i in 1:steps){
        iAMAdata[[i,0]]<-as.tclObj(paste("Step", i),drop=T)
        iAMAdata[[i,1]]<-iAMAschedule[i,1]
        iAMAdata[[i,2]]<-iAMAschedule[i,2]
      }
      valChar<-function(S){
        actcol<-as.numeric(tclvalue(tkindex(iAMAtable,"active","col")))
        actrow<-as.numeric(tclvalue(tkindex(iAMAtable,"active","row")))
        if (length(grep("\n",S))>0){ #
          return(tcl("expr", FALSE))
        } else {
          iAMACellChk(S, actrow, actcol)
          return(tcl("expr", TRUE))# but other kinds of space are not--> error-checking if good value is added
        }
      }

      iAMAtable <<- tcltk2::tk2table(scexAMAiFrame,
        rows=dim(iAMAschedule)[1]+1,
        cols=3,
        selectmode="extended",
        variable=iAMAdata,
        titlerows=1,
        titlecols=1,
        background='white',
        resizeborders="none",
        multiline=T,
        rowseparator='\n',
        colseparator='\t',
        wrap = 1,
        validate= 1,
        justify="center",
        vcmd=function(S) valChar(S),
        selecttitle = 1
      )
      tcl(iAMAtable, "tag", "configure", "readonly", state='disabled', bg = colors()[415])
      tcl(iAMAtable, "tag", "configure", "error",    bg=colors()[652]) # cells with error have yellow background
      tcl(iAMAtable, "tag", "configure", "cellok",   bg='white')
      tcl(iAMAtable, "tag", "configure", "active",   fg='black',relief='groove')
      tcl(iAMAtable, "tag", "configure", "ignore",   state='disabled',bg=colors()[365],fg=colors()[42])
      tcl(iAMAtable, "tag", "celltag", "readonly", cellind(steps,1))
      #tkbind(iAMAtable,"<Button-1>", function(){
      #  if (tclvalue(tcl(iAMAtable,"cget","-state"))=="disabled")  tkconfigure(iAMAtable,selectbackground='yellow')
      #})
      colwidths<-c(6,6,12)
      for(i in 1:3) {
        tcl(iAMAtable, "width", i - 1, colwidths[i])
      }
      tcl(iAMAtable,"height",0, 3)
      tkbind(iAMAtable,"<Return>",function(){
        if (tclvalue(tkindex(iAMAtable,"active","row"))==tclvalue(tkindex(iAMAtable,"end","row"))){
          tkevent.generate(iAMAtable,"<Control-KeyPress-a>")
          return(T)
        }
        if(tclvalue(tkindex(iAMAtable,"active","row"))    == tclvalue(tkindex(iAMAtable,"end","row"))){
          if (tclvalue(tkindex(iAMAtable,"active","col")) != tclvalue(tkindex(iAMAtable,"end","col"))){
            tkevent.generate(iAMAtable,"<KeyPress-Right>")
          }
        } else {
          tkevent.generate(iAMAtable,"<KeyPress-Down>")
        }
      })
      tkbind(iAMAtable,"<Control-KeyPress-a>", function(){
      # add a row
        if (!toR(tkvars$iAMA)) return(F)
        tkinsert(iAMAtable, "rows", tclvalue(tkindex(iAMAtable, "active", "row")), -1)
        for (i in 1:2){
          tcl(iAMAtable,"tag","celltag","error", as.tclObj(paste0(tclvalue(tkindex(iAMAtable, "active", "row")),',', i), drop = T))
        }
        for (i in 1:toR(tkindex(iAMAtable, "end", "row"))) iAMAdata[[i,0]]<-as.tclObj(paste0("Step ",i), drop = T)
        rowind<-as.numeric(tclvalue(tcl(iAMAtable, "index", "end", "row")))
        tcl(iAMAtable,"activate", as.tclObj(paste0(rowind,',', 1)))
        iAMAdata[[rowind,7]]<- "NA"
        iAMAdata[[rowind,6]]<-"NA"
      })
      movechars<-c('Up','Down','Left','Right','Tab','Return','Shift_L','Alt_L')
      tkbind(iAMAtable,"<Key>",function(K){
        if (K %in% movechars) {
          if (substr(tclvalue(tcl(iAMAtable,"curvalue")),1,1)=='.') tcl(iAMAtable,"curvalue",paste0(0,tclvalue(tcl(iAMAtable,"curvalue"))))
        }
      })
      tkbind(iAMAtable,"<Tab>",function(){
        if(tclvalue(tkindex(iAMAtable,"active","col"))!=tclvalue(tkindex(iAMAtable,"end","col"))){
          tkevent.generate(iAMAtable,"<Right>")
        } else {
          tkevent.generate(iAMAtable,"<KeyPress-Down>")
        }
      })
      tkbind(iAMAtable,"<1>",function(){ # mouse button is pressed...check for errors
        if (substr(tclvalue(tcl(iAMAtable,"curvalue")),1,1)=='.') tcl(iAMAtable,"curvalue",paste0(0,tclvalue(tcl(iAMAtable,"curvalue"))))
      })
      tkbind(iAMAtable,"<Control-KeyPress-d>", function(){
      # delete current row (if there is only one row, just erase the data)
        actrow<-tclvalue(tkindex(iAMAtable, "active", "row"))
        if (!(as.numeric(actrow) %in% c(0,toR(tkindex(iAMAtable, "end", "row"))))){
          tkdelete(iAMAtable, "rows", actrow, 1)
          for (i in 1:toR(tkindex(iAMAtable, "end", "row"))) iAMAdata[[i,0]]<-as.tclObj(paste0("Step ", i), drop = T)
        }
      })
      tkevent.delete("<<Paste>>","<Control-v>") # hijack the normal activity of the ctrl + v paste command
      junk<-NA
      tkbind(iAMAtable,"<Control-Key-v>", function(){
        if (tclvalue(tkvars$ICEoption)=="e") return(F)
      # after pasting data from the clipboard, read the prior table data into a buffer
      # if there are errors in the new data, give error message and replace
        if (as.numeric(tclvalue(tcl(iAMAtable,"index","active","col")))>5)   return(F)
      #  junk<<-tkXselection.get(selection='CLIPBOARD') # grab dsta from the clipboard for preliminary analysis before pasting to the table
      #  junk<<-gsub('[^\t^\n]','',tclvalue(junk)) #remove everything but tabs and carriage returns for easy parsing of the data
        junk<-readClipboard()
        ncols<<-max(nchar(gsub('[^\t]','',junk)))
        if (ncols+as.numeric(tclvalue(tcl(iAMAtable,"index","active","col")))>5) {
          tkmessageBox(message="Error: cannot paste over summary columns",icon='error')
          return(F)
        }
        niAMA<-as.numeric(tclvalue(tkindex(iAMAtable,"end","row")))
        createTmpiAMA(niAMA)
        tkconfigure(iAMAtable,rows=max(niAMA+1,as.numeric(tclvalue(tcl(iAMAtable,"index","active","row")))+length(junk)),cols=8)
        tkevent.generate(iAMAtable,"<<Paste>>")
      })
      ########
      # labels and entry widgets for variables
      wid<-5
      tkgrid(nyr.lbl, sticky = 'w')
      tkgrid(nyr.edit, sticky='e', padx=8, row = 0, column = 1)
      tkgrid(Tau.lbl, sticky='w')
      tkgrid(Tau.edit,sticky='e',row=1, column=1,padx=8)
      tkgrid(nsim.lbl, nsim.edit)

      tkgrid(lambda.lbl, sticky='w', padx = c(10,0))
      tkgrid(lambda.edit, sticky='w', row=0, column = 1)
      tkgrid(rhoinf.lbl, sticky='w', padx = c(10,0))
      tkgrid(rhoinf.edit, sticky='w', row = 1, column = 1 )
      tkgrid(yrs.lbl, sticky='w', padx = c(10,0))
      tkgrid(yrs.edit, sticky='w', row=2, column =1)
      tkgrid(g1.lbl, g1.edit, sticky='w', padx=c(10,0))
      tkgrid(g1ci.lbl, sticky='w', padx=c(25,3), row=0, column = 2)
      tkgrid(g1lwr.edit, row=0, column = 3)
      tkgrid(g1upr.edit, row=0, column = 4)
      tkgrid(g2.lbl, g2.edit, sticky='w', padx=c(10,0))
      tkgrid(g2ci.lbl, sticky='w', padx=c(25,3), row=0, column =2)
      tkgrid(g2lwr.edit, row=0, column = 3)
      tkgrid(g2upr.edit, row=0, column = 4)
      tkgrid(gpost.box, padx = c(10,0))
      tkgrid(minigframe, padx = c(15,0), sticky='w', row = 0, column = 1)

      tkgrid(ltT.box, lta.lbl, lta.edit, pady=8)
      tkgrid(rT.box, pady=8)
      tkgrid(rta.lbl, padx=c(6,0), row=0, column = 1)
      tkgrid(rta.edit, row=0, column = 2)
      tkgrid(rhorev.lbl, padx=c(6,0), row=0, column = 3)
      tkgrid(rhorev.edit, padx=c(0,8), row=0, column = 4)
      tkgrid(stT.box, pady=8, padx=10, sticky='w')
      tkgrid(sta.lbl, row=0, column = 1, sticky='e')
      tkgrid(sta.edit, row=0, column = 2, sticky='w', padx=4)
      tkgrid(sty.lbl, row = 0, column = 3, sticky='w')
      tkgrid(sty.edit, row = 0, column = 4, sticky='e',padx=c(0,8))
      tkgrid(stTonyesRadio, row=1, column=0, sticky='w', padx=c(0,10), columnspan = 3)
      tkgrid(stTonnoRadio, row = 1, column = 2, sticky='e', columnspan=2)
      if (toR(tkvars$stTon)=='different') tkgrid(tau.edit, row = 1, column = 4)

      tkgrid(iAMA.box)
      tkgrid(iAMAtable,pady=8)
      #steps.lbl<-tklabel(scexLTFrame, text="Number of steps")
      #steps.edit<-tkentry(scexLTFrame, textvariable = tkvars$steps, justify='right', width=wid)

      tkgrid(scexAMAiFrame, columnspan=5, pady=10)
      tkgrid(scexLTFrame,sticky='w', pady=c(10,0))
      tkgrid(scexRTFrame,sticky='w', row=0, column =1, pady=c(10,0))
      tkgrid(scexSTFrame,sticky='w', columnspan=5, pady=15)
      tkgrid(scexg1Frame, columnspan = 2, padx=c(15,0), pady = c(10,0), sticky='w')
      tkgrid(scexg2Frame, columnspan = 2, padx=c(15,0), pady=c(10,0), sticky='w')
      tkgrid(scexg3Frame, columnspan = 2, padx=c(15,0), pady=10, sticky='w')
      tkgrid(scexGenFrame, sticky='nw',padx=10,pady=10)
      tkgrid(scexFieldFrame,padx = 10, sticky='n')
      tkgrid(scexGovFrame,columnspan=4, row=0, column=1, rowspan=4, sticky='nw', pady=10)
      tkgrid(buttFrame, row=2, column=0)
      if (toR(tkvars$ltT)) {
        tkconfigure(lta.lbl, state='normal')
        tkconfigure(lta.edit, state='normal')
      } else {
        tkconfigure(lta.lbl, state='disabled')
        tkconfigure(lta.edit, state='disabled')
      }
      if (toR(tkvars$stT)){
        tkconfigure(sta.lbl, state='normal')
        tkconfigure(sta.edit, state='normal')
        tkconfigure(sty.lbl,state='normal')
        tkconfigure(sty.edit,state='normal')
        tkconfigure(iAMA.box, state='normal')
      } else {
        tkconfigure(sta.lbl, state='disabled')
        tkconfigure(sta.edit, state='disabled')
        tkconfigure(sty.lbl,state='disabled')
        tkconfigure(sty.edit,state='disabled')
        tkconfigure(iAMA.box, state='disabled')
      #  tkconfigure(iAMAtable,state='disabled',bg=colors()[350],fg=colors()[350])
        tkgrid.remove(iAMAtable)
      }
      if (toR(tkvars$rT)){
        tkconfigure(rta.lbl, state='normal')
        tkconfigure(rta.edit, state='normal')
        tkconfigure(rhorev.lbl, state='normal')
        tkconfigure(rhorev.edit, state='normal')
      } else {
        tkconfigure(rta.lbl, state='disabled')
        tkconfigure(rta.edit, state='disabled')
        tkconfigure(rhorev.lbl, state='disabled')
        tkconfigure(rhorev.edit, state='disabled')
      }
      if (!toR(tkvars$gpost)){
        tkconfigure(g3i.intenseRadio, state = 'disabled')
        tkconfigure(g3i.nonintenseRadio, state = 'disabled')
      } else {
        tkconfigure(g3i.intenseRadio, state = 'normal')
        tkconfigure(g3i.nonintenseRadio, state = 'normal')
      }
      tkbind(ltT.box,"<Button-1>",function(){
        if (tclvalue(tkvars$ltT)=="0"){
          tkconfigure(lta.lbl, state='normal')
          tkconfigure(lta.edit, state='normal')
        } else {
          tkconfigure(lta.lbl, state='disabled')
          tkconfigure(lta.edit, state='disabled')
        }
      })
      tkbind(stT.box,"<Button-1>",function(){
        if (tclvalue(tkvars$stT)=="0"){
          tkconfigure(sta.lbl, state='normal')
          tkconfigure(sta.edit, state='normal')
          tkconfigure(sty.lbl,state='normal')
          tkconfigure(sty.edit,state='normal')
          tkconfigure(iAMA.box,state='normal')
          if (toR(tkvars$iAMA)) {
            tkgrid(iAMAtable)
          }
        } else {
          tkconfigure(sta.lbl, state='disabled')
          tkconfigure(sta.edit, state='disabled')
          tkconfigure(sty.lbl,state='disabled')
          tkconfigure(sty.edit,state='disabled')
          tkconfigure(iAMA.box,state='disabled')
      #    tkconfigure(iAMAtable,state='disabled',bg=colors()[350],fg=colors()[350])
          tkgrid.remove(iAMAtable)
        }
      })
      tkbind(rT.box,"<Button-1>",function(){
        if (tclvalue(tkvars$rT)=="0"){
          tkconfigure(rta.lbl, state='normal')
          tkconfigure(rta.edit, state='normal')
          tkconfigure(rhorev.lbl, state='normal')
          tkconfigure(rhorev.edit, state='normal')
        } else {
          tkconfigure(rta.lbl, state='disabled')
          tkconfigure(rta.edit, state='disabled')
          tkconfigure(rhorev.lbl, state='disabled')
          tkconfigure(rhorev.edit, state='disabled')
        }
      })
      tkbind(iAMA.box,"<Button-1>",function(){
        if (tclvalue(tcl(iAMA.box,"cget","-state"))=='disabled') return(F)
        if (tclvalue(tkvars$iAMA)=="1"){
      #    tkconfigure(iAMAtable, state='disabled', bg=colors()[350], fg=colors()[350])
          tkgrid.remove(iAMAtable)
      #    tkconfigure(iAMAtable, rows = 1)
        } else {
          tkgrid(iAMAtable)
      #    tkconfigure(iAMAtable, rows = length(iAMAdata)/2-1)
      #    tkconfigure(iAMAtable, state='normal', bg='white', fg='black')
      #    tcl(iAMAtable, "tag", "celltag", "readonly", cellind(toR(tkindex(iAMAtable, "end", "row")),1))
        }
      })

      tkbind(nyr.edit, "<KeyRelease>", function(){
        v<-suppressWarnings(as.numeric(toR(tkvars$nyr)))
        if (val.numeric(v) && val.integer(v) && val.gt0(v)){
          tkconfigure(nyr.edit, bg='white')
          nyrok<<-T
        } else {
          tkconfigure(nyr.edit, bg=colors()[652])
          nyrok<<-F
        }
      })
      tkbind(Tau.edit, "<KeyRelease>", function(){
        v<-suppressWarnings(as.numeric(toR(tkvars$Tau)))
        if (val.numeric(v) && val.gt0(v)){
          tkconfigure(Tau.edit, bg='white')
          Tauok<<-T
        } else {
          tkconfigure(Tau.edit, bg=colors()[652])
          Tauok<<-F
        }
      })
      tkbind(nsim.edit, "<KeyRelease>", function(){
        v<-suppressWarnings(as.numeric(toR(tkvars$nsim)))
        if (val.numeric(v) && val.integer(v) && val.gt0(v)){
          tkconfigure(nsim.edit, bg='white')
          nsimok<<-T
        } else {
          tkconfigure(nsim.edit, bg=colors()[652])
          nsimok<<-F
        }
      })
      tkbind(lambda.edit, "<KeyRelease>", function(){
        v<-suppressWarnings(as.numeric(toR(tkvars$lambda)))
        if (val.numeric(v) && val.gte0(v)){
          tkconfigure(lambda.edit, bg='white')
          lambdaok<<-T
        } else {
          tkconfigure(lambda.edit, bg=colors()[652])
          lambdaok<<-F
        }
      })
      tkbind(rhoinf.edit, "<KeyRelease>", function(){
        iAMAdata[[toR(tkindex(iAMAtable,"end", "row")),1]] <- tclvalue(tkvars$rhoinf)
        v<-suppressWarnings(as.numeric(toR(tkvars$rhoinf)))
        if (val.numeric(v) && val.gte0(v)){
          tkconfigure(rhoinf.edit, bg='white')
          rhoinfok<<-T
        } else {
          tkconfigure(rhoinf.edit, bg=colors()[652])
          rhoinfok<<-F
        }
      })
      tkbind(yrs.edit, "<KeyRelease>", function(){
        v<-suppressWarnings(as.numeric(toR(tkvars$yrs)))
        if (val.numeric(v) && val.integer(v) && val.gte0(v)){
          tkconfigure(yrs.edit, bg='white')
          yrsok<<-T
        } else {
          tkconfigure(yrs.edit, bg=colors()[652])
          yrsok<<-F
        }
      })
      tkbind(g1.edit, "<KeyRelease>", function(){
        v<-suppressWarnings(as.numeric(toR(tkvars$g1)))
        if (val.numeric(v) && val.gt0(v) && val.lt1(v)){
          tkconfigure(g1.edit, bg='white')
          g1ok<<-T
        } else {
          tkconfigure(g1.edit, bg=colors()[652])
          g1ok<<-F
        }
      })
      tkbind(g1lwr.edit, "<KeyRelease>", function(){
        v<-suppressWarnings(as.numeric(toR(tkvars$g1lwr)))
        if (val.numeric(v) && val.gt0(v) && val.lt1(v)){
          tkconfigure(g1lwr.edit, bg='white')
          g1lwrok<<-T
        } else {
          tkconfigure(g1lwr.edit, bg=colors()[652])
          g1lwrok<<-F
        }
      })
      tkbind(g1upr.edit, "<KeyRelease>", function(){
        v<-suppressWarnings(as.numeric(toR(tkvars$g1upr)))
        if (val.numeric(v) && val.gt0(v) && val.lt1(v)){
          tkconfigure(g1upr.edit, bg='white')
          g1uprok<<-T
        } else {
          tkconfigure(g1upr.edit, bg=colors()[652])
          g1uprok<<-F
        }
      })
      tkbind(g2.edit, "<KeyRelease>", function(){
        v<-suppressWarnings(as.numeric(toR(tkvars$g2)))
        if (val.numeric(v) && val.gt0(v) && val.lt1(v)){
          tkconfigure(g2.edit, bg='white')
          g2ok<<-T
        } else {
          tkconfigure(g2.edit, bg=colors()[652])
          g2ok<<-F
        }
      })
      tkbind(g2lwr.edit, "<KeyRelease>", function(){
        v<-suppressWarnings(as.numeric(toR(tkvars$g2lwr)))
        if (val.numeric(v) && val.gt0(v) && val.lt1(v)){
          tkconfigure(g2lwr.edit, bg='white')
          g2lwrok<<-T
        } else {
          tkconfigure(g2lwr.edit, bg=colors()[652])
          g2lwrok<<-F
        }
      })
      tkbind(g2upr.edit, "<KeyRelease>", function(){
        v<-suppressWarnings(as.numeric(toR(tkvars$g2upr)))
        if (val.numeric(v) && val.gt0(v) && val.lt1(v)){
          tkconfigure(g2upr.edit, bg='white')
          g2uprok<<-T
        } else {
          tkconfigure(g2upr.edit, bg=colors()[652])
          g2uprok<<-F
        }
      })
      tkbind(rta.edit, "<KeyRelease>", function(){
        v<-suppressWarnings(as.numeric(toR(tkvars$rta)))
        if (val.numeric(v) && val.gt0(v) && val.lt1(v)){
          tkconfigure(rta.edit, bg='white')
          rtaok<<-T
        } else {
          tkconfigure(rta.edit, bg=colors()[652])
          rtaok<<-F
        }
      })
      tkbind(rhorev.edit, "<KeyRelease>", function(){
        v<-suppressWarnings(as.numeric(toR(tkvars$rhorev)))
        if (val.numeric(v) && val.gt0(v)){
          tkconfigure(rhorev.edit, bg='white')
          rhorevok<<-T
        } else {
          tkconfigure(rhorev.edit, bg=colors()[652])
          rhorevok<<-F
        }
      })
      tkbind(lta.edit, "<KeyRelease>", function(){
        v<-suppressWarnings(as.numeric(toR(tkvars$lta)))
        if (val.numeric(v) && val.gte0(v) && val.lte1(v)){
          tkconfigure(lta.edit, bg='white')
          ltaok<<-T
        } else {
          tkconfigure(lta.edit, bg=colors()[652])
          ltaok<<-F
        }
      })
      tkbind(sty.edit, "<KeyRelease>", function(){
        v<-suppressWarnings(as.numeric(toR(tkvars$sty)))
        if (val.numeric(v) && val.integer(v) && val.gt0(v)){
          tkconfigure(sty.edit, bg='white')
          styok<<-T
        } else {
          tkconfigure(sty.edit, bg=colors()[652])
          styok<<-F
        }
      })
      tkbind(sta.edit, "<KeyRelease>", function(){
        v<-suppressWarnings(as.numeric(toR(tkvars$sta)))
        if (val.numeric(v) && val.gte0(v) && val.lte1(v)){
          tkconfigure(sta.edit, bg='white')
          staok<<-T
        } else {
          tkconfigure(sta.edit, bg=colors()[652])
          staok<<-F
        }
      })
      tkbind(gpost.box, "<Button-1>", function(){
        if (toR(tkvars$gpost)){
          tkconfigure(g3i.intenseRadio, state = 'disabled')
          tkconfigure(g3i.nonintenseRadio, state = 'disabled')
        } else {
          tkconfigure(g3i.intenseRadio, state = 'normal')
          tkconfigure(g3i.nonintenseRadio, state = 'normal')
        }
      })
#      tcl("wm", "attributes", scexModule, topmost=T)
      tkbind(scexModule,"<Button-1>", function() tcl("wm", "attributes", scexModule, topmost=F))
      tkwm.protocol(scexModule, "WM_DELETE_WINDOW", function(){
        tkdestroy(scexModule)
        scexsave(.Rvar$scexPrevious)
        tkwm.deiconify(.Rvar$EoA)
      })
    },
    allscexok = function(){
      tableok <<- F
      nyrok <<- T
      Tauok <<- T
      nsimok <<- T
      lambdaok <<- T
      rhoinfok <<- T
      yrsok <<- T
      g1ok <<- T
      g1lwrok <<- T
      g1uprok <<- T
      g2ok <<- T
      g2lwrok <<- T
      g2uprok <<- T
      rtaok <<- T
      rhorevok <<- T
      ltaok <<- T
      styok <<- T
      staok <<- T
    },
    feedR.scex = function(tkvars){
      scextmp<-list()
      for (nm in names(scexDefault)) scextmp[[nm]] <-toR(tkvars[[nm]])
      scextmp$iAMAschedule<-array(dim=c(as.numeric(tclvalue(tkindex(iAMAtable,"end","row"))),as.numeric(tclvalue(tkindex(iAMAtable,"end","col")))))
      for (rowi in 1:dim(scextmp$iAMAschedule)[1]){
        for (coli in 1:dim(scextmp$iAMAschedule)[2]){
          scextmp$iAMAschedule[rowi,coli]<-toR(iAMAdata[[rowi,coli]])
        }
      }
      return(scextmp)
    },
    iAMAtablechk = function(){
      rowm<-as.numeric(tclvalue(tkindex(iAMAtable,"end","row")))
      colm<-as.numeric(tclvalue(tkindex(iAMAtable,"end","col")))
      tableok<<-T
      for (i in 1:rowm){
        for (j in 1:colm){
          tableok<<-tableok*iAMACellChk(tclvalue(iAMAdata[[i,j]]),i,j)
        }
      }
      if (!tableok) return(F)
      return(T)
    #  if (toR(iAMAdata[[rowm,colm-1]]) != 0){
    #
    #  }
    },
    iAMACellChk = function(val, rowm, colm){
      if (colm == toR(tkindex(iAMAtable,"end", "row")) && rowm == 1) return(T)
      if (val==''){
        setF(iAMAtable, rowm,colm)
        return(F)
      }
      val<-suppressWarnings(as.numeric(val))
      if (length(val)==0 || is.na(val)){
        setF(iAMAtable, rowm, colm)
        return(F)
      }
      if (colm == 1){ # then val
        if (val < 0){
          setF(iAMAtable, rowm, colm)
          return(F)
        } else {
          setT(iAMAtable, rowm, colm)  # tag the cell
          return(T)
        }
      }
      if (colm == 2) { # years of monitoring must be positive integer
        if (val < 0 || val != round(val)){
          setF(iAMAtable, rowm, colm)
          return(F)
        } else {
          setT(iAMAtable, rowm, colm)  # tag the cell
          return(T)
        }
      }
      if (colm==2) { # integer years of intensive search
        if (val < 0 || val!= round(val)) {
          setF(iAMAtable, rowm, colm) # negative or not an integer
          return(F)
        }
        setT(iAMAtable, rowm, colm)
        return(T)
      }
    }
  )
)
scexsave<-function(scex) save(scex, file = paste0(.Rvar$datadir,"/scexPrevious.Rdata"))

scex1Form<-R6::R6Class("scex1Form",
  portable = FALSE,
  public = list(
    # tcl variables:
    scex1Module = NA, scexinData = NA, scexinTable = NA,
    scexoutSuperLabel = NA, scexoutSuperTable = NA, scexoutData = NA, scexoutTable = NA,
    # R variables:
    wid = 5, # width of (many of) the data entry boxes
    iAMAschedule = NA, iAMAtable = NA, iAMAdata = NA,
    scextmp = list(), steps = NA,
    # data checks...
    tableok = T, nyrok = T, Tauok = T, nsimok = T, lambdaok = T, rhoinfok = T,
    yrsok = T, g1ok = T, g1lwrok = T, g1uprok = T, g2ok = T, g2lwrok = T, g2uprok = T,
    rtaok = T, rhorevok = T, ltaok = T, styok = T, staok = T,
    viewone = 0,

    # initial values
    initialize = function(scexdat) {      # build the form with initialization of input parameters
      tcltk::tcl("option","add","*tearOff",0)
      scex1Module <<- tktoplevel(); tkgrab.set(scex1Module);  tkfocus(scex1Module)
      tkwm.title(scex1Module,paste0("EoA, v", .Rvar$VER, " - Scenario Explorer, One Example"))
      tkwm.resizable(scex1Module,0,0)
      iAMAschedule <<- scexdat$iAMAschedule
      steps<<-ifelse(scexdat$iAMA, dim(scexdat$iAMAschedule)[1], 0)
      scexinData <<- tclArray()
      with (scexdat, {
        scexinData[[0,0]]<-as.tclObj("Input parameters",drop=T)
        scexinData[[1,0]]<-as.tclObj('Permit years',drop=T); scexinData[[1,1]]<-nyr
        scexinData[[2,0]]<-as.tclObj('Total permitted take (\u03a4)',drop=T); scexinData[[2,1]]<-Tau
        scexinData[[3,0]]<-as.tclObj('Baseline fatality rate (\u03bb)',drop=T); scexinData[[3,1]]<-lambda
        scexinData[[4,0]]<-as.tclObj('g (intensive)',drop=T); scexinData[[4,1]]<-as.tclObj(paste0(g1," (", g1lwr, ", ", g1upr,")"), drop=T)
        scexinData[[5,0]]<-as.tclObj('g (non-intensive)',drop=T); scexinData[[5,1]]<-as.tclObj(paste0(g2," (", g2lwr, ", ", g2upr,")"), drop=T)
        scexinData[[6,0]]<-as.tclObj('long-term trigger (\u03b1)',drop=T); scexinData[[6,1]]<-lta
        scexinData[[7,0]]<-as.tclObj('long-term trigger (\u03c1\u221e)',drop=T); scexinData[[7,1]]<-rhoinf
        scexinData[[8,0]]<-as.tclObj('short-term trigger (\u03c4)',drop=T); scexinData[[8,1]]<-ifelse(stTon=='same',round(Tau/nyr,2),tau)
        scexinData[[9,0]]<-as.tclObj('short-term trigger (\u03b1)',drop=T); scexinData[[9,1]]<-sta
        scexinData[[10,0]]<-as.tclObj('short-term trigger (term)',drop=T); scexinData[[10,1]]<-sty
        if (iAMA) {
          scexinData[[11,0]]<-as.tclObj('Incremental AMAs:',drop=T)
          scexinData[[11,1]]<-as.tclObj('\u03c1, monitoring',drop=T)
          steps<-dim(scexdat$iAMAschedule)[1] #as.numeric(tclvalue(tkindex(iAMAtable,'end','row')))
          for (i in 1:steps){
            scexinData[[11+i,0]]<-as.tclObj(paste0("  step ", i),drop=T)
    #        scexinData[[10+i,1]]<-as.tclObj(paste0(toR(iAMAdata[[i,1]]),", ", toR(iAMAdata[[i,2]])),drop=T)
            scexinData[[11+i,1]]<-as.tclObj(paste0(iAMAschedule[i,1],", ", iAMAschedule[i,2]),drop=T)
          }
        } else {
          steps<-0
          scexinData[[11,0]]<-as.tclObj('Incremental AMAs:',drop=T)
          scexinData[[11,1]]<- 'NA'
        }
        if (scexdat$rT){
          scexinData[[12+steps, 0]]<-as.tclObj('Reversion trigger (\u03b1)',drop=T); scexinData[[12+steps,1]]<-rta
          scexinData[[12+steps+1, 0]]<-as.tclObj('Reversion trigger (\u03c1\u2080)',drop=T); scexinData[[12+steps+1,1]]<-rhorev
        } else {
          scexinData[[12+steps, 0]]<-as.tclObj('Reversion trigger (\u03b1)',drop=T); scexinData[[12+steps,1]]<-'NA'
          scexinData[[12+steps+1, 0]]<-as.tclObj('Reversion trigger (\u03c1\u2080)',drop=T); scexinData[[12+steps+1,1]]<-'NA'
        }
        scexinTable <<- tcltk2::tk2table(scex1Module,
          rows = steps+14,
          cols = 2,
          selectmode = "extended",
          variable = scexinData,
          titlecols = 1,
          titlerows = 1,
          resizeborders = "none",
          multiline = 1,
          rowseparator='\n',
          colseparator='\t',
          selecttitle = 1
        )
        tcl(scexinTable,"tag","configure","names",anchor="w")
        tcl(scexinTable,"tag", "col", "names",0)
        tcl(scexinTable, "tag", "configure", "readonly",state='disabled', bg = colors()[415])
        tcl(scexinTable,"tag", "col", "readonly", 1)
        tcl(scexinTable, "width", 0, 25)
        tcl(scexinTable, "width", 1, 20)
        ###################### output table --->
        scexoutSuperLabel<<-tclArray()
        scexoutSuperLabel[[0,0]]<-"Output:"
        scexoutSuperLabel[[0,1]]<-"Annual"
        scexoutSuperLabel[[0,2]]<-'Cumulative'
        scexoutSuperLabel[[0,3]]<-as.tclObj('Long-term',drop=T)
        scexoutSuperLabel[[0,4]]<-as.tclObj('Reversion',drop=T)
        scexoutSuperLabel[[0,5]]<-as.tclObj('Short-term',drop=T)
        scexoutSuperTable<<-tcltk2::tk2table(scex1Module,
          rows = 1,
          cols = 6,
          selectmode = "extended",
          variable = scexoutSuperLabel,
          titlerows = 1,
          resizeborders = "none",
          rowseparator='\n',
          colseparator='\t',
          selecttitle = 1
        )

        scexoutData<<-tclArray()
        scexoutData[[0,0]]<-'Year'
        scexoutData[[0,1]]<-'M'
        scexoutData[[0,2]]<-'X'
        scexoutData[[0,3]]<-'g'
        scexoutData[[0,4]]<-'\u03c1'
        scexoutData[[0,5]]<-'\u03bb'
        scexoutData[[0,6]]<-'M'
        scexoutData[[0,7]]<-'X'
        scexoutData[[0,8]]<-'g'
        scexoutData[[0,9]]<-as.tclObj('M*',drop=T)
        scexoutData[[0,10]]<-as.tclObj('p(\u03bb < \u03c4\u03c1\u2080)',drop=T)
        scexoutData[[0,11]]<-"M"
        scexoutData[[0,12]]<-"X"
        scexoutData[[0,13]]<-"g"
        scexoutData[[0,14]]<-as.tclObj("p(\u03bb > \u03c4)",drop=T)
        scexoutTable<<-tcltk2::tk2table(scex1Module,
          rows = nyr+1,
          cols = 15,
          selectmode = "extended",
          variable = scexoutData,
          titlecols = 1,
          titlerows = 1,
          resizeborders = "none",
          multiline = T,
          rowseparator='\n',
          colseparator='\t',
          selecttitle = 1
        )
        colwidths<-round(c(5,4,4,6,5,5,4,4,6,10,10,4,4,6,10)*.Rvar$charSizeAdjust)
        for (i in 1:15) tcl(scexoutTable, "width", i-1, colwidths[i])
        tcl(scexoutSuperTable,"width",0,sum(colwidths[1:2]))
        tcl(scexoutSuperTable,"width",1,sum(colwidths[3:6]))
        tcl(scexoutSuperTable,"width",2,sum(colwidths[7:9]))
        tcl(scexoutSuperTable,"width",3,sum(colwidths[10]))
        tcl(scexoutSuperTable,"width",4,sum(colwidths[11]))
        tcl(scexoutSuperTable,"width",5,sum(colwidths[12:15]))
        tcl(scexoutSuperTable, "tag", "celltag", "names", cellind(0,0))
        tcl(scexoutTable, "tag", "configure", "exceeds", fg = 'red')
        tcl(scexoutTable, "tag", "configure", "triggered", bg = 'red', fg = 'white')
        tcl(scexoutTable, "tag", "configure", "reverted", bg = colors()[256], fg = colors()[125])
        tcl(scexoutTable, "tag", "configure", "annual", state='disabled', bg = colors()[15])
        tcl(scexoutTable, "tag", "configure", "cumulative", state='disabled', bg = colors()[415])
        tcl(scexoutTable, "tag", "configure", "longterm", state='disabled', bg = colors()[15])
        tcl(scexoutTable, "tag", "configure", "reversion", state='disabled', bg = colors()[415])
        tcl(scexoutTable, "tag", "configure", "shortterm", state='disabled', bg = colors()[15])
        tcl(scexoutTable, "tag", "configure", "flash", bg = colors()[124])
        tkconfigure(scexoutTable, flashmode = T, flashtime = 1)
        for (i in 1:5) tcl(scexoutTable,"tag", "col", "annual", i)
        for (i in 6:8) tcl(scexoutTable,"tag", "col", "cumulative", i)
        tcl(scexoutTable,"tag", "col", "longterm", 9)
        tcl(scexoutTable,"tag", "col", "reversion", 10)
        for (i in 11:15) tcl(scexoutTable,"tag", "col", "shortterm", i)
        tkconfigure(scexoutTable,state='disabled',bg=colors()[415])
        scex1Buttons<-tkframe(scex1Module)
        scex1Another<-tkbutton(scex1Buttons, text = "Generate\nAnother", width=10, height=2, command=function(){
          scexCalc(.Rvar$scexPrevious, viewone = 1, self)
          tcl(scexoutTable,"selection", "clear", "all")
        })
        scex1EditParms<-tkbutton(scex1Buttons, text="Edit/Close",width=10, height=2,  command=function() {
          tkdestroy(scex1Module)
          tkwm.deiconify(.Rvar$scexWindow$scexModule)
        })
        scex1Save<-tkbutton(scex1Buttons, text="Save to file\n(.txt)",width=10, height=2, command=function(){
          # write the given data set to .txt file
          # no error-checking...
          filename <- tclvalue(tkgetSaveFile(filetypes = "{{Text} {.txt}}",defaultextension = ".txt", initialfile = '.txt', title = "Save"))
          tmp<-unlist(strsplit(filename,'/')); pathname<-paste(tmp[-length(tmp)],collapse='/')
          if (nchar(pathname)>0) .Rvar$csvpath<-pathname
          if (filename == "") return(FALSE)
          cat("Output\n", file = filename)
          cat(sprintf("    |\t%s\t\t\t\t\t%s\t\t\t\t%s\n", "Annual", "|\tCumulative", "|\tShort-term"), file = filename, append = T)
          cat(sprintf("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n",
            "Year|", "\tM", "\tX", "\tg", "\trho", "\tlambda", "\t|\tM", "\tX", "\tg", "\tM*", "\trT(a*)", "\t|\tM", "\tX", "\tg", "\tstT(a*)"),
             file = filename, append = T)
          for (i in 1:toR(tkindex(scexoutTable,"end", "row"))){
            cat(sprintf("%4s|\t%s\t%s\t%s\t%s\t%s\t|\t%s\t%s\t%s\t%s\t%s\t|\t%s\t%s\t%s\t%s\n",
              toR(scexoutData[[i,0]]),
              toR(scexoutData[[i,1]]),
              toR(scexoutData[[i,2]]),
              toR(scexoutData[[i,3]]),
              toR(scexoutData[[i,4]]),
              toR(scexoutData[[i,5]]),
              toR(scexoutData[[i,6]]),
              toR(scexoutData[[i,7]]),
              toR(scexoutData[[i,8]]),
              toR(scexoutData[[i,9]]),
              toR(scexoutData[[i,10]]),
              toR(scexoutData[[i,11]]),
              toR(scexoutData[[i,12]]),
              toR(scexoutData[[i,13]]),
              toR(scexoutData[[i,14]])),
              file = filename, append = T)
          }
          cat(paste0(rep("=",100),collapse=''), file = filename, append = T)
          cat("\n\nInput parameters\n", file = filename, append = T)
          cat(sprintf("%-20s%i\n","Years",nyr), file = filename, append = T)
          cat(sprintf("%-20s%.3g\n","Tau",Tau), file = filename, append = T)
          cat(sprintf("%-20s%.3g\n","lambda",lambda), file = filename, append = T)
          cat(sprintf("%-20s%.3f [%.3f, %.3f]\n","g (intensive)", g1, g1lwr, g1upr), file = filename, append = T)
          cat(sprintf("%-20s%.3f [%.3f, %.3f]\n","g (non-intensive)", g2, g2lwr, g2upr), file = filename, append = T)
          cat(sprintf("%-20s%.3f\n","long-term alpha", lta), file = filename, append = T)
          cat(sprintf("%-20s%.3g\n","rho_inf", rhoinf), file = filename, append = T)
          cat(sprintf("%-20s%.3g\n","short-term tau", ifelse(stTon=='same', round(Tau/nyr,2), tau)), file = filename, append = T)
          cat(sprintf("%-20s%.3g\n","short-term alpha", sta), file = filename, append = T)
          cat(sprintf("%-20s%.3g\n","short-term term", sty), file = filename, append = T)
          if (iAMA) {
            cat(sprintf("%-20s\n","Incremental AMAs"), file = filename, append = T)
            for (i in 1:steps){
              cat(sprintf("%-20s%.3g, %i\n", paste0("   step ", i), iAMAschedule[i,1], iAMAschedule[i,2]), file = filename, append = T)
            }
          }
          cat(sprintf("%-20s%.3g\n","reversion alpha", rta), file = filename, append = T)
          cat(sprintf("%-20s%.3g\n","reversion rho", rhorev), file = filename, append = T)
          return(TRUE)
        })
        tkgrid(scex1EditParms, scex1Save, scex1Another)
        tkgrid(scexinTable, rowspan = 2, columnspan = 3, sticky='n', padx = c(15,10), pady=c(40,0))
        tkgrid(scexoutSuperTable,row=0, column=3, pady=c(15,0))
        tkgrid(scexoutTable, row = 1,column=3)
        tkgrid(scex1Buttons, row = 2, column = 0, sticky='s',pady=10)
        tkwm.protocol(scex1Module, "WM_DELETE_WINDOW", function(){
            tkdestroy(scex1Module)
          # reactivate the scex window...BUT do not attempt to do so unless it exists!
            # also, do not allow destruction of the parent window unless it takes the progeny with it
            tkwm.deiconify(.Rvar$scexWindow$scexModule)
        })
      })
    }
  )
)
