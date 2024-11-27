myForm<-R6::R6Class("myForm",
  portable = FALSE,
  public = list(
    # tcl variables:
    myModule = NA, tkvars = list(), yearTable = NA, yearData = NA, projData = NA, projTable = NA, yscr = NA, tkprojYr = NA,
    yearProjections.lbl = NA, ptabFrame = NA, sides = NA,
    # text boxes
    aR.edit = NA, rho.edit = NA, aL.edit = NA, styr.edit = NA, aCI.edit = NA,
    tau.edit = NA, projYr.edit = NA, nyr.edit = NA, crlev.edit = NA,
    prho.edit = NA, gupr.edit = NA, glwr.edit = NA, g.edit = NA, Tau.edit = NA,
    g.lbl = NA, gci.lbl = NA, prho.lbl = NA, projection.lbl = NA, nyr.lbl = NA, Tau.lbl = NA,
    # radio buttons:
    myMRadio = NA, myCradio = NA, myTradio = NA,  myPradio = NA,
    projIradio = NA, projCradio = NA, projVradio = NA,
    myLRadio = NA, myCIradio = NA, myLradio = NA, myRradio = NA,

    # R variables:
    prior_f = NA, prior_M = NA, wid = 6,
    # data checks...
    gok = T, gleg = T, glwrleg = T, guprleg = T, prok = T,
    # initial values
    initialize = function(mydat) {
      # create the required tcl variables from the data input
      for (i in 1:length(myVar)){
        ind <- which(names(mydat)==myVar[i])
        tkvars[[myVar[i]]] <<- tclVar(as.tclObj(ifelse(length(ind) > 0, mydat[[ind]],""),drop=T))
      }
      for (li in 1:length(myArray)){
        ind <- which(names(mydat)==myArray[li])
        tmp<-tclArray()
        if (length(ind) == 0){ # missing data
          tmp[[0]] <- as.tclObj('', drop=T)
        } else if (length(mydat[[ind]]) == 1) { # scalar
          tmp[[0]] <- as.tclObj(ifelse(!is.na(mydat[[ind]]), mydat[[ind]], ''), drop=T)
        } else if (is.null(dim(mydat[[ind]]))) { # vector
          for (j in 1:length(mydat[[ind]])){
            tmp[[j-1]] <- as.tclObj(mydat[[ind]][j],drop=T)
          }
        } else { # matrix
          for (rowi in 1:dim(mydat[[ind]])[1])
            for (coli in 1:dim(mydat[[ind]])[2])
              tmp[[rowi-1, coli-1]] <- as.tclObj(mydat[[ind]][rowi, coli], drop = T)
        }
        tkvars[[myArray[li]]] <<- tmp
      }
      prior_f <<- mydat$prior_f
      prior_M <<- mydat$prior_M
      tkprojYr <<- length(toR(tkvars$projyears))
      ### build the components of the form
      ## toplevel and menus
      tcltk::tcl("option","add","*tearOff",0)

      myModule<<-tktoplevel()
      tkgrab.set(myModule);  tkfocus(myModule)
      tablesFrame<-tkframe(myModule); RHS<-tkframe(myModule)
      tkwm.title(myModule,paste0("EoA, v", .Rvar$VER, " - Multiple Years Module"))
      tkwm.resizable(myModule,0,0)
      #myFrame<-tkframe(myModule)
      nyear<-length(mydat$years)
      tmpyear<-list()
      ###
      MYtopMenu <- tkmenu(myModule); tkconfigure(myModule,menu=MYtopMenu)
      MYhelpMenu <- tkmenu(MYtopMenu,activebackground=colors()[125],activeforeground=colors()[109])
      MYeditMenu <- tkmenu(MYtopMenu,activebackground=colors()[125],activeforeground=colors()[109])
      tkadd(MYeditMenu,"command",label="Restore defaults",command=function(){
        .Rvar$myPrevious <- myDefault
        suppressWarnings(file.remove(paste0(.Rvar$datadir,"/myPrevious.Rdata")))
        tkdestroy(myModule)
        initialize(.Rvar$myPrevious)
      })
      tkadd(MYeditMenu,"command",label="Restore previous",command=function(){
        tkdestroy(myModule)
        initialize(.Rvar$myPrevious)
      })
      tkadd(MYeditMenu,"command",label="Save to file (.rds)",command=function() {
        nyear<-as.numeric(tclvalue(tkindex(yearTable,"end","row")))
        junk<-list()
        junk$years<-numeric(nyear)
        junk$X<-numeric(nyear)
        junk$Bb<-numeric(nyear)
        junk$Ba<-numeric(nyear)
        junk$rel_wt<-numeric(nyear)
        for (i in 1:nyear){
          years<-tclvalue(yearData[[i,1]])
          junk$years[i]<-ifelse(length(years) == 0 || is.na(years),NA, years)
          rel_wt<-suppressWarnings(as.numeric(tclvalue(yearData[[i,2]])))
          junk$rel_wt[i]<-ifelse(length(rel_wt) == 0 || is.na(rel_wt),NA, rel_wt)
          X<-suppressWarnings(as.numeric(tclvalue(yearData[[i,3]])))
          junk$X[i]<-ifelse(length(X) == 0 || is.na(X),NA, X)
          Ba<-suppressWarnings(as.numeric(tclvalue(yearData[[i,4]])))
          junk$Ba[i]<-ifelse(length(Ba) == 0 || is.na(Ba),NA, Ba)
          Bb<-suppressWarnings(as.numeric(tclvalue(yearData[[i,5]])))
          junk$Bb[i]<-ifelse(length(Bb) == 0 || is.na(Bb),NA, Bb)
        }
        junk$option<-tclvalue(tkvars$option)
        junk$Mtype<-tclvalue(tkvars$Mtype)
        junk$Ltype<-tclvalue(tkvars$Ltype)
        junk$Ptype<-tclvalue(tkvars$Ptype)
        crlev<-suppressWarnings(as.numeric(tclvalue(tkvars$crlev)))
        junk$crlev<-ifelse(length(crlev) == 0 || is.na(crlev), NA, crlev)
        nyr<-suppressWarnings(as.numeric(tclvalue(tkvars$nyr)))
        junk$nyr<-ifelse(length(nyr) == 0 || is.na(nyr), NA, nyr)
        Tau<-suppressWarnings(as.numeric(tclvalue(tkvars$Tau)))
        junk$Tau<-ifelse(length(Tau) == 0 || is.na(Tau), NA, Tau)
        g<-suppressWarnings(as.numeric(tclvalue(tkvars$g)))
        junk$g<-ifelse(length(g) == 0 || is.na(g), NA, g)
        glwr<-suppressWarnings(as.numeric(tclvalue(tkvars$glwr)))
        junk$glwr<-ifelse(length(glwr) == 0 || is.na(glwr), NA, glwr)
        gupr<-suppressWarnings(as.numeric(tclvalue(tkvars$gupr)))
        junk$gupr<-ifelse(length(gupr) == 0 || is.na(gupr), NA, gupr)
        prho<-suppressWarnings(as.numeric(tclvalue(tkvars$prho)))
        junk$prho<-ifelse(length(prho) == 0 || is.na(prho), NA, prho)
        aL<-suppressWarnings(as.numeric(tclvalue(tkvars$aL)))
        junk$aL<-ifelse(length(aL) == 0 || is.na(aL), NA, aL)
        styr<-suppressWarnings(as.numeric(tclvalue(tkvars$styr)))
        junk$styr<-ifelse(length(styr) == 0 || is.na(styr), NA, styr)
        tau<-suppressWarnings(as.numeric(tclvalue(tkvars$tau)))
        junk$tau<-ifelse(length(tau) == 0 || is.na(tau), NA, tau)
        aR<-suppressWarnings(as.numeric(tclvalue(tkvars$aR)))
        junk$aR<-ifelse(length(aR) == 0 || is.na(aR), NA, aR)
        rho<-suppressWarnings(as.numeric(tclvalue(tkvars$rho)))
        junk$rho<-ifelse(length(rho) == 0 || is.na(rho), NA, rho)
        aCI<-suppressWarnings(as.numeric(tclvalue(tkvars$aCI)))
        junk$aCI<-ifelse(length(aCI) == 0 || is.na(aCI), NA, aCI)
        if (tclvalue(tkvars$option)=="M" && tclvalue(tkvars$Mtype) == "P" && tclvalue(tkvars$Ptype) == "V"){
          nproj<-as.numeric(tclvalue(tcl(projTable, 'cget', '-rows')))-1
          junk$projyears<-numeric(nproj)
          junk$projg<-numeric(nproj)
          junk$projrho<-numeric(nproj)
          junk$projglwr<-numeric(nproj)
          junk$projgupr<-numeric(nproj)
          for (i in 1:nproj){
            junk$projyears[i]<-tclvalue(projData[[i, 1]])
            xxx<-suppressWarnings(as.numeric(tclvalue(projData[[i, 2]])))
            junk$projrho[i]<-ifelse(!is.na(xxx), xxx, tclvalue(projData[[i, 2]]))
            xxx<-suppressWarnings(as.numeric(tclvalue(projData[[i, 3]])))
            junk$projg[i]<-ifelse(!is.na(xxx), xxx, tclvalue(projData[[i, 3]]))
            xxx<-suppressWarnings(as.numeric(tclvalue(projData[[i, 4]])))
            junk$projglwr[i]<-ifelse(!is.na(xxx), xxx, tclvalue(projData[[i, 4]]))
            xxx<-suppressWarnings(as.numeric(tclvalue(projData[[i, 5]])))
            junk$projgupr[i]<-ifelse(!is.na(xxx), xxx, tclvalue(projData[[i, 5]]))
          }
        }
        saveparm(junk)
        tkwm.title(myModule,paste0(.Rvar$dataFileTitle, " - EoA, v", .Rvar$VER, " - Multiple Years Module"))
      })
      tkadd(MYeditMenu,"command",label="Read from file (.rds)",command=MYreadparm)
      tkadd(MYtopMenu, "cascade", label="Edit", menu=MYeditMenu)
      tkadd(MYhelpMenu, "command", label="Parameter Conversion", command=conversionsCalculator)
      tkadd(MYhelpMenu, "command", label="About", command = function() tkmessageBox(title='Evidence of Absence (EoA)',message=about_text))
      tkadd(MYtopMenu, "cascade", label="Help", menu=MYhelpMenu)
      #tkwm.protocol(myModule,"WM_DELETE_WINDOW",function(){ # same as the following line but a different format?
      tkbind(myModule,"<Destroy>",function() { # red X kills the window
        ### assumption is that no bad data have been previously saved to the .Rvar$myPrevious dataframe
        ## [previous error-checking should assure that is true]
        tkevent.add("<<Paste>>","<Control-v>")
      #  .Rvar$myPrevious<<-myProvisional
        tkdestroy(myModule)
        tkwm.deiconify(.Rvar$EoA)
      })
      ###########
      #### create table for holding, displaying, editing data
      ## table creation assumes valid data in .Rvar$myPrevious
      yearData<<-tclArray()
#      yearData<-tclArray()
      columnNames<-c("","Year", "\u03c1", "X", "Ba", "Bb", "g\u0302", "95% CI")
      for (i in 1:length(columnNames)) yearData[[0,i-1]]<<-strsplit(columnNames[i]," ",fixed=T)[[1]]
      for (i in 1:nyear){
        yearData[[i,0]]<<-as.tclObj('',drop=T)
        yearData[[i,1]]<<- mydat$years[i]
        yearData[[i,2]]<<- mydat$rel_wt[i]
        yearData[[i,3]]<<- mydat$X[i]
        yearData[[i,4]]<<- signif(mydat$Ba[i],4)
        yearData[[i,5]]<<- signif(mydat$Bb[i],4)
        yearData[[i,6]]<<- round(mydat$Ba[i]/(mydat$Ba[i]+mydat$Bb[i]),3)
        yearData[[i,7]]<<- as.tclObj(paste0("[", signif(qbeta(0.025,mydat$Ba[i], mydat$Bb[i]),3), ", ", signif(qbeta(0.975,mydat$Ba[i], mydat$Bb[i]),3),"]"),drop=T)
      }
      valChar<-function(S){
        actcol<-as.numeric(tclvalue(tkindex(yearTable,"active","col")))
        actrow<-as.numeric(tclvalue(tkindex(yearTable,"active","row")))
        if (length(grep("\n",S))>0){ #
          return(tcl("expr", FALSE))
        } else {
          yearCellChk(S, actrow, actcol)
          return(tcl("expr", TRUE))# but other kinds of space are not--> error-checking if good value is added
        }
      }
      yearTable<<-tcltk2::tk2table(tablesFrame,
#      yearTable<-tk2table(tablesFrame,
        rows=nyear+1,
        cols=8,
        selectmode="extended",
        variable=yearData,
        titlerows = 1,
        titlecols = 1,
        background='white',
        resizeborders="none",
        multiline = 1,
        wrap = 1,
        rowseparator='\n',
        colseparator='\t',
        validate = 1,
        vcmd=function(S) valChar(S),
        selecttitle = 1
      )
      tkevent.add("<<Paste>>","<Control-v>")
      tcl(yearTable,"tag","configure", "error",bg=colors()[652]) # cells with error have yellow background
      tcl(yearTable,"tag","configure", "cellok",bg='white')
      tcl(yearTable,"tag","configure", "active",fg='black',relief='groove')
      tcl(yearTable, "tag","configure", "readonly",state='disabled', bg = colors()[443])
      tcl(yearTable,"tag","configure","ignore",state='disabled',bg=colors()[365],fg=colors()[42])
      tcl(yearTable,"tag","col", "readonly", 6)
      tcl(yearTable,"tag","col", "readonly", 7)

      colwidths<-c(1,15,7,7,10,10,8,15)
      for(i in 1:8) {
        tcl(yearTable, "width", i - 1, colwidths[i])
      }
      tkbind(yearTable,"<Return>",function(){
        if (tclvalue(tkindex(yearTable,"active","row"))==tclvalue(tkindex(yearTable,"end","row"))){
          tkevent.generate(yearTable,"<Control-KeyPress-a>")
          return(T)
        }
        if(tclvalue(tkindex(yearTable,"active","row"))==tclvalue(tkindex(yearTable,"end","row"))){
          if (tclvalue(tkindex(yearTable,"active","col"))!=tclvalue(tkindex(yearTable,"end","col"))){
            tkevent.generate(yearTable,"<KeyPress-Right>")
          }
        } else {
          tkevent.generate(yearTable,"<KeyPress-Down>")
        }
      })
      tkbind(yearTable,"<Tab>",function(){
        if(tclvalue(tkindex(yearTable,"active","col"))!=tclvalue(tkindex(yearTable,"end","col"))){
          tkevent.generate(yearTable,"<KeyPress-Right>")
        } else {
          tkevent.generate(yearTable,"<KeyPress-Down>")
        }
      })
      tkbind(yearTable,"<1>",function(){ # mouse button is pressed...check for errors
        if (substr(tclvalue(tcl(yearTable,"curvalue")),1,1)=='.') tcl(yearTable,"curvalue",paste0(0,tclvalue(tcl(yearTable,"curvalue"))))
      })
      tkbind(yearTable,"<Control-KeyPress-d>", function(){
#        if (tclvalue(tkvars$option) == "M" && tclvalue(tkvars$Mtype)=="P" && tclvalue(tkvars$Ptype) == "V") return(F) # do not allow rows to be deleted when custom operations for future years is enabled
        if (tclvalue(tkvars$option) == "M" && tclvalue(tkvars$Mtype)=="P") return(F) # do not allow rows to be deleted when custom operations for future years is enabled
      # delete current row (if there is only one row, just erase the data)
        actrow<-tclvalue(tkindex(yearTable, "active", "row"))
        tkdelete(yearTable, "rows", actrow, 1)
      })
      tkevent.delete("<<Paste>>","<Control-v>") # hijack the normal activity of the ctrl + v paste command
      junk<-NA
      tkbind(yearTable,"<Control-Key-v>", function(){
      # after pasting data from the clipboard, read the prior table data into a buffer
      # if there are errors in the new data, give error message and replace
        if (as.numeric(tclvalue(tcl(yearTable,"index","active","col")))>5)   return(F)
      #  junk<<-tkXselection.get(selection='CLIPBOARD') # grab dsta from the clipboard for preliminary analysis before pasting to the table
      #  junk<<-gsub('[^\t^\n]','',tclvalue(junk)) #remove everything but tabs and carriage returns for easy parsing of the data
        junk<-readClipboard()
        ncols<-max(nchar(gsub('[^\t]','',junk)))
        if (ncols+as.numeric(tclvalue(tcl(yearTable,"index","active","col")))>5) {
          tkmessageBox(message="Error: cannot paste over summary columns",icon='error')
          return(F)
        }
        yrs<-as.numeric(tclvalue(tkindex(yearTable,"end","row")))
        nyr<-nyrchk()
        if (nyr && tclvalue(tkvars$option)=="M" && tclvalue(tkvars$Mtype)=="P" && tclvalue(tkvars$Ptype)=="V") {
        # then check to see if there would be too many rows
          if (as.numeric(tclvalue(tcl(yearTable,"index","active","row")))+length(junk)-1 > yrs){
            tkmessageBox(message="Error: cannot add new lines of data when 'Projection of future mortality and estimates' option is selected",icon="error")
            return(F)
          }
        }
        tkconfigure(yearTable,rows=max(yrs+1,as.numeric(tclvalue(tcl(yearTable,"index","active","row")))+length(junk)))
        tkevent.generate(yearTable,"<<Paste>>")
        # check table for errors and missing values...
        yearTableChk()
      })
      tkbind(yearTable,"<Control-KeyPress-a>", function(){
      # add a row
#        if (tclvalue(tkvars$option) == "M" && tclvalue(tkvars$Mtype)=="P" && tclvalue(tkvars$Ptype) == "V") return(F) # do not allow rows to be deleted when custom operations for future years is enabled
#        if (tclvalue(tkvars$Mtype)=="P" && tclvalue(tkvars$Ptype)=="V") return(F) # do not allow rows to be deleted when custom operations for future years is enabled
        if (tclvalue(tkvars$option) == "M" && tclvalue(tkvars$Mtype)=="P") return(F) # do not allow rows to be deleted when custom operations for future years is enabled
        tkinsert(yearTable, "rows", "end", 1)
        rowind<-as.numeric(tclvalue(tcl(yearTable,"index","end","row")))
        for (i in 1:5){
          tcl(yearTable,"tag","celltag","error", as.tclObj(paste0(rowind,',', i)))
        }
        tcl(yearTable,"activate", as.tclObj(paste0(rowind,',', 1)))
        yearData[[rowind,7]]<<- "NA"
        yearData[[rowind,6]]<<-"NA"
        yearData[[as.numeric(tclvalue(tcl(yearTable,"index","end","row"))),1]]<<-1+as.numeric(tclvalue(yearData[[as.numeric(tclvalue(tcl(yearTable,"index","end","row")))-1,1]]))
      })
      movechars<-c('Up','Down','Left','Right','Tab','Return','Shift_L','Alt_L')
      tkbind(yearTable,"<Key>",function(K){
        if (K %in% movechars) {
          if (substr(tclvalue(tcl(yearTable,"curvalue")),1,1)=='.') tcl(yearTable,"curvalue",paste0(0,tclvalue(tcl(yearTable,"curvalue"))))
        }
      })
      ##########
      # the following vars are a little tricky to handle because they are not stored but change with the table
      ny<-length(mydat$years)

      #### presentation
      myOptionsFrame<-ttklabelframe(RHS,text="   Options",padding=10)
      myActionsFrame<-ttklabelframe(RHS,text="   Actions",padding=10)
      ### within the options frame
      myMFrame<-ttklabelframe(myOptionsFrame, text="   Fatalities")
      myLFrame<-ttklabelframe(myOptionsFrame, text="   Average Rate")
      ## within in the fatalities frame (myMFrame)
      mframe<-tkframe(myMFrame)
      myMRadio<<-tkradiobutton(mframe, variable = tkvars$option, text = "Estimate M", value = 'M')
      crlev.lbl<-tklabel(mframe,  text="    Credibility level (1 - \u03b1)") ; crlev.edit<<-tkentry(mframe, textvariable = tkvars$crlev, bg='white', width = wid, justify = 'right')
      tkgrid(myMRadio,crlev.lbl,crlev.edit)
      myCradio<<-tkradiobutton(myMFrame, variable=tkvars$Mtype,text="Total mortality",value="C", command = function(){
        tkconfigure(yearTable,state='normal')
        tkconfigure(nyr.edit,state='normal')
        tkconfigure(projIradio, state ='disabled')
        tkconfigure(projCradio, state ='disabled')
        tkconfigure(projVradio, state ='disabled')
        tkconfigure(g.edit, state ='disabled')
        tkconfigure(g.lbl, state='disabled')
        tkconfigure(gci.lbl, state='disabled')
        tkconfigure(glwr.edit, state ='disabled')
        tkconfigure(gupr.edit, state ='disabled')
        tkconfigure(prho.lbl, state='disabled')
        tkconfigure(prho.edit, state='disabled')
        tkgrid.remove(ptabFrame,yearProjections.lbl)
        tkconfigure(projection.lbl,state='disabled')
        tkconfigure(nyr.lbl,state='disabled')
        tkconfigure(Tau.lbl,state='disabled')
#        tkconfigure(nyr.edit,state='disabled')
        tkconfigure(Tau.edit,state='disabled')
      })
      sidesFrame<-tkframe(myMFrame)
      sides <<- tclVar(1)
      oneside <- tkradiobutton(sidesFrame, text = "One-sided CI (M*)", variable = sides, value = 1)
      twoside <- tkradiobutton(sidesFrame, text = "Two-sided CI", variable = sides, value = 2)
      tkgrid(oneside, sticky='w')#, padx=c(10,0))
      tkgrid(twoside, sticky='w')#, padx=c(10,0))
      tkgrid(mframe, row = 0, column = 0, sticky= 'w', columnspan = 2)
      tkgrid(sidesFrame, row = 1, column = 1, sticky='w')
      projectFrame<-ttklabelframe(myMFrame, text="   Project parameters")
      # within the project frame
      projparms<-tkframe(projectFrame)
      nyr.lbl<<-tklabel(projparms, text="Total years in project") ; nyr.edit<<-tkentry(projparms, textvariable=tkvars$nyr, bg='white', width=wid, justify='right')
      Tau.lbl<<-tklabel(projparms, text="Mortality threshold (\u03a4)") ; Tau.edit<<-tkentry(projparms,textvariable=tkvars$Tau,bg='white', width=wid,justify='right')
      tkgrid(nyr.lbl, nyr.edit,sticky='w')
      tkgrid(Tau.lbl, Tau.edit,sticky='w')
      myTradio<<-tkradiobutton(projectFrame, variable=tkvars$Mtype,text="Track past mortality",value="T", command=function()tkconfigure(yearTable,state='normal'))

      projectionFrame<-tkframe(projectFrame)
      myPradio<<-tkradiobutton(projectionFrame, variable=tkvars$Mtype,text="Projection of future mortality and estimates", value="P", command = function() myPdo())
      projection.lbl<<-tklabel(projectionFrame,text='Future monitoring and operations')
      projIradio<<-tkradiobutton(projectionFrame,variable=tkvars$Ptype, value="I", text = "g and \u03c1 unchanged from most recent year", command = function(){
        myProjFill()
        projTableChk()
      })
      projCradio<<-tkradiobutton(projectionFrame,variable=tkvars$Ptype, value="C", text = "g and \u03c1 constant, different from most recent year", command = function(){
        myProjFill()
        projTableChk()
      })
      gFrame<-tkframe(projectionFrame)
      g.lbl<<-tklabel(gFrame,text = "  g"); g.edit<<-tkentry(gFrame,textvariable=tkvars$g,bg='white',justify='right',width=wid)
      gci.lbl<<-tklabel(gFrame,text = "    95% CI:")
      glwr.edit<<-tkentry(gFrame,textvariable=tkvars$glwr,bg='white',justify='right',width=wid)
      gupr.edit<<-tkentry(gFrame,textvariable=tkvars$gupr,bg='white',justify='right',width=wid)
      prho.lbl<<-tklabel(gFrame,text=" \u03c1")
      prho.edit<<-tkentry(gFrame,textvariable=tkvars$prho,bg='white',justify='right',width=wid)

      tkgrid(g.lbl,padx=5)
      tkgrid(g.edit, row = 0, column = 1, sticky='w')
      tkgrid(gci.lbl,row = 0, column = 2,padx=5)
      tkgrid(glwr.edit, row = 0, column = 3,sticky='w')
      tkgrid(gupr.edit,row = 0, column = 4,sticky='w')
      tkgrid(prho.lbl,row=0,column=5,sticky='w',padx=c(10,3))
      tkgrid(prho.edit, row=0,column=6,sticky='w')
      projVradio<<-tkradiobutton(projectionFrame,variable=tkvars$Ptype, value="V", text = "g and \u03c1 vary among future years", command = function(){
        myProjFill(mydat)
        projTableChk()
      })
      projYr.lbl<-tklabel(projectionFrame,text = "  Years to project")
      projYr.edit<<-tkentry(projectionFrame,textvariable=tkprojYr,bg='white',justify='right',width=wid)
      tkgrid(myPradio,sticky='w',columnspan=6)
      tkgrid(projection.lbl,sticky='w',padx=20,columnspan=6)
      tkgrid(projIradio,sticky='w',padx=24,columnspan=6)
      tkgrid(projCradio,sticky='w',padx=24,columnspan=6)
      tkgrid(gFrame,padx=30)
      tkgrid(projVradio,sticky='w',padx=24,columnspan=6)
      # assemble the projectFrame
      tkgrid(projparms,sticky='w')
      tkgrid(myTradio,sticky='w')
      tkgrid(myPradio,sticky='w')
      tkgrid(projectionFrame)
      tkgrid(nyr.lbl,nyr.edit,sticky='w')
      tkgrid(Tau.lbl,Tau.edit)
      # assemble the fatalities frame
#      tkgrid(mframe,sticky='w')
      #tkgrid(myMRadio,sticky='w')
      #tkgrid(aM.lbl,sticky='w',row=0,column=1)
      #tkgrid(aM.edit,row=0,column=2,sticky='w')
      tkgrid(myCradio, row=1, column=0, sticky='w', padx=c(15,3))
      tkgrid(projectFrame,columnspan=3,sticky='w',padx=10)

      ## within the rate frame (myLFrame)
      myLRadio<<-tkradiobutton(myLFrame,variable=tkvars$option,text="Estimate average annual fatality rate (\u03bb)", value="L", command = function() tkconfigure(yearTable,state='normal'))
      tau.lbl<-tklabel(myLFrame, text="Annual rate theshold (\u03c4)"); tau.edit<<-tkentry(myLFrame, textvariable=tkvars$tau, bg='white',width=wid,justify='right')
      myCIradio<<-tkradiobutton(myLFrame, variable=tkvars$Ltype,text="Credibility level for CI (1-\u03b1)", value="ciT", justify='left') # CI...need alpha value for two-sided CI
      aCI.edit<<-tkentry(myLFrame, textvariable=tkvars$aCI, bg='white',width=wid,justify='right')
      myLradio<<-tkradiobutton(myLFrame, variable=tkvars$Ltype,text="Short-term rate (\u03bb > \u03c4)", value="stT", justify="left")
      styr.lbl<-tklabel(myLFrame, text="   Term:") ; styr.edit<<-tkentry(myLFrame, textvariable=tkvars$styr,bg='white',width=wid,justify='right')
      aL.lbl<-tklabel(myLFrame, text="\u03b1") ; aL.edit<<-tkentry(myLFrame, textvariable=tkvars$aL,bg='white',width=wid,justify='right')
      myRradio<<-tkradiobutton(myLFrame, variable=tkvars$Ltype,text="Reversion test (\u03bb < \u03c1 \u03c4)", value="rT", justify= 'left')
      rho.lbl<-tklabel(myLFrame, text="\u03c1") ; rho.edit<<-tkentry(myLFrame, textvariable=tkvars$rho,bg='white',width=wid,justify='right')
      aR.lbl<-tklabel(myLFrame, text="\u03b1") ; aR.edit<<-tkentry(myLFrame, textvariable=tkvars$aR,bg='white',width=wid,justify='right')
      # assemble the rate frame (myLFrame)
      tkgrid(myLRadio,sticky='w',columnspan=3)
      tkgrid(tau.lbl)
      tkgrid(tau.edit, row = 1 , column = 1, columnspan = 2, sticky='e')
      tkgrid(myCIradio, sticky='w',row=2, padx=15, columnspan=2)
      tkgrid(aCI.edit,sticky='e', row=2, column = 1, columnspan = 2)
      tkgrid(myLradio, sticky='w',padx=15,columnspan=2)
      tkgrid(styr.lbl, sticky='e',padx=3, row = 3,column=2)
      tkgrid(styr.edit, sticky='e',padx=3, row = 3,column=3)
      tkgrid(aL.lbl,sticky='e',padx=3, row = 3,column=4)
      tkgrid(aL.edit,sticky='e',padx=3, row = 3,column=5)
      tkgrid(myRradio,sticky='w',padx=15,columnspan=2)
      tkgrid(rho.lbl,sticky='e',padx=3, row = 4,column=2)
      tkgrid(rho.edit,sticky='e',padx=3, row = 4,column=3)
      tkgrid(aR.lbl,sticky='e',padx=3, row = 4,column=4)
      tkgrid(aR.edit,sticky='e',padx=3, row = 4,column=5)

      ## within the actions frame
      ## actions buttons
      myButtWid<-14
      myOnCalc<-tkbutton(myActionsFrame, text='Calculate', width = myButtWid, command=function(){
        junk<-feedMY() # this brings the form data from tcl to R after error-checking
        if (is.list(junk)) {
      # set options
          if (junk$option=="M" && junk$Mtype=="P" && junk$Ptype == "V"){# then check whether g_lwr < ghat < g_upr


            for (i in 1:length(junk$projg)){
              if (junk$projg[i] <= junk$projglwr[i]){
                setFpt(i , 4)
                tkmessageBox(message="Error in data. Cannot calculate...")
                return(F)
              }
              if (junk$projg[i] >= junk$projgupr[i]){
                setFpt(i , 5)
                tkmessageBox(message="Error in data. Cannot calculate...")
                return(F)
              }
            }
          }

          .Rvar$myPrevious$years<-junk$years
          .Rvar$myPrevious$X<-junk$X
          .Rvar$myPrevious$Ba<-junk$Ba
          .Rvar$myPrevious$Bb<-junk$Bb
          .Rvar$myPrevious$rel_wt<-junk$rel_wt
          .Rvar$myPrevious$option<-junk$option
          if (junk$option=="M"){
            .Rvar$myPrevious$Mtype<-junk$Mtype
            .Rvar$myPrevious$crlev<-junk$crlev
            if (junk$Mtype %in% c("T","P")){
              .Rvar$myPrevious$nyr<-junk$nyr
              .Rvar$myPrevious$Tau<-junk$Tau
            }
            if (junk$Mtype=="P"){
              .Rvar$myPrevious$Ptype<-junk$Ptype
              if (junk$Ptype=="C"){
                .Rvar$myPrevious$g<-junk$g
                .Rvar$myPrevious$glwr<-junk$glwr
                .Rvar$myPrevious$gupr<-junk$gupr
                .Rvar$myPrevious$prho<-junk$prho
              }
              if (junk$Ptype=="V"){
                .Rvar$myPrevious$projyears<-junk$projyears
                .Rvar$myPrevious$projrho<-junk$projrho
                .Rvar$myPrevious$projg<-junk$projg
                .Rvar$myPrevious$projglwr<-junk$projglwr
                .Rvar$myPrevious$projgupr<-junk$projgupr
              }
            }
          } else {
            .Rvar$myPrevious$Ltype<-junk$Ltype
            if (junk$Ltype %in% c('stT','rT')){
              .Rvar$myPrevious$tau<-junk$tau
              if (junk$Ltype=='stT'){
                .Rvar$myPrevious$styr<-junk$styr
                .Rvar$myPrevious$aL<-junk$aL
              }
              if (junk$Ltype=='rT'){
                .Rvar$myPrevious$rho<-junk$rho
                .Rvar$myPrevious$aR<-junk$aR
              }
            } else {
              .Rvar$myPrevious$aCI<-junk$aCI
            }
          }
        } else {
          return(F)
        }
        myCalc(.Rvar$myPrevious, tclvalue(sides))
      })
      myOnClose<-tkbutton(myActionsFrame, text="Close", width=myButtWid, command=function(){
        mysave(.Rvar$myPrevious)
        graphics.off()
        tkdestroy(myModule)
       tkwm.deiconify(.Rvar$EoA)
      })
      tkgrid(myOnCalc, myOnClose)

      ## assemble the full options frame
      tkgrid(myMFrame,sticky='w')
      tkgrid(myLFrame,sticky='w',pady=10)
      ## create the form
      yearData.lbl<-tklabel(tablesFrame,text="Past monitoring and operations data")
      yearProjections.lbl<<-tklabel(tablesFrame,text="Future monitoring and operations parameters")
      tkgrid(yearData.lbl,sticky='nw')
      tkgrid(yearTable,sticky='nw',padx=15, columnspan=2)
      tkgrid(tablesFrame,sticky='n',padx=15,pady=20)
      tkgrid(myOptionsFrame,column=1,row=0,sticky='nw',rowspan=2)
      tkgrid(myActionsFrame, column=1, pady=15)
      tkgrid(RHS,row=0,column=1,sticky='n')
      ## initial configuration the text entry fields and radio buttons
#      setstates()


      ################## radio buttons ###########################
      tkbind(myMRadio, "<Button-1>", function() {
#          tkconfigure(yearTable,state='normal')
          tkgrid.remove(ptabFrame,yearProjections.lbl)
          tkconfigure(nyr.edit,state ='normal')
          tkconfigure(Tau.edit,state ='normal')
          tkconfigure(crlev.edit, state ='normal')
          tkconfigure(myCradio,state ='normal')
          tkconfigure(myTradio,state ='normal')
          tkconfigure(myPradio,state ='normal')
          tkconfigure(tau.lbl,state='disabled')
          tkconfigure(styr.lbl,state='disabled')
          tkconfigure(aL.lbl,state='disabled')
          tkconfigure(rho.lbl,state='disabled')
          tkconfigure(aR.lbl,state='disabled')
          if (tclvalue(tkvars$Mtype) == "P"){
            tkconfigure(projIradio, state ='normal')
            tkconfigure(projCradio, state ='normal')
            tkconfigure(projVradio, state ='normal')
            gstate<-ifelse(tclvalue(tkvars$Ptype)=="C",'normal','disabled')#tkvars$Ptype<-tclVar("I") # I, C, V:
            tkconfigure(g.edit, state =gstate)
            tkconfigure(g.lbl, state=gstate)
            tkconfigure(gci.lbl, state=gstate)
            tkconfigure(glwr.edit, state =gstate)
            tkconfigure(gupr.edit, state =gstate)
            tkconfigure(g.edit, state ='normal')
            tkconfigure(g.lbl, state='normal')
            tkconfigure(gci.lbl, state='normal')
            tkconfigure(glwr.edit, state ='normal')
            tkconfigure(gupr.edit, state ='normal')
            tkconfigure(prho.lbl, state='normal')
            tkconfigure(prho.edit, state='normal')
            tkconfigure(projection.lbl,state='normal')
            tkconfigure(nyr.lbl,state='normal')
            tkconfigure(Tau.lbl,state='normal')
            tkconfigure(nyr.edit,state='normal')
            tkconfigure(Tau.edit,state='normal')
#            if (tclvalue(tkvars$Ptype)=="V") {
              gridProjTable(as.numeric(tclvalue(tkvars$nyr))-as.numeric(tclvalue(tcl(yearTable,'index','end','row'))))
#            }
          } else {
            tkconfigure(projIradio, state ='disabled')
            tkconfigure(projCradio, state ='disabled')
            tkconfigure(projVradio, state ='disabled')
            tkconfigure(g.edit, state ='disabled')
            tkconfigure(g.lbl, state='disabled')
            tkconfigure(gci.lbl, state='disabled')
            tkconfigure(glwr.edit, state ='disabled')
            tkconfigure(gupr.edit, state ='disabled')
            tkconfigure(prho.lbl, state='disabled')
            tkconfigure(prho.edit, state='disabled')
            if (tclvalue(tkvars$Mtype)=="T"){
              tkconfigure(nyr.lbl,state='normal')
              tkconfigure(Tau.lbl,state='normal')
              tkconfigure(nyr.edit,state='normal')
              tkconfigure(Tau.edit,state='normal')
            } else {
              tkconfigure(nyr.lbl,state='disabled')
              tkconfigure(Tau.lbl,state='disabled')
#              tkconfigure(nyr.edit,state='disabled')
              tkconfigure(Tau.edit,state='disabled')
            }
          }

          tkconfigure(myCIradio,state='disabled')
          tkconfigure(myLradio,state='disabled')
          tkconfigure(myRradio,state='disabled')

          tkconfigure(aL.edit,  state = 'disabled')
          tkconfigure(styr.edit, state = 'disabled')
          tkconfigure(tau.edit, state ='disabled')
          tkconfigure(tau.lbl, state ='disabled')
          tkconfigure(aR.edit, state ='disabled')
          tkconfigure(rho.edit, state ='disabled')
          tkconfigure(aCI.edit, state ='disabled')
      })
#      tkbind(myCradio,"<Button-1>",function(){
#print("myCradio clicked..."); flush.console()
#        if (tclvalue(tkvars$option) != "M") return(F)
#        tkconfigure(yearTable,state='normal')
#print("...yearTable normalized"); flush.console()
#        tkconfigure(nyr.edit,state='normal')
#        tkconfigure(projIradio, state ='disabled')
#        tkconfigure(projCradio, state ='disabled')
#        tkconfigure(projVradio, state ='disabled')
#        tkconfigure(g.edit, state ='disabled')
#        tkconfigure(g.lbl, state='disabled')
#        tkconfigure(gci.lbl, state='disabled')
#        tkconfigure(glwr.edit, state ='disabled')
#        tkconfigure(gupr.edit, state ='disabled')
#        tkconfigure(prho.lbl, state='disabled')
#        tkconfigure(prho.edit, state='disabled')
#        tkgrid.remove(ptabFrame,yearProjections.lbl)
#        tkconfigure(projection.lbl,state='disabled')
#        tkconfigure(nyr.lbl,state='disabled')
#        tkconfigure(Tau.lbl,state='disabled')
##        tkconfigure(nyr.edit,state='disabled')
#        tkconfigure(Tau.edit,state='disabled')
#      })
      tkbind(myTradio,"<Button-1>",function(){
        if (tclvalue(tkvars$option) != "M") return(F)
#        tkconfigure(yearTable,state='normal')
        tkconfigure(nyr.edit,state='normal')
        tkconfigure(projIradio, state ='disabled')
        tkconfigure(projCradio, state ='disabled')
        tkconfigure(projVradio, state ='disabled')
        tkconfigure(g.edit, state ='disabled')
        tkconfigure(g.lbl, state='disabled')
        tkconfigure(gci.lbl, state='disabled')
        tkconfigure(glwr.edit, state ='disabled')
        tkconfigure(gupr.edit, state ='disabled')
        tkconfigure(prho.lbl, state='disabled')
        tkconfigure(prho.edit, state='disabled')
        tkgrid.remove(ptabFrame,yearProjections.lbl)
        tkconfigure(projection.lbl,state='disabled')
        tkconfigure(nyr.lbl,state='normal')
        tkconfigure(Tau.lbl,state='normal')
        tkconfigure(nyr.edit,state='normal')
        tkconfigure(Tau.edit,state='normal')
      })
      tkbind(projIradio,"<Button-1>",function(){
        if (tclvalue(tkvars$option) != "M" | tclvalue(tkvars$Mtype) != "P") return(F)#tkvars$Mtype<-tclVar(.Rvar$myPrevious$Mtype)
        tkconfigure(nyr.edit,state='normal')
        tkconfigure(g.edit, state ='disabled')
        tkconfigure(g.lbl, state='disabled')
        tkconfigure(gci.lbl, state='disabled')
        tkconfigure(glwr.edit, state ='disabled')
        tkconfigure(gupr.edit, state ='disabled')
        tkconfigure(prho.lbl, state='disabled')
        tkconfigure(prho.edit, state='disabled')
        tkconfigure(nyr.edit,state='normal')
        tkconfigure(projTable, state='disabled')
      })
      tkbind(projCradio,"<Button-1>",function(){
        if (tclvalue(tkvars$option) != "M" | tclvalue(tkvars$Mtype) != "P") return(F)
        tkconfigure(nyr.edit,state='normal')
        tkconfigure(g.edit, state ='normal')
        tkconfigure(g.lbl, state='normal')
        tkconfigure(gci.lbl, state='normal')
        tkconfigure(glwr.edit, state ='normal')
        tkconfigure(gupr.edit, state ='normal')
        tkconfigure(nyr.edit,state='normal')
        tkconfigure(prho.lbl, state='normal')
        tkconfigure(prho.edit, state='normal')
        tkconfigure(projTable, state='disabled')
      })
      tkbind(projVradio,"<Button-1>",function(){
        if (tclvalue(tcl(projVradio,"cget","-state"))=="disabled") return(F)
        tkconfigure(g.edit, state ='disabled')
        tkconfigure(g.lbl, state='disabled')
        tkconfigure(gci.lbl, state='disabled')
        tkconfigure(glwr.edit, state ='disabled')
        tkconfigure(gupr.edit, state ='disabled')
        tkconfigure(prho.lbl, state='disabled')
        tkconfigure(prho.edit, state='disabled')
#        tkconfigure(nyr.edit,state='disabled')
        tkconfigure(projTable, state='normal')
#        rown<-as.numeric(tclvalue(tkvars$nyr))-as.numeric(tclvalue(tcl(yearTable,'index','end','row')))
#        gridProjTable(rown)
      })
      tkbind(myLRadio,"<Button-1>", function() {
          tkgrid.remove(ptabFrame,yearProjections.lbl)
          tkconfigure(nyr.edit,state='normal')
          tkconfigure(nyr.edit,state='disabled')
          tkconfigure(Tau.edit,state='disabled')
          tkconfigure(crlev.edit, state='disabled')
          tkconfigure(myCradio,state='disabled')
          tkconfigure(myTradio,state='disabled')
          tkconfigure(myPradio,state='disabled')
          tkconfigure(projIradio, state ='disabled')
          tkconfigure(projCradio, state ='disabled')
          tkconfigure(projVradio, state ='disabled')
          tkconfigure(g.edit, state ='disabled')
          tkconfigure(g.lbl, state='disabled')
          tkconfigure(gci.lbl, state='disabled')
          tkconfigure(glwr.edit, state ='disabled')
          tkconfigure(gupr.edit, state ='disabled')
          tkconfigure(prho.lbl, state='disabled')
          tkconfigure(prho.edit, state='disabled')
          tkconfigure(nyr.lbl, state='disabled')
          tkconfigure(Tau.lbl, state='disabled')
          tkconfigure(projection.lbl, state='disabled')

          tkconfigure(myCIradio,state='normal')
          tkconfigure(myLradio,state='normal')
          tkconfigure(myRradio,state='normal')
        if (tclvalue(tkvars$Ltype)=="stT") {
          tkconfigure(aCI.edit,state='disabled')
          tkconfigure(rho.edit,state='disabled')
          tkconfigure(aL.edit,state='normal')
          tkconfigure(styr.edit,state='normal')
          tkconfigure(aR.edit,state='disabled')
          tkconfigure(aL.lbl,state="normal")
          tkconfigure(styr.lbl,state="normal")
          tkconfigure(rho.lbl,state='disabled')
          tkconfigure(aR.lbl,state='disabled')
          tkconfigure(tau.edit,state="normal")
          tkconfigure(tau.lbl,state="normal")

        } else if (tclvalue(tkvars$Ltype)=="ciT"){
          tkconfigure(aCI.edit,state='normal')
          tkconfigure(rho.edit,state='disabled')
          tkconfigure(aL.edit,state='disabled')
          tkconfigure(styr.edit,state='disabled')
          tkconfigure(aR.edit,state='disabled')
          tkconfigure(tau.edit,state="disabled")
          tkconfigure(aL.lbl,state="disabled")
          tkconfigure(styr.lbl,state="disabled")
          tkconfigure(rho.lbl,state='disabled')
          tkconfigure(aR.lbl,state='disabled')
          tkconfigure(tau.lbl,state="disabled")

        } else if (tclvalue(tkvars$Ltype)=="rT"){
          tkconfigure(aCI.edit,state='disabled')
          tkconfigure(rho.edit,state='normal')
          tkconfigure(aL.edit,state='disabled')
          tkconfigure(styr.edit,state='disabled')
          tkconfigure(aR.edit,state='normal')
          tkconfigure(tau.edit,state="normal")
          tkconfigure(aL.lbl,state="disabled")
          tkconfigure(styr.lbl,state="disabled")
          tkconfigure(rho.lbl,state='normal')
          tkconfigure(aR.lbl,state='normal')
          tkconfigure(tau.lbl,state="normal")
        }
      })
      tkbind(myCIradio,"<Button-1>",function(){
        if (tclvalue(tkvars$option)=="M") return(F)
          tkconfigure(aCI.edit,state='normal')
          tkconfigure(rho.edit,state='disabled')
          tkconfigure(aL.edit,state='disabled')
          tkconfigure(styr.edit,state='disabled')
          tkconfigure(aR.edit,state='disabled')
          tkconfigure(tau.edit,state="disabled")
          tkconfigure(aL.lbl,state="disabled")
          tkconfigure(styr.lbl,state="disabled")
          tkconfigure(rho.lbl,state='disabled')
          tkconfigure(aR.lbl,state='disabled')
          tkconfigure(tau.lbl,state="disabled")
      })
      tkbind(myLradio,"<Button-1>", function(){
        if (tclvalue(tkvars$option)=="M") return(F)
          tkconfigure(aCI.edit,state='disabled')
          tkconfigure(rho.edit,state='disabled')
          tkconfigure(aL.edit,state='normal')
          tkconfigure(styr.edit,state='normal')
          tkconfigure(aR.edit,state='disabled')
          tkconfigure(rho.lbl,state='disabled')
          tkconfigure(aR.lbl,state='disabled')
          tkconfigure(tau.edit,state="normal")
          tkconfigure(aL.lbl,state="normal")
          tkconfigure(styr.lbl,state="normal")
          tkconfigure(tau.lbl,state="normal")
      })
      tkbind(myRradio,"<Button-1>",function(){
        if (tclvalue(tkvars$option)=="M") return(F)
          tkconfigure(aCI.edit,state='disabled')
          tkconfigure(rho.edit,state='normal')
          tkconfigure(aL.edit,state='disabled')
          tkconfigure(styr.edit,state='disabled')
          tkconfigure(aR.edit,state='normal')
          tkconfigure(tau.edit,state="normal")
          tkconfigure(rho.lbl,state='normal')
          tkconfigure(aR.lbl,state="normal")
          tkconfigure(aL.lbl,state="disabled")
          tkconfigure(styr.lbl,state="disabled")
          tkconfigure(tau.lbl,state="normal")
      })
      #################################################################
      #################################################################
      #################################################################
      ## the custom editing of projection table:
      #1. error-checking on the data table...do not allow projection until past data is entered correctly! This statement applies to all the projection graphs...
      #2. lock past data table...applies to all
      #3. create an editable projection parameters table (if one doesn't already exist)
      #4. grid the projTable
      projData<<-tclArray()
      myProjFill(mydat)
      valCharproj<-function(S){
        actcol<-as.numeric(tclvalue(tkindex(projTable,"active","col")))
        actrow<-as.numeric(tclvalue(tkindex(projTable,"active","row")))
        if (length(grep("\n",S))>0){ #
          return(tcl("expr", FALSE))
        } else if (length(grep("\t", S)) > 0) {
          return(tcl("expr", FALSE))
        } else {
          projCellChk(S, actrow, actcol)
          return(tcl("expr", TRUE))# but other kinds of space are not--> error-checking if good value is added
        }
      }
      priorYrs<-as.numeric(tclvalue(tcl(yearTable,"index","end","row")))
      nyrTot<-try(as.numeric(tclvalue(tkvars$nyr)),silent=T)
      ptabFrame <<- tkframe(tablesFrame)
      projTable <<- tcltk2::tk2table(ptabFrame,
        rows=nyrTot-priorYrs+1,
        cols=6,
        selectmode="extended",
        variable=projData,
        titlerows="1",
        titlecols="1",
        background=colors()[4],
        resizeborders="none",
        multiline = 1,
        wrap = 1,
#        multiline=F,
        rowseparator='\n',
        colseparator='\t',
        validate=T,
        vcmd=function(S) valCharproj(S),
        yscrollcommand = function(...) tkset(yscr,...),
        selecttitle = 1
      )
      yscr <<- tcltk2::tk2scrollbar(ptabFrame, orient = "vertical", command = function(...) tkyview(projTable, ...))
      tcl(projTable,"tag","configure", "error",bg=colors()[652]) # cells with error have yellow background
      tcl(projTable,"tag","configure", "cellok",bg=colors()[4]) # light blue background on future parameter cells
      tcl(projTable,"tag","configure", "active",fg='black',relief='groove')
      if (tclvalue(tkvars$Ptype) %in% c("C","I")) tkconfigure(projTable, state='disabled')
      colwidths<-c(1,15,7,10,10,10)
      for(i in 1:length(colwidths)) {
        tcl(projTable, "width", i - 1, colwidths[i])
      }
      tkbind(projTable,"<Return>",function(){
        if(tclvalue(tkindex(projTable,"active","row"))==tclvalue(tkindex(projTable,"end","row"))){
          if (tclvalue(tkindex(projTable,"active","col"))!=tclvalue(tkindex(projTable,"end","col"))){
            tkevent.generate(projTable,"<KeyPress-Right>")
          }
        } else {
          tkevent.generate(projTable,"<KeyPress-Down>")
        }
      })
      tkbind(projTable,"<Tab>",function(){
        if(tclvalue(tkindex(projTable,"active","col"))!=tclvalue(tkindex(projTable,"end","col"))){
          tkevent.generate(projTable,"<Right>")
        } else {
          tkevent.generate(projTable,"<KeyPress-Down>")
        }
      })
      tkbind(projTable,"<1>",function(){ # mouse button is pressed...check for errors
        if (substr(tclvalue(tcl(projTable,"curvalue")),1,1)=='.') tcl(projTable,"curvalue",paste0(0,tclvalue(tcl(projTable,"curvalue"))))
      })
      tkevent.delete("<<Paste>>","<Control-v>") # hijack the normal activity of the ctrl + v paste command
      junk<-NA
      tkbind(projTable,"<Control-Key-v>", function(){ # paste but don't add rows or columns
      # after pasting data from the clipboard, read the prior table data into a buffer
      # if there are errors in the new data, give error message and replace
        junk<-readClipboard()
        nrows<-length(junk)
        extracols<-max(nchar(gsub('[^\t]','',junk)))
        if (extracols+as.numeric(tclvalue(tcl(projTable,"index","active","col"))) > 5) {
          tkmessageBox(message="Error in shape of data paste ",icon='error')
          return(F)
        }
        if (nrows+as.numeric(tclvalue(tcl(projTable,"index","active","row")))-1 > as.numeric(tclvalue(tkvars$nyr))-as.numeric(tclvalue(tcl(yearTable,"index","end","row")))){
          tkmessageBox(message="Error in shape of data paste ",icon='error')
          return(F)
        }
        tkevent.generate(projTable,"<<Paste>>")
        projTableChk()
      })
      movechars<-c('Up','Down','Left','Right','Tab','Return','Shift_L','Alt_L')
      tkbind(projTable,"<Key>",function(K){
        if (K %in% movechars) {
          if (substr(tclvalue(tcl(projTable,"curvalue")),1,1)=='.') tcl(projTable,"curvalue",paste0(0,tclvalue(tcl(projTable,"curvalue"))))
        }
      })
      if (tclvalue(tkvars$option)=="M" & tclvalue(tkvars$Mtype)=="P"){
        tkgrid(yearProjections.lbl,sticky='sw',pady=c(20,2))
        tkgrid(projTable,yscr)
#        tkgrid.configure(projTable)
        tkgrid.configure(yscr, sticky='nsw')
        tkgrid(ptabFrame, row = 3, padx = 15, sticky='sw')
      }
      # error-check on nyr:
      tkbind(nyr.edit,"<KeyRelease>", function(){
        val<-suppressWarnings(as.numeric(toR(tkvars$nyr)))
        if (!pintchk(val) || val < as.numeric(tclvalue(tcl(yearTable,"index","end","row")))){
          tkconfigure(nyr.edit,bg=colors()[652])
          return(F)
        }
        tkconfigure(nyr.edit,bg='white')
        tkconfigure(projTable, rows = toR(tkvars$nyr)-as.numeric(tclvalue(tcl(yearTable, 'cget', '-rows')))+2)
        if (tclvalue(tkvars$Mtype) == "P"){
          myProjFill()
          projTableChk()
        }
#        if (tclvalue(tkvars$Mtype)=="P") tkconfigure(projVradio,state='normal')
        return(val)
      })
      tkbind(crlev.edit,"<KeyRelease>", function() {
        bgcolr<-ifelse(pchk(suppressWarnings(as.numeric(toR(tkvars$crlev)))), 'white','yellow')
        tkconfigure(crlev.edit,bg=bgcolr)
      })
      tkbind(aL.edit,"<KeyRelease>", function() {
        bgcolr<-ifelse(pchk(suppressWarnings(as.numeric(toR(tkvars$aL)))), 'white','yellow')
        tkconfigure(aL.edit,bg=bgcolr)
      })
      tkbind(aCI.edit,"<KeyRelease>", function() {
        bgcolr<-ifelse(pchk(suppressWarnings(as.numeric(toR(tkvars$aCI)))), 'white','yellow')
        tkconfigure(aCI.edit,bg=bgcolr)
      })
      tkbind(aR.edit,"<KeyRelease>", function() {
        bgcolr<-ifelse(pchk(suppressWarnings(as.numeric(toR(tkvars$aR)))), 'white','yellow')
        tkconfigure(aR.edit,bg=bgcolr)
      })
      tkbind(g.edit,"<KeyRelease>", function() {
        gleg<<-pchk(suppressWarnings(as.numeric(toR(tkvars$g))))
        bgcolr<-ifelse(gleg, 'white','yellow')
        tkconfigure(g.edit,bg=bgcolr)
        if (!gleg) return(F)
        for (i in 1:as.numeric(tclvalue(tkindex(projTable, "end", "row")))) projData[[i,3]]<-tclvalue(tkvars$g)
        gok<- gleg && glwrleg && guprleg && toR(tkvars$gupr) > toR(tkvars$g) && toR(tkvars$glwr) < toR(tkvars$g)
        gok<<-gok
        if (gok){
          for (i in 1:as.numeric(tclvalue(tkindex(projTable, "end", "row")))) projData[[i,4]]<-tclvalue(tkvars$glwr)
          for (i in 1:as.numeric(tclvalue(tkindex(projTable, "end", "row")))) projData[[i,5]]<-tclvalue(tkvars$gupr)
          tkconfigure(glwr.edit, bg='white')
          tkconfigure(gupr.edit, bg='white')
        } else {
          bgcolr<-ifelse(glwrleg && toR(tkvars$glwr) < toR(tkvars$g), 'white','yellow')
          tkconfigure(glwr.edit, bg=bgcolr)
          bgcolr<-ifelse(guprleg && toR(tkvars$gupr) > toR(tkvars$g), 'white','yellow')
          tkconfigure(gupr.edit, bg=bgcolr)
        }
        return(T)
      })
      tkbind(glwr.edit,"<KeyRelease>",function(){
        glwrleg<<-pchk(suppressWarnings(as.numeric(toR(tkvars$glwr))))
        bgcolr<-ifelse(glwrleg, 'white','yellow')
        tkconfigure(glwr.edit,bg=bgcolr)
        if (!glwrleg) return(F)
        gok<- gleg && glwrleg && guprleg && (toR(tkvars$gupr) > toR(tkvars$g)) && (toR(tkvars$glwr) < toR(tkvars$g))
        gok<<-gok
        if (gok){
          for (i in 1:as.numeric(tclvalue(tkindex(projTable, "end", "row")))) projData[[i,4]]<-tclvalue(tkvars$glwr)
          for (i in 1:as.numeric(tclvalue(tkindex(projTable, "end", "row")))) projData[[i,5]]<-tclvalue(tkvars$gupr)
          tkconfigure(glwr.edit, bg='white')
          tkconfigure(gupr.edit, bg='white')
        } else {
          bgcolr<-ifelse(glwrleg && toR(tkvars$glwr) < toR(tkvars$g), 'white','yellow')
          tkconfigure(glwr.edit, bg=bgcolr)
        }
        return(T)
      })
      tkbind(gupr.edit,"<KeyRelease>",function(){
        guprleg<<-pchk(suppressWarnings(as.numeric(toR(tkvars$gupr))))
        bgcolr<-ifelse(guprleg, 'white','yellow')
        tkconfigure(gupr.edit,bg=bgcolr)
        if (!guprleg) return(F)
        gok<- gleg && glwrleg && guprleg && (toR(tkvars$gupr) > toR(tkvars$g)) && (toR(tkvars$glwr) < toR(tkvars$g))
        gok<<-gok
        if (gok){
          for (i in 1:as.numeric(tclvalue(tkindex(projTable, "end", "row")))) projData[[i,4]]<-tclvalue(tkvars$glwr)
          for (i in 1:as.numeric(tclvalue(tkindex(projTable, "end", "row")))) projData[[i,5]]<-tclvalue(tkvars$gupr)
          tkconfigure(glwr.edit, bg='white')
          tkconfigure(gupr.edit, bg='white')
        } else {
          bgcolr<-ifelse(guprleg && toR(tkvars$gupr) > toR(tkvars$g), 'white','yellow')
          tkconfigure(gupr.edit, bg=bgcolr)
        }
        return(T)
      })
      tkbind(projYr.edit, "<KeyRelease>", function(){
        pyr <- suppressWarnings(as.numeric(toR(tkvars$styr)))
        bgcolr<-ifelse(!pintchk(n) || n > as.numeric(tclvalue(tcl(yearTable,"index","end","row"))), 'yellow','white')
        tkconfigure(styr.edit,bg=bgcolr)
      })
      tkbind(styr.edit,"<KeyRelease>", function() {
        n<-suppressWarnings(as.numeric(toR(tkvars$styr)))
        bgcolr<-ifelse(!pintchk(n) || n > as.numeric(tclvalue(tcl(yearTable,"index","end","row"))), 'yellow','white')
        tkconfigure(styr.edit,bg=bgcolr)
      })
      tkbind(Tau.edit,"<KeyRelease>", function() {
        bgcolr<-ifelse(gt0chk(suppressWarnings(as.numeric(toR(tkvars$Tau)))), 'white','yellow')
        tkconfigure(Tau.edit,bg=bgcolr)
      })
      tkbind(tau.edit,"<KeyRelease>", function() {
        bgcolr<-ifelse(gt0chk(suppressWarnings(as.numeric(toR(tkvars$tau)))), 'white','yellow')
        tkconfigure(tau.edit,bg=bgcolr)
      })
      tkbind(rho.edit,"<KeyRelease>", function() {
        bgcolr<-ifelse(gt0chk(suppressWarnings(as.numeric(toR(tkvars$rho)))), 'white','yellow')
        tkconfigure(rho.edit,bg=bgcolr)
      })
      tkbind(prho.edit,"<KeyRelease>", function() {
        bgcolr<-ifelse(gt0chk(suppressWarnings(as.numeric(toR(tkvars$prho)))), 'white','yellow')
        if (bgcolr=='white') for (i in 1:as.numeric(tclvalue(tkindex(projTable, "end", "row")))) projData[[i,2]]<-tclvalue(tkvars$prho)
        tkconfigure(prho.edit,bg=bgcolr)
      })

      rown<-as.numeric(tclvalue(tcl(projTable,"index","end","row")))
      if (tclvalue(tkvars$option)=="M" & tclvalue(tkvars$Mtype)=="P" & tclvalue(tkvars$Ptype)=="V") gridProjTable(rown)
      tkwm.protocol(myModule, "WM_DELETE_WINDOW", function(){
        tkdestroy(myModule)
        mysave(.Rvar$myPrevious)
        tkwm.deiconify(.Rvar$EoA)
      })
      setstates = function(){
        if (tclvalue(tkvars$option)=="M"){
          tkconfigure(nyr.edit,state ='normal')
          tkconfigure(Tau.edit,state ='normal')
          tkconfigure(crlev.edit, state ='normal')
          tkconfigure(myCradio,state ='normal')
          tkconfigure(myTradio,state ='normal')
          tkconfigure(myPradio,state ='normal')

          tkconfigure(myCIradio,state='disabled')
          tkconfigure(myLradio,state='disabled')
          tkconfigure(myRradio,state='disabled')
          tkconfigure(tau.lbl,state='disabled')
          tkconfigure(tau.edit,state='disabled')
          tkconfigure(styr.lbl, state='disabled')
          tkconfigure(aL.lbl, state='disabled')
          tkconfigure(rho.lbl, state='disabled')
          tkconfigure(aR.lbl,  state='disabled')


          if (tclvalue(tkvars$Mtype) == "P"){
            tkconfigure(projIradio, state ='normal')
            tkconfigure(projCradio, state ='normal')
            tkconfigure(projVradio, state ='normal')
            tkconfigure(g.edit, state ='normal')
            tkconfigure(g.lbl, state='normal')
            tkconfigure(gci.lbl, state='normal')
            tkconfigure(glwr.edit, state ='normal')
            tkconfigure(gupr.edit, state ='normal')
            tkconfigure(prho.lbl, state='normal')
            tkconfigure(prho.edit, state='normal')
          }  else {
            tkconfigure(projIradio, state ='disabled')
            tkconfigure(projCradio, state ='disabled')
            tkconfigure(projVradio, state ='disabled')
            tkconfigure(g.edit, state ='disabled')
            tkconfigure(g.lbl, state='disabled')
            tkconfigure(gci.lbl, state='disabled')
            tkconfigure(glwr.edit, state ='disabled')
            tkconfigure(gupr.edit, state ='disabled')
            tkconfigure(prho.lbl, state='disabled')
            tkconfigure(prho.edit, state='disabled')
          }

          tkconfigure(aL.edit,  state = 'disabled')
          tkconfigure(styr.edit, state = 'disabled')
          tkconfigure(tau.edit, state ='disabled')
          tkconfigure(tau.lbl, state ='disabled')
          tkconfigure(aR.edit, state ='disabled')
          tkconfigure(rho.edit, state ='disabled')
          tkconfigure(aCI.edit, state ='disabled')
        } else if (tclvalue(tkvars$option)=="L"){
          tkconfigure(nyr.edit,state='disabled')
          tkconfigure(Tau.edit,state='disabled')
          tkconfigure(crlev.edit, state='disabled')
          tkconfigure(myCradio,state='disabled')
          tkconfigure(myTradio,state='disabled')
          tkconfigure(myPradio,state='disabled')
          tkconfigure(projIradio, state ='disabled')
          tkconfigure(projCradio, state ='disabled')
          tkconfigure(projVradio, state ='disabled')
          tkconfigure(g.edit, state ='disabled')
          tkconfigure(g.lbl, state='disabled')
          tkconfigure(gci.lbl, state='disabled')
          tkconfigure(glwr.edit, state ='disabled')
          tkconfigure(gupr.edit, state ='disabled')
          tkconfigure(prho.lbl, state='disabled')
          tkconfigure(prho.edit, state='disabled')
          tkconfigure(nyr.lbl, state='disabled')
          tkconfigure(Tau.lbl, state='disabled')
          tkconfigure(projection.lbl, state='disabled')

          tkconfigure(myCIradio,state = 'normal')
          tkconfigure(myLradio, state = 'normal')
          tkconfigure(myRradio, state = 'normal')
          tkconfigure(tau.edit, state = 'normal')
          tkconfigure(tau.lbl, state = 'normal')
          if (tclvalue(tkvars$Ltype) == "stT") {# short-term trigger
            tkconfigure(aCI.edit,state='disabled')
            tkconfigure(rho.edit,state='disabled')
            tkconfigure(aL.edit,state='normal')
            tkconfigure(styr.edit,state='normal')
            tkconfigure(aR.edit,state='disabled')
            tkconfigure(tau.lbl,state='normal')
            tkconfigure(tau.edit,state='normal')
            tkconfigure(styr.lbl,state='normal')
            tkconfigure(aL.lbl,state='normal')
            tkconfigure(rho.lbl,state='disabled')
            tkconfigure(aR.lbl,state='disabled')
          } else if (tclvalue(tkvars$Ltype) == "ciT"){ # confidence interval
            tkconfigure(aCI.edit,state='normal')
            tkconfigure(rho.edit,state='disabled')
            tkconfigure(aL.edit,state='disabled')
            tkconfigure(styr.edit,state='disabled')
            tkconfigure(aR.edit,state='disabled')
            tkconfigure(tau.lbl,state='disabled')
            tkconfigure(tau.edit,state='disabled')
            tkconfigure(styr.lbl,state='disabled')
            tkconfigure(aL.lbl,state='disabled')
            tkconfigure(rho.lbl,state='disabled')
            tkconfigure(aR.lbl,state='disabled')
          } else if (tclvalue(tkvars$Ltype) == "rT"){ # reversion trigger
            tkconfigure(aCI.edit,state='disabled')
            tkconfigure(rho.edit,state='normal')
            tkconfigure(aL.edit,state='disabled')
            tkconfigure(styr.edit,state='disabled')
            tkconfigure(aR.edit,state='normal')
            tkconfigure(tau.lbl,state='normal')
            tkconfigure(tau.edit,state='normal')
            tkconfigure(styr.lbl,state='disabled')
            tkconfigure(aL.lbl,state='disabled')
            tkconfigure(rho.lbl,state='normal')
            tkconfigure(aR.lbl,state='normal')
          }
        }
      }
      setstates()
      projTableChk()
    },
    MYreadparm = function() {
        filename <- tclvalue(tkgetOpenFile(filetypes = "{{R data files} {.rds}}",defaultextension = ".rds", initialfile = '.rds', initialdir = .Rvar$csvpath, title = "Read"))
        tmp<-unlist(strsplit(filename,'/')); pathname<-paste(tmp[-length(tmp)],collapse='/')
        if (filename == "") return(FALSE)
        parms <- tryCatch(readRDS(filename), error=function(){tkmessageBox(icon='error',message=paste0("Error: Unable to read file:\n\n",filename)); return(F)})
        .Rvar$csvpath<-pathname
        tkdestroy(myModule)
        initialize(parms)
        tkwm.title(myModule,paste0(tmp[length(tmp)], " - EoA, v",.Rvar$VER," - Multiple Year Module"))
      },
    myPdo = function(){
     # draw projection table [clear projection table if "Total mortality" or "Track past mortality" button is selected]
      tkevent.generate(yearTable,"<KeyPress-Right>")
      tkevent.generate(yearTable,"<KeyPress-Left>")
      if (!yearTableChk()){
        tkmessageBox(message = "Error in past monitoring and operations data. Cannot make projections from bad data.")
        tclvalue(tkvars$Mtype)<<-"T"
        tkconfigure(yearTable, state = 'normal')
        return(F)
      }
      val<-suppressWarnings(as.numeric(toR(tkvars$nyr)))
      if (!pintchk(val)){
        tkmessageBox(message = "Error. Invalid number of years in project (", val,").")
        tclvalue(tkvars$Mtype)<<-"T"
        tkconfigure(yearTable, state = 'normal')
        return(F)
      }
      if (val <  as.numeric(tclvalue(tcl(yearTable,"index","end","row")))){
        tkmessageBox(message = "Error. Total years in project < years monitored...no future to project.")
        tclvalue(tkvars$Mtype)<<-"T"
        tkconfigure(yearTable, state = 'normal')
        return(F)
      }
      if (val ==  as.numeric(tclvalue(tcl(yearTable,"index","end","row")))){
        tkmessageBox(message = "Error. Total years in project = years monitored...no future to project.")
        tclvalue(tkvars$Mtype)<<-"T"
        tkconfigure(yearTable, state = 'normal')
        return(F)
      }
      tkconfigure(yearTable, state='disabled')
      myProjFill()
      if (pintchk(as.numeric(toR(tkvars$nyr)))){
        gridProjTable(toR(tkvars$nyr)-toR(tcl(yearTable, 'cget', '-rows'))+1)
      }
      tkconfigure(projIradio, state ='normal')
      tkconfigure(projCradio, state ='normal')
      if (!nyrchk()) tkconfigure(projVradio,state="disabled") else  tkconfigure(projVradio, state ='normal')
      gstate<-ifelse(tclvalue(tkvars$Ptype)=="C",'normal','disabled')#tkvars$Ptype<-tclVar("I") # I, C, V:
      tkconfigure(g.edit, state =gstate)
      tkconfigure(g.lbl, state=gstate)
      tkconfigure(gci.lbl, state=gstate)
      tkconfigure(glwr.edit, state =gstate)
      tkconfigure(gupr.edit, state =gstate)
      tkconfigure(prho.lbl, state=gstate)
      tkconfigure(prho.edit, state=gstate)
      tkconfigure(projection.lbl,state='normal')
      tkconfigure(nyr.lbl,state='normal')
      tkconfigure(Tau.lbl,state='normal')
      tkconfigure(nyr.edit,state='normal')
      tkconfigure(Tau.edit,state='normal')
      return(T)
    },
    myProjFill = function(mydat=.Rvar$myPrevious){
      columnNames<-c("","Year", "\u03c1", "g\u0302", "g_lwr", "g_upr")    # default is to fill out the table with blanks
      for (i in 1:length(columnNames)) projData[[0,i-1]]<<-strsplit(columnNames[i]," ",fixed=T)[[1]]
      priorYrs<-as.numeric(tclvalue(tcl(yearTable,"index","end","row")))
      nyrTot<-try(as.numeric(tclvalue(tkvars$nyr)),silent=T)
      if (class(nyrTot)=="try-error" || is.na(nyrTot)){
        tkmessageBox(message="Error. Years in project undefined. Cannot calculate projections for indefinite times.",icon="error")
        return(F)
      }
      if (nyrTot <= priorYrs){
        if (tclvalue(tkvars$Ptype) == "P" || nyrTot<priorYrs){
          tkmessageBox(message="Error. No future years in project. Cannot calculate projections for years that do not exist.",icon="error")
          return(F)
        }
      }
      if (tclvalue(tkvars$Ptype) == "I"){
        if (yearCellChk(tclvalue(self$yearData[[priorYrs,2]]), priorYrs, 2)) {
          prho <-toR(yearData[[priorYrs,2]])
        } else {
          prho <- ''
        }
        if (yearCellChk(tclvalue(self$yearData[[priorYrs,4]]), priorYrs, 4) &&
            yearCellChk(tclvalue(self$yearData[[priorYrs,5]]), priorYrs, 5)){
            atmp<-toR(yearData[[priorYrs,4]]); btmp<-toR(yearData[[priorYrs,5]])
            ghat <- round(atmp/(atmp+btmp), 3)
            glwr <- round(qbeta(0.025, atmp, btmp),4)
            gupr <- round(qbeta(0.975, atmp, btmp),4)
        } else {
          ghat<-''
          glwr<-''
          gupr<-''
        }
      } else if (tclvalue(tkvars$Ptype) == "C") {
        ghat <- ifelse(gok, toR(tkvars$g), '')
        glwr <- ifelse(gok && glwrleg && toR(tkvars$glwr)<toR(tkvars$g), toR(tkvars$glwr), '')
        gupr <- ifelse(gok && guprleg && toR(tkvars$gupr)>toR(tkvars$g), toR(tkvars$gupr), '')
        prho <- ifelse(!is.na(suppressWarnings(as.numeric(tclvalue(tkvars$prho)))) && toR(tkvars$prho) >= 0, toR(tkvars$prho), '')
      }
      for (i in 1:(nyrTot-priorYrs)){# table fills out to nyears from the end of the prior data
        projData[[i,0]]<<- as.tclObj('',drop=T)
        projData[[i,1]]<<- i# toR(tkvars$projyears[[i-1]])#as.tclObj(ifelse(!is.na(.Rvar$myPrevious$projyears[i]), .Rvar$myPrevious$projyears[i], i),drop=T)
        if (tclvalue(tkvars$Ptype) %in% c("C", "I")){
          projData[[i,2]] <<- prho
          projData[[i,3]] <<- ghat
          projData[[i,4]] <<- glwr
          projData[[i,5]] <<- gupr
        } else if (tclvalue(tkvars$Ptype) == "V") {
          projData[[i,2]]<<- as.tclObj(ifelse(!is.null(mydat$projrho[i]) && !is.na(mydat$projrho[i]),mydat$projrho[i], ''),drop=T)
          projData[[i,3]]<<- as.tclObj(ifelse(!is.null(mydat$projrho[i]) && !is.na(mydat$projg[i]),mydat$projg[i], ''),drop=T)
          projData[[i,4]]<<- as.tclObj(ifelse(!is.null(mydat$projrho[i]) && !is.na(mydat$projglwr[i]),mydat$projglwr[i], ''),drop=T)
          projData[[i,5]]<<- as.tclObj(ifelse(!is.null(mydat$projrho[i]) && !is.na(mydat$projgupr[i]),mydat$projgupr[i], ''),drop=T)
#          projData[[i,2]]<<- yearData[[priorYrs,2]]#toR(tkvars$projrho[[i-1]])#as.tclObj(ifelse(!is.na(.Rvar$myPrevious$projrho[i]),.Rvar$myPrevious$projrho[i], ''),drop=T)
#          projData[[i,3]]<<- toR(tkvars$projg[[i-1]])#as.tclObj(ifelse(!is.na(.Rvar$myPrevious$projg[i]),.Rvar$myPrevious$projg[i], ''),drop=T)
#          projData[[i,4]]<<- toR(tkvars$projglwr[[i-1]])#as.tclObj(ifelse(!is.na(.Rvar$myPrevious$projglwr[i]),.Rvar$myPrevious$projglwr[i], ''),drop=T)
#          projData[[i,5]]<<- toR(tkvars$projgupr[[i-1]])#as.tclObj(ifelse(!is.na(.Rvar$myPrevious$projgupr[i]),.Rvar$myPrevious$projgupr[i], ''),drop=T)
        } else {
          tkmessageBox("You've got to be kidding!\nTragic error...total termination imminent")
          q(save = 'no')
        }
      }
    },
    gridProjTable = function(rown){
      tkconfigure(self$projTable, rows=rown + 1)
      tkgrid(yearProjections.lbl,sticky='nw', pady=c(20,3))
      tkgrid(projTable, yscr)
#      tkgrid.configure(projTable)
      tkgrid.configure(yscr, sticky='nsw')
      tkgrid(ptabFrame, row = 3, padx = 15, sticky='sw')
    },
    feedMY = function(){ # error-checking and loading data from tcl form into R list
     # check the data in the yearTable...must be correct
      errmsg<-"Error in data. Aborting calculations..."
      mydat<-list()
      # yearTable data:
      if (!yearTableChk()) {
        tkmessageBox(message=errmsg)
        return(F)
      }
      yrs<-as.numeric(tclvalue(tcl(yearTable,"index","end","row")))
      mydat$years  <- numeric(yrs)
      mydat$rel_wt <- numeric(yrs)
      mydat$X      <- numeric(yrs)
      mydat$Ba     <- numeric(yrs)
      mydat$Bb     <- numeric(yrs)
      for (i in 1:yrs){
        mydat$years[i] <- tclvalue(yearData[[i,1]])
        mydat$rel_wt[i]<- as.numeric(tclvalue(yearData[[i,2]]))
        mydat$X[i]     <- as.numeric(tclvalue(yearData[[i,3]]))
        mydat$Ba[i]    <- as.numeric(tclvalue(yearData[[i,4]]))
        mydat$Bb[i]    <- as.numeric(tclvalue(yearData[[i,5]]))
      }
    ### options
      mydat$option<-tclvalue(tkvars$option) # M, L
      # error-checking on M parameters: tkvars$aM, tkvars$nyr, tkvars$Tau
      if (mydat$option=="M"){
        mydat$Mtype<-tclvalue(tkvars$Mtype) # C, T, P
        val<-try(tclvalue(tkvars$crlev),silent=T)
        if (class(val)=="try-error") {
          tkconfigure(crlev.edit,bg=colors()[652])
          tkmessageBox(message=errmsg,icon='error')
          return(F)
        }
        val<-suppressWarnings(as.numeric(val))
        if (is.na(val)) {
          tkconfigure(crlev.edit,bg=colors()[652])
          tkmessageBox(message=errmsg,icon='error')
          return(F)
        }
        if (val <=0 | val >= 1){
          tkconfigure(crlev.edit,bg=colors()[652])
          tkmessageBox(message=errmsg,icon='error')
          return(F)
        }
        mydat$crlev<-val
        if (tclvalue(tkvars$Mtype) %in% c("T","P")){
          #check nyr
          val<-try(tclvalue(tkvars$nyr),silent=T)
          if (class(val)=="try-error"){
            tkconfigure(nyr.edit,bg=colors()[652])
            tkmessageBox(message=errmsg,icon='error')
            return(F)
          }
          val<-suppressWarnings(as.numeric(val))
          if (is.na(val)){
            tkconfigure(nyr.edit,bg=colors()[652])
            tkmessageBox(message=errmsg,icon='error')
            return(F)
          }
          if (val < as.numeric(tclvalue(tcl(yearTable,'index','end','row')))){
            tkconfigure(nyr.edit,bg=colors()[652])
            tkmessageBox(message=errmsg,icon='error')
            return(F)
          }

          if (tclvalue(tkvars$Mtype)=="P" && val == as.numeric(tclvalue(tcl(yearTable,'index','end','row')))){
            tkconfigure(nyr.edit,bg=colors()[652])
            tkmessageBox(message=errmsg,icon='error')
            return(F)
          }
          if (val != round(val)){
            tkconfigure(nyr.edit,bg=colors()[652])
            tkmessageBox(message=errmsg,icon='error')
            return(F)
          }
          mydat$nyr<-val
          #check Tau
          val<-try(tclvalue(tkvars$Tau),silent=T)
          if (class(val)=="try-error") {
            tkconfigure(Tau.edit,bg=colors()[652])
            tkmessageBox(message=errmsg,icon='error')
            return(F)
          }
          val<-suppressWarnings(as.numeric(val))
          if (is.na(val)){
            tkconfigure(Tau.edit,bg=colors()[652])
            tkmessageBox(message=errmsg,icon='error')
            return(F)
          }
          if (val <= 0){
            tkconfigure(Tau.edit,bg=colors()[652])
            tkmessageBox(message=errmsg,icon='error')
            return(F)
          }
          mydat$Tau<-val
      #   if projection and C, check g, glwr, gupr
          if (tclvalue(tkvars$Mtype) == "P"){
            mydat$Ptype<-tclvalue(tkvars$Ptype) # I, C, V
            if (tclvalue(tkvars$Ptype)=="C"){
            # check g
              val<-try(tclvalue(tkvars$g),silent=T)
              if (class(val)=="try-error") {
                tkconfigure(g.edit,bg=colors()[652])
                tkmessageBox(message=errmsg,icon='error')
                return(F)
              }
              val<-suppressWarnings(as.numeric(val))
              if (is.na(val)) {
                tkconfigure(g.edit,bg=colors()[652])
                tkmessageBox(message=errmsg,icon='error')
                return(F)
              }
              if (val <= 0 | val > 1){
                tkconfigure(g.edit,bg=colors()[652])
                tkmessageBox(message=errmsg,icon='error')
                return(F)
              }
              g<-val
              mydat$g<-val
            # check glwr
              val<-try(tclvalue(tkvars$glwr),silent=T)
              if (class(val)=="try-error"){
                tkconfigure(glwr.edit,bg=colors()[652])
                tkmessageBox(message=errmsg,icon='error')
                return(F)
              }
              val<-suppressWarnings(as.numeric(val))
              if (is.na(val)){
                tkconfigure(glwr.edit,bg=colors()[652])
                tkmessageBox(message=errmsg,icon='error')
                return(F)
              }
              if (val <= 0 | val >= 1){
                tkconfigure(glwr.edit,bg=colors()[652])
                tkmessageBox(message=errmsg,icon='error')
                return(F)
              }
              if (val >= g) {
                tkconfigure(glwr.edit,bg=colors()[652])
                tkmessageBox(message=errmsg,icon='error')
                return(F)
              }
              mydat$glwr<-val
            # check gupr
              val<-try(tclvalue(tkvars$gupr),silent=T)
              if (class(val)=="try-error"){
                tkconfigure(gupr.edit,bg=colors()[652])
                tkmessageBox(message=errmsg,icon='error')
                return(F)
              }
              val<-suppressWarnings(as.numeric(val))
              if (is.na(val)){
                tkconfigure(gupr.edit,bg=colors()[652])
                tkmessageBox(message=errmsg,icon='error')
                return(F)
              }
              if (val <= 0 | val >= 1){
                tkconfigure(gupr.edit,bg=colors()[652])
                tkmessageBox(message=errmsg,icon='error')
                return(F)
              }
              if (val < g){
                tkconfigure(gupr.edit,bg=colors()[652])
                tkmessageBox(message=errmsg,icon='error')
                return(F)
              }
              mydat$gupr<-val

              val<-try(tclvalue(tkvars$prho),silent=T)
              if (class(val)=="try-error"){
                tkconfigure(gupr.edit,bg=colors()[652])
                tkmessageBox(message=errmsg,icon='error')
                return(F)
              }
              val<-suppressWarnings(as.numeric(val))
              if (is.na(val)){
                tkconfigure(tkvars$prho,bg=colors()[652])
                tkmessageBox(message=errmsg,icon='error')
                return(F)
              }
              if (val < 0){
                tkconfigure(gupr.edit,bg=colors()[652])
                tkmessageBox(message=errmsg,icon='error')
                return(F)
              }
              mydat$prho<-val
            }
      #     if projection and V, check projTable
            if (tclvalue(tkvars$Ptype)=="V"){
              if (!projTableChk()){
                tkmessageBox(message=errmsg,icon='error')
                return(F)
              }
              nyrs<-as.numeric(tclvalue(tcl(projTable,"index","end","row")))
              mydat$projyears<-numeric(nyrs)
              mydat$projrho<-numeric(nyrs)
              mydat$projg<-numeric(nyrs)
              mydat$projglwr<-numeric(nyrs)
              mydat$projgupr<-numeric(nyrs)
              for (i in 1:nyrs){
                mydat$projyears[i]<-toR(projData[[i,1]])
                mydat$projrho[i]<-toR(projData[[i,2]])
                mydat$projg[i]<-toR(projData[[i,3]])
                mydat$projglwr[i]<-toR(projData[[i,4]])
                mydat$projgupr[i]<-toR(projData[[i,5]])
              }
            }
          }
        } #    if estimateL
      } else if (tclvalue(tkvars$option) == "L") {
          mydat$Ltype<-tclvalue(tkvars$Ltype) # ciT, stT, rT
      # if calculate CI only
        if (tclvalue(tkvars$Ltype)=='ciT'){
          val<-try(tclvalue(tkvars$aCI),silent=T)
          if (class(val)=="try-error"){
            tkconfigure(aCI.edit,bg=colors()[652])
            tkmessageBox(message=errmsg,icon='error')
            return(F)
          }
          val<-suppressWarnings(as.numeric(val))
          if (is.na(val)){
            tkconfigure(aCI.edit,bg=colors()[652])
            tkmessageBox(message=errmsg,icon='error')
            return(F)
          }
          if (val <= 0 | val >= 1){
            tkconfigure(aCI.edit,bg=colors()[652])
            tkmessageBox(message=errmsg,icon='error')
            return(F)
          }
          mydat$aCI<-val
        # check ci bounds
        } else if (tclvalue(tkvars$Ltype) %in% c('stT','rT')){
        # check tau
          val<-try(tclvalue(tkvars$tau),silent=T)
          if (class(val)=="try-error"){
            tkconfigure(tau.edit,bg=colors()[652])
            tkmessageBox(message=errmsg,icon='error')
            return(F)
          }
          val<-suppressWarnings(as.numeric(val))
          if (is.na(val)){
            tkconfigure(tau.edit,bg=colors()[652])
            tkmessageBox(message=errmsg,icon='error')
            return(F)
          }
          if (val <= 0){
            tkconfigure(tau.edit,bg=colors()[652])
            tkmessageBox(message=errmsg,icon='error')
            return(F)
          }
          mydat$tau<-val
          if (tclvalue(tkvars$Ltype) == 'stT'){
            # check yrs
            val<-try(tclvalue(tkvars$styr),silent=T)
            if (class(val)=="try-error"){
              tkconfigure(styr.edit,bg=colors()[652])
              tkmessageBox(message=errmsg,icon='error')
              return(F)
            }
            val<-suppressWarnings(as.numeric(val))
            if (is.na(val)){
              tkconfigure(styr.edit,bg=colors()[652])
              tkmessageBox(message=errmsg,icon='error')
              return(F)
            }
            if (val <= 0){
              tkconfigure(styr.edit,bg=colors()[652])
              tkmessageBox(message=errmsg,icon='error')
              return(F)
            }
            if (val != round(val)){
              tkconfigure(styr.edit,bg=colors()[652])
              tkmessageBox(message=errmsg,icon='error')
              return(F)
            }
            mydat$styr<-val
            # check aL
            val<-try(tclvalue(tkvars$aL),silent=T)
            if (class(val)=="try-error"){
              tkconfigure(aL.edit,bg=colors()[652])
              tkmessageBox(message=errmsg,icon='error')
              return(F)
            }
            val<-suppressWarnings(as.numeric(val))
            if (is.na(val)){
              tkconfigure(aL.edit,bg=colors()[652])
              tkmessageBox(message=errmsg,icon='error')
              return(F)
            }
            if (val <= 0 | val >= 1){
              tkconfigure(aL.edit,bg=colors()[652])
              tkmessageBox(message=errmsg,icon='error')
              return(F)
            }
            mydat$aL<-val
          } else if (tclvalue(tkvars$Ltype) == 'rT'){
            # check aR
            val<-try(tclvalue(tkvars$aR),silent=T)
            if (class(val)=="try-error"){
              tkconfigure(aR.edit,bg=colors()[652])
              tkmessageBox(message=errmsg,icon='error')
              return(F)
            }
            val<-suppressWarnings(as.numeric(val))
            if (is.na(val)){
              tkconfigure(aR.edit,bg=colors()[652])
              tkmessageBox(message=errmsg,icon='error')
              return(F)
            }
            if (val <= 0 | val >= 1){
              tkconfigure(aR.edit,bg=colors()[652])
              tkmessageBox(message=errmsg,icon='error')
              return(F)
            }
            mydat$aR<-val
            # check rho
            val<-try(tclvalue(tkvars$rho),silent=T)
            if (class(val)=="try-error"){
              tkconfigure(rho.edit,bg=colors()[652])
              tkmessageBox(message=errmsg,icon='error')
              return(F)
            }
            val<-suppressWarnings(as.numeric(val))
            if (is.na(val)){
              tkconfigure(rho.edit,bg=colors()[652])
              tkmessageBox(message=errmsg,icon='error')
              return(F)
            }
            if (val <= 0){
              tkconfigure(rho.edit,bg=colors()[652])
              tkmessageBox(message=errmsg,icon='error')
              return(F)
            }
            mydat$rho<-val
          }
        }
      }
      mydat
    },
    setT = function(rowm,colm) tcl(yearTable,"tag", "celltag", "cellok",as.tclObj(paste0(rowm,',', colm))),
    setF = function(rowm,colm) tcl(yearTable,"tag", "celltag", "error",as.tclObj(paste0(rowm,',', colm))),
    updateClassgLab = function(rowm,pBa,pBb){
      self$yearData[[rowm,6]]<<-signif(pBa/(pBa+pBb),4)
      self$yearData[[rowm,7]]<<-as.tclObj(paste0("[",signif(qbeta(0.025,pBa,pBb),3),", ",signif(qbeta(0.975,pBa,pBb),3),"]"),drop=T)
    },
    allokmy = function(){
        gok<<-T; gleg<<-T; glwrleg<<-T; guprleg<<-T
   },
   myReload = function(mydat){ # reload data from R list into tk (with error-check)
    ### check whether data are loadable (i.e., radio button values are all defined properly)
      if (!(mydat$option %in% c("M", "L")) ||
          !(mydat$Mtype %in% c("C", "T", "P")) ||
          !(mydat$Ltype %in% c("ciT", "stT", "rT")) ||
          !(mydat$Ptype %in% c("I", "C", "V"))
      ) {tkmessageBox(message="Error in data. Cannot load."); return(F)  }
    ### if data are loadable, then load
      ## radio buttons
      tclvalue(tkvars$option)<-mydat$option
      tclvalue(tkvars$Mtype)<-mydat$Mtype
      tclvalue(tkvars$Ltype)<-mydat$Ltype
      tclvalue(tkvars$Ptype)<-mydat$Ptype
      # number of years of past data
      nyear<-max(c(length(mydat$years),
                     length(mydat$X),
                     length(mydat$Ba),
                     length(mydat$Bb),
                     length(mydat$rel_wt)))
      ## edit boxes:
      crlev<-mydat$crlev[1]
      if (length(crlev) == 0){
        tclvalue(tkvars$crlev)<-''
        bgcolr<-'yellow'
      } else if (is.na(crlev)) {
        tclvalue(tkvars$crlev)<-"NA"
        bgcolr<-'yellow'
      } else {
        tclvalue(tkvars$crlev)<-crlev
        bgcolr<-ifelse(pchk(suppressWarnings(as.numeric(crlev))),"white","yellow")
      }
      tkconfigure(crlev.edit,bg=bgcolr)
      aR<-mydat$aR[1]
      if (length(aR) == 0){
        tclvalue(tkvars$aR)<-''
        bgcolr<-'yellow'
      } else if (is.na(aR)) {
        tclvalue(tkvars$aR)<-"NA"
        bgcolr<-'yellow'
      } else {
        tclvalue(tkvars$aR)<-aR
        bgcolr<-ifelse(pchk(suppressWarnings(as.numeric(aR))),"white","yellow")
      }
      tkconfigure(aR.edit,bg=bgcolr)
      aCI<-mydat$aCI[1]
      if (length(aCI) == 0){
        tclvalue(tkvars$aCI)<-''
        bgcolr<-'yellow'
      } else if (is.na(aCI)) {
        tclvalue(tkvars$aCI)<-"NA"
        bgcolr<-'yellow'
      } else {
        tclvalue(tkvars$aCI)<-aCI
        bgcolr<-ifelse(pchk(suppressWarnings(as.numeric(aCI))),"white","yellow")
      }
      tkconfigure(aCI.edit,bg=bgcolr)
      aL<-mydat$aL[1]
      if (length(aL) == 0){
        tclvalue(tkvars$aL)<-''
        bgcolr<-'yellow'
      } else if (is.na(aL)) {
        tclvalue(tkvars$aL)<-"NA"
        bgcolr<-'yellow'
      } else {
        tclvalue(tkvars$aL)<-aL
        bgcolr<-ifelse(pchk(suppressWarnings(as.numeric(aL))),"white","yellow")
      }
      tkconfigure(aL.edit,bg=bgcolr)
      nyr<-mydat$nyr[1]
      if (length(nyr) == 0){
        tclvalue(tkvars$nyr)<-''
        bgcolr<-'yellow'
      } else if (is.na(nyr)) {
        tclvalue(tkvars$nyr)<-"NA"
        bgcolr<-'yellow'
      } else {
        tclvalue(tkvars$nyr)<-nyr
        nok<-pintchk(suppressWarnings(as.numeric(nyr)))
        if (nok) nyr<-as.numeric(nyr)
        if (!(nok && nyr >= nyear)) nok<F
        bgcolr<-ifelse(nok,"white","yellow")
      }
      tkconfigure(nyr.edit,bg=bgcolr)
      Tau<-mydat$Tau[1]
      if (length(Tau) == 0){
        tclvalue(tkvars$Tau)<-''
        bgcolr<-'yellow'
      } else if (is.na(Tau)) {
        tclvalue(tkvars$Tau)<-"NA"
        bgcolr<-'yellow'
      } else {
        tclvalue(tkvars$Tau)<-Tau
        bgcolr<-ifelse(gt0chk(suppressWarnings(as.numeric(Tau))),"white","yellow")
      }
      tkconfigure(Tau.edit, bg=bgcolr)
      tau<-mydat$tau[1]
      if (length(tau) == 0){
        tclvalue(tkvars$tau)<-''
        bgcolr<-'yellow'
      } else if (is.na(tau)) {
        tclvalue(tkvars$tau)<-"NA"
        bgcolr<-'yellow'
      } else {
        tclvalue(tkvars$tau)<-tau
        bgcolr<-ifelse(gt0chk(suppressWarnings(as.numeric(tau))),"white","yellow")
      }
      tkconfigure(tau.edit, bg=bgcolr)
      styr<-mydat$styr[1]
      if (length(styr) == 0){
        tclvalue(tkvars$styr)<-''
        bgcolr<-'yellow'
      } else if (is.na(styr)) {
        tclvalue(tkvars$styr)<-"NA"
        bgcolr<-'yellow'
      } else {
        tclvalue(tkvars$styr)<-styr
        bgcolr<-ifelse(pintchk(suppressWarnings(as.numeric(styr))),"white","yellow")
      }
      tkconfigure(styr.edit, bg=bgcolr)
      rho<-mydat$rho[1]
      if (length(rho) == 0){
        tclvalue(tkvars$rho)<-''
        bgcolr<-'yellow'
      } else if (is.na(rho)) {
        tclvalue(tkvars$rho)<-"NA"
        bgcolr<-'yellow'
      } else {
        tclvalue(tkvars$rho)<-rho
        bgcolr<-ifelse(gt0chk(suppressWarnings(as.numeric(rho))),"white","yellow")
      }
      tkconfigure(rho.edit, bg=bgcolr)
      prho<-mydat$prho[1]
      if (length(prho) == 0){
        tclvalue(tkvars$prho)<-''
        bgcolr<-'yellow'
      } else if (is.na(prho)) {
        tclvalue(tkvars$prho)<-"NA"
        bgcolr<-'yellow'
      } else {
        tclvalue(tkvars$prho)<-prho
        bgcolr<-ifelse(gt0chk(suppressWarnings(as.numeric(prho))),"white","yellow")
      }
      tkconfigure(prho.edit, bg=bgcolr)
      allokmy()
      g<-mydat$g[1]
      if (length(g) == 0){
        tclvalue(tkvars$g)<-''
        gleg<<-F
        gok<<-F
        bgcolr<-'yellow'
      } else if (is.na(g)) {
        tclvalue(tkvars$g)<-"NA"
        gleg<<-F
        gok<<-F
        bgcolr<-'yellow'
      } else {
        tclvalue(tkvars$g)<-g
        gleg<<-pchk(suppressWarnings(as.numeric(g)))
        if (!gleg) gok<<-F
        bgcolr<-ifelse(gleg,"white","yellow")
      }
      tkconfigure(g.edit, bg=bgcolr)
      glwr<-mydat$glwr[1]
      if (length(glwr) == 0){
        tclvalue(tkvars$glwr)<-''
        glwrleg<<-F
        gok<<-F
        bgcolr<-'yellow'
      } else if (is.na(glwr)) {
        tclvalue(tkvars$glwr)<-"NA"
        glwrleg<<-F
        gok<<-F
        bgcolr<-'yellow'
      } else {
        tclvalue(tkvars$glwr)<-glwr
        glwrleg<<-pchk(suppressWarnings(as.numeric(glwr)))
        if (!glwrleg) {
          gok<<-F
        } else {
          bgcolr<-ifelse(glwrleg && gleg && as.numeric(glwr) < as.numeric(g),"white","yellow")
        }
      }
      tkconfigure(glwr.edit, bg=bgcolr)
      gupr<-mydat$gupr[1]
      if (length(gupr) == 0){
        tclvalue(tkvars$gupr)<-''
        guprleg<<-F
        gok<<-F
        bgcolr<-'yellow'
      } else if (is.na(gupr)) {
        tclvalue(tkvars$gupr)<-"NA"
        guprleg<<-F
        gok<<-F
        bgcolr<-'yellow'
      } else {
        tclvalue(tkvars$gupr)<-gupr
        guprleg<<-pchk(suppressWarnings(as.numeric(gupr)))
        if (!guprleg) {
          gok<<-F
        } else {
          bgcolr<-ifelse(guprleg && gleg && as.numeric(gupr) > as.numeric(g),"white","yellow")
        }
      }
      tkconfigure(gupr.edit, bg=bgcolr)

      ## table of past data:
      # clear old data
      tkdelete(yearTable,"rows",1,as.numeric(tclvalue(tcl(yearTable,"index","end","row"))))
      # enter new data
      for (i in 1:nyear){
        yearData[[i,0]]<- as.tclObj('',drop=T)
        years<-suppressWarnings(as.numeric(mydat$years[i]))
        yearData[[i,1]]<- ifelse(length(years) == 1 && !is.na(years), years, "NA")
        rel_wt<-suppressWarnings(as.numeric(mydat$rel_wt[i]))
        yearData[[i,2]]<- ifelse(length(rel_wt) == 1 && !is.na(rel_wt), rel_wt, "NA")
        X<-suppressWarnings(as.numeric(mydat$X[i]))
        yearData[[i,3]]<- ifelse(length(X) ==1 && !is.na(X), X, "NA")
        Ba<-suppressWarnings(as.numeric(mydat$Ba[i]))
        Bb<-suppressWarnings(as.numeric(mydat$Bb[i]))
        Baok<-length(Ba)==1 && !is.na(Ba) && Ba>0
        Bbok<-length(Bb)==1 && !is.na(Bb) && Bb>0
        yearData[[i,4]]<- ifelse(Baok, signif(Ba,4), "NA")
        yearData[[i,5]]<- ifelse(Bbok, signif(Bb,4), "NA")
        yearData[[i,6]]<- ifelse(Baok && Bbok, round(Ba/(Ba+Bb),3), "NA")
        if (Baok && Bbok){
          yearData[[i,7]]<-as.tclObj(paste0("[", signif(qbeta(0.025, Ba, Bb),3), ", ", signif(qbeta(0.975,Ba, Bb),3),"]"),drop=T)
        } else {
          yearData[[i,7]]<-"NA"
        }
      }
      tkconfigure(yearTable, rows=nyear+1)
      yearTableChk()
      setstates()
    },
   yearCellChk = function(val, rowm, colm){
  ## function for checking value in a particular cell and tags the cell as error or normal
    if (length(val) == 0) return(F)
    if (colm==1){ # year name: can be anything (but not \n or \t...these values are impossible to enter by hand or copy-paste but may be possible by reading file)
      if (length(grep("[\n\t]",val))>0) {
        setF(rowm, colm)
        return(F)
      } else {
        setT(rowm, colm)
        return(T)
      }
    } else {
      val<-suppressWarnings(as.numeric(val))
      if (is.na(val)){
        setF(rowm, colm)
        return(F)
      }
    }
    if (is.na(val)){ # then the value entered is non-numeric
      setF(rowm,colm)
      if (colm %in% 4:5){
        self$yearData[[rowm,6]]<<-"NA"
        self$yearData[[rowm,7]]<<-"NA"
      }
      return(F)
    }
    if (val==''){
      setF(rowm,colm)
      if (colm %in% 4:5){
        self$yearData[[rowm,6]]<<-"NA"
        self$yearData[[rowm,7]]<<-"NA"
      }
      return(F)
    }
    if (colm==2) { # rel_wt: must be positive
      if (val<=0){
        setF(rowm, colm)
  #      yearData[[rowm,6]]<-"NA" # don't need to set g data to NA unless Ba or Bb are bad
  #      yearData[[rowm,7]]<-"NA"
        return(F)
      } else {
        setT(rowm, colm)  # tag the cell
        if (tclvalue(tcl(yearTable,"tag","includes","error",paste0(rowm,',',4)))=="0" &
            tclvalue(tcl(yearTable,"tag","includes","error",paste0(rowm,',',5)))=="0"){
          pBa<-as.numeric(tclvalue(yearData[[rowm,4]])); pBb<-as.numeric(tclvalue(yearData[[rowm,5]]))
          self$yearData[[rowm,6]]<<-signif(pBa/(pBa+pBb),4)
          self$yearData[[rowm,7]]<<-as.tclObj(paste0("[",signif(qbeta(0.025,pBa,pBb),3),", ",signif(qbeta(0.975,pBa,pBb),3),"]"),drop=T)
        }
        return(T)
      }
    }
    if (colm==3) { #X: must be non-negative integer
      if (val<0) {
        setF(rowm, colm) # negative
        return(F)
      }
      if (val-round(val)!=0){
        setF(rowm, colm) # not an integer
        return(F)
      }
  #    yearData[[rowm,colm]]<<-as.tclObj(val)
      setT(rowm, colm)
      if (tclvalue(tcl(yearTable,"tag","includes","error",paste0(rowm,',',4)))=="0" &
          tclvalue(tcl(yearTable,"tag","includes","error",paste0(rowm,',',5)))=="0"){
        pBa<-as.numeric(tclvalue(yearData[[rowm,4]])); pBb<-as.numeric(tclvalue(yearData[[rowm,5]]))
        self$yearData[[rowm,6]]<<-signif(pBa/(pBa+pBb),4)
        self$yearData[[rowm,7]]<<-as.tclObj(paste0("[",signif(qbeta(0.025,pBa,pBb),3),", ",signif(qbeta(0.975,pBa,pBb),3),"]"),drop=T)
      }
      return(T)
    }
    if (colm %in% 4:5) { #Ba and Bb: must positive
      if (val<=0){
        setF(rowm, colm)
        self$yearData[[rowm,6]]<<-"NA"
        self$yearData[[rowm,7]]<<-"NA"
        return(F)
      } else {
        setT(rowm, colm) # the new value is OK
        if (colm == 4){
          if (tclvalue(tcl(yearTable,"tag","includes","error",paste0(rowm,',',5)))=="0"){
            pBa<-val; pBb<-as.numeric(tclvalue(yearData[[rowm,5]]))
            updateClassgLab(rowm,pBa,pBb)
          }
        } else {
          if (tclvalue(tcl(yearTable,"tag","includes","error",paste0(rowm,',',4)))=="0"){
            pBa<-as.numeric(tclvalue(yearData[[rowm,4]])); pBb<-val
            updateClassgLab(rowm,pBa,pBb)
          }
        }
        return(T)
      }
    }
  },
  projTableChk = function(){
    prok<-T
    rown<-as.numeric(tclvalue(tkindex(projTable,"end","row")))
    for (i in 1:rown){
      for (j in 1:5){
        val<-try(tclvalue(projData[[i,j]]),silent=T)
        if (class(val)=='try-error') val<-''
        prok<-prok*projCellChk(val, i, j)
      }
    }
    prok<<-prok
    prok
  },
  yearTableChk = function(){
    #get number of rows and then check values, cell by cell
    tableOK<- T
    rown<-as.numeric(tclvalue(tcl(yearTable,"index","end","row")))
    for (i in 1:rown){
      for (j in 1:5){
        val<-try(tclvalue(yearData[[i,j]]),silent=T)
        if (class(val)=='try-error') val<-''
        tableOK<-tableOK*yearCellChk(val, i, j)
      }
    }
    tableOK
  },
  nyrchk = function(){
    val<-suppressWarnings(as.numeric(toR(tkvars$nyr)))
    if (length(val) != 1 || is.na(val) || nchar(val) == 0 || val<=0 || val != round(val)) return(F)
    if (val < as.numeric(tclvalue(tkindex(yearTable,"end","row")))) return(F)
    return(val)
  },
  projCellChk = function(val, rowm, colm){
    if (length(val) == 0) return(F)
    if (colm == 1){ # year name: can be anything (but not \n or \t...these values are impossible to enter by hand or copy-paste but may be possible by reading file)
      if (length(grep("[\n\t]",val))>0) {
        setFpt(rowm, colm)
        return(F)
      } else {
        setTpt(rowm, colm)
        return(T)
      }
    }
    val<-suppressWarnings(as.numeric(val))
    if (is.na(val)){ # then the value entered is non-numeric
      setFpt(rowm,colm)
      return(F)
    }
    if (val==''){
      setFpt(rowm,colm)
      return(F)
    }
    if (colm==2) { # rel_wt: must be positive
      if (val<=0){
        setFpt(rowm, colm)
        return(F)
      } else {
        setTpt(rowm, colm)  # tag the cell
        return(T)
      }
    }
    if (colm %in% 3:5) { #ghat, glwr, gupr: must be in (0, 1)... no live error-checking for ghat to be within CI or CI/mean combo compatible with beta distribution
      if (val <= 0 | val >= 1) {
        setFpt(rowm, colm) # not in interval
        return(F)
      }
      setTpt(rowm, colm)
      return(T)
    }
  },
  setTpt = function(rowm,colm) tcl(projTable,"tag", "celltag", "cellok",as.tclObj(paste0(rowm,',', colm))),
  setFpt = function(rowm,colm) tcl(projTable,"tag", "celltag", "error",as.tclObj(paste0(rowm,',', colm)))
 )
)
mysave<-function(my){
  save(my, file = paste0(.Rvar$datadir,"/myPrevious.Rdata"))
}


