symcForm <- R6::R6Class("symcForm",
  portable = FALSE,
  public = list(
    # tcl variables:
    symcModule = NA, cprFrame = NA, tkvars = list(), tksysc = list(), classData = NA, classTable = NA, tk.dwp = NA, tk.cname = NA,
    syCompClass = NA, start.edit = NA, nsearch.lbl = NA, nsearch.edit = NA,
    priorViewButton = NA, priorEditButton = NA, prior.lbox = NA, prior.lbl = NA, sides = NA,
    # R variables:
    prior_f = NA, prior_M = NA, myr = NA, custom.prior = NA,
    objective.prior = NA, arrcomponents = NA, persistence_distn = NA, ncols = NA,
    prdat = NA, geo1 = NA, arrfig.mini = NA, Ir = NA,

    # data checks...
    ArrConfirmed = T, spanok = T, prok = T, dateset = T, Udo = T, P1do = T, P2do = T,
    Xok = T, aok = T, SEnok = T, SExok = T, SExleg = T, kok = T, Iok = T, nsearchok = T,
    nsearchleg = T, startok = T, pdaok = T, pdbok = T, bminok = T, bminleg = T, bmaxok = T,
    bmaxleg = T,

    # initial values
    initialize = function(symcdat) {
      # create the required tcl variables from the data input
#      allsymcok() # assuming a valid data set (should be true if loading symcDefault or symcPrevious)
                  # error-checking happens later
      # store the initial parameters into .Rvar$symcdat
      # when arrival and prior parameters are updated, make changes in .Rvar$symcdat
      # when calculating or saving, write arrival and prior parameters from .Rvar$symcdat into .Rvar$symcPrevious
      .Rvar$symcdat<-symcdat
      for (i in 1:length(symcVar)){
        ind <- which(names(symcdat) == symcVar[i])
        tkvars[[symcVar[i]]] <<- tclVar(as.tclObj(ifelse(length(ind) > 0, symcdat[[ind]],""),drop=T))
      }
      narrvar<-length(symcArray)
      if (narrvar > 0){
        for (li in 1:narrvar){
          ind <- which(names(symcdat)==symcArray[li])
          tmp<-tclArray()
          if (length(ind) == 0){ # missing data
            tmp[[0]] <- as.tclObj('', drop=T)
          } else if (length(symcdat[[ind]]) == 1) { # scalar
            tmp[[0]] <- as.tclObj(ifelse(!is.na(symcdat[[ind]]), symcdat[[ind]], ''), drop=T)
          } else if (is.null(dim(symcdat[[ind]]))) { # vector
            for (j in 1:length(symcdat[[ind]])) tmp[[j-1]] <- as.tclObj(symcdat[[ind]][j],drop=T)
          } else { # matrix
            for (rowi in 1:dim(symcdat[[ind]])[1])
              for (coli in 1:dim(symcdat[[ind]])[2])
                tmp[[rowi-1, coli-1]] <- as.tclObj(symcdat[[ind]][rowi, coli], drop = T)
          }
          tkvars[[symcArray[li]]] <<- tmp
        }
      }
      nclass<-length(tkvars$class) # one for 'not searched'; new blank row is added automatically upon opening
      tmpclass<-list()

      classData <<- tclArray()
      columnNames<-c("","Class", "dwp", "X", "Ba", "Bb", "g\u0302", "95% CI")
      for (i in 1:length(columnNames)) classData[[0,i-1]]<-strsplit(columnNames[i]," ",fixed=T)[[1]]
      # unsearched area
      classData[[1,0]]<<-as.tclObj('',drop=T)
      classData[[1,1]]<<- "unsearched"
      classData[[1,2]]<<-tclvalue(tkvars$rel_wt[[0]])
      classData[[1,3]]<<- 0
      classData[[1,4]]<<- "---"
      classData[[1,5]]<<- "---"
      classData[[1,6]]<<- 0
      classData[[1,7]]<<- as.tclObj("[0, 0]",drop=T)
      for (i in 2:nclass){
        classData[[i,0]]<<-as.tclObj('',drop=T)
        classData[[i,1]]<<-tclvalue(tkvars$class[[i-1]])
        classData[[i,2]]<<-tclvalue(tkvars$rel_wt[[i-1]])
        classData[[i,3]]<<-tclvalue(tkvars$Xmc[[i-1]])
        .Rvar$Ba<-as.numeric(tclvalue(tkvars$Ba[[i-1]]))
        .Rvar$Bb<-as.numeric(tclvalue(tkvars$Bb[[i-1]]))
        classData[[i,4]]<<-signif(.Rvar$Ba,5)
        classData[[i,5]]<<-signif(.Rvar$Bb,5)
        classData[[i,6]]<<-signif(.Rvar$Ba/(.Rvar$Ba+.Rvar$Bb), 3)
        classData[[i,7]]<<-as.tclObj(paste0("[", round(qbeta(0.025,.Rvar$Ba,.Rvar$Bb),3),", ",round(qbeta(0.975, .Rvar$Ba,.Rvar$Bb),3),"]"),drop=T)
      }
      if (!basicMode){
        prior_f<<-symcdat$prior_f
        prior_M<<-symcdat$prior_M
        custom.prior<<-symcdat$custom.prior
        objective.prior<<-symcdat$objective.prior
        arrcomponents<<-symcdat$arrcomponents
      } else {
        prior_f<<-"Objective"
      }
      symcModule <<- tktoplevel()
      tkgrab.set(symcModule);  tkfocus(symcModule)
      tkwm.title(symcModule,paste0("EoA, v", .Rvar$VER, " - Multiple Class Module"))
      tkwm.resizable(symcModule,0,0)
      tkwm.deiconify(symcModule)
      tkwm.withdraw(.Rvar$EoA)
      symcFrame<-tkframe(symcModule)
      if (!basicMode) ArrConfirmed<-F
      ############
      symctopMenu <- tkmenu(symcModule); tkconfigure(symcModule,menu=symctopMenu)
      symchelpMenu <- tkmenu(symctopMenu,activebackground=colors()[125],activeforeground=colors()[109])
      symceditMenu <- tkmenu(symctopMenu,activebackground=colors()[125],activeforeground=colors()[109])
      tkadd(symchelpMenu, "command", label="Parameter Conversion", command=conversionsCalculator)
      tkadd(symchelpMenu, "command", label="About", command=function()tkmessageBox(title='Evidence of Absence (EoA)',message=about_text))
      tkadd(symceditMenu,"command",label="Restore defaults",command=function(){
        .Rvar$symcPrevious <- symcDefault
        suppressWarnings(file.remove(paste0(.Rvar$datadir,"/symcPrevious.Rdata")))
        tkdestroy(symcModule)
        initialize(symcDefault)
      })
      tkadd(symceditMenu,"command",label="Restore previous",command=function(){
        tkdestroy(symcModule)
        initialize(.Rvar$symcPrevious)
      })
      tkadd(symceditMenu,"command",label="Save to file (.Rdata)",command=function() {
        if (!exists('csvpath', env = .Rvar)) assign('csvpath', getwd(), env = .Rvar)
        filename <- tclvalue(tkgetSaveFile(filetypes = "{{R images} {.Rdata}}", defaultextension = ".Rdata", initialfile = '.Rdata', title = "Save", initialdir = .Rvar$csvpath))
        tmp<-unlist(strsplit(filename,'/')); pathname<-paste(tmp[-length(tmp)],collapse='/')
        if (nchar(pathname)>0) .Rvar$csvpath <- pathname
        if (filename == "") return(FALSE)
        symcProvisional<-list()
        for (nm in names(symcDefault)) symcProvisional[[nm]]<-toR(tkvars[[nm]])
        n_class <- toR(tkindex(classTable,"end","row"))
        symcProvisional$class<-numeric(n_class)
        symcProvisional$Xmc<-numeric(n_class)
        symcProvisional$rel_wt<-numeric(n_class)
        symcProvisional$Ba<-numeric(n_class)
        symcProvisional$Bb<-numeric(n_class)
        for (i in 1:n_class){
          symcProvisional$class[i]<-ifelse(is.na(toR(classData[[i,1]])), "NA", toR(classData[[i,1]]))
          symcProvisional$rel_wt[i]<-ifelse(is.na(toR(classData[[i,2]])), "NA", toR(classData[[i,2]]))
          symcProvisional$Xmc[i]<-ifelse(is.na(toR(classData[[i,3]])), "NA", toR(classData[[i,3]]))
          symcProvisional$Ba[i]<-ifelse(is.na(toR(classData[[i,4]])), "NA", toR(classData[[i,4]]))
          symcProvisional$Bb[i]<-ifelse(is.na(toR(classData[[i,5]])), "NA", toR(classData[[i,5]]))
        }
        save(symcProvisional, file = filename)
        tkwm.title(symcModule,paste0(tmp[length(tmp)], " - EoA, v", .Rvar$VER, " - Multiple Class Module"))
      })
      tkadd(symceditMenu,"command",label="Read from file (.Rdata)",command=function() {
        if (!exists('csvpath', env = .Rvar)) assign('csvpath', getwd(), env = .Rvar)
        filename <- tclvalue(tkgetOpenFile(filetypes = "{{R images} {.Rdata}}", defaultextension = ".Rdata", initialfile = '.Rdata', title = "Read", initialdir = .Rvar$csvpath))
        tmp<-unlist(strsplit(filename,'/')); pathname<-paste(tmp[-length(tmp)],collapse='/')
        if (nchar(pathname)>0) .Rvar$csvpath <- pathname
        if (filename == "") return(FALSE)
        symcdat<-get(load(filename))
        tkdestroy(symcModule)
        initialize(symcdat)
        symcFormChk()
        tkwm.title(symcModule,paste0(tmp[length(tmp)], " - EoA, v", .Rvar$VER, " - Multiple Class Module"))

#        for (i in 1:length(symcProvisional$Xmc)){
#          for (j in 1:5) {
#            classCellChk(tclvalue(classData[[i,j]]), i, j)
#          }
#        }
      })
      tkadd(symctopMenu, "cascade", label="Edit", menu=symceditMenu)
      tkadd(symctopMenu, "cascade", label="Help", menu=symchelpMenu)

      tkwm.protocol(symcModule,"WM_DELETE_WINDOW",function(){ # same as the following line but a different format?
        while(1){
          if (sink.number() == 0) break else sink()
        }
        tkdestroy(symcModule)
        tkwm.deiconify(.Rvar$EoA)
      })
      tkbind(symcModule,"<Destroy>",function() { # red X kills the window
        ### assumption is that no bad data have been previously saved to the symcProvisional dataframe
        ## [previous error-checking should assure that is true]
        tkevent.add("<<Paste>>","<Control-v>")
        while(1){
          if (sink.number() == 0) break else sink()
        }
        tkdestroy(symcModule)
        tkwm.deiconify(.Rvar$EoA)
      })
 ###########
      setT<-function(tabdat,rowm,colm) tcl(tabdat,"tag", "celltag", "cellok",as.tclObj(paste0(rowm,',', colm)))  #
      setF<-function(tabdat,rowm,colm) tcl(tabdat,"tag", "celltag", "error",as.tclObj(paste0(rowm,',', colm)))
      adat<-.Rvar$symcPrevious
      valChar<-function(S){
        actcol<-as.numeric(tclvalue(tkindex(classTable,"active","col")))
        actrow<-as.numeric(tclvalue(tkindex(classTable,"active","row")))
        if (length(grep("\n",S))>0){ #
          return(tcl("expr", FALSE))
        } else {
          classCellChk(S, actrow, actcol)
          return(tcl("expr", TRUE))# but other kinds of space are not--> error-checking if good value is added
        }
      }
      tableFrame<-tkframe(symcFrame)
      classTable <<- tcltk2::tk2table(tableFrame,
        rows=nclass+1,
        cols=8,
        selectmode="extended",
        variable=classData,
        titlerows="1",
        titlecols="1",
        background='white',
        resizeborders="none",
        multiline=F,
        rowseparator='\n',
        colseparator='\t',
        validate=T,
        vcmd=function(S) valChar(S),
        selecttitle = 1
      )

      #### step 1: create table for holding, displaying, editing data
      tcl(classTable, "tag", "configure", "error",bg=colors()[652]) # cells with error have yellow background
      tcl(classTable, "tag", "configure", "adderr", bg=colors()[400]) # cells that have valid format but incompatible with other values
      tcl(classTable, "tag", "configure", "cellok",bg='white')
      tcl(classTable, "tag", "configure", "active",fg='black',relief='groove')
      tcl(classTable, "tag", "configure", "readonly",state='disabled', bg = colors()[415])
      tcl(classTable, "tag", "configure", "ignore",state='disabled',bg=colors()[365],fg=colors()[42])
      tcl(classTable, "tag", "configure", "lock", state='disabled')
      tcl(classTable, "tag", "configure", "unlock", state='normal')
      tcl(classTable, "tag", "col", "readonly", 6)
      tcl(classTable, "tag", "col", "readonly", 7)

      tcl(classTable,"tag", "celltag", "readonly",as.tclObj(paste0(1,',',1)))
      tcl(classTable,"tag", "celltag", "readonly",as.tclObj(paste0(1,',',3)))
      tcl(classTable,"tag", "celltag", "readonly",as.tclObj(paste0(1,',',4)))
      tcl(classTable,"tag", "celltag", "readonly",as.tclObj(paste0(1,',',5)))
      if (tclvalue(tkvars$ICEoption)=="e"){
        tcl(classTable, "tag", "col", "lock", 4)
        tcl(classTable, "tag", "col", "lock", 5)
      }
      colwidths<-c(1,15,7,7,10,10,8,15)
      for(i in 1:8) {
        tcl(classTable, "width", i - 1, colwidths[i])
      }
      tkgrid(classTable,sticky='w', row=1)
      tkbind(classTable,"<Return>",function(){
        if (tclvalue(tkindex(classTable,"active","row"))==tclvalue(tkindex(classTable,"end","row"))){
          tkevent.generate(classTable,"<Control-KeyPress-a>")
          return(T)
        }
        if(tclvalue(tkindex(classTable,"active","row"))==tclvalue(tkindex(classTable,"end","row"))){
          if (tclvalue(tkindex(classTable,"active","col"))!=tclvalue(tkindex(classTable,"end","col"))){
            tkevent.generate(classTable,"<KeyPress-Right>")
          }
        } else {
          tkevent.generate(classTable,"<KeyPress-Down>")
        }
      })
      tkbind(classTable,"<Tab>",function(){
        if(tclvalue(tkindex(classTable,"active","col"))!=tclvalue(tkindex(classTable,"end","col"))){
          tkevent.generate(classTable,"<Right>")
        } else {
          tkevent.generate(classTable,"<KeyPress-Down>")
        }
      })

      tkbind(classTable,"<1>",function(){ # mouse button is pressed...check for errors
        if (substr(tclvalue(tcl(classTable,"curvalue")),1,1)=='.') tcl(classTable,"curvalue",paste0(0,tclvalue(tcl(classTable,"curvalue"))))
      })
      tkbind(classTable,"<Control-KeyPress-d>", function(){
      # delete current row (if there is only one row, just erase the data)
        actrow<-tclvalue(tkindex(classTable, "active", "row"))
        if (as.numeric(actrow)>1) tkdelete(classTable, "rows", actrow, 1)
      })
      tkevent.delete("<<Paste>>","<Control-v>") # hijack the normal activity of the ctrl + v paste command
      junk<-NA
      tkbind(classTable,"<Control-Key-v>", function(){
        if (tclvalue(tkvars$ICEoption)=="e") return(F)
      # after pasting data from the clipboard, read the prior table data into a buffer
      # if there are errors in the new data, give error message and replace
        if (as.numeric(tclvalue(tcl(classTable,"index","active","col")))>5)   return(F)
      #  junk<<-tkXselection.get(selection='CLIPBOARD') # grab dsta from the clipboard for preliminary analysis before pasting to the table
      #  junk<<-gsub('[^\t^\n]','',tclvalue(junk)) #remove everything but tabs and carriage returns for easy parsing of the data
        junk<-readClipboard()
        ncols<<-max(nchar(gsub('[^\t]','',junk)))
        if (ncols+as.numeric(tclvalue(tcl(classTable,"index","active","col")))>5) {
          tkmessageBox(message="Error: cannot paste over summary columns",icon='error')
          return(F)
        }
        nclass<-as.numeric(tclvalue(tkindex(classTable,"end","row")))
        createTmpClass(nclass)
        tkconfigure(classTable,rows=max(nclass+1,as.numeric(tclvalue(tcl(classTable,"index","active","row")))+length(junk)),cols=8)
        tkevent.generate(classTable,"<<Paste>>")
      })
      createTmpClass<-function(nclass){
        if (nclass>0){
          tmpclass$class<-numeric(nclass)
          tmpclass$Xmc<-numeric(nclass)
          tmpclass$Ba<-numeric(nclass)
          tmpclass$Bb<-numeric(nclass)
          tmpclass$rel_wt<-numeric(nclass)
          tmpclass$ghat<-numeric(nclass)
          tmpclass$cis<-numeric(nclass)
          tmpclass$data<-T
          for (i in 1:nclass){
            tmpclass$class[i]<-ifelse(is.null(classData[[i,1]]),'', tclvalue(classData[[i,1]]))
            tmpclass$rel_wt[i]<-ifelse(is.null(classData[[i,2]]),'', tclvalue(classData[[i,2]]))
            tmpclass$Xmc[i]<-ifelse(is.null(classData[[i,3]]),'', tclvalue(classData[[i,3]]))
            tmpclass$Ba[i]<-ifelse(is.null(classData[[i,4]]),'', tclvalue(classData[[i,4]]))
            tmpclass$Bb[i]<-ifelse(is.null(classData[[i,5]]),'', tclvalue(classData[[i,5]]))
            tmpclass$ghat[i]<-ifelse(is.null(classData[[i,6]]),'', tclvalue(classData[[i,6]]))
            tmpclass$cis[i]<-ifelse(is.null(classData[[i,7]]),'', tclvalue(classData[[i,7]]))
          }
        } else {
          tmpclass$data<-F
        }
        tmpclass$est<-tclvalue(tkvars$option)
        tmpclass$crlev<-as.numeric(tclvalue(tkvars$crlev))
      }
      restoreTmpClass<-function(){
        # restore values of table prior to changes
        if (!tmpclass$data)return(F)
        for (i in 1:length(tmpclass$class)){
          classData[[i,1]]<-as.tclObj(tmpclass$class[i])
          classData[[i,2]]<-as.tclObj(tmpclass$rel_wt[i])
          classData[[i,3]]<-as.tclObj(tmpclass$Xmc[i])
          classData[[i,4]]<-as.tclObj(tmpclass$Ba[i])
          classData[[i,5]]<-as.tclObj(tmpclass$Bb[i])
          classData[[i,6]]<-as.tclObj(tmpclass$ghat[i])
          classData[[i,7]]<-as.tclObj(tmpclass$cis[i])
        }
      }
      tkbind(classTable,"<Control-KeyPress-a>", function(){
      # add a row
        if (tclvalue(tkvars$ICEoption)=="e") return(F)
        tkinsert(classTable, "rows", "end", 1)
        rowind<-as.numeric(tclvalue(tcl(classTable,"index","end","row")))
        for (i in 1:5){
          tcl(classTable,"tag","celltag","error", as.tclObj(paste0(rowind,',', i)))
        }
        tcl(classTable,"activate", as.tclObj(paste0(rowind,',', 1)))
        classData[[rowind,7]]<- "NA"
        classData[[rowind,6]]<-"NA"
      })
      movechars<-c('Up','Down','Left','Right','Tab','Return','Shift_L','Alt_L')
      tkbind(classTable,"<Key>",function(K){
        if (K %in% movechars) {
          if (substr(tclvalue(tcl(classTable,"curvalue")),1,1)=='.') tcl(classTable,"curvalue",paste0(0,tclvalue(tcl(classTable,"curvalue"))))
        }
      })
      ######
      # prior frame is added if we need to estimate M
      # deactivated if M is not to be estimated
      ## prior frame
      parent<-symcFrame
      arrfFrame<-ttklabelframe(parent,text="Arrival Function")
      editArrivalsbutton<-tkbutton(arrfFrame,text="Edit")

        tksysc$firstsearch<<-tclVar(.Rvar$syscPrevious$firstsearch)
      if (!basicMode){
        if (.Rvar$symcPrevious$prior_f=="Custom") .Rvar$symcPrevious$prior_M<-.Rvar$symcPrevious$custom.prior
        if (.Rvar$symcPrevious$prior_f=="Objective") .Rvar$symcPrevious$prior_M<-.Rvar$symcPrevious$objective.prior
        priorFrame<-ttklabelframe(symcFrame,text="Prior Distribution")
        prior.lbox<<-tklistbox(priorFrame,height=3,width=12,selectmode="single",exportselection="0",bg='white');
        ttPriortmp1<-tkframe(priorFrame)
        priorEditButton<<-tkbutton(ttPriortmp1,text="Edit",width = 5, command=function(){
          prdat<<-.Rvar$symcPrevious
          cpr_form()
        })
        priorViewButton<<-tkbutton(ttPriortmp1,text="View",width=5,command=function() plotPrior(prior_M, prior_f))
        prnames<-c("Objective","Custom")
        for (i in 1:length(prnames)) tkinsert(prior.lbox,i-1,prnames[i])
        if (.Rvar$symcPrevious$prior_f=="Objective"){
          tmplbl<-"          ...           "
        } else if (.Rvar$symcPrevious$prior_f=="Custom") {
          ans<-try(suppressWarnings(min(which(cumsum(.Rvar$symcPrevious$custom.prior)>=0.95))-1),silent=T)
        }
        tkselection.set(prior.lbox,which(prnames==.Rvar$symcPrevious$prior_f)-1)
        tkactivate(prior.lbox,which(prnames==.Rvar$symcPrevious$prior_f)-1)
        prior.lbl<<-tklabel(priorFrame,text=tmplbl)
        if (.Rvar$symcPrevious$prior_f!="Objective") {
          tkconfigure(priorEditButton, state='normal')
        } else {
          tkconfigure(priorEditButton, state='disabled')
        }
        #priorFrame and arrfFrame are gridded together on the same row after arrfFrame has been initialized
        ## arrival frame has to be gridded

  #      parent$name<<-'symc'
        ######
        tksysc$firstsearch<<-tclVar(.Rvar$syscPrevious$firstsearch)
        tksysc$samtype<<-tclVar(.Rvar$syscPrevious$samtype)
        tksysc$span<<-tclVar(ifelse(.Rvar$syscPrevious$samtype=="Formula", .Rvar$syscPrevious$nsearch*.Rvar$syscPrevious$Isam, max(.Rvar$syscPrevious$days)))
        ######
        # this is a variable initialized to a dummy value to allow buildArr.r to graph a uniform arrival
        # if the uniform is selected, then all the classes must have the same start date and span

        tclvalue(tksysc$samtype)<-"Other" # needed for proper construction of the "sysc add" frame; set back to proper value once the arrivals frame is drawn
        arrUniformRadio<-tkradiobutton(arrfFrame,variable=tkvars$arrfun,text="Uniform",value="Uniform")
        arrCompoundRadio<-tkradiobutton(arrfFrame,variable=tkvars$arrfun,text="Compound",value="Compound")
        buttwid<-6
        viewArrivalsbutton<-tkbutton(arrfFrame,text="View",width=buttwid,command=function() {
          arrcomponents<-toR(tkvars$arrcomponents)
          plot_arrivals(tclvalue(tkvars$arrfun),
            toR(tkvars$arrstart), arrcomponents,
            toR(tkvars$lwr.u), toR(tkvars$upr.u), toR(tkvars$wt.u),
            toR(tkvars$lwr.p1),toR(tkvars$upr.p1),toR(tkvars$wt.p1), toR(tkvars$a.p1), toR(tkvars$b.p1),
            toR(tkvars$lwr.p2),toR(tkvars$upr.p2),toR(tkvars$wt.p2), toR(tkvars$a.p2), toR(tkvars$b.p2))
        })
        if (tclvalue(tkvars$arrfun)=="Uniform") tkconfigure(editArrivalsbutton, state="disabled")
        tkbind(editArrivalsbutton,"<ButtonPress-1>",function(X,Y){  # edit arrivals button just brings up a window to edit the arrival variables
          if (tclvalue(tcl(editArrivalsbutton,"cget","-state"))=="disabled") return(F)
          geo1<<-as.numeric(c(X,Y))
          if (tclvalue(tkvars$arrfun)=="Uniform") {
#            source('defineArrivalSeasonU.Rea')
          } else {
#           source('arrival_fun.Rea')
            arrival_fun()
            tkrplot::tkrreplot(arrfig.mini)
          }
        })
        plotarr.mini <- function() { # not currently scaled to endpoints
          par(mar=c(0,0,0,0))
          if (tclvalue(tkvars$arrfun)=="Uniform"){
            plot(0,0,type='n',xlim=c(0,364),ylim=c(0,1),yaxs='i',xaxs='i')
            x0<-as.numeric(as.Date(tclvalue(tksysc$firstsearch))-as.Date(paste0(format(as.Date(tclvalue(tksysc$firstsearch)),"%Y"),"-01-01")))
            if (tclvalue(tksysc$samtype)=="Formula"){
              span <- toR(tkvars$nsearch) * toR(tkvars$Isam)
            } else if (tclvalue(tksysc$samtype)=="Custom"){
              span <- max(toR(tkvars$days))
            } else {
              span<-toR(tksysc$span)
            }
            x1<-x0 + span
            lht<-.3
            polygon(c(x0,x0,x1,x1),c(0,lht, lht, 0))
          } else {
              arrcomponents<-toR(tkvars$arrcomponents)
            if (sum(arrcomponents)==0){
              tkmessageBox(message="No model componenents have been defined.\n Click 'Edit' to build a model or select 'Uniform' for simple model.")
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
        }
        arrfig.mini <<- tkrplot::tkrplot(arrfFrame, fun = plotarr.mini, hscale=.25, vscale=.08)
        # construct arrival frame
        tkbind(arrUniformRadio,"<Button-1>", function() {
          if (tclvalue(tkcget(arrCompoundRadio,"-state"))=="disabled") return(F)
          tclvalue(tkvars$arrfun)<-"Uniform"
          tkrplot::tkrreplot(arrfig.mini)
          tkconfigure(editArrivalsbutton,state='disabled')
          tkgrid(searchFrame)
#          tkgrid(firstsearch.lbl)
#          tkgrid(firstsearch.edit)
#          tkgrid(span.lbl)
#          tkgrid(span.edit)
        })
        tkbind(arrCompoundRadio,"<Button-1>", function() {
          if (tclvalue(tkcget(arrUniformRadio,"-state"))=="disabled") return(F)
          tclvalue(tkvars$arrfun)<-"Compound"
          tkrplot::tkrreplot(arrfig.mini)
          tkconfigure(editArrivalsbutton,state='normal')
          tkgrid.remove(searchFrame)
#          tkgrid.remove(firstsearch.lbl)
#          tkgrid.remove(firstsearch.edit)
#          tkgrid.remove(span.lbl)
#          tkgrid.remove(span.edit)
        })

        tkgrid(arrUniformRadio,sticky='w')
        tkgrid(arrCompoundRadio,sticky='w')
        tkgrid(arrfig.mini, column=1,row=0,rowspan=2,padx=10)
        tkgrid(editArrivalsbutton,column=2,row=0,rowspan=2, padx=10)
      }
        arrConfirmButton<-tkbutton(arrfFrame, text="Confirm", width = 7, command=function(){
          if (!basicMode){
            ArrConfirmed<<-T
            tkconfigure(arrConfirmButton,state="disabled")
            tkconfigure(editArrivalsbutton,state="disabled")
            tkconfigure(arrUniformRadio,state="disabled")
            tkconfigure(arrCompoundRadio,state="disabled")
            tkconfigure(firstsearch.lbl,state='disabled')
            tkconfigure(firstsearch.edit,state='disabled')
            tkconfigure(span.lbl,state='disabled')
            tkconfigure(span.edit,state='disabled')
            tkconfigure(hRadio,state='disabled')
          }
        ### clear table
#          }
        rowm<-as.numeric(tclvalue(tcl(classTable,"index","end","row")))
        if (rowm > 1) tkdelete(classTable,"rows",2,rowm)
        ncl<-as.numeric(tclvalue(tcl(classTable,"index","end","row")))
        if (ncl == 1){
          # write form data into .Rvar$syscPrevious
          .Rvar$symcdat$firstsearch <-toR(tksysc$firstsearch)
          .Rvar$symcdat$arrfun<-toR(tkvars$arrfun)
          .Rvar$symcdat$arrcomponents<-toR(tkvars$arrcomponents)
          .Rvar$symcdat$lwr.u<-toR(tkvars$lwr.u)
          .Rvar$symcdat$upr.u<-toR(tkvars$upr.u)
          .Rvar$symcdat$wt.u<-toR(tkvars$wt.u)
          .Rvar$symcdat$lwr.p1<-toR(tkvars$lwr.p1)
          .Rvar$symcdat$upr.p1<-toR(tkvars$upr.p1)
          .Rvar$symcdat$wt.p1<-toR(tkvars$wt.p1)
          .Rvar$symcdat$a.p1<-toR(tkvars$a.p1)
          .Rvar$symcdat$b.p1<-toR(tkvars$b.p1)
          .Rvar$symcdat$lwr.p2<-toR(tkvars$lwr.p2)
          .Rvar$symcdat$upr.p2<-toR(tkvars$upr.p2)
          .Rvar$symcdat$wt.p2<-toR(tkvars$wt.p2)
          .Rvar$symcdat$a.p2<-toR(tkvars$a.p2)
          .Rvar$symcdat$b.p2<-toR(tkvars$b.p2)
        }
        .Rvar$syscdat<-.Rvar$syscPrevious
          for (nm in names(.Rvar$symcdat)){
            if (nm %in% names(syscDefault)){
              .Rvar$syscdat[[nm]]<-.Rvar$symcdat[[nm]]
            }
          }
          classData[[1,2]]<-as.tclObj('',drop=T)
          tcl(classTable,"tag","celltag","error", as.tclObj(paste0(1, ',', 2)))
          suppressWarnings(try(file.remove(paste0(.Rvar$datadir,"/sysc.tmp")),silent=T))
        })
        tkconfigure(editArrivalsbutton, width = 7, pady = 0)
        tkgrid(editArrivalsbutton, row=0, column=2,sticky='n')
        tkgrid(arrConfirmButton, row=1, column=2)
#      }
        if (tclvalue(tkvars$ICEoption)=='e' & ! basicMode) tkgrid(arrfFrame, columnspan=2, sticky='nw')
        tclvalue(tksysc$samtype)<-.Rvar$syscPrevious$samtype
        ### search frame for entering date of first search and span (visible only when "Uniform" arrivals is selected)
        searchFrame<-tkframe(symcFrame)
        if (!basicMode){
          firstsearch.lbl<-tklabel(searchFrame,text="Date of first search\n  (yyyy-mm-dd)", anchor='w')
          firstsearch.edit<-tkentry(searchFrame,width=10, textvariable=tksysc$firstsearch, bg='white', justify = 'left')
          span.lbl<-tklabel(searchFrame,text="Span of monitoring season",anchor='e')
          span.edit<-tkentry(searchFrame,width=4, textvariable=tksysc$span, justify='right',bg='white')
          tkgrid(firstsearch.lbl, sticky='w')
          tkgrid(firstsearch.edit, padx=c(10,0), row=0, column = 1, sticky='e')
          tkgrid(span.lbl, columnspan=2, sticky='w')
          tkgrid(span.edit, column=1, row=1, sticky='e')
      }
      # slider for dates for arrival function
      # label for start and end dates of uniform arrivals (period of inference)
      ##########
      symcOptionsFrame<-ttklabelframe(symcFrame,text="   Options",padding=10)
      symcActionsFrame<-ttklabelframe(tableFrame,text="   Actions",padding=10)
      overall.lbl<-tklabel(symcOptionsFrame, text='Overall',fg=colors()[645])
      gRadio<-tkradiobutton(symcOptionsFrame,variable=tkvars$option,text="Estimate overall detection probability (g)",value="g")
      MRadio<-tkradiobutton(symcOptionsFrame,variable=tkvars$option,text="Estimate total mortality (M)",value="M")
      aFrame<-tkframe(symcOptionsFrame)
      crlev.lbl<-tklabel(aFrame,text="Credibility level (1 - \u03b1)");   crlev.edit<-tkentry(aFrame,width=4, textvariable=tkvars$crlev, justify='right',bg='white')
      sidesFrame<-tkframe(symcOptionsFrame)
      sides <<- tclVar(1)
      oneside <- tkradiobutton(sidesFrame, text = "One-sided CI (M*)", variable = sides, value = 1)
      twoside <- tkradiobutton(sidesFrame, text = "Two-sided CI", variable = sides, value = 2)
      tkgrid(oneside, sticky='w', padx=5)
      tkgrid(twoside, sticky='w', padx=5)
      tkgrid(crlev.lbl, padx=3, row = 0, column = 0)
      tkgrid(crlev.edit, padx = 3, row = 0, column = 1)
      individual.lbl<-tklabel(symcOptionsFrame,text='Individual classes',fg=colors()[645])
      eRadio<-tkradiobutton(symcOptionsFrame, variable=tkvars$ICEoption, text="Calculate g parameters from monitoring data", value="e", command=function(){
        if (!basicMode){
          tkgrid(arrfFrame, row = 1, column = 4, columnspan = 2, sticky='nw')
          tkgrid(searchFrame, row = 2, column = 4, sticky='w')
        } else {
          tkinvoke(arrConfirmButton)
        }
        tkconfigure(optionsCaveat.lbl,text=eCaveat)
        if (!ArrConfirmed & !basicMode){
          tkconfigure(arrUniformRadio, state='normal')
          tkconfigure(arrCompoundRadio, state='normal')
          if (tclvalue(tkvars$arrfun)=="Compound"){
            tkconfigure(editArrivalsbutton,state="normal")
          } else {
            tkconfigure(firstsearch.lbl,state='normal')
            tkconfigure(firstsearch.edit,state='normal')
            tkconfigure(span.lbl,state='normal')
            tkconfigure(span.edit,state='normal')
          }
          tkconfigure(arrConfirmButton, state='normal')
        }
        tcl(classTable, "tag", "col", "lock", 4)
        tcl(classTable, "tag", "col", "lock", 5)
        rowm<-as.numeric(tclvalue(tcl(classTable,"index","end","row")))
        if (rowm > 1) tkdelete(classTable,"rows",2,rowm)
        classData[[1,2]]<-as.tclObj('',drop=T)
        tcl(classTable,"tag","celltag","error", as.tclObj(paste0(1, ',', 2)))
      })
      hRadio<-tkradiobutton(symcOptionsFrame, variable=tkvars$ICEoption, text="Enter g parameters manually", value="h", command=function(){
        tkconfigure(optionsCaveat.lbl, text=hCaveat)
        if (!basicMode){
          tkgrid.remove(arrfFrame)
          tkgrid.remove(searchFrame)
        }
        tcl(classTable, "tag", "col", "unlock", 4)
        tcl(classTable, "tag", "col", "unlock", 5)
      })
#      hCaveat<-"Caveat:
#      Every class is assumed to have identical period of inference.
#      Either monitoring seasons are fully aligned, or
#      g's assume the same arrival function and apply to entire year."
#      eCaveat<-"Caveat:
#      All classes share the same arrival function.
#      Set arrival parameters and confirm before adding new classes.\n"

      hCaveat<-''
      eCaveat<-""
      optionsCaveat.lbl<-tklabel(symcOptionsFrame,text=ifelse(tclvalue(tkvars$ICEoption)=='e',eCaveat, hCaveat),fg=colors()[125],justify='left')
      tkconfigure(optionsCaveat.lbl,width=53)
      # construct options frame
      tkgrid(overall.lbl, sticky='w')
      tkgrid(MRadio, sticky='w', padx=c(8,0))
      tkgrid(aFrame, row=2, column = 0, padx=c(10,0))
      tkgrid(sidesFrame, row = 2, column = 1)
      tkgrid(gRadio,sticky='w', padx=8, sticky='w', columnspan=4)
      tkgrid(individual.lbl, sticky='w')
      tkgrid(eRadio,padx=8, sticky='w', columnspan=3)
      tkgrid(hRadio, padx = 8, sticky='w')
      tkgrid(optionsCaveat.lbl,columnspan=3)

      tkbind(gRadio,"<Button-1>", function() {
        tkconfigure(crlev.edit,state='disabled')
        tkconfigure(crlev.lbl,state='disabled')
        if (!basicMode) tkconfigure(prior.lbox,state='disabled',bg=colors()[247])
      })
      tkbind(MRadio,"<Button-1>", function() {
        tkconfigure(crlev.edit,state='normal')
        tkconfigure(crlev.lbl,state='normal')
        if (!basicMode) tkconfigure(prior.lbox,state='normal',bg='white')
      })
      tkbind(crlev.edit,"<KeyRelease>",function(){
        testcrlev<-suppressWarnings(as.numeric(tclvalue(tkvars$crlev)))
        if (length(testcrlev)==0){
          tkconfigure(crlev.edit, bg=colors()[652])
          return(F)
        } else if (is.na(testcrlev)) {
          tkconfigure(crlev.edit, bg=colors()[652])
          return(F)
        } else if (testcrlev<=0 | testcrlev>=1) {
          tkconfigure(crlev.edit, bg=colors()[652])
          return(F)
        } else {
          tkconfigure(crlev.edit, bg='white')
        }
        return(T)
      })
      if (!basicMode){
        tkbind(firstsearch.edit, "<KeyRelease>", function(){
          val<-tclvalue(tksysc$firstsearch)
          if (!startchk(val)) {
            tkconfigure(firstsearch.edit, bg=colors()[652])
            return(F)
          } else {
            tkconfigure(firstsearch.edit, bg='white')
            return(T)
          }
        })
        tkbind(span.edit,"<KeyRelease>", function(){
          val<-suppressWarnings(as.numeric(tclvalue(tksysc$span)))
          if (length(val) != 1 || is.na(val) || val<=0){
            tkconfigure(span.edit, bg=colors()[652])
            return(F)
          } else {
            tkconfigure(span.edit, bg='white')
          }
          return(T)
        })
        tkbind(prior.lbox,"<<ListboxSelect>>", function() {
          v<-tclvalue(tkget(prior.lbox,tkcurselection(prior.lbox)))
          if (v=='Objective') {
            tkconfigure(prior.lbl,text="          ...           ")
            tkconfigure(priorEditButton, state='disabled')
            tkconfigure(priorViewButton, state='normal')
            prok<<-T
            prior_f<<-"Objective"
            prior_M<<-objective.prior
          }
          if (v=='Custom') {
            pr95<-ifelse(sum(is.na(custom.prior)), NA, min(which(cumsum(custom.prior)>=0.95))-1)
            tkconfigure(prior.lbl,text=paste0("95th percentile = ", pr95))
            tkconfigure(priorEditButton, state='normal')
            buttstate<-ifelse(sum(is.na(custom.prior)),"disabled","normal")
            tkconfigure(priorViewButton,state=buttstate)
            if (buttstate=='disabled') prok<<-F else prok<<-T
            prior_f<<-"Custom"
            prior_M<<-custom.prior
          }
        })
      }
      tkgrid(symcFrame)
      # construct prior frame
      if (!basicMode){
        tkgrid(priorEditButton); tkgrid(priorViewButton)
        tkgrid(prior.lbox,ttPriortmp1,sticky='w')
        tkgrid(prior.lbl,columnspan=2)
      }
      # construct actions frame
      add.button<-tkbutton(symcActionsFrame,text="Add class",width=8,command=function(){
        if (toR(tkvars$ICEoption)=="e"){# "e" is estimation via entering parameters
          if (!basicMode && tclvalue(tkcget(arrConfirmButton,'-state'))=="normal") {
           # if confirmed, disable the arrConfirm and edit buttons
           if (tclvalue(tkmessageBox(message="Search classes share a common temporal arrival function.\nArrival function must be defined and confirmed before estimating g's for individual classes.\n\nConfirm arrival function now?", type="yesno"))=="yes"){
             tkinvoke(arrConfirmButton)
           } else {
            return(F)
           }
          } else if (basicMode) {
            #tkinvoke(arrConfirmButton)
          }
          syCompClass <<- syForm$new(.Rvar$syscdat, partial = T)
        } else {
          tkinsert(classTable, "rows", "end", 1)
          rowind<-as.numeric(tclvalue(tcl(classTable,"index","end","row")))
          for (i in 1:5){
            tcl(classTable,"tag","celltag","error", as.tclObj(paste0(rowind,',', i)))
          }
          tcl(classTable,"activate", as.tclObj(paste0(rowind,',', 1)))
          classData[[rowind,7]]<- "NA"
          classData[[rowind,6]]<-"NA"
        }
      })
      calc.button<-tkbutton(symcActionsFrame,text="Calculate",width=8,command=function() {
        # error-checking
        classData[[1,1]]<-classData[[1,1]]
        if (!symcFormChk()){
          if (as.numeric(tclvalue(tcl(classTable,"tag","includes","adderr","1,2")))){
            Sys.sleep(0.5)
            tcl(classTable, "tag", "col", "cellok", 2)
          }
          return(F)
        }
        # write form data to .Rvar$symcPrevious list
        nclass<-as.numeric(tclvalue(tkindex(classTable,"end","row")))
        .Rvar$symcdat$class<-numeric(nclass)
        .Rvar$symcdat$rel_wt<-numeric(nclass)
        .Rvar$symcdat$Xmc<-numeric(nclass)
        .Rvar$symcdat$Ba<-numeric(nclass)
        .Rvar$symcdat$Bb<-numeric(nclass)
        .Rvar$symcdat$class[1]<-"unsearched"
        .Rvar$symcdat$rel_wt[1]<-toR(classData[[1,2]])
        .Rvar$symcdat$Xmc[1]<-0
        .Rvar$symcdat$Ba[1]<-0.0001
        .Rvar$symcdat$Bb[1]<-1000
        for (i in 2:nclass){
          .Rvar$symcdat$class[i]<-toR(classData[[i,1]])
          .Rvar$symcdat$rel_wt[i]<-toR(classData[[i,2]])
          .Rvar$symcdat$Xmc[i]<-toR(classData[[i,3]])
          .Rvar$symcdat$Ba[i]<-toR(classData[[i,4]])
          .Rvar$symcdat$Bb[i]<-toR(classData[[i,5]])
        }
        .Rvar$symcdat$option<-tclvalue(tkvars$option)
        if (.Rvar$symcdat$option=="M") { # should these not be updated here but, rather, when the prior is updated?
          .Rvar$symcdat$crlev<-toR(tkvars$crlev)
          .Rvar$symcdat$prior_f<-prior_f
          .Rvar$symcdat$prior_M<-prior_M
          if (prior_f=="Objective"){
            .Rvar$symcPrevious$objective.prior<-prior_M
          } else if (prior_f=="Custom"){
            .Rvar$symcPrevious$custom.prior<-prior_M
          }
        } else {
          .Rvar$symcdat$crlev<-0.5
        }
        .Rvar$symcPrevious$ICEoption<-toR(tkvars$ICEoption)
        .Rvar$symcPrevious$arrfun<-tclvalue(tkvars$arrfun)
        if (.Rvar$symcPrevious$ICEoption=='e' && tclvalue(tkvars$arrfun)=="Compound"){
          .Rvar$symcPrevious$arrstart<-as.numeric(tclvalue(tkvars$arrstart))
          .Rvar$symcPrevious$arrcomponents<-toR(tkvars$arrcomponents)
          .Rvar$symcPrevious$lwr.u<-toR(tkvars$lwr.u)
          .Rvar$symcPrevious$upr.u<-toR(tkvars$lwr.u)
          .Rvar$symcPrevious$wt.u<-toR(tkvars$lwr.u)
          .Rvar$symcPrevious$lwr.p1<-toR(tkvars$lwr.u)
          .Rvar$symcPrevious$upr.p1<-toR(tkvars$lwr.u)
          .Rvar$symcPrevious$wt.p1<-toR(tkvars$lwr.u)
          .Rvar$symcPrevious$a.p1<-toR(tkvars$lwr.u)
          .Rvar$symcPrevious$b.p1<-toR(tkvars$lwr.u)
          .Rvar$symcPrevious$lwr.p2<-toR(tkvars$lwr.u)
          .Rvar$symcPrevious$upr.p2<-toR(tkvars$lwr.u)
          .Rvar$symcPrevious$wt.p2<-toR(tkvars$lwr.u)
          .Rvar$symcPrevious$a.p2<-toR(tkvars$lwr.u)
          .Rvar$symcPrevious$b.p2<-toR(tkvars$lwr.u)
        }
        # calculate
        symcCalc()
      })
      reset.button<-tkbutton(symcActionsFrame, text="Clear", width=8, command=function(){
        # unlock the arrivals function
        ArrConfirmed<<-F
        tkconfigure(hRadio,state='normal')
        if (!basicMode & tclvalue(tkvars$ICEoption)=="e"){
          tkconfigure(arrConfirmButton,state="normal")
          if (tclvalue(tkvars$arrfun)=="Compound") {
            tkconfigure(editArrivalsbutton,state="normal")
          } else {
            tkconfigure(firstsearch.lbl,state='normal')
            tkconfigure(firstsearch.edit,state='normal')
            tkconfigure(span.lbl,state='normal')
            tkconfigure(span.edit,state='normal')
          }
          tkconfigure(arrUniformRadio,state="normal")
          tkconfigure(arrCompoundRadio,state="normal")
        }
        # clear the table
        rowm<-as.numeric(tclvalue(tcl(classTable,"index","end","row")))
        if (rowm > 1) tkdelete(classTable,"rows",2,rowm)
        classData[[1,2]]<-as.tclObj('',drop=T)
        tcl(classTable,"tag","celltag","error", as.tclObj(paste0(1, ',', 2)))
      })
      symcClose<-tkbutton(symcActionsFrame, text = 'Close', command = function(){
        graphics.off()
        tkdestroy(symcModule)
        symcsave(.Rvar$symcPrevious, .Rvar$syscPrevious)
        tkwm.deiconify(.Rvar$EoA)
      })
      tkgrid(add.button, calc.button, reset.button, symcClose)
      tkgrid(symcActionsFrame, row = 0, column = 0, sticky='nw')
      #tkgrid(addbyhand.button)
      #tkgrid(comb.button)
      if (!basicMode & tclvalue(tkvars$ICEoption)=='h'){
        tkconfigure(arrConfirmButton,state="disabled")
        tkconfigure(editArrivalsbutton,state="disabled")
        tkconfigure(arrUniformRadio,state="disabled")
        tkconfigure(arrCompoundRadio,state="disabled")
        tkconfigure(firstsearch.lbl,state='disabled')
        tkconfigure(firstsearch.edit,state='disabled')
        tkconfigure(span.lbl,state='disabled')
        tkconfigure(span.edit,state='disabled')
      }
      # put form together
      tkgrid(symcOptionsFrame, columnspan = 3, rowspan = 3, row = 0, sticky='nw')
      if (!basicMode) {
        tkgrid(priorFrame, row = 0, column = 4, sticky='nw')
        if (tclvalue(tkvars$ICEoption) == 'e'){
          tkgrid(arrfFrame, row = 1, column = 4, columnspan = 2, sticky='nw')
          tkgrid(searchFrame, row = 2, column = 4, sticky='w')
        }
      }
      tkgrid(tableFrame, columnspan=6, sticky='nw', row=0, column = 4)
      if (tclvalue(tkvars$arrfun) == 'Compound'){
        tkgrid.remove(firstsearch.lbl)
        tkgrid.remove(firstsearch.edit)
        tkgrid.remove(span.lbl)
        tkgrid.remove(span.edit)
      }
      tkwm.protocol(symcModule, "WM_DELETE_WINDOW", function(){
        graphics.off()
        tkdestroy(symcModule)
        symcsave(.Rvar$symcPrevious, .Rvar$syscPrevious)
        tkwm.deiconify(.Rvar$EoA)
      })
    },
    updateClassgLab = function(rowm,pBa,pBb){
      classData[[rowm,6]]<-signif(pBa/(pBa+pBb),4)
      classData[[rowm,7]]<-as.tclObj(paste0("[",signif(qbeta(0.025,pBa,pBb),3),", ",signif(qbeta(0.975,pBa,pBb),3),"]"),drop=T)
    },
    classCellChk = function(val, rowm, colm){
      if (val==''){
        setF(classTable, rowm,colm)
        if (colm %in% 4:5){
          classData[[rowm,6]]<-"NA"
          classData[[rowm,7]]<-"NA"
        }
        return(F)
      }
      if (colm==1){ # class name: can be anything (but not \n or \t...these values are impossible to enter by hand or copy-paste but may be possible by reading file)
        if (length(grep("[\n\t]",val))>0) {
          setF(classTable, rowm, colm)
          return(F)
        } else {
          setT(classTable, rowm, colm)
          return(T)
        }
      } else {
        val<-suppressWarnings(as.numeric(val))
        if (length(val)==0 || is.na(val)){
          setF(classTable, rowm, colm)
          if (colm %in% 4:5){
            classData[[rowm,6]]<-"NA"
            classData[[rowm,7]]<-"NA"
          }
          return(F)
        }
      }
      if (rowm==1 & colm == 2){ # then val is for dwp unsearched. must be between 0 and 1; ghat and CIs should not be updated
        if (val < 0 | val >= 1){
          setF(classTable, rowm, colm)
          return(F)
        } else {
          setT(classTable, rowm, colm)  # tag the cell
          return(T)
        }
      }
      if (colm==2) { # rel_wt: must be positive and less than 1
        if (val<0 | val>1){
          setF(classTable, rowm, colm)
          return(F)
        } else {
          setT(classTable, rowm, colm)  # tag the cell
          return(T)
        }
      }
      if (colm==3) { #X: must be non-negative integer
        if (val<0 || val!= round(val)) {
          setF(classTable, rowm, colm) # negative or not an integer
          return(F)
        }
        setT(classTable, rowm, colm)
        return(T)
      }
      if (colm %in% 4:5) { #Ba and Bb: must positive...and g labels need to be updated
        if (val<=0){
          setF(classTable, rowm, colm)
          classData[[rowm,6]]<-"NA"
          classData[[rowm,7]]<-"NA"
          return(F)
        } else {
          setT(classTable, rowm, colm) # the new value is OK
          if (colm == 4){
            if (tclvalue(tcl(classTable,"tag","includes","error",paste0(rowm,',',5)))=="0"){
              pBa<-val; pBb<-as.numeric(tclvalue(classData[[rowm,5]]))
              updateClassgLab(rowm,pBa,pBb)
            }
          } else {
            if (tclvalue(tcl(classTable,"tag","includes","error",paste0(rowm,',',4)))=="0"){
              pBa<-as.numeric(tclvalue(classData[[rowm,4]])); pBb<-val
              updateClassgLab(rowm,pBa,pBb)
            }
          }
          return(T)
        }
      }
    },
    symcCalc = function(){
      nclass<-as.numeric(tclvalue(tkindex(classTable,"end","row")))
      coverage<-1-toR(classData[[1,2]])
      symcdat<-data.frame(array(dim=c(nclass,5)))
      names(symcdat)<-c('class','rel_wt','Xmc','Ba','Bb')
      symcdat$class[1]<-'unsearched'
      symcdat$rel_wt[1]<-1-coverage
      symcdat$Xmc[1]<-0
      symcdat$Ba[1]<-0.001
      symcdat$Bb[1]<-1000
      for (i in 2:nclass)
        for (j in 1:5)
          symcdat[i,j]<-toR(classData[[i,j]])

      ### The function calculates an overall average g along with estimated variance.
      ### Posterior is betabinomial on g and is calculated based on a uniform prior.
      x<-sum(symcdat$Xmc)
      pBa0<-symcdat$Ba; pBb0<-symcdat$Bb
      a<-symcdat$rel_wt/sum(symcdat$rel_wt)
      # p is assumed to be distributed as a beta RV with mean = pmean and variance = ((maxp-minp)/4)^2.
      # This determines beta parameters as follows:
      if (tclvalue(tkvars$option)=="M") alpha<-1-toR(tkvars$crlev)
      mu<-pBa0/(pBa0+pBb0)
      sig2<-pBa0*pBb0/((pBa0+pBb0)^2*(pBa0+pBb0+1))
      Eg<-sum(mu*a) # average, expected G is the weighted average of the observed g's
#      lambda<-max(sum(X), 0.5)/Eg
      Vg<-sum(sig2*a^2)#+(sum(a*(mu^2+sig2))-Eg^2)/lambda
      pBa<-Eg^2/Vg*(1-Eg)-Eg; pBb<-pBa*(1/Eg-1) # these are the two shape parameters for a beta distribution underlying the beta-binomial for X | M
      gmin<-qbeta(.025,shape1=pBa,shape2=pBb); gmax<-qbeta(.975,shape1=pBa,shape2=pBb)
      ### M needs to be calculated up to a reasonable maximum, but the final distribution is clipped at P(M > m) < 0.0001
      if (tclvalue(tkvars$option)=="M"){
        mmax<-ifelse (Vg < 0.00001, fmmax(x,Eg),fmmax.ab(x,pBa,pBb))
        M<-x:mmax
        if (Vg>0.00001){
        #      pBa<-Eg^2/Vg*(1-Eg)-Eg; pBb<-pBa*(1/Eg-1) # these are the two shape parameters for the beta distribution underlying the beta-binomial for X | M
          pXgM<-VGAM::dbetabinom.ab(x,size=M,shape1=pBa,shape2=pBb) # the probabilities of X for M = 0:mmax
        } else {
          pXgM<-dbinom(x,size=M,prob=Eg) # the probabilities of X for M = x:mmax
        }
        if (prior_f=="Objective"){
          pM<-diff((sqrt(c(x-1,x:mmax)+1)))
        } else {
          pM<-numeric(mmax+1)
          pM[1:length(prior_M)]<-prior_M
          pM<-pM[(x+1):length(pM)]
        }
        pMgX<-pXgM*pM; pMgX<-pMgX/sum(pMgX) # posterior distribution for M (ignoring M < x, which has probability = zero)
        M<-c(rep(0,x), M)
        pMgX<-c(rep(0,x),pMgX)
        cs<-cumsum(pMgX)
        pMem<-pMgX[1:min(which(cs>.999))]
        pMem<-pMem/sum(pMem)
        cs<-cumsum(pMem)
        if (tclvalue(sides)=='1') py<-c(1,1-cumsum(pMem[-length(pMem)])) else py <- pMem
        m<-1:length(py)-1
        cut.off<-min(which(py<=alpha))-1
        if (cut.off > 100) {
           wid<-1
        } else if (cut.off < 50){
          wid<-3
        } else {
           wid<-2
        }
        plot(m,py,
          xlab='m',
          ylab=ifelse(tclvalue(sides) == '1', paste0("Prob( M \u2265 m | X = ", x,")"),paste0("Prob( M = m | X = ", x,")")),
          cex.axis=.9*.Rvar$charSizeAdjust,
          cex.lab=.9*.Rvar$charSizeAdjust,
          xlim = c(ifelse(tclvalue(sides)=='1', 0, m[min(which(cumsum(py)>0.001))]),max(m))
        )
        col.use=rep("black", length(py))
        if (tclvalue(sides)=='1'){
           cut.off.cr<-min(which(cs>=1-alpha))
           col.use[m[1]:(cut.off.cr)]='red'
           for(i in 1:length(py)){
             lines(x=rep(m[i],2), y= c(0,py[i]), col=col.use[i], lwd=wid)
           }
        } else {
           lwrbnd<-min(which(cs > alpha/2))-1
           lwrArea<-ifelse(lwrbnd == 0, 0, cs[lwrbnd])
           uprbnd<-min(which(cs > 1-alpha + lwrArea))-1
           reds<-lwrbnd:uprbnd+1
           col.use[reds]=colors()[90]
          for(i in 1:length(py)){
            lines(x=rep(m[i],2), y= c(0,py[i]), col=col.use[i], lwd=wid)
          }
        }
        top<-par('usr')[4]/1.04
        if(.Rvar$platform == "windows") bringToTop()
        ## estimate annual rate
#        x<-sum(X)
        mmax<-ifelse (Vg< 0.00001, fmmax(x,Eg),fmmax.ab(x,pBa,pBb))
        ctprob<-0.0001 # A Jeffreys prior is used for the distribution of lambda. What is a reasonable maximum (and minimum) for lambda?

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
          meanL<-integrate(Vectorize(cpLgX,"L"),lower=0.0001,upper=Lmax,x=x,pBa,pBb,mmax=mmax)$val/sum(symcdat$rel_wt)#styr
          CI<-c(optimize(f=optqL.ab,interval=c(0.00001,Lmax),x=x,pBa,pBb,mmax=mmax,p=.025)$minimum, optimize(f=optqL.ab,interval=c(0.00001,Lmax),x=x,pBa,pBb,mmax=mmax,p=.975)$minimum)/sum(symcdat$rel_wt)#styr
        } else {
          cpLgX<-function(L,x,g,mmax) 1-pLgX(L,x,g,mmax)
          meanL<-integrate(Vectorize(cpLgX,"L"),lower=0.0001,upper=Lmax,x=x,g=Eg,mmax=mmax)$val/sum(symcdat$rel_wt)#styr
          CI<-c(optimize(f=optqL,interval=c(0.00001,Lmax),x=x,g=Eg,mmax=mmax,p=.025)$minimum, optimize(f=optqL,interval=c(0.00001,Lmax),x=x,g=Eg,mmax=mmax,p=.975)$minimum)/sum(symcdat$rel_wt)#styr
        }
        text(max(m),top,paste("Estimated annual fatality rate: \u03bb\u0303 = ",signif(meanL,3), "\n 95% CI = [ ",signif(CI[1],3),", ",signif(CI[2],4),"]",sep=""),adj=1,cex=.85*.Rvar$charSizeAdjust)
        text(max(m),0.9*top,paste("Overall detection probability: g = ",signif(Eg,4), "\n 95% CI = [ ",signif(gmin,4),", ",signif(gmax,4),"]",sep=""),adj=1,cex=.85*.Rvar$charSizeAdjust)
        if (tclvalue(sides)=='1') {
          text(max(m),.8*top,paste("M* = ",m[cut.off], ", i.e., P(M \u2264 ",  m[cut.off],") \u2265 ", (1-alpha)*100, "%",sep=''),adj=1,cex=.85*.Rvar$charSizeAdjust)
        } else {
          text(max(m),.8*top,paste0(100*(1-alpha),"% CI for M = [",m[min(reds)],", ",  m[max(reds)],"]  "),adj=1,cex=.85*.Rvar$charSizeAdjust)
        }
        mtext(side=3,paste('Credibility level (1 - \u03b1) = ',1-alpha),adj=0,family='serif', cex=.Rvar$charSizeAdjust)
        title(paste("Posterior Distribution of Total Fatalities over ",nclass-1," classes",sep=""))
        M<-1:length(pMem)-1
        CpMem<-1-cumsum(pMem)
      }
      ############ write results to a file:
      while(1){ if (sink.number()==0) break else sink() }
        sink(paste0(.Rvar$datadir,"/output"))
        cat("Summary statistics for multiple class estimate\n")
        cat(paste0(unlist(rep("=",80)),collapse=''))
        cat("\nInput: Detection probability, by search class\n")
        cat(paste0("  Search coverage = ",1-toR(classData[[1,2]]),"\n\n"))
        blanks<-paste(rep(" ",max(sapply(symcdat$class,FUN=nchar))-3),collapse='')
        cat(paste0("  Class", blanks, "   DWP     X    Ba     Bb   ghat    95% CI\n"))
        cat(sprintf("  %-*s %5.3g %5i %6s %6s %5s [%5s, %5s]\n",5 + nchar(blanks),
          symcdat$class[1],
          symcdat$rel_wt[1],
          0,
          "  --- ",
          "  --- ",
          "  0  ",
          "0",
          "0")
        )
        for (i in 2:length(symcdat$class)){
          cat(sprintf("  %-*s %5.3g %5i %6.4g %6.4g %5.3f [%5.3f, %5.3f]\n",5+nchar(blanks),
            symcdat$class[i],
            symcdat$rel_wt[i],
            symcdat$Xmc[i],
            symcdat$Ba[i],
            symcdat$Bb[i],
            symcdat$Ba[i]/(symcdat$Ba[i]+symcdat$Bb[i]),
            qbeta(0.025,symcdat$Ba[i],symcdat$Bb[i]),
            qbeta(0.975,symcdat$Ba[i],symcdat$Bb[i]))
          )
        }
        cat(paste0(unlist(rep("=",80)),collapse=''))
        Bab<-c(pBa, pBb)
        if (tclvalue(tkvars$arrfun)=="Uniform" && tclvalue(tkvars$ICEoption)=='e'){
          cat("\nResults for full site over the monitored period (assuming uniform arrivals)\n")
#          Bab<-.Rvar$syresult$Bab
#          pBa<-Bab[1]; pBb<-Bab[2]
        } else {
          cat(paste0("\nResults for full site\n"))
          cat(paste0(unlist(rep("_",80)),collapse=''))
          cat('\n')
#          Bab<-.Rvar$syresult$BabAnn
        }
        cat("\nDetection probability\n")
        cat(paste0('  Estimated g = ', round(Eg,3),", 95% CI = [", round(gmin,3),", ",round(gmax,3),"]\n",sep=''))
        cat(paste0("  Fitted beta distribution parameters for estimated g: Ba = ", round(pBa,4),", Bb = ",round(pBb,4),"\n"))
        cat("\nMortality\n")
        if (tclvalue(tkvars$option)=="M"){
          if (tclvalue(sides)=='1'){
            cat(paste0("  M* = ",  m[cut.off], " for credibility 1 - alpha = ", 1-alpha, ", i.e., P(M <= ", m[cut.off],") >= ",(1-alpha)*100),"%\n",sep='')
          } else {
            cat("   ",paste0(100*(1-alpha),"% CI for M = [",m[min(reds)],", ",  m[max(reds)],"]\n"))
          }
          cat(paste("  Estimated annual fatality rate: lambda = ",signif(meanL,3), ", 95% CI = [ ",signif(CI[1],3),", ",signif(CI[2],4),"]\n",sep=""))
        }
        ## results for each class:
        if (tclvalue(tkvars$ICEoption) == 'e'){
          cat("\nCarcass persistence:\n")
          arrivals.lbl<-paste0("   ",tclvalue(tkvars$arrfun)," arrivals")
          if (tclvalue(tkvars$arrfun)=="Uniform"){
            arrivals.lbl<-paste0(arrivals.lbl,'\n')
          } else {
            with (.Rvar$symcPrevious,{
              arrivals.lbl<-paste0(arrivals.lbl,' with ', sum(arrcomponents),' components:\n')  # this needs to be filled in...
              if (arrcomponents[1]) arrivals.lbl<-paste0(arrivals.lbl,'     uniform arrivals (', round(wt.u/sum(c(wt.u,wt.p1,wt.p2)*arrcomponents),3),') from ', format(as.Date(arrstart+lwr.u,origin="1970-01-01"),"%b %d"),
                ' to ', format(as.Date(arrstart+upr.u,origin="1970-01-01"),"%b %d"),'\n')
              if (arrcomponents[2]) arrivals.lbl<-paste0(arrivals.lbl, '     beta pulse (',round(wt.p1/sum(c(wt.u,wt.p1,wt.p2)*arrcomponents),3),')  ',
                format(as.Date(arrstart+lwr.p1,origin="1970-01-01"),"%b %d"), ' to ',
                format(as.Date(arrstart+upr.p1,origin="1970-01-01"),"%b %d"), ' with a = ',signif(a.p1,4),' and b = ',signif(b.p1,4),'\n')
              if (arrcomponents[3]) arrivals.lbl<-paste0(arrivals.lbl, '     beta pulse (',round(wt.p2/sum(c(wt.u,wt.p1,wt.p2)*arrcomponents),3),') from ',
                format(as.Date(arrstart+lwr.p2,origin="1970-01-01"),"%b %d"), ' to ',
                format(as.Date(arrstart+upr.p2,origin="1970-01-01"),"%b %d"), ' with a = ',signif(a.p2,4),' and b = ',signif(b.p2,4),'\n')
            })
          }
          cat(arrivals.lbl)
          cat(paste0(unlist(rep("_",80)),collapse=''))
          cat('\n')
          cat(paste0(unlist(rep("=",80)),collapse=''))
          conn<- suppressWarnings(try(file(paste0(.Rvar$datadir,"/sysc.tmp"),open='r'),silent=T))
          if (class(conn)[1]!="try-error"){
            cat("\nClass specific data:\n")
            classInput<-readLines(conn)
            for (i in 1:length(classInput)){
              if (i %% 4 == 1) cat(paste0(symcdat$class[i],"\n"))
              cat(classInput[i])
            }
          }
          try(close(conn),silent=T)
        }
        cat(paste0("\nTest of assumed relative weights (rho)\n"))
        rhosumry<-rhotest(mydat=list(X=symcdat$Xmc[-1], Ba=symcdat$Ba[-1], Bb=symcdat$Bb[-1], rel_wt=symcdat$rel_wt[-1], crlev = .Rvar$symcdat$crlev))
        blanks<-paste(rep(" ",max(sapply(symcdat$class,FUN=nchar))-3),collapse='')
        cat(paste0("  Class", blanks, "   Assumed   Fitted (95% CI)\n"))
        cat(sprintf("  %-*s    %5.3f      NA\n",5 + nchar(blanks),
          symcdat$class[1],
          1-sum(rhosumry$rho.ass)
        ))
        for (yi in 1:length(rhosumry$rho.ass)) cat(sprintf("  %-*s    %-5.3f    [%5.3f, %5.3f]\n", 5 + nchar(blanks), symcdat$class[yi + 1], rhosumry$rho.ass[yi], rhosumry$rhoqtls[1,yi], rhosumry$rhoqtls[5,yi]))
        cat(paste0("  p = ", round(rhosumry$pval, 5), ' for likelihood ratio test of H0: assumed rho = true rho\n'))
#        cat(paste0("    Quick test of relative bias: ", round(rhosumry$quicktest, 3),'\n'))
        if (tclvalue(tkvars$option)=="M"){
          cat("\nMortality rates (lambda) by class\n")
          cat(paste0("  Class", blanks, "   Median        IQR           95% CI\n"))
          cat(sprintf("  %-*s% 7s %12s  %12s", 5 + nchar(blanks),
            symcdat$class[1],
            " --- ",
            " --- ",
            "   --- "
          ))
          cat('\n')
          for (cli in 2:length(symcdat$class)){
            qtls<-postL.sumry(symcdat$Xmc[cli], symcdat$Ba[cli], symcdat$Bb[cli])
            cat(sprintf("  %-*s     %-7.2f [%5.2f, %5.2f]  [%5.2f, %5.2f]", 5 + nchar(blanks),
              symcdat$class[cli],
              qtls[3],
              qtls[2],
              qtls[4],
              qtls[1],
              qtls[5]
            ))
            cat('\n')
          }
          cat('\nPosterior distribution of M\n') # this belongs at the bottom b/c it takes up so much space
          cat('m     p(M = m) p(M > m)\n')
          CpMem<-1-cumsum(pMem)
          for (i in 1:length(pMem)) {
            cat(sprintf("%-5.0f  %6.4f   %6.4f\n",i-1, round(pMem[i],4), round(CpMem[i],4)))
          }
        }
        sink()
        if (tclvalue(tkvars$option) == 'g'){
          tit<-"Estimated detection probability (g) for multiple classes"
        } else {
          tit<-'Estimated mortality (M) & detection probability (g) for multiple classes'
        }
        file.show(paste0(.Rvar$datadir,"/output"),delete.file=T,title=tit)
    },
    symcFormChk = function(){
    # data outside the table:
    # if estimation option is "M" (rather than "g"), check alpha
      a<-1-toR(tkvars$crlev)
      if (tclvalue(tkvars$option) == "M" && (!val.numeric(a) || a<=0 || a>=1)){
        tkmessageBox(message="Error in data (alpha). Cannot calculate.")
        return(F)
      }
    # check firstsearch and span
#      if (tclvalue(tkvars$ICEoption) == "e" && tclvalue(tkvars$arrfun) == "Uniform"){
#        v<-tclvalue(tksysc$firstsearch)
#        if (!startchk(v)){ # first search is bad
#          tkmessageBox(message="Error in data (date of first search). Cannot calculate.")
#          return(F)
#        }
#        span<-toR(tksysc$span)
#        if (!val.numeric(span) || span <= 0){
#          tkmessageBox(message="Error in data (span). Cannot calculate.")
#          return(F)
#        }
#      }
    # check prior
      if (tclvalue(tkvars$option)=="M" & !basicMode){
        if (!exists("prior_M") || is.na(sum(prior_M)) || sum(prior_M<0)>0 || abs(sum(prior_M)-1)>0.00001) {
          tkmessageBox(message="Error in prior distribution. Cannot estimate M.")
          return(F)
        }
      }
    # cell check all cells except first row
      nclass<-as.numeric(tclvalue(tkindex(classTable,"end","row")))
      if (nclass<2){
        tkmessageBox(message="Missing data. Cannot calculate.")
        return(F)
      }
      if (is.null(classData[[1, 2]]) || !classCellChk(tclvalue(classData[[1, 2]]), 1, 2)){
        setF(classTable, 1, 2)
        tkmessageBox(message="Error in data. Cannot calculate.")
        return(F)
      }
      for (i in 2:nclass){
        for (j in 1:5){
          if (is.null(classData[[i, j]]) || !classCellChk(tclvalue(classData[[i, j]]), i, j)){
            setF(classTable, i, j)
            tkmessageBox(message="Error in data. Cannot calculate.")
            return(F)
          }
        }
      }
      a<-numeric(nclass)
      for (i in 1:nclass) a[i]<-toR(classData[[i,2]])
      if (is.na(sum(a)) || abs(sum(a)-round(sum(a)))>0.001){
        tcl(classTable, "tag", "col", "adderr", 2)
        tkmessageBox(message="Error. Weights (DWP) must sum to 1.")
        tcl("after",500,{})
        tcl(classTable, "tag", "col", "cellok", 2)
        return(F)
      }
      return(T)
    },
    arrival_fun = function(){
      .Rvar$arrProcess <- tktoplevel()
      tktitle(.Rvar$arrProcess) <- "Arrival process"
      tkgrab.set(.Rvar$arrProcess);  tkfocus(.Rvar$arrProcess)
      tkconfigure(.Rvar$arrProcess,width=1000)
      tkwm.resizable(.Rvar$arrProcess,0,0)
      tkwm.deiconify(.Rvar$arrProcess)
      arrcomponents<-toR(tkvars$arrcomponents)
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
      a.p1<-as.numeric(tclvalue(tkvars$a.p1))
      b.p1<-as.numeric(tclvalue(tkvars$b.p1))
      mu.p1<-tclVar(a.p1/(a.p1+b.p1))
      s2.p1<-tclVar(log(12*a.p1*b.p1/((a.p1+b.p1)^2*(a.p1+b.p1+1))))
      a.p2<-as.numeric(tclvalue(tkvars$a.p2)); b.p2<-as.numeric(tclvalue(tkvars$b.p2))
      mu.p2<-tclVar(a.p2/(a.p2+b.p2))
      s2.p2<-tclVar(log(12*a.p2*b.p2/((a.p2+b.p2)^2*(a.p2+b.p2+1))))

      arrparms<-tclArray()
      columnNames<-c("","Population", "wt", "start", "end", "Ba", "Bb")
      for (i in 1:length(columnNames)) arrparms[[0,i-1]]<-strsplit(columnNames[i]," ",fixed=T)[[1]]
      # unsearched area
      arrparms[[1,0]]<-as.tclObj('',drop=T)
      arrparms[[1,1]]<-"resident"; arrparms[[1,2]]<-tclvalue(tkvars$wt.u) ; arrparms[[1,3]]<-tclvalue(tkvars$lwr.u) ; arrparms[[1,4]]<-tclvalue(tkvars$upr.u); arrparms[[1,5]]<-1   ; arrparms[[1,6]]<-1
      arrparms[[2,1]]<-"pulse1";   arrparms[[2,2]]<-round(as.numeric(tclvalue(tkvars$wt.p1)),3); arrparms[[2,3]]<-round(as.numeric(tclvalue(tkvars$lwr.p1)),3); arrparms[[2,4]]<-round(as.numeric(tclvalue(tkvars$upr.p1)),3);
      arrparms[[2,5]]<-signif(a.p1,4); arrparms[[2,6]]<-signif(b.p1,4)
      arrparms[[3,1]]<-"pulse2";   arrparms[[3,2]]<-round(as.numeric(tclvalue(tkvars$wt.p2)),3); arrparms[[3,3]]<-round(as.numeric(tclvalue(tkvars$lwr.p2)),3); arrparms[[3,4]]<-round(as.numeric(tclvalue(tkvars$upr.p2)),3)
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
        fg='white',
        selecttitle = 1
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

      ats0 <- c(0,31,59,90,120,151,181,212,243,273,304,334)
      figlab<-format(as.Date(ats0,origin="1970-01-01"),"%b-%d")
      arrstart <- tclVar("0")
      ats <- (ats0-as.numeric(tclvalue(arrstart)))%%365
      startDateFrame<-tkframe(.Rvar$arrProcess)
      startlbl <- tklabel(.Rvar$arrProcess,
        text = paste0("Starting date: ",format(as.Date(as.numeric(tclvalue(arrstart)),
        origin="1970-01-01"),"%b-%d")),
        width=20,
        anchor="w")
      # A function that changes the label
      startDateOnSlide <- function(...) {
        tkconfigure(startlbl, text = paste0("Starting date: ",format(as.Date(as.numeric(tclvalue(arrstart)),origin="1970-01-01"),"%b-%d")))
        ats<<-(ats0-as.numeric(tclvalue(arrstart)))%%365
        tkrplot::tkrreplot(arrfig)
      }
      # Add the slider
      startDateSlider <- tkscale(startDateFrame, from = 0, to = 364, variable = arrstart, orient = "horizontal", length = scwid,
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
          mtext(side=2,line=.7,"                                  ", cex=.Rvar$charSizeAdjust)
          lines(-rep(as.numeric(tclvalue(arrstart)),2)%%365,par('usr')[3:4],lty=3)
          box()
          ats<-(ats0-as.numeric(tclvalue(arrstart)))%%365
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
        lwr.u<-as.numeric(tclvalue(tkvars$lwr.u));   upr.u<-as.numeric(tclvalue(tkvars$upr.u))
        lwr.p1<-as.numeric(tclvalue(tkvars$lwr.p1)); upr.p1<-as.numeric(tclvalue(tkvars$upr.p1))
        lwr.p2<-as.numeric(tclvalue(tkvars$lwr.p2)); upr.p2<-as.numeric(tclvalue(tkvars$upr.p2))
        mu.p1 <- as.numeric(tclvalue(mu.p1)); s2.p1<- exp(as.numeric(tclvalue(s2.p1))+log(1/12))
        a.p1 <<- mu.p1^2/s2.p1*(1-mu.p1)-mu.p1; b.p1<<-a.p1*(1/mu.p1-1)
        arrparms[[2,5]]<<-signif(a.p1,4); arrparms[[2,6]]<<-signif(b.p1,4)
        mu.p2 <- as.numeric(tclvalue(mu.p2)); s2.p2<- exp(as.numeric(tclvalue(s2.p2))+log(1/12))
        a.p2 <<- mu.p2^2/s2.p2*(1-mu.p2)-mu.p2; b.p2<<-a.p2*(1/mu.p2-1)
        arrparms[[3,5]]<<-signif(a.p2,4); arrparms[[3,6]]<<-signif(b.p2,4)
        wt.u<-as.numeric(tclvalue(tkvars$wt.u))*Udo
        wt.p1<-as.numeric(tclvalue(tkvars$wt.p1))*P1do
        wt.p2<-as.numeric(tclvalue(tkvars$wt.p2))*P2do
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
      arrfig <- tkrplot::tkrplot(.Rvar$arrProcess, fun = plotarr, hscale=2.5, vscale=1.4)
      onChange <- function(...) {
        tkrplot::tkrreplot(arrfig)
      }
      # radio buttons for defining model form
      arrcomponents<-toR(tkvars$arrcomponents)
      modelDefinition<-tkframe(.Rvar$arrProcess) # frame for holding the model definition options
      Ulab<-tklabel(modelDefinition,text="uniform")
      Uradio<-tclVar(ifelse(arrcomponents[1],"yes","no"))
      UyesRadio<-tkradiobutton(modelDefinition,variable=Uradio,value="yes", command=function(){
        if (!dateset) return(F)
        tkconfigure(Uwt,state='normal',troughcolor=uclr)
        tkconfigure(Ulwr,state='normal',troughcolor=uclr)
        tkconfigure(Uupr,state='normal',troughcolor=uclr)
        tcl(parmTable,"tag", "rowtag", "U","1")
        # redraw the figure with the uniform added back in
        Udo<<-T
        setWts()
        tkrplot::tkrreplot(arrfig)
      })
      UnoRadio<-tkradiobutton(modelDefinition,variable=Uradio,value="no",command=function(){
        if (!dateset) return(F)
        tkconfigure(Uwt,state='disabled',troughcolor=disclr)
        tkconfigure(Ulwr,state='disabled',troughcolor=disclr)
        tkconfigure(Uupr,state='disabled',troughcolor=disclr)
        tcl(parmTable,"tag", "rowtag", "hide","1")
        # redraw the figure with the uniform taken out
        Udo<<-F
        setWts()
        tkrplot::tkrreplot(arrfig)
      })
      P1lab<-tklabel(modelDefinition,text="pulse1")
      P1radio<-tclVar(ifelse(arrcomponents[2],"yes","no"))
      P1yesRadio<-tkradiobutton(modelDefinition,variable=P1radio,value="yes", command=function(){
        if (!dateset) return(F)
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
      P1noRadio<-tkradiobutton(modelDefinition,variable=P1radio,value="no", command = function(){
        if (!dateset) return(F)
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
      P2lab<-tklabel(modelDefinition,text="pulse2")
      P2radio<-tclVar(ifelse(arrcomponents[3],"yes","no"))
      P2yesRadio<-tkradiobutton(modelDefinition,variable=P2radio,value="yes", command=function(){
        if (!dateset) return(F)
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
      P2noRadio<-tkradiobutton(modelDefinition,variable=P2radio,value="no", command=function(){
        if (!dateset) return(F)
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
      tkgrid(Ulab,UyesRadio,UnoRadio,rowspan=2)
      tkgrid(P1lab,P1yesRadio,P1noRadio,rowspan=2)
      tkgrid(P2lab,P2yesRadio,P2noRadio,rowspan=2)

      processBounds<-tkframe(.Rvar$arrProcess)

      Ulwr <- tkscale(processBounds, from = 0, to = 365, variable = tkvars$lwr.u, orient = "horizontal",
        length = scwid,
        width=barwidth,
        command = function(...) {
          if (as.numeric(tclvalue(tkvars$lwr.u))>=as.numeric(tclvalue(tkvars$upr.u))){
            tclvalue(tkvars$lwr.u)<-as.numeric(tclvalue(tkvars$upr.u))-1
            arrparms[[1,3]]<<-round(as.numeric(tclvalue(tkvars$lwr.u)))
            return(F)
          }
          arrparms[[1,3]]<<-round(as.numeric(tclvalue(tkvars$lwr.u)))
          onChange()
        },
        resolution = 1,
        sliderrelief=relief,
        sliderlength=8,
        troughcolor=uclr,
        showvalue=F,
        borderwidth=wbord
      )
      Uupr <- tkscale(processBounds, from = 0, to = 365, variable = tkvars$upr.u, orient = "horizontal",
        length = scwid,
        width=barwidth,
        command = function(...) {
          if (as.numeric(tclvalue(tkvars$lwr.u))>=as.numeric(tclvalue(tkvars$upr.u))){
            tclvalue(tkvars$upr.u)<-as.numeric(tclvalue(tkvars$lwr.u))+1
            arrparms[[1,4]]<<-round(as.numeric(tclvalue(tkvars$lwr.u)))
            return(F)
          }
          arrparms[[1,4]]<<-round(as.numeric(tclvalue(tkvars$upr.u)))
          onChange()
        },
        resolution = 1,
        sliderrelief=relief,
        sliderlength=8,
        troughcolor=uclr,
        showvalue=F,
        borderwidth=wbord
      )
      P1lwr <- tkscale(processBounds, from = 0, to = 365, variable = tkvars$lwr.p1, orient = "horizontal",
        length = scwid,
        width=barwidth,
        command = function(...) {
          if (as.numeric(tclvalue(tkvars$lwr.p1))>=as.numeric(tclvalue(tkvars$upr.p1))){
            tclvalue(tkvars$lwr.p1)<-as.numeric(tclvalue(tkvars$upr.p1))-1
            arrparms[[2,3]]<<-round(as.numeric(tclvalue(tkvars$lwr.p1)))
            return(F)
          }
          arrparms[[2,3]]<<-round(as.numeric(tclvalue(tkvars$lwr.p1)))
          onChange()
        },
        resolution = 1,
        sliderrelief=relief,
        sliderlength=8,
        troughcolor=p1clr,
        showvalue=F,
        borderwidth=wbord
      )
      P1upr <- tkscale(processBounds, from = 0, to = 365, variable = tkvars$upr.p1, orient = "horizontal",
        length = scwid,
        width=barwidth,
        command = function(...) {
          if (as.numeric(tclvalue(tkvars$lwr.p1))>=as.numeric(tclvalue(tkvars$upr.p1))){
            tclvalue(tkvars$upr.p1)<-as.numeric(tclvalue(tkvars$lwr.p1))+1
            arrparms[[2,4]]<<-round(as.numeric(tclvalue(tkvars$upr.p1)))
            return(F)
          }
          arrparms[[2,4]]<<-round(as.numeric(tclvalue(tkvars$upr.p1)))
          onChange()
        },
        resolution = 1,
        sliderrelief=relief,
        sliderlength=8,
        troughcolor=p1clr,
        showvalue=F,
        borderwidth=wbord
      )
      P2lwr <- tkscale(processBounds, from = 0, to = 365, variable = tkvars$lwr.p2, orient = "horizontal",
        length = scwid,
        width=barwidth,
        command = function(...) {
          if (as.numeric(tclvalue(tkvars$lwr.p2))>=as.numeric(tclvalue(tkvars$upr.p2))){
            tclvalue(tkvars$lwr.p2)<-as.numeric(tclvalue(tkvars$upr.p2))-1
            arrparms[[3,3]]<<-round(as.numeric(tclvalue(tkvars$lwr.p2)))
            return(F)
          }
          arrparms[[3,3]]<<-round(as.numeric(tclvalue(tkvars$lwr.p2)))
          onChange()
        },
        resolution = 1,
        sliderrelief=relief,
        sliderlength=8,
        troughcolor=p2clr,
        showvalue=F,
        borderwidth=wbord
      )
      P2upr <- tkscale(processBounds, from = 0, to = 365, variable = tkvars$upr.p2, orient = "horizontal",
        length = scwid,
        width=barwidth,
        command = function(...) {
          if (as.numeric(tclvalue(tkvars$lwr.p2))>=as.numeric(tclvalue(tkvars$upr.p2))){
            tclvalue(tkvars$upr.p2)<-as.numeric(tclvalue(tkvars$lwr.p2))+1
            arrparms[[3,4]]<<-round(as.numeric(tclvalue(tkvars$upr.p2)))
            return(F)
          }
          arrparms[[3,4]]<<-round(as.numeric(tclvalue(tkvars$upr.p2)))
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
        ww<-as.numeric(tclvalue(tkvars$wt.u))
        wo<-as.numeric(tclvalue(tkvars$wt.p1))*P1do+as.numeric(tclvalue(tkvars$wt.p2))*P2do
        arrparms[[1,2]]<<-round(ifelse(wo==0,1,ww/(wo+ww)),2)
        ww<-as.numeric(tclvalue(tkvars$wt.p1))
        wo<-as.numeric(tclvalue(tkvars$wt.u))*Udo+as.numeric(tclvalue(tkvars$wt.p2))*P2do
        arrparms[[2,2]]<<-round(ifelse(wo==0,1,ww/(wo+ww)),2)
        ww<-as.numeric(tclvalue(tkvars$wt.p2))
        wo<-as.numeric(tclvalue(tkvars$wt.u))*Udo+as.numeric(tclvalue(tkvars$wt.p1))*P1do
        arrparms[[3,2]]<<-round(ifelse(wo==0,1,ww/(wo+ww)),2)
      }
      Uwt<-tkscale(processWts, from = 1, to = 0, variable = tkvars$wt.u, orient = "vertical",
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
      P1wt<-tkscale(processWts, from = 1, to = 0, variable = tkvars$wt.p1, orient = "vertical",
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
      P2wt<-tkscale(processWts, from = 1, to = 0, variable = tkvars$wt.p2, orient = "vertical",
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
        tclvalue(arrstart)<-0
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
        tclvalue(arrstart)<-adat$arrstart
        tclvalue(tkvars$lwr.u)<-adat$lwr.u
        tclvalue(tkvars$upr.u)<-adat$upr.u
        tclvalue(tkvars$lwr.p1)<-adat$lwr.p1
        tclvalue(tkvars$upr.p1)<-adat$upr.p1
        tclvalue(tkvars$lwr.p2)<-adat$lwr.p2
        tclvalue(tkvars$upr.p2)<-adat$upr.p2
        a.p1<-adat$a.p1; b.p1<-adat$b.p1
        tclvalue(mu.p1)<-a.p1/(a.p1+b.p1)
        tclvalue(s2.p1)<-log(12*a.p1*b.p1/((a.p1+b.p1)^2*(a.p1+b.p1+1)))
        a.p2<-adat$a.p2; b.p2<-adat$b.p2
        tclvalue(mu.p2)<-a.p2/(a.p2+b.p2)
        tclvalue(s2.p2)<-log(12*a.p2*b.p2/((a.p2+b.p2)^2*(a.p2+b.p2+1)))
        tclvalue(tkvars$wt.u)<-adat$wt.u
        tclvalue(tkvars$wt.p1)<-adat$wt.p1
        tclvalue(tkvars$wt.p2)<-adat$wt.p2
        arrparms[[1,0]]<-as.tclObj('',drop=T)
        arrparms[[1,1]]<-"resident"; arrparms[[1,2]]<-tclvalue(tkvars$wt.u) ; arrparms[[1,3]]<-tclvalue(tkvars$lwr.u) ; arrparms[[1,4]]<-tclvalue(tkvars$upr.u); arrparms[[1,5]]<-1   ; arrparms[[1,6]]<-1
        arrparms[[2,1]]<-"pulse1";   arrparms[[2,2]]<-round(as.numeric(tclvalue(tkvars$wt.p1)),3); arrparms[[2,3]]<-round(as.numeric(tclvalue(tkvars$lwr.p1)),3); arrparms[[2,4]]<-round(as.numeric(tclvalue(tkvars$upr.p1)),3);
        arrparms[[2,5]]<-signif(a.p1,4); arrparms[[2,6]]<-signif(b.p1,4)
        arrparms[[3,1]]<-"pulse2";   arrparms[[3,2]]<-round(as.numeric(tclvalue(tkvars$wt.p2)),3); arrparms[[3,3]]<-round(as.numeric(tclvalue(tkvars$lwr.p2)),3); arrparms[[3,4]]<-round(as.numeric(tclvalue(tkvars$upr.p2)),3)
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
        tkvars$arrcomponents[[0]]<-1
        tkvars$arrcomponents[[1]]<-1
        tkvars$arrcomponents[[2]]<-1
        tkrplot::tkrreplot(arrfig)
      })
      arrSave<-tkbutton(arrButtFrame,text="Save", width=15,command=function(){
        tclvalue(tkvars$arrfun)<-"Compound"
#        tclvalue(tkvars$arrcomponents[[0]])<-Udo; tclvalue(tkvars$arrcomponents[[1]])<-P1do; tclvalue(tkvars$arrcomponents[[2]])<-P2do
        tkvars$arrcomponents[[0]]<-ifelse(tclvalue(Uradio)=='yes', 1, 0)
        tkvars$arrcomponents[[1]]<-ifelse(tclvalue(P1radio)=='yes', 1, 0)
        tkvars$arrcomponents[[2]]<-ifelse(tclvalue(P2radio)=='yes', 1, 0)
        mu<-as.numeric(tclvalue(mu.p1))
        s2<-exp(as.numeric(tclvalue(s2.p1))+log(1/12))
        tclvalue(tkvars$a.p1)<-mu^2/s2*(1-mu)-mu
        tclvalue(tkvars$b.p1)<-as.numeric(tclvalue(tkvars$a.p1))*(1/mu-1)
        mu<-as.numeric(tclvalue(mu.p2))
        s2<-exp(as.numeric(tclvalue(s2.p2))+log(1/12))
        tclvalue(tkvars$a.p2)<-mu^2/s2*(1-mu)-mu
        tclvalue(tkvars$b.p2)<-as.numeric(tclvalue(tkvars$a.p2))*(1/mu-1)
        tkdestroy(.Rvar$arrProcess)
        tkrplot::tkrreplot(arrfig.mini)
      })
      arrCancel<-tkbutton(arrButtFrame,text='Cancel', width=15,command=function()tkdestroy(.Rvar$arrProcess))
      tkgrid(arrRestore)
      tkgrid(arrSave)
      tkgrid(arrCancel)
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
#      if (parent$name != "design") tkgrid(startlbl,arrReset,sticky='w') else tkgrid(startlbl, sticky='w')
      tkgrid(parmTable,columnspan=2,sticky='sw',pady=10,padx=10)
      tkgrid(arrButtFrame, column=1,row=1)
      tkgrid(modelDefinition,processBounds)
      tkgrid(processWts, arrfig, processVariances, sticky='ne',padx=10)
      tkgrid(processMeans, column=1,sticky='n')
      blankframe<-tklabel(.Rvar$arrProcess,text=as.tclObj('        ',drop=T))
      tkgrid(blankframe)
    },
    addsymc = function(){
      # update the table
      tkinsert(classTable, "rows", "end", 1)
      rowind<-as.numeric(tclvalue(tcl(classTable,"index","end","row")))
      classData[[rowind, 2]] <<- .Rvar$syresult$a #value is entered by hand from the symc form
      classData[[rowind, 3]] <<- .Rvar$syresult$X
      if (.Rvar$syscPrevious$arrfun=='Uniform' & !basicMode){
        Bab<-.Rvar$syresult$Bab
       } else {
        Bab<-.Rvar$syresult$BabAnn
      }
      classData[[rowind, 4]] <<- signif(Bab[1],5)
      classData[[rowind, 5]] <<- signif(Bab[2],5)
      classData[[rowind, 6]] <<- round(Bab[1]/sum(Bab),3)
      classData[[rowind, 7]] <<- as.tclObj(paste0("[", round(qbeta(0.025,Bab[1],Bab[2]),3),", ", round(qbeta(0.975,Bab[1],Bab[2]),3),"]"),drop=T)
      tcl(classTable,"activate", as.tclObj(paste0(rowind,',', 1)))
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
        custom.prior<<-toR(tkvars$custom.prior)
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
        tclvalue(tkvars$prior_f)<<-"Objective"
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
    startchk = function(startparm){
      if (length(startparm)!=1 || nchar(startparm)==0 || class(try(as.Date(startparm),silent=T))=="try-error"){ # text box is empty OR not a date
        startok<<-F
        return(F)
      } else {
        startok<<-T
    #    singleYearProvisional$firstsearch<<-startparm
      }
      return(T)
    }
  )
)
symcsave<-function(symc, sysc){
  save(symc, file = paste0(.Rvar$datadir,"/symcPrevious.Rdata"))
  save(sysc,     file = paste0(.Rvar$datadir,"/syscrevious.Rdata"))
}
