.Rvar <- new.env(parent=emptyenv()) # an environment for R variables to live in
.Rvar$platform <- ifelse(.Platform$OS.type == 'windows', 'windows', ifelse(Sys.info()['sysname'] == "Darwin",  'mac', 'linux'))

tcltk::tcl("option","add","*tearOff",0)
.onLoad<-function(libname, pkgname){
  list.of.packages <- c("R6", "matrixStats", "rjags", "tensorA", "actuar","VGAM", "tcltk2", "tkrplot", "gsl")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  .Rvar$datadir<-file.path(system.file(package = "eoa"), "data")
  if (file.exists(paste0(.Rvar$datadir,"/dtPrevious.Rdata"))){
    .Rvar$dtPrevious<-get(load(paste0(.Rvar$datadir,"/dtPrevious.Rdata")))
  } else {
    .Rvar$dtPrevious<-dtDefault
    dtPrevious<-dtDefault
    save(dtPrevious,file = paste0(.Rvar$datadir,"/dtPrevious.Rdata"))
  }
  if (file.exists(paste0(.Rvar$datadir,"/singleYearPrevious.Rdata"))){
#    .Rvar$singleYearPrevious<-get(load(paste0(.Rvar$datadir,"/singleYearPrevious.Rdata")))
    load(paste0(.Rvar$datadir,"/singleYearPrevious.Rdata"), envir = .Rvar)
  } else {
    .Rvar$singleYearPrevious<-singleYearDefault
    singleYearPrevious<-singleYearDefault
    save(singleYearPrevious,file = paste0(.Rvar$datadir,"/singleYearPrevious.Rdata"))
  }
  if (file.exists(paste0(.Rvar$datadir,"/scexPrevious.Rdata"))){
    .Rvar$scexPrevious<-get(load(paste0(.Rvar$datadir,"/scexPrevious.Rdata")))
  } else {
    .Rvar$scexPrevious<-scexDefault
    scexPrevious<-scexDefault
    save(scexPrevious,file = paste0(.Rvar$datadir,"/scexPrevious.Rdata"))
  }
  if (file.exists(paste0(.Rvar$datadir,"/symcPrevious.Rdata"))){
    .Rvar$symcPrevious<-get(load(paste0(.Rvar$datadir,"/symcPrevious.Rdata")))
  } else {
    .Rvar$symcPrevious<-symcDefault
    symcPrevious<-symcDefault
    save(symcPrevious,file = paste0(.Rvar$datadir,"/symcPrevious.Rdata"))
  }
  if (file.exists(paste0(.Rvar$datadir,"/myPrevious.Rdata"))){
    .Rvar$myPrevious<-get(load(paste0(.Rvar$datadir,"/myPrevious.Rdata")))
  } else {
    .Rvar$myPrevious<-myDefault
    myPrevious<-myDefault
    save(myPrevious,file = paste0(.Rvar$datadir,"/myPrevious.Rdata"))
  }
  if (file.exists(paste0(.Rvar$datadir,"/syscPrevious.Rdata"))){
    .Rvar$syscPrevious<-get(load(paste0(.Rvar$datadir,"/syscPrevious.Rdata")))
  } else {
    .Rvar$syscPrevious<-syscDefault
    syscPrevious<-syscDefault
    save(syscPrevious,file = paste0(.Rvar$datadir,"/syscPrevious.Rdata"))
  }
  if (file.exists(paste0(.Rvar$datadir,"/CPdataPrevious.Rdata"))){
    .Rvar$CPdataPrevious<-get(load(paste0(.Rvar$datadir,"/CPdataPrevious.Rdata")))
  } else {
    .Rvar$CPdataPrevious<-CPdataDefault
    CPdataPrevious<-CPdataDefault
    save(CPdataPrevious,file = paste0(.Rvar$datadir,"/CPdataPrevious.Rdata"))
  }
  if (file.exists(paste0(.Rvar$datadir,"/pkdatPrevious.Rdata"))){
    .Rvar$pkdatPrevious<-get(load(paste0(.Rvar$datadir,"/pkdatPrevious.Rdata")))
  } else {
    .Rvar$pkdatPrevious<-pkdatDefault
    pkdatPrevious<-pkdatDefault
    save(pkdatPrevious,file = paste0(.Rvar$datadir,"/pkdatPrevious.Rdata"))
  }
  if (file.exists(paste0(.Rvar$datadir,"/csvpath.Rdata"))){
    .Rvar$csvpath<-get(load(paste0(.Rvar$datadir,"/csvpath.Rdata")))
  } else {
    assign("csvpath", getwd(), env = .Rvar)
  }
}
.onAttach<-function(libname, pkgname){
  if (file.exists(file.path(system.file(package = "eoa"), "data/csvpath.Rdata"))){
    load(file = file.path(system.file(package = "eoa"), "data/csvpath.Rdata"), envir = as.environment("package:eoa"))
  } else {
    assign("csvpath", value = getwd(), envir = as.environment("package:eoa"))
  }
  packageStartupMessage("\nEnter eoa() to begin...")
}
eoa<-function(){
  .Rvar$VER <- packageDescription("eoa", field = "Version")
  .Rvar$DATE <- packageDescription("eoa", field = "Date")
  about_text<-paste0('Evidence of Absence, v', .Rvar$VER, ' (', .Rvar$DATE,
    ')\n  by D. Dalthorp, M. Huso, and D. Dail
      Email ddalthorp@usgs.gov with questions or problems.')
  graphics.off()
  if (.Rvar$platform == 'windows') windows.options(reset = T)
  if (.Rvar$platform == 'mac') quartz.options(reset = T)
  if (.Rvar$platform == 'linux') X11.options(reset = T)

  .Rvar$charSizeAdjust<-par('din')[2]/7
  graphics.off()
#  .Rvar$charSizeAdjust<-(96/(par('cra')/par('cin'))[2])#/(par('din')[2]/7)
#  loadRconsole(paste0(.Rvar$datadir,"/console.cfg"))
  tmppath<-suppressWarnings(try(load("csvpath.Rdata"),silent=T))
  .Rvar$csvpath<-ifelse(class(tmppath)=="try-error", getwd(), csvpath)
  if (.Rvar$platform == 'windows') loadRconsole(paste0(.Rvar$datadir,"/console.cfg"))
  tcltk::tclRequire("Tktable")
# close all previously open eoa windows...
# sy windows...
  try(tkdestroy(.Rvar$syWindow$cdFrame), silent = T)
  try(tkdestroy(.Rvar$syWindow$cprFrame), silent = T)
  try(tkdestroy(.Rvar$syWindow$fcpModule), silent = T)
  try(tkdestroy(.Rvar$syWindow$pkModule), silent = T)
  try(tkdestroy(.Rvar$arrProcess), silent=T)
  try(tkdestroy(.Rvar$syWindow$syModule), silent = T)
# symc windows
  try(tkdestroy(.Rvar$symcWindow$cprFrame), silent = T)
  try(tkdestroy(.Rvar$symcWindow$syCompClass$fcpModule), silent = T)
  try(tkdestroy(.Rvar$symcWindow$syCompClass$pkModule), silent = T)
  try(tkdestroy(.Rvar$symcWindow$syCompClass$syModule), silent = T)
  try(tkdestroy(.Rvar$symcWindow$symcModule), silent = T)
# my windows
  try(tkdestroy(.Rvar$myWindow$myModule), silent = T)
# dt windows
  try(tkdestroy(.Rvar$viewMvg), silent = T)
  try(tkdestory(.Rvar$dtWindowdtModule), silent = T)
# scex windows
  try(tkdestroy(.Rvar$scexWindow$scexModule), silent=T)
  try(tkdestroy(.Rvar$scex1mod$scex1Module), silent = T)
# eoa window
  try(tkdestroy(.Rvar$EoA), silent=T)

#  tcltk::tkmessageBox(message=
#"This software is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. The software is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the software.\n")
  .Rvar$EoA <- tktoplevel()
  tkwm.resizable(.Rvar$EoA,FALSE,FALSE)
  tkwm.minsize(.Rvar$EoA)
  topMenu <- tkmenu(.Rvar$EoA); tkconfigure(.Rvar$EoA,menu=topMenu)
  helpMenu <- tkmenu(topMenu,activebackground=colors()[125],activeforeground=colors()[109])
  tkadd(helpMenu,"command", label="About",
    command = function() tkmessageBox(title=paste0('Evidence of Absence, v', .Rvar$VER), message=about_text))
  tkadd(topMenu,"cascade",label="Help",menu=helpMenu)
  titletext<-paste("EoA - Evidence of Absence, v",.Rvar$VER,sep='')
  tkwm.title(.Rvar$EoA,titletext)  # Give the window a title
  ### .Rvar$EoA: .Rvar$EoA to hold buttons for opening modules
  # text for welcoming and brief instructions in column 1 on the front page
  welcome.column1<-tkframe(.Rvar$EoA);
  welcome.column1lbl<-tklabel(welcome.column1,text='Control Panel',justify='center',padx=17,width=20); tkgrid(welcome.column1lbl)
  # column of buttons for opening other modules
  welcome.buttonwidth<-25
  welcome.column2<-tkframe(.Rvar$EoA)
  welcome.singleYearbutton<-tkbutton(welcome.column2,text='Single Class', width=welcome.buttonwidth, command = function() {
    load(paste0(.Rvar$datadir,'/singleYearPrevious.Rdata'), envir = .Rvar)
    .Rvar$syWindow<-syForm$new(.Rvar$singleYearPrevious)
    tkwm.state(.Rvar$EoA, "withdraw")
  })
  welcome.mclassbutton<-tkbutton(welcome.column2,text='Multiple Classes', width=welcome.buttonwidth, command = function() {
    load(paste0(.Rvar$datadir,'/symcPrevious.Rdata'), envir = .Rvar)
    .Rvar$symcWindow<-symcForm$new(.Rvar$symcPrevious)
    tkwm.state(.Rvar$EoA, "withdraw")
  });
  welcome.myearbutton<-tkbutton(welcome.column2,text='Multiple Years', width=welcome.buttonwidth, command = function() {
    .Rvar$myWindow<-myForm$new(.Rvar$myPrevious)
    tkwm.state(.Rvar$EoA, "withdraw")
  });
  welcome.designTradeoffsbutton<-tkbutton(welcome.column2,text='Design Tradeoffs', width=welcome.buttonwidth, command = function(){
    .Rvar$dtWindow<-dtForm$new(.Rvar$dtPrevious)
    tkwm.state(.Rvar$EoA, "withdraw")
  })
  welcome.scenariosbutton<-tkbutton(welcome.column2,text='Long-term Scenario Explorer', width=welcome.buttonwidth, command = function(){
    .Rvar$scexWindow<-scexForm$new(.Rvar$scexPrevious)
    tkwm.state(.Rvar$EoA, "withdraw")
  });
  welcome.defaultsbutton<-tkbutton(welcome.column2,text='Restore Default Data Sets', width=welcome.buttonwidth,  bg = colors()[238], fg = colors()[477],
    command = function(){
    .Rvar$dtPrevious<-dtDefault
    dtPrevious<-dtDefault
    save(dtPrevious,file = paste0(.Rvar$datadir,"/dtPrevious.Rdata"))
    .Rvar$singleYearPrevious<-singleYearDefault
    singleYearPrevious<-singleYearDefault
    save(singleYearPrevious,file = paste0(.Rvar$datadir,"/singleYearPrevious.Rdata"))
    .Rvar$scexPrevious<-scexDefault
    scexPrevious<-scexDefault
    save(scexPrevious,file = paste0(.Rvar$datadir,"/scexPrevious.Rdata"))
    .Rvar$symcPrevious<-symcDefault
    symcPrevious<-symcDefault
    save(symcPrevious,file = paste0(.Rvar$datadir,"/symcPrevious.Rdata"))
    .Rvar$myPrevious<-myDefault
    myPrevious<-myDefault
    save(myPrevious,file = paste0(.Rvar$datadir,"/myPrevious.Rdata"))
    .Rvar$syscPrevious<-syscDefault
    syscPrevious<-syscDefault
    save(syscPrevious,file = paste0(.Rvar$datadir,"/syscPrevious.Rdata"))
    .Rvar$CPdataPrevious<-CPdataDefault
    CPdataPrevious<-CPdataDefault
    save(CPdataPrevious,file = paste0(.Rvar$datadir,"/CPdataPrevious.Rdata"))
    .Rvar$pkdatPrevious<-pkdatDefault
    pkdatPrevious<-pkdatDefault
    save(pkdatPrevious,file = paste0(.Rvar$datadir,"/pkdatPrevious.Rdata"))
    if (file.exists(paste0(.Rvar$datadir,"/csvpath.Rdata"))){
      .Rvar$csvpath<-get(load(paste0(.Rvar$datadir,"/csvpath.Rdata")))
    } else {
      assign("csvpath", getwd(), env = .Rvar)
    }
    tkconfigure(welcome.defaultsbutton, bg = colors()[81], text = "Restoring...", fg = colors()[1])
    tclAfter(500, function() tkconfigure(welcome.defaultsbutton, text='Restore Default Data Sets', bg = "SystemButtonFace", fg = "SystemWindowText") )
  });
  tkbind(.Rvar$EoA,"<Destroy>",function(){
    if (length(ls(env=.GlobalEnv))==0) quit(save="no")
  })

  tkgrid(welcome.singleYearbutton)
  tkgrid(welcome.mclassbutton)
  tkgrid(welcome.myearbutton)
  tkgrid(welcome.designTradeoffsbutton)
  tkgrid(welcome.scenariosbutton)
  tkgrid(welcome.defaultsbutton, pady = c(15,5))
  tkgrid(welcome.column1,welcome.column2)
  tkwm.deiconify(.Rvar$EoA);  tkgrab.set(.Rvar$EoA);  tkfocus(.Rvar$EoA)
}
