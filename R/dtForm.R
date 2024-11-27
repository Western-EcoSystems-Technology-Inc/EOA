dtForm<-R6::R6Class("dtForm",
  portable = FALSE,
  public = list(
    # global tcl variables
    dtModule = NA,
    arrfig.mini = NA,
    pdb.lbl = NA, pda.lbl = NA,
    tkvars = list(),

    # R variables:
    persistence_distn = NA, arrcomponents = NA, name = "design",

    # data checks...
    target_gok = T, grincok = T, tarrok = T,
    fok = T, fleg = T, fminleg = T, fmaxleg = T, # SE
    schedok = T, Ileg = T, Iminleg = T, Imaxleg = T, spanleg = T,  # I
    kok = T, startok = T, tauok = T, alphaok = T, Xok = T,
    phiok = T, phileg = T, phiminleg = T, phimaxleg = T, # coverage
    lamok = T,  pdaok = T, pdbok = T,  # persistence
    arrok = T,  #arrivals

    # initial values
    initialize = function(dtdat) {
      # set data set values to T b/c dtDefault, dtPrevious have been error-checked previously
      # if other saved data set is loaded (.rds), then subsequent error-check sets flags
       allokdt = function(){
          target_gok <<- T; grincok <<- T;  tarrok<<-T;
          fok <<- T; fleg <<- T; fminleg <<- T; fmaxleg <<- T; # SE
          schedok <<- T; Ileg <<- T; Iminleg <<- T; Imaxleg <<- T; spanleg <<- T;  # I
          kok <<- T;    startok <<- T;
          phiok <<- T; phileg <<- T; phiminleg <<- T; phimaxleg <<- T; # coverage
          lamok <<- T;  pdaok <<- T; pdbok <<- T;  # persistence
          arrok <<- T;  #arrivals
        }
        allokdt()
      # create the required tcl variables from the data input
      for (i in 1:length(dtVar)){
        ind <- which(names(dtdat)==dtVar[i])
        tkvars[[dtVar[i]]] <<- tclVar(as.tclObj(ifelse(length(ind) > 0, dtdat[[ind]],""),drop=T))
      }
      for (li in 1:length(dtArray)){
        ind <- which(names(dtdat)==dtArray[li])
        tmp<-tclArray()
        if (length(ind) == 0){ # missing data
          tmp[[0]] <- as.tclObj('', drop=T)
        } else if (length(dtdat[[ind]]) == 1) { # scalar
          tmp[[0]] <- as.tclObj(ifelse(!is.na(dtdat[[ind]]), dtdat[[ind]], ''), drop=T)
        } else if (is.null(dim(dtdat[[ind]]))) { # vector
          for (j in 1:length(dtdat[[ind]])) tmp[[j-1]] <- as.tclObj(dtdat[[ind]][j],drop=T)
        } else { # matrix
          for (rowi in 1:dim(dtdat[[ind]])[1])
            for (coli in 1:dim(dtdat[[ind]])[2])
              tmp[[rowi-1, coli-1]] <- as.tclObj(dtdat[[ind]][rowi, coli], drop = T)
        }
        tkvars[[dtArray[li]]] <<- tmp
      }
      persistence_distn <<- dtdat$persistence_distn
      arrcomponents <<- dtdat$arrcomponents
        tarrchk= function(tarr){
         if (length(tarr) != 1 || nchar(tarr)==0 || is.na(tarr) || tarr<=0 || tarr>1) {
            tkconfigure(tarr.edit, bg=errcolor)
            tkconfigure(drawButton,state='disabled')
            tarrok<<-F
            return(F)
          } else {
            tkconfigure(tarr.edit, bg='white')
            tarrok<<-T
            if (checkokd()) {
              tkconfigure(drawButton,state='normal')
            }
          }
          return(T)
        }
        phichk = function(phi){
          if (length(phi) != 1 || nchar(phi)==0 || is.na(phi) || phi<=0 || phi>1) {
            tkconfigure(phi.edit, bg=errcolor)
            tkconfigure(drawButton,state='disabled')
            phiok<<-F
            return(F)
          } else {
            tkconfigure(phi.edit, bg='white')
            phiok<<-T
            if (checkokd()) {
              tkconfigure(drawButton,state='normal')
            }
          }
          return(T)
        }

        phiminchk = function(phimin){
          if (length(phimin) != 1 || nchar(phimin)==0 || is.na(phimin) || phimin<=0 || phimin > 1) {
            tkconfigure(phimin.edit, bg=errcolor)
            tkconfigure(drawButton,state='disabled')
            if (phimaxleg) tkconfigure(phimax.edit,bg='white')
            phiok<<-F
            phiminleg<<-F
            return(F)
          } else {
            tkconfigure(phimin.edit, bg='white')
            phiminleg<<-T
            if (phimaxleg) { # if phiminleg, then new phimin might change color of phimax.edit (provided it has a legitimate value)
              if (toR(tkvars$phimin)<toR(tkvars$phimax)){
                phiok<<-T # both min and max have legitimate values and min < max
                tkconfigure(phimin.edit,bg='white')
                tkconfigure(phimax.edit,bg='white')
              } else {
                phiok<<-F
                tkconfigure(phimin.edit,bg=compcolor)
                tkconfigure(phimax.edit,bg=compcolor)
                tkconfigure(drawButton,state='disabled')
                return(F)
              }
            }
            if (checkokd()) tkconfigure(drawButton,state='normal')
          }
          return(T)
        }
        phimaxchk = function(phimax){
          if (nchar(phimax)==0 | is.na(phimax) | phimax<=0 | phimax>1 ) {
            tkconfigure(phimax.edit, bg=errcolor)
            tkconfigure(drawButton,state='disabled')
            if (phiminleg) tkconfigure(phimin.edit,bg='white')
            phiok<<-F
            phimaxleg<<-F
            return(F)
          } else {
            tkconfigure(phimax.edit, bg='white')
            phimaxleg<<-T
            if (phiminleg) {
              if (toR(tkvars$phimin)<toR(tkvars$phimax)){
                phiok<<-T
                tkconfigure(phimax.edit,bg='white')
                tkconfigure(phimin.edit,bg='white')
              } else {
                phiok<<-F
                tkconfigure(phimin.edit,bg=compcolor)
                tkconfigure(phimax.edit,bg=compcolor)
                tkconfigure(drawButton,state='disabled')
                return(F)
              }
            }
            if (checkokd()) tkconfigure(drawButton,state='normal')
          }
          return(T)
        }
        ###
        fchk = function(f){
          if (nchar(f)==0 | is.na(f) | f<=0 | f>1) {
            tkconfigure(f.edit, bg=errcolor)
            tkconfigure(drawButton,state='disabled')
            fok<<-F
            return(F)
          } else {
            tkconfigure(f.edit, bg='white')
            fok<<-T
            if (checkokd()) {
              tkconfigure(drawButton,state='normal')
            }
          }
          return(T)
        }
        fminchk = function(fmin){
          if (nchar(fmin)==0 | is.na(fmin) | fmin<=0 | fmin > 1) {
            tkconfigure(fmin.edit, bg=errcolor)
            tkconfigure(drawButton,state='disabled')
            if (fmaxleg) tkconfigure(fmax.edit,bg='white')
            fok<<-F
            fminleg<<-F
            return(F)
          } else {
            tkconfigure(fmin.edit, bg='white')
            fminleg<<-T
            if (fmaxleg) { # if fminleg, then new fmin might change color of fmax.edit (provided it has a legitimate value)
              if (toR(tkvars$fmin)<toR(tkvars$fmax)){
                fok<<-T # both min and max have legitimate values and min < max
                tkconfigure(fmin.edit,bg='white')
                tkconfigure(fmax.edit,bg='white')
              } else {
                fok<<-F
                tkconfigure(fmin.edit,bg=compcolor)
                tkconfigure(fmax.edit,bg=compcolor)
                tkconfigure(drawButton,state='disabled')
                return(F)
              }
            }
            if (checkokd()) tkconfigure(drawButton,state='normal')
          }
          return(T)
        }
        fmaxchk = function(fmax){
          if (nchar(fmax)==0 | is.na(fmax) | fmax<=0 | fmax>1 ) {
            tkconfigure(fmax.edit, bg=errcolor)
            tkconfigure(drawButton,state='disabled')
            if (fminleg) tkconfigure(fmin.edit,bg='white')
            fok<<-F
            fmaxleg<<-F
            return(F)
          } else {
            tkconfigure(fmax.edit, bg='white')
            fmaxleg<<-T
            if (fminleg) {
              if (toR(tkvars$fmin)<toR(tkvars$fmax)){
                fok<<-T
                tkconfigure(fmax.edit,bg='white')
                tkconfigure(fmin.edit,bg='white')
              } else {
                fok<<-F
                tkconfigure(fmin.edit,bg=compcolor)
                tkconfigure(fmax.edit,bg=compcolor)
                tkconfigure(drawButton,state='disabled')
                return(F)
              }
            }
            if (checkokd()) tkconfigure(drawButton,state='normal')
          }
          return(T)
        }
        ###
        dIchk = function(I){
          I<-suppressWarnings(as.numeric(I))
          if (length(I) != 1 || nchar(I)==0 || is.na(I) || I<=0) {
            tkconfigure(Isam.edit, bg=errcolor)
            tkconfigure(drawButton,state='disabled')
            schedok<<-F
            Ileg<<-F
            if (spanleg) tkconfigure(span.edit,bg='white') else tkconfigure(span.edit,bg=errcolor)
            tkconfigure(rCP.lbl1,text=paste0("meanCP = ", NA, ", r = ", NA, " for Ir = ", NA))
            return(F)
          } else {
            Ileg<<-T
            dschedchk()
          }
          if (pdaok & pdbok){
            rCP<-rCPgab(toR(tkvars$persistence_distn),toR(tkvars$pda0),toR(tkvars$pdb0),toR(tkvars$Isam))
            tkconfigure(rCP.lbl1,text=paste0("meanCP = ",signif(rCP[1],3),", r = ",signif(rCP[2],3), " for Ir = ",toR(tkvars$Isam)))
          }
          return(T)
        }
        dIminchk = function(Imin){
          # Imin must pass three tests: (1) positive number, (2) Imin < Imax, (3) span%%Imin == 0
          Imin<-suppressWarnings(as.numeric(Imin))
          if (length(Imin) != 1 || nchar(Imin)==0 || is.na(Imin) || Imin<=0) { # fails test1 --> Imin not a positive number
            Iminleg<<-F
            schedok<<-F
            tkconfigure(Imin.edit, bg=errcolor)
            tkconfigure(drawButton,state='disabled')
            if (Imaxleg & spanleg){
              if(toR(tkvars$span)%%toR(tkvars$Imax)==0){
                 tkconfigure(Imax.edit,bg='white')
                 tkconfigure(span.edit,bg='white')
              } else {
                 tkconfigure(Imax.edit,bg=compcolor)
                 tkconfigure(span.edit,bg=compcolor)
              }
            } else if (!Imaxleg & !spanleg) {
                tkconfigure(Imax.edit,bg=errcolor)
                tkconfigure(span.edit,bg=errcolor)
            } else if (Imaxleg & !spanleg) {
                tkconfigure(Imax.edit,bg='white')
                tkconfigure(span.edit,bg=errcolor)
            } else if (!Imaxleg & spanleg) {
                tkconfigure(Imax.edit,bg=errcolor)
                tkconfigure(span.edit,bg='white')
            }
            return(F)
          } else { # Imin is a positive number
            Iminleg<<-T
            tkconfigure(Imin.edit, bg='white')
            dschedchk()
          }
          return(T)
        }
        dImaxchk = function(Imax){
          # Imin must pass three tests: (1) positive number, (2) Imin < Imax, (3) span%%Imin == 0
          Imax<-suppressWarnings(as.numeric(Imax))
          if (length(Imax) != 1 || nchar(Imax)==0 || is.na(Imax) || Imax<=0) { # fails test1 --> Imax not a positive number
            Imaxleg<<-F
            schedok<<-F
            tkconfigure(Imax.edit, bg=errcolor)
            tkconfigure(drawButton,state='disabled')
            if (Iminleg & spanleg){
              if(toR(tkvars$span)%%toR(tkvars$Imin)==0){
                 tkconfigure(Imin.edit,bg='white')
                 tkconfigure(span.edit,bg='white')
              } else {
                 tkconfigure(Imin.edit,bg=compcolor)
                 tkconfigure(span.edit,bg=compcolor)
              }
            } else if (!Iminleg & !spanleg) {
                tkconfigure(Imin.edit,bg=errcolor)
                tkconfigure(span.edit,bg=errcolor)
            } else if (Iminleg & !spanleg) {
                tkconfigure(Imin.edit,bg='white')
                tkconfigure(span.edit,bg=errcolor)
            } else if (!Iminleg & spanleg) {
                tkconfigure(Imin.edit,bg=errcolor)
                tkconfigure(span.edit,bg='white')
            }
            return(F)
          } else { # Imax is a positive number
            Imaxleg<<-T
            tkconfigure(Imax.edit, bg='white')
            dschedchk()
          }
          return(T)
        }
        dschedchk = function(){
          Ifix<-toR(tkvars$Ifix)
          span<-suppressWarnings(as.numeric(toR(tkvars$span)))
          tkconfigure(drawButton,state='disabled') # will turn back on if all data pass error-checking
          if (Ifix) {
            Isam<-toR(tkvars$Isam)
          } else {
            Imin<-toR(tkvars$Imin)
            Imax<-toR(tkvars$Imax)
          }
          if (length(span) != 1 || nchar(span)==0 || is.na(span) || span<=0 ) {
            spanleg<<-F
            schedok<<-F
            tkconfigure(span.edit, bg=errcolor)
            if (Ifix==1){
              if(Ileg) tkconfigure(Isam.edit,bg='white') else tkconfigure(Isam.edit,bg=errcolor)
            }
            if (Ifix==0){
              if (Iminleg & Imaxleg){
                if (Imin < Imax){
                  tkconfigure(Imin.edit,bg='white')
                  tkconfigure(Imax.edit,bg='white')
                } else {
                  tkconfigure(Imin.edit,bg=compcolor)
                  tkconfigure(Imax.edit,bg=compcolor)
                }
              } else {
                if (Iminleg) tkconfigure(Imin.edit,bg='white') else tkconfigure(Imin.edit,bg=errcolor)
                if (Imaxleg) tkconfigure(Imax.edit,bg='white') else tkconfigure(Imax.edit,bg=errcolor)
              }
            }
            return(F)
          } else {
        #    designProvisional$span<<-span
            spanleg<<-T
            if (Ifix==1){
              if (Ileg) {
                if (span%%Isam==0){
                  schedok<<-T
                  tkconfigure(span.edit,bg='white')
                  tkconfigure(Isam.edit,bg='white')
                } else {
                  schedok<<-F
                  tkconfigure(Isam.edit,bg=compcolor)
                  tkconfigure(span.edit,bg=compcolor)
                  return(F)
                }
              } else {
                  tkconfigure(Isam.edit,bg=errcolor)
                  tkconfigure(span.edit,bg='white')
                  return(F)
              }
            } else {# Ifix == 0
              if (Iminleg & Imaxleg) {
                if (span%%Imin==0 & span%%Imax==0 & Imin<Imax){
                  tkconfigure(span.edit,bg='white')
                  tkconfigure(Imin.edit,bg='white')
                  tkconfigure(Imax.edit,bg='white')
                } else if (span%%Imin==0 & span%%Imax!=0 & Imin<Imax) {
                  tkconfigure(span.edit,bg=compcolor)
                  tkconfigure(Imin.edit,bg='white')
                  tkconfigure(Imax.edit,bg=compcolor)
                  return(F)
                } else if (span%%Imin!=0 & span%%Imax==0 & Imin<Imax) {
                  tkconfigure(span.edit,bg=compcolor)
                  tkconfigure(Imin.edit,bg=compcolor)
                  tkconfigure(Imax.edit,bg='white')
                  return(F)
                } else if (span%%Imin!=0 & span%%Imax!=0 & Imin<Imax) {
                  tkconfigure(span.edit,bg=compcolor)
                  tkconfigure(Imin.edit,bg=compcolor)
                  tkconfigure(Imax.edit,bg=compcolor)
                  return(F)
                } else if (span%%Imin==0 & span%%Imax==0 & Imin>=Imax){
                  tkconfigure(span.edit,bg='white')
                  tkconfigure(Imin.edit,bg=compcolor)
                  tkconfigure(Imax.edit,bg=compcolor)
                  return(F)
                } else {
                  tkconfigure(span.edit,bg=compcolor)
                  tkconfigure(Imin.edit,bg=compcolor)
                  tkconfigure(Imax.edit,bg=compcolor)
                  return(F)
                }
              } else if (!Iminleg & !Imaxleg) {
                tkconfigure(span.edit,bg='white')
                tkconfigure(Imin.edit,bg=errcolor)
                tkconfigure(Imax.edit,bg=errcolor)
                return(F)
              } else if (!Iminleg) {
                tkconfigure(Imin.edit,bg=errcolor)
                if (span%%Imax==0){
                  tkconfigure(span.edit,bg='white')
                  tkconfigure(Imax.edit,bg='white')
                } else {
                  tkconfigure(span.edit,bg=compcolor)
                  tkconfigure(Imax.edit,bg=compcolor)
                }
                return(F)
              } else if (!Imaxleg) {
                tkconfigure(Imax.edit,bg=errcolor)
                if (span%%Imin==0){
                  tkconfigure(span.edit,bg='white')
                  tkconfigure(Imin.edit,bg='white')
                } else {
                  tkconfigure(span.edit,bg=compcolor)
                  tkconfigure(Imin.edit,bg=compcolor)
                }
                return(F)
              }
            }
            schedok<<-T
            if (checkokd()) tkconfigure(drawButton,state='normal')
          }
          return(T)
        }
        startchk = function(startparm){
          if (length(startparm) != 1 || nchar(startparm)==0 || class(try(as.Date(startparm),silent=T))=="try-error"){ # text box is empty OR not a date
            startok<<-F
            tkconfigure(firstsearch.edit,bg=errcolor)
            tkconfigure(drawButton,state='disabled')
            return(F)
          } else {
            tkconfigure(firstsearch.edit,bg='white')
            startok<<-T
            if (checkokd()) tkconfigure(drawButton,state='normal')
          }
          return(T)
        }

        dgchk = function(g){
          if (length(g) != 1 || nchar(g) == 0 || is.na(g) || g<=0 || g>1) {
            tkconfigure(target_g.edit, bg=errcolor)
            tkconfigure(drawButton,state='disabled')
            target_gok<<-F
            return(F)
          } else {
            tkconfigure(target_g.edit, bg='white')
            target_gok<<-T
            if (checkokd()) {
              tkconfigure(drawButton,state='normal')
            }
          }
          return(T)
        }
        dkchk = function(k){
          if (nchar(k)==0 || is.na(k) || k<0 || k>1) {
            tkconfigure(k.edit, bg=errcolor)
            tkconfigure(drawButton,state='disabled')
            kok<<-F
            return(F)
          } else {
            tkconfigure(k.edit, bg='white')
            kok<<-T
            if (checkokd()) {
              tkconfigure(drawButton,state='normal')
            }
          }
          return(T)
        }
        dpdachk = function(pda) {          # a.edit is the X box itself
          if (nchar(pda)==0 || is.na(pda) || pda<=0) {
            tkconfigure(pda0.edit, bg=colors()[652])
            tkconfigure(drawButton,state='disabled')
            tkconfigure(pdView,state='disabled')
            pdaok<<-F
            return(F)
          } else {
            tkconfigure(pda0.edit, bg='white')
            pdaok<<-T
            if (pdaok & pdbok & Ileg){
              rCP<-rCPgab(toR(tkvars$persistence_distn),toR(tkvars$pda0), toR(tkvars$pdb0), toR(tkvars$Isam))
              tkconfigure(rCP.lbl1,text=paste0("meanCP = ",signif(rCP[1],3),", r = ",signif(rCP[2],3), " for Ir = ",toR(tkvars$Isam)))
              tkconfigure(pdView,state='normal')
            }
            if (checkokd()) tkconfigure(drawButton,state='normal')
          }
          return(T)
        }
        dpdbchk = function(pdb) {          # a.edit is the X box itself
          if (nchar(pdb)==0 || is.na(pdb) || (pdb<=0 & toR(tkvars$persistence_distn)!="Lognormal")) {
            tkconfigure(pdb0.edit, bg=colors()[652])
            tkconfigure(drawButton,state='disabled')
            tkconfigure(pdView,state='disabled')
            pdbok<<-F
            return(F)
          } else {
            tkconfigure(pdb0.edit, bg='white')
            pdbok<<-T
            if (toR(tkvars$persistence_distn)!="Exponential"){
              if (pdaok & pdbok & Ileg){
                rCP<-rCPgab(toR(tkvars$persistence_distn),toR(tkvars$pda0), toR(tkvars$pdb0), toR(tkvars$Isam))
                tkconfigure(rCP.lbl1,text=paste0("meanCP = ",signif(rCP[1],3),", r = ",signif(rCP[2],3), " for Ir = ",toR(tkvars$Isam)))
                tkconfigure(pdView,state='normal')
              }
              if (checkokd()) tkconfigure(drawButton,state='normal')
            } else {
              if (pdaok & pdbok & Ileg){
                rCP<-rCPgab(toR(tkvars$persistence_distn),toR(tkvars$pda0), toR(tkvars$pdb0), toR(tkvars$Isam))
                tkconfigure(rCP.lbl1,text=paste0("meanCP = ",signif(rCP[1],3),", r = ",signif(rCP[2],3), " for Ir = ",toR(tkvars$Isam)))
                tclvalue(tkvars$pda0)<<-1/rCP[1]
                tkconfigure(pdView,state='normal')
              }
              if (checkokd()) tkconfigure(drawButton,state='normal')
             }
          }
          return(T)
        }
        dtauchk = function(tau){
          if (length(tau) != 1 || nchar(tau) == 0 || is.na(tau) || tau <= 0) {
            tkconfigure(tau.edit, bg=errcolor)
            tkconfigure(drawButton,state='disabled')
            tauok<<-F
            return(F)
          } else {
            tkconfigure(tau.edit, bg='white')
            tauok<<-T
            if (checkokd()) {
              tkconfigure(drawButton,state='normal')
            }
          }
          return(T)
        }
        dalphachk = function(alpha){
          if (length(alpha) != 1 || nchar(alpha) == 0 || is.na(alpha) || alpha <= 0 || alpha>=1) {
            tkconfigure(crlev.edit, bg=errcolor)
            tkconfigure(drawButton,state='disabled')
            alphaok<<-F
            return(F)
          } else {
            tkconfigure(crlev.edit, bg='white')
            alphaok<<-T
            if (checkokd()) {
              tkconfigure(drawButton,state='normal')
            }
          }
          return(T)
        }
        dXchk = function(X){
          if (length(X) != 1 || nchar(X) == 0 || is.na(X) || X < 0 || round(X) != X) {
            tkconfigure(X.edit, bg=errcolor)
            tkconfigure(drawButton,state='disabled')
            Xok<<-F
            return(F)
          } else {
            tkconfigure(X.edit, bg='white')
            Xok<<-T
            if (checkokd()) {
              tkconfigure(drawButton,state='normal')
            }
          }
          return(T)
        }
        checkokd = function(){
          pd<-toR(tkvars$persistence_distn)
          (ifelse(tclvalue(tkvars$dfoption)=='g', target_gok, Xok && alphaok && tauok) & grincok & kok & tarrok & ((pd=="Exponential" & lamok) | (pd!="Exponential" & pdaok & pdbok)) & fok & phiok & schedok & startok)
        }
#      tkvars[['distr']]<<- tclVar() #text variable for field CP distribution shape and scale
      ### R variables stored in the dtForm list
      ### build the components of the form
      ## toplevel and menus
      buttwid<-6
      tcltk::tcl("option","add","*tearOff",0)
      dtModule<<-tktoplevel()
      tkgrab.set(dtModule)
      tkfocus(dtModule)
      tkwm.title(dtModule,paste0("EoA, v",.Rvar$VER," - Design Tradeoffs Module"))
      tkwm.resizable(dtModule,0,0)
      tkwm.deiconify(dtModule)

      #  tkconfigure(dtModule)
      DTtopMenu  <- tkmenu(dtModule)
      DTeditMenu <- tkmenu(DTtopMenu)
      tkadd(DTeditMenu,"command",label="Restore defaults",command=function(){
        assign("dtPrevious", dtDefault, env = .Rvar)
        suppressWarnings(file.remove(paste0(.Rvar$datadir,"/dtPrevious.Rdata")))
        tkdestroy(dtModule)
        initialize(dtDefault)
      })
      tkadd(DTeditMenu,"command",label="Restore previous",command=function(){
        tkdestroy(dtModule)
        initialize(.Rvar$dtPrevious)
      })
      tkadd(DTeditMenu,"command",label="Save to file (.rds)",command=function() {
        dtProvisional<-list()
        for (nm in names(tkvars)) dtProvisional[[nm]]<-toR(tkvars[[nm]])
        dtProvisional$persistence_distn<-toR(tkvars$persistence_distn)
        dtProvisional$arrcomponents<-arrcomponents
        saveparm(dtProvisional)
        tkwm.title(dtModule,paste0(.Rvar$dataFileTitle, " - EoA, v",.Rvar$VER," - Design Tradeoffs Module"))
      })
      tkadd(DTeditMenu,"command",label="Read from file (.rds)",command=function(){
        filename <- tclvalue(tkgetOpenFile(filetypes = "{{R data files} {.rds}}",defaultextension = ".rds", initialfile = '.rds', title = "Read"))
        tmp<-unlist(strsplit(filename,'/')); pathname<-paste(tmp[-length(tmp)],collapse='/')
        if (filename == "") return(FALSE)
        parms <- tryCatch(readRDS(filename), error=function(){tkmessageBox(icon='error',message=paste0("Error: Unable to read file:\n\n",filename)); return(F)})
        assign("csvpath", pathname, env = .Rvar)
      # check whether search schedule type, persistence distribution, prior, and arrival function are defined properly
        if (class(parms) != "try-error" && DTparmLoadable(parms)){ # if loadable, then redraw form
          tkdestroy(dtModule)
          initialize(parms) # loads the values to the form
          # error-check data...
        }
        tkwm.title(dtModule,paste0(tmp[length(tmp)], " - EoA, v",.Rvar$VER," - Design Tradeoffs Module"))
      })
      tkconfigure(DTeditMenu, activebackground=colors()[125], activeforeground=colors()[109])
      tkadd(DTtopMenu, "cascade", label="Edit", menu=DTeditMenu)
      tkconfigure(dtModule, menu=DTtopMenu)

      #### read data from R into form variables
      monitoringFrame<-ttklabelframe(dtModule, text = "     Monitoring Parmeters", padding = 10)

      objectivesFrame<-ttklabelframe(dtModule, text = "     Design Objectives", padding = 10)
      df_g<-tkradiobutton(objectivesFrame,variable = tkvars$dfoption, value='g', text = "Design for detection probability", command=function(){
        for (obj in unlist(strsplit(tclvalue(tkwinfo("children",designforaFrame)), ' '))) tkconfigure(obj,state='disabled')
        for (obj in unlist(strsplit(tclvalue(tkwinfo("children",designforgFrame)), ' '))) tkconfigure(obj,state='normal')
      })
      df_a<-tkradiobutton(objectivesFrame,variable = tkvars$dfoption, value='a', text = "Design for credibility level", command=function(){
        for (obj in unlist(strsplit(tclvalue(tkwinfo("children",designforgFrame)), ' '))) tkconfigure(obj,state='disabled')
        for (obj in unlist(strsplit(tclvalue(tkwinfo("children",designforaFrame)), ' '))) tkconfigure(obj,state='normal')
      })
      designforgFrame<-tkframe(objectivesFrame)
      target_g.lbl<-tklabel(designforgFrame,text="Target overall detection probability (g): ")
      target_g.edit<-tkentry(designforgFrame, width=5, textvariable=tkvars$target_g,justify='right',bg='white')
      MvgButton<-tkbutton(designforgFrame, text = "Explore M* vs. g", command = dtmvg)
      tkgrid(target_g.lbl, target_g.edit, sticky='w')
      tkgrid(MvgButton, row = 0, column = 2, padx = c(15,0))
      designforaFrame<-tkframe(objectivesFrame)
      tau.lbl<-tklabel(designforaFrame, text = "Threshold (\u03c4) ")
      tau.edit<-tkentry(designforaFrame, width = 5, textvariable=tkvars$tau, justify= 'right', bg= 'white')
      crlev.lbl<-tklabel(designforaFrame, text = "1 - \u03b1 ")
      crlev.edit<-tkentry(designforaFrame, width = 5, textvariable=tkvars$crlev, justify= 'right', bg= 'white')
      X.lbl<-tklabel(designforaFrame, text ='Carcass count (X) ')
      X.edit<-tkentry(designforaFrame, width = 5, textvariable = tkvars$X, justify = 'right', bg='white')
      tkgrid(tau.lbl, tau.edit)
      tkgrid(crlev.lbl, row = 0, column = 2, padx = c(10, 0))
      tkgrid(crlev.edit, row = 0, column = 3, padx = c(0,10))
      tkgrid(X.lbl, row = 0, column = 4)
      tkgrid(X.edit, row= 0, column =5)
      tkgrid(df_g, sticky='w')
      tkgrid(designforgFrame, sticky='w', padx = 15, pady = c(0,10))
      tkgrid(df_a, sticky='w')
      tkgrid(designforaFrame, sticky='w', padx = 15, pady = c(0,10))
      if (tclvalue(tkvars$dfoption)=='g'){
        for (obj in unlist(strsplit(tclvalue(tkwinfo("children",designforaFrame)), ' '))) tkconfigure(obj,state='disabled')
        for (obj in unlist(strsplit(tclvalue(tkwinfo("children",designforgFrame)), ' '))) tkconfigure(obj,state='normal')
      } else {
        for (obj in unlist(strsplit(tclvalue(tkwinfo("children",designforgFrame)), ' '))) tkconfigure(obj,state='disabled')
        for (obj in unlist(strsplit(tclvalue(tkwinfo("children",designforaFrame)), ' '))) tkconfigure(obj,state='normal')
      }
      SE.lbl<-tklabel(monitoringFrame,text="Searcher efficiency (p)")
      fdatFrame<-tkframe(monitoringFrame)
      f.lbl<-tklabel(fdatFrame,text="p: "); f.edit<-tkentry(fdatFrame,width=5,textvariable=tkvars$f,justify='right',bg='white')
      fmin.lbl<-tklabel(fdatFrame,text="min(p): "); fmin.edit<-tkentry(fdatFrame,width=5,textvariable=tkvars$fmin,justify='right',bg='white')
      fmax.lbl<-tklabel(fdatFrame,text="max(p): "); fmax.edit<-tkentry(fdatFrame,width=5,textvariable=tkvars$fmax,justify='right',bg='white')
      PHI.lbl<-tklabel(monitoringFrame,text="Spatial coverage (a)")
      phidatFrame<-tkframe(monitoringFrame)
      phi.lbl<-tklabel(phidatFrame,text="a: "); phi.edit<-tkentry(phidatFrame,width=5,textvariable=tkvars$phi,justify='right',bg='white')
      phimin.lbl<-tklabel(phidatFrame,text="min(a): "); phimin.edit<-tkentry(phidatFrame,width=5,textvariable=tkvars$phimin,justify='right',bg='white')
      phimax.lbl<-tklabel(phidatFrame,text="max(a): "); phimax.edit<-tkentry(phidatFrame,width=5,textvariable=tkvars$phimax,justify='right',bg='white')
      SI.lbl<-tklabel(monitoringFrame,text="Search interval (I)")
      IdatFrame<-tkframe(monitoringFrame)
      I.lbl<-tklabel(IdatFrame,text="I: "); Isam.edit<-tkentry(IdatFrame,width=5,textvariable=tkvars$Isam,justify='right',bg='white')
      Imin.lbl<-tklabel(IdatFrame,text="min(I): "); Imin.edit<-tkentry(IdatFrame,width=5,textvariable=tkvars$Imin,justify='right',bg='white')
      Imax.lbl<-tklabel(IdatFrame,text="max(I): "); Imax.edit<-tkentry(IdatFrame,width=5,textvariable=tkvars$Imax,justify='right',bg='white')
      span.lbl<-tklabel(IdatFrame,text="Span: "); span.edit<-tkentry(IdatFrame,width=5,textvariable=tkvars$span,justify='right',bg='white')
      fixed.lbl<-tklabel(monitoringFrame,text="Fixed"); vary.lbl<-tklabel(monitoringFrame,text="Variable")
      tkvars$ffix<-tclVar(dtdat$ffix)
      ffixedRadio<-tkradiobutton(monitoringFrame,variable=tkvars$ffix,value="1")
      fvaryRadio<-tkradiobutton(monitoringFrame,variable=tkvars$ffix,value="0")
      tkvars$phifix<-tclVar(dtdat$phifix)
      phifixedRadio<-tkradiobutton(monitoringFrame,variable=tkvars$phifix,value="1")
      phivaryRadio<-tkradiobutton(monitoringFrame,variable=tkvars$phifix,value="0")
      tkvars$Ifix<-tclVar(dtdat$Ifix)
      IfixedRadio<-tkradiobutton(monitoringFrame,variable=tkvars$Ifix,value="1")
      IvaryRadio<-tkradiobutton(monitoringFrame,variable=tkvars$Ifix,value="0")

      tkgrid(fixed.lbl,row=1,column=1,sticky='e');  tkgrid(vary.lbl,row=1,column=2)
      tkgrid(SE.lbl,sticky='w');  tkgrid(ffixedRadio, row=2, column=1, sticky='e'); tkgrid(fvaryRadio,row=2,column=2)
      tkgrid(f.lbl,f.edit,fmin.lbl,fmin.edit,fmax.lbl,fmax.edit)
      if (toR(tkvars$ffix)==1){
        tkconfigure(f.edit,state='normal')
        tkconfigure(fmin.edit,state='disabled')
        tkconfigure(fmax.edit,state='disabled')
      } else {
        tkconfigure(f.edit,state='disabled')
        tkconfigure(fmin.edit,state='normal')
        tkconfigure(fmax.edit,state='normal')
      }
      tkgrid(fdatFrame,row=2,column=3,sticky='w')

      tkgrid(PHI.lbl,sticky='w');  tkgrid(phifixedRadio,row=3,column=1, sticky='e'); tkgrid(phivaryRadio,row=3,column=2)
      tkgrid(phi.lbl,phi.edit,phimin.lbl,phimin.edit,phimax.lbl,phimax.edit)
      if (toR(tkvars$phifix)==1){
        tkconfigure(phi.edit,state='normal')
        tkconfigure(phimin.edit,state='disabled')
        tkconfigure(phimax.edit,state='disabled')
      } else {
        tkconfigure(phi.edit,state='disabled')
        tkconfigure(phimin.edit,state='normal')
        tkconfigure(phimax.edit,state='normal')
      }
      tkgrid(phidatFrame,row=3,column=3,sticky='w')

      start.lbl<-tklabel(IdatFrame,text="Date of first search\n   (yyyy-mm-dd)",anchor='w')
      firstsearch.edit<-tkentry(IdatFrame,width=10,textvariable=tkvars$firstsearch,justify='left',bg='white')
      tkgrid(start.lbl, columnspan = 3, sticky='w')
      tkgrid(firstsearch.edit, column=3, row=0)
      tkgrid(span.lbl, row=0, column=4)
      tkgrid(span.edit, row=0, column =5)
      tkgrid(SI.lbl,sticky='w');  tkgrid(IfixedRadio,row=4, column=1, sticky='e'); tkgrid(IvaryRadio,row=4,column=2)
      tkgrid(I.lbl,Isam.edit,Imin.lbl,Imin.edit,Imax.lbl,Imax.edit)
      if (toR(tkvars$Ifix)==1){
        tkconfigure(Isam.edit,state='normal')
        tkconfigure(Imin.edit,state='disabled')
        tkconfigure(Imax.edit,state='disabled')
      } else {
        tkconfigure(Isam.edit,state='disabled')
        tkconfigure(Imin.edit,state='normal')
        tkconfigure(Imax.edit,state='normal')
      }
      tkgrid(IdatFrame,row=4,column=3,sticky='w',pady=10)
      k.lbl1<-tklabel(monitoringFrame,text='Factor by which searcher efficiency\n      changes with each search (k):',justify='left')
      k.edit<-tkentry(monitoringFrame,width=6,textvariable=tkvars$k,justify='right',bg='white')
      tkgrid(k.lbl1,k.edit, padx=5, pady=5, columnspan=3, sticky='w')

      tkgrid(objectivesFrame, sticky='nw', pady=10, padx=c(5,0))
#      tkgrid(monitoringFrame, columnspan=2, sticky='w')
      tkgrid(monitoringFrame, column=1, row=0, rowspan = 2, sticky='nw', pady=10, padx=c(0,5))


      ## frame for background parameters
      bvarFrame<-ttklabelframe(dtModule,text="     Background Parameters",padding=10)
      persistenceFrame<-ttklabelframe(bvarFrame,text="   Persistence Distribution",padding=10)
      pdparmFrame<-ttklabelframe(persistenceFrame,text='Parameters')
      shapescaleFrame<-tkframe(pdparmFrame)
      pda.lbl<-tklabel(shapescaleFrame,text="shape (\u03b1) "); pda0.edit<-tkentry(pdparmFrame,width=6,textvariable=tkvars$pda0,justify='right',bg='white')
      pdb.lbl<-tklabel(shapescaleFrame,text="scale (\u03b2) "); pdb0.edit<-tkentry(pdparmFrame,width=6,textvariable=tkvars$pdb0,justify='right',bg='white')
      if (pdaok & pdbok){
        rCP<-rCPgab(toR(tkvars$persistence_distn),toR(tkvars$pda0),toR(tkvars$pdb0),toR(tkvars$Isam))
        rCP.lbl1<-tklabel(pdparmFrame,text=paste0("meanCP = ",signif(rCP[1],3),", r = ",signif(rCP[2],3), " for Ir = ",toR(tkvars$Isam)),width=30)
      } else {
        rCP.lbl1<-tklabel(pdparmFrame,text=paste0("meanCP = ",NA,", r = ",NA, " for Ir = ",NA),width=30)
      }
      tkgrid(pda.lbl,pda0.edit)
      tkgrid(pdb.lbl,pdb0.edit)
      tkgrid(shapescaleFrame)
      tkgrid(rCP.lbl1)
      pdlboxFrame<-tkframe(persistenceFrame)
      pdView<-tkbutton(pdlboxFrame,text="View",command=function() plotPersDist(
        toR(tkvars$persistence_distn),
        toR(tkvars$Isam),
        toR(tkvars$pda0),
        toR(tkvars$pdb0)
      ))

      pd.lbox<-tklistbox(pdlboxFrame,height=4,width=12,selectmode="single",exportselection="0",bg='white')
      pdnames<-c("Exponential","Weibull","Log-Logistic","Lognormal")
      for (i in 1:length(pdnames)) tkinsert(pd.lbox,"end",pdnames[i])
      tkselection.set(pd.lbox,which(pdnames==toR(tkvars$persistence_distn))-1)
      tkactivate(pd.lbox,which(pdnames==toR(tkvars$persistence_distn))-1)
      tkgrid(pd.lbox); tkgrid(pdView)
      tkgrid(pdlboxFrame,pdparmFrame)

      parent<-bvarFrame
      parent$name<-"dt"
      tkvars$samtype<-tclVar("dt")
      if (!basicMode){
        arrfFrame<-ttklabelframe(bvarFrame,text="Arrival Function")
        arrUniformRadio<-tkradiobutton(arrfFrame, variable=tkvars$arrfun, text="Uniform",  value="Uniform")
        arrCompoundRadio<-tkradiobutton(arrfFrame,variable=tkvars$arrfun, text="Compound", value="Compound")
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
          dttmp<-
          arrivalFun(self)
          tkrplot::tkrreplot(arrfig.mini)
        })

        arrfig.mini <<- tkrplot::tkrplot(arrfFrame, fun = plotarr.mini, hscale=.25, vscale=.08)

        # construct arrival frame
        tkbind(arrUniformRadio,"<Button-1>", function() {
          if (tclvalue(tkcget(arrUniformRadio,"-state"))=="disabled") return(F)
          tclvalue(tkvars$arrfun)<-"Uniform"
          tkrplot::tkrreplot(arrfig.mini)
          tkconfigure(editArrivalsbutton,state='disabled')
        })
        tkbind(arrCompoundRadio,"<Button-1>", function() {
          if (tclvalue(tkcget(arrUniformRadio,"-state"))=="disabled") return(F)
          tclvalue(tkvars$arrfun)<-"Compound"
          tkrplot::tkrreplot(arrfig.mini)
          tkconfigure(editArrivalsbutton,state='normal')
        })
        tkgrid(arrUniformRadio,sticky='w')
        tkgrid(arrCompoundRadio,sticky='w')
        tkgrid(arrfig.mini, column=1,row=0,rowspan=2,padx=10)
        tkgrid(editArrivalsbutton,column=2,row=0,rowspan=2, padx=10)
      } else {
        tarrFrame<-tkframe(monitoringFrame)
        tarr.lbl<-tklabel(tarrFrame, text = "Temporal coverage (v) ")
        tarr.edit<-tkentry(tarrFrame, width=5, textvariable=tkvars$tarr, justify='right', bg='white')
        tkgrid(tarr.lbl, tarr.edit)
        tkgrid(tarrFrame, row=5, column =3, sticky='e')
      }
      tkgrid(persistenceFrame,sticky='nw',columnspan=2, padx=10)
      if (!basicMode) tkgrid(arrfFrame,sticky='nw', row = 1, column=0, padx=10, pady=10)
      tkgrid(bvarFrame, sticky='nw', row = 1, pady=c(10,0), padx=5)#, row = 1, column = 12)

      #  butFrame<-tkframe(dtModule)
      wid<-8
      buttonFrame<-tkframe(dtModule)
      drawButton<-tkbutton(buttonFrame,text="Calculate", command=function() { # no error-check?
        for (nm in names(dtDefault)) .Rvar$dtPrevious[[nm]] <- toR(tkvars[[nm]])
        drawdt(.Rvar$dtPrevious)
      }, width=wid)
      closeButton<-tkbutton(buttonFrame,text="Close",command=function() {
        graphics.off()
        dtsave(.Rvar$dtPrevious)
        tkdestroy(dtModule)
        tkwm.deiconify(.Rvar$EoA)
      }, width=wid)
      tkwm.protocol(dtModule, "WM_DELETE_WINDOW", function(){
        tkdestroy(dtModule)
        dtsave(.Rvar$dtPrevious)
        tkwm.deiconify(.Rvar$EoA)
      })

      #cancelButton<-tkbutton(buttonFrame,text="Cancel",command=canceldt, width=wid)
      tkgrid(drawButton, closeButton)
      tkgrid(buttonFrame, row = 1, column = 1, sticky='e', padx = c(10,20))
      # error-checking
      # bind radio buttons
      tkbind(phifixedRadio,"<Button-1>", function() {
        tkconfigure(phi.edit,state='normal')
        tkconfigure(phimin.edit,state='disabled')
        tkconfigure(phimax.edit,state='disabled')
        phichk(toR(tkvars$phi))
      })
      tkbind(phivaryRadio,"<Button-1>", function() {
        tkconfigure(phi.edit,state='disabled')
        tkconfigure(phimin.edit,state='normal')
        tkconfigure(phimax.edit,state='normal')
        phiminchk(toR(tkvars$phimin))
        phimaxchk(toR(tkvars$phimax))
      })
      tkbind(ffixedRadio,"<Button-1>", function() {
        tkconfigure(f.edit,state='normal')
        tkconfigure(fmin.edit,state='disabled')
        tkconfigure(fmax.edit,state='disabled')
      })
      tkbind(fvaryRadio,"<Button-1>", function() {
        tkconfigure(f.edit,state='disabled')
        tkconfigure(fmin.edit,state='normal')
        tkconfigure(fmax.edit,state='normal')
      })
      tkbind(IfixedRadio,"<Button-1>", function() {
        tkconfigure(Isam.edit,state='normal')
        tkconfigure(Imin.edit,state='disabled')
        tkconfigure(Imax.edit,state='disabled')
#        Ifix<<-T
      #  if (dIchk(toR(tkvars$Isam))){
          dschedchk()
       # }
      })
      tkbind(IvaryRadio,"<Button-1>", function() {
        tkconfigure(Isam.edit,state='disabled')
        tkconfigure(Imin.edit,state='normal')
        tkconfigure(Imax.edit,state='normal')
#        Ifix<<-F
        dschedchk()
      })
      #tkbind(arrUniformRadio,"<Button-1>", function() {dtProvisional$arrfun<<-"Uniform"; dupdatearrrad("Uniform")})
      #tkbind(arrCompoundRadio,"<Button-1>", function()  {dtProvisional$arrfun<<-"Compound"; dupdatearrrad("Compound")})
      if (!basicMode){
        tkbind(arrUniformRadio,"<Button-1>", function() {
          tclvalue(tkvars$arrfun)<-"Uniform"
          tkrplot::tkrreplot(arrfig.mini)
          tkconfigure(editArrivalsbutton,state="disabled")
        })
        tkbind(arrCompoundRadio,"<Button-1>", function() {
          tclvalue(tkvars$arrfun)<-"Compound"
          tkrplot::tkrreplot(arrfig.mini)
          tkconfigure(editArrivalsbutton,state="normal")
        })
      }
      # bind data entry textboxes
      tkbind(tarr.edit,"<KeyRelease>", function() {tarrchk(suppressWarnings(as.numeric(toR(tkvars$tarr))))})
      tkbind(phi.edit,"<KeyRelease>", function() {phichk(suppressWarnings(as.numeric(toR(tkvars$phi))))})
      tkbind(phimin.edit,"<KeyRelease>", function() {phiminchk(suppressWarnings(as.numeric(toR(tkvars$phimin))))})
      tkbind(phimax.edit,"<KeyRelease>", function() {phimaxchk(suppressWarnings(as.numeric(toR(tkvars$phimax))))})
      tkbind(f.edit,"<KeyRelease>", function() {fchk(suppressWarnings(as.numeric(toR(tkvars$f))))})
      tkbind(fmin.edit,"<KeyRelease>", function() {fminchk(suppressWarnings(as.numeric(toR(tkvars$fmin))))})
      tkbind(fmax.edit,"<KeyRelease>", function() {fmaxchk(suppressWarnings(as.numeric(toR(tkvars$fmax))))})
      tkbind(Isam.edit,"<KeyRelease>", function() {Ifix<<-toR(tkvars$Ifix); dIchk(suppressWarnings(as.numeric(toR(tkvars$Isam))))})
      tkbind(Imin.edit,"<KeyRelease>", function() {Ifix<<-toR(tkvars$Ifix); dIminchk(suppressWarnings(as.numeric(toR(tkvars$Imin))))})
      tkbind(Imax.edit,"<KeyRelease>", function() {Ifix<<-toR(tkvars$Ifix); dImaxchk(suppressWarnings(as.numeric(toR(tkvars$Imax))))})
      tkbind(span.edit,"<KeyRelease>", function() {Ifix<<-toR(tkvars$Ifix); dschedchk()})
      tkbind(firstsearch.edit,"<KeyRelease>", function() {firstsearch<-tclvalue(tkvars$firstsearch); startchk(firstsearch)})
      tkbind(target_g.edit,"<KeyRelease>", function() {dgchk(suppressWarnings(as.numeric(toR(tkvars$target_g))))})
      tkbind(k.edit,"<KeyRelease>", function() {dkchk(suppressWarnings(as.numeric(toR(tkvars$k))))})
      tkbind(pda0.edit,"<KeyRelease>", function() {dpdachk(suppressWarnings(as.numeric(toR(tkvars$pda0))))})
      tkbind(pdb0.edit,"<KeyRelease>", function() {dpdbchk(suppressWarnings(as.numeric(toR(tkvars$pdb0))))})
      tkbind(X.edit, "<KeyRelease>", function(){dXchk(suppressWarnings(as.numeric(toR(tkvars$X))))})
      tkbind(tau.edit, "<KeyRelease>", function(){dtauchk(suppressWarnings(as.numeric(toR(tkvars$tau))))})
      tkbind(crlev.edit, "<KeyRelease>", function(){dalphachk(suppressWarnings(as.numeric(toR(tkvars$crlev))))})
       # note: chk names for I vars get a 'd' prefix to distinguish them from check functions for singleYear module
      tkbind(pd.lbox,"<<ListboxSelect>>", function() {
        v<-tclvalue(tkget(pd.lbox,tkcurselection(pd.lbox)))
        tclvalue(tkvars$persistence_distn)<<-v
        if (v=='Exponential') {
          tkconfigure(pdb.lbl,text="meanCP ")
          tkconfigure(pda.lbl,text='rate');
          tclvalue(tkvars$pda0)<-signif(1/toR(tkvars$pdb.e),5); tkconfigure(pda0.edit,state='disabled')
          tclvalue(tkvars$pdb0)<-toR(tkvars$pdb.e)
        } else if (v=='Weibull' |  v=='Log-Logistic' | v=='Lognormal') {
          if (v=='Weibull'){
            tclvalue(tkvars$pda0)<-tclvalue(tkvars$pda.w)
            tclvalue(tkvars$pdb0)<-tclvalue(tkvars$pdb.w)
          } else if (v=='Log-Logistic') {
            tclvalue(tkvars$pda0)<-tclvalue(tkvars$pda.ll)
            tclvalue(tkvars$pdb0)<-tclvalue(tkvars$pdb.ll)
          } else if (v=='Lognormal') {
            tclvalue(tkvars$pda0)<-tclvalue(tkvars$pda.ln)
            tclvalue(tkvars$pdb0)<-tclvalue(tkvars$pdb.ln)
          }
          tkconfigure(pda.lbl,text='shape (\u03b1) '); tkconfigure(pda0.edit,state='normal')
          tkconfigure(pdb.lbl,text='scale (\u03b2) ')
        }
          pda<-as.numeric(tclvalue(tkvars$pda0))
          pdb<-as.numeric(tclvalue(tkvars$pdb0))
          Ir<-as.numeric(tclvalue(tkvars$Isam))
          rCP<-rCPgab(v, pda, pdb, Ir)
          meanCP<<-rCP[1]; rhat<-rCP[2]
          tkconfigure(rCP.lbl1,text=paste0("meanCP = ",signif(rCP[1],4),", r = ",signif(rCP[2],4), " for Ir = ", Ir))
          # when listbox selection is changed, then the parameter set that is loaded is proper, so...
          #  -- all the backgrounds are whited
          #  -- the OK and leg flags are set to true
          #  -- the "view" button is set to true
          #  -- the full suite of flags is checked for possible activation of the 'calculate' buttons
          tkconfigure(pda0.edit, bg='white'); pdaok<<-T
          tkconfigure(pdb0.edit, bg='white'); pdbok<<-T
          tkconfigure(pdView,state='normal')
          if (checkokd()) tkconfigure(drawButton,state='normal')
      })

        # destructor
      #  tkbind(dtModule,"<Destroy>",command=function() {tkgrid(welcomeFrame); tkmessageBox(message='xxx')})
    },
    plotarr.mini = function() { # not currently scaled to endpoints
      par(mar=c(0,0,0,0))
      if (tclvalue(tkvars$arrfun)=="Uniform"){
        plot(0,0,type='n',xlim=c(0,364),ylim=c(0,1),yaxs='i',xaxs='i')
        x0<-as.numeric(as.Date(tclvalue(tkvars$firstsearch))-as.Date(paste0(format(as.Date(tclvalue(tkvars$firstsearch)),"%Y"),"-01-01")))
        span<-toR(tkvars$span)
        x1<-x0 + span
        lht<-.3
        polygon(c(x0,x0,x1,x1),c(0,lht, lht, 0))
      } else {
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
  )
)
dtsave<-function(dt){
  save(dt, file = paste0(.Rvar$datadir,"/dtPrevious.Rdata"))
}
