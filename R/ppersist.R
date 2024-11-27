#' @export
ppersist <- function(
    persistence_distn,
    t_arrive0, 
    t_arrive1, 
    t_search, 
    pdb, 
    pda = NULL
){
  
  if (persistence_distn == "Weibull") {
    return (t((pgamma(outer(1/pdb, t_search-t_arrive0)^pda, 1/pda)-pgamma(outer(1/pdb, t_search-t_arrive1)^pda, 1/pda))*gamma(1+1/pda)*outer(pdb, 1/(t_arrive1-t_arrive0))))
  }
  if (persistence_distn == "Exponential") {
    return ((exp(outer(t_arrive1-t_search,1/pdb))-exp(outer(t_arrive0-t_search,1/pdb)))/(outer(t_arrive1-t_arrive0, 1/pdb)))
  }
  if (persistence_distn == "Lognormal") {
    #     part1<-pnorm(-outer(log(t_search-t_arrive1),-pdb,FUN="+")/sqrt(pda))*(t_search-t_arrive1)+pnorm(outer(log(t_search-t_arrive1),-pdb,FUN="+")/sqrt(pda)-sqrt(pda))*rep(exp(0.5*pda+pdb),length(t_arrive0))
    #     part0<-pnorm(-outer(log(t_search-t_arrive0),-pdb,FUN="+")/sqrt(pda))*(t_search-t_arrive0)+pnorm(outer(log(t_search-t_arrive0),-pdb,FUN="+")/sqrt(pda)-sqrt(pda))*rep(exp(0.5*pda+pdb),length(t_arrive0))
    root_pda <- sqrt(pda)
    exp_value <- exp((pda/2) + pdb)
    tt <- t_search-t_arrive0
    part0 <- t(pnorm(outer(pdb, -log(tt), "+")/root_pda))*tt+t(pnorm(outer(-pdb, log(tt), "+")/root_pda-root_pda)*exp_value)
    tt <- t_search-t_arrive1
    part1 <- t(pnorm(outer(pdb, -log(tt), "+")/root_pda))*tt+t(pnorm(outer(-pdb, log(tt), "+")/root_pda-root_pda)*exp_value)
    return( -(part1 - part0) / (t_arrive1 - t_arrive0) )
  }
  if (persistence_distn == "Log-Logistic") {
    return(
      matrix(
        Vectorize(
          function(
            t_arrive0, 
            t_arrive1, 
            t_search, 
            pda, 
            pdb
          ){
            
            t1 <- t_search - t_arrive1
            
            t0 <- t_search - t_arrive0
            
            part1 <- ifelse(
              t1 == 0, 0, t1 / (1 + (t1 / pdb)^pda) * gsl::hyperg_2F1(
                1, 1, 1 + 1 / pda, 1/ (1 + (t1 / pdb)^(-pda))
              )
            )
            
            part0 <- t0 / (1 + (t0 / pdb)^pda) * gsl::hyperg_2F1(
              1, 1, 1 + 1 / pda, 1/(1 + (t0 / pdb)^(-pda))
            ) 
            
            return( -(part1 - part0) / (t_arrive1 - t_arrive0) )
          },
      vectorize.args = c("pdb", "pda"))(t_arrive0, t_arrive1, t_search, pda, pdb),
      ncol = length(pda)
      )
    )
  }
}
