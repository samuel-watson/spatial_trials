## reanalysis functions

fn <- function(x,l,kappa,nu,del_e, del_i, b){
  b*((1-((sign(x)/l)*log(exp(sign(x)*l*(x + del_i)/(del_e + del_i)) + exp(l*(sign(x)+1)/2)))^kappa)^nu) 
}

fn0 <- function(x,int,b,del_e, del_i){
  int + b*((1-((sign(x)/-50)*log(exp(sign(x)*-50*(x + del_i)/(del_e + del_i)) + exp(-50*(sign(x)+1)/2)))^4)^8) 
}

genrep0 <- function(dfanal,f1,L){
  dfanal$ystar <- f1 + L%*%(rnorm(length(f1))) 
  fitn <- tryCatch(nls(ystar ~ fn0(distance, int,b, del_e, del_i),data = dfanal, 
                       start = list(int = 0, b = -0.2, del_e = 0.2, del_i = 0.1),
                       lower = c(-10,-10,0.01,0.0), upper = c(10,10,3.0,1.5), algorithm = "port"), error= function(i)return(NA))
  if(is(fitn,"nls")){
    np <- fitn$m$getPars()
  } else {
    np <- rep(NA, 3)
  }
  return(np)
}

fn2a <- function(x,d1,d2,d3,d4,int,b,del_e, del_i, b_d1, b_d2,b_d3,b_d4){
  int + d1*b_d1 + d2*b_d2 + d3*b_d3 + d4*b_d4 + b*((1-((sign(x)/-50)*log(exp(sign(x)*-50*(x + del_i)/(del_e + del_i)) + exp(-50*(sign(x)+1)/2)))^4)^8) 
}

genrep2a <- function(dfanal,f1,L){
  dfanal$ystar <- f1 + L%*%(rnorm(length(f1))) 
  fitn <- tryCatch(nls(ystar ~ fn2a(distance,distance_all1,distance_all2,distance_all3,distance_all4, int,b, del_e, del_i,b_1, b_2,b_3,b_4),data = dfanal, 
                       start = list(int = 0, b = -0.2, del_e = 0.2, del_i = 0.1, b_1 = 0, b_2 = 0,b_3 = 0, b_4 = 0),
                       lower = c(-10,-10,0.01,0.0,rep(-10,4)), upper = c(10,10,4.0,1.5,rep(10,4)), algorithm = "port"), error= function(i)return(NA))
  if(is(fitn,"nls")){
    np <- fitn$m$getPars()
  } else {
    np <- rep(NA, 8)
  }
  return(np)
}

fn2b <- function(x,d1,d2,int,b,del_e, del_i, b_d1, b_d2){
  int + d1*b_d1 + d2*b_d2 + b*((1-((sign(x)/-50)*log(exp(sign(x)*-50*(x + del_i)/(del_e + del_i)) + exp(-50*(sign(x)+1)/2)))^4)^8) 
}

genrep2b <- function(dfanal,f1,L){
  dfanal$ystar <- f1 + L%*%(rnorm(length(f1))) 
  fitn <- tryCatch(nls(ystar ~ fn2b(distance,distance_all,distance_all_sq,int,b, del_e, del_i,b_1, b_2),data = dfanal, 
                       start = list(int = 0, b = -0.2, del_e = 0.2, del_i = 0.1, b_1 = 0, b_2 = 0),
                       lower = c(-10,-10,0.01,0.0,rep(-10,2)), upper = c(10,10,1.0,0.5,rep(10,2)), algorithm = "port"), error= function(i)return(NA))
  if(is(fitn,"nls")){
    np <- fitn$m$getPars()
  } else {
    np <- rep(NA, 6)
  }
  return(np)
}

genrep2c <- function(dfanal,f1,L){
  dfanal$ystar <- f1 + L%*%(rnorm(length(f1))) 
  fitn <- tryCatch(nls(ystar ~ fn2b(distance,distance_all,distance_all_sq,int,b, del_e, 0,b_1, b_2),data = dfanal, 
                       start = list(int = 0, b = -0.2, del_e = 0.2,  b_1 = 0, b_2 = 0),
                       lower = c(-10,-10,0.01,rep(-10,2)), upper = c(10,10,1.0,rep(10,2)), algorithm = "port"), error= function(i)return(NA))
  if(is(fitn,"nls")){
    np <- fitn$m$getPars()
  } else {
    np <- rep(NA, 5)
  }
  return(np)
}
