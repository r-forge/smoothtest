#rZIP<-function(n,p,lambda) {
#  if(p>=0) {
#    n0<-rbinom(1,size=n,prob=p)
#    n1<-n-n0
#    r<-c(rep(0,n0),rpois(n1,lambda=lambda))
#  }
#  else {
#    n0<-rbinom(1,size=n,prob=dZIP(0,p=p,lambda=lambda))
#    if(is.na(n0)) n0<-0 # NA may be produced by a too
#    # negative p caused by the estimation
#    r<-rep(0,n0)
#    while(length(r)<n) {
#      s<-rpois(1,lambda=lambda)
#      if(s!=0) r<-c(r,s)
#    }
#  }
#  return(r) 
#}

#rlaplace<-function(n,location=0,scale=1) { 
#  U<-runif(n)
#  T<-rep(0,n)
#  optfn.lapl<-function(x,y) {
#    (cdflaplace(x,location,scale)-y)^2
#  }
#  for (i in 1:n) {
#    tmp<-optim(location,optfn.lapl,y=U[i])
#    T[i]<-tmp$par
#  }
#  return(T)
#}

#cdflaplace<-function(x,location=0,scale=1) {
#  (1+sign(x-location)*(1-exp(-abs(x-location)/scale)))/2
#}

rZTP<-function(n,lambda) {
 qpois(runif(n, dpois(0, lambda), 1), lambda) 
}

rsamplecont<-function(n,f) { 
  U<-runif(n)
  T<-rep(0,n)
  cdf<-function(x){
    temp<-try(integrate(f,-Inf,x)$value,silent=TRUE)
    if (inherits(temp,"try-error")) {temp<-100}
    return(temp)
  }
  optfn<-function(x,y) {
    (cdf(x)-y)
  }
  for (i in 1:n) {
    tmp<-nls.lm(0,optfn,y=U[i])
    T[i]<-tmp$par
  }
  return(T)
}

#rsamplecont<-function(n,f) { 
#  U<-runif(n)
#  T<-rep(0,n)
#  cdf<-function(x){
#    integrate(f,-Inf,x)$value
#  }
#  optfn<-function(x,y) {
#    (cdf(x)-y)^2
#  }
#  for (i in 1:n) {
#    tmp<-optim(0,optfn,y=U[i])
#    T[i]<-tmp$par
#  }
#  return(T)
#}

# het genereren van discrete random samples gaat wel zeer traag
#rlogar<-function(n,Beta=1/2) { 
#  U<-runif(n)
#  T<-rep(0,n)
#  g<-function(x){ 
#    y<-cdflogar(x,Beta)
#    return(y)
#  }
#  tmp<-1:1000
#  K<-rep(FALSE,1000)
#  for (i in 1:n) { 
#    for (j in tmp) { 
#       K[j]<-g(j)>=U[i]
#    }
#    T[i]<-min(tmp[K])
#  }
#  return(T)
#}
  
#cdflogar<-function(x,Beta=1/2) {
#  t<-1:x
#  sum(-Beta^t/(t*log(1-Beta)))  
#}

# of identiek
#cdflogar2<-function(x,Beta=1/2) {
#  p<-Beta
#  1+beta_x(p,x+1,0)/log(1-p) 
#}

#beta_x<-function(x,a,b){ 
# g<-function(t) {
#    t^(a-1)*(1-t)^(b-1)
# }
# integrate(g,0,x)$value
#}

#rsampledisc<-function(n,f) { 
#  U<-runif(n)
#  T<-rep(0,n)
#  cdf<-function(x){
#    t<-0:x
#    sum(f(t))
#  }
#  tmp<-0:1000
#  K<-rep(FALSE,1001)
#  for (i in 1:n) { 
#    for (j in tmp) { 
#       K[j+1]<-cdf(j)>=U[i]
#    }
#    T[i]<-min(tmp[K])
#  }
#  return(T)
#}

#veel snellere versie

rsampledisc<-function(n,f) { 
  U<-runif(n)
  T<-rep(0,n)
  U<-sort(U)
  quant<-0
  cdf<-f(quant)
  for (i in 1:n) {
    while (cdf<U[i]) {
	quant<-quant+1
	cdf<-cdf+f(quant) #cdf is on every moment cdf(quant)
    }
    T[i]<-quant
  }
  return(T)
}

#nog iets sneller, maar minder goed wegens beperking 0:1000 (bvb)
#rsampledisc<-function(n,f) {
#	x<-0:1000
#	fx<-f(x) 
#	rsample<-sample(x,size=n,replace=TRUE,prob=fx)
#	return(rsample)
#}
