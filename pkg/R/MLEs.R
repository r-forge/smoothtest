unif.MLE<-function(sample){
  MIN<-min(sample)
  MAX<-max(sample)
  return(c(MIN,MAX))
}

unif.MME<-function(sample){
  m<-mean(sample)
  s<-sum((sample-m)^2)/length(sample)
  pars<-c(m-sqrt(3*s),m+sqrt(3*s))
  return(pars)
}

norm.MLE<-function(sample){
  MEAN<-mean(sample)
  SD<-sqrt(sum((sample-MEAN)^2)/length(sample))
  return(c(MEAN,SD))
}

exp.MLE<-function(sample){
  LAMBDA<-1/mean(sample)
  return(LAMBDA)
}

logis.MLE<-function(sample){
 n<-length(sample)
 tmp<-nls.lm(c(mean(sample),sqrt(3)*sqrt(sum((sample-mean(sample))^2)/n)/pi),logis.optfn,x=sample,n=n)
 pars<-tmp$par
 return(pars)
}

logis.optfn<-function(pars,x,n){
 t2<-(pars[2]-(1/n)*sum((x-pars[1])*(1-exp(-(x-pars[1])/pars[2]))/(1+exp(-(x-pars[1])/pars[2]))))
 t1<- (n/2-sum(exp(-(x-pars[1])/pars[2])/(1+exp(-(x-pars[1])/pars[2]))))
 return(c(t1,t2))
}

logis.MME<-function(sample){
  MU<-mean(sample)
  SIGMA<-sqrt(3)*sqrt(sum((sample-MU)^2)/length(sample))/pi
  return(c(MU,SIGMA))
}

ZIP.MLE<-function(sample) {
    x.bar<-mean(sample)
    n<-length(sample)
    n0<-length(sample[sample==0])
    tmp<-optim(c(n0/n,x.bar),ZIP.optfn,x.bar=x.bar,n0=n0,n=n)
    pars<-tmp$par
    # names(pars)<-c("p","lambda")
    pars<-pars[2:1]
    return(pars)
    # returns c(lambda,p)
}

ZIP.optfn<-function(pars,x.bar,n,n0) {
  # pars[1]: p ; pars[2]=lambda
  (1-x.bar/pars[2]-pars[1])^2+(pars[2]/(1-exp(-pars[2]))-x.bar/(1-n0/n))^2
}

ZIP.MME<-function(sample) {
  m<-mean(sample)
  s<-sum((sample-m)^2)/length(sample)
  pars<-c((s-m)/(s-m+m^2),(s-m+m^2)/m)
  # names(pars)<-c("p","lambda")
  pars<-pars[2:1]
  return(pars)
  # returns c(lambda,p)
}

logar.MLE<-function(sample){
  x.bar<-mean(sample)
  tmp<-optimize(logar.optfn,x.bar=x.bar,interval=c(0,1),tol=.Machine$double.eps^10)
  pars<-tmp$minimum
  return(pars)
}

logar.optfn<-function(Beta,x.bar){
  (x.bar+Beta/((1-Beta)*log(1-Beta)))^2
}

pois.MLE<-function(sample){
  lambda<-mean(sample)
  return(lambda)
}

negb.MLE<-function(sample){
  r<-negb.MLEr(sample)
  p<-r/(r+mean(sample))
  return(c(r,p))
}

negb.MLEr<-function(sample){
  r<-negb.MME(sample)[2]
  tmp<-optimize(negb.optfn,x=sample,interval=c(max(r-10,0),r+10),tol=.Machine$double.eps^10)
  pars<-tmp$minimum  
  return(pars)
}

negb.optfn<-function(r,x){
  n<-length(x)
  x.bar<-mean(x)
  k<-max(x)
  tmp<-rep(0,k)
  N<-rep(0,k)
  for (i in 1:k){
	h<-0
	for (j in 0:(i-1)){
	  h<-h+j/(1+j/r)
	}
	tmp[i]<-h
	N[i]<-sum(x==i)
  }
  (n*r^2*log(1+x.bar/r)-x.bar*n*(x.bar+r)/(1+x.bar/r)+sum(N*tmp))^2
}

negb.MME<-function(sample){
  m<-mean(sample)
  s<-sum((sample-m)^2)/length(sample)
  p<-m/s
  r<-m*p/(1-p)
  return(c(r,p))
}

laplace.MLE<-function(sample){
  a<-median(sample)
  b<-sum(abs(sample-a))/length(sample)
  return(c(a,b))
}

laplace.MME<-function(sample){
  a<-mean(sample)
  b<-sqrt(sum((sample-a)^2)/(2*length(sample)))
  return(c(a,b))
}

extrval.MLE<-function(sample){
 n<-length(sample)
 b<-sqrt(6)*sqrt(sum((sample-mean(sample))^2)/n)/pi
 tmp<-nls.lm(c(mean(sample)+digamma(1)*b,b),extrval.optfn,x=sample,n=n)
 pars<-tmp$par
 return(pars)
}

extrval.optfn<-function(pars,x,n){
 t1<-n-sum(exp((pars[1]-x)/pars[2]))
 t2<-sum((pars[1]-x)*exp((pars[1]-x)/pars[2]))+sum(x)-n*pars[1]-n*pars[2]
 return(c(t1,t2))
}

extrval.MME<-function(sample){
  b<-sqrt(6)*sqrt(sum((sample-mean(sample))^2)/length(sample))/pi
  a<-mean(sample)+digamma(1)*b
  return(c(a,b))
}

geom.MLE<-function(sample){
  p<-1/(1+mean(sample))
  return(p)
}

ZTP.MLE<-function(sample){
  x.bar<-mean(sample)
  tmp<-nls.lm(x.bar,ZTP.optfn,x.bar=x.bar)
  pars<-tmp$par
  return(pars)
}

ZTP.optfn<-function(lambda,x.bar){
  y<-Inf
  if (lambda>0) {
   y<-x.bar-lambda/(1-exp(-lambda))
  }
  return(y)
}

genpareto.MLE<-function(sample) {
	n<-length(sample)
	par_start<-genpareto.MME(sample)
	tmp <- try(nls.lm(par_start[2]/par_start[1], genpareto.optfn, sample = sample, n = n),silent=TRUE)
	if (!inherits(tmp,"try-error")) {
	  theta <- tmp$par
	  k <- -sum(log(1-theta*sample))/n
	  s <- k/theta
	  pars <- c(s, k)
	  if (is.nan(pars[1])||is.nan(pars[2])||is.na(pars[1])||is.na(pars[2])||is.infinite(pars[1])||is.infinite(pars[2])) {
		pars<-genpareto.MLE2(sample)
	  }
	  else {
		W1<-(1-k)/s^2*sum(sample/(1-k*sample/s)) - n/s
		W2<- -sum(log(1-k*sample/s))/k^2 - (1/k-1)/s*sum(sample/(1-k*sample/s))
		if (!is.na((W1^2+W2^2)>1)) {
			if ((W1^2+W2^2)>1) {pars<-genpareto.MLE2(sample)}
		}
	      else {pars<-genpareto.MLE2(sample)}
	  }
	}
	else {pars<-genpareto.MLE2(sample)}
	s <- pars[1]
	k <- pars[2]
	return(c(s, k))
}

genpareto.optfn<-function(pars,sample,n) { #pars is theta
	if (pars<(1/max(sample))) {
		k = -1/n*sum(log(1-pars*sample))	
		temp<-n/pars-(1/k-1)*sum(sample/(1-pars*sample))
	}
	else {temp<-Inf}
	return(temp)
}

genpareto.MLE2<-function(sample) { #methode van Choulakian (met herschaling)
	n<-length(sample)
	xbar<-mean(sample)
	x<-sample/xbar
	thetapos<-seq(-4*xbar,1/max(x),length.out=4000)
	start<- -n-sum(log(1-thetapos[1]*x))-n*log(-sum(log(1-thetapos[1]*x))/(n*thetapos[1]))
	thetadef<- thetapos[1]
	thetapos<-thetapos[-c(1,4000)]
	for (theta in thetapos) {
		tmp<--n-sum(log(1-theta*x))-n*log(-sum(log(1-theta*x))/(n*theta))
		if (tmp>start) {	
			start<-tmp
			thetadef<-theta
		}
	}
	theta<-thetadef
	k <- -sum(log(1-theta*x))/n
	s <- k/theta*xbar
	return(c(s,k))
}

genpareto.MME<-function(sample) { #only valid if shape=k > -1/2
	MEAN<-mean(sample)
	VAR<-sum((sample-MEAN)^2)/length(sample)
	k = (MEAN^2/VAR - 1)/2
	s = MEAN*(MEAN^2/VAR + 1)/2
	return(c(s,k))
}

betab.MLE<-function(sample,ntrials) {
	n<-length(sample)
	tmp<-nls.lm(betab.MME(sample,ntrials),betab.optfn,x=sample,ntrials=ntrials,n=n)	
	pars<-tmp$par
	return(pars) #these pars are the estimations of alpha and beta
}

betab.optfn<-function(pars,x,ntrials,n){
	a<-pars[1] #this is alpha
	b<-pars[2] #this is beta
	t1<-sum(digamma(x+a))-n*digamma(a+ntrials+b)-n*digamma(a)+n*digamma(a+b)
	t2<-sum(digamma(ntrials-x+b))-n*digamma(a+ntrials+b)-n*digamma(b)+n*digamma(a+b)
	return(c(t1,t2))
}

betab.MME<-function(sample,ntrials) { #the number of trials is known (fixed)
	n<-ntrials
	MEAN<-mean(sample)
	VAR<-sum((sample-MEAN)^2)/length(sample)
	alpha<-(n-MEAN-VAR/MEAN)*MEAN/(n*VAR/MEAN+MEAN-n)
	beta<-(n-MEAN-VAR/MEAN)*(n-MEAN)/(n*VAR/MEAN+MEAN-n)
	return(c(alpha,beta))
}

gamma.MLE<-function(sample) {
	xbar<-mean(sample)
	if (min(sample)<=0) {stop("There are negative sample observations so the MLE does not exist. ")}
	s<-log(xbar)-mean(log(sample))
	startpar<-(3-s+sqrt((s-3)^2+24*s))/(12*s) #alternatief is via gamma.MME startwaarde voor a bekomen
	tmp<-nls.lm(startpar, gamma.optfn, sample = sample)
	a <- tmp$par
	b <- xbar/a
	return(c(a,b)) # shape a and scale b
}

gamma.optfn<-function(pars,sample) { #pars is shape a
	if (pars>0) {
	temp<-log(pars)-digamma(pars)-log(mean(sample))+mean(log(sample))
	}
	else {temp<-Inf}
	return(temp)
}

gamma.MME<-function(sample) {
	MEAN <- mean(sample)
	VAR <- sum((sample - MEAN)^2)/length(sample)
	a<-MEAN^2/VAR
	b<-VAR/MEAN
	return(c(a,b)) # shape a and scale b
}

############################
#alternative functions (are less quick, but maybe better defined with respect to parameter bounderies)
############################

#logis.MLE<-function(sample) {
#	fit <- vglm(sample ~ 1, logistic2)
#	MU<-coef(fit)[1]
#	SIGMA<-loge(coef(fit)[2], inverse=TRUE)
#	return(c(MU,SIGMA))
#}

#ZIP.MLE<-function(sample) {
#	fit <- vglm(sample ~ 1, zipoisson)
#	phi<-logit(coef(fit)[1], inverse=TRUE) 
#	lambda<-loge(coef(fit)[2], inverse=TRUE)
#	return(c(lambda,phi))
#}

#logar.MLE<-function(sample) {
#	fit <- vglm(sample ~ 1, logff)
#	p<-logit(coef(fit), inverse=TRUE)
#	return(p)
#}

#negb.MLE<-function(sample) {
#	fit <- vglm(sample ~ 1, negbinomial)
#	mu<-loge(coef(fit)[1], inverse=TRUE) #equals E[X]=r(1-p)/p
#	r<-loge(coef(fit)[2], inverse=TRUE)
#	p<-r/(r+mu)
#	return(c(r,p))
#}

#extrval.MLE<-function(sample) {
#	fit <- vglm(sample ~ 1, gumbel)
#	a<-coef(fit)[1]
#	b<-loge(coef(fit)[2], inverse=TRUE)
#	return(c(a,b))
#}

#genpareto.MLE<-function(sample) {
#	fit <- vglm(sample ~ 1, gpd) #unfortunately offset=1 does not seem to work like mentioned in the help file
#	s<-loge(coef(fit)[1], inverse=TRUE)
#	k<-logoff(coef(fit)[2], earg = list(offset=0.5), inverse=TRUE)
#	return(c(s,-k)) #returns c(s,k)
#}

#betab.MLE<-function(sample,ntrials) {
#	n<-ntrials
#	fit <- vglm(cbind(sample,n-sample) ~ 1, betabin.ab)
#	alpha<-loge(coef(fit)[1], inverse=TRUE)
#	beta<-loge(coef(fit)[2], inverse=TRUE)
#	return(c(alpha,beta))
#}

#genpareto.MLE<-function(sample) {
#	n<-length(sample)
#	par_start<-genpareto.MME(sample)
#	tmp <- try(nls.lm(par_start[2]/par_start[1], genpareto.optfn, sample = sample, n = n),silent=TRUE)
#	if (!inherits(tmp,"try-error")) {
#	theta <- tmp$par
#	k <- -sum(log(1-theta*sample))/n
#	s <- k/theta
#	return(c(s, k))
#	}
#	else {stop("The MLE does not exist.")}
#}

#genpareto.optfn<-function(pars,sample,n) { #pars is theta
#	if (pars<(1/max(sample))) {
#		k = -1/n*sum(log(1-pars*sample))	
#		temp<-n/pars-(1/k-1)*sum(sample/(1-pars*sample))
#	}
#	else {temp<-Inf}
#	return(temp)
#}


