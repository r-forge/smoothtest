#moments is not allowed (bootstrap impossible), f is
#subset horizon is automated
#to do: improve graphs?

ddsmooth.test<-function(sample,max.order,horizon="order",criterion="AIC",distr="norm",method="NONE",
		pars=c(0,1),B=200,f=NA,typedistr="cont",plot=FALSE,plot.range=NULL,ntrials=NA,...) {
  # horizon = "order" or "subset"; when subset, the models are automatically constructed
  # moreover, max.order always refers to the real order! 
  # criterion = "AIC", "BIC", "MISE"

  if (!is.numeric(sample)) stop("The given i.i.d. sample has to be a numeric vector.")
  if ( ((trunc(max.order)-max.order)!=0) || (max.order<1) ) stop("The maximum order of the data driven smooth test has to be a strict positive integer.")
  if (!(method=="MLE"||method=="MME"||method=="NONE")) stop("The nuisance parameter estimation method is not 1 of 3 given possibilities.")
  if (!(criterion=="AIC"||criterion=="BIC"||criterion=="MISE")) stop("The criterion is not 1 of 3 given possibilities.")
  if (distr=="betab" && is.na(ntrials)) stop("The number of trials has to be given through \"ntrials\".")

  n<-length(sample)
  p.value<-NULL
  min.order<-NULL

  # minimal order and nuisance par estimation if method!="NONE"
  if (method=="MLE") {
    min.order<-switch(distr,
  	    norm=3,
	    exp=2,
	    logar=2,
	    pois=2,
	    geom=2,
	    ZTP=2,
   	    logis=1,
    	    negb=2,
     	    laplace=1,
	    ZIP=2,
  	    extrval=1,
	    unif=1,	
	    genpareto=1,
	    betab=1,
	    gamma=2)
    pars<-switch(distr,
  	    norm=norm.MLE(sample),
	    exp=exp.MLE(sample),
	    logar=logar.MLE(sample),
	    pois=pois.MLE(sample),
	    geom=geom.MLE(sample),
	    ZTP=ZTP.MLE(sample),
   	    logis=logis.MLE(sample),
    	    negb=negb.MLE(sample),
     	    laplace=laplace.MLE(sample),
	    ZIP=ZIP.MLE(sample),
  	    extrval=extrval.MLE(sample),
	    unif=unif.MLE(sample),
	    genpareto=genpareto.MLE(sample),
	    betab=betab.MLE(sample,ntrials),
	    gamma=gamma.MLE(sample)) # the given pars are only used when method=="NONE"
  }
  else if (method=="MME") {
    min.order<-switch(distr,
  	    norm=3,
	    exp=2,
	    logar=2,
	    pois=2,
	    geom=2,
	    ZTP=2,
   	    logis=3,
    	    negb=3,
     	    laplace=3,
	    ZIP=3,
  	    extrval=3,
	    unif=3,	
	    genpareto=3,
	    betab=3,
	    gamma=3) 
    pars<-switch(distr,
  	    norm=norm.MLE(sample),
	    exp=exp.MLE(sample),
	    logar=logar.MLE(sample),
	    pois=pois.MLE(sample),
	    geom=geom.MLE(sample),
	    ZTP=ZTP.MLE(sample),
   	    logis=logis.MME(sample),
    	    negb=negb.MME(sample),
     	    laplace=laplace.MME(sample),
	    ZIP=ZIP.MME(sample),
  	    extrval=extrval.MME(sample),
	    unif=unif.MME(sample),
	    genpareto=genpareto.MME(sample),
	    betab=betab.MME(sample,ntrials),
	    gamma=gamma.MME(sample))
  }
  else { #in this case automatically method="NONE"
    min.order<-switch(distr,
  	    norm=1,
	    exp=1,
	    logar=1,
	    pois=1,
	    geom=1,
	    ZTP=1,
   	    logis=1,
    	    negb=1,
     	    laplace=1,
	    ZIP=1,
  	    extrval=1,
	    unif=1,
	    genpareto=1,
	    betab=1,
	    gamma=1) 
  }
  if (is.null(min.order)) { #in this case the null distr is not specified by distr
	if (!is.function(f)) stop("The null distribution is not (correctly) specified.") 
					#i.e. nor through distr (1st choice), nor f (2nd choice)
	if (!(typedistr=="cont"||typedistr=="disc")) stop("The type of the null distribution \"f\" is not correctly specified.")
	if (method!="NONE") stop("The data driven smooth test for testing a distribution given by \"f\", 
					can only be performed without nuisance parameter estimation.")
	ind.f<-1 #to indicate that the null distribution is given by \"f\" and \"typedistr\"
	min.order<-1	
  }
  else {
	ind.f<-0
      if (max.order<min.order) stop("The maximum order is too small to perform a meaningfull 
						data driven smooth test for this null distribution.")
      if (criterion=="MISE" && method=="MLE" && !(distr == "norm" || distr == "exp" || distr == "logar" || distr == "pois" 
     	   || distr == "geom" || distr == "ZTP")) stop("The MISE criterion together with MLE is not implemented for this distribution.")
      if (criterion=="MISE" && method=="MME" && (distr == "ZIP" || distr == "negb" || distr == "genpareto" 
         || distr == "betab" || distr == "gamma")) stop("The MISE criterion (together with MME) is not implemented for this distribution.")
      if (criterion=="MISE" && method=="MME" && (distr == "logis" || distr == "laplace" || distr == "extrval" || distr == "unif") 
  	   && (max.order>8)) stop("The MISE criterion together with MME is only implemented until the 8th order for this distribution.")
	if (method=="NONE") { 
	      if (!is.numeric(pars)) stop("The nuisance parameter(s) have to be given properly,
						since no estimation method is chosen.")
     	  	if (is.na(pars[1])||is.nan(pars[1])||is.infinite(pars[1])) stop("The nuisance parameter(s) have to be given properly,
													since no estimation method is chosen.")
		if (length(pars)==2) {
		    if (is.na(pars[2])||is.nan(pars[2])||is.infinite(pars[2])) stop("The nuisance parameter(s) have to be given properly,
													since no estimation method is chosen.")
		}
	}
	else {
    		if (is.na(pars[1]) || is.nan(pars[1]) || is.infinite(pars[1])) stop("The nuisance parameter(s) estimation failed.")
		if (length(pars)==2) {
		    if (is.na(pars[2]) || is.nan(pars[2]) || is.infinite(pars[2])) stop("The nuisance parameter(s) estimation failed.")
		}
		if (distr=="genpareto") {warning("The data driven smooth test for the generalized Pareto distribution could give an error. 
						 If this is the case, use of the alternative estimation method is suggested.")}
		if (distr=="gamma" && method=="MME") {warning("The MME based data driven smooth test for the gamma distribution is unreliable. 
						 Therefore the MLE based data driven smooth test is recommended.")}
	}
  }

  # construction of the horizon
  if(horizon=="order") {
    n.S<-max.order-min.order+1
    S<-list()
    S.df<-rep(NA,n.S)
    for(i in 1:n.S) {
      S[[i]]<-(min.order:max.order)[1:i]
      S.df[i]<-length(S[[i]])
    }
  }
  if(horizon=="subset") {
    models<-create.model(max.order-min.order+1)
    # max.order<-max(unlist(models))+min.order-1
    n.S<-length(models)
    S.df<-rep(NA,n.S)
    S<-models
    for(i in 1:n.S) {
      S[[i]]<-models[[i]]+min.order-1  # because models always start at component one
      S.df[i]<-length(models[[i]])
    }
  }

  # computing the score test statistics for models in horizon S
  if(criterion=="MISE") {
     AllStats<-MISE.sequence(sample,S=S,min.order=min.order,
	distr=distr,method=method,pars=pars,f=f,typedistr=typedistr,ntrials=ntrials)
     mises<-AllStats$mises
  }
  else {
     AllStats<-partialsum.sequence(sample,S=S,min.order=min.order,
	distr=distr,method=method,pars=pars,f=f,typedistr=typedistr,ntrials=ntrials)
  }
  stats<-AllStats$stats
  thetas<-AllStats$theta
  selector<-switch(criterion,
     AIC=stats-2*S.df,
     BIC=stats-S.df*log(n),
     MISE= -mises)
  
  modelerrors<-sum(is.na(stats))
  if (modelerrors!=0) {
	warning("There is (are) ", modelerrors, " constructed model(s) for which the test statistic is undefined.
		   This (these) model(s) is (are) removed from the considered set of models.")
  }
  model.dd<-which.max(selector) # no model will be chosen only if every component of selector is NA
  stat.dd<-NA
  try(stat.dd<-stats[model.dd],silent=TRUE) # only if no model is chosen this is a problem and stat.dd will be NA
  if (is.na(stat.dd)) stop("The data driven smooth test is not defined.")
  
  S.dd<-S[[model.dd]]
  theta<-thetas[S.dd-min.order+1]

  # bootstrap null distribution

  #if(!method=="NONE") {
  #   pars<-smooth.test(sample,distr=distr,method=method,B=NULL,order=min.order+1)$par.est
  #}

  if(!is.null(B)) {
    nd<-rep(NA,B)
    for(i in 1:B) {
     if (ind.f==0) {
      x<-switch(distr,
           logis=rlogis(n, location=pars[1], scale=pars[2]),
           negb=rnbinom(n, size=pars[1], prob=pars[2]),
           laplace=rlaplace(n, location=pars[1], scale=pars[2]),
           ZIP=rzipois(n, lambda=pars[1], phi=pars[2]),
           extrval=rgumbel(n, location=pars[1], scale=pars[2]),
           norm=rnorm(n, mean=pars[1], sd=pars[2]),
           exp=rexp(n, rate=pars[1]),
	     logar=rlog(n, prob=pars[1]),	
           pois=rpois(n, lambda=pars[1]),
	     geom=rgeom(n, prob=pars[1]),
	     ZTP=rZTP(n, lambda=pars[1]),
	     unif=runif(n, min=pars[1], max=pars[2]),
	     genpareto=rgpd(n, scale=pars[1], shape=-pars[2]),
	     betab=rbetabin.ab(n, size=ntrials, shape1=pars[1], shape2=pars[2]),
	     gamma=rgamma(n, shape=pars[1], scale=pars[2]))
     }
     else {
      x<-switch(typedistr,
			cont=rsamplecont(n, f),
			disc=rsampledisc(n, f))
     }
     if(criterion=="MISE") {
         AllStats<-MISE.sequence(x,S=S,min.order=min.order,
		distr=distr,method=method,pars=pars,f=f,typedistr=typedistr,ntrials=ntrials)
         mises<-AllStats$mises
     }
     else {
         AllStats<-partialsum.sequence(x,S=S,min.order=min.order,
		distr=distr,method=method,pars=pars,f=f,typedistr=typedistr,ntrials=ntrials)
     }
     stats<-AllStats$stats
     selector<-switch(criterion,
        AIC=stats-2*S.df,
        BIC=stats-S.df*log(n),
        MISE=-mises)

     modelerrors<-sum(is.na(stats))
     if (modelerrors!=0) {
	  warning("There is (are) ", modelerrors, " constructed model(s) for which the test statistic in a bootstrap sample
		     is undefined. This (these) model(s) is (are) removed from the considered set of models.")
     }
     model.dd<-which.max(selector) # no model will be chosen only if every component of selector is NA
     try(nd[i]<-stats[model.dd],silent=TRUE) # only if no model is chosen this is a problem and nd[i] will be NA
    }
    p.value<-mean(stat.dd<=nd[!is.na(nd)]) # this way the bootstrap samples with NA test stat are ignored
    if (length(nd[is.na(nd)])!=0) warning(length(nd[is.na(nd)]), " bootstrap samples had to be ignored, 
							since their test statistic was undefined.")  
  }  

  carrier<-NULL
  dimproved<-NULL

  if(plot) {
    # nonpar density estimate (improved smooth)
    if (ind.f==0) {
     carrier<-function(x) {
       switch(distr,
           logis=dlogis(x, location=pars[1], scale=pars[2]),
           negb=dnbinom(x, size=pars[1], prob=pars[2]),
           laplace=dlaplace(x, location=pars[1], scale=pars[2]),
           ZIP=dzipois(x, lambda=pars[1], phi=pars[2]),
           extrval=dgumbel(x, location=pars[1], scale=pars[2]),
           norm=dnorm(x, mean=pars[1], sd=pars[2]),
           exp=dexp(x, rate=pars[1]),
	     logar=dlog(x, prob=pars[1]),	
           pois=dpois(x, lambda=pars[1]),
	     geom=dgeom(x, prob=pars[1]),
	     ZTP=dZTP(x, lambda=pars[1]),
	     unif=dunif(x, min=pars[1], max=pars[2]),
	     genpareto=dgpd(x, scale=pars[1], shape=-pars[2]),
	     betab=dbetabin.ab(x, size=ntrials, shape1=pars[1], shape2=pars[2]),
	     gamma=dgamma(x, shape=pars[1], scale=pars[2]))
     }
    }
    else {
     carrier<-f
    } 
    if (ind.f==0) {
     correction<-function(x) {
      y<-1
      for(i in 1:length(S.dd)) {
         y<-y+theta[i]*orth.poly(S.dd[i], distr = distr, pars = pars, ntrials = ntrials)(x)
      }
      return(y)
     }
    }
    else {
     correction<-function(x) {
      y<-1
      for(i in 1:length(S.dd)) {
         y<-y+theta[i]*orth.poly(S.dd[i], distr = "otherwise", f = f, typedistr = typedistr)(x)
      }
      return(y)
     }
    }
    dimproved<-function(x) {
      carrier(x)*correction(x)
    }

    if(distr=="negb"||distr=="ZIP"||distr=="pois"||distr=="logar"||distr=="geom"
     ||distr=="ZTP"||distr=="betab"||(ind.f==1&&typedistr=="disc")) { # all discrete distributions
      if(is.null(plot.range)) plot.range<-range(sample)
      plot.range[2]<-plot.range[2]+1
      x.range<-plot.range[1]:plot.range[2] 
      N<-length(x.range)
      dgajek<-gajek.d(dimproved,carrier,plot.range[1],plot.range[2],N=N)$fc  # I think that N may be limited for discrete distriubutions
      sample.tab<-table(sample)
      #sq<-min(sample):max(sample)
      sq<-x.range
      w<-rep(0,length(sq))
      w[sq%in%as.numeric(names(sample.tab))]<-sample.tab
      barplot(w/length(sample),space=0,col=0,names.arg=sq,
         main=paste("Hypothesised and improved density estimates of ",deparse(substitute(sample))),
         xlab=deparse(substitute(sample)),xlim=plot.range,...)

#      hist(sample,prob=T,
#          main=paste("Hypothesised and improved density estimates of ",deparse(substitute(sample))),
#          xlab=deparse(substitute(sample)),xlim=plot.range) # ToDo: improve graph; it is actually not exactly a boxplot that has to be plotted
      points(x.range+0.5,dgajek(x.range),type="p",pch=15,col="blue")
      points(x.range+0.5,carrier(x.range),type="p",pch=1)
    }
    if(distr=="norm"||distr=="logis"||distr=="laplace"||distr=="extrval"||distr=="exp"
     ||distr=="unif"||distr=="genpareto"||distr=="gamma"||(ind.f==1&&typedistr=="cont")) { # all continuous distributions
      if(is.null(plot.range)) plot.range<-range(sample)
      x.range<-seq(plot.range[1],plot.range[2],length.out=1000)
      N<-1000
      dgajek<-gajek(dimproved,carrier,lower=plot.range[1],upper=plot.range[2],N=N)$fc  
      hist.density<-hist(sample,plot=F)$density
      hist.breaks<-hist(sample,prob=T,
          main=paste("Hypothesised and improved density estimates of ",deparse(substitute(sample))),
          xlab=deparse(substitute(sample)),ylim=c(0,max(dgajek(x.range),carrier(x.range),hist.density)),...)$breaks #,xlim=plot.range
      #width<-hist.breaks[2]-hist.breaks[1]
	#x.range<-seq(min(plot.range[1],hist.breaks),max(plot.range[2],hist.breaks),length.out=1000)
      lines(x.range,dgajek(x.range),type="l",lty=1,col="blue")
      lines(x.range,carrier(x.range),type="l",lty=2)
      #lines(x.range,dgajek(x.range)/width,type="l",lty=1)
      #lines(x.range,carrier(x.range)/width,type="l",lty=2)
      #lines(x.range,dgajek(x.range)/(plot.range[2]-plot.range[1]),type="l",lty=1)
      #lines(x.range,carrier(x.range)/(plot.range[2]-plot.range[1]),type="l",lty=2)
    }
  }

  cat("Data-Driven Smooth goodness-of-fit test\n")
  cat("Null hypothesis:", distr, "against", max.order, 
      "th order alternative\n")
  if (distr=="betab") {cat("The number of trials are assumed known and fixed as", ntrials, "\n")}
  cat("Nuisance parameter estimation:", method, "\n")
  if (method!="NONE") {
	par_est <- pars
      names(par_est)<-switch(distr,
  	    norm=c("mean","sd"),
	    exp="rate",
	    logar="prob",
	    pois="lambda",
	    geom="prob",
	    ZTP="lambda",
   	    logis=c("location","scale"),
    	    negb=c("size","prob"),
     	    laplace=c("location","scale"),
	    ZIP=c("lambda","phi"),
  	    extrval=c("location","scale"),
	    unif=c("min","max"),
	    genpareto=c("scale","-shape"),
	    betab=c("shape1","shape2"),
	    gamma=c("shape","scale"))
	cat("Parameter estimates:", par_est, " (",names(par_est),")\n")
  } 
  else {
	par_est <- "no parameter estimation needed"
  }
  cat("Horizon:",horizon,"\n")
  cat("Selection criterion:",criterion,"\n\n")
  cat("Data-Driven Smooth test statistic S_k =", stat.dd, " p-value =", p.value,"\n")
  cat("Selected model:",S.dd,"\n\n")
  if (is.null(B)) {
      cat("No p-values calculated; only the bootstrap is allowed here\n\n")
  }
  else {
      cat("All p-values are obtained by the bootstrap with", B, "runs\n\n")
  }
  invisible(list(stat=stat.dd,model=S.dd,p.value=p.value,par.est=par_est,dnull=carrier,dimp=dimproved))
}


## auxiliary functions

# to automate the modelbuilding when horizon="subset" is chosen

create.model<-function(k) {
  models<-list()
  for (mod.length in 1:k) {
	mod.extra<-combn(k,mod.length,simplify=F)
	for (i in 1:choose(k,mod.length)) {		#or 1:length(mod.extra)
		models[[length(models)+1]]<-mod.extra[[i]]
	}
  }
  return(models)
}

partialsum.sequence<-function(sample,S,min.order,...) {
  # returns all score statistics within the horizon S
  # min.order depends on method of estimation and distribution; 
  #   it is determined in ddsmooth.test
  max.order<-max(unlist(S))
  st<-smooth.test(sample,order=max.order,B=NULL,output=F,...)
  Sigma<-st$Sigma
  U<-matrix(st$comp,ncol=1) #$comp are now the unscaled components
  theta<-U/sqrt(length(sample))
  stats<-rep(NA,length(S))
  for(i in 1:length(S)) {
    ind<-S[[i]]-min.order+1
    temp<-try(t(U[ind])%*%solve(Sigma[ind,ind])%*%U[ind],silent=TRUE)
    if (!inherits(temp,"try-error")) { #det(Sigma[ind,ind]) \neq 0
	stats[i]<-temp
    }
    else {
	stats[i]<-NA
    }
  }
  return(list(stats=stats,theta=theta))
}

MISE.sequence<-function(sample,S,min.order,...) {
  # returns all theta, and the MISEs within the horizon S
  # min.order depends on method of estimation and distribution; 
  #   it is determined in ddsmooth.test
  max.order<-max(unlist(S))
  st<-smooth.test(sample,rescale=T,order=max.order,B=NULL,output=F,...) 
  Sigma<-st$Sigma
  EVar<-st$EVar/length(sample)
  theta<-st$comp/sqrt(length(sample))
  theta2<-theta^2
  U<-matrix(st$comp,ncol=1)
  stats<-rep(NA,length(S))
  mises<-stats
  for(i in 1:length(S)) {
    ind<-S[[i]]-min.order+1 # indices corresponding to S
    nind<-!(1:(max.order-min.order+1))%in%ind # indices corresponding to \bar{S}
    temp<-try(t(U[ind])%*%solve(Sigma[ind,ind])%*%U[ind],silent=TRUE)
    if (!inherits(temp,"try-error")) { #det(Sigma[ind,ind]) \neq 0
	stats[i]<-temp
	mises[i]<-sum(EVar[ind])+sum(theta2[nind])-sum(EVar[nind])
    }
    else {
	stats[i]<-NA
	mises[i]<-NA # to ensure that this model is removed from the set of models so it is surely not selected 
    }
  }
  return(list(mises=mises,stats=stats,theta=theta))
}

gajek<-function(f,g,lower=0,upper=1,N=100) {
 sx<-seq(lower,upper,1/N)
 f.sx<-f(sx)
 eps.s<-seq(0,max(f.sx),1/N)
 integrals<-c()
 for(eps in eps.s) {
   fc<-function(x) {
     d<-f(x)-eps*g(x)
     d[d<0]<-0
     return(d)  
   }
   i<-integrate(fc,lower=lower,upper=upper)$value
   integrals<-c(integrals,i)
 }
 errors<-abs(integrals-1)
 error<-min(errors)

 int<-integrals[errors==error]
 epsilon<-eps.s[errors==error]
 fc<-function(x) {
     d<-f(x)-epsilon*g(x)
     d[d<0]<-0
     return(d)  
 }
 r<-list(error=error,integral=int,epsilon=epsilon,fc=fc)
 return(r)
}

gajek.d<-function(f,g,lower=0,upper=1,N=100) {
 # gajek correction for discrete distributions
 sx<-lower:upper
 f.sx<-f(sx)
 eps.s<-seq(0,max(f.sx),1/N)
 integrals<-c()
 for(eps in eps.s) {
   fc<-function(x) {
     d<-f(x)-eps*g(x)
     d[d<0]<-0
     return(d)  
   }
   i<-sum(fc(lower:upper))
   integrals<-c(integrals,i)
 }
 errors<-abs(integrals-1)
 error<-min(errors)

 int<-integrals[errors==error]
 epsilon<-eps.s[errors==error]
 fc<-function(x) {
     d<-f(x)-epsilon*g(x)
     d[d<0]<-0
     return(d)  
 }
 r<-list(error=error,integral=int,epsilon=epsilon,fc=fc)
 return(r)
}

