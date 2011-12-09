# the function now also invisibles the HK-variance-estimator/NVE 
# for all MME I've implemented a Cholesky decomposition (chol=T) 
# note that $comp are now the unscaled components (nor th. nor emp.) and also not quadratic
# change in test.stat: if "pars" is given these are used as the MLEs or MMEs when this method is chosen
# but if "pars" is NA the estimates are made according to the chosen method of estimation

test.stat<-function (degree, distr = "unif", pars = c(0, 1), sample, method = "NONE",
    f = NA, moments = NA, typedistr = "cont", chol = FALSE, ntrials=NA) 
{
    n <- length(sample)
    if (method == "NONE") {
        q <- 0
    }
    else { # method will for sure be "MME" or "MLE" and distr is from the list
        if (distr == "unif" || distr == "norm" || distr == "logis" || 
            distr == "ZIP" || distr == "negb" || distr == "laplace" || 
            distr == "extrval" || distr == "genpareto" || distr == "betab" || distr == "gamma") {
            q <- 2
        }
        # if (distr == "exp" || distr == "logar" || distr == "pois" || 
        #    distr == "geom" || distr == "ZTP") {
	  else {
            q <- 1
        }
	  if (method == "MLE" && is.na(pars[1])) { # else pars are already the MLEs
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
	    gamma=gamma.MLE(sample))
	  }
	  else if (method == "MME" && is.na(pars[1])) { # else pars are already the MMEs
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
    }
    # now pars are the (estimated) nuisance parameters and together with the method they determine the varcov matrix
    if (method == "NONE") { # else: distr is from the list (stops ensure this)
        C <- diag(rep(1, degree))
    }
    else if (distr == "norm" || distr == "exp" || distr == "logar" || 
        distr == "pois" || distr == "geom" || distr == "ZTP") {
        C <- diag(rep(1, degree - q))
    }
    else if (degree<5 && method=="MME" && distr!="negb") {
	  C<-switch(distr,
  	    logis=logis.covarMME(),
    	    laplace=laplace.covarMME(),
	    ZIP=ZIP.covarMME(lambda=pars[1],p=pars[2]),
  	    extrval=extrval.covarMME(),
	    genpareto=genpareto.covarMME(k=pars[2]),
	    betab=betab.covarMME(a=pars[1],b=pars[2],n=ntrials),
	    gamma=gamma.covarMME(a=pars[1]),
	    unif=unif.covarMME())
	  C<-C[1:(degree-2),1:(degree-2)]
    }
    else if (degree<5 && method=="MLE" && (distr=="logis" || distr=="laplace" || distr=="extrval")) {
	  C<-switch(distr,
  	    logis=logis.covarMLE(),
    	    laplace=laplace.covarMLE(),
  	    extrval=extrval.covarMLE())
	  C<-C[1:degree,1:degree]
    }
    else if (degree<5 && method=="MLE" && (distr=="ZIP" || distr=="gamma")) { # the 1st component disappears (also when distr="negb")
	  C<-switch(distr,
  	    ZIP=ZIP.covarMLE(lambda=pars[1],p=pars[2]),
  	    gamma=gamma.covarMLE(a=pars[1]))
	  C<-C[1:(degree-1),1:(degree-1)]
    }
    else if (distr == "logis") {
        MU <- 0
        SIGMA <- 1
        f <- function(x) {
            y <- dlogis(x, location = MU, scale = SIGMA)
            invisible(y)
        }
        S1 <- function(x) {
            y <- (1 - exp(-(x - MU)/SIGMA))/(SIGMA * (1 + exp(-(x - 
                MU)/SIGMA)))
            invisible(y)
        }
        S2 <- function(x) {
            y <- ((x - MU) * (1 - exp(-(x - MU)/SIGMA)) - SIGMA * 
                (1 + exp(-(x - MU)/SIGMA)))/(SIGMA^2 * (1 + exp(-(x - 
                MU)/SIGMA)))
            invisible(y)
        }
        A1 <- rep(0, degree)
        A2 <- rep(0, degree)
        for (i in 1:degree) {
            h <- orth.poly(i, distr = "logis", pars = c(MU, SIGMA))
            g <- function(x) {
                y <- h(x) * S1(x) * f(x)
                y[is.nan(y)] <- 0
                invisible(y)
            }
            templist <- integrate(g, -Inf, Inf)
            A1[i] <- templist$value
            g <- function(x) {
                y <- h(x) * S2(x) * f(x)
                y[is.nan(y)] <- 0
                invisible(y)
            }
            templist <- integrate(g, -Inf, Inf)
            A2[i] <- templist$value
        }
        A <- cbind(A1, A2)
        B <- matrix(0, 2, 2)
        g <- function(x) {
            y <- S1(x)^2 * f(x)
            y[is.nan(y)] <- 0
            invisible(y)
        }
        templist <- integrate(g, -Inf, Inf)
        B[1, 1] <- templist$value
        g <- function(x) {
            y <- S1(x) * S2(x) * f(x)
            y[is.nan(y)] <- 0
            invisible(y)
        }
        templist <- integrate(g, -Inf, Inf)
        B[1, 2] <- templist$value
        B[2, 1] <- B[1, 2]
        g <- function(x) {
            y <- S2(x)^2 * f(x)
            y[is.nan(y)] <- 0
            invisible(y)
        }
        templist <- integrate(g, -Inf, Inf)
        B[2, 2] <- templist$value
        if (method == "MLE") {
            I <- diag(rep(1, degree))
            C <- I - A %*% solve(B) %*% t(A)
        }
        if (method == "MME") {
            I <- diag(rep(1, degree - 2))
            K <- A[3:degree, ] %*% solve(A[1:2, ])
            C <- I + K %*% t(K)
        }
    }
    else if (distr == "laplace") {
        a <- 0
        b <- 1
        f <- function(x) {
            y <- exp(-abs(x - a)/b)/(2 * b)
            invisible(y)
        }
        S1 <- function(x) {
            y <- sign(x - a)/b
            invisible(y)
        }
        S2 <- function(x) {
            y <- (abs(x - a) - b)/b^2
            invisible(y)
        }
        A1 <- rep(0, degree)
        A2 <- rep(0, degree)
        for (i in 1:degree) {
            h <- orth.poly(i, distr = "laplace", pars = c(a,b))
            g <- function(x) {
                y <- h(x) * S1(x) * f(x)
                y[is.nan(y)] <- 0
                invisible(y)
            }
            templist <- integrate(g, -Inf, Inf)
            A1[i] <- templist$value
            g <- function(x) {
                y <- h(x) * S2(x) * f(x)
                y[is.nan(y)] <- 0
                invisible(y)
            }
            templist <- integrate(g, -Inf, Inf)
            A2[i] <- templist$value
        }
        A <- cbind(A1, A2)
        B <- matrix(0, 2, 2)
        g <- function(x) {
            y <- S1(x)^2 * f(x)
            y[is.nan(y)] <- 0
            invisible(y)
        }
        templist <- integrate(g, -Inf, Inf)
        B[1, 1] <- templist$value
        g <- function(x) {
            y <- S1(x) * S2(x) * f(x)
            y[is.nan(y)] <- 0
            invisible(y)
        }
        templist <- integrate(g, -Inf, Inf)
        B[1, 2] <- templist$value
        B[2, 1] <- B[1, 2]
        g <- function(x) {
            y <- S2(x)^2 * f(x)
            y[is.nan(y)] <- 0
            invisible(y)
        }
        templist <- integrate(g, -Inf, Inf)
        B[2, 2] <- templist$value
        if (method == "MLE") {
            I <- diag(rep(1, degree))
            C <- I - A %*% solve(B) %*% t(A)
        }
        if (method == "MME") {
            I <- diag(rep(1, degree - 2))
            K <- A[3:degree, ] %*% solve(A[1:2, ])
            C <- I + K %*% t(K)
        }
    }
    else if (distr == "extrval") {
        a <- 0
        b <- 1
        f <- function(x) {
            y <- exp((a - x)/b - exp((a - x)/b))/b
            invisible(y)
        }
        S1 <- function(x) {
            y <- 1/b - 1/b * exp((a - x)/b)
            invisible(y)
        }
        S2 <- function(x) {
            y <- (-a + x + exp((a - x)/b) * a - exp((a - x)/b) * 
                x - b)/b^2
            invisible(y)
        }
        A1 <- rep(0, degree)
        A2 <- rep(0, degree)
        for (i in 1:degree) {
            h <- orth.poly(i, distr = "extrval", pars = c(a,b))
            g <- function(x) {
                y <- h(x) * S1(x) * f(x)
                y[is.nan(y)] <- 0
                invisible(y)
            }
            templist <- integrate(g, -Inf, Inf)
            A1[i] <- templist$value
            g <- function(x) {
                y <- h(x) * S2(x) * f(x)
                y[is.nan(y)] <- 0
                invisible(y)
            }
            templist <- integrate(g, -Inf, Inf)
            A2[i] <- templist$value
        }
        A <- cbind(A1, A2)
        B <- matrix(0, 2, 2)
        g <- function(x) {
            y <- S1(x)^2 * f(x)
            y[is.nan(y)] <- 0
            invisible(y)
        }
        templist <- integrate(g, -Inf, Inf)
        B[1, 1] <- templist$value
        g <- function(x) {
            y <- S1(x) * S2(x) * f(x)
            y[is.nan(y)] <- 0
            invisible(y)
        }
        templist <- integrate(g, -Inf, Inf)
        B[1, 2] <- templist$value
        B[2, 1] <- B[1, 2]
        g <- function(x) {
            y <- S2(x)^2 * f(x)
            y[is.nan(y)] <- 0
            invisible(y)
        }
        templist <- integrate(g, -Inf, Inf)
        B[2, 2] <- templist$value
        if (method == "MLE") {
            I <- diag(rep(1, degree))
            C <- I - A %*% solve(B) %*% t(A)
        }
        if (method == "MME") {
            I <- diag(rep(1, degree - 2))
            K <- A[3:degree, ] %*% solve(A[1:2, ])
            C <- I + K %*% t(K)
        }
    }
    else if (distr == "unif") {
	  if (method == "MLE") {
		C<-diag(rep(1, degree)) # this seams to be true by simulation
	  }
	  if (method == "MME") {
            a <- 0
            b <- 1
            f <- function(x) {
                y <- dunif(x, min=a, max=b)
                invisible(y)
            }
            S1 <- function(x) {
                y <- 1/(b-a)	#this is not entirely correct, but will be fixed while integrating
                invisible(y)
            }
            S2 <- function(x) {
                y <- -1/(b-a)
                invisible(y)
            }
            A1 <- rep(0, degree)
            A2 <- rep(0, degree)
            for (i in 1:degree) {
                h <- orth.poly(i, distr = "unif", pars = c(a,b))
                g <- function(x) {
                    y <- h(x) * S1(x) * f(x)
                    y[is.nan(y)] <- 0
                    invisible(y)
                }
                templist <- integrate(g, a, b)
                A1[i] <- templist$value - h(a)*f(a) #this is the correction
                g <- function(x) {
                    y <- h(x) * S2(x) * f(x)
                    y[is.nan(y)] <- 0
                    invisible(y)
                }
                templist <- integrate(g, a, b)
		    A2[i] <- templist$value + h(b)*f(b) #this is the correction           
            }
            A <- cbind(A1, A2)
            I <- diag(rep(1, degree - 2))
            K <- A[3:degree, ] %*% solve(A[1:2, ])
            C <- I + K %*% t(K)
        }
    }
    else if (method == "MLE") {
        if (distr == "ZIP") {
            lambda <- pars[1]
            p <- pars[2]
            f <- function(x) {
                y <- dpois(x, lambda) * (1 - p)
                y[x == 0] <- p + (1 - p) * exp(-lambda)
                invisible(y)
            }
            S1 <- function(x) {
                y <- -1/(1 - p) + 0 * x
                y[x == 0] <- (1 - exp(-lambda))/(p + (1 - p) * 
                  exp(-lambda))
                invisible(y)
            }
            S2 <- function(x) {
                y <- -1 + x/lambda
                y[x == 0] <- -((1 - p) * exp(-lambda))/(p + (1 - 
                  p) * exp(-lambda))
                invisible(y)
            }
            A1 <- rep(0, degree)
            A2 <- rep(0, degree)
            x <- 0:500
            for (i in 2:degree) {
                h <- orth.poly(i, distr = "ZIP", pars = c(lambda,p))
                A1[i] <- sum(h(x) * S1(x) * f(x))
                A2[i] <- sum(h(x) * S2(x) * f(x))
            }
            A <- cbind(A1, A2)
            B <- matrix(0, 2, 2)
            B[1, 1] <- sum(S1(x)^2 * f(x))
            B[1, 2] <- sum(S1(x) * S2(x) * f(x))
            B[2, 1] <- B[1, 2]
            B[2, 2] <- sum(S2(x)^2 * f(x))
            I <- diag(rep(1, degree - 1))
            C <- I - A[2:degree,] %*% solve(B) %*% t(A[2:degree,])
        }
        else if (distr == "negb") {
            r <- pars[1]
            p <- pars[2]
            f <- function(x) {
                y <- dnbinom(x, size = r, prob = p)
                invisible(y)
            }
            S1 <- function(x) {
                y <- (-r + r * p + x * p)/(p * (p - 1))
                invisible(y)
            }
            S2 <- function(x) {
                y <- digamma(r + x) - digamma(r) + log(p)
                invisible(y)
            }
            A1 <- rep(0, degree)
            A2 <- rep(0, degree)
            x <- 0:500
            for (i in 2:degree) {
                h <- orth.poly(i, distr = "negb", pars = c(r,p))
                A1[i] <- sum(h(x) * S1(x) * f(x))
                A2[i] <- sum(h(x) * S2(x) * f(x))
            }
            A <- cbind(A1, A2)
            B <- matrix(0, 2, 2)
            B[1, 1] <- sum(S1(x)^2 * f(x))
            B[1, 2] <- sum(S1(x) * S2(x) * f(x))
            B[2, 1] <- B[1, 2]
            B[2, 2] <- sum(S2(x)^2 * f(x))
            I <- diag(rep(1, degree-1))
            C <- I - A[2:degree,] %*% solve(B) %*% t(A[2:degree,])
        }
        else if (distr == "genpareto") {
		s <- pars[1]
            k <- pars[2]
		if (!(k>-1/(2*degree))) stop("The smooth test cannot be performed, 
					since the estimated moments are not finite.")
            f <- function(x) {
                y <- dgpd(x, location=0, scale=s, shape=-k)
                invisible(y)
            }
            S1 <- function(x) {
                y <- -1/s+(1/k-1)*k*x/s^2/(1-k*x/s)	#this is not entirely correct
                invisible(y)
            }
            S2 <- function(x) {
                y <- -1/k^2*log(1-k*x/s)+(1/k-1)*(-x)/s/(1-k*x/s) #this is not entirely correct
                invisible(y)
            }
            A1 <- rep(0, degree)
            A2 <- rep(0, degree)
            for (i in 1:degree) {
                h <- orth.poly(i, distr = "genpareto", pars = c(s,k))
                g <- function(x) {
                    y <- h(x) * S1(x) * f(x)
                    y[is.nan(y)] <- 0
                    invisible(y)
                }
		    if (k>0) {templist <- integrate(g, 0, s/k)} else {templist <- integrate(g, 0, Inf)}
                A1[i] <- templist$value #there is no real correction
                g <- function(x) {
                    y <- h(x) * S2(x) * f(x)
                    y[is.nan(y)] <- 0
                    invisible(y)
                }
                if (k>0) {templist <- integrate(g, 0, s/k)} else {templist <- integrate(g, 0, Inf)}
		    A2[i] <- templist$value #there is no real correction           
            }
            A <- cbind(A1, A2)
	      B <- matrix(0, 2, 2)
            g <- function(x) {
                y <- S1(x)^2 * f(x)
                y[is.nan(y)] <- 0
                invisible(y)
            }
            if (k>0) {templist <- integrate(g, 0, s/k)} else {templist <- integrate(g, 0, Inf)}
            B[1, 1] <- templist$value 
            g <- function(x) {
                y <- S1(x) * S2(x) * f(x)
                y[is.nan(y)] <- 0
                invisible(y)
            }
            if (k>0) {templist <- integrate(g, 0, s/k)} else {templist <- integrate(g, 0, Inf)}
            B[1, 2] <- templist$value 
            B[2, 1] <- B[1, 2]
            g <- function(x) {
                y <- S2(x)^2 * f(x)
                y[is.nan(y)] <- 0
                invisible(y)
            }
            if (k>0) {templist <- integrate(g, 0, s/k)} else {templist <- integrate(g, 0, Inf)}
            B[2, 2] <- templist$value 
            I <- diag(rep(1, degree))
            C <- I - A %*% solve(B) %*% t(A)
        }
	  else if (distr == "betab") {
		a <- pars[1]
            b <- pars[2]
            f <- function(x) {
                y <- dbetabin.ab(x, size=ntrials, shape1=a, shape2=b)
                invisible(y)
            }
            S1 <- function(x) {
                y <- digamma(x+a)-digamma(a+ntrials+b)-digamma(a)+digamma(a+b)
                invisible(y)
            }
            S2 <- function(x) {
                y <- digamma(ntrials-x+b)-digamma(a+ntrials+b)-digamma(b)+digamma(a+b)
                invisible(y)
            }
            A1 <- rep(0, degree)
            A2 <- rep(0, degree)
            x <- 0:ntrials
            for (i in 1:degree) {
                h <- orth.poly(i, distr = "betab", pars = c(a,b), ntrials=ntrials)
                A1[i] <- sum(h(x) * S1(x) * f(x))
                A2[i] <- sum(h(x) * S2(x) * f(x))
            }
            A <- cbind(A1, A2)
            B <- matrix(0, 2, 2)
            B[1, 1] <- sum(S1(x)^2 * f(x))
            B[1, 2] <- sum(S1(x) * S2(x) * f(x))
            B[2, 1] <- B[1, 2]
            B[2, 2] <- sum(S2(x)^2 * f(x))
            I <- diag(rep(1, degree))
            C <- I - A %*% solve(B) %*% t(A)
	  }
        else if (distr == "gamma") {
		a <- pars[1] # shape
            b <- pars[2] # scale
		f <- function(x) {
                y <- dgamma(x, shape=a, scale=b)
                invisible(y)
            }
            S1 <- function(x) {
                y <- log(x)-digamma(a)-log(b)
                invisible(y)
            }
            S2 <- function(x) {
                y <- -(-x+b*a)/b^2
                invisible(y)
            }
            A1 <- rep(0, degree)
            A2 <- rep(0, degree)
            for (i in 2:degree) {
                h <- orth.poly(i, distr = "gamma", pars = c(a,b))
                g <- function(x) {
                    y <- h(x) * S1(x) * f(x)
                    y[is.nan(y)] <- 0
                    invisible(y)
                }
		    templist <- integrate(g, 0, Inf)
                A1[i] <- templist$value 
                g <- function(x) {
                    y <- h(x) * S2(x) * f(x)
                    y[is.nan(y)] <- 0
                    invisible(y)
                }
                templist <- integrate(g, 0, Inf)
		    A2[i] <- templist$value            
            }
            A <- cbind(A1, A2)
	      B <- matrix(0, 2, 2)
            g <- function(x) {
                y <- S1(x)^2 * f(x)
                y[is.nan(y)] <- 0
                invisible(y)
            }
            templist <- integrate(g, 0, Inf)
            B[1, 1] <- templist$value 
            g <- function(x) {
                y <- S1(x) * S2(x) * f(x)
                y[is.nan(y)] <- 0
                invisible(y)
            }
            templist <- integrate(g, 0, Inf)
            B[1, 2] <- templist$value 
            B[2, 1] <- B[1, 2]
            g <- function(x) {
                y <- S2(x)^2 * f(x)
                y[is.nan(y)] <- 0
                invisible(y)
            }
            templist <- integrate(g, 0, Inf)
            B[2, 2] <- templist$value 
            I <- diag(rep(1, degree - 1))
            C <- I - A[2:degree,] %*% solve(B) %*% t(A[2:degree,])
        }
    }
    else if (method == "MME") {
        if (distr == "ZIP") {
		lambda <- pars[1]
            p <- pars[2]
            f <- function(x) {
                y <- dpois(x, lambda) * (1 - p)
                y[x == 0] <- p + (1 - p) * exp(-lambda)
                invisible(y)
            }
            S1 <- function(x) {
                y <- -1/(1 - p) + 0 * x
                y[x == 0] <- (1 - exp(-lambda))/(p + (1 - p) * 
                  exp(-lambda))
                invisible(y)
            }
            S2 <- function(x) {
                y <- -1 + x/lambda
                y[x == 0] <- -((1 - p) * exp(-lambda))/(p + (1 - 
                  p) * exp(-lambda))
                invisible(y)
            }
            A1 <- rep(0, degree)
            A2 <- rep(0, degree)
            x <- 0:500
            for (i in 1:degree) {
                h <- orth.poly(i, distr = "ZIP", pars = c(lambda,p))
                A1[i] <- sum(h(x) * S1(x) * f(x))
                A2[i] <- sum(h(x) * S2(x) * f(x))
            }
            A <- cbind(A1, A2)
            I <- diag(rep(1, degree - 2))
            K <- A[3:degree, ] %*% solve(A[1:2, ])
            C <- I + K %*% t(K)
        }
        else if (distr == "negb") {
		r <- pars[1]
            p <- pars[2]
            f <- function(x) {
                y <- dnbinom(x, size = r, prob = p)
                invisible(y)
            }
            S1 <- function(x) {
                y <- (-r + r * p + x * p)/(p * (p - 1))
                invisible(y)
            }
            S2 <- function(x) {
                y <- digamma(r + x) - digamma(r) + log(p)
                invisible(y)
            }
            A1 <- rep(0, degree)
            A2 <- rep(0, degree)
            x <- 0:500
            for (i in 1:degree) {
                h <- orth.poly(i, distr = "negb", pars = c(r,p))
                A1[i] <- sum(h(x) * S1(x) * f(x))
                A2[i] <- sum(h(x) * S2(x) * f(x))
            }
            A <- cbind(A1, A2)
            I <- diag(rep(1, degree - 2))
            K <- A[3:degree, ] %*% solve(A[1:2, ])
            C <- I + K %*% t(K)
        }
        else if (distr == "genpareto") {
		s <- pars[1]
            k <- pars[2]
		if (!(k>-1/(2*degree))) stop("The smooth test cannot be performed, 
					since the estimated moments are not finite.")
            f <- function(x) {
                y <- dgpd(x, location=0, scale=s, shape=-k)
                invisible(y)
            }
            S1 <- function(x) {
                y <- -1/s+(1/k-1)*k*x/s^2/(1-k*x/s)	#this is not entirely correct
                invisible(y)
            }
            S2 <- function(x) {
                y <- -1/k^2*log(1-k*x/s)+(1/k-1)*(-x)/s/(1-k*x/s) #this is not entirely correct
                invisible(y)
            }
            A1 <- rep(0, degree)
            A2 <- rep(0, degree)
            for (i in 1:degree) {
                h <- orth.poly(i, distr = "genpareto", pars = c(s,k))
                g <- function(x) {
                    y <- h(x) * S1(x) * f(x)
                    y[is.nan(y)] <- 0
                    invisible(y)
                }
                if (k>0) {templist <- integrate(g, 0, s/k)} else {templist <- integrate(g, 0, Inf)}
		    #if (k>0) {A1[i] <- templist$value + h(s/k)*f(s/k)} else {A1[i] <- templist$value}
                #there is no real correction since f(s/k)=0
		    A1[i] <- templist$value
                g <- function(x) {
                    y <- h(x) * S2(x) * f(x)
                    y[is.nan(y)] <- 0
                    invisible(y)
                }
                if (k>0) {templist <- integrate(g, 0, s/k)} else {templist <- integrate(g, 0, Inf)}
		    A2[i] <- templist$value #there is no real correction           
            }
            A <- cbind(A1, A2)
            I <- diag(rep(1, degree - 2))
            K <- A[3:degree, ] %*% solve(A[1:2, ])
            C <- I + K %*% t(K)
        }
        else if (distr == "betab") {
		a <- pars[1]
            b <- pars[2]
            f <- function(x) {
                y <- dbetabin.ab(x, size=ntrials, shape1=a, shape2=b)
                invisible(y)
            }
            S1 <- function(x) {
                y <- digamma(x+a)-digamma(a+ntrials+b)-digamma(a)+digamma(a+b)
                invisible(y)
            }
            S2 <- function(x) {
                y <- digamma(ntrials-x+b)-digamma(a+ntrials+b)-digamma(b)+digamma(a+b)
                invisible(y)
            }
            A1 <- rep(0, degree)
            A2 <- rep(0, degree)
            x <- 0:ntrials 
            for (i in 1:degree) {
                h <- orth.poly(i, distr = "betab", pars = c(a,b), ntrials=ntrials)
                A1[i] <- sum(h(x) * S1(x) * f(x))
                A2[i] <- sum(h(x) * S2(x) * f(x))
            }
            A <- cbind(A1, A2)
            I <- diag(rep(1, degree - 2))
            K <- A[3:degree, ] %*% solve(A[1:2, ])
            C <- I + K %*% t(K)
	  }
        else if (distr == "gamma") {
		a <- pars[1] # shape
            b <- pars[2] # scale
		f <- function(x) {
                y <- dgamma(x, shape=a, scale=b)
                invisible(y)
            }
            S1 <- function(x) {
                y <- log(x)-digamma(a)-log(b)
                invisible(y)
            }
            S2 <- function(x) {
                y <- -(-x+b*a)/b^2
                invisible(y)
            }
            A1 <- rep(0, degree)
            A2 <- rep(0, degree)
            for (i in 1:degree) {
                h <- orth.poly(i, distr = "gamma", pars = c(a,b))
                g <- function(x) {
                    y <- h(x) * S1(x) * f(x)
                    y[is.nan(y)] <- 0
                    invisible(y)
                }
		    templist <- integrate(g, 0, Inf)
		    A1[i] <- templist$value
                g <- function(x) {
                    y <- h(x) * S2(x) * f(x)
                    y[is.nan(y)] <- 0
                    invisible(y)
                }
                templist <- integrate(g, 0, Inf)
		    A2[i] <- templist$value            
            }
            A <- cbind(A1, A2)
            I <- diag(rep(1, degree - 2))
            K <- A[3:degree, ] %*% solve(A[1:2, ])
            C <- I + K %*% t(K)
        }
    }
    if (method == "MLE" && distr == "logis") {
        U <- rep(0, degree)
        EVar<-NA
        MU <- pars[1]
        SIGMA <- pars[2]
        for (i in 1:degree) {
            h <- orth.poly(i, distr = distr, pars = c(MU, SIGMA))
            U[i] <- sum(h(sample))/sqrt(n)
        }
	  comp <- U
        dim(U) <- c(degree, 1)
        Tstat <- t(U) %*% solve(C) %*% U
        Tstat <- c(U/sqrt(diag(C)), Tstat)
    }
    else if (method == "MME" && distr == "logis") {
        U <- rep(0, degree - 2)
        EVar <- U
	  MU <- pars[1]
        SIGMA <- pars[2]
        for (i in 3:degree) {
            h <- orth.poly(i, distr = distr, pars = c(MU, SIGMA))
            U[i - 2] <- sum(h(sample))/sqrt(n)
            EVar[i-2]<-NVE.logis(sample,order=i,parest=c(MU,SIGMA))
        }
        comp <- U
        dim(U) <- c(degree - 2, 1)
        Tstat <- t(U) %*% solve(C) %*% U
        if(chol) {
          tmp<-solve(C)
          d<-dim(tmp)[1]
          Sigma.inv.rev<-tmp[d:1,d:1]
          U.rev<-U[d:1,]
          M<-chol(Sigma.inv.rev)
          U<-M%*%U.rev
          Tstat <- c(U[d:1], Tstat)
        }
        else {
          Tstat <- c(U/sqrt(diag(C)), Tstat)
        }
    }
    else if (method == "MLE" && distr == "laplace") {
        U <- rep(0, degree)
        EVar<-NA
	  a <- pars[1]
        b <- pars[2]
        for (i in 1:degree) {
            h <- orth.poly(i, distr = distr, pars = c(a, b))
            U[i] <- sum(h(sample))/sqrt(n)
        }
        comp <- U
        dim(U) <- c(degree, 1)
        Tstat <- t(U) %*% solve(C) %*% U
        Tstat <- c(U/sqrt(diag(C)), Tstat)
    }
    else if (method == "MME" && distr == "laplace") {
        U <- rep(0, degree - 2)
        EVar<-U
	  a <- pars[1]
        b <- pars[2]
        for (i in 3:degree) {
            h <- orth.poly(i, distr = distr, pars = c(a, b))
            U[i - 2] <- sum(h(sample))/sqrt(n)
            EVar[i-2]<-NVE.laplace(sample,order=i,parest=c(a,b))
        }
        comp <- U
        dim(U) <- c(degree - 2, 1)
        Tstat <- t(U) %*% solve(C) %*% U
        if(chol) {
          tmp<-solve(C)
          d<-dim(tmp)[1]
          Sigma.inv.rev<-tmp[d:1,d:1]
          U.rev<-U[d:1,]
          M<-chol(Sigma.inv.rev)
          U<-M%*%U.rev
          Tstat <- c(U[d:1], Tstat)
        }
        else {
          Tstat <- c(U/sqrt(diag(C)), Tstat)
        }
    }
    else if (method == "MLE" && distr == "extrval") {
        U <- rep(0, degree)
        EVar<-NA
	  a <- pars[1]
        b <- pars[2]
        for (i in 1:degree) {
            h <- orth.poly(i, distr = distr, pars = c(a, b))
            U[i] <- sum(h(sample))/sqrt(n)
        }
        comp <- U
        dim(U) <- c(degree, 1)
        Tstat <- t(U) %*% solve(C) %*% U
        Tstat <- c(U/sqrt(diag(C)), Tstat)
    }
    else if (method == "MME" && distr == "extrval") {
        U <- rep(0, degree - 2)
        EVar<-U
        a <- pars[1]
        b <- pars[2]
        for (i in 3:degree) {
            h <- orth.poly(i, distr = distr, pars = c(a, b))
            U[i - 2] <- sum(h(sample))/sqrt(n)
            EVar[i-2]<-NVE.extrval(sample,order=i,parest=c(a,b))
        }
        comp <- U
        dim(U) <- c(degree - 2, 1)
        Tstat <- t(U) %*% solve(C) %*% U
        if(chol) {
          tmp<-solve(C)
          d<-dim(tmp)[1]
          Sigma.inv.rev<-tmp[d:1,d:1]
          U.rev<-U[d:1,]
          M<-chol(Sigma.inv.rev)
          U<-M%*%U.rev
          Tstat <- c(U[d:1], Tstat)
        }
        else {
          Tstat <- c(U/sqrt(diag(C)), Tstat)
        }
    }
    else if (method == "MLE" && distr == "ZIP") {
        U <- rep(0, degree - 1)
        EVar<-NA
        for (i in 2:degree) {
            h <- orth.poly(i, distr = distr, pars = pars)
            U[i - 1] <- sum(h(sample))/sqrt(n)
        }
        comp <- U
        dim(U) <- c(degree - 1, 1)
        Tstat <- t(U) %*% solve(C) %*% U
        Tstat <- c(U/sqrt(diag(C)), Tstat)
    }
    else if (method == "MME" && distr == "ZIP") {
        U <- rep(0, degree - 2)
        EVar<-NA
        for (i in 3:degree) {
            h <- orth.poly(i, distr = distr, pars = pars)
            U[i - 2] <- sum(h(sample))/sqrt(n)
        }
        comp <- U
        dim(U) <- c(degree - 2, 1)
        Tstat <- t(U) %*% solve(C) %*% U
        if(chol) {
          tmp<-solve(C)
          d<-dim(tmp)[1]
          Sigma.inv.rev<-tmp[d:1,d:1]
          U.rev<-U[d:1,]
          M<-chol(Sigma.inv.rev)
          U<-M%*%U.rev
          Tstat <- c(U[d:1], Tstat)
        }
        else {
          Tstat <- c(U/sqrt(diag(C)), Tstat)
        }
    }
    else if (method == "MLE" && distr == "negb") {
        U <- rep(0, degree-1)
        EVar<-NA
        for (i in 2:degree) {
            h <- orth.poly(i, distr = distr, pars = pars)
            U[i-1] <- sum(h(sample))/sqrt(n)
        }
        comp <- U
        dim(U) <- c(degree-1, 1)
        Tstat <- t(U) %*% solve(C) %*% U
        Tstat <- c(U/sqrt(diag(C)), Tstat)
    }
    else if (method == "MME" && distr == "negb") {
        U <- rep(0, degree - 2)
        EVar<-NA
        for (i in 3:degree) {
            h <- orth.poly(i, distr = distr, pars = pars)
            U[i - 2] <- sum(h(sample))/sqrt(n)
        }
        comp <- U
        dim(U) <- c(degree - 2, 1)
        Tstat <- t(U) %*% solve(C) %*% U
        if(chol) {
          tmp<-solve(C)
          d<-dim(tmp)[1]
          Sigma.inv.rev<-tmp[d:1,d:1]
          U.rev<-U[d:1,]
          M<-chol(Sigma.inv.rev)
          U<-M%*%U.rev
          Tstat <- c(U[d:1], Tstat)
        }
        else {
          Tstat <- c(U/sqrt(diag(C)), Tstat)
        }
    }
    else if (method == "MLE" && distr == "unif") {
        U <- rep(0, degree)
        EVar<-NA
	  MIN <- pars[1]
        MAX <- pars[2]
        for (i in 1:degree) {
            h <- orth.poly(i, distr = distr, pars = c(MIN, MAX))
            U[i] <- sum(h(sample))/sqrt(n)
        }
	  comp <- U
        dim(U) <- c(degree, 1)
        Tstat <- t(U) %*% solve(C) %*% U
        Tstat <- c(U/sqrt(diag(C)), Tstat)
    }
    else if (method == "MME" && distr == "unif") {
        U <- rep(0, degree - 2)
        EVar <- U
	  MIN <- pars[1]
        MAX <- pars[2]
        for (i in 3:degree) {
            h <- orth.poly(i, distr = distr, pars = c(MIN, MAX))
            U[i - 2] <- sum(h(sample))/sqrt(n)
            EVar[i-2]<-NVE.unif(sample,order=i,parest=c(MIN, MAX))
        }
        comp <- U
        dim(U) <- c(degree - 2, 1)
        Tstat <- t(U) %*% solve(C) %*% U
        if(chol) {
          tmp<-solve(C)
          d<-dim(tmp)[1]
          Sigma.inv.rev<-tmp[d:1,d:1]
          U.rev<-U[d:1,]
          M<-chol(Sigma.inv.rev)
          U<-M%*%U.rev
          Tstat <- c(U[d:1], Tstat)
        }
        else {
          Tstat <- c(U/sqrt(diag(C)), Tstat)
        }
    }
    else if (method == "MLE" && distr == "genpareto") {
        U <- rep(0, degree)
        EVar<-NA
        for (i in 1:degree) {
            h <- orth.poly(i, distr = distr, pars = pars)
            U[i] <- sum(h(sample))/sqrt(n)
        }
        comp <- U
        dim(U) <- c(degree, 1)
        Tstat <- t(U) %*% solve(C) %*% U
        Tstat <- c(U/sqrt(diag(C)), Tstat)
    }
    else if (method == "MME" && distr == "genpareto") {
        U <- rep(0, degree - 2)
        EVar<-NA
        for (i in 3:degree) {
            h <- orth.poly(i, distr = distr, pars = pars)
            U[i - 2] <- sum(h(sample))/sqrt(n)
        }
        comp <- U
        dim(U) <- c(degree - 2, 1)
        Tstat <- t(U) %*% solve(C) %*% U
        if(chol) {
          tmp<-solve(C)
          d<-dim(tmp)[1]
          Sigma.inv.rev<-tmp[d:1,d:1]
          U.rev<-U[d:1,]
          M<-chol(Sigma.inv.rev)
          U<-M%*%U.rev
          Tstat <- c(U[d:1], Tstat)
        }
        else {
          Tstat <- c(U/sqrt(diag(C)), Tstat)
        }
    }
    else if (method == "MLE" && distr == "betab") {
        U <- rep(0, degree)
        EVar<-NA
        for (i in 1:degree) {
            h <- orth.poly(i, distr = distr, pars = pars, ntrials=ntrials)
            U[i] <- sum(h(sample))/sqrt(n)
        }
        comp <- U
        dim(U) <- c(degree, 1)
        Tstat <- t(U) %*% solve(C) %*% U
        Tstat <- c(U/sqrt(diag(C)), Tstat)
    }
    else if (method == "MME" && distr == "betab") {
        U <- rep(0, degree - 2)
        EVar<-NA
        for (i in 3:degree) {
            h <- orth.poly(i, distr = distr, pars = pars, ntrials=ntrials)
            U[i - 2] <- sum(h(sample))/sqrt(n)
        }
        comp <- U
        dim(U) <- c(degree - 2, 1)
        Tstat <- t(U) %*% solve(C) %*% U
        if(chol) {
          tmp<-solve(C)
          d<-dim(tmp)[1]
          Sigma.inv.rev<-tmp[d:1,d:1]
          U.rev<-U[d:1,]
          M<-chol(Sigma.inv.rev)
          U<-M%*%U.rev
          Tstat <- c(U[d:1], Tstat)
        }
        else {
          Tstat <- c(U/sqrt(diag(C)), Tstat)
        }
    }
    else if (method == "MLE" && distr == "gamma") {
        U <- rep(0, degree - 1)
        EVar<-NA
        for (i in 2:degree) {
            h <- orth.poly(i, distr = distr, pars = pars)
            U[i - 1] <- sum(h(sample))/sqrt(n)
        }
        comp <- U
        dim(U) <- c(degree - 1, 1)
        Tstat <- t(U) %*% solve(C) %*% U
        Tstat <- c(U/sqrt(diag(C)), Tstat)
    }
    else if (method == "MME" && distr == "gamma") {
        U <- rep(0, degree - 2)
        EVar<-NA
        for (i in 3:degree) {
            h <- orth.poly(i, distr = distr, pars = pars)
            U[i - 2] <- sum(h(sample))/sqrt(n)
        }
        comp <- U
        dim(U) <- c(degree - 2, 1)
        Tstat <- t(U) %*% solve(C) %*% U
        if(chol) {
          tmp<-solve(C)
          d<-dim(tmp)[1]
          Sigma.inv.rev<-tmp[d:1,d:1]
          U.rev<-U[d:1,]
          M<-chol(Sigma.inv.rev)
          U<-M%*%U.rev
          Tstat <- c(U[d:1], Tstat)
        }
        else {
          Tstat <- c(U/sqrt(diag(C)), Tstat)
        }
    }
    else { #only MME=MLE or no parameter estimation remains
        U <- rep(0, degree - q)
        EVar<-U
        for (i in (q + 1):degree) {
            h <- orth.poly(i, distr, pars, f=f, moments=moments, 
				   typedistr=typedistr, ntrials=ntrials) # default method="NONE" is used
            U[i - q] <- sum(h(sample))/sqrt(n)
            EVar[i-q]<-var(h(sample))
        }
        comp <- U 
        dim(U) <- c(degree - q, 1)
        Tstat <- t(U) %*% solve(C) %*% U
        Tstat <- c(U/sqrt(diag(C)), Tstat)
    }
    invisible(list(Tstat = Tstat, comp = comp, Sigma=C, EVar=EVar))
}

