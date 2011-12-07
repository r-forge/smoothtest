orthpol_norm<-function(k,mu,sigma) {
	#location-scale family!
	if (k==0) {
        h<-function(x){
          y<-1
          return(y)
        }
      }
      else if (k==1) {
        h<-function(x){
          z<-(x-mu)/sigma
          y<-z
          return(y)
        }
      }
	else if (k==2) {
        h<-function(x){
          z<-(x-mu)/sigma
          y<-1/2*2^(1/2)*z^2-1/2*2^(1/2)
          return(y)
        }
      }
	else if (k==3) {
        h<-function(x){
          z<-(x-mu)/sigma
          y<-1/6*3^(1/2)*2^(1/2)*z^3-1/2*3^(1/2)*2^(1/2)*z
          return(y)
        }
      }
	else if (k==4) {
        h<-function(x){
          z<-(x-mu)/sigma
          y<-1/12*6^(1/2)*z^4-1/2*6^(1/2)*z^2+1/4*6^(1/2)
          return(y)
        }
      }
	else if (k==5) {
        h<-function(x){
          z<-(x-mu)/sigma
          y<-1/60*15^(1/2)*2^(1/2)*z^5-1/6*15^(1/2)*2^(1/2)*z^3+1/4*15^(1/2)*2^(1/2)*z
          return(y)
        }
      }
	else if (k==6) {
        h<-function(x){
          z<-(x-mu)/sigma
          y<-1/60*5^(1/2)*z^6-1/4*5^(1/2)*z^4+3/4*5^(1/2)*z^2-1/4*5^(1/2)
          return(y)
        }
      }
	else if (k==7) {
        h<-function(x){
          z<-(x-mu)/sigma
          y<-1/840*70^(1/2)*2^(1/2)*z^7-1/40*70^(1/2)*2^(1/2)*z^5+1/8*70^(1/2)*2^(1/2)*z^3-1/8*70^(1/2)*2^(1/2)*z
          return(y)
        }
      }
	else if (k==8) {
        h<-function(x){
          z<-(x-mu)/sigma
          y<-1/1680*70^(1/2)*z^8-1/60*70^(1/2)*z^6+1/8*70^(1/2)*z^4-1/4*70^(1/2)*z^2+1/16*70^(1/2)
          return(y)
        }
      }
	else if (k==9) {
        h<-function(x){
          z<-(x-mu)/sigma
          y<-1/5040*35^(1/2)*2^(1/2)*z^9-1/140*35^(1/2)*2^(1/2)*z^7+3/40*35^(1/2)*2^(1/2)*z^5-1/4*35^(1/2)*2^(1/2)*z^3+3/16*35^(1/2)*2^(1/2)*z
          return(y)
        }
      }
	else if (k==10) {
        h<-function(x){
          z<-(x-mu)/sigma
          y<-1/5040*7^(1/2)*z^10-1/112*7^(1/2)*z^8+1/8*7^(1/2)*z^6-5/8*7^(1/2)*z^4+15/16*7^(1/2)*z^2-3/16*7^(1/2)
          return(y)
        }
      }
	return(h)
}

orthpol_exp<-function(k,sigma) {
	#scale family!
	if (k==0) {
        h<-function(x){
          y<-1
          return(y)
        }
      }
      else if (k==1) {
        h<-function(x){
          z<-x/sigma
          y<--(-z+1)
          return(y)
        }
      }
	else if (k==2) {
        h<-function(x){
          z<-x/sigma
          y<-1/2*z^2-2*z+1
          return(y)
        }
      }
	else if (k==3) {
        h<-function(x){
          z<-x/sigma
          y<--(-1/6*z^3+3/2*z^2-3*z+1)
          return(y)
        }
      }
	else if (k==4) {
        h<-function(x){
          z<-x/sigma
          y<-1/24*z^4-2/3*z^3+3*z^2-4*z+1
          return(y)
        }
      }
	else if (k==5) {
        h<-function(x){
          z<-x/sigma
          y<--(-1/120*z^5+5/24*z^4-5/3*z^3+5*z^2-5*z+1)
          return(y)
        }
      }
	else if (k==6) {
        h<-function(x){
          z<-x/sigma
          y<-1/720*z^6-1/20*z^5+5/8*z^4-10/3*z^3+15/2*z^2-6*z+1
          return(y)
        }
      }
	else if (k==7) {
        h<-function(x){
          z<-x/sigma
          y<--(-1/5040*z^7+7/720*z^6-7/40*z^5+35/24*z^4-35/6*z^3+21/2*z^2-7*z+1)
          return(y)
        }
      }
	else if (k==8) {
        h<-function(x){
          z<-x/sigma
          y<-1/40320*z^8-1/630*z^7+7/180*z^6-7/15*z^5+35/12*z^4-28/3*z^3+14*z^2-8*z+1
          return(y)
        }
      }
	else if (k==9) {
        h<-function(x){
          z<-x/sigma
          y<--(-1/362880*z^9+1/4480*z^8-1/140*z^7+7/60*z^6-21/20*z^5+21/4*z^4-14*z^3+18*z^2-9*z+1)
          return(y)
        }
      }
	else if (k==10) {
        h<-function(x){
          z<-x/sigma
          y<-1/3628800*z^10-1/36288*z^9+1/896*z^8-1/42*z^7+7/24*z^6-21/10*z^5+35/4*z^4-20*z^3+45/2*z^2-10*z+1
          return(y)
        }
      }
	return(h)
}

orthpol_pois<-function(k,lambda) {
      n<-k
	h<-function(x){
	  y<-0
        for (i in 0:n) {
	    y<-y + ((-1)^(n-i)*factorial(i)*lambda^(-i)*choose(n,i)*choose(x,i))
	  }
	  y<-y*sqrt(lambda^n/factorial(n))
	  return(y)
	}
	return(h)
}

orthpol_unif<-function(k,a,b) {
      n<-k
	h<-function(x){
	  z<-(x-a)/(b-a)
	  y<-0
        for (i in 0:n) {
	    y<-y + ((-1)^(n)*(-z)^(i)*choose(n,i)*choose(n+i,i))
	  }
	  y<-y*sqrt(2*n+1)
	  return(y)
	}
	return(h)
}

orthpol_logis<-function(k,mu,sigma) {
	#location-scale family!
	if (k==0) {
        h<-function(x){
          y<-1
          return(y)
        }
      }
	else if (k==1) {
        h<-function(x){
          z<-(x-mu)/sigma
          y<-sqrt(3)/pi*z
          return(y)
        }
      }
	else if (k==2) {
        h<-function(x){
          z<-(x-mu)/sigma
          y<-3/sqrt(3.2)/pi^2*(z^2-pi^2/3)
          return(y)
        }
      }
	else if (k==3) {
        h<-function(x){
          z<-(x-mu)/sigma
          y<-5*sqrt(7)/12/pi^3*(z^3-4.2*pi^2/3*z)
          return(y)
        }
      }
	else if (k==4) {
        h<-function(x){
          z<-(x-mu)/sigma
          y<-35/64/pi^4*(z^4-26/7*pi^2*z^2+27/35*pi^4)
          return(y)
        }
      }
	return(h)
}

orthpol_ZIP<-function(k,p,lambda) {
	if (k==0) {
        h<-function(x){
          y<-1
          return(y)
        }
      }
	else if (k==1) {
        h<-function(x){
          y<-(x+(-1+p)*lambda)/(-lambda^2*p^2+lambda-lambda*p+lambda^2*p)^(1/2)
          return(y)
        }
      }
	else if (k==2) {
        h<-function(x){
          y<-(x^2*(lambda*p+1)-x*(lambda^2*p+lambda*p+2*lambda+1)-lambda^2*p+lambda^2)/((1-p)*(2*lambda*p+lambda^2*p+2)*lambda^2*(lambda*p+1))^(1/2)
          return(y)
        }
      }
	else if (k==3) {
        h<-function(x){
          y<-1/2*lambda^3*2^(1/2)*((2*lambda*p+2+2*p^2+lambda^2*p+2*lambda*p^3-2*lambda^2*p^2-4*p-4*lambda*p^2+lambda^2*p^3)*x^3+(6*lambda*p+12*p+4*lambda^3*p^2-6*p^2-9*lambda^2*p-6+18*lambda^2*p^2+6*lambda*p^2-6*lambda-9*lambda^2*p^3-2*lambda^3*p-2*lambda^3*p^3-6*lambda*p^3)*x^2+(p*lambda^4-12*lambda^3*p^2-2*lambda^4*p^2+6*lambda+6*lambda^3*p+4+4*p^2-2*lambda*p^2+8*lambda^2*p^3+lambda^4*p^3+6*lambda^2+4*lambda*p^3+6*lambda^3*p^3-4*lambda^2*p-10*lambda^2*p^2-8*p-8*lambda*p)*x-2*lambda^3-6*lambda^3*p^2+6*lambda^3*p+2*lambda^3*p^3)/(-(-1+p)^5*lambda^9*(2*lambda*p+lambda^2*p+2)*(3*lambda^2*p+lambda^3*p+6*lambda*p+6))^(1/2)
          return(y)
        }
      }
	else if (k==4) {
        h<-function(x){
          y<-1/6*((3*lambda^2*p+lambda^3*p+6*lambda*p+6)*x^4+(-42*lambda^2*p-24*lambda-3*p*lambda^4-18*lambda^3*p-36*lambda*p-36)*x^3+(72*lambda+66+83*lambda^3*p+66*lambda*p+36*lambda^2+3*p*lambda^5+24*p*lambda^4+105*lambda^2*p)*x^2+(-66*lambda^2*p-24*lambda^3-p*lambda^6-36*lambda*p-9*p*lambda^5-48*lambda-36*lambda^2-36-39*p*lambda^4-66*lambda^3*p)*x-6*p*lambda^4+6*lambda^4)*6^(1/2)/((1-p)*lambda^4*(4*lambda^3*p+24*lambda*p+p*lambda^4+12*lambda^2*p+24)*(3*lambda^2*p+lambda^3*p+6*lambda*p+6))^(1/2)
          return(y)
        }
      }
	return(h)
}

orthpol_laplace<-function(k,alpha,beta) {
	#location-scale family!
	if (k==0) {
        h<-function(x){
          y<-1
          return(y)
        }
      }
	else if (k==1) {
        h<-function(x){
          z<-(x-alpha)/beta
          y<-z/sqrt(2)
          return(y)
        }
      }
	else if (k==2) {
        h<-function(x){
          z<-(x-alpha)/beta
          y<-(z^2-2)/(2*sqrt(5))
          return(y)
        }
      }
	else if (k==3) {
        h<-function(x){
          z<-(x-alpha)/beta
          y<-(z^3-12*z)/(12*sqrt(3))
          return(y)
        }
      }
	else if (k==4) {
        h<-function(x){
          z<-(x-alpha)/beta
          y<-1/3576*sqrt(745)*(z^4-168/5*z^2+216/5)
          return(y)
        }
      }
	return(h)
}

orthpol_logar<-function(k,p) {
	if (k==0) {
        h<-function(x){
          y<-1
          return(y)
        }
      }
	else if (k==1) {
        h<-function(x){
          y<-(x*ln(1-p)*p-x*ln(1-p)-p)/ln(1-p)/(p-1)/(-(p+ln(1-p))*p/ln(1-p)^2/(p-1)^2)^(1/2)
          return(y)
        }
      }
	else if (k==2) {
        h<-function(x){
          y<-(-p^2+x^2*ln(1-p)-x*ln(1-p)+x*ln(1-p)*p^2+x*p^2-x*p+x^2*p^3+x^2*ln(1-p)*p^2-2*x^2*p^2-2*x^2*ln(1-p)*p+x^2*p)/(p-1)^2/(p+ln(1-p))/(-p^2/ln(1-p)*(2*p+p^2+2*ln(1-p))/(p-1)^4/(p+ln(1-p)))^(1/2)
          return(y)
        }
      }
	else if (k==3) {
        h<-function(x){
          y<-1/2*p^3/ln(1-p)^3*(p^5*x^2+7*p^4*x^2-6*x^2*ln(1-p)*p^2-6*x^2*ln(1-p)*p-4*x*p+6*ln(1-p)*p^3*x^2+6*p^3*x-2*ln(1-p)*x^3-2*x*p^2+4*x*ln(1-p)*p^3-4*x*ln(1-p)-4*p^3-3*x^3*p^3+5*x^3*p^2-2*x^3*p-x^3*p^4+p^5*x^3-6*ln(1-p)*x^3*p^2-11*x^2*p^3+6*ln(1-p)*x^3*p-3*x^2*p^2+6*x^2*p+2*ln(1-p)*x^3*p^3+6*x^2*ln(1-p))/(-p^9/ln(1-p)^7*(2*p+p^2+2*ln(1-p))*(6*p+2*p^3+3*p^2+6*ln(1-p))/(p-1)^18)^(1/2)/(p-1)^9
          return(y)
        }
      }
	else if (k==4) {
        h<-function(x){
          y<-1/12*(26*p^3*x^4-21*p^2*x^4+6*p*x^4-14*p^4*x^4+36*ln(1-p)*p^2*x^4-24*ln(1-p)*p*x^4+6*ln(1-p)*p^4*x^4-24*ln(1-p)*p^3*x^4-5*p^6*x^4+2*p^7*x^4+6*ln(1-p)*x^4+6*p^5*x^4+36*ln(1-p)*p^4*x^3+72*ln(1-p)*p*x^3-72*ln(1-p)*p^3*x^3-36*p*x^3+24*p^3*x^3+54*p^2*x^3-18*p^5*x^3-36*ln(1-p)*x^3-42*p^4*x^3+6*p^7*x^3+12*p^6*x^3+84*p^5*x^2+66*ln(1-p)*x^2-118*p^4*x^2-48*ln(1-p)*p*x^2+66*ln(1-p)*p^4*x^2+17*p^6*x^2-36*ln(1-p)*p^2*x^2-48*ln(1-p)*p^3*x^2+4*p^7*x^2-15*p^2*x^2+66*p*x^2-38*p^3*x^2+36*ln(1-p)*p^4*x-36*ln(1-p)*x-18*p^2*x+66*p^4*x-12*p^3*x-36*p*x-36*p^4)*2^(1/2)/(-p^4*(12*ln(1-p)+12*p+4*p^3+3*p^4+6*p^2)/(6*p+2*p^3+3*p^2+6*ln(1-p))/(p-1)^8/ln(1-p))^(1/2)/(p-1)^4/(6*p+2*p^3+3*p^2+6*ln(1-p))
          return(y)
        }
      }
	return(h)
}
	       
orthpol_negb<-function(k,p,r) {
	if (k==0) {
        h<-function(x){
          y<-1
          return(y)
        }
      }
	else if (k==1) {
        h<-function(x){
          y<-(p*x-r+p*r)/p/(-(-1+p)*r/p^2)^(1/2)
          return(y)
        }
      }
	else if (k==2) {
        h<-function(x){
          y<-1/2*(2*p^2*r*x+r^2*p^2+p^2*x+p^2*x^2+p^2*r-2*p*x-2*p*r*x-2*p*r-2*r^2*p+r+r^2)/p^2*2^(1/2)/((-1+p)^2*(1+r)*r/p^4)^(1/2)
          return(y)
        }
      }
	else if (k==3) {
        h<-function(x){
          y<--1/6*6^(1/2)/(-r^5*(-1+p)^9*(1+r)^3*(2+r)/p^18)^(1/2)*r^2*(2*p^6*x+30*p^4*x-24*p^5*r^3-40*p^3*r-100*p^3*r^2+2*r+24*p^2*x+30*p^4*r+r^4-6*r^4*p+15*r^4*p^2-20*r^4*p^3+15*r^4*p^4+4*r^3-80*r^3*p^3+5*r^2+66*p^2*r*x-15*p*r*x-116*p^3*r*x+2*r*p^6+5*r^2*p^6+60*r^3*p^2-24*r^3*p+75*r^2*p^4+60*r^3*p^4+30*p^2*r+75*r^2*p^2-30*r^2*p-12*p*r+9*p^2*r*x^2-33*p^3*r*x^2-12*p*r^2*x+57*p^2*r^2*x-108*p^3*r^2*x+4*r^3*p^6+r^4*p^6-p^3*x^3+3*p^4*r*x^3-12*p^3*r^2*x^2+6*p^2*x^2-21*p^3*x^2-38*p^3*x+18*p^4*r^2*x^2+15*p^2*r^3*x-30*p^3*r^3*x+30*p^4*r^3*x+3*p^2*r^2*x^2-3*p*r^3*x+45*p^4*r*x^2-12*p^5*r-30*p^5*r^2-6*p^5*r^4-3*p^5*x^3+3*p^4*x^3-27*p^5*x^2*r+102*p^4*r^2*x-48*p^5*x*r^2-3*p^5*r*x^3+p^6*x^3-12*p^5*r^2*x^2-15*p^5*r^3*x+6*p^6*x^2*r+9*p^6*x*r^2+p^6*r*x^3+3*p^6*r^2*x^2+3*p^6*r^3*x-r*p^3*x^3+102*p^4*r*x-45*p^5*x*r-15*p^5*x^2+27*p^4*x^2+3*p^6*x^2-12*p^5*x-6*p*x+8*p^6*x*r)/p^9
          return(y)
        }
      }
	else if (k==4) {
        h<-function(x){
          y<-1/12/p^4*(p^4*x^4+6*p^4*x^3-12*p^3*x^3-4*r*p^3*x^3+4*p^4*r*x^3+30*p^2*r*x^2+6*p^2*r^2*x^2-48*p^3*r*x^2+11*p^4*x^2-36*p^3*x^2+18*p^4*r*x^2-12*p^3*r^2*x^2+36*p^2*x^2+6*p^4*r^2*x^2-4*p*r^3*x-24*p^3*x+12*p^2*r^3*x+6*p^4*x+18*p^4*r^2*x+102*p^2*r*x-44*p*r*x-80*p^3*r*x-24*p*x+4*p^4*r^3*x-24*p*r^2*x+66*p^2*r^2*x-12*p^3*r^3*x+36*p^2*x+22*p^4*r*x-60*p^3*r^2*x+r^4*p^4-4*r^4*p^3-24*r^3*p+6*r^3+6*r^3*p^4+36*p^2*r+66*r^2*p^2-44*r^2*p-24*p*r-24*r^3*p^3+11*r^2-24*p^3*r-44*p^3*r^2+6*r+36*r^3*p^2+6*p^4*r+11*r^2*p^4+r^4-4*r^4*p+6*r^4*p^2)*6^(1/2)/((-1+p)^4*r*(r^3+11*r+6*r^2+6)/p^8)^(1/2)
          return(y)
        }
      }
	return(h)
}

orthpol_extrval<-function(k,a,b) {
	#location-scale family!
	if (k==0) {
        h<-function(x){
          y<-1
          return(y)
        }
      }
	else if (k==1) {
        h<-function(x){
          z<-(x-a)/b
          y<-0.7796968012*z-0.4500532075
          return(y)
        }
      }
	else if (k==2) {
        h<-function(x){
          z<-(x-a)/b
          y<-.3451996482*(z-.5772156649)^2-.5045182395*z-.2766148303
          return(y)
        }
      }
	else if (k==3) {
        h<-function(x){
          z<-(x-a)/b
          y<-.1060499473*(z-.5772156649)^3-.4944037009*(z-.5772156649)^2-.2194200910*z+.6849580626
          return(y)
        }
      }
	else if (k==4) {
        h<-function(x){
          z<-(x-a)/b
          y<-.2493263979e-1*(z-.5772156649)^4-.2416834756*(z-.5772156649)^3+.2690771450*(z-.5772156649)^2+.7769092062*z-.6743236152
          return(y)
        }
      }
	return(h)
}

orthpol_geom<-function(k,p) {
	if (k==0) {
        h<-function(x){
          y<-1
          return(y)
        }
      }
	else if (k==1) {
        h<-function(x){
          y<-(x*p-1+p)/(-(-1+p))^(1/2)
          return(y)
        }
      }
	else if (k==2) {
        h<-function(x){
          y<-1/2*(3*p^2*x+p^2*x^2+2*p^2-4*p*x-4*p+2)/((-1+p)^2)^(1/2)
          return(y)
        }
      }
	else if (k==3) {
        h<-function(x){
          y<--1/6/(-(-1+p)^9)^(1/2)*(6-36*p+90*p^2-120*p^3+90*p^4+9*x^2*p^2-36*p^5+6*p^6-3*p^5*x^3+45*p^4*x^2-60*p^5*x+11*p^6*x+6*p^6*x^2+132*p^4*x-27*p^5*x^2+p^6*x^3-33*x^2*p^3+81*x*p^2-146*x*p^3-x^3*p^3+3*x^3*p^4-18*x*p)
          return(y)
        }
      }
	else if (k==4) {
        h<-function(x){
          y<-1/24*(p^4*x^4+10*p^4*x^3-16*p^3*x^3-96*p^3*x^2+35*p^4*x^2+72*p^2*x^2+216*p^2*x+50*p^4*x-176*p^3*x-96*p*x+24*p^4+144*p^2-96*p^3-96*p+24)/((-1+p)^4)^(1/2)
          return(y)
        }
      }
	return(h)
}

orthpol_ZTP<-function(k,lambda) {
	if (k==0) {
        h<-function(x){
          y<-1
          return(y)
        }
      }
	else if (k==1) {
        h<-function(x){
          y<--(-x*exp(lambda)+x+lambda*exp(lambda))/(exp(lambda)-1)/(-lambda*exp(lambda)*(-exp(lambda)+1+lambda)/(exp(lambda)-1)^2)^(1/2)
          return(y)
        }
      }
	else if (k==2) {
        h<-function(x){
          y<--(lambda^2*exp(lambda)-2*x*lambda*exp(lambda)-x*exp(lambda)+x^2*exp(lambda)+x-x^2-x^2*lambda+lambda^2*x+3*lambda*x)/((-2*exp(lambda)+lambda^2+2+2*lambda)*lambda^2*exp(lambda)/(-exp(lambda)+1+lambda)/(exp(lambda)-1))^(1/2)/(-exp(lambda)+1+lambda)
          return(y)
        }
      }
	else if (k==3) {
        h<-function(x){
          y<--1/2*(-4*x*exp(lambda)+6*x^2*exp(lambda)+2*lambda^3*exp(lambda)-6*x*lambda*exp(lambda)-6*x*lambda^2*exp(lambda)-2*x^3*exp(lambda)+6*x^2*lambda*exp(lambda)+10*lambda*x+6*lambda^3*x-9*lambda^2*x^2-2*lambda^3*x^2+4*x+lambda^2*x^3-12*x^2*lambda-6*x^2+lambda^4*x+2*lambda*x^3+14*lambda^2*x+2*x^3)*exp(2*lambda)*lambda^3*2^(1/2)/(exp(lambda)-1)^3/(lambda^9*exp(lambda)^5*(-2*exp(lambda)+lambda^2+2+2*lambda)*(-6*exp(lambda)+6+3*lambda^2+6*lambda+lambda^3)/(exp(lambda)-1)^7)^(1/2)
          return(y)
        }
      }
	else if (k==4) {
        h<-function(x){
          y<-1/6*(lambda^3*x^4+6*lambda*x^4+3*lambda^2*x^4+6*x^4-6*exp(lambda)*x^4-36*x^3+24*lambda*exp(lambda)*x^3-3*lambda^4*x^3-42*lambda^2*x^3-60*lambda*x^3-18*lambda^3*x^3+36*exp(lambda)*x^3+3*lambda^5*x^2+24*lambda^4*x^2-36*lambda^2*exp(lambda)*x^2+66*x^2-72*exp(lambda)*lambda*x^2+141*lambda^2*x^2+83*lambda^3*x^2+138*lambda*x^2-66*exp(lambda)*x^2-90*lambda^3*x+36*exp(lambda)*lambda^2*x-9*lambda^5*x+24*lambda^3*exp(lambda)*x-lambda^6*x-84*lambda*x-36*x-39*lambda^4*x+36*exp(lambda)*x-102*lambda^2*x+48*exp(lambda)*lambda*x-6*lambda^4*exp(lambda))*6^(1/2)/(-6*exp(lambda)+6+lambda^3+6*lambda+3*lambda^2)/((-24*exp(lambda)+24+12*lambda^2+24*lambda+4*lambda^3+lambda^4)*lambda^4*exp(lambda)/(exp(lambda)-1)/(-6*exp(lambda)+6+lambda^3+6*lambda+3*lambda^2))^(1/2)
          return(y)
        }
      }
	return(h)
}

orthpol_genpareto<-function(k.degree,s,k) { #s is scale and k is shape
	if (k.degree==0) {
        h<-function(x){
          y<-1
          return(y)
        }
      }
	else if (k.degree==1) {
        h<-function(x){
          y<- -(-x-x*k+s)/(1+k)/(s^2/(1+k)^2/(1+2*k))^(1/2)
          return(y)
        }
      }
	else if (k.degree==2) {
        h<-function(x){
          y<- 1/2*(6*x^2*k^2+5*x^2*k-8*x*k*s+x^2+2*s^2-4*x*s)/(1+2*k)/(1+3*k)/(s^4/(1+4*k)/(1+3*k)^2/(1+2*k)^2)^(1/2)    
          return(y)
        }
      }
	else if (k.degree==3) {
        h<-function(x){
          y<--1/6/(s^6/(1+6*k)/(1+5*k)^2/(1+4*k)^4/(1+3*k)^6/(1+2*k)^6/(1+k)^4)^(1/2)*(-60*x^3*k^3+108*x^2*k^2*s-47*x^3*k^2-12*x^3*k-54*x*k*s^2+63*x^2*k*s+9*x^2*s-x^3+6*s^3-18*x*s^2)/(1+5*k)/(1+4*k)^2/(1+3*k)^3/(1+2*k)^3/(1+k)^2  
          return(y)
        }
      }
	else if (k.degree==4) {
        h<-function(x){
          y<- 1/24*(840*k^4*x^4+x^4+638*k^3*x^4+179*k^2*x^4+22*k*x^4-1184*k^2*s*x^3-1920*k^3*s*x^3-240*k*s*x^3-16*s*x^3+1440*k^2*s^2*x^2+72*s^2*x^2+648*k*s^2*x^2-384*k*s^3*x-96*s^3*x+24*s^4)/(1+7*k)/(1+6*k)/(1+5*k)/(1+4*k)/(s^8/(1+8*k)/(1+7*k)^2/(1+6*k)^2/(1+5*k)^2/(1+4*k)^2)^(1/2)
          return(y)
        }
      }
	return(h)
}

orthpol_betab<-function(k,ntrials,a,b) {
	n<-ntrials
	if (k==0) {
        h<-function(x){
          y<-1
          return(y)
        }
      }
	else if (k==1) {
        h<-function(x){
          y<--(-a*x-x*b+a*n)/(a+b)/(n/(a+b)^2/(a+b+1)*a*b*(a+b+n))^(1/2)
          return(y)
        }
      }
	else if (k==2) {
        h<-function(x){
          y<-1/2*(-a^2*n+a^2*n^2+a^2*x^2-2*a^2*x*n+a^2*x-4*a*x*n+a*n^2+3*a*x^2+2*a*x^2*b-2*x*b*a*n+a*x-a*n+x^2*b^2-2*x*n*b+2*x^2-2*x*n+3*x^2*b-x*b-b^2*x)/(a+b+1)/(a+b+2)*2^(1/2)/(b*(1+b)*(1+a)*a*(a*n-a+n*b-b+n^2-1)*n*(a+b+n)/(a+b+2)^2/(a+b)/(a+b+3)/(a+b+1)^2)^(1/2)
          return(y)
        }
      }
	else if (k==3) {
        h<-function(x){
          y<--1/6*6^(1/2)/((1+b)^3*(b+2)*b^5*(a+2)*(1+a)^3*a^5*(n^2+a*n+n*b-2*a-2*b-4)*(a*n-a+n*b-b+n^2-1)^3*(a+b+n)^5*n^5/(a+b+5)/(a+b+4)^2/(a+b+3)^4/(a+b+2)^6/(a+b)^7/(a+b+1)^7)^(1/2)*(a+b+n)^2*b^2*a^2*n^2*(25*a^2*n^2-18*n*a^3+16*a*b*n^2+10*a*b^2*n^2+50*a^2*n^2*b-10*a^2*n*b^2-26*a^3*n*b-4*n*a*b^2-24*n*b*a^2+5*a^5*n^2-4*a^5*n^3-8*a^3*n*b^2-12*a^4*n*b+5*a^4*n^2*b^2+5*a^5*n^2*b-4*a^5*n^3*b-2*a^5*n*b-2*b^2*a^4*n+28*b*a^4*n^2+b^2*a^4*n^4+5*a^2*b^2*n^4-4*a^4*b^2*n^3+2*a^4*n^4*b+a^4*n^4+4*a^3*b^2*n^4-3*a^3*b*n^4-2*b^5*x*n+20*a^3*b^2*n^2-19*a^4*n^3*b-16*a^3*b^2*n^3-3*b^5*x^2+51*a^2*x^2-51*x^2*b^2+40*a^3*b^2*x^3-20*a^2*b^2*n^3-8*a*b^2*n^3+2*a*n^3-3*a^2*n^3+2*x*b^5+18*x*b^3+37*a^3*n^2-10*n^2*x^3*b^3+36*n^3*x^2*b^2+26*b*x*n^2+25*a^2*b^2*n^2+9*n^2*b^4*x^2-8*a^2*b*n^4-4*a*b*n^4-12*a^4*n^3*x-10*a^4*n*x^3+21*a^4*n^2*x^2-7*a^3*n^4-13*a^2*n^4+12*b^4*x^2*n-76*a^2*x*n-20*b^3*x*n^2+2*b^2*x*n^2+b^5*x^3-32*a^3*b*n^3+57*a^3*b*n^2+150*a*n^2*x^2*b-100*a*n*x^3*b-18*x^2*b+14*b^2*x-16*a^3*n^3-4*a*n-36*x^2*n-10*n*b^4*x^3-23*a^2*b*n^3+6*a*n^2-14*a^2*n-204*a^2*b^2*x^2*n-15*a^2*n^3*b^3*x+85*a^2*x^3+51*a^3*x^2+18*a^3*x-6*a*b*n^3+66*n^3*x^2*b-35*n^2*x^3*b^2+2*a^5*x*b+9*a^5*x*n^2+4*a*x+4*x*b-51*b^3*x^2-21*b^4*x^2-3*a^2*b^4*x*n^2+23*a^4*n^2-2*a^5*n+11*b^4*x^3-3*a^5*n^3*x*b-8*a^5*x*n*b+3*b^5*x^2*n-10*a^4*n-15*a^4*n^3-12*n^4*x+4*a^3*n^5*b-50*n*b^2*x^3+a^4*n^5*b+14*b^3*x*n-9*a^3*b^3*x^2*n-6*n^4*b^2*x+15*a*b^4*x^2*n-11*a*b^4*x*n^2-a*n^2*x^3*b^4-3*a*b^5*x^2+3*a^5*n^2*x^2*b-a^5*n*x^3*b+4*a^3*n^5+a^4*n^5+5*a^2*n^5+a^5*b*n^4+15*a*b^4*x^3-27*a*b^4*x^2+12*a*x*b^4-8*b^4*x*n^2+163*a^3*b*x*n^2+30*a^3*n^3*x^2*b+2*n^5*a+a^5*n^4+12*x*n^2+2*a*n^5*b-35*n*x^3*b^3+54*n^2*x^2*b^3+2*a*n^4*b^2-a*n*b^5*x^3-129*a^2*x*n*b-129*a^2*n*x^3*b+174*a^2*n^2*x^2*b-50*n^2*x^3*b-13*a*n^2*x^3*b^3+12*a*n^2*b^4*x^2-20*a*b^3*x*n^2-14*a*n*b^4*x^3+10*x*b^4+5*a^2*n^5*b-2*b^4*x*n+109*a*b*x*n^2+3*a*b^5*x^2*n+141*a*n^3*x^2*b+74*b*x^3+85*b^2*x^3-129*a*n*b^2*x^3-94*a*n^2*x^3*b+24*x^3-8*a*b*n+26*b^2*x*n-51*a*n^4*x*b-3*a^3*n^4*b^2*x-6*a^3*n*x^3*b^3+9*a^3*n^2*x^2*b^3+14*a^2*x-78*a*n^3*x*b^2+45*x^3*b^3-24*a*x*n-6*a*n^4-3*a^3*n^3*b^3*x-93*a^3*b^2*x^2*n-3*a^2*n^2*x^3*b^3-n^2*x^3*b^4+213*a*n^2*x^2*b^2-291*a*b*x^2*n+a^5*x^3+185*a^2*x^3*b-126*a^3*x^2*n+185*a*x^3*b^2+97*a^3*x*n^2+4*a^4*b^2*x^3-24*n^2*x^3-24*n^3*x*b-36*n^3*x*b^2+2*a^5*x+3*a^5*x^2+6*n^3*x^2*b^3-15*a^4*b^2*x^2*n+8*a^2*x*b^3-18*n^4*x*b-a*b^4*x*n+9*a*n^3*x^2*b^3-27*a*b^3*x^2*n+79*a*x^3*b^3-6*a^2*b^4*x^2+2*a^2*x*b^4+40*a^2*x^3*b^3-345*a^2*x^2*b*n+36*n^3*x^2+3*a^3*b^3*x*n^2+45*a^3*x^3+9*a^5*b*x*n^2+18*a*x^2-12*n^3*b^3*x+24*a^3*x^2*b^2+8*a^3*b^2*x+79*a^3*b*x^3+3*a^2*n^2*b^4*x^2+3*a^2*b^3*x*n^2-4*a^2*n*b^4*x^3+a*b^5*x^3+6*a^4*x^2*b^2+2*a^4*b^2*x+3*a^2*n^3*x^2*b^3-27*a^2*b^3*x^2*n+4*a^2*b^3*x*n-12*a^2*n^4*b^2*x-37*a^2*n*x^3*b^3+48*a^2*n^2*x^2*b^3+a^2*b^4*x*n+12*a*n^3*x+64*a^4*b*x*n^2+3*a^4*n^3*x^2*b-4*a^4*n*b^2*x^3-88*a^3*x*n-9*a^3*n^3*x+27*a^4*x^2*b-45*a^4*x^2*n+12*a^4*x*b-24*a*n*x^3+18*a*n^2*x^2-3*a^4*n^4*x+49*a^4*x*n^2-a^4*n^2*x^3+3*a^4*n^3*x^2+51*a^3*n^2*x^2-35*a^3*n*x^3+15*a^4*b*x^3-a^4*n^2*x^3*b-11*a^4*b^2*x*n-3*a^4*n^4*x*b+9*a^4*n^2*x^2*b^2			+a^5*b*x^3-32*x*b*a*n-171*a^2*x^2*n+69*a^2*x^2*b+28*a^2*x*b+102*a^2*n^3*x^2*b-112*a^2*n*b^2*x^3-56*a^2*n^2*x^3*b-50*a^2*b^2*x*n-66*a^2*n^3*x*b^2-51*a^2*n^4*x*b-39*a^2*n^4*x+89*a^2*x*n^2-6*a^5*b*x^2*n+69*a^2*n^3*x^2-35*a^2*n^2*x^3-69*a*x^2*b^2+99*n^2*x^2*b^2+28*a*b^2*x+194*a*b*x^3-60*a^4*b*x^2*n+6*a^3*n^3*x^2*b^2+69*a^3*b^2*x*n^2-3*a^3*n^2*x^3*b^2-18*a^4*n^3*x*b-6*a^4*n^3*x*b^2+11*a^4*x^3+6*a^3*x^3*b^3+12*a^2*n^3*x-75*a*b^3*x^2-37*a^3*n*b^2*x^3-13*a^3*n^2*x^3*b-44*a^3*b^2*x*n-30*a^3*n^3*x*b^2-21*a^3*n^4*x*b-24*a^2*b^3*x^2+136*a^2*b^2*x^3+36*a^2*n^3*x^2*b^2+108*a^2*b^2*x*n^2-24*a^2*n^2*x^3*b^2+12*x*n*b-55*a^4*x*n*b+30*a^4*n^2*x^2*b-14*a^4*n*x^3*b-39*a^2*n^3*x*b+2*a*x*b^5-n*b^5*x^3-50*a^2*n*x^3+51*a^2*n^2*x^2-2*a*b^5*x*n+3*a^5*x^2*b-6*a^5*x^2*n-24*a*n^3*b^3*x+74*a*x^3+75*a^3*x^2*b+26*a^3*x*b-18*a^3*n^4*x+63*a^3*n^2*x^2*b^2-210*a^3*b*x^2*n-102*b*x^2*n-84*b^2*x^2*n-9*b^3*x^2*n-10*a^3*n^2*x^3+24*a^3*n^3*x^2-44*a^4*x*n+4*a^2*b^4*x^3-36*a^3*n^3*x*b+54*n^2*x^2*b+10*a^4*x+21*a^4*x^2+191*x*b*a^2*n^2+20*a^2*b^2*x+26*a*x*b^3-120*a*x^2*n+168*a^2*n^2*x^2*b^2+66*a*n^3*x^2*b^2+56*a*b^2*x*n^2-56*a*n^2*x^3*b^2+18*a*b^3*x*n-15*a*n^4*b^2*x-66*a*n*x^3*b^3+93*a*n^2*x^2*b^3-24*n*x^3*b-132*a^3*x*n*b+105*a^3*n^2*x^2*b-66*a^3*n*x^3*b+16*a*x*b+3*a^2*b^4*x^2*n-36*a*n^4*x+44*a*x*n^2-8*a^5*x*n-50*a*n^2*x^3+9*a*b^2*x*n-210*x^2*b^2*a*n-42*a*n^3*x*b-3*a^5*n^3*x-a^5*n*x^3+3*a^5*n^2*x^2+84*a*n^3*x^2+15*a^4*b^2*x*n^2)/(a^6+6*b*a^5+9*a^5+15*b^2*a^4+26*a^4+45*a^4*b+90*a^3*b^2+20*a^3*b^3+24*a^3+104*b*a^3+156*b^2*a^2+72*b*a^2+90*a^2*b^3+15*a^2*b^4+72*a*b^2+104*a*b^3+6*a*b^5+45*b^4*a+b^6+26*b^4+9*b^5+24*b^3)/(a+b+2)^2/(a+b+3)/(a+b+1)^3
          return(y)
        }
      }
	else if (k==4) {
        h<-function(x){
          y<-1/12*(342*b*x^4+119*b^2*x^4+4*a^3*b*x^4+18*a^3*x^4+18*b^3*x^4+54*a^2*b*x^4+54*a*b^2*x^4+a^4*x^4+b^4*x^4+238*a*b*x^4+342*a*x^4+119*a^2*x^4+6*a^2*b^2*x^4+4*a*b^3*x^4+360*x^4+6*a^4*x^3-12*b^3*n*x^3-60*a^3*n*x^3-282*b^2*x^3+72*a^2*b*x^3+12*a^3*b*x^3-12*a^3*n*b*x^3-564*n*b*x^3-72*b^3*x^3-332*a^2*n*x^3+282*a^2*x^3-4*a*b^3*n*x^3-6*b^4*x^3-720*n*x^3+72*a^3*x^3-476*a*n*b*x^3-72*a*b^2*x^3-12*a*b^3*x^3-4*a^4*n*x^3-84*a*b^2*n*x^3+360*a*x^3-144*b^2*n*x^3-132*a^2*n*b*x^3-360*b*x^3-12*a^2*b^2*n*x^3-804*a*n*x^3+229*a^2*x^2+18*a*b^2*x^2+72*a^3*n^2*x^2-180*a^3*n*x^2+198*a*x^2-648*a*n*x^2+8*a^3*b*x^2-594*a^2*n*x^2+198*b*x^2+432*b*n*x^2+36*b^3*n*x^2+12*a*b^3*n*x^2+102*a^2*n^2*b*x^2-18*a^4*n*x^2+8*a*b^3*x^2+229*b^2*x^2+90*b^3*x^2+12*a^3*n^2*b*x^2+318*a^2*n^2*x^2+72*x^2+36*b^2*n^2*x^2+432*n^2*x^2+282*a*n^2*b*x^2+6*a^2*b^2*n*x^2+11*b^4*x^2-6*a^2*b^2*x^2-24*a^3*b*n*x^2+252*n^2*b*x^2+252*b^2*n*x^2+90*a^3*x^2+26*a*b*x^2+6*a^4*n^2*x^2-114*a^2*b*n*x^2+18*a^2*b*x^2+11*a^4*x^2+30*a*b^2*n^2*x^2+102*b^2*a*n*x^2+612*a*n^2*x^2+6*a^2*b^2*n^2*x^2+18*a*b*n*x^2-24*a^2*n^3*b*x-6*b^4*x-132*n*b*x-30*a*b^2*n^2*x-66*b^2*x-156*a^3*n*x-18*a*b^2*n*x-36*b^2*n^2*x+324*a*n^2*x+6*a^2*b^2*n*x-108*b^2*n*x-6*a^2*b^2*n^2*x-36*a^3*n^3*x-24*b^3*n*x-36*b*x-24*n^3*b*x+18*a^4*n^2*x-4*a^4*n^3*x+42*b*a^2*n^2*x+66*a^2*x+12*a^3*b*n^2*x-18*a^2*n*b*x-264*a*n*x-72*n*x-116*a^2*n^3*x-22*a^4*n*x+36*a*x-8*a^3*n*b*x-44*a*n^3*b*x-72*n^3*x-156*a*n^3*x-18*a*b*n^2*x-26*b*a*n*x+144*a^3*n^2*x-4*a^3*n^3*b*x+6*a^4*x-8*a*b^3*n*x-36*b^3*x+378*a^2*n^2*x+36*a^3*x-350*a^2*n*x-108*b*n^2*x+11*a^4*n^2+11*a^2*n^4+a^4*n^4-36*a*n^3-6*a^4*n-6*a^4*n^3+6*a*n^4+66*a*n^2+66*a^3*n^2-36*a^3*n+6*a^3*n^4-36*a^3*n^3+121*a^2*n^2-66*a^2*n-36*a*n-66*a^2*n^3)*6^(1/2)/(b*(6+6*b^2+b^3+11*b)*n*(11*n-6*n^2+n^3-6)*a*(36*a+36*n+36*b+4*a^6*n+4*a^3*n*b^3+12*a^4*b^2*n+198*a*n+12*a^5*n*b+36*n^3+30*a^3*n^3+386*a^2*n+80*a^2*n^3+6*a^3*b^2*n^2+24*a^2*b^3*n+12*a^4*b*n^2+108*a^4*b*n+240*b*n^2*a^2+132*a^2+42*a^5*n+108*n^2*b+90*b^2*a^3*n+6*a^5*n^2+6*b^4+458*b*n*a+24*b^3*n+36*n^2*b^2+24*n^3*b+6*n^4+193*a^3+270*a*n*b^2+386*b*a^2+229*a*n^2+198*a*b+66*b^2+185*a^3*n^2+360*a^3*n+270*b*n^2*a+144*a^4+108*n*b^2+132*b*n+229*a*b^2+90*a*n^3+66*n^2+44*a*n^3*b+66*a*n^2*b^2+600*a^2*n*b+300*a^2*n^2+370*b*n*a^3+4*b*n^3*a^3+90*b*n^2*a^3+58*a^5+24*b*n^3*a^2+12*a^6+44*a*b^3*n+240*a^2*n*b^2+36*a^2*n^2*b^2+6*a^2*n^4+4*a^4*n^3+a^3*n^4+54*a^4*n^2+174*a^4*n+300*a^2*b^2+360*a^3*b+36*b^3+a^7+6*a^5*b^2+90*a*b^3+11*a*n^4+4*a^4*b^3+80*a^2*b^3+54*a^4*b^2+42*a^5*b+185*a^3*b^2+11*a*b^4+a^3*b^4+174*a^4*b+30*a^3*b^3+6*b^4*a^2+4*a^6*b)/(a^2+2*a+2*a*b+2*b+b^2)/(b+7+a)/(a+b+1)/(a+b+3)^2/(a+b+4)^2/(a+b+5)^2/(b+6+a)^2)^(1/2)/(a+b+3)/(a^2+2*a*b+9*a+20+9*b+b^2)/(b+6+a)
          return(y)
        }
      }
	return(h)
}

orthpol_gamma<-function(k,a,b) { # shape a and scale b
	if (k==0) {
        h<-function(x){
          y<-1
          return(y)
        }
      }
	else if (k==1) {
        h<-function(x){
          y<-(x-b*a)/b/a^(1/2)
          return(y)
        }
      }
	else if (k==2) {
        h<-function(x){
          y<-1/2*(x^2-2*x*b*a+b^2*a^2-2*b*x+b^2*a)*2^(1/2)/b^2/a^(1/2)/(a+1)^(1/2)
          return(y)
        }
      }
	else if (k==3) {
        h<-function(x){
          y<-1/6*(x^3-6*b*x^2-3*b*a*x^2+9*b^2*a*x+3*b^2*a^2*x+6*b^2*x-2*b^3*a-3*b^3*a^2-b^3*a^3)*6^(1/2)/(a+2)^(1/2)/(a+1)^(1/2)/a^(1/2)/b^3
          return(y)
        }
      }
	else if (k==4) {
        h<-function(x){
          y<-1/12*(x^4-12*b*x^3-4*b*a*x^3+6*b^2*a^2*x^2+36*b^2*x^2+30*b^2*a*x^2-4*b^3*a^3*x-24*b^3*x-44*b^3*a*x-24*b^3*a^2*x+6*b^4*a^3+11*b^4*a^2+6*b^4*a+b^4*a^4)*6^(1/2)/b^4/a^(1/2)/(a^3+6*a^2+11*a+6)^(1/2)
          return(y)
        }
      }
	return(h)
}


