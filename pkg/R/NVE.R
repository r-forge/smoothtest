########NVE.logis#########

NVE.logis<-function(sample,order,parest=NA){
  #choose order between 3 and 8#
 if (order<9 && order>2) {
  n<-length(sample)
  if (is.na(parest[1])) {
  tmp<-logis.MME(sample)
  beta<-tmp[1]
  sigma<-tmp[2]
  }
  else {
  beta<-parest[1]
  sigma<-parest[2]
  }
  z<-(sample-beta)/sigma
  var_est<-var.logis2(order,z,sigma)
 }
 else {var_est<-NA}
 return(var_est)
}

var.logis2<-function(r,z,sigma) {
  n<-length(z)
  Eh.Dmu<-switch(r-2,
	mean(h3.Dmu.logis2(z,sigma)),
	mean(h4.Dmu.logis2(z,sigma)),
	mean(h5.Dmu.logis2(z,sigma)),
	mean(h6.Dmu.logis2(z,sigma)),
	mean(h7.Dmu.logis2(z,sigma)),
	mean(h8.Dmu.logis2(z,sigma)))
  Eh.Dsigma<-switch(r-2,
	mean(h3.Dsigma.logis2(z,sigma)),
	mean(h4.Dsigma.logis2(z,sigma)),
	mean(h5.Dsigma.logis2(z,sigma)),
	mean(h6.Dsigma.logis2(z,sigma)),
	mean(h7.Dsigma.logis2(z,sigma)),
	mean(h8.Dsigma.logis2(z,sigma)))
  v<-mean((g_r.logis2(r,z,sigma,Eh.Dmu,Eh.Dsigma))^2)-(mean(g_r.logis2(r,z,sigma,Eh.Dmu,Eh.Dsigma)))^2
  v<-v*n/(n-1)
  return(v)
}

g_r.logis2<-function(r,z,sigma,Eh.Dmu,Eh.Dsigma) {
  h<-orth.poly(r,"logis",c(0,1))
  tmp<-h(z)+b.mu.logis2(z,sigma)*Eh.Dmu+b.sigma.logis2(z,sigma)*Eh.Dsigma
  return(tmp)
}

b.mu.logis2<-function(z,sigma) {sigma*z}
b.sigma.logis2<-function(z,sigma) {sigma/2*(3*z^2/pi^2-1)}

h3.Dmu.logis2<-function(z,sigma) -5*sqrt(7)/(12*pi^3)*(3*z^2-4.2*pi^2/3)/sigma
h3.Dsigma.logis2<-function(z,sigma) -5*sqrt(7)/(12*pi^3)*(3*z^2-4.2*pi^2/3)*(z/sigma)

h4.Dmu.logis2<-function(z,sigma) {-35/64/pi^4*(4*z^3-26/7*pi^2*2*z)/sigma}
h4.Dsigma.logis2<-function(z,sigma) {-35/64/pi^4*(4*z^3-26/7*pi^2*2*z)*z/sigma}

h5.Dmu.logis2<-function(z,sigma) {-(0.0007112402*5*z^4 - 0.0545973510*3*z^2 + 0.4475789427)/sigma}
h5.Dsigma.logis2<-function(z,sigma) {-(0.0007112402*5*z^4 - 0.0545973510*3*z^2 + 0.4475789427)*z/sigma}

h6.Dmu.logis2<-function(z,sigma) {-(7.520244e-05*6*z^5 - 1.045853e-02*4*z^3 + 2.190961e-01*2*z)/sigma}
h6.Dsigma.logis2<-function(z,sigma) {-(7.520244e-05*6*z^5 - 1.045853e-02*4*z^3 + 2.190961e-01*2*z)*z/sigma}

h7.Dmu.logis2<-function(z,sigma) {-(6.821868e-06*7*z^6 - 1.558928e-03*5*z^4 + 6.671607e-02*3*z^2 - 4.159354e-01)/sigma}
h7.Dsigma.logis2<-function(z,sigma) {-(6.821868e-06*7*z^6 - 1.558928e-03*5*z^4 + 6.671607e-02*3*z^2 - 4.159354e-01)*z/sigma}

h8.Dmu.logis2<-function(z,sigma) {-(5.418057e-07*8*z^7 - 1.896547e-04*6*z^5 + 1.445543e-02*4*z^3 - 2.248586e-01*2*z)/sigma}
h8.Dsigma.logis2<-function(z,sigma) {-(5.418057e-07*8*z^7 - 1.896547e-04*6*z^5 + 1.445543e-02*4*z^3 - 2.248586e-01*2*z)*z/sigma}

##########NVE.extrval####################

NVE.extrval<-function(sample,order,parest=NA){
  #choose order between 3 and 8#
 if (order<9 && order>2) {
  n<-length(sample)
  if (is.na(parest[1])) {
  tmp<-extrval.MME(sample)
  a<-tmp[1]
  b<-tmp[2]
  }
  else {
  a<-parest[1]
  b<-parest[2]
  }
  z<-(sample-a)/b
  var_est<-var.extrval2(order,z,b)
 }
 else {var_est<-NA}
 return(var_est)
}

var.extrval2<-function(r,z,b) {
  n<-length(z)
  Eh.Da<-switch(r-2,
	mean(h3.Da.extrval2(z,b)),
	mean(h4.Da.extrval2(z,b)),
	mean(h5.Da.extrval2(z,b)),
	mean(h6.Da.extrval2(z,b)),
	mean(h7.Da.extrval2(z,b)),
	mean(h8.Da.extrval2(z,b)))
  Eh.Db<-switch(r-2,
	mean(h3.Db.extrval2(z,b)),
	mean(h4.Db.extrval2(z,b)),
	mean(h5.Db.extrval2(z,b)),
	mean(h6.Db.extrval2(z,b)),
	mean(h7.Db.extrval2(z,b)),
	mean(h8.Db.extrval2(z,b)))
  v<-mean((g_r.extrval2(r,z,b,Eh.Da,Eh.Db))^2)-(mean(g_r.extrval2(r,z,b,Eh.Da,Eh.Db)))^2
  v<-v*n/(n-1)
  return(v)
}

g_r.extrval2<-function(r,z,b,Eh.Da,Eh.Db) {
  h<-orth.poly(r,"extrval",c(0,1))
  tmp<-h(z)+b.a.extrval2(z,b)*Eh.Da+b.b.extrval2(z,b)*Eh.Db
  return(tmp)
}

gam<--digamma(1)
b.a.extrval2<-function(z,b) {b*(-3/pi^2*gam*(z-gam)^2+z-gam/2)}
b.b.extrval2<-function(z,b) {b*(3/pi^2*(z-gam)^2-1/2)}

h3.Da.extrval2<-function(z,b) -(0.1060499473*3*(z - 0.5772156649)^2-0.4944037009*2*(z - 0.5772156649)-0.219420091)/b
h3.Db.extrval2<-function(z,b) -(0.1060499473*3*(z - 0.5772156649)^2-0.4944037009*2*(z - 0.5772156649)-0.219420091)*(z/b)

h4.Da.extrval2<-function(z,b) {-(0.02493263979*4*(z - 0.5772156649)^3-0.2416834756*3*(z - 0.5772156649)^2+0.269077145*2*(z - 0.5772156649)+0.7769092062)/b}
h4.Db.extrval2<-function(z,b) {-(0.02493263979*4*(z - 0.5772156649)^3-0.2416834756*3*(z - 0.5772156649)^2+0.269077145*2*(z - 0.5772156649)+0.7769092062)*z/b}

h5.Da.extrval2<-function(z,b) {-(-0.765066730-0.462335263*2*z+0.480821661*3*z^2-0.092490112*4*z^3+0.004746706*5*z^4)/b}
h5.Db.extrval2<-function(z,b) {-(-0.765066730-0.462335263*2*z+0.480821661*3*z^2-0.092490112*4*z^3+0.004746706*5*z^4)*z/b}

h6.Da.extrval2<-function(z,b) {-(1.0327007209-0.0480989854*2*z-0.5440920752*3*z^2+0.1947556972*4*z^3-0.0219166650*5*z^4+0.0007592173*6*z^5)/b}
h6.Db.extrval2<-function(z,b) {-(1.0327007209-0.0480989854*2*z-0.5440920752*3*z^2+0.1947556972*4*z^3-0.0219166650*5*z^4+0.0007592173*6*z^5)*z/b}

h7.Da.extrval2<-function(z,b) {-(-0.9614699721+0.6159510209*2*z+0.4265597373*3*z^2-0.2987314444*4*z^3+0.0569116488*5*z^4-0.0042051687*6*z^5+0.0001046942*7*z^6)/b}
h7.Db.extrval2<-function(z,b) {-(-0.9614699721+0.6159510209*2*z+0.4265597373*3*z^2-0.2987314444*4*z^3+0.0569116488*5*z^4-0.0042051687*6*z^5+0.0001046942*7*z^6)*z/b}

h8.Da.extrval2<-function(z,b) {-(6.127549e-01-1.061618e+00*2*z-1.309572e-01*3*z^2+3.594015e-01*4*z^3-1.081713e-01*5*z^4+1.297077e-02*6*z^5-6.774490e-04*7*z^6+1.268790e-05*8*z^7)/b}
h8.Db.extrval2<-function(z,b) {-(6.127549e-01-1.061618e+00*2*z-1.309572e-01*3*z^2+3.594015e-01*4*z^3-1.081713e-01*5*z^4+1.297077e-02*6*z^5-6.774490e-04*7*z^6+1.268790e-05*8*z^7)*z/b}

###############NVE.laplace#######################

NVE.laplace<-function(sample,order,parest=NA){
  #choose order between 3 and 8#
 if (order<9 && order>2) {
  n<-length(sample)
  if (is.na(parest[1])) {
  tmp<-laplace.MME(sample)
  a<-tmp[1]
  b<-tmp[2]
  }
  else {
  a<-parest[1]
  b<-parest[2]
  }
  z<-(sample-a)/b
  var_est<-var.laplace2(order,z,b)
 }
 else {var_est<-NA}
 return(var_est)
}

var.laplace2<-function(r,z,b) {
  n<-length(z)
  Eh.Da<-switch(r-2,
	mean(h3.Da.laplace2(z,b)),
	mean(h4.Da.laplace2(z,b)),
	mean(h5.Da.laplace2(z,b)),
	mean(h6.Da.laplace2(z,b)),
	mean(h7.Da.laplace2(z,b)),
	mean(h8.Da.laplace2(z,b)))
  Eh.Db<-switch(r-2,
	mean(h3.Db.laplace2(z,b)),
	mean(h4.Db.laplace2(z,b)),
	mean(h5.Db.laplace2(z,b)),
	mean(h6.Db.laplace2(z,b)),
	mean(h7.Db.laplace2(z,b)),
	mean(h8.Db.laplace2(z,b)))
  v<-mean((g_r.laplace2(r,z,b,Eh.Da,Eh.Db))^2)-(mean(g_r.laplace2(r,z,b,Eh.Da,Eh.Db)))^2
  v<-v*n/(n-1)
  return(v)
}

g_r.laplace2<-function(r,z,b,Eh.Da,Eh.Db) {
  h<-orth.poly(r,"laplace",c(0,1))
  tmp<-h(z)+b.a.laplace2(z,b)*Eh.Da+b.b.laplace2(z,b)*Eh.Db
  return(tmp)
}

b.a.laplace2<-function(z,b) {b*z}
b.b.laplace2<-function(z,b) {b/4*(z^2-2)}

h3.Da.laplace2<-function(z,b) {-(0.04811252*3*z^2-0.57735027)/b}
h3.Db.laplace2<-function(z,b) {-(0.04811252*3*z^2-0.57735027)*z/b}

h4.Da.laplace2<-function(z,b) {-(0.007632743*4*z^3-0.256460157*2*z)/b}
h4.Db.laplace2<-function(z,b) {-(0.007632743*4*z^3-0.256460157*2*z)*z/b}

h5.Da.laplace2<-function(z,b) {-(0.0009775774*5*z^4-0.0716890063*3*z^2+0.5083402265)/b}
h5.Db.laplace2<-function(z,b) {-(0.0009775774*5*z^4-0.0716890063*3*z^2+0.5083402265)*z/b}

h6.Da.laplace2<-function(z,b) {-(0.0001035120*6*z^5-0.0139011776*4*z^3+0.2658522066*2*z)/b}
h6.Db.laplace2<-function(z,b) {-(0.0001035120*6*z^5-0.0139011776*4*z^3+0.2658522066*2*z)*z/b}

h7.Da.laplace2<-function(z,b) {-(9.446357e-06*7*z^6-2.111131e-03*5*z^4+8.604678e-02*3*z^2-4.629928e-01)/b}
h7.Db.laplace2<-function(z,b) {-(9.446357e-06*7*z^6-2.111131e-03*5*z^4+8.604678e-02*3*z^2-4.629928e-01)*z/b}

h8.Da.laplace2<-function(z,b) {-(7.506133e-07*8*z^7-2.578817e-04*6*z^5+1.894136e-02*4*z^3-2.682724e-01*2*z)/b}
h8.Db.laplace2<-function(z,b) {-(7.506133e-07*8*z^7-2.578817e-04*6*z^5+1.894136e-02*4*z^3-2.682724e-01*2*z)*z/b}

###############NVE.unif#######################

NVE.unif<-function(sample,order,parest=NA){
  #choose order between 3 and 8#
 if (order<9 && order>2) {
  n<-length(sample)
  if (is.na(parest[1])) {
  tmp<-unif.MME(sample)
  a<-tmp[1]
  b<-tmp[2]
  }
  else {
  a<-parest[1]
  b<-parest[2]
  }
  sigma<-b-a
  z<-(sample-a)/sigma # also location-scale family but scale sigma depends on location a!!
  var_est<-var.unif2(order,z,sigma)
 }
 else {var_est<-NA}
 return(var_est)
}

var.unif2<-function(r,z,sigma) {
  n<-length(z)
  Eh.Da<-switch(r-2,
	mean(h3.Da.unif2(z,sigma)),
	mean(h4.Da.unif2(z,sigma)),
	mean(h5.Da.unif2(z,sigma)),
	mean(h6.Da.unif2(z,sigma)),
	mean(h7.Da.unif2(z,sigma)),
	mean(h8.Da.unif2(z,sigma)))
  Eh.Db<-switch(r-2,
	mean(h3.Db.unif2(z,sigma)),
	mean(h4.Db.unif2(z,sigma)),
	mean(h5.Db.unif2(z,sigma)),
	mean(h6.Db.unif2(z,sigma)),
	mean(h7.Db.unif2(z,sigma)),
	mean(h8.Db.unif2(z,sigma)))
  v<-mean((g_r.unif2(r,z,sigma,Eh.Da,Eh.Db))^2)-(mean(g_r.unif2(r,z,sigma,Eh.Da,Eh.Db)))^2
  v<-v*n/(n-1)
  return(v)
}

g_r.unif2<-function(r,z,sigma,Eh.Da,Eh.Db) {
  h<-orth.poly(r,"unif",c(0,1))
  tmp<-h(z)+b.a.unif2(z,sigma)*Eh.Da+b.b.unif2(z,sigma)*Eh.Db
  return(tmp)
}

b.a.unif2<-function(z,sigma) {sigma*(4*z-1-3*z^2)}
b.b.unif2<-function(z,sigma) {sigma*(-2*z+3*z^2)}

h3.Da.unif2<-function(z,sigma) {-(7^(1/2)*(12-60*z+60*z^2))*(1-z)/sigma}
h3.Db.unif2<-function(z,sigma) {-(7^(1/2)*(12-60*z+60*z^2))*z/sigma}

h4.Da.unif2<-function(z,sigma) {-(-60+540*z-1260*z^2+840*z^3)*(1-z)/sigma}
h4.Db.unif2<-function(z,sigma) {-(-60+540*z-1260*z^2+840*z^3)*z/sigma}

h5.Da.unif2<-function(z,sigma) {-(11^(1/2)*(30-420*z+1680*z^2-2520*z^3+1260*z^4))*(1-z)/sigma}
h5.Db.unif2<-function(z,sigma) {-(11^(1/2)*(30-420*z+1680*z^2-2520*z^3+1260*z^4))*z/sigma}

h6.Da.unif2<-function(z,sigma) {-(13^(1/2)*(-42+840*z-5040*z^2+12600*z^3-13860*z^4+5544*z^5))*(1-z)/sigma}
h6.Db.unif2<-function(z,sigma) {-(13^(1/2)*(-42+840*z-5040*z^2+12600*z^3-13860*z^4+5544*z^5))*z/sigma}

h7.Da.unif2<-function(z,sigma) {-(15^(1/2)*(56-1512*z+12600*z^2-46200*z^3+83160*z^4-72072*z^5+24024*z^6))*(1-z)/sigma}
h7.Db.unif2<-function(z,sigma) {-(15^(1/2)*(56-1512*z+12600*z^2-46200*z^3+83160*z^4-72072*z^5+24024*z^6))*z/sigma}

h8.Da.unif2<-function(z,sigma) {-(17^(1/2)*(-72+2520*z-27720*z^2+138600*z^3-360360*z^4+504504*z^5-360360*z^6+102960*z^7))*(1-z)/sigma}
h8.Db.unif2<-function(z,sigma) {-(17^(1/2)*(-72+2520*z-27720*z^2+138600*z^3-360360*z^4+504504*z^5-360360*z^6+102960*z^7))*z/sigma}

###############NVE.ZIP#######################

#NVE.ZIP<-function(sample,order,parest=NA){
#  #choose order between 3 and 4#
# if (order<5 && order>2) {
#  n<-length(sample)
#  if (is.na(parest[1])) {
#  tmp<-ZIP.MME(sample)
#  lambda<-tmp[1]
#  p<-tmp[2]
#  }
#  else {
#  lambda<-parest[1]
#  p<-parest[2]
#  }
#  var_est<-var.ZIP2(order,sample,lambda,p)
# }
# else {var_est<-NA}
# return(var_est)
#}

#var.ZIP2<-function(r,x,lambda,p) {
#  n<-length(x)
#  Eh.Da<-switch(r-2,
#	mean(h3.Da.ZIP2(x,lambda,p)),
#	mean(h4.Da.ZIP2(x,lambda,p))) #a is short for lambda
#  Eh.Db<-switch(r-2,
#	mean(h3.Db.ZIP2(x,lambda,p)),
#	mean(h4.Db.ZIP2(x,lambda,p))) #b is short for p
#  v<-mean((g_r.ZIP2(r,x,lambda,p,Eh.Da,Eh.Db))^2)-(mean(g_r.ZIP2(r,x,lambda,p,Eh.Da,Eh.Db)))^2
#  v<-v*n/(n-1)
#  return(v)
#}

#g_r.ZIP2<-function(r,x,lambda,p,Eh.Da,Eh.Db) {
#  h<-orth.poly(r,"ZIP",c(lambda,p))
#  tmp<-h(x)+b.a.ZIP2(x,lambda,p)*Eh.Da+b.b.ZIP2(x,lambda,p)*Eh.Db
#  return(tmp)
#}

#b.a.ZIP2<-function(x,lambda,p) {}
#b.b.ZIP2<-function(x,lambda,p) {}

#h3.Da.ZIP2<-function(x,lambda,p) {-1/4*(-24*lambda^5*p^5+48*lambda^5*p^4-48*lambda^5*p^2-216*lambda^3*p+72*lambda^3*p^5*x^3-72*lambda*x^2-72*lambda^2*x+72*x^3+216*lambda^2*p^4*x^3+288*lambda^4*p^5*x+234*lambda^4*x^2*p^4-396*lambda^4*x^2*p^5-118*lambda^4*p^4*x^3-648*lambda^2*x^2*p^4+216*lambda^3*p^2-336*lambda^3*p^3*x^3+360*lambda^3*p^4*x+144*lambda^3*p^5*x-432*lambda^3*x^2*p^4-216*lambda^3*x^2*p^5+72*lambda^3*p^4*x^3-72*lambda^3*p^3+27*lambda^6*p^3*x^3+7*lambda^8*p^5*x+8*lambda^8*x^2*p^4-4*lambda^8*x^2*p^5-236*lambda^4*p^4*x-171*lambda^6*x^2*p^5-54*lambda^6*p^4*x^3+27*lambda^6*p^5*x^3-96*lambda^4*p^4+96*lambda^4*p-288*lambda^4*p^2+288*lambda^4*p^3+10*lambda^7*p^5+4*lambda^8*p^5-30*lambda^7*p^4+4*lambda^7*p^3*x^3-46*lambda^7*p^3*x^2+92*lambda^7*x^2*p^4-46*lambda^7*x^2*p^5-8*lambda^7*p^4*x^3+4*lambda^7*p^5*x^3-12*lambda^8*p^4+12*lambda^8*p^3+30*lambda^7*p^3-10*lambda^7*p^2-4*lambda^8*p^2-216*x^2+24*lambda^5*p+72*lambda*x+72*lambda^3+20*lambda^7*p^2*x-36*lambda^6*p^2*x^2+20*lambda^5*p^2*x^3+38*lambda^7*p^5*x-56*lambda^7*p^4*x-2*lambda^7*p^3*x+306*lambda^6*p^4*x^2-99*lambda^6*p^3*x^2+34*lambda^5*p^3*x^3-504*lambda^2*p^2*x+576*lambda^2*p*x-504*lambda*p*x^2+144*x+288*lambda*p*x+168*lambda^3*p^2*x^3-136*lambda^5*p^3*x+142*lambda^5*p^2*x+24*lambda^5*p*x-462*lambda^4*p^2*x^2-48*lambda^4*p*x^2+24*lambda^3*p*x^3+98*lambda^4*p^2*x^3-72*lambda^6*p^3*x+90*lambda^6*p^2*x-210*lambda^5*p^2*x^2-250*lambda^5*p^4*x+672*lambda^4*p^3*x^2+108*lambda^6*p^5*x+474*lambda^5*p^4*x^2-88*lambda^4*p^3*x^3-126*lambda^6*p^4*x+78*lambda^5*p^3*x^2-332*lambda^4*p^3*x-432*lambda^3*p^2*x^2+220*lambda^4*p^2*x+60*lambda^4*p*x-216*lambda^3*p*x^2+552*lambda^3*p^2*x+220*lambda^5*p^5*x-342*lambda^5*x^2*p^5-128*lambda^5*p^4*x^3+74*lambda^5*p^5*x^3+108*lambda^4*p^5*x^3+1296*lambda^3*p^3*x^2+7*lambda^8*p^3*x+48*lambda^3*p*x+432*lambda^2*p^4*x-324*lambda^2*p^3*x^3+72*p^2*x^3-432*lambda^2*p^3*x+1224*lambda*p^2*x^2-144*p*x^3-432*lambda*p^2*x^3-1104*lambda^3*p^3*x+432*lambda^2*p^2*x^2-540*lambda^2*p*x^2+216*lambda*p*x^3-648*lambda*p^3*x^2+108*lambda^2*p*x^3+216*lambda*p^3*x^3+756*lambda^2*p^3*x^2+432*lambda*p^3*x-216*p^2*x^2-792*lambda*p^2*x+432*p*x^2+144*p^2*x-288*p*x-4*lambda^8*p^3*x^2-14*lambda^8*p^4*x)*lambda^2*2^(1/2)/(-lambda^9*(lambda^2*p+2*lambda*p+2)*(6*lambda*p+3*lambda^2*p+lambda^3*p+6)*(-1+p)^5)^(1/2)/(6*lambda*p+3*lambda^2*p+lambda^3*p+6)/(lambda^2*p+2*lambda*p+2)}
#h3.Db.ZIP2<-function(x,lambda,p) {-1/4*(24*lambda^5*p^4-48*lambda^4+96*lambda^5*p^2+48*lambda^3*p-24*lambda^5-4*lambda^6+72*lambda*x^2-72*lambda^2*x-24*x^3-8*lambda^6*x^2*p-lambda^9*p^3*x-180*lambda^4*x^2*p^4+36*lambda^4*p^4*x^3-24*lambda^3*p^2+48*lambda^3*p^3*x^3+48*lambda^3*p^4*x-72*lambda^3*x^2*p^4+24*lambda^3*p^4*x^3-96*lambda^5*p^3-7*lambda^6*p^3*x^3-40*lambda^6*p-36*lambda^2*x^3*p^2-2*lambda^8*x^2*p^4-4*lambda^8*p+144*lambda^4*p^4*x+7*lambda^6*p^4*x^3+96*lambda^4*p-48*lambda^4*p^2+10*lambda^7*p^4-lambda^7*p^3*x^3+19*lambda^7*p^3*x^2-19*lambda^7*x^2*p^4+lambda^7*p^4*x^3+2*lambda^8*p^4-8*lambda^8*p^3-40*lambda^7*p^3+50*lambda^7*p^2+10*lambda^8*p^2+72*x^2+lambda^9*p^4*x+2*lambda^5*x^3*p+4*lambda^3*x^3+60*lambda^4*x-12*lambda^4*x^2-72*lambda*x+12*lambda^5*x-24*lambda^3-36*lambda^3*x^2-28*lambda^7*p^2*x+26*lambda^6*p^2*x^2-8*lambda^5*p^2*x^3+50*lambda^7*p^4*x-32*lambda^7*p^3*x-75*lambda^6*p^4*x^2+57*lambda^6*p^3*x^2-16*lambda^5*p^3*x^3+144*lambda^2*p^2*x-216*lambda^2*p*x+144*lambda*p*x^2-48*x-72*lambda*p*x-64*lambda^3*p^2*x^3+130*lambda^5*p^3*x-394*lambda^5*p^2*x+64*lambda^5*p*x+318*lambda^4*p^2*x^2-12*lambda^3*p*x^3-34*lambda^4*p^2*x^3-20*lambda^6*p^3*x-156*lambda^6*p^2*x+138*lambda^5*p^2*x^2+188*lambda^5*p^4*x-126*lambda^4*p^3*x^2-162*lambda^5*p^4*x^2-6*lambda^4*p^3*x^3+128*lambda^6*p^4*x+54*lambda^5*p^3*x^2+348*lambda^4*p^3*x+288*lambda^3*p^2*x^2-416*lambda^4*p^2*x-136*lambda^4*p*x+180*lambda^3*p*x^2-30*lambda^5*p*x^2+48*lambda^6*p*x+4*lambda^4*p*x^3+116*lambda^6*p^2-96*lambda^6*p^3+24*lambda^6*p^4+16*lambda^3*p^2*x+80*lambda^3*x+22*lambda^5*p^4*x^3-360*lambda^3*p^3*x^2-11*lambda^8*p^3*x+10*lambda^7*x*p-456*lambda^3*p*x-20*lambda^7*p+72*lambda^2*p^3*x^3+144*lambda^2*p^3*x-216*lambda*p^2*x^2+24*p*x^3+72*lambda*p^2*x^3+312*lambda^3*p^3*x-108*lambda^2*p^2*x^2+324*lambda^2*p*x^2-72*lambda*p*x^3-36*lambda^2*p*x^3-216*lambda^2*p^3*x^2+144*lambda*p^2*x-72*p*x^2+48*p*x+2*lambda^8*p^3*x^2+11*lambda^8*p^4*x)*lambda^3*2^(1/2)/(-lambda^9*(lambda^2*p+2*lambda*p+2)*(6*lambda*p+3*lambda^2*p+lambda^3*p+6)*(-1+p)^5)^(1/2)/(lambda^2*p+2*lambda*p+2)/(6*lambda*p+3*lambda^2*p+lambda^3*p+6)}

#h4.Da.ZIP2<-function(x,lambda,p) {-1/12*(-3456*lambda^4+5184*lambda^5*p^2+20736*lambda*x^2-20736*x^3+864*lambda^6*x^2*p+3456*x^4*lambda^3*p^3+288*x^4*lambda^4*p+10368*lambda^2*p^2*x^4-1248*lambda^9*p^3*x+49*lambda^9*p^3*x^4+270*lambda^8*p^3*x^4+5*lambda^10*p^3*x^4-20736*lambda^3*p^3*x^3-23400*lambda^6*p^3*x^3-1728*lambda^6*p-62208*lambda^2*x^3*p^2-8568*lambda^6*x^3*p^2-3708*lambda^8*x^3*p^3-216*lambda^8*x^3*p^2-1812*lambda^7*x^3*p^2+3456*lambda^4*p-10908*lambda^7*p^3*x^3+28206*lambda^7*p^3*x^2+5184*lambda^2*p*x^4-1944*lambda^8*x*p^2+864*lambda^7*p^3-576*lambda^7*p^2+38016*x^2-864*lambda^5*x^3*p-132*lambda^10*p^3*x^3-888*lambda^9*p^3*x^3-5184*lambda^5*p-13824*lambda*x-9*lambda^11*p^3*x^3-84*lambda^10*p^3+84*lambda^10*p^2+180*lambda^9*p^2-270*lambda^10*p^3*x+5*lambda^12*p^3*x-27*lambda^11*p^3*x-180*lambda^9*p^3-4236*lambda^7*p^2*x+27108*lambda^6*p^2*x^2-26568*lambda^5*p^2*x^3-12204*lambda^7*p^3*x+54180*lambda^6*p^3*x^2-36288*lambda^5*p^3*x^3-62208*lambda^2*p^2*x-72576*lambda^2*p*x+114048*lambda*p*x^2-20736*x-62208*lambda*p*x-576*lambda^9*p^2*x-82944*lambda^3*p^2*x^3-46656*lambda^5*p^3*x-26568*lambda^5*p^2*x+864*lambda^5*p*x+128736*lambda^4*p^2*x^2+16992*lambda^4*p*x^2-20736*lambda^3*p*x^3+1728*lambda^3*p*x^4-57024*lambda^4*p^2*x^3-28008*lambda^6*p^3*x-8712*lambda^6*p^2*x+66420*lambda^5*p^2*x^2+77760*lambda^4*p^3*x^2-38016*lambda^4*p^3*x^3+78624*lambda^5*p^3*x^2-44928*lambda^4*p^3*x+176256*lambda^3*p^2*x^2-77760*lambda^4*p^2*x+5760*lambda^4*p*x+50112*lambda^3*p*x^2+5184*lambda^5*p*x^2-864*lambda^6*p*x-5760*lambda^4*p*x^3+1728*lambda^6*p^3+lambda^13*p^3*x-96*lambda^10*p^2*x-103680*lambda^3*p^2*x-6912*lambda*x^3+6912*lambda^3*x+8274*x^2*lambda^7*p^2-18*lambda^11*p^3+18*lambda^11*p^2+38016*lambda^3*p^3*x^2+3275*lambda^9*p^3*x^2-4284*lambda^8*p^3*x-288*lambda^7*x*p+72*lambda^11*p^3*x^2+661*lambda^10*p^3*x^2+1728*lambda^8*p^2*x^2-31104*lambda^3*p*x+3456*x^4-288*lambda^7*p+4320*x^4*lambda^5*p^3+5184*x^4*lambda^4*p^3+6048*x^4*lambda^4*p^2+2268*x^4*lambda^5*p^2+954*x^4*lambda^7*p^3+78*x^4*lambda^7*p^2+540*x^4*lambda^6*p^2+2412*x^4*lambda^6*p^3+10368*x^4*lambda^3*p^2-20736*lambda^3*p^3*x+114048*lambda^2*p^2*x^2+119232*lambda^2*p*x^2-62208*lambda*p*x^3-51840*lambda^2*p*x^3+216*lambda^9*p^2*x^2+3*lambda^12*p^3*x^2+11178*lambda^8*p^3*x^2+10368*lambda*p*x^4)*6^(1/2)/(-(-1+p)*lambda^4*(lambda^4*p+4*lambda^3*p+12*lambda^2*p+24*lambda*p+24)*(6*lambda*p+3*lambda^2*p+lambda^3*p+6))^(1/2)/(6*lambda*p+3*lambda^2*p+lambda^3*p+6)/(lambda^4*p+4*lambda^3*p+12*lambda^2*p+24*lambda*p+24)/lambda}
#h4.Db.ZIP2<-function(x,lambda,p) {-1/12*6^(1/2)*(-864*lambda^4-1728*lambda^5-864*lambda^6+10368*lambda*x^2-5184*lambda^2*x-5184*x^3-1926*lambda^6*x^2*p+864*x^4*lambda^3*p^3+72*x^4*lambda^4*p+2592*lambda^2*p^2*x^4-2598*lambda^9*p^3*x+10*lambda^9*p^3*x^4+57*lambda^8*p^3*x^4-84*lambda^10*p-12*lambda^11*p-360*lambda^9*p+lambda^10*p^3*x^4-5184*lambda^3*p^3*x^3-7632*lambda^6*p^3*x^3-864*lambda^6*p-15552*lambda^2*x^3*p^2-2916*lambda^6*x^3*p^2-1122*lambda^8*x^3*p^3-72*lambda^8*x^3*p^2-612*lambda^7*x^3*p^2-972*lambda^8*p+864*lambda^4*p-3456*lambda^7*p^3*x^3+14472*lambda^7*p^3*x^2+1296*lambda^2*p*x^4-36*lambda^8-3060*lambda^8*x*p^2-504*lambda^8*p^3-864*lambda^7*p^3+2592*lambda^7*p^2+1512*lambda^8*p^2+9504*x^2+144*lambda^7*x-1296*lambda^5*x^2+216*lambda^5*x^3*p+1080*lambda^6*x+504*lambda^4*x^3+3384*lambda^4*x-39*lambda^10*p^3*x^3-258*lambda^9*p^3*x^3+1230*lambda^8*x*p+1728*lambda^5*p+306*lambda^9*p*x-1218*x^2*lambda^7*p-2124*lambda^4*x^2-6912*lambda*x+2880*lambda^5*x+42*lambda^10*p*x-3*lambda^11*p^3*x^3+144*lambda^5*x^3-18*x^4*lambda^6*p-6*x^4*lambda^7*p-36*x^4*lambda^5*p-42*lambda^10*p^3+126*lambda^10*p^2+540*lambda^9*p^2-360*lambda^8*p*x^2-693*lambda^10*p^3*x-16*lambda^12*p^3*x-132*lambda^11*p^3*x-180*lambda^9*p^3-10404*lambda^7*p^2*x+16938*lambda^6*p^2*x^2-9288*lambda^5*p^2*x^3-14688*lambda^7*p^3*x+26208*lambda^6*p^3*x^2-11664*lambda^5*p^3*x^3-15552*lambda^2*p^2*x-28512*lambda^2*p*x+28512*lambda*p*x^2-5184*x-15552*lambda*p*x-612*lambda^9*p^2*x-25920*lambda^3*p^2*x^3-22032*lambda^5*p^3*x-41256*lambda^5*p^2*x-6264*lambda^5*p*x+63288*lambda^4*p^2*x^2+14616*lambda^4*p*x^2-7776*lambda^3*p*x^3+432*lambda^3*p*x^4-19440*lambda^4*p^2*x^3-21744*lambda^6*p^3*x-24516*lambda^6*p^2*x-216*lambda^6*x^2+39636*lambda^5*p^2*x^2-36*x^4*lambda^4+24624*lambda^4*p^3*x^2-11232*lambda^4*p^3*x^3+32616*lambda^5*p^3*x^2-14688*lambda^4*p^3*x+59616*lambda^3*p^2*x^2-45360*lambda^4*p^2*x-23616*lambda^4*p*x+35856*lambda^3*p*x^2+1764*lambda^5*p*x^2+1980*lambda^6*p*x-2304*lambda^4*p*x^3+2592*lambda^6*p^2-864*lambda^6*p^3-lambda^13*p^3*x-72*lambda^10*p^2*x-288*lambda^7+156*lambda^7*x^3*p+30*lambda^8*x^3*p+396*lambda^6*x^3*p-36288*lambda^3*p^2*x-3456*lambda*x^3+5184*lambda^2*x^2-3456*lambda^3*x-54*lambda^9*p*x^2+4950*x^2*lambda^7*p^2-6*lambda^11*p^3+18*lambda^11*p^2+9504*lambda^3*p^3*x^2+1658*lambda^9*p^3*x^2-7194*lambda^8*p^3*x+2796*lambda^7*x*p+45*lambda^11*p^3*x^2+341*lambda^10*p^3*x^2+972*lambda^8*p^2*x^2-28512*lambda^3*p*x+864*x^4-1440*lambda^7*p+1080*x^4*lambda^5*p^3+1296*x^4*lambda^4*p^3+1512*x^4*lambda^4*p^2+540*x^4*lambda^5*p^2+216*x^4*lambda^7*p^3+18*x^4*lambda^7*p^2+126*x^4*lambda^6*p^2+576*x^4*lambda^6*p^3+2592*x^4*lambda^3*p^2-5184*lambda^3*p^3*x+28512*lambda^2*p^2*x^2+45360*lambda^2*p*x^2-15552*lambda*p*x^3-18144*lambda^2*p*x^3+108*lambda^9*p^2*x^2+3*lambda^12*p^3*x^2+5739*lambda^8*p^3*x^2+2592*lambda*p*x^4)/(-1+p)/(lambda^4*p+4*lambda^3*p+12*lambda^2*p+24*lambda*p+24)/(6*lambda*p+3*lambda^2*p+lambda^3*p+6)/(-(-1+p)*lambda^4*(lambda^4*p+4*lambda^3*p+12*lambda^2*p+24*lambda*p+24)*(6*lambda*p+3*lambda^2*p+lambda^3*p+6))^(1/2)}

