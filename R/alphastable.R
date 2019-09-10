urstab<-function(n,alpha,beta,sigma,mu,param)
{
	stopifnot(0<alpha,alpha<=2,length(alpha)==1,-1<=beta,beta<=1,length(beta)==1,0<=sigma,length(param)==1, param %in% 0:1)
	theta<-runif(n,-pi/2,pi/2)
	theta0<-atan(beta*tan(pi*alpha/2))/alpha
	x<-c()
	w<-rexp(n,1)
		if (param==0)
		{
			if (alpha==1)
			{
				x<-sigma*2/pi*((pi/2+beta*theta)*tan(theta)-beta*log((pi/2*w*cos(theta))/(pi/2+beta*theta)))+2/pi*beta*sigma*log(sigma)+mu
			}
			else
			{
			x<-sigma*sin(alpha*(theta0+theta))/(cos(alpha*theta0)*cos(theta))^(1/alpha)*(cos(alpha*theta0+(alpha-1)*theta)/w)^((1-alpha)/alpha)+mu
			}
		}
		else
		{
			if(alpha!= 1)
			{
				x<-x-beta*sigma*tan(pi*alpha/2)
			}
			else
			{
				x<-x-2/pi*beta*sigma*log(sigma)
			}
		}
	return(x)
}
urstab.trunc<-function(n,alpha,beta,sigma,mu,a,b,param)
{
	stopifnot(0<alpha,alpha<=2,length(alpha)==1,-1<beta,beta<1,length(beta)==1,0<sigma,length(sigma)==1,a<b)
	y<-c()
		if (alpha==1)
		{
			if (param==0)
			{
				y<-c()
                        for (i in 1:n)
				{
				w<-rexp(1)
				l<-uniroot(function(theta){2/pi*((pi/2+beta*theta)*tan(theta)-beta*log((pi/2*cos(theta))/(pi/2+beta*theta)))-(a-mu)/sigma-2/pi*beta*log(w)},lower=-pi/2,upper=pi/2)$root
				u<-uniroot(function(theta){2/pi*((pi/2+beta*theta)*tan(theta)-beta*log((pi/2*cos(theta))/(pi/2+beta*theta)))-(b-mu)/sigma-2/pi*beta*log(w)},lower=-pi/2,upper=pi/2)$root
                        tc<-runif(1,l,u)
				y[i]<-sigma*2/pi*((pi/2+beta*tc)*tan(tc)-beta*log((pi/2*w*cos(tc))/(pi/2+beta*tc)))+mu
				}
			}
			else
			{
				for (i in 1:n)
				{
				w<-rexp(1)
				l<-uniroot(function(theta){2/pi*((pi/2+beta*theta)*tan(theta)-beta*log((pi/2*cos(theta))/(pi/2+beta*theta)))-(a-mu-beta*2/pi*sigma*log(sigma))/sigma-2/pi*beta*log(w)},lower=-pi/2,upper=pi/2)$root
				u<-uniroot(function(theta){2/pi*((pi/2+beta*theta)*tan(theta)-beta*log((pi/2*cos(theta))/(pi/2+beta*theta)))-(b-mu-beta*2/pi*sigma*log(sigma))/sigma-2/pi*beta*log(w)},lower=-pi/2,upper=pi/2)$root
                        tc<-runif(1,l,u)
				y[i]<-sigma*2/pi*((pi/2+beta*tc)*tan(tc)-beta*log((pi/2*w*cos(tc))/(pi/2+beta*tc)))+mu+beta*2/pi*sigma*log(sigma)
				}
			}
		}
		else
		{
				theta0<-atan(beta*tan(pi*alpha/2))/alpha
			if (param==1)
			{
					for(i in 1:n)
					{
						w<-rexp(1)
						l<-uniroot(function(theta){sin(alpha*(theta0+theta))/(cos(alpha*theta0)*cos(theta))^(1/alpha)*(cos(alpha*theta0+(alpha-1)*theta))^((1-alpha)/alpha)-(a-mu)/(sigma*w^((alpha-1)/alpha))},lower=-pi/2,upper=pi/2)$root
						u<-uniroot(function(theta){sin(alpha*(theta0+theta))/(cos(alpha*theta0)*cos(theta))^(1/alpha)*(cos(alpha*theta0+(alpha-1)*theta))^((1-alpha)/alpha)-(b-mu)/(sigma*w^((alpha-1)/alpha))},lower=-pi/2,upper=pi/2)$root
						uu<-runif(1,l,u)
						y[i]<-sigma*sin(alpha*(theta0+uu))/(cos(alpha*theta0)*cos(uu))^(1/alpha)*(cos(alpha*theta0+(alpha-1)*uu)/w)^((1-alpha)/alpha)+mu
					}
			}
			else
			{
					for(i in 1:n)
					{
						w<-rexp(1)
						l<-uniroot(function(theta){sin(alpha*(theta0+theta))/(cos(alpha*theta0)*cos(theta))^(1/alpha)*(cos(alpha*theta0+(alpha-1)*theta))^((1-alpha)/alpha)-(a-mu+sigma*beta*tan(pi*alpha/2))/(sigma*w^((alpha-1)/alpha))},lower=-pi/2,upper=pi/2)$root
						u<-uniroot(function(theta){sin(alpha*(theta0+theta))/(cos(alpha*theta0)*cos(theta))^(1/alpha)*(cos(alpha*theta0+(alpha-1)*theta))^((1-alpha)/alpha)-(b-mu+sigma*beta*tan(pi*alpha/2))/(sigma*w^((alpha-1)/alpha))},lower=-pi/2,upper=pi/2)$root
						uu<-runif(1,l,u)
						y[i]<-sigma*sin(alpha*(theta0+uu))/(cos(alpha*theta0)*cos(uu))^(1/alpha)*(cos(alpha*theta0+(alpha-1)*uu)/w)^((1-alpha)/alpha)+mu-sigma*beta*tan(pi*alpha/2)
					}
			}
		}
	return(y)
}
mrstab.elliptical<-function(n,alpha,Sigma,Mu)
{
	stopifnot(0<alpha,alpha<=2,length(alpha)==1,dim(Sigma)[1]==dim(Sigma)[2],length(Mu)==dim(Sigma)[1])
	d<-dim(Sigma)[1]
	x<-matrix(0, nrow=n, ncol=d)
		for(i in 1:n)
		{
			x[i,]<-suppressWarnings(Mu+sqrt(rstable(1,alpha/2,1,cos(pi*alpha/4)^(2/alpha),0,1))*rmvnorm(1,c(rep(0,d)),Sigma))
		}
	return(x)
}
mrstab<-function(n,m,alpha,Gamma,Mu)
{
	stopifnot(0<alpha,alpha<=2,length(alpha)==1,length(Gamma)==m,length(Mu)==2)
	x<-matrix(0,nrow=n,ncol=2)
	S<-L<-matrix(2*m,nrow=2,ncol=m)
		for (j in 1:m)
			{
				S[1,j]<-cos(2*(j-1)*pi/m)
				S[2,j]<-sin(2*(j-1)*pi/m)
			}
					for (i in 1:n)
					{
							for (j in 1:m)
							{
							L[,j]<-(Gamma[j])^(1/alpha)*(rstable(1,alpha,1,1,0,1)*S[,j])
							}
						x[i,]<-apply(L,1,sum)+Mu
					}
	return(x)
}
udstab<-function(x,alpha,beta,sigma,mu,param)
{
	stopifnot(0<alpha,alpha<=2,length(alpha)==1,-1<=beta,beta<=1,length(beta)==1,0<=sigma,length(param)==1, param %in% 0:1)
	k<-seq(1,150)
	xi<--beta*sigma*tan(pi*alpha/2)
	eta<--beta*tan(pi*alpha/2)
	r<-(1+eta^2)^(1/(2*alpha))
	i<-150
	if (x==mu)
	{
	pdf<-suppressWarnings(dstable(x,alpha,beta,sigma,mu,param))
	}
	else
	{
		if(alpha==1)
		{
			pdf<-suppressWarnings(dstable(x,1,beta,sigma,mu,param))
		}
		else
		{
			if (param==0)
			{
				L1<--sigma*r*(alpha*exp(lgamma(alpha*i+alpha)-lgamma(alpha*i+1)))^(1/alpha)+mu+xi-alpha/2
				L2<-sigma*r*(alpha*exp(lgamma(alpha*i+alpha)-lgamma(alpha*i+1)))^(1/alpha)+mu+xi+alpha/2
				U1<--sigma*r*alpha*exp(lgamma(i/alpha+1)-lgamma(i/alpha+1/alpha))+mu+xi+2*alpha
				U2<-sigma*r*alpha*exp(lgamma(i/alpha+1)-lgamma(i/alpha+1/alpha))+mu+xi-2*alpha
					if(x<L1 || x>L2)
					{
						pdf<-1/(pi*abs(x-mu-xi))*sum((-1)^(k-1)*exp(lgamma(alpha*k+1)-lgamma(k+1))*(abs(x-mu-xi)/(sigma*r))^(-alpha*k)*sin(k*pi/2*(alpha+2/pi*atan(beta*tan(pi*alpha/2))*sign(x-mu-xi))))
					}
					else
					{
						if(x<U2 && x>U1)
						{
							pdf<-1/(pi*abs(x-mu-xi))*sum((-1)^(k-1)*exp(lgamma(k/alpha+1)-lgamma(k+1))*(abs(x-mu-xi)/(sigma*r))^k*sin(k*pi*(alpha+2/pi*atan(beta*tan(pi*alpha/2))*sign(x-mu-xi))/(2*alpha)))
						}
						else
						{
							pdf<-suppressWarnings(dstable(x,alpha,beta,sigma,mu,0))
						}
					}
			}
			else
			{
				L1<--sigma*r*(alpha*exp(lgamma(alpha*i+alpha)-lgamma(alpha*i+1)))^(1/alpha)+mu-alpha
				L2<-sigma*r*(alpha*exp(lgamma(alpha*i+alpha)-lgamma(alpha*i+1)))^(1/alpha)+mu+alpha
				U1<--sigma*r*alpha*exp(lgamma(i/alpha+1)-lgamma(i/alpha+1/alpha))+mu+2*alpha
				U2<-sigma*r*alpha*exp(lgamma(i/alpha+1)-lgamma(i/alpha+1/alpha))+mu-2*alpha
					if(x<L1 || x>L2)
					{
						pdf<-1/(pi*abs(x-mu))*sum((-1)^(k-1)*exp(lgamma(alpha*k+1)-lgamma(k+1))*(abs(x-mu)/(sigma*r))^(-alpha*k)*sin(k*pi/2*(alpha+2/pi*atan(beta*tan(pi*alpha/2))*sign(x-mu))))
					}
					else
					{
						if(x<U2 && x>U1)
						{
							pdf<-1/(pi*abs(x-mu))*sum((-1)^(k-1)*exp(lgamma(k/alpha+1)-lgamma(k+1))*(abs(x-mu)/(sigma*r))^k*sin(k*pi*(alpha+2/pi*atan(beta*tan(pi*alpha/2))*sign(x-mu))/(2*alpha)))
						}
						else
						{
							pdf<-suppressWarnings(dstable(x,alpha,beta,sigma,mu,1))
						}
					}
			}
		}
	}
	return(pdf)
}
mdstab.elliptical<-function(x,alpha,Sigma,Mu)
{
	stopifnot(0<alpha,alpha<=2,length(alpha)==1,length(Mu)==length(x),dim(Sigma)[1]==dim(Sigma)[2],length(Mu)==dim(Sigma)[1])
	if(dim(Sigma)[1]!=dim(Sigma)[2]){message("matrix Sigma must be square")}
	if(length(Mu)!=dim(Sigma)[1]){message("matrix Sigma and Mu must be of the same dimensions")}
	if(length(Mu)!=length(x)){message("vector x and Mu must be of the same dimensions")}
	if(min(eigen(Sigma)$values)<0){message("matrix Sigma is not positive definite")}
	d<-length(Mu)
	dd<-(x-Mu)%*%solve(Sigma)%*%cbind(x-Mu)/2
	k<-150
	j<-seq(1,150)
	r<-2*(exp(lgamma(alpha*k/2+alpha/2+1)+lgamma(alpha*k/2+alpha/2+d/2)-lgamma(alpha*k/2+1)-lgamma(alpha*k/2+d/2))/(k+1))^(2/alpha)
		if (dd>r)
		{
			pdf<-.5*suppressWarnings(1/((2*pi)^(d/2)*pi*sqrt(det(Sigma)))*sum(2^(alpha*j/2+d/2)*(-1)^(j-1)*dd^(-alpha*j/2-d/2)*exp(lgamma(alpha*j/2+1)+lgamma(alpha*j/2+d/2)-lgamma(j+1))*sin(j*pi*alpha*.5)))
		}
		else
		{
			rr<-rstable(5000,alpha/2,1,(cos(pi*alpha/4))^(2/alpha),0,1)
			pdf<-suppressWarnings(mean(exp(-dd/(4*rr))/((4*pi*rr)^(d/2)*sqrt(det(Sigma)))))
		}
	return(pdf)
}
upstab<-function(x,alpha,beta,sigma,mu,param)
{
	stopifnot(0<alpha,alpha<=2,length(alpha)==1,-1<=beta,beta<=1,length(beta)==1,0<=sigma,length(param)==1, param %in% 0:1)
	k<-seq(1,150)
	xi<--beta*sigma*tan(pi*alpha/2)
	eta<--beta*tan(pi*alpha/2)
	r<-(1+eta^2)^(1/(2*alpha))
	i<-150
	if (x==mu)
	{
	cdf<-suppressWarnings(pstable(x,alpha,beta,sigma,mu,param))
	}
	else
	{
		if(alpha==1)
		{
			cdf<-suppressWarnings(pstable(x,1,beta,sigma,mu,param))
		}
		else
		{
			if (param==0)
			{
				L1<--sigma*r*(alpha*exp(lgamma(alpha*i+alpha)-lgamma(alpha*i+1)))^(1/alpha)+mu+xi-alpha/2
				L2<-sigma*r*(alpha*exp(lgamma(alpha*i+alpha)-lgamma(alpha*i+1)))^(1/alpha)+mu+xi+alpha/2
				U1<--sigma*r*alpha*exp(lgamma(i/alpha+1)-lgamma(i/alpha+1/alpha))+mu+xi+2*alpha
				U2<-sigma*r*alpha*exp(lgamma(i/alpha+1)-lgamma(i/alpha+1/alpha))+mu+xi-2*alpha
					if(x<L1 || x>L2)
					{
							cdf<-(1+sign(x-mu-xi))/2+sign(x-mu-xi)/(pi)*sum((-1)^(k)*exp(lgamma(alpha*k+1)-lgamma(k+1))*(abs(x-mu-xi)/(sigma*r))^(-alpha*k)/(alpha*k)*sin(k*pi/2*(alpha+2/pi*atan(beta*tan(pi*alpha/2))*sign(x-mu-xi))))
					}
					else
					{
						if(x<U2 && x>U1)
						{
                        			cdf<-(1/2-atan(beta*tan(pi*alpha/2))/(alpha*pi))-sign(x-mu-xi)/pi*sum((-1)^(k)*exp(lgamma(k/alpha+1)-lgamma(k+1))*(abs(x-mu-xi)/(sigma*r))^(k)/(k)*sin(k*pi*(alpha+2/pi*atan(beta*tan(pi*alpha/2))*sign(x-mu-xi))/(2*alpha)))
						}
						else
						{
							cdf<-suppressWarnings(pstable(x,alpha,beta,sigma,mu,0))
						}
					}
			}
			else
			{
				L1<--sigma*r*(alpha*exp(lgamma(alpha*i+alpha)-lgamma(alpha*i+1)))^(1/alpha)+mu-alpha
				L2<-sigma*r*(alpha*exp(lgamma(alpha*i+alpha)-lgamma(alpha*i+1)))^(1/alpha)+mu+alpha
				U1<--sigma*r*alpha*exp(lgamma(i/alpha+1)-lgamma(i/alpha+1/alpha))+mu+2*alpha
				U2<-sigma*r*alpha*exp(lgamma(i/alpha+1)-lgamma(i/alpha+1/alpha))+mu-2*alpha
					if(x<L1 || x>L2)
					{
							cdf<-(1+sign(x-mu))/2+sign(x-mu)/(pi)*sum((-1)^(k)*exp(lgamma(alpha*k+1)-lgamma(k+1))*(abs(x-mu)/(sigma*r))^(-alpha*k)/(alpha*k)*sin(k*pi/2*(alpha+2/pi*atan(beta*tan(pi*alpha/2))*sign(x-mu))))
					}
					else
					{
						if(x<U2 && x>U1)
						{
                      				cdf<-(1/2-atan(beta*tan(pi*alpha/2))/(alpha*pi))-sign(x-mu)/pi*sum((-1)^(k)*exp(lgamma(k/alpha+1)-lgamma(k+1))*(abs(x-mu)/(sigma*r))^(k)/(k)*sin(k*pi*(alpha+2/pi*atan(beta*tan(pi*alpha/2))*sign(x-mu))/(2*alpha)))
						}
						else
						{
							cdf<-suppressWarnings(pstable(x,alpha,beta,sigma,mu,1))
						}
			}
		}
	}
}
	return(cdf)
}
mfitstab.ustat<-function(u,m,method=method)
{
	stopifnot(length(u[1,])>1,length(u[,1])>2,m %in% 2:40, method %in% 1:2)
	S<-matrix(2*m,nrow=2,ncol=m)
      	for (j in 1:m)
    		{
    			S[1,j]<-cos(2*(j-1)*pi/m)
    			S[2,j]<-sin(2*(j-1)*pi/m)
    		}
	T<-t(S)
	u<-t(u)
	n<-length(u[1,])
	mass<-W<-V<-c()
	L<-matrix(m*m,nrow=m,ncol=m)
				if(method=="1")
					{
						s<-0
							for (i in 1:(n-1))
      		    				{
                      					for (j in (i+1):n)
								{
      								s1<-sqrt(sum(u[,i]^2))
      								s2<-sqrt(sum(u[,j]^2))
      								s3<-sqrt(sum((u[,i]+u[,j])^2))
      								s<-s+(log(s3)-1/2*(log(s1)+log(s2)))/log(2)
      							}
      						}
   	   					alpha.hat<-n*(n-1)/(2*s)
         									for (i in 1:m)
         									{
	   										for (j in 1:m)
	   										{
	   										q<-T[i,1]*S[1,j]+T[i,2]*S[2,j]
	   										L[i,j]<--(abs(q))^(alpha.hat)*(1-sign(q)*tan(pi*alpha.hat/2))
	   										}
	   									}
   												for (i in 1:m)
           											{
	      											hh<-T[i,1]*u[1,]+T[i,2]*u[2,]
	      											V[i]<-complex(real=mean(cos(hh)),imaginary=mean(sin(hh)))
	      											W[i]<-Re(log(V[i]))+Im(log(V[i]))
	      										}
		    											for (i in 1:m)
		    											{
		    												mass[i]<-((nnls(-L,-W)[1])$x[i])
		    											}
					}
					else
					{
						d1<-0;
						d2<-0;
							for (j2 in 1:n)
            					{
								d1[j2]<-u[1,j2]
								d2[j2]<-u[2,j2]
							}
						y<-0
						j1=1:(n/2)
						y[j1]<-log(((u[1,2*j1]+u[1,2*j1-1])^2+(u[2,2*j1]+u[2,2*j1-1])^2)^.5)
						xx<-0
						j1=1:n
						xx[j1]<-log((u[1,j1]^2+u[2,j1]^2)^.5)
   						alpha.hat<-log(2)/(mean(y)-mean(xx))
         							for (i in 1:m)
         							{
	   								for (j in 1:m)
	   								{
	   									q<-T[i,1]*S[1,j]+T[i,2]*S[2,j]
	   									L[i,j]<--(abs(q))^(alpha.hat)*(1-sign(q)*tan(pi*alpha.hat/2))
	   								}
	   							}
   											for (i in 1:m)
            									{
	      										hh<-T[i,1]*u[1,]+T[i,2]*u[2,]
	      										V[i]<-complex(real=mean(cos(hh)),imaginary=mean(sin(hh)))
	      										W[i]<-Re(log(V[i]))+Im(log(V[i]))
	      									}
		    											for (i in 1:m)
		    											{
		    												mass[i]<-((nnls(-L,-W)[1])$x[i])
		    											}
					}
	return(list(alpha=alpha.hat,mass=mass))
}
mfitstab.elliptical<-function(yy,alpha0,Sigma0,Mu0){
s0<-Sigma0
a0<-alpha0
m0<-Mu0
b<-1;
N<-2000;
n<-length(yy[,1]);
d<-length(yy[1,]);
mm<-150;
nn<-120;
jj<-1;
Sigma.matrix<-array(0,dim=c(d,ncol=d,mm))
mu.matrix<-matrix(0,ncol=d,nrow=mm)
alpha.matrix<-c()
Sigma.matrix[,,1]<-s0
alpha.matrix[1]<-a0
mu.matrix[1,]<-m0
sm<-c()
ww<-matrix(0,ncol=d,nrow=n)
	for (r in 2:mm)
	{
		ee<-sqrt(rexp(n,1))
		k<-round(min(160,160/a0))
		vv<-0.5+(exp(lgamma(a0*k/2+a0/2+1)+lgamma(a0*k/2+a0/2+d/2+1)-lgamma(a0*k/2+1)-lgamma(a0*k/2+d/2+1))/((k+1)))^(2/a0)
			for (i in 1:n)
			{
			dis<-(yy[i,]-m0)%*%solve(s0)%*%cbind(yy[i,]-m0)
				if (dis>vv)
				{
					s1<-s2<-0;
						for (j in 1:k)
						{
							s1<-s1+(-1)^(j)*dis^(-a0*j/2-d/2-1)*exp(lgamma(a0*j/2+1)+lgamma(a0*j/2+d/2+1)-lgamma(j+1))*sin(j*pi*a0*.75)
							s2<-s2+(-1)^(j)*dis^(-a0*j/2-d/2)*exp(lgamma(a0*j/2+1)+lgamma(a0*j/2+d/2)-lgamma(j+1))*sin(j*pi*a0*.75)
						}
					sm[i]<-s1/s2
				}
				else
				{
					rr<-(rstable(N,a0/2,1,(cos(pi*a0/4))^(2/a0),0,1))
					ss1<-sum(rr^(-d/2-1)*exp(-.5*dis/rr),na.rm=TRUE)
					ss2<-sum(rr^(-d/2)*exp(-.5*dis/rr),na.rm=TRUE)
					sm[i]<-ss1/ss2
				}
			}
		ee<-sqrt(rexp(n,1))
					for (j in 1:d)
					{
						mu.matrix[r,j]<-suppressWarnings(sum(yy[,j]*sm,na.rm=TRUE)/sum(sm,na.rm=TRUE))
						ww[,j]<-yy[,j]-m0[j]
					}
						m0<-mu.matrix[r,]
						y<-ww/ee
						Z<-c()
							for (i in 1:n)
							{
								num<-y[i,]%*%solve(s0)%*%y[i,]
								up<-d^(d/2)*exp(-d/2)/num^(d/2)
								j<-1
									while (j<2)
									{
										w<-rweibull(1,a0,1)
										ex<-exp(-.5*num*w^2)
											if (runif(1)<w^d*ex/up)
											{
												Z[i]<-w
												j<-j+1
											}
									}
							}
												f<-function(p)
												{
													sum(-log(p[1])-(p[1]-1)*log(Z)+Z^p[1])
												}
													a0<-suppressWarnings(nlm(f,p<-c(a0),hessian=TRUE)$estimate)
													sum<-0
													for (i in 1:n)
													{
														sum<-sum+cbind(y[i,])%*%rbind(y[i,])*Z[i]^2
													}
	s0<-sum/n
	Sigma.matrix[,,r]<-s0
														if (a0>2)
														{
															a0<-1.99
														}
	alpha.matrix[r]<-a0
	}
		a1<-matrix(0,nrow=(mm-nn+1),ncol=1)
		a1<-alpha.matrix[(mm-nn):mm]
		s1<-matrix(0,nrow=d,ncol=d)
				for (i in 1:d)
				{
					for (j in 1:d)
					{
						s1[i,j]<-mean(Sigma.matrix[i,j,(mm-nn):(mm)])
					}
				}
		Sigma<-s1
		alpha=mean(a1)
		mu=apply(mu.matrix[(mm-nn):mm,],2,mean)
		suppressWarnings(return(list(alpha=alpha,Sigma=Sigma,Mu=mu)))
}
ufitstab.sym<-function(yy,alpha0,sigma0,mu0)
{
	n<-length(yy)
 	m<-120
 	N<-2000
 	alphahat<-c()
 	sigmahat<-c()
 	muhat<-c()
 	m0<-mu0
 	a<-matrix(m*3,nrow=m,ncol=3)
 	a[1,1]<-alpha0
 	a[1,2]<-sigma0
 	a[1,3]<-m0
 	a0<-alpha0
 	s0<-sigma0
 	ss<-sm<-Z1<-Z<-c()
		for (r in 1:m-1)
		{
 			k<-round(min(165,165/a0))
 			ee<-sqrt(rexp(n,1))
 			vv<-1+2*s0*(exp(lgamma(a0*k/2+a0/2+1)+lgamma(a0*k/2+a0/2+1/2)-lgamma(a0*k/2+1)-lgamma(a0*k/2+1/2))/((k+1)))^(1/a0)
 				for (i in 1:n)
				{
 					d<-abs(yy[i]-m0)
 						if (d>vv)
						{
 							s1<-s2<-0
 								for (j in 1:k)
								{
 									s1<-s1+(2*s0)^(a0*j+2)/(abs(d)^(a0*j+3))*(-1)^(j-1)*exp(lgamma(a0*j/2+1)+lgamma(a0*j/2+3/2)-lgamma(j+1))*sin(j*pi*a0*.75)/(pi^1.5)
									s2<-s2+(2*s0)^(a0*j)/(abs(d)^(a0*j+1))*(-1)^(j-1)*exp(lgamma(a0*j/2+1)+lgamma(a0*j/2+1/2)-lgamma(j+1))*sin(j*pi*a0*.75)/(pi^1.5)
								}
 							sm[i]<-s1/s2
 						}
						else
						{
 							rr<-rstable(N,a0/2,1,(cos(pi*a0/4))^(2/a0),0,1)
 							sm[i]<-sum(1/(rr^(1.5))*exp(-(d^2/(2*sqrt(rr)*s0)^2)),na.rm=TRUE)/sum(1/(rr^(0.5))*exp(-(d^2/(2*sqrt(rr)*s0)^2)),na.rm=TRUE)
						}
				}
 			m0<-sum(yy*sm,na.rm=TRUE)/sum(sm,na.rm=TRUE)
 			y<-(yy-m0)/ee
 				for (i in 1:n)
				{
 					y0<-y[i]
 					j<-1
 						while (j<2)
						{
 							tt<-rweibull(1,a0,1)
 							ra<-exp(-.5)/(sqrt(2*pi)*abs(y0))
 							u<-runif(1)
 								if (u<dnorm(y0,0,sqrt(2)*s0/tt)/ra)
 								{
									Z1[j]<-tt
 									j<-j+1
								}
						}
 					Z[i]<-Z1
				}
 			f<-function(p){sum(-log(p[1])-(p[1]-1)*log(Z)+Z^p[1])}
 			out<-suppressWarnings(nlm(f, p<-c(a0), hessian=FALSE))
 			a0<-out$estimate[]
 			s0<-sqrt(sum(y^2*Z^2)/(2*n))
 			a[r+1,3]<-m0
 			a[r+1,2]<-s0
 			if (a0>2){a0<-1.95}
 			a[r+1,1]<-a0
		}
	return((list(alpha=mean(a[(m-30):m,1]),sigma=mean(a[(m-30):m,2]),mu=mean(a[(m-30):m,3]))))
}
ufitstab.sym.mix<-function(yy,k,omega0,alpha0,sigma0,mu0)
{
	n<-length(yy)
	m<-150
	N<-4000
	mu.matrix<-matrix(m*k,ncol=k,nrow=m)
	sigma.matrix<-matrix(m*k,ncol=k,nrow=m)
	alpha.matrix<-matrix(m*k,ncol=k,nrow=m)
	p.matrix<-matrix(m*k,ncol=k,nrow=m)
	tau.matrix<-matrix(n*k,ncol=k,nrow=n)
	d<-matrix(n*k,ncol=k,nrow=n)
	sm<-matrix(n*k,ncol=k,nrow=n)
	ss<-matrix(n*k,ncol=k,nrow=n)
	clustering<-rep(0,length(yy))
	vv<-c()
	mu.matrix[1,]<-mu0
	p.matrix[1,]<-omega0
	alpha.matrix[1,]<-alpha0
	sigma.matrix[1,]<-sigma0
	p0<-p.matrix[1,]
	a0<-alpha.matrix[1,]
	s0<-sigma.matrix[1,]
	m0<-mu.matrix[1,]
	a11<-matrix(0,ncol=5,nrow=k)
	a12<-matrix(0,ncol=5,nrow=k)
		for (r in 2:m)
		{
			for (j in 1:n)
			{
				for (ii in 1:k)
				{
					kk<-round(min(168,168/a0[ii]))
					vv[ii]<-2+2*s0[ii]*(exp(lgamma(a0[ii]*kk/2+a0[ii]/2+1)+lgamma(a0[ii]*kk/2+a0[ii]/2+1/2)-lgamma(a0[ii]*kk/2+1)-lgamma(a0[ii]*kk/2+1/2))/((kk+1)))^(1/a0[ii])
					d<-abs(yy[j]-m0[ii])
						if (d>vv[ii])
						{
							s.1<-s.2<-0
								for (jj in 1:kk)
								{
									s.1<-s.1+(2*s0[ii])^(a0[ii]*jj+2)/(abs(d)^(a0[ii]*jj+3))*(-1)^(jj-1)*exp(lgamma(a0[ii]*jj/2+1)+lgamma(a0[ii]*jj/2+3/2)-lgamma(jj+1))*sin(jj*pi*a0[ii]*.75)/(pi^1.5)
									s.2<-s.2+(2*s0[ii])^(a0[ii]*jj+0)/(abs(d)^(a0[ii]*jj+1))*(-1)^(jj-1)*exp(lgamma(a0[ii]*jj/2+1)+lgamma(a0[ii]*jj/2+1/2)-lgamma(jj+1))*sin(jj*pi*a0[ii]*.75)/(pi^1.5)
								}
							sm[j,ii]<-s.1/s.2
						}
						else
						{
							rr<-rstable(N,a0[ii]/2,1,(cos(pi*a0[ii]/4))^(2/a0[ii]),0,1)
							ss1<-sum(1/rr^(1.5)*exp(-(yy[j]-m0[ii])^2/(2*sqrt(rr)*s0[ii])^2),na.rm=TRUE)
							ss2<-sum(1/rr^(0.5)*exp(-(yy[j]-m0[ii])^2/(2*sqrt(rr)*s0[ii])^2),na.rm=TRUE)
							sm[j,ii]<-ss1/ss2
						}
					s.pdf<-0
									for (mm in 1:k)
									{
										s.pdf<-s.pdf+p0[mm]*dstable(yy[j],a0[mm],0,s0[mm],m0[mm],1)
									}
					tau.matrix[j,ii]<-p0[ii]*dstable(yy[j],a0[ii],0,s0[ii],m0[ii],1)/s.pdf
				}
			}
										for (ii in 1:k)
										{
											mu.matrix[r,ii]<-sum(yy*sm[,ii]*tau.matrix[,ii],na.rm=TRUE)/sum(sm[,ii]*tau.matrix[,ii],na.rm=TRUE)
											m0[ii]<-mu.matrix[r,ii]
											p0[ii]<-sum(tau.matrix[,ii])/n
											p.matrix[r,ii]<-p0[ii]
										}
					z<-matrix(0,ncol=k,nrow=n)
											for (j in 1:n)
											{
												max<-tau.matrix[j,1]
												tt<-1
													for (ii in 2:k)
													{
														if (tau.matrix[j,ii]> max)
														{
															max<-tau.matrix[j,ii]
															tt<-ii
														}
													}
												z[j,tt]<-1
											}
															for (bb in 1:k)
															{
																for (rrr in 1:5)
																{
																	n00<-length(yy[z[,bb]==1])
																	y00<-(yy[z[,bb]==1]-m0[bb])/sqrt(rexp(n00,1))
																	Z<-c()
																		for (i in 1:n00)
																		{
																			up<-exp(-.5)/(sqrt(2*pi)*abs(y00[i]))
																			j<-1
																				while (j<2)
																				{
																					w<-rweibull(1,a0[bb],1)
																					ex<-dnorm(y00[i],0,sqrt(2)*s0[bb]/w)
																					if (runif(1)<ex/up)
																					{
																						Z[i]<-w
																						j<-j+1
																					}
																					if (ex==0)
																					{
																						Z[i]<-sqrt(2)*s0[bb]/abs(y00[i])
																						j<-j+1
																					}
																				}
																		}

																	f<-function(v){sum(-log(v[1])-(v[1]-1)*log(Z)+Z^v[1])}
																	out<-suppressWarnings(nlm(f,v<-c(a0[bb]),hessian=FALSE))
																	a11[bb,rrr]<-out$estimate[]
																	if (a11[bb,rrr]>2)
																	{
																		a11[bb,rrr]<-1.99
																	}
																		a12[bb,rrr]<-sqrt(sum(y00^2*Z^2,na.rm=TRUE)/(2*n00))
																}
																alpha.matrix[r,bb]<-mean(a11[bb,])
																sigma.matrix[r,bb]<-mean(a12[bb,])
																a0[bb]<-alpha.matrix[r,bb]
																s0[bb]<-sigma.matrix[r,bb]
															}
		}
	for (i in 1:length(yy)){clustering[i]<-which(z[i,]==1)[1]}
	return(list(omega=apply(p.matrix[(m-50):m,],2,mean),alpha=apply(alpha.matrix[(m-50):m,],2,mean),sigma=apply(sigma.matrix[(m-50):m,],2,mean),mu=apply(mu.matrix[(m-50):m,],2,mean),cluster=clustering))
}
ufitstab.cauchy<-function(y,beta0,sigma0,mu0,param)
{
	stopifnot(-1<=beta0,beta0<=1,length(beta0)==1,0<=sigma0,length(mu0)==1,length(param)==1,param %in% 0:1)
	n<-length(y)
	ep11<-ep1<-ep2<-ep22<-ep222<-ep11<-ep0<-m<-k<-t1<-t2<-c()
	t2[1]<-sigma0*beta0
	t1[1]<-max(sigma0*(1-abs(beta0)),.001)
	m[1]<-mu0
	nn<-1000;mm<-950
		for(j in 1:nn)
		{
				for(i in 1:n)
				{
					p2<-rstable(3000,1,1,1,0,1)
					k<-(y[i]-m[j]-t2[j]*p2)/t1[j]
					tt<-0;r<-0;
					dy<-2^(tt/2)*gamma(tt/2+1)/(pi*t1[j])*mean(p2^r/(1+k^2)^(tt/2+1))
					tt<-2;r<-2;
					ep22[i]<-2^(tt/2)*gamma(tt/2+1)/(pi*t1[j])*mean(p2^r/(1+k^2)^(tt/2+1))/dy
					tt<-2;r<-1;
					ep2[i]<-2^(tt/2)*gamma(tt/2+1)/(pi*t1[j])*mean(p2^r/(1+k^2)^(tt/2+1))/dy
					tt<-2;r<-0;
					ep1[i]<-2^(tt/2)*gamma(tt/2+1)/(pi*t1[j])*mean(p2^r/(1+k^2)^(tt/2+1))/dy
				}
			m[j]<-(sum(y*ep1,na.rm=TRUE)-t2[j]*sum(ep2,na.rm=TRUE))/sum(ep1,na.rm=TRUE)
			m[j+1]<-m[j]
			t1[j]<-sqrt((sum((y-m[j])^2*ep1)+t2[j]^2*sum(ep22)-2*t2[j]*sum((y-m[j])*ep2))/n)
			t1[j+1]<-t1[j]
			t2[j]<-sum((y-m[j])*ep2)/sum(ep22)
			t2[j+1]<-t2[j]
		}
	beta.hat<-uniroot(function(p) mean(t2[mm:nn])/mean(t1[mm:nn])-p/(1-abs(p)),c(-.999999,.999999))$root
	sigma.hat<-min(c(mean(t2[mm:nn])/beta.hat, mean(t1[mm:nn])/(1-abs(beta.hat))))
	mu.hat<-mean(m[mm:nn])
		if (param==0)
		{
		return(list(beta=beta.hat,sigma=sigma.hat,mu=mu.hat))
		}
		else
		{
		return(list(beta=beta.hat,sigma=sigma.hat,mu=(mu.hat-2/pi*beta.hat*sigma.hat*log(sigma.hat))))
		}
}
ufitstab.cauchy.mix<-function(y,k,omega0,beta0,sigma0,mu0)
{
	stopifnot(-1<=beta0,beta0<=1,length(beta0)==k,0<=sigma0,length(mu0)==k,length(sigma0)==k,sum(omega0)==1,0<omega0,omega0<1)
	  n <- length(y)
    MM <- 1300
    NN <- 1500
    m <- 1500
    N <- 2000
    estim.matrix <- array(0, dim = c(4, ncol = k, m))
    mu.matrix <- matrix(m * k, ncol = k, nrow = m)
    t1.matrix <- matrix(m * k, ncol = k, nrow = m)
    t2.matrix <- matrix(m * k, ncol = k, nrow = m)
    p.matrix <- matrix(m * k, ncol = k, nrow = m)
    tau.matrix <- matrix(m * k, ncol = k, nrow = n)
    e1ij <- matrix(m * k, ncol = k, nrow = n)
    e2ij <- matrix(m * k, ncol = k, nrow = n)
    e3ij <- matrix(m * k, ncol = k, nrow = n)
    e4ij <- matrix(m * k, ncol = k, nrow = n)
    mu.matrix[1, ] <- mu0
    p.matrix[1, ] <- omega0
    t2.matrix[1, ] <- sigma0 * beta0
    t1.matrix[1, ] <- sigma0 * (1 - abs(beta0))
    p0 <- p.matrix[1, ]
    t1 <- t1.matrix[1, ]
    t2 <- t2.matrix[1, ]
    m0 <- mu.matrix[1, ]
    b0 <- s0 <- dy <- c()
    clustering<-rep(0,n)
		for (r in 2:m)
		{
			for (i in 1:n)
			{
				for (bb in 1:k)
				{
					p2<-rstable(N,1,1,1,0,1)
					t1[bb]<-ifelse (abs(t1[bb])< 0.000001,t1[bb]<-.000001,t1[bb])
					kk<-(y[i]-m0[bb]-t2[bb]*p2)/t1[bb]
					tt<-0;rr<-0;
					dy[bb]<-2^(tt/2)*gamma(tt/2+1)/(pi*t1[bb])*mean(p2^rr/(1+kk^2)^(tt/2+1),na.rm=TRUE)
					tt<-2;rr<-2;
					e4ij[i,bb]<-2^(tt/2)*gamma(tt/2+1)/(pi*t1[bb])*mean(p2^rr/(1+kk^2)^(tt/2+1),na.rm=TRUE)/dy[bb]
					tt<-2;rr<-1;
					e3ij[i,bb]<-2^(tt/2)*gamma(tt/2+1)/(pi*t1[bb])*mean(p2^rr/(1+kk^2)^(tt/2+1),na.rm=TRUE)/dy[bb]
					tt<-2;rr<-0;
					e2ij[i,bb]<-2^(tt/2)*gamma(tt/2+1)/(pi*t1[bb])*mean(p2^rr/(1+kk^2)^(tt/2+1),na.rm=TRUE)/dy[bb]
				}
						for (aa in 1:k)
						{
							e1ij[i,aa]<-p0[aa]*dy[aa]/sum(p0*dy,na.rm=TRUE)
						}
			}
								for (ii in 1:k)
								{
									mu.matrix[r,ii]<-(sum(y*e2ij[,ii]*e1ij[,ii],na.rm=TRUE)-t2[ii]*sum(e3ij[,ii]*e1ij[,ii],na.rm=TRUE))/
									sum(e1ij[,ii]*e2ij[,ii],na.rm=TRUE)
									m0[ii]<-mu.matrix[r,ii]
									t1.matrix[r,ii]<-sqrt((sum((y-m0[ii])^2*e2ij[,ii]*e1ij[,ii],na.rm=TRUE)-2*t2[ii]*
									sum((y-m0[ii])*e3ij[,ii]*e1ij[,ii],na.rm=TRUE)+t2[ii]^2*sum(e4ij[,ii]*e1ij[,ii],na.rm=TRUE))/sum(e1ij[,ii],na.rm=TRUE))
									t1[ii]<-t1.matrix[r,ii]
									t2.matrix[r,ii]<-sum((y-m0[ii])*e3ij[,ii]*e1ij[,ii],na.rm=TRUE)/sum(e4ij[,ii]*e1ij[,ii],na.rm=TRUE)
									t2[ii]<-t2.matrix[r,ii]
									p.matrix[r,ii]<-sum(e1ij[,ii])/n
									p0[ii]<-p.matrix[r,ii]
								}
										z<-matrix(0,ncol=k,nrow=n)
											for (j in 1:n)
											{
												max<-e1ij[j,1]
												uu<-1
													for (ii in 2:k)
													{
														if (e1ij[j,ii]> max)
														{
															max<-e1ij[j,ii]
															uu<-ii
														}
													}
												z[j,uu]<-1
											}
																for (aa in 1:k)
																{
																	b0[aa]<-suppressWarnings(uniroot(function(p) t2[aa]/t1[aa]-p/(1-abs(p)),c(-.9999999,.9999999))$root)
																	s0[aa]<-t1[aa]/(1-abs(b0[aa]))
																}
			estim.matrix[1,,r]<-p0
			estim.matrix[2,,r]<-b0
			estim.matrix[3,,r]<-s0
			estim.matrix[4,,r]<-m0
		}
	estim.matrix[1,,1]<-omega0
	estim.matrix[2,,1]<-beta0
	estim.matrix[3,,1]<-sigma0
	estim.matrix[4,,1]<-mu0
	for (i in 1:length(y)){clustering[i]<-which(z[i,]==1)[1]}
	return(list(omega=apply(estim.matrix[1,,(MM:NN)],1,mean),beta=apply(estim.matrix[2,,(MM:NN)],1,mean),sigma=apply(estim.matrix[3,,(MM:NN)],1,mean),mu=apply(estim.matrix[4,,(MM:NN)],1,mean),cluster=clustering))
}
ufitstab.skew<-function(y,alpha0,beta0,sigma0,mu0,param)
{
	stopifnot(length(y)>=4,0<alpha0,alpha0<=2,alpha0!=1,length(alpha0)==1,-1<=beta0,beta0<=1,length(beta0)==1,0<=sigma0,length(sigma0)==1,length(param)==1, param %in% 0:1)
	n<-length(y)
	M<-100;N0<-100;N<-120;
	sss<-m<-s<-a<-b<-ep11<-ep1<-ep2<-ep3<-a.estim<-c()
	m[1]<-mu0;s[1]<-sigma0;b[1]<-beta0;a[1]<-alpha0
		for(j in 1:N)
		{
			pi1<-matrix(suppressWarnings(rstable(M^2,a[j]/2,1,(cos(pi*a[j]/4))^(2/a[j]),0,1)),M,M)
			pi2<-matrix(suppressWarnings(rstable(M^2,a[j],1,1,0,1)),M,M)
				for(i in 1:n)
				{
					ss<-dnorm(y[i],m[j]-s[j]*b[j]*tan(pi*a[j]/2)+sign(b[j])*abs(b[j])^(1/a[j])*s[j]*pi2,sd=sqrt(2*pi1)*s[j]*(1-abs(b[j]))^(1/a[j]))
					ep11[i]<-mean(ss,na.rm=TRUE)
					dy<-ep11[i]
					ep1[i]<-mean(ss/pi1,na.rm=TRUE)/dy
					ep2[i]<-mean(ss*pi2/pi1,na.rm=TRUE)/dy
					ep3[i]<-mean((ss*pi2^2/pi1),na.rm=TRUE)/dy
				}
						m[j+1]<-(sum((y+s[j]*b[j]*tan(pi*a[j]/2))*ep1,na.rm=TRUE)-s[j]*sign(b[j])*abs(b[j])^(1/a[j])*sum(ep2,na.rm=TRUE))/sum(ep1,na.rm=TRUE)
						fs<-function(p){.5*sum((y-m[j])^2*ep1)/(abs(p)^3*(1-abs(b[j]))^(2/a[j]))+.5*b[j]*(tan(pi*a[j]/2))*sum((y-m[j])*ep1)/(p^2*(1-abs(b[j]))^(2/a[j]))-.5*sign(b[j])*abs(b[j])^(1/a[j])*sum((y-m[j])*ep2)/(p^2*(1-abs(b[j]))^(2/a[j]))-n/abs(p)}
						s[j+1]<-suppressWarnings(uniroot(fs,c(0.000000001,10000000))$root)
						fb<-function(p){.25*sum((y-m[j])^2*ep1)/(s[j]^2*(1-abs(p))^(2/a[j]))+.25*p^2*(tan(pi*a[j]/2))^2*sum(ep1)/((1-abs(p))^(2/a[j]))+.25*abs(p)^(2/a[j])*sum(ep3)/((1-abs(p))^(2/a[j]))+.5*p*(tan(pi*a[j]/2))*sum((y-m[j])*ep1)/(s[j]*(1-abs(p))^(2/a[j]))-.5*abs(p)^(1+1/a[j])*(tan(pi*a[j]/2))*sum(ep2)/((1-abs(p))^(2/a[j]))-.5*sign(p)*abs(p)^(1/a[j])*sum((y-m[j])*ep2)/(s[j]*(1-abs(p))^(2/a[j]))+sum(1/a[j]*(log(1-abs(p))))}
						b[j+1]<-suppressWarnings(optimize(fb,lower=-.999999,upper=.999999))$minimum[[1]]
						st<-suppressWarnings(rstable(n,a[j],1,1,0,1))
						sss[j]<-s[j]*(1+abs(b[j]))^(1/a[j])
						yy<-(y-m[j]-s[j]*sign(b[j])*(abs(b[j]))^(1/a[j])*st)/sss[j]
				for (ii in 1:20)
				{
 					nn<-length(yy)
 					Z1<-c()
 					Z<-c()
	 				yyy<-(yy)/(sqrt(rexp(nn,1)))
 						for (i in 1:nn)
						{
 							y0<-yyy[i]
 							jj<-1
 								while (jj<2)
								{
 									tt<-rweibull(1,a[j],1)
 									ra<-exp(-.5)/(sqrt(2*pi)*abs(y0))
	 								u<-runif(1)
 										if (u<dnorm(y0,0,sqrt(2)/tt)/ra)
 										{
											Z1[jj]<-tt
 											jj<-jj+1
										}
								}
 							Z[i]<-Z1
						}
 					f<-function(p){sum(-log(p[1])-(p[1]-1)*log(Z)+Z^p[1])}
 					a.hat<-suppressWarnings(nlm(f, p<-c(a[j])))$estimate[]
 					if (a.hat>1.99){a.hat<-1.98}
 					a.estim[ii]<-a.hat
				}
			a[j+1]<-mean(a.estim[10:20])
		}
	alpha.hat<-mean(a[N0:N])
	mu.hat<-mean(m[N0:N])
	sigma.hat<-mean(s[N0:N])
	w<-function(p){-sum(log(dstable(y,alpha.hat,p,sigma.hat,mu.hat,0)))}
	beta.hat<-suppressWarnings(optimize(w,c(-1,1)))$minimum
		if(param==0)
		{
		return(list(alpha=alpha.hat,beta=beta.hat,sigma=sigma.hat,mu=mu.hat))
		}
		else
		{
		return(list(alpha=alpha.hat,beta=beta.hat,sigma=sigma.hat,mu=mu.hat-beta.hat*sigma.hat*tan(pi*alpha.hat/2)))
		}
}
ufitstab.ustat<-function(x)
{
	n<-length(x)
	s1<-s2<-0
		for (i in 1:(n-1))
		{
			for (j in (i+1):n)
			{
				s1<-s1+(log(abs(x[i]+x[j]))-(log(abs(x[i]))+log(abs(x[j])))/2)/log(2)
				s2<-s2+(1+0.57721566/log(2))*(log(abs(x[i]))+log(abs(x[j])))/2-0.57721566/log(2)*log(abs(x[i]+x[j]))+0.57721566
			}
		}
	return(list(alpha=n*(n-1)/(2*s1),sigma=exp(2*s2/(n*(n-1)))))
}
