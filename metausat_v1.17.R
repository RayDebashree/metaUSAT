###------------- Code for metaUSAT for testing association of multiple continuous phenotypes with SNPs using GWAS summary statistics-----------------
#
# Ref 1: "Methods for Meta-analysis of Multiple Traits using GWAS Summary Statistics". Genetic Epidemiology, 2017. DOI: 10.1002/gepi.22105
# Ref 2: "USAT: A Unified Score-based Association Test for Multiple Phenotype-Genotype Analysis". Genetic Epidemiology, 40(1):20-34, 2016.
#
# metaUSAT: A linear combination of SSU and MANOVA score-type statistics using GWAS summary statistics (univariate Wald test statistics as summary statistics)
# T_MANOVA ~ chisq(k), where k = no of traits, (no. of SNPs=1)
# T_SSU ~ a*chisq(d) + b
# T_rho = rho*T_MANOVA + (1-rho)*T_SSU
# optimal rho found using grid search over [0,1]
# T_metaUSAT = min_rho p_rho, where p_rho is pvalue of T_rho
# calculates the appropriate p-value of T_metaUSAT using a one-dimensional numerical integration

##--------------------------------------------- Version 1.17 (dated December 11, 2017) --------------------------------------------------
# Corresponding Author: Debashree Ray, Ph.D. <dray@jhu.edu>

#### load the necessary libraries
library(CompQuadForm)
library(survey)	
library(minqa)
library(psych)

message("=====================================")
message("     metaUSAT v1.17 is loaded")
message("=====================================")
message("If you use this software, please cite:")
message("Ray et al.(2017) Methods for Meta-analysis of Multiple Traits using")
message("    GWAS Summary Statistics. Genetic Epidemiology, DOI: 10.1002/gepi.22105")
message("-------------------------------------")
message("metaUSAT is based on USAT:")
message("Ray et al.(2016) USAT: A Unified Score-based Association Test for Multiple")
message("    Phenotype-Genotype Analysis. Genetic Epidemiology, 40(1):20-34")
message("-------------------------------------")


############################################
eps<-1e-13	# a threshold for checking non-positive values
# function for Moore-Penrose Inverse
mpinv <- function(A, eps = 1e-13, power=-1) {
        s<-eigen(A)
        e<-s$values
        V<-s$vectors
        e[e>eps] <- 1/(e[e > eps])^abs(power)
        return(V%*%diag(e)%*%t(V))
}

##----------------------- Begin: Functions borrowed from Dr. Baolin Wu's website -----------------------
################ functions required for calculating higher order cumulants in order to do higher order moment matching (Wu and Pankow 2015)
#####
cum2mnc = function(kappa){
### convert cumulants to non-central moments
###    recursive formula produces as many cumulants as moments
###    References: Kenneth Lange: Numerical Analysis for Statisticians, 2nd ed. Page 15
  N = length(kappa)+1
  mc = rep(0, N); mc[1] = 1
  for(k in 1:(N-1)){
    mc[k+1] = sum(choose(k-1, 0:(k-1))*kappa[k:1]*mc[1:k])
  }
  return(mc[-1])
}
mnc2mc = function(mnc){
### convert non-central to central moments, uses recursive formula
  N = length(mnc)
  mc = rep(0,N); mc[1] = 0
  s1 = rep(c(1,-1), N)
  mnc = c(1,mnc)
  for(k in 1:(N-1)){
    mc[k+1] = sum( choose(k+1, 0:(k+1))*s1[(k+2):1]*mnc[1:(k+2)]*mnc[2]^((k+1):0) )
  }
  return(mc)
}
#### non-central chi-square cumulants
chisq.cum = function(k, lam, N){
### k: DF; lam: ncp 
  ik = 1:N
  2^(ik-1)*gamma(ik)*(k+ik*lam)
}
## 1-DF chisq mix cumulants
chi1sqm.cum = function(lam, N){
### lam: weight coef
  ik = 1:N
  a1 = 2^(ik-1)*gamma(ik)
  cl = rep(0, N)
  for(i in 1:N) cl[i] = a1[i]*sum(lam^i)
  cl
}

## match higher moments
wu.lambda = function(lam, N=12){	# Nth cumulant (Wu Pankow suggests using N=12)
  cl = chi1sqm.cum(lam, N)
  muQ = cl[1]; sigmaQ = sqrt(cl[2])
  a1 = mnc2mc(cum2mnc(cl))
  a1 = a1/sqrt(a1[2])^(1:N)  
  f1 = function(xpar){
    k = exp(xpar[1])
    v = xpar[2]
    a2 = mnc2mc(cum2mnc(chisq.cum(k,v,N)))
    a2 = a2/sqrt(a2[2])^(1:N)  
    (a1[N-1]-a2[N-1])^2 + (a1[N]-a2[N])^2
  }
  tmp = bobyqa(c(0,1), f1, lower=c(-Inf,0),upper=c(Inf,Inf))
  xpar = tmp$par
  l = exp(xpar[1])
  d = xpar[2]
  if(f1(c(xpar[1],0))<=tmp$fval){
    d=0
    f.1 = function(xpar) f1(c(xpar,0))
    l = exp(bobyqa(xpar[1], f.1)$par)
  }
  muX = l+d; sigmaX = sqrt(chisq.cum(l,d,N=2)[2])
  list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
}
wu.pval = function(Q.all, lambda, N=12){
  param = wu.lambda(lambda,N)
  Q.Norm = (Q.all - param$muQ)/param$sigmaQ
  Q.Norm1 = Q.Norm*param$sigmaX + param$muX
  pchisq(Q.Norm1, df = param$l,ncp=param$d, lower.tail=FALSE)
}
##----------------------- End: Functions borrowed from Dr. Baolin Wu's website -----------------------

###### Functions for estimating correlation matrices
cor.pearson<-function(Z.matrix, P.matrix, p.threshold=1e-5)
{
  # estimating correlation
  row.exclude<-which( apply(P.matrix, MARGIN = 1, function(x) any(x < p.threshold)) == TRUE )
  Z.matrix<-Z.matrix[-row.exclude,]
  R<-cor(Z.matrix)
  return(R)
}

cor.tetrachor<-function(Z.matrix)
{
  # estimating correlation
  Z.matrix<-as.matrix(Z.matrix)
  Z.matrix[which(Z.matrix<0)]<-0
  Z.matrix[which(Z.matrix>0)]<-1
  R<-tetrachoric(Z.matrix)$"rho"
  return(R)
}

############################################
############################################

metausat<-function(Z, R, weights=1, metamanova=FALSE, AbsTol=.Machine$double.eps^0.8) 
{		# Z: vector of Z-scores
		# R: estimated correlation matrix of the Z-scores
		# weights: vector of weights for the Z-scores (Default 1)
		# metamanova: if TRUE, metaMANOVA statistic and p-value are output as well

   ## Checks
   k<-length(Z)		# no. of studies*no. of traits in each study
   if(k<2) stop("This method is meant for meta-analysis of multiple traits from 1 or more studies or for meta-analysis of single trait from multiple studies. Information on <2 traits/studies provided.")
   if(nrow(R)!=k | ncol(R)!=k) stop("Order of covariance matrix does not match with Z-vector.")
   if(weights!=1){
     if(length(weights)!=k) stop("Length of vector 'weights' must be same as length of Z-vector.")
     Z<-Z*weights
     W<-diag(weights)
     R<-W%*%R%*%W
   }
   Z<-matrix(Z, nrow=k, ncol=1)

   ############ MetaSSU #####################
   #SumSqU:
     #it's possible Z=0 and Cov(Z)=0:
     if (all(abs(Z)<1e-20)) { pssu<-1; stop("Check Z vector: all entries very close to 0.") }else{
       Tssu<- t(Z) %*% Z
       ##distr of Tssu is sum of c_r Chisq_1 (c_r is rth eigen value of CovZ):
         eigR<-eigen(R, only.values=TRUE)$values
         dfman<-sum(eigR>eps)           # number of positive eigen values (rank of the covariance matrix)
       ###approximate the distri by alpha Chisq_d + beta:
         alpha1<-as.double( sum(eigR^3)/sum(eigR^2) )
         beta1<-as.double( sum(eigR) - (sum(eigR^2))^2/(sum(eigR^3)) )
         d1<-as.double( (sum(eigR^2))^3/(sum(eigR^3))^2 )
         pssu<-as.double(pchisq((Tssu-beta1)/alpha1, d1, lower.tail=FALSE))    #p-value for SSU test
    }

   ############ metaMANOVA: MANOVA type standard chi-squared test for meta-analysis #####################
   Tman <- t(Z) %*% mpinv(R) %*% Z
   pman <- pchisq(Tman, df=dfman, lower.tail=FALSE)

   #################################### metaUSAT: combined statistic  ###################################
     Ts = Tssu
     Tc = Tman
     pval.rho<-q.rho<-To.rho<-NULL
     rho.set<-seq(0,1,0.1)	# set of rho values (weights) to consider in the linear combination
     for(rho in rho.set)
     {
       	To<-(1-rho)*Ts + rho*Tc
       	To.rho<-c(To.rho,To)            # value of unified statistics for different values of rho
       	eig<- (1-rho)*eigR+rho
      	me <- sum(eig>eps)
        dav.qf<-davies(q=To, lambda=eig[1:me], h=rep(1,me),acc=AbsTol,lim=1e+4)
        if(dav.qf$ifault!=0 | dav.qf$Qq<0 | dav.qf$Qq==0) pval.qf<-pchisqsum(x=To, df=rep(1,me), a=eig[1:me], lower.tail=FALSE, method="saddlepoint") else pval.qf<-dav.qf$Qq
       	pval.rho<-c(pval.rho, pval.qf)
     }
     pval.opt<-min(pval.rho)
     rho.opt<-0+(which(pval.rho==pval.opt)-1)*0.1

     # calculate the p-value using a one-dimensional numerical integration
     # pval.opt is our test statistic T, whose p-value needs to be computed
        # first find the (1-T)th quantiles of the linear combination of chi-sq distributions for every rho
        T<-pval.opt
        n.r<-length(rho.set)
        qminrho<-rep(0,n.r)
        c1 <- rep(0, 4)
        for (i in 1:n.r) {
               	rho <- rho.set[i]
               	eig<- (1-rho)*eigR+rho
		if(T>1e-4) {
               	   c1[1] <- sum(eig)
               	   c1[2] <- sum(eig^2)
               	   c1[3] <- sum(eig^3)
                   c1[4] <- sum(eig^4)
               	   muQ <- c1[1]
               	   varQ <- 2*c1[2]
               	   s1 <- c1[3]/c1[2]^1.5
                   s2 <- c1[4]/c1[2]^2
               	   if(s1^2>s2){
                     a=1/(s1-sqrt(s1^2-s2))
                     ncp.d=s1*a^3-a^2
                     df.l=a^2-2*ncp.d
               	   }else{
                     a=1/s1 ; ncp.d=0 ; df.l=c1[2]^3/c1[3]^2
               	   }
                       	mu.chi<-df.l+ncp.d
                       	sig.chi<-sqrt(2)*a
                       	q.org <- qchisq(T, df=df.l, ncp=ncp.d, lower.tail=FALSE)
                       	q.q <- ( (q.org-mu.chi)/sig.chi )*sqrt(varQ) + muQ
	 	}else {
		  wu.out<-wu.lambda(lam=eig, N=6)
			mu.chi<-wu.out$muX
			sig.chi<-wu.out$sigmaX
                       	q.org <- qchisq(T, df=wu.out$l, ncp=wu.out$d, lower.tail=FALSE)
                       	q.q <- ( (q.org-mu.chi)/sig.chi )*(wu.out$sigmaQ) + wu.out$muQ
		}
               	qminrho[i] <- q.q
        }

     # the required pvalue is 1-P(chi2_k<quant|Ts)*f(Ts), where quant=min_{rho}((qminrho-(1-rho)*Ts)/rho) and
     # f(Ts) is the density of Ts, where Ts~alpha1*chi2_(d1)+beta1
     # i.e., pvalue is P(chi2_k>=quant|Ts)*f(Ts)
                f.Liu.chi.upp<-function (x)
                {
                    temp <- (qminrho - rho.set*x)/(1 - rho.set)
                    temp.min <- min(temp)
                    re <- pchisq((temp.min-beta1)/alpha1, df=d1, ncp=0, lower.tail=FALSE) * dchisq(x, df=dfman)   # using the scaled and shifted chi-sq distn of SSU
                    return(re)
                }
       integfunc<-integrate(Vectorize(f.Liu.chi.upp),0,Inf,abs.tol=AbsTol,stop.on.error=FALSE)
       if(integfunc$message!="OK"){
	  message(paste(integfunc$message," at AbsTol=",AbsTol,"; NA assigned to p.metausat",sep=""))
	  pval.T<-NA
       }else pval.T<-as.double(integfunc$value)

 if(metamanova)
   return(list(T.metamanova=drop(Tman), p.metamanova=drop(pman), T.metausat=pval.opt, omg.opt=rho.opt[1], p.metausat=pval.T, AbsTol=AbsTol, error.msg=integfunc$message))
 else
   return(list(T.metausat=pval.opt, omg.opt=rho.opt[1], p.metausat=pval.T, AbsTol=AbsTol, error.msg=integfunc$message))

}

