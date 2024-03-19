#' Two-parameter Maxwell-Boltzmann Distribution
#'
#' Density, distribution function, quantile function, and random generation for the two-parameter Maxwell-Boltzmann distribution. This probability distribution has a location parameter \eqn{\mu}, and a scale parameter \eqn{\sigma}.

#' @rdname infMax
#' @param x,q numeric vector of qunatiles.
#' @param mu vector of location parameter. Default is 0.
#' @param sig vector of scale parameter. Default is 1.

#' @details
#' If location and scale parameter are not specified, above functions will assume default value as 0 and 1, respectively.
#'
#' The probability density function (pdf) of two parameter Maxwell-Boltzmann distribution (popularly known as Maxwell distribution) is given by
#' \deqn{f(x|\mu,\sigma) = \frac{4}{\sigma\Gamma(1/2)}\left(\frac{x-\mu}{\sigma}\right)^2\exp\left\{-\left(\frac{x-\mu}{\sigma}\right)^2\right\}, \quad x>\mu, \ \sigma >0}
#' Where \eqn{\mu} and \eqn{\sigma} are the location and scale parameter respectively.
#'
#' See \href{ https://doi.org/10.1007/s42519-022-00270-y}{Statistical Intervals for Maxwell Distributions} for more details.
#' @returns \code{dMaxwell} gives the density, \code{pMaxwell} gives the cumulative distribution function, \code{qMaxwell} gives quantile functions, and \code{rMaxwell} generates random numbers from the two-parameter Maxwell distribution.
#'
#' The length of the result is determined by \code{n} for \code{rMaxwell},  and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#'The numerical arguments other than \code{n} are recycles to the length of the result. Only the first element of the logical vector is being used.
#'
#' In case of \code{pMaxwell} and and \code{qMaxwell}, if there is no specification about the logical argument \code{lower.tail}, they assume it's \code{True} by default, which means probabilities are \eqn{P[X \leq x]}.
#' @references Chowdhury, F. and Krishnamoorthy, K. (2022). Statistical intervals for Maxwell distributions. Journal of Statistical Theory and Practice. 16, 45.
#'
#' @keywords dMaxwell
#' @export
#' @examples
#' #For calculating densities
#' dMaxwell(2,1,2)
#' dMaxwell(5, 0, 5)
#'
#' # In order to plot the pdf with a fixed location parameter (say 0) with varying scale parameter
#' x<-seq(0,5,length=100)
#' sigmas<-c(.5,1,2)
#' plot(x,dMaxwell(x,0,.5), type="l",col="blue",xlab="x",ylab="Probability density f(x)")
#' lines(x,dMaxwell(x,0,1),col="red")
#' lines(x,dMaxwell(x,0,2),col="green")
#' legend("topright", legend = paste("\u03C3 =", sigmas), col = c("blue", "red", "green"), lty = 1, lwd = 2)
#'
dMaxwell<-function(x,mu=0,sig=1)
{
  if(max(length(x),length(mu),length(sig))==length(x)){
    x<-x;mu<-rep_len(mu,length(x));sig<-rep_len(sig,length(x))

  }
  if(max(length(x),length(mu),length(sig))==length(mu)){
    x<-rep_len(x,length(mu));mu<-mu;sig<-rep_len(sig,length(mu))
  }
  if(max(length(x),length(mu),length(sig))==length(sig)){
    x<-rep_len(x,length(sig));mu<-rep_len(mu,length(sig));sig<-sig
  }
  x<-pmax(mu,x)
  sigma<-c()
  for(i in 1:length(sig))
  {
    sigma[i]<-ifelse(sig[i]<=0,NaN,sig[i])
  }
  d<-4/(sigma*sqrt(pi))*((x-mu)/sigma)^2*exp(-((x-mu)/sigma)^2)
  return(d)

}



#' @rdname infMax
#' @keywords pMaxwell
#' @export
#' @examples
#' # For calculating cumulative probabilities
#' pMaxwell(2,1,2)
#' pMaxwell(2,0,2,lower.tail=F)
#'
pMaxwell<-function(q, mu=0, sig=1, lower.tail=TRUE)
{
  if(max(length(q),length(mu),length(sig))==length(q)){
    q<-q;mu<-rep_len(mu,length(q));sig<-rep_len(sig,length(q))

  }
  if(max(length(q),length(mu),length(sig))==length(mu)){
    q<-rep_len(q,length(mu));mu<-mu;sig<-rep_len(sig,length(mu))
  }
  if(max(length(q),length(mu),length(sig))==length(sig)){
    q<-rep_len(q,length(sig));mu<-rep_len(mu,length(sig));sig<-sig
  }
  q<-pmax(mu,q)
  sigma<-c()
  for(i in 1:length(sig))
  {
    sigma[i]<-ifelse(sig[i]<=0,NaN,sig[i])
  }
  z<-(q-mu)^2/sigma^2
  if(lower.tail == TRUE){
    return(pgamma(z,1.5))
  }
  else{
    return(pgamma(z,1.5,lower.tail = F))
  }

}




#' @rdname infMax
#' @param p vector of probabilities.
#' @param lower.tail logical; If TRUE(default), probabilities are \eqn{P[ X \leq x]}, otherwise, \eqn{P[X \gt x]}.
#' @keywords qMaxwell
#' @export
#' @examples
#' # For calculating quantiles
#' qMaxwell(.5,1,2)
#' qMaxwell(.5,0,2,lower.tail=F)
#'
qMaxwell<-function(p,mu=0,sig=1, lower.tail=TRUE)
{

  if(max(length(p),length(mu),length(sig))==length(p)){
    p<-p;mu<-rep_len(mu,length(p));sig<-rep_len(sig,length(p))

  }
  if(max(length(p),length(mu),length(sig))==length(mu)){
    p<-rep_len(p,length(mu));mu<-mu;sig<-rep_len(sig,length(mu))
  }
  if(max(length(p),length(mu),length(sig))==length(sig)){
    p<-rep_len(p,length(sig));mu<-rep_len(mu,length(sig));sig<-sig
  }
  pr<-c()
  sigma<-c()
  for(i in 1:length(p))
  {
    pr[i]<-ifelse(p[i]<0 | p[i]>1,NaN,p[i])
    sigma[i]<-ifelse(sig[i]<=0,NaN,sig[i])
  }
  if(lower.tail==TRUE){
    x<-mu+(sigma*sqrt(qgamma(pr,shape = 1.5,scale = 1)))
    return(x)
  }
  else{
    x<-mu+(sigma*sqrt(qgamma(pr,shape = 1.5,scale = 1, lower.tail = F)))
    return(x)
  }

}

#' @rdname infMax
#' @param n number of observations; If \code{length(n)>1}, the length is taken to be the number required.

#' @keywords rMaxwell
#' @export
#' @examples
#' # For generating random numbers
#' rMaxwell(5,1,2)
#' rMaxwell(25,0,2)
#'

#' @export
rMaxwell<-function(n, mu=0,sig=1)
{

  n<-ifelse(length(n)>1,length(n),n)
  if(n<=0){
    return(warning("invalid arguments"))
  }
  if(length(mu)<n | length(mu)>n){
    mu<-rep_len(mu,length.out = n)
  }
  if(length(sig)<n | length(sig)>n){
    sig<-rep_len(sig,length.out = n)
  }
  sigma<-c()
  for(i in 1:length(sig))
  {
    sigma[i]<-ifelse(sig[i]<=0,NaN,sig[i])
  }

  z<-rgamma(n, 3/2)
  x<-mu +sqrt(z)*sigma
  return(x)

}



#' Moment Estimates (MEs) calculating function of two-parameter Maxwell distributions
#'
#' This function allows you to calculate the MEs of location & scale parameter respectively
#' @param x a numeric vector of data input.

#' @keywords MEs
#' @export
#' @examples
#'  # MEs(x)

MEs<-function(x)
{
  n<-length(x)
  s<-sqrt((1/n)*sum((x-mean(x))^2))
  sig<-s*sqrt((2*pi)/((3*pi)-8))
  mu<-mean(x)-((2/sqrt(pi))*sig)
  return(c(mu,sig))
}


#' Maximum likelihood estimates (MLEs) calculating function of two-parameter Maxwell distributions
#'
#' This function allows you to calculate the MLEs of location & scale parameter respectively
#' @param x a numeric vector of data input.

#' @keywords MLEs
#' @export
#' @examples
#'  # MLEs(x)
MLEs<-function(x){
  n<-length(x); Pf<-1-(1-.9999)^(1/n)
  qp<-sqrt(qgamma(Pf, 3/2))
  sigsq.x<-(n-1)*var(x)/n
  xb<-mean(x); s<-sqrt(sigsq.x); xmin<-min(x)
  sigh<-s*sqrt(2*pi/(3*pi-8))
  tconst<-sigh*qp
  fn<-function(v){
    gt<-n*(xb-v)-.6666667*(sigsq.x+(xb-v)^2)*sum(1/(x-v))
  }
  muh<-uniroot(fn, c(xmin-tconst, xmin))[[1]]
  sigh<-sqrt(2*(sum((x-muh)^2))/3/n)
  return(c(muh,sigh))
}

#' Q-Q plot function for two-parameter Maxwell distributions
#'
#' This Q-Q plot function will calculate theoretical Maxwell quantiles based on moment estimates and plot it against sample quantiles in order check goodness of fit.
#' @param x a numeric vector of data input, should be with length more than one.

#' @keywords Q_Q_plot
#' @export
#' @examples
#'  # Q_Q_plot(x)
Q_Q_plot<-function(x)
{
  j<-seq(1,length(x),1)
  const<-(1.5-(4/pi))
  sigma_hat<-sqrt(var(x)/const)
  mu_hat<-mean(x)-((2*sigma_hat)/sqrt(pi))
  prop_data<-(j-0.5)/length(x)
  q_j<-mu_hat+(sqrt(qgamma(prop_data,1.5))*sigma_hat)
  x_j<-sort(x)
  plot(q_j,x_j, col="red",xlab = "Theoretical Maxwell Quantile", ylab = "Sample quantile")
  abline(0,1)
}

#' Confidence intervals (MLEs based) for two-parameter Maxwell distributions
#'
#' This function computes the required percentiles to compute confidence interval for the mean (MLEs based) of Maxwell distributions.

#'
#' @param cl confidence level as in decimal
#' @param n a numeric vector, the number of observation in each sample
#' @param nr numeric vector length one, number of simulation runs, defaults to 100,000.
#' @keywords perc.ci.mean.MLEs
#' @export
#' @examples
#'  # perc.ci.mean.MLEs(.95, 10,10^5)

perc.ci.mean.MLEs<-function(cl,n,nr=10^5){
  if(0<=cl & cl<=1){
  al<-(1-cl)/2; cs<-2/sqrt(pi)
  x<-matrix(rMaxwell(nr*n,0,1),nr,n)
  ml<-apply(x, 1, function(x) MLEs(x))
  muh<-ml[1,]; sigh<-ml[2,]
  piv<-(cs-muh)/sigh
  crt<-quantile(piv, c(al,1-al))
  print(crt,3)
  }
  else{print("Invalid confidence level input: Input as a two decimal point")}
}



#' Confidence intervals (MEs based) for two-parameter Maxwell distributions
#'
#' This function computes the required percentiles to compute confidence interval for the mean (moment estimator based) of Maxwell distributions.

#'
#' @param cl confidence level as in decimal, numeric value between 0 and 1.
#' @param n a numeric vector, the number of observation in each sample
#' @param nr numeric vector length one, number of simulation runs, defaults to 100,000.
#' @keywords perc.ci.mean
#' @export
#' @examples
#'  # perc.ci.mean.MEs(.95,5,10^4)


perc.ci.mean.MEs<-function(cl, n,nr=10^5){
  if(0<=cl & cl<=1){
    al<-(1-cl)/2; cs<-2/sqrt(pi)
    x<-matrix(rMaxwell(nr*n,0,1),nr,n)
    ml<-apply(x, 1, function(x) MEs(x))
    muh<-ml[1,]; sigh<-ml[2,]
    piv<-(cs-muh)/sigh
    crt<-quantile(piv, c(al,1-al))
    print(crt,3)
  }
  else{print("Invalid confidence level input: Input as a two decimal point")}
}

#' Prediction interval (PIs) of two-parameter Maxwell distributions
#'
#' This function computes the required percentiles to compute prediction intervals of Maxwell distributions.

#' @param nr numeric vector length one, number of simulation runs, defaults to 100,000.
#' @param n a numeric vector, the number of observation in each sample
#' @param m a numeric vector, the mean of future sample size
#' @param cl confidence level as in decimal, numeric value between 0 and 1.


#' @keywords PI.fac
#' @export
#' @examples
#'  # PI.fac(10^5, n = 45, m = 5, cl = .95)

PI.fac<-function(nr, n, m, cl){
  if(0<=cl & cl<=1){
  x<-matrix(rMaxwell(nr*n,0,1),nr,n)
  ml<-apply(x, 1, function(x) MLEs(x))
  muh<-ml[1,]; sigh<-ml[2,]
  y<-matrix(rMaxwell(nr*m,0,1),nr,m)
  yb<-apply(y, 1, function(x) mean(x))
  piv<-(yb-muh)/sigh
  crt<-unname(quantile(piv,c(.025,.975)))
  print(crt,3)
  }
  else{print("Invalid confidence level input: Input as a two decimal point")}
}

#' One-sided tolerance interval of two-parameter Maxwell distributions
#'
#' This function computes (p, cl) one-sided tolerance factors in order to calculate one sided tolerance intervals.

#' @param nr numeric vector length one, number of simulation runs, defaults to 100,000.
#' @param n a numeric vector, the number of observation in each sample
#' @param p a numeric vector of probabilities
#' @param cl confidence level as in decimal, numeric value between 0 and 1.


#' @keywords fac.os.TLs
#' @export
#' @examples
#' # fac.os.TLs(10^5, n = 10, p = .9, cl = .95)
#'
fac.os.TLs<-function(nr=10^5, n, p, cl){
  if(0<=p & p<=1 & 0<=cl & cl<=1){
  x<-matrix(rMaxwell(nr*n,0,1),nr,n)
  ml<-apply(x, 1, function(x) MLEs(x))
  muh<-ml[1,]; sigh<-ml[2,]
  UppP<-sqrt(qgamma(p,3/2))
  piv<-(UppP-muh)/sigh
  crtU<-quantile(piv, cl)
  LowP<-sqrt(qgamma(1-p,3/2))
  piv<-(LowP-muh)/sigh
  crtL<-quantile(piv, 1-cl)
  print(c(crtL, crtU),3)
  }
  else{print("Invalid input of either probability or confidence level: input must be as a decimal")}
}


#' One-sided tolerance interval (lower and upper) of two-parameter Maxwell distributions
#'
#' This function computes (p, cl) one-sided lower tolerance limit (LTL) and one sided upper tolerance limit (UTL).
#' @param x numeric vector of data input
#' @param nr numeric vector length one, number of simulation runs, defaults to 100,000.
#' @param p a numeric vector of probabilities
#' @param cl confidence level as in decimal, numeric value between 0 and 1.
#' @param tails logical; If TRUE (default), function calculate lower tolerance limit. Otherwise, function calculate upper tolerance limit.


#' @keywords one.sided.MLE.TL
#' @export
#' @examples
#' # x<-c(135, 98, 114, 137, 138, 144, 99, 93, 115, 106, 132, 122, 94, 98, 127,
#' # 122, 102, 133, 114, 120, 93, 126, 119, 104, 119, 114, 125, 107, 98, 117, 111,
#' # 106, 108, 127, 126, 135, 112, 94, 127, 99, 120, 120, 121, 122, 96, 109, 123, 105)
#' # lower tolerance limit
#' # one.sided.MLE.TL(x,10^5,.9,.95,T)
#' # upper tolerance limit
#' # one.sided.MLE.TL(x,10^5,.9,.95,F)
#'

one.sided.MLE.TL<-function(x, nr, p, cl, tails=T){
  if(0<=p & p<=1 & 0<=cl & cl<=1){
    xmle<-MLEs(x)
    n<-length(x)
    y<-matrix(rMaxwell(nr*n,0,1),nr,n)
    ml<-apply(y, 1, function(y) MLEs(y))
    muh<-ml[1,]; sigh<-ml[2,]
    UppP<-sqrt(qgamma(p,3/2))
    piv<-(UppP-muh)/sigh
    crtU<-as.numeric(quantile(piv, cl))
    LowP<-sqrt(qgamma(1-p,3/2))
    piv<-(LowP-muh)/sigh
    crtL<-as.numeric(quantile(piv, 1-cl))
    if(tails== T){
      LTL<-xmle[1]+(crtL*xmle[2])
      return(LTL)}
    if(tails== F){
      UTL<-xmle[1]+(crtU*xmle[2])
      return(UTL)}
     else{print("tails is a logical: either input 'T for LTL' or 'F for UTL' ")}

  }
  else{print("Invalid input of either probability or confidence level: input must be as a decimal")}
}

#' Two-sided or equal-tailed tolerance factor
#'
#' This function computes (p, cl) two-sided or equal-tailed tolerance factor
#' @param nr numeric vector length one, number of simulation runs, defaults to 100,000.
#' @param n numeric vector, number of observation in each sample
#' @param p a numeric vector of probabilities
#' @param cl confidence level as in decimal, numeric value between 0 and 1.
#' @param tails logical; If TRUE, function calculate Equal-tailed tolerance factors. Otherwise, function calculate two-sided tolerance factor.


#' @keywords TFs.Maxwell
#' @export
#' @examples
#' For two-sided tolerence interval
#' # TFs.Maxwell(10^5, n = 20, p = .9, .95, tails = F)
#' For equal-tailed tolerance interval
#' # TFs.Maxwell(10^5, n = 20, p = .9, .95, tails = T)
#'
TFs.Maxwell = function(nr=10^5,n,p,cl,tails){
  if(0<=p & p<=1 & 0<=cl & cl<=1){
  ql = sqrt(qgamma((1-p)/2,1.5))
  qu = sqrt(qgamma((1+p)/2,1.5))
  x = matrix(rMaxwell(nr*n,0,1), nr, n)
  mles = apply(x, 1, function(x) MLEs(x))
  muh = mles[1,]; sigh = mles[2,]
  # one-sided factor
  TF = function(nr, n, muh, sigh, p, cl){
    wu = (qu-muh)/sigh; wl = (ql-muh)/sigh
    UTF = quantile(wu, cl); LTF = quantile(wl, 1-cl)
    return(c(LTF,UTF))
  }

  fn = function(gamt){
    fac = TF(nr,n, muh, sigh, (1+p)/2, (1+gamt)/2)
    Low = muh + fac[1]*sigh; Upp = muh + fac[2]*sigh
    if(tails == T){
      cont = Low <= ql & qu <= Upp}
    else{
      cont = pMaxwell(Upp,0,1)-pMaxwell(Low,0,1)}
    covr = mean(cont >= p)
    return(covr-cl)
  }
  # bisection method
  xl = .5; xr = cl; k = 1
  repeat{
    fl = fn(xl); fr = fn(xr)
    xm = (xl+xr)/2
    fm = fn(xm)
    if(abs(fm) < 1.0e-5 | k > 50){break}
    if(fl*fm > 0){xl = xm}
    else{xr = xm}
    k=k+1
  }
  fac = unname(TF(nr,n, muh, sigh, (1+p)/2, (1+xm)/2))
  print(fac,3)
  }
  else{print("Invalid input of either probability or confidence level: input must be as a decimal")}
}

#' Two-sided or equal-tailed tolerance Interval
#'
#' This function computes (p, cl) two-sided or equal-tailed tolerance Interval
#' @param x a numeric vercor of data input
#' @param nr numeric vector length one, number of simulation runs, defaults to 100,000.
#' @param p a numeric vector of probabilities
#' @param cl confidence level as in decimal, numeric value between 0 and 1.
#' @param tails logical; If TRUE, function calculate Equal-tailed tolerance factors. Otherwise, function calculate two-sided tolerance factor.


#' @keywords TI.MLE
#' @export
#' @examples
#' #' # x<-c(135, 98, 114, 137, 138, 144, 99, 93, 115, 106, 132, 122, 94, 98, 127,
#' # 122, 102, 133, 114, 120, 93, 126, 119, 104, 119, 114, 125, 107, 98, 117, 111,
#' # 106, 108, 127, 126, 135, 112, 94, 127, 99, 120, 120, 121, 122, 96, 109, 123, 105)
#'
#' For equal-tailed tolerance interval:
#' # TI.MLE(x,10^5,.9,.95,T)
#' For two-sided tolerance interval:
#' # TI.MLE(x,10^5,.9,.95,F)
#'
TI.MLE<-function(x,nr, p, cl, tails)
{
  if(0<=p & p<=1 & 0<=cl & cl<=1){
    xmle<-MLEs(x)
    if(tails==T){
      Equal.TF<-TFs.Maxwell(nr,length(x),p,cl,T)
      return(c(xmle[1]+(Equal.TF[1]*xmle[2]), xmle[1]+(Equal.TF[2]*xmle[2])))
    }
    else{
      Two.TF<-TFs.Maxwell(nr,length(x),p,cl,F)
      return(c(xmle[1]+(Two.TF[1]*xmle[2]), xmle[1]+(Two.TF[2]*xmle[2])))
    }
  }
  else{print("Invalid input of either probability or confidence level: input must be as a decimal")}
}

