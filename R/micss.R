##
## Modified Iterative Cummulative Sum of Squares Algorithm Package
##
## J.L. Carrion-i-Silvestre and A. Sanso (2023): Generalized Extreme Value
##    Approximation to the CUMSUMQ Test for Constant Unconditional Variance in Heavy-Tailed Time Series.
##
## 30 / 7 / 2023
##

##
## Test statistic
##

kappa_test <- function(e,sig.lev=0.05,alpha=NULL,kmax=NULL){
  t <- length(e)
  if(t<26){
    kt <- list(kappa=0,tb=0,cv=999)
  } else {
    dk <- numeric(t)
    a <- e-mean(e)
    s2 <- sum(a^2)/t
    for (k in 1:t) {dk[k] <- abs(sum(a[1:k]^2)-k*s2)}
    tb <- which.max(dk)
    a2 <- a^2-s2
    c <- sqrt(lrv.spc.bartlett(a2,kmax))
    kap <- sqrt(1/t)*dk[tb]/c
    if (is.null(alpha)) {alpha=4}
    cv <- cv.kappa(t=t,alpha=alpha,sig.lev=sig.lev)
    p.val <- p.val.kappa(x=kap,t=t,alpha=alpha)
  }
  return(list(kappa=kap,tb=tb,cv=cv,p.val=p.val))
}

#' @title Critical Values
#'
#' @description Computes critical values for the \link{kappa_test} statistic
#'
#' @param t Sample size.
#' @param alpha Value of the tail index between 2 and 4.
#' @param sig.lev Significance level.
#' @return Critical value.
#'
#' @seealso \link{p.val.kappa}
#'
#' @references
#' J.L. Carrion-i-Silvestre & A. Sansó (2023): Generalized Extreme Value Approximation to the CUMSUMQ Test for Constant Unconditional Variance in Heavy-Tailed Time Series.
#'
#' @keywords internal

cv.kappa <- function(t,alpha,sig.lev){
  loc <- sr_loc(t,alpha)
  sc <- sr_scale(t,alpha)
  sh <- sr_shape(t,alpha)
  p <- 1-sig.lev
  if (sh==0) cv <- loc-sc*log(-log(p))
  else cv <- loc+sc/sh*((-log(p))^(-sh)-1)
  return(cv)
}

#' @title P-values
#'
#' @description Computes p-values for the \link{kappa_test} statistic
#'
#' @param x Value of the statistic.
#' @param t Sample size.
#' @param alpha Value of the tail index between 2 and 4.
#' @return p-value.
#'
#' @seealso \link{cv.kappa}
#'
#' @references
#' J.L. Carrion-i-Silvestre & A. Sansó (2023): Generalized Extreme Value Approximation to the CUMSUMQ Test for Constant Unconditional Variance in Heavy-Tailed Time Series.
#'
#' @keywords internal
p.val.kappa <- function(x,t,alpha){
  loc <- sr_loc(t,alpha)
  sc <- sr_scale(t,alpha)
  sh <- sr_shape(t,alpha)
  if (sh==0) p.val <- exp(-exp(-(x-loc)/sc))
  else p.val <- 1-exp(-(1+sh*((x-loc)/sc))^(-1/sh))
  return(p.val)
}

#' @title sr_loc
#'
#' @description Computes the location parameter for the GEV approximation
#' to the distribution of \link{kappa_test} statistic
#'
#' @param t Sample size
#' @param alpha Value of the tail index between 2 and 4
#' @param lmax Maximum lag to be used for the long-run estimation of the fourth order moment of the
#' innovations. If not specified it is generated automatically depending on the sample size.
#' @return Location parameter
#'
#' @details used internally by \link{cv.kappa} and \link{p.val.kappa}
#'
#' @references
#' J.L. Carrion-i-Silvestre & A. Sansó (2023): Generalized Extreme Value Approximation to the CUMSUMQ Test for Constant Unconditional Variance in Heavy-Tailed Time Series.
#'
#' @keywords internal
sr_loc <- function(t,alpha,lmax=NULL){
  if (is.null(lmax)) lmax <- floor(12*(t/100)^(1/4))
  b <- c(7.428131e-01,-5.386952e-02,1.906445e-02,-1.358935e-03,
         3.594350,-1.216379e+02,1.258597e+03,-3.416918e-01,
         6.286310e-01,-3.392324e-01,3.162517e+01)
  reg <- c(1,alpha,alpha^2,alpha^3,1/t,1/t^2,1/t^3,1/lmax,
           alpha/t,alpha^2/t,alpha/t^2)
  return(sum(b*reg))
}

#' @title sr_scale
#'
#' @description Computes the scale parameter for the GEV approximation
#' to the distribution of \link{kappa_test} statistic
#'
#' @param t Sample size
#' @param alpha Value of the tail index between 2 and 4
#' @param lmax Maximum lag to be used for the long-run estimation of the fourth order moment of the
#' innovations. If not specified it is generated automatically depending on the sample size.
#' @return Scale parameter
#'
#' @details Used internally by \link{cv.kappa} and \link{p.val.kappa}
#'
#' @references
#' J.L. Carrion-i-Silvestre & A. Sansó (2023): Generalized Extreme Value Approximation to the CUMSUMQ Test for Constant Unconditional Variance in Heavy-Tailed Time Series.
#'
#' @keywords internal
sr_scale <- function(t,alpha,lmax=NULL){
  if (is.null(lmax)) lmax <- floor(12*(t/100)^(1/4))
  b <- c(8.508460e-02,1.018279e-01,-3.171846e-02,3.800244e-03,
         1.684208e+01,-4.080742e+02,3.683425e+03,-5.104544e-01,
         -3.426550,7.149008e+01)
  reg <- c(1,alpha,alpha^2,alpha^3,1/t,1/t^2,1/t^3,1/lmax,
           alpha/t,alpha/t^2)
  return(sum(b*reg))
}

#' @title sr_shape
#'
#' @description Computes the scale parameter for the GEV approximation
#' to the distribution of \link{kappa_test} statistic
#'
#' @param t Sample size
#' @param alpha Value of the tail index between 2 and 4
#' @param lmax Maximum lag to be used for the long-run estimation of the fourth order moment of the
#' innovations. If not specified it is generated automatically depending on the sample size.
#' @return Shape parameter
#'
#' @details Used internally by \link{cv.kappa} and \link{p.val.kappa}
#'
#' @references
#' J.L. Carrion-i-Silvestre & A. Sansó (2023): Generalized Extreme Value Approximation to the CUMSUMQ Test for Constant Unconditional Variance in Heavy-Tailed Time Series.
#'
#' @keywords internal
sr_shape <- function(t,alpha,lmax=NULL){
  if (is.null(lmax)) lmax <- floor(12*(t/100)^(1/4))
  b <- c(8.640358e-01,-7.356562e-01,2.516790e-01,-3.017798e-02,
         3.076360e+01,-1.264638e+03,1.620025e+04,-9.650881e-01,
         -4.956778,2.510715e+02,-2.373360e+03)
  reg <- c(1,alpha,alpha^2,alpha^3,1/t,1/t^2,1/t^3,1/lmax,
           alpha/t,alpha/t^2,alpha/t^3)
  return(sum(b*reg))
}

##
## MICSS algorithm
##
micss <- function(e,sig.lev=0.05,kmax=NULL,alpha=NULL,
                  tail.est="NR",k=0.1){
  if (is.null(alpha)) { # estimates tail index with absolut values
    alpha.est <- estimate.alpha(x=abs(e),sig.lev=sig.lev,tail.est=tail.est,k=k)
  } else {alpha.est <- list(alpha=alpha, sd.alpha=0, tail.est="fixed")}
  if (alpha.est$alpha<=2) { # No finite variance
    icss.results <- list(nb=0,tb=0,e=e,v=Inf)
  } else {
    icss.results <- icss(e,sig.lev=sig.lev,kmax=kmax,
                         alpha=alpha.est$alpha)
  }
  bbm <- list(icss=icss.results,alpha=alpha.est)
  class(bbm) <- "micss"
  return(bbm)
}

#' @title print.micss
#'
#' @description Prints the output of \link{micss}.
#'
#' @param x An object with the output of the \link{micss} algorithm.
#' @param ... Further arguments passed to or from other methods.
#' @return{No return value. It prints the output of \link{micss}}
#'
#' @examples
#' set.seed(2)
#' e <- c(stats::rnorm(200),3*stats::rnorm(200))
#' o <- micss(e)
#' print.micss(o)
print.micss <- function(x,...){
  if (!is(x,"micss")) {
    stop("Wrong object. Enter the output of micss()")
  }
  aa <- x$alpha
  cat("Tail estimator: ", aa$tail.est,".\n", sep="")
  cat("Estimated alpha = ", aa$alpha.fit, "; used alpha = ",
      aa$alpha, "; std error = ", aa$sd.alpha,
      "\n", sep="")
  cat("Test Ho: (alpha>=4) = ", aa$t4, "; p-value = ",
      stats::pnorm(aa$t4),"\n", sep="")
  cat("Test Ho: (alpha<=2) = ", aa$t2, "; p-value = ",
      1-stats::pnorm(aa$t2),"\n\n", sep="")
  bb <- x$icss
  nb <- bb$nb
  if (length(bb$v)==1){
    if (bb$v==Inf) cat("WARNING: No finite variance\n\n")
  } else {
    cat("Number of breaks: ", nb, " (", nb+1," different periods)\n", sep="")
    if (nb==0)  {
      cat("\nMax value of kappa: ", bb$kappa.max,
          "; p-value = ", bb$p.val,"\n", sep="")
    } else print(taula.micss(bb,alpha=aa$alpha))
  }
}

#' @title taula.micss
#'
#' @description When there are breaks in the unconditional variance, it
#' creates a matrix with the different periods, estimated variances, values of
#'  \link{kappa_test} and \link{p.val.kappa}.
#'
#' @param x An object with the output of \link{icss}.
#' @param alpha Value of the tail index between 2 and 4
#' @return Matrix with the different periods, estimated variances, values of
#'  \link{kappa_test} and \link{p.val.kappa}.
#'
#' @details Used internally by \link{print.micss}.
#'
#' @keywords internal
taula.micss <- function(x,alpha){
  e <- x$e
  t <- length(e)
  nb <- x$nb
  if (nb>0) {
    ta <- taula.icss(x)
    kappas <- rep(NA,nb+1)
    pval <- rep(NA,nb+1)
    for (i in 1:nb){
      kap <- kappa_test(e[ta[i,1]:ta[i+1,2]],alpha=alpha)
      kappas[i] <- kap$kappa
      pval[i] <- kap$p.val
    }
    taula <- cbind(ta,kappas,pval)
    colnames(taula) <- c("Start","End","Var", "Kappa", "p-value")
  } else {taula <- NA}
  return(taula)
}

##
## ICSS algorithm
##
icss <- function(e,sig.lev=0.05,kmax=NULL,alpha=NULL){
  bb <- steps1.to.2b(e=e,sig.lev=sig.lev,kmax=kmax,alpha=alpha)
  if (bb$nb>1){
    bb <- step2c(e=e,bb=bb,sig.lev=sig.lev,kmax=kmax,alpha=alpha)
    bb <- step3(e=e,bb=bb,sig.lev=sig.lev,kmax=kmax,alpha=alpha)
  }
  if (bb$nb==0) bb <- list(nb=bb$nb,tb=bb$tb,e=e,
                           kappa.max=bb$kappa.max,
                           p.val=bb$p.val)
  else bb <- list(nb=bb$nb,tb=bb$tb,e=e)
  bb$v <- var.est.icss(bb)
  class(bb)<-"icss"
  return(bb)
}

#' @title steps1.to.2b
#'
#' @description Computes steps 1 to 2b of the \link{icss} Algorithm.
#'
#' @param e Stationary variable on which the constancy of unconditional variance is tested.
#' @param sig.lev Significance level.
#' @param kmax Maximum lag to be used for the long-run estimation of the fourth order moment of the innovations.
#' @param alpha Tail index.
#' @return
#' \itemize{
#'   \item \code{nb}: Number of breaks.
#'   \item \code{tb}: Time of the breaks.
#'   \item \code{kappa.max}: Maximum value of the \link{kappa_test} if there is no break.
#'   \item \code{p.val}: p-value.
#' }
#'
#' @details Used internally by \link{icss}.
#'
#' @references
#' C. Inclan & G.C. Tiao (1994): Use of Cumulative Sums of Squares for Retrospective Detection of Changes of Variance. Journal of the American Statistical Association 89, 913-923.
#'
#' @keywords internal
steps1.to.2b <- function(e,sig.lev,kmax,alpha){
  if(length(e)<26){
    nb <- 0
    tb <- 0
  } else {
    kt <- kappa_test(e=e,sig.lev=sig.lev,kmax=kmax,alpha=alpha)  # step1
    if (kt$kappa > kt$cv) {         # there is a break
      t2 <- kt$tb
      tb.f <- step2a(e[1:t2],sig.lev=sig.lev,kmax=kmax,alpha=alpha)
      tb.l <- step2b(e[(t2+1):length(e)],sig.lev=sig.lev,kmax=kmax,
                     alpha=alpha) + t2-1
      if (is.null(kmax)) {kmax <- floor(12*(length(e)/100)^(1/4))}
      if ((tb.l-tb.f)<(kmax+3)){    # there is only one break
        nb <- 1
        tb <- tb.f
      } else {             # more than one break
        nb <- 2
        tb <- c(tb.f,tb.l)
      }
    } else {              # No breaks found
      nb <- 0
      tb <- 0
      kappa.max <-kt$kappa
    }
  }
  if (nb==0) return(list(nb=nb,tb=tb,kappa.max=kappa.max,p.val=kt$p.val))
  else return(list(nb=nb,tb=tb))
}

#' @title step2a
#'
#' @description Computes step 2a of the \link{icss} Algorithm.
#'
#' @param e Stationary variable on which the constancy of unconditional variance is tested.
#' @param sig.lev Significance level.
#' @param kmax Maximum lag to be used for the long-run estimation of the fourth order moment of the innovations.
#' @param alpha Tail index.
#' @return
#' \itemize{
#'   \item \code{tb}: Time of the break.
#' }
#'
#' @details Used internally by \link{icss}.
#'
#' @references
#' C. Inclan & G.C. Tiao (1994): Use of Cumulative Sums of Squares for Retrospective Detection of Changes of Variance. Journal of the American Statistical Association 89, 913-923.
#'
#' @keywords internal
step2a <- function(e,sig.lev,kmax,alpha){
  tb <- length(e)
  a <- e
  repeat{
    kt <- kappa_test(a,sig.lev=sig.lev,kmax=kmax,alpha=alpha)
    if (kt$kappa > kt$cv) {  # there is a break
      tb <- kt$tb
      a <- a[1:tb]
    } else break
  }
  return(tb)
}

#' @title step2b
#'
#' @description Computes step 2b of the \link{icss} Algorithm.
#'
#' @param e Stationary variable on which the constancy of unconditional variance is tested.
#' @param sig.lev Significance level.
#' @param kmax Maximum lag to be used for the long-run estimation of the fourth order moment of the innovations.
#' @param alpha Tail index.
#' @return
#' \itemize{
#'   \item \code{tb}: Time of the break.
#' }
#'
#' @details Used internally by \link{icss}.
#'
#' @references
#' C. Inclan & G.C. Tiao (1994): Use of Cumulative Sums of Squares for Retrospective Detection of Changes of Variance. Journal of the American Statistical Association 89, 913-923.
#'
#' @keywords internal
step2b <- function(e,sig.lev,kmax,alpha){
  a <- e
  t <- length(e)
  repeat{
    kt <- kappa_test(a,sig.lev=sig.lev,kmax=kmax,alpha=alpha)
    if (kt$kappa > kt$cv) {  # there is a break
      tb <- kt$tb
      a <- a[(tb+1):length(a)]
    } else break
  }
  return(tb=t-length(a)+1) # OK
}

#' @title step2c
#'
#' @description Computes step 2c of the \link{icss} Algorithm.
#'
#' @param e Stationary variable on which the constancy of unconditional variance is tested.
#' @param bb The output of function \link{steps1.to.2b}.
#' @param sig.lev Significance level.
#' @param kmax Maximum lag to be used for the long-run estimation of the fourth order moment of the innovations.
#' @param alpha Tail index.
#' @return
#' \itemize{
#'   \item \code{nb}: Number of breaks.
#'   \item \code{tb}: Time of the breaks.
#' }
#'
#' @details Used internally by \link{icss}.
#'
#' @references
#' C. Inclan & G.C. Tiao (1994): Use of Cumulative Sums of Squares for Retrospective Detection of Changes of Variance. Journal of the American Statistical Association 89, 913-923.
#'
#' @keywords internal
step2c <- function(e,bb,sig.lev,kmax,alpha){
  nb <- bb$nb
  tb <- bb$tb
  i <- 1      # iteration
  a <- e[(tb[i]+1):tb[i+1]]
  repeat{
    if (length(a)<26) {break} # minimal sample
    bb <- steps1.to.2b(a,sig.lev=sig.lev,kmax=kmax,alpha=alpha)
    if (bb$nb==0) {       # no new breaks
      break
    } else if (bb$nb==1){ # only one new break
      tb <- c(tb,bb$tb+tb[i])
      tb <- sort(tb)
      nb <- nb+1
      break
    } else { # two new breaks
      tb <- c(tb,(bb$tb+tb[i]))
      tb <- sort(tb)
      i <- i+1
      a <- e[(tb[i]+1):tb[i+1]]
      nb <- nb+2
    }
  }
  return(list(nb=nb,tb=tb))
}

#' @title step3
#'
#' @description Computes step 3 of the \link{icss} Algorithm.
#'
#' @param e Stationary variable on which the constancy of unconditional variance is tested.
#' @param bb The output of \link{step2c}, A list with the number of break and the time of the breaks.
#' @param sig.lev Significance level.
#' @param kmax Maximum lag to be used for the long-run estimation of the fourth order moment of the innovations.
#' @param alpha Tail index.
#' @return
#' \itemize{
#'   \item \code{nb}: Number of breaks.
#'   \item \code{tb}: Time of the breaks.
#' }
#'
#' @details Used internally by \link{icss}.
#'
#' @references
#' C. Inclan & G.C. Tiao (1994): Use of Cumulative Sums of Squares for Retrospective Detection of Changes of Variance. Journal of the American Statistical Association 89, 913-923.
#'
#' @keywords internal
step3 <- function(e,bb,sig.lev,kmax,alpha){
  b0 <- c(0,bb$tb,length(e))
  iter <- 1
  max.iter <- 10
  repeat{
    nb <- (length(b0)-2)
    if (nb>1){
      if (iter == max.iter) {break}
      b1 <- 0
      for (i in 1:nb){   # checking the break points
        a <- e[(b0[i]+1):b0[i+2]]
        kt <- kappa_test(a,sig.lev=sig.lev,kmax=kmax,alpha=alpha)
        if (kt$kappa > kt$cv) {b1 <- c(b1,(b0[i]+kt$tb))}
      }
      b1 <- c(b1,length(e))
      b1 <- unique(b1)
      if ((length(b1)==length(b0))){
        if ((max(abs(b1-b0))<=2)) {break}
        else {
          b0 <- b1
          iter <- iter+1
        }
      } else {
        b0 <- b1
        iter <- iter+1
      }
    } else {
      b1 <- b0
      break
    }
  }
  nb <- (length(b1)-2)
  if (nb>0) {tb <- b1[2:(length(b1)-1)]}
  else tb <- 0
  return(list(nb=nb,tb=tb))
}

#' @title plot.icss
#'
#' @description Plots the output of the ICSS algorithm.
#'
#' @param x An object with the output of \link{icss} or \link{micss}.
#' @param type Type of graphic. 3 possibilities: "std", which is the default,
#' plots the absolute value of the variable and the standard deviation; "var"
#' plots the squares of the variable and the variance; "res.std" plots
#' the standardized residuals.
#' @param title Title of the graphic.
#' @param ... Further arguments passed to or from other methods.
#' @return{No return value. It generates a plot the output of \link{micss} or \link{icss}}
#'
#' @examples
#' set.seed(2)
#' e <- c(stats::rnorm(200),3*stats::rnorm(200))
#' o <- micss(e)
#' plot.icss(o,title="Example of the MICSS algorithm")
plot.icss <- function(x,type="std",title=NULL,...){
  if (is(x,"icss")) {kk <- x}
  else if (is(x,"micss")) {kk <- x$icss}
  else {
    stop("Wrong object. Enter the output of micss() or icss()")
  }
  e <- kk$e
  t <- length(e)
  if (is.null(kk$v)) v <- var.est.icss(kk)
  else v <- kk$v
  if (type=="std"){
    plot(x=1:t, y=abs(e), type = "l", xlab="t",
         ylab="Absolute residuals and std.dev", main=title)
    graphics::lines(sqrt(v),col="red")
  } else if (type=="var"){
    plot(x=1:t, y=e^2, type = "l", xlab="t",
         ylab="Squared residuals", main=title)
    graphics::lines(v,col="red")
  } else if (type=="res.std"){
    plot(x=1:t, y=e/sqrt(v), type = "l", xlab="t",
         ylab="Standarized residuals", main=title)
  }
}


#' @title taula.icss
#'
#' @description When there are breaks in the unconditional variance, it
#' creates a matrix with the different periods and estimated variances.
#'
#' @param x An object with the output of \link{icss}.
#' @return Matrix with the different periods and estimated variances.
#'
#' @details Used internally by \link{print.icss}. If there are no break returns NA value.
#'
#' @keywords internal
taula.icss <- function(x){
  e <- x$e
  t <- length(e)
  nb <- x$nb
  if (nb>0) {
    taula <- matrix(NA, nrow=(nb+1), ncol=3)
    colnames(taula) <- c("Start","End","Var")
    rownames(taula) <- seq(1,(nb+1),1)
    b <- c(0,x$tb,t)
    for (i in 1:(nb+1)){
      taula[i,1] <- b[i]+1
      taula[i,2] <- b[i+1]
      taula[i,3] <- stats::var(e[(b[i]+1):b[i+1]])
    }
  } else {taula <- NA}
  return(taula)
}

#' @title print.icss
#'
#' @description Prints the output of \link{icss}.
#'
#' @param x An object with the output of the \link{icss} algorithm.
#' @param ... Further arguments passed to or from other methods.
#' @return{No return value. It prints the output of \link{icss}}
#'
#' @details Used internally by \link{icss}.
#'
#' @examples
#' set.seed(2)
#' e <- c(stats::rnorm(200),3*stats::rnorm(200))
#' o <- icss(e)
#' print.icss(o)
print.icss <- function(x,...){
  if (!is(x,"icss")) {
    stop("Wrong object. Enter the output of icss()")
  }
  nb <- x$nb
  cat("Number of breaks: ", nb, " (", nb+1," different periods)\n", sep="")
  if (nb>0) print(taula.icss(x))
}

#' @title var.est.icss
#'
#' @description Computes the variance of each period according to the breaks found with \link{icss} or \link{micss}.
#'
#' @param x An object with the output of \link{icss} or \link{micss}.
#' @return A vector with the variances.
#'
#' @keywords internal
var.est.icss <- function(x){
  e <- x$e
  t <- length(e)
  if (x$nb==0) {v <- rep(stats::var(e),t)}
  else {
    taula <- taula.icss(x)
    v <- rep(taula[1,3],taula[1,2])
    for (i in 2:nrow(taula)){
      v <- c(v,rep(taula[i,3],(taula[i,2]-taula[i,1]+1)))
    }
  }
  return(v)
}


##
## Procedures to estimate and test the index of tail thickness
##

#' @title estimate.alpha
#'
#' @description Computes the estimator of the tail index (alpha) using Hill (1975) or Nicolau & Rodrigues (2019) estimators and
#' tests both the hypothesis that alpha is bigger or equal to 4, and that alpha is lower or equal to 2.
#'
#' @param x A numeric vector.
#' @param sig.lev Significance level. The default value is 0.05.
#' @param tail.est Estimator of the tail index. The default value is "Hill", which uses Hill's (1975) estimator. "NR" uses Nicolau & Rodrigues (2019) estimator.
#' @param k Fraction of the upper tail to be used to estimate of the tail index.
#' @return
#' \itemize{
#'   \item \code{alpha}: Value of the tail index (alpha) to be used in \link{micss}.
#'   \item \code{alpha.fit}: Estimated tail index.
#'   \item \code{t4}: Test of the null hypothesis alpha>=4 against alpha<4.
#'   \item \code{t2}: Test of the null hypothesis alpha<=2 against alpha>2.
#' }
#'
#' @details Used internally by \link{micss}.
#'
#' @references
#' B. Hill (1975): A Simple General Approach to Inference About the Tail of a Distribution. The Annals of Mathematical Statistics 3, 1163-1174.
#'
#' J. Nicolau and P.M.M. Rodrigues (2019): A new regression-based tail index estimator. The Review of Economics and Statistics 101, 667-680.
#'
#' @keywords internal
estimate.alpha<-function(x,sig.lev=0.05,tail.est="Hill",k=0.1){
  if (tail.est=="NR") {
    alpha.fit <- alpha_nr(x,k)
    t4 <- ((alpha.fit$alpha-4)/alpha.fit$sd.alpha) # test Ho: alpha=4, Ha: alpha<4
    t2 <- ((alpha.fit$alpha-2)/alpha.fit$sd.alpha) # test Ho: alpha=2, Ha: alpha>2
    if (t4>stats::qnorm(sig.lev)) {alpha <- 4}
    else {
      if (t2<stats::qnorm(1-sig.lev)) {alpha=2}
      else {alpha=alpha.fit$alpha}
    }
  }
  else {         # "Hill"
    alpha.fit <- alpha_hill(x,k)
    t4 <- sqrt(alpha.fit$s)*(alpha.fit$alpha/4-1) # test Ho: alpha=4, Ha: alpha<4
    t2 <- sqrt(alpha.fit$s)*(alpha.fit$alpha/2-1) # test Ho: alpha=2, Ha: alpha>2
    if (t4>stats::qnorm(sig.lev)) {alpha <- 4}
    else {
      if (t2<stats::qnorm(1-sig.lev)) {alpha=2}
      else {alpha=alpha.fit$alpha}
    }
  }
  return(list(alpha=alpha, alpha.fit=alpha.fit$alpha, t4=t4, t2=t2,
              sd.alpha=alpha.fit$sd.alpha, tail.est=tail.est))
}

alpha_hill <- function(x,k){
  y <- stats::na.omit(x)
  x0 <- stats::quantile(y,1-k, na.rm = TRUE)
  s <- length(y[(y>x0)])
  alpha <- s/sum(log(y[(y>x0)]/x0))
  sd.alpha <- sqrt(alpha^2/s)
  return(list(alpha=alpha,sd.alpha=sd.alpha,s=s))
}

alpha_nr <- function(y,k){
  y <- stats::na.omit(y)
  t <- length(y)
  ah <- alpha_hill(y,k)$alpha
  m <- floor(t*k)+10
  u <- seq(1/m,(m-1)/m,1/m)
  x0 <- stats::quantile(y,1-k, na.rm = TRUE)
  x <- (1-u)^(-1/ah)*x0
  x <- x[(x<=max(y))]
  m <- length(x)
  pr <- rep(0,m)
  for(i in 1:m){
    pr[i] <- length(y[(y>x[i])])/t
  }
  yy <- log(pr)
  z <- -log(x)
  lm.fit <- stats::lm(yy~z)
  b <- stats::coef(lm.fit)
  sd.alpha <- sqrt(2*b[2]^2/m)
  return(list(alpha=b[2],sd.alpha=sd.alpha))
}


##
## Whitening and long-run variance
##

whitening <- function(y,kmax=NULL){
  t <- length(y)
  min_BIC <- abs(log(stats::var(y)))*1000
  res <- e <- y-mean(y)
  rho <- 0
  k <- 0
  if (is.null(kmax)) {kmax <- floor(12*(t/100)^(1/4))}
  if (kmax==0) {kmax <- 1}
  for (i in 1:kmax){
    temp <- e
    for (j in 1:i) {temp <- cbind(temp, dplyr::lag(e,j))}
    temp <- temp[(i+1):nrow(temp),]
    X <- as.matrix(temp[,(2:ncol(temp))])
    y <- as.matrix(temp[,1])
    lm.fit <- stats::lm(y~0+X, na.action=na.omit)
    rho_temp <- stats::coef(lm.fit)
    res_temp <- stats::resid(lm.fit)
    BIC <- log(sum(res_temp^2)/(t-kmax))+(i*log(t-kmax)/(t-kmax))
    if (BIC < min_BIC){
      min_BIC <- BIC
      k <- i
      rho <- rho_temp
      res <- res_temp  # Whited
    }
  }
  return(list(e=res,rho=rho,lag=k))
}

#' @title lrv.spc.bartlett
#'
#' @description Estimation of the long-run variance using the Barlett window.
#'
#' @param x Stationary variable. A numeric vector.
#' @param kmax Maximum lag to be used for the long-run estimation of the variance.
#' @return Estimation of the long-run variance.
#'
#' @details Estimates the log-run fourth order moment when x are the squares of a variable.
#'
#' @references
#' D. Sul, P.C.B. Phillips & C.Y. Choi (2005): Prewhitening Bias in HAC Estimation, Oxford Bulletin of Economics and Statistics 67, 517-546.
#'
#' D.W.K. Andrews & J.C. Monahan (1992): An Improved Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimator. Econometrica 60, 953-966.
#'
#' @examples
#' lrv.spc.bartlett(rnorm(100))
lrv.spc.bartlett <- function(x,kmax=NULL){
  t <- length(x)
  w <- whitening(x,kmax=kmax)
  res <- w$e
  rho <- w$rho
  k <- w$lag
  temp <- cbind(res,dplyr::lag(res,1))[2:length(res),]
  X <- as.matrix(temp[,2])
  y <- as.matrix(temp[,1])
  a <- coef(stats::lm(y~0+X)) # AR(1) approximation as in Andrews & Monahan (1992, pag. 958)
  l <- 1.1447*(4*a^2*t/((1+a)^2*(1-a)^2))^(1/3) # Obtaining the bandwidth for the spectral window
  l <- min(trunc(l),(length(res)-1)) # Truncate the estimated bandwidth
  lrv <- sum(res^2)/t # Short-run variance
  if (is.na(l)) {l <- 0}
  if (l>0){
    for (i in 1:l){     # Long-run variance
      w <- (1-i/(l+1))  # Bartlett kernel
      lrv <- lrv+2*w*sum(res[1:(length(res)-i)]*res[(1+i):length(res)])/t
    }
  }
  lrv_recolored <- lrv/(1-sum(rho))^2
  lrv <- min(lrv_recolored,(t*lrv)) # Sul, Phillips & Choi (2003) boundary rule
  return(lrv)
}
