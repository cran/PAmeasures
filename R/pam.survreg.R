

#' @title Prediction Accuracy Measures for Parametric Survival Regression Models
#'
#' @description This function calculates a pair of measures, R-Squared and L-Squared, for parametric survival regression models. R-squared is an extension of the classical R2 statistic for a linear model, quantifying the amount of variability in the response that is explained by a corrected prediction based on the original prediction function. L-squared is the proportion of the prediction error of the original prediction function that is explained by the corrected prediction function, quantifying the distance between the corrected and uncorrected predictions. When used together, they give a complete summary of the predictive power of a prediction function.
#' @export
#' @param fit.survreg object inheriting from class survreg representing a fitted parametric survival
#' regression model. Specifying x = TRUE and y=TRUE are required in the call to survreg( )
#'to include the design matrix and the response vector in the object fit.
#' @return  A list containing two components: R-squared and L-squared
#' @examples
#' library(survival)
#' library(PAmeasures)
#'
#'# Use Mayo Clinic Primary Biliary Cirrhosis Data
#'data(pbc)
#'
#' head(pbc)
#'
#'# Fit an exponential model with bilirubin
#'fit1 <- survreg(Surv(time, status==2) ~ bili, data = pbc,dist="exponential",x=TRUE,y=TRUE)
#'
#'# R.squared and L.squared of exponential model
#'pam.survreg(fit1)
#'
#'# Fit a lognormal model with standardised blood clotting time
#'fit2 <- survreg(Surv(time, status==2) ~ protime, data = pbc,dist="lognormal",x=TRUE,y=TRUE)
#'
#'# R.squared and L.squared of lognormal model
#'pam.survreg(fit2)
#'
#'# Fit a weibull model with bilirubin and standardised blood clotting time
#'fit3 <- survreg(Surv(time, status==2) ~ bili + protime, data = pbc,dist="weibull",x=TRUE,y=TRUE)
#'
#'# R.squared and L.squared of weibull model
#'pam.survreg(fit3)






pam.survreg<-function(fit.survreg){


  x.matrix.unsorted<-fit.survreg$x
  if(fit.survreg$dist == "exponential" || fit.survreg$dist == "weibull" || fit.survreg$dist == "lognormal" || fit.survreg$dist == "loglogistic"  ){y.unsorted <-exp(fit.survreg$y[,1])} else{y.unsorted <-fit.survreg$y[,1]}

  censor.unsorted<-fit.survreg$y[,2]


  nsize<-length(y.unsorted)
  y<-sort(y.unsorted)
  delta<-censor.unsorted[order(y.unsorted)]


  p<-dim(as.matrix(x.matrix.unsorted))[2]

  if(p==1){x.matrix<-as.matrix(x.matrix.unsorted[order(y.unsorted)])
  }else{x.matrix<-x.matrix.unsorted[order(y.unsorted),]}

  nsize<-length(y )


  #KM estimate for censoring distribution
  fit.km.censoring <- survfit(Surv(y, 1-delta) ~ 1  )
  #sum.km.censoring<-summary(fit.km.censoring,  censored=TRUE)
  sum.km.censoring<-summary(fit.km.censoring,times=y,extend=TRUE)
  km.censoring<-sum.km.censoring$surv

  km.censoring.minus<-c(1,km.censoring[-length(km.censoring)])
  ratio.km<- delta/km.censoring.minus

  ratio.km[is.nan(ratio.km)]<-0

  weight.km<-ratio.km/(sum(ratio.km))


  if(fit.survreg$dist == "exponential"){t.predicted <-predict(fit.survreg,newdata=data.frame(x.matrix), type="response")*gamma(2)
  }else if(fit.survreg$dist == "weibull"){t.predicted <-predict(fit.survreg,newdata=data.frame(x.matrix), type="response")*gamma(1+fit.survreg$scale)
  }else if(fit.survreg$dist == "lognormal"){t.predicted <-predict(fit.survreg,newdata=data.frame(x.matrix), type="response")*exp((fit.survreg$scale)^2/2)
  }else if(fit.survreg$dist == "loglogistic"){t.predicted <-predict(fit.survreg,newdata=data.frame(x.matrix), type="response")*gamma(1+fit.survreg$scale)*gamma(1-fit.survreg$scale)
  }else{t.predicted <-predict(fit.survreg,newdata=data.frame(x.matrix), type="response")}



  wls.fitted<-tryCatch(lm(y~t.predicted,weights=weight.km),error=function(e){return(c(NA,NA))})
  calibrate.fitted<-tryCatch(predict(wls.fitted),error=function(e){return(c(NA,NA))})



  num.rho2<-sum(weight.km*(calibrate.fitted-sum(weight.km*y))^2)
  denom.rho2<-sum(weight.km*(y-sum(weight.km*y))^2)

  R2 <-format(round(num.rho2/denom.rho2,digits=4) ,nsmall=4)


  num.L2<- sum(weight.km*(y-calibrate.fitted)^2)
  denom.L2<- sum(weight.km*(y-t.predicted)^2)
  L2 <-format(round(num.L2/denom.L2,digits=4),nsmall=4)

  return(list(R.squared=R2,L.squared=L2))
}
