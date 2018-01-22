#' @title Prediction Accuracy Measures for Regression Models of Right-Censored Data
#'
#' @description This function calculates a pair of measures, R-Squared and L-Squared, for any regression models of right-censored data. R-squared is an extension of the classical R2 statistic for a linear model, quantifying the amount of variability in the response that is explained by a corrected prediction based on the original prediction function. L-squared is the proportion of the prediction error of the original prediction function that is explained by the corrected prediction function, quantifying the distance between the corrected and uncorrected predictions. When used together, they give a complete summary of the predictive power of a prediction function.
#' @export
#' @param y A numeric vector containing the response values.
#' @param y.predict A numeric vector containing the predicted response values from a fitted model.
#' @param delta A numeric vector indicating the status of the event, normally 0=alive, 1=dead.
#' @return  A list containing two components: R-squared and L-squared
#' @examples
#' library(survival)
#' library(PAmeasures)
#'
#'# Use Mayo Clinic Primary Biliary Cirrhosis Data
#'data(pbc)
#'
#'# Fit an exponential model with bilirubin
#'fit1 <- survreg(Surv(time, status==2) ~ bili, data = pbc,dist="exponential" )
#'
#'# Obtain predicted response from the fitted exponential model
#'predict.time<-predict(fit1,type="response")
#'
#'# Recode status at endpoint, 0 for censored, 1 for dead
#'delta.pbc<- as.numeric(pbc$status == 2)
#'
#'# R.squared and L.squared of log-linear model
#'pam.censor(pbc$time, predict.time, delta.pbc)





pam.censor<-function(y,y.predict,delta){



  y.sorted<-sort(y)
  delta.sorted<-delta[order(y)]
  y.predict.sorted<-y.predict[order(y)]

  #KM estimate for censoring distribution
  fit.km.censoring <- survfit(Surv(y.sorted, 1-delta.sorted) ~ 1  )
  #sum.km.censoring<-summary(fit.km.censoring,  censored=TRUE)
  sum.km.censoring<-summary(fit.km.censoring,times=y.sorted,extend=TRUE)
  km.censoring<-sum.km.censoring$surv

  km.censoring.minus<-c(1,km.censoring[-length(km.censoring)])
  ratio.km<- delta.sorted/km.censoring.minus

  ratio.km[is.nan(ratio.km)]<-0

  weight.km<-ratio.km/(sum(ratio.km))






  wls.fitted<-tryCatch(lm(y.sorted~y.predict.sorted,weights=weight.km),error=function(e){return(c(NA,NA))})
  calibrate.fitted<-tryCatch(predict(wls.fitted),error=function(e){return(c(NA,NA))})



  num.rho2<-sum(weight.km*(calibrate.fitted-sum(weight.km*y.sorted))^2)
  denom.rho2<-sum(weight.km*(y.sorted-sum(weight.km*y.sorted))^2)

  R2 <-format(round(num.rho2/denom.rho2,digits = 4) ,nsmall=4)


  num.L2<- sum(weight.km*(y.sorted-calibrate.fitted)^2)
  denom.L2<- sum(weight.km*(y.sorted-y.predict.sorted)^2)
  L2 <-format(round(num.L2/denom.L2,digits = 4),nsmall=4)

  return(list(R.squared=R2,L.squared=L2))
}
