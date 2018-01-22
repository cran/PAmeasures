#' @title Prediction Accuracy Measures for Cox proportional hazards model
#'
#' @description This function calculates a pair of measures, R-Squared and L-Squared, for Cox proportional hazards model. R-squared is an extension of the classical R2 statistic for a linear model, quantifying the amount of variability in the response that is explained by a corrected prediction based on the original prediction function. L-squared is the proportion of the prediction error of the original prediction function that is explained by the corrected prediction function, quantifying the distance between the corrected and uncorrected predictions. When used together, they give a complete summary of the predictive power of a prediction function.
#' @export
#' @param fit.cox object inheriting from class coxph representing a fitted Cox proportional hazards
#'regression model. Specifying x = TRUE and y=TRUE are required in the call to coxph( )
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
#'
#'# Fit a univariate Cox PH model with standardised blood clotting time
#'fit1 <- coxph(Surv(time, status==2) ~ protime, data = pbc,x=TRUE,y=TRUE)
#'
#'# R.squared and L.squared of Cox PH model
#'pam.coxph(fit1)
#'
#'# Fit a multiple Cox PH model with bilirubin and standardised blood clotting time
#'fit2 <- coxph(Surv(time, status==2) ~ bili + protime, data = pbc,x=TRUE,y=TRUE)
#'
#'# R.squared and L.squared of Cox PH model
#' pam.coxph(fit2)


pam.coxph<-function(fit.cox){


x.matrix.unsorted<-fit.cox$x
y.unsorted<-fit.cox$y[,1]
censor.unsorted<-fit.cox$y[,2]

my.beta<-fit.cox$coeff


nsize<-length(y.unsorted)
y<-sort(y.unsorted)
delta<-censor.unsorted[order(y.unsorted)]


p<-dim(as.matrix(x.matrix.unsorted))[2]

if(p==1){x.matrix<-as.matrix(x.matrix.unsorted[order(y.unsorted)])
}else{x.matrix<-x.matrix.unsorted[order(y.unsorted),]}


y.length<-length(y)
yi.matrix<- matrix(rep(y,each=y.length),nrow=y.length)
yj.matrix<- t(yi.matrix)



#KM estimate for censoring distribution
fit.km.censoring <- survfit(Surv(y, 1-delta) ~ 1  )
#sum.km.censoring<-summary(fit.km.censoring,  censored=TRUE)
sum.km.censoring<-summary(fit.km.censoring,times=y,extend=TRUE)
km.censoring<-sum.km.censoring$surv

km.censoring.minus<-c(1,km.censoring[-length(km.censoring)])
ratio.km<- delta/km.censoring.minus

ratio.km[is.nan(ratio.km)]<-0

weight.km<-ratio.km/(sum(ratio.km))




###################################################################
#calculate the predicted survival times based on the cox model
####################################################################
#each column is the risk set for the i-th group
R<-((yj.matrix>=yi.matrix)*1)


#cumulative baseline hazard estimate
my.Lambda<-R%*%(( delta )/t(t(exp(x.matrix%*%my.beta)) %*%R))


rm(R)

my.power<-matrix(rep(t(exp(x.matrix%*%my.beta)),each=y.length),nrow=y.length)



my.factor<-(max(my.Lambda)/100)

my.Lambda2<- t(matrix(rep(exp(-my.Lambda/my.factor),each=y.length),nrow=y.length))



#survival estimate,each column is for the i-th subject, evaluated at all sorted time-points

 S.hat.x<-my.Lambda2^(my.factor*my.power)



rm(my.Lambda)
rm(my.power)
rm(my.Lambda2)


t1<-y
t2<-c(0, t1[1:length(t1)-1])
delta.t<- t1-t2

#predicted survival time from Cox PH model
t.predicted<-colSums(delta.t*S.hat.x)




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
