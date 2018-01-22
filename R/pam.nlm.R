#' @title Prediction Accuracy Measures for Nonlinear Regression Models.
#'
#' @description This function calculates a pair of measures, R-Squared and L-Squared, for any nonlinear regression model. R-squared is an extension of the classical R2 statistic for a linear model, quantifying the amount of variability in the response that is explained by a corrected prediction based on the original prediction function. L-squared is the proportion of the prediction error of the original prediction function that is explained by the corrected prediction function, quantifying the distance between the corrected and uncorrected predictions. When used together, they give a complete summary of the predictive power of a prediction function.
#' @export
#' @param y A numeric vector containing the response values.
#' @param y.predict A numeric vector containing the predicted response values from a fitted model.
#' @return  A list containing two components: R-squared and L-squared
#' @examples
#' library(PAmeasures)
#'
#' data(moore)
#'
#' head(moore)
#'
#' # Transistor count
#' count <- moore$count
#'
#' time<-moore$time
#'
#' # Fit a log-linear model
#' moore.glm= glm(log2(count) ~ time, family=gaussian(link = "identity") )
#'
#' # Obtain predicted transistor count
#' count.predict<-2^(predict(moore.glm,newdata = data.frame(X = time),type = "response" ))
#'
#' # R.squared and L.squared of log-linear model
#' pam.nlm(count, count.predict)

pam.nlm<-function(y,y.predict){

nsize<-length(y)


cal.lm<-lm(y~y.predict)

cal.predict<-predict(cal.lm)

R2<-(sum((cal.predict-sum(y)/nsize)^2))/(sum((y-sum(y)/nsize)^2))

R2 <-format(round(R2,digits=4) ,nsmall=4)

L2<-(sum((y-cal.predict)^2))/(sum((y-y.predict)^2))
L2 <-format(round(L2,digits=4),nsmall=4)

return(list(R.squared=R2,L.squared=L2))
}
