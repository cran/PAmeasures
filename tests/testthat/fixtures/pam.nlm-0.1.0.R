# PAmeasures 0.1.0 pam.nlm reference, Xiaoyan Wang and Gang Li, GPL-3.
# Source: CRAN Archive/PAmeasures/PAmeasures_0.1.0.tar.gz.
# Function body preserved; line endings normalized.
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
