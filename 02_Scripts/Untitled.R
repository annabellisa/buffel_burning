
# make logistic looking plot for fig 1:


x<-1:100
par(mar=c(5,5,1,1))
plot(1:100,dlogis(rev(x), location = 0, scale = 20, log = FALSE), type="l", lwd=5, xlab="", ylab="", xaxt="n", yaxt="n", bty="l", col="cornflowerblue")
box(bty="l", lwd=5)




