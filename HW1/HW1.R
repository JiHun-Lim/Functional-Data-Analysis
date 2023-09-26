library(fds)
library(fda)
library(ggplot2)
library(fields)
library(expm)
# rm(list=ls())

time = pinchtime
basis <- create.bspline.basis(c(0,0.3),nbasis=15, norder=4)
plot(basis)
pinch.F<-Data2fd(time, pinch, basis)
plot(pinch.F)

plot(pinch.F, col = "grey")
plot(mean.fd(pinch.F), lwd =4, add = TRUE)
plot(std.fd(pinch.F), lwd =4, add = TRUE)

pinch_var<-var.fd(pinch.F)
pts<-seq(from=0, to=0.3, length = 50)
pinch_mat = eval.bifd(pts, pts, pinch_var)
persp(pts, pts, pinch_mat)
contour(pts, pts, pinch_mat)

pinch.pca = pca.fd(pinch.F, nharm=4)
pinch.pca$varprop

##########################################################################################

yield = FedYieldcurve
terms = yield$x
plot(terms, yield$y[,1], pch=15, ylab="Yield", ylim=c(0,16))
points(terms, yield$y[,330], pch=16)


yield_data = yield$y
basis <- create.bspline.basis(c(0,330), nbasis=4)
yield.F <- Data2fd(terms, yield_data, basis)
plot(yield.F)

plot(yield.F, col = "grey")
plot(mean.fd(yield.F), lwd =4, add = TRUE)

yield.pca = pca.fd(yield.F, nharm=1)
yield.pca$varprop
plot(yield.pca$harmonics)

yield = FedYieldcurve
terms = yield$x
yield_data = yield$y[,1]
basis <- create.bspline.basis(c(0,330), nbasis=4)
yield_smooth <- smooth.basis(terms, yield_data, basis)
plot(terms, yield$y[,1], pch=15, ylab="Yield", ylim=c(0,16))
plot(yield_smooth, add = TRUE)



yield = FedYieldcurve
terms = yield$x
yield_data = yield$y[,1]
basis <- create.bspline.basis(c(0,330), nbasis=4)
yield_smooth <- smooth.basis(terms, yield_data, basis)
plot(terms, yield$y[,1], pch=15, ylab="Yield", ylim=c(0,16))
plot(yield_smooth, add = TRUE)
basis <- create.bspline.basis(c(0,330), nbasis=6)
yield_par <- fdPar(basis, Lfdobj =2, lambda = 1)
yield_smooth <- smooth.basis(terms, yield_data, yield_par)
plot(yield_smooth, col = "red", add = TRUE)

yield = FedYieldcurve
terms = yield$x
yield_data = yield$y[,1]
basis <- create.bspline.basis(c(0,330), nbasis=4)
yield_smooth <- smooth.basis(terms, yield_data, basis)
plot(terms, yield$y[,1], pch=15, ylab="Yield", ylim=c(0,16))
plot(yield_smooth, add = TRUE)
basis <- create.bspline.basis(c(0,330), nbasis=6)
yield_par <- fdPar(basis, Lfdobj =2, lambda = 0.1)
yield_smooth <- smooth.basis(terms, yield_data, yield_par)
yield_smooth$gcv
plot(yield_smooth, col = "red", add = TRUE)
basis <- create.bspline.basis(c(0,330), nbasis=6)
yield_par <- fdPar(basis, Lfdobj =2, lambda = 0.2)
yield_smooth <- smooth.basis(terms, yield_data, yield_par)
yield_smooth$gcv
plot(yield_smooth, col = "blue", add = TRUE)
basis <- create.bspline.basis(c(0,330), nbasis=6)
yield_par <- fdPar(basis, Lfdobj =2, lambda = 0.3)
yield_smooth <- smooth.basis(terms, yield_data, yield_par)
yield_smooth$gcv
plot(yield_smooth, col = "green", add = TRUE)
basis <- create.bspline.basis(c(0,330), nbasis=6)
yield_par <- fdPar(basis, Lfdobj =2, lambda = 0.5)
yield_smooth <- smooth.basis(terms, yield_data, yield_par)
yield_smooth$gcv
plot(yield_smooth, col = "purple", add = TRUE)
basis <- create.bspline.basis(c(0,330), nbasis=6)
yield_par <- fdPar(basis, Lfdobj =2, lambda = 1)
yield_smooth <- smooth.basis(terms, yield_data, yield_par)
yield_smooth$gcv
plot(yield_smooth, col = "orange", add = TRUE)
basis <- create.bspline.basis(c(0,330), nbasis=6)
yield_par <- fdPar(basis, Lfdobj =2, lambda = 2)
yield_smooth <- smooth.basis(terms, yield_data, yield_par)
yield_smooth$gcv
plot(yield_smooth, col = "skyblue", add = TRUE)


########################################################################################



set.seed(201609)
m<-50; times<-seq(0,1,length=m)
range<-1; nu1=1/2; sig2<-1; nu2=3/2
Matern(.5,range=range,nu=nu1)

d_mat<-abs(outer(times,times,"-"))
C_1<-apply(d_mat,c(1,2),FUN=Matern,range=range,nu=nu1)
C_1<-C_1*sig2
C_1_sq<-sqrtm(C_1)
C_2<-apply(d_mat,c(1,2),FUN=Matern,range=range,nu=nu2)
C_2<-C_2*sig2
C_2_sq<-sqrtm(C_2)
Z<-rnorm(m)
X1<-C_1_sq%*%Z; X2<-C_2_sq%*%Z

par(mar=c(3,3,1,1),mfrow=c(1,2))
plot(times,X1,type="l")
plot(times,X2,type="l")

par(mar=c(3,3,1,1),mfrow=c(1,2))
Xd1<-(tail(X1,-1) - head(X1,-1))*m
Xd2<-(tail(X2,-1) - head(X2,-1))*m
plot(times[-1],Xd1,type="l")
plot(times[-1],Xd2,type="l")












