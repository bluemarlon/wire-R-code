library(TDA)

# ex.2.1
set.seed(1); PC <- circleUnif(n = 60, r = 1)
plot(PC, main = "(a)")
(pers.diag.1 <- ripsDiag(X=PC, maxdimension = 1, maxscale = max(dist(PC))))
plot(x=1,y=1,type="n",ylim=c(0,2),xlim=c(0,2),ylab="death",xlab="birth",main="(b)")
abline(v = 0, lty = 2)
plot(PC, pch = 16, cex = 5, col = "blue", main = "(c)")
death.time = sort(pers.diag.1$diagram[pers.diag.1$diagram[, 1]==0, 3])
plot(x = rep(0, sum(death.time<=0.1)), y = death.time[which(death.time<=0.1)],
     ylim = c(0, 2), xlim = c(0, 2), ylab = "death")
abline(v = 0, lty = 2); abline(h = 0.1, lty = 2)
plot(PC, pch = 16, cex = 12, col = "blue", main = "(e)")
plot(x = rep(0, sum(death.time<=0.32)), y = death.time[which(death.time<=0.32)],
     ylim = c(0, 2), xlim = c(0, 2), ylab = "death", xlab = "birth", main = "(f)")
abline(h = pers.diag.1$diagram[pers.diag.1$diagram[, 1]==1, 2], lty = 2)
abline(v = 0, lty = 2); abline(v = 0.32, col = "red", lty = 1)
plot(PC, pch = 16, cex = 40, col = "blue", main = "(g)")
plot(pers.diag.1$diagram, main = "(h)")

# ex.2.2
set.seed(1); PC <- circleUnif(n = 60, r = 1)
dist.PC <- dist(PC)
(pers.diag.2=ripsDiag(X=dist.PC,dist="arbitrary",maxdimension=1,maxscale=max(dist.PC))$diagram)

# ex.2.3
set.seed(1); PC <- circleUnif(n = 60, r = 1)
dist.PC <- dist(PC)
(pers.diag.2=ripsDiag(X=dist.PC,dist="arbitrary",maxdimension=1,maxscale=max(dist.PC))$diagram)
set.seed(2); PC2 <- circleUnif(n = 60, r = 1)
pers.diag.3 <- ripsDiag(X=PC2, maxdimension = 1,maxscale = max(dist(PC2)))$diagram
wasserstein(pers.diag.2, pers.diag.3, p=1, dimension = 0)
bottleneck(pers.diag.2, pers.diag.3, dimension = 0)
dist.PC.man <- dist(PC, method = "manhattan"); max.dist=max(dist.PC.man)
pers.diag.4=ripsDiag(X=dist.PC.man,dist="arbitrary",maxdimension=1,
                     maxscale=max.dist)$diagram
wasserstein(pers.diag.2, pers.diag.4, p=1, dimension = 0)
bottleneck(pers.diag.2, pers.diag.4, dimension = 0)

# ex.2.4
T = 480
per1=12;ts1 = cos(1:T*2*pi/per1);d=2;
ts.ex = ts1
tau <- which(abs(acf(ts.ex, plot = F)$acf) < 2/sqrt(T))[1]-1
PC=t(purrr::map_dfc(1:(T-(d-1)*tau+1),~ts.ex[seq(from=.x, by=tau, length.out=d)]))
diag=ripsDiag(PC, maxdimension=1, maxscale=max(dist(PC)))
ts.plot(ts.ex);plot(PC,xlab ="x1",ylab="x2",main="PC");plot(diag$diagram)

# ex.2.5 diag for plot
library(tidyverse)
T = 480;per1=12;ts1 = cos(1:T*2*pi/per1);
x.ts = ts1; d=15; N=201; T1 = 216;
x.ts <- pracma::movavg(x.ts, 5, type = "s") #step 0
sp.ts <- stats::spline(1:T*2*pi/T, x.ts, n=T1)$y #step 1.1
PC <- plyr::ldply(map(1:N, ~sp.ts[.x:(.x+d-1)]))#step 1.2
X.PC=t(apply(PC,1,FUN=function(x){(x-mean(x))/sqrt(sum((x-mean(x))^2))})) #step2
diag <- ripsDiag(X=X.PC, maxdimension = 1, maxscale = sqrt(3))#step 3

# ex.3.1
funval = c(1, 0.5, 1, 1.5, 0.5, 0, 1, 1, 0.5, 1)
(pers.diag.4 <- gridDiag(FUNvalues = funval, sublevel = TRUE))
plot(funval, x=1:10, type = "l",yaxt='n',ylim = c(0,2), cex.axis=1.4, xlab = "z",
     ylab= "y", cex.lab= 1.3, cex=1.2, lwd=1.5, lty=1, pch=1, bty='n', main="(a)")
points(x=6, y=0, pch=16, col="blue", type = "p")
ticks<-c(0, 0.5, 1, 1.5, 2); axis(2,at=ticks,labels=ticks)
abline(a=0, b=0, lty=2, lwd=1,pch=1)
plot(rep(0,11),x=0:10/5,type = "l",cex.axis=1.4,xaxt='n',xlab= "birth",yaxt='n',
     ylim=c(0,2),ylab="death",cex.lab=1.3,cex=1.2,lwd=1.5,lty=2,pch=1,bty='n',main="(b)")
axis(1,at=ticks,labels=ticks); axis(2,at=ticks,labels=ticks); abline(v=0, lty=2)
plot(funval, x=1:10, type = "l",yaxt='n',ylim = c(0,2),cex.axis=1.4,xlab = "z",
     ylab = "y",cex.lab= 1.3, cex=1.2,lwd=1.5, lty=1, pch=1, bty='n', main = "(c)")
axis(2,at=ticks,labels=ticks); abline(a=0.5, b=0, lty=2, lwd=1,pch=1)
points(x=6, y=0, pch=16, col="green4", type = "p",xlab = "z", ylab= "y")
points(x=c(2,5,9), y=rep(0.5,3), pch=16, col="blue", type = "p",xlab="z",ylab="y")
segments(x0=6, y0=0, x1=5, y1=0.5, lty = 1, pch=1, lwd=2.5, col = "blue")

# ex.3.2
set.seed(1); PC <- circleUnif(n = 60, r = 1)
m0=0.05; by <- 0.065; Xlim <- range(PC[,1]); Ylim <- range(PC[,2])
(pers.diag.5 = gridDiag(X=PC, FUN=dtm, lim=cbind(Xlim, Ylim), by=by, m0=m0) )
par(mfrow = c(1,2))
Xseq <- seq(from = Xlim[1], to = Xlim[2], by = by)
Yseq <- seq(from = Ylim[1], to = Ylim[2], by = by)
Grid <- expand.grid(Xseq, Yseq); DTM = dtm(X = PC, Grid = Grid, m0 = m0)
persp(x = Xseq, y = Yseq,z = matrix(DTM, nrow = length(Xseq), ncol = length(Yseq)),
      xlab = "", ylab = "", zlab = "", theta = -20, phi = 35, scale = FALSE,
      expand = 2, col = "red", border = NA, ltheta = 50, shade = 0.5,main = "(a)")
plot(pers.diag.5[["diagram"]], main = "(b)")

# ex.3.3
T = 480
per1=12;ts1 = cos(1:T*2*pi/per1);d=2;
x.1 = 1:length(ts1); ts1 = lm(ts1~x.1)$residuals; ts1 = ts1/sd(ts1) #Step 1
spc.t1=spec.pgram(ts1, kernel=kernel("modified.daniell", c(1)), plot = F)#Step 2
#Step 3
PD.t1=gridDiag(FUNvalues=spc.t1$spec,location = FALSE, sublevel=TRUE)$diagram

# ex.4.1
set.seed(2); PC2 <- circleUnif(n = 60, r = 1)
pers.diag.3 <- ripsDiag(X=PC2, maxdimension = 1,maxscale = max(dist(PC2)))
Diag = pers.diag.3$diagram; Land <- c(); k=1; threshold=1; tot_sum=1
while(threshold>0.01){
  Land <- cbind(Land, landscape(Diag = Diag, dimension = 0, KK = k,
                                tseq = seq(min(Diag[,2:3]), max(Diag[,2:3]), length=500)))
  pre_sum <- sum(abs(Land[, ncol(Land)]))
  threshold = abs(tot_sum - pre_sum)
  tot_sum = pre_sum
}
