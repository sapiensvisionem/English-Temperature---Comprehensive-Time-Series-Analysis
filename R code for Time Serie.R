# R Code

data <- read.table("C:/Users/jihun/Desktop/stat 361 data file.txt", quote="\"", comment.char="")
data.soi <- read.table("C:/Users/jihun/Desktop/soi index.txt", quote="\"", comment.char="")
# winter, spring, summer, autumn
Year = data$V1
DJF = data$V2
MAM = data$V3
JJA = data$V4
SON = data$V5

# exclude the first point
par(mfrow=c(2,2))
plot(Year[-1],DJF[-1], xlab="Year", ylab="Winter Temp", main="Winter")
plot(Year, MAM, xlab="Year", ylab="Spring Temp", main="Spring")
plot(Year, JJA, xlab="Year", ylab="Summer Temp", main="Summer")
plot(Year[-359], SON[-359], xlab="Year", ylab="Fall Temp", main="Fall")

# filtered by moving average
par(mfrow=c(2,4))
v1 = filter(MAM, sides=2, rep(1/10,10))
plot.ts(MAM, xlab="rescaled year", ylab= "Spring Temp",main="Spring time series")
plot.ts(v1, xlab="rescaled year", ylab="Moving Average", main="moving average")

v2 = filter(JJA, sides=2, rep(1/10,10))
plot.ts(JJA, xlab="rescaled year", ylab= "Summer Temp",main="Summer time series")
plot.ts(v2, xlab="rescaled year", ylab="Moving Average", main="moving average")

v3 = filter(SON[-359], sides=2, rep(1/10,10))
plot.ts(SON[-359], xlab="rescaled year", ylab= "Fall Temp",main="Fall time series")
plot.ts(v3, xlab="rescaled year", ylab="Moving Average", main="moving average")


v4 = filter(DJF[-1], sides=2, rep(1/10,10))
plot.ts(DJF[-1], xlab="rescaled year", ylab= "Winter Temp",main="Winter time series")
plot.ts(v4, xlab="rescaled year", ylab="Moving Average", main="moving average")

# ACF
par(mfrow=c(2,2))
acf(DJF[-1],lag.max=50,main="Winter")
acf(MAM,lag.max=50,main="Spring")
acf(JJA,lag.max=50,main="Summer")
acf(SON[-359],lag.max=50,main="Fall")

# plot with fitted linear trend line
par(mfrow=c(2,2))

plot(MAM, type="o", ylab="Global Temperature", main="Spring")
abline(lm(MAM ~time(MAM)))

plot(JJA, type="o", ylab="Global Temperature", main="Summer")
abline(lm(JJA ~time(JJA)))

plot(SON[-359], type="o", ylab="Global Temperature", main="Fall")
abline(lm(SON[-359] ~time(SON[-359])))

plot(DJF[-1], type="o", ylab="Global Temperature", main="Winter")
abline(lm(DJF[-1] ~time(DJF[-1])))

par(mfrow=c(1,1))
plot(time(SON[-359]), resid(lm(SON[-359] ~ time(SON[-359]))), type="o", main="detrended", xlab ="scaled time", ylab="spring")
plot(time(SON[-359:-358]), diff(SON[-359]), type="o", main="differenced", xlab ="scaled time", ylab="autumn")

df <- diff(SON[-359])
par(mfrow=c(2,1))
acf(df, main="ACF of Differenced")
pacf(df, main="PACF of Differenced")

library(astsa)
sarima(df, 0,0,1)
sarima(SON[-359], 0, 1, 1)

# broken stick linear regression
I = ifelse(Year > 2000, 1, 0)
fit1 = lm(MAM ~ Year)
fit2 = lm(MAM ~ Year + Year:I)
anova(fit1,fit2)

summary(lm(JJA ~ time(JJA)))


# Periodogram
n=length(MAM)
I = abs(fft(MAM))^2/n
P = (4/n)*I[1:(n/2)]
f = 0:((n/2)-1)/n
par(mfrow=c(1,1))
plot(f, P, type="l", xlab="Frequency", ylab="Scaled Periodogram")

par(mfrow=c(2,2))
plot(predict(arima(df[-342:-358], order = c(0,0,1)), n.ahead = 17)$pred[-1],c(2001:2016))
plot(c(2001:2016),df[343:358])
pre = predict(arima(df[-342:-358], order = c(0,0,1)), n.ahead = 17)$pred[-1]
length(pre)
plot(c(2001:2016),pre)
plot(c(2001:2016),df[343:358])


# Trend plus periodic (global)
par(mfrow=c(1,1))
t = Year - mean(Year)
t2 = t^2; t3 = t^3
cs=cos(2*pi*t); sn=sin(2*pi*t)
reg1 = lm(MAM ~ t + t2 + t3)
reg2 = lm(MAM ~ t + t2 + t3 + cs + sn)
plot(MAM, type="p", ylab="Spring Temperature")
lines(fitted(reg1)); lines(fitted(reg2))

# Kernel smoother
plot(MAM, type="p",ylab="Spring temperature", main="Kernel Smoother")
lines(ksmooth(time(MAM), MAM, "normal", bandwidth=20))
lines(ksmooth(time(MAM), MAM, "normal", bandwidth=100))

# Nearest Neighbor and Lowess
par(mfrow=c(2,1))
plot(MAM, type="p", ylab="Spring Temperature", main="nearest neighbor")
lines(supsmu(time(MAM), MAM, span=.2))
lines(supsmu(time(MAM), MAM, span=.005))
plot(MAM, type="p", ylab="Spring Temperature", main="Lowess")
lines(lowess(MAM, f=0.01)); lines(lowess(MAM, f=0.2))

# Smoothing spline
par(mfrow=c(2,1))
plot(time(MAM),MAM, main="Smoothing spline")
lines(smooth.spline(MAM, spar=1))
plot(time(MAM),MAM, main="Smoothing spline")
lines(smooth.spline(MAM, spar=1/2))

# periodogram of detrended
par(mfrow=c(1,1))
fit <- lm(MAM ~ time(MAM))
res <- resid(fit)
n = length(res)
ord = 1:((n-1)/2)
freq = (ord)/n
per = as.vector(abs(fft(res))^2/n)
# P = Mod(2*fft(res)/100)^2 ; Fr = 0:99/100
plot(freq,per[ord],type="l",main="scaled periodogram")

# PACf vs acf
par(mfrow=c(2,1))
pacf(MAM)
acf(MAM)

# AR(2) regression
regr = ar.ols(MAM[-(350:359)], order=2, demean=FALSE, intercept=TRUE)
fore = predict(regr, n.ahead=10)
par(mfrow=c(1,1))
ts.plot(fore$pred)
ts.plot(MAM[350:359], fore$pred, ylab="MAM temp")
U=fore$pred + fore$se
L = fore$pred - fore$se
xx=c(time(U), rev(time(U))); yy = c(L, rev(U))
polyhon(xx, yy, border=8, col=gray(0.6, alpha= 0.2))
lines(fore$pred)

# Yule Walker Estimation
MAM.yw = ar.yw(MAM, order=2)
MAM.yw$x.mean # 8.169 mean estimate
MAM.yw$ar # (0.2145, 0.2081) coefficient estimates
sqrt(diag(MAM.yw$asy.var.coef)) # 0.05183, 0.05183 stabdard errors
MAM.yw$var.pred # 0.7096 error variance estimate

MAM.pr = predict(MAM.yw, n.ahead=24)
plot(c(2017:2040),MAM.pr$pre)
# ts.plot(MAM, MAM.pr$pred)

# log-difference like GNP plot
ldMAM = log(diff(MAM)) # diff(MAM) has zeroes'
acf(ldMAM, na.action = na.pass)
pacf(ldMAM)
# AR (1)
sarima(ldMAM, 1, 0, 0)
# MA (2)
sarima(ldMAM, 0,0, 2)
ARMAtoMA(ar = 0.35, ma=0, 10)


# diagnostics of ARIMA model
sarima(gnpgr, 1, 0, 0)
# (xt - one-step-ahead prediction) / 
# normality of residuals - must be standard normal
# zero correlation
# correlations shouldn't be too large - Ljung-Box test
# if there is autocorrelation, run WLS

# Regression with Lagged Variables

# spectrum
arma.spec(asdlog="no", main="MAM")
