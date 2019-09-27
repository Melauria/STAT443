######################################################################
############################### DATA #################################
######################################################################
library(forecast)
data = read.csv("cleaned_data.csv", header = T)

# deal with data
y = ts(log(data$value))
train_y = ts(y[1:288])
holdo_y = ts(y[289:323])

df = diff(y)
train_df = ts(df[1:288])
holdo_df = ts(df[289:322])


nominal_y = ts(data$nominal)
real_y = ts(data$real)
nominal_df = diff(nominal_y)
real_df = diff(real_y)

# plots
par(mfrow=c(2,3))
plot(y, ylab = "log(value)", main = "log of traveller flow volumn to BC per month 
     from 1992-01 to 2018-11")
acf(y, main = "ACF")
pacf(y, main = "PACF")
plot(df, ylab = "differenced log(value)", main = "differenced series of log flow volumn")
acf(df, main = "ACF")
pacf(df, main = "PACF")

######################################################################
############################## MODELS ################################
######################################################################

# persistence
# rmse = 0.01891757
ntrain = length(train_y)
nholdo = length(holdo_y)
mse_p = 0
fc_p = train_y[ntrain]
zt_p = holdo_y[1] 
fcerror_p = zt_p-fc_p
mse_p = mse_p + fcerror_p^2
for (i in c(2:nholdo)){
  fc_p = holdo_y[i-1] 
  zt_p = holdo_y[i] 
  fcerror_p = zt_p-fc_p
  mse_p = fcerror_p^2 + mse_p
}
rmse_p = sqrt(mse_p/nholdo)
rmse_p

# average
# rmse = 0.1628434
ntrain = length(train_y)
nholdo = length(holdo_y)
mse_avg = 0
cumsum = sum(train_y)
fc_avg = cumsum/ntrain
zt_avg = holdo_y[1] 
fcerror_avg = zt_avg-fc_avg
mse_avg = mse_avg+(fcerror_avg)^2 
for (i in c(2:nholdo)){
  cumsum <- holdo_y[i-1] + cumsum
  fc_avg <- cumsum/(ntrain+i-1)
  zt_avg <- holdo_y[i] 
  fcerror_avg <- zt_avg-fc_avg
  mse_avg <- mse_avg + fcerror_avg^2
}
rmse_avg <- sqrt(mse_avg/nholdo)
rmse_avg

# simple exponential smoothing
# rmse = 0.01929983
esfit = HoltWinters(train_y, beta=F, gamma=F)
nholdo = length(holdo_y)
alpha = esfit$alpha 
lev = esfit$coef
sse = 0
fc = lev
zt = holdo_y[1]
newv = lev
fcerror = zt-fc
sse = sse + fcerror^2
for(i in 2:nholdo){
  newv = alpha*holdo_y[i-1]+(1-alpha)*newv
  fc = newv
  zt = holdo_y[i]
  fcerror = zt-fc
  sse = sse + fcerror^2
}
rmse = sqrt(sse/nholdo)
rmse

# holt linear exponential smoothing
# rmse = 0.01857013
hlfit = HoltWinters(train_y, gamma = F) 
nholdo = length(holdo_y)
alpha = hlfit$alpha 
beta = hlfit$beta
lev = hlfit$coef[1]
trend = hlfit$coef[2]
sse = 0
fc = lev+trend
zt = holdo_y[1]
newfc = fc
fcerror = zt-fc
sse = sse + fcerror^2
vprev = lev
bprev = trend
for (i in 2:nholdo){
  vnew = alpha*holdo_y[i-1]+(1-alpha)*fc
  bnew = beta*(vnew-vprev)+(1-beta)*bprev
  fc = vnew + bnew
  zt = holdo_y[i]
  fcerror = zt-fc
  sse = sse + fcerror^2
  vprev = vnew
  bprev = bnew
}
rmse = sqrt(sse/nholdo)
rmse

# MA(2)
# rmse = 0.01686141
ma2 = Arima(train_df, order = c(0,0,2),include.mean=T, method = "CSS")
ma2.fit = Arima(holdo_df, model = ma2)
accuracy(ma2.fit)

# ARMA(1,1)
# rmse = 0.01716506
arma11 = Arima(train_df, order = c(1,0,1),include.mean=T, method = "CSS")
arma11.fit = Arima(holdo_df, model = arma11)
accuracy(arma11.fit)

# ARIMA(1,1,1)
# rmse = 0.01780088
arima111 = Arima(train_y, order = c(1,1,1),include.mean=T, method = "CSS")
arima111.fit = Arima(holdo_y, model = arima111)
accuracy(arima111.fit)

# ARIMA(0,1,2)
# rmse = 0.01743835
arima012 = Arima(train_y, order = c(0,1,2),include.mean=T, method = "CSS")
arima012.fit = Arima(holdo_y, model = arima012)
accuracy(arima012.fit)

# ARMAX(0,0,2)
# rmse = 0.01714317
armax2_nominal = Arima(train_df, order = c(0,0,2),include.mean=T, method = "CSS", xreg=nominal_df[1:288])
armax2_nominal.fit = Arima(holdo_df, model = armax2_nominal, xreg=nominal_df[289:322])
accuracy(armax2_nominal.fit)
# rmse = 0.01709704
armax2_real = Arima(train_df, order = c(0,0,2),include.mean=T, method = "CSS", xreg=real_df[1:288])
armax2_real.fit = Arima(holdo_df, model = armax2_real, xreg=real_df[289:322])
accuracy(armax2_real.fit)
# rmse = 0.0171441
armax2_both = Arima(train_df, order = c(0,0,2),include.mean=T, method = "CSS", 
               xreg=cbind(real_df[1:288],nominal_df[1:288]))
armax2_both.fit = Arima(holdo_df, model = armax2_both, xreg=cbind(real_df[289:322],nominal_df[289:322]))
accuracy(armax2_both.fit)

# ARMAX(0,0,1)
# rmse = 0.01765109
armax1_nominal = Arima(train_df, order = c(0,0,1),include.mean=T, method = "CSS", xreg=nominal_df[1:288])
armax1_nominal.fit = Arima(holdo_df, model = armax1_nominal, xreg=nominal_df[289:322])
accuracy(armax1_nominal.fit)
# rmse = 0.01758477
armax1_real = Arima(train_df, order = c(0,0,1),include.mean=T, method = "CSS", xreg=real_df[1:288])
armax1_real.fit = Arima(holdo_df, model = armax1_real, xreg=real_df[289:322])
accuracy(armax1_real.fit)
# rmse = 0.01788879
armax1_both = Arima(train_df, order = c(0,0,1),include.mean=T, method = "CSS", 
               xreg=cbind(real_df[1:288],nominal_df[1:288]))
armax1_both.fit = Arima(holdo_df, model = armax1_both, xreg=cbind(real_df[289:322],nominal_df[289:322]))
accuracy(armax1_both.fit)

# ARIMAX(0,1,1)
# rmse = 0.01793455
arimax011_nominal = Arima(train_y, order = c(0,1,1),include.mean=T, method = "CSS", xreg=nominal_y[1:288])
arimax011_nominal.fit = Arima(holdo_y, model = arimax011_nominal, xreg=nominal_y[289:323])
accuracy(arimax011_nominal.fit)
# rmse = 0.01787375
arimax011_real = Arima(train_y, order = c(0,1,1),include.mean=T, method = "CSS", xreg=real_y[1:288])
arimax011_real.fit = Arima(holdo_y, model = arimax011_real, xreg=real_y[289:323])
accuracy(arimax011_real.fit)
# rmse = 0.01818749
arimax011_both = Arima(train_y, order = c(0,1,1),include.mean=T, method = "CSS", 
                       xreg=cbind(nominal_y[1:288],real_y[1:288]))
arimax011_both.fit = Arima(holdo_y, model = arimax011_both, xreg=cbind(nominal_y[289:323],real_y[289:323]))
accuracy(arimax011_both.fit)

# ARIMAX(0,1,2)
# rmse = 0.01769919
arimax012_nominal = Arima(train_y, order = c(0,1,2),include.mean=T, method = "CSS", xreg=nominal_y[1:288])
arimax012_nominal.fit = Arima(holdo_y, model = arimax012_nominal, xreg=nominal_y[289:323])
accuracy(arimax012_nominal.fit)
# rmse = 0.01765058
arimax012_real = Arima(train_y, order = c(0,1,2),include.mean=T, method = "CSS", xreg=real_y[1:288])
arimax012_real.fit = Arima(holdo_y, model = arimax012_real, xreg=real_y[289:323])
accuracy(arimax012_real.fit)
# rmse = 0.01774462
arimax012_both = Arima(train_y, order = c(0,1,2),include.mean=T, method = "CSS",
                xreg=cbind(nominal_y[1:288],real_y[1:288]))
arimax012_both.fit = Arima(holdo_y, model = arimax012_both, xreg=cbind(nominal_y[289:323],real_y[289:323]))
accuracy(arimax012_both.fit)

# final prediction using MA(2)
# MA(2)
# rmse = 0.01686141
ma2final = Arima(df, order = c(0,0,2),include.mean=T, method = "CSS")
pred = predict(ma2final, n.ahead=12)$pred
diff.pred = c(y[length(y)])
for (i in c(1:12)) {
  print(pred[i])
  print(diff.pred[i])
  diff.pred[i+1] = pred[i] + diff.pred[i]
}
pred=exp(diff.pred)
round(pred)