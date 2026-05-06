setwd("~/OneDrive/Desktop/IU Class Work/Spring 2026/STAT-S 650 (Time Series Analysis)")
TSdata = read.csv("seattle-weather.csv")

#1.Describe the data set, the variables, and the basic hypotheses you wish to test.

#View(TSdata)

N = nrow(TSdata)
year =  2012 + (1:N)/365

#1461 time points, one per day for four years, (from 2012 - 2015)
N
#Data over 4 years
N/365.25

#No NAs/missing values
sum(is.na(TSdata))

#Wind is my output variable
wind = TSdata$wind
#Max temperature is my input variable
temp_max_C = TSdata$temp_max
#Convert temperature from Celsius to Fahrenheit.
temp_max = (temp_max_C * (9/5)) + 32
date = 1:length(wind)

#Distribution of temperature
par(mfrow=c(2,1))
hist(temp_max,
     main = 'Distribution of Max Temp',
     col = "blue",
     xlab = 'Max Temp (Fahrenheit)',
     prob = TRUE)
#Distribution of Wind Speed
hist(wind,
     main = 'Distribution of Wind Speed',
     col = "red",
     xlab = 'Wind Speed (mph)',
     prob = TRUE)


#Hypotheses: 
  #1) Both max temperature and wind speed exhibit strong seasonal trends.
  #2) There is an inverse relationship between max temperature and wind speed.
  #3) Past values of max temperature and wind speed significantly predict current wind speed.


#2. Plot the time series data.

#Time series of max temp and wind speed
par(mfrow=c(2,1))
M = cbind(year,temp_max)
plot(M,type="l", main = 'Max Temperature over Time', col="blue", xlab="Year",ylab="Max Temperature")
M = cbind(year,wind)
plot(M,type="l", main = 'Wind Speed over Time', col="red", xlab="Year",ylab="Wind Speed")


# 31-day moving average (centered)
wind_smooth = stats::filter(wind, rep(1/31, 31), sides = 2)
temp_smooth = stats::filter(temp_max, rep(1/31, 31), sides = 2)
par(mfrow=c(2,1))

#Time series of max temp and wind speed with 31 day smoothing line
plot(year, temp_max,
     type="l",
     col="gray",
     main="Max Temperature (Raw vs 31-day Smoothed)",
     xlab="Year", ylab="Temperature (F)")
lines(year, temp_smooth, col="blue", lwd=2)

# Wind
plot(year, wind,
     type="l",
     col="gray",
     main="Wind Speed (Raw vs 31-day Smoothed)",
     xlab="Year", ylab="Wind Speed")
lines(year, wind_smooth, col="red", lwd=2)


#Periodogram before detrending
par(mfrow=c(2,1))
m=0
k = kernel("daniell",m)
H1 = spectrum(temp_max,k, log="no",plot=TRUE, main='Periodogram for Max Temp')
H2 = spectrum(wind,k, log="no",plot=TRUE, main='Periodogram for Wind Speed')


#Need to detrend seasonality to see autocorrelations and cross correlations without seasonality
f = 1/365
Xsin = sin(2*pi*f*date)
Xcos = cos(2*pi*f*date)
f2 = 2/365
Xsin2 = sin(2*pi*f2*date)
Xcos2 = cos(2*pi*f2*date)
#Detrending using two sinusoidal equations with sine and cosine in each
trend_temp <- lm(temp_max ~ Xsin+Xcos +Xsin2 + Xcos2)
trend_wind <- lm(wind ~ Xsin+Xcos +Xsin2 + Xcos2)

#Calculate Residuals for detrended data
resid_temp = residuals(trend_temp)
resid_wind = residuals(trend_wind)

# Plotting the Residuals of detrended data with Year on the X-axis
par(mfrow=c(2,1))
# Plot for Max Temp Residuals
plot(year, resid_temp, type='l', 
     xlab='Year', 
     ylab='Temp Residuals', 
     main='Max Temperature Residuals (Seasonality Removed)', 
     col="blue")
# Plot for Wind Speed Residuals
plot(year, resid_wind, type='l', 
     xlab='Year', 
     ylab='Wind Residuals', 
     main='Wind Speed Residuals (Seasonality Removed)', 
     col="red")


#Periodograms before and after detrending
par(mfrow=c(2,2))
#Periodograms before Detrending
H1 = spectrum(temp_max,k, log="no",plot=TRUE, main='Periodogram for Max Temp')
H2 = spectrum(wind,k, log="no",plot=TRUE, main='Periodogram for Wind Speed')
#Periodograms after Detrending
H3 = spectrum(resid_temp,k, log="no",plot=TRUE, main='Periodogram for Max Temp after Removing Seasonality')
H4 = spectrum(resid_wind,k, log="no",plot=TRUE, main='Periodogram for Wind Speed after Removing Seasonality')


#3. Report and summarize the autocorrelations and the cross correlations.

#Autocorrelations before and after detrending
par(mfrow=c(2,2))
#Autocorrelations before Detrending
acf(temp_max, 365, main="ACF of Max Temp")
acf(wind, 365, main="ACF of Wind Speed")
#Autocorrelations after Detrending
acf(resid_temp, 50, main="ACF of Max Temp after Removing Seasonality")
acf(resid_wind, 50, main="ACF of Wind Speed after Removing Seasonality")

#Cross Correlations before and after detrending
par(mfrow=c(2,1))
ccf(temp_max, wind, 365, main="CCF: Wind Speed (y) vs Max Temp (x)", ylab="CCF")
ccf(resid_temp, resid_wind, 50, main="CCF: Wind Speed (y) vs Max Temp (x) after Removing Seasonality", ylab="CCF")


#4. Report and summaries the spectral analysis of each series.

#Spectral analysis before and after detrending
par(mfrow=c(2,2))
m = 31
k = kernel("daniell",m)
#Spectral Analysis before Detrending
spec.pgram(temp_max, k, taper=0, log="no", main="Spectral Density for Max Temp with WS = 31")
spec.pgram(wind, k, taper=0, log="no", main="Spectral Density for Wind Speed with WS = 31")
#Spectral Analysis after Detrending
spec.pgram(resid_temp, k, taper=0, log="no", main="Spectral Density for Max Temp with WS = 31 (No Seasonality)")
spec.pgram(resid_wind, k, taper=0, log="no", main="Spectral Density for Wind Speed with WS = 31 (No Seasonality)")

#5. Fit a linear dynamic model to the data and compare several competing models.
#E.g., test for treatment effects of a input on the output of the system. Check the residuals to 
#make sure they are white noise.
library(forecast)
library(dplyr)

#Create wind and temp lag predictors
TSdata <- TSdata %>%
  mutate(
    temp_lag1 = lag(temp_max, 1),
    temp_lag2 = lag(temp_max, 2),
    wind_lag1 = lag(wind, 1)
  ) %>%
  na.omit()

TSdata$temp = TSdata$temp_max
weather = TSdata


# 0. Baseline: Just a lm of predicting wind with temperature
model0 <- lm(weather$wind ~ weather$temp)

#All are ARIMA/ARIMAX (1, 0, 0) AR(1) models, chosen by auto.arima for having the lowest BIC

# 1. Wind lag 1
model1 <- auto.arima(weather$wind, ic ='bic')

# 2. Wind lag 1 and Temp
model2 <- auto.arima(weather$wind, xreg = weather$temp, ic ='bic')

# 3. Wind lag 1 and Temp lag 1
model3 <- auto.arima(weather$wind, xreg = weather$temp_lag1, ic ='bic')

#4. Wind lag 1, Temp and Temp lag1
xreg = cbind(temp = weather$temp, temp_lag1 = weather$temp_lag1)
model4 <- auto.arima(weather$wind, xreg = xreg, ic ='bic') 

#5. Wind lag 1, Temp, Temp lag1, and Temp lag2
xreg = cbind(temp = weather$temp, temp_lag1 = weather$temp_lag1, temp_lag2 = weather$temp_lag2)
model5 <- auto.arima(weather$wind, xreg = xreg, ic ='bic') 

AIC(model0, model1, model2, model3, model4, model5)
BIC(model0, model1, model2, model3, model4, model5)

#Checking ARIMAX parameters (p, d, q), chosen for minimizing BIC for the model
#All are ARIMAX (1, 0, 0)/AR(1)
arimaorder(model1)
arimaorder(model2)
arimaorder(model3)
arimaorder(model4)
arimaorder(model5)

#Model coefficients
coef(model0)
coef(model1)
#Our best model (model 2) coefficients
coef(model2)
coef(model3)
coef(model4)
coef(model5)

summary(model2)

#P-values for Coefficients for best model(model 2)
# 1. Extract coefficients
coefs2 <- coef(model2)
# 2. Extract the diagonal of the variance-covariance matrix to get Standard Errors
se2 <- sqrt(diag(vcov(model2)))
# 3. Calculate the z-statistics
z_values2 <- coefs2 / se2
# 4. Calculate two-tailed p-values
p_values2 <- 2 * (1 - pnorm(abs(z_values2)))
# 5. Combine into a readable table
model2_stats <- cbind(Coefficients = coefs2, 
                      Std.Error = se2, 
                      z_stat = z_values2, 
                      p_value = p_values2)
print(round(model2_stats, 5))


#Chi-square test comparing nested models

#Chi-square comparing models 1 and 2
LRT_stat <- 2 * (as.numeric(logLik(model2)) - as.numeric(logLik(model1)))
p_val <- pchisq(LRT_stat, df = 1, lower.tail = FALSE)
cat("LRT Statistic:", LRT_stat, "\n")
cat("P-value:", p_val)
#P-value is 0.0001, adding current temp is imporving the model.

#Chi-square comparing models 2 and 4
LRT_stat <- 2 * (as.numeric(logLik(model4)) - as.numeric(logLik(model2)))
p_val <- pchisq(LRT_stat, df = 1, lower.tail = FALSE)
cat("LRT Statistic:", LRT_stat, "\n")
cat("P-value:", p_val)
#P-value is 0.87, adding temp lag 1 isn't adding anything.

#Chi-square comparing models 4 and 5
LRT_stat <- 2 * (as.numeric(logLik(model5)) - as.numeric(logLik(model4)))
p_val <- pchisq(LRT_stat, df = 1, lower.tail = FALSE)
cat("LRT Statistic:", LRT_stat, "\n")
cat("P-value:", p_val)
#P-value is 0.58, adding temp lag 2 isn't adding anything.

#This is consistent with our findings that model 2 is our best model.


#Residuals for best model (Model 2)

#ACF residuals
par(mfrow=c(1,1))
resids_best_model = residuals(model2)
acf(resids_best_model, 30, main ='ACF of Wind Residuals for Best Model (Model 2)')
acf(resids_best_model, 365, main ='ACF of Wind Residuals for Best Model (Model 2)')

#Spectral analysis of residuals (window size of 31 days)
m=31
k = kernel("daniell",m)
#Spectral analysis of residuals from model 2
H2 = spec.pgram(resids_best_model, k, taper=0, log="no")
nspec = H2$spec/sum(H2$spec)
plot(H2$freq,nspec,type="l",xlab="freq",xlim=c(0,.5),ylim=c(0,.01), main='Spectral Analysis of Residuals for Best Model (Model 2)')

checkresiduals(model2)

#Lijung Box-test has p-value of 0.3948. This is much larger
#than 0.05. There is no evidence of leftover autocorrelation
#in residuals.

#Time series of residuals
par(mfrow = c(1,1))
plot(residuals(model2),
     type = "l",
     main = "Residuals from Best Model (Model 2)",
     ylab = "Residuals")

#ACF of residuals
acf(residuals(model2),
    main = "ACF of Residuals from Best Model (Model 2)")

#Histogram of residuals
hist(residuals(model2),
     main = "Histogram of Residuals from Best Model (Model 2)",
     xlab = "Residuals")


#Adding spectral analysis, time series, and ACF residuals to one plot

m=31
k = kernel("daniell",m)
#Spectral analysis of residuals from model 2
H2 = spec.pgram(resids_best_model, k, taper=0, log="no")
par(mfrow = c(3,1))
nspec = H2$spec/sum(H2$spec)
plot(H2$freq,nspec,type="l",xlab="freq",xlim=c(0,.5),ylim=c(0,.01), main='Spectral Analysis of Residuals for Best Model (Model 2)')

#Time Series of residuals
resids_ts <- ts(resids_best_model, start = c(2012, 3), frequency = 365)
plot(resids_ts,
     type = "l",
     col = "red",
     xaxt = "n",  # suppress default x-axis
     main = "Residual Time Series (Model 2)",
     ylab = "Residuals",
     xlab = "Year")
#Add custom integer year ticks
years <- floor(time(resids_ts))
year_ticks <- unique(years)
axis(1,
     at = year_ticks,
     labels = year_ticks)

#ACF of residuals
acf(residuals(model2),
    main = "ACF of Residuals from Best Model (Model 2)")


#VAR
library(vars)
var_data <- na.omit(weather[, c("wind", "temp")])

#Select the lags for the VAR model with the lowest BIC
lag_select <- VARselect(var_data, lag.max = 10, type = "const")
#Model Recommends kag 1, has the lowest BIC (SC(n) = 1)
lag_select$selectio

#VAR model
var_model_bic <- VAR(var_data, p = 1, type = "const")
#VAR model equation and coefficients
summary(var_model_bic)

#Compare performance of model2 and best VAR model

# Check RMSE for Model 2 (1.300)
accuracy(model2)
#Check RMSE for VAR model predicting wind
# For the VAR, you have to extract the wind-specific residuals
var_resids_wind <- residuals(var_model_bic)[, "wind"]
rmse_var <- sqrt(mean(var_resids_wind^2))
#VAR RMSE 1.301
print(rmse_var)
#ARIMAX Model 2 is equal or slightly better than VAR model

#Does Temp Granger-cause wind?
causality(var_model_bic, cause = "temp")$Granger
# Does Wind Granger-cause Temp?
causality(var_model_bic, cause = "wind")$Granger
#Both Granger-cause each other and are statistically significant.


#6. Write your conclusions from the modeling analyses and relate these to your initial hypotheses.

#Hypotheses: 
#1) Both max temperature and wind speed exhibit strong seasonal trends.
#2) There is an inverse relationship between max temperature and wind speed.
#3) Past values of max temperature and wind speed significantly predict current wind speed.

#Conclusions:
#1) I was partially correct; the max temperature does exhibit strong seasonality. 
# Wind speed exhibits weak seasonality. This can be seen in the time series and spectral 
# analysis plots before detrending.

#2) At face value, I was correct; max temperature and wind speed have a strong inverse 
# relationship. However, this relationship is so strong due to seasonality, and once the 
# data has been detrended, the relationship is much weaker. After detrending, there is still a 
# weak inverse relationship. This can be seen when looking at the detrended CCF and the negative 
# coefficient value for temperature when predicting wind speed in our best model.


#3) I was incorrect. Wind speed lag 1 did end up being a significant predictor in our best model, 
# but past values of temperature were not statistically significant for predicting wind speed, 
# only current values of temperature. No values of wind speed or temperature were statistically 
# significant beyond lag 1. Our best model for wind speed only had wind speed lag 1 and 
# current temperature as predictors.

# Beyond our initial hypotheses, we also found that the VAR model that has the lowest BIC with 
# respect to our predictors only used 1 lag and had approximately equal or slightly worse performance 
# than our best ARIMAX model. Additionally, wind speed Granger-causes temperature, and temperature 
# Granger-causes wind speed. There is bidirectional Granger-causality between the two variables.
