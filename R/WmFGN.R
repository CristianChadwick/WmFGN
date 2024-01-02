# The codes of this paper where developed in MATLAB
# and originally commented in Spanish. The most important
# lines have been translated into English, and to R,
# to improve the accessibility of a broader readership

# Code developed by Cristian Chadwick
# email: cristian.chadwick@uai.cl
# email: cchadwi1@uc.cl
# date: 10-21-2023

# The "rm()" command erases the elements of the Environment, in this case the
# entire list 
rm(list=ls())
# The command "graphics.off()" closes the open plots.
graphics.off()

# Set your pathway to the file where the “historicos_biobio_CEN.csv” file is stored
setwd("C:/FilePathWay…")

# library() allows activating packages, the “openxlsx” will allow to read csv or Excel
# files
library(openxlsx)  

# Read the "historicos_biobio_CEN.csv" file with the command "read.csv()"
# this allows importing the data, with the data frame format.
Qimp <- read.csv("historicos_biobio_CEN.csv")

# It is always good to verify that the file was properly imported to R, for 
# this the View() is quite useful
View(Qimp)

# We will work with hydrological years, reason why we use data from April to
# march. To do this we have to extract the rows "1:(dim(Qimp)[1]-1)", 
# while in the columns form the forth onward: “(4:dim(Qimp)[2])”
Q_dat<-as.matrix(Qimp[1:(dim(Qimp)[1]-1),4:dim(Qimp)[2]])

# "n" is the number of years of data that we have (59 years)
n <- dim(Q_dat)[1]/12
# E is the number of stations (9 in this case)
E <- dim(Q_dat)[2]  

# Let’s generate a time (Tiempo) vector to be able to plot the time series
Tiempo<-Qimp[1:(dim(Qimp)[1]-1),2]+(Qimp[1:(dim(Qimp)[1]-1),3]-1)/12

# Plotting the time series (just for the first station Q_dat[,1], for other 
# one has to change the number):
plot(Tiempo,Q_dat[,1],type = "l",xlab = "Años",ylab = "Q [m^3/s]")

# The time series of streamflows, must be converted from a 2D matrix
# into a 3D matrix, passing from the dimensions time x stations to 
# Months x Years x Stations [12 x 59 x 9]
Q_3D <- array(Q_dat,c(12,n,E))
# Then one can permute the order of the dimensions with “aperm” to:
# Years x Months x Stations
Q_3D <- aperm(Q_3D,c(2,1,3))

# After this, estimate the mean and the standard deviation for each month
# and station:
mu_x <- apply(Q_3D,2:3,mean)
sig_x <- apply(Q_3D,2:3,sd)
# Estimate the parameters of the log-normal distributions 
# with the method of moments, one distribution per month 
# and station:
sig_y <- sqrt(log(1+(sig_x/mu_x)^2))
mu_y <- log(mu_x)-sig_y^2/2

# The streamflow values (from Q_3D) are normalized and stored in the 
# matrix Q_3D_N
Q_3D_N <- array(0,dim(Q_3D))
for (i in 1:n){
  # Normalization and standardization 
  Q_3D_N[i,,] <- (log(Q_3D[i,,])-mu_y)/sig_y
}
# Before, normalizing it is important to verify that they have a normal
# distribution (done before in the selection of the streamflows).

# Defining the number of years to be generated
ng <- 10000

# A random uniform 3D matrix is generated (X_3D), which will be used for
# mFGNS, SmFGN and WmFGN
X_3D <- rnorm((ng+1)*12*E,mean = 0,sd = 1)
# the data is reshaped into a 3D matrix with dimensions 
# ng+1 x 12 x E 
# that correspond to the Number of years +1 x Months x Stations 
X_3D <- array(X_3D,c((ng+1),12,E))

#  For the multi-site mFGN method, one has to resample the historical data,
# for this we have to generate random numbers between 1 and “n” (number
# of observed years). The numbers to be generated need to add one (which is
# lost in the process). 
AUX <- sample.int(n,(ng+1)*12,replace = T)
AUX <- array(AUX,c((ng+1),12))
# The X random matrix (is generated) from the Q_3D_N matrix,
# by resampling its values.
X_3D_Sem <- array(0,dim(X_3D))
for (i in 1:(ng+1)){
  for (j in 1:12){
    X_3D_Sem[i,j,] <- Q_3D_N[AUX[i,j],j,]
  }
}

# The matrix where the results from the mFGN will be stored is
# generated
mFGN <- array(0,c((ng),12,E))
# mFGN Model
for (k in 1:E){
  # The data from the “k” station are extracted from Q_3D_N
  Q_aux <- Q_3D_N[,,k]
  # The transformation that permutes the second half of the year and
  # the first half is performed (respecting the years). 
  Q_aux_p <- cbind(Q_aux[1:(n-1),(12/2+1):12],Q_aux[(1+1):n,1:(12/2)])
  # The Cholesky factorization is estimated from the correlation matrix from
  # "Q_aux" and "Q_aux_p"
  Q <- chol(cor(Q_aux))
  Q_p <- chol(cor(Q_aux_p))
  # The random number (by resampling the historical numbers), of each station
  # number “k” are extracted
  X_aux <- X_3D_Sem[,,k]
  # The first and second half of the year of “X” are permuted:
  X_aux_p <- cbind(X_aux[1:(ng+1-1),(12/2+1):12],X_aux[(1+1):(ng+1),1:(12/2)])
  # The mFGN model is generated:
  Z <- X_aux%*%Q
  Z_p <- X_aux_p%*%Q_p
  # Both halves of the mFGN model are binded:
  Zc <- cbind(Z_p[,(12/2+1):12],Z[(1+1):(ng+1),(12/2+1):12])
  # The results are finally stored in a 3D matrix:
  mFGN[,,k] <- Zc
}
# Similar process, but for the mFGNS model is started:
mFGNS <- array(0,c((ng),12,E))
# Modelo mFGNS
for (k in 1:E){
  # The streamflow data of station “k” is extracted:
  Q_aux <- Q_3D_N[,,k]
  # The transformation that permutes the second half of the year and
  # the first half is performed (respecting the years). 
  Q_aux_p <- cbind(Q_aux[1:(n-1),(12/2+1):12],Q_aux[(1+1):n,1:(12/2)])
  # The Cholesky factorization is estimated from the correlation matrix from
  # "Q_aux" and "Q_aux_p"
  Q <- chol(cor(Q_aux))
  Q_p <- chol(cor(Q_aux_p))
  # The random numbers of the parametrized synthetic matrix are extracted 
  # for the station “k”
  X_aux <- X_3D[,,k]
  # The first and second half of the year of “X” are permuted:
  X_aux_p <- cbind(X_aux[1:(ng+1-1),(12/2+1):12],X_aux[(1+1):(ng+1),1:(12/2)])
  # The mFGN model is applied:
  Z <- X_aux%*%Q
  Z_p <- X_aux_p%*%Q_p
  # Both halves of the mFGN model are binded:
  Zc <- cbind(Z_p[,(12/2+1):12],Z[(1+1):(ng+1),(12/2+1):12])
  # The results are finally stored in a 3D matrix:
  mFGNS[,,k] <- Zc
}

# The spatial correlation is added in the following “for” cycle
# for each month:
for (j in 1:12){
  # The streamflow per month is extracted:
  Q_aux_E <- Q_3D_N[,j,]
  # The Cholesky factorization is estimated for the spatial correlation
  Q_E <- chol(cor(Q_aux_E))
  # The spatial correlation is added to the "mFGNS" model
  mFGNS[,j,] <- mFGNS[,j,]%*%Q_E
}

# The matrix for the "SmFGN" model is initialized:
SmFGN <- array(0,c((ng),12,E))
# A matrix "X_3D_AUX" with the same dimensions than "X_3D" is generated, but this 
# will incorporate the spatial correlation
X_3D_AUX <- array(0,dim(X_3D))
# SmFGN model
for (j in 1:12){
  # The streamflow of each month is extracted
  Q_aux_E <- Q_3D_N[,j,]
  # The Cholesky factorization is estimated for the spatial correlation
  Q_E <- chol(cor(Q_aux_E))
  # The spatial correlation is added to the “X_3D” matrix to obtain:
  # the “X_3D_AUX” matrix
  X_3D_AUX[,j,] <- X_3D[,j,]%*%Q_E
}

# In a second step the stations are extracted to incorporate the temporal correlation
for (k in 1:E){ 
  # The streamflow data of station “k” is extracted:
  Q_aux <- Q_3D_N[,,k] 
  # The transformation that permutes the second half of the year and
  # the first half is performed (respecting the years). 
  Q_aux_p <- cbind(Q_aux[1:(n-1),(12/2+1):12],Q_aux[(1+1):n,1:(12/2)]) 
  # The Cholesky factorization is estimated from the correlation matrix from
  # "Q_aux" and "Q_aux_p"
  Q <- chol(cor(Q_aux))
  Q_p <- chol(cor(Q_aux_p))
  # The random numbers of the parametrized synthetic matrix are extracted 
  # for the station “k”. Note that the random numbers already have the spatial
  # correlation
  X_aux <- X_3D_AUX[,,k] 
  # The first and second half of the year of “X” are permuted:
  X_aux_p <- cbind(X_aux[1:(ng+1-1),(12/2+1):12],X_aux[(1+1):(ng+1),1:(12/2)])
  # The mFGN model is applied:
  Z <- X_aux%*%Q
  Z_p <- X_aux_p%*%Q_p
  # Both halves of the mFGN model are binded:
  Zc <- cbind(Z_p[,(12/2+1):12],Z[(1+1):(ng+1),(12/2+1):12]) 
  # The results are finally stored in a 3D matrix:
  SmFGN[,,k] <- Zc
}

# Finally, the optimization procedure is performed, for that 
# Alpha values between 0 and 1 are generated:
Alpha <- seq(from=0, to=1, by=0.01)


# The vectors that will store the error are generated:
Err_Temp_Med<- array(0, c(length(Alpha)))
Err_Temp_Max<- array(0, c(length(Alpha)))
Err_Spa_Med<- array(0, c(length(Alpha)))
Err_Spa_Max<- array(0, c(length(Alpha)))

for (i in 1:length(Alpha)){
  # The WmFGN model is generated as by adding the “mFGNS” and “SmFGN” 
  # weighted by the Alpha values
  WmFGN <- (1-Alpha[i])*mFGNS + Alpha[i]*SmFGN
  
  # Using a “for” loop one can estimate the temporal correlation matrix and
  # its errors.
  
  # The temporal error matrix is initialized
  E_corr_Temp <-array(0,c(12,12,E)) 
  for (k in 1:E){
    # The dimension of the station “k” is extracted:
    Q_aux <- Q_3D_N[,,k]
    # The correlation matrix is estimated
    QQ_Obs <- cor(Q_aux)
    # The same procedure is repeated for the WmFGN
    Q_aux <- WmFGN[,,k]
    QQ_WmFGN <- cor(Q_aux)
    
    # The absolute difference are estimated (between observed and generated data)
    E_corr_Temp[,,k]<-abs(QQ_Obs-QQ_WmFGN)
  }
  # The mean and maximum error are estimated
  Err_Temp_Med[i]<- mean(E_corr_Temp)
  Err_Temp_Max[i]<- max(E_corr_Temp)
  
  # The spatial error matrix is initialized
  E_corr_Spa <-array(0,c(E,E,12)) 
  for (j in 1:12){
    # The dimension of the month “j” is extracted:
    Q_aux <- Q_3D_N[,j,]
    # The correlation matrix is estimated
    QQ_Obs <- cor(Q_aux)
    # The same procedure is repeated for the WmFGN
    Q_aux <- WmFGN[,j,]
    QQ_WmFGN <- cor(Q_aux)
    
    # The absolute difference are estimated (between observed and generated data)
    E_corr_Spa[,,j]<-abs(QQ_Obs-QQ_WmFGN)
  }
  # The mean and maximum error are estimated
  Err_Spa_Med[i]<- mean(E_corr_Spa)
  Err_Spa_Max[i]<- max(E_corr_Spa)
}

# Similar procedure is performed to estimate the error in the mFGN

# The temporal error matrix is initialized
E_corr_Temp_mFGN <-array(0,c(12,12,E)) 
for (k in 1:E){
  # The dimension of the station “k” is extracted:
  Q_aux <- Q_3D_N[,,k]    
  # The correlation matrix is estimated
  QQ_Obs <- cor(Q_aux)
  # The same procedure is repeated for the mFGN
  Q_aux <- mFGN[,,k]
  QQ_mFGN <- cor(Q_aux)
  
  # The absolute difference are estimated (between observed and generated data)
  E_corr_Temp_mFGN[,,k]<-abs(QQ_Obs-QQ_mFGN)
}
# The mean and maximum error are estimated
Err_Temp_Med_mFGN<- mean(E_corr_Temp_mFGN)
Err_Temp_Max_mFGN<- max(E_corr_Temp_mFGN)

# The spatial error matrix is initialized
E_corr_Spa_mFGN <-array(0,c(E,E,12)) 
for (j in 1:12){
  # The dimension of the month “j” is extracted:
  Q_aux <- Q_3D_N[,j,]
  # The correlation matrix is estimated
  QQ_Obs <- cor(Q_aux)
  # The same procedure is repeated for the mFGN
  Q_aux <- mFGN[,j,]
  QQ_mFGN <- cor(Q_aux)
  
  # The absolute difference are estimated (between observed and generated data)
  E_corr_Spa_mFGN[,,j]<-abs(QQ_Obs-QQ_mFGN)
}
# The mean and maximum error are estimated
Err_Spa_Med_mFGN<- mean(E_corr_Spa_mFGN)
Err_Spa_Max_mFGN<- max(E_corr_Spa_mFGN)


# The Lambda values are generated
Lambda<- Alpha
# The variables A1 and B1, are generated as the minimum and its position
# of the optimum function 1 (of the mean error). The position (B1) allows 
# indicating which Alpha value
# is related to each Lambda
A1<- array(0,c(length(Lambda),1))
B1<- array(0,c(length(Lambda),1))
# A2 and B2 are the same as A1 and B2, but for the objective function 2,
# hence, they work with the maximum errors.
A2<- array(0,c(length(Lambda),1))
B2<- array(0,c(length(Lambda),1))
# The errors associated to mFGN are stored in the following variables
Obj1_mFGN<- array(0,c(length(Lambda),1))
Obj2_mFGN<- array(0,c(length(Lambda),1))

# A “for” loop is generated to go through Lambda values
for(k in 1:length(Lambda)){
  # Objective function 1 values are estimated. Note that the variables 
  # Err_Temp_Med and Err_Spa_Med already consider all the possible Alpha
  # values, hence we have to find the best one
  Obj1<- Lambda[k]*Err_Temp_Med+(1-Lambda[k])*Err_Spa_Med
  # The minimum values are found
  A1[k]<- min(Obj1)
  # The values of Alpha, related to the minimum A1[k] are found. This are
  # stored in B1[k]. Hence, Alpha[B1[k]] is the optimal parameter associated
  # to Lambda [k]
  B1[k]<- which(Obj1==A1[k])
  
  # The same procedure is performed for the objective function 2 (maximums):
  Obj2<- Lambda[k]*Err_Temp_Max+(1-Lambda[k])*Err_Spa_Max
  A2[k]<- min(Obj2)
  B2[k]<- which(Obj2==A2[k])
  
  # Similar procedure than that of Objective Function 1 and 2, but for 
  # the mFGN method
  Obj1_mFGN[k]<- Lambda[k]*Err_Temp_Med_mFGN+(1-Lambda[k])*Err_Spa_Med_mFGN
  Obj2_mFGN[k]<- Lambda[k]*Err_Temp_Max_mFGN+(1-Lambda[k])*Err_Spa_Max_mFGN      
}

# Several plots of results:
# Optimum Alpha (function 1), against Lambda
plot(Lambda,Alpha[B1],type="l",xlim=c(0,1),ylim=c(0,1),col=rgb(0.25,0.85,0.2),lwd=2)
legend(x=0.6,y=0.2,legend = "Mean Error",lt=1,col=rgb(0.25,0.85,0.2),lwd=2)

# Optimum Alpha (function 2), against Lambda
plot(Lambda,Alpha[B2],type="l",xlim=c(0,1),ylim=c(0,1),col=rgb(0.75,0.80,0.5),lwd=2)
legend(x=0.1,y=0.8,legend = "Max Error",lt=1,col=rgb(0.75,0.80,0.5),lwd=2)


# Objective function 1 is plotted against Lambda
plot(Lambda,A1,type="l",xlim=c(0,1),ylim=c(0,0.08),ylab = "Objective Function 1",col=rgb(0.25,0.35,0.85),lwd=2)

matplot(Lambda,Obj1_mFGN,type="l",lt=2,add=T,col=rgb(0.85,0.25,0.25),lwd=2)
legend(x=0.5,y=0.08,legend = c("WmFGN Mean Error","mFGN Mean Error"),lt=c(1,2),col=c(rgb(0.25,0.35,0.85),rgb(0.85,0.25,0.25)),lwd=2)


# Objective function 2 is plotted against Lambda
plot(Lambda,A2,type="l",xlim=c(0,1),ylim=c(0,0.35),ylab = "Objective Function 2",col=rgb(0.25,0.35,0.85),lwd=2)

matplot(Lambda,Obj2_mFGN,type="l",lt=2,add=T,col=rgb(0.85,0.25,0.25),lwd=2)
legend(x=0.5,y=0.08,legend = c("WmFGN Max Error","mFGN Max Error"),lt=c(1,2),col=c(rgb(0.25,0.35,0.85),rgb(0.85,0.25,0.25)),lwd=2)

# Mean error is plotted against Alpha
plot(Alpha,Err_Spa_Med,type="l",xlim=c(0,1),ylim=c(0,0.13),ylab = "Mean Error",col=rgb(0.25,0.35,0.85),lwd=2)
matplot(Alpha,Err_Temp_Med,type="l",lt=2,add=T,col=rgb(0.25,0.35,0.85),lwd=2)

matplot(Alpha,array(1,c(length(Alpha),1))%*%Err_Spa_Med_mFGN,type="l",add=T,col=rgb(0.85,0.25,0.25),lwd=2)
matplot(Alpha,array(1,c(length(Alpha),1))%*%Err_Temp_Med_mFGN,type="l",lt=2,add=T,col=rgb(0.85,0.25,0.25),lwd=2)
legend(x=0.5,y=0.13,legend = c("WmFGN Temp Err","WmFGN Spa Err","mFGN Temp Err","mFGN Spa Err"),lt=c(2,1,2,1),col=c(rgb(0.25,0.35,0.85),rgb(0.25,0.35,0.85),rgb(0.85,0.25,0.25),rgb(0.85,0.25,0.25)),lwd=2)


# Maximum error is plotted against Alpha
plot(Alpha,Err_Spa_Max,type="l",xlim=c(0,1),ylim=c(0,0.35),ylab = "Mean Error",col=rgb(0.25,0.35,0.85),lwd=2)
matplot(Alpha,Err_Temp_Max,type="l",lt=2,add=T,col=rgb(0.25,0.35,0.85),lwd=2)

matplot(Alpha,array(1,c(length(Alpha),1))%*%Err_Spa_Max_mFGN,type="l",add=T,col=rgb(0.85,0.25,0.25),lwd=2)
matplot(Alpha,array(1,c(length(Alpha),1))%*%Err_Temp_Max_mFGN,type="l",lt=2,add=T,col=rgb(0.85,0.25,0.25),lwd=2)
legend(x=0.5,y=0.13,legend = c("WmFGN Temp Err","WmFGN Spa Err","mFGN Temp Err","mFGN Spa Err"),lt=c(2,1,2,1),col=c(rgb(0.25,0.35,0.85),rgb(0.25,0.35,0.85),rgb(0.85,0.25,0.25),rgb(0.85,0.25,0.25)),lwd=2)

# Pareto Curve of the mean errors 
plot(Err_Spa_Med,Err_Temp_Med,type="l",lty=4,xlim=c(0,0.08),ylim=c(0,0.08),ylab = "Mean Temporal Error",xlab = "Mean Spatial Error",col=rgb(0.25,0.35,0.85),lwd=2)
matplot(Err_Spa_Med_mFGN,Err_Temp_Med_mFGN,pch=1,add=T,col=rgb(0.85,0.25,0.25),lwd=2)
legend(x=0.04,y=0.08,legend = c("WmFGN","mFGN"),lt=c(4,NA),pch=c(NA,1),col=c(rgb(0.25,0.35,0.85),rgb(0.85,0.25,0.25)),lwd=2)

# Pareto Curve of the maximum errors 
plot(Err_Spa_Max,Err_Temp_Max,type="l",lty=4,xlim=c(0,0.35),ylim=c(0,0.35),ylab = "Max Temporal Error",xlab = "Max Spatial Error",col=rgb(0.25,0.35,0.85),lwd=2)
matplot(Err_Spa_Max_mFGN,Err_Temp_Max_mFGN,pch=1,add=T,col=rgb(0.85,0.25,0.25),lwd=2)
legend(x=0.01,y=0.1,legend = c("WmFGN","mFGN"),lt=c(4,NA),pch=c(NA,1),col=c(rgb(0.25,0.35,0.85),rgb(0.85,0.25,0.25)),lwd=2)


# One selects a Lambda value:
lamb_sel <- 0.5
# The Alpha associated to that Lambda (using the mean error function 1),
# is found.
Alph_Sel<- Alpha[B1[which(Lambda==lamb_sel)]]

# The WmFGN is estimated (associated to the selected Lambda value)
WmFGN <- (1-Alph_Sel)*mFGNS + Alph_Sel*SmFGN

# One can transform the data back to streamflows by:
New_Q_3D <- array(0,dim(WmFGN))
for (i in 1:ng){
  New_Q_3D[i,,] <- exp(WmFGN[i,,]*sig_y+mu_y)
}

# Verifying that the mean and the standard deviation are properly simulated 
N_mu_x <- apply(New_Q_3D,2:3,mean)
N_mu_x-mu_x


N_sig_x <- apply(New_Q_3D,2:3,sd)
N_sig_x-sig_x

Mes <- 1:12

# Select a station from 1 to 9 to plot:
Est <- 2
# Plot the mean values of the station 
plot(Mes,mu_x[,Est],type = 'l',xlab = 'Months of a Hydrological Year',ylab = 'Mean Q [m^3/s]',ylim=c(0, max(mu_x[,Est])*1.5),col=rgb(0.27,0.37,0.82),lwd=2)
matplot(Mes,N_mu_x[,Est],type = 'l',add=T,col=rgb(0.77,0.37,0.22),lty = 2,lwd=2.5)
legend(x=9,y=max(mu_x[,Est])*1.45,legend = c('Obs','WmFGN'),col=c(rgb(0.27,0.37,0.82),rgb(0.77,0.37,0.22)),lwd=c(2,2.5),lty = c(1,2))


# Plot the standard deviation values of the station 
plot(Mes,sig_x[,Est],type = 'l',xlab = 'Months of a Hydrological Year',ylab = 'Std. Q [m^3/s]',ylim=c(0, max(sig_x[,Est])*1.5),col=rgb(0.27,0.37,0.82),lwd=2)
matplot(Mes,N_sig_x[,Est],type = 'l',add=T,col=rgb(0.77,0.37,0.22),lty = 2,lwd=2.5)
legend(x=9,y=max(sig_x[,Est])*1.45,legend = c('Obs','WmFGN'),col=c(rgb(0.27,0.37,0.82),rgb(0.77,0.37,0.22)),lwd=c(2,2.5),lty = c(1,2))



