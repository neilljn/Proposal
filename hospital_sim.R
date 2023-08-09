### now trying to simulate hospital data: SI model, discrete time, patients entering and exiting, some colonised on admission, inaccurate testing

# function to do the simulation
hospital_sim <- function(beta,col_prop,sens,spec,admit_rate,dis_rate,test_rate,T_max,seed){
  
  # list of inputs:
  # beta - transmission rate
  # col_prop - proportion of being colonised on admission
  # sens - test sensitivity
  # spec - test specificity
  # admit_rate - rate of admittance to the hospital (overall)
  # dis_rate - rate of discharge from the hospital (per person)
  # test_rate - rate of testing per patient
  # T_max - length of time to run the simulation for
  # seed - value of the random seed
  
  # setting the seed
  set.seed(seed)
  
  # creating empty variables for the current susceptible and infectious status of each patient
  S <- c()
  I <- c()
  
  # creating empty variables for storage
  S_list <- list()
  I_list <- list()
  sum_S <- c()
  sum_I <- c()
  sum_I_admit <- rep(0,times=T_max)
  admit_times <- c()
  dis_times <- c()
  col_times <- c()
  col_admit <- c()
  test_times <- c()
  test_ids <- c()
  test_results <- c()
  I_star <- c()
  
  # initial number of patients
  n <- 0
  
  # initial number of tests
  m <- 0
  
  # simulating for T_max time
  for(t in 1:T_max){
    
    # discharging patients from the previous day
    if(n!=0){
      for(i in 1:n){
        dis<- rbinom(1,S[i]+I[i],1-exp(-1*dis_rate)) # determining if a discharge would happen
        if(dis==1){
          dis_times[i] <- t-1 # recording the time of discharge as the previous day
          S[i] <- 0 
          I[i] <- 0
        }
      }
    }
    
    # generating the number of arrivals from a poisson distribution
    n_new <- rpois(1,admit_rate)
    
    # for each new arrival, assign either susceptible or infective
    if(n_new!=0){
      for(i in (n+1):(n+n_new)){
        admit_times[i] <- t #saving the admission time for this patient
        col_arrival <- rbinom(1,1,col_prop) #determining if the new patient is colonised on admission
        if(col_arrival==1){
          S[i] <- 0
          I[i] <- 1
          col_times[i] <- Inf # this individual cannot be colonised in the hospital, since they are colonised on admission
          col_admit[i] <- 1 # recording being colonised on admission for that individual
          sum_I_admit[t] <- sum_I_admit[t] + 1 # recording being colonised on admission for that time
        }
        else{
          S[i] <- 1
          I[i] <- 0
          col_admit[i] <- 0 # recording not being colonised on admission
        }
      }
    }
    
    # updating the total number of patients ever admitted
    n <- n + n_new
    
    # tracking the number of infectives when transmission occurs
    I_star[t] <- sum(I)
    
    # creating a vector for transitions
    if(n!=0){
      SI <- rep(0,times=n)
    }
    
    # simulating each individual's transmission
    if(n!=0){
      for(i in 1:n){
        
        # simulating transitions for each individual
        SI[i] <- rbinom(1,S[i],1-exp(-1*beta*sum(I)))
        
        # storing times of transitions
        if(SI[i]==1){
          col_times[i] <- t
        }
        
      }
    }
    
    # applying the transitions
    if(n!=0){
      S <- S - SI
      I <- I + SI
    }
    
    # storing values of S, I, and R
    S_list[[t]] <- S
    I_list[[t]] <- I
    sum_S[t] <- sum(S)
    sum_I[t] <- sum(I)
    
    # possibly testing on each individual currently admitted
    if(n!=0){
      for(i in 1:n){
        test <- rbinom(1,S[i]+I[i],1-exp(-1*test_rate)) # determining if a test will happen on this patient
        if(test==1){
          m <- m+1 # increase number of tests
          test_times[m] <- t  # recording time of current test
          test_ids[m] <- i # recording patient of current test
          prob_pos <- (sens*I[i])+((1-spec)*S[i]) # probability of positive test result
          test_results[m] <- rbinom(1,1,prob_pos) # recording outcome of current test
        }
      }
    }
    
  }
  
  # if colonisation time is NA, turn it into infinity
  if(n!=0){
    for(i in 1:n){
      if(is.na(col_times[i])){
        col_times[i] <- Inf
      }
    }
  }
  
  # if discharge time is NA, turn it into T_max
  if(n!=0){
    for(i in 1:n){
      if(is.na(dis_times[i])){
        dis_times[i] <- T_max
      }
    }
  }
  
  # completing colonisation and discharge time vectors
  n_col <- length(col_times)
  n_dis <- length(dis_times)
  if(n_col != n){
    for(i in (n_col+1):n){
      col_times[i] <- Inf
    }
  }
  if(n_dis != n){
    for(i in (n_dis+1):n){
      dis_times[i] <- T_max
    }
  }
  
  # returning outputs
  return(list(S=S_list,I=I_list,sum_S=sum_S,sum_I=sum_I,sum_I_admit=sum_I_admit,I_star=I_star,admit_times=admit_times,dis_times=dis_times,col_times=col_times,col_admit=col_admit,test_times=test_times,test_ids=test_ids,test_results=test_results,patient_total=n,test_total=m))
  
  # list of outputs:
  # S - list of susceptible status for each individual over time
  # I - list of infective status for each individual over time
  # sum_S - total number of susceptible individuals over time
  # sum_I - total number of infective individuals over time
  # sum_I_admit - total number of infective individuals on admission over time
  # I_star - total number of infective individuals over time at time of infection
  # admit_times - the time of admission for each individual
  # dis_times - the time of discharge for each individual
  # col_times - the time of colonisation for each individual
  # col_admit - whether that individual is colonised on admission
  # test_times - the time of each test
  # test_ids - the individual number associated with this test
  # test_results - the result of this test
  # patient_total - number of patients admitted
  # test_total - number of tests carried out
  
}

# trying out the simulation
data_hospital <- hospital_sim(1,0.1,0.95,0.95,1,0.5,1,100,1)

# plotting the results of the simulation (base plot)
col_admit_count <- data_hospital$sum_I_admit
col_admit_count[col_admit_count == 0] <- NA
plot(1:100,data_hospital$sum_S,type="l",ylim=c(0,max(c(data_hospital$sum_S,data_hospital$sum_I))),xlab="Day",ylab="Number of Individuals",main="Hospital Model"); lines(1:100,data_hospital$sum_I,col="red"); points(col_admit_count,pch=19,col="darkred",lwd=2.5); legend("topright",legend=c("S", "I", "I (admission)"),col=c("black","red","darkred"), lty= c(1,1,NA), pch=c(NA,NA,19), lwd=c(1,1,2.5))

# plotting the results of the simulation (ggplot)
plot_S <- as.data.frame(cbind(1:100,data_hospital$sum_S))
plot_I <- as.data.frame(cbind(1:100,data_hospital$sum_I))
plot_I_admit <- as.data.frame(cbind(1:100,col_admit_count))
sim_plot <- ggplot() + geom_line(aes(x=V1,y=V2,col),plot_S) + geom_line(aes(x=V1,y=V2),plot_I,col="red") + geom_point(aes(x=V1,y=col_admit_count),plot_I_admit,col="darkred",size=2.25)
sim_plot <- sim_plot + xlab("Day") + ylab("Individuals") + ggtitle("Hospital Simulation") + theme_classic()
print(sim_plot)

# generating data for the MCMC
data_MRSA <- hospital_sim(0.1,0.05,0.8,1,3,1/3,0.5,500,1)

# finding the number of true postives and false negatives
TP <- 0
FN <- 0
for(j in 1:(data_MRSA$test_total)){ # checking over every test
  if(data_MRSA$test_times[j] >= data_MRSA$col_times[data_MRSA$test_ids[j]]){ # if the time of the test is after the colonisation time, the patient is colonised
    if(data_MRSA$test_results[j]==T){ # if the test is positive, we have a true positive
      TP <- TP+1
    }
    else{ # if the test is negative, we have a false negative
      FN <- FN+1
    }
  }
  if((data_MRSA$test_times[j] >= data_MRSA$admit_times[data_MRSA$test_ids[j]])&(data_MRSA$col_admit[data_MRSA$test_ids[j]]==1)){ # if the time of the test is after the colonisation time, the patient is colonised
    if(data_MRSA$test_results[j]==T){ # if the test is positive, we have a true positive
      TP <- TP+1
    }
    else{ # if the test is negative, we have a false negative
      FN <- FN+1
    }
  }
}
print(c(TP,FN)) # number of true positives and false negatives

# summary stats
data_MRSA$patient_total # total number of patients
mean(data_MRSA$dis_times-data_MRSA$admit_times) #average time in hospital
sum(is.finite(data_MRSA$col_times)) # patients infected in the hospital
sum(data_MRSA$col_admit) # patients infected on admission
data_MRSA$patient_total - sum(is.finite(data_MRSA$col_times)) - sum(data_MRSA$col_admit) # patients never admitted
data_MRSA$test_total # total number of tests
data_MRSA$test_total - TP - FN # number of true negatives