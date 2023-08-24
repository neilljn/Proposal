### coding up an MCMC algorithm myself (from Theodore Kypraios et al 2013)

# libraries
library(beepr)
library(pheatmap)
library(ggplot2)
library(coda)

# function for how much to scale lambda by
scale <- function(k){
  return(1-1.001^(-(k+2000)))
}

# function for the conditional log likelihood of beta
log_like_beta <- function(beta,t_admit,t_dis,t_col,phi,c,n,sum_I){
  x <- rep(0,times=n)
  for(i in 1:n){ # summing over each patient
    if(phi[i]==0){ # non-zero only if susceptible on admission
      if(t_col[i]==Inf){
        for(t in t_admit[i]:t_dis[i]){
          x[i] <- x[i] - beta*sum_I[t] # if never colonised during stay, this is the likelihood of being never colonised
        }
      }
      else{
        if(t_admit[i]!=t_col[i]){
          for(t in t_admit[i]:(t_col[i]-1)){
            x[i] <- x[i] - beta*sum_I[t] # the likelihood of not being colonised, up to the point where there were
          }
        }
        x[i] <- x[i] + log(1-exp(-1*beta*sum_I[t_col[i]])) # the likelihood of being colonised at that point in time
        if(is.nan(log(1-exp(-1*beta*sum_I[t_col[i]])))==1){
          print(c(i,t_col[i],sum_I[t_col[i]],t_admit[i],t_dis[i]))
        }
      }
    }
  }
  return(sum(x)-(c*beta))
}

# function for the log likelihood
log_like <- function(t_admit,t_dis,t_test,id_test,res_test,z,p,beta,t_col,phi,a_z,b_z,a_p,b_p,c,n,TP,FN,sum_I){
  term_z <- ((a_z+TP-1)*log(z)) + ((b_z+FN-1)*log(1-z))
  term_p <- ((a_p-1+sum(phi))*log(p)) + ((b_p - 1 + n - sum(phi))*log(1-p))
  term_beta <- log_like_beta(beta,t_admit,t_dis,t_col,phi,c,n,sum_I)
  return(term_z+term_p+term_beta)
}

# function to run the MCMC
MRSA_MCMC <- function(t_admit,t_dis,t_test,id_test,res_test,z_0,p_0,beta_0,t_col_0,phi_0,lambda_0,a_z,b_z,a_p,b_p,c,w,M,M_repeat,seed,end_noise,print_it){
  
  # list of inputs:
  # t_admit - vector of admission times
  # t_dis - vector of discharge times
  # t_test - vector of test times
  # id_test - the patient id for each test
  # res_test - the result of each test
  # z_0 - initial value of z
  # p_0 - initial value of p
  # beta_0 - initial value of lambda
  # t_col_0 - initial vector of colonisation times
  # phi_0 - initial vector of admission colonisation states
  # lambda_0 - initial tuning parameter for the RWM for beta
  # a_z - hyperparameter for the prior of z
  # b_z - hyperparameter for the prior of z
  # a_p - hyperparameter for the prior of p
  # b_p - hyperparameter for the prior of p
  # c - hyperparameter for the prior of beta
  # w - some sort of efficiency value? (between 0 and 1)
  # M - number of MCMC iterations
  # M_repeat - number of updates to the latent variables within each MCMC iteration
  # seed - value of the random seed
  # end_noise - whether a noise is made at the end
  # print_it - whether the iteration number is printed
  
  # setting the seed
  set.seed(seed)
  
  # total number of patients
  n <- length(t_admit)
  
  # total number of tests
  n_tests <- length(t_test)
  
  # end time - note that time is indexed from 1
  T_max <- max(data_MRSA$dis_times)
  
  # making current parameter, latent variable, and tuning parameter values equal to the initial values
  z <- z_0
  p <- p_0
  beta <- beta_0
  t_col <- t_col_0
  phi <- phi_0
  lambda <- lambda_0
  
  # creating empty vectors and matrices for storage
  z_vec <- rep(0,times=M+1)
  p_vec <- rep(0,times=M+1)
  beta_vec <- rep(0,times=M+1)
  t_col_matrix <- matrix(0,nrow=M+1,ncol=n)
  phi_matrix <- matrix(0,nrow=M+1,ncol=n)
  sum_I_matrix <- matrix(0,nrow=M+1,ncol=T_max)
  lambda_vec <- rep(0,times=M+1)
  
  # adding the initial values to the storage vectors/matrices
  z_vec[1] <- z
  p_vec[1] <- p
  beta_vec[1] <- beta
  t_col_matrix[1,] <- t_col
  phi_matrix[1,] <- phi
  lambda_vec[1] <- lambda
  
  # initial number of accepted RWM steps for updating beta
  acc_beta <- 0
  
  # initial number of accepted latent variable moves
  acc_latent <- rep(0,times=12)
  M_latent <- rep(0,times=12)
  
  # finding the time of each patients first positive test result
  t_first_pos <- rep(Inf,times=n)
  for(k in 1:n_tests){ # checking over each test
    if((res_test[k] == 1) & (t_first_pos[id_test[k]] == Inf)){ # checking if the test is positive and that individual has no previous positive test
      t_first_pos[id_test[k]] <- t_test[k] # setting the time of the first positive test to the time of the current test
    }
  }
  
  # determining if a patient never has any positive test results
  always_negative <- 1-is.finite(t_first_pos)
  
  # calculating whether each individual ever gets MRSA using the initial latent variables
  col_ever <- rep(0,times=n)
  for(i in 1:n){
    if(phi[i]==1){
      col_ever[i] <- 1
    }
    if(t_col[i]<Inf){
      col_ever[i] <- 1
    }
  }
  col_ever_total <- sum(col_ever)
  
  # calculating the initial number of infected individuals at each point in time
  sum_I <- rep(0,times=T_max)
  for(j in 1:T_max){ # checking over every point in time
    for(i in 1:n){ # checking over every patient
      if((j > t_col[i]) & (j <= t_dis[i])){ # if the time is between colonisation and discharge, there is an infective patient
        sum_I[j] <- sum_I[j]+1
      }
      if((j >= t_admit[i]) & (j <= t_dis[i]) & (phi[i]==1)){ # if the time is between admission and discharge, and they were infected on admission, there is an infective patient
        sum_I[j] <- sum_I[j]+1
      }
    }
  }
  sum_I_matrix[1,] <- sum_I
  
  # calculating the initial number of true positives and false negatives
  TP <- 0
  FN <- 0
  for(j in 1:n_tests){ # checking over every test
    if(t_test[j] >= t_col[id_test[j]]){ # if the time of the test is after the colonisation time, the patient is colonised
      if(res_test[j]==T){ # if the test is positive, we have a true positive
        TP <- TP+1
      }
      else{ # if the test is negative, we have a false negative
        FN <- FN+1
      }
    }
    if((t_test[j] >= t_admit[id_test[j]])&(phi[id_test[j]]==1)){ # if the time of the test is after the colonisation time, the patient is colonised
      if(res_test[j]==T){ # if the test is positive, we have a true positive
        TP <- TP+1
      }
      else{ # if the test is negative, we have a false negative
        FN <- FN+1
      }
    }
  }
  
  # running the macroreplications of the MCMC
  for(k in 1:M){
    
    # printing the current iteration number
    if(print_it==1){
      print(k)
    }
    
    # Gibbs step to update z
    z <- rbeta(1,a_z+TP,b_z+FN)
    z_vec[k+1] <- z
    
    # Gibbs step to update p
    p <- rbeta(1,a_p+sum(phi),b_p+n-sum(phi))
    p_vec[k+1] <- p
    
    # RWM step to update beta
    beta_new <- rnorm(1,beta,lambda) # generating new beta from a random walk
    if(beta_new>0){ # there will only ever be an acceptance if beta is positive
      log_alpha <- min(0,log_like_beta(beta_new,t_admit,t_dis,t_col,phi,c,n,sum_I)-log_like_beta(beta,t_admit,t_dis,t_col,phi,c,n,sum_I)) # log of the acceptance probability
      u <- runif(1)
      if(log(u) < log_alpha){ #accepting with probability alpha
        beta <- beta_new
        lambda <- lambda*(scale(k)^(-2))
        acc_beta <- acc_beta + 1
      }
      else{
        lambda <- lambda*scale(k)
      }
    }
    else{
      lambda <- lambda*scale(k)
    }
    beta_vec[k+1] <- beta
    lambda_vec[k+1] <- lambda
    
    # updating t_col and phi, M_repeat times
    for(k_latent in 1:M_repeat){

      # updating some varibles
      col_ever <- phi | is.finite(t_col)
      col_ever_prime <- col_ever & always_negative
      m <- sum(col_ever)
      m_prime <- sum(col_ever_prime)

      # choosing which M-H step to update t_col and phi
      if(m==0){
        choice <- 3
      }
      else if((m==n)&(m_prime==0)){
        choice <- 1
      }
      else if((m==n)&(m_prime!=0)){
        choice <- sample(1:2,1)
      }
      else if((m!=n)&(m_prime==0)){
        choice <- sample(c(1,3),1)
      }
      else{
        choice <- sample(1:3,1)
      }
      
      # creating variables for the proposed values of the latent variables and sum_I
      t_col_new <- t_col
      phi_new <- phi
      sum_I_new <- sum_I

      # the probability ratio is reset to 0 in case of errors
      q_move <- 0

      # option 1: moving a colonisation time
      if(choice==1){

        # choosing a patient who at some point has MRSA
        patients_MRSA <- which(col_ever %in% 1) # the index of all patients who at some point have MRSA
        if(length(patients_MRSA)==1){
          patient_move <- patients_MRSA[1]
        }
        else{
          patient_move <- sample(patients_MRSA,1) # randomly choosing one of these patients
        }

        # the latest time they can be infected
        l <- min(t_dis[patient_move],t_first_pos[patient_move])

        # with probability w, they are colonised on admission
        u <- runif(1)
        if(u<w){

          # case where they are already colonised on admission
          if(phi[patient_move]==1){
            case <- 5
            q_move <- 1
          }

          # case where they were previously colonised during the hospital stay
          if(phi[patient_move]==0){
            case <- 6
            phi_new[patient_move] <- 1 # now colonised on admission
            t_col_new[patient_move] <- Inf # since colonised on admission, no colonisation time
            q_move <- (1-w)/(w*(l - t_admit[patient_move] + 1))
            for(t in t_admit[patient_move]:t_col[patient_move]){  # updating sum I
              sum_I_new[t] <- sum_I_new[t] + 1
            }
          }

        }

        # otherwise, they have a colonisation time uniformly before their first positive test (or discharge)
        if(u>=w){

          # case where they were previously colonised on admission
          if(phi[patient_move]==1){
            case <- 7
            phi_new[patient_move] <- 0 # no longer colonised on admission
            if(t_admit[patient_move]==l){
              t_col_new[patient_move] <- l # only possible new colonisation time
            }
            else{
              t_col_new[patient_move] <- sample(t_admit[patient_move]:l,1) # random new colonisation time
            }
            q_move <- (w*(l - t_admit[patient_move] + 1))/(1-w)
            for(t in t_admit[patient_move]:t_col_new[patient_move]){ # updating sum I
              sum_I_new[t] <- sum_I_new[t] - 1
            }
          }

          # case where they were already colonised during the hospital stay
          if(phi[patient_move]==0){
            case <- 8
            if(t_admit[patient_move]==l){
              t_col_new[patient_move] <- l # only possible new colonisation time
            }
            else{
              t_col_new[patient_move] <- sample(t_admit[patient_move]:l,1) # random new colonisation time
            }
            q_move <- 1
            if(t_col_new[patient_move]>t_col[patient_move]){ # updating sum I if new colonisation time is later
              for(t in (t_col[patient_move]+1):t_col_new[patient_move]){
                sum_I_new[t] <- sum_I_new[t] - 1
              }
            }
            if(t_col_new[patient_move]<t_col[patient_move]){ # updating sum I if new colonisation time is earlier
              for(t in (t_col_new[patient_move]+1):t_col[patient_move]){
                sum_I_new[t] <- sum_I_new[t] + 1
              }
            }
          }
        }

      }

      # option 2: removing a colonisation time
      if(choice==2){

        # choosing a patient who at some point has MRSA but never tested positive
        patients_MRSA <- which(col_ever_prime %in% 1) # the index of all patients who at some point have MRSA
        if(length(patients_MRSA)==1){
          patient_remove <- patients_MRSA[1]
        }
        else{
          patient_remove <- sample(patients_MRSA,1) # randomly choosing one of these patients
        }

        # the latest time they can be infected
        l <- min(t_dis[patient_remove],t_first_pos[patient_remove])

        # if they were colonised on admission, we make them susceptible on admission (and never colonised)
        if(phi[patient_remove]==1){
          case <- 9
          phi_new[patient_remove] <- 0 # no longer colonised on admission
          q_move <- (w*m_prime)/(n - m + 1)
          for(t in t_admit[patient_remove]:t_dis[patient_remove]){ # updating sum_I
            sum_I_new[t] <- sum_I_new[t] - 1
          }
        }

        # if they were colonised after admission, we remove their colonisation time
        if(phi[patient_remove]==0){
          case <- 10
          t_col_new[patient_remove] <- Inf # no longer colonised at any point
          q_move <- ((1-w)*m_prime)/((l - t_admit[patient_remove] + 1)*(n - m + 1))
          if(t_col[patient_remove]!=t_dis[patient_remove]){
            for(t in (t_col[patient_remove]+1):t_dis[patient_remove]){ # updating sum_I
              sum_I_new[t] <- sum_I_new[t] - 1
            }
          }
        }

      }

      # option 3: adding a colonisation time
      if(choice==3){

        # choosing a patient who never has MRSA
        patients_no_MRSA <- which(col_ever %in% 0) # the index of all patients who never have MRSA
        if(length(patients_no_MRSA)==1){
          patient_add <- patients_no_MRSA[1]
        }
        else{
          patient_add <- sample(patients_no_MRSA,1) # randomly choosing one of these patients
        }

        # the latest time they can be infected
        l <- min(t_dis[patient_add],t_first_pos[patient_add])

        # with probability w, they become colonised on admission
        u <- runif(1)
        if(u<w){
          case <- 11
          phi_new[patient_add] <- 1 # now colonised on admission
          q_move <- (n - m)/(w*(m_prime + 1))
          for(t in t_admit[patient_add]:t_dis[patient_add]){ # updating sum_I
            sum_I_new[t] <- sum_I_new[t] + 1
          }
        }

        # otherwise, they now have a colonisation time uniformly before their first positive test (or discharge)
        if(u>=w){
          case <- 12
          if(t_admit[patient_add]==l){
            t_col_new[patient_add] <- l # only possible new colonisation time
          }
          else{
            t_col_new[patient_add] <- sample(t_admit[patient_add]:l,1) # random new colonisation time
          }
          q_move <- ((n - m)*(l - t_admit[patient_add] + 1))/((1-w)*(m_prime + 1))
          if(t_col_new[patient_add]!=t_dis[patient_add]){
            for(t in (t_col_new[patient_add]+1):t_dis[patient_add]){ # updating sum_I
              sum_I_new[t] <- sum_I_new[t] + 1
            }
          }
        }

      }

      # calculating the number of true positives and false negatives under the proposed latent variables
      TP_new <- 0
      FN_new <- 0
      for(j in 1:n_tests){ # checking over every test
        if(t_test[j] >= t_col_new[id_test[j]]){ # if the time of the test is after the colonisation time, the patient is colonised
          if(res_test[j]==T){ # if the test is positive, we have a true positive
            TP_new <- TP_new+1
          }
          else{ # if the test is negative, we have a false negative
            FN_new <- FN_new+1
          }
        }
        if((t_test[j] >= t_admit[id_test[j]])&(phi_new[id_test[j]]==1)){ # if the time of the test is after the colonisation time, the patient is colonised
          if(res_test[j]==T){ # if the test is positive, we have a true positive
            TP_new <- TP_new+1
          }
          else{ # if the test is negative, we have a false negative
            FN_new <- FN_new+1
          }
        }
      }

      # calculating the log-likelihood under both the current and proposed latent variables
      ll_current <- log_like(t_admit,t_dis,t_test,id_test,res_test,z,p,beta,t_col,phi,a_z,b_z,a_p,b_p,c,n,TP,FN,sum_I)
      ll_new <- log_like(t_admit,t_dis,t_test,id_test,res_test,z,p,beta,t_col_new,phi_new,a_z,b_z,a_p,b_p,c,n,TP_new,FN_new,sum_I_new)

      # performing the M-H step
      log_alpha <- min(0,log(q_move)+ll_new-ll_current)
      u_MH <- runif(1)
      if(log(u_MH)<log_alpha){ # accepting the move with probability alpha
        t_col <- t_col_new # updating latent variables
        phi <- phi_new
        sum_I <- sum_I_new # updating dependents of latent variables
        TP <- TP_new
        FN <- FN_new
        acc_latent[1] <- acc_latent[1] + 1 # updating acceptance trackers
        acc_latent[choice+1] <- acc_latent[choice+1]+ 1
        acc_latent[case] <- acc_latent[case] + 1
      }

      # updating the number of steps of each type
      M_latent[1] <- M_latent[1] + 1
      M_latent[choice+1] <- M_latent[choice+1]+ 1
      M_latent[case] <- M_latent[case] + 1

    }

    # storing latent variables
    t_col_matrix[k+1,] <- t_col
    phi_matrix[k+1,] <- phi
    sum_I_matrix[k+1,] <- sum_I
    
  }
  
  # making a noise since the code has finished (if enabled)
  if(end_noise==1){
    beep(5)
  }
  
  # returning outputs
  return(list(z=z_vec,p=p_vec,beta=beta_vec,t_col=t_col_matrix,phi=phi_matrix,sum_I=sum_I_matrix,lambda=lambda_vec,accept_beta=acc_beta/M,accept_latent=acc_latent/M_latent))
  
  # list of outputs
  # z - vector of generated values for the sensitivity
  # p - vector of generated values for the proportion of individuals who are colonised on admission
  # beta - vector of generated values for the transmission rate
  # t_col - matrix of generated vectors for colonisation times
  # phi - matrix of generated vectors for being colonised on admission
  # sum_I - matrix of generate vectors for the number of individuals infected at each point in time
  # lambda - vector of the changing values of the tuning parameter
  # accept_beta - acceptance rate for the RWM step (updating beta)
  # accept_latent - acceptance rate for the latent variable step 
  
}

# testing the MCMC (parameters only)
M <- 10000
MRSA_inference <- MRSA_MCMC(data_MRSA$admit_times,data_MRSA$dis_times,data_MRSA$test_times,data_MRSA$test_ids,data_MRSA$test_results,0.8,0.05,0.1,data_MRSA$col_times,data_MRSA$col_admit,0.01,1,1,1,1,0.001,0.3,M,50,1,0,1)

# looking at the acceptance rate of beta
MRSA_inference$accept_beta

# looking at the acceptance rate of the latent variables
MRSA_inference$accept_latent[1] # all changes
MRSA_inference$accept_latent[2:4] # move / remove / add
MRSA_inference$accept_latent[5:8] # subcases of move
MRSA_inference$accept_latent[9:10] # subcases of remove
MRSA_inference$accept_latent[11:12] # subcases of add

# trace plots for the parameters
plot(1:(M+1),MRSA_inference$z,type="l",main="z"); abline(h=0.8,col="red")
plot(1:(M+1),MRSA_inference$p,type="l",main="p"); abline(h=0.05,col="red") 
plot(1:(M+1),MRSA_inference$beta,type="l",main=expression(beta));  abline(h=0.1,col="red")

# parameter averages
mean(MRSA_inference$z)
mean(MRSA_inference$p)
mean(MRSA_inference$beta)

# relative difference between true value and MCMC average
(0.8-mean(MRSA_inference$z))/0.8
(0.05-mean(MRSA_inference$p))/0.1
(0.1-mean(MRSA_inference$beta))/0.1

# effective sample sizes
effectiveSize(MRSA_inference$z)
effectiveSize(MRSA_inference$p)
effectiveSize(MRSA_inference$beta)

# trace plots for the latent variables
col_times <- MRSA_inference$t_col
col_times[col_times == Inf] <- NA
pheatmap(col_times[,1:200],cluster_rows = FALSE, cluster_cols = FALSE, main="colonisation times (first 200 patients)")
pheatmap((MRSA_inference$phi)[,1:200],cluster_rows = FALSE, cluster_cols = FALSE, main=expression(paste(phi, " (first 200 patients)")))
plot(rowSums(MRSA_inference$phi),type="l",main=expression(paste(Sigma,phi)))
mean(rowSums(MRSA_inference$phi)); sum(data_MRSA$col_admit)

# trace plot for the RWM tuning parameter
plot(1:(M+1),MRSA_inference$lambda,type="l",main=expression(lambda))
MRSA_inference$lambda[10001]

# investigating correlation between parameters / latent variables
plot(MRSA_inference$z,MRSA_inference$p)
plot(MRSA_inference$z,MRSA_inference$beta)
plot(MRSA_inference$p,MRSA_inference$beta)
plot(MRSA_inference$p,rowSums(MRSA_inference$phi))
plot(MRSA_inference$beta,rowSums(is.finite(MRSA_inference$t_col)))
plot(rowSums(MRSA_inference$phi),rowSums(is.finite(MRSA_inference$t_col)))
cor(MRSA_inference$z,MRSA_inference$p)
cor(MRSA_inference$z,MRSA_inference$beta)
cor(MRSA_inference$p,MRSA_inference$beta)
cor(MRSA_inference$p,rowSums(MRSA_inference$phi))
cor(MRSA_inference$beta,rowSums(is.finite(MRSA_inference$t_col)))
cor(rowSums(MRSA_inference$phi),rowSums(is.finite(MRSA_inference$t_col)))

# histograms for the parameters
hist(MRSA_inference$z,main="z", freq=F); abline(v=0.8,col="red")
hist(MRSA_inference$p,main="p", freq=F); abline(v=0.05,col="red")
hist(MRSA_inference$beta,main=expression(beta),freq=F); abline(v=0.1,col="red")

# number of infective individuals over time
plot(1:500, apply(MRSA_inference$sum_I, 2, quantile, probs=0.975),type="n")
polygon(c(1:500, 500:1), c(apply(MRSA_inference$sum_I, 2, quantile, probs=0.975), rev(apply(MRSA_inference$sum_I, 2, quantile, probs=0.025))), col = "lightgreen",border="lightgreen")
lines(apply(MRSA_inference$sum_I, 2, median),type="l",col="blue")
lines(data_MRSA$I_star, col="red")

# histograms for select individuals
# 20: infected on arrival (6-10)
# 43: never infected (15-20)
# 66: infected in the hospital (25-30, 26)
index_20 <- c("arrival",6:10,"never")
index_43 <- c("arrival",15:20,"never")
index_66 <- c("arrival",25:30,"never")
densities_20 <- rep(0,times=7)
densities_43 <- rep(0,times=8)
densities_66 <- rep(0,times=8)
for(k in 1:(M+1)){
  if(is.finite(MRSA_inference$t_col[k,20])){
    time <- MRSA_inference$t_col[k,20]
    densities_20[time-4] <- densities_20[time-4]+1
  }
  else if(MRSA_inference$phi[k,20]==1){
    densities_20[1] <- densities_20[1]+1
  }
  else{
    densities_20[7] <- densities_20[7]+1
  }
}
for(k in 1:(M+1)){
  if(is.finite(MRSA_inference$t_col[k,43])){
    time <- MRSA_inference$t_col[k,43]
    densities_43[time-13] <- densities_43[time-13]+1
  }
  else if(MRSA_inference$phi[k,43]==1){
    densities_43[1] <- densities_43[1]+1
  }
  else{
    densities_43[8] <- densities_43[8]+1
  }
}
for(k in 1:(M+1)){
  if(is.finite(MRSA_inference$t_col[k,66])){
    time <- MRSA_inference$t_col[k,66]
    densities_66[time-23] <- densities_66[time-23]+1
  }
  else if(MRSA_inference$phi[k,20]==1){
    densities_66[1] <- densities_66[1]+1
  }
  else{
    densities_66[8] <- densities_66[8]+1
  }
}
densities_20 <- densities_20/(M+1)
densities_43 <- densities_43/(M+1)
densities_66 <- densities_66/(M+1)
barplot(densities_20,space=0,names.arg=index_20,col=c("red","grey","grey","grey","grey","grey","grey"))
barplot(densities_43,space=0,names.arg=index_43,col=c("grey","grey","grey","grey","grey","grey","grey","red"))
barplot(densities_66,space=0,names.arg=index_66,col=c("grey","grey","red","grey","grey","grey","grey","grey"))

# trace plots for the parameters in ggplot
df_z <- as.data.frame(cbind(1:(M+1),MRSA_inference$z))
trace_z <- ggplot() + geom_line(aes(x=V1,y=V2),df_z) + geom_hline(yintercept=0.8, color = "red", linewidth=1.2)
trace_z <- trace_z + xlab("Iteration") + ylab("z") + ggtitle("Trace Plot of z") + theme_classic()
print(trace_z)
df_p <- as.data.frame(cbind(1:(M+1),MRSA_inference$p))
trace_p <- ggplot() + geom_line(aes(x=V1,y=V2),df_p) + geom_hline(yintercept=0.05, color = "red", linewidth=1.2)
trace_p <- trace_p + xlab("Iteration") + ylab("p") + ggtitle("Trace Plot of p") + theme_classic()
print(trace_p)
df_beta <- as.data.frame(cbind(1:(M+1),MRSA_inference$beta))
trace_beta <- ggplot() + geom_line(aes(x=V1,y=V2),df_beta) + geom_hline(yintercept=0.1, color = "red", linewidth=1.2)
trace_beta <- trace_beta + xlab("Iteration") + ylab(expression(beta)) + ggtitle(expression(paste("Trace Plot of ", beta))) + theme_classic()
print(trace_beta)

# histograms for the parameters in ggplot
df_z_hist <- as.data.frame(MRSA_inference$z)
hist_z <- ggplot() + geom_histogram(aes(x=MRSA_inference$z, y = ..density..),df_z_hist,bins=15) + geom_vline(xintercept=0.8, color = "red", linewidth=1.2)
hist_z <- hist_z + xlab("z") + ylab("Density") + ggtitle("Histogram of z") + theme_classic()
print(hist_z)
df_p_hist <- as.data.frame(MRSA_inference$p)
hist_p <- ggplot() + geom_histogram(aes(x=MRSA_inference$p, y = ..density..),df_p_hist,bins=15) + geom_vline(xintercept=0.05, color = "red", linewidth=1.2)
hist_p <- hist_p + xlab("p") + ylab("Density") + ggtitle("Histogram of p") + theme_classic()
print(hist_p)
df_beta_hist <- as.data.frame(MRSA_inference$beta)
hist_beta <- ggplot() + geom_histogram(aes(x=MRSA_inference$beta, y = ..density..),df_beta_hist,bins=15) + geom_vline(xintercept=0.1, color = "red", linewidth=1.2)
hist_beta <- hist_beta + xlab(expression(beta)) + ylab("Density") + ggtitle(expression(paste("Histogram of ", beta))) + theme_classic()
print(hist_beta)

# plot of the scale function and trace of lambda in ggplot
df_scale <- as.data.frame(cbind(1:(M+1),scale(1:(M+1))))
scale_plot <- ggplot() + geom_line(aes(x=V1,y=V2),df_scale) + geom_hline(yintercept=1, color = "blue", linewidth=0.6)
scale_plot <- scale_plot + xlab("Iteration") + ylab("Scaling Factor") + ggtitle("Scaling Factor") + theme_classic()
print(scale_plot)
df_lambda <- as.data.frame(cbind(1:(M+1),MRSA_inference$lambda))
trace_lambda <- ggplot() + geom_line(aes(x=V1,y=V2),df_lambda)
trace_lambda <- trace_lambda + xlab("Iteration") + ylab(expression(lambda)) + ggtitle(expression(paste("Trace Plot of ", lambda))) + theme_classic()
print(trace_lambda)

# correlation plots in ggplot
df_corr_z_p <- as.data.frame(cbind(MRSA_inference$z,MRSA_inference$p))
corr_z_p <- ggplot() + geom_point(aes(x=V1,y=V2),df_corr_z_p)
corr_z_p <- corr_z_p + xlab("z") + ylab("p") + ggtitle("Correlation between z and p") + theme_classic()
print(corr_z_p)
df_corr_z_beta <- as.data.frame(cbind(MRSA_inference$z,MRSA_inference$beta))
corr_z_beta <- ggplot() + geom_point(aes(x=V1,y=V2),df_corr_z_beta)
corr_z_beta <- corr_z_beta + xlab("z") + ylab(expression(beta)) + ggtitle(expression(paste("Correlation between z and ", beta))) + theme_classic()
print(corr_z_beta)
df_corr_p_beta <- as.data.frame(cbind(MRSA_inference$p,MRSA_inference$beta))
corr_p_beta <- ggplot() + geom_point(aes(x=V1,y=V2),df_corr_p_beta)
corr_p_beta <- corr_p_beta + xlab("p") + ylab(expression(beta)) + ggtitle(expression(paste("Correlation between p and ", beta))) + theme_classic()
print(corr_p_beta)
df_corr_p_phi <- as.data.frame(cbind(MRSA_inference$p,rowSums(MRSA_inference$phi)))
corr_p_phi <- ggplot() + geom_point(aes(x=V1,y=V2),df_corr_p_phi)
corr_p_phi <- corr_p_phi + xlab("p") + ylab(expression(paste(Sigma,phi))) + ggtitle(expression(paste("Correlation between p and ", Sigma,phi))) + theme_classic()
print(corr_p_phi)
df_corr_beta_t <- as.data.frame(cbind(MRSA_inference$beta,rowSums(is.finite(MRSA_inference$t_col))))
corr_beta_t <- ggplot() + geom_point(aes(x=V1,y=V2),df_corr_beta_t)
corr_beta_t <- corr_beta_t + xlab(expression(beta)) + ylab("Number of Finite Infection Times") + ggtitle(expression(paste("Correlation between ", beta, " and the Number of Finite Infection Times"))) + theme_classic()
print(corr_beta_t)
df_corr_phi_t <- as.data.frame(cbind(rowSums(MRSA_inference$phi),rowSums(is.finite(MRSA_inference$t_col))))
corr_phi_t <- ggplot() + geom_point(aes(x=V1,y=V2),df_corr_phi_t)
corr_phi_t <- corr_phi_t + xlab(expression(paste(Sigma,phi))) + ylab("Number of Finite Infection Times") + ggtitle(expression(paste("Correlation between ", Sigma, phi, " and the Number of Finite Infection Times"))) + theme_classic()
print(corr_phi_t)

# plots for latent variables in ggplot
df_phi <- as.data.frame(cbind(1:(M+1),rowSums(MRSA_inference$phi)))
trace_phi <- ggplot() + geom_line(aes(x=V1,y=V2),df_phi) + geom_hline(yintercept=sum(data_MRSA$col_admit), color = "red", linewidth=1.2)
trace_phi <- trace_phi + xlab("Iteration") + ylab(expression(paste(Sigma,phi))) + ggtitle(expression(paste("Trace Plot of ", Sigma, phi))) + theme_classic()
print(trace_phi)
df_phi_hist <- as.data.frame(rowSums(MRSA_inference$phi))
hist_phi <- ggplot() + geom_histogram(aes(x=rowSums(MRSA_inference$phi), y = ..density..),df_phi_hist,bins=15) + geom_vline(xintercept=sum(data_MRSA$col_admit), color = "red", linewidth=1.2)
hist_phi <- hist_phi + xlab(expression(paste(Sigma,phi))) + ylab("Density") + ggtitle(expression(paste("Histogram of ", Sigma, phi))) + theme_classic()
print(hist_phi)
df_t <- as.data.frame(cbind(1:(M+1),rowSums(is.finite(MRSA_inference$t_col))))
trace_t <- ggplot() + geom_line(aes(x=V1,y=V2),df_t) + geom_hline(yintercept=814, color = "red", linewidth=1.2)
trace_t <- trace_t + xlab("Iteration") + ylab("Number of Finite Infection Times") + ggtitle("Trace Plot of the Number of Finite Infection Times") + theme_classic()
print(trace_t)
df_t_hist <- as.data.frame(rowSums(is.finite(MRSA_inference$t_col)))
hist_t <- ggplot() + geom_histogram(aes(x=rowSums(is.finite(MRSA_inference$t_col)), y = ..density..),df_t_hist,bins=15) + geom_vline(xintercept=814, color = "red", linewidth=1.2)
hist_t <- hist_t + xlab("Number of Finite Infection Times") + ylab("Density") + ggtitle("Histogram of the Number of Finite Infection Times") + theme_classic()
print(hist_t)

# number of infective individuals over time in ggplot
df_infectives <- as.data.frame(cbind(1:500,apply(MRSA_inference$sum_I, 2, quantile, probs=c(0.025)),apply(MRSA_inference$sum_I, 2, quantile, probs=c(0.5)),apply(MRSA_inference$sum_I, 2, quantile, probs=c(0.975)),data_MRSA$I_star))
infectives1 <- ggplot(df_infectives,aes(x=V1)) + xlim(1, 100) + geom_ribbon(aes(ymin=V2,ymax=V4,fill="lightgreen")) + geom_line(aes(y=V3,col="blue")) + geom_line(aes(y=V5,col="red")) + scale_fill_identity(name = '', guide = 'legend',labels = c('95% Credible Interval')) + scale_colour_manual(name = '', values =c('blue'='blue','red'='red'), labels = c('Median','True Value'))
infectives1 <- infectives1 + xlab("") + ylab("") + ggtitle("Number of Infective Individuals Over Time") + theme_classic() + theme(legend.position=c(0.85,0.95),legend.direction="horizontal")
print(infectives1)
infectives2 <- ggplot(df_infectives,aes(x=V1)) + xlim(101, 200) + geom_ribbon(aes(ymin=V2,ymax=V4),fill="lightgreen") + geom_line(aes(y=V3),col="blue") + geom_line(aes(y=V5),col="red")
infectives2 <- infectives2 + xlab("") + ylab("") + ggtitle("") + theme_classic()
print(infectives2)
infectives3 <- ggplot(df_infectives,aes(x=V1)) + xlim(201, 300) + geom_ribbon(aes(ymin=V2,ymax=V4),fill="lightgreen") + geom_line(aes(y=V3),col="blue") + geom_line(aes(y=V5),col="red")
infectives3 <- infectives3 + xlab("") + ylab("Infective Individuals") + ggtitle("") + theme_classic()
print(infectives3)
infectives4 <- ggplot(df_infectives,aes(x=V1)) + xlim(301, 400) + geom_ribbon(aes(ymin=V2,ymax=V4),fill="lightgreen") + geom_line(aes(y=V3),col="blue") + geom_line(aes(y=V5),col="red")
infectives4 <- infectives4 + xlab("") + ylab("") + ggtitle("") + theme_classic()
print(infectives4)
infectives5 <- ggplot(df_infectives,aes(x=V1)) + xlim(401, 500) + geom_ribbon(aes(ymin=V2,ymax=V4),fill="lightgreen") + geom_line(aes(y=V3),col="blue") + geom_line(aes(y=V5),col="red")
infectives5 <- infectives5 + xlab("Day") + ylab("") + ggtitle("") + theme_classic() 
print(infectives5)

# histograms for select individuals in ggplot
index_20 <- factor(index_20,levels=index_20)
df_20 <- data.frame(index=index_20,dens=densities_20)
hist_20 <- ggplot(df_20, aes(x=index, y=dens)) + geom_bar(stat="identity",fill=c("red","#595959","#595959","#595959","#595959","#595959","#595959"),width = 1)
hist_20 <- hist_20 + xlab("Infection Time") + ylab("Density") + ggtitle("Histogram of Infection Times for Individual 20") + theme_classic()
print(hist_20)
index_43 <- factor(index_43,levels=index_43)
df_43 <- data.frame(index=index_43,dens=densities_43)
hist_43 <- ggplot(df_43, aes(x=index, y=dens)) + geom_bar(stat="identity",fill=c("#595959","#595959","#595959","#595959","#595959","#595959","#595959","red"),width = 1)
hist_43 <- hist_43 + xlab("Infection Time") + ylab("Density") + ggtitle("Histogram of Infection Times for Individual 43") + theme_classic()
print(hist_43)
index_66 <- factor(index_66,levels=index_66)
df_66 <- data.frame(index=index_66,dens=densities_66)
hist_66 <- ggplot(df_66, aes(x=index, y=dens)) + geom_bar(stat="identity",fill=c("#595959","#595959","red","#595959","#595959","#595959","#595959","#595959"),width = 1)
hist_66 <- hist_66 + xlab("Infection Time") + ylab("Density") + ggtitle("Histogram of Infection Times for Individual 66") + theme_classic()
print(hist_66)