############################# Pareto optimisation ####################################
############################# Libraries ##############################################
library(deSolve)
library(fields)
library(colorspace)
############################# Functions ##############################################
Tolerance_resistance_conventional_model <- function(time, y, params){
  #Initial conditions
  SU <- y[1]
  EU <- y[2]
  IU <- y[3]
  ST <- y[4]
  ET <- y[5]
  IT <- y[6]
  SR <- y[7]
  ER <- y[8]
  IR <- y[9]
  
  #Parameters
  beta_R <- params$beta*params$delta_beta_R
  beta_T <-params$beta*params$delta_beta_T
  beta_U <- params$beta
  epsilon_R <- params$epsilon*params$delta_epsilon_R
  epsilon_T <- params$epsilon*params$delta_epsilon_T
  epsilon_U <- params$epsilon
  eta_R <- params$eta
  eta_T <- params$eta
  eta_U <- params$eta
  gamma <- params$gamma
  mu_R <- params$mu_R
  mu_T <- params$mu_T
  mu_U <- params$mu_U 
  Y_R <- params$Y*params$delta_Y_R
  Y_T <- params$Y*params$delta_Y_T
  Y_U <- params$Y
  phi_R <- params$phi_R
  phi_T <- params$phi_T
  phi_rog <- params$phi_Q
  L_R <- params$L*params$delta_L_R
  L_T <-params$L*params$delta_L_T
  L_U <- params$L
  sigma <- params$delta_sigma_R
  #Outcome profits
  P_SR <- Y_R - phi_R
  P_ER <- Y_R - phi_R
  P_IRH <- Y_R - phi_R - L_R
  P_IRR <- Y_R - phi_R - phi_rog*L_R
  
  P_ST <- Y_T - phi_T
  P_ET <- Y_T - phi_T
  P_ITH <- Y_T - phi_T - L_T
  P_ITR <- Y_T - phi_T - phi_rog*L_T
  
  P_SU <- Y_U
  P_EU <- Y_U
  P_IUH <- Y_U - L_U
  P_IUR <- Y_U - phi_rog*L_U
  
  #Expected profits by strategy
  q_R <- (beta_R*(sigma*IR + IT + IU))/(beta_R*(sigma*IR + IT + IU) + gamma)
  q_ER <- q_R*(gamma/(gamma + epsilon_R))
  q_IR <- q_R*(epsilon_R/(gamma + epsilon_R))
  q_IRH <- q_IR*(gamma/(gamma + mu_R))
  q_IRR <- q_IR*(mu_R/(gamma + mu_R))
  
  q_T <- (beta_T*(sigma*IR + IT + IU))/(beta_T*(sigma*IR + IT + IU) + gamma)
  q_ET <- q_T*(gamma/(gamma + epsilon_T))
  q_IT <- q_T*(epsilon_T/(gamma + epsilon_T))
  q_ITH <- q_IT*(gamma/(gamma + mu_T))
  q_ITR <- q_IT*(mu_T/(gamma + mu_T))
  
  
  q_U <- (beta_U*(sigma*IR + IT + IU))/(beta_U*(sigma*IR + IT + IU) + gamma)
  q_EU <- q_U*(gamma/(gamma + epsilon_U))
  q_IU <- q_U*(epsilon_U/(gamma + epsilon_U))
  q_IUH <- q_IU*(gamma/(gamma + mu_U))
  q_IUR <- q_IU*(mu_U/(gamma + mu_U))
  
  P_R <- (1 - q_R)*P_SR + q_ER*P_ER + q_IRH*P_IRH + q_IRR*P_IRR
  P_T <- (1 - q_T)*P_ST + q_ET*P_ET + q_ITH*P_ITH + q_ITR*P_ITR
  P_U <- (1 - q_U)*P_SU + q_EU*P_EU + q_IUH*P_IUH + q_IUR*P_IUR
  
  #Totals
  R <- SR + ER + IR
  Tol <- ST + ET + IT
  U <- SU + EU + IU
  #Switching terms
  if(P_R > P_T | (P_T == P_R & R > Tol)){ ## U INTO R
    z_SUR <- max(0, 1-exp(-eta_U*(P_R - P_SU)))
    z_EUR <- max(0, 1-exp(-eta_U*(P_R - P_EU)))
    z_IUHR <- max(0, 1-exp(-eta_U*(P_R - P_IUH)))
    z_IURR <- max(0, 1-exp(-eta_U*(P_R - P_IUR)))
    z_SUT <- 0
    z_EUT <- 0
    z_IUHT <- 0
    z_IURT <- 0
  } else if(P_R < P_T | (P_T == P_R & R < Tol)){ ## U INTO T
    z_SUR <-0
    z_EUR <- 0
    z_IUHR <- 0
    z_IURR <- 0
    z_SUT <- max(0, 1-exp(-eta_U*(P_T - P_SU)))
    z_EUT <- max(0, 1-exp(-eta_U*(P_T - P_EU)))
    z_IUHT <- max(0, 1-exp(-eta_U*(P_T - P_IUH)))
    z_IURT <- max(0, 1-exp(-eta_U*(P_T - P_IUR)))
  }
  if(P_R > P_U | (P_U == P_R & R > U)){ ## T INTO R
    z_STU <- 0
    z_ETU <- 0
    z_ITHU <- 0
    z_ITRU <- 0
    z_STR <- max(0, 1-exp(-eta_T*(P_R - P_ST)))
    z_ETR <- max(0, 1-exp(-eta_T*(P_R - P_ET)))
    z_ITHR <- max(0, 1-exp(-eta_T*(P_R - P_ITH)))
    z_ITRR <- max(0, 1-exp(-eta_T*(P_R - P_ITR)))
  } else if(P_R < P_U | (P_U == P_R & R < U)){ ## T INTO U
    z_STU <- max(0, 1-exp(-eta_T*(P_U - P_ST)))
    z_ETU <- max(0, 1-exp(-eta_T*(P_U - P_ET)))
    z_ITHU <- max(0, 1-exp(-eta_T*(P_U - P_ITH)))
    z_ITRU <- max(0, 1-exp(-eta_T*(P_U - P_ITR)))
    z_STR <- 0
    z_ETR <- 0
    z_ITHR <- 0
    z_ITRR <- 0
  }
  if(P_T > P_U | (P_U == P_T & Tol > U)){ ## R INTO T
    z_SRU <- 0
    z_ERU <- 0
    z_IRHU <- 0
    z_IRRU <- 0
    z_SRT <- max(0, 1-exp(-eta_R*(P_T - P_SR)))
    z_ERT <- max(0, 1-exp(-eta_R*(P_T - P_ER)))
    z_IRHT <- max(0, 1-exp(-eta_R*(P_T - P_IRH)))
    z_IRRT <- max(0, 1-exp(-eta_R*(P_T - P_IRR)))
  } else if(P_T < P_U | (P_U == P_T & Tol < U)){ ## R INTO U
    z_SRU <- max(0, 1-exp(-eta_R*(P_U - P_ST)))
    z_ERU <- max(0, 1-exp(-eta_R*(P_U - P_ET)))
    z_IRHU <- max(0, 1-exp(-eta_R*(P_U - P_ITH)))
    z_IRRU <- max(0, 1-exp(-eta_R*(P_U - P_ITR)))
    z_SRT <- 0
    z_ERT <- 0
    z_IRHT <- 0
    z_IRRT <- 0
  }
  
  z_SU <- max(z_SUT, z_SUR)
  z_EU <- max(z_EUT, z_EUR)
  z_IUH <- max(z_IUHT, z_IUHR)
  z_IUR <- max(z_IURT, z_IURR)
  z_ST <- max(z_STU, z_STR)
  z_ET <- max(z_ETU, z_ETR)
  z_ITH <- max(z_ITHU, z_ITHR)
  z_ITR <- max(z_ITRU, z_ITRR)
  z_SR <- max(z_SRU, z_SRT)
  z_ER <- max(z_ERU, z_ERT)
  z_IRH <- max(z_IRHU, z_IRHT)
  z_IRR <- max(z_IRRU, z_IRRT)
  
  ## Model
  
  dSU <- gamma*((1-z_SU)*SU + (1- z_EU)*EU + (1 - z_IUH)*IU + z_STU*ST + z_ETU*ET + z_ITHU*IT + z_SRU*SR + z_ERU*ER + z_IRHU*IR) - beta_U*SU*(IU + IT + sigma*IR) + mu_U*(1 - z_IUR)*IU + mu_T*z_ITRU*IT + mu_R*z_IRRU*IR - gamma*SU
  dEU <-  beta_U*SU*(IU + IT + sigma*IR) - epsilon_U*EU - gamma*EU
  dIU <- epsilon_U*EU - mu_U*IU - gamma*IU
  
  dST <- gamma*((1-z_ST)*ST + (1- z_ET)*ET + (1 - z_ITH)*IT + z_SUT*SU + z_EUT*EU + z_IUHT*IU + z_SRT*SR + z_ERT*ER + z_IRHT*IR) - beta_T*ST*(IU + IT + sigma*IR) + mu_T*(1 - z_ITR)*IT + mu_U*z_IURT*IU + mu_R*z_IRRT*IR - gamma*ST
  dET <-  beta_T*ST*(IU + IT + sigma*IR) - epsilon_T*ET - gamma*ET
  dIT <- epsilon_T*ET - mu_T*IT - gamma*IT
  
  dSR <- gamma*((1-z_SR)*SR + (1- z_ER)*ER + (1 - z_IRH)*IR + z_SUR*SU + z_EUR*EU + z_IUHR*IU + z_STR*ST + z_ETR*ET + z_ITHR*IT) - beta_R*SR*(IU + IT + sigma*IR) + mu_R*(1 - z_IRR)*IR + mu_U*z_IURR*IU + mu_T*z_ITRR*IT - gamma*SR
  dER <-  beta_R*SR*(IU + IT + sigma*IR) - epsilon_R*ER - gamma*ER
  dIR <- epsilon_R*ER - mu_R*IR - gamma*IR
  
  return(list(c(dSU, dEU, dIU, dST, dET, dIT, dSR, dER, dIR)))
}
readParams <- function(N, beta,delta_beta_T,delta_beta_R,gamma,Y,delta_Y_T,delta_Y_R, phi_T,phi_R,phi_Q,L, delta_L_T,delta_L_R,mu_U,mu_T,mu_R,eta,delta_sigma_T,delta_sigma_R,epsilon,delta_epsilon_T, delta_epsilon_R)
{
  retval <- list(N = N, 
                 beta = beta/N,
                 delta_beta_T = delta_beta_T,
                 delta_beta_R = delta_beta_R,
                 gamma = gamma, 
                 Y = Y, 
                 delta_Y_T = delta_Y_T, 
                 delta_Y_R = delta_Y_R, 
                 phi_T= phi_T,
                 phi_R= phi_R,
                 phi_Q = phi_Q, 
                 L = L, 
                 delta_L_T = delta_L_T,
                 delta_L_R = delta_L_R,
                 mu_U = mu_U, 
                 mu_T = mu_T,
                 mu_R = mu_R,
                 eta = eta,
                 delta_sigma_T = delta_sigma_T,
                 delta_sigma_R = delta_sigma_R,
                 epsilon = epsilon,
                 delta_epsilon_T = delta_epsilon_T,
                 delta_epsilon_R = delta_epsilon_R)
  
  return(retval)
}

##### First, need two-way parameter scan over phi_T and phi_R #######################
##### Parameters are "ineffective tolerance, ineffective resistance", so delta_L_T = 0.5 and delta_beta_R = 0.5 
parameters <- readParams(N = 1, beta = 0.055, delta_beta_T = 1, delta_beta_R = 0.5, gamma = 1/120, Y = 1, delta_Y_T = 1, delta_Y_R = 1, phi_T = 0.1,
                         phi_R = 0.1, phi_Q = 0.7, L = 0.6, delta_L_T = 0.5, delta_L_R = 1,mu_U = 1/(((1/1) - 1/2)*120), mu_T =1/(((1/0.1) - 1/2)*120), mu_R =1/(((1/1) - 1/2)*120),
                         eta = 10, delta_sigma_T= 1, delta_sigma_R = 0.5, epsilon = 1/41, delta_epsilon_T = 1, delta_epsilon_R = 1/2)
### Values to scan over
phi_T_seq <- seq(0, 0.4, 0.005)
phi_R_seq <- seq(0, 0.4, 0.005)

### Initial conditions - S_U, E_U, I_U, S_T, E_T, I_T, S_R, E_R, I_R
yini <- c((1 - 0.01 - 0.1 - 0.1), 0, 0.01, 0.1, 0, 0, 0.1, 0, 0 )

### Length of epidemic 
times <- seq(0, 100*120, 100)

profit_matrix_bad_t_bad_r <- matrix(NA, nrow = length(phi_T_seq), ncol = length(phi_R_seq))
infectious_matrix_bad_t_bad_r <- matrix(NA, nrow = length(phi_T_seq), ncol = length(phi_R_seq))
cost_matrix_bad_t_bad_r <- matrix(NA, nrow = length(phi_T_seq), ncol = length(phi_R_seq))
equilibrium_matrix_bad_t_bad_r <- matrix(NA, nrow = length(phi_T_seq), ncol = length(phi_R_seq))
tolerance_matrix_bad_t_bad_r <- matrix(NA, nrow = length(phi_T_seq), ncol = length(phi_R_seq))
resistance_matrix_bad_t_bad_r <- matrix(NA, nrow = length(phi_T_seq), ncol = length(phi_R_seq))
average_infectious_matrix_bad_t_bad_r <- matrix(NA, nrow = length(phi_T_seq), ncol = length(phi_R_seq))
average_profit_matrix_bad_t_bad_r <- matrix(NA, nrow = length(phi_T_seq), ncol = length(phi_R_seq))
average_cost_matrix_bad_t_bad_r <- matrix(NA, nrow = length(phi_T_seq), ncol = length(phi_R_seq))
average_tolerance_matrix_bad_t_bad_r <- matrix(NA, nrow = length(phi_T_seq), ncol = length(phi_R_seq))
average_resistance_matrix_bad_t_bad_r <- matrix(NA, nrow = length(phi_T_seq), ncol = length(phi_R_seq))
average_equilibrium_matrix_bad_t_bad_r <- matrix(NA, nrow = length(phi_T_seq), ncol = length(phi_R_seq))

z <- 1
for(phiT in phi_T_seq){
  j <- 1
  parameters$phi_T <- phiT
  
  for(phiR in phi_R_seq){
    
    parameters$phi_R <- phiR
    out <- ode(func =Tolerance_resistance_conventional_model,  y = yini, time = times, parms = parameters)
    SU <- mean(tail(out,50)[,2])
    ST <- mean(tail(out,50)[,5])
    SR <- mean(tail(out,50)[,8])
    EU <- mean(tail(out,50)[,3])
    ET <- mean(tail(out,50)[,6])
    ER <- mean(tail(out,50)[,9])
    IU <- mean(tail(out,50)[,4])
    IT <- mean(tail(out,50)[,7])
    IR <- mean(tail(out,50)[,10])
    
    infectious_matrix_bad_t_bad_r[z,j] <- IU + IT + IR
    #profit <- P_U*(SU + EU + IU) + P_T*(ST + ET + IT) + P_R*(SR + ER + IR)
    P_SR <- parameters$Y_R - parameters$phi_R
    P_ER <- parameters$Y_R - parameters$phi_R
    P_IRH <- parameters$Y_R - parameters$phi_R - parameters$L_R
    P_IRR <- parameters$Y_R - parameters$phi_R - parameters$phi_rog*parameters$L_R
    
    P_ST <- parameters$Y_T - parameters$phi_T
    P_ET <- parameters$Y_T - parameters$phi_T
    P_ITH <- parameters$Y_T - parameters$phi_T - parameters$L_T
    P_ITR <- parameters$Y_T - parameters$phi_T - parameters$phi_rog*parameters$L_T
    
    tolerance_matrix_bad_t_bad_r[z,j] <- ST + ET + IT
    resistance_matrix_bad_t_bad_r[z,j] <- SR + ER + IR
    profit_matrix_bad_t_bad_r[z,j] <- SU*P_SU + EU*P_EU + (parameters$mu_U/(parameters$gamma + parameters$mu_U))*IU*P_IUR + (parameters$gamma/(parameters$gamma + parameters$mu_U))*IU*P_IUH + ST*P_ST + ET*P_ET + (parameters$mu_T/(parameters$gamma + parameters$mu_T))*IT*P_ITR + (parameters$gamma/(parameters$gamma + parameters$mu_T))*IT*P_ITH + SR*P_SR + ER*P_ER + (parameters$mu_R/(parameters$gamma + parameters$mu_R))*IR*P_IRR + (parameters$gamma/(parameters$gamma + parameters$mu_R))*IR*P_IRH
    cost_matrix_bad_t_bad_r[z,j] <- (0.4 - phiR)*(SR + ER + IR) + (0.4 - phiT)*(ST + ET + IT)
    U = SU + EU + IU
    Tol = ST + ET + IT
    R = SR + ER + IR
    C = Tol + R
    I = IU + IT + IR
    if(round(C, 4) == 0 & round(I, 4) == 0){
      equilibrium_matrix_bad_t_bad_r[z,j] <- 1
    } else if(round(C, 4) == 0 & round(I, 5) > 0){
      equilibrium_matrix_bad_t_bad_r[z,j] <- 2
    } else if(round(Tol, 4) == 1 & round(I, 5) > 0){
      equilibrium_matrix_bad_t_bad_r[z,j] <- 3
    }else if(round(R, 5) == 1 & round(I, 5) > 0){
      equilibrium_matrix_bad_t_bad_r[z,j] <- 4
    }else if(round(Tol, 5)> 0 & round(R, 5)> 0 & round(C, 5) < 1 & round(I, 5) > 0){
      equilibrium_matrix_bad_t_bad_r[z,j] <-5
    }else if(round(Tol, 5) > 0 & round(U, 5) > 0 & round(I, 5) > 0){
      equilibrium_matrix_bad_t_bad_r[z,j] <- 6
    }else if(round(R, 5) > 0 & round(U, 5) > 0 & round(I, 5) > 0){
      equilibrium_matrix_bad_t_bad_r[z,j] <- 7
    }else if(round(R, 5) > 0 & round(Tol, 5) > 0 & round(I, 5) > 0){
      equilibrium_matrix_bad_t_bad_r[z,j] <- 8
    }
    
    SU <- out[,2]
    ST <- out[,5]
    SR <- out[,8]
    EU <- out[,3]
    ET <- out[,6]
    ER <- out[,9]
    IU <- out[,4]
    IT <- out[,7]
    IR <- out[,10]
    U = mean(SU + EU + IU)
    Tol = mean(ST + ET + IT)
    R = mean(SR + ER + IR)
    C = mean(Tol + R)
    I = mean(IU + IT + IR)
    profit_long <- SU*P_SU + EU*P_EU + (parameters$mu_U/(parameters$gamma + parameters$mu_U))*IU*P_IUR + (parameters$gamma/(parameters$gamma + parameters$mu_U))*IU*P_IUH + ST*P_ST + ET*P_ET + (parameters$mu_T/(parameters$gamma + parameters$mu_T))*IT*P_ITR + (parameters$gamma/(parameters$gamma + parameters$mu_T))*IT*P_ITH + SR*P_SR + ER*P_ER + (parameters$mu_R/(parameters$gamma + parameters$mu_R))*IR*P_IRR + (parameters$gamma/(parameters$gamma + parameters$mu_R))*IR*P_IRH
    average_profit_matrix_bad_t_bad_r[z,j] <- mean(profit_long)
    average_cost_matrix_bad_t_bad_r[z,j] <- mean((0.4 - phiR)*(SR + ER + IR) + (0.4 - phiT)*(ST + ET + IT))
    average_infectious_matrix_bad_t_bad_r[z,j] <- mean(IU + IT + IR)
    average_tolerance_matrix_bad_t_bad_r[z,j] <- mean(ST + ET + IT)
    average_resistance_matrix_bad_t_bad_r[z,j] <- mean(SR + ER + IR)
    if(round(C, 4) == 0 & round(I, 4) == 0){
      average_equilibrium_matrix_bad_t_bad_r[z,j] <- 1
    } else if(round(C, 4) == 0 & round(I, 5) > 0){
      average_equilibrium_matrix_bad_t_bad_r[z,j] <- 2
    } else if(round(Tol, 4) == 1 & round(I, 5) > 0){
      average_equilibrium_matrix_bad_t_bad_r[z,j] <- 3
    }else if(round(R, 5) == 1 & round(I, 5) > 0){
      average_equilibrium_matrix_bad_t_bad_r[z,j] <- 4
    }else if(round(Tol, 5)> 0 & round(R, 5)> 0 & round(C, 5) < 1 & round(I, 5) > 0){
      average_equilibrium_matrix_bad_t_bad_r[z,j] <-5
    }else if(round(Tol, 5) > 0 & round(U, 5) > 0 & round(I, 5) > 0){
      average_equilibrium_matrix_bad_t_bad_r[z,j] <- 6
    }else if(round(R, 5) > 0 & round(U, 5) > 0 & round(I, 5) > 0){
      average_equilibrium_matrix_bad_t_bad_r[z,j] <- 7
    }else if(round(R, 5) > 0 & round(Tol, 5) > 0 & round(I, 5) > 0){
      average_equilibrium_matrix_bad_t_bad_r[z,j] <- 8
    }
    j <- j + 1
    
  }
  z <- z + 1
}


###############  Pareto front using get_frontier() #############
par(mfrow = c(1,1), mar = c(5.2, 6.2, 4.1, 3.1))

############## Plotting the cost-profit combinations 

cost_seq <- sort(as.vector(unique(cost_matrix_bad_t_bad_r)))
cost_seq <- unique(cost_seq)
best_profit_dataframe_bad_t_bad_r <- data.frame(cost_seq, NA)
colnames(best_profit_dataframe_bad_t_bad_r) <- c("Cost", "Profit")
plot(NA, xlim = c(max(cost_seq), min(cost_seq)), ylim = c(min(profit_seq), max(profit_seq)), xlab = "", ylab = "")

i <- 1
for(x in cost_seq){
  which_x <- unique(profit_matrix_bad_t_bad_r[which(cost_matrix_bad_t_bad_r == x)])
  #print(which_x)
  points(rep(x, length(which_x)), which_x, pch = 19, col = "grey78", cex = 1)
  which_x <- max(profit_matrix_bad_t_bad_r[which(cost_matrix_bad_t_bad_r == x)])
  best_profit_dataframe_bad_t_bad_r[i,2] <- which_x
  i <- i + 1
}

#subsidies_bad_t_bad_r <- subsidies_bad_t_bad_r[order(-1*subsidies_bad_t_bad_r$Cost), ]
KM_frontier <- get_frontier(best_profit_dataframe_bad_t_bad_r, Cost, Profit, quadrant = "top.left",decreasing = FALSE)
pal_frontier <- sequential_hcl( length(KM_frontier$Cost), "Sunset")

points(KM_frontier$Cost, KM_frontier$Profit, col = pal_frontier, cex = .8, pch = 19)
points(KM_frontier$Cost, KM_frontier$Profit, col = pal_frontier, cex = 1.1)
lines(KM_frontier$Cost, KM_frontier$Profit, col = "gray48", lwd = 2, lty = 2)
title(ylab = expression(paste("Profit of growers (P)")), cex.lab = 2, line = 4)
title(xlab = expression(paste("Cost to planner (",tau, ")")), cex.lab = 2, line = 2.75)
legend("topright", legend = c("Pareto front"), col = "gray48", lwd = 2, lty = 2, cex = 1.5)

########### Gini coefficients #############
subsidies_bad_t_bad_r$profits_rescaled <- (max(subsidies_bad_t_bad_r$Profit) - subsidies_bad_t_bad_r$Profit)/(max(subsidies_bad_t_bad_r$Profit) - min(subsidies_bad_t_bad_r$Profit))
subsidies_bad_t_bad_r$costs_rescaled <- (min(subsidies_bad_t_bad_r$Cost) - subsidies_bad_t_bad_r$Cost)/(min(subsidies_bad_t_bad_r$Cost) - max(subsidies_bad_t_bad_r$Cost))
subsidies_bad_t_bad_r$Gini <- (abs(subsidies_bad_t_bad_r$costs_rescaled  - subsidies_bad_t_bad_r$profits_rescaled) + abs(subsidies_bad_t_bad_r$profits_rescaled - subsidies_bad_t_bad_r$costs_rescaled ))/(2*2*(subsidies_bad_t_bad_r$profits_rescaled + subsidies_bad_t_bad_r$costs_rescaled ))

which(subsidies_bad_t_bad_r$Gini == min(subsidies_bad_t_bad_r$Gini))

text(subsidies_bad_t_bad_r$Cost[c(1,24,46)] + 0.025, subsidies_bad_t_bad_r$Profit[c(1,24, 46)], round(subsidies_bad_t_bad_r$Gini[c(1,24,46)], 3), cex = 1.4)

(abs(subsidies_bad_t_bad_r$Cost[1] - subsidies_bad_t_bad_r$Profit[1]) + abs(subsidies_bad_t_bad_r$Profit[1] - subsidies_bad_t_bad_r$Cost[1]))/(2*2*(subsidies_bad_t_bad_r$Profit[1] + subsidies_bad_t_bad_r$Cost[1]))


(abs(0 - 1)+ abs(1 - 0) + 0 + 0)/(2*2*(0+ 1))
(abs(1 - 1)+ abs(1 - 1) + 0 + 0)/(2*2*(1+ 1))

(abs(0.4 - 1)+ abs(1 - 0.4) + 0 + 0)/(2*2*(0.4+ 1))