################### Base model ################
############### Three-strategy decision model with grower vs alternative comparison ####################
library(deSolve)
############### Here growers will first evaluate which alternative strategy is higher in payoff ###########
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

par(mar = c(5.1, 6.1, 3.2, 3.2))

### Initial conditions - S_U, E_U, I_U, S_T, E_T, I_T, S_R, E_R, I_R
yini <- c((1 - 0.01 - 0.1 - 0.1), 0, 0.01, 0.1, 0, 0, 0.1, 0, 0 )

### Length of epidemic 
times <- seq(0, 100*120, 100)

### Default parameters
parameters <- readParams(N = 1, beta = 0.055, delta_beta_T = 1, delta_beta_R = 0.5, gamma = 1/120, Y = 1, delta_Y_T = 1, delta_Y_R = 1, phi_T = 0.1,
                         phi_R = 0.1, phi_Q = 0.7, L = 0.6, delta_L_T = 0.1, delta_L_R = 1,mu_U = 1/(((1/1) - 1/2)*120), mu_T =1/(((1/0.1) - 1/2)*120), mu_R =1/(((1/1) - 1/2)*120),
                         eta = 10, delta_sigma_T= 1, delta_sigma_R = 0.5, epsilon = 1/41, delta_epsilon_T = 1, delta_epsilon_R = 1/2)


### Simulation and plotting
out <- ode(func =Tolerance_resistance_conventional_model, y = yini, time = times, parms = parameters)
plot(out[,1], out[,2], lwd = 4, ty = "l", ylim = c(0,1), col = "#274c77", main = "", ylab = expression("Proportion of growers"), xlab = expression("Time (seasons)"), xaxt = "n", cex.lab = 1.9, cex.axis = 1.6)

axis(side = 1, at = seq(0, max(times ), max(times /5)),labels = seq(0, max(times /120), max(times /120)/5), cex.axis = 1.6)
lines(out[,1], out[,3], lwd = 4, col = "#6096ba")
lines(out[,1], out[,4], lwd = 4, col = "#a3cef1")
lines(out[,1], out[,5], lwd = 4, col = "#436436")
lines(out[,1], out[,6], lwd = 4, lty = 1, col = "#80B06D")
lines(out[,1], out[,7], lwd = 4, lty = 1, col = "#AFE996")
lines(out[,1], out[,8], lwd = 4, lty = 1, col = "#e2711d")
lines(out[,1], out[,9], lwd = 4, lty = 1, col = "#ff9505")
lines(out[,1], out[,10], lwd = 4, lty = 1, col = "#ffb627")
legend("topright", legend = c(expression(paste("S"[U]*"")),expression(paste("S"[T]*"")), expression(paste("S"[R]*"")),expression(paste("E"[U]*"")),expression(paste("E"[T]*"")), expression(paste("E"[R]*"")),expression(paste("I"[U]*"")),expression(paste("I"[T]*"")), expression(paste("I"[R]*""))),lty = 1, lwd = 3, col = c("#274c77", "#436436", "#e2711d", "#6096ba","#80B06D", "#ff9505", "#A3CEF1","#AFE996","#ffb627"), cex = 1.4, ncol = 3)
