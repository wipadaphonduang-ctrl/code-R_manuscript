library(sn)
library(VGAM)
library(pwr)
library(openxlsx)
library(MASS)
library(DescTools)
library(BSDA)
library(PairedData)
library(car)
library(e1071)

#Welch base on rank test
WBR.test <- function(x1, x2, alpha = 0.05) {
  n1 <- length(x1)
  n2 <- length(x2)
  
  pooled_ranks <- rank(c(x1, x2), ties.method = "average")
  
  R1 <- pooled_ranks[seq_len(n1)]
  R2 <- pooled_ranks[n1 + seq_len(n2)]
  
  fit <- stats::t.test(R1, R2, alternative = "two.sided", var.equal = FALSE)
  
  list(
    WBR.reject = as.integer(fit$p.value < alpha),
    R1 = R1,
    R2 = R2,
    WBR = unname(fit$statistic),
    df = unname(fit$parameter),
    p.value = fit$p.value
  )
}

# Yeo–Johnson transformation with pooled parameter 
yeo_johnson_transform_pooled <- function(x1, x2) {
  
  apply_yj <- function(x, lambda) {
    transformed_x <- numeric(length(x))
    for (i in seq_along(x)) {
      if (lambda != 0 && x[i] >= 0) {
        transformed_x[i] <- (((x[i] + 1)^lambda) - 1) / lambda
      } else if (lambda == 0 && x[i] >= 0) {
        transformed_x[i] <- log(x[i] + 1)
      } else if (lambda != 2 && x[i] < 0) {
        transformed_x[i] <- -(((-x[i] + 1)^(2 - lambda)) - 1) / (2 - lambda)
      } else if (lambda == 2 && x[i] < 0) {
        transformed_x[i] <- -log(-x[i] + 1)
      }
    }
    return(transformed_x)
  }
  
  pooled_x      <- c(x1, x2)
  yeo_result    <- car::powerTransform(pooled_x, family = "yjPower")
  lambda_common <- as.numeric(yeo_result$lambda)
  
  return(list(
    transformed_x1 = apply_yj(x1, lambda_common),
    transformed_x2 = apply_yj(x2, lambda_common),
    lambda_common  = lambda_common
  ))
}

# Imposed-shift script
# Condition: with outliers in both groups
power_results <- data.frame()
set.seed(131169L)

mu    <- 1
sigma <- 1      
alpha <- 0.05
M     <- 5000

for (n1 in c(5, 7, 8, 10, 12, 15, 20, 25, 50, 100, 200, 500)) {
  for (m in c(1, 9)) {
    for (lambda in c(0, 9)) {
      #n2   <- n1
      n2 <- round(1.2*n1)
      
      dalta <- lambda/(sqrt(1+(lambda^2)))
      mean <- mu+(sigma*dalta*sqrt(2/pi))
      
      power_t_value          <- 0
      power_Welch_value      <- 0
      power_Z_value          <- 0
      power_Yuen_value       <- 0
      
      power_MWU_value        <- 0
      power_WBR_value        <- 0
      power_VW_value         <- 0
      
      power_Yeo_t_value      <- 0
      power_Yeo_Welch_value  <- 0
      power_Yeo_Z_value      <- 0
      power_Yeo_Yuen_value   <- 0
      
      for (i in 1:M) {
        # Generate data
        x1 <- (rsn(n = n1, xi = mu, omega = sigma, alpha = lambda, tau = 0, dp = NULL)-mean)
        x2 <- m *(rsn(n = n2, xi = mu, omega = sigma, alpha = lambda, tau = 0, dp = NULL)-mean)-0.5
        
        #generate outlier                                                                                                      
        k1 <- ceiling(n1*0.15)
        k2 <- ceiling(n2*0.15)
        
        x_1 <- sort(x1)
        x_2 <- sort(x2)
        
        x3 <- head(x_1, length(x_1)-k1)
        x4 <- head(x_2, length(x_2)-k2)
        
        Q3_1             <- quantile(x1,0.75)
        IQR_1            <- IQR(x1)
        innerfence3_1    <- Q3_1 + (3*IQR_1)
        outerfence4.5_1  <- Q3_1 + (4.5*IQR_1)            
        outlier_1        <- sample((seq(innerfence3_1,outerfence4.5_1,0.001)),k1)           
        y1               <- c(x3,outlier_1)
        
        Q3_2             <- quantile(x2,0.75)
        IQR_2            <- IQR(x2)
        innerfence3_2    <- Q3_2 + (3*IQR_2)
        outerfence4.5_2  <- Q3_2 + (4.5*IQR_2)            
        outlier_2        <- sample((seq(innerfence3_2,outerfence4.5_2,0.001)),k2)           
        y2               <- c(x4,outlier_2)
        
        #Parametric 
        #Student’s t-test 
        t_test        <- t.test(y1, y2, var.equal = TRUE)
        t_reject      <- ifelse(t_test$p.value < alpha, 1, 0)
        power_t_value <- power_t_value + t_reject
        
        #Welch’s t-test 
        Welch_test        <- t.test(y1, y2, var.equal = FALSE)
        Welch_reject      <- ifelse(Welch_test$p.value < alpha, 1, 0)
        power_Welch_value <- power_Welch_value + Welch_reject
        
        #Z-test 
        S1 <- sd(y1)
        S2 <- sd(y2)
        
        Z_test        <- z.test(y1, y2, sigma.x = S1, sigma.y = S2, alternative = "two.sided")
        Z_reject      <- ifelse(Z_test$p.value < alpha, 1, 0)
        power_Z_value <- power_Z_value + Z_reject
        
        #Yuen-Welch (YW) test 
        Yuen_test        <- yuen.t.test(y1, y2, tr = 0.15, alternative = "two.sided")
        Yuen_reject      <- ifelse(Yuen_test$p.value < alpha, 1, 0)
        power_Yuen_value <- power_Yuen_value + Yuen_reject
        
        #Non-parametric 
        #Mann-Whitney U (MWU) test
        MWU_test        <- wilcox.test(y1, y2, alternative = "two.sided")
        MWU_reject      <- ifelse(MWU_test$p.value < alpha, 1, 0)
        power_MWU_value <- power_MWU_value + MWU_reject
        
        #Welch test based on ranks (WBR) 
        WBR_test        <- WBR.test(y1, y2, alpha = alpha)
        WBR_reject      <- WBR_test$WBR.reject
        power_WBR_value <- power_WBR_value + WBR_reject
        
        #Van der Waerden (VW) test 
        VW_test            <- VanWaerdenTest(list(y1, y2))
        VW_reject          <- ifelse(VW_test$p.value < alpha, 1, 0)
        power_VW_value     <- power_VW_value + VW_reject
        
        #Transformed using one pooled parameter
        pooled_transform <- yeo_johnson_transform_pooled(y1, y2)
        
        transformed_x1 <- list(transformed_x = pooled_transform$transformed_x1, lambda_Yeo = pooled_transform$lambda_common)
        transformed_x2 <- list(transformed_x = pooled_transform$transformed_x2, lambda_Yeo = pooled_transform$lambda_common)
        
        lambda_Yeo1 <- transformed_x1$lambda_Yeo
        lambda_Yeo2 <- transformed_x2$lambda_Yeo

        stopifnot(identical(lambda_Yeo1, lambda_Yeo2))
        
        #Yeo t-test
        Yeo_t_test <- t.test(transformed_x1$transformed_x, transformed_x2$transformed_x, var.equal = TRUE)
        Yeo_t_reject <- ifelse(Yeo_t_test$p.value < alpha, 1, 0)
        power_Yeo_t_value <- power_Yeo_t_value + Yeo_t_reject
        
        #Yeo Welch t-test
        Yeo_Welch_test <- t.test(transformed_x1$transformed_x, transformed_x2$transformed_x, var.equal = FALSE)
        Yeo_Welch_reject <- ifelse(Yeo_Welch_test$p.value < alpha, 1, 0)
        power_Yeo_Welch_value <- power_Yeo_Welch_value + Yeo_Welch_reject
        
        #Yeo Z-test
        Yeo_S1 <- sd(transformed_x1$transformed_x)
        Yeo_S2 <- sd(transformed_x2$transformed_x)
        
        Yeo_Z_test <- z.test(transformed_x1$transformed_x, transformed_x2$transformed_x, sigma.x = Yeo_S1, sigma.y = Yeo_S2, alternative = "two.sided")
        Yeo_Z_reject <- ifelse(Yeo_Z_test$p.value < alpha, 1, 0)
        power_Yeo_Z_value <- power_Yeo_Z_value + Yeo_Z_reject
        
        #Yeo Yuen Welch
        Yeo_Yuen_test <- yuen.t.test(transformed_x1$transformed_x, transformed_x2$transformed_x, tr = 0.15, alternative = "two.sided")
        Yeo_Yuen_reject <- ifelse(Yeo_Yuen_test$p.value < alpha, 1, 0)
        power_Yeo_Yuen_value <- power_Yeo_Yuen_value + Yeo_Yuen_reject
        
      }
      
      power_t     <- power_t_value/M
      power_Welch <- power_Welch_value/M
      power_Z     <- power_Z_value/M
      power_Yuen  <- power_Yuen_value/M
      
      power_MWU   <- power_MWU_value/M
      power_WBR   <- power_WBR_value/M
      power_VW    <- power_VW_value/M
      
      power_Yeo_t     <- power_Yeo_t_value/M
      power_Yeo_Welch <- power_Yeo_Welch_value/M
      power_Yeo_Z     <- power_Yeo_Z_value/M
      power_Yeo_Yuen  <- power_Yeo_Yuen_value/M
      
      
      
      power_results <- rbind(power_results, c(mu ,sigma, m, lambda, n1, n2, power_t, power_Welch, power_Z, power_Yuen, power_MWU, power_WBR, power_VW, power_Yeo_t, power_Yeo_Welch, power_Yeo_Z, power_Yeo_Yuen))
    }
  }
}
colnames(power_results) <- c("mu", "sigma", "m", "lambda", "n1", "n2", "t-test", "Welch t-test", "Z-test","YW", "MWU", "WBR", "VW", "Yeo t", "Yeo Welch", "Yeo Z", "Yeo YW")
power_results

# Save the final data frame to a CSV file
write.csv(power_results, file = "D:/Thesis/Thesis_amm/thesis and peper/JSCT_manuscript/Result_151169/power_outlier2_n1!=n2.csv", row.names = FALSE)

        