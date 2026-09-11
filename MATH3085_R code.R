#install.packages("survival")

# Q1(a):
library(survival)
two_times <- c(85, 32, 38, 45, 4, 84, 49, 180, 87, 75, 102, 39, 12, 11, 80, 35, 6)
two_status <- c(1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)  
triple_times <- c(22, 2, 48, 85, 160, 238, 56, 94, 51, 12, 171, 80, 180, 4, 90, 180, 3)
triple_status <- c(1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1)  
data_q1 <- data.frame(
  time = c(two_times, triple_times),
  status = c(two_status, triple_status),
  treatment = factor(c(rep("Two-drug", length(two_times)), rep("Triple-drug", length(triple_times))),
                     levels = c("Two-drug", "Triple-drug"))
)
data_q1$treatment <- relevel(data_q1$treatment, ref = "Two-drug")
km_fit <- survfit(Surv(time, status) ~ treatment, data = data_q1)
plot(km_fit, 
     lty = c(1, 2),  
     xlab = "Time (days)", 
     ylab = "Survival Probability", 
     main = "Kaplan-Meier Survival Curves by Treatment")
legend("topright", 
       legend = c("Two-drug", "Triple-drug"), 
       lty = c(1, 2))

# Q1(b)
summary_120 <- summary(km_fit, times = 120)
two_drug_ci <- summary_120$lower[1] 
two_drug_survival <- summary_120$surv[1]  
two_drug_upper <- summary_120$upper[1]  
triple_drug_ci <- summary_120$lower[2]
triple_drug_survival <- summary_120$surv[2]
triple_drug_upper <- summary_120$upper[2]
cat("Two-drug: 95% CI for S(120) = [", two_drug_ci, ", ", two_drug_upper, "]\n")
cat("Triple-drug: 95% CI for S(120) = [", triple_drug_ci, ", ", triple_drug_upper, "]\n")

# Q1(c)
cox_model <- coxph(Surv(time, status) ~ treatment, data = data_q1)
summary(cox_model)
plot(km_fit, lty = c(1, 2), col = "black", xlab = "Time (days)", 
     ylab = "Survival Probability", main = "KM and PH Survival Curves by Treatment")
ph_two <- survfit(cox_model, newdata = data.frame(treatment = "Two-drug"))
lines(ph_two, lty = 1, col = "red", conf.int = FALSE)
ph_triple <- survfit(cox_model, newdata = data.frame(treatment = "Triple-drug"))
lines(ph_triple, lty = 2, col = "red", conf.int = FALSE)
legend("topright", legend = c("KM Two-drug", "KM Triple-drug", "PH Two-drug", "PH Triple-drug"), 
       lty = c(1, 2, 1, 2), col = c("black", "black", "red", "red"))

# Q1(d)
logmlog <- function(x) { log(-log(x)) }
plot(km_fit, fun = logmlog, conf.int = FALSE,
     xlab = "Time (days)", ylab = "log(-log S(t))",
     main = "Complementary Log-Log Plot for PH Assumption",
     lty = c(1, 2))
legend("bottomright", legend = c("Two-drug", "Triple-drug"), lty = c(1, 2))






# Q2(a)
drs <- read.csv("DRS.csv")
full_model_2a <- coxph(Surv(time, status) ~ treatment + age + riskgp, data = drs)
summary(full_model_2a)
back_model_1_2a <- coxph(Surv(time, status) ~ treatment + riskgp, data = drs)
summary(back_model_1_2a)
back_model_2_2a <- coxph(Surv(time, status) ~ age + treatment, data = drs)
summary(back_model_2_2a)
final_model_2a <- step(full_model_2a, direction = "backward")
summary(final_model_2a)

# Q2(b)
event_times <- drs$time[drs$status == 1]
counts <- table(event_times)
tied_list <- counts[counts > 1]
print(tied_list)
mod_exact_2b <- coxph(final_model_2a$formula, data = drs, ties = "exact")
mod_efron_2b <- coxph(final_model_2a$formula, data = drs, ties = "efron")
summary(mod_exact_2b)
summary(mod_efron_2b)

# Q2(c)
drs$treatment <- factor(drs$treatment, levels = c("laser", "placebo"))  
aft_model <- survreg(Surv(time, status) ~ treatment + riskgp, data = drs, dist = "weibull")
summary(aft_model)  
confint(aft_model)[c("treatmentplacebo", "riskgp"), ]  

scaled <- drs$time / exp(predict(aft_model, type = "lp"))
S <- survfit(Surv(scaled, drs$status) ~ 1)
plot(log(S$time), log(-log(S$surv)), 
     xlab = "log t", ylab = "log(-log S(t))", 
     main = "Check for Weibull Baseline Distribution")
abline(lm(log(-log(S$surv)) ~ log(S$time)))