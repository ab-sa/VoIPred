library(glmnet)
library(VoIPred)

settings <- list()
settings$master_formula <- day30 ~ sex + age + dia + miloc + pmi + htn + smk + kill + tx - 1
settings$default_th <- 0.1
settings$n_sim <- 500 #if 0 wont do this part
settings$subsample <- 1000
settings$auc_n_sim <- 0   #If set to 0, it will not calculate AUC with optimism correction with the same n_sim.
settings$sample_size_n_sim_outer <- 0 #if set to 0 will not do
settings$sample_size_n_sim_inner <- 100 #Voi calculations for each point witin each iteration
settings$sample_sizes <- c(250, 500, 1000, 2000, 4000, 8000, 16000, 32000, Inf)
# settings$sampleInfo_size <- c(100, 200, 400, 800, 1500)
settings$sampleInfo_size <- c(100, 500, 5000, 10000)

set.seed(1234)
results <- list()

data("gusto")
gusto$kill <- (as.numeric(gusto$Killip)>1)*1
master_formula <- settings$master_formula
sample <- gusto[sample(1:dim(gusto)[1],settings$subsample,F),]
x <- model.matrix(master_formula, sample)
y <- sample$day30

cv_reg <- cv.glmnet(x, y, family="binomial")
reg <- glmnet(x, y, family = "binomial", lambda = cv_reg$lambda.min)
results$reg_obj <- reg

# results$res0 <- list()
# results$res1 <- list()


results$res0 <- voi.glmnet(reg_obj = reg, x, y, n_sim = settings$n_sim, Bayesian_bootstrap = F,
                           sampleInfo_size = settings$sampleInfo_size,  type.measure = "default")
results$res1 <- voi.glmnet(reg, x, y, n_sim = settings$n_sim, Bayesian_bootstrap = T,
                           sampleInfo_size = settings$sampleInfo_size)

names(results$res0)


res_def <- voi.glmnet(reg_obj = reg, x, y, n_sim = settings$n_sim, Bayesian_bootstrap = T,
                      sampleInfo_size = settings$sampleInfo_size, type.measure = "default", LASSO = "both")
res_def_PBs <- voi.glmnet(reg_obj = reg, x, y, n_sim = settings$n_sim, Bayesian_bootstrap = F, ParamBootS = T,
                          sampleInfo_size = settings$sampleInfo_size, type.measure = "default", LASSO = "both")
# res_def_VS <- voi.glmnet(reg_obj = reg, x, y, n_sim = settings$n_sim, Bayesian_bootstrap = F,
#                          sampleInfo_size = settings$sampleInfo_size, type.measure = "default", LASSO = "VS")
# res_def_sample <- voi.glmnet(reg_obj = reg, x, y, n_sim = settings$n_sim, Bayesian_bootstrap = F,
#                               sampleInfo_size = settings$sampleInfo_size, type.measure = "default", LASSO = "sample")
res_def_none <- voi.glmnet(reg_obj = reg, x, y, n_sim = settings$n_sim, Bayesian_bootstrap = F,
                              sampleInfo_size = settings$sampleInfo_size, type.measure = "default", LASSO = "none")
res_def_none_PBs <- voi.glmnet(reg_obj = reg, x, y, n_sim = settings$n_sim, Bayesian_bootstrap = F, ParamBootS = T,
                               sampleInfo_size = settings$sampleInfo_size, type.measure = "default", LASSO = "none")


res_mae <- voi.glmnet(reg_obj = reg, x, y, n_sim = settings$n_sim, Bayesian_bootstrap = F,
                      sampleInfo_size = settings$sampleInfo_size, type.measure = "mae")
res_class <- voi.glmnet(reg_obj = reg, x, y, n_sim = settings$n_sim, Bayesian_bootstrap = F,
                        sampleInfo_size = settings$sampleInfo_size, type.measure = "class", LASSO = "both")
res_dev <- voi.glmnet(reg_obj = reg, x, y, n_sim = settings$n_sim, Bayesian_bootstrap = F,
                      sampleInfo_size = settings$sampleInfo_size, type.measure = "deviance")


plot(res_def$lambda, res_def$voi, type = "l", xlab = "Threshold", ylab = "NB diff",
     main = paste0("# of sim: ", settings$n_sim))
lines(res_def$lambda, res_def$vosi[1, ], col = "red")
lines(res_def$lambda, res_def$vosi[2, ], col = "blue")
lines(res_def$lambda, res_def$vosi[3, ], col = "orange")
lines(res_def$lambda, res_def$vosi[4, ], col = "green")
legend("topright", c("EVPI", "EVSI-100", "EVSI-500", "EVSI-5000", "EVSI-10000"),
       lty = rep(1, 5), col = c("black",  "red", "blue", "orange", "green"))

plot(res_def_PBs$lambda, res_def_PBs$voi, type = "l", xlab = "Threshold", ylab = "NB diff")
lines(res_def_PBs$lambda, res_def_PBs$vosi[1, ], col = "red")
lines(res_def_PBs$lambda, res_def_PBs$vosi[2, ], col = "blue")
lines(res_def_PBs$lambda, res_def_PBs$vosi[3, ], col = "orange")
lines(res_def_PBs$lambda, res_def_PBs$vosi[4, ], col = "green")
legend("topright", c("EVPI", "EVSI-100", "EVSI-500", "EVSI-5000", "EVSI-10000"),
       lty = rep(1, 5), col = c("black",  "red", "blue", "orange", "green"))



plot(res_def_none$lambda, res_def_none$voi, type = "l", xlab = "Threshold", ylab = "NB diff")
lines(res_def_none$lambda, res_def_none$vosi[1, ], col = "red")
lines(res_def_none$lambda, res_def_none$vosi[2, ], col = "blue")
lines(res_def_none$lambda, res_def_none$vosi[3, ], col = "orange")
lines(res_def_none$lambda, res_def_none$vosi[4, ], col = "green")
legend("topright", c("EVPI", "EVSI-100", "EVSI-500", "EVSI-5000", "EVSI-10000"),
       lty = rep(1, 5), col = c("black",  "red", "blue", "orange", "green"))

plot(res_def_none_PBs$lambda, res_def_none_PBs$voi, type = "l", xlab = "Threshold", ylab = "NB diff")
lines(res_def_none_PBs$lambda, res_def_none_PBs$vosi[1, ], col = "red")
lines(res_def_none_PBs$lambda, res_def_none_PBs$vosi[2, ], col = "blue")
lines(res_def_none_PBs$lambda, res_def_none_PBs$vosi[3, ], col = "orange")
lines(res_def_none_PBs$lambda, res_def_none_PBs$vosi[4, ], col = "green")
legend("topright", c("EVPI", "EVSI-100", "EVSI-500", "EVSI-5000", "EVSI-10000"),
       lty = rep(1, 5), col = c("black",  "red", "blue", "orange", "green"))











plot(res_def_VS$lambda, res_def_VS$voi, type = "l", xlab = "Threshold", ylab = "NB diff")
lines(res_def_VS$lambda, res_def_VS$vosi[1, ], col = "red")
lines(res_def_VS$lambda, res_def_VS$vosi[2, ], col = "blue")
lines(res_def_VS$lambda, res_def_VS$vosi[3, ], col = "orange")
lines(res_def_VS$lambda, res_def_VS$vosi[4, ], col = "green")
legend("topright", c("EVPI", "EVSI-100", "EVSI-500", "EVSI-5000", "EVSI-10000"),
       lty = rep(1, 5), col = c("black",  "red", "blue", "orange", "green"))

plot(res_def_sample$lambda, res_def_sample$voi, type = "l", xlab = "Threshold", ylab = "NB diff")
lines(res_def_sample$lambda, res_def_sample$vosi[1, ], col = "red")
lines(res_def_sample$lambda, res_def_sample$vosi[2, ], col = "blue")
lines(res_def_sample$lambda, res_def_sample$vosi[3, ], col = "orange")
lines(res_def_sample$lambda, res_def_sample$vosi[4, ], col = "green")
legend("topright", c("EVPI", "EVSI-100", "EVSI-500", "EVSI-5000", "EVSI-10000"),
       lty = rep(1, 5), col = c("black",  "red", "blue", "orange", "green"))

plot(res_def_none$lambda, res_def_none$voi, type = "l", xlab = "Threshold", ylab = "NB diff")
lines(res_def_none$lambda, res_def_none$vosi[1, ], col = "red")
lines(res_def_none$lambda, res_def_none$vosi[2, ], col = "blue")
lines(res_def_none$lambda, res_def_none$vosi[3, ], col = "orange")
lines(res_def_none$lambda, res_def_none$vosi[4, ], col = "green")
legend("topright", c("EVPI", "EVSI-100", "EVSI-500", "EVSI-5000", "EVSI-10000"),
       lty = rep(1, 5), col = c("black",  "red", "blue", "orange", "green"))


cbind(res_def_none$lambda,
      res_def_none$voi,
      t(res_def_none$vosi))

NBs <-
  data.frame(lambda = res_def_none$lambda,
             all = res_def_none$NB_all,
             model =res_def_none$NB_model,
             max = res_def_none$NB_max,
             bs2= res_def_none$NB_max_merged[3, ])

mean(NBs$max >= NBs$bs2)
mean(NBs$max >= NBs$model)

which(NBs$max < NBs$bs2)

NBs[c(1 : 10) , ]


plot(res_mae$lambda, res_mae$voi, type = "l")
lines(res_mae$lambda, res_mae$vosi[1, ], col = "red")
lines(res_mae$lambda, res_mae$vosi[2, ], col = "blue")

plot(res_class$lambda, res_class$voi, type = "l")
lines(res_class$lambda, res_class$vosi[1, ], col = "red")
lines(res_class$lambda, res_class$vosi[2, ], col = "blue")





plot(results$res0$lambda, results$res0$voi, type = "l")
lines(results$res0$lambda, results$res0$vosi[1, ], col = "red")
lines(results$res0$lambda, results$res0$vosi[1, ], col = "blue")


plot(x = res_def$NB_max,
     y = res_def$NB_max_merged[1, ])
abline(a  = 0, b = 1)


plot(results$res0$voi, results$res0$vosi[1, ])
abline(a  = 0, b = 1)
points(results$res0$voi, results$res0$vosi[2, ], col = "red")
points(results$res0$voi, results$res0$vosi[3, ], col = "blue")
points(results$res0$voi, results$res0$vosi[4, ], col = "green")
points(results$res0$voi, results$res0$vosi[5, ], col = "orange")

plot(results$res1[ , c(2, 9)])
abline(a  = 0, b = 1)


