#' @import mvtnorm
#' @import glmnet

aux <- new.env()

#' @export
get_aux<-function()
{
  return(aux)
}


bootstrap <- function (n, Bayesian=F)
{
  if(Bayesian)
  {
    u <- c(0, sort(runif(n - 1)), 1)
    return((u[-1] - u[-length(u)]) * n)
  }
  else
  {
    u <- rmultinom(1, n, rep(1 / n, n))
    return(u[ , 1])
  }
}



# bootstrap: 0=no (parametric); 1=ordinary; 2=Bayesian
#' @export
voi.glm<-function(reg_obj, n_sim=1000, bootstrap=0, lambdas=(1:99)/100)
{
  sample_size <- dim(reg_obj$data)[1]
  mu <- coefficients(reg_obj)
  sigma <- vcov(reg_obj)

  pi<-predict(reg_obj,type="response")
  aux$pi <- pi

  NB_test <-  NB_all <-  NB_max  <- rep(0, length(lambdas))

  for(i in 1:n_sim)
  {
    if(bootstrap>0)
    {
      ws <- bootstrap(sample_size, bootstrap-1)
      bs_data <- cbind(reg_obj$data, ws=ws)
      bs_reg <- glm(data=bs_data, formula = reg_obj$formula, family=reg_obj$family, weights = ws)
      p <- as.vector(predict(bs_reg, newdata=reg_obj$data, type="response"))
    }
    else
    {
      new_betas <- rmvnorm(1, mean=mu, sigma=sigma)
      p <- as.vector(1/(1+exp(-(new_betas %*% t(model.matrix(reg_obj))))))
    }

    if(i==1) plot(pi,p)

    for(j in 1:length(lambdas))
    {
      NB_test[j] <- NB_test[j] + mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * (pi > lambdas[j]))

      NB_all[j] <- NB_all[j] + mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * 1)

      NB_max[j] <- NB_max[j] + mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * (p > lambdas[j]))
    }
    #cat('.')
  }

  voi <- (NB_max-pmax(0,NB_test,NB_all))/n_sim

  res <-cbind(lambda=lambdas, voi=voi, NB_all=NB_all/n_sim, NB_test=NB_test/n_sim, NB_max=NB_max/n_sim)

  return(res)
}





#' @export
#' @param reg_obj: any object that you can aplpy predict with new data to get predictions
#' @param x: The model matrix of predictors
#' @param y: The vector of responses
voi.glmnet <- function(reg_obj, x, y, n_sim=1000, lambdas=(1:99)/100, Bayesian_bootstrap=F, empirical=F,
                       sampleInfo_size = NULL, type.measure = "default", LASSO = "both", ParamBootS = FALSE)
{
  aux$coeffs <- t(as.matrix(coefficients(reg_obj)))
  aux$x <- x
  aux$y <- y
  aux$reg_obj <- reg_obj

  sample_size <- dim(x)[1]

  pi <- predict(reg_obj, type="response", newx=x)
  aux$pi <- pi

  NB_model <- NB_all <- NB_max <- rep(0, length(lambdas))

  dc_model <- NB_model
  dc_all <- NB_model

  optimism <- NB_model

  aux$bs_coeffs <- matrix(NA,nrow=n_sim, ncol=dim(coefficients(reg_obj)))
  colnames(aux$bs_coeffs) <- rownames(coefficients(reg_obj))

  if (! is.null(sampleInfo_size)) {
    NB_max_merged <- matrix(0, nrow = length(sampleInfo_size), ncol = length(lambdas))
    aux$bs_coeffs_bs2 <- array(NA, dim = c(length(sampleInfo_size),
                                           n_sim, nrow(coefficients(reg_obj))))
    dimnames(aux$bs_coeffs_bs2)[[1]] <- sampleInfo_size
    dimnames(aux$bs_coeffs_bs2)[[2]] <- c(1 : n_sim)
    dimnames(aux$bs_coeffs_bs2)[[3]] <- rownames(coefficients(reg_obj))

    p_bs2 <- list()
    pi_bs2 <- list()
  }


  for (j in 1:length(lambdas))
  {
    dc_model[j] <- dc_model[j] + mean((y - (1 - y) * lambdas[j] / (1 - lambdas[j])) * (pi > lambdas[j]))
    dc_all[j] <- dc_all[j] + mean((y - (1 - y) * lambdas[j] / (1 - lambdas[j])) * 1)
  }

  pb <- txtProgressBar(min = 0, max = n_sim, initial = 0)
  # i <- 1
  for(i in 1 : n_sim)
  {
    # if (n_sim <= 100) cat(".") else if (trunc(i / trunc(n_sim / 100)) == i / trunc(n_sim / 100)) cat(".")
    setTxtProgressBar(pb, i)

    if (ParamBootS) {
      y <- rbinom(n = length(y), size = 1, prob = pi)
    }

    weights <- bootstrap(sample_size, Bayesian = Bayesian_bootstrap)
    # set.seed(1)
    if (LASSO %in% c("both", "VS")) {
      tmp <- cv.glmnet(x=x, y=y, family="binomial", weights = weights, type.measure = type.measure)
      bs_reg <- glmnet(x=x, y=y, family="binomial", lambda = tmp$lambda.min, weights = weights)
    }
    if (LASSO %in% c("sample", "none")) {
      bs_reg <- glmnet(x=x, y=y, family="binomial", lambda = 0, weights = weights)
    }

    if (LASSO == "VS") {
      bs_reg_noShrink <- glmnet(x = x[ , rownames(bs_reg$beta)[bs_reg$beta[, 1] != 0]], y = y,
                                family="binomial", lambda = 0, weights = weights)
      aux$bs_coeffs[i,] <- 0
      aux$bs_coeffs[i , c("(Intercept)", rownames(bs_reg$beta)[bs_reg$beta[, 1] != 0])] <- t(as.matrix(coefficients(bs_reg_noShrink)))
      p <- as.vector(predict(bs_reg_noShrink, newx=x[ , rownames(bs_reg$beta)[bs_reg$beta[, 1] != 0]], type="response"))
    }
    else {
      aux$bs_coeffs[i,] <- t(as.matrix(coefficients(bs_reg)))
      p <- as.vector(predict(bs_reg, newx=x, type="response"))
    }


    if (! is.null(sampleInfo_size)) {

      # j = 4
      for (j in 1 : length(sampleInfo_size)) {

        bs2_ind <- sample(c(1 : nrow(x)), size = sampleInfo_size[j], replace = TRUE, prob = weights)
        x_merged <- rbind(x, x[bs2_ind, ])
        y_merged <- c(y, y[bs2_ind])
        p_bs2[[j]] <- c(p, p[bs2_ind])

        sample_size_bs2 <- dim(x_merged)[1]
        # set.seed(1)
        if (LASSO %in% c("both", "sample", "VS")) {
          tmp_bs2 <- cv.glmnet(x = x_merged, y = y_merged, family = "binomial", type.measure = type.measure)
          # weights = as.vector(weights_bs2))
          bs_reg_bs2 <- glmnet(x = x_merged, y = y_merged, family = "binomial",
                               lambda = tmp_bs2$lambda.min)
          # weights = as.vector(weights_bs2))

          cbind(pop = coefficients(bs_reg),
                sample = coefficients(bs_reg_bs2))

        }
        else {
          bs_reg_bs2 <- glmnet(x = x_merged, y = y_merged, family = "binomial", lambda = 0)
        }

        aux$bs_coeffs_bs2[j, i, ] <- t(as.matrix(coefficients(bs_reg_bs2)))
        pi_bs2[[j]] <- as.vector(predict(bs_reg_bs2, newx = x_merged, type = "response"))

      }

    }


    if(i==1) plot(pi,p)


    ## (p - lambda) / (1 - lambda)


    # j = 1
    for(j in 1:length(lambdas))
    {
      NB_model[j] <- NB_model[j] + mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * (pi[ , 1] > lambdas[j]))
      NB_all[j] <- NB_all[j] + mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * 1)
      NB_max[j] <- NB_max[j] + mean((p - (1 - p) * lambdas[j] / (1 - lambdas[j])) * (p > lambdas[j]))

      if (! is.null(sampleInfo_size)) {

        for (k in 1 : length(sampleInfo_size)) {

          NB_max_merged[k, j] <- NB_max_merged[k, j] +
            mean((p_bs2[[k]] - (1 - p_bs2[[k]]) * lambdas[j] / (1 - lambdas[j])) * (pi_bs2[[k]] > lambdas[j]))

        }

      }

      dc_model_int <- sum((y + (1 - y) * lambdas[j] / (1 - lambdas[j])) * (p > lambdas[j]) * weights) / sum(weights)
      dc_model_ext <- mean((y + (1 - y) * lambdas[j] / (1 - lambdas[j])) * (p > lambdas[j]))
      optimism[j] <- optimism[j] + dc_model_int - dc_model_ext
    }
  }

  NB_model <- NB_model / n_sim
  NB_all <- NB_all / n_sim
  NB_max <- NB_max / n_sim
  optimism <- optimism / n_sim

  voi <- (NB_max-pmax(0, NB_model, NB_all))

  res <- list(lambda=lambdas, voi=voi, NB_all=NB_all, NB_model=NB_model, NB_max=NB_max, dc_model=dc_model, dc_all=dc_all, optimism=optimism)

  if (! is.null(sampleInfo_size)) {
    NB_max_merged <- NB_max_merged / n_sim
    vosi <- t(apply(NB_max_merged, 1, function(x) (x - pmax(0, NB_model, NB_all))))
    res$vosi <- vosi
    res$NB_max_merged <- NB_max_merged
  }

  return(res)
}









#' @export
process_results <- function(res, graphs=c("voi","summit","dc"))
{
  out <- list()
  out$inb_current <- mean(pmax(0,res[,'NB_model'],res[,'NB_all'])-pmax(0,res[,'NB_all']))
  out$inb_perfect <- mean(res[,'NB_max']-pmax(0,res[,'NB_all']))
  out$voi_r <- out$inb_perfect/out$inb_current

  if(!is.na(match("voi",graphs)))
  {
    plot(res[,'lambda'], res[,'NB_max']-pmax(0,res[,'NB_model'],res[,'NB_all']),type='l', lwd=2, col="red", xlab="Threshold", ylab="EVPI")
  }

  if(!is.na(match("summit",graphs)))
  {
    plot(res[,'lambda'],pmax(0,res[,'NB_model'],res[,'NB_all'])-pmax(0,res[,'NB_all']),type='l', xlab='Threshold', ylab='Incremental net benefit', lwd=2)
    lines(res[,'lambda'],res[,'NB_max']-pmax(0,res[,'NB_all']),type='l',col="red", lwd=2)
  }

  if(!is.na(match("dc",graphs)))
  {
    max_y <- max(res[,'NB_all'],res[,'NB_model'])
    plot(res[,'lambda'],res[,'NB_model'],type='l', xlab='Threshold', ylab='Net benefit', lwd=2, xlim=c(0,1), ylim=c(0,max_y), col="red")
    lines(res[,'lambda'],res[,'lambda']*0,type='l', lwd=1, col="gray")
    lines(res[,'lambda'],res[,'NB_all'],type='l', lwd=1, col="black")
    lines(res[,'lambda'],res[,'NB_max'],type='l',col="blue", lwd=2)
  }

  return(out)
}




#' @export
decision_curve <- function(y, pi, lambdas=(1:99)/100)
{
  NB_model <- rep(0, length(lambdas))
  NB_all <- NB_model

  for(j in 1:length(lambdas))
  {
    NB_model[j] <- mean((y - (1 - y) * lambdas[j] / (1 - lambdas[j])) * (pi > lambdas[j]))
    NB_all[j] <- mean((y - (1 - y) * lambdas[j] / (1 - lambdas[j])) * 1)
  }

  return(cbind(lambda=lambdas, NB_none=NB_all*0, NB_model=NB_model, NB_all=NB_all))
}



#' @export
plot_decision_curve <- function(dc_data)
{
  max_y <-max(dc_data[,c('NB_none','NB_model','NB_all')])
  plot(dc_data[,'lambda'],dc_data[,'NB_none'],type='l',ylim=c(0,max_y),xlim=c(0,1),col='gray', xlab="Threshold", ylab="Net benefit")
  lines(dc_data[,'lambda'],dc_data[,'NB_all'],type='l',ylim=c(0,max_y),col='black')
  lines(dc_data[,'lambda'],dc_data[,'NB_model'],type='l',ylim=c(0,max_y),col='red',lw=2)
}
