#*********************************************************************
# This is free software and comes with ABSOLUTELY NO WARRANTY.
# You are welcome to redistribute it under GPL conditions.
# Copyright (C) Andrea Dal Pozzolo, March 2015
# License: GPL (>= 2)                                           
#*********************************************************************

# #' @title warpingWrapper: Wrapper to warpingUnder function
# #'
# #' @description
# #' Calls the function \code{\link{warpingUnder}} using a classisification algorithm and an unbalanced dataset
# #'
# #' @param algo classification algorithm supported by mlr package
# #' @param dataset dataset used for the analysis
# #' @param nCV number of repetation of the Cross Validation (CV)
# #' @param B  number of models to create for each fold of the CV
# #' @param nFold number of folds of the CV
# #' @param ncore number of cores to use in multicore computation
# #' @param ... extra parameters passed to mlr train function
# #' 
# #' @return The function returns a list: 
# #' \item{results}{results of undersampling}
# #' \item{probs}{probabilities}
# #' \item{mres}{mres}
# #' \item{meltResALL}{meltResALL}
# #' 
# #' @examples 
# #' library(warping)
# #' warpingWrapper("randomForest", "pima", nCV=50, nFold=3)
# #' 
# #' 
# #' @export
# warpingWrapper <- function(algo="randomForest", dataset="pima", nCV=10, B=1, nFold=100, ncore=1, ...){
#   stopifnot(length(algo)==1, length(dataset)==1, nCV > 0)
#   
#   urlDataset <- "http://www.ulb.ac.be/di/map/adalpozz/imbalanced-datasets/"
#   load(url(paste0(urlDataset, dataset,".Rdata")))
#   
#   cat("dataset:", dataset, "\t algo:", algo, "\n")
#   d <- ncol(data)
#   colnames(data) <- c(paste0("V", 1:(d-1)), "Y")
#   
#   results <- warpingUnder(Y ~., data, algo, task.id=paste(algo, dataset, sep="_"), 
#                           positive=1, costs=NULL, verbose=TRUE, dirPlot="plot", 
#                           ncore, nCV, B, nFold, ...)     
#   
#   retunr(results)  
# }



#' @title Brier Score
#'
#' @description
#' Function used to compute Brier Score, a measure of probability calibration.
#'
#' @param phat probability of a instance to belogn to the class positive
#' @param truth class of the instance
#' @param positive value of the positive (minority) class
#' 
#' @return Brier Score statistics 
#' 
#' @examples 
#' truth <- c(1,1,1,0,0,0,1,0,1,0)
#' phat <- c(0.9,0.8,0.4,0.5,0.3,0.2,0.8,0.3,0.8,0.3)
#' BS(phat, truth, 1)
#' 
#' @export
BS <- function(phat, truth, positive) {
  stopifnot(length(phat) == length(truth), positive %in% truth)
  y <- as.numeric(truth == positive)
  mean((y - phat)^2)
}  



#***********************
# The following code is used to asses the wrapping effect 
# of undersampling on the posterior probability
# use a large number of folds to make sure we use almost all the dataset for training
#***********************
#' @title warpingUnder function
#'
#' @description
#' Function used to asses the wrapping effect of undersampling on the posterior probability.
#'  use a large number of folds to make sure we use almost all the dataset for training
#'
#' @param formula A formula of the form y ~ x1 + x2 + ...
#' @param data Data frame from which variables specified in formula are preferentially to be taken
#' @param algo classification algorithm supported by mlr package
#' @param task.id name of the task
#' @param positive value of the positive (minority) class
#' @param costs cost matrix
#' @param nCV number of repetation of the Cross Validation (CV)
#' @param B  number of models to create for each fold of the CV
#' @param nFold number of folds of the CV
#' @param ncore number of cores to use in multicore computation
#' @param dirPlot directory where to save plots (set dirPlot=NA to avoid plots)
#' @param verbose print extra information (logical variable)
#' @param ... extra parameters passed to mlr train function
#' 
#' @return The function returns a list: 
#' \item{results}{results of undersampling}
#' \item{probs}{probabilities}
#' \item{mres}{mres}
#' \item{meltResALL}{meltResALL}
#' 
#' @examples 
#' library(mlbench) 
#' data(Ionosphere)
#' library(warping)
#' res <- warpingUnder(Class ~., Ionosphere, "randomForest", task.id="rf_Ionosphere", positive="bad", nCV=3, B=1, nFold=5)
#' 
#' @export
warpingUnder <- function(formula, data, algo, task.id="cv", positive=1, costs=NULL, nCV=10, B=1, nFold=100, ncore = 1, dirPlot=NA, verbose=TRUE, ...){
  
  
  stopifnot(class(formula) == "formula", 
            NCOL(data)>1, NROW(data)>1, 
            is.logical(verbose), ncore > 0, 
            nCV >= 0, B > 0, nFold > 1)
  
  #derivative of ps w.r.t. p, where ps (p) is the posterior with (without) undersampling
  dps <- function(p, beta) beta / (p +beta * (1-p))^2
  
  #compute the sd of the posterior probabilities for each instance across all the repeated CV
  probSd <- function(pred.data, positive) {
    #library(plyr)
    #rename the column to be used with ddply
    colnames(pred.data)[which(colnames(pred.data)==paste("prob", positive, sep = "."))] <- "prob.1"
    #take the variance of the probability for each instance over all the CV.
    dt <- ddply(pred.data, .(id, truth), summarize, sd.prob = sd(prob.1))
    dt$sd.prob
  }
  
  if (!is.na(dirPlot))
    if (!file.exists(dirPlot))
      dir.create(dirPlot)  
  
  
  target <- as.character(formula[[2]])
  tgt <- which(names(data) == target)
  if(length(tgt) == 0)
    stop("target variable not defined")
  
  #define task
  type <- ifelse(is.factor(data[ ,tgt]), "classif", "regr")
  if(type == "classif")
    task <- makeClassifTask(id=task.id, data=data, target=target, positive=positive)
  else
    stop("target variable must be a factor for classification")
  
  L <- task$task.desc$class.levels
  negative <- setdiff(L, positive)
  if(length(L) > 2)
    stop("only binary classification supported yet")
  
  type <- task$task.desc$type
  cv.lrn <- paste(type, algo, sep=".")
  ## Define the learner: set it up for predicting probabilities
  lrn <- makeLearner(cv.lrn, predict.type = "prob", ...)
  if(verbose)
    print(lrn)
  
  task <- removeConstantFeatures(task)
  
  Y <- data[ ,tgt]
  N <- length(Y)
  N.pos <- sum(Y == positive)
  
  #define a cost matrix if not given using class priors
  if(is.null(costs))
    costs <- matrix(c(cTP=0, cFP=N.pos/N, cFN=(N-N.pos)/N, cTN=0), 2)
  else
    stopifnot(is.matrix(costs), ncol(costs) == 2, nrow(costs) == 2)
  colnames(costs) <- rownames(costs) <- task$task.desc$class.levels
  if(verbose){
    cat("cost matrix \n")
    print(costs)
  }
  
  ## Calculate the theoretical threshold for the positive class
  th <- costs[2,1]/(costs[2,1] + costs[1,2])  #FP cost / (FP cost + FN cost)
  
  metrics <- list(gmean, f1, ppv, tpr) #mlr::auc
  submetrics <- list(gmean, f1)
  
  if(nCV > 1) #change the aggration method in the case of repeated CV
    metrics <- lapply(metrics, function(x)  setAggregation(x, testgroup.mean))  
  
  
  #predict once with all samples
  mod <- train(lrn, task)
  #test on the training set itself
  pred <- predict(mod, newdata=data)
  pred.all <- setThreshold(pred, th)
  res.all <- mlr::performance(pred.all, metrics, task)
  prob.all <- getProbabilities(pred)
  
  
  ## define a 10-fold cross-validation that will be used for all tasks (unbalanced, undersampling, oversampling ..)
  if(nCV == 1)
    rdesc <- makeResampleDesc("CV", iters = nFold, stratify = TRUE)
  if(nCV > 1)
    rdesc <- makeResampleDesc("RepCV", folds = nFold, reps = nCV, stratify = TRUE)
  #   if(nCV == 0)
  #     rdesc <- makeResampleDesc("LOO")
  
  # Create a resample instance based an a task in order to calculate the performance of all learners on the same instances
  rin <- makeResampleInstance(rdesc, task = task)
  
  if(ncore > 1){
    #library("parallelMap")
    parallelStart("multicore", ncore) #Starting parallelization in mode=multicore with cpus=ncore.
  }
  
  #cross-validation without resampling (dataset unbalanced)
  cv <- resample(lrn, task, rin, measures = metrics, models=FALSE, show.info=FALSE)
  # cv <- crossval(lrn, task, iters=10L, stratify=TRUE, models=FALSE, show.info=FALSE)
  
  if(ncore > 1)
    parallelStop() #Stop parallelization
  
  predcv <- cv$pred
  
  #extract the posterior probability of the positive class
  probs <- data.frame(id=factor(predcv$data$id), class=predcv$data$truth, iter=predcv$data$iter, prob.unbal=getProbabilities(predcv))
  pred.th <- setThreshold(predcv, th)
  res.th <- mlr::performance(pred.th, metrics, task)
  if(verbose){
    cat("\n Performance with theoretical threshold", th, "\n")
    print(res.th)
  }
  
  #take the sd of the probability for each instance over all the CV.
  sdUnbal <- probSd(pred.th$data, positive)
  
  #extract 10 values of beta
  betas.under <- betasUnder(data[ ,tgt], positive, 5, "perc") 
  # repeat undersampling for different values of beta
  probs.under <- res.under <- res.ucal <- NULL
  preds.under <- preds.ucal <- list()
  for(b in betas.under){
    dirPlot.b <- ifelse(b == betas.under[1], dirPlot, NA)
    res.b <- undersampling(task, lrn, b, metrics, rin, positive, negative, verbose, dirPlot.b, ncore, B=B)
    preds.under[[length(preds.under) + 1]] <- res.b$pred.under$data
    preds.ucal[[length(preds.ucal) + 1]] <- res.b$pred.cal$data
    prob.b <- res.b$prob[ ,-1]
    colnames(prob.b) <- paste(colnames(prob.b), b, sep="_")
    probs <- cbind(probs, prob.b)
    res.under <- rbind(res.under, res.b$under)
    res.ucal <- rbind(res.ucal, res.b$cal)
  }
  row.names(res.under) <- paste("res.under.b", betas.under, sep="-")
  row.names(res.ucal) <- paste("res.ucal.b", betas.under, sep="-")  
  names(preds.under) <- paste("pred.under.b", betas.under, sep="-")  
  names(preds.ucal) <- paste("pred.ucal.b", betas.under, sep="-")  
  sdUnder <- sapply(preds.under, function(x) probSd(x, positive))
  sdUcal <- sapply(preds.ucal, function(x) probSd(x, positive))
  
  sdRatio <- sdUnder/sdUnbal
  #compute the probability of having dps > sdRatio using the probability obtained without CV
  #higher <- matrix(NA, nrow=length(prob.all), ncol=length(betas.under))
  dpsMat <- matrix(NA, nrow=length(prob.all), ncol=length(betas.under))
  for(b in 1:length(betas.under)){
    dpsMat[ ,b] <- dps(prob.all, betas.under[b])
  }
  higher <- dpsMat >= sdRatio
  probHigher <- apply(higher, 2, function(x) mean(x, na.rm=TRUE))
  cat("\n probability of the derivative being higher \n")
  names(probHigher) <- sub("pred.under.b", "beta", names(probHigher))
  print(probHigher)
  #write.csv(probHigher, file="probHigher.csv")
  
  resBoth <- resLow <- resHigh <- NULL
  for(b in 1:length(betas.under)){
    higher.b <- higher[ ,b]
    preds.under.b <- preds.under[[b]]
    id.prob <- which(colnames(preds.under.b) == paste("prob", positive, sep="."))
    y <- factor(preds.under.b$truth == positive, levels=c(TRUE, FALSE), labels=c(1, 0))
    
    #to do compute the performances separately for each iteration
    
    id.high <- which(higher.b)
    if(length(id.high) > 0){
      pred.under.high <- subset(preds.under.b, id %in% id.high)
      pred.unbal.high <- subset(pred.th$data, id %in% id.high)
      
      dd.under <- res.b$pred.under
      dd.under$data <- pred.under.high
      res.under.high <- mlr::performance(dd.under, metrics, task.under)
      dd.unbal <- pred.th
      dd.under$data <- pred.unbal.high
      res.unbal.high <- mlr::performance(dd.under, metrics, task)
      res.high <- data.frame(rbind(res.unbal.high, res.under.high))
      res.high$type <- c("unbal", "under")
      
      #       y.high <- factor(pred.under.high$truth == positive, levels=c(TRUE, FALSE), labels=c(1, 0))
      #       res.under.high <- probMetrics(pred.under.high[ ,id.prob], y.high, verbose=FALSE)
      #       res.unbal.high <- probMetrics(pred.unbal.high[ ,id.prob], y.high, verbose=FALSE)
      #       res.high <- rbind(data.frame(res.unbal.high, type="unbal"), data.frame(res.under.high, type="under"))
      
      res.high$beta <- betas.under[b]
      resHigh <- rbind(resHigh, res.high)      
    }
    
    id.low <- which(!higher.b)
    if(length(id.low) > 0){
      pred.under.low <- subset(preds.under.b, id %in% id.low)
      pred.unbal.low <- subset(pred.th$data, id %in% id.low) 
      
      dd.under <- res.b$pred.under
      dd.under$data <- pred.under.low
      res.under.low <- mlr::performance(dd.under, metrics, task.under)
      dd.unbal <- pred.th
      dd.under$data <- pred.unbal.low
      res.unbal.low <- mlr::performance(dd.under, metrics, task)
      res.low <- data.frame(rbind(res.unbal.low, res.under.low))
      res.low$type <- c("unbal", "under")
      
      #       y.low <- factor(pred.under.low$truth == positive, levels=c(TRUE, FALSE), labels=c(1, 0))
      #       res.under.low <- probMetrics(pred.under.low[ ,id.prob], y.low, verbose=FALSE)
      #       res.unbal.low <- probMetrics(pred.unbal.low[ ,id.prob], y.low, verbose=FALSE)
      #       res.low <- rbind(data.frame(res.unbal.low, type="unbal"), data.frame(res.under.low, type="under"))
      
      res.low$beta <- betas.under[b]
      resLow <- rbind(resLow, res.low)
    }
    
    dd.under <- res.b$pred.under
    dd.under$data <- preds.under.b
    res.under.both <- mlr::performance(dd.under, metrics, task.under)
    res.unbal.both <- mlr::performance(pred.th, metrics, task)
    res.both <- data.frame(rbind(res.unbal.both, res.under.both))
    res.both$type <- c("unbal", "under")
    
    #     res.under.both <- probMetrics(preds.under.b[ ,id.prob], y, verbose=FALSE)
    #     res.unbal.both <- probMetrics(pred.th$data[ ,id.prob], y, verbose=FALSE)
    #     res.both <- rbind(data.frame(res.unbal.both, type="unbal"), data.frame(res.under.both, type="under"))
    
    res.both$beta <- betas.under[b]
    resBoth <- rbind(resBoth, res.both)
  }
  resALL <- rbind(data.frame(resHigh, high=TRUE), data.frame(resLow, high=FALSE), data.frame(resBoth, high="ALL"))
  resALL$high <- factor(resALL$high)
  #write.csv(resALL, "resALL.csv")
  
  #library(reshape2)
  meltResALL <- reshape2::melt(resALL, id.vars = c("type", "beta", "high"), variable.name = "metric")
  #meltResALL <- melt(resALL, .(type, beta, high), variable.name="metric")
  #save(meltResALL, file="meltResALL.Rdata")
  
  if (!is.na(dirPlot)){
    #plot the density of the posterior probability
    dd <- predcv$data
    colnames(dd)[which(colnames(dd)==paste("prob", positive, sep = "."))] <- "prob.1"
    pp <- as.numeric(quantile(prob.all, c(0.25, 0.5, 0.75)))
    if(all(pp < 0.001))
      pp <- round(pp, 3)
    #library(ggplot2)
    d <- ggplot(dd, aes(x=prob.1, fill=truth)) + geom_density()
    d <- d + labs(list(title = paste0(task.id, '\n'), x = "\n p(+|x)"))
    for(q in pp)
      d <- d + geom_vline(xintercept = q, colour="blue", linetype = "longdash")
    ggsave(filename=paste0(dirPlot, "/dens_beta_1.pdf"), plot=d)
    ggsave(filename=paste0(dirPlot, "/dens_beta_1.svg"), plot=d)
    
    #plot the derivative and sdRatioMean
    sdRatioMean <- apply(sdRatio, 2, function(x) mean(x[is.finite(x)]))
    vardd <- NULL
    for(z in 1:length(pp))
      vardd <- rbind(vardd, data.frame(beta=betas.under, sdRatioMean, p=pp[z], dps=dps(pp[z], betas.under)))
    vardd$p <- factor(vardd$p)
    p <- ggplot(data = vardd, aes(x=beta, y = dps, colour=p))
    p <- p + geom_line() + ggtitle(paste0(task.id, "\n"))
    p <- p + geom_hline(yintercept = 1, colour="blue", linetype = "longdash")
    p <- p + geom_line(aes(x=beta, y = sdRatioMean), colour="black", linetype = "longdash")
    p <- p + ylim(0, max(vardd$dps, vardd$sdRatioMean)) + theme_bw()
    ggsave(filename=paste0(dirPlot, "/dps_",task.id,".pdf"), plot=p)
    ggsave(filename=paste0(dirPlot, "/dps_",task.id,".svg"), plot=p)
    #   library(gridExtra)
    #   pdf(paste0(dirPlot, "/dps_p_",task.id,".pdf"))
    #   grid.arrange(arrangeGrob(p, d, ncol=2))
    #   dev.off()
  }
  
  
  #library(reshape2)
  m.var <- setdiff(colnames(probs), c("class", "id", "iter"))
  mdd <- reshape2::melt(probs, measure.vars = m.var, variable.name = "prob")
  mdd$beta <- as.numeric(sapply(strsplit(as.character(mdd$prob), "_", fixed=TRUE), function(x) x[2]))
  mdd$type <- sapply(strsplit(as.character(mdd$prob), ".", fixed=TRUE), function(x) x[2])
  mdd$type <- factor(sapply(strsplit(mdd$type, "_", fixed=TRUE), function(x) x[1]))
  mdd[which(mdd$type=="unbal"), 'beta'] <- 1
  mdd$lambda <- with(mdd, (sqrt(beta)-beta)/(1-beta))
  
  brier.all <- BS(getProbabilities(pred.all), pred.all$data$truth, positive)
  res.all <- c(res.all, beta=1, pi=N.pos/N, brier=brier.all)
  brier.th <- BS(getProbabilities(pred.th), pred.th$data$truth, positive)
  res.th <- c(res.th, beta=1, pi=N.pos/N, brier=brier.th) 
  
  results <- rbind(res.all, res.th, res.under, res.ucal)
  m <- colnames(results)
  m <- sub("tpr", "recall", m)
  m <- sub("ppv", "precision", m)
  colnames(results) <- m
  #save(results, file="results.Rdata")
  #write.csv(results, file=paste0("results_", task.id, ".csv"))
  
  df <- data.frame(results, prob=row.names(results))
  mres <- reshape2::melt(df, id.vars=c("prob", "beta"), variable.name = "metric")
  mres$type <- sapply(strsplit(as.character(mres$prob), ".", fixed=TRUE), function(x) x[2])
  mres$type <- factor(sapply(strsplit(mres$type, "_", fixed=TRUE), function(x) x[1]))
  #save(mres, file="mres.Rdata")
  
  
  l <- list(results=results, probs=probs, mres=mres, meltResALL=meltResALL)
  return(l)
}






