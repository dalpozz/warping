#*********************************************************************
# This is free software and comes with ABSOLUTELY NO WARRANTY.
# You are welcome to redistribute it under GPL conditions.
# Copyright (C) Andrea Dal Pozzolo, March 2015
# License: GPL (>= 2)                                           
#*********************************************************************



#*********************************************************************
# undersampling the majority class with rate = beta
# beta is the probability of selecting a negative instance during undersampling
# the probability estimate is calibrated to adjust it for class prior change
#*********************************************************************
#' @title Undersampling
#'
#' @description
#' Undersampling the majority class with rate = beta, where beta is the probability of selecting a negative instance during undersampling.
#' The probability estimate is calibrated to adjust it for class prior change
#'
#' @param task classification task to perform (see mlr package)
#' @param lrn learner to use with the task (see mlr package)
#' @param beta probability of selecting a negative instance during undersampling
#' @param metrics metrics to asess classification performances
#' @param rdesc resampling methods to use (see mlr package)
#' @param positive value of the positive (minority) class
#' @param negative value of the negative (majority) class
#' @param verbose print extra information (logical variable)
#' @param ncore number of cores to use in multicore computation
#' @param dirPlot directory where to save plots
#' @param type way to apply undersampling: within or before the CV
#' @param B  number of models to create for each fold of the CV
#' 
#' 
#' @return The function returns a list: 
#' \item{under}{results with undersampling}
#' \item{cal}{results after calibration}
#' \item{pred.under}{predictions with undersampling}
#' \item{pred.cal}{predictions after calibration}
#' \item{prob}{posterior probability of the classifier for the positive class}
#' 
#' @examples #NOT FINISHED YET
#' task <- makeClassifTask(id="class", data=pima, target=Y, positive=1)
#' res <- undersampling(task, lrn, beta, metrics, positive=1, negative=0)
#' 
#' 
#' @export
undersampling <- function(task, lrn, beta, metrics, rdesc=NULL, positive, negative, verbose=TRUE, dirPlot="plot", ncore=1, type="within", B=1){
  
  stopifnot(beta > 0, beta <= 1, is.logical(verbose), ncore > 0, type %in% c("within", "before"), B > 0)
  
  if (!is.na(dirPlot))
    if(!file.exists(dirPlot))
      dir.create(dirPlot)  
  
  submetrics <- list(gmean, f1)
  
  task.id <- task$task.desc$id  
  # undersampling the majority class (negative) class with rate beta
  task.under <- undersample(task, rate = beta)
  #extract the proportion of positives after undersampling
  y <- getTaskTargets(task.under)
  pi.under <- sum(y==positive)/length(y)
  if(pi.under > 0.5)
    cat(" -- WARNING more positive examples than negatives -- \n")
  
  
  if(type == "before"){ #undersampling before the CV
    # the undersampling is done outsite the CV
    lrn.under <- lrn
    task.under <- task.under    
  } else { #undersampling within the CV
    if (B==1) #use a wrapper to perform undersampling before learning a model
      lrn.under <- makeUndersampleWrapper(lrn, usw.rate=beta)
    else{
      #undersampling + bagging within the CV
      lrn <- setPredictType(lrn, "response")
      lrn.under <- makeUnderBaggingWrapper(lrn, ubw.rate=beta, ubw.iters=B,  ubw.maxcl="all")
      lrn.under <- setPredictType(lrn.under, "prob")
    }
    task.under <- task
  }
  
  if(ncore > 1){
    library("parallelMap")
    parallelStart("multicore", ncore) #Starting parallelization in mode=multicore with cpus=ncore.
  }
  
  if (is.null(rdesc))
    cv.under <- crossval(lrn.under, task.under, iters=10L, stratify=TRUE, models=FALSE, show.info=FALSE)
  else
    cv.under <- resample(lrn.under, task.under, rdesc, measures = metrics, models=FALSE, show.info=FALSE)
  
  if(ncore > 1){
    #Stop parallelization
    parallelStop()
  }
  
  #extract predictions from CV
  pred.under <- cv.under$pred
  
  #use the proportion of positives as threshold
  pred.under <- setThreshold(pred.under, pi.under)
  res.under <- mlr::performance(pred.under, metrics, task.under)
  res.under <- c(res.under, beta=beta, pi=pi.under)
  if(verbose){
    cat("\n Performance with undersampling repeated", B, "time using rate", beta, "and threshold", pi.under, "\n")
    print(res.under)
  }
  
  #   pref(getProbabilities(pred.under), pred.under$data$response, text=paste0("under_", beta), dirPlot)
  #   p <- plotThreshVsPerf(pred = pred.under, measures = submetrics, mark.th = pi.under)
  #   ggsave(filename=paste0(dirPlot, "/Perf_", task.id, "_beta_", beta, "_under.pdf"), plot=p, width=10, height=5)
  
  
  #extract posterior probabilities of the positive class after undersampling.
  prob.under <- getProbabilities(pred.under) 
  #calibrate the probability for the change in priors with the testing set (see Elkan paper)
  #we assume that testing and unbalanced training sets have the same priors
  prob.cal <- beta*prob.under / (beta*prob.under - prob.under  +1)
  
  #adjust the threhsold to keep performances of undersampling
  th.cal <- pi.under*beta / (pi.under*(beta - 1 ) + 1)
  pred.cal <- pred.under
  cl.pos <- paste("prob", positive, sep = ".")
  cl.neg <- paste("prob", negative, sep = ".")
  pred.cal$data[ ,which(colnames(pred.cal$data) == cl.pos)] <- prob.cal
  pred.cal$data[ ,which(colnames(pred.cal$data) == cl.neg)] <- 1-prob.cal
  response <- factor(prob.cal >= th.cal, levels=c(TRUE, FALSE), labels=c(positive, negative))
  pred.cal$data$response <- response
  pred.cal$threshold[which(names(pred.cal$threshold) == positive)] <- th.cal
  pred.cal$threshold[which(names(pred.cal$threshold) == negative)] <- 1-th.cal
  
  res.cal <- mlr::performance(pred.cal, metrics, task.under)
  res.cal <- c(res.cal, beta=beta, pi=pi.under)
  if(verbose){
    cat("\n Performance with undersampling after calibration using threshold", th.cal, "\n")
    print(res.cal)
  }
  
  #   pref(prob.cal, response, text=paste0("under_cal_", beta), dirPlot) 
  #   p <- plotThreshVsPerf(pred = pred.cal, measures = submetrics, mark.th = th.cal)
  #   ggsave(filename=paste0("plotThreshVsPerf_", task.id, "_under_cal_", beta, ".pdf"), plot=p, width=10, height=5)
  
  
  #probabiltiy value for which derivative of pSamp wrt p is = 1
  lambda <- (sqrt(beta)-beta)/(1-beta)
  
  probClassUnder <- data.frame(class=pred.under$data$truth, phat.under=prob.under, phat.cal=prob.cal)
  #   library(caret)
  #   calUnder <- calibration(class ~ phat.under + phat.cal, data = probClassUnder) #cuts = 13
  #   p <- plot(calUnder, type = "l", auto.key = list(columns = 3, lines = TRUE, points = FALSE))
  #   pdf(file=paste0("calibration_", beta, ".pdf"), width = 7, height = 7)
  #   print(p)
  #   dev.off()
  
  
  dd <- probClassUnder[order(probClassUnder$phat.under, decreasing=T), ]
  dd$x <- 1:nrow(dd)
  library(reshape2)
  mdd <- reshape2::melt(dd, measure.vars = c("phat.under", "phat.cal"), variable.name = "prob")  
  tit <- paste0("beta = ", beta, " pi.under = ", round(pi.under, 3), "\n")
  
  if(!is.na(dirPlot)){
    #save(mdd, tit, file=paste0(dirPlot, "/mdd_beta_",beta,".Rdata"))  
    library(ggplot2)
    #plot the posterior probability from larger to smaller
    #     p <- ggplot(data = mdd, aes(x=x, y = value, colour=prob))
    #     p <- p + geom_line() + labs(list(title = tit, y = "p(+|x) \n"))
    #     ggsave(filename=paste0(dirPlot, "/phat_beta_",beta,".pdf"), plot=p)
    h <- ggplot(mdd, aes(x=value, fill=prob)) + geom_bar(position="dodge")
    h <- h + facet_grid(.~ class) + theme_bw() + labs(list(title = tit, x = "\n p(+|x)"))
    h <- h + geom_vline(xintercept = lambda, colour="blue", linetype = "longdash")
    ggsave(filename=paste0(dirPlot, "/hist_beta_",beta,".pdf"), plot=h)
    #ggsave(filename=paste0(dirPlot, "/hist_beta_",beta,".svg"), plot=h)
    d <- ggplot(mdd, aes(x=value, fill=prob)) + geom_density()
    d <- d + facet_grid(.~ class) + theme_bw() + labs(list(title = tit, x = "\n p(+|x)"))
    d <- d + geom_vline(xintercept = lambda, colour="blue", linetype = "longdash")
    ggsave(filename=paste0(dirPlot, "/dens_beta_",beta,".pdf"), plot=d)
    #ggsave(filename=paste0(dirPlot, "/dens_beta_",beta,".svg"), plot=d)
  }
  
  return(list(under=res.under, cal=res.cal, pred.under=pred.under, pred.cal=pred.cal, prob=probClassUnder))
  
}


#*********************************************************************
#defines the possible levels of undersampling given the class proporiton
#N: number of values of beta to return (excluding beta=1)
#method: defines the way to compute beta
#*********************************************************************
#' @title betasUnder
#'
#' @description
#' Defines the possible levels of undersampling given the class proporiton
#'
#' @param y response variable
#' @param positive value of the positive (minority) class
#' @param N number of values of beta to return 
#' @param method method to compute beta: perc or prob
#' 
#' 
#' @return values of beta for a give response variable
#' 
#' @examples 
#' y <- rep(c(1, 0, 0, 0, 0, 0), 100)
#' betasUnder(y, positive=1, N = 10, method="perc")
#' @export
betasUnder <- function(y, positive=1, N = 10, method="perc"){
  
  type <- match.arg(method, c("prob", "perc"))
  N.neg <- length(which(y != positive))
  N.pos <- length(which(y == positive))
  stopifnot(N.pos <= N.neg)
  
  if(type == "prob"){
    #compute betas to have 10% 20% .. 90% of negative instances removed in undersampling
    minBeta <- N.pos/N.neg
    #betas <- unique(c(seq(minBeta, 1, by=0.1), 1))
    betas <- seq(minBeta, 1, by = (1-minBeta)/N)
    betas <- round(betas, 2)
    #make sure beta is not zero
    betas[1] <- ifelse(betas[1] == 0, minBeta, betas[1])
    betas <- sort(betas, decreasing=TRUE)
  }
  
  if(type == "perc"){
    #compute betas to have 10% 20% .. 50% of positive instances after undersampling
    pPosMin <- N.pos / (N.pos+N.neg)
    pPosValues <- seq(pPosMin, 0.5, by = (0.5-pPosMin)/N)
    betas <- (N.pos - N.pos*pPosValues)/pPosValues/N.neg  
  }
  
  #remove beta=1
  betas <- betas[-1] 
  
  #cat("\n undersampling betas with method", method, ":\n", betas, "\n")
  
  return(betas)
}





