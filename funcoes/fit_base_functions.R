fit_density_forest <- function(xTrain,yTrain,xValidation,yValidation,nCores=1)
{
  fit=fitFlexCoDE(xTrain=xTrain,zTrain=yTrain,
                  xValidation=xValidation,zValidation=yValidation,nIMax = 35,
                  regressionFunction = regressionFunction.Forest,
                  regressionFunction.extra = list(nCores=nCores))
  return(fit)
}

fit_quantile_forest <- function(xTrain,yTrain,xValidation,yValidation)
{
  
  fit <- quantregForest(x=rbind(xTrain,xValidation),
                        y=c(yTrain[,1],yValidation[,1]),nthreads =1)
  return(fit)
}


fit_regression_forest <- function(xTrain,yTrain,xValidation,yValidation)
{
  fit <- randomForest(x=rbind(xTrain,xValidation),
                      y=c(yTrain[,1],yValidation[,1]))
  #class(fit) <- "forest"
  return(fit)
}

fit_regression_mean_error_forest <- function(xTrain,yTrain,xValidation,yValidation)
{
  fit_mean <- randomForest(x=xTrain,
                           y=yTrain[,1])
  pred_val<- predict(fit_mean, xValidation)
  fit_error <- randomForest(x=xValidation,
                            y=abs(yValidation[,1]-pred_val))
  fit <- list(fit_mean=fit_mean,
              fit_error=fit_error)
  class(fit) <- "forest_weighted"
  return(fit)
}


















    

