
# returns prediction bands for each element of xTest
# xTrain and yTrain are is I2, hold out data not used
#                          to fit the density cde_fit
quantile_split_prediction_bands <- function(quantile_fit,
                                            xTrain,yTrain,
                                            xTest,
                                            alpha=0.1,
                                            yTest=NULL)
{
  pred_test <- predict(quantile_fit,xTest,what=c(alpha/2,1-alpha/2))
  pred_train <- predict(quantile_fit,xTrain,what=c(alpha/2,1-alpha/2))
  E <- apply(cbind(pred_train[,1]-yTrain,yTrain-pred_train[,2]),1,max)
  Q <- quantile(E,probs = 1-alpha)
  
  limits <- cbind(pred_test[,1]-Q,pred_test[,2]+Q) 
  return(list(limits=limits))
  
}


# Returns fit_dist_split with indicator of if each yTest  belongs to prediction band
#        (used for diagnostics)
quantile_split_prediction_bands_evalY <- function(fit_quantile_split,
                                                  yTest)
{
  fit_quantile_split$yTest_covered <- fit_quantile_split$limits[,1]<=yTest & fit_quantile_split$limits[,2]>=yTest
  fit_quantile_split$prediction_bands_size <- fit_quantile_split$limits[,2]-fit_quantile_split$limits[,1]
  fit_quantile_split$limits <- NULL
  return(fit_quantile_split)
}
