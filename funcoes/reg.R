# returns prediction bands for each element of xTest
# xTrain and yTrain are is I2, hold out data not used
#                          to fit the density cde_fit
# t_grid is a grid of values for the f(.|x)
reg_split_prediction_bands <- function(reg_fit,
                                       xTrain,yTrain,
                                       xTest,
                                       alpha=0.1,
                                       y_grid)
{
  yTrain=yTrain[,1]
  
  pred_test <- predict(reg_fit,xTest)
  pred_train <- predict(reg_fit,xTrain)
  # observed densities:
  conformity_score_train <- -(pred_train-yTrain)^2
  
  #prediction_bands <- list()
  prediction_bands_which_belong <- list()
  
  ths <- quantile(conformity_score_train,probs=alpha)
  for(ii in 1:nrow(xTest))
  {
    prediction_bands_which_belong[[ii]] <- -(y_grid-pred_test[ii])^2>=ths
  }
  
  return(list(prediction_bands=prediction_bands_which_belong,
              y_grid=y_grid,pred_test=pred_test,ths=ths))
  
}


# Returns fit_reg_split with indicator of if each yTest  belongs to prediction band
#        (used for diagnostics)
reg_split_prediction_bands_evalY <- function(fit_reg_split,
                                             yTest)
{
  conformity_score_test <- -(fit_reg_split$pred_test-yTest)^2
  yTest_covered <- conformity_score_test >= fit_reg_split$ths
  fit_reg_split$yTest_covered <- yTest_covered
  fit_reg_split$prediction_bands_size <- sapply(fit_reg_split$prediction_bands,
                                                function(x)
                                                  mean(x))*diff(range(fit_reg_split$y_grid))
  fit_reg_split$prediction_bands <- NULL
  fit_reg_split$pred_test <- NULL
  return(fit_reg_split)
}