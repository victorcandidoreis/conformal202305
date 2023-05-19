# returns prediction bands for each element of xTest
# xTrain and yTrain are is I2, hold out data not used
#                          to fit the density cde_fit
# t_grid is a grid of values for the f(.|x)
reg_weighted_split_prediction_bands <- function(reg_fit_mean,
                                                reg_fit_error,
                                                xTrain,yTrain,
                                                xTest,
                                                alpha=0.1,
                                                y_grid)
{
  yTrain=yTrain[,1]
  
  pred_train_mean <- predict(reg_fit_mean,xTrain)
  pred_test_mean <- predict(reg_fit_mean,xTest)
  pred_train_error <- predict(reg_fit_error,xTrain)
  pred_test_error <- predict(reg_fit_error,xTest)
  
  # observed densities:
  conformity_score_train <- -abs(pred_train_mean-yTrain)/pred_train_error
  
  #prediction_bands <- list()
  prediction_bands_which_belong <- list()
  
  ths <- quantile(conformity_score_train,probs=alpha)
  for(ii in 1:nrow(xTest))
  {
    prediction_bands_which_belong[[ii]] <- -abs(y_grid-pred_test_mean[ii])/pred_test_error[ii] >=ths
  }
  return(list(prediction_bands=prediction_bands_which_belong,
              y_grid=y_grid,pred_test_mean=pred_test_mean,
              pred_test_error=pred_test_error,
              ths=ths))
  
  
}

reg_weighted_split_prediction_bands_evalY <- function(fit_reg_weighted_split,
                                                      yTest)
{
  conformity_score_test <- -abs(fit_reg_weighted_split$pred_test_mean-yTest)/fit_reg_weighted_split$pred_test_error
  yTest_covered <- conformity_score_test >= fit_reg_weighted_split$ths
  fit_reg_weighted_split$yTest_covered <- yTest_covered
  fit_reg_weighted_split$prediction_bands_size <- sapply(fit_reg_weighted_split$prediction_bands,
                                                         function(x)
                                                           mean(x))*diff(range(fit_reg_weighted_split$y_grid))
  fit_reg_weighted_split$prediction_bands <- NULL
  fit_reg_weighted_split$pred_test_mean <- NULL
  fit_reg_weighted_split$pred_test_error <- NULL
  
  return(fit_reg_weighted_split)
}

