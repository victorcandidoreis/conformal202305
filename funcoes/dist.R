cum_dist <- function(y_grid,cde_estimates,y_values)
{
  which_closest <- FNN::get.knnx(data = y_grid,
                                 query=y_values,k=1)$nn.index
  apply(as.matrix(1:nrow(cde_estimates)),1,function(xx){
    return(sum(cde_estimates[xx,1:which_closest[xx]])*diff(y_grid)[1])
  })
}

# returns prediction bands for each element of xTest
# xTrain and yTrain are is I2, hold out data not used
#                          to fit the density cde_fit
dist_split_prediction_bands <- function(cde_fit,
                                        xTrain,yTrain,
                                        xTest,
                                        alpha=0.1,
                                        yTest=NULL, 
                                        median=FALSE)
{
  pred_test <- predict(cde_fit,xTest)
  pred_train <- predict(cde_fit,xTrain)
  
  cum_dist_evaluated_train <- cum_dist(pred_train$z,pred_train$CDE,yTrain)
  if(median==FALSE)
  {
    ths <-  quantile(cum_dist_evaluated_train,
                     probs = c(alpha/2,1-alpha/2))
  } else {
    ths <- quantile(abs(cum_dist_evaluated_train-0.5),probs = 1-alpha)
  }
  prediction_bands_which_belong <- list()
  FTest <- matrix(NA,nrow(xTest),length(pred_train$z))
  for (ii in 1:nrow(xTest)){
    FTest[ii,] <- cumsum(pred_test$CDE[ii,])*diff(pred_train$z)[1]
    if(median==FALSE)
    {
      prediction_bands_which_belong[[ii]] <- FTest[ii,]>=ths[1]&FTest[ii,]<=ths[2]  
    } else {
      prediction_bands_which_belong[[ii]] <- abs(FTest[ii,]-0.5)<=ths
    }
    
    
  }
  
  
  return(list(prediction_bands=prediction_bands_which_belong,
              y_grid=pred_train$z,FTest=FTest,ths=ths,median=median))
  
}



# Returns fit_dist_split with indicator of if each yTest  belongs to prediction band
#        (used for diagnostics)
dist_split_prediction_bands_evalY <- function(fit_dist_split,
                                              yTest)
{
  yTest_covered <- rep(NA,length(yTest))
  ths <- fit_dist_split$ths
  for (ii in 1:length(yTest)){
    which_closest <- which.min(abs(fit_dist_split$y_grid-yTest[ii]))
    if(fit_dist_split$median==FALSE)
    {
      yTest_covered[ii] <- fit_dist_split$FTest[ii,which_closest]>=ths[1]&fit_dist_split$FTest[ii,which_closest]<=ths[2]  
    } else {
      yTest_covered[ii] <- abs(fit_dist_split$FTest[ii,which_closest]-0.5)<=ths  
    }
  }
  fit_dist_split$yTest_covered <- yTest_covered
  fit_dist_split$prediction_bands_size <- sapply(fit_dist_split$prediction_bands,
                                                 function(x)
                                                   mean(x))*diff(range(fit_dist_split$y_grid))
  fit_dist_split$prediction_bands <- NULL
  fit_dist_split$FTest <- NULL
  return(fit_dist_split)
}