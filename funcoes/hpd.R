
which_neighbors <-  function(xTrain,xTest,k){
  return(FNN::get.knnx(data=xTrain,query = xTest,k=k)$nn.index)
}


findThresholdHPD=function(binSize,estimates,confidence)
{
  estimates=as.vector(estimates)
  maxDensity=max(estimates)
  minDensity=min(estimates)
  newCut=(maxDensity+minDensity)/2
  eps=1
  ii=1
  while(ii<=1000)
  {
    prob=sum(binSize*estimates*(estimates>newCut))
    eps=abs(confidence-prob)
    if(eps<0.0000001) break; # level found
    if(confidence>prob) maxDensity=newCut
    if(confidence<prob) minDensity=newCut
    newCut=(maxDensity+minDensity)/2
    ii=ii+1
  }
  return(newCut)
}



#                          to fit the density cde_fit
# t_grid is a grid of values for the f(.|x)
hpd_split_prediction_bands <- function(cde_fit,
                                       xTrain,yTrain,
                                       xTest,
                                       alpha=0.1)
{
  pred_test <- predict(cde_fit,xTest)
  pred_train <- predict(cde_fit,xTrain)
  # observed densities:
  which_select <- cbind(1:length(yTrain),
                        which_neighbors(as.matrix(pred_train$z),
                                        as.matrix(yTrain),k=1))
  
  
  which_smaller <- apply(pred_train$CDE<=pred_train$CDE[which_select],1,which)
  conformity_score_train <- rep(NA,nrow(xTrain))
  for(ii in 1:nrow(xTrain))
  {
    conformity_score_train[ii] <- sum(pred_train$CDE[ii,which_smaller[[ii]]])
  }
  band <- diff(pred_train$z)[1]
  conformity_score_train <- conformity_score_train*band
  
  #plot(pred_train$z,pred_train$CDE[2,])
  #abline(v=yTrain[2])
  #points(y=rep(0.05,length(which_smaller[[2]])),x=pred_train$z[which_smaller[[2]]])
  
  
  th <- quantile(conformity_score_train,probs=alpha)
  prediction_bands_which_belong <- list()
  for(ii in 1:nrow(xTest))
  {
    th_hpd <- findThresholdHPD(band,pred_test$CDE[ii,],1-th)
    prediction_bands_which_belong[[ii]] <- pred_test$CDE[ii,]>=th_hpd
  }
  return(list(prediction_bands=prediction_bands_which_belong,
              y_grid=pred_test$z,th=th,pred_test=pred_test))
}



hpd_split_prediction_bands_evalY <- function(fit_hpd_split,yTest)
{
  which_select <- cbind(1:length(yTest),
                        which_neighbors(as.matrix(fit_hpd_split$pred_test$z),
                                        as.matrix(yTest),k=1))
  which_smaller <- apply(fit_hpd_split$pred_test$CDE<=fit_hpd_split$pred_test$CDE[which_select],1,which)
  conformity_score_test <- rep(NA,length(yTest))
  for(ii in 1:length(yTest))
  {
    conformity_score_test[ii] <- sum(fit_hpd_split$pred_test$CDE[ii,which_smaller[[ii]]])
  }
  band <- diff(fit_hpd_split$pred_test$z)[1]
  conformity_score_test <- conformity_score_test*band
  
  
  yTest_covered <- conformity_score_test >= fit_hpd_split$th
  fit_hpd_split$yTest_covered <- yTest_covered
  fit_hpd_split$pred_test <- NULL
  fit_hpd_split$prediction_bands_size <- sapply(fit_hpd_split$prediction_bands,function(x)mean(x))*diff(range(fit_hpd_split$y_grid))
  fit_hpd_split$prediction_bands <- NULL
  return(fit_hpd_split)
}


