fit.orchestrator=function(type,
                          func=NULL,
                          fit=NULL,
                          fit2=NULL,
                          xTrain,
                          yTrain,
                          xValidationTest,
                          yValidationTest=NULL,
                          alpha=NULL,
                          y_grid=NULL){
  
  if(type=="base"){
    return(func(xTrain,
                yTrain,
                xValidationTest,
                yValidationTest))
  }
  
  if(type=="method"){
    if(!is.null(y_grid)){
      if(!is.null(fit$fit_error)){
        return(func(fit$fit_mean,
                    fit$fit_error,
                    xTrain,
                    yTrain,
                    xValidationTest,
                    alpha,
                    y_grid))
      } else {
        return(func(fit,
                    xTrain,
                    yTrain,
                    xValidationTest,
                    alpha,
                    y_grid))
      }
    } else {
      return(func(fit,
                  xTrain,
                  yTrain,
                  xValidationTest,
                  alpha))
    }
  }
}

region.orchestrator=function(method.name,
                             method.fit){
  if(method.name=="quantile"){
    return(matrix(method.fit$limits,1,2))
  }
  
  if(method.name=="dist" || method.name=="reg_n_weighted" || method.name=="reg_weighted"){
    return(method.fit$y_grid[method.fit$prediction_bands[[1]]]
           %>% (function(a) matrix(c(min(a),max(a)),1,2)))
  }
  
  if(method.name=="hpd_flex"){
    
    aux=c(F,method.fit$prediction_bands[[1]]) %>% rle()
    
    cum.aux=cumsum(aux$lengths)
    
    indext=which(aux$values==T)
    
    cum.auxt=cum.aux[indext]
    
    cum.auxf=cum.aux[-(indext)]
    
    ref=length(cum.auxt)
    
    start.v=cum.auxf[1:ref]
    
    end.v=cum.auxt[1:ref]-1
    
    return(data.frame(start=method.fit$y_grid[start.v],
                      end=method.fit$y_grid[end.v]))
    
  }
  
}

