#
#setwd("~/conformal/comparacao_geral/bimodal/")

library(parallel)

library(foreach)

library(doParallel)


cl = makeCluster(3)

registerDoParallel(cl)

h.start = Sys.time()

saveRDS(h.start,"inicio.RDS")

n.sim = 5e2


# n.cond=1e2
# 
# #Fixar covariaveis para averiguar a cobertura marginal
# 
# source("funcoes/generate_bimodal.R")
# 
# set.seed(10)
# 
# covs.cond=generate_bimodal(n.cond,20)
# 
# covs.cond=as.data.frame(covs.cond[,-1])
# 
# colnames(covs.cond)=paste0("x",1:ncol(covs.cond))
# 
# saveRDS(covs.cond,"dados_condicional.RDS")


process=foreach(i = 1:n.sim) %dopar% {
  
  # source("conformal_f.R")
  # #
  # source("generate_bimodal.R")
  
  #
  arquivos=list.files("funcoes")
  
  for(jj in 1:length(arquivos)){
    source(paste0("funcoes/",arquivos[jj]))
  }
  
  
  n.df = 2500
  #
  df=generate_bimodal(n.df,20)
  #
  new.yx=generate_bimodal(1,20)
  #
  ygrid = seq(-6,10,length.out = 1e2)
  
  dados.condicional=readRDS("dados_condicional.RDS")
  
  #----------------------------------------------------------------
  

  base.fits=list(hpd_flex=fit_density_forest,
                 dist=fit_density_forest,
                 quantile=fit_quantile_forest,
                 reg_n_weighted=fit_regression_forest,
                 reg_weighted=fit_regression_mean_error_forest)
  
  method.fits=list(hpd_flex=hpd_split_prediction_bands,
                   dist=dist_split_prediction_bands,
                   quantile=quantile_split_prediction_bands,
                   reg_n_weighted=reg_split_prediction_bands,
                   reg_weighted=reg_weighted_split_prediction_bands)
  
  
  conformal.methods=c("hpd_leave",names(base.fits),"hpd_forest")
  
  for(ii in 1:length(conformal.methods)){
    h.p.start=Sys.time()
    if(conformal.methods[ii]=="hpd_forest"){
      
      conf.sup=conformal.supplies(dataset = df
                                  ,y.name = "y"
                                  ,est.func = dens.est
                                  ,score.func = int.score
                                  ,y.grid = ygrid
                                  ,cov.names = paste0("x",1:(ncol(df)-1))
                                  ,split.obj=T)
      
      
      saveRDS(Sys.time()-h.p.start,
              paste0("tempo_treino_",conformal.methods[ii],"_",i,".RDS"))
      h.p.start=Sys.time()
      
      
      region=conformal.regions(covariates = as.data.frame(t(new.yx[1,-1])) #new.yx[1,-1]
                               ,conf.suppls = conf.sup
                               ,score.compare.func = score.compare
                               ,coverage.level = 0.025)
      
      saveRDS(Sys.time()-h.p.start,
              paste0("tempo_predicao_",conformal.methods[ii],"_",i,".RDS"))
      
      saveRDS(list(new.yx=new.yx,region=region),
              paste0("infos_marginal_",conformal.methods[ii],"_",i,".RDS"))
      
      
      infos.condicional=list()
      
      gc()
      
      for(j in 1:nrow(dados.condicional)){
        
        new.x=dados.condicional[j,]
        
        colnames(new.x)=paste0("x",1:ncol(dados.condicional))
        
        infos.condicional[[j]]=list(new.x=new.x,region=conformal.regions(covariates = new.x
                                                                         ,conf.suppls = conf.sup
                                                                         ,score.compare.func = score.compare
                                                                        ,coverage.level = 0.025))
      }
      
      saveRDS(infos.condicional,
              paste0("infos_condicional_",conformal.methods[ii],"_",i,".RDS"))
      
      
    } 
    else if (conformal.methods[ii]=="hpd_leave"){
      
      conf.sup=conformal.supplies(dataset = df
                                  ,y.name = "y"
                                  ,est.func = dens.est
                                  ,score.func = int.score
                                  ,y.grid = ygrid
                                  ,cov.names = paste0("x",1:(ncol(df)-1))
                                  ,split.obj=F)
      
      
      saveRDS(Sys.time()-h.p.start,
              paste0("tempo_treino_",conformal.methods[ii],"_",i,".RDS"))
      h.p.start=Sys.time()
      
      
      region=conformal.regions(covariates = as.data.frame(t(new.yx[1,-1])) #new.yx[1,-1] 
                               ,conf.suppls = conf.sup
                               ,score.compare.func = score.compare
                               ,coverage.level = 0.05)
      saveRDS(Sys.time()-h.p.start,
              paste0("tempo_predicao_",conformal.methods[ii],"_",i,".RDS"))
      
      saveRDS(list(new.yx=new.yx,region=region),
              paste0("infos_marginal_",conformal.methods[ii],"_",i,".RDS"))
      
      
      infos.condicional=list()
      
      for(j in 1:nrow(dados.condicional)){
        
        new.x=dados.condicional[j,]
        
        colnames(new.x)=paste0("x",1:ncol(dados.condicional))
        
        infos.condicional[[j]]=list(new.x=new.x,region=conformal.regions(covariates = new.x
                                                                         ,conf.suppls = conf.sup
                                                                         ,score.compare.func = score.compare
                                                                         ,coverage.level = 0.05))
      }
      
      saveRDS(infos.condicional,
              paste0("infos_condicional_",conformal.methods[ii],"_",i,".RDS"))
      
      
    }
    else{
      
      which_train <- sample(1:(nrow(df)/2),nrow(df)*0.7/2)
      
      which_valid = (1:(nrow(df)/2))[!((1:(nrow(df)/2)) %in% which_train)]
      
      base.fit=fit.orchestrator(type="base",
                                func=base.fits[[conformal.methods[ii]]],
                                fit=NULL,
                                xTrain=df[which_train,-1,drop=FALSE],
                                yTrain=df[which_train,1,drop=FALSE],
                                xValidationTest=df[which_valid,-1,drop=FALSE],
                                yValidationTest=df[which_valid,1,drop=FALSE],
                                alpha=NULL,
                                y_grid=NULL)
      
      y.grid.method=NULL
      
      if(conformal.methods[ii]=="reg_n_weighted" || 
         conformal.methods[ii]=="reg_weighted"){
        y.grid.method=ygrid
      }
      
      saveRDS(Sys.time()-h.p.start,
              paste0("tempo_treino_",conformal.methods[ii],"_",i,".RDS"))
      h.p.start=Sys.time()
      
      
      
      method.fit=fit.orchestrator(type="method",
                                  func=method.fits[[conformal.methods[ii]]],
                                  fit=base.fit,
                                  xTrain=df[(nrow(df)/2+1):nrow(df),-1,drop=FALSE],
                                  yTrain=df[(nrow(df)/2+1):nrow(df),1,drop=FALSE],
                                  xValidationTest=new.yx[,-1,drop=FALSE],
                                  yValidationTest=NULL,
                                  alpha=0.05,
                                  y_grid=y.grid.method)
      
      saveRDS(Sys.time()-h.p.start,
              paste0("tempo_predicao_",conformal.methods[ii],"_",i,".RDS"))
      
      region=region.orchestrator(conformal.methods[ii],
                                 method.fit)
      
      
      saveRDS(list(new.yx=new.yx,region=region),
              paste0("infos_marginal_",conformal.methods[ii],"_",i,".RDS"))
      
      
      infos.condicional=list()
      
      for(j in 1:nrow(dados.condicional)){
        
        new.x=dados.condicional[j,]
        
        colnames(new.x)=paste0("x",1:ncol(dados.condicional))
        
        method.fit=fit.orchestrator(type="method",
                                    func=method.fits[[conformal.methods[ii]]],
                                    fit=base.fit,
                                    xTrain=df[(nrow(df)/2+1):nrow(df),-1,drop=FALSE],
                                    yTrain=df[(nrow(df)/2+1):nrow(df),1,drop=FALSE],
                                    xValidationTest=new.x,
                                    yValidationTest=NULL,
                                    alpha=0.05,
                                    y_grid=y.grid.method)
        
        region=region.orchestrator(conformal.methods[ii],
                                   method.fit)
        
        infos.condicional[[j]]=list(new.x=new.x,region=region)
      }
      
      saveRDS(infos.condicional,
              paste0("infos_condicional_",conformal.methods[ii],"_",i,".RDS"))
    }
  }
  
  
  rm(list=ls())
  
  gc()
  
}

h.end = Sys.time()

stopCluster(cl)

saveRDS(h.end-h.start,"tempo.RDS")