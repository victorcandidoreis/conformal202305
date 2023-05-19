library(ranger)

library(pracma)

library(magrittr)

conformal.supplies = function(dataset
              ,y.name
              ,est.func
              ,score.func
              ,y.grid
              ,cov.names
              ,...){
  
  y.column = which(dataset  %>%  colnames() == y.name)
  
  base.est = dataset %>% est.func(ycol = y.column
  ,ygrid = y.grid
  ,xnames = cov.names
  ,...)
  
  if(base.est$split.object %>% length() == 1){
    
    eval.scores = sapply(1:base.est$n,
           function(a){
             dens.x=base.est$def.estimation.x(dataset[a,-y.column]
                                       ,index.train=a)
             
             score.supp = list(
               cut.off = dens.x(dataset[a,y.column])
               ,densities = list(dens.x(y.grid)))
             
             return(score.supp %>% score.func(ygrid = y.grid,...))
           }
    )
  } else {
    
    eval.scores = sapply(1:length(base.est$split.object),
           function(a){
              dens.x=base.est$def.estimation.x(dataset[
                base.est$split.object[a],-y.column]
               ,index.train=base.est$split.object[a])
                           
              score.supp = list(
                 cut.off = dens.x(dataset[base.est$split.object[a],y.column])
                ,densities = list(dens.x(y.grid)))
              
              return(score.supp %>% score.func(ygrid = y.grid,...))
          }
    )
    
  }

  
  return(list(base.estimation = base.est$def.estimation.x
      ,evaluated.scores = eval.scores
      ,score.function = score.func
      ,split.object = base.est$split.object
      ,n = base.est$n
      ,covariates.names = cov.names
      ,ygrid = y.grid))
}

conformal.regions = function(covariates
             ,coverage.level = 0.025
             ,score.compare.func
             ,conf.suppls
             ,...){
  
  regions = covariates %>% 
    score.compare.func(alpha = coverage.level
    ,suppls = conf.suppls
    ,...)
  
  return(regions)
  
}

dens.est = function(df
    ,ycol
    ,ygrid
    ,xnames
    ,ntrees=250
    ,split.obj = F
    ,...){

  n = df %>% nrow()
  
  trees.without.obs = NULL
  
  model.form = colnames(df)[ycol] %>% paste0("~.") %>% 
    formula()
  
  if(split.obj){
    
    split.obj = sample(n,n/2)
    
    forest = model.form %>% ranger(data = df[-split.obj,]
          ,num.trees = ntrees
          ,keep.inbag = T)
    
  } else {
    
    forest = model.form %>% ranger(data = df
    ,num.trees = ntrees
    ,keep.inbag = T)
   
    trees.without.obs = lapply(1:n,
      function(obs.i){
        
        return(sapply(forest$inbag.counts,
          function(counts){
            
            return(counts[obs.i]==0)
            
          }
        ) %>% which())
        
      }
    )
     
  }
  
  term.nodes.lab.data = predict(forest,data=df
            ,type="terminalNodes")$predictions
  
  dens.est.x=function(covs=NULL
      ,index.train=NULL
      ,model = forest
      ,d.f = df
      ,y.col = ycol
      ,split.ob = split.obj
      ,ter.nod.l.d = term.nodes.lab.data
      ,tr.without.o = trees.without.obs
      ,y.grid=ygrid
      ,nt = n
      ,x.names=xnames
      ,...){
    
    
    if(covs %>% dim() %>% is.null()){
      
      covs = covs %>% t() %>%  as.data.frame()
      
    }
    
    if(!is.null(x.names)){
      
      colnames(covs) = x.names
      
    }
    
    if(split.ob %>% length() == 1){
      
      if(index.train %>% is.null()){
        
        term.nodes.new.data = predict(forest, data = covs
  ,type = "terminalNodes")$predictions
      
        mixture.weights = lapply(1:nt,
          function(i){
            
            mixture.weights=sapply(tr.without.o[[i]],
              function(tree.j)
              {
 flag.same.node=ter.nod.l.d[-i,
   tree.j]==term.nodes.new.data[tree.j]
       
 return(flag.same.node/sum(flag.same.node))
              }
            )
          }
        )
        
        y = lapply(1:nt,
          function(i){
           
            return(d.f[-i,y.col])

          }
        )
        
      } else {
        
        term.nodes.new.data = ter.nod.l.d[index.train,]
        
        tr.without.i = tr.without.o[[index.train]]
        
        mixture.weights = sapply(tr.without.i,
          function(tree.j)
          {
            flag.same.node=ter.nod.l.d[-index.train,
              tree.j] == term.nodes.new.data[tree.j]
          
            return(flag.same.node/sum(flag.same.node))
          }
        )
        
        mixture.weights = list(mixture.weights)
        
        y=list(d.f[-index.train,y.col])
        
      }
      
    } else {
      
      if(index.train %>% is.null()){
        
        term.nodes.new.data = predict(forest, data = covs
  ,type = "terminalNodes")$predictions
        
      } else {
        
        term.nodes.new.data = ter.nod.l.d[index.train,]
        
      }
      
      mixture.weights = sapply(1:ncol(ter.nod.l.d),
        function(tree.j)
        {
          flag.same.node=ter.nod.l.d[-split.obj,
            tree.j] == term.nodes.new.data[tree.j]
   
          return(flag.same.node/sum(flag.same.node))
        }
      )
      
      mixture.weights = list(mixture.weights)
      
      y=list(d.f[-split.ob,y.col])
      
    }
    
    mixture.weights = lapply(1:length(mixture.weights),
      function(a){
 
        a = mixture.weights[[a]] %>% rowMeans() %>% sqrt()
 
        a = a/sum(a)
 
        return(a)
 
      }
    )
    
    bandwidths = lapply(1:length(mixture.weights),
     function(a){
       
       return(RFCDE:::select_bandwidth(
          as.matrix(y[[a]])
         ,length(mixture.weights[[a]])*mixture.weights[[a]]
         ,"plugin"))
       
     }
    )
    
    dens.given.x = function(responses
          ,mixt.ws=mixture.weights
          ,y.train=y
          ,bands=bandwidths
          ,obs.i=NULL){
      
      if(obs.i %>% is.null()){
        
        return(RFCDE:::kde_estimate(
          as.matrix(y[[1]])
         ,responses
         ,length(mixture.weights[[1]])*mixture.weights[[1]]
         ,bands[[1]]))
        
      } else {
        
        return(RFCDE:::kde_estimate(
          as.matrix(y[[obs.i]])
          ,responses
          ,length(mixture.weights[[obs.i]])*
    mixture.weights[[obs.i]]
          ,bands[[obs.i]]))
        
      }
      
    }
    
    return(dens.given.x)
    
  }
  
  nt=n
  
  return(list(def.estimation.x = dens.est.x
      ,split.object = split.obj
      ,n = nt))
  
}


int.score = function(score.supplies
   ,ygrid,...){
  
  sapply(1:length(score.supplies$cut.off),
   function(a){

    index.cond=score.supplies$densities[[a]]<=
      score.supplies$cut.off[a]
    
    return(1-trapz(ygrid[index.cond],
              score.supplies$densities[[a]][index.cond]))
    
   }
  )
  
}

score.compare = function(covs
         ,suppls
         ,alpha
         ,...){
  
  intervals=lapply(1:nrow(covs),
   function(a){
     
     if(suppls$split.object %>% length() == 1){
      
      dens.given.x = suppls$base.estimation(covs[a,])
     
      densities = lapply(1:suppls$n,
                         function(i){
                           
                           suppls$ygrid %>% dens.given.x(obs.i=i)
                           
                         }
      )
      
      score.matrix = sapply(1:length(suppls$ygrid),
       function(b){
        
         # densities = lapply(1:suppls$n,
         #                    function(i){
         #                      
         #                      suppls$ygrid %>% dens.given.x(obs.i=i)
         #                      
         #                    }
         # ) 
         
        scores = sapply(1:suppls$n,
          function(a){

              score.supp = list(
              cut.off =densities[[a]][b] #suppls$ygrid[b] %>% dens.given.x(obs.i=a)
              ,densities = densities[a])
                               
              return(score.supp %>% suppls$score.function(
                ygrid=suppls$ygrid,...))
         }
        )
        
        return(scores)
          
       }
      )
      
      region.b = score.matrix %>% sweep(
           MARGIN = 1
          ,STATS = suppls$evaluated.scores
          ,FUN = ">"
          ,check.margin = F)
      
      region.b = region.b %>% colSums() < 
        ((1-alpha)*(suppls$n+1))
      
      
       
     } else {
       
       dens.given.x = suppls$base.estimation(covs[a,])
       
       densities = suppls$ygrid %>% dens.given.x()
 
 
      score.vector = sapply(1:length(suppls$ygrid),
            function(a){
              score.supp = list(
                cut.off = densities[a]
                ,densities = list(densities))
                
              return(score.supp %>% suppls$score.function(
                ygrid=suppls$ygrid,...))
            })
        
      region.b = suppls$evaluated.scores %>% 
            quantile(probs = (1-2*alpha)) >= score.vector
         
       
     }
     
     aux=c(F,region.b) %>% rle()
     
     cum.aux=cumsum(aux$lengths)
     
     indext=which(aux$values==T)
     
     cum.auxt=cum.aux[indext]
     
     cum.auxf=cum.aux[-(indext)]
     
     ref=length(cum.auxt)
     
     start.v=cum.auxf[1:ref]
     
     end.v=cum.auxt[1:ref]-1
     
     return(data.frame(start=suppls$ygrid[start.v],
        end=suppls$ygrid[end.v]))

   }
  )
  
  if( intervals %>% length() == 1){
    
    return(intervals[[1]])
    
  } else {
    
    return(intervals)
    
  }
  
}
