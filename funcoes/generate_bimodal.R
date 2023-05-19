generate_bimodal <- function(n,d,x=NULL)
{
  if(is.null(x))
  {
    x=matrix(runif(n*d,-1.5,1.5),n,d)
  }
  f=(x[,1]-1)^2*(x[,1]+1)
  g=rep(0,n)
  g[x[,1]> -0.5]=2*sqrt(x[x[,1]> -0.5,1]+0.5)
  s=1/4+abs(x[,1])
  # response
  y=ifelse(runif(n)>0.5,f-g,f+g)+rnorm(n,0,sqrt(s))
  
  df=cbind(y,x)
  
  colnames(df)=c("y",paste0("x",1:ncol(x)))
  
  return(df)
}
