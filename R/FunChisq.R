# FunChisq.R
#
# YZ, HZ, MS
# Modified: Feb 8, 2015

fun.chisq.test <- function (x, method="default")
{
  #Check numbers in x are integers
  if(sum(x%%1!=0)>=1)stop("Numbers in the contingency table should be integers!", call. = TRUE)
  ####
  
  DNAME <- deparse(substitute(x))
  
  row.chisq.sum <- sum(apply(x, 1,
                             function(v){
                               if(sum(v)>0) sum((v - mean(v))^2/mean(v))
                               else 0
                             }))
  
  pooled <- apply(x, 2, sum)
  fun.chisq <- row.chisq.sum - sum((pooled - mean(pooled))^2/mean(pooled))
  df <- nrow(x) * (ncol(x) - 1) - (ncol(x) - 1)
    
  if(method=="normalized") {

    normalized <- as.numeric((fun.chisq-df)/sqrt(2*df))
    names(normalized) <- "statistic"
    names(df) <- "parameter"
    return(structure(list(statistic = normalized, parameter = df, 
                          p.value = pnorm( normalized, lower.tail=F ), 
                          data.name = DNAME,
                          method = "Normalized functional chi-square test for functional dependencies"), 
                     class = "htest"))
    
  } else if(method=="default") {
    
    names(fun.chisq) <- "statistic"
    names(df) <- "parameter"
    p.value <- pchisq(fun.chisq, df = df, lower.tail=FALSE)
    return(structure(list( statistic=fun.chisq, parameter=df, p.value=p.value, data.name= DNAME,
                           method = "Functional chi-square test for functional dependencies"),
                     class = "htest"))
    
  } else if(method=="exact") {
    ####
    #Hua added, Nov 13, 2014
    #Exact functional test
    if((sum(x) <= 200 || sum(x)/nrow(x)/ncol(x) <=5) && nrow(x)<=5 && ncol(x)<=5){
      DNAME <- deparse(substitute(x))
      p.value <- exact.functional.test(x)
      names(fun.chisq) <- "statistic"
      return(structure(list(statistic = fun.chisq, p.value = p.value, data.name= DNAME, 
                            method="Exact functional test for functional dependencies"), 
                       class = "htest"))
    }else{
      return(fun.chisq.test(x, "default"))
    }
    
    ####
  } else {
    stop("method argument invalid!")
  }
}

cp.fun.chisq.test <- function(x, method="default")
{
  DNAME <- deparse(substitute(x))
  if(mode(x)!="list" || length(x)<2 )
  {
    stop("only accept list of 2 or more matrices as input!")
  }
  
  finalStat <- 0
  finalDf <- 0
  
  for(i in 1:nrow(x[[1]]))
  {
    oneT <- c() # one table for each row
    for(j in 1:length(x))
    {
      oneT <- rbind(oneT, x[[j]][i,])
    }
    oneresult <- fun.chisq.test(oneT)
    finalStat <- finalStat + oneresult$statistic
    finalDf <- finalDf + oneresult$parameter
  }
  
  if(method=="normalized") {
    
    finalStat <- as.numeric((finalStat-finalDf)/sqrt(2*finalDf))
    names(finalStat) <- "statistic"
    names(finalDf) <- "parameter"
    return(structure(list(statistic = finalStat, parameter = finalDf, 
                          p.value = pnorm( finalStat, lower.tail=F ), 
                          method = "Nomalized comparative functional chi-square test for heterogeneity",
                          data.name= DNAME), 
                     class = "htest"))
    
  } else if(method=="default") {
    
    names(finalStat) <- "statistic"
    names(finalDf) <- "parameter"
    p.value <- pchisq(finalStat, df = finalDf, lower.tail=FALSE)
    return(structure(list( statistic=finalStat, parameter=finalDf, p.value=p.value, 
                           method = "Comparative functional chi-square test for heterogeneity",
                           data.name= DNAME), 
                     class = "htest"))
    
  } else {
    stop("method can be either \"default\" or \"normalized\" or \"exact\"")
  }  
}

exact.functional.test <- function(x){
  res <- .Call("ExactFunctionalTest", x, PACKAGE="FunChisq")
  return (as.double(res))
}
####