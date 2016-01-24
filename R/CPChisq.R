# CPChisq.R -- comparative chi-square test for heterogeneous associations
#
# MS
# Created: July 3, 2015 generalized from two contingency tables in 
#                       "hetero-chisq-steel.R"

cp.chisq.test <- function(x, method="default")
{
  DNAME <- deparse(substitute(x))
  if(mode(x)!="list" || length(x)<2 )
  {
    stop("only accept list of 2 or more matrices as input!")
  }
  
  pooled <- x[[1]]
  for(i in 2:length(x)) {
    pooled <- pooled + x[[i]]
  }

  # Row sum of matrix pooled
  pooled.row.sum <- apply(pooled, 1, sum)
  
  # Column sum of matrix pooled
  pooled.col.sum <- apply(pooled, 2, sum)
  
  #Sum of matrix pooled
  pooled.sum <- sum(pooled.col.sum)

  pooled.chisq <- 0
  
  if(pooled.sum > 0) {
    #Expected value in pooled
    pooled.expected <- pooled.row.sum %*% t(pooled.col.sum) / pooled.sum
    
    #Chi-square value of pooled
    pooled.chisq <- sum((pooled - pooled.expected)^2 / pooled.expected, na.rm = TRUE)
    
    total.chisq <- 0 
    
    for(i in seq_along(length(x))) {
      #Expected values in matrix i
      expected <- pooled.expected / pooled.sum * sum(x[[i]])
      
      #Chi-square value of matrix i
      chisq <- sum((x[[i]] - expected)^2 / expected, na.rm = TRUE)
      
      total.chisq <- total.chisq + chisq
    }

    hetero.chisq  <- total.chisq - pooled.chisq
    df <- (length(x) - 1) * (sum(pooled.row.sum != 0) - 1) * 
      (sum(pooled.col.sum != 0) - 1)
    
  } else {
    
    hetero.chisq <- 0
    df <- 0
    
  }

  if(method=="normalized") {

    finalStat <- ifelse(df > 0, ( hetero.chisq - df) / sqrt( 2 * df ), 0)
    names(finalStat) <- "statistic"
    
    finalDf <- df
    names(finalDf) <- "parameter"
    
    p.val <- ifelse(df > 0, pnorm( finalStat, lower.tail=FALSE ), 1)
    names(p.val) <- "p.value"
    
    return(structure(list(statistic = finalStat, parameter = finalDf, p.value = p.val, 
                          method = "Nomalized comparative chi-square test for heterogeneity",
                          data.name= DNAME),
                     class = "htest"))
    
  } else if(method=="default") {
    
    finalStat <- hetero.chisq
    names(finalStat) <- "statistic"
    
    finalDf <- df
    names(finalDf) <- "parameter"
    
    p.val <- pchisq(finalStat, df = finalDf, lower.tail=FALSE)
    names(p.val) <- "p.value"
    
    return(structure(list( statistic=finalStat, parameter=finalDf, p.value=p.val, 
                           method = "Comparative chi-square test for heterogeneity",
                           data.name= DNAME), 
                     class = "htest"))
    
  } else {
    stop("method can only be \"default\", \"normalized\", or \"exact\".\n")
  }  
}
