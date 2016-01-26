# FunChisq.R -- statistical tests for nonparametric functional dependencies
#
# YZ, HZ, MS
# Modified: Feb 8, 2015; June 25, 2015
# Jan 22, 2016 MS:
#    1. introduced the type argument to indicate either functional or non-constant
#       functional chi-square
#    2. introduced the log.p argument to obtain the log of p-value
# Jan 23, 2016 MS:
#    Added a return of an estimate of functional index between 0 and 1. It is asymetrical and
#      different from Cramer's V.

fun.chisq.test <- function (x, method="default", type="non-constant", log.p=FALSE)
{
  if(!is.matrix(x) && !is.data.frame(x)) stop("input x must be matrix or data frame\n")

  row.chisq.sum <- sum(apply(x, 1,
                             function(v){
                               row.sum <- sum(v)
                               expected <- row.sum / length(v)
                               if(row.sum>0) sum( (v - expected)^2 / expected)
                               else 0
                             }))

  col.sum <- apply(x, 2, sum)

  if(type == "non-constant") { # non-constant functional chi-square
    fun.chisq <- row.chisq.sum - sum((col.sum - mean(col.sum))^2/mean(col.sum))
    df <- nrow(x) * (ncol(x) - 1) - (ncol(x) - 1)
  } else if(type == "all") { # functional chi-square
    fun.chisq <- row.chisq.sum
    df <- nrow(x) * (ncol(x) - 1)
  } else {
    stop("ERROR: unknown type functional chi-square!\n")
  }

  if(ncol(x) > 1) {
    estimate <- sqrt(abs(fun.chisq) / sum(col.sum) / (ncol(x) - 1))
  } else {
    estimate <- 0
  }
  names(estimate) <- "estimate"

  DNAME <- deparse(substitute(x))

  if(method=="normalized") {
    method.text <- "Normalized functional chi-square test"
    normalized <- as.numeric((fun.chisq-df)/sqrt(2*df))
    names(normalized) <- "statistic"
    names(df) <- "parameter"
    p.value <- pnorm( normalized, lower.tail=FALSE, log.p=log.p )
    return(structure(list(statistic = normalized, parameter = df, p.value = p.value,
                          estimate = estimate, data.name = DNAME, method = method.text),
                     class = "htest"))

  } else if(method=="default") {
    method.text <- "Functional chi-square test"
    names(fun.chisq) <- "statistic"
    names(df) <- "parameter"
    p.value <- pchisq( fun.chisq, df = df, lower.tail=FALSE, log.p=log.p )
    return(structure(list( statistic=fun.chisq, parameter=df, p.value=p.value,
                           estimate = estimate, data.name=DNAME,
                           method = method.text),
                     class = "htest"))

  } else if(method=="exact") {

    method.text <- "Exact functional test"
    if(sum(x%%1!=0)>=1) { # Check whether numbers in x are all integers
      stop("ERROR: Exact test requires integers in the contingency table!", call. = TRUE)
    }
    ####

    ####
    #Hua added, Nov 13, 2014
    #Exact functional test
    if((sum(x) <= 200 || sum(x)/nrow(x)/ncol(x) <=5)
       && nrow(x)<=5 && ncol(x)<=5) {
      p.value <- exact.functional.test(x)
      if(log.p) p.value <- log(p.value)
      names(fun.chisq) <- "statistic"
      return(structure(list(statistic = fun.chisq, p.value = p.value, estimate = estimate,
                            data.name = DNAME, method = method.text),
                       class = "htest"))
    } else {
      return(fun.chisq.test(x, method="default", type=type, log.p=log.p))
    }

    ####
  } else {
    stop("method argument invalid!")
  }
}

cp.fun.chisq.test <- function(x, method="default")
{
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

  DNAME <- deparse(substitute(x))

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
    stop("method can only be \"default\", \"normalized\", or \"exact\".\n")
  }
}

exact.functional.test <- function(x){
  res <- .Call("ExactFunctionalTest", x, PACKAGE="FunChisq")
  return (as.double(res))
}
####
