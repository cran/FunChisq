fun.chisq.test <- function (x, method="default")
{
 DNAME <- deparse(substitute(x))
 row.chisq.sum <- sum(apply(x, 1,
 		      function(v){
 			if(sum(v)>0) sum((v - mean(v))^2/mean(v))
 			else 0
 			}))
 
 pooled <- apply(x, 2, sum)
 fun.chisq <- row.chisq.sum - sum((pooled - mean(pooled))^2/mean(pooled))
 df <- nrow(x) * (ncol(x) - 1) - (ncol(x) - 1)
 normalized <- as.numeric((fun.chisq-df)/sqrt(2*df))
 
 if(method=="normalized")
 {
   names(normalized) <- "statistic"
   names(df) <- "parameter"
   return(structure(list(statistic = normalized, parameter = df, p.value = pnorm( normalized, lower.tail=F ), data.name= DNAME), class = "htest"))
 }else if(method=="default")
 {
   names(fun.chisq) <- "statistic"
   names(df) <- "parameter"
   p.value <- pchisq(fun.chisq, df = df, lower.tail=FALSE)
   return(structure(list( statistic=fun.chisq, parameter=df, p.value=p.value, data.name= DNAME) ,class = "htest"))
 }
}
