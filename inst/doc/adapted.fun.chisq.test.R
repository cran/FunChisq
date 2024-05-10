## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----fig.show="hold", out.width="45%", eval=FALSE, echo=FALSE-----------------
#  require(FunChisq)
#  require(infotheo)
#  require(DescTools)
#  require(dqrng)
#  
#  # set seeds
#  set.seed(123)
#  dqset.seed(123)
#  
#  # simulate a noisy `many.to.one' function (Y = f(X) and X != f(Y))
#  tab = simulate_tables(n=1000, nrow=3, ncol=4, type="many.to.one", noise=0.3)$noise.list[[1]]
#  
#  # calculate correct and incorrect direction scores for AdpFunChisq
#  # The P-value of AdpFunChisq is lower the better.
#  correct_afc = fun.chisq.test(tab, method="adapted")
#  incorrect_afc = fun.chisq.test(t(tab), method="adapted")
#  
#  # calculate correct and incorrect direction scores for Conditional Entropy
#  # The score of Conditional Entropy is lower the better.
#  vars = Untable(tab)
#  correct_ce = condentropy( as.numeric(vars[,2]), as.numeric(vars[,1]) )
#  incorrect_ce = condentropy( as.numeric(vars[,1]), as.numeric(vars[,2]) )
#  
#  # plot table with AdpFunChisq and Conditional Entropy Scores
#  # correct direction
#  plot_table(tab, xlab = 'Y', ylab = 'X', value.cex = 0.9, cex.main=0.9, col="seagreen",
#             main=paste0("AdpFunChisq P = ",format(correct_afc$p.value, digits = 2),
#                         "\nConditional Entropy = ",format(correct_ce, digits = 2)))
#  # incorrect direction
#  plot_table(t(tab), xlab = 'X', ylab = 'Y', value.cex = 0.9, cex.main=0.9, col="firebrick1",
#             main=paste0("AdpFunChisq P = ",format(incorrect_afc$p.value, digits = 2),
#                         "\nConditional Entropy = ",format(incorrect_ce, digits = 2)))

## ----out.width="45%", message=FALSE-------------------------------------------
require(FunChisq)
require(infotheo)
require(DescTools)

compare.methods <- function(func) 
{
  plot_table(func, xlab = 'Y', ylab = 'X', value.cex = 0.9, col="seagreen", 
             main = expression("Function:"~italic(Y%~~%f(X))))
  
  plot_table(t(func), xlab = 'X', ylab = 'Y', value.cex = 0.9, col="firebrick1", 
             main = expression("Inverse non-function:"~ italic(X != f^-1*(Y))))

  sample.matrix <- Untable(func)
  x <- as.numeric(sample.matrix[, 1])
  y <- as.numeric(sample.matrix[, 2])
  
  H.y.given.x <- condentropy( y, x )
  H.x.given.y <- condentropy( x, y )
  
  if(H.y.given.x < H.x.given.y) {
    cat("Conditional entropy picked CORRECT direction X to Y: H(Y|X) =", 
        format(H.y.given.x, digits = 2), "< H(X|Y) =", 
        format(H.x.given.y, digits = 2), "\n")
  } else {
    cat("Conditional entropy picked *WRONG* direction Y to X: H(Y|X) =", 
        format(H.y.given.x, digits = 2), ">= H(X|Y) =", 
        format(H.x.given.y, digits = 2), "\n")
  }

  func.fc.pval <- fun.chisq.test(func)$p.value
  ifunc.fc.pval <- fun.chisq.test(t(func))$p.value
  
  if(func.fc.pval < ifunc.fc.pval) {
    cat("Original FunChisq picked CORRECT direction X to Y: p(X->Y) =", 
        format(func.fc.pval, digits = 2), "< p(Y->X) =", 
        format(ifunc.fc.pval, digits = 2), "\n")
  } else {
    cat("Original FunChisq picked *WRONG* direction Y to X: p(X->Y) =", 
        format(func.fc.pval, digits = 2), ">= p(Y->X) =", 
        format(ifunc.fc.pval, digits = 2), "\n")
  }
  
  func.afc.pval <- fun.chisq.test(func, method="adapted")$p.value
  ifunc.afc.pval <- fun.chisq.test(t(func), method="adapted")$p.value

  if(func.afc.pval < ifunc.afc.pval) {
    cat("Adapted FunChisq picked CORRECT direction X to Y: p(X->Y) =", 
        format(func.afc.pval, digits = 2), "< p(Y->X) =", 
        format(ifunc.afc.pval, digits = 2), "\n")
  } else {
    cat("Adapted FunChisq picked *WRONG* direction Y to X: p(X->Y) =", 
        format(func.afc.pval, digits = 2), ">= p(Y->X) =",
        format(ifunc.afc.pval, digits = 2), "\n")
  }
}

## ----out.width="45%"----------------------------------------------------------
func <- matrix(c(
  1, 5, 1,
  1, 5, 1,
  6, 1, 1
), nrow=3, byrow=TRUE)
compare.methods(func)

## ----out.width="45%"----------------------------------------------------------
func <- matrix(c(
  1, 1, 6, 1,
  1, 1, 6, 1,
  1, 6, 1, 1
), nrow=3, byrow=TRUE)
compare.methods(func)

## ----out.width="45%"----------------------------------------------------------
func <- matrix(c(
  1, 5,
  1, 6,
  4, 1
), nrow=3, byrow=TRUE)
compare.methods(func)

## ----fig.show="hold", out.width="45%"-----------------------------------------
# approximately functional pattern
func = matrix(c(
  8, 0, 1, 0, 
  1, 8, 1, 1, 
  8, 1, 1, 0
), nrow=3, ncol=4, byrow=TRUE)

# empirically independent (constant) pattern
const = matrix(c(
  9, 0, 0, 0, 
  11, 0, 0, 0, 
  10, 0, 0, 0
), nrow=3, ncol=4, byrow=TRUE)

# calculate the AdpFunChisq P-value of 'func' and 'const' tables, lower the better.
func_afc = fun.chisq.test(func, method="adapted")
const_afc = fun.chisq.test(const, method="adapted")

# calculate the Conditional Entropy scores of 'func' and 'const' tables, lower the better.
sample.matrix = Untable(func)
x <- as.numeric(sample.matrix[,1])
y <- as.numeric(sample.matrix[,2])
func_H = condentropy(y, x)

sample.matrix = Untable(const)
x <- as.numeric(sample.matrix[,1])
y <- as.numeric(sample.matrix[,2])
const_H = condentropy(y, x)

# plot table with AdpFunChisq and Conditional Entropy Scores
# functional pattern
plot_table(
  func, xlab = 'Y', ylab = 'X', value.cex = 0.9, cex.main=0.9, col="seagreen",
  main=paste0("AdpFunChisq P = ", format(func_afc$p.value, digits = 2),
              "\nConditional Entropy = ", format(func_H, digits = 2)))

# constant pattern
plot_table(
  const, xlab = 'Y', ylab = 'X', value.cex = 0.9, cex.main=0.9, col="firebrick1",
  main=paste0("AdpFunChisq P = ", format(const_afc$p.value, digits = 2),
              "\nConditional Entropy = ", format(const_H, digits = 2)))

