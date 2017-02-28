## Hua created at Mar 3 2014.
## Adding noise to contingency table by using House noise model.
## House noise model was designed by Yang Zhang and Joe Song.
## Note: The adding noise code only supports on parent and one child, does not support combinatorial parents.
##
## Hua modified at Feb 14, 2017.
## Updates:
## 1. Parameter "all.tables" must be A list of tables or one table (matrix/data frame).
## 2. u is the noise level between [0, 1].
## 3. Allow adding noise to one margin only (X/Y).
##    Added parameter margin with values between [0,1,2]. Default is 0.
##       0: add noise to both row and column margins of a table;
##       1: add noise to row margin. Row sums are fixed;
##       2: add noise to column margin. Column sums are fixed.
##

add.house.noise <- function(tables, u, margin=0) {
  if(is.matrix(tables) || is.data.frame(tables)){
    tables.noised <- add.noise.one.table(tables, u, margin)
  }else{
    tables.noised <- lapply(1:length(tables),
                            function(k){
                              t <- tables[[k]]
                              t <- add.noise.one.table(t, u, margin)
                              return (t)
                            }
    )
  }
  return (tables.noised)
}

add.noise.one.table <- function(one.table, u, margin=0){
  # margin: 0, both X and Y; 1, change along X (in a column); 2, change along Y (in a row).
  t <- one.table
  table.result <- NULL

  if(margin == 0){
    t.row <- nrow(t)
    t.col <- ncol(t)

    tables <- lapply(1:(t.row * t.col), function(x) {
      row <- ceiling(x/t.col)
      col <- (x-1)%%t.col+1

      t.tmp <- lapply(1:t[row, col], function(y){
        t.one.sample <- matrix(0, nrow=t.row, ncol=t.col)
        if(t[row, col]==0)return (t.one.sample)

        ##########
        #Add noise
        row.new <- row
        col.new <- col
        if(t.row!=1){
          row.new <- house.noise.model.to.node(t.row, u, row)#Index of rows, jump in columns
        }else{
          row.new <- row
        }
        if(t.col!=1){
          col.new <- house.noise.model.to.node(t.col, u, col)#Index of columns, jump in rows
        }else{
          col.new <- col
        }

        t.one.sample[row.new, col.new] <- 1
        ##########

        return (t.one.sample)
      })
      return (Reduce('+', t.tmp))
    })
    table.result <- Reduce('+', tables)
  } else if(margin == 1) {
      table.result <- t(apply(t, margin, function(x){
        t.tmp <- as.matrix(x)
        add.noise.one.table(t.tmp, u, 0)
      }))
  } else if(margin == 2){
     table.result <- apply(t, margin, function(x){
       t.tmp <- as.matrix(x)
       add.noise.one.table(t.tmp, u, 0)
     })
  }

  # else if(margin == 1){
  #   table.result <- apply(t, 2,function(x){
  #     t.tmp <- as.matrix(x)
  #     add.noise.one.table(t.tmp, u, 0)
  #   })
  # }else if(margin == 2){
  #   table.result <- t(apply(t, 1, function(x){
  #     t.tmp <- as.matrix(x)
  #     add.noise.one.table(t.tmp, u, 0)
  #   }))
  # }

  else{
    stop("Wrong margin values! ([0,1,2])")
  }
  return (table.result)
}

house.noise.model.to.node <- function(baseNum, u, B){#u: mnoise level, B: current base level
  #Refer to the C++ code in GLN, the code is translated to here in R
  prvector <- matrix (0, nrow=baseNum, ncol=baseNum)
  tempsum <- sapply(1:baseNum, simplify="array", USE.NAMES=FALSE, function(x){
    a <- c(1:baseNum)
    a <- abs(a-x)
    return (sum(a))
  })
  prvector <- t(apply(matrix(1:baseNum, ncol=1), 1, function(x){
    a <- c(1:baseNum)
    res <- (1 - abs(x - a)/tempsum[x]) * u/(baseNum-1)
    res[x] <- res[x]+1-u
    return (res)
  }))
  prvector <- prvector * (1-u) + u/baseNum

  valuevector <- matrix (0, nrow=baseNum, ncol=baseNum)
  valuevector <- apply(matrix(1:baseNum, nrow=1), 2, function(x){
    if(x==1){
      return (prvector[,1:x])
    }else{
      return (rowSums(prvector[,1:x]))
    }

  })
  valuevector[,baseNum] <- 1

  randomNumber <- runif(1, 0, 1)
  valuevector.B <- valuevector[B,]
  index <- c(1:length(valuevector.B))
  return (index[randomNumber<=valuevector.B][1])
}


# ####Test: random k 5*5 test cases.
# house.noise.model.test.case <- function(k=3){
#   # A list of k tables
#   t <- lapply(1:k, function(x){matrix(sample.int(5, 25, TRUE), 5, 5)})
#   t.XY <- add.house.noise(t, 0.5, 0)
#   t.X <- add.house.noise(t, 0.5, 2)
#   t.Y <- add.house.noise(t, 0.5, 1)
#
#   t.res <- lapply(1:length(t), function(x){
#     return(
#       sum(t[[x]]) == sum(t.XY[[x]]) &
#         all(colSums(t[[x]]) == colSums(t.X[[x]])) &
#         all(rowSums(t[[x]]) == rowSums(t.Y[[x]]))
#     )
#   })
#   t.res <- unlist(t.res)
#
#
#   # One table
#   t <- matrix(sample.int(5, 25, TRUE), 5, 5)
#   t.XY <- add.house.noise(t, 0.5, 0)
#   t.X <- add.house.noise(t, 0.5, 2)
#   t.Y <- add.house.noise(t, 0.5, 1)
#
#   t.res <- c(t.res,
#              sum(t) == sum(t.XY) &
#                all(colSums(t) == colSums(t.X)) &
#                all(rowSums(t) == rowSums(t.Y)))
#
#   t.res <- unlist(t.res)
#   if(all(t.res)){
#     message("All test cases passed!")
#   }else{
#     message(paste(length(t.res)-sum(t.res), '/', length(t.res), " test cases failed!", sep=''))
#   }
# }
