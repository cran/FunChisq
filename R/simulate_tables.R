# simulate_tables : Generates simulated contingency tables (with or without noise) of a given "type".
#
# Parameters details
# n              : Expected sample size, should be atleast nrow for "functional", "many.to.one", "discontinuous"
#                  and nrow * ncol for "independent" and "dependent.non.functional" [DEFAULT = 100].
# nrow           : Expected number of rows in the table, should be greater than 2 (greater than 3 for type = "many.to.one") [DEFAULT = 3].
# ncol           : Expected number of columns in the table , should be greater than 2 [DEFAULT = 3].
# types include  : functional: y is a function of x but x may or may not be a function of y [DEFAULT] ,
#                  many.to.one: y is a function of x and x is not a function of y,
#                  dependent.non.functional: Non functional table with statistical dependency,
#                : discontinuous: y is a duscontinuous function of x but x may or may not be function of y,
#                  independent: Independent tables (null population).
# n.tables       : number of tables to be generated.[DEFAULT = 1].
# row.marginal   : A vector of row probabilities [DEFAULT = Equally likely].
# col.marginal   : A vector of column probabilities (Only required for type = "independent") [DEFAULT = Equally likely].
# noise          : Factor of noise to be added to the table between 0 and 1 (uses House noise model, 0 means no noise) [DEFAULT = 0].
#
# Created by     : Ruby Sharma
# Date           : October 16 2016
# Last Modified  : April 27 2017
# Version        : 0.0.2
# Updates        : Added pattern type "discontinuous"
#                : Added noise.model parameter


simulate_tables <- function(n=100, nrow=3, ncol=3,
                            type = c("functional",  #=F
                                     "many.to.one", #= MO
                                     "discontinuous",#=DS
                                     "independent", #=I
                                     "dependent.non.functional" #=NF
                            ), noise.model = c("house","candle"),
                            noise=0.0, n.tables=1,
                            row.marginal=rep(1/nrow, nrow),
                            col.marginal= rep(1/ncol, ncol))
{
  type <- match.arg(type)
  noise.model <- match.arg(noise.model)

  prelim.check(nrow, ncol, type, n, noise, row.marginal, col.marginal)

  if(ncol<2 || nrow<2)
    stop("Number of rows and/or columns should atleast be 2")

  if(length(row.marginal)<nrow)
    stop("Row probabilites for all rows expected")

  if(length(which(row.marginal!=0))<2)
    stop("Atleast two non-zero row probabilites expected")

  row.marginal=row.marginal/sum(row.marginal)
  col.marginal=col.marginal/sum(col.marginal)


  #Intialization
  pattern.list = list()
  sample.list = list()
  noise.list = list()
  p.value.list = list()

  for(i in 1:n.tables)
  {
    alltables = table.generate(nrow, ncol, type, n, noise,
                               noise.model, row.marginal, col.marginal)
    pattern.list[[i]] = alltables$pattern.table
    sample.list[[i]] = alltables$sampled.table
    noise.list[[i]] = alltables$noise.table
    p.value.list[[i]] = alltables$p.value
  }


  #return a list of pattern table, sampled contingency table and noise table.

  list(pattern.list = pattern.list, sample.list = sample.list,
       noise.list = noise.list, pvalue.list = p.value.list)
}


#generating pattern table, sampled contingency table and noise table
table.generate=function(nrow, ncol, type, n, noise, noise.model, row.marginal, col.marginal)
{

  if(type=="dependent.non.functional"){

    if(n<(nrow*ncol))
      stop(paste("For dependent.non.functional, n should be greater than or equal to",(nrow*ncol)))

    sam.val.row = sample.in.rows(n, row.marginal, type, ncol)
    pattern.table = nonfunctional.table(nrow, ncol, row.marginal, sam.val.row)
    prob.table = non.functional.prob(nrow, ncol, pattern.table)
    sampled.table = dis.sample.prob(nrow, ncol, sam.val.row, prob.table)
    sampled.table = is_dependent(n, nrow, ncol, row.marginal, col.marginal, sampled.table, sam.val.row)
    p.val = chisq.test.pval(sampled.table)

  } else if (type=="many.to.one") {

    if(nrow<3)
      stop("For many.to.one, number of rows should be at least be 3!")

    if(length(which(row.marginal!=0))<3)
      stop("For many.to.one, at least three non-zero row probabilities are expected!")

    if(n<nrow)
      stop(paste("n should be greater than or equal to", nrow))

    sam.val.row = sample.in.rows(n, row.marginal, type, ncol)
    pattern.table = many.to.one.table(nrow, ncol, row.marginal, sam.val.row)
    sampled.table = dis.sample.prob(nrow, ncol, sam.val.row, pattern.table)
    p.val=FunChisq::fun.chisq.test(sampled.table)$p.value

  } else if (type=="independent") {

    if(length(col.marginal)<ncol)
      stop(paste("For independent, exactly",ncol,"column probabilites expected"))

    if(n<(nrow*ncol))
      stop(paste("For independent, n should be greater than or equal to",(nrow*ncol)))

    sam.val.row = sample.in.rows(n, row.marginal, type, ncol)
    pattern.table = independent.table(nrow, ncol, row.marginal, col.marginal, sam.val.row)
    prob.table = indep.prob.table(nrow, ncol, row.marginal, col.marginal)
    sampled.table = dis.sample.prob(nrow, ncol, sam.val.row, prob.table)
    p.val = chisq.test.pval(sampled.table)

  } else if(type=="discontinuous"){
    if(n<nrow)
      stop(paste("n should be greater than or equal to", nrow))

    sam.val.row = sample.in.rows(n, row.marginal, type, ncol)
    pattern.table = discontinuous.table(nrow, ncol, row.marginal, sam.val.row)
    sampled.table = dis.sample.prob(nrow, ncol, sam.val.row, pattern.table)
    p.val=FunChisq::fun.chisq.test(sampled.table)$p.value

  }else {

    if(n<nrow)
      stop(paste("n should be greater than or equal to",nrow))

    sam.val.row = sample.in.rows(n, row.marginal, type, ncol)
    pattern.table = functional.table(nrow, ncol, row.marginal, sam.val.row)
    sampled.table = dis.sample.prob(nrow, ncol, sam.val.row, pattern.table)
    p.val=FunChisq::fun.chisq.test(sampled.table)$p.value

  }

  if(noise==0){
    noise.table = sampled.table
  } else {
    noise.table = add.noise(tables = sampled.table, u = noise, noise.model = noise.model, margin = 1)
  }


  # return singular tables
  list(pattern.table = pattern.table, sampled.table = sampled.table, noise.table = noise.table, p.value = p.val)
}

# distributing samples across rows guided by row probabilities
sample.in.rows=function(n, row.marginal, type, ncol)
{
  non.zero.rows.length = length(which(row.marginal!=0))
  if(type!="independent") {

    n = n-non.zero.rows.length
    sam.val = rmultinom(1, n, row.marginal)
    ind = which(row.marginal!=0)
    sam.val[ind] = sam.val[ind]+1


  } else {

    n = n-(non.zero.rows.length*ncol)
    sam.val = rmultinom(1, n, row.marginal)
    ind = which(row.marginal!=0)
    sam.val[ind] = sam.val[ind]+ncol

  }

  return(sam.val)
}


# distribute samples
dis.sample.prob=function(nrow, ncol, sam.val.row, table)
{
  sampled.table=as.data.frame(matrix(0, nrow = nrow, ncol = ncol))
  for(i in 1:nrow)
  {
    # determine if all columns in the ith row of the supplied table are zero
    all.zero.column = all(table[i,]==0)
    if(!all.zero.column)
    {
      non.zero.columns = which(table[i,]!=0)
      size = sam.val.row[i]-length(non.zero.columns)
      sam.val.cell = rmultinom(1, size, table[i,])
      sam.val.cell[non.zero.columns] = sam.val.cell[non.zero.columns]+1
      sampled.table[i,] = sam.val.cell
    }
  }

  return(sampled.table)
}


# generating pattern table for "functional"
functional.table=function(nrow, ncol, row.marginal, sam.val.row)
{
  pattern.table = as.data.frame(matrix(0,ncol = ncol,nrow = nrow))

  for(i in 1:nrow)
  {
    if(sam.val.row[i]!=0) {

      index = sample(1:ncol,1)
      pattern.table[i,index] = 1

    }
  }

  # check for constant functions
  pattern.table=not_constant(ncol, pattern.table)

  return(pattern.table)
}


# generating pattern table for "many.to.one"
many.to.one.table=function(nrow, ncol, row.marginal, sam.val.row)
{

  pattern.table = as.data.frame(matrix(0,ncol = ncol,nrow = nrow))
  # get the functional table
  pattern.table = functional.table(nrow, ncol, row.marginal, sam.val.row)
  # check for non-monotonicity
  pattern.table = is_many.to.one(pattern.table)

  return(pattern.table)
}


# generating pattern table for "independent"
independent.table=function(nrow, ncol, row.marginal, col.marginal, sam.val.row)
{

  pattern.table = as.data.frame(matrix(0,ncol = ncol,nrow = nrow))
  prob.table = indep.prob.table(nrow, ncol, row.marginal, col.marginal)
  indexes = non.zero.index(prob.table)
  rows = indexes$rows
  cols = indexes$cols

  for(i in 1:length(rows))
    pattern.table[rows[i],cols[i]] = 1

  return(pattern.table)
}


# generating pattern table for "dependent.non.functional"
nonfunctional.table=function(nrow, ncol, row.marginal, sam.val.row)
{

  pattern.table = as.data.frame(matrix(0,ncol = ncol,nrow = nrow))

  for(i in 1:nrow)
  {
    if(sam.val.row[i]!=0) {

      if(sam.val.row[i]<ncol){
        index = sample(1:ncol,sample(1:sam.val.row[i],1))
      } else {
        index = sample(1:ncol,sample(1:ncol,1))
      }

      pattern.table[i,index] = 1
    }
  }
  pattern.table = make.non.functional(pattern.table, ncol, sam.val.row, nonfunc = "notf.x")
  pattern.table = make.non.functional(t(pattern.table), ncol, sam.val.row, nonfunc = "notf.y")
  pattern.table = t(pattern.table)
  return(pattern.table)
}

# generating pattern table for "discontinuous function"
discontinuous.table=function(nrow, ncol, row.marginal, sam.val.row)
{
  pattern.table = as.data.frame(matrix(0,ncol = ncol,nrow = nrow))

  sample.from =seq_along(1:ncol(pattern.table))
  prev.col.ind = NULL

  for(i in 1:nrow)
  {
    if(sam.val.row[i]!=0) {

      index = sample(sample.from,1)
      sample.from =seq_along(1:ncol(pattern.table))
      prev.col.ind = index
      pattern.table[i,index] = 1
      sample.from = sample.from[-which(sample.from==prev.col.ind)]


    }
  }
  # check for constant functions
  pattern.table=not_constant(ncol, pattern.table)
  return(pattern.table)
}

# introducting atleast two entries in one row for nonfunctional table
make.non.functional=function(pattern.table, ncol, sam.val.row, nonfunc)
{
  indexes = non.zero.index(pattern.table)
  rows = indexes$rows
  cols = indexes$cols


  if(nonfunc == "notf.x")
  {
    if(anyDuplicated(rows)==0) {
      only.row = which(sam.val.row>1)
      chng.from = c(only.row)
      if(length(chng.from)==1)
      {
        chng.row.index = chng.from
      }else
      {
        chng.row.index = sample(chng.from, 1)
      }

      prev.col.index = cols[which(rows==chng.row.index)]
      chng.col.index = sample.sec.col.ind(prev.col.index,ncol)
      pattern.table[chng.row.index,chng.col.index] = 1
    }

  }

  if(nonfunc == "notf.y")
  {
    nr = nrow(pattern.table)
    nc = ncol(pattern.table)
    wth.more.samp <- c()
    wth.less.samp <- c()
    for(i in 1:nc)
    {
      sample.lim <- length(cols[cols[]==i])
      if(sam.val.row[i]>sample.lim)
      {
        wth.more.samp = c(i, wth.more.samp)
      }else{
        wth.less.samp = c(i, wth.less.samp )
      }

    }
    if(anyDuplicated(rows)==0) {


      chng.from = rows
      if(length(chng.from)==1)
      {
        chng.row.index = chng.from
      }else
      {
        chng.row.index = sample(chng.from, 1)
      }

      prev.col.index = cols[which(rows==chng.row.index)]
      chng.col.index = sample.sec.col.ind_notf.y(prev.col.index, nc, wth.more.samp)
      if(chng.col.index == 0)
      {

        if(length(rows[cols[]==wth.less.samp])==1)
        {
          chng.row.index = rows[cols[]==wth.less.samp]
        }else{
          chng.row.index = sample(rows[cols[]==wth.less.samp], 1)
        }
        chng.col.index = wth.more.samp
      }
      pattern.table[chng.row.index,chng.col.index] = 1
    }

  }


  return(pattern.table)
}


# generating probability table for "independent"
indep.prob.table=function(nrow, ncol, row.marginal, col.marginal)
{
  prob.table = as.data.frame(matrix(0,ncol = ncol,nrow = nrow))

  # multiplying row probability and column probability in in each cell of prob.table
  for(i in 1:nrow){
    for(j in 1:ncol)
      prob.table[i,j] = row.marginal[i]*col.marginal[j]
  }

  return(prob.table)
}


# generating probability table for "dependent.non.functional"
non.functional.prob=function(nrow, ncol, pattern.table)
{

  prob.table = as.data.frame(matrix(0,ncol = ncol,nrow = nrow))

  indexes = non.zero.index(pattern.table)
  rows = indexes$rows
  cols = indexes$cols

  for(i in 1:nrow)
  {
    if(is.element(i,rows))
    {
      col.ind = which(rows==i)
      freq.no = length(col.ind)
      prob.ele.row = 1/freq.no
      prob.table[i,cols[col.ind]] = prob.ele.row
    }
  }

  return(prob.table)
}


# check for constant function, if found then change one column index
not_constant=function(ncol, pattern.table)
{

  indexes = non.zero.index(pattern.table)
  rows = indexes$rows
  cols = indexes$cols

  if(length(unique(cols))==1) {

    index = cols[1]
    chng.col.index = sample.sec.col.ind(index, ncol)
    chng.row.index = sample(rows, 1)
    pattern.table[chng.row.index,index] = 0
    pattern.table[chng.row.index,chng.col.index] = 1
  }

  return(pattern.table)
}


# check monotonicity, if found then atleast two rows would share samples in the same column
is_many.to.one=function(pattern.table)
{

  indexes = non.zero.index(pattern.table)
  rows = indexes$rows
  cols = indexes$cols

  if(length(unique(cols))==length(rows)){

    id = sample(rows,2)
    pattern.table[id[2],cols[which(rows==id[2])]] = 0
    pattern.table[id[2],cols[which(rows==id[1])]] = 1
  }

  return(pattern.table)
}
#check whether the dependent.non.functional table is dependent, if not make dependent
is_dependent = function(n, nrow, ncol, row.marginal,  col.marginal, sampled.table, sam.val.row)
{

  sam.val.indep =  sample.in.rows(n, row.marginal, type ="independent", ncol)
  expec.prob.table = as.data.frame(matrix(0,ncol = ncol,nrow = nrow))
  expec.prob.table = indep.prob.table(nrow, ncol, row.marginal, col.marginal)
  indep.sampled.table = dis.sample.prob(nrow, ncol, sam.val.indep, expec.prob.table)
  difference = indep.sampled.table - sampled.table
  if(all(difference==0))
  {
    indexes = non.zero.index(sampled.table)
    rows = indexes$rows
    cols = indexes$cols
    sel.row <- sample(rows, 1)
    sel.from <- cols[which(rows==sel.row)]
    if(length(sel.from)==1)
    {
      sel.colj1 = sel.from
    }else
    {
      sel.colj1 <- sample(sel.from,1)
    }
    sel.from <- 1:ncol
    sel.colj2 <-sample(sel.from[-which(sel.from==sel.colj1)],1)
    sampled.table[sel.row,sel.colj2] <- sampled.table[sel.row,sel.colj2] + sampled.table[sel.row,sel.colj1]
    sampled.table[sel.row,sel.colj1] <- 0
  }

  return(sampled.table)



}



# sorting row and column indexes on the basis of row
sort.index=function(rows, cols)
{
  row.col.ind = as.data.frame(matrix(nrow = length(rows),ncol =2))
  row.col.ind[,1] = rows
  row.col.ind[,2] = cols
  row.col.ind = row.col.ind[order(row.col.ind[,1],row.col.ind[,2]),]
  rows = row.col.ind[,1]
  cols = row.col.ind[,2]

  list(row.in = rows,col.in = cols)
}


# retrieving non zero indexes from the table
non.zero.index=function(table)
{
  rows = row(table)[which(!table == 0)]
  cols = col(table)[which(!table == 0)]
  sorted = sort.index(rows,cols)
  rows = sorted$row.in
  cols = sorted$col.in

  list(rows = rows,cols = cols)
}


# sampling secondary column index
sample.sec.col.ind=function(index, ncol)
{

  if(index==1)
    vec.to.sample = c(2:ncol)

  if(index==ncol)
    vec.to.sample = c(1:(ncol-1))

  if(index>1 && index<ncol)
    vec.to.sample = c(1:(index-1),(index+1):ncol)

  if(length(vec.to.sample)==1)
  {
    chng.col.index=vec.to.sample
  }else{
    chng.col.index = sample(vec.to.sample,1)
  }


  return(chng.col.index)
}

prelim.check = function(nrow, ncol, type, n, noise, row.marginal, col.marginal)
{
  if(class(nrow)!="numeric" && class(nrow)!="integer")
    stop("nrow should be numeric")

  if(class(ncol)!="numeric" && class(ncol)!="integer")
    stop("ncol should be numeric")

  if(class(type)!="character")
    stop("type should be character")

  if(type!="functional" && type!="many.to.one" && type!="dependent.non.functional" &&
     type!="independent"&& type!="discontinuous")
    stop("type can only be functional, many.to.one, dependent.non.functional, independent or discontinuous")

  if(class(n)!="numeric" && class(n)!="integer")
    stop("n should be numeric")

  if(class(noise)!="numeric" && class(noise)!="integer")
    stop("noise should be numeric")

  if(class(row.marginal)!="numeric" && class(row.marginal)!="integer")
    stop("row.marginal should be numeric")

  if(class(col.marginal)!="numeric" && class(col.marginal)!="integer")
    stop("col.marginal should be numeric")
}

chisq.test.pval <- function(table)
{
  # identify non-zero rows and non-zero columns
  non.zero.rows <- apply(table, 1, function(row) { 0 != sum(abs(row)) } )

  non.zero.cols <- apply(table, 2, function(col) { 0 != sum(abs(col)) } )

  # perform Pearson chi-square test
  chisq <- suppressWarnings(chisq.test(table[non.zero.rows, non.zero.cols])$statistic)

  # compute p-value using the orgional table size
  pval <- pchisq( chisq, prod(dim(table) - 1), lower.tail = FALSE)
  names(pval) <- NULL
  return(pval)
}

#sampling second column for making x!=f(y) for dependent.non.functional
sample.sec.col.ind_notf.y <- function(index ,ncol, only.col)
{


  if(index==1)
    vec.to.sample = c(2:ncol)

  if(index==ncol)
    vec.to.sample = c(1:(ncol-1))

  if(index>1 && index<ncol)
    vec.to.sample = c(1:(index-1),(index+1):ncol)

  vec.to.sample =  vec.to.sample[vec.to.sample %in% only.col]

  if(length(vec.to.sample)==0)
  {
    chng.col.index = 0
  }else if(length(vec.to.sample)==1)
  {
    chng.col.index = vec.to.sample
  }else{
    chng.col.index = sample(vec.to.sample,1)
  }
  return(chng.col.index)
}


