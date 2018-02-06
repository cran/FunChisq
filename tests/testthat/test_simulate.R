# tests simulate_tables.R
# Created by: Sajal Kumar
# Modified by : Ruby Sharma
# Date : February 27 2017

library(testthat)
library(FunChisq)


#Attributes to test
#1) y = f(x)
#2) Not a constant function - All populated samples are not supposed to be on the same column.

Test_Functional_table = function(iter)
{

  func.flag = FALSE

  for(i in 1:iter)
  {
    Get.Stats = Construct_Table("functional")

    conti.table = Get.Stats$conti.table
    noise.table = Get.Stats$noise.table

    # check if all samples are populated

    check.all.samples = All.sample.check(conti.table,Get.Stats$sample.size)

    if(check.all.samples$flag)
    {
      func.flag = TRUE
      failure.summary = check.all.samples$failure.table
      break
    }

    # Check for non zero cells

    check.non.zero = No.Non.Zero.Check(conti.table)

    if(check.non.zero$flag)
    {
      func.flag = TRUE
      failure.summary = check.non.zero$failure.table
      break
    }

    # Check y = f(x)

    check.functional = Functional.check(conti.table)

    if(check.functional$flag)
    {
      func.flag = TRUE
      failure.summary = check.functional$failure.table
      break
    }

    # check constant functions

    check.constant = Constant.check(conti.table)

    if(check.constant$flag)
    {
      func.flag = TRUE
      failure.summary = check.constant$failure.table
      break
    }

    check.margin = margin.check(conti.table, noise.table)
    if(check.margin$flag)
    {
      func.flag = TRUE
      failure.summary = check.margin$failure.table
      break
    }


  }

  #if any table was flagged, the test failed.
  expect_identical(func.flag, FALSE)
}


#Attributes to test
#1) y = f(x)
#2) Not a constant function - All populated samples are not supposed to be on the same column.
Test_Functional_Discontinuous_table = function(iter)
{

  func.flag = FALSE

  for(i in 1:iter)
  {
    Get.Stats = Construct_Table("discontinuous")

    conti.table = Get.Stats$conti.table

    # check if all samples are populated

    check.all.samples = All.sample.check(conti.table,Get.Stats$sample.size)

    if(check.all.samples$flag)
    {
      func.flag = TRUE
      failure.summary = check.all.samples$failure.table
      break
    }

    # Check for non zero cells

    check.non.zero = No.Non.Zero.Check(conti.table)

    if(check.non.zero$flag)
    {
      func.flag = TRUE
      failure.summary = check.non.zero$failure.table
      break
    }

    # Check y = f(x)

    check.functional = Functional.check(conti.table)

    if(check.functional$flag)
    {
      func.flag = TRUE
      failure.summary = check.functional$failure.table
      break
    }

    # check constant functions

    check.constant = Constant.check(conti.table)

    if(check.functional$flag)
    {
      func.flag = TRUE
      failure.summary = check.constant$failure.table
      break
    }
  }

  #if any table was flagged, the test failed.
  expect_identical(func.flag, FALSE)
}

#Attributes to test
#1) y = f(x)
#2) x ! = f(y)
#2) Not a constant function - All populated samples are not supposed to be on the same column.

Test_Functional_Many_to_one_table = function(iter)
{

  non.mono.func.flag = FALSE

  for(i in 1:iter)
  {

    Get.Stats = Construct_Table("many.to.one")

    conti.table = Get.Stats$conti.table

    # check if all samples are populated

    check.all.samples = All.sample.check(conti.table,Get.Stats$sample.size)

    if(check.all.samples$flag)
    {
      non.mono.func.flag = TRUE
      failure.summary = check.all.samples$failure.table
      break
    }

    # Check for non zero cells

    check.non.zero = No.Non.Zero.Check(conti.table)

    if(check.non.zero$flag)
    {
      non.mono.func.flag = TRUE
      failure.summary = check.non.zero$failure.table
      break
    }

    # Check y = f(x)

    check.functional = Functional.check(conti.table)

    if(check.functional$flag)
    {
      non.mono.func.flag = TRUE
      failure.summary = check.functional$failure.table
      break
    }

    # check x != f(y)

    check.non.functional.invert = Non.functional.check(t(conti.table))

    if(check.non.functional.invert$flag)
    {
      non.mono.func.flag = TRUE
      failure.summary = conti.table
      break
    }

    # check constant functions

    check.constant = Constant.check(conti.table)

    if(check.constant$flag)
    {
      non.mono.func.flag = TRUE
      failure.summary = check.constant$failure.table
      break
    }
  }

  #if any table was flagged, the test failed.

  expect_identical(non.mono.func.flag, FALSE)
}


#Attributes to test
#1) y ! = f(x)

Test_Non_Functional_table = function(iter)
{

  non.func.flag = FALSE

  for(i in 1:iter)
  {

    Get.Stats = Construct_Table("dependent.non.functional")

    conti.table = Get.Stats$conti.table

    # check if all samples are populated

    check.all.samples = All.sample.check(conti.table,Get.Stats$sample.size)

    if(check.all.samples$flag)
    {
      non.func.flag = TRUE
      failure.summary = check.all.samples$failure.table
      break
    }

    # Check for non zero cells

    check.non.zero = No.Non.Zero.Check(conti.table)

    if(check.non.zero$flag)
    {
      non.func.flag = TRUE
      failure.summary = check.non.zero$failure.table
      break
    }

    # check y != f(x)

    check.non.functional = Non.functional.check(conti.table)

    if(check.non.functional$flag)
    {
      non.func.flag = TRUE
      failure.summary = conti.table
      break
    }
  }

  #if any table was flagged, the test failed.
  expect_identical(non.func.flag, FALSE)
}

#Attributes to test
#1) No non zero cell
#2) No dependency
#3) y != f(x) AND x != f(y)

Test_Independent_table = function(iter)
{

  ind.flag = FALSE

  for(i in 1:iter)
  {

    Get.Stats = Construct_Table("independent")

    conti.table = Get.Stats$conti.table

    # check if all samples are populated

    check.all.samples = All.sample.check(conti.table,Get.Stats$sample.size)

    if(check.all.samples$flag)
    {
      ind.flag = TRUE
      failure.summary = check.all.samples$failure.table
      break
    }

    # check for all non zero cells

    check.all.non.zero = All.Non.Zero.Check(conti.table)

    if(check.all.non.zero$flag)
    {
      ind.flag = TRUE
      failure.summary = check.all.non.zero$failure.table
      break
    }

    # check for no dependency (x and y are independent)

    check.dependency = Dependency.check(conti.table)

    if(!check.dependency$flag)
    {
      ind.flag = TRUE
      failure.summary = check.dependency$failure.table
      break
    }

    # check y != f(x)

    check.non.functional = Non.functional.check(conti.table)

    if(check.non.functional$flag)
    {
      non.func.flag = TRUE
      failure.summary = conti.table
      break
    }

    # check x != f(y)

    check.non.functional.invert = Non.functional.check(t(conti.table))

    if(check.non.functional.invert$flag)
    {
      ind.flag = TRUE
      failure.summary = conti.table
      break
    }
  }

  #if any table was flagged, the test failed.
  expect_identical(ind.flag, FALSE)
}


#Setup random parameters and construct a table using simulate.tables
Construct_Table = function(type)
{
  row_col.equality = sample(c(TRUE,FALSE),1)
  if(row_col.equality) {
    nrows = sample(c(2:15),1)
    ncols = nrows
  } else {
    nrows = sample(c(2:15),1)
    ncols = sample(c(2:15),1)
  }

  if(type=="many.to.one" && nrows == 2)
    nrows = 3

  sample.size = ((sample(c(1:100),1) * sample(c(1:100),1))) + nrows
  row.marginal = rep(0,nrows)

  if(type=="independent" && sample.size < (nrows * ncols))
    sample.size = nrows * ncols

  if(type=="dependent.non.functional" && sample.size < (nrows * ncols))
    sample.size = nrows * ncols

  row.marginal.set = sample(c(TRUE,FALSE),1)
  if(row.marginal.set){
    row.marginal = runif(nrows)
    row.marginal = row.marginal/sum(row.marginal)
    conti.table = simulate_tables(n=sample.size, nrow = nrows, ncol = ncols, type = type, n.tables = 1, row.marginal = row.marginal, noise = 0.2)
  } else {
    conti.table = simulate_tables(n=sample.size, nrow = nrows, ncol = ncols, type = type, n.tables = 1, noise =0.2)
  }

  list(conti.table = conti.table$sample.list[[1]], nrows = nrows, ncols = ncols, sample.size = sample.size, row.marginal = row.marginal, noise.table = conti.table$noise.list[[1]])
}

# The table should contain all the samples specified in sample size
All.sample.check = function(conti.table, sample.size)
{
  flag=FALSE
  failure.summary = NULL

  if(sum(conti.table)!=sample.size){
    flag = TRUE
    failure.summary = conti.table
  }

  list(flag=flag,failure.table = failure.summary)
}

# If there are no non-zero elements in the table then something's wrong.
No.Non.Zero.Check = function(conti.table)
{
  index = which(conti.table!=0, arr.ind=TRUE)

  flag=FALSE
  failure.summary = NULL

  if(length(index)==0){
    flag = TRUE
    failure.summary = conti.table
  }

  list(flag=flag, failure.table = failure.summary)
}

# If there are any zero element in the table then something's wrong.
All.Non.Zero.Check = function(conti.table)
{
  index = which(conti.table==0, arr.ind=TRUE)

  flag=FALSE
  failure.summary = NULL

  if(length(index)!=0){
    flag = TRUE
    failure.summary = conti.table
  }

  list(flag=flag, failure.table = failure.summary)
}

# Testing y = f(x)
#If there exists a row such that it contains samples for more than one y then the table is invalid
Functional.check = function(conti.table)
{
  index = which(conti.table!=0, arr.ind=TRUE)

  flag=FALSE
  failure.summary = NULL

  if(length(index[,1])!=length(unique(index[,1]))){
    flag = TRUE
    failure.summary = conti.table
  }

  list(flag=flag, failure.table = failure.summary)
}

# Testing y != f(X)
Non.functional.check = function(conti.table)
{

  status = Functional.check(conti.table)

  flag=FALSE
  failure.summary = NULL

  if(!status$flag)
    flag = TRUE

  list(flag=flag, failure.table = failure.summary)
}

# If there exists a column (or y) such that it contains all the samples then the table is invalid
Constant.check = function(conti.table)
{
  index = which(conti.table!=0, arr.ind=TRUE)

  flag=FALSE
  failure.summary = NULL

  if(length(unique(index[,2]))==1){
    flag=TRUE
    failure.summary = conti.table
  }

  list(flag=flag, failure.table = failure.summary)
}


#If there exists a zero cell after removing all zero rows and columns then the table is invalid
Dependency.check = function(conti.table)
{
  index.zero = which(conti.table==0, arr.ind=TRUE)
  index.non.zero = which(conti.table!=0, arr.ind=TRUE)

  flag=FALSE
  failure.summary = NULL

  if(length(index.zero) == 0)
  {
    flag = TRUE
    failure.summary = conti.table

  } else {

    index.zero.copy = index.zero
    for(j in 1:length(index.zero.copy[,1]))
    {

      if(length(which(index.zero[,1]==index.zero.copy[j,1]))==ncol(conti.table)){
        index.zero = index.zero[-which(index.zero[,1]==index.zero.copy[j,1]),,drop=FALSE]
      }

      if(length(which(index.zero[,2]==index.zero.copy[j,2]))==nrow(conti.table)){
        index.zero = index.zero[-which(index.zero[,2]==index.zero.copy[j,2]),,drop=FALSE]
      }
    }

    if(nrow(index.zero)==0){
      flag = TRUE
      failure.summary = conti.table
    }
  }

  list(flag=flag, failure.table = failure.summary)
}

margin.check = function(conti.table, noise.table)
{
  flag=FALSE
  failure.summary = NULL


    if(any(rowSums(conti.table)!= rowSums(noise.table)))
    {
      flag = TRUE
      failure.summary = noise.table
    }
  list(flag=flag, failure.table = failure.summary)

}

test_that("Testing the Simulate_tables()", {
  Test_Functional_table(2)
  Test_Independent_table(2)
  Test_Non_Functional_table(2)
  Test_Functional_Many_to_one_table(2)
  Test_Functional_Discontinuous_table(2)
})
