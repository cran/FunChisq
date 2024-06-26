# NEWS

## Version 2.5.4

2024-05-10

1. Created version 2.5.4 from 2.5.3.2.


## Version 2.5.3.2

2024-05-10

  1. Add dqRNGkind("xoroshiro128+") to test code for incoming 'dqrng' package's change of default random number generator.
  2. Update the test cases to be consistent with testthat (>= 3.0.0).


2023-09-09

Major changes were made for the exact function test (EFT). The EFT DQP code has been rewritten to fix bugs that might cause crash for certain inputs:

1. Implemented the DQP algorithm using hash table instead of vector.
2. Used row sums as the key to build the search tree. It fixed one bug in DQP.
3. Renamed the variables to avoid reversing rows and columns in the input table.
4. Renamed the header and cpp files of the DQP algorithm to fastEnu.h and fastEnu.cpp.
5. Edited the main.cpp file to implement the EFT using revised code.

2023-09-09

1. Created version 2.5.3.2 from 2.5.3.1.

## Version 2.5.3.1

2023-06-21

1. Created version 2.5.3.1 from 2.5.3.
2. Added seed functions to reproduce simulated tables in the test_simulate.R

## Version 2.5.3

2023-05-25

1. Created version 2.5.3 from 2.5.2.3.
2. Modified README.md and FunChisq-package.Rd.

## Version 2.5.2.3 (not publicly released)

2023-05-25

1. Renamed sort.index() to sort_index() in simulate_tables.R to avoid conflict with S3 methods of sort().

2023-05-22

1. Changed NEWS from plain text file to markdown file.
2. Edited NEWS.

2023-05-20

1. Created 2.5.2.3 from 2.5.2.2
2. Fixed a bug in the adapted Funchisq test. The R code cannot get column number from an integer. Fixed by changing the integer to matrix data type.

## Version 2.5.2.2 (not publicly released)

2023-05-19

1. Created 2.5.2.2 from 2.5.2.1
2. Modified test_simulate.R to examine the pattern tables for the maintenance of table type property.
3. Modified the simulate_table.R, not_constant() function by commenting out random table generation part (it was eliminating the table type property).

## Version 2.5.2.1 (not publicly released)

2023-05-19

1. Dropped the C++11 specification.
2. Added adapted the functional chi-squared test reference (Kumar and Song, 2022) in DESCRIPTION and fun.chisq.test() manual.

2022-01-18

1. Updated the README.md file.

2021-07-18

1. Created 2.5.2.1 from 2.5.2.
2. Modified simulate_tables() function for usage of seq(n) to seq(length=n), which gives different results when n is 0.
3. Revised the manual for simulate_tables() function regarding marginal distributions.

## Version 2.5.2

2021-05-19

1. Created 2.5.2 from 2.5.1.2

## Version 2.5.1.2 (not publicly released)

2021-05-17

1. Updated the is_many.to.one() function in simulate_tables()
to satisfy the non-monotonic property of functional patterns within the rows and columns having non-zero marginals.
2. Added a new method option "adapted" to fun.chisq.test.
3. Added a vignette 'adapted.fun.chisq.test.Rmd' for adapted functional chi-squared test.
4. Mentioned "adapted functional chi-squared test" in DESCRIPTION.
5. Imported R packages 'dqrng', 'DescTools' and 'infotheo' in DESCRIPTION.
6. Added three testthat test cases for adapted functional chi-squared test in test_AdpFunChisq.R.
7. Updated fun.chisq.test.Rd with details about the method option "adapted".
8. Added imports from 'dqrng', 'infotheo' and 'DescTools' in NAMESPACE.

2021-05-11

1. Created 2.5.1.2 from 2.5.1.1
2. Updated the not_constant() function in simulate_tables() to satisfy the non-constant property of functional patterns within the rows and columns having non-zero marginals.

## Version 2.5.1.1 (not publicly released)

2021-04-18

1. Removed the 'LazyData: TRUE' option in the package DESCRIPTION as the package does not have a 'data' directory.

2021-03-26

1. Further improved numerical stability for calculating FunChisq statistic and function index (estimate) when alternative = "non-constant" (default) for fun.chisq.test().
2. Added two more testthat cases for in test_FunChisq.R

2021-03-24

1. Fixed a bug that if the input has all zeros, the
fun.chisq.test() function would stop with an error. A test case is added to the testthat examples.
2. Fixed a bug that if an input table has small non-integer values, round-off errors can result in calculating function index (estimate), now the issue has been fixed. A test case is added to the testthat examples.

2020-11-27

1. Fixed the p.value entry to be unamed in the return object of fun.chisq.test(). It was incorrectly named "statistic" in some test options. The change was made inside file "R/FunChisq.R".

2020-07-18

1. Updated REFERENCES.bib and CITATION
2. Updated README.md

2020-06-14

1. Created 2.5.1.1 from 2.5.1

## Version 2.5.1

2020-06-13

1. Created 2.5.1 from 2.5.0.1

## Version 2.5.0.1 (not publicly released)

2020-06-13

1. Updated the README file to include badges about the package.
2. Replaced \dontrun by \donttest for examples in manuals.

2020-05-31

1. API change to simulate_tables(). For functional patterns, if both row and column marginal distributions are specified, each marginal is observed half of the time with the other marginal not strictly observed. In the past, only the row marginal distribution was observed.
2. The simulate_tables() function now can generate functional patterns for a given column marginal distribution.
3. simulate_tables(): Default row and column marginals are changed to NULL at the parameter.
Any NULL is replaced by uniform distribution if not specified by user.
4. simulate_tables(): Introduced decide.marginal() to select the
marginal type if both are provided or none of them are provided.
5. simulate_tables(): Changed the name of functional.table() function to disc.functional.pattern().
6. simulate_tables(): Functional 0/1 pattern table is no longer dependent on either row and column marginal.
7. Added new testthat tests for simulate_tables() to check for marginal distribution of simulated tables.
8. Added three more examples of functional tables in the simulate_tables() manual with user specified column marginals.
9. Updated simulate_tables() manual.

2020-04-24

1. Created Version 2.5.0.1 from 2.5.0.
2. Updated manuals of most functions in the package.

## Version 2.5.0

2020-04-23

1. API change: simulate_tables() now includes a margin parameter to specify how noise model will be applied. The default noise for functional patterns has now been changed to apply along both row and column. In previous versions, the noise was only applied along columns within each row such that row sums do not change.
2. Test_Functional_table() has been modified to not check for equal row sums before and after adding noise, which is true only if the noise is applied along the column.
3. Updated references, CITATION, and manuals to reflect a newly accepted manuscript (Nguyen et al. 2020).

2020-03-15

1. Fixed a NOTE in 2.4.9.2 Rdpack was not imported in NAMESPACE.
2. Updated references in FunChisq-package.Rd and REFERENCES.bib.

2020-03-14

1. Created version 2.5.0 from 2.4.9.2.

## Version 2.4.9.2

2020-03-13

1. Updated CITATION, DESCRIPTION, README.md
2. Introduced the Rdpack to manage references

2020-02-13

1. Added ORCID for author Sajal Kumar.

2019-11-10

1. Updated the REAMDE.md file of the package.

2019-11-04

2. Fixed "page" entry to "pages" entry in the CITATION file.
3. Deprecated 'cp.chisq.test'. The function has been moved to R package 'DiffXTables'.

2019-09-28

1. Created version 2.4.9.2.
2. Fixed signed / unsigned mismatch in trimTable.cpp.

## Version 2.4.9.1

2019-09-20

1. Created version 2.4.9.1.
2. Fixed a bug in trimming all-zero rows or columns from the input table in EFT DQP and DP. Added a new test case x16 for this bug.
3. Eliminated the use of static array in the EFT DP and DQP code.

## Version 2.4.9 (not publicly released)

2019-09-08

1. Created version 2.4.9.
2. Revised README.md.
3. Fixed a bug for EFT DQP. Included the table that caused the bug as test case table x15.

## Version 2.4.8-1

2019-09-08

1. Revised 'DESCRITION'.
2. Revised 'README.md'.

2019-09-07

1. Fixed a bug in the exact.dqp option of the exact functional test. Added a new test case that had given this bug in the test functions.
2. Re-organized the test cases into more R files based on 'testthat'.
3. Added 'README.md'.
4. Included references into DESCRIPTION.

2019-09-06

1. Default index.kind is changed from uncoditional to function index in function test.interactions(). (API change)
2. Edited documentation.
3. Changed the package title to "Model-Free Functional Chi-Squared and Exact Tests"

2019-09-05

1. Created version 2.4.8-1.
2. Updated package documentation (DESCRIPTION, CITATION). Now included a reference to Hien Nguyen's unpublished dissertation.
3. Updated manuals for the package.
4. Updated vignette for fun.chisq.test.Rmd.
5. Commented out std::cout usage in Node::show() function defined in Node.cpp.
6. Changed to use chi-squared test instead of chi-square test.
7. Created a new manual (EFT.Rd) for exact functional test functions EFTDP() and EFTDQP().
8. Updated the table and sample size restrictions to the exact test options in function fun.chisq.test().

## Version 2.4.8 (not publicly released)

2018-12-25

1. Provide new method options "exact.qp", "exact.dp", and "exact.dqp" for fun.chisq.test().

2018-12-06

1. Created version 2.4.8 to merge the exact tests in 2.4.6 and 2.4.7.

## Version 2.4.7 (not publicly released)

2018-07-24

1. Improved the runtime of the exact functional test by applying the quadratic lower bound and upper bounds for each node at the step of network construction during branch-and-bound. This sped up the test thousands of times from the quadratic-programming-only or dynamic-programming-only option on all contingency tables.

## Version 2.4.6 (not publicly released)

2018-06-17

1. Added an alternative implementation of the exact functional test using branch-and-bound with dynamic programming. This reduced runtime on some contigency tables with large table size or sample size by thousands of times.
2. Added a new method option (method = "exact.dp") to function fun.chisq.test() in "R/FunChisq.R" for the above new exact functional test implementation.

## Version 2.4.5-3

2018-12-04

1. Updated simulate_tables() function to be more accurate when generating small tables where row and column variables are independent. Specifically, zero counts are now allowed.
2. Tests for simulate_tables() in test_simulate.R is also updated.
3. Introduced simulate_independent_tables() from the preivous code, where noise is applied along both rows and columns (previously only along the rows). The funciton is in simulate_tables.R.
4. Modified prelim.check() in simulate_tables.R.
5. The manual page for simulate_tables() is updated
accordingly.
6. Updated the examples in the manual of plot_tables().

2018-11-01
1. Added an argument in plot_table() to highlight row maxima with a box.

2018-10-24
1. Created version 2.4.5-3.
2. Make "conditional function index" to be the default, previously "unconditional".
3. Added a new vignette "Examples of discrete patterns".
4. Set the default function index.kind to "conditional". The "unconditional" option for function index will phase out.
5. Updated vignettes.

## Version 2.4.5-2

2018-10-22

1. Reintroduced a test source code file test_FunChisq.R back, which was accidentally deleted from the previous version.

2018-10-17

1. Added a plot_table() function to visualize a contingency table.

2018-10-16

1. Created version 2.4.5-2.
2. Included conditional FunChisq R code.

## Version 2.4.5-1

2018-06-16

1. In the exact functional test C++ code, improved the numerical precision of the FunChisq statistic using an equivalent mathematical form already used in the R version. This fixed bugs in the exact functional test when the table is of certain dimension and sample size.
2. Added a few more test cases for the exact functional test.
3. Increased the maximum table size from 5x5 to 10x10 for the exact functional test.
4. Added a warning message when asymptotic test is used in place of the exact test to avoid long computational time.
5. Updated the DOI of the reference to exact functional test.

## Version 2.4.5

2018-02-08

1. Added a parameter "exact.mode.bound" in the fun.chisq.test() function to switch ON/OFF the fast branch-and-bound algorithm in the exact functional test.
2. Updated testthat cases for exact functional test.

## Version 2.4.4

2018-02-05

1. Added a new reference (Zhong and Song, 2018) for exact functional test.
2. Changed the package title to "Chi-Square and Exact Tests for Model-Free Functional Dependency"

2017-10-26

1. Included a reference for the simulate_tables() function that describes the strategies used by the function.

2017-06-11

1. Updated the functional indices to be more accurate when the number of row is less than the number of columns and also when conditioned on the column sum.

2017-05-29

1. Added "simulate.p.value" method option in fun.chisq.test() to calculate p-value with Monte carlo simulated distribution.

2017-05-16

1. Updated the vignette for choosing quantities for
functional dependencies.

2017-05-05

1. Fixed a bug in simulating 'discontinuous' contingency tables when the number of column is two.
2. Updated manuals.

## Version 2.4.3

2017-04-27

1. Fixed a bug in specifying the margin to apply noise in simulate_tables().
2. Added a noise.model parameter in simulate_tables() to specify either the "house" or "candle" noise model for applying noise to contingency tables representing ordinal or categorical data.

2017-04-25

1. Revised package description and manuals.

## Version 2.4.2

2017-04-20

1. Renamed pattern type "nonmonotonic" to "many.to.one" in simulate_tables(). The latter is mathematically correct.
2. Added a "discontinuous" pattern type in simulate_tables().
3. Added a new candle noise model for categorical variables.
4. Sped up add.noise function code by taking advantage of vectorized multinomial distribution R function.

## Version 2.4.1

2017-04-02

1. Revised the manual for function simulate_tables()

2017-02-28

1. Suppressed warning messages when calling chisq.test() to compute chi-squares in function simulate_tables().
2. Edited the manual for simulate_tables().
3. Edited other manuals for improved consistency with R package reference format.

## Version 2.4.0

2017-02-26

1. Added a new R function simulate_tables() with supporting test functions.
2. Added a new R function add.house.noise() with supporting test functions.
3. Introduced the use of 'R_registerRoutines' and 'R_useDynamicSymbols'.

## Version 2.3.4

2016-10-01

1. Added keywords to the package manual.

## Version 2.3.3

2016-09-02

1. Updated the vignettes.

## Version 2.3.2

2016-08-29

1. Updated the reference section to fix the year of Pearson's chi-square test paper which was published in 1900 not 1990.

2016-05-03

1. Now check all options for the method argument so that only valid methods are allowed.

## Version 2.3.1

2016-05-01

1. Expanded test.interactions() to test many-to-one combinatorial interactions in C++ via Rcpp.

## Version 2.3.0 (not deposited to CRAN)

2016-04-30
1. Data frame input is now converted to numeric matrix before exact functional test.
2. Added function test.interactions() to test pairwise (one-to-one) interactions, implemented in C++ via Rcpp for computational efficiency.

## Version 2.2.4

2016-04-21

1. Added vignette "Which statistic to use for functional dependency from fun.chisq.test()?"
2. Updated references and description.

## Version 2.2.3

2016-03-31

1. Use [[Rcpp::export]] to automatically generate R interface.

2016-03-15

1. Revised the code to remove dependency on RcppClassic.

2016-03-12

1. Revised help documentation and included new references.

## Version 2.2.2

2016-02-07

1. Handled a special case for the normalized FunChisq when the degrees of freedom are zero.

2016-02-05

1. Renamed the "type" argument to "alternative" in
fun.chisq.test().
2. Revised values of the "method" argument in fun.chisq.test(). Previous values ("default" and "normalized") are still supported but obsolete.
3. Added the "log.p" argument in both cp.fun.chisq.test() and cp.chisq.test().
4. Updated the documentation.

2016-01-30

1. Added a new argument "index.kind" to fun.chisq.test() to specify the function index kind: "unconditional" or "conditional" on a given marginal of Y.

2016-01-26

1. Fixed a bug in cp.chisq.test().

## Version 2.2.1

2016-01-25

1. Fixed a bug introduced in the previous version.
2. Updated package description to be more reflective of recent additions.
3. Included missing references.

## Version 2.2.0

2016-01-23

1. Improved code efficiency for fun.chisq.test().
2. Now fun.chisq.test() returns an estimate of function index between 0 and 1, in analogy to Cramer's V, but asymmetrical.
3. Added additional test examples for fun.chisq.test().

2016-01-22

1. Added type argument to fun.chisq.test() to specify functional or non-constant functional chi-squares.
2. Added log.p argument to fun.chisq.test() to obtain log of p-value to improve accuracy when sample size is large and p-value is close to zero.

2015-07-03 (Unpublished version 2.1.1 )

1. Added cp.chisq.test() for comparative chi-square test, not considering functional dependencies.
2. Added test examples for cp.chisq.test().

## Version 2.1.0

2015-06-29

1. Substantially reduced the run time of exact functional test by designing a better branch and bound strategy.
2. Enabled C++11 compling by adding two files (Makevars and Makevars.win).
3. Fixed a bug occured under WIN32 using long double by adding define.h to specify double precision for WIN32 and long double precision for WIN64.
4. Increased float comparison precision by considering a tolerance.
5. Imported pnorm() and pchisq() from the "stats" package.

## Version 2.0.2

2015-03-03

1. Updated documentation for the package and its functions.
2. Fixed a testing issue under r-release-linux-ix86 flavor.

## Version 2.0.1

2015-02-23

1. Removed assert() function from StatDistribution.cpp.
2. Removed Makevar and Makevar.win files.
3. Removed all iostream and sstream and cassert.
4. Removed TransitionTableIO.cpp.
5. In the testthat examples, signif(x,8) is used to compare 8 significant digits of the results.
6. Added a normalized functional chi-square test example.

## Version 2.0.0

2015-02-12

1. Added a new exact functional test as a method option to fun.chisq.test(). The exact functional test is an exact version of FunChisq test, it is more precise to detect functional dependencies in small sample-sized contingency tables.
2. Added a new comparative functional chi-square test for detecting heterogeneity in functional dependencies among contingency tables.
3. Revised the examples and documentation to improve usability.
4. Added automated test cases into the package.
5. Started the NEWS file.

## Version 1.0

2014-03-08

1. The first release of this package implements the functional chi-square test and a normalized version.
