# r_function_of_ncmm_and_ncrc
The R function for NC-MM and NC-RC methods for causal effect estimation based on proxy variable by controlling unmeasured confounders.

We provides three functions: NCreplicate(), NCvalidation(), and NCcorrected(), where the third function is built on the basis of the first two functions. 

When the basic assumptions of methods D-1, R-1, and C-1 hold, the NCcorrected() function, based on the proposed NC-MM and NC-RC methods, can be used respectively for the detection, reduction, or correction of unmeasured confounding when one or more of exposure, negative control exposure, and observable confounding factors are measured with traditional additive non-differential measurement errors. The function and related parameters are described as follows:

NCcorrected(mismeasuredvariable, precisievarable, dataused, datavalidate, datatype, aim, boottime)

mismeasuredvariable: A vector composed of variables for which only observed values with measurement errors exist, such as mismeasuredvariable=c("X", "C").

precisievatable: A vector composed of variables that can be directly and accurately measured, with the outcome variable "Y" included, such as precisievatable = c("Y", "M").

datatype: Specifies whether to estimate variance based on replicated or internally validated data sets, and can take values of "replicate" or "validation". If datatype = replicate, the main data set needs to be provided through the dataused parameter, for example, dataused = dataA, where dataA is a dataframe satisfying:

(a) The aliases for negative control exposure, exposure, and outcome variables should be named as "C", "X", and "Y" respectively, while observable confounding factors should be named as "M" or other single uppercase letters.

(b) The repeated measurements of variables with measurement errors need to be distinguished by numbers such as "1, 2, 3, ...", and Arabic numerals should not appear in the names of variables that can be accurately measured.

(c) The order of variables without accurate measurements in dataA should match the order of variables in the mismeasuredvariable parameter, such as dataA = data.frame(X1, X2, C1, C2, M, Y).

(d) datavalidate = NA.

If datatype = "validation", the main data set needs to be provided through the dataused parameter, such as dataused = dataB, where the measurements of variables that cannot be accurately obtained are named with the number "1". At the same time, the internally validated data set needs to be provided through the datavalidate parameter, such as datavalidate = dataC, where the observed variables with measurement errors are named with the number "2", for example: dataB = data.frame(X1, C1, Y, M) and dataC = data.frame(vX, vC, X2, C2). It is important to note that the order of variable names in the "mismeasuredvariable" and "precisievarable" parameters must match the order of variables in dataB and dataC.

aim: Specifies the purpose of the study, which can be "correction", "detection", or "reduction", corresponding to methods D-1, R-1, and C-1, respectively.

boottime: Specifies the number of repetitions for Bootstrap when estimating variances using NC-RC and NC-MM.

The output of the function includes the effect estimates, standard errors, and 95\% confidence intervals for the NC-MM and NC-RC methods.
