# CHANGES IN enpls VERSION 5.6 (2016-11-25)

## NEW FEATURES

* New argument `cvfolds` now available in all applicable functions for finer control of cross-validation folds in automatic parameter selection of each PLS/SPLS model.

## IMPROVEMENTS

* Critical implementation improvements for processing the case where argument `maxcomp = NULL` (maximum number of components not specified explicitly) in `enpls.` functions. Now it should correctly determine the maximum number of components to use, considering both cross-validation and special cases such as n < p. Thanks to Dr. You-Wu Lin for the feedback.

# CHANGES IN enpls VERSION 5.0 (2016-10-20)

## NEW FEATURES

* New functions `enpls.ad()` and `enspls.ad()` for PLS and sparse PLS model applicability domain evaluation.
* New plot functions `plot.enpls.ad()` and `plot.enpls.ad()` for exploring model applicability domain evaluation results with traditional static plot support and interactive plot support.
* New argument `alpha` available for setting transparency level (to reduce overplotting) in `plot.cv.enpls`, `plot.enpls.od`, `plot.cv.enspls`, and `plot.enspls.od`.

# CHANGES IN enpls VERSION 4.5 (2016-09-15)

## IMPROVEMENTS

* Reduced memory footprints for `enpls.fit()` and `enspls.fit()`.

## NEW FEATURES

* New functions `enpls.rmse()`, `enpls.mae()`, and `enpls.rmsle()` for computing RMSE, MAE, and RMSLE.

# CHANGES IN enpls VERSION 4.0 (2016-08-28)

## IMPROVEMENTS

* General improvements on function documentation.
* Changing the argument name `MCtimes` to `reptimes`.
* Changing the option name `"bootstrap"` to `"boot"`.

# CHANGES IN enpls VERSION 3.0 (2016-06-22)

## NEW FEATURES

* Add sparse partial least squares regression.

## IMPROVEMENTS

* Improvements on documentation; rewritten vignette.

# CHANGES IN enpls VERSION 2.0 (2016-06-19)

## NEW FEATURES

* New project website: [http://enpls.org](http://enpls.org)

## IMPROVEMENTS

* General improvements on plotting functions.

## BUG FIXES

* Fixed major bugs in `cv.enpls` and plotting functions.
* Many other bug fixes.

# CHANGES IN enpls VERSION 1.1 (2015-11-26)

##  BUG FIXES

* Fixed the bugs in automatic component number selection which could make
the intercept-only model rank best. Thanks for the test and patch from
Max Kuhn <<max.kuhn@pfizer.com>>.

# CHANGES IN enpls VERSION 1.0 (2014-10-03)

## NEW FEATURES

* initial release
