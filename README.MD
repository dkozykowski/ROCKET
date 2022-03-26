# ROCKET Algorithm implementation
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Implementation of the *RandOm Convolutional KErnel Transform* [ROCKET] published in [arXiv:1910.13051](https://arxiv.org/pdf/1910.13051.pdf) in R language.

This implementation works for all multi-dimensional time series. It also works with time series with partially missing data, since data approximation has been implemented. 

This project is a part of Warsaw University of Technology Machine Learning Course.

## Overview

ROCKET is an exceptionally fast and accurate time series classification algorithm. 
The algorithm randomly generates a great variety of convolutional kernels and extracts two features for each convolution: the maximum and the proportion of positive values.

- ROCKET doesn’t use a hidden layer or any non-linearities
- Features produced by ROCKET are independent of each other
- ROCKET works with any kind of classifier
- ROCKET uses a very large number of kernels
- In CNN, a group of kernels tend to share the same size, dilation and padding. ROCKET has all 5 parameters randomized.
- In CNN, Dilation increases exponentially with depth; ROCKET has random dilation values
- CNNs only have average/max pooling. ROCKET has a unique pooling called as ppv which has proven to provide much better classification accuracy on time series.

## Important notes

 Time series used as example are provided by [timeseriesclassification.com](https://timeseriesclassification.com/dataset.php).
