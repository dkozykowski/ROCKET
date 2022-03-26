# load required libraries
library(reticulate)

datasets <- import("sktime.datasets")
np <- import("numpy")
utils <- import("sktime.utils.validation.panel")


# define global settings
seed = 5
num_kernels = 50 # increase this number for better accuracy


# set libraries seeds
set.seed(seed)
np$random$seed(as.integer(seed))


# define global settings
seed = 5
num_kernels = 50 # increase this number for better accuracy


# define functions to edit global settings
set_kernel_number <- function(num) {
  num_kernels = num
}

set_seed <- function(new_seed) {
  seed = new_seed
  set.seed(seed)
  np$random$seed(as.integer(seed))
}


# fill missing time series values using approximation
fix_length <- function(x, expected_length) {
  current_length = length(x)
  
  if (current_length != expected_length) {
    x_args = strtoi(stringr::str_trim(names(x)), 10)
    y_args = as.vector(x)
    return(approx(x_args, y_args, 0:(expected_length - 1))$y)
  }
  
  return(x)
}


# standardize time series length using linear approximation for missing values
fill_missing_data <- function(x_train) {
  time_series_count = length(x_train[[1]])
  time_series_dims = length(x_train[1,])
  
  for (i in 1 : time_series_count)
    for (j in 1 : time_series_dims)
      x_train[i,][[j]][[1]] = fix_length(x_train[i,][[j]][[1]], original_length)
  
  return(x_train)
}


# define rocket functions
generate_kernels <- function(X) {
  X = utils$check_X(X, coerce_to_numpy=TRUE)
  
  num_columns = dim(X)[2]
  num_timepoints = dim(X)[3]
  
  lengths = array(as.integer(sample(c(7, 9, 11), num_kernels, replace = TRUE)))
  
  limit = pmin(num_columns, lengths)
  
  num_channel_indices = as.integer(2 ** np$random$uniform(0, log2(limit + 1)))
  
  channel_indices = as.integer(rep(0, sum(num_channel_indices)))
  
  weights = as.double(rep(0, sum(lengths * num_channel_indices)))
  
  biases = array(as.double(rep(0, num_kernels)))
  dilations = array(as.integer(rep(0, num_kernels)))
  paddings = array(as.integer(rep(0, num_kernels)))
  
  a1 = 1  # for weights
  a2 = 1  # for channel_indices
  
  for (i in 1 : num_kernels) {
    temp_length = lengths[i]
    temp_num_channel_indices = num_channel_indices[i]
    
    temp_weights = as.double(np$random$normal(0, 1, temp_num_channel_indices * temp_length))
    
    b1 = a1 + (temp_num_channel_indices * temp_length) - 1 
    b2 = a2 + temp_num_channel_indices - 1
    
    a3 = 1 # for weights (per channel)
    for (j in 1 : temp_num_channel_indices) {
      b3 = a3 + temp_length - 1
      temp_weights[a3 : b3] = temp_weights[a3 : b3] - mean(temp_weights[a3 : b3])
      a3 = b3 + 1
    }
    
    
    weights[a1 : b1] = temp_weights
    
    channel_indices[a2 : b2] = sample(0 : (num_columns - 1), temp_num_channel_indices)
    
    biases[i] = np$random$uniform(-1, 1)
    
    dilation = 2 ** np$random$uniform(0, log2((num_timepoints - 1) / (temp_length - 1)))
    dilation = as.integer(dilation)
    dilations[i] = dilation
    
    if (sample(0 : 1, 1) == 1)
      paddings[i] = as.integer(((temp_length - 1) * dilation) / 2)
    else paddings[i] = 0
    
    a1 = b1 + 1
    a2 = b2 + 1
  } 
  
  channel_indices = channel_indices + 1
  
  return(list(weights,
              lengths,
              biases,
              dilations,
              paddings,
              num_channel_indices,
              channel_indices))
}

apply_kernel_multivariate <- function(X, 
                                      weights, 
                                      length, 
                                      bias, 
                                      dilation, 
                                      padding, 
                                      num_channel_indices, 
                                      channel_indices) {
  
  num_columns = dim(X)[1]
  num_timepoints = dim(X)[2]
  
  output_length = (num_timepoints + (2 * padding)) - ((length - 1) * dilation)
  
  temp_ppv = 0
  temp_max = -Inf
  
  end = (num_timepoints + padding) - ((length - 1) * dilation) - 1
  
  for (i in (-padding + 1) : end) 
  {
    temp_sum = bias
    
    index = i
    
    for (j in 1 : length) {
      if (index > 0 && index <= num_timepoints) {
        for (k in 1 : num_channel_indices) {
          temp_sum = temp_sum + weights[k, j] * X[channel_indices[k], index]
        }
        
      }
      index = index + dilation
    }
    
    if (temp_sum > temp_max)
      temp_max = temp_sum
    
    if (temp_sum > 0)
      temp_ppv = temp_ppv + 1 
  }
  return(c(as.double(temp_ppv / output_length), as.double(temp_max)))
}

apply_kernel_univariate <- function(X, weights, length, bias, dilation, padding) {
  num_timepoints = length(X)
  
  output_length = (num_timepoints + (2 * padding)) - ((length - 1) * dilation)
  
  temp_ppv = 0
  temp_max = -Inf
  
  end = (num_timepoints + padding) - ((length - 1) * dilation)
  
  for (i in (-padding + 1) : end) 
  {
    temp_sum = bias
    
    index = i
    
    for (j in 1 : length) {
      if (index > 0 && index <= num_timepoints)
        temp_sum = temp_sum + weights[j] * X[index]
      
      index = index + dilation
      
    }
    
    if (temp_sum > temp_max)
      temp_max = temp_sum
    
    if (temp_sum > 0)
      temp_ppv = temp_ppv + 1
  }
  return(c(as.double(temp_ppv / output_length), as.double(temp_max)))
}

apply_kernels <- function(x, kernels) {
  
  weights = kernels[[1]]
  lengths = kernels[[2]]
  biases = kernels[[3]]
  dilations = kernels[[4]]
  paddings = kernels[[5]]
  num_channel_indices = kernels[[6]]
  channel_indices = kernels[[7]]
  
  
  x = utils$check_X(x, coerce_to_numpy=TRUE)
  
  num_instances = dim(x)[1]
  num_columns = dim(x)[2]
  
  for(i in 1 : num_instances) {
    for (o in 1 : num_columns) {
      x[i,o,] = (x[i,o,] - mean(x[i,o,])) / sd(x[i,o,]) + 1e-8 
    }
  }
  
  num_kernels = length(lengths)
  
  new_x = matrix(0, num_instances, num_kernels * 2)  # 2 features per kernel
  
  for (i in 1 : num_instances) {
    a1 = 1  # for weights
    a2 = 1  # for channel_indices
    a3 = 1  # for features
    
    for (j in 1 : num_kernels)
    {
      b1 = a1 + num_channel_indices[j] * lengths[j] - 1
      b2 = a2 + num_channel_indices[j] - 1
      b3 = a3 + 2 - 1
      
      if (num_channel_indices[j] == 1)
      {
        
        temp_result = apply_kernel_univariate(
          x[i, channel_indices[a2],],
          weights[a1:b1],
          lengths[j],
          biases[j],
          dilations[j],
          paddings[j]
        ) 
        new_x[i, a3] = temp_result[1]
        new_x[i, a3 + 1] = temp_result[2]
      }
      else
      {
        temp_weights = matrix(weights[a1 : b1], num_channel_indices[j], lengths[j])
        
        temp_result = apply_kernel_multivariate(
          x[i,,],
          temp_weights,
          lengths[j],
          biases[j],
          dilations[j],
          paddings[j],
          num_channel_indices[j],
          channel_indices[a2:b2]
        )
        new_x[i, a3] = temp_result[1]
        new_x[i, a3 + 1] = temp_result[2]
      }
      
      a1 = b1 + 1
      a2 = b2 + 1
      a3 = b3 + 1       
    }
  }
  
  return(as.data.frame(new_x))
}