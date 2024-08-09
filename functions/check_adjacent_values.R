check_adjacent_values <- function(matrix) {
  # Get the dimensions of the matrix
  nrows <- nrow(matrix)
  ncols <- ncol(matrix)
  
  # Initialize the output matrix with zeros
  result <- matrix(0, nrow = nrows, ncol = ncols)
  
  # Loop through each non-edge entry of the matrix
  for (i in 2:(nrows - 1)) {
    for (j in 2:(ncols - 1)) {
      # Get the center value
      center <- matrix[i, j]
      
      # Check the adjacent values
      top <- matrix[i - 1, j]
      bottom <- matrix[i + 1, j]
      left <- matrix[i, j - 1]
      right <- matrix[i, j + 1]
      
      # If all adjacent values are the same as the center, set result to 1
      if (center == top && center == bottom && center == left && center == right) {
        result[i, j] <- 0
      } else {
        result[i, j] <- 1
      }
    }
  }
  
  return(result)
}
