most_frequent_value <- function(numeric_matrix) {
  # Initialize a vector to store the results
  result <- vector("numeric", nrow(numeric_matrix))
  
  # Loop through each row of the numeric_matrix
  for (i in 1:nrow(numeric_matrix)) {
    # Get the current row
    row <- numeric_matrix[i, ]
    
    # Create a frequency table for the current row
    freq_table <- table(row)
    
    # Find the maximum frequency
    max_freq <- max(freq_table)
    
    # Get the values with the maximum frequency
    max_freq_values <- as.numeric(names(freq_table[freq_table == max_freq]))
    
    # If there's a tie, choose the smallest column index
    if (length(max_freq_values) > 1) {
      for (j in 1:ncol(numeric_matrix)) {
        if (row[j] %in% max_freq_values) {
          result[i] <- row[j]
          break
        }
      }
    } else {
      result[i] <- max_freq_values
    }
  }
  
  return(result)
}

## Matrix to test
#test_matrix <- matrix(
# c(c(3,3,1,1,2,2),
#   c(3,1,1,3,1,3),
#   c(6,5,4,3,2,1),
#   c(0,1,3,3,2,2),
#   c(3,3,1,2,2,2),
#   c(1,1,1,3,3,3)),byrow = T,ncol = 6)
#
#most_frequent_value(test_matrix)
