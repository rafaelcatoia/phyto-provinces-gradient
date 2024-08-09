# Load necessary library
#library(dplyr)

# Define the function
stack_dataframes <- function(df_list) {
  # Check if the input is a list
  # if (!is.list(df_list)) {
  #   stop("The input must be a list.")
  # }
  
  # Initialize an empty list to store modified dataframes
  modified_dfs <- list()
  
  # Loop through each dataframe in the list
  for (i in seq_along(df_list)) {
    # Check if the current element is a dataframe
    # if (!is.data.frame(df_list[[i]])) {
    #   stop("All elements of the list must be dataframes.")
    # }
    
    # Determine the index name
    index_name <- if (!is.null(names(df_list))) names(df_list)[i] else as.character(i)
    
    # Add a new column with the name of the list index
    df_list[[i]]$index_name <- index_name
    
    # Append the modified dataframe to the list
    modified_dfs[[i]] <- df_list[[i]]
  }
  
  # Combine all modified dataframes into one
  combined_df <- dplyr::bind_rows(modified_dfs)
  
  return(combined_df)
}

# # Example usage
# df1 <- data.frame(a = 1:3, b = 4:6)
# df2 <- data.frame(a = 7:9, b = 10:12)
# df_list <- list(df1 = df1, df2 = df2)
# 
# stacked_df <- stack_dataframes(df_list)
# print(stacked_df)
# class(stacked_df)
# 
# # Example with unnamed list
# df_list2 <- list(df1, df2)
# stacked_df2 <- stack_dataframes(df_list2)
# print(stacked_df2)# 