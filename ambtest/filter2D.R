filter2D <- function(img, filter) {
  # Algoritm for 2D Convolution
  filter_center_index_y <- median(1:dim(filter)[1])
  filter_max_index_y <- dim(filter)[1]
  filter_center_index_x <- median(1:dim(filter)[2])
  filter_max_index_x <- dim(filter)[2]
  
  # For each position in the picture, 2D convolution is done by 
  # calculating a score for all overlapping values within the two matrices
  x_min <- 1
  x_max <- dim(img)[2]
  y_min <- 1
  y_max <- dim(img)[1]
  
  df <- NULL
  for (x_val in c(x_min:x_max)){
    for (y_val in c(y_min:y_max)){
      # Distanced from cell
      img_dist_left <- x_val-1
      img_dist_right <- x_max-x_val
      img_dist_up <- y_val-1
      img_dist_down <- y_max-y_val
      
      # Overlapping filter cells
      filter_x_start <- filter_center_index_x-img_dist_left
      if (filter_x_start < 1) {
        filter_x_start <- 1
      }
      filter_x_end <- filter_center_index_x+img_dist_right
      if (filter_x_end > filter_max_index_x) {
        filter_x_end <- filter_max_index_x
      }
      filter_y_start <- filter_center_index_y-img_dist_up
      if (filter_y_start < 1) {
        filter_y_start <- 1
      }
      filter_y_end <- filter_center_index_y+img_dist_down
      if (filter_y_end > filter_max_index_y) {
        filter_y_end <- filter_max_index_y
      }
      
      # Part of filter that overlaps
      filter_overlap_matrix <- filter[filter_y_start:filter_y_end, filter_x_start:filter_x_end]
      
      # Overlapped image cells
      image_x_start <- x_val-filter_center_index_x+1
      if (image_x_start < 1) {
        image_x_start <- 1
      }
      image_x_end <- x_val+filter_max_index_x-filter_center_index_x
      if (image_x_end > x_max) {
        image_x_end <- x_max
      }
      image_y_start <- y_val-filter_center_index_y+1
      if (image_y_start < 1) {
        image_y_start <- 1
      }
      image_y_end <- y_val+filter_max_index_y-filter_center_index_y
      if (image_y_end > y_max) {
        image_y_end <- y_max
      }
      
      # Part of image that is overlapped
      image_overlap_matrix <- img[image_y_start:image_y_end, image_x_start:image_x_end]
      
      # Calculating the cell value
      cell_value <- sum(filter_overlap_matrix*image_overlap_matrix)
      df = rbind(df,data.frame(x_val,y_val, cell_value))
    }
  }
  
  # Axis labels
  x_axis <- c(x_min:x_max)
  y_axis <- c(y_min:y_max)
  
  # Populating matrix
  filter_matrix <- matrix(df[,3], nrow = x_max, ncol = y_max, dimnames = list(x_axis, y_axis))
  
  # New - Determine rows and columns of matrix as well as the filter kernel
  nrow_filter <- nrow(filter)
  ncol_filter <- ncol(filter)
  nrows <- nrow(filter_matrix)
  ncols <- ncol(filter_matrix)
  
  # New - Figure out where to cut off
  row_cutoff <- floor(nrow_filter/2)
  col_cutoff <- floor(ncol_filter/2)
  
  # New - Remove out borders
  filter_matrix <- filter_matrix[((1+row_cutoff):(nrows-row_cutoff)), ((1+col_cutoff):(ncols-col_cutoff))]
  
  # Finally return matrix
  return(filter_matrix)   
}
