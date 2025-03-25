#' Create the 'mob_in' object.
#'
#' The 'mob_in' object will be passed on for analyses of biodiversity across
#' scales.
#'
#' @param comm community matrix in which rows are samples (e.g., plots) and
#'   columns are species.
#' @param plot_attr matrix which includes the environmental attributes and
#'   spatial coordinates of the plots. Environmental attributes are mandatory,
#'   while spatial coordinates are optional.
#' @param coord_names character vector with the names of the columns of
#'   \code{plot_attr} that specify the coordinates of the samples. Defaults to
#'   NULL (no coordinates). When providing coordinate names, the order the names
#'   are provided matters when working with latitude-longitude coordinates
#'   (i.e., argument \code{latlong = TRUE}, and it is expected that the column
#'   specifying the x-coordinate or the longitude is provided first, y-coordinate
#'    or latitude provided second. To provide coordinate names use the following
#'    syntax: \code{coord_names = c('longitude_col_name','latitude_col_name')}
#' @param binary Boolean, defaults to FALSE. Whether the plot by species matrix
#'   "comm" is in abundances or presence/absence.
#' @param latlong Boolean, defaults to FALSE. Whether the coordinates are
#'   latitude-longitudes. If TRUE, distance calculations by downstream functions
#'   are based upon great circle distances rather than Euclidean distances. Note
#'   latitude-longitudes should be in decimal degree.
#'
#' @return a "mob_in" object with four attributes. "comm" is the plot by species
#'   matrix. "env" is the environmental attribute matrix, without the spatial
#'   coordinates. "spat" contains the spatial coordinates (1-D or 2-D). "tests"
#'   specifies whether each of the three tests in the biodiversity analyses is
#'   allowed by data.
#'
#' @author Dan McGlinn and Xiao Xiao
#' @export
#' @examples
#'  data(inv_comm)
#'  data(inv_plot_attr)
#'  inv_mob_in = make_mob_in(inv_comm, inv_plot_attr, coord_names = c('x', 'y'))
#'  
#'  
#'  

#Inbar's code

red_bin_mob_in <- make_mob_in(sp_matric, env_red, # define the sp matrix and the variables   # no wornings apear - v
                              coord_names = c('lon_x', 'lat_y'))   # define coor


#harmonising names
comm <- sp_matric
plot_attr <- env_red
coord_names = c('lon_x', 'lat_y')
binary = FALSE
latlong = FALSE

make_mob_in = function(comm, plot_attr, coord_names = NULL, binary = FALSE,
                       latlong = FALSE) {
  # possibly make group_var and ref_group mandatory arguments
  out = list(tests = list(N = TRUE, SAD = TRUE, agg = TRUE))
  # carry out some basic checks
  if (nrow(comm) < 5) {
    warning("Number of plots in community is less than five therefore only individual rarefaction will be computed")
    out$tests$N = FALSE
    out$tests$agg = FALSE
  }
  
  if (nrow(comm) != nrow(plot_attr)) {
    stop("Number of plots in community does not equal number of plots in plot attribute table")
    
  }
  
  if (is.null(coord_names) == FALSE) {
    spat_cols = sapply(coord_names, function(x) which(x == names(plot_attr)))
    
    if (length(spat_cols) == 1 & latlong == TRUE) {
      stop("Both latitude and longitude have to be specified")
    }
  }
  
  if (any(row.names(comm) != row.names(plot_attr)))
    warning("Row names of community and plot attributes tables do not match
                which may indicate different identities or orderings of samples")
  
  if (binary)  {
    warning("Only spatially-explicit sampled based forms of rarefaction can be computed on binary data")
    out$tests$SAD = FALSE
    out$tests$N = FALSE
  }
  else {
    if (max(comm) == 1)
      warning("Maximum abundance is 1 which suggests data is binary, change the binary argument to TRUE")
  }
  
  if (any(colSums(comm) == 0)) {
    warning("Some species have zero occurrences and will be dropped from the community table")
    comm = comm[ , colSums(comm) != 0]
  }
  
  out$comm = data.frame(comm)
  if (is.null(coord_names) == FALSE) {
    if (length(spat_cols) > 0) {
      out$env = data.frame(plot_attr[ , -spat_cols])
      colnames(out$env) = colnames(plot_attr)[-spat_cols]
      out$spat = data.frame(plot_attr[ , spat_cols])
    }
  }
  else {
    warning("Note: 'coord_names' was not supplied and therefore spatial aggregation will not be examined in downstream analyses")
    out$tests$agg = FALSE
    out$env = data.frame(plot_attr)
    out$spat = NULL
  }
  
  out$latlong = latlong
  class(out) = 'mob_in'
  return(out)
}