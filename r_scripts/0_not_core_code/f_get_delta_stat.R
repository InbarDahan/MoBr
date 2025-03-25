get_delta_stats(mob_in = red_bin_mob_in, env_var = 'depth',  # define the env gradient
                group_var = 'depth_group', stats = c('betas', 'r'), # define the site\grouping var
                type = 'continuous', spat_algo = 'kNCN', # define the type of var # define the spatial rarefaction method
                n_perm = 199, overall_p = FALSE) #, density_stat = c("mean", "max", "min")) # 

# mine: 
mob_in = red_bin_mob_in
env_var = 'depth'
group_var = 'depth_group'
ref_level = NULL
tests = c('SAD', 'N', 'agg')
spat_algo = 'kNCN'
type = 'continuous'
stats = c('betas', 'r')
inds = NULL
log_scale = FALSE
min_plots = NULL
density_stat = c('mean', 'max', 'min')
n_perm=199
overall_p = FALSE





#' Conduct the MoB tests on drivers of biodiversity across scales.
#' 
#' There are three tests, on effects of 1. the shape of the SAD, 2.
#' treatment/group-level density, 3. degree of aggregation. The user can
#' specifically to conduct one or more of these tests.
#' 
#' @param mob_in an object of class mob_in created by make_mob_in()
#' @param env_var a character string specifying the environmental variable in
#'   \code{mob_in$env} to be used for explaining the change in richness
#' @param group_var an optional character string 
#'   in \code{mob_in$env} which defines how samples are pooled. If not provided
#'   then each unique value of the argument \code{env_var} is used define the
#'   groups. 
#' @param ref_level a character string used to define the reference level of
#'   \code{env_var} to which all other groups are compared with. Only makes sense
#'   if \code{env_var} is a factor (i.e. \code{type == 'discrete'})
#' @param tests specifies which one or more of the three tests ('SAD', N',
#'   'agg') are to be performed. Default is to include all three tests.
#' @param spat_algo character string that can be either: \code{'kNN'} or
#'   \code{'kNCN'} for k-nearest neighbor and k-nearest centroid neighbor
#'   sampling respectively. It defaults to k-nearest neighbor which is a more
#'   computationally efficient algorithm that closely approximates the
#'   potentially more correct k-NCN algo (see Details of ?rarefaction).
#' @param type "discrete" or "continuous". If "discrete", pair-wise comparisons
#'   are conducted between all other groups and the reference group. If
#'   "continuous", a correlation analysis is conducted between the response
#'   variables and env_var.
#' @param stats a vector of character strings that specifies what statistics to
#'   summarize effect sizes with. Options include: \code{c('betas', 'r2',
#'   'r2adj', 'f', 'p')} for the beta-coefficients, r-squared, adjusted
#'   r-squared, F-statistic, and p-value respectively. The default value of
#'   \code{NULL} will result in only betas being calculated when \code{type ==
#'   'discrete'} and all possible stats being computed when \code{type ==
#'   'continuous'}. Note that for a discrete analysis all non-betas stats are
#'   meaningless because the model has zero degrees of freedom in this context.
#' @param inds effort size at which the individual-based rarefaction curves are
#'   to be evaluated, and to which the sample-based rarefaction curves are to be
#'   interpolated. It can take three types of values, a single integer, a vector
#'   of integers, and NULL. If \code{inds = NULL} (the default), the curves are
#'   evaluated at every possible effort size, from 1 to the total number of
#'   individuals within the group (slow). If inds is a single integer, it is
#'   taken as the number of points at which the curves are evaluated; the
#'   positions of the points are determined by the "log_scale" argument. If inds
#'   is a vector of integers, it is taken as the exact points at which the
#'   curves are evaluated.
#' @param log_scale if "inds" is given a single integer, "log_scale" determines
#'   the position of the points. If log_scale is TRUE, the points are equally
#'   spaced on logarithmic scale. If it is FALSE (default), the points are
#'   equally spaced on arithmetic scale.
#' @param min_plots minimal number of plots for test 'agg', where plots are
#'   randomized within groups as null test. If it is given a value, all groups
#'   with fewer plots than min_plot are removed for this test. If it is NULL
#'   (default), all groups are kept. Warnings are issued if 1. there is only one
#'   group left and "type" is discrete, or 2. there are less than three groups
#'   left and "type" is continuous, or 3. reference group ("ref_group") is
#'   removed and "type" is discrete. In these three scenarios, the function will
#'   terminate. A different warning is issued if any of the remaining groups
#'   have less than five plots (which have less than 120 permutations), but the 
#'   test will be carried out.
#' @param density_stat reference density used in converting number of plots to
#'   numbers of individuals, a step in test "N". It can take one of the
#'   three values: "mean", "max", or "min". If it is "mean", the average
#'   plot-level abundance across plots (all plots when "type" is "continuous,
#'   all plots within the two groups for each pair-wise comparison when "type"
#'   is "discrete") are used. If it is "min" or "max", the minimum/maximum
#'   plot-level density is used.
#' @param n_perm number of iterations to run for null tests, defaults to 1000.
#' @param overall_p Boolean defaults to FALSE specifies if overall across scale 
#'  p-values for the null tests. This should be interpreted with caution because
#'  the overall p-values depend on scales of measurement yet do not explicitly 
#'  reflect significance at any particular scale. 
#' @return a "mob_out" object with attributes
#' @author Dan McGlinn and Xiao Xiao
#' @import dplyr
#' @import purrr
#' @importFrom stats sd
#' @export
#' @seealso	\code{\link{rarefaction}}
#' @examples
#' data(inv_comm)
#' data(inv_plot_attr)
#' inv_mob_in = make_mob_in(inv_comm, inv_plot_attr, coord_names = c('x', 'y'))
#' inv_mob_out = get_delta_stats(inv_mob_in, 'group', ref_level='uninvaded',
#'                            type='discrete', log_scale=TRUE, n_perm=3)
#' plot(inv_mob_out)
get_delta_stats = function(mob_in, env_var, group_var=NULL, ref_level = NULL, 
                           tests = c('SAD', 'N', 'agg'), spat_algo = NULL,
                           type = c('continuous', 'discrete'),
                           stats = NULL, inds = NULL,
                           log_scale = FALSE, min_plots = NULL,
                           density_stat = c('mean', 'max', 'min'),
                           n_perm=1000, overall_p = FALSE) {
  # perform preliminary checks and variable assignments
  if (!methods::is(mob_in, "mob_in"))
    stop('mob_in must be output of function make_mob_in (i.e., of class mob_in')
  if (!(env_var %in% names(mob_in$env)))
    stop(paste(env_var, ' is not one of the columns in mob_in$env.'))
  if (!is.null(group_var))
    if  (!(group_var %in% names(mob_in$env)))
      stop(paste(group_var, ' is not one of the columns in mob_in$env.')) 
  tests = match.arg(tests, several.ok = TRUE)
  test_status = tests %in% names(unlist(mob_in$tests)) 
  approved_tests = tests[test_status]
  if (length(approved_tests) < length(tests)) {
    tests_string = paste(approved_tests, collapse = ' and ')
    warning(paste('Based upon the attributes of the community object only the following tests will be performed:',
                  tests_string))
    tests = approved_tests
  }
  type = match.arg(type)
  density_stat = match.arg(density_stat)
  
  env = mob_in$env[ , env_var]
  # if group_var is NULL then set all samples to same group (??)
  if (is.null(group_var))
    groups = env
  else {
    groups = mob_in$env[ , group_var]
    # check that for the defined groups all samples have same environmental value
    if (any(tapply(env, groups, stats::sd) > 0)) {
      # bc all env values not the same for a group then compute mean value
      message("Computed average environmental value for each group")
      env = tapply(env, groups, mean)
    }
  }    
  if (type == 'discrete') {
    if (!methods::is(env, 'factor')) {
      warning(paste("Converting", env_var, "to a factor with the default contrasts because the argument type = 'discrete'."))
      env = as.factor(env)
    }
    if (!is.null(ref_level)) { # need to ensure that contrasts on the reference level set
      env_levels = levels(env) 
      if (ref_level %in% env_levels) {
        if (env_levels[1] != ref_level)
          env = factor(env, levels = c(ref_level, env_levels[env_levels != ref_level]))
      } else
        stop(paste(ref_level, "is not in", env_var))
    }    
  } else if (type == 'continuous') {
    if (!is.numeric(env)) {
      warning(paste("Converting", env_var, "to numeric because the argument type = 'continuous'"))
      env = as.numeric(as.character(env))
    }
    if (!is.null(ref_level))
      stop('Defining a reference level (i.e., ref_level) only makes sense when doing a discrete analysis (i.e., type = "discrete")')
  }
  #TODO It needs to be clear which beta coefficients apply to 
  # which factor level - this is likely most easily accomplished by appending
  # a variable name to the beta column or adding an additional column
  #
  #if (is.null(env_var)){
  #    env_levels = as.numeric(names(sad_groups))
  #} else {
  #    env_levels = tapply(mob_in$env[, env_var],
  #                        list(groups), mean)
  #}
  N_max = min(tapply(rowSums(mob_in$comm), groups, sum))
  inds = get_inds(N_max, inds, log_scale)
  ind_dens = get_ind_dens(mob_in$comm, density_stat)
  n_plots = min(tapply(mob_in$comm[ , 1], groups, length))
  
  out = list()
  out$env_var = env_var
  if (!is.null(group_var))
    out$group_var = group_var
  out$type = type
  out$tests = tests
  out$log_scale = log_scale
  out$density_stat = list(density_stat = density_stat,
                          ind_dens = ind_dens)
  out = append(out, 
               get_results(mob_in, env, groups, tests, inds, ind_dens, n_plots,
                           type, stats, spat_algo))
  
  null_results = run_null_models(mob_in, env, groups, tests, inds, ind_dens,
                                 n_plots, type, stats, spat_algo,
                                 n_perm, overall_p)
  # merge the null_results into the model data.frame
  out$S_df = left_join(out$S_df, null_results$S_df, 
                       by = c("env", "test", "sample", "effort"))
  out$mod_df = left_join(out$mod_df, null_results$mod_df, 
                         by = c("test", "sample", "effort", "index"))
  if (overall_p)
    out$p = attr(null_results, "p")
  class(out) = 'mob_out'
  return(out)
}
