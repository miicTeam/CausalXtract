#' MIIC, causal network learning algorithm including latent variables
#'
#' @description MIIC (Multivariate Information based Inductive Causation) combines
#' constraint-based and information-theoretic approaches to disentangle direct
#' from indirect effects amongst correlated variables, including cause-effect
#' relationships and the effect of unobserved latent causes.
#'
#' @details In regular mode, starting from a complete graph, the method iteratively removes
#' dispensable edges, by uncovering significant information contributions from
#' indirect paths, and assesses edge-specific confidences from randomization of
#' available data. The remaining edges are then oriented based on the signature
#' of causality in observational data.
#'
#' In temporal mode (when any \emph{tau} >= 1), miic reorganizes the dataset
#' using the \emph{tau} and \emph{delta_tau} parameters to transform the timesteps
#' into lagged samples. As starting point, a lagged graph is created with
#' only edges having at least one node laying on the last timestep.
#' Then, miic standard algorithm is applied to remove dispensable edges.
#' The remaining edges are then oriented by using the temporality and the
#' signature of causality in observational data.
#'
#' The method relies on an information theoretic based (conditional) independence
#' test which is described in (Verny \emph{et al.}, PLoS Comp. Bio. 2017),
#' (Cabeli \emph{et al.}, PLoS Comp. Bio. 2020). It deals with both categorical
#' and continuous variables by performing optimal context-dependent discretization.
#' As such, the input data frame may contain both numerical columns which will be
#' treated as continuous, or character / factor columns which will be treated
#' as categorical. For further details on the optimal discretization method and
#' the conditional independence test, see the function discretizeMutual.
#' The user may also choose to run miic with scheme presented in
#' (Li \emph{et al.}, NeurIPS 2019) to improve the end result's interpretability
#' by ensuring consistent separating set during the skeleton iterations.
#'
#' @seealso \code{\link{discretizeMutual}} for optimal discretization and
#' (conditional) independence test.
#'
#' @references
#' \itemize{
#' \item Verny et al., \emph{PLoS Comp. Bio. 2017.}  https://doi.org/10.1371/journal.pcbi.1005662
#' \item Cabeli et al., \emph{PLoS Comp. Bio. 2020.}  https://doi.org/10.1371/journal.pcbi.1007866
#' \item Li et al., \emph{NeurIPS 2019} http://papers.nips.cc/paper/9573-constraint-based-causal-structure-learning-with-consistent-separating-sets.pdf
#' }
#'
#' @param input_data [a data frame]
#' A n*d data frame (n samples, d variables) that contains the observational data.
#'
#' In regular mode, each column corresponds to one variable and each row is a sample that gives the
#' values for all the observed variables. The column names correspond to the
#' names of the observed variables. Numeric columns will be treated as continuous
#' values, factors and character as categorical.
#'
#' In temporal mode (when any \emph{tau} >= 1), the expected dataframe layout
#' is variables as columns and timeseries/timesteps as rows.
#' The timestep information must be supplied in the first column and,
#' for each timeseries, be consecutive (increment of 1) and in ascending order.
#' Multiple trajectories can be provided, miic will consider that a new trajectory
#' starts each time a smaller timestep than the one of the previous row is encountered.
#'
#' @param black_box [a data frame]
#' An optional E*2 data frame containing E pairs of variables that will be considered
#' as independent during the network reconstruction. In practice, these edges will not
#' be included in the skeleton initialization and cannot be part of the final result.
#' Variable names must correspond to the \emph{input_data} data frame.
#'
#' @param n_eff [a positive integer]
#' In regular mode, the n samples given in the \emph{input_data} data frame are
#' expected to be independent. In case of correlated samples such as in
#' Monte Carlo sampling approaches, the effective number of independent samples
#' \emph{n_eff} can be estimated using the decay of the autocorrelation function
#' (Verny \emph{et al.}, PLoS Comp. Bio. 2017). This \emph{effective} number \emph{n_eff}
#' of \emph{independent} samples can be provided using this parameter.
#'
#' @param cplx [a string; \emph{c("nml", "mdl")}]
#' In practice, the finite size of the input
#' dataset requires that the 2-point and 3-point information measures should be
#' \emph{shifted} by a \emph{complexity} term. The finite size corrections can be
#' based on the Minimal Description Length (MDL) criterion (set the option with "mdl").
#' In practice, the MDL complexity criterion tends to underestimate the relevance of
#' edges connecting variables with many different categories, leading to the removal of
#' false negative edges. To avoid such biases with finite datasets, the (universal)
#' Normalized Maximum Likelihood (NML) criterion can be used (set the option with "nml").
#' The default is "nml" (see Affeldt \emph{et al.}, UAI 2015).
#'
#' @param latent [a string; \emph{c("no", "yes", "orientation")}]
#' When set to "yes", the network reconstruction is taking into account hidden (latent)
#' variables. When set to "orientation", latent variables are not considered during the skeleton
#' reconstruction but allows bi-directed edges during the orientation. Dependence
#' between two observed variables due to a latent variable is indicated with a '6' in
#' the adjacency matrix and in the network edges.summary and by a bi-directed edge
#' in the (partially) oriented graph.
#'
#' @param orientation [a boolean value]
#' The miic network skeleton can be partially directed
#' by orienting and propagating edge directions, based on the sign and magnitude
#' of the conditional 3-point information of unshielded triples. The propagation
#' procedure relyes on probabilities; for more details, see Verny \emph{et al.}, PLoS Comp. Bio. 2017).
#' If set to FALSE the orientation step is not performed.
#'
#' @param ori_proba_ratio [a floating point between 0 and 1] The threshold when
#' deducing the type of an edge tip (head/tail) from the probability of
#' orientation. For a given edge tip, denote by p the probability of it being a
#' head, the orientation is accepted if (1 - p) / p < ori_proba_ratio. 0 means
#' reject all orientations, 1 means accept all orientations.
#'
#' @param ori_consensus_ratio [a floating point between 0 and 1] The threshold
#' when deducing the type of an consensus edge tip (head/tail) from the average
#' probability of orientation. For a given edge tip, denote by p the probability
#' of it being a head, the orientation is accepted if (1 - p) / p <
#' ori_consensus_ratio. 0 means reject all orientations, 1 means accept all
#' orientations.
#'
#' @param propagation [a boolean value]
#' If set to FALSE, the skeleton is partially oriented with only the
#' v-structure orientations. Otherwise, the v-structure orientations are
#' propagated to downstream undirected edges in unshielded triples following
#' the orientation method
#'
#' @param state_order [a data frame] An optional data frame providing extra
#' information for variables. It must have d rows where d is the number of input
#' variables, and the following structure (named columns):
#'
#' "var_names" (required) contains the name of each variable as specified
#' by colnames(input_data).
#'
#' "var_type" (optional) contains a binary value that specifies if each
#' variable is to be considered as discrete (0) or continuous (1).
#'
#' "levels_increasing_order" (optional) contains a single character string
#' with all of the unique levels of the ordinal variable in increasing order,
#' delimited by comma ','. It will be used during the post-processing to compute
#' the sign of an edge using Spearman's rank correlation. If a variable is
#' continuous or is categorical but not ordinal, this column should be NA.
#'
#' "is_contextual" (optional) contains a binary value that specifies if a
#' variable is to be considered as a contextual variable (1) or not (0).
#' Contextual variables cannot be the child node of any other variable (cannot
#' have edge with arrowhead pointing to them).
#'
#' "tau" (optional) contains an integer value that specifies the number of
#' layers to be considered for the variable.\cr
#' Note that if "tau" column is present in the \emph{state_order} and valid
#' (at least one tau is >= 1), it will induce the switch from regular to
#' temporal mode.\cr
#' Note also that the tau values presents in the state_oder will overwrite the
#' function parameter.
#'
#' "delta_tau" (optional) contains an integer value that specifies the number
#' of timesteps between each layer for the variable.
#' Note that if "delta_tau" column is present in the  \emph{state_order},
#' its values will overwrite the function parameter.
#'
#' "movavg" (optional) contains an integer value that specifies the size of
#' the moving average window to be applied to the variable.
#' Note that if "movavg" column is present in the \emph{state_order},
#' its values will overwrite the function parameter.
#'
#' @param true_edges [a data frame]
#' An optional E*2 data frame containing the E edges of the true graph for
#' computing performance after the run.
#'
#' @param n_shuffles [a positive integer] The number of shufflings of
#' the original dataset in order to evaluate the edge specific confidence
#' ratio of all inferred edges.
#'
#' @param conf_threshold [a positive floating point] The threshold used
#' to filter the less probable edges following the skeleton step. See Verny
#' \emph{et al.}, PLoS Comp. Bio. 2017.
#'
#' @param sample_weights [a numeric vector]
#' An optional vector containing the weight of each observation.
#'
#' @param test_mar [a boolean value]
#' If set to TRUE, distributions with missing values will be tested with Kullback-Leibler
#' divergence : conditioning variables for the given link \eqn{X\rightarrow Y}\eqn{Z} will be
#' considered only if the divergence between the full distribution and the non-missing
#' distribution \eqn{KL(P(X,Y) | P(X,Y)_{!NA})} is low enough (with \eqn{P(X,Y)_{!NA}} as
#' the joint distribution of \eqn{X} and \eqn{Y} on samples which are not missing on Z.
#' This is a way to ensure that data are missing at random for the considered
#' interaction and to avoid selection bias. Set to TRUE by default
#'
#' @param consistent [a string; \emph{c("no", "orientation", "skeleton")}]
#' if "orientation": iterate over skeleton and orientation steps to ensure
#' consistency of the network;
#' if "skeleton": iterate over skeleton step to get a consistent skeleton, then
#' orient edges and discard inconsistent orientations to ensure consistency of
#' the network. See (Li \emph{et al.}, NeurIPS 2019) for details.
#'
#' @param max_iteration [a positive integer] When the \emph{consistent} parameter
#' is set to "skeleton" or "orientation", the maximum number of iterations
#' allowed when trying to find a consistent graph. Set to 100 by default.
#'
#' @param consensus_threshold [a floating point between 0.5 and 1.0] When the
#' \emph{consistent} parameter is set to "skeleton" or "orientation", and when
#' the result graph is inconsistent, or is a union of more than one inconsistent
#' graphs, a consensus graph will be produced based on a pool of graphs. If the
#' result graph is inconsistent, then the pool is made of [max_iteration] graphs
#' from the iterations, otherwise it is made of those graphs in the union. In
#' the consensus graph, an edge is present when the proportion of non-zero
#' status in the pool is above the threshold. For example, if the pool contains
#' [A, B, B, 0, 0], where "A", "B" are different status of the edge and "0"
#' indicates the absence of the edge. Then the edge is set to connected ("1") if
#' the proportion of non-zero status (0.6 in the example) is equal to or higher
#' than [consensus_threshold]. (When set to connected, the orientation of the
#' edge will be further determined by the average probability of orientation.)
#' Set to 0.8 by default.
#'
#' @param tau [an int > 0 or a list of int] Optional, -1 by default.\cr
#' Tau(s) define(s) the number of layers that will be considered for the
#' variables. The layers will be distant of \emph{delta_tau} timesteps.\cr
#' \itemize{
#' \item When an integer is supplied, it must be strictly positive and all the
#'       variables will be lagged using this value as the number of layers.
#' \item When a list is supplied, the size of the list must match the number of
#'       variables and each variable (so excluding the timesteps column)
#'       will be lagged by the corresponding number of layers.
#'       Variables associated with a \emph{tau} values < 1  will not be lagged.
#' }
#' If supplied and valid (at least one tau is >= 1), the \emph{tau} parameter
#' induces the switch from regular to temporal mode.
#'
#' @param delta_tau [an integer or a list of int] Optional, 1 by default.\cr
#' Used only in temporal mode. the delta_tau(s) define(s) the number of timesteps
#' between each layer.\cr
#' i.e.: on 1000 timesteps with  \emph{tau} = 2 and \emph{delta_tau} = 7,
#' the timesteps kept for the samples conversion will be 1, 8, 15
#' for the first sample, the next sample will use 2, 9, 16 and so on.\cr
#' \itemize{
#' \item When an integer is supplied (integer > 0). All the variables will be
#'       lagged using this value as the number of timesteps between layers.
#' \item When a list is supplied, the size of the list must match the number of
#'       variables (so excluding the timesteps column) and each variable will
#'       be lagged with corresponding \emph{delta_tau} timesteps between layers.
#'       If a variable is not lagged ( \emph{tau} < 1), the \emph{delta_tau}
#'       is ignored.
#' }
#'
#' @param movavg [an integer or a list of int] Optional, -1 by default.\cr
#' Used only in temporal mode.\cr
#' \itemize{
#' \item When an integer is supplied (integer > 1), a moving average
#'       operation is applied to all integer and numeric variables.
#' \item When a list is supplied, the size of the list must match the number of
#'       variables. A moving average will be applied to each variable having
#'       a movavg value greater than 1
#' }
#'
#' @param keep_max_data [a boolean value] Optional, FALSE by default.\cr
#' Used only in temporal mode. If TRUE, rows where some NAs have been
#' introduced during the lagging will be kept whilst they will be dropped
#' if FALSE.
#'
#' @param verbose [a boolean value] If TRUE, debugging output is printed.
#'
#' @param n_threads [a positive integer]
#' When set greater than 1, n_threads parallel threads will be used for computation. Make sure
#' your compiler is compatible with openmp if you wish to use multithreading.
#'
#' @param negative_info [a boolean value] For test purpose only. FALSE by
#' default. If TRUE, negative shifted mutual information is allowed during the
#' computation when mutual information is inferior to the complexity term. For
#' small dateset with complicated structures, e.g., discrete variables with many
#' levels, allowing for negative shifted mutual information may help identifying
#' weak v-structures related to those discrete variables, as the negative
#' three-point information in those cases will come from the difference between
#' two negative shifted mutual information terms (expected to be negative due to
#' the small sample size). However, under this setting, a v-structure (X -> Z <-
#' Y) in the final graph does not necessarily imply that X is dependent on Y
#' conditioning on Z, As a consequence, the interpretability of the final graph
#' is hindered. In practice, it's advised to keep this parameter as FALSE.
#'
#' @return A \emph{miic-like} object that contains:
#' \itemize{
#'  \item{all.edges.summary:}{ a data frame with information about the relationship between
#'  each pair of variables
#'  \itemize{
#'  \item \emph{x:} X node
#'  \item \emph{y:} Y node
#'  \item \emph{type:} contains 'N' if the edge has
#'  been removed or 'P' for retained edges. If a true edges file is given,
#'  'P' becomes 'TP' (True Positive) or 'FP' (False Positive), while
#'  'N' becomes 'TN' (True Negative) or 'FN' (False Negative).
#'  \item \emph{ai:} the contributing nodes found by the method which participate in
#'  the mutual information between \emph{x} and \emph{y}, and possibly separate them.
#'  \item \emph{info:} provides the pairwise mutual information times \emph{Nxyi} for
#'  the pair (\emph{x}, \emph{y}).
#'  \item \emph{info_cond:} provides the conditional mutual information times \emph{Nxy_ai} for
#'  the pair (\emph{x}, \emph{y}) when conditioned on the collected nodes \emph{ai}. It is
#'  equal to the \emph{info} column when \emph{ai} is an empty set.
#'  \item \emph{cplx:} gives the computed complexity between the (\emph{x}, \emph{y})
#'  variables taking into account the contributing nodes \emph{ai}. Edges that have
#'  have more conditional information \emph{info_cond} than \emph{cplx} are retained in the
#'  final graph.
#'  \item \emph{Nxy_ai:} gives the number of complete samples on which the information and
#'  the  complexity have been computed. If the input dataset has no missing value, the
#'  number of samples is the same for all pairs and corresponds to the total
#'  number of samples.
#'  \item \emph{info_shifted:} represents the \emph{info} - \emph{cplx} value.
#'  It is a way to quantify the strength of the edge (\emph{x}, \emph{y}).
#'  \item \emph{confidenceRatio:} this column is present if the confidence cut
#'  is > 0 and it represents the ratio between the probability to reject
#'  the edge (\emph{x}, \emph{y}) in the dataset versus the mean probability
#'  to do the same in multiple (user defined) number of randomized datasets.
#'  \item \emph{infOrt:} the orientation of the edge (\emph{x}, \emph{y}). It is
#'  the same value as in the adjacency matrix at row \emph{x} and column \emph{y} : 1 for
#'  unoriented, 2 for an edge from X to Y, -2 from Y to X and 6 for bidirectional.
#'  \item \emph{trueOrt:} the orientation of the edge (\emph{x}, \emph{y}) present
#'  in the true edges file if provided.
#'  \item \emph{isOrtOk:} information about the consistency of the inferred graph’s
#'  orientations with a reference graph is given (i.e. if true edges file is provided).
#'  Y: the orientation is consistent; N: the orientation is not consistent with
#'  the PAG (Partial Ancestor Graph) derived from the given true graph.
#'  \item \emph{sign:} the sign of the partial correlation between variables
#'  \emph{x} and \emph{y}, conditioned on the contributing nodes \emph{ai}.
#'  \item \emph{partial_correlation:} value of the partial correlation for the
#'  edge (\emph{x}, \emph{y}) conditioned on the contributing nodes \emph{ai}.
#'  \item \emph{isCausal:} details about the nature of the arrow tip for a directed
#'  edge. A directed edge in a causal graph does not necessarily imply causation but it
#'  does imply that the cause-effect relationship is not the other way around. An arrow-tip
#'  which is itself downstream of another directed edge suggests stronger causal sense and is
#'  marked by a 'Y', or 'N' otherwise.
#'  \item \emph{proba:} probabilities for the inferred orientation, derived from the three-point
#'  mutual information (cf Affeldt & Isambert, UAI 2015 proceedings) and noted as p(x->y);p(x<-y).
#'  }
#'  }
#'
#'  \item{retained.edges.summary:} {a data frame in the format of all.edges.summary containing only the inferred edges.}
#'
#'  \item{orientations.prob:} {this data frame lists the orientation probabilities of the two edges of all unshielded triples
#'  of the reconstructed network with the structure: node1 -- mid-node -- node2:}
#'  \itemize{
#'  \item node1: node at the end of the unshielded triplet
#'  \item p1: probability of the arrowhead node1 <- mid-node
#'  \item p2: probability of the arrowhead node1 -> mid-node
#'  \item mid-node: node at the center of the unshielded triplet
#'  \item p3: probability of the arrowhead mid-node <- node2
#'  \item p4: probability of the arrowhead mid-node -> node2
#'  \item node2: node at the end of the unshielded triplet
#'  \item NI3: 3 point (conditional) mutual information * N
#'  }
#'
#'
#'  \item {AdjMatrix:} the adjacency matrix is a square matrix used to represent
#'  the inferred graph. The entries of the matrix indicate whether pairs of
#'  vertices are adjacent or not in the graph. The matrix can be read as a
#'  (row, column) set of couples where the row represents the source node and
#'  the column the target node. Since miic can reconstruct mixed networks
#'  (including directed, undirected and bidirected edges), we will have a
#'  different digit for each case:
#'  \itemize{
#'  \item 1: (\emph{x}, \emph{y}) edge is undirected
#'  \item 2: (\emph{x}, \emph{y}) edge is directed as \emph{x} -> \emph{y}
#'  \item -2: (\emph{x}, \emph{y}) edge is directed as \emph{x} <- \emph{y}
#'  \item 6: (\emph{x}, \emph{y}) edge is bidirected
#'  }
#' }
#' @export
#' @useDynLib miic
#' @import Rcpp
#'
#' @examples
#' library(miic)
#'
#' # EXAMPLE HEMATOPOIESIS
#' data(hematoData)
#'
#' # execute MIIC (reconstruct graph)
#' miic.res <- miic(
#'   input_data = hematoData[1:1000,], latent = "yes",
#'   n_shuffles = 10, conf_threshold = 0.001
#' )
#'
#' # plot graph
#' if(require(igraph)) {
#'  plot(miic.res, method="igraph")
#' }
#'
#' \donttest{
#' # write graph to graphml format. Note that to correctly visualize
#' # the network we created the miic style for Cytoscape (http://www.cytoscape.org/).
#'
#' miic.write.network.cytoscape(g = miic.res, file = file.path(tempdir(), "temp"))
#'
#' # EXAMPLE CANCER
#' data(cosmicCancer)
#' data(cosmicCancer_stateOrder)
#' # execute MIIC (reconstruct graph)
#' miic.res <- miic(
#'   input_data = cosmicCancer, state_order = cosmicCancer_stateOrder, latent = "yes",
#'   n_shuffles = 100, conf_threshold = 0.001
#' )
#'
#' # plot graph
#' if(require(igraph)) {
#'  plot(miic.res)
#' }
#'
#' # write graph to graphml format. Note that to correctly visualize
#' # the network we created the miic style for Cytoscape (http://www.cytoscape.org/).
#' miic.write.network.cytoscape(g = miic.res, file = file.path(tempdir(), "temp"))
#'
#' # EXAMPLE OHNOLOGS
#' data(ohno)
#' data(ohno_stateOrder)
#' # execute MIIC (reconstruct graph)
#' miic.res <- miic(
#'   input_data = ohno, latent = "yes", state_order = ohno_stateOrder,
#'   n_shuffles = 100, conf_threshold = 0.001
#' )
#'
#' # plot graph
#' if(require(igraph)) {
#'  plot(miic.res)
#' }
#'
#' # write graph to graphml format. Note that to correctly visualize
#' # the network we created the miic style for Cytoscape (http://www.cytoscape.org/).
#' miic.write.network.cytoscape(g = miic.res, file = file.path(tempdir(), "temp"))
#'
#' # EXAMPLE COVID CASES (timeseries demo)
#' data(covidCases)
#' # execute MIIC (reconstruct graph in temporal mode)
#' tmiic.res <- miic(input_data = covidCases, tau = 2, movavg = 14)
#'
#' # to plot the default graph (compact)
#' if(require(igraph)) {
#'  plot(tmiic.res)
#' }
#'
#' # to plot the raw temporal network Using igraph
#' if(require(igraph)) {
#'   plot(tmiic.res, display="raw")
#' }
#'
#' # to plot the full temporal network Using igraph
#' if(require(igraph)) {
#'   plot(tmiic.res, display="lagged")
#' }
#'
#' }
#'
miic <- function(input_data,
                 state_order = NULL,
                 true_edges = NULL,
                 black_box = NULL,
                 n_threads = 1,
                 cplx = c("nml", "mdl"),
                 orientation = TRUE,
                 ori_proba_ratio = 1,
                 ori_consensus_ratio = NULL,
                 propagation = TRUE,
                 latent = c("no", "yes", "orientation"),
                 n_eff = -1,
                 n_shuffles = 0,
                 conf_threshold = 0,
                 sample_weights = NULL,
                 test_mar = TRUE,
                 consistent = c("no", "orientation", "skeleton"),
                 max_iteration = 100,
                 consensus_threshold = 0.8,
                 tau = -1,
                 delta_tau = 1,
                 movavg = -1,
                 keep_max_data = FALSE,
                 negative_info = FALSE,
                 verbose = FALSE) {
  res <- NULL

  if (is.null(input_data)) {
    stop("The input data file is required")
  }
  if (!is.data.frame(input_data)) {
    stop("The input data is not a dataframe")
  }
  if (!is.null(state_order)) {
    err_code <- checkStateOrder(state_order, input_data)
    if (err_code != "0") {
      warning(paste(errorCodeToString(err_code),
          "state_order file will be ignored.", sep=", "), call.=FALSE)
      state_order <- NULL
    } else if (is.null(state_order$var_names)) {
      warning(paste("Named column var_names is required for state_order,",
                    "state_order file will be ignored."), call.=FALSE)
      state_order <- NULL
    }
  }
  if (!is.null(true_edges)) {
    err_code <- checkTrueEdges(true_edges)
    if (err_code != "0") {
      warning(paste(errorCodeToString(err_code),
          "true_edges file will be ignored.", sep=", "), call.=FALSE)
      true_edges <- NULL
    }
  }
  #
  # Check if we use regular or temporal mode
  #
  miic_mode = "regular"
  if (any(tau >= 1))
      miic_mode <- "temporal"
  if (!is.null(state_order)) {
    if (!is.null(state_order$tau)) {
      if (any (as.numeric (state_order$tau) >= 1) )
        miic_mode <- "temporal"
    }
  }

  if (miic_mode != "temporal")
    tau <- c(-1)
  else {
    #
    # In temporal mode, convert data into lagged samples
    # and lag other inputs accordingly
    #
    cat ("Using temporal mode of miic (CausalXtract)\n")
    nodes_not_lagged <- colnames(input_data)[-1]
    struct_ret <- tmiic.lag_inputs (input_data, tau, delta_tau=delta_tau,
                                    movavg=movavg, state_order=state_order,
                                    true_edges=true_edges, black_box=black_box,
                                    keep_max_data=keep_max_data)
    input_data <- struct_ret$input_data
    tau <- struct_ret$tau
    delta_tau <- struct_ret$delta_tau
    state_order <- struct_ret$state_order
    true_edges <- struct_ret$true_edges
    black_box <- struct_ret$black_box
  }

  # Remove rows with only NAs
  input_data <- input_data[rowSums(is.na(input_data)) != ncol(input_data), ]
  if (length(input_data) == 0) {
    stop("The input data is empty or contains only NAs")
  }

  cplx <- tryCatch(
    {
      match.arg(cplx)
    },
    error = function(e) {
      if (grepl("object .* not found", e$message)) {
        message(e, "")
        return("")
      }
      return(toString(cplx))
    }
  )
  cplx <- match.arg(cplx)

  latent <- tryCatch(
    {
      match.arg(latent)
    },
    error = function(e) {
      if (grepl("object .* not found", e$message)) {
        message(e, "")
        return("")
      }
      return(toString(latent))
    }
  )
  latent <- match.arg(latent)

  consistent <- tryCatch(
    {
      match.arg(consistent)
    },
    error = function(e) {
      if (grepl("object .* not found", e$message)) {
        message(e, "")
        return("")
      }
      return(toString(consistent))
    }
  )
  consistent <- match.arg(consistent)

  if (n_eff > nrow(input_data)) {
    stop(
      paste0(
        "The number of effective samples cannot be greater than the ",
        "number of samples."
      )
    )
  }

  if (propagation != TRUE && propagation != FALSE) {
    stop("The propagation type is not correct. Allowed types are TRUE or FALSE")
  }

  if (orientation != TRUE && orientation != FALSE) {
    stop("The orientation type is not correct. Allowed types are TRUE or FALSE")
  }

  if (is.null(ori_consensus_ratio)) {
    ori_consensus_ratio <- ori_proba_ratio
  }

  if (verbose) {
    cat("START miic...\n")
  }

  is_contextual <- rep(0, ncol(input_data))
  names(is_contextual) <- colnames(input_data)
  is_continuous <- sapply(input_data, is.numeric)
  # Parse "state_order" file
  if (!is.null(state_order)) {
    mismatch <- is.na(match(state_order$var_names, colnames(input_data)))
    if (any(mismatch)) {
      var_str <- paste(state_order$var_names[mismatch], collapse=", ")
      warning(paste("Variable(s)", var_str, "specified in state_order file",
          "will not match any name in input_data, and will be ignored."),
          call.=FALSE)
      state_order <- state_order[!mismatch, ]
    }
    not_found <- is.na(match(colnames(input_data), state_order$var_names))
    if (any(not_found)) {
      var_str <- paste(colnames(input_data)[not_found], collapse=", ")
      warning(paste("Variable(s)", var_str, "in input_data not found",
          "in state_order file."), call.=FALSE)
    }
    for (row in 1:nrow(state_order)) {
      col <- as.character(state_order[row, "var_names"])
      if (!is.null(state_order$var_type)) {
        if (state_order[row, "var_type"] == 0) {
          input_data[, col] <- factor(input_data[, col])
          is_continuous[[col]] <- FALSE
        }
      }
      if (!is.null(state_order$is_contextual)) {
        if (state_order[row, "is_contextual"] == 1) {
          is_contextual[[col]] <- 1
        }
      }
      if (!is.null(state_order$levels_increasing_order)) {
        order_string <- state_order[row, "levels_increasing_order"]
        if (!is.na(order_string)) {
          if (is_continuous[[col]] == TRUE) {
            warning(paste(col, "is considered as a continuous variable,",
                "the provided orders will be ignored."), call.=FALSE)
          } else {
            orders <- unlist(strsplit(as.character(order_string), ","))
            values <- unique(input_data[[col]][!is.na(input_data[[col]])])
            absent <- is.na(match(values, orders))
            if (any(absent)) {
              var_str <- ""
              if (length(values[absent]) > 10) {
                # Only show the first 3 if the list is too long
                var_str <- paste(paste(values[absent][1:3], collapse=", "),
                    "and", length(values[absent]) - 3, "more")
              } else {
                var_str <- paste(values[absent], collapse=", ")
              }
              warning(paste("Variable", col, "has value(s)", var_str,
                  "that will not match the provided orders,",
                  "the provided orders will be ignored."), call.=FALSE)
            }
          }
        }
      }
    }
  }
  # Check the number of unique values of continuous and discrete variables
  for (col in colnames(input_data)) {
    unique_values <- length(unique(input_data[[col]][!is.na(input_data[[col]])]))
    if (is_continuous[[col]] &&
      (unique_values <= 40) && (nrow(input_data) > 40)) {
      if (is_continuous[[col]] &&
        (unique_values <= 2) && (nrow(input_data) > 40)) {
        # Less than 3 unique values does not make sense for a continuous variable
        stop(
          paste0(
            "Numerical variable ",
            col,
            " only has ",
            unique_values,
            " non-NA unique values. Is this a factor?"
          )
        )
      }
      # Less than 40 unique variables can be discretized but may not be truly
      # continuous
      warning(
        paste0(
          "Numerical variable ",
          col,
          " is treated as continuous but only has ",
          unique_values,
          " unique values."
        ), call.=FALSE
      )
    }
    if ((!is_continuous[[col]]) &&
      (unique_values >= 40) && (nrow(input_data) > 40)) {
      warning(paste0(
        col,
        " is treated as discrete but has many levels (",
        unique_values,
        ")."
      ), call.=FALSE)
    }
  }

  err_code <- checkInput(input_data, "miic")
  if (err_code != "0") {
    warning(errorCodeToString(err_code), call.=FALSE)
  } else {
    if (verbose) {
      cat("\t# -> START reconstruction...\n")
    }
    res <-
      miic.reconstruct(
        input_data = input_data,
        n_threads = n_threads,
        cplx = cplx,
        latent = latent,
        n_eff = n_eff,
        black_box = black_box,
        n_shuffles = n_shuffles,
        orientation = orientation,
        ori_proba_ratio = ori_proba_ratio,
        propagation = propagation,
        conf_threshold = conf_threshold,
        verbose = verbose,
        is_contextual = is_contextual,
        is_continuous = is_continuous,
        sample_weights = sample_weights,
        test_mar = test_mar,
        consistent = consistent,
        max_iteration = max_iteration,
        tau = tau,
        delta_tau = delta_tau,
        negative_info = negative_info
      )
    if (res$interrupted) {
      stop("Interupted by user")
    }
    if (verbose) {
      cat("\t# -> END reconstruction...\n\t# --------\n")
    }

    res$all.edges.summary <- summarizeResults(
      observations = input_data,
      results = res,
      true_edges = true_edges,
      state_order = state_order,
      consensus_threshold = consensus_threshold,
      ori_consensus_ratio = ori_consensus_ratio,
      latent = latent != "no",
      propagation = propagation,
      verbose = verbose
    )
  }

  if (miic_mode != "temporal") {
    class(res) <- "miic"
  }
  else {
    class(res) <- "tmiic"
    res$tmiic_specific <- list (graph_type="raw",
                                nodes_not_lagged=nodes_not_lagged,
                                is_contextual=is_contextual,
                                tau=tau,
                                delta_tau=delta_tau)
  }
  return(res)
}


#' Basic plot function of a miic network inference result
#'
#' @description This function calls \code{\link{miic.export}} to build a
#' plottable object from the result returned by \code{\link{miic}} and plot it.
#'
#' @details See the documentation of \code{\link{miic.export}} for further
#' details.
#'
#' @param x [a miic graph object]
#' The graph object returned by \code{\link{miic}}.
#' @param method A string representing the plotting method. Default to "igraph".
#' Currently only "igraph" is supported.
#' @param pcor_palette Optional. The color palette used to represent the partial
#' correlations (the color of the edges). See \code{\link{miic.export}} for details.
#' @param \dots Additional plotting parameters. See the corresponding plot function
#' for the complete list.
#' For igraph, see \code{\link[igraph]{igraph.plotting}}.
#'
#' @export
#'
#' @seealso \code{\link{miic.export}} for generic exports,
#' \code{\link{getIgraph}} for igraph export,
#' \code{\link[igraph]{igraph.plotting}}
#'
plot.miic = function(x, method = 'igraph', pcor_palette = NULL, ...) {
  if (class(x) != "miic"){
    stop("Not a miic object.")
  }
  if (method == 'igraph'){
    if (base::requireNamespace("igraph", quietly = TRUE)) {
      igraph_obj = miic.export (x, 'igraph', pcor_palette = pcor_palette)
      igraph::plot.igraph (igraph_obj, ...)
    } else {
      stop("Package 'igraph' is required.")
    }
  } else {
    stop("Method not supported. See ?miic.export for supported methods.")
  }
}
