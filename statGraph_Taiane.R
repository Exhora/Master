# -----------------------------------------------------------------------------
# Statistical methods to analyze graphs
# -----------------------------------------------------------------------------

# Helper functions ------------------------------------------------------------

# Given a partition x[1]...x[n] and y[i] = f(x[i]), returns the trapezoid sum
# approximation for int_{x[1]}^{x[n]}{f(x)dx}
trapezoidSum <- function (x, y) {
    n <- length(x)
    delta <- (x[2] - x[1])
    area <- sum(y[2:(n-1)])
    area <- (area + (y[1] + y[n])/2)*delta
    return(area)
}

# Returns the kernel bandwidth for a given sample x
kernelBandwidth <- function(x) {
    n <- length(x)
    # Sturges' criterion
    nbins <- ceiling(log2(n) + 1)
    return(abs(max(x) - min(x))/nbins)
}

# Returns the density function for a given sample x at n points in the interval
# [form, to]
gaussianDensity <- function(x, from=NULL, to=NULL, bandwidth="Silverman", 
                            npoints=1024) {
    if (bandwidth == "Sturges")
        bw <- kernelBandwidth(x)
    else if (bandwidth == "Silverman")
        bw <- bw.nrd0(x)
    if (bw == 0)
       return(NA)
    if (is.null(from) || is.null(to))
        f <- density(x, bw=bw, n=npoints)
        #f <- density(x, n=npoints)

    else
        f <- density(x, bw=bw, from=from, to=to, n=npoints)
        #f <- density(x, from=from, to=to, n=npoints)

    area <- trapezoidSum(f$x, f$y)
    return(list("x"=f$x, "y"=f$y/area))
}

# Returns the spectral density for a given adjacency matrix A
spectralDensity <- function(A, from=NULL, to=NULL, bandwidth="Silverman", 
                            npoints=1024) {
    eigenvalues <- (as.numeric(eigen(A, only.values = TRUE)$values)/
                    sqrt(nrow(A)))
    return(gaussianDensity(eigenvalues, from, to, bandwidth=bandwidth, 
                           npoints=npoints))
}

# Returns the spectral densities for given adjacency matrices A1 and A2 at the
# same points
spectralDensities <- function(A1, A2, bandwidth="Silverman",
                              npoints=1024) {
    n1 <- nrow(A1)
    n2 <- nrow(A2)
    e1 <- (as.numeric(eigen(A1, only.values = TRUE)$values)/sqrt(n1))
    e2 <- (as.numeric(eigen(A2, only.values = TRUE)$values)/sqrt(n2))
    #b1 <- kernelBandwidth(e1)
    #b2 <- kernelBandwidth(e2)
    #from <- min(min(e1) - 3*b1, min(e2) - 3*b2)
    #to <- max(max(e1) + 3*b1, max(e2) + 3*b2)
    from <- min(e1, e2)
    to <- max(e1, e2)
    f1 <- gaussianDensity(e1, from=from, to=to, bandwidth=bandwidth, 
                          npoints=npoints)
    f2 <- gaussianDensity(e2, from=from, to=to, bandwidth=bandwidth, 
                          npoints=npoints)
    if (sum(is.na(f1)) > 0 || sum(is.na(f2)) > 0)
        return(NA)
    return(list("f1"=f1, "f2"=f2))
}

# Returns the spectral densities for a list of adjacency matrices at the
# same points
nSpectralDensities <- function (adjacencyMatrices, from=NULL, to=NULL, 
                                bandwidth="Silverman", npoints=1024) {
    if (class(adjacencyMatrices) == "matrix")
        return(spectralDensity(adjacencyMatrices, from, to, bandwidth, npoints))    
    ngraphs <- length(adjacencyMatrices)
    n <- ncol(adjacencyMatrices[[1]])
    spectra <- matrix(NA, n, ngraphs)
    for (i in 1:ngraphs) {
        A <- adjacencyMatrices[[i]]
        eigenvalues <- (as.numeric(eigen(A, only.values = TRUE)$values)/
                       sqrt(nrow(A)))
        spectra[,i] <- eigenvalues
    }
    densities <- matrix(NA, npoints, ngraphs)
    minimum <- min(spectra)
    maximum <- max(spectra)
    if (!is.null(from) && !is.null(to)) {
        minimum <- min(minimum, from)
        maximum <- max(maximum, to)
    }
    for (i in 1:ngraphs) {
        f <- gaussianDensity(spectra[,i], bandwidth=bandwidth,
                             from=minimum, to=maximum,
                             npoints=npoints)

        densities[,i] <- f$y
        x <- f$x
    }
    return(list("x"=x, "densities"=densities))
}

# Estimates the spectral density of a graph model
modelSpectralDensity <- function(fun, n, p, K=0, ngraphs=50, from=NULL, to=NULL, 
                                 bandwidth="Silverman", npoints=1024) {
    spectra <- matrix(NA, n, ngraphs)
    for (i in 1:ngraphs) {
	if (K == 0){
        A <- fun(n, p)
	}else{
        A <- fun(n, p, K)
	}
        eigenvalues <- (as.numeric(eigen(A, only.values = TRUE)$values)/
                       sqrt(nrow(A)))
        spectra[,i] <- eigenvalues
    }
    densities <- matrix(NA, npoints, ngraphs)
    minimum <- min(spectra)
    maximum <- max(spectra)
    if (!is.null(from) && !is.null(to)) {
        minimum <- min(minimum, from)
        maximum <- max(maximum, to)
    }
    for (i in 1:ngraphs) {
        f <- gaussianDensity(spectra[,i], from=minimum, to=maximum, 
                             bandwidth=bandwidth, npoints=npoints)

        densities[,i] <- f$y
        x <- f$x
    }
    return(list("x" = x, "y" = rowMeans(densities)))
}

matchFunction <- function(name) {
    return(match.fun(name))
}

# Graph models  ----------------------------------------------------------------

# Erdos-Renyi graph 
ER <- function(n, p) {
    return(as.matrix(igraph::get.adjacency(igraph::erdos.renyi.game(n, p, 
                                                                   type="gnp",
                                                            directed = FALSE))))
}

# Geometric graph
GRG <- function(n, r) {
    return(as.matrix(igraph::get.adjacency(igraph::grg.game(n, r))))
}

# Barabasi-Albert graph
BA <- function(n, ps) {
    return(as.matrix(igraph::get.adjacency(igraph::barabasi.game(n, power = ps,
                                                 directed = FALSE))))
}

# Watts-Strougatz graph
WS <- function(n, pr, K=4) {
    return(as.matrix(igraph::get.adjacency(igraph::watts.strogatz.game(1, n, K, 
                                                                       pr, ))))
}

# K-regular game
KR <- function(n, k) {
    return(as.matrix(get.adjacency(igraph::k.regular.game(n, k))))
}

# Entropy ---------------------------------------------------------------------

#' Graph spectral entropy
#'
#' Given an adjacency matrix A, it returns the spectral entropy of the 
#' corresponding graph.
#' @param A a symmetric matrix of nonnegative real values. For unweighted graphs 
#' it contains only 0s and 1s and corresponds to the adjacency matrix. For 
#' weighted graphs, it may contain any nonnegative real value and corresponds to
#' the weighted adjancency matrix.
#' @param bandwidth string indicating which criterion should be used
#' to choose the bandwidth for the spectral density estimation. The available 
#' criteria are "Silverman" (default) and "Sturges".
#'
#' @return a real number corresponding to the graph spectral entropy.  
#'
#' @references
#' Daniel Yasumasa Takahashi, João Ricardo Sato, Carlos Eduardo Ferreira, and 
#' André Fujita. Discriminating Different Classes of Biological Networks by 
#' Analyzing the Graph Spectra Distribution. PLoS ONE 7, no. 12 
#' (December 2012): e49949. doi:10.1371/journal.pone.0049949.
#' 
#' @examples
#' library(igraph)
#' G <- erdos.renyi.game(100, p=0.5)
#' A <- as.matrix(get.adjacency(G))
#' entropy <- graphEntropy(A)
#' entropy
#'
#' @export
graphEntropy <- function(A, bandwidth="Silverman") {
    f <- spectralDensity(A, bandwidth=bandwidth)
    if (sum(is.na(f)) > 0)
        return(NA)
    y <- f$y
    n <- length(y)
    i <- which(y != 0)
    y[i] <- y[i]*log(y[i])
    return(-trapezoidSum(f$x, y))
}

# Kullback-Leibler divergence between graphs -----------------------------------

# Given two spectral densities, returns the Kullback-Leibler divergence between
# the corresponding graphs
KL <- function(f1, f2) {
    y <- f1$y
    n <- length(y)
    for (i in 1:n) {
        if (y[i] != 0 && f2$y[i] == 0)
            return(Inf)
        if (y[i] != 0)
            y[i] <- y[i]*log(y[i]/f2$y[i])
    }
    return(trapezoidSum(f1$x, y))
}

#' Graph Information Criterion (GIC)
#'
#' Returns the Kullback-Leibler divergence between an observed graph and a given 
#' graph model. 
#'
#' @param adjacencyMatrices a matrix or a list of symmetric matrix of nonnegative 
#' real values. For unweighted graphs each matrix contains only 0s and 1s and 
#' corresponds to the adjacency matrix. For weighted graphs, it may contain any 
#' nonnegative real value and corresponds to the weighted adjancency matrix.
#' @param model either a list, a string or a function describing a graph model: 
#' 
#' A list that represents the spectral density of a model. It contains the 
#' components "x" and "y", which are numeric vectors of the same size. The "x" 
#' component contains the points at which the density was computed and the "y" 
#' component contains the observed density.
#'
#' A string that indicates one of the following models: "ER" (Erdos-Renyi random 
#' graph), "GRG" (geometric random graph), "KR" (k regular 
#' random graph), "WS" (Watts-Strogatz model), and "BA" (Barabási-Albert model).
#' When the argument model is a string, the user must also provide the model 
#' parameter by using the argument p.  
#'
#' A function that returns a graph according to a graph model. It must have two 
#' arguments: n (graph size) and p (a model parameter). The model parameter 
#' will be provided by the argument p of the GIC function.
#' @param p the model parameter. The user must provide a valid parameter if the
#' argument "model" is a string or a function. For the "ER", "GRG", "KR", "WS", 
#' and "BA" models, the parameter p corresponds to the probability to connect a 
#' pair of vertices, the radius used to contruct the geometric graph in a unit 
#' square, the degree k of a regular graph, the probability to reconnect a 
#' vertex, and the scale exponent, respectively.
#' @param bandwidth string indicating which criterion should be used
#' to choose the bandwidth for the spectral density estimation. The available 
#' criteria are "Silverman" (default) and "Sturges".
#'
#' @return A real number corresponding to the Kullback-Leibler divergence 
#' between the observed graph and a given model.
#'
#' @references
#' Daniel Yasumasa Takahashi, João Ricardo Sato, Carlos Eduardo Ferreira, and 
#' André Fujita. Discriminating Different Classes of Biological Networks by 
#' Analyzing the adjacencyMatrices Spectra Distribution. PLoS ONE 7, no. 12 
#' (December 2012): e49949. doi:10.1371/journal.pone.0049949.
#' 
#' @examples
#' library(igraph)
#' A <- as.matrix(get.adjacency(erdos.renyi.game(100, p=0.5)))
#' # Using a string to indicate the graph model
#' result1 <- GIC(A, "ER", 0.5)
#' result1
#'
#' # Using a function to describe the graph model
#' # Erdos-Renyi graph
#' model <- function(n, p) {
#'    return(as.matrix(get.adjacency(erdos.renyi.game(n, p))))
#' }
#' result2 <- GIC(A, model, 0.5)
#' result2
#' @export
GIC <- function(adjacencyMatrices, model, p=NULL, bandwidth="Silverman", K=0) {
   	npoints <- 1024    
    if (class(model) == "list") {
        f2 <- model
        f1 <- nSpectralDensities(adjacencyMatrices, from=min(f2$x), to=max(f2$x), 
                              bandwidth=bandwidth, npoints=length(f2$x))
        if (class(adjacencyMatrices) == "list")
            f1$y <- rowMeans(f1$densities) 
    }
    else {
        fun <- model
        if (class(model) == "character")
            fun <- matchFunction(model)
        if (class(adjacencyMatrices) == "matrix") {
            spectra  <- (as.numeric(eigen(adjacencyMatrices, only.values=TRUE)$
                         values)/sqrt(nrow(adjacencyMatrices)))
            n <- ncol(adjacencyMatrices)        
        }
        else {
            n <- ncol(adjacencyMatrices[[1]])
            ngraphs <- length(adjacencyMatrices)
            spectra <- matrix(NA, n, ngraphs)
            for (i in 1:ngraphs) {
                A <- adjacencyMatrices[[i]]
                eigenvalues <- (as.numeric(eigen(A, only.values = TRUE)$values)/
                               sqrt(nrow(A)))
                spectra[,i] <- eigenvalues
            }
        }
	
	f2 <- modelSpectralDensity(fun, n, p, K, from=min(spectra), 
                                   to=max(spectra), bandwidth=bandwidth,
                                   npoints=1024)
	
	
	if (class(adjacencyMatrices) == "matrix") {
            f1 <- gaussianDensity(spectra, from=min(f2$x), to=max(f2$x), 
                                 bandwidth=bandwidth, npoints=1024)
        }
        else {
            densities <- matrix(NA, npoints, ngraphs)
            for (i in 1:ngraphs) {
                f <- gaussianDensity(spectra[,i], from=min(f2$x), to=max(f2$x),
                                     bandwidth=bandwidth, npoints=1024)

                densities[,i] <- f$y
                x <- f$x
            }
            f1 <- list()
            f1$x <- f$x
            f1$y <- rowMeans(densities)
        }
    }
    if (sum(is.na(f1)) > 0)
        return(NA)
    if (sum(is.na(f2)) > 0)
        return(NA)
    return(KL(f1, f2))
}

#' Graph parameter estimator
#'
#' Estimates the model parameter that best approximates the model to the 
#' observed graph. It can take a long time to execute, specially for the k 
#' regular graph model.
#'
#' @param A a symmetric matrix of nonnegative real values. For unweighted graphs 
#' it contains only 0s and 1s and corresponds to the adjacency matrix. For 
#' weighted graphs, it may contain any nonnegative real value and corresponds to
#' the weighted adjancency matrix.
#' @param model either a string or a function describing a graph model: 
#' 
#' A string that indicates one of the following models: "ER" (Erdos-Renyi random 
#' graph), "GRG" (geometric random graph), "KR" (k regular 
#' random graph), "WS" (Watts-Strogatz model), and "BA" (Barabási-Albert model).
#'
#' A function that returns a graph according to a graph model. It must have two 
#' arguments: n (graph size) and p (a model parameter).
#' @param parameters numeric vector containing the values that will be 
#' considerated for the model parameter. For the "ER", "GRG", "KR", "WS", and 
#' "BA" models, the parameter corresponds to the probability to connect a pair 
#' of vertices, the radius used to contruct the geometric graph in a unit 
#' square, the degree k of a regular graph, the probability to reconnect a 
#' vertex, and the scale exponent, respectively. For those models there are 
#' default parameter values that will be considered, and then it is optional to 
#' use this argument. However if the user provides a custom function to generate 
#' graphs, then it is also necessary to provide the values that will be 
#' considered for the estimation. 
#' @param bandwidth string indicating which criterion should be used
#' to choose the bandwidth for the spectral density estimation. The available 
#' criteria are "Silverman" (default) and "Sturges".
#'
#' @return A list containing:
#' \item{p}{the parameter estimate. For the "ER", "GRG", "KR", "WS", and "BA"
#' models, the parameter corresponds to the probability to connect a pair of 
#' vertices, the radius used to contruct the geometric graph in a unit square,
#' the degree k of a regular graph, the probability to reconnect a vertex, and 
#' the scale exponent, respectively.}  
#' \item{KL}{the Kullback-Leibler divergence between the observed graph and the 
#' model obtained by the estimator.}
#'
#' @references
#' Daniel Yasumasa Takahashi, João Ricardo Sato, Carlos Eduardo Ferreira, and 
#' André Fujita. Discriminating Different Classes of Biological Networks by 
#' Analyzing the adjacencyMatrices Spectra Distribution. PLoS ONE 7, no. 12 
#' (December 2012): e49949. doi:10.1371/journal.pone.0049949.
#' 
#' @examples
#' library(igraph)
#' A <- as.matrix(get.adjacency(erdos.renyi.game(100, p=0.5)))
#'
#' # Using a string to indicate the graph model
#' result1 <- graphParamEstimator(A, "ER")
#' result1
#'
#' # Using a function to describe the graph model
#' # Erdos-Renyi graph
#' model <- function(n, p) {
#'    return(as.matrix(get.adjacency(erdos.renyi.game(n, p))))
#' }
#' result2 <- graphParamEstimator(A, model,  seq(0, 1, 0.01))
#' result2
#'
#' @export
graphParamEstimator <- function(A, model, parameters=NULL, 
                                  bandwidth="Silverman") {
    n <- ncol(A)
    if (class(model) == "function" && is.null(parameters)) {
        stop("It is necessary to enter the paremeters that will be evaluated.")
    }
    if (is.null(parameters)) {
            parameters <- seq(0, 1, 0.01)
            if (model == "GRG")
                parameters <- seq(0, sqrt(2), 0.01)
            if (model == "BA")
                parameters <- seq(1, 4, 0.01)
            if (model == "KR")
                parameters <- as.integer(seq(0, 1, 0.01)*n)
    }
    K <- 0
	if(model == "WS"){
	arestas <- 0
	#Calculo do K para WS
	if(class(A)=="list"){
	for(i in 1:length(A)){
	    arestas <- arestas + (sum(A[[i]])/2)
	}
	K <- round(arestas/(length(A)*n))
	}else{
		K <- round(sum(A)/(n*2))
	}
	}

    pmin <- -1
    klmin <- Inf
    kmin <- 0
    for (p in parameters) {
       kl <- GIC(A, model, p, bandwidth, K) 
       if (kl < klmin) {
           klmin <- kl
           pmin <- p
	   if(K!=0){
		kmin <- K
	   }
       }
    }
    return(list("p"=pmin, "KL"=klmin, "K"= kmin))
}

#' Graph model selection
#'
#' Selects the graph model that best approximates the observed graph. It can 
#' take a long time to execute, specially for the k regular graph model.
#'
#' @param A a symmetric matrix of nonnegative real values. For unweighted graphs 
#' it contains only 0s and 1s and corresponds to the adjacency matrix. For 
#' weighted graphs, it may contain any nonnegative real value and corresponds to
#' the weighted adjancency matrix.
#' @param models either a character vector or a list of functions describing 
#' graph models: 
#' 
#' A character vector cotaining some of the following models: "ER" (Erdos-Renyi 
#' random graph), "GRG" (geometric random graph), "KR" (k regular random graph), 
#' "WS" (Watts-Strogatz model), and "BA" (Barabási-Albert model).
#'
#' A list of functions. Each function returns a graph according to a graph model
#' and has two arguments: n (graph size) and p (a model parameter), in this 
#' order.
#' 
#' If the argument "models" is NULL, then the "ER", "WS", and "BA" models will 
#' be considered.
#' @param parameters list of numeric vectors. Each vector contains the values 
#' that will be considerated for the parameter of the corresponding model.
#' For the "ER", "GRG", "KR", "WS", and "BA" models, the parameter corresponds 
#' to the probability to connect a pair of vertices, the radius used to contruct 
#' the geometric graph in a unit square, the degree k of a regular graph, the 
#' probability to reconnect a vertex, and the scale exponent, respectively. For 
#' those models there are default parameter values that will be considered, and 
#' then it is optional to use this argument. However if the user provides custom 
#' functions to generate graphs, then it is also necessary to provide the values 
#' that will be considered for the parameter estimation. 
#' @param bandwidth string indicating which criterion should be used
#' to choose the bandwidth for the spectral density estimation. The available 
#' criteria are "Silverman" (default) and "Sturges".
#' @return A matrix in which each row corresponds to a model, the column "p"
#' corresponds to the parameter estimate, and the column "KL" corresponds to
#' the Kullback-Leibler divergence between the observed graph and the model in 
#' the row. 
#' @references
#' Daniel Yasumasa Takahashi, João Ricardo Sato, Carlos Eduardo Ferreira, and 
#' André Fujita. Discriminating Different Classes of Biological Networks by 
#' Analyzing the adjacencyMatrices Spectra Distribution. PLoS ONE 7, no. 12 
#' (December 2012): e49949. doi:10.1371/journal.pone.0049949.
#' 
#' @examples
#' library(igraph)
#' A <- as.matrix(get.adjacency(erdos.renyi.game(100, p=0.5)))
#'
#' # Using the default graph models (ER, WS and BA)
#' result1 <- graphModelSelection(A)
#' result1
#'
#' # Using strings to indicate the graph models
#' result2 <- graphModelSelection(A, models=c("ER", "WS", "BA"))
#' result2
#'
#' # Using functions to describe the graph models
#' # Erdos-Renyi graph
#' model1 <- function(n, p) {
#'    return(as.matrix(get.adjacency(erdos.renyi.game(n, p))))
#' }
#' # Watts-Strougatz graph
#' model2 <- function(n, pr, K=8) {
#'     return(as.matrix(igraph::get.adjacency(igraph::watts.strogatz.game(1, n, 
#'                                              K, pr))))
#'}
#' # Barabasi-Albert graph
#' model3 <- function(n, ps) {
#'     return(as.matrix(igraph::get.adjacency(igraph::barabasi.game(n, power=ps,
#'                                                 directed = FALSE))))
#'}
#' parameters <- list(seq(0, 1, 0.01), seq(0, 1, 0.01), seq(1, 4, 0.01))
#' result3 <- graphModelSelection(A, list(model1, model2, model3), parameters)
#' result3
#' @export
graphModelSelection <- function(A, models=NULL, parameters=NULL, 
                                  bandwidth="Silverman") {
    n <- ncol(A)
    if (class(models) == "list" && is.null(parameters)) {
        stop("It is necessary to enter the paremeters that will be evaluated.")
    }
    if (is.null(models)) {
        models <- c("ER", "WS", "BA")
    }
    results <- matrix(NA, length(models), 3)
    colnames(results) <- c("p", "KL", "K")
    if (class(models) == "character") {
        rownames(results) <- models
    }
    p <- NULL
    for (i in 1:length(models)) {
        if (!is.null(parameters))
            p <- parameters[[i]]
        r <- graphParamEstimator(A, models[[i]], p, bandwidth) 
        results[i, "p"] <- r$p
        results[i, "KL"] <- r$KL
	results[i, "K"] <- r$K
    }
    return(results)
}

# Jensen-Shannon divergence ---------------------------------------------------

# Given two spectral densities, returns the Jensen-Shannon divergence between
# the corresponding adjacencyMatrices
#' @export
JS <- function(f1, f2) {
    fm <- f1
    fm$y <- (f1$y + f2$y)/2
    return((KL(f1, fm) + KL(f2, fm))/2)
}

# Given two adjacency matrices, returns the Jensen-Shannon divergence between
# the corresponding adjacencyMatrices
#' @export
jensenShannon <- function(A1, A2, bandwidth="Silverman") {
    f <- spectralDensities(A1, A2, bandwidth=bandwidth)
    if (sum(is.na(f)) > 0)
        return(NA)
    f1 <- f$f1
    f2 <- f$f2
    fm <- f1
    fm$y <- (f1$y + f2$y)/2
    return((KL(f1, fm) + KL(f2, fm))/2)
}

#' Bootstrap test for the graph Jensen-Shannon (JS) divergence.
#'
#' Given two lists of graphs, x and y, it tests H0: "JS divergence between x and  
#' y is 0" against H1: "JS divergence between x and y is larger than 0".
#'
#' @param x a list of symmetric matrices with nonnegative real 
#' entries. For unweighted graphs, each matrix contains only 0s and 1s and 
#' corresponds to an adjacency matrix. For weighted graphs, each matrix may 
#' contain any nonnegative real value and corresponds to the weighted adjancency 
#' matrix.
#' @param y a list of symmetric matrices with nonnegative real 
#' entries. For unweighted graphs, each matrix contains only 0s and 1s and 
#' corresponds to an adjacency matrix. For weighted graphs, each matrix may 
#' contain any nonnegative real value and corresponds to the weighted adjancency 
#' matrix.
#' @param numBoot integer indicating the number of bootstrap resamplings.
#' @param bandwidth character string indicating which criterion should be used
#' to choose the bandwidth for the spectral density estimation. The available 
#' criteria are "Silverman" (default) and "Sturges".
#'
#' @return A list containing:
#' \item{JS}{the Jensen-Shannon divergence between the two graph sets.}  
#' \item{pvalue}{the p-value of the test.}
#'
#' @references
#' Daniel Yasumasa Takahashi, João Ricardo Sato, Carlos Eduardo Ferreira, and 
#' André Fujita. Discriminating Different Classes of Biological Networks by 
#' Analyzing the adjacencyMatrices Spectra Distribution. PLoS ONE 7, no. 12 
#' (December 2012): e49949. doi:10.1371/journal.pone.0049949.
#' 
#' @examples
#' library(igraph)
#' x <- y <- list()
#' for (i in 1:40)
#'    x[[i]] <- as.matrix(get.adjacency(erdos.renyi.game(50, p=0.5)))
#' for (i in 1:45)
#'    y[[i]] <- as.matrix(get.adjacency(erdos.renyi.game(50, p=0.2)))
#' 
#' result <- graphTest(x, y)
#' result
#'
#' @export
graphTest <- function(x, y, numBoot=1000, bandwidth="Silverman") {
    adjacencyMatrices <- append(x, y)
    labels <- c(rep(0, length(x)), rep(1, length(y)))
    f <- nSpectralDensities(adjacencyMatrices)
    densities <- f$densities
    x <- f$x
    y1 <- rowMeans(densities[, labels==0])
    y2 <- rowMeans(densities[, labels==1])
    n1 <- length(which(labels==0))
    n2 <- length(which(labels==1))
    results <- vector(len=numBoot)
    ngraphs <- length(adjacencyMatrices)
    result <- JS(list("x"=x, "y"=y1), list("x"=x, "y"=y2))
    for (i in 1:numBoot) {
        b1 <- sample(1:ngraphs, n1, replace=T)
        b2 <- sample(1:ngraphs, n2, replace=T)
        y1 <- rowMeans(densities[, b1])
        y2 <- rowMeans(densities[, b2])
        results[i] <- JS(list("x"=x, "y"=y1), list("x"=x, "y"=y2))
    }
    pvalue <- (sum(results >= result))/numBoot
    return(list("JS"=result, "pvalue"=pvalue))
}
