---
name: GraphicalModels
topic: Graphical Models
maintainer: Soren Hojsgaard
email: sorenh@math.aau.dk
version: 2021-12-27
source: https://github.com/cran-task-views/GraphicalModels/
---

[Wikipedia](http://en.wikipedia.org/wiki/Graphical_model) says: 

*A graphical model or probabilistic graphical model (PGM) or
structured probabilistic model is a probabilistic model for which a
graph expresses the conditional dependence structure between random
variables. They are commonly used in probability theory, statistics —
particularly Bayesian statistics — and machine learning.*


A supplementary view is that graphical models are based on exploiting
conditional independencies for constructing complex stochastic models
with a modular structure. That is, a complex stochastic model is built
up by simpler building blocks. 

This task view is a collection of
packages intended to supply R code to deal with graphical models.

Notice that Structural Equation Models (SEM) packages are in a sense
also graphical models. However, SEM packages are not presented here
but are they have their own section in the `r view("Psychometrics")`
task view.

The packages can be roughly structured into the following topics
(although several of them have functionalities which go across these
categories):

### Representation, manipulation and display of graphs

-   `r pkg("backbone")`: An implementation of methods for extracting
    an unweighted unipartite graph (i.e. a backbone) from an
    unweighted unipartite graph, a weighted unipartite graph, the
    projection of an unweighted bipartite graph, or the projection of
    a weighted bipartite graph.
-   `r pkg("diagram")`: Visualises simple graphs (networks)
    based on a transition matrix, utilities to plot flow diagrams,
    visualising webs, electrical networks, \...
-   `r pkg("DiagrammeR")`: Build graph/network structures using
	functions for stepwise addition and deletion of nodes and edges.
-   `r bioc("graph")`: A package that implements some simple
    graph handling capabilities.
-   `r pkg("gRbase", priority = "core")`: The gRbase package
    provides certain general constructs which are used by other
    graphical modelling packages. This includes 1) the concept of gmData
    (graphical meta data), 2) several graph algorithms 3) facilities for
    table operations, 4) functions for testing for conditional
    independence. gRbase also illustrates how hierarchical log-linear
    models (hllm) may be implemented.
-   `r pkg("igraph")`: Routines for simple graphs and network
    analysis. It can handle large graphs very well and provides
    functions for generating random and regular graphs, graph
    visualization, centrality methods and much more.
-   `r pkg("network")`: Tools to create and modify network
    objects. The network class can represent a range of relational data
    types, and supports arbitrary vertex/edge/graph attributes.
-   `r pkg("qgraph")`: Weighted network visualization and
    analysis, as well as Gaussian graphical model computation. See
    Epskamp et al. (2012)
    [doi:10.18637/jss.v048.i04](https://dx.doi.org/doi:10.18637/jss.v048.i04)
-   `r bioc("Rgraphviz")`: Provides plotting capabilities for
    R graph objects.
-   `r bioc("RBGL")`: A fairly extensive and comprehensive
    interface to the graph algorithms contained in the BOOST library.
    (based on graph objects from the `r bioc("graph")`
    package).

### Classical models - General purpose packages

-   `r pkg("ggm")`: Fitting graphical Gaussian models.
-   `r pkg("gRbase")`: The gRbase package provides certain
    general constructs which are used by other graphical modelling
    packages (in particular by `r pkg("gRain")`). This
    includes 1) the concept of gmData (graphical meta data), 2) several
    graph algorithms 3) facilities for table operations, 4) functions
    for testing for conditional independence. gRbase also illustrates
    how hierarchical log-linear models (hllm) may be implemented. Link:
    [doi:10.18637/jss.v014.i17](https://dx.doi.org/doi:10.18637/jss.v014.i17)
-   `r pkg("sna")`: A range of tools for social network
    analysis, including node and graph-level indices, structural
    distance and covariance methods, structural equivalence detection,
    network regression, random graph generation, and 2D/3D network
    visualization.
-   `r pkg("mgm")`: Estimation of k-Order time-varying Mixed
    Graphical Models and mixed VAR(p) models via elastic-net regularized
    neighborhood regression.

### Miscellaneous: Model search, structure learning, specialized types of models etc.

-   `r pkg("abn")`: Modelling Multivariate Data with Additive Bayesian
	Networks.  An additive Bayesian network model consists of a form
	of a DAG where each node comprises a generalized linear model,
	GLM. Additive Bayesian network models are equivalent to Bayesian
	multivariate regression using graphical modelling, they
	generalises the usual multivariable regression, GLM, to multiple
	dependent variables. 'abn' provides routines to help determine
	optimal Bayesian network models for a given data set, where these
	models are used to identify statistical dependencies in messy,
	complex data.
	
-   `r pkg("BDgraph")`: Bayesian Graph Selection Based on
    Birth-Death MCMC Approach. Bayesian inference for structure learning
    in undirected graphical models. The main target is to uncover
    complicated patterns in multivariate data wherein either continuous
    or discrete variables.
-   `r pkg("bnclassify")`: Learning Discrete Bayesian Network
    Classifiers from Data. State-of-the art algorithms for learning
    discrete Bayesian network classifiers from data, including a number
    of those described in Bielza and Larranaga (2014) with functions for
    prediction, model evaluation and inspection. The package implements
    a number of structure-learning algorithms, with both discriminative
    and generative network scores, and a number of naive Bayes-specific
    parameter estimation methods, such as the Model Averaged Naïve
    Bayes. Prediction with complete data is rather fast, allowing for
    discriminative scores for structure learning with moderately
    sized-data sets (e.g. the splice data set, with 61 variables and
    3190 instances in \~100 seconds on my machine). There is a number of
    utility functionalities, such as cross-validation.
-   `r pkg("bnstruct")`: Bayesian Network Structure Learning
    from Data with Missing Values. Bayesian Network Structure Learning
    from Data with Missing Values. The package implements the
    Silander-Myllymaki complete search, the Max-Min
    Parents-and-Children, the Hill-Climbing, the Max-Min Hill-climbing
    heuristic searches, and the Structural Expectation-Maximization
    algorithm. Available scoring functions are BDeu, AIC, BIC. The
    package also implements methods for generating and using bootstrap
    samples, imputed data, inference.

-   `r pkg("catnet")`: Categorical Bayesian Network Inference
	Structure learning and parameter estimation of discrete Bayesian
	networks using likelihood-based criteria. Exhaustive search for fixed
	node orders and stochastic search of optimal orders via simulated
	annealing algorithm are implemented.

-   `r pkg("deal")`: Learning Bayesian Networks with Mixed Variables
	Bayesian networks with continuous and/or discrete variables can be
	learned and compared from data.

-   `r pkg("FBFsearch")`: Algorithm for searching the space
    of Gaussian directed acyclic graphical models through moment
    fractional Bayes factors

-   `r pkg("GeneNet")`: Modeling and Inferring Gene Networks.
    GeneNet is a package for analyzing gene expression (time series)
    data with focus on the inference of gene networks.

-   `r pkg("gRc")`: Inference in Graphical Gaussian Models with Edge
	and Vertex Symmetries Estimation, model selection and other aspects of
	statistical inference in Graphical Gaussian models with edge and
	vertex symmetries (Graphical Gaussian models with colours).
	
-	`r pkg("gRim")`: Graphical Interaction Models.  Provides the
    following types of models: Models for contingency tables
	(i.e. log-linear models) Graphical Gaussian models for multivariate
	normal data (i.e. covariance selection models) Mixed interaction
	models.

-   `r pkg("huge")`: High-dimensional Undirected Graph
    Estimation.
-   `r pkg("lvnet")`: Latent Variable Network
    Modeling. Estimate, fit and compare Structural Equation Models (SEM)
    and network models (Gaussian Graphical Models; GGM) using OpenMx.
    Allows for two possible generalizations to include GGMs in SEM: GGMs
    can be used between latent variables (latent network modeling; LNM)
    or between residuals (residual network modeling; RNM).

-   `r pkg("MXM")`: Feature Selection (Including Multiple Solutions)
    and Bayesian Networks.

-   `r pkg("networkDynamic")`: Dynamic Extensions for Network
    Objects. Simple interface routines to facilitate the handling of
    network objects with complex intertemporal data. "networkDynamic"
    is a part of the "statnet" suite of packages for network analysis.

-   `r pkg("ndtv")`: Network Dynamic Temporal Visualizations.  Renders
	dynamic network data from 'networkDynamic' objects as movies,
	interactive animations, or other representations of changing
	relational structures and attributes.

-   `r pkg("pcalg")`: Standard and robust estimation of the
    skeleton (ugraph) and the equivalence class of a Directed Acyclic
    Graph (DAG) via the PC-Algorithm. The equivalence class is
    represented by its (unique) Completed Partially Directed Acyclic
    Graph (CPDAG).
-   `r rforge("qp")`: This package is deprecated and it is
    now only a stub for the newer version called qpgraph available
    through the Bioconductor project. The q-order partial correlation
    graph search algorithm, q-partial, or qp, algorithm for short, is a
    robust procedure for structure learning of undirected Gaussian
    graphical Markov models from "small n, large p" data, that is,
    multivariate normal data coming from a number of random variables p
    larger than the number of multidimensional data points n as in the
    case of, e.g., microarray data.
-   `r bioc("qpgraph")`: q-order partial correlation graphs,
    or qp-graphs for short, are undirected Gaussian graphical Markov
    models that represent q-order partial correlations. They are useful
    for learning undirected graphical Gaussian Markov models from data
    sets where the number of random variables p exceeds the available
    sample size n as, for instance, in the case of microarray data where
    they can be employed to reverse engineer a molecular regulatory
    network.
<!-- -   `r pkg("SIN")`: This package provides routines to perform -->
    <!-- SIN model selection as described in Drton and Perlman (2004). The -->
    <!-- selected models are represented in the format of the 'ggm' -->
    <!-- package, which allows in particular parameter estimation in the -->
    <!-- selected model. -->
-   `r pkg("spectralGraphTopology")`: The package provides
    estimators to learn k-component, bipartite, and k-component
    bipartite graphs from data by imposing spectral constraints on the
    eigenvalues and eigenvectors of the Laplacian and adjacency
    matrices. Those estimators leverages spectral properties of the
    graphical models as a prior information, which turn out to play key
    roles in unsupervised machine learning tasks such as community
    detection.

### Bayesian Networks/Probabilistic expert systems

-   `r pkg("bnlearn")`: Bayesian network structure learning
    via constraint-based (also known as 'conditional independence')
    and score-based algorithms. This package implements the Grow-Shrink
    (GS) algorithm, the Incremental Association (IAMB) algorithm, the
    Interleaved-IAMB (Inter-IAMB) algorithm, the Fast-IAMB (Fast-IAMB)
    algorithm, the Max-Min Parents and Children (MMPC) algorithm and the
    Hill-Climbing (HC) greedy search algorithm for both discrete and
    Gaussian networks, along with many score functions and conditional
    independence tests. Some utility functions (model comparison and
    manipulation, random data generation, arc orientation testing) are
    also included.
-   `r pkg("gRain")`: A package for probability propagation
    in graphical independence networks, also known as probabilistic
    expert systems (which includes Bayesian networks as a special case).
    Link:
    [doi:10.18637/jss.v046.i10](https://dx.doi.org/doi:10.18637/jss.v046.i10)
-   `r rforge("RHugin")`: The Hugin Decision Engine (HDE) is
    commercial software produced by HUGIN EXPERT A/S for building and
    making inference from Bayesian belief networks. The RHugin package
    provides a suite of functions allowing the HDE to be controlled from
    within the R environment for statistical computing. The RHugin
    package can thus be used to build Bayesian belief networks, enter
    and propagate evidence, and to retrieve beliefs. Additionally, the
    RHugin package can read and write hkb and NET files, making it easy
    to work simultaneously with both the RHugin package and the Hugin
    GUI. A licensed copy of the HDE (or the trial version) is required
    for the RHugin package to function, hence the target audience for
    the package is Hugin users who would like to take advantage of the
    statistical and programmatic capabilities of R. Notice: RHugin is
    NOT on CRAN. Link: <http://rhugin.r-forge.r-project.org/>
<!-- -   `r pkg("sparsebn")`: Fast methods for learning sparse -->
    <!-- Bayesian networks from high-dimensional data using coordinate -->
    <!-- descent and sparse regularization. Designed to handle mixed -->
    <!-- experimental and observational data with thousands of variables with -->
    <!-- either continuous or discrete observations. -->


### BUGS models

-   `r pkg("bayesmix")`: Bayesian mixture models of
    univariate Gaussian distributions using JAGS.
-   `r pkg("dclone")`: Data Cloning and MCMC Tools for
    Maximum Likelihood Methods. Low level functions for implementing
    maximum likelihood estimating procedures for complex models using
    data cloning and Bayesian Markov chain Monte Carlo methods with
    support for JAGS, WinBUGS and OpenBUGS. Parallel MCMC computation is
    supported and can result in considerable speed-up.
-   `r pkg("boa")`: Bayesian Output Analysis Program
    (BOA) for MCMC. A menu-driven program and library of functions for
    carrying out convergence diagnostics and statistical and graphical
    analysis of Markov chain Monte Carlo sampling output.
-   `r pkg("BRugs")`: R interface to the OpenBUGS MCMC
    software. Fully-interactive R interface to the OpenBUGS software for
    Bayesian analysis using MCMC sampling. Runs natively and stably in
    32-bit R under Windows. Versions running on Linux and on 64-bit R
    under Windows are in "beta" status and less efficient.
-   `r pkg("coda")`: Output analysis and diagnostics
    for MCMC. Output analysis and diagnostics for Markov Chain Monte
    Carlo simulations.
-   `r pkg("ergm")`: Fit, Simulate and Diagnose
    Exponential-Family Models for Networks. An integrated set of tools
    to analyze and simulate networks based on exponential-family random
    graph models (ERGM). "ergm" is a part of the
    [statnet](http://statnet.org/) suite of packages for network
    analysis.
-   `r pkg("R2OpenBUGS")`: Running OpenBUGS from
    R. Using this package, it is possible to call a BUGS model,
    summarize inferences and convergence in a table and graph, and save
    the simulations in arrays for easy access in R.
-   `r pkg("R2WinBUGS")`: Running WinBUGS and OpenBUGS from R
    / S-PLUS. Using this package, it is possible to call a BUGS model,
    summarize inferences and convergence in a table and graph, and save
    the simulations in arrays for easy access in R / S-PLUS. In S-PLUS,
    the openbugs functionality and the windows emulation functionality
    is not yet available.
-   `r pkg("rjags")`: Bayesian graphical models using MCMC. Interface to
    the JAGS MCMC library.



### Links
-   [gR initiative homepage and mailing list](http://www.R-project.org/gR/)
-   [Bioconductor](http://www.Bioconductor.org/)
-   [JAGS](http://mcmc-jags.sourceforge.net/)
-   [OpenBUGS](http://www.openbugs.info/w/)
-   [WinBUGS](http://www.mrc-bsu.cam.ac.uk/software/bugs/)
-   [Bayes Net Toolbox for MATLAB](https://code.google.com/p/bnt/)
-   [TETRAD](http://www.phil.cmu.edu/projects/tetrad/)
-   [Grappa](http://www.stats.bris.ac.uk/~peter/Grappa/)
-   [RHugin](http://rhugin.r-forge.r-project.org/)
