## Correlation Analysis with SPIEC-EASI

#### Tutorial #### 
library(devtools)
# install_github("zdk123/SpiecEasi")
library(SpiecEasi)



data(amgut1.filt)
depths <- rowSums(amgut1.filt)
amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
amgut1.filt.cs <- round(amgut1.filt.n * min(depths))

d <- ncol(amgut1.filt.cs)
n <- nrow(amgut1.filt.cs)
e <- d

# Synthesize the data
set.seed(10010)
graph <- make_graph('cluster', d, e)
Prec  <- graph2prec(graph)
Cor   <- cov2cor(prec2cov(Prec))
X <- synth_comm_from_counts(amgut1.filt.cs, mar=2, distr='zinegbin', Sigma=Cor, n=n)

# the main SPIEC-EASI pipeline: Data transformation, sparse inverse covariance estimation and model selection
se <- spiec.easi(X, method='mb', lambda.min.ratio=1e-2, nlambda=15)

# examine ROC over lambda path and PR over the stars index for the selected graph
huge::huge.roc(se$est$path, graph, verbose=FALSE)
stars.pr(getOptMerge(se), graph, verbose=FALSE)
# stars selected final network under: se.est$refit$stars
# The above example does not cover all possible options and parameters. For example, other generative network models are available, the lambda.min.ratio (the scaling factor that determines the minimum sparsity/lambda parameter) shown here might not be right for your dataset, and its possible that you'll want more repetitions (number of subsamples) for StARS.




# Analysis of American Gut data
# Now let's apply SpiecEasi directly to the American Gut data. Don't forget that the normalization is performed internally in the spiec.easi function. Also, we should use a larger number of stars repetitions for real data. We can pass in arguments to the inner stars selection function as a list via the parameter pulsar.params. If you have more than one processor available, you can also supply a number to ncores. Also, let's compare results from the MB and glasso methods as well as SparCC (correlation).

se.mb.amgut <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-2,
                          nlambda=20, pulsar.params=list(rep.num=50))
se.gl.amgut <- spiec.easi(amgut1.filt, method='glasso', lambda.min.ratio=1e-2,
                          nlambda=20, pulsar.params=list(rep.num=50))
sparcc.amgut <- sparcc(amgut1.filt)
## Define arbitrary threshold for SparCC correlation matrix for the graph
sparcc.graph <- abs(sparcc.amgut$Cor) >= 0.3
diag(sparcc.graph) <- 0
sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
## Create igraph objects
ig.mb     <- adj2igraph(getRefit(se.mb.amgut))
ig.gl     <- adj2igraph(getRefit(se.gl.amgut))
ig.sparcc <- adj2igraph(sparcc.graph)
# Visualize using igraph plotting:
  
library(igraph)
## set size of vertex proportional to clr-mean
vsize    <- rowMeans(clr(amgut1.filt, 1))+6
am.coord <- layout.fruchterman.reingold(ig.mb)

par(mfrow=c(1,3))
plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")
plot(ig.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")


library(Matrix)
secor  <- cov2cor(getOptCov(se.gl.amgut))
sebeta <- symBeta(getOptBeta(se.mb.amgut), mode='maxabs')
elist.gl     <- summary(triu(secor*getRefit(se.gl.amgut), k=1))
elist.mb     <- summary(sebeta)
elist.sparcc <- summary(sparcc.graph*sparcc.amgut$Cor)

hist(elist.sparcc[,3], main='', xlab='edge weights')
hist(elist.mb[,3], add=TRUE, col='forestgreen')
hist(elist.gl[,3], add=TRUE, col='red')


# Lets look at the degree statistics from the networks inferred by each method.

dd.gl     <- degree.distribution(ig.gl)
dd.mb     <- degree.distribution(ig.mb)
dd.sparcc <- degree.distribution(ig.sparcc)

plot(0:(length(dd.sparcc)-1), dd.sparcc, ylim=c(0,.35), type='b',
     ylab="Frequency", xlab="Degree", main="Degree Distributions")
points(0:(length(dd.gl)-1), dd.gl, col="red" , type='b')
points(0:(length(dd.mb)-1), dd.mb, col="forestgreen", type='b')
legend("topright", c("MB", "glasso", "sparcc"),
       col=c("forestgreen", "red", "black"), pch=1, lty=1)


## Working with phyloseq

library(phyloseq)
## Load round 2 of American gut project
data('amgut2.filt.phy')
se.mb.amgut2 <- spiec.easi(amgut2.filt.phy, method='mb', lambda.min.ratio=1e-2,
                           nlambda=20, pulsar.params=list(rep.num=50))
ig2.mb <- adj2igraph(getRefit(se.mb.amgut2),  vertex.attr=list(name=taxa_names(amgut2.filt.phy)))
plot_network(ig2.mb, amgut2.filt.phy, type='taxa', color="Rank3")

# Cross domain interactions
data(hmp2)
se.hmp2 <- spiec.easi(list(hmp216S, hmp2prot), method='mb', nlambda=40,
                      lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05))

dtype <- c(rep(1,ntaxa(hmp216S)), rep(2,ntaxa(hmp2prot)))
plot(adj2igraph(getRefit(se.hmp2)), vertex.color=dtype+1, vertex.size=9)




######## Load Data & functions
source("src/load_phyloseq_obj.R")

dat.plot <- core(dat, detection = 0, prevalence = 0.1)
se.mb.amgut2 <- spiec.easi(dat.plot, method='mb', lambda.min.ratio=1e-2,
                           nlambda=20, pulsar.params=list(rep.num=50))
ig2.mb <- adj2igraph(getRefit(se.mb.amgut2),  vertex.attr=list(name=taxa_names(dat.plot)))
plot_network(ig2.mb, dat.plot, type='taxa', color="Phylum")


# dtype <- c(rep(1,ntaxa(hmp216S)), rep(2,ntaxa(hmp2prot)))
# plot(adj2igraph(getRefit(se.hmp2)), vertex.color=dtype+1, vertex.size=9)




