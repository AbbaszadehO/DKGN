# Data-driven and Knowledge-based Algorithms for Gene Network Reconstruction on High-dimensional Data
An R implementation (dev version)

# PREREQUISITES
* [R (3.6.x)](https://cran.r-project.org/bin/windows/base/) (R 3.6.1 is recommended).
* [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) (>=1.0.2).
* [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/RcppArmadillo.pdf) (>=0.9.800.1.0).
* [fdrtool](https://cran.r-project.org/web/packages/fdrtool/index.html) (>=1.2.15).

# Installing
This version can be directly installed from the R console using:

```
# install.packages("devtools")
devtools::install_github("AbbaszadehO/DKGN")
```

# Running the examples
## Scenario 1 ##

Data-driven algorithm on Network1 

```
library(DKGN)
TFnum = 195
n = nrow(Network1)
p = ncol(Network1)
S = t(Network1) %*% Network1 / n
target = getTarget(S = S, type = "mean")
sigma_hat = shrinkCovariance(S, target = target, n = n, lambda = seq(0.01, 0.99, 0.01))
gamma = getGammamatrix(sigma_hat,confidence=0.95)
omega_hat = sparsePrecision(
  sigma_hat,
  n,
  TFnum,
  gamma = gamma,
  rho = 1.0,
  max_iter = 100,
  tol = 1e-10,
)
partialcorr = pCorr(precision = omega_hat)
colnames(partialcorr) = rownames(partialcorr) = rownames(S)
network = data.frame(as.table(partialcorr[1:TFnum, (TFnum + 1):p]))
colnames(network) = c("node1", "node2", "partialcorr")
network = network[order(abs(network$partialcorr), decreasing = TRUE), ]
write.table(
  x = network,
  file = "net1.txt",
  sep = "\t",
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)
```

Knowledge-based (PM+TZ) algorithm on Network1 when 25% of true edges are available

```
library(DKGN)
set.seed(100)
percentage = 0.25
TFnum = 195
n = nrow(Network1)
p = ncol(Network1)
S = t(Network1) %*% Network1 / n
prior =
  Network1_gold[sample(nrow(Network1_gold), floor(percentage * dim(Network1_gold)[1])), c(1, 2)]
target = getTarget(S = S, type = "prior", prior = prior)
sigma_hat = shrinkCovariance(S,
                             target = target,
                             n = n,
                             lambda = seq(0.01, 0.99, 0.01))
diag(target) = target[target == 0] = 1
target[target!=1] = 0
gamma = getGammamatrix(sigma_hat, confidence = 0.95, prior = target)
omega_hat = sparsePrecision(
  sigma_hat,
  n,
  TFnum,
  gamma,
  rho = 1.0,
  max_iter = 100,
  tol = 1e-10,
)
partialcorr = pCorr(precision = omega_hat)
colnames(partialcorr) = rownames(partialcorr) = rownames(S)
network = data.frame(as.table(partialcorr[1:TFnum, (TFnum + 1):p]))
colnames(network) = c("node1", "node2", "partialcorr")
network = network[order(abs(network$partialcorr), decreasing = TRUE),]
write.table(
  x = network,
  file = "net1+PM+TZ.txt",
  sep = "\t",
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)
```

## Scenario 2 ##

Data-driven algorithm on Network4

```
library(DKGN)
library(fdrtool)
n = nrow(Network4)
p = ncol(Network4)
S = t(Network4) %*% Network4 / n
target = getTarget(S = S, type = "mean")
sigma_hat = shrinkCovariance(S,
                             target = target,
                             n = n,
                             lambda = seq(0.01, 0.99, 0.01))
gamma = getGammamatrix(sigma_hat, confidence=0.95)
omega_hat = sparsePrecision(
  sigma_hat,
  n,
  p,
  gamma = gamma,
  rho = 1.0,
  max_iter = 200,
  tol = 1e-10,
)
partialcorr = pCorr(precision = omega_hat)
colnames(partialcorr) = rownames(partialcorr) = rownames(S)
partialcorr[lower.tri(partialcorr, diag = TRUE)] = NA
network = na.omit(data.frame(as.table(partialcorr)))
colnames(network) = c("node1", "node2", "partialcorr")
fdrout = fdrtool(network$partialcorr,statistic = "correlation",plot =FALSE)
network$q_values = fdrout$qval
network = network[network$q_values<=0.01,]
```

Knowledge-based (PM+TZ) algorithm on Network4 when 25% of true edges are available

```
library(DKGN)
library(fdrtool)
set.seed(100)
percentage = 0.25
n = nrow(Network4)
p = ncol(Network4)
S = t(Network4) %*% Network4 / n
prior =
  Network4_gold[sample(nrow(Network4_gold), floor(percentage * dim(Network4_gold)[1])), c(1, 2)]
target = getTarget(S = S, type = "prior", prior)
sigma_hat = shrinkCovariance(S,
                             target = target,
                             n = n,
                             lambda = seq(0.01, 0.99, 0.01))
diag(target) = target[target == 0] = 1
target[target != 1] = 0
gamma = getGammamatrix(sigma_hat, confidence = 0.95, prior = target)
omega_hat = sparsePrecision(
  sigma_hat,
  n,
  p,
  gamma = gamma,
  rho = 1.0,
  max_iter = 200,
  tol = 1e-10,
)
partialcorr = pCorr(precision = omega_hat)
colnames(partialcorr) = rownames(partialcorr) = rownames(S)
partialcorr[lower.tri(partialcorr, diag = TRUE)] = NA
network = na.omit(data.frame(as.table(partialcorr)))
colnames(network) = c("node1", "node2", "partialcorr")
fdrout = fdrtool(network$partialcorr,statistic = "correlation",plot =FALSE)
network$q_values = fdrout$qval
network = network[network$q_values<=0.01,]
```


## Scenario 3 ##

Knowledge-based (PM+Z) algorithm on Network4 when 60% of edges in the prior are true and 40% false


```
set.seed(100)
data = Network4
n = nrow(data)
p = ncol(data)
S = t(data) %*% data / n
x = 0.6
genenames = colnames(Network4)
rownames(S) = colnames(S) = genenames
string = data.frame(Network4_gold)
gold <- data.frame(string$X1,string$X2,1)
gold_matrix <- matrix(data = 0,nrow = p, ncol = p)
rownames(gold_matrix) <- genenames
colnames(gold_matrix) <- genenames
for(i in 1:length(string$X1)){
  tmp <- gold[i,]
  tmp<-as.matrix(tmp)
  gold_matrix[tmp[1,1],tmp[1,2]] <- 1
  gold_matrix[tmp[1,2],tmp[1,1]] <- 1
}
diag(gold_matrix) <- 0
graph = gold_matrix
true.mat = graph
p.mat <- true.mat
p.mat[] <- 1
t.edge   <- which(true.mat == 1)
t.unedge <- which(true.mat == 0)
p.edge <- sample(t.edge, size=length(t.edge) * x)
p.edge <- sort(c(p.edge, sample(t.unedge,
                                size=length(t.edge) * (1-x))))
p.mat[p.edge] <- runif(length(p.edge), min=0, max=(1-x))
p.mat[p.mat==1] <- runif(length(t.unedge), min=x, max=1)
p.mat<-as.matrix(forceSymmetric(p.mat))
diag(p.mat) <- 1
prior_knowledge = p.mat
target = getTarget(S = S, type = "mean")
sigma_hat = shrinkCovariance(S, target = target, n = n,
                             lambda = seq(0.01, 0.99, 0.01))
gamma = getGammamatrix(S = sigma_hat, confidence = 0.95,
                       prior = prior_knowledge)
omega_hat = sparsePrecision(
  sigma_hat,
  n,
  p,
  gamma = gamma,
  rho = 1,
  max_iter = 200,
  tol = 1e-10
)
partialcorr = pCorr(precision = omega_hat)
colnames(partialcorr) = rownames(partialcorr) = rownames(S)
partialcorr[lower.tri(partialcorr, diag = TRUE)] = NA
network = na.omit(data.frame(as.table(partialcorr)))
colnames(network) = c("node1", "node2", "partialcorr")
fdrout = fdrtool(network$partialcorr,statistic = "correlation",plot =FALSE)
network$q_values = fdrout$qval
network = network[network$q_values <=0.01,]
```


# Authors

* **Omid Abbaszadeh** - - [Scholar](https://scholar.google.com/citations?user=6J2pSg4AAAAJ&hl=en)
* **Ali Azarpeyvand** - - [Scholar](https://scholar.google.com/citations?user=XK7HXzMAAAAJ&hl=en)
* **Alireza Khanteymoori** - - [Scholar](https://scholar.google.com/citations?user=VAyELgcAAAAJ&hl=en)
* **Abbas Bahari** - - [Scholar](https://scholar.google.com/citations?user=9yNynDoAAAAJ&hl=en)

# Reference

Reference will be added after the article is published

# License

This project is licensed under the GNU GENERAL PUBLIC LICENSE - see the [LICENSE](LICENSE) file for details
