---
title: "SignIT Vignette"
author: "Eric Zhao"
date: "September 20, 2017"
output: html_document
---



# SignIT: Mutation Signature Analysis in Individual Tumours with MCMC

Cancer is a disease characterized by the formation of tumours composed of cells which divide aggressively due to the loss of specialized functions and growth inhibition. Cancers arise due to the stepwise accumulation of spontaneous genetic mutations specific to the tumour, also known as somatic mutations. These mutations arise from many causes, including chemical mutagens, intracellular mutational processes, and deficient DNA repair genes. 

While some mutations may contribute directly to tumour invasiveness, the majority are thought to be "passengers," with minimal impact on the evolutionary fitness of cancer cells. Nevertheless, it is possible to glean insights about a cancer's molecular history from these passengers. They provide evidence of the mutational processes which drove the acquisition of mutations in cancer.

Recent research has described characteristic signatures of mutation deciphered from cancer whole genome and exome sequences.

This vignette demonstrates the use of SignIT to decipher mutation signatures in individual tumours and explore the solution space using Markov Chain Monte Carlo techniques.

# SignIT on a Patient's Mutation Catalog

## The Mutation Catalog

A mutation catalog is a Data Frame or Tibble containing two columns named mutation_type and count. The mutation_type column is a vector of unique strings each denoting a mutation type. It is possible to define custom mutation types, so long as they match the mutation types defined in a reference mutation signature table.


```r
data(example_catalog)
```

## Running MCMC Analysis and Exploring Posterior Distributions

Once a mutation catalog is available, obtaining the contributions of each mutation signature (also known as the exposures) is straightforward. The following line of code initializes the Bayesian model which analyzes mutation signatures and stores the output into a list structure.

By default, the analysis runs 4 chains in parallel, with 10 000 iterations per chain. However, for the purposes of this vignette, we will use only 2000 iterations per chain. This will help to lessen the runtime of this and later analyses for demonstration purposes.


```r
stan_model <- get_stan_model()
```

```
## Error in get_stan_model(): could not find function "get_stan_model"
```

```r
patient_signature_exposures <- get_stan_exposures(
  mutation_catalog = example_catalog,
  stan_model = stan_model
)
```

```
## Error in get_stan_exposures(mutation_catalog = example_catalog, stan_model = stan_model): could not find function "get_stan_exposures"
```

A few functions are available to explore the solutions obtained. The easiest way is to plot the posterior distributions for each mutation signature, as follows.


```r
patient_exposure_posterior_plot <- plot_exposure_posteriors(patient_signature_exposures)
```

```
## Error in eval(lhs, parent, parent): object 'patient_signature_exposures' not found
```

```r
print(patient_exposure_posterior_plot)
```

```
## Error in print(patient_exposure_posterior_plot): object 'patient_exposure_posterior_plot' not found
```

Note that all exposure plotting functions automatically shorten signature names by removing the word "Signature". For comparison, we also provide a convenience function which computes a naive non-negative least squares (NNLS) best fit solution.


```r
patient_exposure_nnls_plot <- plot_nnls_solution(patient_signature_exposures)
```

```
## Error in nnls(exposures_mcmc_output$reference_signatures %>% select(-mutation_type) %>% : could not find function "nnls"
```

```r
print(patient_exposure_nnls_plot)
```

```
## Error in print(patient_exposure_nnls_plot): object 'patient_exposure_nnls_plot' not found
```

***Write some observations here about the differences observed.***

***Write functions which provide quantitative summary statistics about the exposures MCMC output***.

# Density-based clustering of MCMC output

To mitigate overfitting mutation signatures in an individual tumour scenario, it is important to consider good alternative solutions which may or may not differ from the ideal point estimate determined by NNLS. Consider the following scenario. A tumour has a high contribution from Signature 3. However, Signature 3 and Signature 8 share some similarities in their mutation profiles. While the ideal solution may have a high Signature 3 and low Signature 8, an alternative solution with high Signature 8 and low Signature 3 may be nearly as good a fit.

In order to account for such solutions, we need to explore the density landscape of the 30-dimensional MCMC solution space for local optima. To do this, SignIT performs density-based clustering using the HDBSCAN algorithm. It then provides tools to explore the best-fit and highest-density regions within each solution cluster.

HDBSCAN accepts a parameter named minPts, which stands for minimum points. The higher this parameter, the fewer clusters will form and the more outlier points there will be. While SignIT uses a default minPts of 40, this parameter can be adjusted depending on the needs of your specific dataset.

The following command performs density clustering. Note that this step can place high demands on processor and memory resources, and may take a significant amount of time. For longer MCMC chains, it is recommended that this step is run on a high performance computing server.


```r
patient_density_clustering <- density_clustering(patient_signature_exposures, minPts = 40)
```

```
## Error in eval(lhs, parent, parent): object 'patient_signature_exposures' not found
```

```r
print(patient_density_clustering)
```

```
## Error in print(patient_density_clustering): object 'patient_density_clustering' not found
```

## Visualizing High-Dimensional Clusters with t-SNE

It is always handy to be visualize the high dimensional solution space to try and understand the structure of data and verify the accuracy of density clustering. Because the MCMC solutions exist in 30-dimensional space, however, it is necessary to dimensionally reduce the solutions down to 2 dimensions for visualization. The t-SNE algorithm is especially well suited to such a task. The following function produces a visualization of the exposures MCMC chain along two projected t-SNE dimensions. Points are coloured according to the cluster they were assigned to by density clustering in 30-dimensional space.

This visualization allows for verification that the data clusters along sensible borders and is not over-clustered or under-clustered.

Please note that this figure can take a long time to generate. We will generate three tsne graphs at different perplexities.


```r
plot_density_tsne(patient_signature_exposures, patient_density_clustering, perplexity = 1)
```

```
## Error in eval(lhs, parent, parent): object 'patient_signature_exposures' not found
```


```r
plot_density_tsne(patient_signature_exposures, patient_density_clustering, perplexity = 10)
```

```
## Error in eval(lhs, parent, parent): object 'patient_signature_exposures' not found
```


```r
plot_density_tsne(patient_signature_exposures, patient_density_clustering, perplexity = 50)
```

```
## Error in eval(lhs, parent, parent): object 'patient_signature_exposures' not found
```


```r
plot_density_tsne(patient_signature_exposures, patient_density_clustering, perplexity = 100)
```

```
## Error in eval(lhs, parent, parent): object 'patient_signature_exposures' not found
```


```r
plot_density_tsne(patient_signature_exposures, patient_density_clustering, perplexity = 200)
```

```
## Error in eval(lhs, parent, parent): object 'patient_signature_exposures' not found
```


```r
plot_density_tsne(patient_signature_exposures, patient_density_clustering, perplexity = 500)
```

```
## Error in eval(lhs, parent, parent): object 'patient_signature_exposures' not found
```


```r
plot_density_tsne(patient_signature_exposures, patient_density_clustering, perplexity = 1000)
```

```
## Error in eval(lhs, parent, parent): object 'patient_signature_exposures' not found
```


```r
plot_density_tsne(patient_signature_exposures, patient_density_clustering, perplexity = 10000)
```

```
## Error in eval(lhs, parent, parent): object 'patient_signature_exposures' not found
```


## Exploring Alternative Mutation Signature Solutions

SignIT provides a number of ways to explore optimal solutions within density clusters. The following method locates the solution within each cluster which minimizes the residual sum of squares, and thus provides the "best fit" to the data. Solutions are visualized as line graphs to show related changes in different signatures.


```r
best_fit_cluster_solutions_plot <- plot_cluster_solutions(patient_signature_exposures, patient_density_clustering, method = 'bestfit')
```

```
## Error in eval(lhs, parent, parent): object 'patient_signature_exposures' not found
```

```r
print(best_fit_cluster_solutions_plot)
```

```
## Error in print(best_fit_cluster_solutions_plot): object 'best_fit_cluster_solutions_plot' not found
```

Alternatively, the following solution identifies the "most dense" point per cluster by computing the mean distance to k neighbouring points. The k-value can be set via a parameter, and defaults to 400.


```r
most_dense_cluster_solutions_plot <- plot_cluster_solutions(patient_signature_exposures, patient_density_clustering, method = 'mostdense', k_value = 400)
```

```
## Error in eval(lhs, parent, parent): object 'patient_signature_exposures' not found
```

```r
print(most_dense_cluster_solutions_plot)
```

```
## Error in print(most_dense_cluster_solutions_plot): object 'most_dense_cluster_solutions_plot' not found
```

To concisely examine the clusters resulting from density clustering, the following function plots a scatterplot showing the minimum residual sum of squares versus the cluster score for each cluster. Better solutions should aim to minimize RSS, but more stable clusters should have a higher cluster score. Ideal solutions should therefore be in the bottom-right of this plot.

***Show naive NNLS solution as part of this plot***

***Replace the points with text labels showing the cluster number?***


```r
cluster_error_plot <- plot_cluster_solution_errors(patient_signature_exposures, patient_density_clustering)
```

```
## Error in eval(lhs, parent, parent): object 'patient_density_clustering' not found
```

```r
print(cluster_error_plot)
```

```
## Error in print(cluster_error_plot): object 'cluster_error_plot' not found
```

For comparison's sake, here are the various plots of solutions from above in a single grid.


```r
plot_grid(
  best_fit_cluster_solutions_plot + theme(legend.position = 'none') + labs(title = 'Solutions of each cluster (by lowest RSS)'),
  most_dense_cluster_solutions_plot + theme(legend.position = 'none') + labs(title = 'Solutions of each cluster (by highest density region)'),
  patient_exposure_nnls_plot + labs(title = 'Naive NNLS solution'),
  patient_exposure_posterior_plot + labs(title = 'Posteriors of the complete MCMC solution'),
  ncol = 1, align = 'v', 
  labels = c('A', 'B', 'C', 'D')
  )
```

```
## Error in plot_grid(best_fit_cluster_solutions_plot + theme(legend.position = "none") + : object 'best_fit_cluster_solutions_plot' not found
```

## Examining Signature Bleed

The previously described reciprocal interaction between Signature 3 and Signature 8 is referred to as "signature bleed". The similarity between these signatures results in instability of possible solutions. Bleed can be identified by an inverse correlation between the two signatures in the MCMC solution chain. The following method visualizes the projection of MCMC solutions onto the Signature 3 - Signature 8 plane.


```r
plot_two_signature_hexplot(patient_signature_exposures, 'Signature 3', 'Signature 8', patient_density_clustering, trendline = FALSE)
```

```
## Error in eval(lhs, parent, parent): object 'patient_signature_exposures' not found
```

By computing pairwise Spearman correlations for each signature pair, we can observe which signatures have anti-correlated exposure posterior distributions. The blue squares below indicate some degree of bleed between signatures.


```r
plot_signature_pairwise_bleed(patient_signature_exposures)
```

```
## Error in eval(lhs, parent, parent): object 'patient_signature_exposures' not found
```

***Create function to example signature bleed across the full structure of MCMC space, perhaps in a multivariate manner.***

# SignIT on a Simulated Mutation Catalog

In order to further illustrate the components of SignIT analysis, this next portion will simulate a catalog of somatic mutations.

In this vignette, we will demonstrate using the 30 reference mutation signatures defined on the [COSMIC](http://cancer.sanger.ac.uk/cosmic/signatures) website. The corresponding matrix of mutation signature probabilities can be downloaded [here](http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt). This dataset is built into SignIT and can be accessed via the get_reference_signatures method.


```r
reference_signatures <- get_reference_signatures()
print(reference_signatures)
```

```
## # A tibble: 96 x 31
##    mutation_type `Signature 1` `Signature 2` `Signature 3` `Signature 4`
##            <chr>         <dbl>         <dbl>         <dbl>         <dbl>
##  1       A[C>A]A   0.011098326  6.827082e-04   0.022172307        0.0365
##  2       A[C>A]C   0.009149341  6.191072e-04   0.017871675        0.0309
##  3       A[C>A]G   0.001490070  9.927896e-05   0.002138340        0.0183
##  4       A[C>A]T   0.006233885  3.238914e-04   0.016265146        0.0243
##  5       A[C>G]A   0.001801068  2.634810e-04   0.024002615        0.0097
##  6       A[C>G]C   0.002580909  2.698660e-04   0.012160304        0.0054
##  7       A[C>G]G   0.000592548  2.192339e-04   0.005275419        0.0031
##  8       A[C>G]T   0.002963986  6.109735e-04   0.023277656        0.0054
##  9       A[C>T]A   0.029514533  7.441557e-03   0.017872171        0.0120
## 10       A[C>T]C   0.014322747  2.726312e-03   0.008896034        0.0075
## # ... with 86 more rows, and 26 more variables: `Signature 5` <dbl>,
## #   `Signature 6` <dbl>, `Signature 7` <dbl>, `Signature 8` <dbl>,
## #   `Signature 9` <dbl>, `Signature 10` <dbl>, `Signature 11` <dbl>,
## #   `Signature 12` <dbl>, `Signature 13` <dbl>, `Signature 14` <dbl>,
## #   `Signature 15` <dbl>, `Signature 16` <dbl>, `Signature 17` <dbl>,
## #   `Signature 18` <dbl>, `Signature 19` <dbl>, `Signature 20` <dbl>,
## #   `Signature 21` <dbl>, `Signature 22` <dbl>, `Signature 23` <dbl>,
## #   `Signature 24` <dbl>, `Signature 25` <dbl>, `Signature 26` <dbl>,
## #   `Signature 27` <dbl>, `Signature 28` <dbl>, `Signature 29` <dbl>,
## #   `Signature 30` <dbl>
```

In order to simulate a mutation catalog, we need to multiply the reference signature matrix with a vector of signature exposures, which indicates the number of mutations contributed by each mutational process. For simplicity, we will simulate a catalog with only contributions from Signature 3.


```r
reference_signature_matrix <- reference_signatures_as_matrix(reference_signatures, mutation_catalog = reference_signatures)

simulated_exposure <- c(
  rep(0, 2),
  1000,
  rep(0, 27)
)

n_mutations <- sum(simulated_exposure)

catalog_probability_vector <- reference_signature_matrix %*% simulated_exposure / n_mutations
```

To perform the actual simulation, we draw mutations from a multinomial distribution.


```r
simulated_catalog <- tibble(
  mutation_type = reference_signatures[['mutation_type']],
  count = rmultinom(
    n = 1,
    size = n_mutations,
    prob = catalog_probability_vector
    ) %>% as.numeric
)

# Randomize the order to simulated catalog to test that functions are robust to reordering
simulated_catalog <- simulated_catalog[sample(1:dim(simulated_catalog)[1], replace = FALSE), ]
print(simulated_catalog)
```

```
## # A tibble: 96 x 2
##    mutation_type count
##            <chr> <dbl>
##  1       C[T>A]A     9
##  2       C[C>A]T    13
##  3       T[T>A]G     7
##  4       A[T>G]A     2
##  5       A[T>A]C     6
##  6       C[C>G]A    18
##  7       C[T>G]C     5
##  8       A[C>T]G     9
##  9       T[T>C]T     6
## 10       T[C>G]C    24
## # ... with 86 more rows
```

get_exposures will automatically account for the ordering of mutation types. Note that we randomized the order of our catalog rows to demonstrate this. Next, we compute signature exposures using SignIT.

Again, we use 2000 iterations instead of the default 10 000 in order to minimize runtimes.


```r
simulated_signature_exposures <- get_stan_exposures(
  mutation_catalog = simulated_catalog %>% mutate(count = as.integer(count)),
  n_chains = 4,
  n_iter = 10000
)
```

```
## Error in get_stan_exposures(mutation_catalog = simulated_catalog %>% mutate(count = as.integer(count)), : could not find function "get_stan_exposures"
```

```r
simulated_density_clustering <- density_clustering(simulated_signature_exposures, minPts = 40)
```

```
## Error in eval(lhs, parent, parent): object 'simulated_signature_exposures' not found
```

We proceed now to generate the plots demonstrated previously. Notice how even in this example with pure simulated Signature 3, the naive NNLS solution reports signature exposures arising from other mutation signatures. It is likely overfitting to the noise introduced by drawing from a multinomial random variable. The SignIT posteriors, however, provide a much more reasonable guess of the true exposure fractions, although they still do underestimate the exposure of Signature 3 (the true value was 1000). Note the slightly overestimated posteriors of other signatures including Signature 8, revealing some signature bleed.


```r
simulated_best_fit_plot <- plot_cluster_solutions(simulated_signature_exposures, simulated_density_clustering, method = 'bestfit')
```

```
## Error in eval(lhs, parent, parent): object 'simulated_signature_exposures' not found
```

```r
simulated_most_densest_region_plot <- plot_cluster_solutions(simulated_signature_exposures, simulated_density_clustering, method = 'mostdense', k_value = 400)
```

```
## Error in eval(lhs, parent, parent): object 'simulated_signature_exposures' not found
```

```r
simulated_exposure_nnls_plot <- plot_nnls_solution(simulated_signature_exposures)
```

```
## Error in nnls(exposures_mcmc_output$reference_signatures %>% select(-mutation_type) %>% : could not find function "nnls"
```

```r
simulated_exposure_posteriors <- plot_exposure_posteriors(simulated_signature_exposures, view = 'violin')
```

```
## Error in eval(lhs, parent, parent): object 'simulated_signature_exposures' not found
```

```r
plot_grid(
  simulated_best_fit_plot + theme(legend.position = 'none') + labs(title = 'Solutions of each cluster (by lowest RSS)'),
  simulated_most_densest_region_plot + theme(legend.position = 'none') + labs(title = 'Solutions of each cluster (by highest density region)'),
  simulated_exposure_nnls_plot + labs(title = 'Naive NNLS solution'),
  simulated_exposure_posteriors + labs(title = 'Posteriors of the complete MCMC solution'),
  ncol = 1, align = 'v', 
  labels = c('A', 'B', 'C', 'D')
  )
```

```
## Error in plot_grid(simulated_best_fit_plot + theme(legend.position = "none") + : object 'simulated_best_fit_plot' not found
```
