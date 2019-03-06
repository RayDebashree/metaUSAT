
### Description
metaUSAT is a data-adaptive statistical approach for testing genetic associations of multiple traits from single/multiple studies using univariate GWAS summary statistics. This multivariate meta-analysis method can appropriately account for overlapping samples (if any) and can potentially test binary and/or continuous traits. The R function `metausat` implements this association test. For details of this statistical method, please refer/cite:

Ray, D., Boehnke, M. "[Methods for meta-analysis of multiple traits using GWAS summary statistics](http://onlinelibrary.wiley.com/doi/10.1002/gepi.22105/full)". *Genetic Epidemiology*, 42(2):134-145, 2018. [PMID: 29226385](https://www.ncbi.nlm.nih.gov/pubmed/29226385)

"The power of multivariate association test depends on a complex interplay of the number of truly associated traits, their correlation structure and the directions of the signals. The underlying association scenario changes from one genetic variant to another, and is not known *a priori* for any real dataset. There is no uniformly most powerful multivariate test, and a particular choice of association test may not be powerful enough to detect true signals. metaUSAT, being data-adaptive in nature, is less affected by the true (unknown) association scenario, and proves to be a robust yet
computationally efficient choice for investigators."

**Key Words:** Cross-phenotype association; GWAS; Meta-analysis; Multiple traits; Multivariate analysis;
Overlapping samples; PheWAS; Pleiotropy; Score test; Summary statistics

### Requirements
R (>= 3.0.1), CompQuadForm, minqa, psych, survey


### Changes
Version 1.17 - December 11, 2017
> First public release of the software.


### Notes
1. The method metaUSAT and its software is designed to test association of multiple traits (categorical and/or continuous) from a single study, or a single trait from multiple studies, or multiple traits from multiple studies. It only requires individual trait summary statistics for testing joint association of traits with a single genetic variant. 
    * Caution: metaUSAT is a single-variant association test; so not expected to work well for rare variants (i.e., genetic variants with very low allele-frequencies).

2. metaUSAT uses the summary statistics for a given genetic variant and the estimated correlation matrix to test association. If one or more studies have overlapping samples/individuals (which may or may not be known), the estimated correlation matrix reflects this overlap, and metaUSAT can appropriately account for that.
    * Caution: the joint analysis framework on which metaUSAT or metaMANOVA depends does not work well when the overlap between studies is large.

3. Since metaUSAT uses only summary statistics, it is assumed that all necessary covariate adjustments were performed when the individual trait summary statistics were obtained.

4. metaUSAT does not require independence of samples. When samples are related (e.g., in family-based GWAS), metaUSAT can use the summary statistics from [EMMAX](https://genome.sph.umich.edu/wiki/EMMAX) (or other univariate mixed model framework) to appropriately test for genetic associations.

5. metaUSAT does not assume homogeneity of trait effects across studies, hence powerful in detecting association if trait effects are heterogeneous across studies. However, if the studies are nearly independent and the trait effects are believed to be homogeneous across studies, you may use meta-analyzed summary statistics for each trait (e.g., Z-statistic output from [METAL](https://genome.sph.umich.edu/wiki/METAL_Documentation)) to perform joint meta-analysis of multiple traits. This may increase the power of the test by reducing the degrees of freedom of the test.

6. If you receive an error message like `the integral is probably divergent`, try reducing the absolute tolerance parameter `AbsTol`.

7. A user kindly pointed out a warning message displayed when `weights` are provided with version 1.17:<br/>
`In if (weights != 1) { :`<br/>
&nbsp;&nbsp;`the condition has length > 1 and only the first element will be used`<br/>
You can ignore this message (as long as the first weight is not 1, which would mean your sample size is 1); the weights you provided will be appropriately used inside the function. If the sample sizes are different across traits/studies, provide square root of sample sizes as `weights`, else the default is `weights=1`.
