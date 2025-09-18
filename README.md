### *SpotSpace*: Methods for Projecting Network Signals in Spatial Transcriptomics
  <!-- badges: start -->
  [![](https://img.shields.io/badge/lifecycle-maturing-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#maturing)
  [![](https://img.shields.io/badge/license-Artistic--2.0-blue.svg)](https://cran.r-project.org/web/licenses/Artistic-2.0)
  [![](https://img.shields.io/badge/doi-10.32614/CRAN.package.PathwaySpace-blue.svg)](https://doi.org/10.32614/CRAN.package.PathwaySpace)
  <!-- badges: end -->
The *SpotSpace* package extends *PathwaySpace* methods alongside *Seurat* workflows, providing tools for signal propagation and visualization in spatial transcriptomics. By integrating signal processing with spatial visualization, *SpotSpace* allows users to project network signals onto spot-level coordinates to explore signal patterns on tissue microenvironments.
### Installation in R (>=4.4)

##### Install dependencies to build the package's vignettes

```r
install.packages("knitr")
install.packages("rmarkdown")
install.packages("Seurat")
install.packages("SeuratObject")
install.packages("patchwork")
install.packages("remotes")
remotes::install_github("satijalab/seurat-data")
```

##### Install the SpotSpace package

```r
install.packages("remotes")
remotes::install_github("sysbiolab/RGraphSpace", build_vignettes=TRUE)
remotes::install_github("sysbiolab/PathwaySpace", build_vignettes=TRUE)
remotes::install_github("sysbiolab/SpotSpace", build_vignettes=TRUE)
```

### Examples

Follow the *SpotSpace* vignette and try to make some *plots*!

```r
library(SpotSpace)
vignette("SpotSpace")
```

### Citation

If you use *SpotSpace*, please cite:

* Tercan & Apolonio *et al.* Protocol for assessing distances in pathway space for classifier feature sets from machine learning methods. *STAR Protocols*, 2025. https://doi.org/10.1016/j.xpro.2025.103681

### Licenses

The *SpotSpace* package is distributed under [Artistic-2.0](https://www.r-project.org/Licenses/Artistic-2.0)
