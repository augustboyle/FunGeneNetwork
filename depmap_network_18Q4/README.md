## Co-functional interaction visualization

This Shiny app visualizes Avana CRISPR screen data version 18Q4 from the Cancer Dependency map following Meyers, et al 2017 ([https://doi.org/10.1038/ng.3984](https://doi.org/10.1038/ng.3984)), as analyzed in the corresponding article in Molecular Systems Biology (Boyle, Pritchard and Greenleaf 2018).

Two visualizations are supported, viewed by clicking on the tabs above:


1. Network visualization
   - Genes are nodes
   - Co-functional interactions are edges
   - Other aesthetics can be mapped to color and size
2. Heatmap
   - Genes are rows
   - Cell lines are columns
   - Cell lines and genes are clustered by variance-normalized gene profiles
   - Other aesthetics can be added in an additional column

Layout options for the network view are taken from [igraph](http://igraph.org/r/doc/layout_.html)

'Mask by community' grays out nodes that belong to different de novo called gene communities in the network view 

'alpha' modifies the edge opacity in the network view.

'seed' allows different, reproducible layouts for the network view.

'min |R|' only draws edges between genes correlated above the threshold in the network view

[STRING](https://string-db.org/) version 10.5 is used to distinguish known (solid) from uncharacterized (dashed) interactions.

You may also upload a file to add your own scores ad hoc for visualization via the 'Upload gene scores' toolbar.
To do this, simply upload a file with 'gene_symbol' and 'score' columns. The 'score' column will be used when the 'Custom' option is selected for the gene score.

You can also add your own files to the script.

Example CRISPR screen data was taken from Jost, et al 20127 ([https://doi.org/10.1016/j.molcel.2017.09.012](https://doi.org/10.1016/j.molcel.2017.09.012))


Running this visualization locally requires the [R programming language](https://www.r-project.org/) and various packages.

After installing R, enter an R session and install the required packages
```
install.packages(c("tidyverse", "shiny", "visNetwork", "shinythemes", "circlize", "igraph"))

source("https://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap")
```
Now, you can run the app from a terminal
```
Rscript --vanilla app.R
```

When successful, R will tell you it is 'listening' at a given port
```
Listening on http://127.0.0.1:xxxx
```
Copy the address into your address bar in a browser like Chrome or Safari. Now you can begin visualizing co-functional interactions by query gene.

