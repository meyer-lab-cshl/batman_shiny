#.libPaths("/opt/R/4.2.1/lib/R/library")
.libPaths(c("/usr/lib64/R/library","/usr/share/R/library"))

install.packages(c('igraph','readxl','tidyverse', 'plotly', 'ggalluvial'), repos='http://cran.us.r-project.org')
