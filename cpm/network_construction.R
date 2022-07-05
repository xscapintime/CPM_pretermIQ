data("neoOpen")

# TMFG network
tmfg <- TMFG(neoOpen)$A

# logo
logo <- LoGo(neoOpen)



#Load qgraph
library(qgraph)
#NEO-PI-3 defined facets
facets <- c(rep("actions", 8), rep("aesthetics", 8), rep("fantasy", 8),
rep("feelings", 8), rep("ideas", 8), rep("values", 8))
#Visualize TMFG
A <- qgraph(tmfg, groups = facets, palette = "ggplot2")
#Visualize LoGo
B <- qgraph(logo, groups = facets, palette = "ggplot2")


#Visualize TMFG and LoGo side-by-side
layout(t(1:2))
Layout <- averageLayout(A, B)
qgraph(A, layout = Layout, esize = 20, title = "TMFG")
qgraph(B, layout = Layout, esize = 20, title = "LoGo")
