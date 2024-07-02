options(tidyverse.quiet = TRUE)
library(tidyverse)
library(edgeR)
library(ggplot2)
library(scales)
library(dplyr)
library(viridis)
library(igraph)
library(hash)
library(mgcv)
library(gridExtra)
library(assertthat)
library(XML)
library(grid)
library(ggnewscale)
library(ragg)
library(truncnorm)
library(edgeR)
library(factoextra)
library(umap)
library(Seurat)
library(reshape2)

#sizes for the big/ web Illustrator panel
AXIS_TITLE = 18
AXIS_TEXT = 12
LEGEND_TITLE = 15
LEGEND_TEXT = 12
PLOT_TITLE = 20
PLOT_CAPTION = 20
theme_beauty <- function(x_axis_label_rotation = 90){ 
  font <- "Helvetica"   #font family
  theme_light() %+replace%  #generally light scheme
    theme(
      
      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      #axis.line = element_line(colour = "black", linetype='solid'),
      axis.ticks = element_line(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      
      #since theme_minimal() already strips axis lines, 
      #we don't need to do that again
      
      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = PLOT_TITLE,                #set font size
        #face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
        family = font,            #font family
        size = PLOT_CAPTION),               #font size
      
      plot.caption = element_text(           #caption
        family = font,            #font family
        size = PLOT_CAPTION,                 #font size
        hjust = 1),               #right align
      
      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = AXIS_TITLE),               #font size
      
      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        size = AXIS_TEXT),                #font size
      
      #margin dismensions: t,r,b,l
      #axis.text.x = element_text(            #margin for axis text:
      #  margin=margin(t=5, b = 10), angle = x_axis_label_rotation),         #t=distance of axis numbering from axis, b=distance of axis label from axis
      #axis.text.y = element_text(            #margin for axis text
      #  margin=margin(r=5,l = 10)),         #r=distance of axis numbering from axis, l=distance of axis label from axis
      
      #legend
      legend.text = element_text(colour="#000000", 
                                 size=LEGEND_TEXT),
      legend.title=element_text(size=LEGEND_TITLE)
    )
}

blue <- "#0138FE"
yellow <- "#FEC701"
black <- "#000000"
#purple <- "#990080"
purple <- "#660066"
grey <- "#A4A4A4"

#darker red and grey that fits to viridis scheme
red <- "#b22222"


calculate_pairwise_correlations<-function(data, order = NA, makeOrder = FALSE)
{
  wide_df <- data %>%
    dplyr::select(sample_id, ab_id, ab_count_normalized) %>%
    pivot_wider(names_from = ab_id, values_from = ab_count_normalized)  
  cor_matrix <- cor(wide_df[,-1], use = "complete.obs")
  
  #order by hieracical clustering
  if(makeOrder)
  {
    hc <- hclust(as.dist(1 - cor_matrix)) 
    ordered_cor_matrix <- cor_matrix[hc$order, hc$order]
    cor_matrix <- ordered_cor_matrix 
  }
  if(length(order) > 1)
  {
    orderList = order
  }
  else
  {
    prots <- unique(data[data$ab_type == "trajectoryProteins",]$ab_id)
    #order according to ab_type (trajectory proteins)
    orderFrame <- data %>%
      dplyr::select(c(ab_id, ab_type)) %>%
      unique()
    orderFrame <- orderFrame %>% arrange(ab_id) %>% arrange(ab_type)
    orderList <- orderFrame$ab_id
  }
  
  cor_matrix <- cor_matrix[orderList, orderList, drop=FALSE]
  
  melted_cor_matrix <- melt(cor_matrix)
  
  melted_cor_matrix_reordered <- melted_cor_matrix %>%
    dplyr::filter(Var1 %in% prots & Var2 %in% prots) %>%
    arrange(desc(Var1)) %>%
    bind_rows(melted_cor_matrix %>% 
    dplyr::filter(!(Var1 %in% prots & Var2 %in% prots)))
  
  
  #melted_cor_matrix$value <- melted_cor_matrix$value - mean(melted_cor_matrix$value)
  # Generate the heatmap
  plot <- (ggplot(data = melted_cor_matrix_reordered, aes(x = Var1, y = Var2, fill = value)) +
             geom_tile() +
             scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                                  midpoint = 0, limit = c(-1, 1), space = "Lab", 
                                  name="Correlation") +
             theme_minimal() + 
             theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                              size = 12, hjust = 1)) +
             coord_fixed())
  plot
  return(plot)
}

#data <- read_tsv("/Users/t.stohn/Desktop/Projects/scNormalization/PRESENTATIONS/BioSB_data/x/0_0/0/TMM.tsv") %>% group_by(sample_id) %>% mutate(libsize = sum(ab_count_normalized))
#mT <- calculate_pairwise_correlations(data)


dataG <- read_tsv("/Users/t.stohn/Desktop/Projects/scNormalization/CorrelationPreservationNorm/test/StrongLinearTrajectory/Groundtruth.tsv") %>% group_by(sample_id) %>% mutate(libsize = sum(ab_count_normalized)) %>%
  ungroup()
order <- unique(dataG$ab_id)


mG <- calculate_pairwise_correlations(dataG, order) + 
  theme_beauty() + theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                         axis.title.x = element_blank(), axis.title.y = element_blank()) 
mG

dataTmm <- read_tsv("/Users/t.stohn/Desktop/Projects/scNormalization/CorrelationPreservationNorm/test/StrongLinearTrajectory/TMM.tsv") %>% group_by(sample_id) %>% mutate(libsize = sum(ab_count_normalized)) %>%
  ungroup()
mT <- calculate_pairwise_correlations(dataTmm, order) + 
  theme_beauty() + theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                         axis.title.x = element_blank(), axis.title.y = element_blank()) 
mT

data <- read_tsv("/Users/t.stohn/Desktop/Projects/scNormalization/CorrelationPreservationNorm/test/corrPres.tsv") 
data <- data %>%pivot_longer(cols = -c("sample_id"),names_to = "ab_id", values_to = "ab_count_normalized")
data <- data %>% group_by(sample_id) %>% mutate(libsize = sum(ab_count_normalized)) %>% ungroup()
metaInfo <- dataTmm %>% dplyr::select(c(sample_id, ab_id, ab_type))
data <- inner_join(data, metaInfo, by = c("sample_id", "ab_id"))
mC <- calculate_pairwise_correlations(data, order)+ 
  theme_beauty() + theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                         axis.title.x = element_blank(), axis.title.y = element_blank()) 
mC

