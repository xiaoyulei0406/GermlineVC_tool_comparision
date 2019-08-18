library(dplyr)
library(ggplot2)
library(reshape2)
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
round_digits <- -2
files <- list.files(pattern = "summary\\.csv$")
print (files)
dlist <- lapply(files, read.csv)
names <- lapply(files, function(x) gsub("happy_", "", gsub(".summary.csv", "", x)))

dnamed <- mapply(cbind, dlist, "Name"=names, SIMPLIFY=F)
merged <- Reduce(function(...) merge(..., all=T), dnamed)
print(merged)
names(merged) <- c( "Type", "Filter", "Total", "True Positives", "False Negatives", "QTotal", "False Positives", "Unknown", "Genotype Error", "Recall", "Precision", "NA", "F1 Score", "T TiTv" , "Q TiTv" , "T Het Hom" , "Q Het Hom", "Name")
#names(merged) <- c( "Type", "Filter", "TRUTH.TOTAL", "TRUTH.TP","TRUTH.FN","QUERY.TOTAL", "QUERY.FP", "QUERY.UNK", "FP.gt", "METRIC.Recall", "METRIC.Precision", "METRIC.Frac_NA", "METRIC.F1_Score", "TRUTH.TOTAL.TiTv_ratio" , "QUERY.TOTAL.TiTv_ratio" , "TRUTH.TOTAL.het_hom_ratio" , "QUERY.TOTAL.het_hom_ratio","Name")

melted <- melt(merged, id.vars=c("Name", "Filter", "Type"))
metrics <- subset(melted, variable%in%c("Recall", "Precision", "F1 Score"))
p1 <- ggplot(metrics, aes(x=Name, y=value, color=Filter)) +
  geom_point(stat="identity", position = position_jitter(w = 0.06, h = 0), size =4) +
  geom_text(aes(label=ifelse(Filter=="PASS", round(value, 3), "")), color="black", size=2.5, hjust=-0.4, vjust=0.5) +
  geom_text(aes(label=ifelse(Filter!="PASS", round(value, 3), "")), color="darkgrey", size=2.5, hjust=1.6, vjust=0.5) +
  facet_grid( variable ~ Type, scales="free_y" ) +
  ylab("Metrics") +
  scale_color_manual(values=c("darkblue", "red")) +
  theme(axis.text.x=element_text(angle=30, hjust = 1))

ggsave(plot=multiplot(p1, cols=1), filename = 'metrics.png', width=10, height=8, units="in")
