plot_counts <- function(geneID, dds, annot,
                        logscale = TRUE,
                        xvar = "daylength", colorvar = "age") {
  
  description <- annot$Description[annot$geneID == geneID] %>%
    str_wrap(width = 50)
  
  d <- plotCounts(dds, gene=geneID,
                  intgroup=c("age", "daylength"), 
                  returnData=TRUE)
  
  p <- ggplot(d,
              aes(x=.data[[xvar]], y=count, color = .data[[colorvar]])) + 
    geom_point(position=position_jitter(w=0.1,h=0)) + 
    labs(title = geneID, subtitle = description) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  if (logscale == TRUE) p <- p + scale_y_log10()
  
  ggsave(paste0("results/DE/geneplots/", geneID, ".png"), p)
  
  return(p)
}
