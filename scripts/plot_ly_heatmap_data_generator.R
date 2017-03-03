# creates plot_ly heatmap for web app for GO enrichment heatmap

source('~/git/unified_gene_expression/figures and tables/GO_heatmap_eye_vs_body__megalong.R')

# manually cluster the rows
m <- wide_data[,2:ncol(wide_data)]
go_order <-hclust(dist(m))$order
# manually cluster the cols
sample_order <-hclust(dist(t(m)))$order

# now actually reorder
m <- m[go_order, sample_order]
# get colnames and rownames as ordered factors
order_columns <- factor(colnames(m), levels = colnames(m))
order_rows <-  factor(row.names(m), levels =row.names(m))


save(m, file='~/git/human_eyeIntegration_web_app/www/go_heatmap.Rdata')


        
# example plot_ly below
options <- list(
  l = 500,
  b = 300,
  tickangle = 90,
  type = 'category')
options2 <- list(
  l = 500,
  b = 300,
  type = 'category')
# example plot
plot_ly(z=as.matrix(m), 
                       y = order_rows,
                       x = order_columns,
                       type='heatmap', zmax=30, zauto=FALSE) %>% 
                 layout(margin=options, xaxis=options, yaxis=options2)
