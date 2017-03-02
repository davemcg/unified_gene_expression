# generate word cloud files for all vs all diff exp for web app

load('~/git/unified_gene_expression/data/go_enrichment_all_vs_all.Rdata') #all_vs_all
library(tm)
library(SnowballC)
library(wordcloud)
library(dplyr)

background_go_words <- function(go_enrichment, ontology){
  all_keywords <-go_enrichment %>% filter(Ontology=='BP') %>% dplyr::select(`GO ID`, Term) %>% unique() %>% .[['Term']]
  docs <- Corpus(VectorSource( wordStem(all_keywords, language = "english")))
  # Convert the text to lower case
  docs <- tm_map(docs, content_transformer(tolower))
  # Remove numbers
  docs <- tm_map(docs, removeNumbers)
  # Remove english common stopwords
  docs <- tm_map(docs, removeWords, stopwords("english"))
  # Remove your own stop word
  # specify your stopwords as a character vector
  # docs <- tm_map(docs, removeWords, c("blabla1", "blabla2")) 
  # Remove punctuations
  docs <- tm_map(docs, removePunctuation)
  # Eliminate extra white spaces
  docs <- tm_map(docs, stripWhitespace)
  
  dtm <- TermDocumentMatrix(docs)
  m <- as.matrix(dtm)
  v <- sort(rowSums(m),decreasing=TRUE)
  all_go <- data.frame(word = names(v),freq=v)
  #head(d, 10)
  all_go$ratio <- all_go$freq/sum(all_go$freq)
  all_go
}
all_BP_go <- background_go_words(all_vs_all_go, 'BP')
all_MF_go <- background_go_words(all_vs_all_go, 'MF')
#Word cloud function
cloud_maker <- function(word_vector, exclusion_words, max.words, all_go, title){
  docs <- Corpus(VectorSource(word_vector))
  docs <- tm_map(docs, content_transformer(tolower))
  docs <- tm_map(docs, removeNumbers)
  docs <- tm_map(docs, stripWhitespace)
  docs <- tm_map(docs, removePunctuation)
  docs <- tm_map(docs, removeWords, exclusion_words) 
  dtm <- TermDocumentMatrix(docs)
  m <- as.matrix(dtm)
  v <- sort(rowSums(m),decreasing=TRUE)
  d <- data.frame(word = names(v),freq=v)
  d$ratio <- d$freq/sum(d$freq)
  calculated <- left_join(d, all_go, by='word') %>% mutate(delta_ratio = ratio.x-ratio.y, transformed_delta =log10(delta_ratio+1)*1e6) %>% filter(transformed_delta>0)
  set.seed(1234)
  
  layout(matrix(c(1, 2), nrow=2), heights=c(1, 4))
  par(mar=rep(0, 4))
  
  png(filename = title, width = 400, height = 400, res = 300, pointsize = 3)
  
  #text(x=0.5, y=0.5, title)
  
  print(wordcloud(words = calculated$word, scale=c(3,1),freq = calculated$transformed_delta, min.freq = 1,
                  max.words=max.words, random.order=FALSE, rot.per=0.35, 
                  colors=c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3"),
                  main="Title"))
  calculated$word %>% data.frame()
  
  dev.off()
}



all_vs_all_go %>%
  group_by(Set, Test) %>% 
  filter(as.numeric(`P value (FDR)`) < 0.01, Ontology=='BP') %>% 
  mutate(Title = paste(Set, Test, sep='__'), Title=paste0(Title, '.png')) %>% 
  do(cloud_maker(.$Term,c('involved','process','regulation','negative','positive'),75,all_BP_go,.$Title))
