#必要なパッケージ
library(GBClust)
"""
事前に以下のインストールが必要
install.packages("devtools")
devtools::install_github("tommasorigon/GBClust") 
"""
library(mcclust)
library(ggplot2)
library(FactoMineR)
library(viridis) 
library(factoextra)
library(ca)
library(readxl)
library(tidyverse)
library(ggrepel)
library(dplyr)
library(tidytext)
library(purrr) 

#多重対応分析
df <- read.delim(file.choose(),
                 header = TRUE, sep = "\t",
                 stringsAsFactors = TRUE)

res.mca <- MCA(df, graph = FALSE)

#図1の図示
grp <- as.factor(df[, "X"])　#Xは調査対象となる形式の列の名称
fviz_mca_biplot(res.mca, repel = TRUE, col.var = "#4D4D4D",label=c("var"),#invisible='ind',
                habillage = grp, addEllipses = TRUE, ellipse.level = 0.9,
                labelsize = 3, pointsize = 3,arrows = c(FALSE, FALSE),geom.var = c("point", "text"))+labs(title = "MCA", x = "Dim.1", y ="Dim.2" )+theme_minimal()

#図2,3の図示
fviz_contrib(res.mca,  
             choice = "var",  
             axes = 1,  
             top = 5)

fviz_contrib(res.mca,  
             choice = "var",  
             axes = 2,  
             top = 5)　

#適切なクラスター数の確認
hcpc_res <- HCPC(re.ca, graph = FALSE)
fviz_cluster(hcpc_res,
             repel = FALSE,            
             show.clust.cent = TRUE, 
             palette = "jco",      
             ggtheme = theme_minimal(),
             main = "Factor map"
)

#図4の図示の準備
coords <- data.frame(
  Dim1 = res.mca$ind$coord[, 1],
  Dim2 = res.mca$ind$coord[, 2]
)

X <- model.matrix(~ . - 1, data = df)
K <- 4　#場合によって異なる

set.seed(123)
init    <- kmeans2(X, k = K, nstart = 20)
gibbs   <- kmeans_gibbs(X, k = K,
                        a_lambda = 0, b_lambda = 0,
                        R = 5000, burn_in = 3000, trace = TRUE)
S       <- mcclust::comp.psm(gibbs$G)
D       <- as.matrix(dist(X))^2
medoids <- comp_medoids(D, init$cluster)

Miscl <- function(S, cluster, medoids) {
  n <- nrow(S)
  sapply(seq_len(n),
         function(i) 1 - S[i, medoids[cluster[i]]])
}
miscl_probs <- Miscl(S, init$cluster, medoids)

results <- cbind(coords,
                 Cluster   = factor(init$cluster),
                 MisclProb = miscl_probs)
         
#図4の図示
results$Row <- seq_len(nrow(results))

results2 <- results %>%
  mutate(
    bin1 = round(Dim1, 2),
    bin2 = round(Dim2, 2)
  ) %>%
  group_by(bin1, bin2) %>%

  mutate(label_flag = (row_number() == 1)) %>%
  ungroup()

p1_one_label <- ggplot(results2, aes(x = Dim1, y = Dim2, color = Cluster)) +
  geom_point(size = 2) +
  geom_text_repel(
    data = filter(results2, label_flag),
    aes(label = Row),
    size         = 2,
    box.padding  = 0.3,
    point.padding= 0.5,
    segment.size = 0.2,
    segment.color= "gray50",
    max.overlaps = Inf
  ) +
  theme_minimal() +
  ggtitle("クラスタごとの散布図（重なり位置は1つだけラベル表示）")

print(p1_one_label)

#図5の図示
K <- 4   # クラスタ数
top_df <- lapply(seq_len(K), function(clu){
  tab2 <- as.data.frame(desc$category[[clu]])
  tab2$modality <- rownames(tab2)
  tab2$v.test   <- as.numeric(as.character(tab2$v.test))
  
  tab2 %>% 
    filter(v.test > 0) %>%            # 正の v.test のみ
    arrange(desc(v.test)) %>%         
    slice_head(n = 5) %>%             
    mutate(Cluster = paste0("Cluster ", clu))
}) |> bind_rows()

ggplot(top_df,
       aes(x = reorder_within(modality, v.test, Cluster),
           y = v.test,
           fill = Cluster)) +          
  geom_col(show.legend = TRUE) +
  coord_flip() +
  facet_wrap(~ Cluster, ncol = 2, scales = "free_y") +
  scale_x_reordered() +
  scale_color_viridis_c()+
  #scale_fill_brewer(palette = "Set2") + 
  labs(
    title = "各クラスタの寄与水準 Top 5（正の v.test のみ）",
    x     = "Variable = Level",
    y     = "v.test"
  ) +
  theme_minimal(base_size = 12)+theme(legend.position = 'none')

#図6の図示
p2_one_label <- ggplot(results2, aes(x = Dim1, y = Dim2, color = MisclProb)) +
  geom_point(size = 2) +
  geom_text_repel(
    data          = filter(results2, label_flag),
    aes(label     = Row),
    size          = 2,
    max.overlaps  = Inf,
    box.padding   = 0.3,
    point.padding = 0.5,
    segment.size  = 0.2,
    segment.color = "gray50"
  ) +
  scale_color_viridis_c() +
  theme_minimal() +
  ggtitle("Misclassification Probability (MCA space)\n重複位置グループから1つだけラベル表示")

print(p2_one_label)

#図7の図示
hist(miscl_probs,
     breaks = 30,
     main   = "Misclassification Probability Distribution",
     xlab   = "Misclassification Probability")

#個々の点の誤分類確率の確認
thr <- 0.01　#閾値を設定

high_uncert_df <- results2 %>% 
  filter(MisclProb >= thr, label_flag) %>%   
  arrange(Row) %>%                          
  select(Row, MisclProb)


high_uncert_list <- results2 %>% 
  filter(MisclProb >= thr) %>%               
  group_by(bin1, bin2) %>% 
  summarise(
    rows = list(Row),                      
    miscl = list(MisclProb),                
    .groups = "drop"
  ) %>% 
  mutate(pair = map2(rows, miscl, ~ tibble(Row = .x, MisclProb = .y))) %>% 
  pull(pair)

high_uncert_list
