# ==================== 1. 加载必要包 ====================
library(tidyverse)
library(ggsci)
library(ggrepel)
library(forcats)
library(tidygraph)
library(ggraph)
library(igraph)
dir()
# ==================== 2. 读取数据 ====================
diff_kegg  <- read_csv("KEGG_diff_results_IT_vs_ST.csv",
                       col_types = cols(
                         KO           = col_character(),
                         genes        = col_character(),
                         logFC        = col_double(),
                         logCPM       = col_double(),
                         F            = col_double(),
                         PValue       = col_double(),
                         FDR          = col_double(),
                         Level3       = col_character(),
                         Significance = col_character()
                       ))

diff_kegg  <- read_csv("KEGG_diff_results_IT-H_vs_ST-H.csv",
                       col_types = cols(
                         KO           = col_character(),
                         genes        = col_character(),
                         logFC        = col_double(),
                         logCPM       = col_double(),
                         F            = col_double(),
                         PValue       = col_double(),
                         FDR          = col_double(),
                         Level3       = col_character(),
                         Significance = col_character()
                       ))

kegg_paths <- read_csv("KEGG_paths.csv",
                       col_types = cols(
                         KO    = col_character(),
                         paths = col_character()
                       ))

# ==================== 3. 统计背景 KO–Level3 对应数 ====================
bg_counts <- kegg_paths %>%
  separate_rows(paths, sep = "\\s*\\|\\s*") %>%
  separate(paths,
           into = c("Level1", "Level2", "Level3"),
           sep    = ";\\s*",
           fill   = "right",
           extra  = "merge") %>%
  filter(!is.na(Level3), Level3 != "") %>%
  distinct(KO, Level3) %>%
  count(Level3, name = "bg_count")

# ==================== 4. 挑选显著差异 KO 并展开 Level3 ====================
sig_kegg <- diff_kegg %>%
  filter(PValue < 0.05, abs(logFC) > 1) %>%
  separate_rows(Level3, sep = "\\s*\\|\\s*") %>%
  mutate(Level3 = str_trim(Level3)) %>%
  filter(Level3 != "") %>%
  distinct(KO, Level3, .keep_all = TRUE)

# ==================== 5. 汇总富集信息 ====================
enriched <- sig_kegg %>%
  group_by(Level3) %>%
  summarise(
    n_sig     = n(),
    avg_logFC = mean(logFC),
    min_FDR   = min(FDR),
    .groups   = "drop"
  ) %>%
  left_join(bg_counts, by = "Level3") %>%
  mutate(
    rich_factor = n_sig / bg_count,
    Group       = if_else(avg_logFC > 0, "IT-enriched", "ST-enriched")
  )

# ==================== 6. 添加 hub KO 加权 ====================
# 从网络分析中获得的枢纽 KO 列表
hub_KO <- c("K04498", "K08738", "K03283", "K03174", "K06237",
            "K00626", "K01593", "K00274", "K01692", "K00128",
            "K05692", "K00718", "K02790", "K02791", "K18933",
            "K00461", "K14004", "K02580", "K14754", "K02660",
            "K03531", "K00850", "K14260")
# 提取包含枢纽 KO 的通路名
pathways_with_hubKO <- sig_kegg %>%
  filter(KO %in% hub_KO) %>%
  pull(Level3) %>%
  unique()
# 计算加权分数：rich_factor + 20% hub_flag
enriched <- enriched %>%
  mutate(
    hub_flag = if_else(Level3 %in% pathways_with_hubKO, 1, 0),
    score    = rich_factor * (1 + 0 * hub_flag)
  )



# ==================== 7. （可选）二级功能标签 ====================
enriched2 <- enriched %>%
  mutate(
    Function2 = case_when(
      str_detect(Level3, regex("cancer|infection|diabetes|myocarditis|leukemia|hepatocellular|cardiomyopathy",
                               ignore_case = TRUE)) ~ "Other",
      str_detect(Level3, regex("metabolism|biosynthesis|degradation|pathway",
                               ignore_case = TRUE)) ~ "Metabolism",
      str_detect(Level3, regex("signaling|signal|kinase|receptor",
                               ignore_case = TRUE)) ~ "Signaling",
      str_detect(Level3, regex("cell cycle|transport|trafficking|phagosome",
                               ignore_case = TRUE)) ~ "Cellular_process",
      str_detect(Level3, regex("toxin|resistance|antibiotic|macrolide",
                               ignore_case = TRUE)) ~ "Resistance",
      TRUE ~ "Other"
    )
  )

# ==================== 8. 各组复合指标排序并取前20 ====================
top_enriched <- enriched2 %>%
  filter(n_sig >= 3) %>%                          # 保证至少3个差异 KO
  group_by(Group) %>%
  slice_max(order_by = score, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Level3 = fct_reorder(Level3, score))

# ==================== 9. 绘制复合指标气泡图 ====================
bubble_plot <- top_enriched %>%
  ggplot(aes(
    x     = rich_factor,
    y     = Level3,
    size  = n_sig,
    color = -log10(min_FDR)
  )) +
  geom_point(alpha = 0.8) +
  scale_color_gradientn(
    colours = pal_nejm()(7),
    name    = "-log10(FDR)"
  ) +
  scale_size_area(
    max_size = 8,
    name     = "n_sig"
  ) +
  labs(
    x        = "Rich factor (n_sig / background KO count)",
    y        = NULL,
    title    = "KEGG Pathway Enrichment Bubble Plot",
    subtitle = "FDR < 0.05 & |log2FC| > 1"
  ) +
  facet_wrap(~Group, ncol = 1, scales = "free_y") +
  theme_bw(base_size = 13, base_family = "Helvetica") +
  theme(
    strip.text       = element_text(face = "bold", size = 12),
    axis.text.y      = element_text(size = 8),
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    plot.margin      = margin(5,5,5,5)
  ) +
  geom_text_repel(
    aes(label = sprintf("%.2f", rich_factor)),
    size        = 2.5,
    show.legend = FALSE,
    seed        = 42
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, max(top_enriched$rich_factor) * 1.1)
  )

# 保存图形
ggsave("KEGG_Enrichment_BubblePlot_weighted-1.pdf",
       bubble_plot,
       width = 7, height = 8,
       device = cairo_pdf)

print(bubble_plot)








# ==================== 1. 加载必要包 ====================
library(tidyverse)
library(ggsci)
library(ggrepel)
library(forcats)
library(tidygraph)
library(ggraph)
library(igraph)

# ==================== 2. 读取数据 ====================
diff_kegg  <- read_csv("KEGG_diff_results_IT_vs_ST.csv",
                       col_types = cols(
                         KO           = col_character(),
                         genes        = col_character(),
                         logFC        = col_double(),
                         logCPM       = col_double(),
                         F            = col_double(),
                         PValue       = col_double(),
                         FDR          = col_double(),
                         Level3       = col_character(),
                         Significance = col_character()
                       ))
kegg_paths <- read_csv("KEGG_paths.csv",
                       col_types = cols(
                         KO    = col_character(),
                         paths = col_character()
                       ))

# ==================== 3. 统计背景 KO–Level3 对应数 ====================
bg_counts <- kegg_paths %>%
  # 拆分多重通路
  separate_rows(paths, sep = "\\s*\\|\\s*") %>%
  mutate(paths = str_trim(paths)) %>%
  # 拆三级注释
  separate(paths,
           into = c("Level1", "Level2", "Level3"),
           sep    = ";\\s*",
           fill   = "right",
           extra  = "merge") %>%
  filter(!is.na(Level3), Level3 != "") %>%
  distinct(KO, Level3) %>%
  count(Level3, name = "bg_count")

# ==================== 4. 挑选显著差异 KO 并展开 Level3 ====================
sig_kegg <- diff_kegg %>%
  filter(FDR < 0.05, abs(logFC) > 1) %>%
  separate_rows(Level3, sep = "\\s*\\|\\s*") %>%
  mutate(Level3 = str_trim(Level3)) %>%
  filter(Level3 != "") %>%
  distinct(KO, Level3, .keep_all = TRUE)

# ==================== 5. 汇总富集信息 ====================
enriched <- sig_kegg %>%
  group_by(Level3) %>%
  summarise(
    n_sig     = n(),
    avg_logFC = mean(logFC),
    min_FDR   = min(FDR),
    .groups   = "drop"
  ) %>%
  left_join(bg_counts, by = "Level3") %>%
  mutate(
    rich_factor = n_sig / bg_count,
    Group       = if_else(avg_logFC > 0, "IT-enriched", "ST-enriched")
  ) %>%
  arrange(Group, min_FDR)

# ==================== 6. （可选）二级功能标签 ====================
enriched2 <- enriched %>%
  mutate(
    Function2 = case_when(
      str_detect(Level3, regex("cancer|infection|diabetes|myocarditis|leukemia|hepatocellular|cardiomyopathy",
                               ignore_case = TRUE)) ~ "Other",
      str_detect(Level3, regex("metabolism|biosynthesis|degradation|pathway",
                               ignore_case = TRUE)) ~ "Metabolism",
      str_detect(Level3, regex("signaling|signal|kinase|receptor",
                               ignore_case = TRUE)) ~ "Signaling",
      str_detect(Level3, regex("cell cycle|transport|trafficking|phagosome",
                               ignore_case = TRUE)) ~ "Cellular_process",
      str_detect(Level3, regex("toxin|resistance|antibiotic|macrolide",
                               ignore_case = TRUE)) ~ "Resistance",
      TRUE ~ "Other"
    )
  ) %>%
  arrange(Group, Function2, min_FDR)

# ==================== 7. 各组取前 20 通路（严格不超过20） ====================
top_enriched <- enriched2 %>%
  group_by(Group) %>%
  slice_min(order_by = min_FDR, n = 20, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Level3 = fct_reorder(Level3, rich_factor))

# ==================== 8. 绘制气泡图 ====================
bubble_plot <- top_enriched %>%
  ggplot(aes(x = rich_factor,
             y = Level3,
             size = n_sig,
             color = -log10(min_FDR))) +
  geom_point(alpha = 0.8) +
  scale_color_gradientn(
    colours = pal_nejm()(7),
    name    = "-log10(FDR)"
  ) +
  scale_size_area(
    max_size = 8,
    name     = "n_sig"
  ) +
  labs(
    x        = "Rich factor (n_sig / background KO count)",
    y        = NULL,
    title    = "KEGG Pathway Enrichment Bubble Plot",
    subtitle = "FDR < 0.05 & |log2FC| > 1"
  ) +
  facet_wrap(~Group, ncol = 1, scales = "free_y") +
  theme_bw(base_size = 13, base_family = "Helvetica") +
  theme(
    strip.text        = element_text(face = "bold", size = 12),
    axis.text.y       = element_text(size = 8),
    legend.position   = "right",
    panel.grid.minor  = element_blank(),
    plot.margin       = margin(5,5,5,5)
  ) +
  geom_text_repel(
    aes(label = sprintf("%.2f", rich_factor)),
    size       = 2.5,
    show.legend = FALSE,
    seed       = 42
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, max(top_enriched$rich_factor) * 1.1)
  )

# 直接保存即可
ggsave("KEGG_Enrichment_BubblePlot.pdf",
       bubble_plot,
       width = 7, height = 8,
       device = cairo_pdf)

# ==================== 9. 构建 IT-enriched 二部网络 ====================
it_edges <- sig_kegg %>%
  filter(logFC > 0) %>%
  distinct(KO, Level3) %>%
  mutate(KO = paste0("KO:", KO))

nodes <- tibble(
  name = c(unique(it_edges$KO), unique(it_edges$Level3)),
  type = c(
    rep("KO", length(unique(it_edges$KO))),
    rep("Pathway", length(unique(it_edges$Level3)))
  )
)
edges <- it_edges %>% rename(from = KO, to = Level3)

graph_it <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE) %>%
  activate(nodes) %>%
  mutate(degree = centrality_degree())

# ==================== 10. 导出并绘制子网络 (Top10 KO) ====================
top_kos <- graph_it %>%
  as_tibble() %>%
  filter(type == "KO") %>%
  slice_max(order_by = degree, n = 10, with_ties = FALSE) %>%
  pull(name)

subgraph <- graph_it %>%
  activate(nodes) %>%
  filter(
    (type == "KO"      & name %in% top_kos) |
      (type == "Pathway" & name %in% edges$to[edges$from %in% top_kos])
  )

# 可视化子网
ggraph(subgraph, layout = "fr") +
  geom_edge_link(alpha = 0.2) +
  geom_node_point(aes(
    size  = if_else(type == "KO", degree, 3),
    color = type
  )) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c(KO = "#D73027", Pathway = "#4575B4")) +
  theme_void()

# 导出 GML
write_graph(as.igraph(subgraph),
            file   = "IT_KO_Pathway_Subnetwork_top10.gml",
            format = "gml")
write_graph(as.igraph(graph_it),
            file   = "IT_KO_Pathway_FullNetwork.gml",
            format = "gml")

#############################################################
it_edges <- sig_kegg %>% filter(logFC > 0) %>% distinct(KO, Level3)
# 添加前缀标识
it_edges <- it_edges %>% mutate(KO = paste0("IT_KO:", KO))

nodes_it <- tibble(
  name = c(unique(it_edges$KO), unique(it_edges$Level3)),
  type = c(rep("KO", length(unique(it_edges$KO))), rep("Pathway", length(unique(it_edges$Level3))))
)
edges_it <- it_edges %>% rename(from = KO, to = Level3)

graph_it <- tbl_graph(nodes = nodes_it, edges = edges_it, directed = FALSE) %>%
  activate(nodes) %>% mutate(degree = centrality_degree())

# 导出并绘制 IT 子网络 (Top10)
top_kos_it <- graph_it %>% as_tibble() %>% filter(type == "KO") %>% slice_max(order_by = degree, n = 10) %>% pull(name)
sub_it <- graph_it %>% activate(nodes) %>% filter((type == "KO" & name %in% top_kos_it) | (type == "Pathway" & name %in% edges_it$to[edges_it$from %in% top_kos_it]))
ggraph(sub_it, layout = "fr") + geom_edge_link(alpha = 0.2) +
  geom_node_point(aes(size = if_else(type == "KO", degree, 3), color = type)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c(KO = "#D73027", Pathway = "#4575B4")) + theme_void()
write_graph(as.igraph(sub_it), file = "IT_KO_Pathway_Subnetwork_top10.gml", format = "gml")





#############################################################
st_edges <- sig_kegg %>% filter(logFC < 0) %>% distinct(KO, Level3)
# 添加前缀标识
st_edges <- st_edges %>% mutate(KO = paste0("ST_KO:", KO))

nodes_st <- tibble(
  name = c(unique(st_edges$KO), unique(st_edges$Level3)),
  type = c(rep("KO", length(unique(st_edges$KO))), rep("Pathway", length(unique(st_edges$Level3))))
)
edges_st <- st_edges %>% rename(from = KO, to = Level3)

graph_st <- tbl_graph(nodes = nodes_st, edges = edges_st, directed = FALSE) %>%
  activate(nodes) %>% mutate(degree = centrality_degree())

# 导出并绘制 ST 子网络 (Top10)
top_kos_st <- graph_st %>% as_tibble() %>% filter(type == "KO") %>% slice_max(order_by = degree, n = 10) %>% pull(name)
sub_st <- graph_st %>% activate(nodes) %>% filter((type == "KO" & name %in% top_kos_st) | (type == "Pathway" & name %in% edges_st$to[edges_st$from %in% top_kos_st]))
ggraph(sub_st, layout = "fr") + geom_edge_link(alpha = 0.2) +
  geom_node_point(aes(size = if_else(type == "KO", degree, 3), color = type)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c(KO = "#FDAE61", Pathway = "#313695")) + theme_void()
write_graph(as.igraph(sub_st), file = "ST_KO_Pathway_Subnetwork_top10.gml", format = "gml")













##############################################################
# ———— 1. 加载必要包 ————
library(tidyverse)
library(ggsci)
library(ggrepel)
library(forcats)

# ———— 2. 读取数据 ————
diff_kegg <- read.csv("KEGG_diff_results_IT_vs_ST.csv", stringsAsFactors = FALSE)
kegg_paths <- read.csv("KEGG_paths.csv", stringsAsFactors = FALSE)

# ———— 3. 统计每条 Level3 通路的背景 KO 数 ————


# 假设 kegg_paths 长这样：
#    KO    paths
# 1  K00001  Metabolism; Carbohydrate metabolism; Glycolysis / Gluconeogenesis | Metabolism; Lipid metabolism; Fatty acid degradation | …

bg_counts <- kegg_paths %>%
  # 1) 先按 “ | ” 拆分每个通路条目到多行
  separate_rows(paths, sep = " \\| ") %>%
  mutate(paths = str_trim(paths)) %>%
  # 2) 再把每行按 “; ” 拆成 Level1, Level2, Level3
  separate(
    col = paths,
    into = c("Level1", "Level2", "Level3"),
    sep = "; ",
    fill = "right",
    extra = "merge"
  ) %>%
  # 3) 只保留确实有三级信息的行
  filter(!is.na(Level3) & Level3 != "") %>%
  # 4) 去重每个 KO–Level3 组合
  distinct(KO, Level3) %>%
  # 5) 统计每个 Level3 的背景 KO 数
  count(Level3, name = "bg_count")

# 查看一下结果
head(bg_counts)

# ———— 4. 挑选显著差异 KO 并展开 Level3 ————
sig_kegg <- diff_kegg %>%
  filter(FDR < 0.05, abs(logFC) > 1) %>%
  separate_rows(Level3, sep = " \\| ") %>%
  mutate(Level3 = str_trim(Level3)) %>%
  filter(Level3 != "")

# ———— 5. 汇总各通路富集信息 ————
enriched <- sig_kegg %>%
  group_by(Level3) %>%
  summarise(
    n_sig = n(),                                      # 差异 KO 数
    avg_logFC = mean(logFC),                         # 平均 logFC
    min_FDR = min(FDR)                               # 最小 FDR
  ) %>%
  left_join(bg_counts, by = "Level3") %>%             # 加入背景基因数
  mutate(
    rich_factor = n_sig / bg_count,                   # Rich factor
    Group = if_else(avg_logFC > 0, "IT-enriched", "ST-enriched")
  ) %>%
  arrange(Group, min_FDR)


# 已有 enriched 表，包含 Level3, n_sig, avg_logFC, min_FDR, bg_count, rich_factor, Group

enriched2 <- enriched %>%
  mutate(
    Function2 = case_when(
      # 先把疾病类注释都归到 Other
      str_detect(Level3,
                 regex("cancer|infection|diabetes|myocarditis|leukemia|hepatocellular|cardiomyopathy",
                       ignore_case = TRUE)) ~ "Other",
      # 然后再打其它二级功能标签
      str_detect(Level3, regex("metabolism|biosynthesis|degradation|pathway",
                               ignore_case = TRUE))      ~ "Metabolism",
      str_detect(Level3, regex("signaling|signal|kinase|receptor",
                               ignore_case = TRUE))      ~ "Signaling",
      str_detect(Level3, regex("cell cycle|transport|trafficking|phagosome",
                               ignore_case = TRUE))      ~ "Cellular_process",
      str_detect(Level3, regex("toxin|resistance|antibiotic|macrolide",
                               ignore_case = TRUE))      ~ "Resistance",
      # 剩余的都归到 Other
      TRUE                                               ~ "Other"
    )
  ) %>%
  arrange(Group, Function2, min_FDR)


# 查看打标签后前几行
enriched2 %>% select(Group, Level3, Function2, n_sig, rich_factor) %>% head()

# ———— 6. 各组取前 20 通路 ————
top_enriched <- enriched2 %>%
  group_by(Group) %>%
  slice_min(order_by = min_FDR, n = 20) %>%
  ungroup()

# ———— 7. 绘制气泡图 ————
p <- top_enriched %>% 
  mutate(Level3 = fct_reorder(Level3, rich_factor)) %>%
  ggplot(aes(rich_factor, Level3, size = n_sig, color = -log10(min_FDR))) +
  geom_point(alpha = 0.8) +
  scale_color_gradientn(
    colours = pal_nejm()(7),
    name    = "-log10(FDR)"
  ) +
  scale_size_area(
    max_size = 8,
    name     = "n_sig"
  ) +
  labs(
    x        = "Rich factor (n_sig / background KO count)",
    y        = NULL,
    title    = "KEGG Pathway Enrichment Bubble Plot",
    subtitle = "FDR < 0.05 & |log2FC| > 1"
  ) +
  facet_wrap(~Group, ncol = 1, scales = "free_y") +
  theme_bw(base_size = 13, base_family = "Helvetica") +
  theme(
    strip.text        = element_text(face = "bold", size = 12),
    axis.text.y       = element_text(size = 8),
    legend.position   = "right"
  ) +
  geom_text_repel(
    aes(label = sprintf("%.2f", rich_factor)),
    size       = 2.5,
    show.legend = FALSE,
    seed       = 42
  )

p <- p + theme(
  panel.grid.minor    = element_blank(),
  axis.title.y        = element_blank(),
  plot.margin         = margin(5,5,5,5)
) + 
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, max(top_enriched$rich_factor) * 1.1)
  )

print(p)
# ———— 8. 保存高质量输出 ————
ggsave("KEGG_Enrichment_BubblePlot-1.pdf",
       p, width = 7, height = 8, device = cairo_pdf)


library(tidygraph)
library(ggraph)
library(igraph)   # 导出 GML 需要 igraph

# 前面步骤同您的代码
it_edges <- sig_kegg %>%
  filter(logFC > 0, FDR < 0.05) %>%
  distinct(KO, Level3)

nodes <- bind_rows(
  tibble(name = unique(it_edges$KO),     type = "KO"),
  tibble(name = unique(it_edges$Level3), type = "Pathway")
)
edges <- it_edges %>% rename(from = KO, to = Level3)
graph_it <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE) %>%
  activate(nodes) %>%
  mutate(degree = centrality_degree())

top_kos <- graph_it %>%
  as_tibble() %>%
  filter(type == "KO") %>%
  slice_max(degree, n = 10) %>%
  pull(name)

subgraph <- graph_it %>%
  activate(nodes) %>%
  filter(
    (type == "KO"      & name %in% top_kos) |
      (type == "Pathway" & name %in% it_edges$Level3[it_edges$KO %in% top_kos])
  )
# 4) 可视化
ggraph(subgraph, layout = "fr") +
  geom_edge_link(alpha = 0.2) +
  geom_node_point(aes(
    size  = if_else(type == "KO", degree, 3),
    color = type
  )) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c(KO = "#D73027", Pathway = "#4575B4")) +
  theme_void()
# ======================
# 导出 GML 文件部分
# ======================
# 1. 转 igraph 对象
sub_ig <- as.igraph(subgraph)
# 2. 写入 .gml 文件
write_graph(sub_ig, file = "IT_KO_Pathway_Subnetwork_top10.gml", format = "gml")
write_graph(as.igraph(graph_it), file = "IT_KO_Pathway_FullNetwork.gml", format = "gml")














###########################################################
library(tidyverse)
library(ggsci)
library(ggrepel)
library(showtext)
showtext_auto()

# 1. Load and preprocess data
diff_kegg <- read.csv("KEGG_diff_results_IT_vs_ST.csv", stringsAsFactors = FALSE)
sig_kegg <- diff_kegg %>%
  filter(FDR < 0.05 & abs(logFC) > 1) %>%
  separate_rows(Level3, sep = "\\|") %>%
  mutate(Level3 = str_trim(Level3))

# 2. Summarize pathway enrichment
enriched_pathways <- sig_kegg %>%
  group_by(Level3) %>%
  summarise(
    KO_count = n(),
    avg_logFC = mean(logFC),
    min_FDR = min(FDR)
  ) %>%
  mutate(
    Group = case_when(
      avg_logFC > 0 ~ "IT-enriched",
      avg_logFC < 0 ~ "ST-enriched",
      TRUE ~ "Mixed"
    )
  ) %>%
  filter(Group != "Mixed") %>% # only keep clear enrichment direction
  ungroup()

# 3. Select top 20 by FDR for each group
top_IT <- enriched_pathways %>%
  filter(Group == "IT-enriched") %>%
  arrange(min_FDR) %>%
  slice_head(n = 20)

top_ST <- enriched_pathways %>%
  filter(Group == "ST-enriched") %>%
  arrange(min_FDR) %>%
  slice_head(n = 20)

top_pathways <- bind_rows(top_IT, top_ST)

# 4. Bubble plot
p <- ggplot(top_pathways, 
            aes(x = avg_logFC, y = reorder(Level3, avg_logFC),
                size = KO_count, color = Group)) +
  geom_point(alpha = 0.85) +
  scale_color_manual(values = c("IT-enriched" = "#E94F37", 
                                "ST-enriched" = "#3E92CC")) +
  scale_size_continuous(name = "Number of significant KOs", range = c(3, 13)) +
  labs(
    x = "Average logFC (IT vs ST)",
    y = "KEGG Pathway",
    title = "Top 20 IT-enriched and ST-enriched KEGG Pathways (by FDR)",
    subtitle = "Each group: top 20 significant pathways (FDR < 0.05 & |logFC| > 1)",
    color = "Group"
  ) +
  geom_text_repel(aes(label = KO_count), size = 3, color = "black", show.legend = FALSE) +
  theme_bw(base_size = 15) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    legend.position = "right"
  )

print(p)

# 5. Save as high-res PDF
ggsave("Top20_KEGG_Pathways_PerGroup_BubblePlot-1.pdf", p, width = 14, height = 12)

















# 加载必要包
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(ggrepel)

# 1. 数据预处理 --------------------------------------------------------------
kegg_data <- read.csv("KEGG_diff_results_IT_vs_ST.csv", stringsAsFactors = FALSE)

# 计算Rich factor (差异基因数/通路总基因数)
# 注：需要补充每个KO对应的背景基因总数，这里用示例数据代替
background_counts <- data.frame(
  KO = unique(kegg_data$KO),
  bg_count = sample(50:200, length(unique(kegg_data$KO)), replace = TRUE)
)

kegg_processed <- kegg_data %>%
  left_join(background_counts, by = "KO") %>%
  separate_rows(Level3, sep = " \\| ") %>%
  mutate(Level3 = str_trim(Level3)) %>%
  filter(Level3 != "") %>%
  group_by(Level3, Significance) %>%
  summarise(
    n_genes = n(),
    rich_factor = n_genes / sum(bg_count),
    min_FDR = min(FDR)
  ) %>%
  ungroup() %>%
  filter(Significance != "Non-significant") %>%
  arrange(min_FDR)

# 2. 分开IT和ST富集通路 -----------------------------------------------------
# IT富集通路
it_enriched <- kegg_processed %>%
  filter(str_detect(Significance, "IT enriched")) %>%
  slice_min(min_FDR, n = 15)  # 取最显著前20

# ST富集通路
st_enriched <- kegg_processed %>%
  filter(str_detect(Significance, "ST enriched")) %>%
  slice_min(min_FDR, n = 15)

# 3. 绘制IT富集气泡图 ------------------------------------------------------
it_bubble <- ggplot(it_enriched, 
                    aes(x = rich_factor, 
                        y = fct_reorder(Level3, rich_factor))) +
  geom_point(aes(size = n_genes, color = -log10(min_FDR))) +
  scale_size_continuous(
    range = c(3, 8),
    name = "Gene count",
    breaks = c(min(it_enriched$n_genes), 
               median(it_enriched$n_genes), 
               max(it_enriched$n_genes))
  ) +
  scale_color_gradientn(
    colours = rev(brewer.pal(11, "Spectral")),
    name = "-log10(FDR)",
    guide = guide_colorbar(reverse = FALSE)
  ) +
  labs(
    x = "Rich factor (IT-enriched)",
    y = "",
    title = "Top 20 IT-enriched KEGG Pathways"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(color = "black", size = 8),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    legend.key.size = unit(0.4, "cm")
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")

print(it_bubble)
# 4. 绘制ST富集气泡图 ------------------------------------------------------
st_bubble <- ggplot(st_enriched, 
                    aes(x = rich_factor, 
                        y = fct_reorder(Level3, rich_factor))) +
  geom_point(aes(size = n_genes, color = -log10(min_FDR))) +
  scale_size_continuous(
    range = c(3, 8),
    name = "Gene count",
    breaks = c(min(st_enriched$n_genes), 
               median(st_enriched$n_genes), 
               max(st_enriched$n_genes))
  ) +
  scale_color_gradientn(
    colours = rev(brewer.pal(11, "Spectral")),
    name = "-log10(FDR)",
    guide = guide_colorbar(reverse = FALSE)
  ) +
  labs(
    x = "Rich factor (ST-enriched)",
    y = "",
    title = "Top 20 ST-enriched KEGG Pathways"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(color = "black", size = 8),
    axis.text.x = element_text(color = "black", size = 8),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    legend.key.size = unit(0.4, "cm")
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")

print(st_bubble)
# 5. 合并图形并保存 --------------------------------------------------------
library(patchwork)

combined_plot <- (it_bubble | st_bubble) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("KEGG_enrichment_bubbles.pdf", combined_plot, 
       width = 14, height = 8, dpi = 300)
ggsave("KEGG_enrichment_bubbles.png", combined_plot, 
       width = 14, height = 8, dpi = 300)













# 加载必要包
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(ggrepel)
library(forcats)
library(ComplexHeatmap)
library(circlize)

# 1. 数据预处理 --------------------------------------------------------------
kegg_data <- read.csv("KEGG_diff_results_IT_vs_ST.csv", stringsAsFactors = FALSE)

# 提取Level3通路信息并整理
kegg_processed <- kegg_data %>%
  separate_rows(Level3, sep = " \\| ") %>%
  mutate(Level3 = str_trim(Level3)) %>%
  filter(Level3 != "") %>%
  group_by(Level3) %>%
  summarise(
    n_genes = n(),
    avg_logFC = mean(logFC),
    min_FDR = min(FDR),
    Significance = first(Significance)
  ) %>%
  mutate(
    pathway_type = case_when(
      str_detect(Significance, "IT enriched") ~ "IT-enriched",
      str_detect(Significance, "ST enriched") ~ "ST-enriched",
      TRUE ~ "Non-significant"
    ),
    log10_FDR = -log10(min_FDR)
  ) %>%
  filter(pathway_type != "Non-significant") %>%
  arrange(desc(avg_logFC))

# 2. 气泡图（Bubble plot）----------------------------------------------------
top_pathways <- 20 # 展示top20通路

bubble_plot <- kegg_processed %>%
  slice_max(order_by = abs(avg_logFC), n = top_pathways) %>%
  ggplot(aes(x = avg_logFC, y = fct_reorder(Level3, avg_logFC))) +
  geom_point(aes(size = n_genes, fill = avg_logFC), shape = 21, color = "black") +
  scale_fill_gradient2(
    low = "#2166AC", 
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Log2 Fold Change"
  ) +
  scale_size_continuous(
    range = c(3, 8),
    name = "Gene count"
  ) +
  geom_text(
    aes(label = ifelse(log10_FDR > 10, 
                       formatC(min_FDR, format = "e", digits = 1), 
                       ""),
        hjust = 0.5, vjust = -1.5, size = 2.5
    )) +
      labs(
        x = "Log2 Fold Change (IT vs ST)",
        y = "KEGG Pathway",
        title = "Differentially Enriched KEGG Pathways"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color = "black", size = 8),
        axis.text.x = element_text(color = "black", size = 8),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        legend.key.size = unit(0.4, "cm")
      ) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 35))
    
    # 保存气泡图
    ggsave("KEGG_bubble_plot.pdf", bubble_plot, 
           width = 8, height = 6, dpi = 300)
    
    # 3. 高级热图（ComplexHeatmap）-----------------------------------------------
    # 准备矩阵数据
    heatmap_data <- kegg_data %>%
      filter(FDR < 0.01 & abs(logFC) > 1.5) %>%
      group_by(Level3) %>%
      summarise(
        mean_logFC = mean(logFC),
        gene_count = n()
      ) %>%
      separate_rows(Level3, sep = " \\| ") %>%
      mutate(Level3 = str_trim(Level3)) %>%
      filter(Level3 != "") %>%
      group_by(Level3) %>%
      summarise(
        mean_logFC = mean(mean_logFC),
        gene_count = sum(gene_count)
      ) %>%
      arrange(desc(abs(mean_logFC))) %>%
      slice_head(n = 25)  # 选择top25通路
    
    # 创建注释条
    pathway_annot <- ifelse(heatmap_data$mean_logFC > 0, "IT-enriched", "ST-enriched")
    names(pathway_annot) <- heatmap_data$Level3
    
    # 颜色映射
    col_fun <- colorRamp2(
      c(min(heatmap_data$mean_logFC), 0, max(heatmap_data$mean_logFC)), 
      c("#4575B4", "white", "#D73027")
    )
    
    # 绘制热图
    kegg_heatmap <- Heatmap(
      matrix = as.matrix(heatmap_data$mean_logFC),
      name = "Log2FC",
      col = col_fun,
      row_labels = str_wrap(heatmap_data$Level3, width = 40),
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 8),
      row_title = "KEGG Pathways",
      column_title = "Log2 Fold Change (IT vs ST)",
      column_title_side = "bottom",
      show_column_names = FALSE,
      row_split = pathway_annot,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      rect_gp = gpar(col = "white", lwd = 0.5),
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = 8),
        labels_gp = gpar(fontsize = 7)
      ),
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(
          round(heatmap_data$gene_count[i], 0), 
          x, y, 
          gp = gpar(fontsize = 6))
      }
    )
    
    # 保存热图
    pdf("KEGG_complex_heatmap.pdf", width = 6, height = 8)
    draw(kegg_heatmap)
    dev.off()
    
    # 4. 通路分类网络图 ---------------------------------------------------------
    library(ggraph)
    library(tidygraph)
    
    # 创建通路-基因关系网络
    network_data <- kegg_data %>%
      filter(FDR < 0.01 & abs(logFC) > 1.5) %>%
      separate_rows(Level3, sep = " \\| ") %>%
      mutate(Level3 = str_trim(Level3)) %>%
      filter(Level3 != "") %>%
      select(KO, Level3, logFC) %>%
      group_by(Level3, KO) %>%
      summarise(logFC = mean(logFC)) %>%
      ungroup()
    
    # 创建网络对象
    kegg_graph <- as_tbl_graph(network_data, directed = FALSE) %>%
      mutate(
        node_type = ifelse(name %in% unique(network_data$Level3), "pathway", "gene"),
        pathway_class = ifelse(node_type == "pathway", 
                               ifelse(name %in% unique(network_data$Level3[network_data$logFC > 0]),
                                      "IT-enriched", "ST-enriched"), NA)
      )
    
    # 绘制网络图
    set.seed(123)
    network_plot <- ggraph(kegg_graph, layout = "fr") + 
      geom_edge_link(alpha = 0.1, width = 0.2) +
      geom_node_point(aes(
        color = pathway_class,
        size = ifelse(node_type == "pathway", 4, 1),
        alpha = ifelse(node_type == "pathway", 1, 0.5)
      )) +
      geom_node_text(
        aes(label = ifelse(node_type == "pathway", name, ""),
            size = ifelse(node_type == "pathway", 3, 0)),
        repel = TRUE,
        max.overlaps = 100,
        segment.color = NA
      ) +
      scale_color_manual(
        values = c("IT-enriched" = "#D73027", "ST-enriched" = "#4575B4"),
        na.value = "grey70"
      ) +
      scale_size_identity() +
      scale_alpha_identity() +
      theme_void() +
      labs(
        title = "KEGG Pathway-Gene Interaction Network",
        color = "Pathway Type"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
      )
    
    # 保存网络图
    ggsave("KEGG_network_plot.pdf", network_plot, 
           width = 10, height = 8, dpi = 300)
    
    # 5. 通路分类柱状图 --------------------------------------------------------
    # 按功能大类分类统计
    kegg_category <- kegg_data %>%
      filter(FDR < 0.01 & abs(logFC) > 1.5) %>%
      mutate(
        Category = case_when(
          str_detect(Level3, "signal|adhesion|receptor") ~ "Cell Signaling",
          str_detect(Level3, "metaboli") ~ "Metabolism",
          str_detect(Level3, "synthesis|biosynthesis") ~ "Biosynthesis",
          str_detect(Level3, "phosphorylation|oxidative") ~ "Energy Production",
          str_detect(Level3, "infection|virus") ~ "Host-Pathogen Interaction",
          TRUE ~ "Other"
        )
      ) %>%
      group_by(Category, Significance) %>%
      summarise(count = n()) %>%
      ungroup() %>%
      complete(Category, Significance, fill = list(count = 0)) %>%
      filter(!is.na(Category))
    
    # 绘制分类柱状图
    category_plot <- ggplot(kegg_category, 
                            aes(x = reorder(Category, count), 
                                y = count, 
                                fill = Significance)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      geom_text(aes(label = count), 
                position = position_dodge(width = 0.8),
                hjust = -0.2, size = 3) +
      scale_fill_manual(
        values = c("IT enriched (FDR<1%, LFC>1.5)" = "#D73027", 
                   "ST enriched (FDR<1%, LFC>1.5)" = "#4575B4"),
        labels = c("IT-enriched", "ST-enriched")
      ) +
      labs(
        x = "Functional Category",
        y = "Number of Significant KEGG Terms",
        title = "Functional Distribution of Enriched KEGG Pathways"
      ) +
      coord_flip() +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
    # 保存分类图
    ggsave("KEGG_category_plot.pdf", category_plot, 
           width = 7, height = 5, dpi = 300)
    