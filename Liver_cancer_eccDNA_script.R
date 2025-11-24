# 加载必要的包
library(openxlsx)
library(ggpubr)
library(ggplot2)
library(stringr)

# 初始化函数
source("0.fun.R")
init()

# 将数据保存为RDS格式（只需运行一次）----
# 癌症样本
cancer_dir <- paste0("/home/liuyuchen/Projects/eccDNA/",
                     "data/cancer/13.eccDNAanno_v2.1")
cancer_files <- list.files(cancer_dir, pattern = "anno\\.xls", full.names = TRUE)

# 正常样本
normal_dir <- "/home/liuyuchen/Projects/eccDNA/data/normal/7.eccDNAanno_v2.1"
normal_files <- list.files(normal_dir, pattern = "anno\\.xls", full.names = TRUE)

# 创建目录保存RDS文件
dir.create("rds_data", showWarnings = FALSE)

# 处理并保存癌症样本
for (file in cancer_files) {
  sample_id <- str_match(file, ".*/(.*)\\.anno\\.xls")[, 2]
  rt <- read.table(file, header = TRUE, sep = "\t", comment.char = "")
  saveRDS(rt, file = paste0("rds_data/", sample_id, ".rds"))
}

# 处理并保存正常样本
for (file in normal_files) {
  sample_id <- str_match(file, ".*/(.*)\\.anno\\.xls")[, 2]
  rt <- read.table(file, header = TRUE, sep = "\t", comment.char = "")
  saveRDS(rt, file = paste0("rds_data/", sample_id, ".rds"))
}


# 改文件名使得更易读
# 设置工作目录（根据实际情况修改）
setwd("/home/liuyuchen/Projects/eccDNA/rds_data")

# 获取所有RDS文件
rds_files <- list.files(pattern = "\\.rds$")

# 计数器
renamed_count <- 0

# 遍历并重命名文件
for (file in rds_files) {
  if (grepl("-N-1-3\\.rds$", file)) {
    # 构建新文件名
    new_name <- gsub("-N-1-3\\.rds$", ".N.rds", file)
    
    # 执行重命名
    file.rename(file, new_name)
    
    # 输出日志
    message(paste("重命名:", file, "->", new_name))
    renamed_count <- renamed_count + 1
  }
}

# 输出总结
message(paste("\n完成! 共重命名", renamed_count, "个文件。"))
message("原始文件列表:")
message(paste(rds_files, collapse = "\n"))

# 显示更新后的文件列表
new_files <- list.files(pattern = "\\.rds$")
message("\n更新后的文件列表:")
message(paste(new_files, collapse = "\n"))

# save.image("/home/liuyuchen/Projects/eccDNA/workspace.RData")





# 从RDS文件读取数据进行分析（后续分析使用这部分）----
setwd("/home/liuyuchen/Projects/eccDNA")


rds_files <- list.files("rds_data", pattern = "\\.rds$", full.names = TRUE)
sample_ids <- gsub("\\.rds$", "", basename(rds_files))

# 初始化统计表格
nCategory <- c("lncRNA", "mRNA", "miRNA", "pseudogene", "intergenic")
Category_stat <- data.frame(matrix(0, ncol = length(sample_ids), nrow = length(nCategory) + 2))
colnames(Category_stat) <- sample_ids
rownames(Category_stat) <- c(nCategory, "multipleType", "Total")

Genic_stat <- data.frame(matrix(0, ncol = length(sample_ids), nrow = 2))
colnames(Genic_stat) <- sample_ids
rownames(Genic_stat) <- c("Genic", "Intergenic")

Subgenic_stat <- data.frame(matrix(0, ncol = length(sample_ids), nrow = 4))
colnames(Subgenic_stat) <- sample_ids
rownames(Subgenic_stat) <- c("lncRNA", "mRNA", "miRNA", "pseudogene")

# 处理每个样本
for (i in seq_along(rds_files)) {
  id <- sample_ids[i]
  message(paste0("Analyzing sample: ", id))
  
  # 从RDS读取数据
  rt <- readRDS(rds_files[i])
  
  # 初始化计数器
  multiple_type_count <- 0
  genic_count <- 0
  intergenic_count <- 0
  subgenic_counts <- setNames(rep(0, 4), c("lncRNA", "mRNA", "miRNA", "pseudogene"))
  
  # 逐行统计
  for (j in 1:nrow(rt)) {
    category <- rt$Category[j]
    
    # 按逗号分割所有类型
    categories <- unlist(strsplit(category, ","))
    categories <- trimws(categories)  # 去除空格
    
    # 检查是否为多类型
    if (length(categories) > 1) {
      multiple_type_count <- multiple_type_count + 1
    }
    
    # 统计每个类型
    for (cat in categories) {
      # 将"others"改为"intergenic"
      if (cat == "others") {
        cat <- "intergenic"
      }
      
      # 统计到对应的分类
      if (cat %in% nCategory) {
        Category_stat[cat, id] <- Category_stat[cat, id] + 1
        
        # 统计Genic/Intergenic
        if (cat == "intergenic") {
          intergenic_count <- intergenic_count + 1
        } else {
          genic_count <- genic_count + 1
          subgenic_counts[cat] <- subgenic_counts[cat] + 1
        }
      }
    }
  }
  
  # 填充统计结果
  Category_stat["multipleType", id] <- multiple_type_count
  Category_stat["Total", id] <- nrow(rt)
  
  Genic_stat["Genic", id] <- genic_count
  Genic_stat["Intergenic", id] <- intergenic_count
  
  for (subcat in rownames(Subgenic_stat)) {
    Subgenic_stat[subcat, id] <- subgenic_counts[subcat]
  }
  
  # 验证统计结果
  total_categories <- sum(Category_stat[nCategory, id])
  if ((total_categories + multiple_type_count - (length(categories) - 1) * multiple_type_count) != nrow(rt)) {
    warning(paste("Sample", id, ": category count validation failed"))
  }
}
# 输出统计结果
print("详细分类统计:")
print(Category_stat)
print("Genic/Intergenic统计:")
print(Genic_stat)
print("Genic子分类统计:")
print(Subgenic_stat)

# 保存结果
write.csv(Category_stat, "Detailed_Category_Statistics.csv")
write.csv(Genic_stat, "Genic_Intergenic_Statistics.csv")
write.csv(Subgenic_stat, "Genic_Subcategory_Statistics.csv")

# 上面的MultipleType有点问题，画图的时候不展示，用以下代码的SingleType 和 MultipleType来画图
Category_type <- data.frame(matrix(0, ncol = length(sample_ids), nrow = 2))
colnames(Category_type) <- sample_ids
rownames(Category_type) <- c("singleType", "multipleType")

# 只统计singleType和multipleType的循环
for (i in seq_along(rds_files)) {
  id <- sample_ids[i]
  message(paste0("Analyzing sample for type statistics: ", id))
  
  # 从RDS读取数据
  rt <- readRDS(rds_files[i])
  
  # 初始化计数器
  single_type_count <- 0
  multiple_type_count <- 0
  
  # 逐行统计
  for (j in 1:nrow(rt)) {
    category <- rt$Category[j]
    
    # 按逗号分割所有类型
    categories <- unlist(strsplit(category, ","))
    categories <- trimws(categories)  # 去除空格
    
    # 统计singleType vs multipleType
    if (length(categories) == 1) {
      single_type_count <- single_type_count + 1
    } else {
      multiple_type_count <- multiple_type_count + 1
    }
  }
  
  # 填充统计结果
  Category_type["singleType", id] <- single_type_count
  Category_type["multipleType", id] <- multiple_type_count
  
  # 验证统计结果
  if ((single_type_count + multiple_type_count) != nrow(rt)) {
    warning(paste("Sample", id, ": singleType + multipleType != total rows"))
  }
}

# 输出Type统计结果
print("Type统计结果 (singleType vs multipleType):")
print(Category_type)

# 保存结果
write.csv(Category_type, "Category_Type_Statistics.csv")

# save.image("/home/liuyuchen/Projects/eccDNA/workspace.RData")

# 可视化部分
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)


# 处理样本名称和分组信息----
# 读取分组信息
group_save <- read.xlsx("group_save.xlsx", sheet = 1)
rownames(group_save) <- group_save$Person

# 创建符合Cell Reports风格的条形图----
# 确定肿瘤和正常样本的分界位置
# 转换T列为二分类：Tumor (T/M) 和 Normal (N)
group_save$Group <- ifelse(group_save$Tumor %in% c("T", "M"), "Tumor", "Normal")

# 准备数据：计算每个样本的Total eccDNA数量

Category_type_df <- as.data.frame(t(Category_type)) %>%
  rownames_to_column("Sample") %>%
  mutate(Total = singleType + multipleType) %>%
  # 使用group_save中的准确分组信息和Number
  left_join(group_save %>% select(Person, Group, Number), by = c("Sample" = "Person")) %>%
  arrange(Group, Number) %>%  # 按Group和Number排序
  mutate(x_position = 1:n())  # 创建x轴位置

# 计算肿瘤样本的结束位置
tumor_end <- sum(Category_type_df$Group == "Tumor")

# 绘制图表
p <- ggplot(Category_type_df, aes(x = Number, y = Total, fill = Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_vline(xintercept = tumor_end + 0.5, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_text(aes(label = format(Total, big.mark = ",")),  # 添加千分位分隔符
            vjust = -0.8,  # 向上移动标签
            hjust = 0.5,   # 水平居中
            size = 3, 
            angle = 90, 
            position = position_nudge(x = 0.1)) +  # 向右微调
  scale_fill_manual(values = c("Tumor" = "#1F77B4", "Normal" = "#AEC7E8")) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE),  # y轴标签添加千分位分隔符
    breaks = 10^(0:5), 
    expand = expansion(mult = c(0, 0.2))) +  # 增加顶部空间
  labs(y = "eccDNA Number (log10 scaled)", x = "Sample", fill = "Group") +  # 修改纵坐标标签
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "top",
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(10, 10, 20, 10)  # 增加顶部边距
  ) +
  annotate("text", x = tumor_end/2, y = max(Category_type_df$Total) * 0.9, 
           label = "Tumor", size = 4, fontface = "bold") +
  annotate("text", x = tumor_end + (nrow(Category_type_df) - tumor_end)/2, y = max(Category_type_df$Total) * 0.9, 
           label = "Normal", size = 4, fontface = "bold") + 
  coord_cartesian(clip = "off")  # 允许标签超出绘图区域

# 保存图表
ggsave("figures/1.eccDNA_distribution.pdf", p, width = 11, height = 6, dpi = 1200)

# 显示图表
print(p)

# 额外输出数据摘要
cat("eccDNA数量统计摘要:\n")
cat("总样本数:", nrow(Category_type_df), "\n")
cat("肿瘤样本数:", tumor_end, "\n")
cat("正常样本数:", nrow(Category_type_df) - tumor_end, "\n")
cat("肿瘤样本eccDNA总数:", sum(Category_type_df$Total[Category_type_df$Group == "Tumor"]), "\n")
cat("正常样本eccDNA总数:", sum(Category_type_df$Total[Category_type_df$Group == "Normal"]), "\n")
cat("每样本平均eccDNA数:\n")
cat("  - 肿瘤样本:", round(mean(Category_type_df$Total[Category_type_df$Group == "Tumor"]), 2), "\n")
cat("  - 正常样本:", round(mean(Category_type_df$Total[Category_type_df$Group == "Normal"]), 2), "\n")





# 堆叠图：Genic vs Intergenic（使用group_save中的准确分组信息）

genic_grouped <- as.data.frame(t(Genic_stat)) %>%
  rownames_to_column("Sample") %>%
  # 使用group_save中的准确分组信息
  left_join(group_save %>% select(Person, Group), by = c("Sample" = "Person")) %>%
  # 按Group汇总，计算总和
  group_by(Group) %>%
  summarise(
    Genic = sum(Genic),
    Intergenic = sum(Intergenic)
  ) %>%
  # 计算每个组内的比例（分别计算）
  rowwise() %>%
  mutate(
    Total = Genic + Intergenic,
    Genic_Prop = Genic / Total,
    Intergenic_Prop = Intergenic / Total
  ) %>%
  # 转换为长格式
  pivot_longer(cols = c("Genic_Prop", "Intergenic_Prop"), 
               names_to = "Category", 
               values_to = "Proportion") %>%
  mutate(Category = ifelse(Category == "Genic_Prop", "Genic", "Intergenic"))

# 创建堆叠图（分别显示各组内部比例）
p1 <- ggplot(genic_grouped, aes(x = Group, y = Proportion, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.5) +
  labs(title = "Genic vs Intergenic Distribution", 
       x = "Sample Group", 
       y = "Proportion") +
  scale_fill_manual(values = c("Genic" = "#1f77b4", "Intergenic" = "#ff7f0e")) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.2),
    panel.grid.minor = element_blank()
  )

# 添加百分比标签
p1_with_labels <- p1 + 
  geom_text(aes(label = scales::percent(Proportion, accuracy = 0.1)), 
            position = position_stack(vjust = 0.5), 
            size = 6,  # 增大标签字体
            color = "white",
            fontface = "bold")

# 保存PDF
ggsave("figures/2.Genic_Intergenic_Distribution_Grouped.pdf", p1_with_labels, 
       width = 8, height = 8, device = "pdf")


# 画图展示每个样本的Intergenic 和Genic 占比

# 计算每个样本的Genic/Intergenic比例
genic_sample_proportions <- as.data.frame(t(Genic_stat)) %>%
  rownames_to_column("Sample") %>%
  # 使用group_save中的准确分组信息
  left_join(group_save %>% select(Person, Group), by = c("Sample" = "Person")) %>%
  # 计算每个样本的比例
  mutate(
    Total = Genic + Intergenic,
    Genic_Prop = Genic / Total,
    Intergenic_Prop = Intergenic / Total
  ) %>%
  select(Sample, Group, Genic_Prop, Intergenic_Prop) %>%
  # 转换为长格式，便于绘图
  pivot_longer(cols = c("Genic_Prop", "Intergenic_Prop"), 
               names_to = "Category", 
               values_to = "Proportion") %>%
  mutate(Category = ifelse(Category == "Genic_Prop", "Genic", "Intergenic"))

# 2. 进行T检验
# Genic比例的T检验
genic_ttest <- t.test(Proportion ~ Group, 
                      data = genic_sample_proportions %>% filter(Category == "Genic"))
genic_pvalue <- genic_ttest$p.value

# Intergenic比例的T检验
intergenic_ttest <- t.test(Proportion ~ Group, 
                          data = genic_sample_proportions %>% filter(Category == "Intergenic"))
intergenic_pvalue <- intergenic_ttest$p.value

# 3. 创建小提琴图
library(ggpubr)

p_violin <- ggplot(genic_sample_proportions, aes(x = Category, y = Proportion, fill = Group)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9), 
              size = 1.5, alpha = 0.6) +
  labs(title = "Genic vs Intergenic Proportion Distribution by Group",
       x = "Variant Category",
       y = "Proportion") +
  scale_fill_manual(values = c("Tumor" = "#e74c3c", "Normal" = "#3498db")) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    axis.line = element_line(color = "black", linewidth = 0.5),  # 添加坐标轴线
    axis.ticks = element_line(color = "black", linewidth = 0.5), # 添加坐标轴刻度
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_blank()       # 去掉面板边框
  )

# 4. 添加T检验结果到图中
p_violin_with_stats <- p_violin +
  stat_compare_means(
    aes(label = paste0("p = ", ..p.format..)),
    method = "t.test",
    label.x = 1.5,
    label.y = 1.05,
    size = 5
  ) 

# 5. 保存PDF
ggsave("figures/2.S.Genic_Intergenic_ViolinPlot_TTest.pdf", p_violin_with_stats, 
       width = 8, height = 10, device = "pdf")

# 6. 显示图表
print(p_violin_with_stats)

# 7. 输出详细的统计结果
cat("Genic/Intergenic比例分布统计:\n\n")

# Genic统计
cat("Genic比例:\n")
genic_summary <- genic_sample_proportions %>% 
  filter(Category == "Genic") %>%
  group_by(Group) %>%
  summarise(
    Mean = mean(Proportion),
    SD = sd(Proportion),
    N = n(),
    Min = min(Proportion),
    Max = max(Proportion)
  )
print(genic_summary)
cat("T检验p值:", format(genic_pvalue, scientific = TRUE, digits = 3), "\n\n")

# Intergenic统计
cat("Intergenic比例:\n")
intergenic_summary <- genic_sample_proportions %>% 
  filter(Category == "Intergenic") %>%
  group_by(Group) %>%
  summarise(
    Mean = mean(Proportion),
    SD = sd(Proportion),
    N = n(),
    Min = min(Proportion),
    Max = max(Proportion)
  )
print(intergenic_summary)
cat("T检验p值:", format(intergenic_pvalue, scientific = TRUE, digits = 3), "\n")




# 2. 饼图部分也使用准确的分组信息
# Genic子分类饼图
subgenic_group <- as.data.frame(t(Subgenic_stat)) %>%
  rownames_to_column("Sample") %>%
  # 使用group_save中的准确分组信息
  left_join(group_save %>% select(Person, Group), by = c("Sample" = "Person")) %>%
  group_by(Group) %>%
  summarise(across(c("lncRNA", "mRNA", "miRNA", "pseudogene"), sum)) %>%
  pivot_longer(cols = -Group, names_to = "Subcategory", values_to = "Count") %>%
  # 计算每个组内的百分比
  group_by(Group) %>%
  mutate(
    Total = sum(Count),
    Percentage = Count / Total * 100,
    Label = sprintf("%.1f%%", Percentage)
  )

# 为每个组创建饼图（增大字体并添加百分比标签）
create_pie_chart <- function(data, group_name) {
  group_data <- data %>% filter(Group == group_name)
  
  # 计算标签位置
  group_data <- group_data %>%
    arrange(desc(Subcategory)) %>%
    mutate(
      ypos = cumsum(Percentage) - 0.5 * Percentage
    )
  
  ggplot(group_data, aes(x = "", y = Percentage, fill = Subcategory)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
    coord_polar("y", start = 0) +
    geom_text(aes(y = ypos, label = Label), 
              color = "black", size = 5, fontface = "bold") +
    labs(title = paste0(group_name, ": Genic Subcategories"), fill = "Subcategory") +
     theme_void(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, 
                              margin = margin(b = 10)),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right"
  ) +
  scale_fill_manual(values = genic_colors) +
  # 添加一些美化设置
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )
}

p2_tumor <- create_pie_chart(subgenic_group, "Tumor")
p2_normal <- create_pie_chart(subgenic_group, "Normal")

# 保存饼图
ggsave("figures/2.Genic_Subcategories_Tumor_PieChart.pdf", p2_tumor, width = 8, height = 8, device = "pdf")
ggsave("figures/2.Genic_Subcategories_Normal_PieChart.pdf", p2_normal, width = 8, height = 8, device = "pdf")

# 显示饼图
print(p2_tumor)
print(p2_normal)

# 输出详细的百分比统计
cat("Genic子分类百分比统计:\n\n")
cat("Tumor组:\n")
tumor_stats <- subgenic_group %>% filter(Group == "Tumor")
for(i in 1:nrow(tumor_stats)) {
  cat(sprintf("  %s: %.1f%% (%d variants)\n", 
              tumor_stats$Subcategory[i], 
              tumor_stats$Percentage[i],
              tumor_stats$Count[i]))
}

cat("\nNormal组:\n")
normal_stats <- subgenic_group %>% filter(Group == "Normal")
for(i in 1:nrow(normal_stats)) {
  cat(sprintf("  %s: %.1f%% (%d variants)\n", 
              normal_stats$Subcategory[i], 
              normal_stats$Percentage[i],
              normal_stats$Count[i]))
}

group_save$Group <- ifelse(group_save$Tumor %in% c("T", "M"), "Tumor", "Normal")

# 准备数据：计算每个样本的Total eccDNA数量

Category_stat_df <- as.data.frame(t(Category_type)) %>%
  rownames_to_column("Sample") %>%
  mutate(Total = singleType + multipleType) %>%
  # 使用group_save中的准确分组信息和Number
  left_join(group_save %>% select(Person, Group, Number), by = c("Sample" = "Person")) %>%
  arrange(Group, Number) %>%  # 按Group和Number排序
  mutate(x_position = 1:n())  # 创建x轴位置



#  eccDNA数量箱型图 Tumor vs Normal 分组
use_colors <- c("Tumor" = "#1F77B4", "Normal" = "#AEC7E8") # 颜色方案
#  设置因子顺序为Tumor在前Normal在后
group_save$Group <- ifelse(group_save$Tumor %in% c("T", "M"), "Tumor", "Normal")
Category_stat_df$Group <- factor(Category_stat_df$Group, levels = c("Tumor", "Normal"))
ggplot(Category_stat_df, aes(x = Group, y = Total, fill = Group)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white") +
  scale_fill_manual(values = use_colors, guide = "none") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5) +
  scale_y_continuous(trans = "log10") +
  xlab("") + ylab("EccDNA Number per Sample") +
  theme_pubr() +
  theme(text = element_text(size = 18))

ggsave("figures/S1.TumorNormal.group.numbercompare.pdf", width = 5, height = 6)

#  eccDNA数量箱型图 Tumor vs Normal 分组
# 用新的代码画上述分组图


# 转换Group为二分类：Stage 3-4 和 Stage 1-2
group_save$Group <- ifelse(group_save$Stage==1|group_save$Stage==2,"Stage 1-2","Stage 3-4")

# 转换统计结果为长格式并合并分组信息
Category_stat_df <- as.data.frame(t(Category_stat))
rownames(Category_stat_df)=Category_stat_df$Person
Category_stat_df <- merge(Category_stat_df, group_save[, c("Person", "Number", "Group")], by = "Person") #Number是每个病人的ID
Category_stat_df=Category_stat_df[!is.na(Category_stat_df$Group),]
# 设置因子水平保持原始顺序
Category_stat_df$Group <- factor(Category_stat_df$Group, levels = c("Stage 3-4", "Stage 1-2"))
use_colors <- c("Stage 3-4" = "#1F77B4", "Stage 1-2" = "#AEC7E8") # 颜色方案
ggplot(Category_stat_df, aes(x = Group, y = Total, fill = Group)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white") +
  scale_fill_manual(values = use_colors, guide = "none") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5) +
  scale_y_continuous(trans = "log10") +
  xlab("") + ylab("EccDNA Number per Sample") +
  theme_pubr() +
  theme(text = element_text(size = 18))

ggsave("figures/S1.Stage.group.numbercompare.pdf", width = 5, height = 6)


# 转换Group为二分类：Metastasis tumor 和 Primary tumor
group_save$Group <- ifelse(group_save$Metastasis=="T","Primary tumor","Metastasis tumor")

# 转换统计结果为长格式并合并分组信息
Category_stat_df <- as.data.frame(t(Category_stat))
rownames(Category_stat_df)=Category_stat_df$Person
Category_stat_df <- merge(Category_stat_df, group_save[, c("Person", "Number", "Group")], by = "Person") #Number是每个病人的ID
Category_stat_df=Category_stat_df[!is.na(Category_stat_df$Group),]
# 设置因子水平保持原始顺序
Category_stat_df$Group <- factor(Category_stat_df$Group, levels = c("Metastasis tumor", "Primary tumor"))
use_colors <- c("Metastasis tumor" = "#1F77B4", "Primary tumor" = "#AEC7E8") # 颜色方案
ggplot(Category_stat_df, aes(x = Group, y = Total, fill = Group)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white") +
  scale_fill_manual(values = use_colors, guide = "none") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5) +
  scale_y_continuous(trans = "log10") +
  xlab("") + ylab("EccDNA Number per Sample") +
  theme_pubr() +
  theme(text = element_text(size = 18))

ggsave("figures/S1.Metastasis.group.numbercompare.pdf", width = 5, height = 6)


# 箱型图展示 singleType 和 MultipleType 两种eccDNA的差异，还有eccDNA包含基因的categories
# 修正prepare_data函数来处理Category_type的数据结构
prepare_data_category_type <- function(data, group_save) {
  # 将数据转置，使singleType和multipleType成为列名
  df_transposed <- as.data.frame(t(data))
  
  # 添加Person列
  df_transposed$Person <- rownames(df_transposed)
  
  # 转换为长格式
  df_long <- df_transposed %>%
    pivot_longer(cols = -Person, names_to = "Category", values_to = "Count")
  
  # 添加分组信息
  df_long <- df_long %>%
    left_join(group_save %>% select(Person, Group), by = "Person") %>%
    filter(!is.na(Group))
  
  # 计算每个样本的总数（用于计算比例）
  total_counts <- df_long %>%
    group_by(Person) %>%
    summarise(Total = sum(Count), .groups = "drop")
  
  # 合并总数信息并计算比例
  df_long <- df_long %>%
    left_join(total_counts, by = "Person") %>%
    mutate(Proportion = ifelse(Total > 0, Count / Total * 100, 0))
  
  return(df_long)
}

# 定义颜色方案
box_colors <- c(Tumor = "#E71B1C, Normal = "#2762AE")

# 1. SingleType vs MultiType 分类数据
df_single_multi <- prepare_data_category_type(Category_type, group_save)

# 检查数据
head(df_single_multi)

# Count箱型图 - SingleType/MultiType
p_values_count <- df_single_multi %>%
  group_by(Category) %>%
  summarise(
    p_value = wilcox.test(Count ~ Group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_signif = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    x_pos = match(Category, unique(df_single_multi$Category)),
    y_pos = max(df_single_multi$Count, na.rm = TRUE) * 1.1
  )
df_single_multi$Group=factor(df_single_multi$Group,levels = c("Tumor","Normal"))
ggplot(df_single_multi, aes(x = Category, y = Count, fill = Group)) +
  geom_boxplot(
    outlier.shape = 21, 
    outlier.fill = "white", 
    position = position_dodge(0.8)
  ) +
  geom_text(
    data = p_values_count,
    aes(x = x_pos, y = y_pos, label = p_signif),
    inherit.aes = FALSE,
    size = 6,
    fontface = "bold"
  ) +
  scale_y_log10(labels = scales::comma) +
  scale_fill_manual(values = box_colors) +
  labs(x = "Functional Category", y = "eccDNA Number (log10 scaled)") +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank()
  )
ggsave("figures/S1.TumorNormal_category_count_single_multi.pdf", width = 5, height = 6)

# Proportion箱型图 - SingleType/MultiType
p_values_prop <- df_single_multi %>%
  group_by(Category) %>%
  summarise(
    p_value = wilcox.test(Proportion ~ Group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_signif = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    x_pos = match(Category, unique(df_single_multi$Category))
  )

# 计算每个Category的最大比例值作为y轴参考
max_proportions <- df_single_multi %>%
  group_by(Category) %>%
  summarise(max_prop = max(Proportion, na.rm = TRUE), .groups = "drop")
p_values_prop$y_pos <- max_proportions$max_prop + 5

ggplot(df_single_multi, aes(x = Category, y = Proportion, fill = Group)) +
  geom_boxplot(
    outlier.shape = 21, 
    outlier.fill = "white", 
    position = position_dodge(0.8)
  ) +
  geom_text(
    data = p_values_prop,
    aes(x = x_pos, y = y_pos, label = p_signif),
    inherit.aes = FALSE,
    size = 6,
    fontface = "bold"
  ) +
  scale_y_sqrt(
    limits = c(0, 100), 
    labels = scales::percent_format(scale = 1)
  ) +
  scale_fill_manual(values = box_colors) +
  labs(x = "Functional Category", y = "Proportion (%)") +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank()
  )
ggsave("figures/S1.TumorNormal_category_proportion_single_multi.pdf", width = 5, height = 6)

# 2. 对于Subgenic_stat数据，使用原来的函数
prepare_data_subgenic <- function(data, categories, group_save) {
  # 将数据转换为长格式
  df_long <- as.data.frame(t(data)) %>%
    rownames_to_column("Person") %>%
    pivot_longer(cols = -Person, names_to = "Category", values_to = "Count") %>%
    filter(Category %in% categories)
  
  # 添加分组信息
  df_long <- df_long %>%
    left_join(group_save %>% select(Person, Group), by = "Person") %>%
    filter(!is.na(Group))
  
  # 计算每个样本的总数（用于计算比例）
  total_counts <- df_long %>%
    group_by(Person) %>%
    summarise(Total = sum(Count), .groups = "drop")
  
  # 合并总数信息并计算比例
  df_long <- df_long %>%
    left_join(total_counts, by = "Person") %>%
    mutate(Proportion = ifelse(Total > 0, Count / Total * 100, 0))
  
  return(df_long)
}

# Genic子分类数据
categories_genic <- c("lncRNA", "mRNA", "miRNA", "pseudogene")
df_genic <- prepare_data_subgenic(Subgenic_stat, categories_genic, group_save)

# 检查数据
head(df_genic)
df_genic$Group=factor(df_genic$Group,levels = c("Tumor","Normal"))
box_colors <- c(Tumor = "#1F77B4", Normal = "#AEC7E8")

# 然后继续绘制Genic子分类的箱型图（使用之前修正后的代码）
# Count箱型图 - Genic子分类
p_values_count_genic <- df_genic %>%
  group_by(Category) %>%
  summarise(
    p_value = wilcox.test(Count ~ Group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_signif = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    x_pos = match(Category, categories_genic),
    y_pos = max(df_genic$Count, na.rm = TRUE) * 1.1
  )
df_genic$Category=factor(df_genic$Category,levels=c("miRNA","mRNA","lncRNA","pseudogene"))
ggplot(df_genic, aes(x = Category, y = Count, fill = Group)) +
  geom_boxplot(
    outlier.shape = 21, 
    outlier.fill = "white", 
    position = position_dodge(0.8)
  ) +
  geom_text(
    data = p_values_count_genic,
    aes(x = x_pos, y = y_pos, label = p_signif),
    inherit.aes = FALSE,
    size = 6,
    fontface = "bold"
  ) +
  scale_y_log10(labels = scales::comma) +
  scale_fill_manual(values = box_colors) +
  labs(x = "Functional Category", y = "eccDNA Number per Sample (log10 scaled)") +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank()
  )
ggsave("figures/S1.TumorNormal_category_count_genic.pdf", width = 9, height = 6)

# Proportion箱型图 - Genic子分类
p_values_prop_genic <- df_genic %>%
  group_by(Category) %>%
  summarise(
    p_value = wilcox.test(Proportion ~ Group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_signif = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    x_pos = match(Category, categories_genic)
  )

max_proportions_genic <- df_genic %>%
  group_by(Category) %>%
  summarise(max_prop = max(Proportion, na.rm = TRUE), .groups = "drop")
p_values_prop_genic$y_pos <- max_proportions_genic$max_prop + 5

ggplot(df_genic, aes(x = Category, y = Proportion, fill = Group)) +
  geom_boxplot(
    outlier.shape = 21, 
    outlier.fill = "white", 
    position = position_dodge(0.8)
  ) +
  geom_text(
    data = p_values_prop_genic,
    aes(x = x_pos, y = y_pos, label = p_signif),
    inherit.aes = FALSE,
    size = 6,
    fontface = "bold"
  ) +
  scale_y_sqrt(
    limits = c(0, 100), 
    labels = scales::percent_format(scale = 1)
  ) +
  scale_fill_manual(values = box_colors) +
  labs(x = "Functional Category", y = "Proportion (%)") +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank()
  )
ggsave("figures/S1.TumorNormal_category_proportion_genic.pdf", width = 9, height = 6)







# eccDNA长度分布分析 ----

# 提取每个样本的eccDNA长度并存储
rds_files <- list.files("rds_data", pattern = "\\.rds$", full.names = TRUE)
sample_ids <- gsub("\\.rds$", "", basename(rds_files))
ecc_length_list <- list()  # 存储每个样本的长度列表

for (i in seq_along(rds_files)) {
  id <- sample_ids[i]
  rt <- readRDS(rds_files[i])
  rt$length <- as.numeric(rt$End) - as.numeric(rt$Start)  # 计算长度
  ecc_length_list[[id]] <- rt$length
}


# Tumor vs Normal 长度分布对比（密度图+累积分布图）
# 整合所有样本长度数据
tumor_samples <- group_save[group_save$Tumor == "T", "Person"]  # 从group_save获取Tumor样本ID
normal_samples <- group_save[group_save$Tumor == "N", "Person"]  # 从group_save获取Normal样本ID

# 提取对应样本的eccDNA长度（确保样本ID匹配）
# 过滤ecc_length_list中属于Tumor组的样本
tumor_len <- unlist(ecc_length_list[names(ecc_length_list) %in% tumor_samples])
# 过滤ecc_length_list中属于Normal组的样本
normal_len <- unlist(ecc_length_list[names(ecc_length_list) %in% normal_samples])
len_df <- data.frame(
  Length = c(tumor_len, normal_len),
  Group = factor(c(rep("Tumor", length(tumor_len)), rep("Normal", length(normal_len))),levels = c("Tumor", "Normal"))
)

# 执行KS检验
ks_test <- ks.test(tumor_len, normal_len)
ks_p_value <- ks_test$p.value
ks_p_label <- ifelse(ks_p_value < 0.001, "P < 0.001", paste("P =", format(ks_p_value, digits = 3)))


use_colors <- c(
  "Tumor" = "#0A4C6A",  # 深青色（深色）
  "Normal" = "#AEC7E8"  # 蓝色（浅色，比Tumor浅一度）
)
# 密度图
# 计算Tumor组密度曲线的峰值位置
tumor_density <- density(tumor_len, from = 0, to = 4050)
# 找到所有局部峰值（density值比左右两侧都大的点）
is_peak <- diff(sign(diff(tumor_density$y))) < 0
peak_indices <- which(is_peak) + 1  # 修正索引

# 选择密度最高的4个峰值
top_peak_indices <- peak_indices[order(tumor_density$y[peak_indices], decreasing = TRUE)[1:4]]
top_peak_positions <- tumor_density$x[top_peak_indices]
top_peak_values <- tumor_density$y[top_peak_indices]
top_peak_positions=round(top_peak_positions,0)
# 密度图（添加峰值标注）
p_density <- ggdensity(len_df, x = "Length", 
                       fill = "Group",
                       palette = use_colors,
                       alpha = 0.6, 
                       xlab = "Fragment length (bp)", 
                       ylab = "Density") +
  scale_y_continuous(limits = c(0, 0.005)) +
  scale_x_continuous(limits = c(0, 4050), 
                     breaks = c(0, 2000, 4000)) +

# 原有p值标注（y值略低于检验方法文字）
  annotate("text", x = 2500, y = 0.004,  
           label = paste0("K-S test, ",ks_p_label), size = 3, fontface = "bold") +
  
# 添加Tumor组峰值的垂直线段和标注
  geom_segment(data = data.frame(x = top_peak_positions, 
                                 y = rep(0, length(top_peak_positions)),
                                 xend = top_peak_positions,
                                 yend = top_peak_values),
               aes(x = x, y = y, xend = xend, yend = yend),
               linetype = "dashed", color = "#0A4C6A", size = 0.5) +
  
  geom_text(data = data.frame(x = top_peak_positions, 
                              y = top_peak_values + 0.0003,  # 稍微高于峰值
                              label = paste0(round(top_peak_positions), " bp")),
            aes(x = x, y = y, label = label),
            color = "#0A4C6A", size = 4, fontface = "bold") +
  
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave("figures/1.TumorNormal.length_density.pdf", p_density, width = 6, height = 4)
ggsave("figures/1.TumorNormal.length_density_rect.pdf", p_density, width = 5, height = 5)
# 累积分布图


p_ecdf <- ggplot(len_df, aes(x = Length, color = Group)) +
  stat_ecdf(geom = "step", size = 1.2) +
  scale_color_manual(values = use_colors) +
  xlim(0, 10000) +
  ylim(0, 1) +
  xlab("Fragment length (bp)") + 
  ylab("Cumulative Frequency") +
  
  # 添加K-S检验p值标注
  annotate("text", x = 5000, y = 0.8,
           label = paste0("K-S test, ", ks_p_label), size = 3, fontface = "bold") +
  
  # 计算并标注累积频率为0.5时的长度（中位数）
  geom_vline(aes(xintercept = median(tumor_len)), linetype = "dashed", color = use_colors["Tumor"], size = 0.8) +
  geom_vline(aes(xintercept = median(normal_len)), linetype = "dashed", color = use_colors["Normal"], size = 0.8) +
  
  # 添加中位数文本标注（稍微偏移以避免重叠）
  annotate("text", x = median(tumor_len), y = 0.45,
           label = paste0("Median: ", round(median(tumor_len)), " bp"),
           color = use_colors["Tumor"], hjust = -0.1, size = 4) +
  annotate("text", x = median(normal_len), y = 0.55,
           label = paste0("Median: ", round(median(normal_len)), " bp"),
           color = use_colors["Normal"], hjust = 1.1, size = 4) +
  
  # 添加水平参考线到中位数
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "gray50", size = 0.5) +
  
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


ggsave("figures/1.TumorNormal.length_ecdf.pdf", p_ecdf,  width = 6, height = 4)
ggsave("figures/1.TumorNormal.length_ecdf.pdf", p_ecdf,  width = 5, height = 5)


# Metastasis tumor vs Primary tumor 长度分布对比（密度图+累积分布图）

metastasis_samples <- group_save[group_save$Metastasis == "M", "Person"]  # 从group_save获取metastasis样本ID
primary_samples <- group_save[group_save$Metastasis == "T", "Person"]  # 从group_save获取Primary样本ID

# 提取对应样本的eccDNA长度（确保样本ID匹配）
# 过滤ecc_length_list中属于Metastasis组的样本
metastasis_len <- unlist(ecc_length_list[names(ecc_length_list) %in% metastasis_samples])
# 过滤ecc_length_list中属于Primary组的样本
primary_len <- unlist(ecc_length_list[names(ecc_length_list) %in% primary_samples])
len_df <- data.frame(
  Length = c(metastasis_len, primary_len),
  Group = factor(c(rep("Metastasis", length(metastasis_len)), rep("Primary", length(primary_len))),levels = c("Metastasis", "Primary"))
)

# 执行KS检验
ks_test <- ks.test(metastasis_len, primary_len)
ks_p_value <- ks_test$p.value
ks_p_label <- ifelse(ks_p_value < 0.001, "P < 0.001", paste("P =", format(ks_p_value, digits = 3)))


use_colors <- c(
  "Metastasis" = "#0A4C6A",  # 深青色（深色）
  "Primary" = "#AEC7E8"  # 蓝色（浅色，比Metastasis浅一度）
)
# 密度图
# 计算Metastasis组密度曲线的峰值位置
metastasis_density <- density(metastasis_len, from = 0, to = 4050)
# 找到所有局部峰值（density值比左右两侧都大的点）
is_peak <- diff(sign(diff(metastasis_density$y))) < 0
peak_indices <- which(is_peak) + 1  # 修正索引

# 选择密度最高的4个峰值
top_peak_indices <- peak_indices[order(metastasis_density$y[peak_indices], decreasing = TRUE)[1:4]]
top_peak_positions <- metastasis_density$x[top_peak_indices]
top_peak_values <- metastasis_density$y[top_peak_indices]
top_peak_positions=round(top_peak_positions,0)
# 密度图
p_density <- ggdensity(len_df, x = "Length", 
                       fill = "Group",
                       palette = use_colors,
                       alpha = 0.6, 
                       xlab = "Fragment length (bp)", 
                       ylab = "Density") +
  scale_y_continuous(limits = c(0, 0.006)) +
  scale_x_continuous(limits = c(0, 4050), 
                     breaks = c(0, 2000, 4000)) +
  
  # 原有p值标注（y值略低于检验方法文字）
  annotate("text", x = 3500, y = 0.004,  
           label = paste0("K-S test, ",ks_p_label), size = 5, fontface = "bold") +
  
  # 添加Metastasis组峰值的垂直线段和标注
  geom_segment(data = data.frame(x = top_peak_positions, 
                                 y = rep(0, length(top_peak_positions)),
                                 xend = top_peak_positions,
                                 yend = top_peak_values),
               aes(x = x, y = y, xend = xend, yend = yend),
               linetype = "dashed", color = "#0A4C6A", size = 0.5) +
  
  geom_text(data = data.frame(x = top_peak_positions, 
                              y = top_peak_values + 0.0003,  # 稍微高于峰值
                              label = paste0(round(top_peak_positions), " bp")),
            aes(x = x, y = y, label = label),
            color = "#0A4C6A", size = 4, fontface = "bold") +
  
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave("figures/1.MetastasisPrimary.length_density.pdf", p_density, width = 6, height = 4)
# 累积分布图


p_ecdf <- ggplot(len_df, aes(x = Length, color = Group)) +
  stat_ecdf(geom = "step", size = 1.2) +
  scale_color_manual(values = use_colors) +
  xlim(0, 10000) +
  ylim(0, 1) +
  xlab("Fragment length (bp)") + 
  ylab("Cumulative Frequency") +
  
  # 添加K-S检验p值标注
  annotate("text", x = 7000, y = 0.8,
           label = paste0("K-S test, ", ks_p_label), size = 5, fontface = "bold") +
  
  # 计算并标注累积频率为0.5时的长度（中位数）
  geom_vline(aes(xintercept = median(metastasis_len)), linetype = "dashed", color = use_colors["Metastasis"], size = 0.8) +
  geom_vline(aes(xintercept = median(primary_len)), linetype = "dashed", color = use_colors["Primary"], size = 0.8) +
  
  # 添加中位数文本标注（稍微偏移以避免重叠）
  annotate("text", x = median(metastasis_len), y = 0.45,
           label = paste0("Median: ", round(median(metastasis_len)), " bp"),
           color = use_colors["Metastasis"], hjust = -0.1, size = 4) +
  annotate("text", x = median(primary_len), y = 0.55,
           label = paste0("Median: ", round(median(primary_len)), " bp"),
           color = use_colors["Primary"], hjust = 1.1, size = 4) +
  
  # 添加水平参考线到中位数
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "gray50", size = 0.5) +
  
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


ggsave("figures/1.MetastasisPrimary.length_ecdf.pdf", p_ecdf,  width = 6, height = 4)




# Stage3-4 vs Stage1-2 长度分布对比（密度图+累积分布图）
group_save$Group <- ifelse(group_save$Stage==1|group_save$Stage==2,"Stage 1-2","Stage 3-4")
stage34_samples <- group_save[group_save$Group  == "Stage 3-4", "Person"]  # 从group_save获取metastasis样本ID
stage12_samples <- group_save[group_save$Group  == "Stage 1-2", "Person"]  # 从group_save获取Primary样本ID

# 提取对应样本的eccDNA长度（确保样本ID匹配）
# 过滤ecc_length_list中属于Metastasis组的样本
stage34_len <- unlist(ecc_length_list[names(ecc_length_list) %in% stage34_samples])
# 过滤ecc_length_list中属于Primary组的样本
stage12_len <- unlist(ecc_length_list[names(ecc_length_list) %in% stage12_samples])
len_df <- data.frame(
  Length = c(stage34_len, stage12_len),
  Group = factor(c(rep("Stage 3-4", length(stage34_len)), rep("Stage 1-2", length(stage12_len))),levels = c("Stage 3-4", "Stage 1-2"))
)

# 执行KS检验
ks_test <- ks.test(stage34_len, stage12_len)
ks_p_value <- ks_test$p.value
ks_p_label <- ifelse(ks_p_value < 0.001, "P < 0.001", paste("P =", format(ks_p_value, digits = 3)))


use_colors <- c(
  `Stage 3-4` = "#0A4C6A",  # 深青色（深色）
  `Stage 1-2` = "#AEC7E8"  # 蓝色（浅色，比Metastasis浅一度）
)
# 密度图
# 计算Metastasis组密度曲线的峰值位置
metastasis_density <- density(stage34_len, from = 0, to = 4050)
# 找到所有局部峰值（density值比左右两侧都大的点）
is_peak <- diff(sign(diff(metastasis_density$y))) < 0
peak_indices <- which(is_peak) + 1  # 修正索引

# 选择密度最高的4个峰值
top_peak_indices <- peak_indices[order(metastasis_density$y[peak_indices], decreasing = TRUE)[1:4]]
top_peak_positions <- metastasis_density$x[top_peak_indices]
top_peak_values <- metastasis_density$y[top_peak_indices]
top_peak_positions=round(top_peak_positions,0)
# 密度图（添加峰值标注）
p_density <- ggdensity(len_df, x = "Length", 
                       fill = "Group",
                       palette = use_colors,
                       alpha = 0.6, 
                       xlab = "Fragment length (bp)", 
                       ylab = "Density") +
  scale_y_continuous(limits = c(0, 0.006)) +
  scale_x_continuous(limits = c(0, 4050), 
                     breaks = c(0, 2000, 4000)) +
  
  # 原有p值标注（y值略低于检验方法文字）
  annotate("text", x = 3500, y = 0.004,  
           label = paste0("K-S test, ",ks_p_label), size = 5, fontface = "bold") +
  
  # 添加Metastasis组峰值的垂直线段和标注
  geom_segment(data = data.frame(x = top_peak_positions, 
                                 y = rep(0, length(top_peak_positions)),
                                 xend = top_peak_positions,
                                 yend = top_peak_values),
               aes(x = x, y = y, xend = xend, yend = yend),
               linetype = "dashed", color = "#0A4C6A", size = 0.5) +
  
  geom_text(data = data.frame(x = top_peak_positions, 
                              y = top_peak_values + 0.0003,  # 稍微高于峰值
                              label = paste0(round(top_peak_positions), " bp")),
            aes(x = x, y = y, label = label),
            color = "#0A4C6A", size = 4, fontface = "bold") +
  
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ggsave("figures/1.Stage.length_density.pdf", p_density, width = 6, height = 4)
# 累积分布图


p_ecdf <- ggplot(len_df, aes(x = Length, color = Group)) +
  stat_ecdf(geom = "step", size = 1.2) +
  scale_color_manual(values = use_colors) +
  xlim(0, 10000) +
  ylim(0, 1) +
  xlab("Fragment length (bp)") + 
  ylab("Cumulative Frequency") +
  
  # 添加K-S检验p值标注
  annotate("text", x = 7000, y = 0.8,
           label = paste0("K-S test, ", ks_p_label), size = 5, fontface = "bold") +
  
  # 计算并标注累积频率为0.5时的长度（中位数）
  geom_vline(aes(xintercept = median(stage34_len)), linetype = "dashed", color = use_colors["Stage 3-4"], size = 0.8) +
  geom_vline(aes(xintercept = median(stage12_len)), linetype = "dashed", color = use_colors["Stage 1-2"], size = 0.8) +
  # 添加中位数文本标注（稍微偏移以避免重叠）
  annotate("text", x = median(stage34_len), y = 0.45,
           label = paste0("Median: ", round(median(stage34_len)), " bp"),
           color = use_colors["Stage 3-4"], hjust = -0.1, size = 4) +
  annotate("text", x = median(stage12_len), y = 0.55,
           label = paste0("Median: ", round(median(stage12_len)), " bp"),
           color = use_colors["Stage 1-2"], hjust = 1.1, size = 4) +
  
  # 添加水平参考线到中位数
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "gray50", size = 0.5) +
  
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


ggsave("figures/1.Stage.length_ecdf.pdf", p_ecdf,  width = 6, height = 4)




# eccDNA长度分段统计（<2k, 2k~8k, 8k~10k, >10k）
length_cut <- function(len) {
  cut(len, 
      breaks = c(-Inf, 2000, 8000, 10000, Inf),
      labels = c("<2k", "2k~8k", "8k~10k", ">10k"))
}

# 统计每个样本的长度分段比例

len_segment_stat <- data.frame()
for (id in sample_ids) {
  len <- as.numeric(ecc_length_list[[id]])
  segment <- length_cut(len)
  
  # 计算各段比例（确保结果正确）
  segment_counts <- table(segment)
  segment_props <- prop.table(segment_counts) * 100  # 转换为百分比
  
  # 创建样本的比例数据框
  temp <- data.frame(
    Person = id,
    LengthSegment = names(segment_props),
    Proportion = as.numeric(segment_props)
  )
  
  len_segment_stat <- rbind(len_segment_stat, temp)
}



# 合并分组信息
group_save$Group <- ifelse(group_save$Tumor %in% c("T", "M"), "Tumor", "Normal")
len_segment_stat <- merge(len_segment_stat, group_save[, c("Person", "Group", "Number")], by = "Person")

# 验证数据格式（确保Proportion是0-100的数值）
str(len_segment_stat)
head(len_segment_stat)

# 设置分组顺序：Metastasis tumor, Primary tumor, Normal
len_segment_stat$Group <- factor(len_segment_stat$Group, 
                                 levels = c("Tumor", "Normal"))

# 设置长度段顺序和颜色
len_segment_stat$LengthSegment <- factor(len_segment_stat$LengthSegment, 
                                         levels = c("<2k", "2k~8k", "8k~10k", ">10k"))

# 指定颜色：绿、深绿、蓝、深蓝
custom_colors <- c(
  "<2k"      = "#B2DF8A",  # 绿色
  "2k~8k"    = "#33A02C",  # 深绿色
  "8k~10k"   = "#A6CEE3",  # 蓝色
  ">10k"     = "#1F78B4"   # 深蓝色
)

# 绘制堆叠柱状图
ggplot(len_segment_stat, aes(x = factor(Number), y = Proportion, fill = LengthSegment)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~Group, scales = "free_x",space="free_x") +  # 按组分开显示
  scale_fill_manual(values = custom_colors) +  # 使用自定义颜色
  xlab("") + 
  ylab("Proportion (%)") +
  theme_pubr() +
  theme(
    axis.text.x = element_blank(),  # 移除x轴文字
    axis.ticks.x = element_blank(),  # 移除x轴刻度线
    legend.position = "top",
    legend.title = element_blank()
  )
ggsave("figures/1.length_segment_stackbar_tumornormal.pdf", width = 8, height = 4)

# --------------------------------------------------------------------------------------------------------


# 小提琴图
library(dplyr)
library(ggplot2)
library(ggpubr)

# 定义长度分段函数（4个类别）
length_cut_extended <- function(len) {
  cut(
    len, 
    breaks = c(-Inf, 2000, 8000, 10000, Inf),
    labels = c("<2k", "2k~8k", "8k~10k", ">10k"),
    right = FALSE
  )
}

# 合并所有样本数据（原始每条eccDNA + 分组信息）
violin_raw <- data.frame()
for (id in sample_ids) {
  len_vec <- as.numeric(ecc_length_list[[id]])
  seg_vec <- length_cut_extended(len_vec)
  
  group_info <- subset(group_save, Person == id, select = c("Person", "Group", "Number"))
  
  temp_df <- data.frame(
    Person = id,
    Length = len_vec,
    LengthCategory = seg_vec,
    Group = group_info$Group,
    Number = group_info$Number
  )
  
  violin_raw <- rbind(violin_raw, temp_df)
}

# 为了加快后续处理，对每段长度按比例抽样（例如保留万分之一）
sample_frac <- 1e-2  # 可调整
violin_sampled <- violin_raw %>%
  group_by(LengthCategory) %>%
  slice_sample(prop = sample_frac, replace = FALSE) %>%
  ungroup()

# 标记肿瘤与正常
violin_sampled$TumorStatus <- ifelse(
  violin_sampled$Group %in% c("Metastasis tumor", "Primary tumor"),
  "Tumor",
  "Normal"
)
violin_sampled$TumorStatus <- factor(violin_sampled$TumorStatus, levels = c("Tumor", "Normal"))
violin_sampled$LengthCategory <- factor(
  violin_sampled$LengthCategory,
  levels =  c("<2k", "2k~8k", "8k~10k", ">10k")
)

# 计算每个样本在每个长度段的比例
violin_data <- violin_sampled %>%
  group_by(Person, LengthCategory, TumorStatus) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Person) %>%
  mutate(Proportion = Count / sum(Count) * 100) %>%
  ungroup()

# Wilcoxon 检验
p_values <- violin_data %>%
  group_by(LengthCategory) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(Proportion[TumorStatus == "Tumor"],
                  Proportion[TumorStatus == "Normal"],
                  alternative = "two.sided")$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_label = case_when(
      is.na(p_value) ~ "NA",
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

# 定义颜色
use_colors <- c("Tumor" = "#0A4C6A", "Normal" = "#AEC7E8")

# 绘图
 ggplot(
  violin_data,
  aes(x = LengthCategory, y = Proportion, fill = TumorStatus)
) +
  geom_violin(position = position_dodge(0.9), alpha = 0.7, width = 0.8) +
  geom_text(
    data = p_values,
    aes(x = LengthCategory, y = max(violin_data$Proportion, na.rm = TRUE) + 5, label = p_label),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = use_colors) + scale_y_sqrt(breaks = c(0, 10, 20, 30, 90),
                                                        limits = c(0, max(violin_data$Proportion, na.rm = TRUE) * 1.2))+
  ylab("Proportion (%)") +
  xlab("eccDNA Length Category") +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", size = 0.2)
  )

# 保存图
ggsave("figures/S1.length_category_violin_proportion_plot.pdf", width = 8, height = 5)


# 统计每个样本在不同长度分段的 eccDNA 数量
violin_count_data <- violin_raw %>%
  group_by(Person, LengthCategory, TumorStatus) %>%
  summarise(Count = n(), .groups = "drop")

# 统计每类是否存在显著性差异（Wilcoxon检验）
p_values_count <- violin_count_data %>%
  group_by(LengthCategory) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(Count[TumorStatus == "Tumor"],
                  Count[TumorStatus == "Normal"],
                  alternative = "two.sided")$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_label = case_when(
      is.na(p_value) ~ "NA",
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

# 绘图：Count 小提琴图
ggplot(
  violin_count_data,
  aes(x = LengthCategory, y = Count, fill = TumorStatus)
) +
  geom_boxplot(position = position_dodge(0.9), alpha = 0.7, width = 0.8, outlier.size = 1) +
  geom_text(
    data = p_values_count,
    aes(x = LengthCategory, y = max(violin_count_data$Count, na.rm = TRUE) + 5, label = p_label),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +scale_y_log10()+
  scale_fill_manual(values = use_colors) +
  ylab("eccDNA Number (log10 scaled)") +
  xlab("eccDNA Length Category") +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", size = 0.2)
  )

# 保存图
ggsave("figures/S1.length_category_count_plot.pdf", width = 8, height = 5)


#-------------------

# eccDNA per Mb 作图

#-------------------



#-----------------------------------------
# 基于RDS文件的eccDNA per Mb差异分析和箱型图绘制
#-----------------------------------------
# 加载必要的包
library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)
library(ggpubr)
library(purrr)

# 设置工作目录
setwd("/home/liuyuchen/Projects/eccDNA")

# 确保输出目录存在
dir.create("figures/ChrCompare", recursive = TRUE, showWarnings = FALSE)

# 读取group_save.xlsx获取样本分组信息
group_data <- readxl::read_excel("group_save.xlsx")

# 读取hg18染色体长度信息
hg18_karyotype <- read.table("hg18_karyotype_data.txt", header = FALSE)
colnames(hg18_karyotype) <- c("chr", "id", "start", "end", "label")

# 创建染色体长度的映射
chr_length <- hg18_karyotype$end
names(chr_length) <- hg18_karyotype$chr

# 处理RDS文件并计算每个样本每条染色体的eccDNA per Mb
sample_data_list <- list()
rds_dir <- "rds_data"

# 获取所有RDS文件
rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)

# 处理每个RDS文件
for (rds_file in rds_files) {
  sample_name <- sub("\\.rds$", "", basename(rds_file))
  cat("Processing sample:", sample_name, "\n")
  
  sample_data <- readRDS(rds_file)
  
  # 检查X.Chr列是否存在
  if (!"X.Chr" %in% colnames(sample_data)) {
    cat("Warning: 'X.Chr' column not found in", sample_name, "data. Skipping.\n")
    next
  }
  
  # 统计每条染色体的eccDNA数量
  chr_counts <- table(sample_data$X.Chr)
  
  # 计算每条染色体的eccDNA per Mb
  chr_per_mb <- list()
  for (chromosome in names(chr_counts)) {
    if (chromosome %in% names(chr_length)) {
      mb_length <- chr_length[chromosome] / 1000000
      per_mb <- chr_counts[chromosome] / mb_length
      chr_per_mb[[chromosome]] <- per_mb
    }
  }
  
  # 转换为数据框
  per_mb_df <- data.frame(
    sample = sample_name,
    chr = names(chr_per_mb),
    eccDNA_per_Mb = unlist(chr_per_mb)
  )
  
  sample_data_list[[sample_name]] <- per_mb_df
}

# 合并所有样本数据
all_data <- do.call(rbind, sample_data_list)
group_data=group_data[,c("Person","Tumor","Stage","Metastasis")]
# 去掉性染色体和Mitochondria
all_data <- all_data[!all_data$chr %in% c("chrX", "chrY", "chrM"), ]
# 添加样本分组信息
all_data <- merge(all_data, group_data, by.x = "sample", by.y = "Person", all.x = TRUE)

#提取Tumor列的分组信息
all_data$Group <- ifelse(all_data$Tumor=="T", "Tumor", "Normal")

# 合并样本总数和分组信息


# 计算Group列中各分类的平均eccDNA密度
cat("\n计算Group列中各分类的平均eccDNA密度...\n")
group_avg_density <- all_data %>% group_by(Group) %>% summarize(Average_eccDNA_per_Mb = mean(eccDNA_per_Mb, na.rm = TRUE),standard_eccDNA_per_Mb=sd(eccDNA_per_Mb, na.rm = TRUE))

# 输出结果
print("Group分类的平均eccDNA密度（per Mb）:")
print(group_avg_density)
#  Group  Average_eccDNA_per_Mb
#   <chr>                  <dbl>                  <dbl>
# 1 Normal                  13.2                   5.23
# 2 Tumor                   65.0                 133. 
# 按染色体统计eccDNA密度
cat("\n按染色体统计eccDNA密度...\n")
# 提取染色体数字，创建适当的因子顺序

# 按染色体和分组统计平均eccDNA密度
chr_group_avg_density <- all_data %>% 
  group_by(chr, Group) %>% 
  summarize(Average_eccDNA_per_Mb = sum(eccDNA_per_Mb, na.rm = TRUE)) %>%
  arrange(chr)

# 输出结果
print("按染色体和分组统计的eccDNA数量（per Mb）:")
print(chr_group_avg_density)


# 计算每个染色体的总体平均eccDNA密度（不分组）
chr_group_avg_density1 <- chr_group_avg_density %>% 
  group_by(Group) %>% 
  summarize(Average_eccDNA_per_Mb1 = mean(Average_eccDNA_per_Mb, na.rm = TRUE))

print("按染色体统计的总体平均eccDNA密度（per Mb）:")
print(chr_group_avg_density1)

#   Group  Average_eccDNA_per_Mb_per_chr
#   <chr>                   <dbl>
# 1 Normal                   145.
# 2 Tumor                   1754.
# 这个数据后续可以用于标示在circos图中间，有别于前面只统计总数。


# 画箱型图
#提取Tumor列的分组信息
all_data$Group <- ifelse(all_data$Tumor=="T", "Tumor", "Normal")
all_data$Group=factor(all_data$Group,levels=c("Tumor","Normal"))

# 添加染色体排序
all_data$chr <- factor(all_data$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
                                              "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", 
                                              "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"))

# 准备用于箱型图的数据
# 按染色体和分组计算平均密度
chr_group_stats <- all_data %>% 
  group_by(chr, Group) %>% 
  summarize(mean_eccDNA = mean(eccDNA_per_Mb, na.rm = TRUE),
            sd_eccDNA = sd(eccDNA_per_Mb, na.rm = TRUE),
            n = n())

# 创建箱型图 - 更接近示例图的样式
p <- ggplot(all_data, aes(x = chr, y = log(eccDNA_per_Mb,2), fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) + # 移除离群点显示
  ylim(0, 10) + # 设置纵坐标范围为0-300
  scale_fill_manual(values = c("Normal" = "#2762AE", "Tumor" = "#E71B1C")) + # 使用与示例图一致的颜色
  labs(title = "",
       x = "",
       y = "Log2 Normalized eccDNA per Mb") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15, angle = 90, vjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "white"))

# 添加统计显著性标记（使用Wilcoxon秩和检验）
# 修复：使用更直接的方式计算每个染色体的统计检验
wilcox_results <- data.frame()
for (chr in unique(all_data$chr)) {
  chr_data <- all_data[all_data$chr == chr, ]
  if (length(unique(chr_data$Group)) > 1) {
    test_result <- t.test(log(eccDNA_per_Mb,2) ~ Group, data = chr_data)
    wilcox_results <- rbind(wilcox_results, data.frame(chr = chr, p_value = test_result$p.value))
  }
}

# 添加显著性标记（模仿示例图的样式，将所有标记放在顶部）
# 找到整个数据集的最大y值
max_overall_y <- 8

# 创建所有染色体的显著性标记文本
all_sig_labels <- sapply(unique(all_data$chr), function(chr) {
  if (chr %in% wilcox_results$chr) {
    p_val <- wilcox_results$p_value[wilcox_results$chr == chr]
    # 根据p值确定显著性标记
    if (p_val < 0.0001) {
      return("****")
    } else if (p_val < 0.001) {
      return("***")
    } else if (p_val < 0.01) {
      return("**")
    } else if (p_val < 0.05) {
      return("*")
    }
  }
  return("ns")
})

# 在顶部添加所有显著性标记
p <- p + annotate("text", 
                  x = unique(all_data$chr), 
                  y = rep(max_overall_y, length(unique(all_data$chr))), 
                  label = all_sig_labels, 
                  size = 3, 
                  hjust = 0.5, 
                  vjust = 0)

# 保存箱型图
pdf("figures/ChrCompare/Tumor_vs_Normal_chr_perMb_Boxplot.pdf", width = 10, height = 5)
print(p)
dev.off()




#提取Stage列的分组信息
all_data$Group <- ifelse(all_data$Stage==1|all_data$Stage==2,"Stage 1-2","Stage 3-4")
all_data$Group=factor(all_data$Group,levels=c("Stage 3-4","Stage 1-2"))
all_data_stage=all_data[which(!is.na(all_data$Stage)),]
# 添加染色体排序
all_data_stage$chr <- factor(all_data_stage$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
                                              "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", 
                                              "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"))

# 准备用于箱型图的数据
# 按染色体和分组计算平均密度
chr_group_stats <- all_data_stage %>% 
  group_by(chr, Group) %>% 
  summarize(mean_eccDNA = mean(eccDNA_per_Mb, na.rm = TRUE),
            sd_eccDNA = sd(eccDNA_per_Mb, na.rm = TRUE),
            n = n())

# 创建箱型图 - 更接近示例图的样式
p <- ggplot(all_data_stage, aes(x = chr, y = log(eccDNA_per_Mb,2), fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) + # 移除离群点显示
  ylim(0, 10) +scale_fill_manual(values = c("Stage 1-2" = "#2762AE", "Stage 3-4" = "#E71B1C")) + # 使用与示例图一致的颜色
  labs(title = "",
       x = "",
       y = "Log2 Normalized eccDNA per Mb") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15, angle = 90, vjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "white"))

# 添加统计显著性标记（使用Wilcoxon秩和检验）
# 修复：使用更直接的方式计算每个染色体的统计检验
wilcox_results <- data.frame()
for (chr in unique(all_data_stage$chr)) {
  chr_data <- all_data_stage[all_data_stage$chr == chr, ]
  if (length(unique(chr_data$Group)) > 1) {
    test_result <- t.test(log(eccDNA_per_Mb,2) ~ Group, data = chr_data)
    wilcox_results <- rbind(wilcox_results, data.frame(chr = chr, p_value = test_result$p.value))
  }
}

# 添加显著性标记（模仿示例图的样式，将所有标记放在顶部）
# 找到整个数据集的最大y值
max_overall_y <- 9.5

# 创建所有染色体的显著性标记文本
all_sig_labels <- sapply(unique(all_data_metastasis$chr), function(chr) {
  if (chr %in% wilcox_results$chr) {
    p_val <- wilcox_results$p_value[wilcox_results$chr == chr]
    # 根据p值确定显著性标记
    if (p_val < 0.0001) {
      return("****")
    } else if (p_val < 0.001) {
      return("***")
    } else if (p_val < 0.01) {
      return("**")
    } else if (p_val < 0.05) {
      return("*")
    }
  }
  return("ns") # 不显示非显著结果
})

# 在顶部添加所有显著性标记
p <- p + annotate("text", 
                  x = unique(all_data_metastasis$chr), 
                  y = rep(max_overall_y, length(unique(all_data_metastasis$chr))), 
                  label = all_sig_labels, 
                  size = 3, 
                  hjust = 0.5, 
                  vjust = 0)

# 保存箱型图
pdf("figures/ChrCompare/Stage34_vs_Stage12_chr_perMb_Boxplot.pdf", width = 10, height = 5)
print(p)
dev.off()


#提取Metastasis列的分组信息
all_data$Group <- ifelse(all_data$Metastasis=="T","Primary tumor","Metastasis tumor")
all_data$Group=factor(all_data$Group,levels=c("Metastasis tumor","Primary tumor"))
all_data_metastasis=all_data[which(!is.na(all_data$Metastasis)),]
# 添加染色体排序
all_data_metastasis$chr <- factor(all_data_metastasis$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
                                              "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", 
                                              "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"))

# 准备用于箱型图的数据
# 按染色体和分组计算平均密度
chr_group_stats <- all_data_metastasis %>% 
  group_by(chr, Group) %>% 
  summarize(mean_eccDNA = mean(eccDNA_per_Mb, na.rm = TRUE),
            sd_eccDNA = sd(eccDNA_per_Mb, na.rm = TRUE),
            n = n())

# 创建箱型图 - 更接近示例图的样式
p <- ggplot(all_data_metastasis, aes(x = chr, y = log(eccDNA_per_Mb,2), fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) + # 移除离群点显示
  ylim(0, 12) + 
  scale_fill_manual(values = c(`Metastasis tumor` ="#E71B1C" ,`Primary tumor` = "#2762AE" )) + 
  labs(title = "",
       x = "",
       y = "Log2 Normalized eccDNA per Mb") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15, angle = 90, vjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "white"))

# 添加统计显著性标记（使用Wilcoxon秩和检验）
# 修复：使用更直接的方式计算每个染色体的统计检验
wilcox_results <- data.frame()
for (chr in unique(all_data_metastasis$chr)) {
  chr_data <- all_data_metastasis[all_data_metastasis$chr == chr, ]
  if (length(unique(chr_data$Group)) > 1) {
    test_result <- t.test(log(eccDNA_per_Mb,2) ~ Group, data = chr_data)
    wilcox_results <- rbind(wilcox_results, data.frame(chr = chr, p_value = test_result$p.value))
  }
}

# 添加显著性标记（模仿示例图的样式，将所有标记放在顶部）
# 找到整个数据集的最大y值
max_overall_y <- 10

# 创建所有染色体的显著性标记文本
all_sig_labels <- sapply(unique(all_data$chr), function(chr) {
  if (chr %in% wilcox_results$chr) {
    p_val <- wilcox_results$p_value[wilcox_results$chr == chr]
    # 根据p值确定显著性标记
    if (p_val < 0.0001) {
      return("****")
    } else if (p_val < 0.001) {
      return("***")
    } else if (p_val < 0.01) {
      return("**")
    } else if (p_val < 0.05) {
      return("*")
    }
  }
  return("ns")
})

# 在顶部添加所有显著性标记
p <- p + annotate("text", 
                  x = unique(all_data$chr), 
                  y = rep(max_overall_y, length(unique(all_data$chr))), 
                  label = all_sig_labels, 
                  size = 3, 
                  hjust = 0.5, 
                  vjust = 0)

# 保存箱型图
pdf("figures/ChrCompare/Metastasis_vs_Primary_chr_perMb_Boxplot.pdf", width = 10, height = 5)
print(p)
dev.off()
# save.image("/home/liuyuchen/Projects/eccDNA/workspace.RData")


# 加载必要的包
library(readxl)

#-----------------------------------------
# 生成各分组的bed文件用于circos图
#-----------------------------------------
cat("\n生成各分组的bed文件用于circos图...\n")

# 函数：根据分组信息生成bed文件
generate_group_bed_files <- function() {
  # 确保输出目录存在
  dir.create("figures/ChrCompare", recursive = TRUE, showWarnings = FALSE)
  
  # 处理Tumor分组
  cat("处理Tumor分组...\n")
  tumor_samples <- group_data$Person[group_data$Tumor == "T"]
  normal_samples <- group_data$Person[group_data$Tumor == "N"]
  
  # 生成Tumor bed文件
  process_group_samples(tumor_samples, "eccDNA.Tumor.bed")
  # 生成Normal bed文件
  process_group_samples(normal_samples, "eccDNA.Normal.bed")
  
  # 处理Stage分组
  cat("处理Stage分组...\n")
  stage12_samples <- group_data$Person[group_data$Stage %in% c(1, 2)]
  stage34_samples <- group_data$Person[group_data$Stage %in% c(3, 4)]
  
  # 生成Stage bed文件
  process_group_samples(stage12_samples, "eccDNA.stage12.bed")
  process_group_samples(stage34_samples, "eccDNA.stage34.bed")
  
  # 处理Metastasis分组
  cat("处理Metastasis分组...\n")
  metastasis_samples <- group_data$Person[group_data$Metastasis == "M" & !is.na(group_data$Metastasis) & !is.na(group_data$Person)]
  primary_samples <- group_data$Person[group_data$Metastasis == "T"& !is.na(group_data$Metastasis) & !is.na(group_data$Person)]
 
  # 生成Metastasis bed文件
  process_group_samples(metastasis_samples, "eccDNA.Metastasis.bed")
  process_group_samples(primary_samples, "eccDNA.Primary.bed")
  
  cat("所有bed文件已生成到figures/ChrCompare目录\n")
}

# 函数：处理指定分组的样本并生成bed文件
process_group_samples <- function(samples, output_file) {
  if (length(samples) == 0) {
    cat("警告: 没有找到样本，跳过", output_file, "的生成\n")
    return()
  }
  
  # 合并所有样本的数据
  all_sample_data <- data.frame()
  
  for (sample_name in samples) {
    rds_file <- file.path("rds_data", paste0(sample_name, ".rds"))
    
    if (file.exists(rds_file)) {
      cat("处理样本:", sample_name, "\n")
      sample_data <- readRDS(rds_file)
      
      # 检查必要的列是否存在
      if (all(c("X.Chr", "Start", "End") %in% colnames(sample_data))) {
        # 选择需要的列
        sample_data_subset <- sample_data %>% dplyr::select(X.Chr, Start, End)
        all_sample_data <- rbind(all_sample_data, sample_data_subset)
      } else {
        cat("警告: 样本", sample_name, "缺少必要的列\n")
      }
    } else {
      cat("警告: RDS文件不存在:", rds_file, "\n")
    }
  }
  
  # 如果有数据，保存为bed文件
  if (nrow(all_sample_data) > 0) {
    output_path <- file.path("figures/ChrCompare", output_file)
    write.table(all_sample_data, output_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    cat("已保存:", output_path, " (包含", nrow(all_sample_data), "条记录)\n")
  } else {
    cat("警告: 没有找到有效数据，跳过", output_file, "的生成\n")
  }
}

# 调用函数生成bed文件
generate_group_bed_files()





# 开始生成circos图，使用顺序结构便于后续调整和验证
cat("\n开始生成circos图...\n")

# 确保输出目录存在
output_dir <- "figures/ChrCompare"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 加载必要的包
if (!require(circlize)) install.packages("circlize")
if (!require(tidyverse)) install.packages("tidyverse")
library(circlize)
library(tidyverse)

# 读取染色体信息
cat("读取染色体信息...\n")

# 首先读取hg18_karyotype_data.txt获取染色体基本信息
if (!file.exists("hg18_karyotype_data.txt")) {
  stop("未找到hg18_karyotype_data.txt文件！")
}

# 读取染色体基本信息（5列数据）
chr_basic <- read.table("hg18_karyotype_data.txt", sep = " ", header = FALSE,
                      col.names = c("chr", "id", "start", "end", "label"),
                      stringsAsFactors = FALSE)

cat("成功读取", nrow(chr_basic), "条染色体基本信息\n")

# 然后读取cytoBand.txt获取染色体区带信息
if (!file.exists("cytoBand.txt")) {
  stop("未找到cytoBand.txt文件！")
}

# 读取染色体区带信息（5列数据）
cytoband <- read.table("cytoBand.txt", sep = "\t", header = FALSE,
                      col.names = c("chr", "start", "end", "name", "gieStain"),
                      stringsAsFactors = FALSE)

cat("成功读取", nrow(cytoband), "条染色体区带信息\n")

# 过滤常用染色体
# 确保染色体名称格式一致并只保留1-22号染色体
cytoband <- subset(cytoband, chr %in% paste0("chr", 1:22))  # 只保留1-22号染色体，不包含X/Y染色体

# 检查是否成功读取染色体信息
if (nrow(cytoband) == 0) {
  warning("未找到有效的1-22号染色体信息，请检查cytoBand.txt文件格式")
} else {
  cat("过滤后保留", nrow(cytoband), "条染色体区带信息（1-22号）\n")
}

# 将bed文件转换为密度数据的代码块
convert_bed_to_density <- function(bed_file, n_samples) {
  if (!file.exists(bed_file)) {
    cat("警告: 文件不存在:", bed_file, "\n")
    return(NULL)
  }
  
  cat("处理文件:", bed_file, "\n")
  df <- read.table(bed_file, sep = "\t", header = FALSE,
                  col.names = c("chrom", "start", "end"),
                  stringsAsFactors = FALSE)
  
  # 过滤X染色体
  df <- df %>% filter(chrom %in% paste0("chr", 1:22))
  
  # 计算1Mb窗口的密度
  df %>% 
    mutate(mb_bin = floor(start / 1e6) + 1) %>%
    group_by(chrom, mb_bin) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(
      start = (mb_bin - 1) * 1e6,
      end = mb_bin * 1e6,
      density = count / n_samples
    ) %>%
    dplyr::select(chrom, start, end, density)
}

# 读取group_save.xlsx获取样本分组信息
group_data <- readxl::read_excel("group_save.xlsx")

# 从group_data获取样本数量（添加NA过滤）
tumor_samples <- group_data$Person[group_data$Tumor == "T" & !is.na(group_data$Tumor) & !is.na(group_data$Person)]
normal_samples <- group_data$Person[group_data$Tumor == "N" & !is.na(group_data$Tumor) & !is.na(group_data$Person)]
n_tumor <- length(tumor_samples)
n_normal <- length(normal_samples)

# 转换bed文件为密度数据
tumor_density <- convert_bed_to_density(file.path(output_dir, "eccDNA.Tumor.bed"), n_tumor)
normal_density <- convert_bed_to_density(file.path(output_dir, "eccDNA.Normal.bed"), n_normal)

if (is.null(tumor_density) || is.null(normal_density)) {
  stop("无法处理bed文件，请检查文件是否存在！")
}

# 计算平均密度（用于中间文本）
mean_tumor_density <- round(mean(tumor_density$density), 1)
mean_normal_density <- round(mean(normal_density$density), 1)

cat("Tumor平均密度:", mean_tumor_density, "per Mb\n")
cat("Normal平均密度:", mean_normal_density, "per Mb\n")

# 生成Tumor vs Normal的circos图


## 条带颜色映射（灰阶 + 着丝粒红） -----------------------------------
stain_col <- c("gneg"   = "#FAFAFA",
               "gpos25" = "#D0D0D0",
               "gpos50" = "#bbb9b9",
               "gpos75" = "#8d8989",
               "gpos100"= "#797777",
               "acen"   = "#E41A1C",   # 着丝粒
               "gvar"   = "#E0E0E0",
               "stalk"  = "#E0E0E0")

# 设置pdf输出
pdf(file.path(output_dir, "circos_TumorNormal.pdf"), width = 8, height = 8)

## 4. 全局参数 ----------------------------------------------------------
circos.clear()
circos.par(start.degree = 90, gap.degree = 2.5, track.height = 0.08)

## 5. 初始化染色体 ------------------------------------------------------
circos.initialize(factors = chr_basic$chr, xlim = cbind(0, chr_basic$end))
## 7. 染色体标签（圈外） ------------------------------------------------
circos.track(track.index = 1, ylim = c(0,1), bg.border = NA,
             panel.fun = function(x, y){
               chr <- CELL_META$sector.index
               circos.text(CELL_META$xcenter, -0.4, chr,
                           facing = "outside", cex = 0.9,
                           col = "grey20", font = 2)
             })
## 6. 画 cytoband 轨道 --------------------------------------------------
circos.track(track.index = 2, ylim = c(0,1), bg.border = NA,
             panel.fun = function(x, y) {
               chr <- CELL_META$sector.index
               cb  <- cytoband[cytoband$chr == chr, ]
               for(i in seq_len(nrow(cb))){
                 circos.rect(cb$start[i], 0, cb$end[i], 1,
                             col  = stain_col[cb$gieStain[i]],
                             border = NA)
               }
             })

# 绘制Tumor密度图（更显眼）
max_tumor_density <- max(tumor_density$density) * 1.2
circos.genomicTrackPlotRegion(
  tumor_density,
  ylim = c(0, max_tumor_density),
  track.height = 0.08,  # 增加轨道高度
  bg.border = "black",  # 添加边框使其更显眼
  panel.fun = function(region, value, ...) {
    # 填充区域使其更显眼
    circos.genomicRect(region, value, ytop = value$density, ybottom = 0,
                      col = rgb(231, 27, 28, maxColorValue = 255, alpha = 200),
                      border = rgb(231, 27, 28, maxColorValue = 255))
    # 添加线条增强视觉效果
    circos.genomicLines(region, value, col = rgb(231, 27, 28, maxColorValue = 255),
                      lwd = 2)  # 增加线宽使其更显眼
  }
)

# 绘制Normal密度图（更显眼）
max_normal_density <- max(normal_density$density) * 1.2
circos.genomicTrackPlotRegion(
  normal_density,
  ylim = c(0, max_normal_density),
  track.height = 0.08,  # 增加轨道高度
  bg.border = "black",  # 添加边框使其更显眼
  panel.fun = function(region, value, ...) {
    # 填充区域使其更显眼
    circos.genomicRect(region, value, ytop = value$density, ybottom = 0,
                      col = rgb(39, 98, 174, maxColorValue = 255, alpha = 200),
                      border = rgb(39, 98, 174, maxColorValue = 255))
    # 添加线条增强视觉效果
    circos.genomicLines(region, value, col = rgb(39, 98, 174, maxColorValue = 255),
                      lwd = 2)  # 增加线宽使其更显眼
  }
)

# 在图中间添加文字
text(0, 0, 
     paste0("Mean eccDNA Density per Mb\nTumor: ", mean_tumor_density, "\nNormal: ", mean_normal_density),
     cex = 1.5,  # 增加字体大小使其更显眼
     font = 2,   # 粗体
     col = "black")

# 添加图例
legend("bottomright", 
       legend = c("Normal", "Tumor"), 
       fill = c(rgb(39, 98, 174, maxColorValue = 255), rgb(231, 27, 28, maxColorValue = 255)),
       border = "black",
       bty = "n", 
       cex = 1.5)  # 增加图例大小

circos.clear()
dev.off()





# 函数：从bed文件中筛选特定长度的eccDNA并进行统计分析
extract_large_eccdna <- function(bed_file, length_thresholds, group_name) {
  if (!file.exists(bed_file)) {
    cat("警告: 文件不存在:", bed_file, "\n")
    return(list())
  }
  
  cat("提取并统计大跨度eccDNA:", bed_file, "\n")
  df <- read.table(bed_file, sep = "\t", header = FALSE,
                  col.names = c("chrom", "start", "end"),
                  stringsAsFactors = FALSE)
  
  # 过滤X染色体并计算长度
  df <- df %>% 
    filter(chrom %in% paste0("chr", 1:22)) %>%
    mutate(length = end - start,
           group = group_name)
  
  # 按长度阈值分组
  result_list <- list()
  for (threshold in length_thresholds) {
    key <- paste0("gt_", threshold)
    result_list[[key]] <- df %>% 
      filter(length > threshold) %>%
      dplyr::select(chrom, start, end, length, group)
    cat("  ", group_name, " - 长度>", threshold, "bp的eccDNA数量:", nrow(result_list[[key]]), "\n")
  }
  
  # 返回完整数据用于统计
  result_list$all_data <- df
  
  return(result_list)
}

# 函数：进行Tumor vs Normal的长度分布统计和差异分析
analyze_length_differences <- function(tumor_data, normal_data) {
  cat("\n进行Tumor vs Normal长度分布差异分析...\n")
  
  # 合并数据
  all_data <- rbind(tumor_data$all_data, normal_data$all_data)
  
  # 基本统计
  tumor_stats <- tumor_data$all_data %>% 
    summarise(mean_length = mean(length),
              median_length = median(length),
              max_length = max(length),
              n_total = n())
  
  normal_stats <- normal_data$all_data %>% 
    summarise(mean_length = mean(length),
              median_length = median(length),
              max_length = max(length),
              n_total = n())
  
  cat("Tumor统计: 平均长度=", round(tumor_stats$mean_length, 1), 
      ", 中位数=", tumor_stats$median_length, 
      ", 最大长度=", tumor_stats$max_length, 
      ", 总数=", tumor_stats$n_total, "\n")
  cat("Normal统计: 平均长度=", round(normal_stats$mean_length, 1), 
      ", 中位数=", normal_stats$median_length, 
      ", 最大长度=", normal_stats$max_length, 
      ", 总数=", normal_stats$n_total, "\n")
  
  # 进行Wilcoxon检验
  if (nrow(tumor_data$all_data) > 0 && nrow(normal_data$all_data) > 0) {
    wilcox_test <- wilcox.test(tumor_data$all_data$length, normal_data$all_data$length)
    cat("Wilcoxon检验p值:", wilcox_test$p.value, "\n")
    
    # 计算不同长度区间的比例差异
    length_categories <- c(1000, 5000, 10000, 50000, 100000, 500000, 1000000)
    
    cat("\n长度区间分布差异:\n")
    for (i in 1:(length(length_categories) - 1)) {
      lower <- length_categories[i]
      upper <- length_categories[i + 1]
      
      tumor_count <- sum(tumor_data$all_data$length > lower & tumor_data$all_data$length <= upper)
      normal_count <- sum(normal_data$all_data$length > lower & normal_data$all_data$length <= upper)
      
      tumor_prop <- tumor_count / tumor_stats$n_total
      normal_prop <- normal_count / normal_stats$n_total
      
      cat(sprintf("  %d-%d bp: Tumor=%.2f%%, Normal=%.2f%%, 差异=%.2f%%\n", 
                  lower, upper, tumor_prop*100, normal_prop*100, (tumor_prop-normal_prop)*100))
    }
    
    # 找出差异最显著的长度区间
    diff_thresholds <- c(10000, 50000, 100000, 500000)
    best_threshold <- NULL
    max_diff <- 0
    
    for (threshold in diff_thresholds) {
      tumor_gt <- sum(tumor_data$all_data$length > threshold) / tumor_stats$n_total
      normal_gt <- sum(normal_data$all_data$length > threshold) / normal_stats$n_total
      diff <- abs(tumor_gt - normal_gt)
      
      if (diff > max_diff && tumor_gt > 0 && normal_gt > 0) {
        max_diff <- diff
        best_threshold <- threshold
      }
    }
    
    if (!is.null(best_threshold)) {
      cat("\n推荐使用的长度阈值:", best_threshold, "bp (差异最大: ", round(max_diff*100, 1), "%)\n")
      return(best_threshold)
    }
  }
  
  # 如果没有找到合适的阈值，返回默认值
  return(100000)
}

# 基于实际数据的智能阈值选择
length_thresholds <- c(50000000)  # 基于实际分布调整 50Mb

# 提取Tumor和Normal中的大跨度eccDNA
tumor_large_eccdna <- extract_large_eccdna(file.path(output_dir, "eccDNA.Tumor.bed"), length_thresholds, "Tumor")
normal_large_eccdna <- extract_large_eccdna(file.path(output_dir, "eccDNA.Normal.bed"), length_thresholds, "Normal")



chr_basic=chr_basic %>% dplyr::filter(chr %in% paste0("chr", 1:22))

# 生成包含大跨度eccDNA弧线的circos图
pdf(file.path(output_dir, "circos_TumorNormal_large_eccdna.pdf"), width = 10, height = 10)

# 全局参数设置
circos.clear()
circos.par(start.degree = 90, gap.degree = 2.5, track.height = 0.08)

# 初始化染色体
circos.initialize(factors = chr_basic$chr, xlim = cbind(0, chr_basic$end))

# 染色体标签
circos.track(track.index = 1, ylim = c(0,1), bg.border = NA,
             panel.fun = function(x, y){
               chr <- CELL_META$sector.index
               circos.text(CELL_META$xcenter, -0.4, chr,
                           facing = "outside", cex = 1.2,
                           col = "grey20", font = 2)
             })

# 画cytoband轨道
circos.track(track.index = 2, ylim = c(0,1), bg.border = NA,
             panel.fun = function(x, y) {
               chr <- CELL_META$sector.index
               cb <- cytoband[cytoband$chr == chr, ]
               for(i in seq_len(nrow(cb))){
                 circos.rect(cb$start[i], 0, cb$end[i], 1,
                             col = stain_col[cb$gieStain[i]],
                             border = NA)
               }
             })

# 绘制Tumor密度图
max_tumor_density <- max(tumor_density$density) * 1.2
circos.genomicTrackPlotRegion(
  tumor_density,
  ylim = c(0, max_tumor_density),
  track.height = 0.08,
  bg.border = "black",
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, ytop = value$density, ybottom = 0,
                      col = rgb(231, 27, 28, maxColorValue = 255, alpha = 200),
                      border = rgb(231, 27, 28, maxColorValue = 255))
    circos.genomicLines(region, value, col = rgb(231, 27, 28, maxColorValue = 255),
                      lwd = 2)
  }
)

# 绘制Normal密度图
max_normal_density <- max(normal_density$density) * 1.2
circos.genomicTrackPlotRegion(
  normal_density,
  ylim = c(0, max_normal_density),
  track.height = 0.08,
  bg.border = "black",
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, ytop = value$density, ybottom = 0,
                      col = rgb(39, 98, 174, maxColorValue = 255, alpha = 200),
                      border = rgb(39, 98, 174, maxColorValue = 255))
    circos.genomicLines(region, value, col = rgb(39, 98, 174, maxColorValue = 255),
                      lwd = 2)
  }
)

# 绘制大跨度eccDNA弧线，基于实际数据调整抽样策略
cat("\n绘制大跨度eccDNA弧线...\n")

# 函数：简单绘制所有>1Mb的eccDNA弧线（不采样）
smart_draw_links <- function(eccdna_data, threshold, color, lwd) {
  key <- paste0("gt_", threshold)

  # 检查数据是否存在
  if (!key %in% names(eccdna_data)) {
    cat("警告: 在eccdna_data中找不到键", key, "\n")
    return()
  }
  
  # 不进行采样，显示所有数据以体现tumor比normal多
  data_to_draw <- eccdna_data[[key]]
  
  cat("  绘制", nrow(data_to_draw), "条弧线\n")
  
  if (nrow(data_to_draw) > 0) {
    for (i in seq_len(nrow(data_to_draw))) {
      chr <- data_to_draw$chrom[i]
      start <- data_to_draw$start[i]
      end <- data_to_draw$end[i]
      
      # 修复：调整h.ratio使弧线更明显，使用正确的RGB格式
      circos.link(chr, start, chr, end,
                 col = color,  # 颜色已经包含了alpha值
                 lwd = lwd,
                 h.ratio = 0.8)  # 减小h.ratio使弧线更明显
    }
  }
}

# 绘制不同长度阈值的弧线，基于实际数据调整
cat("\n=== 弧线绘制策略 ===\n")

# 只保留1mb以上的eccDNA弧线 - 超长eccDNA，重点显示
# 修复：使用正确的RGB格式（alpha在0-1之间）
smart_draw_links(tumor_large_eccdna, "5e+07", 
                rgb(231, 27, 28, maxColorValue = 255), 
                lwd = 0.8)  # 增加线宽使其更明显
smart_draw_links(normal_large_eccdna, "5e+07", 
                rgb(39, 98, 174, maxColorValue = 255), 
                lwd = 0.8)  # 增加线宽使其更明显

# 使用前面已经计算好的平均密度
# tumor_density和normal_density已经在前面计算过了
# n_tumor和n_normal也已经在前面计算过了

# 计算各种长度阈值的数量
count_1mb_tumor <- nrow(tumor_large_eccdna$`gt_5e+07`)
count_1mb_normal <- nrow(normal_large_eccdna$`gt_5e+07`)

# 在图中间添加详细的统计信息
text(0, 0, 
     paste0("Mean eccDNA Density per Mb\n",
            "Tumor: ", mean_tumor_density, "\n",
            "Normal: ", mean_normal_density, "\n"),
     cex = 1.0,
     font = 2,
     col = "black")

# 添加简化的图例
legend("bottomright", 
       legend = c("Normal density", "Tumor density", 
                 paste0("Normal >50Mb (", count_1mb_normal, ") "),
                 paste0("Tumor >50Mb (", count_1mb_tumor, ") ")), 
       fill = c(rgb(39, 98, 174, maxColorValue = 255, alpha = 200),
               rgb(231, 27, 28, maxColorValue = 255, alpha = 200),
               rgb(39, 98, 174, maxColorValue = 255, alpha = 150),
               rgb(231, 27, 28, maxColorValue = 255, alpha = 150)),
       border = rep("black", 1),
       lwd = c(1, 1, 2, 2),
       bty = "n", 
       cex = 0.8)

circos.clear()
dev.off()

#-----------------------------------------------


#-------------------------------------------

# 处理Stage分组
# 基于实际数据的智能阈值选择
length_thresholds <- c(5000000)  # 基于实际分布调整

# 提取Tumor和Normal中的大跨度eccDNA
tumor_large_eccdna <- extract_large_eccdna(file.path(output_dir, "eccDNA.Tumor.bed"), length_thresholds, "Tumor")
normal_large_eccdna <- extract_large_eccdna(file.path(output_dir, "eccDNA.Normal.bed"), length_thresholds, "Normal")
chr_basic=chr_basic %>% dplyr::filter(chr %in% paste0("chr", 1:22))

# 生成包含大跨度eccDNA弧线的circos图 (Tumor vs Normal)
pdf(file.path(output_dir, "circos_TumorNormal_large_eccdna.pdf"), width = 10, height = 10)

# 全局参数设置
circos.clear()
circos.par(start.degree = 90, gap.degree = 2.5, track.height = 0.08)

# 初始化染色体
circos.initialize(factors = chr_basic$chr, xlim = cbind(0, chr_basic$end))

# 染色体标签
circos.track(track.index = 1, ylim = c(0,1), bg.border = NA,
             panel.fun = function(x, y){
               chr <- CELL_META$sector.index
               circos.text(CELL_META$xcenter, -0.4, chr,
                           facing = "outside", cex = 1.2,
                           col = "grey20", font = 2)
             })

# 画cytoband轨道
circos.track(track.index = 2, ylim = c(0,1), bg.border = NA,
             panel.fun = function(x, y) {
               chr <- CELL_META$sector.index
               cb <- cytoband[cytoband$chr == chr, ]
               for(i in seq_len(nrow(cb))){  
                 circos.rect(cb$start[i], 0, cb$end[i], 1,
                             col = stain_col[cb$gieStain[i]],
                             border = NA)
               }
             })

# 绘制Tumor密度图
max_tumor_density <- max(tumor_density$density) * 1.2
circos.genomicTrackPlotRegion(
  tumor_density,
  ylim = c(0, max_tumor_density),
  track.height = 0.08,
  bg.border = "black",
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, ytop = value$density, ybottom = 0,
                      col = rgb(231, 27, 28, maxColorValue = 255, alpha = 200),
                      border = rgb(231, 27, 28, maxColorValue = 255))
    circos.genomicLines(region, value, col = rgb(231, 27, 28, maxColorValue = 255),
                      lwd = 2)
  }
)

# 绘制Normal密度图
max_normal_density <- max(normal_density$density) * 1.2
circos.genomicTrackPlotRegion(
  normal_density,
  ylim = c(0, max_normal_density),
  track.height = 0.08,
  bg.border = "black",
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, ytop = value$density, ybottom = 0,
                      col = rgb(39, 98, 174, maxColorValue = 255, alpha = 200),
                      border = rgb(39, 98, 174, maxColorValue = 255))
    circos.genomicLines(region, value, col = rgb(39, 98, 174, maxColorValue = 255),
                      lwd = 2)
  }
)

# 绘制大跨度eccDNA弧线，基于实际数据调整抽样策略
cat("\n绘制大跨度eccDNA弧线...\n")

# 函数：绘制eccDNA弧线
smart_draw_links <- function(eccdna_data, threshold, color, lwd) {
  key <- paste0("gt_", threshold)

  # 检查数据是否存在
  if (!key %in% names(eccdna_data)) {
    cat("警告: 在eccdna_data中找不到键", key, "\n")
    return()
  }
  
  data_to_draw <- eccdna_data[[key]]
  
  cat("  绘制", nrow(data_to_draw), "条弧线\n")
  
  if (nrow(data_to_draw) > 0) {
    for (i in seq_len(nrow(data_to_draw))) {
      chr <- data_to_draw$chrom[i]
      start <- data_to_draw$start[i]
      end <- data_to_draw$end[i]
      
      circos.link(chr, start, chr, end,
                 col = color,
                 border = NA,
                 lwd = lwd,
                 h.ratio = 0.6)
    }
  }
}

# 绘制Tumor和Normal的弧线
smart_draw_links(tumor_large_eccdna, "5e+06", 
                rgb(231, 27, 28, maxColorValue = 255, alpha = 0.7), 
                lwd = 1.5)
smart_draw_links(normal_large_eccdna, "5e+06", 
                rgb(39, 98, 174, maxColorValue = 255, alpha = 0.7), 
                lwd = 1.5)

# 计算各种长度阈值的数量
count_1mb_tumor <- ifelse("gt_5e+06" %in% names(tumor_large_eccdna), nrow(tumor_large_eccdna$`gt_5e+06`), 0)
count_1mb_normal <- ifelse("gt_5e+06" %in% names(normal_large_eccdna), nrow(normal_large_eccdna$`gt_5e+06`), 0)

# 在图中间添加详细的统计信息
text(0, 0, 
     paste0("Mean eccDNA Density per Mb\n",
            "Tumor: ", mean_tumor_density, "\n",
            "Normal: ", mean_normal_density, "\n",
            "Tumor >1Mb: ", count_1mb_tumor, "\n",
            "Normal >1Mb: ", count_1mb_normal),
     cex = 1.0,
     font = 2,
     col = "black")

# 添加简化的图例
legend("bottomright", 
       legend = c("Normal density", "Tumor density", 
                 paste0("Normal >1Mb (", count_1mb_normal, ")"),
                 paste0("Tumor >1Mb (", count_1mb_tumor, ")")), 
       fill = c(rgb(39, 98, 174, maxColorValue = 255, alpha = 200),
               rgb(231, 27, 28, maxColorValue = 255, alpha = 200),
               rgb(39, 98, 174, maxColorValue = 255, alpha = 150),
               rgb(231, 27, 28, maxColorValue = 255, alpha = 150)),
       border = rep("black", 1),
       lwd = c(1, 1, 2, 2),
       bty = "n", 
       cex = 0.8)

circos.clear()
dev.off()

# 处理Stage分组
stage12_samples <- group_data$Person[group_data$Stage %in% c(1, 2) & !is.na(group_data$Stage) & !is.na(group_data$Person)]
stage34_samples <- group_data$Person[group_data$Stage %in% c(3, 4) & !is.na(group_data$Stage) & !is.na(group_data$Person)]
n_stage12 <- length(stage12_samples)
n_stage34 <- length(stage34_samples)

stage12_density <- convert_bed_to_density(file.path(output_dir, "eccDNA.stage12.bed"), n_stage12)
stage34_density <- convert_bed_to_density(file.path(output_dir, "eccDNA.stage34.bed"), n_stage34)


  mean_stage12_density <- round(mean(stage12_density$density), 1)
  mean_stage34_density <- round(mean(stage34_density$density), 1)
  
  # 提取Stage分组中的大跨度eccDNA
  length_thresholds <- c(50000000)  # 基于实际分布调整 50Mb
  stage12_large_eccdna <- extract_large_eccdna(file.path(output_dir, "eccDNA.stage12.bed"), length_thresholds, "Stage12")
  stage34_large_eccdna <- extract_large_eccdna(file.path(output_dir, "eccDNA.stage34.bed"), length_thresholds, "Stage34")
  
  # 生成Stage分组的circos图
  cat("\n生成Stage分组circos图...\n")
  pdf(file.path(output_dir, "circos_Stage_large_eccdna.pdf"), width = 10, height = 10)
  
  # 全局参数设置
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 2.5, track.height = 0.08)
  
  # 初始化染色体
  circos.initialize(factors = chr_basic$chr, xlim = cbind(0, chr_basic$end))
  
  # 染色体标签
  circos.track(track.index = 1, ylim = c(0,1), bg.border = NA,
               panel.fun = function(x, y){
                 chr <- CELL_META$sector.index
                 circos.text(CELL_META$xcenter, -0.4, chr,
                             facing = "outside", cex = 1.2,
                             col = "grey20", font = 2)
               })
  
  # 画cytoband轨道
  circos.track(track.index = 2, ylim = c(0,1), bg.border = NA,
               panel.fun = function(x, y) {
                 chr <- CELL_META$sector.index
                 cb <- cytoband[cytoband$chr == chr, ]
                 for(i in seq_len(nrow(cb))){  
                   circos.rect(cb$start[i], 0, cb$end[i], 1,
                               col = stain_col[cb$gieStain[i]],
                               border = NA)
                 }
               })
  
  # 定义颜色（符合Cell Report高级配色）
  color_stage12 <- rgb(114, 158, 206, maxColorValue = 255)  # 天蓝色
  color_stage34 <- rgb(244, 109, 67, maxColorValue = 255)   # 橙色
  
   # 绘制Stage 3-4密度图
  max_stage34_density <- max(stage34_density$density) * 1.2
  circos.genomicTrackPlotRegion(
    stage34_density,
    ylim = c(0, max_stage34_density),
    track.height = 0.08,
    bg.border = "black",
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, ytop = value$density, ybottom = 0,
                        col = rgb(244, 109, 67, maxColorValue = 255, alpha = 200),
                        border = color_stage34)
      circos.genomicLines(region, value, col = color_stage34,
                        lwd = 2)
    }
  )
  # 绘制Stage 1-2密度图
  max_stage12_density <- max(stage12_density$density) * 1.2
  circos.genomicTrackPlotRegion(
    stage12_density,
    ylim = c(0, max_stage12_density),
    track.height = 0.08,
    bg.border = "black",
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, ytop = value$density, ybottom = 0,
                        col = rgb(114, 158, 206, maxColorValue = 255, alpha = 200),
                        border = color_stage12)
      circos.genomicLines(region, value, col = color_stage12,
                        lwd = 2)
    }
  )
  
 
  
  # 绘制大跨度eccDNA弧线
  cat("绘制Stage分组大跨度eccDNA弧线...\n")
    smart_draw_links(stage34_large_eccdna, "5e+07", 
                  rgb(244, 109, 67, maxColorValue = 255), 
                  lwd = 0.8)
  smart_draw_links(stage12_large_eccdna, "5e+07", 
                  rgb(114, 158, 206, maxColorValue = 255), 
                  lwd = 0.8)

  
  # 计算各种长度阈值的数量
  count_1mb_stage12 <- ifelse("gt_5e+07" %in% names(stage12_large_eccdna), nrow(stage12_large_eccdna$`gt_5e+07`), 0)
  count_1mb_stage34 <- ifelse("gt_5e+07" %in% names(stage34_large_eccdna), nrow(stage34_large_eccdna$`gt_5e+07`), 0)
  
  # 在图中间添加详细的统计信息
  text(0, 0, 
       paste0("Mean eccDNA Density per Mb\n",
              "Stage 1-2: ", mean_stage12_density, "\n",
              "Stage 3-4: ", mean_stage34_density, "\n"
             ),
       cex = 1.0,
       font = 2,
       col = "black")
  
  # 添加简化的图例
  legend("bottomright", 
         legend = c("Stage 1-2 density", "Stage 3-4 density", 
                   paste0("Stage 1-2 >50Mb (", count_1mb_stage12, ")"),
                   paste0("Stage 3-4 >50Mb (", count_1mb_stage34, ")")), 
         fill = c(rgb(114, 158, 206, maxColorValue = 255, alpha = 200),
                 rgb(244, 109, 67, maxColorValue = 255, alpha = 200),
                 rgb(114, 158, 206, maxColorValue = 255, alpha = 150),
                 rgb(244, 109, 67, maxColorValue = 255, alpha = 150)),
         border = rep("black", 1),
         lwd = c(1, 1, 2, 2),
         bty = "n", 
         cex = 0.8)
  
  circos.clear()
  dev.off()


# 处理Metastasis分组
metastasis_samples <- group_data$Person[group_data$Metastasis == "M" & !is.na(group_data$Metastasis) & !is.na(group_data$Person)]
primary_samples <- group_data$Person[group_data$Metastasis == "T" & !is.na(group_data$Metastasis) & !is.na(group_data$Person)]
n_metastasis <- length(metastasis_samples)
n_primary <- length(primary_samples)

metastasis_density <- convert_bed_to_density(file.path(output_dir, "eccDNA.Metastasis.bed"), n_metastasis)
primary_density <- convert_bed_to_density(file.path(output_dir, "eccDNA.Primary.bed"), n_primary)


  mean_metastasis_density <- round(mean(metastasis_density$density), 1)
  mean_primary_density <- round(mean(primary_density$density), 1)
  
  # 提取Metastasis分组中的大跨度eccDNA
  metastasis_large_eccdna <- extract_large_eccdna(file.path(output_dir, "eccDNA.Metastasis.bed"), length_thresholds, "Metastasis")
  primary_large_eccdna <- extract_large_eccdna(file.path(output_dir, "eccDNA.Primary.bed"), length_thresholds, "Primary")
  
  # 生成Metastasis分组的circos图
  cat("\n生成Metastasis分组circos图...\n")
  pdf(file.path(output_dir, "circos_Metastasis_large_eccdna.pdf"), width = 10, height = 10)
  
  # 全局参数设置
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 2.5, track.height = 0.08)
  
  # 初始化染色体
  circos.initialize(factors = chr_basic$chr, xlim = cbind(0, chr_basic$end))
  
  # 染色体标签
  circos.track(track.index = 1, ylim = c(0,1), bg.border = NA,
               panel.fun = function(x, y){
                 chr <- CELL_META$sector.index
                 circos.text(CELL_META$xcenter, -0.4, chr,
                             facing = "outside", cex = 1.2,
                             col = "grey20", font = 2)
               })
  
  # 画cytoband轨道
  circos.track(track.index = 2, ylim = c(0,1), bg.border = NA,
               panel.fun = function(x, y) {
                 chr <- CELL_META$sector.index
                 cb <- cytoband[cytoband$chr == chr, ]
                 for(i in seq_len(nrow(cb))){  
                   circos.rect(cb$start[i], 0, cb$end[i], 1,
                               col = stain_col[cb$gieStain[i]],
                               border = NA)
                 }
               })
  
  # 定义颜色（符合Cell Report高级配色）
  color_primary <- rgb(65, 171, 93, maxColorValue = 255)      # 深绿色
  color_metastasis <- rgb(153, 102, 255, maxColorValue = 255)  # 紫色
  

   
  # 绘制Metastasis密度图
  max_metastasis_density <- max(metastasis_density$density) * 1.2
  circos.genomicTrackPlotRegion(
    metastasis_density,
    ylim = c(0, max_metastasis_density),
    track.height = 0.08,
    bg.border = "black",
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, ytop = value$density, ybottom = 0,
                        col = rgb(153, 102, 255, maxColorValue = 255, alpha = 200),
                        border = color_metastasis)
      circos.genomicLines(region, value, col = color_metastasis,
                        lwd = 2)
    }
  )
  # 绘制Primary密度图
  max_primary_density <- max(primary_density$density) * 1.2
  circos.genomicTrackPlotRegion(
    primary_density,
    ylim = c(0, max_primary_density),
    track.height = 0.08,
    bg.border = "black",
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, ytop = value$density, ybottom = 0,
                        col = rgb(65, 171, 93, maxColorValue = 255, alpha = 200),
                        border = color_primary)
      circos.genomicLines(region, value, col = color_primary,
                        lwd = 2)
    }
  )
 
  
  # 绘制大跨度eccDNA弧线
  cat("绘制Metastasis分组大跨度eccDNA弧线...\n")
    smart_draw_links(metastasis_large_eccdna, "5e+07", 
                  rgb(153, 102, 255, maxColorValue = 255), 
                  lwd = 0.8)
  
  smart_draw_links(primary_large_eccdna, "5e+07", 
                  rgb(65, 171, 93, maxColorValue = 255), 
                  lwd = 0.8)

  # 计算各种长度阈值的数量
  count_1mb_primary <- ifelse("gt_5e+07" %in% names(primary_large_eccdna), nrow(primary_large_eccdna$`gt_5e+07`), 0)
  count_1mb_metastasis <- ifelse("gt_5e+07" %in% names(metastasis_large_eccdna), nrow(metastasis_large_eccdna$`gt_5e+07`), 0)
  
  # 在图中间添加详细的统计信息
  text(0, 0, 
       paste0("Mean eccDNA Density per Mb\n",
              "Primary: ", mean_primary_density, "\n",
              "Metastasis: ", mean_metastasis_density, "\n",
              "Primary >50Mb: ", count_1mb_primary, "\n",
              "Metastasis >50Mb: ", count_1mb_metastasis),
       cex = 1.0,
       font = 2,
       col = "black")
  
  # 添加简化的图例
  legend("bottomright", 
         legend = c("Primary density", "Metastasis density", 
                   paste0("Primary >50Mb (", count_1mb_primary, ")"),
                   paste0("Metastasis >50Mb (", count_1mb_metastasis, ")")), 
         fill = c(rgb(65, 171, 93, maxColorValue = 255, alpha = 200),
                 rgb(153, 102, 255, maxColorValue = 255, alpha = 200),
                 rgb(65, 171, 93, maxColorValue = 255, alpha = 150),
                 rgb(153, 102, 255, maxColorValue = 255, alpha = 150)),
         border = rep("black", 1),
         lwd = c(1, 1, 2, 2),
         bty = "n", 
         cex = 0.8)
  
  circos.clear()
  dev.off()
