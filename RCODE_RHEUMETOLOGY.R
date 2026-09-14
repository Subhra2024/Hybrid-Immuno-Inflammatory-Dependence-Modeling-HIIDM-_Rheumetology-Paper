# ==============================================================================
# Project: Hybrid Immuno-Inflammatory Dependence Modeling (HIIDM)
# Purpose: Full Implementation - Tables, Advanced Plots & Network Visualization
# ==============================================================================

# ------------------------------
# Load Required Libraries
# ------------------------------
suppressMessages({
  library(tidyverse)
  library(readxl)
  library(mice)
  library(energy)
  library(igraph)
  library(glmnet)
  library(corrplot)
  library(ggsci)
  library(caret)
  library(naniar)
  library(factoextra)
  library(visNetwork)
  library(GGally)
  library(pheatmap)
  library(reshape2)
  library(viridis)
  library(ggraph)
  library(tidygraph)
  library(tibble)
})

# ==============================================================================
# SECTION 1: DATA LOADING & PREPROCESSING
# ==============================================================================
df_raw <- read_excel("C:/Users/SUBHRAJIT SAHA/Downloads/Rheumatic and Autoimmune Disease Dataset.xlsx")
names(df_raw) <- make.names(names(df_raw))

biomarkers <- c("ESR","CRP","RF","Anti.CCP","HLA.B27","ANA","Anti.Ro","Anti.La",
                "Anti.dsDNA","Anti.Sm","C3","C4")
demographics <- c("Age","Gender")
binary_markers <- c("HLA.B27","ANA","Anti.Ro","Anti.La","Anti.dsDNA","Anti.Sm")

df <- df_raw %>%
  mutate(Gender = ifelse(Gender == "Male",1,0)) %>%
  mutate(across(all_of(binary_markers), ~ifelse(.=="Positive",1,0))) %>%
  mutate(Disease = as.factor(Disease))

# Missing data visualization
if(names(dev.cur())!="null device") dev.new()
gg_miss_var(df[,c(demographics,biomarkers)]) +
  theme_minimal() + labs(title="Missing Values Distribution Across Clinical Variables")

# Multiple imputation
imp <- mice(df[,c(demographics,biomarkers)], m=1, method='pmm', printFlag=FALSE)
df_complete <- complete(imp)
df_complete$Disease <- df$Disease

# Baseline table
baseline_summary <- df_complete %>%
  group_by(Disease) %>%
  summarise(
    N_Patients = n(),
    Mean_Age = round(mean(Age),1),
    Pct_Male = round(mean(Gender)*100,1),
    Mean_ESR = round(mean(ESR),1),
    Mean_CRP = round(mean(CRP),1)
  )
View(as.data.frame(baseline_summary), title="Table 1: Baseline Characteristics")
print(baseline_summary)

# Standardization
df_std <- df_complete %>%
  mutate(across(all_of(biomarkers), ~ scale(.) %>% as.numeric))
Z_mat <- as.matrix(df_std[,biomarkers])

# ==============================================================================
# SECTION 2: INFLAMMATORY SEVERITY INDEX (PCA)
# ==============================================================================
inflam_data <- df_std[,c("ESR","CRP")]
pca_inflam <- prcomp(inflam_data, center=FALSE, scale.=FALSE)
df_complete$I_i <- df_std$I_i <- pca_inflam$x[,1]

# PCA Biplot
if(names(dev.cur())!="null device") dev.new()
fviz_pca_biplot(pca_inflam, col.ind=df_complete$Disease,
                palette="npg", addEllipses=TRUE, label="var",
                title="Inflammatory Severity Index (PCA)", legend.title="Disease Group")

# ==============================================================================
# SECTION 3: COMPOSITE DEPENDENCE MATRIX & NETWORK
# ==============================================================================
set.seed(123)
sample_idx <- if(nrow(df_std)>2000) sample(1:nrow(df_std),2000) else 1:nrow(df_std)
Z_sample <- Z_mat[sample_idx, ]
p <- length(biomarkers)
T_mat <- matrix(0,p,p)
colnames(T_mat) <- rownames(T_mat) <- biomarkers

for(j in 1:p){
  for(k in j:p){
    tau <- cor(Z_sample[,j],Z_sample[,k], method="kendall")
    v_dist <- dcor(Z_sample[,j],Z_sample[,k])
    val <- 0.5*abs(tau)+0.5*v_dist
    T_mat[j,k] <- T_mat[k,j] <- val
  }
}

# Heatmap
if(names(dev.cur())!="null device") dev.new()
corrplot(T_mat, method="color", type="upper", order="hclust",
         addCoef.col="black", tl.col="black", tl.srt=45,
         col=colorRampPalette(c("#f7fbff","#08306b"))(200),
         title="\nComposite Dependence Matrix (T)", mar=c(0,0,2,0))

# Adjacency & Laplacian
delta <- quantile(T_mat[upper.tri(T_mat)],0.75)
A_adj <- ifelse(T_mat>delta,1,0); diag(A_adj)<-0
D_mat <- diag(rowSums(A_adj))
L_mat <- D_mat - A_adj

# ==============================================================================
# SECTION 4: AUTOIMMUNE ACTIVITY SCORE & CENTRALITY
# ==============================================================================
g <- graph_from_adjacency_matrix(A_adj, mode="undirected")
ev <- eigen_centrality(g)$vector
theta <- ev / sum(ev)

# Centrality table
centrality_df <- data.frame(Biomarker=names(theta), Centrality_Theta=round(theta,4))
View(centrality_df, title="Table 2: Biomarker Centrality")
print(centrality_df)
# Graph-regularized score
df_complete$A_i <- Z_mat %*% theta
lambda_reg <- 0.1
network_penalty <- rowSums((Z_mat %*% L_mat)*Z_mat)
df_complete$A_tilde_i <- df_complete$A_i - (lambda_reg*network_penalty)

# Network plot
if(names(dev.cur())!="null device") dev.new()
node_colors <- c("ESR"="red","CRP"="red","RF"="blue","Anti.CCP"="blue",
                 "HLA.B27"="blue","ANA"="green","Anti.Ro"="green",
                 "Anti.La"="green","Anti.dsDNA"="green","Anti.Sm"="green",
                 "C3"="orange","C4"="orange")
plot(g, vertex.color=node_colors[names(V(g))], vertex.label.color="black",
     vertex.label.font=2, vertex.size=ev*40, edge.width=2,
     layout=layout_with_fr, main="Biomarker Network (Node Size=Centrality)")
legend("bottomleft", legend=c("Inflammatory","Autoantibodies","ANA","Complement"),
       fill=c("red","blue","green","orange"), bty="n", cex=0.8)

# Autoimmune score distributions
score_plot_data <- df_complete %>%
  select(Disease,A_i,A_tilde_i) %>%
  pivot_longer(cols=c(A_i,A_tilde_i), names_to="Score_Type", values_to="Score_Value")
if(names(dev.cur())!="null device") dev.new()
ggplot(score_plot_data, aes(x=Score_Value, fill=Disease)) +
  geom_density(alpha=0.5) + facet_wrap(~Score_Type, scales="free") +
  theme_minimal() + scale_fill_npg() +
  labs(title="Autoimmune Activity Score Distributions", x="Score Value", y="Density")

# ==============================================================================
# SECTION 5: ELASTIC-NET DISEASE PREDICTION
# ==============================================================================
set.seed(42)
train_idx <- createDataPartition(df_complete$Disease, p=0.8, list=FALSE)
train_data <- df_complete[train_idx,]; test_data <- df_complete[-train_idx,]

X_train <- model.matrix(Disease~Age+Gender+I_i+A_tilde_i,data=train_data)[,-1]
Y_train <- train_data$Disease
X_test <- model.matrix(Disease~Age+Gender+I_i+A_tilde_i,data=test_data)[,-1]
Y_test <- test_data$Disease

cv_fit <- cv.glmnet(X_train,Y_train,family="multinomial", alpha=0.5, type.measure="class")
preds_class <- predict(cv_fit,newx=X_test, s="lambda.min", type="class")
preds_factor <- factor(as.character(preds_class), levels=levels(Y_test))

# Classification metrics
conf_mat <- confusionMatrix(preds_factor,Y_test)
View(as.data.frame(conf_mat$byClass), title="Table 3: Classification Metrics")
print(conf_mat)

# Coefficient importance
coef_list <- coef(cv_fit, s="lambda.min")
imp_df <- bind_rows(lapply(names(coef_list), function(k){
  cf <- as.matrix(coef_list[[k]])[-1,]
  data.frame(Class=k, Feature=names(cf), Value=abs(cf))
})) %>% group_by(Feature) %>% summarise(Overall_Importance=sum(Value))

if(names(dev.cur())!="null device") dev.new()
ggplot(imp_df, aes(x=reorder(Feature,Overall_Importance), y=Overall_Importance, fill=Overall_Importance)) +
  geom_bar(stat="identity") + coord_flip() +
  scale_fill_gradient(low="lightblue", high="navy") +
  theme_minimal() + labs(title="Global Feature Importance", x="Predictors", y="Aggregated Magnitude")

# Confusion matrix heatmap
cm_df <- as.data.frame(conf_mat$table)
if(names(dev.cur())!="null device") dev.new()
ggplot(cm_df, aes(x=Reference, y=Prediction, fill=Freq)) +
  geom_tile() + geom_text(aes(label=Freq), color="white", size=6, fontface="bold") +
  scale_fill_gradient(low="#a6cbf3", high="#0b2a59") +
  theme_minimal() + theme(axis.text.x=element_text(angle=45,hjust=1)) +
  labs(title="Prediction Confusion Matrix", x="True Class", y="Predicted Class")

# ==============================================================================
# SECTION 6: INTERACTIVE NETWORK VISUALIZATION
# ==============================================================================
edges <- as_data_frame(g, what="edges")
nodes <- data.frame(
  id = biomarkers,
  label = biomarkers,
  value = theta*50,
  group = c(rep("Inflammatory",2), rep("Autoantibodies",3), rep("ANA",4), rep("Complement",3))
)
visNetwork(nodes, edges, width="100%") %>%
  visGroups(groupname="Inflammatory", color="red") %>%
  visGroups(groupname="Autoantibodies", color="blue") %>%
  visGroups(groupname="ANA", color="green") %>%
  visGroups(groupname="Complement", color="orange") %>%
  visOptions(highlightNearest=TRUE, nodesIdSelection=TRUE) %>%
  visLayout(randomSeed=123)

# ==============================================================================
# SECTION 7: ADDITIONAL TABLES & ADVANCED PLOTS
# ==============================================================================

# ---- 7a: Biomarker Summary by Disease ----
biomarker_summary <- df_complete %>%
  group_by(Disease) %>%
  summarise(across(all_of(biomarkers),
                   list(mean=~mean(., na.rm=TRUE),
                        sd=~sd(., na.rm=TRUE),
                        median=~median(., na.rm=TRUE),
                        missing=~sum(is.na(.))),
                   .names="{col}_{fn}"))
View(as.data.frame(biomarker_summary), title="Table: Biomarker Summary by Disease")
print(biomarker_summary)
# ---- 7b: Pairwise biomarker relationships ----
if(names(dev.cur())!="null device") dev.new()
ggpairs(df_complete[,c(biomarkers,"Disease")],
        mapping=aes(color=Disease, alpha=0.6)) +
  theme_minimal() + ggtitle("Pairwise Biomarker Relationships by Disease")

# ---- 7c: Composite dependence heatmap ----
if(names(dev.cur())!="null device") dev.new()
pheatmap(T_mat, cluster_rows=TRUE, cluster_cols=TRUE,
         color=colorRampPalette(c("white","navy"))(100),
         main="Hierarchical Clustering Heatmap of Composite Dependence")

# ---- 7d: Autoimmune score vs inflammatory index ----
if(names(dev.cur())!="null device") dev.new()
ggplot(df_complete, aes(x=I_i, y=A_tilde_i, color=Disease,
                        size=apply(Z_mat,1,function(x) max(x*theta)))) +
  geom_point(alpha=0.7) + scale_color_npg() +
  labs(title="Autoimmune Score vs Inflammatory Severity Index",
       x="Inflammatory Index (I_i)", y="Graph-Regularized Autoimmune Score (A~i)") +
  theme_minimal()

# ---- 7e: Elastic-Net coefficient heatmap ----
coef_df <- bind_rows(lapply(names(coef_list), function(k){
  cf <- as.matrix(coef_list[[k]])[-1,]
  data.frame(Class=k, Feature=names(cf), Value=cf)
}))
coef_matrix <- acast(coef_df, Feature ~ Class, value.var="Value")
if(names(dev.cur())!="null device") dev.new()
pheatmap(coef_matrix, cluster_rows=TRUE, cluster_cols=TRUE,
         color=viridis(50), main="Elastic-Net Coefficient Heatmap")

# ---- 7f: Advanced biomarker network (ggraph) ----
edge_list <- which(A_adj==1, arr.ind=TRUE)
edges <- data.frame(
  from=rownames(T_mat)[edge_list[,1]],
  to=colnames(T_mat)[edge_list[,2]],
  weight=T_mat[edge_list]
)
nodes <- data.frame(
  name=biomarkers,
  type=c(rep("Inflammatory",2), rep("Autoantibodies",3), rep("ANA",4), rep("Complement",3)),
  theta=theta
)
g_tbl <- tbl_graph(nodes=nodes, edges=edges, directed=FALSE)
if(names(dev.cur())!="null device") dev.new()
ggraph(g_tbl, layout='fr') +
  geom_edge_link(aes(width=weight), alpha=0.6, color="gray50") +
  geom_node_point(aes(size=theta, color=type)) +
  geom_node_text(aes(label=name), repel=TRUE, size=4) +
  scale_size_continuous(range=c(5,15)) +
  scale_color_manual(values=c("Inflammatory"="red","Autoantibodies"="blue",
                              "ANA"="green","Complement"="orange")) +
  theme_void() + ggtitle("Advanced Biomarker Network (Weighted & Colored)")

# ==============================================================================
# SECTION 8: ADDITIONAL HIIDM TABLES
# ==============================================================================

# 1. Inflammatory Severity Index by Disease
inflammatory_table <- df_complete %>%
  group_by(Disease) %>%
  summarise(
    Mean_ESR = round(mean(ESR, na.rm=TRUE),2),
    SD_ESR   = round(sd(ESR, na.rm=TRUE),2),
    Mean_CRP = round(mean(CRP, na.rm=TRUE),2),
    SD_CRP   = round(sd(CRP, na.rm=TRUE),2),
    Mean_Ii  = round(mean(I_i, na.rm=TRUE),2),
    SD_Ii    = round(sd(I_i, na.rm=TRUE),2)
  )
View(as.data.frame(inflammatory_table), title="Table: Inflammatory Severity Index by Disease")
print(inflammatory_table)
# 2. Composite Dependence Summary
T_summary <- data.frame(
  Biomarker = biomarkers,
  Avg_Dependence = round(colMeans(T_mat), 3),
  Max_Dependence = round(apply(T_mat, 2, max), 3),
  Min_Dependence = round(apply(T_mat, 2, min), 3)
)
View(T_summary, title="Table: Composite Dependence Summary per Biomarker")
print(T_summary)
# 3. Network Degree & Laplacian Contribution
# Per-biomarker Laplacian contribution
laplacian_contrib <- colSums(Z_mat * (Z_mat %*% L_mat))

# Network degree
degree_vec <- rowSums(A_adj)

# Create table
network_table <- data.frame(
  Biomarker = biomarkers,
  Degree = degree_vec,
  Laplacian_Contribution = round(laplacian_contrib, 2),
  Centrality_Theta = round(theta, 4)
)

View(network_table, title="Table: Biomarker Network Metrics")
print(network_table)
# 4. Autoimmune Activity Scores Summary
autoimmune_table <- df_complete %>%
  group_by(Disease) %>%
  summarise(
    Mean_Ai = round(mean(A_i, na.rm=TRUE),2),
    SD_Ai   = round(sd(A_i, na.rm=TRUE),2),
    Mean_Atilde = round(mean(A_tilde_i, na.rm=TRUE),2),
    SD_Atilde   = round(sd(A_tilde_i, na.rm=TRUE),2)
  )
View(as.data.frame(autoimmune_table), title="Table: Autoimmune Activity Scores by Disease")
print(autoimmune_table)
# 5. Predictors Correlation Table
predictor_table <- round(cor(df_complete[, c("Age","Gender","I_i","A_i","A_tilde_i")]), 3)
View(predictor_table, title="Table: Correlation of Predictors Used in Disease Model")
print(predictor_table)
# 6. Elastic-Net Coefficients per Class
elasticnet_coef_table <- bind_rows(lapply(names(coef_list), function(cls){
  cf <- as.matrix(coef_list[[cls]])
  if(nrow(cf) <= 1) return(NULL)  # only intercept, skip
  cf <- cf[-1, , drop=FALSE]      # remove intercept
  data.frame(
    Class = cls,
    Feature = rownames(cf),
    Coefficient = round(cf[,1], 4)
  )
}))
View(elasticnet_coef_table, title="Table: Elastic-Net Coefficients per Disease Class")
print(elasticnet_coef_table)
# 7. Aggregated Feature Importance
feature_importance <- elasticnet_coef_table %>%
  group_by(Feature) %>%
  summarise(Overall_Importance = sum(abs(Coefficient))) %>%
  arrange(desc(Overall_Importance))
View(feature_importance, title="Table: Aggregated Feature Importance Across Classes")
print(feature_importance)
# 8. Confusion Matrix Table
confusion_table <- as.data.frame(conf_mat$table) %>%
  rename(True_Class = Reference, Predicted_Class = Prediction, Count = Freq)
View(confusion_table, title="Table: Confusion Matrix of Predicted vs True Disease Class")
print(confusion_table)
# 9. Per-Class Classification Metrics
class_metrics <- as.data.frame(conf_mat$byClass) %>%
  round(3) %>%
  tibble::rownames_to_column(var = "Metric")
View(class_metrics, title="Table: Per-Class Classification Metrics (Precision, Recall, F1)")
print(class_metrics)
cat("\nAll analyses, tables, and plots (including additional HIIDM tables) are complete.\n")
