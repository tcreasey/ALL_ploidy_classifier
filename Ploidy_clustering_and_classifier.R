library(FSelector)
library(BBmisc)
library(factoextra)
library(ComplexHeatmap)
library(ggplot2)
library(data.table)
library(rpart)
library(rattle)
library(caret)
library(rpart.plot)

# Input dataframe (discovery_cohort.csv) available as supplementary table 4 in accompanying manuscript

d <- read.csv("discovery_cohort.csv")


# Only keep HeH and HoTr cases based on Final diagnosis 
d_HoTr_HeH <- d[d$Final.diagnosis == "HoTr" | d$Final.diagnosis == "HeH", ]
d_HoTr_HeH <- na.omit(d_HoTr_HeH)
colnames(d_HoTr_HeH) <- colnames(d)

# Normalise across patient samples
d_HoTr_HeH_norm <- normalize(as.matrix(d_HoTr_HeH[,c(5:26)]), method = "standardize")
d_HoTr_HeH_norm <- cbind(d_HoTr_HeH[c(1:4)], d_HoTr_HeH_norm)
colnames(d_HoTr_HeH_norm) <- colnames(d_HoTr_HeH)


# Information gain of individual chromosomes for Final diagnosis HoTr vs HeH-------------------------------------
d3 <- subset(d_HoTr_HeH_norm, select = c(4:26))
Info_gain <- as.data.frame(information.gain(Final.diagnosis~., d3))


# PCA -----------------------------------------------------------------------------------------------------------

#PCA based on ploidy status**************************************************************************************

# calculate the PCA (with scaled data)
wdbc.pr.final <- prcomp(d_HoTr_HeH_norm[c(5:26)], center = TRUE, scale = TRUE)
# print the details of PCA
summary(wdbc.pr.final)


lev <- c("Low hypodiploidy","Near triploidy","HoTr by SNP","High hyperdiploidy")

p1 <- fviz_pca_ind(wdbc.pr.final, geom.ind = "point", pointshape = 21, 
                   pointsize = 2, 
                   fill.ind = d_HoTr_HeH_norm$Cytogenetics, 
                   col.var = "black", 
                   addEllipses = FALSE,
                   label = "var",
                   repel = TRUE) + ggtitle(NULL) +
  scale_fill_manual(name = "Cytogenetic subgroup",
                    values = c("High hyperdiploidy" = "dodgerblue2", 
                               "Low hypodiploidy" = "red", 
                               "Near triploidy" = "green3",
                               "HoTr by SNP"="yellow"),
                    breaks = lev) + 
  scale_colour_manual(name = "Cytogenetic subgroup",
                      values = c("High hyperdiploidy" = "dodgerblue2", 
                                 "Low hypodiploidy" = "red", 
                                 "Near triploidy" = "green3",
                                 "HoTr by SNP"="yellow"),
                      breaks = lev)



#Heatmap based on ploidy status----------------------------------------------

Heatmap(as.matrix(d_HoTr_HeH_norm[,5:26]))
annot_df <- data.frame(Cytogenetics = d_HoTr_HeH_norm$Cytogenetics)
col <- list(Cytogenetics = c("High hyperdiploidy" = "dodgerblue2", 
                             "Low hypodiploidy" = "red", 
                             "Near triploidy" = "green3",
                             "HoTr by SNP"="yellow"))

lab_level <- c("Low hypodiploidy","Near triploidy","HoTr By SNP","High hyperdiploidy")

row_ha = rowAnnotation(Cytogenetics = d_HoTr_HeH_norm$Cytogenetics, col = col, gp = gpar(fontsize = 10),
                       annotation_legend_param = list(
                         title = "Cytogenetics",
                         at = c("Low hypodiploidy","Near triploidy","HoTr by SNP","High hyperdiploidy"),
                         labels = c("Low hypodiploidy","Near triploidy","HoTr by SNP","High hyperdiploidy")))

col_ha = HeatmapAnnotation("     Info \n     gain" = anno_barplot(Info_gain$attr_importance))


p2 <- Heatmap(as.matrix(d_HoTr_HeH_norm[,5:26]) ,
              cluster_columns = FALSE, clustering_method_rows = "ward.D2",
              right_annotation = row_ha, 
              bottom_annotation = col_ha,
              row_names_gp = gpar(fontsize = 6),
              name="log2 ratios\n(normalised) ")


### Decision tree using all cases and k-fold cross validation-----------------------------------------------------------------------

d_norm <- normalize(as.matrix(d[,c(5:26)]), method = "standardize")
d_norm <- cbind(d[c(1:4)], d_norm)
colnames(d_norm) <- colnames(d)


classifier.fit <- rpart(Ploidy.genetic.subgroup ~ chr1+chr2+chr3+chr4+chr5+chr6+chr7+chr8+chr9+chr10+chr11+chr12+chr13+chr14+chr15+
                          chr16+chr17+chr18+chr19+chr20+chr21+chr22, data = d_norm, method = "class",
                        control=rpart.control(maxdepth = 2))
prp(classifier.fit, type=5, extra=101, nn=TRUE, fallen.leaves=TRUE, faclen=0, varlen=0, shadow.col="grey", branch.lty = 3,
    box.palette = list("dodgerblue","firebrick1","springgreen3"), legend.x = NA)


# Define training control
set.seed(3)
train.control <- trainControl(method = "cv", number = 10,
                              savePredictions = "final")
# Train the model
model <- train(Ploidy.genetic.subgroup ~ chr1+chr2+chr3+chr4+chr5+chr6+chr7+chr8+chr9+chr10+chr11+chr12+chr13+chr14+chr15+
                 chr16+chr17+chr18+chr19+chr20+chr21+chr22, data = d_norm, 
               method = "rpart2", tuneLength = 2,
               trControl = train.control)
# Summarize the results
print(model)
summary(model)
confusionMatrix(model)
confusionMatrix(model$pred$pred, model$pred$obs)
confusionMatrix(model$pred$pred, model$pred$obs, mode = "prec_recall")

