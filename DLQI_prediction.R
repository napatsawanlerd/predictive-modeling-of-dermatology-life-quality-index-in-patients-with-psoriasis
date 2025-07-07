#---------------------------------------------------------Modeling for Prediction ----------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------------#
library(rpart)
library(ranger)
library(kernlab)
library(gbm)
library(e1071)
library(caret)
library(xgboost)
library(boot)
library(creditmodel)
library(ggplot2)
library(glmnet)
library(Matrix)
library(readxl)
library(tidyr)
library(GGally) 
library(caret)
library(dplyr)
library(Metrics) 
library(stats)
library(randomForest)
library(gridExtra)
library(grid)
library(psych)
library(RColorBrewer)
library(reactable)
library(DT)
library(corrplot)
library(wesanderson)
library(systemfonts)
library(showtext)
library(iml)
library(grid)
library(viridis)
library(patchwork)
library(gridGraphics)
library(pdftools)
library(png)

format_4 <- function(x) { 
  sprintf("%.4f", as.numeric(x)) 
}


font_add(family = "Calisto", regular = "calisto-mt.ttf")
showtext_auto()
#---------------------------------------------Function Graph Predict data / Best Lambda--------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------------#
predict_actual <- function (test_data, pred_data, model) {
  test_data <- as.numeric(test_data)
  pred_data <- as.numeric(pred_data)
  plot(
    test_data,
    type = "l",
    col = "blue",
    lwd = 2,
    xlab = "Index",
    ylab = "DLQI",
    main = paste("Prediction vs Actual for", model)
  )
  lines(
    pred_data,
    col = "red",
    lwd = 2,
    lty = 2
  )
  legend(
    "topright",
    legend = c("Actual", "Predicted"),
    col = c("blue", "red"),
    lty = 1,
    lwd = 2
  )
}

Best_lambda <- function(results, best_lambda, model, x, y) {
  ggplot(results, aes(x = log10(lambda), y = RMSE)) +
    geom_line(color = "#7E38B7", size = 0.8) +
    geom_point(color = "#7E38B7", size = 1.5) +
    geom_vline(xintercept = log10(best_lambda), linetype = "dashed", color = "grey20", size = 0.8) +
    labs(
      #title = bquote("Validation RMSE vs. Log(" ~ lambda ~ ") for L" ~ .(model)),
      x = "Log(Lambda)",
      y = "Validation RMSE"
    ) +
    theme_minimal() +
    theme(
      plot.margin = unit(c(4, 4, 2, 2), "mm"),
      axis.text = element_text(size = 120, family = "Time"),
      axis.title = element_text(size = 120, family = "Time"),
      plot.title = element_text(size = 120, hjust = 0.5, family = "Time"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid = element_blank()
    ) +
    annotate(
      "label", 
      x = log10(best_lambda) - x, 
      y = min(results$RMSE) + y, 
      label = substitute(paste(lambda[x] == v), 
                         list(x = model, 
                              v = format_4(best_lambda))),
      parse = TRUE,
      fill = "white", 
      color = "black", 
      size = 40,  
      label.size = 0.5, 
      label.padding = unit(0.2, "cm"),
      label.r = unit(0.1, "cm"),
      family = "Time",  
      lineheight = 0.1
    )
}

#---------------------------------------------Function Calculate Loss Function--------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------------#
RMSE <- function(y_test, model_pred) {
  rmse <- rmse(y_test, model_pred)
  return(rmse)
}

RMSE_custom <- function(y_test, model_pred, X_model) {
  n <- length(y_test)
  p <- length(X_model)
  sqrt(sum((y_test - model_pred)^2) / (n - p))
}

R_square <- function(y_test, model_pred) {
  (cor(y_test, model_pred)^2)*100
}

adjusted_r_squared <- function(r_squared, y_test, X_model) {
  r_squared <- r_squared/100
  n <- length(y_test)
  k <- length(X_model)
  adj_r2 <- (1 - ((1 - r_squared) * (n - 1) / (n - k - 1)))*100
}

MAPE <- function(y_true, y_pred) {
  n <- length(y_test)
  median_value <- median(y_true[y_true != 0], na.rm = TRUE)
  y_true_adj <- ifelse(y_true == 0, median_value, y_true)   
  
  mape_value <- mean(abs((y_true_adj - y_pred) / y_true_adj), na.rm = TRUE) * 100
  return(mape_value)
}

#----------------------------------- Preparing Data for Prediction Scale (MinMax Scaler) + LOOCV -------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------------#
data <- read_excel("~/Documents/PASI_pj3:2/data/PASIData.xls", sheet = "SumStress_Group")
str(data)  #149 x 26

# Scaling numeric and categorical data
set.seed(845)
names(data)[names(data) == "PASI_W0"] <- "PASI"
numeric_var <- data[, c(2,4,7,26)]
categorical_var <- data[, c(8,3,9,10)]
ordinal_var <- data[, -c(1:10,25,26,27)]
num_sum <- cbind(numeric_var, ordinal_var)
str(num_sum)
var(num_sum)

describe <- describe(num_sum)
describe_table <- data.frame(
  Var = rownames(describe),
  Min = round(describe$min, 4),
  Mean = round(describe$mean, 4),
  Median = round(describe$median, 4),
  Max = round(describe$max, 4),
  SD = round(describe$sd, 4),
  SE = round(describe$se, 4)
)
datatable(describe_table, 
          options = list(pageLength = 5, autoWidth = TRUE), 
          caption = "Table of Describe")

# the result_table
tables <- list()
for (var in names(categorical_var)) {
  tables[[var]] <- table(categorical_var[[var]])
}

result_table <- data.frame(
  Variable = character(0),
  Level = character(0),
  Frequency = integer(0),
  Percentage = numeric(0)
)

for (var_name in names(categorical_var)) {
  table_data <- tables[[var_name]]
  variable_name <- var_name
  total_count <- sum(table_data)
  
  for (level_name in names(table_data)) {
    frequency <- table_data[[level_name]]
    percentage <- (frequency / total_count) * 100
    
    result_table <- rbind(result_table, data.frame(Variable = variable_name, Level = level_name, Frequency = frequency, Percentage = round(percentage, 4)))
  }
}
print(result_table)

# Spearman Correlation
cor_mat <- cor(num_sum, method = "spearman")
corrplot(cor_mat, method = "color", type = "lower",
         tl.col = "black", addCoef.col = "black", tl.pos = "lt",
         tl.cex = 1.3, addCoef.cex = 1.2,
         addCoef.font = 2,
         mar = c(0.2, 0.2, 0.2, 0.2))
cor_heat_plot(cor_mat)

y <- num_sum$DLQI
X_num <- num_sum[, !colnames(num_sum) %in% "DLQI"]
str(X_num)

# Dummy categorical variables
X_cat <- data.frame(lapply(categorical_var, as.factor))
X_cat$Maritial_Status <- relevel(X_cat$Maritial_Status, ref = "1") 
X_trans_cat <- model.matrix(~ . - 1, data = X_cat)
X_trans_cat <- X_trans_cat[, !grepl("Maritial_Status1", colnames(X_trans_cat))]
str(X_trans_cat)

# Scale numeric using Min-Max scaling
preProcessSteps <- preProcess(X_num, method = "range")
X_trans_num <- predict(preProcessSteps, newdata = X_num)

# Combine the Data
complete_data <- data.frame(X_trans_num, X_trans_cat, DLQI = y)
dim(complete_data)
head(complete_data)

# Split Train/Test
set.seed(845)
y <- complete_data$DLQI
X <- complete_data[, !colnames(complete_data) %in% "DLQI"]
train_index <- createDataPartition(y, p = 0.80, list = FALSE)
X_train <- X[train_index, , drop = FALSE]
y_train <- y[train_index]
X_test <- X[-train_index, , drop = FALSE]
y_test <- y[-train_index]

# DLQI <3 replace to median
#median_DLQI <- median(num_sum$DLQI[num_sum$DLQI >= 3], na.rm = TRUE)
#length(num_sum$DLQI[num_sum$DLQI < 3])
#y_train[y_train < 3] <- median_DLQI
#y_test[y_test < 3] <- median_DLQI
train <- cbind(X_train, DLQI = y_train)
# Train control with LOOCV
train_control <- trainControl(method = "LOOCV", returnResamp = "final", savePredictions = TRUE)
train_control
lambda_seq <- 10^seq(3, -3, by = -0.1)
length(y_test)
length(y_train)

#------------------------------Linear Regression (MinMax Scaler numeric & Spilt 80:20 & LOOCV)------------------------------------# 
set.seed(845)
train <- data.frame(X_train, DLQI = y_train)
LR <- train(
  DLQI ~ .,            
  data = train,      
  method = "lm",       
  trControl = train_control 
)
summary(LR$finalModel)
coef_LR <- coef(LR$finalModel)
LR_pred <- predict(LR, newdata = X_test)

# Evaluate Model
predict_actual(y_test, LR_pred, "Linear Regression")
RMSE_LR <- RMSE(y_test , LR_pred)
RMSE_custom_LR <- RMSE_custom(y_test, LR_pred, X_train)
R_squared_LR <- R_square(y_test, LR_pred)
adjR_squared_LR <- adjusted_r_squared(R_squared_LR, y_test, X_train)
MAPE_LR <- MAPE(y_test, LR_pred)

print(paste("RMSE:", RMSE_LR))
print(paste("Modified MAPE (%):", MAPE_LR))
print(paste("R-square (%):", R_squared_LR))
print(paste("Adjusted  R-square (%):", adjR_squared_LR))


#------------------------------Ridge Regression (MinMax Scaler numeric & Spilt 80:20 & LOOCV)------------------------------------# 
# Perform LOOCV to find the best lambda
set.seed(845)
Ridge_model <- train(
  DLQI ~ ., 
  data = train,
  method = "glmnet", 
  trControl = train_control, 
  tuneGrid = expand.grid(alpha = 0, lambda = lambda_seq)
)
best_lambda_Ridge <- Ridge_model$bestTune$lambda 
best_lambda_Ridge 

# Complete Ridge_model using best_lambda 
Best.Ridge_model <- train(
  DLQI ~ .,
  data = train,
  method = "glmnet", 
  trControl = train_control, 
  tuneGrid = expand.grid(alpha = 0, lambda = best_lambda_Ridge),  
)
Best.Ridge_model

coef_Ridge <- coef(Best.Ridge_model$finalModel, s = best_lambda_Ridge)
coef_Ridge
coef_Ridge_matrix <- as.matrix(coef_Ridge)
X_Ridge_features <- rownames(coef_Ridge)
X_Ridge_features <- X_Ridge_features[X_Ridge_features != "(Intercept)"]
Ridge_select <- X[, X_Ridge_features, drop = FALSE]

# Prediction Ridge_model
Ridge_pred <- predict(Best.Ridge_model, newdata = X_test)

# Evaluate Model
predict_actual(y_test, Ridge_pred, "L2")
RMSE_Ridge <- RMSE(y_test , Ridge_pred)
RMSE_custom_Ridge <- RMSE_custom(y_test, Ridge_pred, X_Ridge_features)
R_squared_Ridge <- R_square(y_test, Ridge_pred)
adjR_squared_Ridge <- adjusted_r_squared(R_squared_Ridge, y_test, X_Ridge_features)

MAPE_Ridge <- MAPE(y_test, Ridge_pred)
print(paste("RMSE:", RMSE_Ridge))
print(paste("Modified MAPE (%):", MAPE_Ridge))
print(paste("R-square (%):", R_squared_Ridge))
print(paste("Adjusted  R-square (%):", adjR_squared_Ridge))


#------------------------------LASSO Regression (MinMax Scaler numeric & Spilt 80:20 & LOOCV) ------------------------------------#
# Perform LOOCV find the best lambda
set.seed(845)
LASSO_model <- train(
  DLQI ~ ., 
  data = train,
  method = "glmnet", 
  trControl = train_control,
  tuneGrid = expand.grid(alpha = 1, lambda = lambda_seq)
)
best_lambda_LASSO <- LASSO_model$bestTune$lambda
best_lambda_LASSO

# Complete Ridge_model using best_lambda 
Best.LASSO_model <- train(
  DLQI ~ ., 
  data = train,
  method = "glmnet", 
  trControl = train_control, 
  tuneGrid = expand.grid(alpha = 1, lambda = best_lambda_LASSO),  
)
Best.LASSO_model

png(filename = "~/Documents/SHAP/lambda_L1_coef.png", 
    width = 4, height = 3, units = "in", res = 1000)
par(mar = c(4, 4, 2, 2), family="Time")
a2 <- plot(LASSO_model$finalModel, xvar = "lambda", label = TRUE,
     ylab = "Coefficients",
     xlab = expression(Log(lambda)), 
     cex.axis = 1.1, 
     cex.lab = 1.1, 
     cex.main = 1.1, 
     lwd = 1.1,
     xlim = c(-6, 2))
box(lwd = 1.2)
dev.off()

p2 <- Best_lambda(LASSO_model$results, best_lambda_LASSO, "1", 2, 0.001) +
  scale_x_continuous(limits = c(-3, 3)) +
  scale_y_continuous(limits = c(7.2, 8)) +
  labs(x = expression(Log(lambda))) +
  theme(text = element_text(family = "Time")) 
graph <- "~/Documents/SHAP"
ggsave(
  filename = file.path(graph, "lambda_L1.png"), 
  plot = p2,                 
  width = 4,                                      
  height = 3,                                     
  dpi = 1000                                     
)

coef_LASSO <- coef(Best.LASSO_model$finalModel, s = best_lambda_LASSO)
coef_LASSO

coef_LASSO_matrix <- as.matrix(coef_LASSO)
X_LASSO_model <- rownames(coef_LASSO_matrix)[coef_LASSO_matrix != 0]
length(X_LASSO_model)
X_LASSO_features <- X_LASSO_model
X_LASSO_features <- X_LASSO_features[X_LASSO_features != "(Intercept)"]
LASSO_select <- X[, X_LASSO_features, drop = FALSE]

# Prediction LASSO model
LASSO_pred <- predict(Best.LASSO_model, newdata = X_test)

# Evaluate Model
predict_actual(y_test, LASSO_pred, "L1")
RMSE_LASSO <- RMSE(y_test, LASSO_pred)
RMSE_custom_LASSO <- RMSE_custom(y_test, LASSO_pred, X_LASSO_features)
R_squared_LASSO <- R_square(y_test, LASSO_pred)
MAPE_LASSO <- MAPE(y_test, LASSO_pred)
adjR_squared_LASSO <- adjusted_r_squared(R_squared_LASSO, y_test, X_LASSO_features)

print(paste("RMSE:", RMSE_LASSO))
print(paste("Modified MAPE (%):", MAPE_LASSO))
print(paste("R-square (%):", R_squared_LASSO))
print(paste("Adjusted  R-square (%):", adjR_squared_LASSO))

#------------------------------Adaptive LASSO Regression (MinMax Scaler numeric & Spilt 80:20 & LOOCV) ------------------------------------#
set.seed(845)
weights <- 1 / abs(coef_Ridge[-1])
weights

# Complete AL_model using best_lambda
Best.AL_model <- train(
  DLQI ~ ., 
  data = train,
  method = "glmnet", 
  trControl = train_control, 
  tuneGrid = expand.grid(alpha = 1, lambda = best_lambda_LASSO),
  penalty.factor = weights
)
Best.AL_model
coef_AL <- coef(Best.AL_model$finalModel, s = best_lambda_LASSO)
coef_AL

coef_AL_matrix <- as.matrix(coef_AL)
coef_AL_matrix
X_AL_model <- rownames(coef_AL_matrix)[coef_AL_matrix != 0]
length(X_AL_model)
X_AL_features <- X_AL_model
X_AL_features <- X_AL_features[X_AL_features != "(Intercept)"]
AL_select <- X[, X_AL_features, drop = FALSE]

# Prediction AL_model
AL_pred <- predict(Best.AL_model, newdata = X_test)

# Evaluate Model
predict_actual(y_test, AL_pred, "A-L1")
RMSE_AL <- RMSE(y_test, AL_pred)
RMSE_custom_AL <- RMSE_custom(y_test, AL_pred, X_AL_features)
R_squared_AL <- R_square(y_test, AL_pred)
MAPE_AL <- MAPE(y_test, AL_pred)
adjR_squared_AL <- adjusted_r_squared(R_squared_AL, y_test, X_AL_features)

print(paste("RMSE:", RMSE_AL))
print(paste("Modified MAPE (%):", MAPE_AL))
print(paste("R-square (%):", R_squared_AL))
print(paste("Adjusted  R-square (%):", adjR_squared_AL))


#------------------------------Elastic net Regression (MinMax Scaler numeric & Spilt 80:20 & LOOCV) ------------------------------------#
# Perform LOOCV find the best lambda
set.seed(845)
EN_model <- train(
  DLQI ~ ., 
  data = train,
  method = "glmnet", 
  trControl = train_control, 
  tuneGrid = expand.grid(alpha = 0.5, lambda = lambda_seq)
)
best_lambda_EN <- EN_model$bestTune$lambda
best_lambda_EN

# Complete EN_model using best_lambda
Best.EN_model <- train(
  DLQI ~ ., 
  data = train,
  method = "glmnet", 
  trControl = train_control, 
  tuneGrid = expand.grid(alpha = 0.5, lambda = best_lambda_EN),
)
Best.EN_model
png(filename = "~/Documents/SHAP/lambda_L1L2_coef.png", 
    width = 4, height = 3, units = "in", res = 1000)
par(mar = c(4, 4, 2, 2), family="Time")
a3 <- plot(EN_model$finalModel, xvar = "lambda", label = TRUE,
     ylab = "Coefficients",
     xlab = expression(Log(lambda)), 
     cex.axis = 1.1, 
     cex.lab = 1.1, 
     cex.main = 1.1, 
     lwd = 1.1,
     xlim = c(-6, 2))
box(lwd = 1.2)
dev.off()

p3 <- Best_lambda(EN_model$results, best_lambda_EN, "12", 2, 0.001) +
  scale_x_continuous(limits = c(-3, 3)) +
  scale_y_continuous(limits = c(7.2, 8)) +
  labs(x = expression(Log(lambda))) +
  theme(text = element_text(family = "Time")) 
graph <- "~/Documents/SHAP"
ggsave(
  filename = file.path(graph, "lambda_L1L2.png"), 
  plot = p3,                 
  width = 4,                                      
  height = 3,                                     
  dpi = 1000                                     
)

coef_EN <- coef(Best.EN_model$finalModel, s = best_lambda_EN)
coef_EN
coef_EN_matrix <- as.matrix(coef_EN)
coef_EN_matrix
X_EN_model <- rownames(coef_EN_matrix)[coef_EN_matrix != 0]
length(X_EN_model)
X_EN_features <- X_EN_model
X_EN_features <- X_EN_features[X_EN_features != "(Intercept)"]
EN_select <- X[, X_EN_features, drop = FALSE]

# Prediction EN_model
EN_pred <- predict(Best.EN_model, newdata = X_test)

# Evaluate Model
predict_actual(y_test, EN_pred, "L1L2")
RMSE_EN <- RMSE(y_test, EN_pred)
RMSE_custom_EN <- RMSE_custom(y_test, EN_pred, X_EN_features)
R_squared_EN <- R_square(y_test, EN_pred)
MAPE_EN <- MAPE(y_test, EN_pred)
adjR_squared_EN <- adjusted_r_squared(R_squared_EN, y_test, X_EN_features)

print(paste("RMSE:", RMSE_EN))
print(paste("Modified MAPE (%):", MAPE_EN))
print(paste("R-square (%):", R_squared_EN))
print(paste("Adjusted  R-square (%):", adjR_squared_EN))


#------------------------------ Adaptive Elastic net Regression (MinMax Scaler numeric & Spilt 80:20 & LOOCV) -----------------------------------#
# Complete AEN_model using best_lambda
set.seed(845)
Best.AEN_model <- train(
  DLQI ~ ., 
  data = train,
  method = "glmnet", 
  trControl = train_control, 
  tuneGrid = expand.grid(alpha = 0.5, lambda = best_lambda_EN),
  penalty.factor = weights
)
Best.AEN_model
coef_AEN <- coef(Best.AEN_model$finalModel, s = best_lambda_EN)
coef_AEN

coef_AEN_matrix <- as.matrix(coef_AEN)
coef_AEN_matrix
X_AEN_model <- rownames(coef_AEN_matrix)[coef_AEN_matrix != 0]
length(X_AEN_model)
X_AEN_features <- X_AEN_model
X_AEN_features <- X_AEN_features[X_AEN_features != "(Intercept)"]
AEN_select <- X[, X_AEN_features, drop = FALSE]

# Prediction AEN_model
AEN_pred <- predict(Best.AEN_model, newdata = X_test)

# Evaluate Model
predict_actual(y_test, AEN_pred, "A-L1L2")
RMSE_AEN <- RMSE(y_test, AEN_pred)
RMSE_custom_AEN <- RMSE_custom(y_test, AEN_pred, X_AEN_features)
R_squared_AEN <- R_square(y_test, AEN_pred)
MAPE_AEN <- MAPE(y_test, AEN_pred)
adjR_squared_AEN <- adjusted_r_squared(R_squared_AEN, y_test, X_AEN_features)

print(paste("RMSE:", RMSE_EN))
print(paste("Modified MAPE (%):", MAPE_AEN))
print(paste("R-square (%):", R_squared_AEN))
print(paste("Adjusted  R-square (%):", adjR_squared_AEN))

#----------------------------------- Summerize Selection Data (MinMax Scaler numeric & Spilt 80:20 & LOOCV)  ----------------------------------------#
# Remove Intercept
dim(Ridge_select)
dim(LASSO_select)
dim(AL_select)
dim(EN_select)
dim(AEN_select)

Best_lambda <- function(results, best_lambda, model, x, y) {
  ggplot(results, aes(x = log10(lambda), y = RMSE)) +
    geom_line(color = "#7E38B7", size = 1.3) +
    geom_point(color = "#7E38B7", size = 1.9) +
    geom_vline(xintercept = log10(best_lambda), linetype = "dashed", color = "grey20", size = 0.8) +
    labs(
      #title = bquote("Validation RMSE vs. Log(" ~ lambda ~ ") for L" ~ .(model)),
      x = "Log(Lambda)",
      y = "Validation RMSE"
    ) +
    theme_minimal() +
    theme(
      plot.margin = unit(c(3, 3, 1, 1), "mm"),
      axis.text = element_text(color = "black", size = 23, family = "Time"),
      axis.title = element_text(color = "black", size = 23, family = "Time"),
      plot.title = element_text(color = "black", size = 23, hjust = 0.5, family = "Time"),
      panel.border = element_rect(color = "black", fill = NA, size = 1.3),
      panel.grid = element_blank()
    ) +
    annotate(
      "label", 
      x = log10(best_lambda) - x, 
      y = min(results$RMSE) + y, 
      label = substitute(paste(lambda[x] == v),
                         list(x = model,
                              v = format_4(best_lambda))),
      parse = TRUE,
      fill = "white", 
      color = "black", 
      size = 8,  
      label.size = 0.5, 
      label.padding = unit(0.2, "cm"),
      label.r = unit(0.1, "cm"),
      family = "Time",  
      lineheight = 0.1
    )
}

# Ridge Model Plot
png_file1 <- tempfile(fileext = ".png")
png(png_file1, width = 3600, height = 3000, res = 600)
plot(Ridge_model$finalModel, xvar = "lambda", label = TRUE,
     ylab = "Coefficients",
     xlab = expression(Log(lambda)),
     cex.axis = 1.7,
     cex.lab = 1.7,
     cex.main = 1.7,
     lwd = 1.7,
     xlim = c(-2, 8))
box(lwd = 1.5)
dev.off()
img1 <- png::readPNG(png_file1)

a1 <- grid::rasterGrob(img1, width = unit(1, "npc"), height = unit(1.5, "npc"))
a1_p <- wrap_elements(a1, clip = FALSE) + 
  theme(plot.margin = unit(c(0, 10, 0, 10), "mm"), 
        panel.background = element_blank(),
        plot.background = element_blank())

# LASSO Model Plot
png_file2 <- tempfile(fileext = ".png")
png(png_file2, width = 3600, height = 3000, res = 600)
plot(LASSO_model$finalModel, xvar = "lambda", label = TRUE,
     ylab = "Coefficients",
     xlab = expression(Log(lambda)),
     cex.axis = 1.7,
     cex.lab = 1.7,
     cex.main = 1.7,
     lwd = 1.7,
     xlim = c(-6, 2))
box(lwd = 1.5)
dev.off()
img2 <- png::readPNG(png_file2)

a2 <- grid::rasterGrob(img2, width = unit(1, "npc"), height = unit(1.5, "npc"))
a2_p <- wrap_elements(a2, clip = FALSE) + 
  theme(plot.margin = unit(c(0, 10, 0, 10), "mm"), 
        panel.background = element_blank(),
        plot.background = element_blank())

# Elastic Net Model Plot
png_file3 <- tempfile(fileext = ".png")
png(png_file3, width = 3600, height = 3000, res = 600)
plot(EN_model$finalModel, xvar = "lambda", label = TRUE,
     ylab = "Coefficients",
     xlab = expression(Log(lambda)),
     cex.axis = 1.7,
     cex.lab = 1.7,
     cex.main = 1.7,
     lwd = 1.7,
     xlim = c(-6, 2))
box(lwd = 1.5)
dev.off()
img3 <- png::readPNG(png_file3)

a3 <- grid::rasterGrob(img3, width = unit(1, "npc"), height = unit(1.5, "npc"))
a3_p <- wrap_elements(a3, clip = FALSE) + 
  theme(plot.margin = unit(c(0, 10, 0, 10), "mm"), 
        panel.background = element_blank(),
        plot.background = element_blank())

# Create the Best_lambda plots
p1 <- Best_lambda(Ridge_model$results, best_lambda_Ridge, "2", 3, 0.001) +
  scale_x_continuous(limits = c(-3, 3)) +
  scale_y_continuous(limits = c(7.2, 8)) +
  labs(x = expression(Log(lambda))) +
  theme(text = element_text(family = "Time"),
        plot.margin = unit(c(10, 20, 10, 10), "mm"))

p2 <- Best_lambda(LASSO_model$results, best_lambda_LASSO, "1", 2, 0.001) +
  scale_x_continuous(limits = c(-3, 3)) +
  scale_y_continuous(limits = c(7.2, 8)) +
  labs(x = expression(Log(lambda))) +
  theme(text = element_text(family = "Time"),
        plot.margin = unit(c(10, 20, 10, 10), "mm"))  

p3 <- Best_lambda(EN_model$results, best_lambda_EN, "12", 2, 0.001) +
  scale_x_continuous(limits = c(-3, 3)) +
  scale_y_continuous(limits = c(7.2, 8)) +
  labs(x = expression(Log(lambda))) +
  theme(text = element_text(family = "Time"),
        plot.margin = unit(c(10, 20, 10, 10), "mm"))  

# Combine all plots 
combined_all_plots <- (a1_p | p1) / (a2_p | p2) / (a3_p | p3) +
  plot_layout(heights = c(1, 1, 1))

graph <- "~/Documents/SHAP"
ggsave(
  filename = file.path(graph, "best_lambda2.pdf"),
  plot = combined_all_plots,
  width = 16,
  height = 16, 
  device = "pdf"
)


#------------ Evaluate Table
resultstable_penal <- data.frame(
  Model = c("LR", "L2", "L1", "A-L1","L1L2", "A-L1L2"),
  num_var = c(length(X_train), dim(Ridge_select)[2], dim(LASSO_select)[2], dim(AL_select)[2], dim(EN_select)[2], dim(AEN_select)[2]),
  Lambda = c("-", format_4(best_lambda_Ridge), format_4(best_lambda_LASSO), format_4(best_lambda_LASSO), format_4(best_lambda_EN), format_4(best_lambda_EN)),
  RMSE = c(format_4(RMSE_LR), format_4(RMSE_Ridge), format_4(RMSE_LASSO), format_4(RMSE_AL), format_4(RMSE_EN), format_4(RMSE_AEN)),
  MAPE = c(format_4(MAPE_LR), format_4(MAPE_Ridge), format_4(MAPE_LASSO), format_4(MAPE_AL), format_4(MAPE_EN), format_4(MAPE_AEN)),
  RMSE_custom = c(format_4(RMSE_custom_LR), format_4(RMSE_custom_Ridge), format_4(RMSE_custom_LASSO), format_4(RMSE_custom_AL), format_4(RMSE_custom_EN), format_4(RMSE_custom_AEN)),
  R_squared = c(format_4(R_squared_LR), format_4(R_squared_Ridge), format_4(R_squared_LASSO), format_4(R_squared_AL), format_4(R_squared_EN), format_4(R_squared_AEN)),
  Adj_R_squared = c(format_4(adjR_squared_LR), format_4(adjR_squared_Ridge), format_4(adjR_squared_LASSO), format_4(adjR_squared_AL), format_4(adjR_squared_EN), format_4(adjR_squared_AEN))
)
resultstable_penal
datatable(resultstable_penal, 
          options = list(pageLength = 5, autoWidth = TRUE), 
          caption = "Table of Results (Penalized Regression)")

#------------ Coef Table
s0 <- format_4(as.matrix(coef_LR))
s1 <- format_4(as.matrix(coef_Ridge))
s2 <- format_4(as.matrix(coef_LASSO))
s3 <- format_4(as.matrix(coef_AL))
s4 <- format_4(as.matrix(coef_EN))
s5 <- format_4(as.matrix(coef_AEN))
coef_LR_df <- as.data.frame(s0) 
coef_Ridge_df <- as.data.frame(s1)  
coef_LASSO_df <- as.data.frame(s2)  
coef_AL_df <- as.data.frame(s3)  
coef_EN_df <- as.data.frame(s4)
coef_AEN_df <- as.data.frame(s5)
coef_results <- data.frame(
  Variables = rownames(coef_Ridge),
  LR = coef_LR_df$s0,
  L2 = coef_Ridge_df$s1,
  L1 = coef_LASSO_df$s2,
  A_L1 = coef_AL_df$s3,
  L1L2 = coef_EN_df$s4,
  A_L1L2 = coef_AEN_df$s5
)
coef_results <- coef_results %>% 
  mutate(across(everything(), ~ ifelse(. == "0.0000", "-", .)))
coef_results

datatable(coef_results, 
          options = list(pageLength = 5, autoWidth = TRUE), 
          caption = "Table of Coefficient (Penalized Regression)")

#------------ Feature Selection Graph
model_names <- c("L2", "L1", "A-L1", "L1L2", "A-L1L2")
num_vars <- c(dim(Ridge_select)[2], dim(LASSO_select)[2], dim(AL_select)[2], dim(EN_select)[2], dim(AEN_select)[2])
# save graph
png(filename = "~/Documents/PASI_pj3:2/Penalized_Report/graph_penalized/num_features.png", 
    width = 12, height = 8, units = "in", res = 500)
par(mar = c(5, 6, 4, 2) + 0.1, family="Calisto")
bar_positions <- barplot(num_vars, 
                         names.arg = model_names, 
                         main = "Number of Features Selected (149 data)",   
                         xlab = "Model",                                   
                         ylab = "Number of Features Selected",            
                         col = "#436EEE",                                  
                         ylim = c(0, max(num_vars) + 5),                   
                         family = "Calisto",                             
                         cex.names = 2,                                   
                         cex.axis = 2,                                 
                         cex.lab = 2.3,                          
                         cex.main = 2.5)                                  
text(x = bar_positions, y = num_vars, label = num_vars, pos = 3, cex = 2.5, col = "black", font = 2)
dev.off()

#------------ Evaluate Graph
theme_custom <- theme_minimal() +
  theme(
    axis.title.x = element_text(size = 18, family = "Calisto"),
    axis.title.y = element_text(size = 18, family = "Calisto"),
    plot.title = element_text(size = 14, face = "bold", family = "Calisto"),
    axis.text.x = element_text(size = 18, family = "Calisto"),
    axis.text.y = element_text(size = 18, family = "Calisto"),
    legend.position = "none",
    panel.grid = element_blank()
  )
RMSE_penal_model = c(RMSE_LR, RMSE_Ridge,RMSE_LASSO,RMSE_AL,RMSE_EN,RMSE_AEN)
data_RMSE <- data.frame(
  Model = c("LR", "L2", "L1", "A-L1", "L1L2", "A-L1L2"),
  RMSE_model = c(RMSE_LR, RMSE_Ridge,RMSE_LASSO,RMSE_AL,RMSE_EN,RMSE_AEN)
)
min_rmse <- min(data_RMSE$RMSE_model)
data_RMSE$Model <- factor(data_RMSE$Model, levels = c("LR", "L2", "L1", "A-L1", "L1L2", "A-L1L2"))
p1 <- ggplot(data_RMSE, aes(x = Model, y = RMSE_model, fill = RMSE_model == min_rmse)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_text(aes(label = format_4(RMSE_model)), vjust = -0.5, family = "Calisto") +
  labs(title = "RMSE", x = "Model", y = "RMSE") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "#D3D3D3")) +
  geom_hline(yintercept = min_rmse, linetype = "dashed", color = "blue") +
  theme_custom

RMSE_custom_penal_model = c(RMSE_custom_LR, RMSE_custom_Ridge,RMSE_custom_LASSO,RMSE_custom_AL,RMSE_custom_EN,RMSE_custom_AEN)
data_RMSE_custom <- data.frame(
  Model = c("LR", "L2", "L1", "A-L1", "L1L2", "A-L1L2"),
  RMSE_custom_model = c(RMSE_custom_LR, RMSE_custom_Ridge,RMSE_custom_LASSO,RMSE_custom_AL,RMSE_custom_EN,RMSE_custom_AEN)
)
min_rmse_custom <- min(data_RMSE_custom$RMSE_custom_model)
data_RMSE_custom$Model <- factor(data_RMSE_custom$Model, levels = c("LR", "L2", "L1", "A-L1", "L1L2", "A-L1L2"))
p2 <- ggplot(data_RMSE_custom, aes(x = Model, y = RMSE_custom_model, fill = RMSE_custom_model == min_rmse_custom)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_text(aes(label = format_4(RMSE_custom_model)), vjust = -0.5, family = "Calisto") +
  labs(title = "RMSE", x = "Model", y = "RMSE") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "#D3D3D3")) +
  geom_hline(yintercept = min_rmse_custom, linetype = "dashed", color = "blue") +
  theme_custom

MAPE_penal_model = c(MAPE_LR, MAPE_Ridge, MAPE_LASSO, MAPE_AL, MAPE_EN, MAPE_AEN)
data_MAPE <- data.frame(
  Model = c("LR", "L2", "L1", "A-L1", "L1L2", "A-L1L2"),
  MAPE_model = c(MAPE_LR, MAPE_Ridge, MAPE_LASSO, MAPE_AL, MAPE_EN, MAPE_AEN)
)
min_mape <- min(data_MAPE$MAPE_model)
data_MAPE$Model <- factor(data_MAPE$Model, levels = c("LR", "L2", "L1", "A-L1", "L1L2", "A-L1L2"))
p3 <- ggplot(data_MAPE, aes(x = Model, y = MAPE_model, fill = MAPE_model == min_mape)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_text(aes(label = format_4(MAPE_model)), vjust = -0.5, family = "Calisto") +
  labs(title = "MAPE (%)", x = "Model", y = "MAPE (%)") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "#D3D3D3")) +
  geom_hline(yintercept = min_mape, linetype = "dashed", color = "blue") +
  theme_custom

R_squared_penal_model = c(R_squared_LR, R_squared_Ridge, R_squared_LASSO, R_squared_AL, R_squared_EN, R_squared_AEN)
data_R_squared <- data.frame(
  Model = c("LR", "L2", "L1", "A-L1", "L1L2", "A-L1L2"),
  R_squared_model = c(R_squared_LR, R_squared_Ridge, R_squared_LASSO, R_squared_AL, R_squared_EN, R_squared_AEN)
)
max_r_squared <- max(data_R_squared$R_squared_model)
data_R_squared$Model <- factor(data_R_squared$Model, levels = c("LR", "L2", "L1", "A-L1", "L1L2", "A-L1L2"))
p4 <- ggplot(data_R_squared, aes(x = Model, y = R_squared_model, fill = R_squared_model == max_r_squared)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_text(aes(label = format_4(R_squared_model)), vjust = -0.5, family = "Calisto") +
  labs(title = "R_squared (%)", x = "Model", y = "R_squared (%)") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "#D3D3D3")) +
  geom_hline(yintercept = max_r_squared, linetype = "dashed", color = "blue") +
  theme_custom
grid.arrange(
  p1, p3, p4,
  ncol = 2,
  top = textGrob("Penalized Regression", 
                 gp = gpar(fontsize = 15, fontface = "bold", family = "Calisto"))
)


custom_colors <- c(
  "Test Data" = "grey",       
  "L2" = "#0000CD",       
  "L1" = "#00FFFF",        
  "A-L1" = "#40E0D0",
  "L1L2" = "#FF66FF",  
  "A-L1L2" = "#8470FF" 
)
#------------------------------- Graph Compare Model
pred_penal = list(LR_pred, Ridge_pred, LASSO_pred, AL_pred, EN_pred, AEN_pred)
model_names = c("LR", "L2", "L1", "A-L1", "L1L2", "A-L1L2")
n <- length(y_test)
x <- 1:n
data <- data.frame(
  x = rep(x, length(pred_penal) + 1),
  y = c(y_test, unlist(pred_penal)),
  category = rep(c("Test Data", model_names), each = n)
)
theme_model <- theme(
  panel.border = element_rect(color = "black", fill = NA, size = 1), 
  panel.grid = element_blank(), 
  plot.title = element_text(size = 150, face = "bold", hjust = 0.5, family="Calisto"),
  axis.title.x = element_text(size = 150, family="Calisto"),
  axis.title.y = element_text(size = 150, family="Calisto"),
  axis.text.x = element_text(size = 150, family="Calisto"), 
  axis.text.y = element_text(size = 150, family="Calisto"), 
  legend.text = element_text(size = 90, family="Calisto"), 
  legend.position = c(0.98, 0.03), 
  legend.justification = c("right", "bottom"), 
  legend.background = element_rect(fill = "white", color = "black", size = 0.5),
  legend.box.background = element_rect(color = "black", size = 0.7),
  legend.key.width = unit(2, "cm"), 
  legend.key.height = unit(0.7, "cm") 
)
custom_colors <- c(
  "Test Data" = "grey",
  "LR" = "tomato",
  "L2" = "#61D04F",       
  "L1" =  "#CD0BBC",        
  "A-L1" = "#F5C710",
  "L1L2" = "#8470FF",  
  "A-L1L2" = "#2297E6"
)
data$category <- factor(
  data$category, 
  levels = c("Test Data", "LR", "L2", "L1", "A-L1", "L1L2", "A-L1L2")
)
ggplot(data, aes(x = x, y = y)) +
  geom_line(aes(color = category, linetype = category), size = 1.5) + 
  scale_color_manual(values = custom_colors) + 
  scale_linetype_manual(
    values = c(
      "Test Data" = "dashed", 
      "LR" = "solid",
      "L2" = "solid", 
      "L1" = "solid", 
      "A-L1" = "solid", 
      "L1L2" = "solid", 
      "A-L1L2" = "solid"
    )
  ) +
  labs(
    title = "Comparison of Model Predictions (Penalized Regression)",
    x = "Index",
    y = "Prediction",
    color = NULL,
    linetype = NULL
  ) +
  theme_minimal() +
  theme_model
graph <- "~/Documents/PASI_pj3:2/Penalized_Report/graph_penalized/Prediction_Graph"
ggsave(
  filename = file.path(graph, "pred_penal.png"), 
  plot = last_plot(),                    
  width = 12,                                     
  height = 8,                                     
  dpi = 500                                    
)


#-----------Loop for Actual vs Predict & Actual VS Error----------#
predict_actual_plot <- function(test_data, pred_data, model_name) {
  test_data <- as.numeric(unlist(test_data))
  pred_data <- as.numeric(unlist(pred_data))
  
  if (length(test_data) != length(pred_data)) stop("Vectors must have the same length.")
  median_test <- median(test_data[test_data != 0], na.rm = TRUE)
  test_data <- ifelse(test_data == 0, median_test, test_data) 
  print(paste("Median value:", median_test))
  
  analysis_data <- data.frame(
    Actual = test_data,
    Predicted = pred_data,
    Error = pred_data - test_data,
    Relative_Error = ((pred_data - test_data)/test_data)*100,
    Average_Actaul = (pred_data + test_data)/2
  )
  print(analysis_data)
  
  # Calculate Mean and ±1.96 SD for Bland-Altman
  pearson_corr <- cor(test_data, pred_data, method = "pearson", use = "complete.obs")
  mean_error <- mean(analysis_data$Relative_Error, na.rm = TRUE)
  print(paste("Mean Relative Error:", mean_error))
  sd_error <- sd(analysis_data$Relative_Error, na.rm = TRUE)
  print(paste("SD Relative Error:", sd_error))
  upper_limit <- mean_error + 1.96 * sd_error
  lower_limit <- mean_error - 1.96 * sd_error
  
  # homoscedasticity
  max_limit <- 25
  path_homoscedasticity <- "~/Documents/PASI_pj3:2/Penalized_Report/graph_penalized/Homoscedasticity/penal"
  p_homoscedasticity <- ggplot(analysis_data, aes(x = Actual, y = Predicted)) +
    geom_point(alpha = 0.7, size = 3, color = "#1F77B4") +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") + 
    coord_cartesian(xlim = c(0, max_limit), ylim = c(0, max_limit)) + 
    scale_x_continuous(breaks = seq(0, max_limit, by = 5)) + 
    scale_y_continuous(breaks = seq(0, max_limit, by = 5)) + 
    labs(
      title = paste("Actual vs Predicted (", model_name, ")", sep = ""),
      x = "Actual Y",
      y = "Predicted Y"
    ) +
    theme_minimal(base_family = "Calisto") +
    theme(
      plot.title = element_text(size = 110, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 110),
      axis.title.y = element_text(size = 110),
      axis.text = element_text(size = 110),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1.5)
    )
  
  # Bland-Altman Plot
  path_bland_altman <- "~/Documents/PASI_pj3:2/Penalized_Report/graph_penalized/Bland_altman/penal"
  p_bland_altman <- ggplot(analysis_data, aes(x = Average_Actaul, y = Relative_Error)) +
    geom_point(alpha = 0.7, size = 3, color = "#FF7F0E") +
    geom_hline(yintercept = mean_error, color = "blue", linetype = "solid", size = 1) +
    geom_hline(yintercept = upper_limit, color = "darkgreen", linetype = "dashed", size = 1) +
    geom_hline(yintercept = lower_limit, color = "darkgreen", linetype = "dashed", size = 1) +
    annotate("text", x = max(analysis_data$Average_Actaul), y = mean_error, 
             label = sprintf("Mean: %.2f", mean_error), color = "blue", hjust = 1, vjust = -1, size = 35, family = "Calisto") +
    annotate("text", x = max(analysis_data$Average_Actaul), y = upper_limit, 
             label = sprintf("+1.96SD: %.2f", upper_limit), color = "darkgreen", hjust = 1, vjust = -1, size = 35, family = "Calisto") +
    annotate("text", max(analysis_data$Average_Actaul), y = lower_limit, 
             label = sprintf("-1.96SD: %.2f", lower_limit), color = "darkgreen", hjust = 1, vjust = -1, size = 35, family = "Calisto") +
    labs(
      title = paste("Actual vs Error (", model_name, ")", sep = ""),
      x = "Mean of Actual and Predicted",
      y = "Relative Error %"
    ) +
    theme_minimal(base_family = "Calisto") +
    theme(
      plot.title = element_text(size = 110, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 110),
      axis.title.y = element_text(size = 110),
      axis.text = element_text(size = 110),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1.5)
    )
  
  ggsave(
    filename = file.path(path_homoscedasticity, paste("homoscedasticity_", model_name, ".png", sep = "")),
    plot = p_homoscedasticity,
    width = 5,
    height = 5,
    dpi = 500
  )
  ggsave(
    filename = file.path(path_bland_altman, paste("bland_altman_", model_name, ".png", sep = "")),
    plot = p_bland_altman,
    width = 7,
    height = 5,
    dpi = 500
  )
  
  return(list(homoscedasticity_plot = p_homoscedasticity, 
              bland_altman_plot = p_bland_altman,
              mean_error = format_4(mean_error),
              sd_error = format_4(sd_error),
              upper_limit = format_4(upper_limit),
              lower_limit = format_4(lower_limit),
              pearson_corr = format_4(pearson_corr)))
}

models <- c("LR", "L2", "L1", "A-L1", "L1L2", "A-L1L2")
predictions <- list(LR_pred, Ridge_pred, LASSO_pred, AL_pred, EN_pred, AEN_pred)

results_penal <- list()
for (i in seq_along(models)) {
  result <- predict_actual_plot(y_test, predictions[[i]], models[i])
  
  results_penal[[i]] <- data.frame(
    Model = models[i],
    Relative_Error = result$mean_error,
    SD_Error = result$sd_error,
    Upper_Limit = result$upper_limit,
    Lower_Limit = result$lower_limit,
    Pearson_Correlation = result$pearson_corr
  )
}

table_Penalerror <- do.call(rbind, results_penal)
print(table_Penalerror)
datatable(table_Penalerror, 
          options = list(pageLength = 5, autoWidth = TRUE), 
          caption = "Table of Error (Penalized Regression)")


#------------------------ Random forest + Penalized 5 Models (Std.numeric & Spilt 80:20 & LOOCV model) ---------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
# Grid Search
Random_forest <- function(X_select, X_features, model, n_tree) {
  set.seed(845)
  select_data <- data.frame(X_select, DLQI = y)
  
  train_index <- createDataPartition(select_data$DLQI, p = 0.8, list = FALSE)
  train_data <- select_data[train_index, ]
  test_data <- select_data[-train_index, ]
  
  num_trees = n_tree
  rf_grid <- expand.grid(
    mtry = seq(2, min(ncol(train_data) - 1, 15), by = 1),   
    splitrule = c("extratrees"),      
    min.node.size = c(5, 8, 10, 12, 14, 16, 18, 20)
  )
  
  train_control <- trainControl(
    method = "LOOCV",
    search = "grid"
  )
  
  rf_model <- train(
    DLQI ~ .,
    data = train_data,
    method = "ranger",
    trControl = train_control,
    tuneGrid = rf_grid,
    num.trees = num_trees  
  )
  
  print("Best parameters found in grid search:")
  print(rf_model$bestTune)
  
  final_model <- ranger(
    formula = DLQI ~ .,
    data = train_data,
    mtry = rf_model$bestTune$mtry,
    splitrule = rf_model$bestTune$splitrule,
    min.node.size = rf_model$bestTune$min.node.size,
    num.trees = num_trees
  )
  
  rf_pred_test <- predict(final_model, data = test_data)$predictions
  RMSE_val <- RMSE(test_data$DLQI, rf_pred_test)
  RMSE_custom_val <- RMSE_custom(test_data$DLQI, rf_pred_test, X_features)
  MAPE_val <- MAPE(test_data$DLQI, rf_pred_test)
  R_squared_val <- R_square(test_data$DLQI, rf_pred_test)
  adjR_squared <- adjusted_r_squared(R_squared_val, rf_pred_test, X_features)
  
  print(paste("Final model RMSE:", RMSE_val))
  print(paste("Final model RMSE Custom:", RMSE_custom_val))
  print(paste("Final model MAPE:", MAPE_val))
  print(paste("Final model R-squared:", R_squared_val))
  print(paste("Final model Adj. R-squared:", adjR_squared))
  
  # graph
  #predict_actual(test_data$DLQI, rf_pred_test, "RF")
  
  return(list(
    final_model = final_model,
    best_params = rf_model$bestTune,
    RMSE = RMSE_val,
    RMSE_custom = RMSE_custom_val,
    MAPE = MAPE_val,
    R_squared = R_squared_val,
    rf_pred = rf_pred_test,
    adjR_squared = adjR_squared
  ))
}
RF_Ridge <- Random_forest(Ridge_select, X_Ridge_features, n_tree = 20)
RF_LASSO <- Random_forest(LASSO_select, X_LASSO_features, n_tree = 20)
RF_AL <- Random_forest(AL_select, X_AL_features, n_tree = 250)
RF_EN <- Random_forest(EN_select, X_EN_features, n_tree = 190)
RF_AEN <- Random_forest(AEN_select, X_AEN_features, n_tree = 10)

#------------ Evaluate Table
resultstable_RF <- data.frame(
  Model = c("RF-L2", "RF-L1", "RF-A-L1","RF-L1L2", "RF-A-L1L2"),
  RMSE = c(format_4(RF_Ridge$RMSE), format_4(RF_LASSO$RMSE), format_4(RF_AL$RMSE), format_4(RF_EN$RMSE), format_4(RF_AEN$RMSE)),
  RMSE_custom = c(format_4(RF_Ridge$RMSE_custom), format_4(RF_LASSO$RMSE_custom), format_4(RF_AL$RMSE_custom), format_4(RF_EN$RMSE_custom), format_4(RF_AEN$RMSE_custom)),
  MAPE = c(format_4(RF_Ridge$MAPE), format_4(RF_LASSO$MAPE), format_4(RF_AL$MAPE), format_4(RF_EN$MAPE), format_4(RF_AEN$MAPE)),
  R_squared = c(format_4(RF_Ridge$R_squared), format_4(RF_LASSO$R_squared), format_4(RF_AL$R_squared), format_4(RF_EN$R_squared), format_4(RF_AEN$R_squared)),
  Adj_R_squared = c(format_4(RF_Ridge$adjR_squared), format_4(RF_LASSO$adjR_squared), format_4(RF_AL$adjR_squared), format_4(RF_EN$adjR_squared), format_4(RF_AEN$adjR_squared))
)
resultstable_RF
datatable(resultstable_RF, 
          options = list(pageLength = 5, autoWidth = TRUE), 
          caption = "Table of Results (Random Forest)")

#------------ Evaluate Graph
RMSE_RF_model = c(RF_Ridge$RMSE, RF_LASSO$RMSE, RF_AL$RMSE, RF_EN$RMSE, RF_AEN$RMSE)
data_RMSE <- data.frame(
  Model = c("RF-L2", "RF-L1", "RF-A-L1", "RF-L1L2", "RF-A-L1L2"),
  RMSE_model = c(RF_Ridge$RMSE, RF_LASSO$RMSE, RF_AL$RMSE, RF_EN$RMSE, RF_AEN$RMSE)
)
min_rmse <- min(data_RMSE$RMSE_model)
data_RMSE$Model <- factor(data_RMSE$Model, levels = c("RF-L2", "RF-L1", "RF-A-L1", "RF-L1L2", "RF-A-L1L2"))
p1 <- ggplot(data_RMSE, aes(x = Model, y = RMSE_model, fill = RMSE_model == min_rmse)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_text(aes(label = format_4(RMSE_model)), vjust = -0.5, family = "Calisto") +
  labs(title = "RMSE (n)", x = "Model", y = "RMSE") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "#D3D3D3")) +
  geom_hline(yintercept = min_rmse, linetype = "dashed", color = "blue") +
  theme_custom

RMSE_custom_RF_model = c(RF_Ridge$RMSE_custom, RF_LASSO$RMSE_custom, RF_AL$RMSE_custom, RF_EN$RMSE_custom, RF_AEN$RMSE_custom)
data_RMSE <- data.frame(
  Model = c("RF-L2", "RF-L1", "RF-A-L1", "RF-L1L2", "RF-A-L1L2"),
  RMSE_custom_model = c(RF_Ridge$RMSE_custom, RF_LASSO$RMSE_custom, RF_AL$RMSE_custom, RF_EN$RMSE_custom, RF_AEN$RMSE_custom)
)
min_rmse_custom <- min(data_RMSE$RMSE_custom_model)
data_RMSE$Model <- factor(data_RMSE$Model, levels = c("RF-L2", "RF-L1", "RF-A-L1", "RF-L1L2", "RF-A-L1L2"))
p2 <- ggplot(data_RMSE, aes(x = Model, y = RMSE_custom_model, fill = RMSE_custom_model == min_rmse_custom)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_text(aes(label = format_4(RMSE_custom_model)), vjust = -0.5, family = "Calisto") +
  labs(title = "RMSE (n-p)", x = "Model", y = "RMSE") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "#D3D3D3")) +
  geom_hline(yintercept = min_rmse_custom, linetype = "dashed", color = "blue") +
  theme_custom

MAPE_RF_model = c(RF_Ridge$MAPE, RF_LASSO$MAPE, RF_AL$MAPE, RF_EN$MAPE, RF_AEN$MAPE)
data_MAPE <- data.frame(
  Model = c("RF-L2", "RF-L1", "RF-A-L1", "RF-L1L2", "RF-A-L1L2"),
  MAPE_model = c(RF_Ridge$MAPE, RF_LASSO$MAPE, RF_AL$MAPE, RF_EN$MAPE, RF_AEN$MAPE)
)
min_mape <- min(data_MAPE$MAPE_model)
data_MAPE$Model <- factor(data_MAPE$Model, levels = c("RF-L2", "RF-L1", "RF-A-L1", "RF-L1L2", "RF-A-L1L2"))
p3 <- ggplot(data_MAPE, aes(x = Model, y = MAPE_model, fill = MAPE_model == min_mape)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_text(aes(label = format_4(MAPE_model)), vjust = -0.5, family = "Calisto") +
  labs(title = "MAPE(%)", x = "Model", y = "MAPE (%)") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "#D3D3D3")) +
  geom_hline(yintercept = min_mape, linetype = "dashed", color = "blue") +
  theme_custom

R_squared_RF_model = c(RF_Ridge$R_squared, RF_LASSO$R_squared, RF_AL$R_squared, RF_EN$R_squared, RF_AEN$R_squared)
data_R_squared <- data.frame(
  Model = c("RF-L2", "RF-L1", "RF-A-L1", "RF-L1L2", "RF-A-L1L2"),
  R_squared_model = c(RF_Ridge$R_square, RF_LASSO$R_square, RF_AL$R_square, result_EN$R_square, RF_AEN$R_square)
)
max_r_squared <- max(data_R_squared$R_squared_model)
data_R_squared$Model <- factor(data_R_squared$Model, levels = c("RF-L2", "RF-L1", "RF-A-L1", "RF-L1L2", "RF-A-L1L2"))
p4 <- ggplot(data_R_squared, aes(x = Model, y = R_squared_model, fill = R_squared_model == max_r_squared)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_text(aes(label = format_4(R_squared_model)), vjust = -0.5, family = "Calisto") +
  labs(title = "R_squared (%)", x = "Model", y = "R_squared (%)") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "#D3D3D3")) +
  geom_hline(yintercept = max_r_squared, linetype = "dashed", color = "blue") +
  theme_custom

grid.arrange(
  p1, p3, p4,
  ncol = 2,
  top = textGrob("Random Forest", gp = gpar(fontsize = 15, fontface = "bold", family = "Calisto"))
)


#---------------------------------- Graph Compare Model
pred_RF = list(RF_Ridge$rf_pred, RF_LASSO$rf_pred, RF_AL$rf_pred, RF_EN$rf_pred, RF_AEN$rf_pred)
model_names = c("RF-L2", "RF-L1", "RF-A-L1", "RF-L1L2", "RF-A-L1L2")
n <- length(y_test)
x <- 1:n
data <- data.frame(
  x = rep(x, length(pred_RF) + 1),
  y = c(y_test, unlist(pred_RF)),
  category = rep(c("Test Data", model_names), each = n)
)

custom_colors <- c(
  "Test Data" = "grey",       
  "RF-L2" = "#61D04F",       
  "RF-L1" =  "#CD0BBC",        
  "RF-A-L1" = "#F5C710",
  "RF-L1L2" = "#8470FF",  
  "RF-A-L1L2" = "#2297E6"
)

data$category <- factor(
  data$category, 
  levels = c("Test Data", "RF-L2", "RF-L1", "RF-A-L1", "RF-L1L2", "RF-A-L1L2")
)

ggplot(data, aes(x = x, y = y)) +
  geom_line(aes(color = category, linetype = category), size = 1.5) + 
  scale_color_manual(values = custom_colors) + 
  scale_linetype_manual(
    values = c(
      "Test Data" = "dashed", 
      "RF-L2" = "solid", 
      "RF-L1" = "solid", 
      "RF-A-L1" = "solid", 
      "RF-L1L2" = "solid", 
      "RF-A-L1L2" = "solid"
    )
  ) +
  labs(
    title = "Comparison of Model Predictions (Random Forest)",
    x = "Index",
    y = "Prediction",
    color = NULL,
    linetype = NULL
  ) +
  theme_minimal() +
  theme_model

graph <- "~/Documents/PASI_pj3:2/Penalized_Report/graph_penalized/Prediction_Graph"
ggsave(
  filename = file.path(graph, "pred_RF.png"), 
  plot = last_plot(),                    
  width = 12,                                     
  height = 8,                                     
  dpi = 500                                    
)

#-----------Loop for Actual vs Predict & Actual VS Error----------#
predict_actual_plot_RF <- function(test_data, pred_data, model_name) {
  test_data <- as.numeric(unlist(test_data))
  pred_data <- as.numeric(unlist(pred_data))
  
  if (length(test_data) != length(pred_data)) stop("Vectors must have the same length.")
  median_test <- median(test_data[test_data != 0], na.rm = TRUE)
  test_data <- ifelse(test_data == 0, median_test, test_data) 
  print(paste("Median value:", median_test))
  
  analysis_data <- data.frame(
    Actual = test_data,
    Predicted = pred_data,
    Error = pred_data - test_data,
    Relative_Error = ((pred_data - test_data)/test_data)*100,
    Average_Actaul = (pred_data + test_data)/2
  )
  print(analysis_data)
  
  # Calculate Mean and ±1.96 SD for Bland-Altman
  pearson_corr <- cor(test_data, pred_data, method = "pearson", use = "complete.obs")
  mean_error <- mean(analysis_data$Relative_Error, na.rm = TRUE)
  print(paste("Mean Relative Error:", mean_error))
  sd_error <- sd(analysis_data$Relative_Error, na.rm = TRUE)
  print(paste("SD Relative Error:", sd_error))
  upper_limit <- mean_error + 1.96 * sd_error
  lower_limit <- mean_error - 1.96 * sd_error
  
  # homoscedasticity
  max_limit <- 25
  path_homoscedasticity <- "~/Documents/PASI_pj3:2/Penalized_Report/graph_penalized/Homoscedasticity/RF"
  p_homoscedasticity <- ggplot(analysis_data, aes(x = Actual, y = Predicted)) +
    geom_point(alpha = 0.7, size = 3, color = "#1F77B4") +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") + 
    coord_cartesian(xlim = c(0, max_limit), ylim = c(0, max_limit)) + 
    scale_x_continuous(breaks = seq(0, max_limit, by = 5)) + 
    scale_y_continuous(breaks = seq(0, max_limit, by = 5)) + 
    labs(
      title = paste("Actual vs Predicted (", model_name, ")", sep = ""),
      x = "Actual Y",
      y = "Predicted Y"
    ) +
    theme_minimal(base_family = "Calisto") +
    theme(
      plot.title = element_text(size = 110, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 110),
      axis.title.y = element_text(size = 110),
      axis.text = element_text(size = 110),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1.5)
    )
  
  # Bland-Altman Plot
  path_bland_altman <- "~/Documents/PASI_pj3:2/Penalized_Report/graph_penalized/Bland_altman/RF"
  p_bland_altman <- ggplot(analysis_data, aes(x = Average_Actaul, y = Relative_Error)) +
    geom_point(alpha = 0.7, size = 3, color = "#FF7F0E") +
    geom_hline(yintercept = mean_error, color = "blue", linetype = "solid", size = 1) +
    geom_hline(yintercept = upper_limit, color = "darkgreen", linetype = "dashed", size = 1) +
    geom_hline(yintercept = lower_limit, color = "darkgreen", linetype = "dashed", size = 1) +
    annotate("text", x = max(analysis_data$Average_Actaul), y = mean_error, 
             label = sprintf("Mean: %.2f", mean_error), color = "blue", hjust = 1, vjust = -1, size = 35, family = "Calisto") +
    annotate("text", x = max(analysis_data$Average_Actaul), y = upper_limit, 
             label = sprintf("+1.96SD: %.2f", upper_limit), color = "darkgreen", hjust = 1, vjust = -1, size = 35, family = "Calisto") +
    annotate("text", max(analysis_data$Average_Actaul), y = lower_limit, 
             label = sprintf("-1.96SD: %.2f", lower_limit), color = "darkgreen", hjust = 1, vjust = -1, size = 35, family = "Calisto") +
    labs(
      title = paste("Actual vs Error (", model_name, ")", sep = ""),
      x = "Mean of Actual and Predicted",
      y = "Relative Error %"
    ) +
    theme_minimal(base_family = "Calisto") +
    theme(
      plot.title = element_text(size = 110, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 110),
      axis.title.y = element_text(size = 110),
      axis.text = element_text(size = 110),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1.5)
    )
  
  ggsave(
    filename = file.path(path_homoscedasticity, paste("homoscedasticity_", model_name, ".png", sep = "")),
    plot = p_homoscedasticity,
    width = 5,
    height = 5,
    dpi = 500
  )
  ggsave(
    filename = file.path(path_bland_altman, paste("bland_altman_", model_name, ".png", sep = "")),
    plot = p_bland_altman,
    width = 7,
    height = 5,
    dpi = 500
  )
  
  return(list(homoscedasticity_plot = p_homoscedasticity, 
              bland_altman_plot = p_bland_altman,
              mean_error = format_4(mean_error),
              sd_error = format_4(sd_error),
              upper_limit = format_4(upper_limit),
              lower_limit = format_4(lower_limit),
              pearson_corr = format_4(pearson_corr)))
}
models <- c("RF-L2", "RF-L1", "RF-A-L1", "RF-L1L2", "RF-A-L1L2")
predictions <- list(RF_Ridge$rf_pred, RF_LASSO$rf_pred, RF_AL$rf_pred, RF_EN$rf_pred, RF_AEN$rf_pred)

results_RF <- list()
for (i in seq_along(models)) {
  result <- predict_actual_plot_RF(y_test, predictions[[i]], models[i])
  
  results_RF[[i]] <- data.frame(
    Model = models[i],
    Relative_Error = result$mean_error,
    SD_Error = result$sd_error,
    Upper_Limit = result$upper_limit,
    Lower_Limit = result$lower_limit,
    Pearson_Correlation = result$pearson_corr
  )
}

table_RFerror <- do.call(rbind, results_RF)
print(table_RFerror)
datatable(table_RFerror, 
          options = list(pageLength = 5, autoWidth = TRUE), 
          caption = "Table of Error (Random Forest)")



#------------------------ SVR + Penalized 5 Models ---------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------#

SVR <- function(X_select, X_features, name) {
  set.seed(845)
  select_data <- data.frame(X_select, DLQI = y)
  y <- select_data$DLQI
  X <- select_data[, !colnames(select_data) %in% "DLQI"]
  train_index <- createDataPartition(y, p = 0.8, list = FALSE)
  X_train <- X[train_index, , drop = FALSE]
  y_train <- y[train_index]
  X_test <- X[-train_index, , drop = FALSE]
  y_test <- y[-train_index]
  
  C_values <- c(0.001, 0.01, 0.1, 1, 10, 100)
  sigma_values <- c(0.0001, 0.001, 0.01, 0.1, 1, 2, 5)
  epsilon_values <- c(0.01, 0.05, 0.1, 0.5, 1)
  
  results <- data.frame(
    C = numeric(),
    sigma = numeric(),
    epsilon = numeric(),
    RMSE_grid = numeric(),
    MAPE_grid = numeric(),
    stringsAsFactors = FALSE
  )
  
  best_RMSE <- Inf
  best_params <- list()
  for (C in C_values) {
    for (sigma in sigma_values) {
      for (epsilon in epsilon_values) {
        rmse_values <- c()
        mape_values <- c()
        
        # Leave-One-Out
        for (i in 1:nrow(X_train)) {
          X_train_cv <- X_train[-i, , drop = FALSE]
          y_train_cv <- y_train[-i]
          X_test_cv <- X_train[i, , drop = FALSE]
          y_test_cv <- y_train[i]
          
          svm_model <- svm(
            formula = DLQI ~ .,
            data = data.frame(X_train_cv, DLQI = y_train_cv),
            type = "eps-regression",
            kernel = "radial",
            cost = C,
            gamma = sigma,
            epsilon = epsilon
          )
          
          svm_pred <- predict(svm_model, newdata = data.frame(X_test_cv))
          rmse_values <- c(rmse_values, RMSE(y_test_cv, svm_pred))
          mape_values <- c(mape_values, MAPE(y_test_cv, svm_pred))
        }
        
        RMSE_avg <- mean(rmse_values)
        MAPE_avg <- mean(mape_values, na.rm = TRUE)
        if (RMSE_avg < best_RMSE) {
          best_RMSE <- RMSE_avg
          best_params <- list(C = C, sigma = sigma, epsilon = epsilon,
                              RMSE_grid = RMSE_avg, MAPE_grid = MAPE_avg)
        }
        
        results <- rbind(results, data.frame(
          C = C,
          sigma = sigma,
          epsilon = epsilon,
          RMSE_grid = RMSE_avg,
          MAPE_grid = MAPE_avg
        ))
      }
    }
  }
  
  print("Best parameters found in grid search:")
  print(best_params)
  
  final_svm_model <- svm(
    formula = DLQI ~ .,
    data = data.frame(X_train, DLQI = y_train),
    type = "eps-regression",
    kernel = "sigmoid",
    cost = best_params$C,
    gamma = best_params$sigma,
    epsilon = best_params$epsilon
  )
  
  svm_pred_test <- predict(final_svm_model, newdata = data.frame(X_test))
  predict_actual(y_test, svm_pred_test, "SVR")
  
  RMSE <- RMSE(y_test, svm_pred_test)
  RMSE_custom <- RMSE_custom(y_test, svm_pred_test, X_features)
  MAPE <- MAPE(y_test, svm_pred_test)
  R_squared <- R_square(y_test, svm_pred_test)
  adjR_squared <- adjusted_r_squared(R_squared, svm_pred_test, X_features)
  
  print(paste("Final model RMSE:", RMSE))
  print(paste("Final model RMSE Custom:", RMSE_custom))
  print(paste("Final model MAPE:", MAPE))
  print(paste("Final model R-square:", R_squared))
  print(paste("Final model Adj. R-square:", adjR_squared))
  
  return(list(model = final_svm_model, 
              RMSE = RMSE, MAPE = MAPE, 
              R_squared = R_squared, 
              RMSE_custom = RMSE_custom,
              svm_pred = svm_pred_test,
              adjR_squared = adjR_squared))
}

SVR_Ridge <- SVR(Ridge_select, X_Ridge_features)
SVR_LASSO <- SVR(LASSO_select, X_LASSO_features)
SVR_AL <- SVR(AL_select, X_AL_features)
SVR_EN <- SVR(EN_select, X_EN_features)
SVR_AEN <- SVR(AEN_select, X_AEN_features)

#--------- Evaluate Table
resultstable_SVR <- data.frame(
  Model = c("SVR-L2", "SVR-L1", "SVR-A-L1","SVR-L1L2", "SVR-A-L1L2"),
  RMSE = c(format_4(SVR_Ridge$RMSE), format_4(SVR_LASSO$RMSE), format_4(SVR_AL$RMSE), format_4(SVR_EN$RMSE), format_4(SVR_AEN$RMSE)),
  RMSE_custom = c(format_4(SVR_Ridge$RMSE_custom), format_4(SVR_LASSO$RMSE_custom), format_4(SVR_AL$RMSE_custom), format_4(SVR_EN$RMSE_custom), format_4(SVR_AEN$RMSE_custom)),
  MAPE = c(format_4(SVR_Ridge$MAPE), format_4(SVR_LASSO$MAPE), format_4(SVR_AL$MAPE), format_4(SVR_EN$MAPE), format_4(SVR_AEN$MAPE)),
  R_squared = c(format_4(SVR_Ridge$R_squared), format_4(SVR_LASSO$R_squared), format_4(SVR_AL$R_squared), format_4(SVR_EN$R_squared), format_4(SVR_AEN$R_squared)),
  Adj_R_squared = c(format_4(SVR_Ridge$adjR_squared), format_4(SVR_LASSO$adjR_squared), format_4(SVR_AL$adjR_squared), format_4(SVR_EN$adjR_squared), format_4(SVR_AEN$adjR_squared))
)
resultstable_SVR
datatable(resultstable_SVR, 
          options = list(pageLength = 5, autoWidth = TRUE), 
          caption = "Table of Results (Random Forest)")

#-------- Evaluate Graph
RMSE_SVR_model = c(SVR_Ridge$RMSE, SVR_LASSO$RMSE, SVR_AL$RMSE, SVR_EN$RMSE, SVR_AEN$RMSE)
data_RMSE <- data.frame(
  Model = c("SVR-L2", "SVR-L1", "SVR-A-L1", "SVR-L1L2", "SVR-A-L1L2"),
  RMSE_model = c(SVR_Ridge$RMSE, SVR_LASSO$RMSE, SVR_AL$RMSE, SVR_EN$RMSE, SVR_AEN$RMSE)
)
min_rmse <- min(data_RMSE$RMSE_model)
data_RMSE$Model <- factor(data_RMSE$Model, levels = c("SVR-L2", "SVR-L1", "SVR-A-L1", "SVR-L1L2", "SVR-A-L1L2"))
p1 <- ggplot(data_RMSE, aes(x = Model, y = RMSE_model, fill = RMSE_model == min_rmse)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_text(aes(label = format_4(RMSE_model)), vjust = -0.5, family = "Calisto") +
  labs(title = "RMSE (n)", x = "Model", y = "RMSE") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "#D3D3D3")) +
  geom_hline(yintercept = min_rmse, linetype = "dashed", color = "blue") +
  theme_custom

RMSE_custom_SVR_model = c(SVR_Ridge$RMSE_custom, SVR_LASSO$RMSE_custom, SVR_AL$RMSE_custom, SVR_EN$RMSE_custom, SVR_AEN$RMSE_custom)
data_RMSE <- data.frame(
  Model = c("SVR-L2", "SVR-L1", "SVR-A-L1", "SVR-L1L2", "SVR-A-L1L2"),
  RMSE_custom_model = c(SVR_Ridge$RMSE_custom, SVR_LASSO$RMSE_custom, SVR_AL$RMSE_custom, SVR_EN$RMSE_custom, SVR_AEN$RMSE_custom)
)
min_rmse_custom <- min(data_RMSE$RMSE_custom_model)
data_RMSE$Model <- factor(data_RMSE$Model, levels = c("SVR-L2", "SVR-L1", "SVR-A-L1", "SVR-L1L2", "SVR-A-L1L2"))
p2 <- ggplot(data_RMSE, aes(x = Model, y = RMSE_custom_model, fill = RMSE_custom_model == min_rmse_custom)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_text(aes(label = format_4(RMSE_custom_model)), vjust = -0.5, family = "Calisto") +
  labs(title = "RMSE (n-p)", x = "Model", y = "RMSE") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "#D3D3D3")) +
  geom_hline(yintercept = min_rmse_custom, linetype = "dashed", color = "blue") +
  theme_custom

MAPE_SVR_model = c(SVR_Ridge$MAPE, SVR_LASSO$MAPE, SVR_AL$MAPE, SVR_EN$MAPE, SVR_AEN$MAPE)
data_MAPE <- data.frame(
  Model = c("SVR-L2", "SVR-L1", "SVR-A-L1", "SVR-L1L2", "SVR-A-L1L2"),
  MAPE_model = c(SVR_Ridge$MAPE, SVR_LASSO$MAPE, SVR_AL$MAPE, SVR_EN$MAPE, SVR_AEN$MAPE)
)
min_mape <- min(data_MAPE$MAPE_model)
data_MAPE$Model <- factor(data_MAPE$Model, levels = c("SVR-L2", "SVR-L1", "SVR-A-L1", "SVR-L1L2", "SVR-A-L1L2"))
p3 <- ggplot(data_MAPE, aes(x = Model, y = MAPE_model, fill = MAPE_model == min_mape)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_text(aes(label = format_4(MAPE_model)), vjust = -0.5, family = "Calisto") +
  labs(title = "MAPE(%)", x = "Model", y = "MAPE (%)") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "#D3D3D3")) +
  geom_hline(yintercept = min_mape, linetype = "dashed", color = "blue") +
  theme_custom

R_squared_SVR_model = c(SVR_Ridge$R_squared, SVR_LASSO$R_squared, SVR_AL$R_squared, SVR_EN$R_squared, SVR_AEN$R_squared)
data_R_squared <- data.frame(
  Model = c("SVR-L2", "SVR-L1", "SVR-A-L1", "SVR-L1L2", "SVR-A-L1L2"),
  R_squared_model = c(SVR_Ridge$R_square, SVR_LASSO$R_square, SVR_AL$R_square, SVR_EN$R_square, SVR_AEN$R_square)
)
max_r_squared <- max(data_R_squared$R_squared_model)
data_R_squared$Model <- factor(data_R_squared$Model, levels = c("SVR-L2", "SVR-L1", "SVR-A-L1", "SVR-L1L2", "SVR-A-L1L2"))
p4 <- ggplot(data_R_squared, aes(x = Model, y = R_squared_model, fill = R_squared_model == max_r_squared)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_text(aes(label = format_4(R_squared_model)), vjust = -0.5, family = "Calisto") +
  labs(title = "R_squared (%)", x = "Model", y = "R_squared (%)") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "#D3D3D3")) +
  geom_hline(yintercept = max_r_squared, linetype = "dashed", color = "blue") +
  theme_custom

grid.arrange(
  p1, p3, p4,
  ncol = 2,
  top = textGrob("Support Vector Regression", gp = gpar(fontsize = 15, fontface = "bold", family = "Calisto"))
)


#------------------------------- Graph Compare Model
pred_SVR = list(SVR_Ridge$svm_pred, SVR_LASSO$svm_pred, SVR_AL$svm_pred, SVR_EN$svm_pred, SVR_AEN$svm_pred)
model_names = c("SVR-L2", "SVR-L1", "SVR-A-L1", "SVR-L1L2", "SVR-A-L1L2")
n <- length(y_test)
x <- 1:n
data <- data.frame(
  x = rep(x, length(pred_SVR) + 1),
  y = c(y_test, unlist(pred_SVR)),
  category = rep(c("Test Data", model_names), each = n)
)

custom_colors <- c(
  "Test Data" = "grey",       
  "SVR-L2" = "#61D04F",       
  "SVR-L1" =  "#CD0BBC",        
  "SVR-A-L1" = "#F5C710",
  "SVR-L1L2" = "#8470FF",  
  "SVR-A-L1L2" = "#2297E6"
)


data$category <- factor(
  data$category, 
  levels = c("Test Data", "SVR-L2", "SVR-L1", "SVR-A-L1", "SVR-L1L2", "SVR-A-L1L2")
)

ggplot(data, aes(x = x, y = y)) +
  geom_line(aes(color = category, linetype = category), size = 1.5) + 
  scale_color_manual(values = custom_colors) + 
  scale_linetype_manual(
    values = c(
      "Test Data" = "dashed", 
      "SVR-L2" = "solid", 
      "SVR-L1" = "solid", 
      "SVR-A-L1" = "solid", 
      "SVR-L1L2" = "solid", 
      "SVR-A-L1L2" = "solid"
    )
  ) +
  labs(
    title = "Comparison of Model Predictions (SVR)",
    x = "Index",
    y = "Prediction",
    color = NULL,
    linetype = NULL
  ) +
  theme_minimal() +
  theme_model

graph <- "~/Documents/PASI_pj3:2/Penalized_Report/graph_penalized/Prediction_Graph"
ggsave(
  filename = file.path(graph, "pred_SVR.png"), 
  plot = last_plot(),                    
  width = 12,                                     
  height = 8,                                     
  dpi = 500                                  
)

#-----------Loop for Actual vs Predict & Actual VS Error----------#
predict_actual_plot_SVR <- function(test_data, pred_data, model_name) {
  test_data <- as.numeric(unlist(test_data))
  pred_data <- as.numeric(unlist(pred_data))
  
  if (length(test_data) != length(pred_data)) stop("Vectors must have the same length.")
  median_test <- median(test_data[test_data != 0], na.rm = TRUE)
  test_data <- ifelse(test_data == 0, median_test, test_data) 
  print(paste("Median value:", median_test))
  
  analysis_data <- data.frame(
    Actual = test_data,
    Predicted = pred_data,
    Error = pred_data - test_data,
    Relative_Error = ((pred_data - test_data)/test_data)*100,
    Average_Actaul = (pred_data + test_data)/2
  )
  print(analysis_data)
  
  # Calculate Mean and ±1.96 SD for Bland-Altman
  pearson_corr <- cor(test_data, pred_data, method = "pearson", use = "complete.obs")
  mean_error <- mean(analysis_data$Relative_Error, na.rm = TRUE)
  print(paste("Mean Relative Error:", mean_error))
  sd_error <- sd(analysis_data$Relative_Error, na.rm = TRUE)
  print(paste("SD Relative Error:", sd_error))
  upper_limit <- mean_error + 1.96 * sd_error
  lower_limit <- mean_error - 1.96 * sd_error
  
  # homoscedasticity
  max_limit <- 25
  path_homoscedasticity <- "~/Documents/PASI_pj3:2/Penalized_Report/graph_penalized/Homoscedasticity/SVR"
  p_homoscedasticity <- ggplot(analysis_data, aes(x = Actual, y = Predicted)) +
    geom_point(alpha = 0.7, size = 3, color = "#1F77B4") +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") + 
    coord_cartesian(xlim = c(0, max_limit), ylim = c(0, max_limit)) + 
    scale_x_continuous(breaks = seq(0, max_limit, by = 5)) + 
    scale_y_continuous(breaks = seq(0, max_limit, by = 5)) + 
    labs(
      title = paste("Actual vs Predicted (", model_name, ")", sep = ""),
      x = "Actual Y",
      y = "Predicted Y"
    ) +
    theme_minimal(base_family = "Calisto") +
    theme(
      plot.title = element_text(size = 110, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 110),
      axis.title.y = element_text(size = 110),
      axis.text = element_text(size = 110),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1.5)
    )
  
  # Bland-Altman Plot
  path_bland_altman <- "~/Documents/PASI_pj3:2/Penalized_Report/graph_penalized/Bland_altman/SVR"
  p_bland_altman <- ggplot(analysis_data, aes(x = Average_Actaul, y = Relative_Error)) +
    geom_point(alpha = 0.7, size = 3, color = "#FF7F0E") +
    geom_hline(yintercept = mean_error, color = "blue", linetype = "solid", size = 1) +
    geom_hline(yintercept = upper_limit, color = "darkgreen", linetype = "dashed", size = 1) +
    geom_hline(yintercept = lower_limit, color = "darkgreen", linetype = "dashed", size = 1) +
    annotate("text", x = max(analysis_data$Average_Actaul), y = mean_error, 
             label = sprintf("Mean: %.2f", mean_error), color = "blue", hjust = 1, vjust = -1, size = 35, family = "Calisto") +
    annotate("text", x = max(analysis_data$Average_Actaul), y = upper_limit, 
             label = sprintf("+1.96SD: %.2f", upper_limit), color = "darkgreen", hjust = 1, vjust = -1, size = 35, family = "Calisto") +
    annotate("text", max(analysis_data$Average_Actaul), y = lower_limit, 
             label = sprintf("-1.96SD: %.2f", lower_limit), color = "darkgreen", hjust = 1, vjust = -1, size = 35, family = "Calisto") +
    labs(
      title = paste("Actual vs Error (", model_name, ")", sep = ""),
      x = "Mean of Actual and Predicted",
      y = "Relative Error %"
    ) +
    theme_minimal(base_family = "Calisto") +
    theme(
      plot.title = element_text(size = 110, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 110),
      axis.title.y = element_text(size = 110),
      axis.text = element_text(size = 110),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1.5)
    )
  
  ggsave(
    filename = file.path(path_homoscedasticity, paste("homoscedasticity_", model_name, ".png", sep = "")),
    plot = p_homoscedasticity,
    width = 5,
    height = 5,
    dpi = 500
  )
  ggsave(
    filename = file.path(path_bland_altman, paste("bland_altman_", model_name, ".png", sep = "")),
    plot = p_bland_altman,
    width = 7,
    height = 5,
    dpi = 500
  )
  
  return(list(homoscedasticity_plot = p_homoscedasticity, 
              bland_altman_plot = p_bland_altman,
              mean_error = format_4(mean_error),
              sd_error = format_4(sd_error),
              upper_limit = format_4(upper_limit),
              lower_limit = format_4(lower_limit),
              pearson_corr = format_4(pearson_corr)))
}
models <- c("SVR-L2", "SVR-L1", "SVR-A-L1", "SVR-L1L2", "SVR-A-L1L2")
predictions <- list(SVR_Ridge$svm_pred, SVR_LASSO$svm_pred, SVR_AL$svm_pred, SVR_EN$svm_pred, SVR_AEN$svm_pred)

results_SVR <- list()
for (i in seq_along(models)) {
  result <- predict_actual_plot_SVR(y_test, predictions[[i]], models[i])
  
  results_SVR[[i]] <- data.frame(
    Model = models[i],
    Relative_Error = result$mean_error,
    SD_Error = result$sd_error,
    Upper_Limit = result$upper_limit,
    Lower_Limit = result$lower_limit,
    Pearson_Correlation = result$pearson_corr
  )
}

table_SVRerror <- do.call(rbind, results_SVR)
print(table_SVRerror)
datatable(table_SVRerror, 
          options = list(pageLength = 5, autoWidth = TRUE), 
          caption = "Table of Error (Support Vector Regression)")


#------------------------ Compaire 15 Models (Std.numeric & Spilt 80:20 & LOOCV model) ---------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------#
y_pred_true <- data.frame(
  y_true = y_test,
  pred_LR = format_4(LR_pred),
  pred_L2 = format_4(Ridge_pred),
  pred_L1 = format_4(LASSO_pred),
  `pred_A-L1` = format_4(AL_pred),
  `pred_L1L2` = format_4(EN_pred),
  `pred_A-L1L2` = format_4(AEN_pred),
  `pred_RF-L2` = format_4(RF_Ridge$rf_pred),
  `pred_RF-L1` = format_4(RF_LASSO$rf_pred),
  `pred_RF-A-L1` = format_4(RF_AL$rf_pred),
  `pred_RF-L1L2` = format_4(RF_EN$rf_pred),
  `pred_RF-A-L1L2` = format_4(RF_AEN$rf_pred),
  `pred_SVR-L2` = format_4(SVR_Ridge$svm_pred),
  `pred_SVR-L1` = format_4(SVR_LASSO$svm_pred),
  `pred_SVR-A-L1` = format_4(SVR_AL$svm_pred),
  `pred_SVR-L1L2` = format_4(SVR_EN$svm_pred),
  `pred_SVR-A-L1L2` = format_4(SVR_AEN$svm_pred)
)

datatable(y_pred_true, 
          options = list(pageLength = 5, autoWidth = TRUE), 
          caption = "Table of y true & y-hat (Linear Regression)")


combined_error <- rbind(table_Penalerror, table_RFerror, table_SVRerror)
datatable(combined_error, 
          options = list(pageLength = 5, autoWidth = TRUE), 
          caption = "Table of Combined Error")

resultstable_penal <- resultstable_penal[, !(colnames(resultstable_penal) %in% c("num_var", "Lambda"))]
combined_result <- rbind(resultstable_penal, resultstable_RF, resultstable_SVR)
datatable(combined_result,
          options = list(pageLength = 5, autoWidth = TRUE),
          caption = "Table of Combined Results")


set3_colors <- colorRampPalette(brewer.pal(12, "Set3"))(10)
paired_colors <- brewer.pal(8, "Paired")
all_colors <- c(paired_colors, set3_colors)
custom_colors <- setNames(all_colors, c("LR", "L2", "L1", "A-L1", "L1L2", "A-L1L2", 
                                        "RF-L2", "RF-L1", "RF-A-L1", "RF-L1L2", "RF-A-L1L2", 
                                        "SVR-L2", "SVR-L1", "SVR-A-L1", "SVR-L1L2", "SVR-A-L1L2"))
theme_compare <- theme(
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  axis.title.x = element_text(size = 130, family = "Calisto"), 
  axis.text.x = element_text(size = 90, angle = 60, hjust = 1, family = "Calisto"), 
  plot.title = element_text(size = 120, face = "bold", hjust = 0.5, family = "Calisto"), 
  legend.position = "right",  
  legend.title = element_text(size = 120, face = "bold", family = "Calisto"), 
  legend.text = element_text(size = 120, family = "Calisto"),
  panel.grid = element_blank()
)

# RMSE
data <- data.frame(
  Model = c("LR", "L2", "L1", "A-L1", "L1L2", "A-L1L2", 
            "RF-L2", "RF-L1", "RF-A-L1", "RF-L1L2", "RF-A-L1L2", 
            "SVR-L2", "SVR-L1", "SVR-A-L1", "SVR-L1L2", "SVR-A-L1L2"),
  RMSE = c(RMSE_penal_model, RMSE_RF_model, RMSE_SVR_model)
)

min_rmse <- min(data$RMSE)
data$Model <- factor(data$Model, levels = data$Model[order(data$RMSE)])
ggplot(data, aes(x = Model, y = RMSE, fill = RMSE)) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  geom_text(aes(label = format_4(RMSE)),  # แสดงตัวเลขค่าเดิม
            position = position_dodge(width = 0.8), 
            vjust = -0.5,
            size = 25,
            family = "Calisto") +  
  theme_minimal() +
  scale_fill_gradient(low = "skyblue", high = "darkblue", 
                      name = "Value") +
  labs(title = "RMSE Comparison for 15 Models", 
       x = "Predictive Model", 
       y = NULL) +
  theme_compare
graph <- "~/Documents/PASI_pj3:2/Penalized_Report/graph_penalized/Comparison_Loss"
ggsave(
  filename = file.path(graph, "compare_RMSE.png"), 
  plot = last_plot(),                  
  width = 12,                                     
  height = 8,                                     
  dpi = 500                                  
)


# RMSE Custom
data <- data.frame(
  Model = c("LR", "L2", "L1", "A-L1", "L1L2", "A-L1L2", 
            "RF-L2", "RF-L1", "RF-A-L1", "RF-L1L2", "RF-A-L1L2", 
            "SVR-L2", "SVR-L1", "SVR-A-L1", "SVR-L1L2", "SVR-A-L1L2"),
  RMSE_custom = c(RMSE_custom_penal_model, RMSE_custom_RF_model, RMSE_custom_SVR_model)
)
min_rmse <- min(data$RMSE_custom)
data$Model <- factor(data$Model, levels = data$Model[order(data$RMSE_custom)])
ggplot(data, aes(x = Model, y = RMSE_custom, fill = RMSE_custom)) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  geom_text(aes(label = format_4(RMSE_custom)),  # แสดงตัวเลขค่าเดิม
            position = position_dodge(width = 0.8), 
            vjust = -0.5,
            size = 25,
            family = "Calisto") +  
  theme_minimal() +
  scale_fill_gradient(low = "skyblue", high = "darkblue", 
                      name = "Value") +
  labs(title = "RMSE_custom (%) Comparison for 15 Models (Ascending Order)", 
       x = "Predictive Model", 
       y = NULL) +
  theme_compare
graph <- "~/Documents/PASI_pj3:2/Penalized_Report/graph_penalized/Comparison_Loss"
ggsave(
  filename = file.path(graph, "compare_RMSE_custom.png"), 
  plot = last_plot(),                  
  width = 12,                                     
  height = 8,                                     
  dpi = 500                                  
)


# MAPE
data <- data.frame(
  Model = c("LR", "L2", "L1", "A-L1", "L1L2", "A-L1L2", 
            "RF-L2", "RF-L1", "RF-A-L1", "RF-L1L2", "RF-A-L1L2", 
            "SVR-L2", "SVR-L1", "SVR-A-L1", "SVR-L1L2", "SVR-A-L1L2"),
  MAPE = c(MAPE_penal_model, MAPE_RF_model, MAPE_SVR_model)
)
min_mape <- min(data$MAPE)
data$Model <- factor(data$Model, levels = data$Model[order(data$MAPE)])
ggplot(data, aes(x = Model, y = MAPE, fill = MAPE)) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  geom_text(aes(label = format_4(MAPE)),
            position = position_dodge(width = 0.8), 
            vjust = -0.5,
            size = 25,
            family = "Calisto") +  
  theme_minimal() +
  scale_fill_gradient(low = "skyblue", high = "darkblue", 
                      name = "Value") +
  labs(title = "MAPE(%) Comparison for 15 Models", 
       x = "Predictive Model", 
       y = NULL) +
  theme_compare
graph <- "~/Documents/PASI_pj3:2/Penalized_Report/graph_penalized/Comparison_Loss"
ggsave(
  filename = file.path(graph, "compare_MAPE.png"), 
  plot = last_plot(),                  
  width = 12,                                     
  height = 8,                                     
  dpi = 500                                  
)


# R-squared
data <- data.frame(
  Model = c("LR", "L2", "L1", "A-L1", "L1L2", "A-L1L2", 
            "RF-L2", "RF-L1", "RF-A-L1", "RF-L1L2", "RF-A-L1L2", 
            "SVR-L2", "SVR-L1", "SVR-A-L1", "SVR-L1L2", "SVR-A-L1L2"),
  R_squared = c(R_squared_penal_model, R_squared_RF_model, R_squared_SVR_model)
)

max_R_squared <- max(data$R_squared)
data$Model <- factor(data$Model, levels = data$Model[order(-data$R_squared)])
ggplot(data, aes(x = Model, y = R_squared, fill = R_squared)) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  geom_text(aes(label = format_4(R_squared)),  # แสดงตัวเลขค่าเดิม
            position = position_dodge(width = 0.8), 
            vjust = -0.5,
            size = 25,
            family = "Calisto") +  
  theme_minimal() +
  scale_fill_gradient(low = "skyblue", high = "darkblue", 
                      name = "Value") +
  labs(title = "R_squared (%) Comparison for 15 Models", 
       x = "Predictive Model", 
       y = NULL) +
  theme_compare
graph <- "~/Documents/PASI_pj3:2/Penalized_Report/graph_penalized/Comparison_Loss"
ggsave(
  filename = file.path(graph, "compare_Rsquare.png"), 
  plot = last_plot(),                  
  width = 12,                                     
  height = 8,                                     
  dpi = 500                                  
)

#-----------------------------------------------------------------------------------------------#
#------------------------------ Shap value with test set-----------------------------------------------#
#-----------------------------------------------------------------------------------------------#
Random_forest <- function(X_select, X_features, n_tree) {
  set.seed(845)
  select_data <- data.frame(X_select, DLQI = y)
  
  train_index <- createDataPartition(select_data$DLQI, p = 0.8, list = FALSE)
  train_data <- select_data[train_index, ]
  test_data <- select_data[-train_index, ]
  
  num_trees <- n_tree
  rf_grid <- expand.grid(
    mtry = seq(2, min(ncol(train_data) - 1, 15), by = 1),
    splitrule = c("extratrees"),
    min.node.size = c(5, 8, 10, 12, 14, 16, 18, 20)
  )
  
  train_control <- trainControl(
    method = "LOOCV",
    search = "grid"
  )
  
  rf_model <- train(
    DLQI ~ .,
    data = train_data,
    method = "ranger",
    trControl = train_control,
    tuneGrid = rf_grid,
    num.trees = num_trees
  )
  
  print("Best parameters found in grid search:")
  print(rf_model$bestTune)
  
  final_model <- ranger(
    formula = DLQI ~ .,
    data = train_data,
    mtry = rf_model$bestTune$mtry,
    splitrule = rf_model$bestTune$splitrule,
    min.node.size = rf_model$bestTune$min.node.size,
    num.trees = num_trees
  )
  
  rf_pred_test <- predict(final_model, data = test_data)$predictions
  
  RMSE_val <- RMSE(test_data$DLQI, rf_pred_test)
  RMSE_custom_val <- RMSE_custom(test_data$DLQI, rf_pred_test, X_features)
  MAPE_val <- MAPE(test_data$DLQI, rf_pred_test)
  R_squared_val <- R_square(test_data$DLQI, rf_pred_test)
  adjR_squared <- adjusted_r_squared(R_squared_val, rf_pred_test, X_features)
  
  print(paste("Final model RMSE:", RMSE_val))
  print(paste("Final model RMSE Custom:", RMSE_custom_val))
  print(paste("Final model MAPE:", MAPE_val))
  print(paste("Final model R-squared:", R_squared_val))
  print(paste("Final model Adj. R-squared:", adjR_squared))
  
  ## -------- SHAP Value --------
  X_train_only <- train_data[, setdiff(names(train_data), "DLQI")]
  y_train_only <- train_data$DLQI
  
  predictor <- Predictor$new(
    model = final_model,
    data = X_train_only,
    y = y_train_only
  )
  
  # calculate SHAP in test set
  X_test_only <- test_data[, setdiff(names(test_data), "DLQI")]
  shap_list <- list()
  shap_df_list <- list()
  
  for (i in 1:nrow(X_test_only)) {
    shap <- Shapley$new(predictor, x.interest = X_test_only[i, , drop = FALSE])
    shap_list[[i]] <- shap
    shap_vector <- setNames(shap$results$phi, shap$results$feature)
    shap_df_list[[i]] <- as.data.frame(t(shap_vector))
  }
  
  shap_df <- do.call(rbind, shap_df_list)
  rownames(shap_df) <- paste0("Sample_", 1:nrow(shap_df))
  
  return(list(
    final_model = final_model,
    best_params = rf_model$bestTune,
    RMSE = RMSE_val,
    RMSE_custom = RMSE_custom_val,
    MAPE = MAPE_val,
    R_squared = R_squared_val,
    rf_pred = rf_pred_test,
    adjR_squared = adjR_squared,
    shap_values = shap_list,
    shap_df = shap_df
  ))
}
RF_EN <- Random_forest(EN_select, X_EN_features, n_tree = 190)
RF_EN$shap_df
RF_EN$shap_values

# Bar chart (mean |SHAP value|)
shap_mean_df <- data.frame(
  feature = colnames(RF_EN$shap_df),
  mean_abs_shap = apply(abs(RF_EN$shap_df), 2, mean)
)
shap_mean_df <- shap_mean_df %>%
  arrange(desc(mean_abs_shap)) %>%
  mutate(
    rank = row_number(),
    color = ifelse(rank <= 4, "#40CC89", "#3399CC")  # green for top 4, bluegrey otherwise
  )
shap_mean_df <- shap_mean_df %>%
  mutate(feature = recode(feature,
                          "Stress6" = "STR06",
                          "Stress8" = "STR08",
                          "Stress11" = "STR11",
                          "other.disease2" = "COM2",
                          "Gender2" = "GEN2",
                          "Age" = "AGE",
                          "PASI" = "PASI",
                          "Stress3" = "STR03"
  ))

xmin <- 0.5
xmax <- nrow(shap_mean_df) + 0.5
ymin <- 0
ymax <- max(shap_mean_df$mean_abs_shap) * 1.15
ggplot(shap_mean_df, aes(x = reorder(feature, mean_abs_shap), y = mean_abs_shap, fill = color)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = sprintf("%.3f", mean_abs_shap)), 
            hjust = -0.1, size = 35, family = "Calisto", fontface = "bold") +  
  scale_fill_identity() +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title = "Mean Absolute SHAP Values",
    x = "Feature",
    y = "Mean |SHAP|"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Calisto", face = "bold", size = 130, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 130, color = "black"),
    #panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    axis.text = element_text(color = "black")
  ) +
  annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
           color = "black", fill = NA, size = 1)

graph <- "~/Downloads"
ggsave(
  filename = file.path(graph, "shap.png"), 
  plot = last_plot(),                  
  width = 12,                                     
  height = 8,                                     
  dpi = 500                                  
)

# the Beeswarm plot
y_actual <- data$DLQI
X_actual <- data[, !colnames(data) %in% "DLQI"]
X_train_act <- X_actual[train_index, , drop = FALSE]
y_train_act <- y_actual[train_index]
X_test_act <- X_actual[-train_index, , drop = FALSE]
y_test_act <- y_actual[-train_index]

shap_long <- RF_EN$shap_df %>%
  mutate(sample_id = row_number()) %>%
  pivot_longer(-sample_id, names_to = "feature", values_to = "shap_value")

X_test <- RF_EN$shap_values[[1]]$x.interest
for (i in 2:length(RF_EN$shap_values)) {
  X_test <- rbind(X_test, RF_EN$shap_values[[i]]$x.interest)}

X_test_long <- X_test %>%
  mutate(sample_id = row_number()) %>%
  pivot_longer(-sample_id, names_to = "feature", values_to = "feature_value")

shap_plot_data <- shap_long %>%
  left_join(X_test_long, by = c("sample_id", "feature"))

shap_plot_data <- shap_plot_data %>%
  mutate(feature = recode(feature,
                          "Stress6" = "STR06",
                          "Stress8" = "STR08",
                          "Stress11" = "STR11",
                          "other.disease2" = "COM2",
                          "Gender2" = "GEN2",
                          "Age" = "AGE",
                          "PASI" = "PASI",
                          "Stress3" = "STR03"
  ))

shap_plot_data <- shap_plot_data %>%
  group_by(feature) %>%
  mutate(mean_abs_shap = mean(abs(shap_value))) %>%
  ungroup() %>%
  mutate(feature = reorder(feature, mean_abs_shap))

ggplot(shap_plot_data, aes(x = shap_value, y = feature, color = feature_value)) +
  geom_blank() +
  geom_jitter(height = 0.25, size = 3, alpha = 0.7) +
  geom_vline(xintercept = 0, color = "grey40", size = 0.7) +  
  scale_color_gradientn(
    colors = c("#0000FF", "#00FFFF", "#00FF80", "#00FF00"), name = "Feature value"
  ) +
  scale_x_continuous(breaks = seq(-3, 3, by = 0.5)) + 
  labs(
    title = "SHAP value (impact on model output)",
    x = "SHAP value (impact on model output)",
    y = NULL
  ) +
  theme_minimal() +
  guides(title = NULL, color = guide_colorbar(barheight = unit(15, "cm"), barwidth = unit(0.6, "cm"))) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
    panel.grid.minor.y = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    plot.title = element_text(family = "Calisto", face = "bold", color = "black", size = 130, hjust = 0.5),
    axis.text = element_text(family = "Calisto", color = "black", size = 100),
    axis.title.x = element_text(family = "Calisto", color = "black", size = 100),
    legend.text = element_text(family = "Calisto", color = "black", size = 80),
    legend.title = element_text(family = "Calisto", color = "black", size = 80)
  )

graph <- "~/Downloads"
ggsave(
  filename = file.path(graph, "impact.png"), 
  plot = last_plot(),                  
  width = 12,                                     
  height = 8,                                     
  dpi = 500                                  
)



# PASI Ratio
set.seed(845)         
y_actual <- data$DLQI
X_actual <- data[, !colnames(data) %in% "DLQI"]
X_train_act <- X_actual[train_index, , drop = FALSE]
y_train_act <- y_actual[train_index]
X_test_act  <- X_actual[-train_index, , drop = FALSE] 
y_test_act  <- y_actual[-train_index]

shap_long <- RF_EN$shap_df |>
  mutate(sample_id = row_number()) |>
  pivot_longer(-sample_id,
               names_to  = "feature",
               values_to = "shap_value")

numeric_feats <- c("PASI", "Age") 
X_test_long <- X_test_act |>
  mutate(sample_id = row_number()) |>
  select(sample_id, all_of(numeric_feats)) |>
  pivot_longer(-sample_id,
               names_to  = "feature",
               values_to = "feature_value")

shap_plot_data <- shap_long |>
  filter(feature %in% numeric_feats) |>     
  left_join(X_test_long, by = c("sample_id", "feature")) |>
  group_by(feature) |>
  mutate(mean_abs_shap = mean(abs(shap_value))) |>
  ungroup() |>
  mutate(feature = reorder(feature, mean_abs_shap)) 

shap_plot_data %>% filter(feature == "PASI")
ggplot(shap_plot_data, aes(x = feature_value, y = shap_value, color = shap_value)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_hline(yintercept = 0, color = "gray70", linetype = "dashed", size = 0.8) +
  scale_color_gradientn(colors = c("#0000FF", "#00FFFF", "#00FF80", "#00FF00")) +
  scale_y_continuous(breaks = seq(-2, 2, by = 0.5)) +
  scale_x_continuous(breaks = seq(0, 70, by = 10)) +
  labs(
    title = "SHAP value vs. PASI",
    x = "PASI",
    y = "SHAP value for PASI"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(family = "Calisto", hjust = 0.5, face = "bold", size = 130),
    axis.text = element_text(family = "Calisto", color = "black", size = 120),
    axis.title.x = element_text(family = "Calisto", color = "black", size = 120),
    axis.title.y = element_text(family = "Calisto", color = "black", size = 120),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid = element_blank()
  )

graph <- "~/Downloads"
ggsave(
  filename = file.path(graph, "PASI.png"), 
  plot = last_plot(),                  
  width = 12,                                     
  height = 8,                                     
  dpi = 500                                  
)


# Age Ratio
shap_plot_data %>% filter(feature == "Age")
ggplot(shap_plot_data, aes(x = feature_value, y = shap_value, color = shap_value)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_hline(yintercept = 0, color = "gray70", linetype = "dashed", size = 0.8) +
  scale_color_gradientn(colors = c("#0000FF", "#00FFFF", "#00FF80", "#00FF00")) +
  scale_y_continuous(breaks = seq(-2, 2, by = 0.5)) +
  scale_x_continuous(breaks = seq(0, 70, by = 10)) +
  labs(
    title = "SHAP value vs. AGE",
    x = "AGE",
    y = "SHAP value for AGE"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(family = "Calisto", hjust = 0.5, face = "bold", size = 130),
    axis.text = element_text(family = "Calisto", color = "black", size = 120),
    axis.title.x = element_text(family = "Calisto", color = "black", size = 120),
    axis.title.y = element_text(family = "Calisto", color = "black", size = 120),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid = element_blank()
  )

graph <- "~/Downloads"
ggsave(
  filename = file.path(graph, "Age.png"), 
  plot = last_plot(),                  
  width = 12,                                     
  height = 8,                                     
  dpi = 500                                  
)


shap_plot_data %>% filter(feature == "Age")
ggplot(shap_plot_data, aes(x = feature_value, y = shap_value, color = shap_value)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_hline(yintercept = 0, color = "gray70", linetype = "dashed", size = 0.8) +
  scale_color_gradientn(colors = c("#0000FF", "#00FFFF", "#00FF80", "#00FF00")) +
  scale_y_continuous(breaks = seq(-2, 2, by = 0.5)) +
  scale_x_continuous(breaks = seq(0, 70, by = 10)) +
  labs(
    title = "SHAP value vs. AGE",
    x = "AGE",
    y = "SHAP value for AGE"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(family = "Calisto", hjust = 0.5, face = "bold", size = 130),
    axis.text = element_text(family = "Calisto", color = "black", size = 120),
    axis.title.x = element_text(family = "Calisto", color = "black", size = 120),
    axis.title.y = element_text(family = "Calisto", color = "black", size = 120),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid = element_blank()
  )

graph <- "~/Downloads"
ggsave(
  filename = file.path(graph, "Age.png"), 
  plot = last_plot(),                  
  width = 12,                                     
  height = 8,                                     
  dpi = 500                                  
)


# PASI Ratio
set.seed(845)         
y_actual <- data$DLQI
X_actual <- data[, !colnames(data) %in% "DLQI"]
X_train_act <- X_actual[train_index, , drop = FALSE]
y_train_act <- y_actual[train_index]
X_test_act  <- X_actual[-train_index, , drop = FALSE] 
y_test_act  <- y_actual[-train_index]

shap_long <- RF_EN$shap_df |>
  mutate(sample_id = row_number()) |>
  pivot_longer(-sample_id,
               names_to  = "feature",
               values_to = "shap_value")

numeric_feats <- c("PASI", "Age") 
X_test_long <- X_test_act |>
  mutate(sample_id = row_number()) |>
  select(sample_id, all_of(numeric_feats)) |>
  pivot_longer(-sample_id,
               names_to  = "feature",
               values_to = "feature_value")

shap_plot_data <- shap_long |>
  filter(feature %in% numeric_feats) |>     
  left_join(X_test_long, by = c("sample_id", "feature")) |>
  group_by(feature) |>
  mutate(mean_abs_shap = mean(abs(shap_value))) |>
  ungroup() |>
  mutate(feature = reorder(feature, mean_abs_shap)) 

shap_plot_data %>% filter(feature == "PASI")
ggplot(shap_plot_data, aes(x = feature_value, y = shap_value, color = shap_value)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_hline(yintercept = 0, color = "gray70", linetype = "dashed", size = 0.8) +
  scale_color_gradientn(colors = c("#0000FF", "#00FFFF", "#00FF80", "#00FF00")) +
  scale_y_continuous(breaks = seq(-2, 2, by = 0.5)) +
  scale_x_continuous(breaks = seq(0, 70, by = 10)) +
  labs(
    title = "SHAP value vs. PASI",
    x = "PASI",
    y = "SHAP value for PASI"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(family = "Calisto", hjust = 0.5, face = "bold", size = 130),
    axis.text = element_text(family = "Calisto", color = "black", size = 120),
    axis.title.x = element_text(family = "Calisto", color = "black", size = 120),
    axis.title.y = element_text(family = "Calisto", color = "black", size = 120),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid = element_blank()
  )

graph <- "~/Downloads"
ggsave(
  filename = file.path(graph, "PASI.png"), 
  plot = last_plot(),                  
  width = 12,                                     
  height = 8,                                     
  dpi = 500                                  
)



#-----------------------------------------------------------------------------------------------#
#------------------------------ Shap value with all dataset-------------------------------------
#-----------------------------------------------------------------------------------------------#
Random_forest_SHAP_full <- function(X_select, X_features) {
  set.seed(845)
  all_data <- data.frame(X_select, DLQI = y)
  
  final_model <- ranger(
    formula = DLQI ~ .,
    data = all_data,
    mtry = 3,
    splitrule = "extratrees",
    min.node.size = 10,
    num.trees = 190
  )
  
  pred_all <- predict(final_model, data = all_data)$predictions
  
  RMSE_val        <- RMSE(all_data$DLQI, pred_all)
  RMSE_custom_val <- RMSE_custom(all_data$DLQI, pred_all, X_features)
  MAPE_val        <- MAPE(all_data$DLQI, pred_all)
  R_squared_val   <- R_square(all_data$DLQI, pred_all)
  adjR_squared    <- adjusted_r_squared(R_squared_val, pred_all, X_features)
  
  cat("Internal metrics (บนข้อมูลเต็ม):\n",
      "  RMSE          :", RMSE_val, "\n",
      "  RMSE Custom   :", RMSE_custom_val, "\n",
      "  MAPE          :", MAPE_val, "\n",
      "  R-squared     :", R_squared_val, "\n",
      "  Adj. R-squared:", adjR_squared, "\n")
  ### ------------------------------------------------------------------
  ### SHAP
  X_all <- all_data |> select(-DLQI)
  
  predictor <- Predictor$new(
    model = final_model,
    data  = X_all,
    y     = all_data$DLQI
  )
  
  shap_list    <- vector("list", nrow(X_all))
  shap_df_list <- vector("list", nrow(X_all))
  
  for (i in seq_len(nrow(X_all))) {
    shap_i          <- Shapley$new(predictor, x.interest = X_all[i, , drop = FALSE])
    shap_list[[i]]  <- shap_i
    shap_vec        <- setNames(shap_i$results$phi, shap_i$results$feature)
    shap_df_list[[i]] <- as.data.frame(t(shap_vec))
  }
  shap_df <- bind_rows(shap_df_list) |> `row.names<-`(paste0("Sample_", seq_len(nrow(X_all))))
  
  list(
    final_model  = final_model,
    RMSE         = RMSE_val,
    RMSE_custom  = RMSE_custom_val,
    MAPE         = MAPE_val,
    R_squared    = R_squared_val,
    adjR_squared = adjR_squared,
    predictions  = pred_all,
    shap_values  = shap_list,
    shap_df      = shap_df
  )
}

res <- Random_forest_SHAP_full(EN_select, X_EN_features)
dim(res$shap_df)

# Bar chart (mean |SHAP value|)
shap_mean_df <- data.frame(
  feature = colnames(res$shap_df),
  mean_abs_shap = apply(abs(res$shap_df), 2, mean)
)

shap_mean_df <- shap_mean_df %>%
  arrange(desc(mean_abs_shap)) %>%
  mutate(
    rank = row_number(),
    color = ifelse(rank <= 4, "#B28DFF", "#99CCED")  # green for top 4, bluegrey otherwise
  )

shap_mean_df <- shap_mean_df %>%
  mutate(feature = recode(feature,
                          "Stress6" = "STR06",
                          "Stress8" = "STR08",
                          "Stress11" = "STR11",
                          "other.disease2" = "COM2",
                          "Gender2" = "GEN2",
                          "Age" = "AGE",
                          "PASI" = "PASI",
                          "Stress3" = "STR03"
  ))

xmin <- 0.5
xmax <- nrow(shap_mean_df) + 0.5
ymin <- 0
ymax <- max(shap_mean_df$mean_abs_shap) * 1.15

# Create Bar Chart
p1 <- ggplot(shap_mean_df, aes(x = reorder(feature, mean_abs_shap), y = mean_abs_shap, fill = color)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = sprintf("%.3f", mean_abs_shap)),
            hjust = -0.1, size = 9.5, family = "Calisto", fontface = "bold") +
  scale_fill_identity() +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Feature",
    y = "Mean |SHAP|"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Calisto", face = "bold", size = 28, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 28, color = "black"),
    panel.border = element_blank(),
    axis.title.x = element_text(color = "black", size = 28),
    axis.title.y = element_text(color = "black", size = 28),
    axis.text = element_text(color = "black", size = 28)
  ) +
  annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
           color = "black", fill = NA, size = 1)

# Prepare data for Beeswarm plot
shap_long <- res$shap_df %>%
  mutate(sample_id = row_number()) %>%
  pivot_longer(-sample_id, names_to = "feature", values_to = "shap_value")

X_test <- res$shap_values[[1]]$x.interest
for (i in 2:length(res$shap_values)) {
  X_test <- rbind(X_test, res$shap_values[[i]]$x.interest)
}

X_test_long <- X_test %>%
  mutate(sample_id = row_number()) %>%
  pivot_longer(-sample_id, names_to = "feature", values_to = "feature_value")

shap_plot_data <- shap_long %>%
  left_join(X_test_long, by = c("sample_id", "feature"))

shap_plot_data <- shap_plot_data %>%
  mutate(feature = recode(feature,
                          "Stress6" = "STR06",
                          "Stress8" = "STR08",
                          "Stress11" = "STR11",
                          "other.disease2" = "COM2",
                          "Gender2" = "GEN2",
                          "Age" = "AGE",
                          "PASI" = "PASI",
                          "Stress3" = "STR03"
  ))

shap_plot_data <- shap_plot_data %>%
  group_by(feature) %>%
  mutate(mean_abs_shap = mean(abs(shap_value))) %>%
  ungroup() %>%
  mutate(feature = reorder(feature, mean_abs_shap))

# Create Beeswarm Plot
p2 <- ggplot(shap_plot_data, aes(x = shap_value, y = feature, color = feature_value)) +
  geom_blank() +
  geom_jitter(height = 0.25, size = 3, alpha = 0.7) +
  geom_vline(xintercept = 0, color = "grey40", size = 0.7) +
  scale_color_gradientn(
    colors = c("#C4FEFF", "#99CCED", "#9C89FF", "#7E38B7", "#541675"),
    name = "Feature value"
  ) +
  scale_x_continuous(breaks = seq(-5, 5, by = 1)) +
  labs(
    x = "SHAP Value (Impact on Model Output)",
    y = NULL
  ) +
  theme_minimal() +
  guides(title = NULL,
         color = guide_colorbar(barheight = unit(15, "cm"),
                                barwidth = unit(0.6, "cm"),
                                title.hjust = 0.5,
                                title.vjust = 5)) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
    panel.grid.minor.y = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    plot.title = element_text(family = "Calisto", face = "bold", color = "black", size = 28, hjust = 0.5),
    axis.text = element_text(family = "Calisto", color = "black", size = 28),
    axis.title.x = element_text(family = "Calisto", color = "black", size = 28),
    legend.text = element_text(family = "Calisto", color = "black", size = 28),
    legend.title = element_text(family = "Calisto", color = "black", size = 28, vjust = 0),
    legend.spacing = unit(1, "cm")
  )

# Combine plots side by side using patchwork
combined_plot <- p1 + plot_spacer() + p2 + plot_layout(widths = c(1, 0.05, 1), ncol = 3)
graph <- "~/Documents/SHAP"
ggsave(
  filename = file.path(graph, "Combined_SHAP_plots2.pdf"),
  plot = combined_plot,
  width = 24,
  height = 8,
  device = "pdf"
)

color = c("#99E6C0", "#66DDAA", "#40CC89", "#40BBAA", "#3399CC")
color = c("#C4FEFF", "#99CCED", "#9C89FF", "#7E38B7", "#541675")
color = c("#797EF6", "#4ADEDE", "#1AA7EC", "#1E2F97")
color = c("#85E3FF", "#6EB5FF", "#B28DFF")
color = c("#61F4DE", "#65CBE9", "#6C8DFA")
color = c('#aec7e8', '#6baed6', '#1f77b4')
color = c("#C4FAF8", "#ACE7FF", "#85E3FF", "#B5B9FF", "#A79AFF")
color = c("#ACE7FF","#85E3FF", "#97A2FF", "#B28DFF")

color = c("#C4FEFF", "#99CCED", "#9C89FF", "#7E38B7", "#541675")
color = c("#ACE7FF","#85E3FF", "#97A2FF", "#B28DFF")


# PASI and Age
shap_cont <- function (feature_name, name, color) {
  set.seed(845)         
  y_actual <- data$DLQI
  X_actual <- data[, !colnames(data) %in% "DLQI"]
  
  shap_long <- res$shap_df |>
    mutate(sample_id = row_number()) |>
    pivot_longer(-sample_id,
                 names_to  = "feature",
                 values_to = "shap_value")
  
  numeric_feats <- c("PASI", "Age")
  
  X_test_long <- X_actual |>
    mutate(sample_id = row_number()) |>
    select(sample_id, all_of(numeric_feats)) |>
    pivot_longer(-sample_id,
                 names_to  = "feature",
                 values_to = "feature_value")
  
  shap_plot_data <- shap_long |>
    filter(feature %in% numeric_feats) |>       
    left_join(X_test_long, by = c("sample_id", "feature")) |>
    group_by(feature) |>
    mutate(mean_abs_shap = mean(abs(shap_value))) |>
    ungroup() |>
    mutate(feature = reorder(feature, mean_abs_shap))   
  
  shap_plot_data <- shap_plot_data |>
    left_join(tibble(sample_id = 1:length(y_actual), y_actual = y_actual), by = "sample_id")
  
  ggplot(shap_plot_data %>% filter(feature == feature_name), aes(x = feature_value, y = shap_value, color = y_actual)) +
    geom_point(alpha = 0.7, size = 3) +
    geom_hline(yintercept = 0, color = "gray70", linetype = "dashed", size = 0.8) +
    scale_color_gradientn(
      colors = color, 
      name = "DLQI") +
    scale_y_continuous(breaks = seq(-3, 3, by = 1)) +
    scale_x_continuous(breaks = seq(0, 80, by = 10)) +
    labs(
      x = name,
      y = "SHAP value"
    ) +
    theme_minimal() +
    guides(title = NULL,
           color = guide_colorbar(barheight = unit(8, "cm"),
                                  barwidth = unit(0.4, "cm"),
                                  title.hjust = 0.5,
                                  title.vjust = 5)) +   
    theme(
      plot.title = element_text(family = "Calisto", hjust = 0.5, face = "bold", size = 23),
      axis.text = element_text(family = "Calisto", color = "black", size = 23),
      axis.title.x = element_text(family = "Calisto", color = "black", size = 23),
      axis.title.y = element_text(family = "Calisto", color = "black", size = 23),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid = element_blank(),          
      legend.key.width  = unit(0.5, "cm"),     
      legend.key.height = unit(1.5, "cm"),     
      legend.spacing.y  = unit(2, "cm"),       
      legend.text = element_text(family = "Calisto", color = "black", size = 18),
      legend.title = element_text(family = "Calisto", color = "black", size = 18)
    )
}

# Stress 3, 6, 8, 11
shap_stress <- function (feature, name, color) {
  set.seed(845)   
  y_actual <- data$DLQI
  X_actual <- data[, !colnames(data) %in% "DLQI"]
  feature_to_plot <- feature
  ordinal_feats <- c("Stress3", "Stress6", "Stress8", "Stress11")   
  
  shap_long <- res$shap_df |>
    mutate(sample_id = row_number()) |>
    pivot_longer(-sample_id, names_to = "feature", values_to = "shap_value")
  
  X_test_long <- X_actual |>
    mutate(sample_id = row_number()) |>
    select(sample_id, all_of(ordinal_feats)) |>
    pivot_longer(-sample_id, names_to = "feature", values_to = "feature_value")
  
  shap_plot_data <- shap_long |>
    filter(feature %in% ordinal_feats) |>
    left_join(X_test_long, by = c("sample_id", "feature")) |>
    group_by(feature) |>
    mutate(mean_abs_shap = mean(abs(shap_value))) |>
    ungroup() |>
    mutate(feature = reorder(feature, mean_abs_shap)) |>
    left_join(tibble(sample_id = 1:length(y_actual), y_actual = y_actual), by = "sample_id")
  
  shap_plot_feature <- shap_plot_data %>%
    filter(feature == feature_to_plot) %>%
    mutate(feature_value = factor(feature_value, levels = 0:3))
  
  mean_shap <- shap_plot_feature %>%
    group_by(feature_value) %>%
    summarise(mean_shap = mean(shap_value, na.rm = TRUE))
  
  ggplot() +
    geom_col(data = mean_shap,
             aes(x = feature_value, y = mean_shap),
             fill = "gray90", width = 0.4) +
    geom_jitter(data = shap_plot_feature,
                aes(x = feature_value, y = shap_value, color = y_actual),
                width = 0.1, alpha = 0.7, size = 3) +
    geom_hline(yintercept = 0, color = "gray70", linetype = "dashed", size = 0.8) +
    scale_color_gradientn(
      colors = color,
      name = "DLQI") +
    labs(
      x = name,
      y = "SHAP value"
    ) +
    theme_minimal() +   
    guides(title = NULL,
           color = guide_colorbar(barheight = unit(8, "cm"),
                                  barwidth = unit(0.4, "cm"),
                                  title.hjust = 0.5,
                                  title.vjust = 5)) +   
    theme(
      plot.title = element_text(family = "Calisto", hjust = 0.5, face = "bold", size = 23),
      axis.text = element_text(family = "Calisto", color = "black", size = 23),
      axis.title.x = element_text(family = "Calisto", color = "black", size = 23),
      axis.title.y = element_text(family = "Calisto", color = "black", size = 23),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid = element_blank(),
      legend.key.width  = unit(0.5, "cm"),
      legend.key.height = unit(1.5, "cm"),
      legend.spacing.y  = unit(2, "cm"),
      legend.text = element_text(family = "Calisto", color = "black", size = 18),
      legend.title = element_text(family = "Calisto", color = "black", size = 18)
    )
}

# Gender2 and other.disease2
shap_cat <- function (feature, name, color) {
  set.seed(845)
  complete_data <- data.frame(X_trans_num, X_trans_cat, DLQI = y)
  y <- complete_data$DLQI
  X <- complete_data[, !colnames(complete_data) %in% "DLQI"]
  
  feature_to_plot <- feature
  ordinal_feats <- c("Gender2", "other.disease2")   
  
  shap_long <- res$shap_df |>
    mutate(sample_id = row_number()) |>
    pivot_longer(-sample_id, names_to = "feature", values_to = "shap_value")
  
  X_test_long <- X |>
    mutate(sample_id = row_number()) |>
    select(sample_id, all_of(ordinal_feats)) |>
    pivot_longer(-sample_id, names_to = "feature", values_to = "feature_value")
  
  shap_plot_data <- shap_long |>
    filter(feature %in% ordinal_feats) |>
    left_join(X_test_long, by = c("sample_id", "feature")) |>
    group_by(feature) |>
    mutate(mean_abs_shap = mean(abs(shap_value))) |>
    ungroup() |>
    mutate(feature = reorder(feature, mean_abs_shap)) |>
    left_join(tibble(sample_id = 1:length(y), y = y), by = "sample_id")
  
  shap_plot_feature <- shap_plot_data %>%
    filter(feature == feature_to_plot) %>%
    mutate(feature_value = factor(feature_value, levels = 0:3))
  
  mean_shap <- shap_plot_feature %>%
    group_by(feature_value) %>%
    summarise(mean_shap = mean(shap_value, na.rm = TRUE))
  
  ggplot() +
    geom_col(data = mean_shap,
             aes(x = feature_value, y = mean_shap),
             fill = "gray90", width = 0.4) +
    geom_jitter(data = shap_plot_feature,
                aes(x = feature_value, y = shap_value, color = y),
                width = 0.1, alpha = 0.7, size = 3) +
    geom_hline(yintercept = 0, color = "gray70", linetype = "dashed", size = 0.8) +
    scale_color_gradientn(
      colors = color,
      name = "DLQI") +
    labs(
      x = name,
      y = "SHAP value"
    ) +
    theme_minimal() +   
    guides(title = NULL,
           color = guide_colorbar(barheight = unit(8, "cm"),
                                  barwidth = unit(0.4, "cm"),
                                  title.hjust = 0.5,
                                  title.vjust = 5)) +   
    theme(
      plot.title = element_text(family = "Calisto", hjust = 0.5, face = "bold", size = 23),
      axis.text = element_text(family = "Calisto", color = "black", size = 23),
      axis.title.x = element_text(family = "Calisto", color = "black", size = 23),
      axis.title.y = element_text(family = "Calisto", color = "black", size = 23),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid = element_blank(),
      legend.key.width  = unit(0.5, "cm"),
      legend.key.height = unit(1.5, "cm"),
      legend.spacing.y  = unit(2, "cm"),
      legend.text = element_text(family = "Calisto", color = "black", size = 18),
      legend.title = element_text(family = "Calisto", color = "black", size = 18)
    )
}

# Create all 8 plots
pasi <- shap_cont(feature_name = "PASI", name = "PASI", color)
age <- shap_cont(feature_name = "Age", name = "AGE", color)
stress3 <- shap_stress(feature = "Stress3", name = "STR03", color)
stress6 <- shap_stress(feature = "Stress6", name = "STR06", color)
stress8 <- shap_stress(feature = "Stress8", name = "STR08", color)
stress11 <- shap_stress(feature = "Stress11", name = "STR11", color)
gender <- shap_cat(feature = "Gender2", name = "GEN2", color)
other_dis <- shap_cat(feature = "other.disease2", name = "COM2", color)

# Combine all plots in 2 rows x 4 columns layout
design <- "
ABCD
EFGH
"
combined_all_plots <- stress6 + stress8 + stress11 + stress3 +
  age + pasi + other_dis + gender +
  plot_layout(design = design)

graph <- "~/Documents/SHAP"
ggsave(
  filename = file.path(graph, "All_SHAP_2x4_2.pdf"),
  plot = combined_all_plots,
  width = 32,
  height = 12, 
  device = "pdf"
)

# Combine all plots in 2 rows x 4 columns layout
design <- "
AB
CD
EF
GH
"
combined_all_plots <- stress6 + stress8 + stress11 + stress3 +
  age + pasi + other_dis + gender +
  plot_layout(design = design)

graph <- "~/Documents/SHAP"
ggsave(
  filename = file.path(graph, "All_SHAP_4x2_2.pdf"),
  plot = combined_all_plots,
  width = 16,
  height = 16, 
  device = "pdf"
)

