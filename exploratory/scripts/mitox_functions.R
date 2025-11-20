
################################################################################
# Functions for running RF
#
# to run call kTest followed by grabVals to get AUCPR and values for ROC curves
#
################################################################################

k_validate <- function(seed, neg.outcome, pos.outcome, outcome, factors_list, tree){
  set.seed(seed)


  out <- list()

  out[[1]] <- getROC(seed, neg.outcome, pos.outcome, outcome, factors_list, tree)


  avg_AUCPR <- lapply(out, function(x){ x[[1]]})
  message("\nAUCPR average: ", (rowMeans(as.data.frame(avg_AUCPR))))
  return(out)
}


# Function create and test RF model #
getROC <- function(seed, neg.outcome, pos.outcome, outcome, factors_list, tree){

  set.seed(seed)

  #splitting data into train and test
  split <- sample.split(demographics[[outcome]], SplitRatio = 0.5)

  train  <- subset(demographics, split == TRUE) %>%
    subset(`Patient Id` %in% blood$`Patient Id`) %>%
    subset(`Patient Id` %in% baseline$`Patient Id`)

  test   <- subset(demographics, split == FALSE) %>%
    subset(`Patient Id` %in% blood$`Patient Id`) %>%
    subset(`Patient Id` %in% baseline$`Patient Id`)

  train.outcomes <- demographics %>%
    select(`Patient Id`, irAE, Response) %>%
    subset(`Patient Id` %in% train$`Patient Id`) %>%
    subset(`Patient Id` %in% blood$`Patient Id`) %>%
    subset(`Patient Id` %in% baseline$`Patient Id`)

  train.outcomes <- arrange(train.outcomes, desc(`Patient Id`))


  test.outcomes <- demographics %>%
    select(`Patient Id`, irAE, Response) %>%
    subset(`Patient Id` %in% test$`Patient Id`) %>%
    subset(`Patient Id` %in% blood$`Patient Id`) %>%
    subset(`Patient Id` %in% baseline$`Patient Id`)

  test.outcomes <- arrange(test.outcomes, desc(`Patient Id`))

  #creating objects for the model to be trained and tested on
  factors <- all_factors %>%
    select((factors_list))

  train.seq <- filter(factors,
                      (`Patient Id` %in% train.outcomes$`Patient Id`))

  train.seq <- arrange(train.seq, desc(`Patient Id`))

  test.seq <- filter(factors,
                     (`Patient Id` %in% test.outcomes$`Patient Id`))
  test.seq <- arrange(test.seq, desc(`Patient Id`))


#################### Add lasso for selection to the model ######################
  # Design matrix for LASSO (drop Patient Id)
#   x_train <- as.matrix(train.seq[ , -1, drop = FALSE ])
#   y_train <- as.factor(train.outcomes[[outcome]])
#
#   # glmnet requires numeric 0/1
#   y_numeric <- ifelse(y_train == pos.outcome, 1, 0)
#
#   # Fit cross-validated LASSO
#   lasso_fit <- cv.glmnet(
#     x_train,
#     y_numeric,
#     family = "binomial",
#     alpha = 1
#   )
#
# # Warning in lognet(x, is.sparse, y, weights, offset, alpha, nobs, nvars,  :
# #  one multinomial or binomial class has fewer than 8  observations; dangerous ground
#
# # This is due to class imbalance and is concerning
#
#   # Extract variables with non-zero coefficients
#   coef_mat <- coef(lasso_fit, s = "lambda.min")
#   selected_vars <- rownames(coef_mat)[coef_mat[,1] != 0]
#   selected_vars <- selected_vars[selected_vars != "(Intercept)"]
#
#   if (length(selected_vars) == 0) {
#     stop("LASSO selected zero variables — cannot run Random Forest.")
#   }
#
#   message("\nLASSO selected ", length(selected_vars), " variables.")
#
#   # Restrict train/test to selected predictors only
#   train.seq <- train.seq[ , c("Patient Id", selected_vars), drop = FALSE ]
#   test.seq  <- test.seq [ , c("Patient Id", selected_vars), drop = FALSE ]
################################################################################

  if(all(train.outcomes$`Patient Id` == train.outcomes$`Patient Id`) == FALSE){
    stop("Training Sample_IDs do not match")
  }
  if(all(test.outcomes$`Patient Id` == test.outcomes$`Patient Id`) == FALSE){
    stop("Testing Sample_IDs do not match")
  }
  model.training <- randomForest(x = train.seq[,-1, drop=FALSE], y =
                          as.factor(train.outcomes[[outcome]]), ntree = tree, importance = TRUE) #, mtry = 9

  test.preds <- predict(model.training,
                        test.seq[,-1, drop=FALSE])

  # Prediction confusion matrix
  pred_cm <- table(observed = test.outcomes[[outcome]],
                   predicted = test.preds)

  print(pred_cm)

  prediction_for_roc_curve <- predict(model.training,
                                      test.seq[,-1, drop=FALSE],
                                      type="prob")

  pred <- prediction(prediction_for_roc_curve[,2],
                     test.outcomes[[outcome]],
                     label.ordering =
                       c(neg.outcome, pos.outcome))
  #for making aucpr curves
  scores <- data.frame(prediction_for_roc_curve[,2], test.outcomes[[outcome]])

  scores <- scores %>%
    mutate(score = case_when((test.outcomes..outcome.. == pos.outcome) ~ 1,
                             TRUE ~ 0))
  # print(prediction_for_roc_curve[,2])
  aucpr <- pr.curve(scores.class0=scores[scores$score=="0",]$`prediction_for_roc_curve...2.`,
                    scores.class1=scores[scores$score=="1",]$`prediction_for_roc_curve...2.`,
                    curve=T)

  aucpr.plot <- as.data.frame(aucpr$curve)
  print(ggplot(aucpr.plot, aes(V1, V2))+geom_path()+ylim(0,1))

  perf <- performance(pred, "tpr", "fpr")

  AUCPR_ROCR <- performance(pred, measure = "aucpr")
  print(AUCPR_ROCR)
  AUCPR <- AUCPR_ROCR@y.values[[1]]
  print(AUCPR)

  df <- data.frame(FalsePositive=c(perf@x.values[[1]]),
                   TruePositive=c(perf@y.values[[1]]))
  varimp <- model.training$importance
  out <- list(AUCPR, df, pred_cm, varimp, aucpr.plot)

  return(out)
}

# Currently unused but var importance
# grabImp <- function(output, input, seed_list){
#   datalist <- list()
#
#   for(i in 1:length(seed_list)){
#     tempdata <- data.frame(as.list(input[[i]][[1]][[4]][,4]))
#     datalist[[i]] <- tempdata
#   }
#   output <- do.call(rbind, datalist)
#
#   output
# }

# Grabbing AUROC values from input models #
grabVals <- function(input, seed_list){

  out <- list()

  for (i in seq_along(seed_list)) {

    # Extract getROC() output for this seed
    res <- input[[i]][[1]]

    AUCPR <- res[[1]]    # AUPCR value
    rocdf <- res[[2]]    # ROC curve dataframe

    rocdf$seed  <- seed_list[i]
    rocdf$AUCPR <- AUCPR

    out[[i]] <- rocdf
  }

  # bind into one big data frame
  dplyr::bind_rows(out)
}

# Main function call to generate RF models using 25 seeds #
kTest <- function(seed_list, neg.outcome, pos.outcome, outcome, factors_list, tree){
  out <- list()
  for(i in 1:length(seed_list)){
    out[[i]] <- k_validate(seed = seed_list[i], neg.outcome, pos.outcome, outcome, factors_list, tree)
  }
  # out <- grabVals(out)
  # CW edit to function
  # out <- do.call(rbind.data.frame, out)
  # colnames(out) <- c("AUCPR")
  out
}
