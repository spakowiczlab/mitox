
################################################################################
# Functions for running RF
#
# to run call kTest followed by grabVals to get AUCPR and values for ROC curves
#
################################################################################

# this function is an additional layer we aren't using right now
# could use this to add cross fold validation as it was originally intended for
k_validate <- function(seed, df, neg.outcome, pos.outcome, outcome, factors_list, validation, val_df, tree){
  set.seed(seed)

  out <- list()

  out[[1]] <- getROC(seed, df, neg.outcome, pos.outcome, outcome, factors_list, validation, val_df, tree)

  avg_AUCPR <- lapply(out, function(x){ x[[3]]})
  message("\nAUCPR average: ", (rowMeans(as.data.frame(avg_AUCPR))))
  return(out)
}

# Function create and test RF model #
getROC <- function(seed, df, neg.outcome, pos.outcome, outcome, factors_list, validation, val_df, tree){

  set.seed(seed)

  if(validation == TRUE) {

    train.outcomes <- df %>%
      dplyr::select(`Patient Id`, irAE, Response) %>%
      arrange(desc(`Patient Id`))

    test.outcomes <- val_df %>%
      dplyr::select(`Patient Id`, outcome) %>%
      arrange(desc(`Patient Id`))

    train.seq <- base.micro %>%
      dplyr::select(factors_list) %>%
      arrange(desc(`Patient Id`))

    test.seq <- val_df %>%
      dplyr::select(factors_list) %>%
      arrange(desc(`Patient Id`))

  } else {

    #splitting data into train and test
    split <- sample.split(df[[outcome]], SplitRatio = 0.5)

    train <- subset(df, split == TRUE)
    test <- subset(df, split == FALSE)

    train.outcomes <- df %>%
      dplyr::select(`Patient Id`, irAE, Response) %>%
      subset(`Patient Id` %in% train$`Patient Id`)

    train.outcomes <- arrange(train.outcomes, desc(`Patient Id`))

    test.outcomes <- df %>%
      dplyr::select(`Patient Id`, irAE, Response) %>%
      subset(`Patient Id` %in% test$`Patient Id`)

    test.outcomes <- arrange(test.outcomes, desc(`Patient Id`))

    #creating objects for the model to be trained and tested on
    factors <- df %>%
      dplyr::select((factors_list))

    train.seq <- filter(factors,
                        (`Patient Id` %in% train.outcomes$`Patient Id`))

    train.seq <- arrange(train.seq, desc(`Patient Id`))

    test.seq <- filter(factors,
                       (`Patient Id` %in% test.outcomes$`Patient Id`))
    test.seq <- arrange(test.seq, desc(`Patient Id`))

  }

  ###################### Boruta for selection to the model ######################
  set.seed(seed)

  # Drop Patient Id for modeling
  x_train <- train.seq[ , -1, drop = FALSE ]
  y_train <- as.factor(train.outcomes[[outcome]])

  # Combine for Boruta formula interface
  boruta_data <- cbind(y = y_train, x_train)

  # Run Boruta
  set.seed(seed)
  bor <- Boruta(y ~ ., data = boruta_data, doTrace = 2, ntree = tree)

  # try with more trees
  # better for response but worse for irAE
  # bor <- Boruta(y ~ .,
  #               data = boruta_data,
  #               maxRuns = 500,
  #               doTrace = 1,
  #               ntree = max(tree, 1000))

  # Resolve tentative features
  bor_fixed <- TentativeRoughFix(bor)

  # Get confirmed features only
  selected_vars <- getSelectedAttributes(bor_fixed, withTentative = FALSE)

  min_vars <- 5   # set your desired minimum number of predictors

  if(length(selected_vars) < min_vars){
    warning("Boruta confirmed fewer than min_vars — including tentative features")
    set.seed(seed)
    selected_vars <- getSelectedAttributes(bor_fixed, withTentative = TRUE)
  }

  # still selecting 0 vars so need another fallback
  # rank based RF selection - get top 5
  if(length(selected_vars) < min_vars){
    warning("Boruta still selected fewer than min_vars — filling with top RF importance features")

    # Fit quick RF on all predictors
    set.seed(seed)
    rf_temp <- randomForest(x = x_train, y = y_train, ntree = max(500, tree))

    # Rank features by importance
    imp <- importance(rf_temp)
    top_vars <- rownames(imp)[order(-imp[,1])]

    # Add top features until min_vars is met
    selected_vars <- unique(c(selected_vars, top_vars[1:min_vars]))
  }

  message("\nBoruta selected ", length(selected_vars), " variables.")
  selected_vars <- selected_vars[!is.na(selected_vars)]

  # Restrict train/test to selected features
  train.seq <- train.seq[ , c("Patient Id", selected_vars), drop = FALSE ]
  test.seq  <- test.seq [ , c("Patient Id", selected_vars), drop = FALSE ]
  ###############################################################################

  if(all(train.outcomes$`Patient Id` == train.outcomes$`Patient Id`) == FALSE){
    stop("Training Sample_IDs do not match")
  }
  if(all(test.outcomes$`Patient Id` == test.outcomes$`Patient Id`) == FALSE){
    stop("Testing Sample_IDs do not match")
  }

  set.seed(seed)

  # Train random forest
  model.training <- randomForest(
    x = train.seq[,-1, drop = FALSE],
    y = as.factor(train.outcomes[[outcome]]),
    ntree = tree,
    importance = TRUE
  )

  # Predict class labels
  test.preds <- predict(model.training, test.seq[,-1, drop = FALSE])

  # Confusion matrix
  pred_cm <- table(
    observed = test.outcomes[[outcome]],
    predicted = test.preds
  )
  print(pred_cm)

  #  Predict probabilities
  pred_probs <- predict(model.training, test.seq[,-1, drop = FALSE], type = "prob")

  # ROC / AUROC (for plotting)
  pred_rocr <- prediction(pred_probs[,2], test.outcomes[[outcome]], label.ordering = c(neg.outcome, pos.outcome))

  # Precision-Recall / AUCPR using ROCR
  pr_perf <- performance(pred_rocr, "prec", "rec")

  # create fresh prediction object
  pred_rocr <- prediction(pred_probs[,2], test.outcomes[[outcome]], label.ordering = c(neg.outcome, pos.outcome))
  AUCPR_ROCR <- performance(pred_rocr, measure = "aucpr")
  auc_pr <- AUCPR_ROCR@y.values[[1]]

  pr_df <- data.frame(
    Recall    = unlist(pr_perf@x.values),
    Precision = unlist(pr_perf@y.values)
  )
  pr_df$seed <- seed

  # ROC curve
  pred_rocr <- prediction(pred_probs[,2], test.outcomes[[outcome]], label.ordering = c(neg.outcome, pos.outcome))
  roc_perf <- performance(pred_rocr, "tpr", "fpr")
  pred_rocr <- prediction(pred_probs[,2], test.outcomes[[outcome]], label.ordering = c(neg.outcome, pos.outcome))
  auc_roc  <- performance(pred_rocr, measure = "auc")@y.values[[1]]

  roc_df <- data.frame(
    FalsePositive = unlist(roc_perf@x.values),
    TruePositive  = unlist(roc_perf@y.values)
  )
  roc_df$seed <- seed

  # Variable importance
  varimp <- model.training$importance

  # Output
  out <- list(
    AUROC      = auc_roc,   # scalar
    ROC_df     = roc_df,    # for ROC plotting
    AUCPR      = auc_pr,    # scalar
    PR_df      = pr_df,     # for PR plotting
    ConfMat    = pred_cm,   # confusion matrix
    VarImp     = varimp     # variable importance
  )

  return(out)
}


# Main function call to generate RF models using 25 seeds #
kTest <- function(seed_list, df, neg.outcome, pos.outcome, outcome, factors_list, validation = FALSE, val_df = NULL, tree = 1000){
  out <- list()
  for(i in 1:length(seed_list)){
    out[[i]] <- k_validate(seed = seed_list[i], df, neg.outcome, pos.outcome, outcome, factors_list, validation, val_df, tree)
  }
  # out <- grabVals(out)
  out
}

grabVals <- function(input, seed_list) {

  # Lists to hold per-seed curve data
  roc_list <- list()
  pr_list  <- list()

  for (i in seq_along(seed_list)) {
    res <- input[[i]][[1]]  # getROC output for this seed

    # --- ROC data ---
    rocdf <- res[[2]]           # assuming ROC df is res[[2]]
    rocdf$AUCROC <- res[[1]]
    roc_list[[i]] <- rocdf

    # --- PR data ---
    prdf <- res[[4]]
    prdf$AUCPR <- res[[3]]
    pr_list[[i]] <- prdf
  }

  # Combine all seeds into one data.frame each
  roc_df <- dplyr::bind_rows(roc_list)
  pr_df  <- dplyr::bind_rows(pr_list)

  list(ROC = roc_df, PR = pr_df)
}

# get vars from AUC > 0.8
get_vars_auc_threshold <- function(myList, threshold = 0.8) {
  aucs <- sapply(myList, function(x) x[[1]][[1]])

  # which models pass the threshold?
  idx <- which(aucs > threshold)
  if (length(idx) == 0) return(character(0))
  vars <- unlist(lapply(idx, function(i) rownames(myList[[i]][[1]][[6]])))
  unique(vars)
}

# get variable importance
grabImp <- function(input, seed_list){

  out <- list()

  for (i in seq_along(seed_list)) {

    # Extract getROC() output for this seed
    res <- input[[i]][[1]]

    impdf <- res[[6]] %>%   # Variable importance export
      as.data.frame() %>%
      rownames_to_column(var = "Factors")  #convert row names to column

    impdf$seed  <- seed_list[i]

    out[[i]] <- impdf
  }

  # bind into one big data frame
  dplyr::bind_rows(out)
}
