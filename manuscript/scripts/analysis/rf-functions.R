
################################################################################
# Functions for running RF
#
# to run call kTest followed by grabVals to get AUCPR and values for ROC curves
#
################################################################################

# this function is an additional layer we aren't using right now
# could use this to add cross fold validation as it was originally intended for
k_validate <- function(seed, neg.outcome, pos.outcome, outcome, factors_list, validation, val_df, tree){
  set.seed(seed)

  out <- list()

  out[[1]] <- getROC(seed, neg.outcome, pos.outcome, outcome, factors_list, validation, val_df, tree)

  avg_AUCPR <- lapply(out, function(x){ x[[1]]})
  message("\nAUCPR average: ", (rowMeans(as.data.frame(avg_AUCPR))))
  return(out)
}

# Function create and test RF model #
getROC <- function(seed, neg.outcome, pos.outcome, outcome, factors_list, validation, val_df, tree){

  set.seed(seed)

  if(validation == TRUE) {

    train.outcomes <- demographics %>%
      dplyr::select(`Patient Id`, irAE, Response) %>%
      subset(`Patient Id` %in% baseline$`Patient Id`) %>%
      arrange(desc(`Patient Id`))

    test.outcomes <- val_df %>%
      select(`Patient Id`, outcome) %>%
      arrange(desc(`Patient Id`))

    train.seq <- base.micro %>%
      select(factors_list) %>%
      arrange(desc(`Patient Id`))

    test.seq <- val_df %>%
      select(factors_list) %>%
      arrange(desc(`Patient Id`))

  } else {

    #splitting data into train and test
    split <- sample.split(demographics[[outcome]], SplitRatio = 0.5)

    train <- subset(demographics, split == TRUE) %>%
      subset(`Patient Id` %in% blood$`Patient Id`) %>%
      subset(`Patient Id` %in% baseline$`Patient Id`)

    test <- subset(demographics, split == FALSE) %>%
      subset(`Patient Id` %in% blood$`Patient Id`) %>%
      subset(`Patient Id` %in% baseline$`Patient Id`)

    train.outcomes <- demographics %>%
      dplyr::select(`Patient Id`, irAE, Response) %>%
      subset(`Patient Id` %in% train$`Patient Id`) %>%
      subset(`Patient Id` %in% blood$`Patient Id`) %>%
      subset(`Patient Id` %in% baseline$`Patient Id`)

    train.outcomes <- arrange(train.outcomes, desc(`Patient Id`))

    test.outcomes <- demographics %>%
      dplyr::select(`Patient Id`, irAE, Response) %>%
      subset(`Patient Id` %in% test$`Patient Id`) %>%
      subset(`Patient Id` %in% blood$`Patient Id`) %>%
      subset(`Patient Id` %in% baseline$`Patient Id`)

    test.outcomes <- arrange(test.outcomes, desc(`Patient Id`))

    #creating objects for the model to be trained and tested on
    factors <- all_factors %>%
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

# Grabbing the variables used in each model
# Function to get distinct rownames from all list elements
get_vars <- function(myList) {
  unique(unlist(lapply(myList, function(x) rownames(x[[1]][[4]]))))
}

# get vars from top AUC
get_best_auc_vars <- function(myList) {
  aucs <- sapply(myList, function(x) x[[1]][[1]])
  best_idx <- which.max(aucs)
  rownames(myList[[best_idx]][[1]][[4]])
}

# get vars from AUC > 0.8
get_vars_auc_threshold <- function(myList, threshold = 0.8) {
  aucs <- sapply(myList, function(x) x[[1]][[1]])

  # which models pass the threshold?
  idx <- which(aucs > threshold)
  if (length(idx) == 0) return(character(0))
  vars <- unlist(lapply(idx, function(i) rownames(myList[[i]][[1]][[4]])))
  unique(vars)
}

# Main function call to generate RF models using 25 seeds #
kTest <- function(seed_list, neg.outcome, pos.outcome, outcome, factors_list, validation, val_df, tree){
  out <- list()
  for(i in 1:length(seed_list)){
    out[[i]] <- k_validate(seed = seed_list[i], neg.outcome, pos.outcome, outcome, factors_list, validation, val_df, tree)
  }
  # out <- grabVals(out)
  # CW edit to function
  # out <- do.call(rbind.data.frame, out)
  # colnames(out) <- c("AUCPR")
  out
}

grabImp <- function(input, seed_list){

  out <- list()

  for (i in seq_along(seed_list)) {

    # Extract getROC() output for this seed
    res <- input[[i]][[1]]

    impdf <- res[[4]] %>%   # Variable importance export
      as.data.frame() %>%
      rownames_to_column(var = "Factors")  #convert row names to column

    impdf$seed  <- seed_list[i]

    out[[i]] <- impdf
  }

  # bind into one big data frame
  dplyr::bind_rows(out)
}
