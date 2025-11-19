


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


  #Add lasso for selection to the model


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


# Grabbing AUROC values from input models #
# grabVals <- function(output, input, seed_list){
#   output <- list()
#
#   for(i in 1:length(seed_list)){
#     output <- append(output, input[[i]][[1]][[1]])
#   }
#
#   output <- do.call(rbind.data.frame, output)
#   colnames(output) <- c("AUCPR")
#   output
# }
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

grabImp <- function(output, input, seed_list){
  datalist <- list()

  for(i in 1:length(seed_list)){
    tempdata <- data.frame(as.list(input[[i]][[1]][[4]][,4]))
    datalist[[i]] <- tempdata
  }
  output <- do.call(rbind, datalist)

  output
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


