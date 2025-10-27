
sorting <- function(data, cohort, outcome1, outcome2){
  filter(data,
         Cohort == cohort,
         Reason.for.discontinuation == outcome1 |
           Reason.for.discontinuation == outcome2)
}

met_sorting <- function(data, outcome1, outcome2){
  filter(data,
         Reason.for.discontinuation == outcome1 |
           Reason.for.discontinuation == outcome2)
}


k_validate <- function(seed, neg.outcome, pos.outcome, outcome, tree){
  set.seed(seed)


  out <- list()

  out[[1]] <- getROC(train.seq, train.outcomes, test.seq, test.outcomes, neg.outcome, pos.outcome, outcome, tree)


  avg_AUCPR <- lapply(out, function(x){ x[[1]]})
  message("\nAccuracy average: ", (rowMeans(as.data.frame(avg_AUCPR))))
  return(out)
}


# Function create and test RF model #
getROC <- function(train.seq, train.outcomes, test.seq, test.outcomes, neg.outcome, pos.outcome, outcome, tree){
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

  y <- as.data.frame(aucpr$curve)
  print(ggplot(y, aes(V1, V2))+geom_path()+ylim(0,1))

  perf <- performance(pred, "tpr", "fpr")

  AUCPR_ROCR <- performance(pred, measure = "aucpr")
  print(AUCPR_ROCR)
  AUCPR <- AUCPR_ROCR@y.values[[1]]
  print(AUCPR)

  df <- data.frame(FalsePositive=c(perf@x.values[[1]]),
                   TruePositive=c(perf@y.values[[1]]))
  varimp <- model.training$importance
  out <- list(AUCPR, df, pred_cm, varimp)

  return(out)
}


# Grabbing AUROC values from input models #
grabVals <- function(output, input, seed_list){
  output <- list()

  for(i in 1:length(seed_list)){
    output <- append(output, input[[i]][[1]][[1]])
  }

  output <- do.call(rbind.data.frame, output)
  colnames(output) <- c("AUCPR")
  output
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
kTest <- function(seed_list, neg.outcome, pos.outcome, outcome, tree){
  out <- list()
  for(i in 1:length(seed_list)){
    out[[i]] <- k_validate(seed = seed_list[i], neg.outcome, pos.outcome, outcome, tree)
  }
  # out <- grabVals(out)
  # CW edit to function
  # out <- do.call(rbind.data.frame, out)
  # colnames(out) <- c("AUCPR")
  out
}


fishing <- function(output, data, seed_list){
  output <- list()
  for (i in 1:length(seed_list)) {
    output[[i]] <- fisher.test(data[[i]][[1]][[3]])
    output[[i]] <- output[[i]][[1]][[1]]
  }
  output <- do.call(rbind.data.frame, output)
  colnames(output) <- c("p.value")
  output

}


p.calc <- function(data, random){
  avg <- (data)
  p <- ecdf(random)
  p(avg)

}




add.metrics <- function(dataobject, output, pos.outcome){
  test_conmat2 <- data.frame(run = 1:25000,
                             Accuracy = NA,
                             Specificity = NA,
                             Precision = NA,
                             Recall = NA,
                             F1 = NA)

  for (i in 1:25000) {
    test_conmat <- confusionMatrix(dataobject[[i]][[1]][[3]],
                                   positive = pos.outcome)

    test_conmat2$Accuracy[i] <- test_conmat[[3]][1]
    test_conmat2$Specificity[i] <- test_conmat[[4]][2]
    test_conmat2$Precision[i] <- test_conmat[[4]][5]
    test_conmat2$Recall[i] <- test_conmat[[4]][6]
    test_conmat2$F1[i] <- test_conmat[[4]][7]


  }
  test_conmat2[is.na(test_conmat2)] <- 0
  x = 1
  output <- data.frame(run = 1:1000,
                       Avg.Accuracy = NA,
                       Avg.Specificity = NA,
                       Avg.Precision = NA,
                       Avg.Recall = NA,
                       Avg.F1 = NA)

  for (i in 1:1000) {
    y = x+24

    output$Avg.Accuracy[i] <- mean(test_conmat2$Accuracy[x:y])
    output$Avg.Specificity[i] <- mean(test_conmat2$Specificity[x:y])
    output$Avg.Precision[i] <- mean(test_conmat2$Precision[x:y])
    output$Avg.Recall[i] <- mean(test_conmat2$Recall[x:y])
    output$Avg.F1[i] <- mean(test_conmat2$F1[x:y])

    x = x+25
  }
  output
}


