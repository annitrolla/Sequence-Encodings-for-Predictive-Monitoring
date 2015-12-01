# CHECK WHETHER ALL PACKAGES ARE INSTALLED. IF NOT, INSTALL MISSING
required_packages <- c('data.table','Hmisc','HMM','reshape2','hmm.discnp','plyr','caret','ggplot2',
                       'gridExtra', 'ROCR','pROC','stringr','randomForest')
for(pkg in required_packages){
  if(pkg %in% rownames(installed.packages()) == FALSE) {
    install.packages(pkg)
  }
}
library(data.table)
library(Hmisc)
library(HMM)
library(reshape2)
library(hmm.discnp)
library(plyr)
library(caret)
library(ggplot2)
library(gridExtra)
library(pROC)
library(stringr)
library(randomForest)
library(e1071)
library(gbm)
library(rpart)

#library(dunn.test)


##--------------##
# GENERAL-PURPOSE FUNCTIONS
sink.reset <- function(){
  for(i in seq_len(sink.number())){
    sink(NULL)
  }
}
# for ploting multiple ggplots in the same window
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  plots <- c(list(...), plotlist)  
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))     
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# function for numeration of events in one sequence
seq_order_fun  <- function(x) {
  seq_along(x)
}

# numerating events in sequences in the dataset
seq_ordering <- function(dt, nr_of_events=c(1,2,3)) {
  sequences_all <- ddply(dt, .(sequence_nr), mutate, seq_order = seq_order_fun(sequence_nr))
  sequences_first_x <- subset(sequences_all, seq_order%in%nr_of_events)
  sequences_first_x$seq_order <- as.factor(sequences_first_x$seq_order)
  return(sequences_first_x)
}

# convertion of factors to numerical features
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

# function for merging particular features together for HMM
merge_fun <- function(x, y) {
  columns_x <- grep("sequence_nr|label|odds_rat.", colnames(x))
  columns_y <- grep("sequence_nr|label|odds_rat.", colnames(y))
  df <- merge(x[,columns_x],y[,columns_y], by=c("sequence_nr","label"),all=T)
  return(df)
}

##--------------##
# DATA PROCESSING COMMON FUNCTIONS

# input data pre-processing to extract features
splitting_features <- function(feature_pairs) {  
  feature_pairs <- gsub('}','', feature_pairs)
  pair_parts <- strsplit(feature_pairs,'=')
  feature_name <- unlist(lapply(pair_parts, function(x) feature_name=x[1]))
  feature_value <- unlist(lapply(pair_parts, function(x) feature_name=x[2]))
  features <- cbind.data.frame(feature_name, feature_value)  
  return(features)
}

#input data pre-processing to extract events
event_splitting <- function(events){
  events_list <- strsplit(events, '{', fixed=T)
  event_names <- lapply(events_list, function(x) x[1])
  event_nr <- seq_along(event_names)
  feature_pairs_list <- strsplit(unlist(lapply(events_list, function(x) x[2])), ',', fixed=T)
  features_list <- lapply(feature_pairs_list, splitting_features)   
  features <- mapply(cbind.data.frame, features_list, event_names, event_nr, SIMPLIFY=F)
  return(features)
}

# input data collection in one data.frame, with pre-specified names of columns
dt_collector <- function(dt){
  file <- readLines(dt)
  events_in_sequence <- strsplit(file,';')
  events_in_sequence <- lapply(events_in_sequence, function(x) gsub(', ', '', x))
  labels <- lapply(events_in_sequence, function(x) x[length(x)])
  sequence_nr <- seq_along(events_in_sequence)
  events_in_sequence_list <- lapply(events_in_sequence, function(x) x[-length(x)])
  events_preproc <- lapply(events_in_sequence_list, event_splitting)
  seq <- lapply(events_preproc, rbindlist)
  seq_labels <- mapply(cbind.data.frame, seq, labels, sequence_nr, SIMPLIFY=F)
  all_dt <- rbindlist(seq_labels)
  setnames(all_dt, c("attribute_name", 'attribute_value', 'activity_name', 'activity_nr','label','sequence_nr'))
  return(all_dt)
}

# ordering and deleting artifacts (duplicating values)
data_preprocess <- function(dt) {
  print('Reading the file and pre-processing it..')
  dat <- dt_collector(dt)
  dat <- as.data.frame(dat)
  dat_wo_duplicated <- dat[!duplicated(dat[, c('attribute_name', 'activity_name','activity_nr','label','sequence_nr')],fromLast=T),]
  dat <- dcast(dat_wo_duplicated, sequence_nr + activity_name + activity_nr + label ~ attribute_name, value.var='attribute_value', fill=NaN)
  dat_ordered <- dat[order(dat$sequence_nr, dat$activity_nr),]
  return(dat_ordered)
}

# renaming of colnames and values, removing space between them (R cannot handle this name convention)
rename_values <- function(data, columns_to_revalue=1:length(colnames(data))) {
  data_revalue <- apply(data[,columns_to_revalue], 2, function(x) str_replace_all(string=x, pattern=" ", repl="_"))
  data <- cbind.data.frame(data[,-columns_to_revalue], data_revalue)
  data <- as.data.frame(unclass(data))
  return(data)
}

rename_columns <- function(data) {
  colnames(data) <- str_replace_all(string=colnames(data), pattern=" ", repl="_")
  return(data)
}

# abbreviation of features if they are too long (R can handle up to some nr of chars)
abbrv_strings <- function(data, letters=TRUE, numbers=FALSE, feature='event_name') {
  if (letters == TRUE){
    name <- paste(feature,'_short',sep='')
    data[[name]] <- abbreviate(data[[feature]])
    data[[feature]] <- NULL
  }
  if (numbers == TRUE){
    name <- paste(feature,'_to_nr',sep='')
    unique_feature <- data.frame(feature_unique=unique(data[[feature]]))  
    data[[name]] <- match(data[[feature]], unique_feature$feature_unique)
  }  
  return(data)
}

# test/train division
train_test_division <- function(data, ratio=2/3, seed=sample.int(1000000,1)){
  unique_sequence_nr <- unique(data$sequence_nr)
  train_idx <- unique_sequence_nr[sample(length(unique_sequence_nr), ceiling(length(unique_sequence_nr)*ratio))]
  train_dt <- subset(data, sequence_nr%in%train_idx)
  test_dt <- subset(data, !(sequence_nr%in%train_idx))
  return(list(train_dt, test_dt))
}

##--------------##
# RANDOM FOREST FUNCTIONS
# splitting frequent_events to different factors (R cannot handle more than 33 levels on rf)
features_to_factorize <- function(data, factor_levels=20) {
  names(which(unlist(lapply(data, function(x) length(levels(x))>=factor_levels))==TRUE))
}

frequent_events_to_factors <- function(train_data=train_dt, test_data=test_dt, validation_data=validation_dt, feature, factor_levels=20){
  feature_freq <- as.data.frame(table(train_data[[feature]]))
  feature_freq <- feature_freq[order(feature_freq$Freq, decreasing = T),]
  total_nr_values <- nrow(feature_freq)
  pieces <- seq(1, total_nr_values, by=factor_levels)
  if (total_nr_values > pieces[length(pieces)]){
    pieces <- c(pieces, total_nr_values)
  }
  for (j in 1:(length(pieces)-1)){
    name <- paste('freq_', feature, '_', j, sep='')
    train_data[[name]] <- ifelse(train_data[[feature]] %in% feature_freq$Var1[pieces[j]:(pieces[j+1]-1)], as.character(train_data[[feature]]), 'other')
    test_data[[name]] <- ifelse(test_data[[feature]] %in% train_data[[name]], as.character(test_data[[feature]]), 'other')
    validation_data[[name]] <- ifelse(validation_data[[feature]] %in% train_data[[name]], as.character(validation_data[[feature]]), 'other')
  }
  return(list(train_data,validation_data,test_data))
}

# adjusting dataformat for digestible in modelling (factor levels unification between different sets of data)
data_to_model_with_selected_features <- function(list_data, features=pattern_to_grep){
  dt_bundle <- data.frame()
  for(i in 1:length(list_data)){
    k <- grep(features, colnames(list_data[[i]]))
    current_dt <- list_data[[i]][,k]
    dt <- as.data.frame(unclass(current_dt))
    dt$is_train <- i
    dt_bundle <- rbind.data.frame(dt_bundle,dt)
  }
  dt_bundle <- as.data.frame(unclass(dt_bundle))
  return(dt_bundle)
}

##--------------##
# HMM FUNCTIONS
# forward algorithm with handling of missing values (by default couldn't)
forward_modified <- function(hmm, observation) {
  hmm$transProbs[is.na(hmm$transProbs)] = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  nObservations = length(observation)
  nStates = length(hmm$States)
  f = array(NA, c(nStates, nObservations))
  dimnames(f) = list(states = hmm$States, index = 1:nObservations)
  for (state in hmm$States) {
    f[state, 1] = suppressWarnings(log(hmm$startProbs[state] * hmm$emissionProbs[state, observation[1]]))
  }
  for (k in 2:nObservations) {
    #print(paste("k:",k,sep=''))
    for (state in hmm$States) {
      #print(paste("state:",state,sep=''))
      logsum = -Inf
      for (previousState in hmm$States) {
        temp = f[previousState, k - 1] + log(hmm$transProbs[previousState, 
                                                            state])
        if (is.na(temp)==TRUE) {
          logsum = 0.0001 + log(1 + exp(logsum - 0.0001))
        }
        else if (temp > -Inf) {
          logsum = temp + log(1 + exp(logsum - temp))
        }
        #print(paste("temp:",temp,sep=''))
      }
      #print(paste("logsum:",logsum,sep=''))
      f[state, k] = log(hmm$emissionProbs[state, observation[k]]) + logsum
    }
  }
  return(f)
}

# creation of the lists for collecting HMM parameters
hmm_list = function(data) {
    lst = list()
    if(class(data)=="data.frame"){
      for (i in 1:nrow(data)) { lst[[i]] = data[i,1:ncol(data)] }}
    else{
      #for (i in 1:length(data)) { lst[[i]] = data[i,1] } 
      lst[[1]] <- data
    }
    return(lst)
}

# Log likelihood calculation for the particular model
log_lik_forward = function(data, model) {
  log_lik_vec = rep(NA, nrow(data)) 
  for (i in 1:nrow(data)){
    #print(i)
    fwd <- forward_modified(model, as.character(data[i,]))
    logforward = ifelse(is.finite(fwd)==TRUE, fwd, 0)
    log_lik_vec[i] = log(sum(exp(logforward[,ncol(data)]))) 
    } 
  return(log_lik_vec)
}

# difference in loglikelihoods of two HMM models
difference = function(data_input, model_1, model_2) {
  log_lik_1 = log_lik_forward(data_input,model_1)
  log_lik_2 = log_lik_forward(data_input,model_2)
  log_lik_1 <- ifelse(is.finite(log_lik_1)==FALSE, 0, log_lik_1)
  log_lik_2 <- ifelse(is.finite(log_lik_2)==FALSE, 0, log_lik_2)
  diff = log_lik_1 - log_lik_2  #>0 ->good, <0 bad
  return(list(diff, log_lik_1, log_lik_2))
}

# data modifications specifically for HMM
data_for_hmm_input <- function(data_to_process, type="true", nr_prefixes, event_data) {
  variation <- 'Varies'
  #if type == all, then no subsetting is done (wrt to label) and it can be used for purposes of testing HMM
  if (type == "all") {
    label_subset <- data_to_process   
  } else {
    label_subset  <- subset(data_to_process, type == label)
  }
  sequences_all_label <- label_subset[,c(event_data,"sequence_nr")]
  sequences_all <- ddply(sequences_all_label, .(sequence_nr), mutate, seq_order = seq_order_fun(sequence_nr))
  sequences_all[[event_data]] <- as.character(sequences_all[[event_data]])
  sequences_all <- dcast(sequences_all, sequence_nr~seq_order, value.var=event_data, fill="NA")  
  sequences_label_min <- sequences_all[,c(1:(nr_prefixes+1))] #choose number of prefixes in use
  #check if the sequence is static or not 
  if (length(table(as.vector(t(sequences_label_min[,!colnames(sequences_label_min) %in% 'sequence_nr']))))==1) {
    variation <- 'Static'
  }
  return(list(sequences_label_min, variation))
}

# creating the symbols for HMM emission matrices that match between true and false sequences
symbol_space_function <- function(hmm_true=result_main_true, hmm_false=result_main_false, test_data, opt="def", event_data) {
  #opt="rare" - option with rare event
  #symbols = c(1,2,3)
  ### The problem - we need full space of observations
  symbols_true  <-  attr(hmm_true$y, "uval")
  symbols_false  <-  attr(hmm_false$y, "uval")
  symbols_train  <- union(symbols_true, symbols_false)
  symbols_test <- unlist(unique(as.data.frame(test_data)[event_data]))
  #missing_symbols <- setdiff(test_symbols,symbols)
  symbols <- unique(c(union(symbols_train, symbols_test), "NA")) 
  if(opt=="rare"){
    return(symbols_train)
  }
  else {
    return(symbols)
  }
}

# estimation of emission matrices
emmision_calc <- function(nr_states, symbols, Rho) {
  emissions <- matrix(0.001, nrow=nr_states, ncol=length(symbols), dimnames=list(c(1:nr_states),symbols))
  for(j in 1:nr_states) {
    emissions[j,match(colnames(emissions), rownames(Rho),nomatch=0.001)] <- Rho[,j]
  }
  return(emissions)
}

# initialization of HMM model with fitted parameters
hmm_initializing = function(train_data, test_data, nr_prefixes, nr_states, seed, output="odds ratios", event_data, iter=200) {  
  lst_true_train <- data_for_hmm_input(data_to_process = train_data, type="true", nr_prefixes, event_data)
  lst_false_train <- data_for_hmm_input(data_to_process = train_data, type="false", nr_prefixes, event_data)
  lst_crossval <- data_for_hmm_input(data_to_process = test_data, type = "all", nr_prefixes, event_data)
  sequences_true_train <- lst_true_train[[1]]
  sequences_false_train <- lst_false_train[[1]]
  sequences_crossval <- lst_crossval[[1]]
  
  variation_true_train <- lst_true_train[[2]]
  variation_false_train <- lst_false_train[[2]]

  # in case  a feature is static, built-in package cannot build a model. Overcoming issue by passing
  # directly values to odds ratios
  if (variation_true_train == 'Static' & variation_false_train == 'Static') {
    sequences_crossval$odds_ratio <- 1
  } else if (variation_true_train == 'Static' & variation_false_train == 'Varies') {
      sequences_crossval$odds_ratio <- 0
  } else if (variation_true_train == 'Varies' & variation_false_train == 'Static') {
    sequences_crossval$odds_ratio <- 0
  } else {
    lst_main_true  <-  hmm_list(sequences_true_train[, !names(sequences_true_train)%in% 'sequence_nr'])
    lst_main_false  <-  hmm_list(sequences_false_train[, !names(sequences_false_train)%in% 'sequence_nr'])
    set.seed(seed)
    result_main_true  <-  try(hmm(lst_main_true, K = nr_states, itmax = iter),silent=T)
    result_main_false  <-  try(hmm(lst_main_false, K = nr_states, itmax = iter), silent=T)
    if ("try-error" %in% class(result_main_false)) {
      result_main_false  <-  hmm(lst_main_false, K = nr_states, mixture=T, itmax = iter)  
    } 
    #general initialization
    states  <-  paste("state", 1:nr_states, sep='')
    symbols <- symbol_space_function(hmm_true=result_main_true, hmm_false=result_main_false, test_data, opt="def", event_data)
    emissions_true <- emmision_calc(nr_states, symbols=symbols, Rho=result_main_true$Rho)
    emissions_false <- emmision_calc(nr_states, symbols=symbols, Rho=result_main_false$Rho)  
    initialize_hmm_true = initHMM(States = states, Symbols = symbols, 
                                startProbs = result_main_true$ispd,
                                transProbs=result_main_true$tpm,
                                emissionProbs = emissions_true)

    initialize_hmm_false = initHMM(States = states, Symbols = symbols, 
                                 startProbs = result_main_false$ispd,
                                 transProbs=result_main_false$tpm,
                                 emissionProbs = emissions_false)

  
     hmm_list_models <- list(initialize_hmm_true[[3]],initialize_hmm_true[[4]],#optional [[5]]
                          initialize_hmm_false[[3]],initialize_hmm_false[[4]])

     hmm_list_models_full <- list(initialize_hmm_true, initialize_hmm_false)
     diff_for_cross_val <- difference(sequences_crossval[,!names(sequences_crossval) %in% 'sequence_nr'], initialize_hmm_true, initialize_hmm_false)[[1]]
     sequences_crossval$odds_ratio <- diff_for_cross_val 
    }
    sequences_crossval$label  <- unique(test_data[,c("label","sequence_nr")])$label
    sequences_crossval$seq_check <- unique(test_data[,c("label","sequence_nr")])$sequence_nr
    #more detailed output can be extracted (used it for descriptive purposes)
    if (output == "odds ratios") {
       return(sequences_crossval)
    } else if(output == "hmm models") {
       return(hmm_list_models)
    } else if(output == "hmm full spec") {
       return(hmm_list_models_full)
    } else {
       print("Warning: output is not defined")
    }
}

# HMMs for one particular prefix
hmm_prefix_run <- function(train_data, test_data, nr_prefixes, features, nr_states, seed, iter){
  hmm_results_train <- list()
  hmm_results_test <- list()
  for (variable in features) {
    #print(variable)
    train_res <- hmm_initializing(train_data=train_data, test_data=train_data, nr_prefixes, nr_states, seed, output="odds ratios", event_data=variable, iter)
    test_res <- hmm_initializing(train_data=train_data, test_data=test_data, nr_prefixes, nr_states, seed, output="odds ratios", event_data=variable, iter)
    hmm_results_train[[variable]] <- train_res
    hmm_results_test[[variable]] <- test_res
  }
  dt_hmm_train <- suppressWarnings(Reduce(merge_fun, hmm_results_train))
  dt_hmm_test <- suppressWarnings(Reduce(merge_fun, hmm_results_test))
  or_names <- grep("odds_rat.",colnames(dt_hmm_train))
  colnames(dt_hmm_train)[or_names] <- paste("odds_ratio_", 1:length(or_names), sep='')
  colnames(dt_hmm_test)[or_names] <- paste("odds_ratio_", 1:length(or_names), sep='')
  return(list(dt_hmm_train, dt_hmm_test))
}

# Siloutte search for best nr of states according to diff parameters
best_nr_states <- function(train_data, validation_data, features, nr_prefixes, seed, measure_optimiz='AUC', states_range, iter=200) {  
  AUC_vs_states <- c()
  States <- c()
  for (st in states_range) {
    #print(st)
    dt_current <- hmm_prefix_run(train_data=train_data, test_data=validation_data, nr_prefixes=nr_prefixes, features=features, nr_states=st, seed=seed, iter)
    dt_basic_train_current <- dt_current[[1]]
    dt_basic_validation_current <- dt_current[[2]]
    model <- randomForest(data=dt_basic_train_current[,!names(dt_basic_train_current)%in% 'sequence_nr'], as.factor(label)~.)
    predictions_bin <- predict(model, dt_basic_validation_current[,!names(dt_basic_validation_current)%in% 'sequence_nr']) 
    predictions_p <- predict(model, dt_basic_validation_current[,!names(dt_basic_validation_current)%in% 'sequence_nr'], type='prob')
    roc_crv <- roc(dt_basic_validation_current$label, predictions_p[,'true'])
    auc <- roc_crv$auc[1]
    AUC_vs_states <- c(AUC_vs_states,auc)
    States <- c(States, st)
  }
  if (measure_optimiz == 'AUC') {
    z <- AUC_vs_states
  } else if (measure_optimiz == 'F1') {
    z <- F1_vs_states
  } else {
    z <- Accuracy_vs_states
  }
  optimal_state <- States[which.max(z)]
  return(c(max(z), optimal_state))
}

##----------##
# HMM + RF functions
# fitting best number of states by siloutte search using also random forest information into optimal search
best_nr_states_hmmrf <- function(train_data, validation_data, features, nr_prefixes, seed=seed_nr, 
                                 measure_optimiz='AUC', rf_train, rf_validation, states_range, iter) {  
  AUC_vs_states <- c()
  States <- c()
  for(st in states_range){
    #print(st)
    dt_current_hmm <- hmm_prefix_run(train_data=train_data, test_data=validation_data, nr_prefixes=nr_prefixes,
                                 features=features, nr_states=st, seed=seed, iter)
    dt_basic_train_current_hmm <- dt_current_hmm[[1]]
    dt_basic_validation_current_hmm <- dt_current_hmm[[2]]    
    
    dt_basic_train_current <- merge(rf_train, dt_basic_train_current_hmm, by=c('sequence_nr','label'), all=T)    
    dt_basic_validation_current <- merge(rf_validation, dt_basic_validation_current_hmm, by=c('sequence_nr','label'), all=T)
    
    dt_basic_train_current$is_train <- NULL
    dt_basic_validation_current$is_train <- NULL    
    
    model <- randomForest(data=dt_basic_train_current[,!names(dt_basic_train_current) %in% 'sequence_nr'], as.factor(label)~.)
    predictions_p <- predict(model, dt_basic_validation_current[,!names(dt_basic_validation_current)%in% 'sequence_nr'], type='prob')
    predictions_bin <- predict(model, dt_basic_validation_current[,!names(dt_basic_validation_current)%in% 'sequence_nr'])   
    roc_crv <- roc(dt_basic_validation_current$label,predictions_p[,'true'], ci=TRUE)
    #conf_m <- confusionMatrix(predictions_bin, as.factor(dt_basic_validation_current$label))
    #accuracy <- conf_m$overall[1]
    #f1 <- 2*(conf_m$byClass[1]*conf_m$byClass[3])/(conf_m$byClass[1]+conf_m$byClass[3])
    auc <- roc_crv$auc[1]
    AUC_vs_states <- c(AUC_vs_states,auc)
    #F1_vs_states <- c(F1_vs_states,f1)
    #Accuracy_vs_states <- c(Accuracy_vs_states,accuracy)
    States <- c(States, st)
  }
  if(measure_optimiz=='AUC'){
    z <- AUC_vs_states
  } else if(measure_optimiz=='F1'){
    z <- F1_vs_states
  }else{
    z <- Accuracy_vs_states
  }
  optimal_state <- States[which.max(z)]
  return(c(max(z), optimal_state))
}

##----------##
# HMM+RF freq functions

# calculation of counts of each event(activity)
count_per_sequence <- function(data) {
  counts_ftable <- ftable(data[,c('sequence_nr','activity_name_short')])
  counts_frame <- suppressMessages(dcast(as.data.frame(counts_ftable), 
                                         as.formula(paste(paste(names(attr(counts_ftable, "row.vars")), 
                                                                collapse="+"), "~", paste(names(attr(counts_ftable, "col.vars")))))))
  colnames(counts_frame)[!colnames(counts_frame)%in%'sequence_nr'] <- paste('count',1:(ncol(counts_frame)-1), sep='_')
  return(counts_frame)  
}


count_per_sequence_events <- function(data) {
  counts_ftable <- ftable(data[,c('sequence_nr','activity_name_short')])
  counts_frame <- suppressMessages(dcast(as.data.frame(counts_ftable), 
                                         as.formula(paste(paste(names(attr(counts_ftable, "row.vars")), 
                                                                collapse="+"), "~", paste(names(attr(counts_ftable, "col.vars")))))))
  #colnames(counts_frame)[!colnames(counts_frame)%in%'sequence_nr'] <- paste('count',1:(ncol(counts_frame)-1), sep='_')
  return(counts_frame)  
}

#------------MAIN PROCEDURE------------------#
# general function for rf results
rf_run <- function(prefix_range, dt_train=division[[1]], dt_validation=division[[2]], dt_test=division[[3]], to_exclude){
  if (getwd() != paste(path,'/', output_folder, sep='')) {
    setwd(paste(path,'/', output_folder, sep=''))
  }
  method='RF_index'
  filecon <- paste("results_", method, ".txt",sep='') 
  filecon_general <- paste("general_results_", method,".txt",sep='')
  cat(paste("AUC", 'CI.AUC_L', 'CI.AUC_U', 'prefix', sep =', '),'\n', file=filecon_general, append=FALSE)
  
  features_factorize <- features_to_factorize(dt_train)
  
  for (feat in features_factorize[!features_factorize %in% to_exclude]) { #timestamp excluded here:
    print(feat)
    division_prep <- frequent_events_to_factors(train_data=dt_train, test_data=dt_test, validation_data=dt_validation, feature=feat)
    dt_train <- division_prep[[1]]
    dt_test <- division_prep[[3]]
    dt_validation <- division_prep[[2]]
  }
  colnms_in_dt <- setdiff(names(division_prep[[1]]), features_factorize) #exclude features that we factorized for the model, currently activity is duplicated (-6)
  pattern_to_grep <- paste0(colnms_in_dt, collapse="|")
  dt_basic <- data_to_model_with_selected_features(list_data=division_prep, features=pattern_to_grep)
  dt_basic_train <- subset(dt_basic, is_train==1)
  dt_basic_train$is_train <- NULL
  dt_basic_validation <-  subset(dt_basic, is_train==2)
  dt_basic_validation$is_train <- NULL
  dt_basic_test <-  subset(dt_basic, is_train==3)
  dt_basic_test$is_train <- NULL
  
  pb <- txtProgressBar(min = 0, max = prefix_range[length(prefix_range)], style = 3)
  for (prefix in prefix_range){
    dt_basic_train_current <- seq_ordering(dt_basic_train, nr_of_events=c(1:prefix))
    dt_basic_test_current <- seq_ordering(dt_basic_test, nr_of_events=c(1:prefix))
    dt_basic_train_current$is_train=1
    dt_basic_test_current$is_train=0
    dt <- rbind.data.frame(dt_basic_train_current, dt_basic_test_current)  
    
    feat <- names(dt)[!names(dt) %in% c('sequence_nr','activity_nr','label','seq_order','is_train','Number_of_executions','Age')] #choose features to produce a sequence
    res <- list()
    for (feature in feat){
      #print(feature)
      dt_long <- suppressWarnings(dcast(dt, label+sequence_nr+is_train ~ seq_order, value.var=feature, fill="NA"))
      res[[feature]] <- dt_long
    }
    merged_res = suppressWarnings(Reduce(function(...) merge(...,by=c("sequence_nr","label",'is_train'),all=TRUE), res))
    merged_res[is.na(merged_res)] <- 'NA'
    #merged_res <- merged_res[complete.cases(merged_res),]
    colnames(merged_res) <- c('sequence_nr','label','is_train', paste(rep(feat,each=prefix),'_seq_order_',rep(1:prefix,length(feat)),sep=''))
    merged_res <- as.data.frame(unclass(merged_res))
    merged_res <- merge(merged_res, unique(dt[,c('sequence_nr','is_train','Age','Number_of_executions')]), by=c('sequence_nr','is_train'), all.x=T)
    
    subdf <- lapply(merged_res, function(x) if(is.factor(x) & nlevels(x)==1) x <- NULL else x)
    merged_res <- as.data.frame(subdf[-(which(sapply(subdf,is.null),arr.ind=TRUE))])
    dt_basic_train_current <- subset(merged_res, is_train==1)
    dt_basic_test_current <- subset(merged_res, is_train==0)
    dt_basic_train_current$is_train <- NULL
    dt_basic_test_current$is_train <- NULL
    
    model <- randomForest(data=dt_basic_train_current[,!names(dt_basic_train_current)%in% 'sequence_nr'], label~., na.action=na.omit)

    predictions_bin <- predict(model, dt_basic_test_current[,!names(dt_basic_test_current)%in% 'sequence_nr']) 
    predictions_p <- predict(model, dt_basic_test_current[,!names(dt_basic_test_current)%in% 'sequence_nr'], type='prob')
    roc_crv <- roc(dt_basic_test_current$label, predictions_p[,'true'], ci=TRUE)
    CI_L <- roc_crv$ci[1]
    AUC <- roc_crv$auc[1]
    CI_U <- roc_crv$ci[3]  
    dt_basic_test_current$predictions_true <- predictions_p[,'true']
    prediction_for_sequence <- dt_basic_test_current[c('sequence_nr','label','predictions_true')]
    prediction_for_sequence$method <- method
    prediction_for_sequence$nr_prefixes <- prefix
    prediction_for_sequence <- merge(prediction_for_sequence, sequence_length_test, by='sequence_nr', all.x=TRUE)
    prediction_for_sequence$relative_earliness <- prediction_for_sequence$nr_prefixes/prediction_for_sequence$seq_length
    if (prefix==2) {
      write.table(prediction_for_sequence, filecon, col.names=T, row.names=F, quote=F, sep=',') 
    } else {
      write.table(prediction_for_sequence, filecon, col.names=F, row.names=F, quote=F, sep=',', append=T)
    } 
    cat(paste(AUC, CI_L, CI_U, prefix, sep=', '),'\n', file=filecon_general, append=TRUE)
    setTxtProgressBar(pb, prefix)
  }
  close(pb)
  setwd("../")
}

#-----------ADDITIONAL ALGORITHMS TO COMPARE WITH---------------#
svm_run <- function(prefix_range, dt_train=division[[1]], dt_validation=division[[2]], dt_test=division[[3]], to_exclude){
  if (getwd() != paste(path,'/', output_folder, sep='')) {
    setwd(paste(path,'/', output_folder, sep=''))
  }
  method='SVM_index'
  filecon <- paste("results_", method, ".txt",sep='') 
  filecon_general <- paste("general_results_", method,".txt",sep='')
  cat(paste("AUC", 'CI.AUC_L', 'CI.AUC_U', 'prefix', sep =', '),'\n', file=filecon_general, append=FALSE)
  
  features_factorize <- features_to_factorize(dt_train)
  
  for (feat in features_factorize[!features_factorize %in% to_exclude]) { #timestamp excluded here:
    #print(feat)
    division_prep <- frequent_events_to_factors(train_data=dt_train, test_data=dt_test, validation_data=dt_validation, feature=feat)
    dt_train <- division_prep[[1]]
    dt_test <- division_prep[[3]]
    dt_validation <- division_prep[[2]]
  }
  colnms_in_dt <- setdiff(names(division_prep[[1]]), features_factorize) #exclude features that we factorized for the model, currently activity is duplicated (-6)
  pattern_to_grep <- paste0(colnms_in_dt, collapse="|")
  dt_basic <- data_to_model_with_selected_features(list_data=division_prep, features=pattern_to_grep)
  dt_basic_train <- subset(dt_basic, is_train==1)
  dt_basic_train$is_train <- NULL
  dt_basic_validation <-  subset(dt_basic, is_train==2)
  dt_basic_validation$is_train <- NULL
  dt_basic_test <-  subset(dt_basic, is_train==3)
  dt_basic_test$is_train <- NULL
  
  pb <- txtProgressBar(min = 0, max = prefix_range[length(prefix_range)], style = 3)
  for (prefix in prefix_range){
    dt_basic_train_current <- seq_ordering(dt_basic_train, nr_of_events=c(1:prefix))
    dt_basic_test_current <- seq_ordering(dt_basic_test, nr_of_events=c(1:prefix))
    dt_basic_train_current$is_train=1
    dt_basic_test_current$is_train=0
    dt <- rbind.data.frame(dt_basic_train_current, dt_basic_test_current)  
    
    feat <- names(dt)[! names(dt) %in% c('sequence_nr','activity_nr','label','seq_order','is_train','Number_of_executions','Age')] #choose features to produce a sequence
    res <- list()
    for (feature in feat){
      #print(feature)
      dt_long <- suppressWarnings(dcast(dt, label+sequence_nr+is_train ~ seq_order, value.var=feature, fill="NA"))
      res[[feature]] <- dt_long
    }
    merged_res = suppressWarnings(Reduce(function(...) merge(...,by=c("sequence_nr","label",'is_train'),all=TRUE), res))
    merged_res[is.na(merged_res)] <- 'NA'
    #merged_res <- merged_res[complete.cases(merged_res),]
    colnames(merged_res) <- c('sequence_nr','label','is_train', paste(rep(feat,each=prefix),'_seq_order_',rep(1:prefix,length(feat)),sep=''))
    merged_res <- as.data.frame(unclass(merged_res))
    merged_res <- merge(merged_res, unique(dt[,c('sequence_nr','is_train','Age','Number_of_executions')]), by=c('sequence_nr','is_train'), all.x=T)
    
    subdf <- lapply(merged_res, function(x) if(is.factor(x) & nlevels(x)==1) x <- NULL else x)
    merged_res <- as.data.frame(subdf[-(which(sapply(subdf,is.null),arr.ind=TRUE))])
    dt_basic_train_current <- subset(merged_res, is_train==1)
    dt_basic_test_current <- subset(merged_res, is_train==0)
    dt_basic_train_current$is_train <- NULL
    dt_basic_test_current$is_train <- NULL
    
    model <- svm(data=dt_basic_train_current[,!names(dt_basic_train_current)%in% 'sequence_nr'], label~., probability=TRUE, na.action=na.omit, kernel="radial")
    
    predictions_bin <- predict(model, dt_basic_test_current[,!names(dt_basic_test_current)%in% 'sequence_nr']) 
    predictions_p <- attr(predict(model, dt_basic_test_current[,!names(dt_basic_test_current)%in% 'sequence_nr'], probability=TRUE), "probabilities")

    roc_crv <- roc(dt_basic_test_current$label, predictions_p[,'true'], ci=TRUE)
    CI_L <- roc_crv$ci[1]
    AUC <- roc_crv$auc[1]
    CI_U <- roc_crv$ci[3]  
    dt_basic_test_current$predictions_true <- predictions_p[,'true']
    prediction_for_sequence <- dt_basic_test_current[c('sequence_nr','label','predictions_true')]
    prediction_for_sequence$method <- method
    prediction_for_sequence$nr_prefixes <- prefix
    prediction_for_sequence <- merge(prediction_for_sequence, sequence_length_test, by='sequence_nr', all.x=TRUE)
    prediction_for_sequence$relative_earliness <- prediction_for_sequence$nr_prefixes/prediction_for_sequence$seq_length
    if (prefix==2) {
      write.table(prediction_for_sequence, filecon, col.names=T, row.names=F, quote=F, sep=',') 
    } else {
      write.table(prediction_for_sequence, filecon, col.names=F, row.names=F, quote=F, sep=',', append=T)
    } 
    cat(paste(AUC, CI_L, CI_U, prefix, sep=', '),'\n', file=filecon_general, append=TRUE)
    setTxtProgressBar(pb, prefix)
  }
  close(pb)
  setwd("../")
}

rpart_run <- function(prefix_range, dt_train=division[[1]], dt_validation=division[[2]], dt_test=division[[3]], to_exclude){
  if (getwd() != paste(path,'/', output_folder, sep='')) {
    setwd(paste(path,'/', output_folder, sep=''))
  }
  method='TREE_index'
  filecon <- paste("results_", method, ".txt",sep='') 
  filecon_general <- paste("general_results_", method,".txt",sep='')
  cat(paste("AUC", 'CI.AUC_L', 'CI.AUC_U', 'prefix', sep =', '),'\n', file=filecon_general, append=FALSE)
  
  features_factorize <- features_to_factorize(dt_train)
  
  for (feat in features_factorize[!features_factorize %in% to_exclude]) { #timestamp excluded here:
    #print(feat)
    division_prep <- frequent_events_to_factors(train_data=dt_train, test_data=dt_test, validation_data=dt_validation, feature=feat)
    dt_train <- division_prep[[1]]
    dt_test <- division_prep[[3]]
    dt_validation <- division_prep[[2]]
  }
  colnms_in_dt <- setdiff(names(division_prep[[1]]), features_factorize) #exclude features that we factorized for the model, currently activity is duplicated (-6)
  pattern_to_grep <- paste0(colnms_in_dt, collapse="|")
  dt_basic <- data_to_model_with_selected_features(list_data=division_prep, features=pattern_to_grep)
  dt_basic_train <- subset(dt_basic, is_train==1)
  dt_basic_train$is_train <- NULL
  dt_basic_validation <-  subset(dt_basic, is_train==2)
  dt_basic_validation$is_train <- NULL
  dt_basic_test <-  subset(dt_basic, is_train==3)
  dt_basic_test$is_train <- NULL
  
  pb <- txtProgressBar(min = 0, max = prefix_range[length(prefix_range)], style = 3)
  for (prefix in prefix_range){
    dt_basic_train_current <- seq_ordering(dt_basic_train, nr_of_events=c(1:prefix))
    dt_basic_test_current <- seq_ordering(dt_basic_test, nr_of_events=c(1:prefix))
    dt_basic_train_current$is_train=1
    dt_basic_test_current$is_train=0
    dt <- rbind.data.frame(dt_basic_train_current, dt_basic_test_current)  
    
    feat <- names(dt)[! names(dt) %in% c('sequence_nr','activity_nr','label','seq_order','is_train','Number_of_executions','Age')] #choose features to produce a sequence
    res <- list()
    for (feature in feat){
      #print(feature)
      dt_long <- suppressWarnings(dcast(dt, label+sequence_nr+is_train ~ seq_order, value.var=feature, fill="NA"))
      res[[feature]] <- dt_long
    }
    merged_res = suppressWarnings(Reduce(function(...) merge(...,by=c("sequence_nr","label",'is_train'),all=TRUE), res))
    merged_res[is.na(merged_res)] <- 'NA'
    #merged_res <- merged_res[complete.cases(merged_res),]
    colnames(merged_res) <- c('sequence_nr','label','is_train', paste(rep(feat,each=prefix),'_seq_order_',rep(1:prefix,length(feat)),sep=''))
    merged_res <- as.data.frame(unclass(merged_res))
    merged_res <- merge(merged_res, unique(dt[,c('sequence_nr','is_train','Age','Number_of_executions')]), by=c('sequence_nr','is_train'), all.x=T)    
    subdf <- lapply(merged_res, function(x) if(is.factor(x) & nlevels(x)==1) x <- NULL else x)
    merged_res <- as.data.frame(subdf[-(which(sapply(subdf,is.null),arr.ind=TRUE))])
 
    dt_basic_train_current <- subset(merged_res, is_train==1)
    dt_basic_test_current <- subset(merged_res, is_train==0)
    dt_basic_train_current$is_train <- NULL
    dt_basic_test_current$is_train <- NULL
    
    dt_to_use <- dt_basic_train_current[,!names(dt_basic_train_current)%in% 'sequence_nr']
    #dt_to_use <- droplevels(dt_to_use)
    #subdf <- lapply(dt_to_use, function(x) if(is.factor(x) & nlevels(x)==1) x <- NULL else x)
    #dt_to_use <- as.data.frame(subdf[-(which(sapply(subdf,is.null),arr.ind=TRUE))])
    model <- rpart(data=dt_to_use, label~., na.action=na.omit)
    #sapply(dt_basic_train_current[,!names(dt_basic_train_current)%in% 'sequence_nr'], function(x) is.factor(x) & nlevels(x)==2)
    predictions_bin <- predict(model, dt_basic_test_current[,!names(dt_basic_test_current)%in% 'sequence_nr'], type='class') 
    predictions_p <- predict(model, dt_basic_test_current[,!names(dt_basic_test_current)%in% 'sequence_nr'], type='prob')
    roc_crv <- roc(dt_basic_test_current$label, predictions_p[,'true'], ci=TRUE)
    CI_L <- roc_crv$ci[1]
    AUC <- roc_crv$auc[1]
    CI_U <- roc_crv$ci[3]  
    dt_basic_test_current$predictions_true <- predictions_p[,'true']
    prediction_for_sequence <- dt_basic_test_current[c('sequence_nr','label','predictions_true')]
    prediction_for_sequence$method <- method
    prediction_for_sequence$nr_prefixes <- prefix
    prediction_for_sequence <- merge(prediction_for_sequence, sequence_length_test, by='sequence_nr', all.x=TRUE)
    prediction_for_sequence$relative_earliness <- prediction_for_sequence$nr_prefixes/prediction_for_sequence$seq_length
    if (prefix==2) {
      write.table(prediction_for_sequence, filecon, col.names=T, row.names=F, quote=F, sep=',') 
    } else {
      write.table(prediction_for_sequence, filecon, col.names=F, row.names=F, quote=F, sep=',', append=T)
    } 
    cat(paste(AUC, CI_L, CI_U, prefix, sep=', '),'\n', file=filecon_general, append=TRUE)
    setTxtProgressBar(pb, prefix)
  }
  close(pb)
  setwd("../")
}

gbm_run <- function(prefix_range, dt_train=division[[1]], dt_validation=division[[2]], dt_test=division[[3]], to_exclude){
  if (getwd() != paste(path,'/', output_folder, sep='')) {
    setwd(paste(path,'/', output_folder, sep=''))
  }
  method='GBM_index'
  filecon <- paste("results_", method, ".txt",sep='') 
  filecon_general <- paste("general_results_", method,".txt",sep='')
  cat(paste("AUC", 'CI.AUC_L', 'CI.AUC_U', 'prefix', sep =', '),'\n', file=filecon_general, append=FALSE)
  
  features_factorize <- features_to_factorize(dt_train)
  
  for (feat in features_factorize[!features_factorize %in% to_exclude]) { #timestamp excluded here:
    #print(feat)
    division_prep <- frequent_events_to_factors(train_data=dt_train, test_data=dt_test, validation_data=dt_validation, feature=feat)
    dt_train <- division_prep[[1]]
    dt_test <- division_prep[[3]]
    dt_validation <- division_prep[[2]]
  }
  colnms_in_dt <- setdiff(names(division_prep[[1]]), features_factorize) #exclude features that we factorized for the model, currently activity is duplicated (-6)
  pattern_to_grep <- paste0(colnms_in_dt, collapse="|")
  dt_basic <- data_to_model_with_selected_features(list_data=division_prep, features=pattern_to_grep)
  dt_basic_train <- subset(dt_basic, is_train==1)
  dt_basic_train$is_train <- NULL
  dt_basic_validation <-  subset(dt_basic, is_train==2)
  dt_basic_validation$is_train <- NULL
  dt_basic_test <-  subset(dt_basic, is_train==3)
  dt_basic_test$is_train <- NULL
  
  pb <- txtProgressBar(min = 0, max = prefix_range[length(prefix_range)], style = 3)
  for (prefix in prefix_range){
    dt_basic_train_current <- seq_ordering(dt_basic_train, nr_of_events=c(1:prefix))
    dt_basic_test_current <- seq_ordering(dt_basic_test, nr_of_events=c(1:prefix))
    dt_basic_train_current$is_train=1
    dt_basic_test_current$is_train=0
    dt <- rbind.data.frame(dt_basic_train_current, dt_basic_test_current)  
    
    feat <- names(dt)[! names(dt) %in% c('sequence_nr','activity_nr','label','seq_order','is_train','Number_of_executions','Age')] #choose features to produce a sequence
    res <- list()
    for (feature in feat){
      #print(feature)
      dt_long <- suppressWarnings(dcast(dt, label+sequence_nr+is_train ~ seq_order, value.var=feature, fill="NA"))
      res[[feature]] <- dt_long
    }
    merged_res = suppressWarnings(Reduce(function(...) merge(...,by=c("sequence_nr","label",'is_train'),all=TRUE), res))
    merged_res[is.na(merged_res)] <- 'NA'
    #merged_res <- merged_res[complete.cases(merged_res),]
    colnames(merged_res) <- c('sequence_nr','label','is_train', paste(rep(feat,each=prefix),'_seq_order_',rep(1:prefix,length(feat)),sep=''))
    merged_res <- as.data.frame(unclass(merged_res))
    merged_res <- merge(merged_res, unique(dt[,c('sequence_nr','is_train','Age','Number_of_executions')]), by=c('sequence_nr','is_train'), all.x=T)
    
    subdf <- lapply(merged_res, function(x) if(is.factor(x) & nlevels(x)==1) x <- NULL else x)
    merged_res <- as.data.frame(subdf[-(which(sapply(subdf,is.null),arr.ind=TRUE))])
    dt_basic_train_current <- subset(merged_res, is_train==1)
    dt_basic_test_current <- subset(merged_res, is_train==0)
    dt_basic_train_current$is_train <- NULL
    dt_basic_test_current$is_train <- NULL
    dt_to_use <- dt_basic_train_current[,!names(dt_basic_train_current)%in% 'sequence_nr']
    dt_to_use$label <- ifelse(dt_to_use$label=='true',1,0)
    model <- gbm(data=dt_to_use, label~., distribution='bernoulli')
    
    dt_to_use_test <- dt_basic_test_current[,!names(dt_basic_test_current)%in% 'sequence_nr']
    dt_to_use_test$label <- ifelse(dt_to_use_test$label=='true',1,0)
    predictions_p <- predict(model, dt_to_use_test, n.trees=100, type='response') 
    predictions_bin <- ifelse(predictions_p>=0.5,"true","false")

    roc_crv <- roc(dt_basic_test_current$label, predictions_p, ci=TRUE)
    CI_L <- roc_crv$ci[1]
    AUC <- roc_crv$auc[1]
    CI_U <- roc_crv$ci[3]  
    dt_basic_test_current$predictions_true <- predictions_p
    prediction_for_sequence <- dt_basic_test_current[c('sequence_nr','label','predictions_true')]
    prediction_for_sequence$method <- method
    prediction_for_sequence$nr_prefixes <- prefix
    prediction_for_sequence <- merge(prediction_for_sequence, sequence_length_test, by='sequence_nr', all.x=TRUE)
    prediction_for_sequence$relative_earliness <- prediction_for_sequence$nr_prefixes/prediction_for_sequence$seq_length
    if (prefix==2) {
      write.table(prediction_for_sequence, filecon, col.names=T, row.names=F, quote=F, sep=',') 
    } else {
      write.table(prediction_for_sequence, filecon, col.names=F, row.names=F, quote=F, sep=',', append=T)
    } 
    cat(paste(AUC, CI_L, CI_U, prefix, sep=', '),'\n', file=filecon_general, append=TRUE)
    setTxtProgressBar(pb, prefix)
  }
  close(pb)
  setwd("../")
}

#general function for calculating hmms for range of prefixes and different states
hmm_run <- function(prefix_range, states_range, train_data=division[[1]], validation_data=division[[2]], test_data=division[[3]], features_for_hmm, iter=200){
  if (getwd() != paste(path,'/', output_folder, sep='')) {
    setwd(paste(path,'/', output_folder, sep=''))
  }
  method='HMM_ratios'
  filecon <- paste("results_", method,".txt",sep='') 
  filecon_general <- paste("general_results_", method,".txt",sep='')
  cat(paste("AUC", 'CI.AUC_L', 'CI.AUC_U', 'prefix', sep =', '),'\n', file=filecon_general, append=FALSE)
  pb <- txtProgressBar(min = 0, max = prefix_range[length(prefix_range)], style = 3)
  for (prefix in prefix_range) {
    best_states <- best_nr_states(train_data, validation_data, features=features_for_hmm,
                                  nr_prefixes = prefix, seed=seed_nr, measure_optimiz='AUC', states_range, iter)
    dt_current_hmm <- hmm_prefix_run(train_data, test_data, features=features_for_hmm, nr_prefixes = prefix, nr_states=best_states[2], seed=seed_nr, iter)
    dt_basic_train_current <- dt_current_hmm[[1]]
    dt_basic_test_current <- dt_current_hmm[[2]]
    model <- randomForest(data=dt_basic_train_current[,!names(dt_basic_train_current)%in% 'sequence_nr'], as.factor(label)~.)
    predictions_bin <- predict(model, dt_basic_test_current[,!names(dt_basic_test_current)%in% 'sequence_nr']) 
    predictions_p <- predict(model, dt_basic_test_current[,!names(dt_basic_test_current)%in% 'sequence_nr'], type='prob')
    roc_crv <- roc(dt_basic_test_current$label, predictions_p[,'true'], ci=TRUE)
    auc <- roc_crv$auc[1]    
    dt_basic_test_current$predictions_true <- predictions_p[,'true']
    prediction_for_sequence <- dt_basic_test_current[c('sequence_nr','label','predictions_true')]
    prediction_for_sequence$method <- method
    prediction_for_sequence$nr_prefixes <- prefix
    prediction_for_sequence <- merge(prediction_for_sequence, sequence_length_test, by='sequence_nr', all.x=TRUE)
    prediction_for_sequence$relative_earliness <- prediction_for_sequence$nr_prefixes/prediction_for_sequence$seq_length
    if (prefix == 2) {
      write.table(prediction_for_sequence, filecon, col.names=T, row.names=F, quote=F, sep=',') 
    } else {
      write.table(prediction_for_sequence, filecon, col.names=F, row.names=F, quote=F, sep=',', append=T)
    } 
    cat(paste(roc_crv$auc[1], roc_crv$ci[1], roc_crv$ci[3], prefix, sep=', '),'\n', file=filecon_general, append=TRUE)
    setTxtProgressBar(pb, prefix)
  }
  close(pb)
  setwd("../")
}
# general function that combines hmm and rf together
hmm_rf_run <- function(prefix_range, states_range, train_data, validation_data, test_data, features_for_hmm, to_exclude, iter=200) {
  if (getwd() != paste(path,'/', output_folder, sep='')) {
    setwd(paste(path,'/', output_folder, sep=''))
  }
  method='index_HMM_ratios'
  filecon <- paste("results_", method,".txt",sep='') 
  filecon_general <- paste("general_results_", method,".txt",sep='')
  cat(paste("AUC", 'CI.AUC_L', 'CI.AUC_U', 'prefix', sep =', '),'\n', file=filecon_general, append=FALSE)
  features_factorize <- features_to_factorize(train_data)  
  for (feat in features_factorize[!features_factorize %in% to_exclude]) {
    division_prep <- frequent_events_to_factors(train_data, test_data, validation_data, feature=feat)
    train_data <- division_prep[[1]]
    test_data <- division_prep[[3]]
    validation_data <- division_prep[[2]]
  }
  colnms_in_dt <- setdiff(names(division_prep[[1]]), features_factorize) #exclude features that we factorized for the model, currently activity is duplicated (-6)
  pattern_to_grep <- paste0(colnms_in_dt, collapse="|")
  dt_basic <- data_to_model_with_selected_features(list_data=division_prep, features=pattern_to_grep)
  dt_basic_train <- subset(dt_basic, is_train==1)
  dt_basic_train$is_train <- NULL
  dt_basic_validation <-  subset(dt_basic, is_train==2)
  dt_basic_validation$is_train <- NULL
  dt_basic_test <-  subset(dt_basic, is_train==3)
  dt_basic_test$is_train <- NULL  
  pb <- txtProgressBar(min = 0, max = prefix_range[length(prefix_range)], style = 3)
  for (prefix in prefix_range) {
    dt_basic_train_current <- seq_ordering(dt_basic_train, nr_of_events=c(1:prefix))
    dt_basic_validation_current <- seq_ordering(dt_basic_validation, nr_of_events=c(1:prefix))
    dt_basic_test_current <- seq_ordering(dt_basic_test, nr_of_events=c(1:prefix))
    dt_basic_train_current$is_train <- 1
    dt_basic_validation_current$is_train <- 2
    dt_basic_test_current$is_train <- 3
    dt <- rbind.data.frame(dt_basic_train_current, dt_basic_validation_current, dt_basic_test_current)
    feat_seq <- names(dt)[!names(dt) %in% c('sequence_nr','activity_nr','label','seq_order','is_train','Number_of_executions','Age')] #choose features to produce a sequence
    res <- list()
    for (feature in feat_seq) {
      #print(feature)
      dt_long <- suppressWarnings(dcast(dt, label+sequence_nr+is_train ~ seq_order, value.var=feature, fill="NA"))
      res[[feature]] <- dt_long
    }
    merged_res = suppressWarnings(Reduce(function(...) merge(..., by=c("sequence_nr","label",'is_train'), all=TRUE), res))
    merged_res[is.na(merged_res)] <- 'NA'
    colnames(merged_res) <- c('sequence_nr','label','is_train', paste(rep(feat_seq, each=prefix),'_seq_order_', rep(1:prefix, length(feat_seq)), sep=''))
    merged_res <- as.data.frame(unclass(merged_res))
    merged_res <- merge(merged_res, unique(dt[,c('sequence_nr','is_train','Age','Number_of_executions')]), by=c('sequence_nr','is_train'), all.x=T)
    subdf <- lapply(merged_res, function(x) if(is.factor(x) & nlevels(x)==1) x <- NULL else x)
    merged_res <- as.data.frame(subdf[-(which(sapply(subdf,is.null),arr.ind=TRUE))])
    
    dt_basic_train_current <- subset(merged_res, is_train==1)
    dt_basic_validation_current <- subset(merged_res, is_train==2)
    dt_basic_test_current <- subset(merged_res, is_train==3)
    best_states <- best_nr_states_hmmrf(train_data, validation_data, features=features_for_hmm,
                                        nr_prefixes = prefix, seed=seed_nr, measure_optimiz='AUC',
                                        rf_train=dt_basic_train_current, rf_validation=dt_basic_validation_current, 
                                        states_range, iter)
    
    dt_current_hmm <- hmm_prefix_run(train_data, test_data, features=features_for_hmm, nr_prefixes = prefix,
                                     nr_states=best_states[2], seed=seed_nr, iter)
    
    dt_basic_train_current_hmm <- dt_current_hmm[[1]]
    dt_basic_test_current_hmm <- dt_current_hmm[[2]]
    
    dt_basic_train_current <- merge(dt_basic_train_current, dt_basic_train_current_hmm, by=c('sequence_nr','label'), all=T)
    dt_basic_test_current <- merge(dt_basic_test_current, dt_basic_test_current_hmm, by=c('sequence_nr','label'), all=T)   
    dt_basic_train_current$is_train <- NULL
    dt_basic_test_current$is_train <- NULL
    model <- randomForest(data=dt_basic_train_current[,!names(dt_basic_train_current) %in% 'sequence_nr'], label~.)
    predictions_p <- predict(model, dt_basic_test_current[,!names(dt_basic_test_current) %in% 'sequence_nr'], type='prob')
    predictions_bin <- predict(model, dt_basic_test_current[,!names(dt_basic_test_current) %in% 'sequence_nr'])
    roc_crv <- roc(dt_basic_test_current$label,predictions_p[,'true'], ci=TRUE)
    CI_L <- roc_crv$ci[1]
    AUC <- roc_crv$auc[1]
    CI_U <- roc_crv$ci[3] 
    dt_basic_test_current$predictions_true <- predictions_p[,'true']
    prediction_for_sequence <- dt_basic_test_current[c('sequence_nr','label','predictions_true')]
    prediction_for_sequence$method <- method
    prediction_for_sequence$nr_prefixes <- prefix
    prediction_for_sequence <- merge(prediction_for_sequence, sequence_length_test, by='sequence_nr', all.x=TRUE)
    prediction_for_sequence$relative_earliness <- prediction_for_sequence$nr_prefixes/prediction_for_sequence$seq_length
    if (prefix == 2) {
      write.table(prediction_for_sequence, filecon, col.names=T, row.names=F, quote=F, sep=',') 
    } else {
      write.table(prediction_for_sequence, filecon, col.names=F, row.names=F, quote=F, sep=',', append=T)
    } 
    cat(paste(AUC, CI_L, CI_U, prefix, sep=', '),'\n', file=filecon_general, append=TRUE)
    setTxtProgressBar(pb, prefix)
  }
  close(pb)
  setwd("../")
}

#general function that adds to basic rf also fcounts of events
rf_with_frequencies_run <- function(prefix_range, dt_train=division[[1]], dt_validation=division[[2]], dt_test=division[[3]], to_exclude){
  if (getwd() != paste(path,'/', output_folder, sep='')) {
    setwd(paste(path,'/', output_folder, sep=''))
  }
    method='index_counts'
    filecon <- paste("results_", method, ".txt",sep='') 
    filecon_general <- paste("general_results_", method,".txt",sep='')
    cat(paste("AUC", 'CI.AUC_L', 'CI.AUC_U', 'prefix', sep =', '),'\n', file=filecon_general, append=FALSE)
    #separation of frequent_events to different factors, R cannot handle more than 33 levels in rf
    #define which features to factorize
    features_factorize <- features_to_factorize(dt_train[,!colnames(dt_train) %in% 'activity_name_short'])
    for (feat in features_factorize[!features_factorize %in% to_exclude]) { #timestamp excluded here:
      division_prep <- frequent_events_to_factors(train_data=dt_train, test_data=dt_test, validation_data=dt_validation, feature=feat)
      dt_train <- division_prep[[1]]
      dt_test <- division_prep[[3]]
      dt_validation <- division_prep[[2]]
    }
    colnms_in_dt <- setdiff(names(division_prep[[1]]), features_factorize) #exclude features that we factorized for the model, currently activity is duplicated (-6)
    pattern_to_grep <- paste0(colnms_in_dt, collapse="|")
    dt_basic <- data_to_model_with_selected_features(list_data=division_prep, features=pattern_to_grep)
    dt_basic_train <- subset(dt_basic, is_train==1)
    dt_basic_train$is_train <- NULL
    dt_basic_validation <-  subset(dt_basic, is_train==2)
    dt_basic_validation$is_train <- NULL
    dt_basic_test <-  subset(dt_basic, is_train==3)
    dt_basic_test$is_train <- NULL
    
    pb <- txtProgressBar(min = 0, max = prefix_range[length(prefix_range)], style = 3)
    for (prefix in prefix_range){
      dt_basic_train_current <- seq_ordering(dt_basic_train, nr_of_events=c(1:prefix))
      dt_basic_test_current <- seq_ordering(dt_basic_test, nr_of_events=c(1:prefix))
      dt_basic_train_current$is_train=1
      dt_basic_test_current$is_train=0
      
      activity_counts_train <- count_per_sequence(dt_basic_train_current)
      activity_counts_test <- count_per_sequence(dt_basic_test_current)
      
      dt <- rbind.data.frame(dt_basic_train_current, dt_basic_test_current)  
      dt$activity_name_short <- NULL
      feat <- names(dt)[! names(dt) %in% c('sequence_nr','activity_nr','label','seq_order','is_train', 'Number_of_executions','Age')] #choose features to produce a sequence
      res <- list()
      for (feature in feat){
        #print(feature)
        dt_long <- suppressWarnings(dcast(dt, label+sequence_nr+is_train ~ seq_order, value.var=feature, fill="NA"))
        res[[feature]] <- dt_long
      }
      merged_res = suppressWarnings(Reduce(function(...) merge(...,by=c("sequence_nr","label",'is_train'),all=TRUE), res))
      merged_res[is.na(merged_res)] <- 'NA'
      colnames(merged_res) <- c('sequence_nr','label','is_train', paste(rep(feat,each=prefix),'_seq_order_',rep(1:prefix,length(feat)),sep=''))
      merged_res <- as.data.frame(unclass(merged_res))
      merged_res <- merge(merged_res, unique(dt[,c('sequence_nr','is_train','Age','Number_of_executions')]), by=c('sequence_nr','is_train'), all.x=T)
      subdf <- lapply(merged_res, function(x) if(is.factor(x) & nlevels(x)==1) x <- NULL else x)
      merged_res <- as.data.frame(subdf[-(which(sapply(subdf,is.null),arr.ind=TRUE))])
      dt_basic_train_current <- subset(merged_res, is_train==1)
      dt_basic_train_current <- merge(dt_basic_train_current, activity_counts_train, by='sequence_nr', all.x=TRUE)
      dt_basic_train_current$is_train <- NULL
      
      dt_basic_test_current <- subset(merged_res, is_train==0)
      dt_basic_test_current <- merge(dt_basic_test_current, activity_counts_test, by='sequence_nr', all.x=TRUE)
      dt_basic_test_current$is_train <- NULL
      
      model <- randomForest(data=dt_basic_train_current[,!names(dt_basic_train_current) %in% 'sequence_nr'], label~.)
      predictions_bin <- predict(model, dt_basic_test_current[,!names(dt_basic_test_current) %in% 'sequence_nr']) 
      predictions_p <- predict(model, dt_basic_test_current[,!names(dt_basic_test_current) %in% 'sequence_nr'], type='prob')
      roc_crv <- roc(dt_basic_test_current$label,predictions_p[,'true'], ci=TRUE)
      CI_L <- roc_crv$ci[1]
      AUC <- roc_crv$auc[1]
      CI_U <- roc_crv$ci[3] 
      dt_basic_test_current$predictions_true <- predictions_p[,'true']
      prediction_for_sequence <- dt_basic_test_current[c('sequence_nr','label','predictions_true')]
      prediction_for_sequence$method <- method
      prediction_for_sequence$nr_prefixes <- prefix
      prediction_for_sequence <- merge(prediction_for_sequence, sequence_length_test, by='sequence_nr', all.x=TRUE)
      prediction_for_sequence$relative_earliness <- prediction_for_sequence$nr_prefixes/prediction_for_sequence$seq_length      
      if (prefix==2) {
        write.table(prediction_for_sequence, filecon, col.names=T, row.names=F, quote=F, sep=',') 
      } else {
        write.table(prediction_for_sequence, filecon, col.names=F, row.names=F, quote=F, sep=',', append=T)
      } 
      cat(paste(AUC, CI_L, CI_U, prefix, sep=', '),'\n', file=filecon_general, append=TRUE)
      setTxtProgressBar(pb, prefix)
    }
    close(pb)
    setwd("../")
}

#general function that adds to rf with frequencies also hmms, also fits optimal nr of states
hmm_rf_with_frequencies_run <- function(prefix_range, states_range, train_data, validation_data, test_data, features_for_hmm, to_exclude, iter){
  if (getwd() != paste(path,'/', output_folder, sep='')) {
    setwd(paste(path,'/', output_folder, sep=''))
  }
  method='index_counts_HMM'
  filecon <- paste("results_", method,".txt",sep='') 
  filecon_general <- paste("general_results_", method,".txt",sep='')
  cat(paste("AUC", 'CI.AUC_L', 'CI.AUC_U', 'prefix', sep =', '),'\n', file=filecon_general, append=FALSE)
  features_factorize <- features_to_factorize(train_data[,!colnames(train_data) %in% 'activity_name_short'])
  #features_factorize <- features_to_factorize(train_data)  
  for (feat in features_factorize[!features_factorize %in% to_exclude]) { #timestamp excluded here:
    division_prep <- frequent_events_to_factors(train_data, test_data, validation_data, feature=feat)
    train_data <- division_prep[[1]]
    test_data <- division_prep[[3]]
    validation_data <- division_prep[[2]]
  }
  colnms_in_dt <- setdiff(names(division_prep[[1]]), features_factorize) #exclude features that we factorized for the model, currently activity is duplicated (-6)
  pattern_to_grep <- paste0(colnms_in_dt, collapse="|")
  dt_basic <- data_to_model_with_selected_features(list_data=division_prep, features=pattern_to_grep)
  dt_basic_train <- subset(dt_basic, is_train==1)
  dt_basic_train$is_train <- NULL
  dt_basic_validation <-  subset(dt_basic, is_train==2)
  dt_basic_validation$is_train <- NULL
  dt_basic_test <-  subset(dt_basic, is_train==3)
  dt_basic_test$is_train <- NULL  
  pb <- txtProgressBar(min = 0, max = prefix_range[length(prefix_range)], style = 3)
  for(prefix in prefix_range){
    dt_basic_train_current <- seq_ordering(dt_basic_train, nr_of_events=c(1:prefix))
    dt_basic_validation_current <- seq_ordering(dt_basic_validation, nr_of_events=c(1:prefix))
    dt_basic_test_current <- seq_ordering(dt_basic_test, nr_of_events=c(1:prefix))
    dt_basic_train_current$is_train <- 1
    dt_basic_validation_current$is_train <- 2
    dt_basic_test_current$is_train <- 3
    
    activity_counts_train <- count_per_sequence(dt_basic_train_current)
    activity_counts_validation <- count_per_sequence(dt_basic_validation_current)
    activity_counts_test <- count_per_sequence(dt_basic_test_current)
    
    dt <- rbind.data.frame(dt_basic_train_current, dt_basic_validation_current, dt_basic_test_current)
    dt$activity_name_short <- NULL
    feat_seq <- names(dt)[!names(dt) %in% c('sequence_nr','activity_nr','label','seq_order','is_train','Number_of_executions','Age')] #choose features to produce a sequence
    res <- list()
    for (feature in feat_seq) {
      #print(feature)
      dt_long <- suppressWarnings(dcast(dt, label+sequence_nr+is_train ~ seq_order, value.var=feature, fill="NA"))
      res[[feature]] <- dt_long
    }
    merged_res = suppressWarnings(Reduce(function(...) merge(...,by=c("sequence_nr","label",'is_train'), all=TRUE), res))
    merged_res[is.na(merged_res)] <- 'NA'
    colnames(merged_res) <- c('sequence_nr','label','is_train', paste(rep(feat_seq, each=prefix),'_seq_order_', rep(1:prefix, length(feat_seq)), sep=''))
    merged_res <- as.data.frame(unclass(merged_res))
    merged_res <- merge(merged_res, unique(dt[,c('sequence_nr','is_train','Age','Number_of_executions')]), by=c('sequence_nr','is_train'), all.x=T)
    subdf <- lapply(merged_res, function(x) if(is.factor(x) & nlevels(x)==1) x <- NULL else x)
    merged_res <- as.data.frame(subdf[-(which(sapply(subdf,is.null),arr.ind=TRUE))])
    
    dt_basic_train_current <- subset(merged_res, is_train==1)
    dt_basic_train_current <- merge(dt_basic_train_current, activity_counts_train, by='sequence_nr', all.x=TRUE)
        
    dt_basic_validation_current <- subset(merged_res, is_train==2)
    dt_basic_validation_current <- merge(dt_basic_validation_current, activity_counts_validation, by='sequence_nr', all.x=TRUE)
    
    dt_basic_test_current <- subset(merged_res, is_train==3)
    dt_basic_test_current <- merge(dt_basic_test_current, activity_counts_test, by='sequence_nr', all.x=TRUE)
    
    best_states <- best_nr_states_hmmrf(train_data, validation_data, features=features_for_hmm,
                                        nr_prefixes = prefix, seed=seed_nr, measure_optimiz='AUC',
                                        rf_train=dt_basic_train_current, rf_validation=dt_basic_validation_current, 
                                        states_range, iter)
    
    dt_current_hmm <- hmm_prefix_run(train_data, test_data, features=features_for_hmm, nr_prefixes = prefix,
                                     nr_states=best_states[2], seed=seed_nr, iter)
    
    dt_basic_train_current_hmm <- dt_current_hmm[[1]]
    dt_basic_test_current_hmm <- dt_current_hmm[[2]]
    
    dt_basic_train_current <- merge(dt_basic_train_current, dt_basic_train_current_hmm, by=c('sequence_nr','label'), all=T)
    dt_basic_test_current <- merge(dt_basic_test_current, dt_basic_test_current_hmm, by=c('sequence_nr','label'), all=T)   
    dt_basic_train_current$is_train <- NULL
    dt_basic_test_current$is_train <- NULL
    
    model <- randomForest(data=dt_basic_train_current[,!names(dt_basic_train_current) %in% 'sequence_nr'], label~.)
    predictions_p <- predict(model, dt_basic_test_current[,!names(dt_basic_test_current) %in% 'sequence_nr'], type='prob')
    predictions_bin <- predict(model, dt_basic_test_current[,!names(dt_basic_test_current) %in% 'sequence_nr'])
    roc_crv <- roc(dt_basic_test_current$label,predictions_p[,'true'], ci=TRUE)
    CI_L <- roc_crv$ci[1]
    AUC <- roc_crv$auc[1]
    CI_U <- roc_crv$ci[3]  
    dt_basic_test_current$predictions_true <- predictions_p[,'true']
    prediction_for_sequence <- dt_basic_test_current[c('sequence_nr','label','predictions_true')]
    prediction_for_sequence$method <- method
    prediction_for_sequence$nr_prefixes <- prefix
    prediction_for_sequence <- merge(prediction_for_sequence, sequence_length_test, by='sequence_nr', all.x=TRUE)
    prediction_for_sequence$relative_earliness <- prediction_for_sequence$nr_prefixes/prediction_for_sequence$seq_length    
    if (prefix == 2) {
      write.table(prediction_for_sequence, filecon, col.names=T, row.names=F, quote=F, sep=',') 
    } else {
      write.table(prediction_for_sequence, filecon, col.names=F, row.names=F, quote=F, sep=',', append=T)
    } 
    cat(paste(AUC, CI_L, CI_U, prefix, sep=', '),'\n', file=filecon_general, append=TRUE)
    setTxtProgressBar(pb, prefix)
  }
  close(pb)
  setwd("../")
} 
  

##--------------##
#---BASELINES
indexes_events_run <- function(prefix_range, dt_train=division[[1]], dt_validation=division[[2]], 
                               dt_test=division[[3]], activity_feature="activity_name_short"){
  if (getwd() != paste(path,'/', output_folder, sep='')) {
    setwd(paste(path,'/', output_folder, sep=''))
  }
  method='index_events'
  filecon <- paste("results_", method, ".txt",sep='') 
  filecon_general <- paste("general_results_", method,".txt",sep='')
  cat(paste("AUC", 'CI.AUC_L', 'CI.AUC_U', 'prefix', sep =', '),'\n', file=filecon_general, append=FALSE)
  division_prep <- frequent_events_to_factors(train_data=dt_train, test_data=dt_test, validation_data=dt_validation, feature=activity_feature)
  dt_train <- division_prep[[1]]
  dt_test <- division_prep[[3]]
  dt_validation <- division_prep[[2]]
  
  colnms_in_dt <- grep(paste('activity_nr|sequence_nr|label|', activity_feature, sep=''), names(dt_train),value = TRUE)
  colnms_in_dt <- colnms_in_dt[!colnms_in_dt%in%activity_feature]
  pattern_to_grep <- paste0(colnms_in_dt, collapse="|")
  dt_basic <- data_to_model_with_selected_features(list_data=division_prep, features=pattern_to_grep)
  dt_basic_train <- subset(dt_basic, is_train==1)
  dt_basic_train$is_train <- NULL
  dt_basic_validation <-  subset(dt_basic, is_train==2)
  dt_basic_validation$is_train <- NULL
  dt_basic_test <-  subset(dt_basic, is_train==3)
  dt_basic_test$is_train <- NULL
  
  pb <- txtProgressBar(min = 0, max = prefix_range[length(prefix_range)], style = 3)
  for (prefix in prefix_range){
    dt_basic_train_current <- seq_ordering(dt_basic_train, nr_of_events=c(1:prefix))
    dt_basic_test_current <- seq_ordering(dt_basic_test, nr_of_events=c(1:prefix))
    dt_basic_train_current$is_train=1
    dt_basic_test_current$is_train=0
    dt <- rbind.data.frame(dt_basic_train_current, dt_basic_test_current)  
    
    feat <- names(dt)[! names(dt) %in% c('sequence_nr','activity_nr','label','seq_order','is_train')] #choose features to produce a sequence
    res <- list()
    for (feature in feat){
      #print(feature)
      dt_long <- suppressWarnings(dcast(dt, label+sequence_nr+is_train ~ seq_order, value.var=feature, fill="NA"))
      res[[feature]] <- dt_long
    }
    merged_res = suppressWarnings(Reduce(function(...) merge(...,by=c("sequence_nr","label",'is_train'),all=TRUE), res))
    merged_res[is.na(merged_res)] <- 'NA'
    colnames(merged_res) <- c('sequence_nr','label','is_train', paste(rep(feat,each=prefix),'_seq_order_',rep(1:prefix,length(feat)),sep=''))
    merged_res <- as.data.frame(unclass(merged_res))
    merged_res <- merge(merged_res, unique(dt[,c('sequence_nr','is_train')]), by=c('sequence_nr','is_train'), all.x=T)
    
    subdf <- lapply(merged_res, function(x) if(is.factor(x) & nlevels(x)==1) x <- NULL else x)
    merged_res <- as.data.frame(subdf[-(which(sapply(subdf,is.null),arr.ind=TRUE))])
    dt_basic_train_current <- subset(merged_res, is_train==1)
    dt_basic_test_current <- subset(merged_res, is_train==0)
    dt_basic_train_current$is_train <- NULL
    dt_basic_test_current$is_train <- NULL
    
    model <- randomForest(data=dt_basic_train_current[,!names(dt_basic_train_current)%in% 'sequence_nr'], label~.)
    predictions_bin <- predict(model, dt_basic_test_current[,!names(dt_basic_test_current)%in% 'sequence_nr']) 
    predictions_p <- predict(model, dt_basic_test_current[,!names(dt_basic_test_current)%in% 'sequence_nr'], type='prob')
    roc_crv <- roc(dt_basic_test_current$label, predictions_p[,'true'], ci=TRUE)
    CI_L <- roc_crv$ci[1]
    AUC <- roc_crv$auc[1]
    CI_U <- roc_crv$ci[3]  
    dt_basic_test_current$predictions_true <- predictions_p[,'true']
    prediction_for_sequence <- dt_basic_test_current[c('sequence_nr','label','predictions_true')]
    prediction_for_sequence$method <- method
    prediction_for_sequence$nr_prefixes <- prefix
    prediction_for_sequence <- merge(prediction_for_sequence, sequence_length_test, by='sequence_nr', all.x=TRUE)
    prediction_for_sequence$relative_earliness <- prediction_for_sequence$nr_prefixes/prediction_for_sequence$seq_length
    if (prefix==2) {
      write.table(prediction_for_sequence, filecon, col.names=T, row.names=F, quote=F, sep=',') 
    } else {
      write.table(prediction_for_sequence, filecon, col.names=F, row.names=F, quote=F, sep=',', append=T)
    } 
    cat(paste(AUC, CI_L, CI_U, prefix, sep=', '),'\n', file=filecon_general, append=TRUE)
    setTxtProgressBar(pb, prefix)
  }
  close(pb)
  setwd("../")
}

frequencies_only <- function(prefix_range, dt_train=division[[1]], dt_validation=division[[2]], dt_test=division[[3]]){
  if (getwd() != paste(path,'/', output_folder, sep='')) {
    setwd(paste(path,'/', output_folder, sep=''))
  }
  method='counts_control_flow'
  filecon <- paste("results_", method, ".txt",sep='') 
  filecon_general <- paste("general_results_", method,".txt",sep='')
  cat(paste("AUC", 'CI.AUC_L', 'CI.AUC_U', 'prefix', sep =', '),'\n', file=filecon_general, append=FALSE)
  
  pb <- txtProgressBar(min = 0, max = prefix_range[length(prefix_range)], style = 3)
  for (prefix in prefix_range){
    dt_basic_train_current <- seq_ordering(dt_train, nr_of_events=c(1:prefix))
    dt_basic_test_current <- seq_ordering(dt_test, nr_of_events=c(1:prefix))
    
    activity_counts_train <- count_per_sequence_events(dt_basic_train_current)
    activity_counts_test <- count_per_sequence_events(dt_basic_test_current)
    
    exist_only_test <- setdiff(colnames(activity_counts_test), colnames(activity_counts_train))
    activity_counts_test <- activity_counts_test[!colnames(activity_counts_test) %in% exist_only_test]
    exist_only_train <- setdiff(colnames(activity_counts_train), colnames(activity_counts_test))
    for ( i in exist_only_train){
      activity_counts_test[i] <- 0 
    }
    matching_table <- data.frame(original_names=colnames(activity_counts_train)[!colnames(activity_counts_train)%in%'sequence_nr'], replaced=paste('count',1:(ncol(activity_counts_train)-1), sep='_'))
    colnames(activity_counts_train)[!colnames(activity_counts_train)%in%'sequence_nr'] <- paste('count',1:(ncol(activity_counts_train)-1), sep='_') 
    match_idx <- match(colnames(activity_counts_test)[!colnames(activity_counts_test)%in%'sequence_nr'], matching_table$original_names)
    colnames(activity_counts_test)[!colnames(activity_counts_test)%in%'sequence_nr'] <- as.character(matching_table[match_idx,2])
    
    dt_basic_train_current <- merge(activity_counts_train, by='sequence_nr', unique(dt_basic_train_current[c('sequence_nr','label')]), all.x=TRUE)
    dt_basic_test_current <- merge(activity_counts_test, by='sequence_nr', unique(dt_basic_test_current[c('sequence_nr','label')]), all.x=TRUE)    
    
    model <- randomForest(data=dt_basic_train_current[,!names(dt_basic_train_current) %in% 'sequence_nr'], as.factor(label)~.)
    predictions_bin <- predict(model, dt_basic_test_current[,!names(dt_basic_test_current) %in% 'sequence_nr']) 
    predictions_p <- predict(model, dt_basic_test_current[,!names(dt_basic_test_current) %in% 'sequence_nr'], type='prob')
    roc_crv <- roc(dt_basic_test_current$label,predictions_p[,'true'], ci=TRUE)
    CI_L <- roc_crv$ci[1]
    AUC <- roc_crv$auc[1]
    CI_U <- roc_crv$ci[3] 
    dt_basic_test_current$predictions_true <- predictions_p[,'true']
    prediction_for_sequence <- dt_basic_test_current[c('sequence_nr','label','predictions_true')]
    prediction_for_sequence$method <- method
    prediction_for_sequence$nr_prefixes <- prefix
    prediction_for_sequence <- merge(prediction_for_sequence, sequence_length_test, by='sequence_nr', all.x=TRUE)
    prediction_for_sequence$relative_earliness <- prediction_for_sequence$nr_prefixes/prediction_for_sequence$seq_length      
    if (prefix==2) {
      write.table(prediction_for_sequence, filecon, col.names=T, row.names=F, quote=F, sep=',') 
    } else {
      write.table(prediction_for_sequence, filecon, col.names=F, row.names=F, quote=F, sep=',', append=T)
    } 
    cat(paste(AUC, CI_L, CI_U, prefix, sep=', '),'\n', file=filecon_general, append=TRUE)
    setTxtProgressBar(pb, prefix)
  }
  close(pb)
  setwd("../")
}

events_boolean <- function(prefix_range, dt_train=division[[1]], dt_validation=division[[2]], dt_test=division[[3]]){
  if (getwd() != paste(path,'/', output_folder, sep='')) {
    setwd(paste(path,'/', output_folder, sep=''))
  }
  method='boolean_control_flow'
  filecon <- paste("results_", method, ".txt",sep='') 
  filecon_general <- paste("general_results_", method,".txt",sep='')
  cat(paste("AUC", 'CI.AUC_L', 'CI.AUC_U', 'prefix', sep =', '),'\n', file=filecon_general, append=FALSE)
  
  pb <- txtProgressBar(min = 0, max = prefix_range[length(prefix_range)], style = 3)
  for (prefix in prefix_range){
    dt_basic_train_current <- seq_ordering(dt_train, nr_of_events=c(1:prefix))
    dt_basic_test_current <- seq_ordering(dt_test, nr_of_events=c(1:prefix))
    #dt_basic_train_current <- seq_ordering(dt_train, nr_of_events=c(30:prefix))
    #dt_basic_test_current <- seq_ordering(dt_test, nr_of_events=c(30:prefix))
    
    activity_counts_train <- count_per_sequence_events(dt_basic_train_current)
    activity_counts_test <- count_per_sequence_events(dt_basic_test_current)
    
    exist_only_test <- setdiff(colnames(activity_counts_test), colnames(activity_counts_train))
    activity_counts_test <- activity_counts_test[!colnames(activity_counts_test) %in% exist_only_test]
    exist_only_train <- setdiff(colnames(activity_counts_train), colnames(activity_counts_test))
    for ( i in exist_only_train){
      activity_counts_test[i] <- 0 
    }
    matching_table <- data.frame(original_names=colnames(activity_counts_train)[!colnames(activity_counts_train)%in%'sequence_nr'], replaced=paste('count',1:(ncol(activity_counts_train)-1), sep='_'))
    colnames(activity_counts_train)[!colnames(activity_counts_train)%in%'sequence_nr'] <- paste('count',1:(ncol(activity_counts_train)-1), sep='_') 
    match_idx <- match(colnames(activity_counts_test)[!colnames(activity_counts_test)%in%'sequence_nr'], matching_table$original_names)
    colnames(activity_counts_test)[!colnames(activity_counts_test)%in%'sequence_nr'] <- as.character(matching_table[match_idx,2])
    
    dt_basic_train_current <- merge(activity_counts_train, by='sequence_nr', unique(dt_basic_train_current[c('sequence_nr','label')]), all.x=TRUE)
    dt_basic_test_current <- merge(activity_counts_test, by='sequence_nr', unique(dt_basic_test_current[c('sequence_nr','label')]), all.x=TRUE)    
    
    dt_basic_train_current[grep('count_',colnames(dt_basic_train_current), value=T)] <- (dt_basic_train_current[grep('count_',colnames(dt_basic_train_current), value = TRUE)]>=1)*1
    dt_basic_test_current[grep('count_',colnames(dt_basic_test_current), value=T)] <- (dt_basic_test_current[grep('count_',colnames(dt_basic_test_current), value = TRUE)]>=1)*1

    model <- randomForest(data=dt_basic_train_current[,!names(dt_basic_train_current) %in% 'sequence_nr'], as.factor(label)~.)
    summary(model)
    predictions_bin <- predict(model, dt_basic_test_current[,!names(dt_basic_test_current) %in% 'sequence_nr']) 
    predictions_p <- predict(model, dt_basic_test_current[,!names(dt_basic_test_current) %in% 'sequence_nr'], type='prob')
    roc_crv <- roc(dt_basic_test_current$label,predictions_p[,'true'], ci=TRUE)
    CI_L <- roc_crv$ci[1]
    AUC <- roc_crv$auc[1]
    CI_U <- roc_crv$ci[3] 
    dt_basic_test_current$predictions_true <- predictions_p[,'true']
    prediction_for_sequence <- dt_basic_test_current[c('sequence_nr','label','predictions_true')]
    prediction_for_sequence$method <- method
    prediction_for_sequence$nr_prefixes <- prefix
    prediction_for_sequence <- merge(prediction_for_sequence, sequence_length_test, by='sequence_nr', all.x=TRUE)
    prediction_for_sequence$relative_earliness <- prediction_for_sequence$nr_prefixes/prediction_for_sequence$seq_length      
    if (prefix==2) {
      write.table(prediction_for_sequence, filecon, col.names=T, row.names=F, quote=F, sep=',') 
    } else {
      write.table(prediction_for_sequence, filecon, col.names=F, row.names=F, quote=F, sep=',', append=T)
    } 
    cat(paste(AUC, CI_L, CI_U, prefix, sep=', '),'\n', file=filecon_general, append=TRUE)
    setTxtProgressBar(pb, prefix)
  }
  close(pb)
  setwd("../")
}


indexes_events_latest_data_run <- function(prefix_range, dt_train=division[[1]], dt_validation=division[[2]],to_exclude, 
                                           dt_test=division[[3]], activity_feature="activity_name_short"){
  if (getwd() != paste(path,'/', output_folder, sep='')) {
    setwd(paste(path,'/', output_folder, sep=''))
  }
  method='index_events_latest_data'
  filecon <- paste("results_", method, ".txt",sep='') 
  filecon_general <- paste("general_results_", method,".txt",sep='')
  cat(paste("AUC", 'CI.AUC_L', 'CI.AUC_U', 'prefix', sep =', '),'\n', file=filecon_general, append=FALSE)
  
  features_factorize <- features_to_factorize(dt_train)
  
  for (feat in features_factorize[!features_factorize %in% to_exclude]) { #timestamp excluded here:
    #print(feat)
    division_prep <- frequent_events_to_factors(train_data=dt_train, test_data=dt_test, validation_data=dt_validation, feature=feat)
    dt_train <- division_prep[[1]]
    dt_test <- division_prep[[3]]
    dt_validation <- division_prep[[2]]
  }
  colnms_in_dt <- setdiff(names(division_prep[[1]]), features_factorize) #exclude features that we factorized for the model, currently activity is duplicated (-6)
  pattern_to_grep <- paste0(colnms_in_dt, collapse="|")
  dt_basic <- data_to_model_with_selected_features(list_data=division_prep, features=pattern_to_grep)
  dt_basic_train <- subset(dt_basic, is_train==1)
  dt_basic_train$is_train <- NULL
  dt_basic_validation <-  subset(dt_basic, is_train==2)
  dt_basic_validation$is_train <- NULL
  dt_basic_test <-  subset(dt_basic, is_train==3)
  dt_basic_test$is_train <- NULL
  
  pb <- txtProgressBar(min = 0, max = prefix_range[length(prefix_range)], style = 3)
  for (prefix in prefix_range){
    dt_basic_train_current <- seq_ordering(dt_basic_train, nr_of_events=c(1:prefix))
    dt_basic_test_current <- seq_ordering(dt_basic_test, nr_of_events=c(1:prefix))
    dt_basic_train_current$is_train=1
    dt_basic_test_current$is_train=0
    dt <- rbind.data.frame(dt_basic_train_current, dt_basic_test_current)      
    
    res <- list()
    colnms_events <- grep(paste("freq", activity_feature, sep='_'), names(dt),value = TRUE)
    
    for (feature in colnms_events){
      #print(feature)
      dt_long <- suppressWarnings(dcast(dt, label+sequence_nr+is_train ~ seq_order, value.var=feature, fill="NA"))
      res[[feature]] <- dt_long
    }
    merged_res = suppressWarnings(Reduce(function(...) merge(...,by=c("sequence_nr","label",'is_train'),all=TRUE), res))
    merged_res[is.na(merged_res)] <- 'NA'
    colnames(merged_res) <- c('sequence_nr','label','is_train', paste(rep(colnms_events, each=prefix),'_seq_order_',rep(1:prefix, length(colnms_events)),sep=''))
    merged_res <- as.data.frame(unclass(merged_res))
    merged_res <- merge(merged_res, unique(dt[,c('sequence_nr','is_train')]), by=c('sequence_nr','is_train'), all.x=T)
    subdf <- lapply(merged_res, function(x) if(is.factor(x) & nlevels(x)==1) x <- NULL else x)
    merged_res <- as.data.frame(subdf[-(which(sapply(subdf,is.null), arr.ind=TRUE))])
    
    feat <- names(dt)[! names(dt) %in% c(grep(paste("freq", activity_feature, sep='_'), names(dt),value = TRUE))]
    dt_last_snapshot <- subset(dt[feat], seq_order==prefix)
    
    merged_res <- merge(merged_res, dt_last_snapshot[!(colnames(dt_last_snapshot)%in%c('activity_nr','seq_order'))], by=c('sequence_nr','is_train','label'),all.x=TRUE)
    dt_basic_train_current <- subset(merged_res, is_train==1)
    dt_basic_test_current <- subset(merged_res, is_train==0)
    dt_basic_train_current$is_train <- NULL
    dt_basic_test_current$is_train <- NULL
    
    model <- randomForest(data=dt_basic_train_current[,!names(dt_basic_train_current)%in% 'sequence_nr'], label~., na.action=na.omit)
    predictions_bin <- predict(model, dt_basic_test_current[,!names(dt_basic_test_current)%in% 'sequence_nr']) 
    predictions_p <- predict(model, dt_basic_test_current[,!names(dt_basic_test_current)%in% 'sequence_nr'], type='prob')
    roc_crv <- roc(dt_basic_test_current$label, predictions_p[,'true'], ci=TRUE)
    CI_L <- roc_crv$ci[1]
    AUC <- roc_crv$auc[1]
    CI_U <- roc_crv$ci[3]  
    dt_basic_test_current$predictions_true <- predictions_p[,'true']
    prediction_for_sequence <- dt_basic_test_current[c('sequence_nr','label','predictions_true')]
    prediction_for_sequence$method <- method
    prediction_for_sequence$nr_prefixes <- prefix
    prediction_for_sequence <- merge(prediction_for_sequence, sequence_length_test, by='sequence_nr', all.x=TRUE)
    prediction_for_sequence$relative_earliness <- prediction_for_sequence$nr_prefixes/prediction_for_sequence$seq_length
    if (prefix==2) {
      write.table(prediction_for_sequence, filecon, col.names=T, row.names=F, quote=F, sep=',') 
    } else {
      write.table(prediction_for_sequence, filecon, col.names=F, row.names=F, quote=F, sep=',', append=T)
    } 
    cat(paste(AUC, CI_L, CI_U, prefix, sep=', '),'\n', file=filecon_general, append=TRUE)
    setTxtProgressBar(pb, prefix)
  }
  close(pb)
  setwd("../")
}


##--------------##
#EVALUATION OF RESULTS
collecting_statistics_relative <- function(earliness_range, bootstrap=FALSE) {
  filenames <- list.files()[grep(paste("^result(?=.*\\.txt)",sep=''), list.files(), perl=TRUE)]
  method_names <- gsub("^(results_)+(.*)+(.txt)+$","\\2", filenames)
  auc_per_earliness <- data.frame()
  roc_objects <- list()
  ldf <- lapply(filenames, fread)  
  results <- do.call(rbind.data.frame, ldf)
  for (current_method in method_names) {
    dt_method <- subset(results, method==current_method)
    for(current_relative_earliness in earliness_range){
      dt_earl <- unique(subset(dt_method, relative_earliness<=current_relative_earliness))
      dt_earl2 <- ddply(dt_earl, .(sequence_nr), summarise, max_earl=max(relative_earliness), predict_max=predictions_true[relative_earliness==max_earl][1],
                        label=label[relative_earliness==max_earl][1], method=method[relative_earliness==max_earl][1])
      
      roc_crv <- roc(dt_earl2$label, dt_earl2$predict_max, ci=TRUE, of='auc')  
      roc_objects[[paste(current_method, current_relative_earliness*100, sep='_')]] <- roc_crv  
      if (bootstrap==TRUE){
        CI <- ci.auc(roc_crv, method='bootstrap', progress='none')
        CI_L <- CI[1]
        AUC <- CI[2]
        CI_U <- CI[3]  
      } else{
        CI_L <- roc_crv$ci[1]
        AUC <- roc_crv$auc[1]
        CI_U <- roc_crv$ci[3]  
      }
      auc_res_current <- data.frame(AUC, CI_L, CI_U, current_relative_earliness, current_method)
      colnames(auc_res_current) <- c("AUC","CI_L","CI_U","relative_earliness", "method")
      auc_per_earliness <- rbind.data.frame(auc_per_earliness, auc_res_current)
    }
  }
  return(list(auc_per_earliness, roc_objects))
}
collecting_statistics_absolute <- function(prefix_range, bootstrap=FALSE) {
  filenames <- list.files()[grep(paste("^result(?=.*\\.txt)",sep=''), list.files(), perl=TRUE)]
  method_names <- gsub("^(results_)+(.*)+(.txt)+$","\\2", filenames)
  auc_per_earliness <- data.frame()
  roc_objects <- list()
  ldf <- lapply(filenames, fread)  
  results <- do.call(rbind.data.frame, ldf)
  for (current_method in method_names) {
    dt_method <- subset(results, method==current_method)
    for(current_absolute_earliness in prefix_range){
      print(current_absolute_earliness)
      dt_earl <- unique(subset(dt_method, nr_prefixes == current_absolute_earliness))
      #dt_earl2 <- ddply(dt_earl, .(sequence_nr), summarise, max_earl=max(relative_earliness), predict_max=predictions_true[relative_earliness==max_earl][1],
      #                  label=label[relative_earliness==max_earl][1], method=method[relative_earliness==max_earl][1])
      
      roc_crv <- roc(dt_earl$label, dt_earl$predictions_true, ci=TRUE, of='auc')  
      roc_objects[[paste(current_method, current_absolute_earliness, sep='_')]] <- roc_crv  
      if (bootstrap==TRUE){
        CI <- ci.auc(roc_crv, method='bootstrap', progress='none')
        CI_L <- CI[1]
        AUC <- CI[2]
        CI_U <- CI[3]  
      } else{
        CI_L <- roc_crv$ci[1]
        AUC <- roc_crv$auc[1]
        CI_U <- roc_crv$ci[3]  
      }
      auc_res_current <- data.frame(AUC, CI_L, CI_U, current_absolute_earliness, current_method)
      colnames(auc_res_current) <- c("AUC","CI_L","CI_U","absolute_earliness", "method")
      auc_per_earliness <- rbind.data.frame(auc_per_earliness, auc_res_current)
    }
  }
  return(list(auc_per_earliness, roc_objects))
}


general_results <- function(){
  filenames <- list.files()[grep(paste("(?=.*general)(?=.*\\.txt)",sep=''),list.files(), perl=T)]
  methods_names <- gsub("^(general_results_)+(.*)+(.txt)+$","\\2", filenames)
  ldf <- lapply(filenames, fread)
  names(ldf) <- methods_names
  for (i in 1:length(ldf)) {
    ldf[[i]]$method <- names(ldf)[i] 
  }
  results <- as.data.frame(do.call(rbind.data.frame, ldf), use.names=F)
  results[,4] <- as.factor(str_trim(results[,4]))
  #friedman.test(results$AUC, as.factor(results$method), as.factor(results[,4]))
  dunn_test <- dunn.test(results$AUC, as.factor(results$method), as.factor(results[,4]), kw=T, method='bonferroni', wrap=TRUE)
  return(dunn_test)
}

roc_comparison_relative <- function(roc_objects, earliness_range) {
  filenames <- list.files()[grep(paste("^result(?=.*\\.txt)",sep=''), list.files(), perl=TRUE)]
  method_names <- gsub("^(results_)+(.*)+(.txt)+$","\\2", filenames)
  p_value_matrix_list <- list()
  for(j in earliness_range*100){
    p_value_matrix_list[[as.character(j)]] <- matrix(NA, ncol=length(method_names), nrow=length(method_names), dimnames=list(method_names,method_names))
    for(m1 in method_names){
      for(m2 in method_names){
        if (m1 != m2) {
          p_value_matrix_list[[as.character(j)]][m1,m2] <- roc.test(roc_objects[[paste(m1, j, sep='_')]], roc_objects[[paste(m2, j, sep='_')]])$p.value
        }
      }
    }
  }
  return(p_value_matrix_list)
}

roc_comparison_absolute <- function(roc_objects, prefix_range) {
  filenames <- list.files()[grep(paste("^result(?=.*\\.txt)",sep=''), list.files(), perl=TRUE)]
  method_names <- gsub("^(results_)+(.*)+(.txt)+$","\\2", filenames)
  p_value_matrix_list <- list()
  for(j in prefix_range){
    p_value_matrix_list[[as.character(j)]] <- matrix(NA, ncol=length(method_names), nrow=length(method_names), dimnames=list(method_names,method_names))
    for(m1 in method_names){
      for(m2 in method_names){
        if (m1 != m2) {
          p_value_matrix_list[[as.character(j)]][m1,m2] <- roc.test(roc_objects[[paste(m1, j, sep='_')]], roc_objects[[paste(m2, j, sep='_')]])$p.value
        }
      }
    }
  }
  return(p_value_matrix_list)
}

pairwise_comparison <- function(x){
  dt <- na.omit(melt(x))
  dt$pairwise_methods <- paste(dt$Var1, dt$Var2, sep='-')
  return(dt)
}

pairwise_comparison_function <- function(p_value_list, relative=TRUE){
  pairwise_comparison_list <- lapply(p_value_list, pairwise_comparison)
  if (relative==TRUE){
   for (i in 1:length(pairwise_comparison_list)){
     pairwise_comparison_list[[i]]$relative_earliness <- names(pairwise_comparison_list)[i]
    }
  } else {
     for (i in 1:length(pairwise_comparison_list)){
       pairwise_comparison_list[[i]]$absolute_earliness <- names(pairwise_comparison_list)[i]  
    }
  }
  pairwise_comparison_dt <- do.call(rbind.data.frame,pairwise_comparison_list) 
  colnames(pairwise_comparison_dt)[3] <- 'pvalue'
  return(pairwise_comparison_dt)
}







