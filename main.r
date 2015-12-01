# SET THE PATH
#the path should be the same as with the source file and Data
path <- getwd() # get the path
# creation of a folder, where the results will be saved 
output_folder <- 'results'

# PACKAGES AND FUNCTIONS
source('functions.r') 
 
# GLOBAL VARIABLES
# used for partitioning the data
seed_nr <- 40981

# for prefix range to explore
#prefix_range <- seq(2, 30, by=1) # prefix cannot be less than 2
prefix_range <- seq(2, 10, by=1) # Note: prefix cannot be less than 2

# for hidden states to fit
states_range <- 2 # note: states cannot be less than 2
# dynamic features for fitting hmm
features_for_hmm <- c('activity_name_short', 'Activity_code', 'Diagnosis_code','Diagnosis_Treatment_Combination_ID',
                      'Number_of_executions', 'Producer_code', 'Section', "Specialism_code", "Treatment_code",
                      "Diagnosis_short", "group_short")

iter <- 50 # number of iterations for hmm

# optional:
# for cross-validation ## 5-fold CV  repeated 3 times, slow!!
# fitControl <- trainControl(method = "cv", number = 2, verboseIter = TRUE)

# LOADING THE DATA
# Read and pre-process the data (takes some time)
dt <- 'Data/f1/BPI2011_80.txt' # path to a datafile
dat <- suppressWarnings(data_preprocess(dt))

# the same in case of a separate file for test data
test <- 'Data/f1/BPI2011_20.txt'
test_dat <- suppressWarnings(data_preprocess(test))

# if column names and values have spaces between the words, remove them:
names(test_dat)
dat_preprocessed <- rename_columns(data=dat) 
test_preprocessed <- rename_columns(data=test_dat)

# Abbreviate features for easy handling (R is not fine with too long variable names)
# choose features to abbreviate (if abbreviated, original feature is replaced with shortened version)

features_abbreviate <- c('activity_name','Diagnosis','group')
for(feat in features_abbreviate){
  dat_preprocessed <- abbrv_strings(data=dat_preprocessed, letters=TRUE, numbers=FALSE, feature=feat)
  test_preprocessed <- abbrv_strings(data=test_preprocessed, letters=TRUE, numbers=FALSE, feature=feat)
}

# record separately length of the sequences for the calculation of earliness as % of the length of a sequence
sequence_length_train <- as.data.frame(table(dat_preprocessed$sequence_nr))
colnames(sequence_length_train) <- c('sequence_nr', 'seq_length')
sequence_length_test <- as.data.frame(table(test_preprocessed$sequence_nr))
colnames(sequence_length_test) <- c('sequence_nr', 'seq_length')

# if there is a need to work with the full dataset, split it on train/validation/test. Validation is required for hidden state fit, 
# if there is a separate test set, split train set on just two parts and add test set as third

# division split 60% train, 20% validation, 20% test
division_full <- train_test_division(data=dat_preprocessed, ratio=0.6, seed=seed_nr)
division <- list(train=division_full[[1]], validation=division_full[[2]], test=test_preprocessed)

# for continius features change the type:
division <- lapply(division, function(x) as.data.frame(unclass(x)))
for (j in 1:length(division)) {
  division[[j]]$Age <- as.numeric(division[[j]]$Age)
  division[[j]]$Number_of_executions <- as.numeric(division[[j]]$Number_of_executions)
}

#create the folder 'output' for the results in the current path
system(paste('mkdir', output_folder, sep=' '))

# RANDOM FOREST
# features to exclude from the mode:
to_exclude <- c('Start_date','End_date','time.timestamp')

# RF index and comparison with other models
setwd(paste(path, output_folder, sep='/'))
rf_run(prefix_range, dt_train=division[[1]], dt_validation=division[[2]], dt_test=division[[3]], to_exclude)
svm_run(prefix_range, dt_train=division[[1]], dt_validation=division[[2]], dt_test=division[[3]], to_exclude)
rpart_run(prefix_range, dt_train=division[[1]], dt_validation=division[[2]], dt_test=division[[3]], to_exclude)
gbm_run(prefix_range, dt_train=division[[1]], dt_validation=division[[2]], dt_test=division[[3]], to_exclude)

# HMM 
# choose features that we use for HMM calculation (sequences)
# features_for_hmm <- c('activity_name_short', 'org.resource', grep('^Attribute.', colnames(division[[1]]), fixed=FALSE, value=TRUE))
hmm_run(prefix_range, states_range, train_data=division[[1]], validation_data=division[[2]], test_data=division[[3]], features_for_hmm, iter)

# HMM + RF
hmm_rf_run(prefix_range, states_range, train_data=division[[1]], validation_data=division[[2]], test_data=division[[3]], features_for_hmm, to_exclude, iter)

# RF with frequencies
rf_with_frequencies_run(prefix_range, dt_train=division[[1]], dt_validation=division[[2]], dt_test=division[[3]], to_exclude)

# HMM + RF with frequencies
hmm_rf_with_frequencies_run(prefix_range, states_range, train_data=division[[1]], validation_data=division[[2]], test_data=division[[3]], features_for_hmm, to_exclude, iter)

# BASELINES
indexes_events_run(prefix_range, dt_train=division[[1]], dt_validation=division[[2]], dt_test=division[[3]], activity_feature="activity_name_short")
frequencies_only(prefix_range, dt_train=division[[1]], dt_validation=division[[2]], dt_test=division[[3]])

events_boolean(prefix_range, dt_train=rbind.data.frame(division[[1]],division[[2]]), dt_validation=division[[2]], dt_test=division[[3]])

indexes_events_latest_data_run(prefix_range, dt_train=division[[1]], dt_validation=division[[2]], to_exclude, 
                               dt_test=division[[3]], activity_feature="activity_name_short")



