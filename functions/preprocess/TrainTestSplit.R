
train_test_split = function(df, tr_prp, yname=NULL, seed=NULL){
#tr_prp is a value between 0 and 1 that is 
#the proportion of rows of df that goes to the training set

#yname is the name of the target variable
#if not provide, then the train test split will not be stratified
  
  if (!is.null(seed)){
    set.seed(seed)
  }
  N = nrow(df)
  
  stratified=FALSE
  if (!is.null(yname)){
    library(caret)
    stratified = TRUE
  }
  
  if (stratified){
    y = unlist(df[,yname])
    id_train = caret::createDataPartition(y, 1, p = tr_prp)[[1]]
  }else{
    id_train = sample(1:N, floor(tr_prp*N))
  }
  
  id_test = (1:N)[-id_train]
  
  train = df[id_train,]
  test = df[-id_train,]
  
  return(list(train=train, test=test, id_train=id_train, id_test = id_test))
}
