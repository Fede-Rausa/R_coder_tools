
##find the idexes of the columns of a data frame
dummy_id = function(df){
  #returns the indexes of the columns of df that are dummy variables
  C = ncol(df)
  cond = numeric(C)
  for (i in 1:C){
    col = df[,i]
    n = length(unique(col[!is.na(col)]))
    cond[i] = (n==2)
  }
  return(which(as.logical(cond)))
}


## similar to the select function of dplyr
## varnames are the names of the columns that should be extracted
## if exclude is TRUE then is returned the dataframe withouth that columns
select_cols = function(df, varnames, exclude=F){
  if (exclude){
    bool = !(colnames(df) %in% varnames)
  }else{
    bool = (colnames(df) %in% varnames)
  }
  return(df[,bool])
}



##get the types of the columns of a data frame
col_types = function(df){
  df = as.data.frame(df)
  C = ncol(df)
  mat = matrix(NA, ncol=2, nrow=C)
  for (i in 1:C){
    mat[i,1] = colnames(df)[i]
    mat[i,2] = class(df[,i])
  }
  colnames(mat) = c('var', 'class')
  return(mat)
}


#gives a table of descriptives for each column of the data frame
describe_data = function(df, all_info=F, max_unique=10){
  df = as.data.frame(df)
  C = ncol(df)
  N = nrow(df)
  mat_type = col_types(df)
  filtro_num = (mat_type[,2] == 'numeric') 
  
  
  n_num = sum(filtro_num)
  n_char = C - n_num
  
  num_df = df[, filtro_num]
  char_df = df[, !filtro_num]
  
  mat_num = matrix(NA, ncol = 10, nrow = n_num)
  for (i in 1:n_num){
    na_count = sum(is.na(num_df[,i]))
    na_ratio = paste0(round(na_count/N * 100, 2), '%')
    Var = round(var(num_df[,i], na.rm=T), 2)
    sdev = round(sqrt(Var), 2)
    su = round(unlist(summary(num_df[,i])[1:6]), 2)
    su[4] = round(su[4], 2)
    mat_num[i,] = c(su ,Var, sdev , na_count, na_ratio)
  }
  rownames(mat_num) = colnames(num_df)
  colnames(mat_num) = c(names(summary(1:5)),'var', 'sdev' , 'NA count', 'NA ratio')
  
  if (n_char != 0){
    mat_char = matrix(NA, ncol=4, nrow=n_char)
    for (i in 1:n_char){
      na_count = sum(is.na(char_df[,i]))
      na_ratio = paste0(round(na_count/N * 100, 2), '%')
      uni_count = length(unique(char_df[,i]))
      uni_ratio = paste0(round(uni_count/N * 100, 2), '%')
      mat_char[i,] = c(uni_count, uni_ratio, na_count, na_ratio)
    }
    rownames(mat_char) = colnames(char_df)
    colnames(mat_char) = c('unique values', 'unique ratio', 'NA count', 'NA ratio')
    
    if (all_info){
      uni_list = list()
      id = 1
      bools_id = (as.numeric(mat_char[,1]) <= max_unique)
      for (i in 1:n_char){
        if (bools_id[i]){
          var = char_df[,i]
          uni_list[[id]] = table(var[!is.na(var)])
          id = id + 1
        }
      }
      names(uni_list) = rownames(mat_char)[bools_id]
      return(list(dim=dim(df), cn = colnames(df), num=mat_num, cat=mat_char, uni=uni_list))
    }else{
      return(list(num=mat_num, cat=mat_char))
    }
    
  }else{
    return(list(num=mat_num))
  }
  
}
