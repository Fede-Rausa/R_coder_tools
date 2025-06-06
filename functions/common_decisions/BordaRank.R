borda_method = function(score_mat){
    #score_mat is a matrix or a dataframe
    #on the rows you have different items, 
    #on the columns their scores, all better when higher
    C = ncol(score_mat)
    N = nrow(score_mat)
    r = numeric(N)
    for (i in 1:C){
        r = r + rank(score_mat[,i])
    }
    r = r/C #get the mean rank
    
    if (!is.null(rownames(score_mat))){
        names(r) = rownames(score_mat)
    }   
    return(r[order(r, decreasing=T)])
}


