##### given a colname from a df from poly() it gives it the complete name

polystr2name = function(polystr, nomi, n){
    vet = unlist(strsplit(polystr, '\\.')[[1]])
    st = ''
    for (i in 1:n){
        if (as.numeric(vet[i]) != 0){
            st = paste0(st, nomi[i], '_',vet[i] ,'_')
        }
    }
return(st)
}

################uses poly but in an interpretable way

polydf = function(df, degree=2, variables=NULL){
    if (!is.null(variables)){
        df = df[,variables]
    }
    mat = as.matrix(df)
    nomi = colnames(df)
    n = length(nomi)

    
    pdf = poly(mat, degree=degree, raw=T)
    n0 = ncol(pdf)
    nomi0 = colnames(pdf)
    
    polystr2name = function(polystr, nomi, n){
        vet = unlist(strsplit(polystr, '\\.')[[1]])
        st = ''
        for (i in 1:n){
            if (as.numeric(vet[i]) != 0){
                st = paste0(st, nomi[i], '_',vet[i] ,'_')
            }
        }
        return(st)
    }
    
    
    new = numeric(n0)
    for (j in 1:n0){
        new[j] = polystr2name(nomi0[j], nomi, n)
    }
    colnames(pdf) = new
    return(pdf)
}




