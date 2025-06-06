clamp = function(vet, min, max){
    vet = ifelse(vet<min, min, vet)
    vet = ifelse(vet>max, max, vet)
    return(vet)
}

invert_dict = function(vet){
    if (is.null(names(vet))){
        names(vet) = 1:length(vet)
    }
    new = names(vet)
    names(new) = vet
    return(new)
}

extract_regex = function(pattern, string){
    return(regmatches(string, gregexec(pattern, string))[[1]][1,])
}
