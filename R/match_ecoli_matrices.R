# Takes outputs from get_ecoli_matrix for 2014 and 2015 to create matrices
# that only contain shared strains
match_ecoli_matrices = function(matrix.2014, matrix.2015){
  
  matrix.2014.sub <- matrix.2014[rownames(matrix.2014) %in% rownames(matrix.2015),
                           colnames(matrix.2014) %in% colnames(matrix.2015)]
  
  matrix.2015.sub <- matrix.2015[rownames(matrix.2015) %in% rownames(matrix.2014),
                           colnames(matrix.2015) %in% colnames(matrix.2014)]
  
  return(list(matrix.2014.sub = matrix.2014.sub,
         matrix.2015.sub = matrix.2015.sub))
}
