# Simple script to perform cross-validation with kin.blup in rrBLUP

# Note: the '...' means you can pass more parameters, which will be passed on to kin.blup() unchanged
cross.validate=function(data, pheno, fold=10, seed=1, ...){
  cat("Performing",fold,"fold cross-validation on",pheno,"with seed",seed,"\n")
  set.seed(seed)
  sets = cut(seq(1:nrow(data)), breaks = fold, labels = F)
  #sets = rep(1:fold, nrow(data)/fold + 1)	# make a list of sets to divide data into
  #sets = sets[1:nrow(data)] # Trim off excess
  sets =  sample(sets, replace=F)	# Randomize
  
  # Go through and calculate
  accuracy=list()
  for(s in sort(unique(sets))){
	tomask=sets==s                  # Identify which values to mask
	mydata=data                     # Copy values over
	mydata$MASKED = data[[pheno]]   # Make new phenotype so don't touch original	
	mydata$MASKED[tomask] = NA      # Mask phenotypes
	
	# Make model
	model = kin.blup(data=mydata, pheno="MASKED", ...)
	# Tabulate accuracy
	accuracy[[s]] = cor(data[[pheno]][tomask], model$pred[tomask])
  }
  
  # Convert list back to vector and return
  accuracy=unlist(accuracy)	
  return(accuracy)
}


multi.validate=function(times, seed=1, ...){
  cat("Peforming",times,"cross-validations\n")
  accuracies = sapply(1:times, function(x){
	cors=cross.validate(seed=x, ...)
	return(mean(cors))
  })
  return(accuracies)
}  