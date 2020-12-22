##### FUNCTIONS FOR REPORTING ######

## reports simmary statistic of posterior parameter (mean, SD, 2.5%, 50%, 97.5%)
as.report <- function(variable, na.rm = TRUE, vector.names = TRUE){
  
  q <- quantile(variable, probs = c(.025, .5, .975), na.rm = na.rm)
  
  y <- round(c(mean(variable, na.rm = na.rm),
               sd(variable, na.rm = na.rm),
               q[[1]],
               q[[2]],
               q[[3]]), digits = 2)
  
  if(vector.names == TRUE){
  
  names(y) <- c("Mean", "SD", "2.5%", "50%", "97.5%") 
  
  } else {
    
    NULL
    
  }
  
  return(y)
  
}

## creates new variable which is the rowmeans of the dataframe input (default name is "mean")
apply.mean <- function(dataframe){
  
  result <- apply(dataframe, 1, mean)
 
  return(result)
   
}


## writing the function for bayes.p

## computes bayes p value of posterior (default round to 2 dp)
bayes.p <- function(variable, round.dp = 2){
  
  y <- round((length(variable[variable > 0]))/(length(variable)),
             digits = round.dp)
  
  return(y)
}


warp.bf <- function(brmhypothesis){
  if(brmhypothesis$hypothesis$Evid.Ratio == Inf){
    return(nrow(brmhypothesis$samples)) 
  } else if (brmhypothesis$hypothesis$Evid.Ratio == 0) {
    
    return(1/nrow(brmhypothesis$samples))
  } else {
    
    return(brmhypothesis$hypothesis$Evid.Ratio)
    
  }
  
}

report.bf <- function(bf, sf = 3, dp = 2, N = NULL){
  
  if (is.na(bf)){
    
    return(NaN)
    
  } else if (bf == Inf){
    
    return(N)
    
  } else if (bf == 0){
    
    return(1/N)
    
  } else if (bf<=999 | bf>= .001){
    
    return(signif(bf, digits = sf)) ## report 3 siginificant figures
    
  } else {
    
    return(format(bf, scientific = T, digits = sf))
  } 
  
  
  
}

report.rope <- function(posterior, inside = c(-.1, .1), direction = TRUE, dp = 2){
  
  require(bayestestR)
  
  post.mean <- mean(posterior)
  
  if(direction == FALSE){
    
    result <- rope(posterior, ci = c(0.95), range(-.1, .1))
    return(round(result[[4]][1]/100, dp))
    
  } else if (direction == TRUE & post.mean < 0){
    
    result <- rope(posterior, ci = c(0.95), range(-.1, Inf))
    return(round((100 - result[[4]][1])/100, dp))
    
  } else if (direction ==  TRUE & post.mean > 0){
    
    result <- rope(posterior, ci = c(0.95), range(-Inf, .1))
    return(round((100 - result[[4]][1])/100, dp))
    
  } else {
    
    print("ERROR")
    
  }
  
} 


report.p <- function(p){
  
  if (p >= .01 & p <= 1 & p > 0){
    
    return(paste("=", round(p, 2)))
    
  } else if (p < .01 & p >= .001 & p <= 1 & p > 0) {
    
    return(paste("=", round(p, 3)))
    
  } else if (p < .001 & p <= 1 & p > 0){
    
    return("< .001")
    
  } else if (p > 1 & p > 0) {
    
    print ("ERROR: p cannot be greater than 1")
    
  } else if (p < 0) {
    
    print("ERROR: p cannot be less than 0")
  } else {
    
    print("ERROR: invalid input. p has to be a number between 0 and 1")
    
  }
  
}

extract.p <- function(freq.model, iv = NULL, exact.match = FALSE){
  
  
  freq.model.sum <- summary(freq.model)
  
  if (length(iv) == 0 & nrow(anova(freq.model)) == 2){
    
    return(report.p(anova(freq.model)$`Pr(>F)`[1]))
    
  } else if (length(iv) == 0 & nrow(anova(freq.model)) > 2) {
    
    return(anova(freq.model))
    
  } else if (length(iv) > 0 & exact.match == TRUE) {
    
    iv.index <- colnames(freq.model.sum$cov.unscaled) %in% iv
    
    iv.p <- report.p(freq.model.sum[["coefficients"]][iv.index, "Pr(>|t|)"])
    
    return(iv.p)
    
    
 }  else  {
    
    iv.index <- grep(iv, colnames(freq.model.sum$cov.unscaled), fixed = T)
    
    iv.p <- report.p(freq.model.sum[["coefficients"]][iv.index, "Pr(>|t|)"])
    
    return(iv.p)
    
  }
}



# test <- data.frame(id = 1:10, value = rnorm(10, 1, 2))
# 
# as.report(test$value)
# 
# apply.mean(test)
# 
# as.report(apply.mean(test))
# 
# bayes.p(test$value)



##### FUNCTIONS FOR DOING BAYESIAN #####


## computes precision from sd
to.precise <- function(x){
  y <- 1/(x^2)
  return(y)
}


## computes sd from precision
to.sd <- function(x){
  y <- sqrt((1/x))
  return(y)
}


###### FUNCTIONS FOR RECODING ######

recode.multiple <- function(dataframe, columns, from, to = 1:length(from), fun = as.numeric){
  
  for(i in 1:length(from)){ ## this does the heavy lifting of recoding
    dataframe[, columns][dataframe[, columns] == from[i]] <- to[i]
  }
  
  if (is.null(fun)){ ## this does the coercion when the fun argument is specified
    NULL             ## note that the default argument is as.numeric
  } else{
    dataframe[, columns] <- lapply(dataframe[, columns], fun)
  }
  
  return(dataframe)
}


library(stringr)
## function to clean text
Clean_String <- function(string){
  # Remove everything that is not a number or letter (may want to keep more 
  # stuff in your actual analyses). 
  temp <- string
  
  temp <- stringr::str_replace_all(temp,"[|]", " ")
  # Shrink down to just one white space
  temp <- stringr::str_replace_all(temp,"[\\s]+", " ")
  # Split it
  temp <- stringr::str_split(temp, " ")[[1]]
  # Get rid of trailing "" if necessary
  indexes <- which(temp == "")
  if(length(indexes) > 0){
    temp <- temp[-indexes]
  } 
  return(temp)
}


## keeps creates a column (newname) which keeps count of number of non-missing values 
## in the columns provided in list
keep.count <- function(df, columns, newname){
  df[, newname] <- NA
  count <- df[,columns]
  count.nomiss <- !is.na(count)
  df[, newname] <- apply(count.nomiss, 1, sum)
  
  return(df)
}


collapse.category <- function(dataframe, return.conflicts = TRUE) {
  
  
  ## identifying conflicting columns (NAs non-inclusive)
  df_temp <- dataframe
  df_conflicts <- c() ## initialise vector to store conflicts
  
  for(i in 1:nrow(df_temp)){
    
    df_non_duplicates <- NA # initialize on each row
    
    ## remove duplicates from these rows
    df_non_duplicates <- df_temp[i , ][!duplicated(t(df_temp[i , ]))]
    ## remove missing data
    df_non_duplicates <- df_non_duplicates[!is.na(df_non_duplicates)]
    
    ## are there conflicting values???
    if (length(df_non_duplicates) > 1){ ## assuming there are no conflicts, this length should be 1
      
      df_conflicts <- c(df_conflicts, i)
      
    } else {
      
      NULL
    }
    
    
  }
  
  #print(df_conflicts)
  
  ## reports if there are conflicts
  if (length(df_conflicts) > 0){
    print("WARNING: There are non-NA conflicts in the category columns")
    
  } else {
    
    print("CONGRATS: There are no conflicts in the category columns")
  }
  
  result <- rep(NA, nrow(dataframe))
  
  #For each specified category column
  for (i in 1:nrow(dataframe)){
    for (j in 1:ncol(dataframe)){
      #if the row of the specified column is not missing
      if(!is.na(dataframe[i,j])){
        #then slot that value into the new column
        result[i] <- dataframe[i,j]
      } else {
        #otherwise leave as it is (NA)
        NULL
      }
      
    }
  }
  
  if (length(df_conflicts) > 0 & return.conflicts==TRUE){
    
    print("returning LIST of conflicts. Output object CANNOT be new column of dataframe")
    conflict.dataframe = dataframe[df_conflicts, ]
    conflict.dataframe$conflict.rows <- df_conflicts
    
    return(list(result = result, conflict.dataframe = conflict.dataframe))
    
  } else {
    return(result)
  }
  
}

## note that order.by is probably a list
maintain.order <- function(df, order.by, vars.to.order, new.category.name = "time", exact.match = TRUE, starts.with = TRUE, check = TRUE) {
  
  
  ## dealing with different list structures
  if(!is.list(df[ ,order.by][[1]])) {
    
    order.col <- df[ ,order.by]
    
  } else if(!is.list(df[ ,order.by][[1]][[1]])) {
    
    order.col <- df[ ,order.by][[1]]
    
  } else if(!is.list(df[ ,order.by][[1]][[1]][[1]])) {
    
    order.col <- df[ ,order.by][[1]][[1]]
    
  } else {
    
    print("ERROR: The list structure for order.by is not appropriate")
    
  } ## order.col fixes the structure of the list
  

  if (exact.match == TRUE) { ## if vars.to.order names matches order.by exactly
    
    ## indexing variables to order by list to order by
    list.index <- list()
    for (i in 1:nrow(df)){
      
      temp.index <- which(order.col[[i]] %in% vars.to.order)
      
      list.index <- c(list.index, list(order.col[[i]][temp.index]))
      
    } ## this produces a list of the correct order for each participant
    
    
    order.length <- c()
    for (i in 1:nrow(df)) {
      
      temp <- length(list.index[[i]])
      order.length <- c(order.length, temp)
      
    }
    
  # max.order.length <- max(order.length) ## this gets maximum number of variables across participants
  #   
  #   
  #   
  #   ## initializing new columns
  #   for (i in 1:max.order.length) {
  #     
  #     ## initialize new column
  #     df[, paste(new.category.name, i, sep = ".")] <- NA
  #     
  #   }
    
    
    
    
    for (i in 1:nrow(df)) {
      
      for (j in 1:order.length[i]) {
        
        
        df[i, paste(new.category.name, j, sep = ".")] <- df[i , list.index[[i]][j]]
        
        
      }
      
    }
    
    if (check == TRUE){
      print(list.index) ## returns order list for each participants for checking
    } else {
      NULL
    }
  
  return(df)
    
  } else if (exact.match == FALSE & starts.with == TRUE) { #if order.by variable pattern in vars.to.order
    
    
    ## getting variable names of vars.to.order via list.index
    
    list.index <- list()
    
    for (i in 1:nrow(df)) {
      
      temp.index <- c()
      
      for (j in 1:length(order.col[[i]])) {
        
        index <- grep(paste0("^", order.col[[i]][j]), vars.to.order)   
        
        if (length(index) == 0){
          
          NULL
          
        } else {
          
          temp.index <- c(temp.index, vars.to.order[index]) ## getting ordered variables for each participant
          
        }
        
      }
      
      list.index <- c(list.index, list(temp.index)) ## getting list of order variables
      
    }
    if (check == TRUE){
      print(list.index)
    } else {
      NULL
    }
    
    order.length <- c()
    for (i in 1:nrow(df)) {
      
      temp <- length(list.index[[i]])
      order.length <- c(order.length, temp)
      
    }
    
 #   max.order.length <- max(order.length) ## this gets maximum number of variables across participants
    
    
    for (i in 1:nrow(df)) {
      
      for (j in 1:order.length[i]) {
        
        
        df[i, paste(new.category.name, j, sep = ".")] <- df[i , list.index[[i]][j]]
        
        
      }
      
    }
    
    return(df)
    
    
    
  } else if (exact.match == FALSE) { #if order.by variable pattern in vars.to.order
    

    
    ## getting variable names of vars.to.order via list.index
    
    list.index <- list()
    
    for (i in 1:nrow(df)) {
      
      temp.index <- c()
      
      for (j in 1:length(order.col[[i]])) {
      
      index <- grep(order.col[[i]][j], vars.to.order)   
      
      if (length(index) == 0){
        
        NULL
        
      } else {
        
      temp.index <- c(temp.index, vars.to.order[index]) ## getting ordered variables for each participant
      
      }
      
      }
      
      list.index <- c(list.index, list(temp.index)) ## getting list of order variables
      
    }
    if (check == TRUE){
    print(list.index)
    } else {
      NULL
    }
    
    order.length <- c()
    for (i in 1:nrow(df)) {
      
      temp <- length(list.index[[i]])
      order.length <- c(order.length, temp)
      
    }
    
 #   max.order.length <- max(order.length) ## this gets maximum number of variables across participants
    
    
    for (i in 1:nrow(df)) {
      
      for (j in 1:order.length[i]) {
        
        
        df[i, paste(new.category.name, j, sep = ".")] <- df[i , list.index[[i]][j]]
        
        
      }
      
    }
    
    return(df)
    
    
    
  } else {
    
    NULL
    
  }
  
}


df.Rep <- function(.data_Frame, .search_Columns, .search_Value, .sub_Value){
  .data_Frame[, .search_Columns] <- ifelse(.data_Frame[, .search_Columns]==.search_Value,.sub_Value/.search_Value,1) * .data_Frame[, .search_Columns]
  return(.data_Frame)
} 

##### FUNCTIONS FOR DEALING WITH LARGE DATASETS #####

reduce.listwise.deletion <- function(dataframe, target.columns = colnames(dataframe), drop){
  
  check.dataframe <- is.na(dataframe[, target.columns])
  
  for (i in 1:nrow(check.dataframe)){
    
    if (sum(check.dataframe[i, ]) <=  drop){ 
      
      NULL
      
    } else if (sum(check.dataframe[i, ]) >  drop) {
      
      check.dataframe[i, ] <- FALSE
      
    } else {
      
      print("ERROR")
      
    }
    
  }
  
  sum.dataframe <- apply(check.dataframe, 2, sum) # sum each column

  return(sum.dataframe)
  
}

cor.with <- function(dataframe, var, with = colnames(dataframe),  group = NULL){
  
  require(tidyverse)
  
  if (length(with) != length(colnames(dataframe))) {
    
    with <- c(var, with) ## adds variable for users
    
  } else {
    
    NULL
    
  }
  
  
  if (length(group) == 0 ) {
    
    result <- cor(dataframe[, with], use = "complete.obs") [ , var]
    
    print(paste("n =", nrow(na.omit(dataframe[, with]))))
                
    return(result)
    
  } else {
    
    df.grp <- as.vector(unlist(unique(dataframe[, group])))
    
    cor.matrix <- c()
    
    while (length(df.grp) > 0) { # for each factor in group
      
      subset.rows <- which(dataframe[, group] == df.grp[1])
      ## subset out group
      temp.df <- dataframe[subset.rows, ]
      
      temp.row <- c() ## initializing
      
      temp.row <- cor(temp.df[, with], use = "complete.obs") [ , var]
      
      
      
      ## add temp.row by row
      cor.matrix <- rbind(cor.matrix, temp.row)
      
      ## closing while loop
      df.grp <- df.grp[-1]
      
    }
    
    ## give row names
    rownames(cor.matrix) <- as.vector(unlist(unique(dataframe[, group])))
    
    ## give col names
    colnames(cor.matrix) <- with
    
    print("MESSAGE: Complete Observations are used for each Group")
    
    return(cor.matrix)
    
    
  }
  
  
  
}

regress.with <- function(dataframe, var, with, control = NULL,  group = NULL){
  
  require(tidyverse)
  
  if (length(control) == 0){
    ## formulae
    formulae.str <- paste(var, "~")
    formulae.str <- paste(formulae.str, with)
    
    formulae <- as.formula(formulae.str)  
    
    
  } else {
    
    formulae.str <- paste(var, "~")
    formulae.str <- paste(formulae.str, with)
    
    for (i in 1:length(control)) {
      formulae.str <- paste(formulae.str, "+")
      formulae.str <- paste(formulae.str, control[i])
    }
    
    formulae <- as.formula(formulae.str)  
  }
  
  if (length(group) == 0 ) {
    
    result <- c()
    
    for (i in 1:length(with)) {
      model <- lm(formulae, data = dataframe)
      
      value <- model$coefficients[2]
      result <- c(result, value)
    }
    
    names(result) <- with
    
    return(result)
    
  } else {
    
    df.grp <- as.vector(unlist(unique(dataframe[, group])))
    
    
    cor.matrix <- c()
    
    while (length(df.grp) > 0) { # for each factor in group
      
      subset.rows <- which(dataframe[, group] == df.grp[1])
      ## subset out group
      temp.df <- dataframe[subset.rows, ]
      
      temp.row <- c() ## initializing
      
      for (i in 1:length(with)) {
        model <- lm(formulae, data = temp.df)
        
        value <- model$coefficients[2]
        
        temp.row <- c(temp.row, value)
        
      }
      
      ## add temp.row by row
      cor.matrix <- rbind(cor.matrix, temp.row)
      
      ## closing while loop
      df.grp <- df.grp[-1]
      
    }
    
    # give row names
    rownames(cor.matrix) <- as.vector(unlist(unique(dataframe[, group])))
    
    ## give col names
    colnames(cor.matrix) <- with
    
    return(cor.matrix)
    
    
  }
  
  
  
}

abstract.regression <- function(form, data) {
  
  locate.element <- function(form){ # unless they are within {}
    
    between <- str_locate_all(form, "[[:space:]]|[[+]]|[[-]]|[[*]]|[[~]]")[[1]][,1]
    
    betw.brack <- str_locate_all(form, "(?<=\\{).+?(?=\\})")[[1]]
    
    if (length(betw.brack > 0)){
      
      ignore <- c()
      
      for(i in 1:nrow(betw.brack)) {
        
        ignore <- c(ignore, betw.brack[i, 1]:betw.brack[i, 2])
        
      }
      
      between <- between[!between %in% ignore]
      
    } else {NULL}
    
    length <- 1:nchar(form)
    text <- length[!length %in% between]
    
    return(unname(text))
    
    
  }
  
  locate.operator <- function(form){ # unless they are within {}
    
    between <- str_locate_all(form, "[[+]]|[[-]]|[[*]]|[[~]]")[[1]][,1]
    
    betw.brack <- str_locate_all(form, "(?<=\\{).+?(?=\\})")[[1]]
    
    if (length(betw.brack > 0)){
      
      ignore <- c()
      
      for(i in 1:nrow(betw.brack)) {
        
        ignore <- c(ignore, betw.brack[i, 1]:betw.brack[i, 2])
        
      }
      
      between <- between[!between %in% ignore]
      
    } else {NULL}
    
    return(unname(between))
    
    
  }
  
  break.formulae <- function(form){
    
    require(stringr)
    
    ## break text
    form <- gsub(" ", "", form)
    
    text <- locate.element(form)
    
    group.text <- list()
    
    for (i in 1:length(text)){
      
      if( i == 1 ) { # store first char
        
        group.text[[1]] <- text[1]
        count <- 1
        
      } else if ( text[i] == text[i-1] + 1 ) {
        
        group.text[[count]] <- c(group.text[[count]], text[i])
        
      } else if ( text[i] != text[i-1] + 1 ) {
        
        group.text[[count + 1]] <- text[i]
        count <- count + 1
        
      }
      
      
    }
    
    
    elements <- c()
    
    for (i in 1:length(group.text)){
      
      elements <- c(elements, substr(form, min(group.text[[i]]), max(group.text[[i]])))
      
    }
    
    return (elements )
    
    
    
  }
  
  simple.operation <- function(expr) {
    
    expr <- gsub(" ", "", expr)
    
    operator <- substr(expr, locate.operator(expr), locate.operator(expr)) 
    
    first <- break.formulae(expr)[1]
    second <- break.formulae(expr)[2]
    
    
    if(operator == "+"){
      
      as.numeric(first) + as.numeric(second)
      
    } else if (operator == "-"){
      
      as.numeric(first) - as.numeric(second)
      
    } else if (operator == "*"){
      
      as.numeric(first) * as.numeric(second)
      
    } else{
      
      NULL
    }
    
    
  }
  
  get.variable <- function(form, data){
    
    template <- break.formulae(form)
    
    # change {} to .
    for (i in 1:length(template)){
      
      betw.brack <- str_locate_all(template[i], "(?<=\\{).+?(?=\\})")[[1]]
      
      # replacing {} with .
      if (length(betw.brack) > 0){ 
        
        expr <- paste0("{", substr(template[i], betw.brack[1,1], betw.brack[1, 2]), "}")
        
        template[i] <- gsub(
          pattern = expr,
          replacement = ".", template[[i]], fixed = T)
        
        
      } else {
        
        NULL
        
      }
      
      
    }
    
    # replacing X with .
    for (i in 1:length(template)){
      
      template[[i]] <- gsub("X", ".", template[[i]] )
      
      
    }
    
    # getting rid of redundant templates
    template <- as.vector(unlist(unique(template)))
    
    # initialising list
    variable.list <- vector(mode = "list", length = length(template))
    
    names(variable.list) <- template
    
    require(gtools)
    
    for (i in 1:length(template)){
      
      variable.list[[i]] <- mixedsort(colnames(data)[grep(names(variable.list)[i], colnames(data))])
      
      
    }
    
    return(as.vector(unlist(variable.list)) )
    
  }
  
  
  get.formulae <- function(form, variables, data){
    
    formulae.list <- list()
    
    for ( i in 1:length(variables)){
      
      betw.brack <- str_locate_all(form, "(?<=\\{).+?(?=\\})")[[1]]
      
      # dealing with {}
      if (length(betw.brack) > 0){
        
        expr <- c() # initialise
        
        for (j in 1:nrow(betw.brack)){ # for each {}
          
          expr <- c(expr, substr(form, betw.brack[j,1], betw.brack[j, 2]))
          
        }
        
        expr <- gsub ("X", i, expr)
        
        for (j in  1:length(expr)){ # for each expression
          
          expr[j] <- simple.operation(expr[j] )
          
          
        } 
        
        
        ### sub expr back to {}
        temp.form <- form
        for (j in 1:nrow(betw.brack)){ # for each {}
          
          temp.expr <- paste0("{", substr(form, betw.brack[j,1], betw.brack[j, 2]), "}")
          temp.form <- gsub(temp.expr, expr[j], temp.form, fixed = T)
          
        }
        
        ### sub X to i
        temp.form <- gsub("X", i, temp.form)
        
        temp.form.elements <- as.vector(unlist(break.formulae(temp.form)))
        
        check <- sum(temp.form.elements %in% colnames(data))
        
        if ( check == length(temp.form.elements)){
          
          formulae.list <- c(formulae.list, as.formula(temp.form))
          
        } else {NULL}
        
        
      } else { # if no {}
        
        ### sub X to i
        temp.form <- form
        
        temp.form <- gsub("X", i, temp.form)
        
        temp.form.elements <- as.vector(unlist(break.formulae(temp.form)))
        
        check <- sum(temp.form.elements %in% colnames(data))
        
        if ( check == length(temp.form.elements)){
          
          formulae.list <- c(formulae.list, as.formula(temp.form))
          
        } else {NULL}
        
      }
      
    }
    
    return(formulae.list)
    
  }
  
  #######
  
  elements <- break.formulae(form)
  elements # get constituent parts
  
  variables <- get.variable(form, data) # variables from data
  
  forms <- get.formulae(form, variables, data) # possible formulaes for modelling
  
  ######
  
  model.list <- list()
  
  for (i in 1:length(forms)){ # for each model
    
    model.list[[i]] <- summary(lm(forms[[i]], data = data))
    
    
  }
  
  for (i in 1:length(forms)){ # for each model
    
    if (i == 1) {
      
      coefs <-  model.list[[i]]$coefficients[2:nrow(model.list[[i]]$coefficients), c("Estimate", "Std. Error", "Pr(>|t|)")]
      
    } else {
      
      coefs <- rbind (coefs, 
                      model.list[[i]]$coefficients[2:nrow(model.list[[i]]$coefficients),c("Estimate", "Std. Error", "Pr(>|t|)")])
      
      
    }
    
  }
  
  print(forms) # print formulaes for checking
  return(coefs)
  
}

extract.df.list <- function(input, parameter, max.int, na.num = length(max.int)){ # give parameter number safer
  
  locate.use <- function(coefs, max.int){
    
    if (is.vector(coefs)){
      
      n <- 1
      
    } else {
      
      n <- nrow(coefs)
    }
    
    
    if (n >= max(max.int)){
      
      return(max.int)
      
    } else{ # if less rows than provided interval
      
      while(n < max(max.int)){
        
        max.int <- max.int[-length(max.int)]
        
      }
      
      return(max.int)
      
    }
    
  }
  
  add.na <- function(coefs, na.num){
    
    if (length(coefs) < na.num){
      
      coefs <- c(coefs, rep(NA, na.num - length(coefs) ))
      
    } else {coefs <- coefs}
    
    
  }
  
  output <- NULL
  
  for (i in 1:length(input)) {
    
    if (is.vector(input[[i]])){
      
      temp <- input[[i]][parameter]
      temp <- add.na(temp, na.num)
      
      output <- rbind(output, temp)
      
      
    } else if (i == 1) {
      
      output <- input[[i]][ locate.use(input[[i]], max.int ), parameter  ]
      output <- add.na(output, na.num)
      
    } else {
      
      temp <- input[[i]][ locate.use(input[[i]], max.int ), parameter  ]
      temp <- add.na(temp, na.num)
      
      output <- rbind(output, temp)
      
      
    }
    
  }
  
  return(output)
  
}


##### TIDYVERSE SPECIFIC FUNCTIONS ####

add.count <- function (x, rename = NULL) {
  
  require(tidyverse)
  
  count.data <- x %>%
    count()
  
  if (length(rename) == 0) {
    
    NULL
    
  } else {
    
    colnames(count.data)[grep("n", colnames(count.data)) ] <- rename
    
    
  }
  
  y <- full_join(x, count.data)
  return(y)
  
}


gather.keep <- function(df, gather.col, keep = NULL, key = "key", value= "value") {
  
  require(tidyverse)
  
  if (is.null(keep)) {
    
    result <- gather(df[, gather.col], key = key, value = value)
    return(result)
    
  } else {
    
    result <- gather(df[, gather.col], key = key, value = value) 
    
    no.var <- length(gather.col) 
    
    grow <- df[, c(gather.col, keep)]
    temp <- grow
    for (i in 1:(no.var - 1)){
    temp <- rbind(temp, grow) ##multiplying temp to match rows for result
    
    }
    
    result <- cbind(result, temp[, keep])
    
    colnames(result) <- c(key, value, keep)
    
    return(result)
  }
  
  
}



##### FUNCTIONS FOR SPECIFIC COMPUTATIONS ####

percentile.filter.lower <- function(data, percentile) {
  
  result <- c()
  
  for (i in 1:nrow(data)) {
    
    cutoff <- quantile(unlist(data[i, ]), probs = percentile, na.rm = T)
    
    result <- c(result, mean(unlist(data[i, ])[unlist(data[i, ]) < cutoff], na.rm = T))
    
    
  }
  
  return(result)
  
  
  
}

rowSD <- function(data) {
  
  result <- c()
  
  for (i in 1:nrow(data)){
    
    result <- c(result, sd(data[i, ], na.rm = T))
    
    
    
    
  }
  
  return(result)
  
}

compute.roc <- function(data, positive, negative){
  
  true.positive <- rep(0, nrow(data))
  false.positive <- rep(0, nrow(data))
  false.negative <- rep(0, nrow(data))
  true.negative <- rep(0, nrow(data))
  
  for(i in 1:nrow(data)){
    
    for(j in 1:length(positive)){
      
      if (is.na(data[i, positive[j]])){
        
        NULL
        
      } else if (data[i , positive[j]] == 1) { # if participant respond positive for positive questions
        
        true.positive[i] <-  true.positive[i] + 1 # add count for true positive
        
      } else if (data[i , positive[j]] == 0) { # if participant respond negative for positive questions
        
        false.positive[i] <- false.positive[i] + 1
        
      }
      
    }
    
    for(j in 1:length(negative)){
      
      if(is.na(data[i , negative[j]])){
        
        NULL
        
      } else if (data[i , negative[j]] == 1) { # if participant respond positive for negative questions
        
        false.negative[i] <-  false.negative[i] + 1 # add count for true positive
        
      } else if (data[i , negative[j]] == 0) { # if participant respond negative for negative questions
        
        true.negative[i] <- true.negative[i] + 1
        
      }
      
    }
    
  }
  
  result <- data.frame(true.positive,
                       false.positive,
                       false.negative,
                       true.negative)
  
  return(result)
  
  
}

remove.outlier <- function(data, sd.num = 2, na.rm = T, print.outlier = TRUE){
  
  lower <- mean(data, na.rm = na.rm) - sd.num * sd(data, na.rm = na.rm)
  upper <- mean(data, na.rm = na.rm) + sd.num * sd(data, na.rm = na.rm)
  
  if(print.outlier == T){
    
    print(paste(length(data[data <= lower | data >= upper]), "outlier(s) removed"))
  } else {NULL}
  
  data[data <= lower | data >= upper] <- NA
  return(data)
}
##### FUNCTIONS FOR DATA WRANGLING #######

multi_spread <- function(df, id, key, var){
  require(gtools)
  require(tidyverse)
  
  df.list <- lapply(var, 
                    function(x){
                      df.temp <- spread(cellphone[, c(id, key, x)], key, x)
                      cname <- mixedsort(unique(as.vector(unlist(df[,key]))))
                      df.temp <- df.temp[, c(id, cname)]
                      colnames(df.temp)[2:ncol(df.temp)]  <- paste(x, colnames(df.temp[,cname]), sep= "_")
                      return(df.temp)
                    })
  
  df.output <- df.list[[1]]
  for(i in 1:(length(var)-1)){
    df.output <- full_join(df.output, df.list[[i+1]])
    
  }
  
  return(as.data.frame(df.output))        
  
}