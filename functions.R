################################################################################################
################################### Naive boostrap resampling ##################################
################################################################################################

bootstrap.sample <- function(data) {
  res <- data[sample(nrow(data), replace = TRUE), ]
  res
}

################################################################################################
########## Descriptive table of categorical covariates as a function of other variable #########
################################################################################################

desc_cat <- function(datos, variables, covariate, siteid, sitios_total = NULL, digits = 2, 
                     pvalue = T, pdigits = 3, n = T){
  
  if(!is.factor(datos[, siteid])){
    datos[, siteid] <- as.factor(datos[, siteid])
  }
  
  variable <- NULL
  ma <- NULL
  N_t <- NULL
  for(j in 1:length(variables)){
    response <- variables[j]
    data <- na.omit(datos[, c(response, covariate, siteid)])
    tabla <- aggregate(as.factor(data[, response]), list(data[, covariate], data[, siteid]),
                       table)
    table_n <- cbind(tabla[, 1:2], tabla$x)
    names(table_n)[3:dim(table_n)[2]] <- paste(colnames(tabla)[3], 
                                               names(table_n)[3:dim(table_n)[2]], sep = ".")
    rm(tabla)
    table_p <- aggregate(as.factor(data[, response]), list(data[, covariate], data[, siteid]),
                         function(x){
                           prop.table(table(x))
                         })
    sitios <- levels(data[, siteid])
    cat <- length(unique(na.omit(as.numeric(data[, response]))))
    
    if(pvalue){
      P <- as.numeric()
      if(cat == 2){
        for(i in 1:length(sitios)){
          test <- try(fisher.test(data[data[, siteid] == sitios[i], response],
                                  data[data[, siteid] == sitios[i], covariate]), silent = T)
          if(class(test) == "try-error"){
            P[i] <- NA
          }else{
            P[i] <- format(round(test$p.value, pdigits), nsmall = pdigits)
          }
        }
      }else{
        for(i in 1:length(sitios)){
          test <- try(chisq.test(data[data[, siteid] == sitios[i], response],
                                 data[data[, siteid] == sitios[i], covariate], correct = T),
                      silent = T)
          if(class(test) == "try-error"){
            P[i] <- NA
          }else{
            P[i] <- format(round(test$p.value, pdigits), nsmall = pdigits)
          }
        }
      }
    }
    
    if(is.null(sitios_total)){
      sitios_total <- sitios
    }
    matriz <- matrix("--", ncol = 2*length(sitios_total), nrow = cat)
    res <- cbind(table_n, table_p$x)
    
    Group.2 <- unique(res$Group.2)
    for(i in length(Group.2):1){
      if(dim(res[res$Group.2 == Group.2[i], ])[1] == 1){
        row <- as.numeric(row.names(res[res$Group.2 == Group.2[i], ]))
        fila <- c(1, Group.2[i], rep(0 , 2*cat))
        if(row == dim(res)[1]){
          res <- rbind(res[1:row, ], fila)
        }else{
          res <- rbind(res[1:row, ], fila, res[(row + 1):dim(res)[1], ])
        }
      }
    }
    
    l <- length(sitios_total)
    pos <- NULL
    while(length(sitios) != length(sitios_total)){
      posi <- sum(sitios == sitios_total)
      if(posi == (l - 1)){
        sitios <- c(sitios, sitios_total[l])
        fila1 <- c(0,sitios_total[l], rep(0 , 2*cat))
        fila2 <- c(1, sitios_total[l], rep(0 , 2*cat))
        res <- rbind(res, fila1, fila2)
        P <- c(P, NA)
      }else{
        sitios <- c(sitios[1:posi], sitios_total[posi + 1], sitios[(posi + 1):length(sitios)])
        fila1 <- c(0, sitios_total[posi + 1], rep(0 , 2*cat))
        fila2 <- c(1, sitios_total[posi + 1], rep(0 , 2*cat))
        res <- rbind(res[1:(2*posi), ], fila1, fila2, res[(2*posi+ 1):dim(res)[1],])
        P <- c(P[1:posi], NA, P[(posi +1):length(P)])
      }
      pos <- c(pos, posi)
    }
    
    if(!n){
      N <- numeric()
    }
    for(i in 1:length(sitios)){
      for(k in 1:cat){
        if(n){
          matriz[k, (2*i - 1):(2*i)] <- 
            apply(res[res$Group.2 == sitios[i], ][, -(1:2)], 1, 
                  function(x){
                    paste0(x[k], "/", sum(x[1:(cat)]), " (",
                           format(round(100*x[cat + k], digits), nsmall = digits), "%", ")")
                  })
        }else{
          matriz[k, (2*i - 1):(2*i)] <- 
            apply(res[res$Group.2 == sitios[i], ][, -(1:2)], 1, 
                  function(x){
                    paste0(x[k], " (", format(round(100*x[cat + k], digits), nsmall = digits), 
                           "%", ")")
                  })
          N[(2*i - 1):(2*i)] <- apply(res[res$Group.2 == sitios[i], ][, -(1:2)], 1, 
                                         function(x){
                                           paste0("(n = ", sum(x[1:(cat)]), ")")
                                         })
        }
      }
    }
    
    variable <- c(variable,
                  paste0(response, "_", sort(unique(data[, response]))))
    if(pvalue){
      matriz <- rbind(matriz, paste0("P = ", rep(P, each = 2)))
      variable <- c(variable, " ")
    }
    ma <- rbind(ma, matriz)
    if(!n){
      N_t <- rbind(N_t, N)
    }
  }
  
  ma <- cbind(variable, ma)
  
  if(n){
    return(ma)
  }else{
    return(list(ma, N_t))
  }
}

################################################################################################
########## Descriptive table of continuous covariates as a function of other variable ##########
################################################################################################

desc_cont <- function(variables, covariate, siteid, datos, median = T, digits = 1, pvalue = T, 
                      digits_p = 3, n = T){
  
  m <- NULL
  N_t <- NULL
  
  for(j in 1:length(variables)){
    response <- variables[j]
    data <- na.omit(datos[, c(response, covariate, siteid)])
    
    if(median == T){
      table_q <- aggregate(data[, response], list(data[, covariate], data[, siteid]), 
                           quantile, na.rm = T, prob = c(0.50, 0.25, 0.75))
    }else{
      table_q1 <- aggregate(data[, response], list(data[, covariate], data[, siteid]), 
                            mean, na.rm = T)
      table_q2 <- aggregate(data[, response], list(data[, covariate], data[, siteid]), 
                            sd, na.rm = T)
      table_q <- cbind(table_q1, table_q2[, 3])
    }
    
    table_n <- aggregate(data[, response], list(data[, covariate], data[, siteid]), length)
    
    table_q <- cbind(table_q, table_n[, 3])
    
    sitios <- sort(unique(data[, siteid]))
    
    if(pvalue == T){
      P <- as.numeric()
      for(i in 1:length(sitios)){
        test <- try(wilcox.test(data[data[, siteid] == sitios[i], response] ~
                                  data[data[, siteid] == sitios[i],covariate]), silent = T)
        if(class(test) == "try-error"){
          P[i] = NA
        }else{
          P[i] = format(round(test$p.value, digits_p), nsmall = digits_p)
        }
      }
    }
    
    matriz <- matrix("--", ncol = 2*length(sitios), nrow = 1)
    
    Group.2 <-  unique(table_q$Group.2)
    for(i in length(Group.2):1){
      if(dim(table_q[table_q$Group.2 == Group.2[i], ])[1] == 1){
        row <- as.numeric(row.names(table_q[table_q$Group.2 == Group.2[i], ]))
        fila <- c(1, Group.2[i], rep("NA" , 3))
        if(row == dim(table_q)[1]){
          table_q <- rbind(table_q[1:row, ], fila)
        }else{ 
          table_q <- rbind(table_q[1:row, ], fila, table_q[(row + 1):dim(table_q)[1],])
        }
      }
    }
    
    level <- sort(unique(table_q$Group.1))
    if(n == F){
      N <- numeric(dim(table_q)[1])
    }
    for(i in 1:length(sitios)){
      if(n == T){
        matriz[, (2*i - 1):(2*i)] <- 
          c(apply(table_q[table_q$Group.2 == sitios[i] & table_q$Group.1 == level[1],], 1, 
                  function(x){
                    if(median == T){
                      paste0(format(round(as.numeric(x[3]), digits), nsmall = digits), " [",
                             format(round(as.numeric(x[4]), digits), nsmall = digits), "-",
                             format(round(as.numeric(x[5]), digits), nsmall = digits), "] \n", 
                             "(n = ", as.numeric(x[6]), ")")
                    }else{
                      paste0(format(round(as.numeric(x[3]), digits), nsmall = digits), " (",
                             format(round(as.numeric(x[4]), digits), nsmall = digits), ") \n", 
                             "(n = ", as.numeric(x[5]), ")")
                    }
                  }), 
            apply(table_q[table_q$Group.2 == sitios[i] & table_q$Group.1 == level[2],], 1, 
                  function(x){
                    if(median == T){
                      paste0(format(round(as.numeric(x[3]), digits), nsmall = digits), " [",
                             format(round(as.numeric(x[4]), digits), nsmall = digits), "-", 
                             format(round(as.numeric(x[5]), digits), nsmall = digits), "] \n", 
                             "(n = ", as.numeric(x[6]), ")")
                    }else{
                      paste0(format(round(as.numeric(x[3]), digits), nsmall = digits), " (",
                             format(round(as.numeric(x[4]), digits), nsmall = digits), ") \n", 
                             "(n = ", as.numeric(x[5]), ")")
                    }
                  }))
      }else{
        matriz[, (2*i - 1):(2*i)] <- 
          c(apply(table_q[table_q$Group.2 == sitios[i] & table_q$Group.1 == level[1],], 1, 
                  function(x){
                    if(median == T){
                      paste0(format(round(as.numeric(x[3]), digits), nsmall = digits), " [",
                             format(round(as.numeric(x[4]), digits), nsmall = digits), "-",
                             format(round(as.numeric(x[5]), digits), nsmall = digits),"]")
                    }else{
                      paste0(format(round(as.numeric(x[3]), digits), nsmall = digits), " (",
                             format(round(as.numeric(x[4]), digits), nsmall = digits), ")")
                    }
                  }), 
            apply(table_q[table_q$Group.2 == sitios[i] & table_q$Group.1 == level[2],], 1, 
                  function(x){
                    if(median == T){
                      paste0(format(round(as.numeric(x[3]), digits), nsmall = digits), " [", 
                             format(round(as.numeric(x[4]), digits), nsmall = digits), "-", 
                             format(round(as.numeric(x[5]), digits), nsmall = digits), "]")
                    }else{
                      paste0(format(round(as.numeric(x[3]), digits), nsmall = digits), " (", 
                             format(round(as.numeric(x[4]), digits), nsmall = digits), ")")
                    }
                  }))
        N[(2*i - 1):(2*i)] <- 
          c(apply(table_q[table_q$Group.2 == sitios[i] & table_q$Group.1 == level[1],], 1, 
                  function(x){
                    if(median == T){
                      paste0("(n = ", as.numeric(x[6]), ")")
                    }else{
                      paste0("(n = ", as.numeric(x[5]), ")")
                    }
                  }), 
            apply(table_q[table_q$Group.2 == sitios[i] & table_q$Group.1 == level[2],], 1, 
                  function(x){
                    if(median == T){
                      paste0("(n = ", as.numeric(x[6]), ")")
                    }else{
                      paste0("(n = ", as.numeric(x[5]), ")")
                    }
                  }))
      }
    }
    
    for(i in 1:dim(matriz)[2]){
      if(substr(matriz[1, i], 1, 2) == "NA"){
        matriz[1, i] <- "--"
      }
    }
    
    if(pvalue == T){
      ma <- rbind(matriz, paste0("P = ", rep(P, each = 2)))
      ma <- cbind(c(response, " "), ma)
    }else{
      ma <- cbind(response, matriz)
    }
    m <- rbind(m, ma)
    if(!n){
      N_t <- rbind(N_t, N)
    }
  }
  
  if(n == T){
    return(m)
  }else{
    return(list(m, N_t))
  }
}

################################################################################################
######################################## Interpolation  ########################################
################################################################################################

interpolate <- function( x, y, target ){
  ooo <- order(x)
  x <- x[ooo]
  y <- y[ooo]
  n <- length(target)
  m <- length(x)
  iii <- c( 1:length(x) )
  output <- rep( NA, n )
  for( j in 1:n ){
    index.low <- max( iii[ x <= target[j] ] )
    index.high <- min( iii[  target[j] <= x ] )
    if( target[j] > max(x) ){
      index.low <- m
      index.high <- m
    }
    if( target[j] < min(x) ){
      index.low <- 1
      index.high <- 1
    }
    #
    if( index.low == index.high ){
      output[j] <- y[ index.low ]
    }else{
      dx <- ( target[j] - x[index.low] )/( x[index.high]-x[index.low] )
      dy <- ( y[ index.high ] - y [ index.low ] )
      output[j] <- y[ index.low ] + dx * dy
    }
  }
  output
}

################################################################################################
################################ ggplot for continuous variables ###############################
################################################################################################

plot_cont <- function(datos, var = NULL, nvar = NULL, hist = F, labs = NULL, breaks = NULL, 
                      xlim = NULL, ylim = NULL, english = F){
  # datos is the data.frame with continuous covariates
  # var is the vector with the grouping variable. For instance, var = datos$var
  # nvar is the name that we want to use for the grouping variable
  # hist is only available when there is not grouping variable
  # labs, breaks, xlim and ylim are lists and they are needed for the histograms 
  # labs are the label of the represented variables
  
  require("ggplot2")
  a <- dim(datos)[2]
  nm <- names(datos)
  plots <- list()  # new empty list
  if(is.null(var) == F){
    require("reshape2")
    my_data <- melt(datos)
    my_data[, 3]<- rep(var, times = a)
    names(my_data)[3] <- nvar
    p <- ggplot(my_data, aes_string(fill = nvar)) + 
      geom_boxplot(aes(x = variable, y = value)) +
      facet_wrap(~ variable, scales = "free") + theme_gray() +
      stat_summary(geom = "point", fun = "mean", col = "blue", 
                   position = position_dodge(width = 0.75), aes(x = variable, y = value)) 
    plots[[1]] <- p + scale_color_brewer(palette = "Dark2") +
      theme_bw()
  }else{
    if (hist == T){
      for(k in 1:a){
        plots[[(2*k-1)]] <- ggplot(data = datos, aes_string(nm[k])) +
          geom_histogram(breaks = breaks[[k]], col = "black", fill = "white", 
                         aes(y = ..density..)) +
          labs(x = labs[[k]], y = ifelse(english == T, "Density", "Densidad")) + 
          ylim(ylim[[k]]) + geom_density(colour = "blue") +
          theme_bw()
        plots[[(2*k)]] <- ggplot(datos, aes_string(y = nm[k])) + 
          geom_boxplot(aes(x = "")) + theme_gray() + ylim(xlim[[k]]) + 
          stat_summary(geom = "point", fun = "mean", col = "blue",
                       position = position_dodge(width = 0.75), aes(x = "")) + 
          scale_color_brewer(palette = "Dark2") + labs(y = labs[[k]], x = "") +
          theme_bw()
      }
    }else{
      require("reshape2")
      my_data <- melt(datos)
      p <- ggplot(my_data, aes(x = variable, y = value)) + geom_boxplot() +
        facet_wrap(~ variable, scales = "free") + theme_gray() +
        stat_summary(geom = "point", fun = "mean", col = "blue",
                     position = position_dodge(width = 0.75)) 
      plots[[1]] <- p + scale_color_brewer(palette = "Dark2")+
        theme_bw()
    }
  } 
  return(plots)
}
