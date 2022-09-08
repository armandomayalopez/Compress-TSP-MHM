###################################################
# path_data: path to original dataset
# path_micro: path to microaggregate dataset
# k: elements in cluster
# sep: csv input file separator
# header: TRUE or FALSE if input file has headers
###################################################
library("sda")

iloss_value <- function(path_data, path_micro, k, sep, header) {
  datos_std<- FALSE
  #Dataset original
  x <- read.table(path_data, header = header, sep = sep)
  x <- as.matrix(x)
  dimnames(x) = list(1:nrow(x),1:ncol(x))
  
  #Dataset microagregado
  y <- read.table(path_micro, header = header, sep = sep)
  y <- as.matrix(y)
  dimnames(y) = list(1:nrow(y),1:ncol(y))

  #y <- y[,-14] #correccion ficheros 
  
  #creamos matriz de grupos
  grupos <- matrix(data = NA, ncol=((2*k)-1), byrow=T)
  grupos <- grupos[-1,]
  
  #buscamos elementos iguales en y para obtener grupos
  
  columnas <- ncol(y)
  
  for(i in 1:nrow(x)){
  
    punto <- as.numeric(unique(which(apply(y, 1, function(q) identical(q[1:columnas], y[i,])))))

    #completamos con NA si el grupo es menor de 2k-1
    while(length(punto) < ((2*k)-1)){
      punto <- c(punto,NA)
    }
    
    grupos <- rbind(grupos,punto)

  }
  
  #preparamos datos para calcular SST y SSE
  grupos <- unique(grupos)
  dimnames(grupos) <- list(1:nrow(grupos),1:ncol(grupos))
  
  datos <- as.matrix(x)
  
  matrix_centroides <- matrix(0,ncol=ncol(datos),byrow=T)
  matrix_centroides <- matrix_centroides[-1,]
  
  if(!datos_std){
    
    datos_standard <- x 
    
    for (i in 1:ncol(x)){
      media <- mean(x[,i])
      desviacion_tipica <- sqrt(var(x[,i]))
      datos_standard[,i] <- (x[,i]-media)/desviacion_tipica 	
    }
    
    dimnames(datos_standard) = list(1:nrow(datos_standard),1:ncol(datos_standard))
    x <- datos_standard
    remove(datos_standard)
    
    datos_standard <- y 
    
    for (i in 1:ncol(y)){
      media <- mean(y[,i])
      desviacion_tipica <- sqrt(var(y[,i]))
      datos_standard[,i] <- (y[,i]-media)/desviacion_tipica 	
    }
    
    dimnames(datos_standard) = list(1:nrow(datos_standard),1:ncol(datos_standard))
    y <- datos_standard
    remove(datos_standard)
    
  }
  
  
  
  
  datos <- as.matrix(x)
  
  #matrix_centroides <- matrix(0,ncol=ncol(datos),byrow=T)
  #matrix_centroides <- matrix_centroides[-1,]
  
  #for(i in 1:nrow(grupos)){
  
  # matrix_centroides <- rbind(matrix_centroides, as.numeric(y[grupos[i,1],]))	
  
  #}
  
  #dimnames(matrix_centroides) = list(1:nrow(matrix_centroides),1:ncol(matrix_centroides))
  
  matrix_centroides <- matrix(0,ncol=ncol(datos),byrow=T)
  matrix_centroides <- matrix_centroides[-1,]
  
  
  contador <- 1
  while(contador<=nrow(grupos)){
    
    particion <- rbind(datos[as.integer(grupos[contador,!is.na(grupos[contador,1:((2*k)-1)])]),])#cogemos los puntos que nos indica grupo_final, entre 2k-1 y k
    particion <- as.matrix(particion)
    centroide_grupo <- as.numeric(centroids(particion,rownames(particion), verbose=F)$means[,"(pooled)"])#calculamos el centrodie   
    matrix_centroides <- rbind(matrix_centroides, centroide_grupo)		
    contador <- (contador+1)
    
  }
  
  
  dimnames(matrix_centroides) = list(1:nrow(matrix_centroides),colnames(datos)) #colocamos los nombres
  
  #Calculamos Centroide del dataset
  datos <- as.matrix(x)
  centroide_dataset <- as.numeric(centroids(datos,rownames(datos),lambda.var=0,lambda.freqs=0, verbose=F)$means[,"(pooled)"])#calculamos el centrodie
  #centroide_dataset <- centroid(datos)
  
  #Inicializamos SST a 0
  sst <- 0
  datos <- as.matrix(x)
  #Calculamos SST
  for(i in 1:nrow(datos)){
    
    sst <-  sst + (datos[i,]-centroide_dataset)%*%(datos[i,]-centroide_dataset)
    
  } 
  
  #print(paste0("SST: ",sst))
  
  #Inicializamos SSE
  sse <- 0	
  datos <- as.matrix(x)
  
  for(i in 1:nrow(grupos)){
    
    for(j in 1:((2*k)-1)){
      
      if(!is.na(grupos[i,j])){
        
        sse <- sse + (datos[grupos[i,j],]- matrix_centroides[i,])%*%(datos[grupos[i,j],]- matrix_centroides[i,])
        
      }
    }	
  }
  
  #print(paste0("SSE: ", sse))
  
  #Calculamos Iloss
  iloss <- round((sse/sst)*100, digits = 6)
  
  #print(paste0("Iloss: ", iloss))
  
  
  
  return(c(k,sse,sst,iloss))
  #return(iloss)
  #return(grupos)


}


grupos_microagregados <- function(source_data, micro_data, k) {
  
  
  #iloss_out -> true:return iloss/false: return groups
  datos_std=FALSE
  
  max_k <- as.integer((2*k)-1)
  #print(max_k)
  #Dataset original
  x <- as.matrix(source_data)
  dimnames(x) = list(1:nrow(x),1:ncol(x))
  
  #Dataset microagregado
  y <- as.matrix(micro_data)
  dimnames(y) = list(1:nrow(y),1:ncol(y))
  

  #creamos matriz de grupos
  grupos <- matrix(data = NA, ncol=((2*k)-1), byrow=T)
  grupos <- grupos[-1,]
  
  #buscamos elementos iguales en y para obtener grupos
  
  columnas <- ncol(y)
  
  for(i in 1:nrow(x)){
    
    punto <- as.numeric(unique(which(apply(y, 1, function(q) identical(q[1:columnas], y[i,])))))
    
    #completamos con NA si el grupo es menor de 2k-1
    if(length(punto) > max_k){
      #print("aqui")
      #print(length(punto))
      #print(punto)
      puntos <- matrix(data = punto, ncol=k, byrow=T)
      relleno <- matrix(data= NA, ncol = (k-1),nrow = dim(puntos)[1], byrow=T)
     
      #print(relleno)
      puntos <- cbind(puntos,relleno)
      #print(puntos)
      grupos <- rbind(grupos,puntos)
      #for(i_puntos in 1:nrow(puntos)){
          
      #  puntos_row <- c(puntos[i_puntos,])
        
       # while(length(puntos_row) < max_k){
        #  puntos_row <- c(puntos_row,NA)
        #}
        
        #print(puntos_row)
        #grupos <- rbind(grupos,puntos_row)
        
      #}
      
     
    }else{
      
      while(length(punto) < max_k){
        punto <- c(punto,NA)
      }
      grupos <- rbind(grupos,punto)
    }
    
    
    
    
    
   
    
  }
  
  #preparamos datos para calcular SST y SSE
  grupos <- unique(grupos)
  dimnames(grupos) <- list(1:nrow(grupos),1:ncol(grupos))
  

  return(grupos)
  

}