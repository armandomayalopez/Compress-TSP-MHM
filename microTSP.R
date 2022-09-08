# Declaramos librerias

library("TSP")
library("sda")
library("igraph")
library("stringr")
concorde_path("/path/to/concorde")


####################################
############Fuciones################
########################################################################

obtener_sse <- function(x, y, lista_ordenada, datos_std) {
  
  grupo <- c((x+1):y)
  
  #print(grupo)
  
  elementos <- lista_ordenada[grupo[1:length(grupo)]]
  
  #print(elementos)
  
  datos_grupo <- as.matrix(datos_std[elementos,])
  
  #print(datos_grupo)
  
  centroide_dataset <- as.numeric(centroids(datos_grupo,rownames(datos_grupo), verbose=F)$means[,"(pooled)"])#calculamos el centrodie
  
  #Inicializamos SSE
  sse <- 0	
  #datos <- as.matrix(x)
  
  for(i in 1:length(elementos)){
    
    sse <- sse + (datos_grupo[i,]- centroide_dataset)%*%(datos_grupo[i,]- centroide_dataset)
    
  }
  
  return(sse)
  
  
}

########################################################################

obtener_grupos <- function(entrada, lista_ordenada, k_value) {
  
  
  grupos = matrix(NA,nrow(entrada),2*k_value - 1)
  #View(grupos)
  for(i in 1:nrow(entrada)){
    
    grupo <- c((entrada[i,1]+1):entrada[i,2])
    #print(grupo)
    elementos <- lista_ordenada[grupo[1:length(grupo)]]
    
    while(length(elementos) < ((2*k_value)-1)){
      
      elementos <- c(elementos,NA)
      
      
    }
    #print(elementos)
    grupos[i,] <- elementos
    
  }
  
  
  return(grupos)
  #print(grupo)
  
  #elementos <- lista_ordenada[grupo[1:length(grupo)]]
  
  
}

########################################################################

microConcorde <- function(nombre_fichero_datos,k, metodo, separador, cabeceras){
  #for(k in 5:5){
  
  
  #Dataset de entrada
  separador_fichero_in = separador #fichero datos
  cabeceras_fichero_in = cabeceras #fichero datos
  
  #Leemos datos
  datos <- read.table(nombre_fichero_datos, header = cabeceras_fichero_in, sep = separador_fichero_in)
  datos_standard <- datos 
  
  #Iniciamos contador
  time0 <- Sys.time() 
  #Estandarizamos
  for (i in 1:ncol(datos)){
    media <- mean(datos[,i])
    desviacion_tipica <- sqrt(var(datos[,i]))
    datos_standard[,i] <- (datos[,i]-media)/desviacion_tipica 	
  }
  
  remove(datos)
  dimnames(datos_standard) = list(1:nrow(datos_standard),1:ncol(datos_standard))
  
  distancias<-dist(datos_standard)
  
  if(metodo == "concorde"){
    
    tsp <- as.TSP(distancias)
    tsp_dummy <- insert_dummy(tsp, label = "dummy")
    tour_solution <- solve_TSP(tsp_dummy, method = "concorde", control = list(precision = 6, verbose = FALSE))
    remove(tsp)
    
  }else{
    atsp <- as.ATSP(distancias)
    tsp_dummy <- insert_dummy(atsp, label = "dummy")
    tour_solution <- solve_TSP(tsp_dummy, method = metodo)
    remove(atsp)
  }
  
  path_solution <- cut_tour(tour_solution, cut = "dummy")
  elementos_ordenados <- as.numeric(labels(path_solution))
  
  
  ##########################################
  ################## MHM ###################
  ##########################################
  
  numRecords <- nrow(datos_standard)

  #elementos_ordenados <- as.numeric(labels(path_solution))

  #Creamos combination      
  a <- t(combn(c(0:numRecords), 2))
  
  b <- matrix(0,(numRecords + 1),(numRecords + 1)) #en R las matrices empiezan en 1
  
  dimnames(b) = list(1:nrow(b), 1:ncol(b))
  
  
  
  #seleccionamos pares entre k,2*k-1
  flag <- 0
  for(i in 1:nrow(a)){
    
    if( ( (a[i,1] + k ) <= a[i,2]) && ( a[i,2] < (a[i,1] + 2*k)) ) {
      
      if(flag == 0){
        
        out <- rbind(a[i,])
        flag <- 1
        
        
      }else{
        
        out <- rbind(out, a[i,])
      }
      
      
    }
    
  }
  
  #eliminamos a, comentar si se quiere mantener la tabla intermedia  
  remove(a)
  
  #añadimos los arcos segun MHM
  for(i in 1:nrow(out)){
    
    b[out[i,1]+1,out[i,2]+1] <- obtener_sse(out[i,1], out[i,2], elementos_ordenados, datos_standard) # el más 1 R comienza las matrices en 1
    
  }  
  
  #convertimos la matriz b en un grafo pesado, calculamos el camino mas corto
  grafo <- graph_from_adjacency_matrix(b, weighted=TRUE)
  path <- get.shortest.paths(grafo,1,(numRecords + 1)) # el temita
  camino <- as.numeric(path$vpath[[1]])
  
  #R no tiene indice 0 en las columnas, realizamos arreglo
  camino <- camino - 1
  
  #Convertimos el camino de un dataframe de salida a una lista de elementos
  for(i in 1:(length(camino)-1)){
    if(i == 1){  
      salida <- c(camino[i],camino[i+1])
    }else{
      
      salida <- rbind(salida, c(camino[i],camino[i+1]))
      
    }  
  }
  
  salida <- as.matrix(salida)
  dimnames(salida) = list(1:nrow(salida),1:ncol(salida))
  grupos_salida <- obtener_grupos(salida, elementos_ordenados, k)
  remove(salida)
  
  #preparamos datos para calcular SST y SSE
  grupos <- unique(grupos_salida)
  dimnames(grupos) <- list(1:nrow(grupos),1:ncol(grupos))
  
  
  
  #preparamos datos para calcular SST y SSE
  grupos <- unique(grupos_salida)
  dimnames(grupos) <- list(1:nrow(grupos),1:ncol(grupos))
  
  #write.table(grupos, file = paste0(root_path, nombre_dataset,"_resultados/", nombre_dataset, "_gruposk_",k,".csv") , sep = ";", row.names = F, col.names = F)
  
  #Creamos matriz de centroides de los grupos
  
  datos <- as.matrix(datos_standard)
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
  
  dimnames(matrix_centroides) = list(1:nrow(matrix_centroides),1:ncol(matrix_centroides))
  
  #Calculamos Centroide del dataset
  datos <- as.matrix(datos_standard)
  centroide_dataset <- as.numeric(centroids(datos,rownames(datos), verbose=F)$means[,"(pooled)"])#calculamos el centrodie
  
  #Detenemos contador
  time1 <- Sys.time() 
  
  #Inicializamos SST a 0
  sst <- 0
  #datos <- as.matrix(datos_standard)
  #Calculamos SST
  for(i in 1:nrow(datos)){
    
    sst <-  sst + (datos[i,]-centroide_dataset)%*%(datos[i,]-centroide_dataset)
    
  } 
  
  
  #Inicializamos SSE
  sse <- 0	
  #datos <- as.matrix(datos_standard)
  
  for(i in 1:nrow(grupos)){
    
    for(j in 1:((2*k)-1)){
      
      if(!is.na(grupos[i,j])){
        
        sse <- sse + (datos[grupos[i,j],]- matrix_centroides[i,])%*%(datos[grupos[i,j],]- matrix_centroides[i,])
        
      }
    }	
  }
  
  
  
  #Calculamos Iloss
  iloss <- round((sse/sst)*100, digits = 4)
  sse <- round(sse, digits = 4)
  sst <- round(sst, digits = 4)
  
  #Compute time
  
  tiempo <- round(as.numeric(difftime(time1, time0, units = "secs")), digits = 2)
  print(paste0("MicroTSP --> K_value: ", k, " SSE: ", sse , " SST: ", sst , " Iloss: ", iloss, " Tiempo: ", tiempo))

  retornar <- rbind(c(k,sse,sst,iloss,tiempo))
  retornar <- as.matrix(retornar)
  dimnames(retornar) = list(1:nrow(retornar),1:ncol(retornar))
      
  return(retornar)
  
  
  
  
}# fin funcion

########################################################################

microConcordeAll <- function(nombre_fichero_datos, metodo, separador, cabeceras){
  
  #Dataset de entrada
  separador_fichero_in = separador #fichero datos
  cabeceras_fichero_in = cabeceras #fichero datos
  
  #Leemos datos
  datos <- read.table(nombre_fichero_datos, header = cabeceras_fichero_in, sep = separador_fichero_in)
  datos_standard <- datos 
  

  #Iniciamos contador
  timea <- Sys.time() 
  #Estandarizamos
  for (i in 1:ncol(datos)){
    media <- mean(datos[,i])
    desviacion_tipica <- sqrt(var(datos[,i]))
    datos_standard[,i] <- (datos[,i]-media)/desviacion_tipica 	
  }
  
  remove(datos)
  
  dimnames(datos_standard) = list(1:nrow(datos_standard),1:ncol(datos_standard))

  
  distancias_camino <- as.matrix(dist(as.matrix(datos_standard)))
  


  distancias<-dist(datos_standard)

  
  if(metodo == "concorde"){
    
    tsp <- as.TSP(distancias)
    tsp_dummy <- insert_dummy(tsp, label = "dummy")
    tour_solution <- solve_TSP(tsp_dummy, method = "concorde", control = list(precision = 6, verbose = FALSE))
    remove(tsp)
    
  }else{
    atsp <- as.ATSP(distancias)
    tsp_dummy <- insert_dummy(atsp, label = "dummy")
    tour_solution <- solve_TSP(tsp_dummy, method = metodo)
    remove(atsp)
  }


  path_solution <- cut_tour(tour_solution, cut = "dummy")
  
  elementos_ordenados <- as.numeric(labels(path_solution))
  timeb <- Sys.time() 

  #Longitud camino
  camino <- elementos_ordenados
  
  longitud <- 0
  
  
  for(i in 1:(length(camino)-1)) {
    
    longitud <- longitud + as.numeric(distancias_camino[as.numeric(camino[i]),as.numeric(camino[i+1])])
    
  }  
    
  ##########################################
  ################## MHM ###################
  ##########################################
  
  numRecords <- nrow(datos_standard)
  flag_retornar <- 0
  #elementos_ordenados <- as.numeric(labels(path_solution))
  for(k in 3:10){
    time0 <- Sys.time() 
    #Creamos combination      
    a <- t(combn(c(0:numRecords), 2))
    
    b <- matrix(0,(numRecords + 1),(numRecords + 1)) #en R las matrices empiezan en 1
    
    dimnames(b) = list(1:nrow(b), 1:ncol(b))
    
    
    
    #seleccionamos pares entre k,2*k-1
    flag <- 0
    for(i in 1:nrow(a)){
      
      if( ( (a[i,1] + k ) <= a[i,2]) && ( a[i,2] < (a[i,1] + 2*k)) ) {
        
        if(flag == 0){
          
          out <- rbind(a[i,])
          flag <- 1
          
          
        }else{
          
          out <- rbind(out, a[i,])
        }
        
        
      }
      
    }
    
    time1 <- Sys.time() 
    
    #eliminamos a, comentar si se quiere mantener la tabla intermedia  
    remove(a)
    timeArcos1 <- Sys.time()
    #añadimos los arcos segun MHM
    for(i in 1:nrow(out)){
      
      b[out[i,1]+1,out[i,2]+1] <- obtener_sse(out[i,1], out[i,2], elementos_ordenados, datos_standard) # el más 1 R comienza las matrices en 1
      
    }  
    
    timeArcos2 <- Sys.time()
    
    
    timeD1 <- Sys.time()
    #convertimos la matriz b en un grafo pesado, calculamos el camino mas corto
    grafo <- graph_from_adjacency_matrix(b, weighted=TRUE)
    path <- get.shortest.paths(grafo,1,(numRecords + 1)) # el temita
    camino <- as.numeric(path$vpath[[1]])
    timeD2 <- Sys.time()
    #R no tiene indice 0 en las columnas, realizamos arreglo
    camino <- camino - 1
    
    #Convertimos el camino de un dataframe de salida a una lista de elementos
    for(i in 1:(length(camino)-1)){
      if(i == 1){  
        salida <- c(camino[i],camino[i+1])
      }else{
        
        salida <- rbind(salida, c(camino[i],camino[i+1]))
        
      }  
    }
    
    salida <- as.matrix(salida)
    dimnames(salida) = list(1:nrow(salida),1:ncol(salida))
    
    
    
    grupos_salida <- obtener_grupos(salida, elementos_ordenados, k)
    remove(salida)
    
    #preparamos datos para calcular SST y SSE
    grupos <- unique(grupos_salida)
    dimnames(grupos) <- list(1:nrow(grupos),1:ncol(grupos))

    ##Cardinalidad media
    cardinalidad <- length(grupos[1,!is.na(grupos[1,])])
    tamaño_grupos <- length(grupos[,1])
    
    for(i in 2:tamaño_grupos){
      
      cardinalidad <- cardinalidad + length(grupos[i,!is.na(grupos[i,])])
      
    }
    
    cardinalidad_media <- cardinalidad/tamaño_grupos
    cardinalidad_media <- round(cardinalidad_media, digits = 4)
    
    
    
    
        
    #write.table(grupos, file = paste0(root_path, nombre_dataset,"_resultados/", nombre_dataset, "_gruposk_",k,".csv") , sep = ";", row.names = F, col.names = F)
    
    #Creamos matriz de centroides de los grupos
    
    datos <- as.matrix(datos_standard)
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
    
    dimnames(matrix_centroides) = list(1:nrow(matrix_centroides),1:ncol(matrix_centroides))
    
    #Calculamos Centroide del dataset
    datos <- as.matrix(datos_standard)
    centroide_dataset <- as.numeric(centroids(datos,rownames(datos), verbose=F)$means[,"(pooled)"])#calculamos el centrodie
    
    #Detenemos contador
    
    
    #Inicializamos SST a 0
    sst <- 0
    #datos <- as.matrix(datos_standard)
    #Calculamos SST
    for(i in 1:nrow(datos)){
      
      sst <-  sst + (datos[i,]-centroide_dataset)%*%(datos[i,]-centroide_dataset)
      
    } 
    
    
    #Inicializamos SSE
    sse <- 0	
    #datos <- as.matrix(datos_standard)
    
    for(i in 1:nrow(grupos)){
      
      for(j in 1:((2*k)-1)){
        
        if(!is.na(grupos[i,j])){
          
          sse <- sse + (datos[grupos[i,j],]- matrix_centroides[i,])%*%(datos[grupos[i,j],]- matrix_centroides[i,])
          
        }
      }	
    }
    
    
    
    #Calculamos Iloss
    iloss <- round((sse/sst)*100, digits = 4)
    sse <- round(sse, digits = 4)
    sst <- round(sst, digits = 4)
    
    #Compute time
    tiempo0 <- round(as.numeric(difftime(timeb, timea, units = "secs")), digits = 2)
    tiempoArcos <- round(as.numeric(difftime(timeArcos2, timeArcos1, units = "secs")), digits = 2)#MHM Arcos
    tiempo <- round(as.numeric(difftime(time1, time0, units = "secs")), digits = 2)
    tiempoDijkstra <- round(as.numeric(difftime(timeD2, timeD1, units = "secs")), digits = 2)#Dijskstra
    
    print(paste0("MicroTSP ; K_value ;", k, "; SSE ;", sse , "; SST ;", sst , "; Iloss ;", iloss, "; Tiempo_MHM_Pares ;", tiempo,"; Tiempo_MHM_Arcoss ;", tiempoArcos, "; Tiempo_MHM_Dijkstra ;", tiempoDijkstra,"; Tiempo_MDAV ; NA ; Tiempo_Grupos ; NA ", "; TiempoTSP ;", tiempo0, "; Longitud ;", longitud, "; K_media ;", cardinalidad_media))
    
    
    if(flag_retornar == 0){
      
      retornar <- rbind(c(k,sse,sst,iloss,tiempo))
      retornar <- as.matrix(retornar)
      dimnames(retornar) = list(1:nrow(retornar),1:ncol(retornar))
      
      flag_retornar = 1  
      
    }else{
      
      retornar <- rbind(retornar, c(k,sse,sst,iloss,tiempo))
      retornar <- as.matrix(retornar)
      dimnames(retornar) = list(1:nrow(retornar),1:ncol(retornar))
    }
    
  }#k
  return(retornar)
  
  
  
  
}# fin funcion


########################################################################
microConcordeDatos <- function(datos, k, metodo){
  #for(k in 5:5){
  
  
  #Dataset de entrada
  #separador_fichero_in = separador #fichero datos
  #cabeceras_fichero_in = cabeceras #fichero datos
  
  #Leemos datos
  #datos <- read.table(nombre_fichero_datos, header = cabeceras_fichero_in, sep = separador_fichero_in)
  datos_standard <- datos 
  
  #Iniciamos contador
  time0 <- Sys.time() 
  #Estandarizamos
  for (i in 1:ncol(datos)){
    media <- mean(datos[,i])
    desviacion_tipica <- sqrt(var(datos[,i]))
    datos_standard[,i] <- (datos[,i]-media)/desviacion_tipica 	
  }
  
  remove(datos)
  dimnames(datos_standard) = list(1:nrow(datos_standard),1:ncol(datos_standard))
  
  distancias<-dist(datos_standard)
  
  if(metodo == "concorde"){
    
    tsp <- as.TSP(distancias)
    tsp_dummy <- insert_dummy(tsp, label = "dummy")
    tour_solution <- solve_TSP(tsp_dummy, method = "concorde", control = list(precision = 6, verbose = FALSE))
    remove(tsp)
    
  }else{
    atsp <- as.ATSP(distancias)
    tsp_dummy <- insert_dummy(atsp, label = "dummy")
    tour_solution <- solve_TSP(tsp_dummy, method = metodo)
    remove(atsp)
  }
  
  path_solution <- cut_tour(tour_solution, cut = "dummy")
  elementos_ordenados <- as.numeric(labels(path_solution))
  
  
  ##########################################
  ################## MHM ###################
  ##########################################
  
  numRecords <- nrow(datos_standard)
  
  #elementos_ordenados <- as.numeric(labels(path_solution))
  
  #Creamos combination      
  a <- t(combn(c(0:numRecords), 2))
  
  b <- matrix(0,(numRecords + 1),(numRecords + 1)) #en R las matrices empiezan en 1
  
  dimnames(b) = list(1:nrow(b), 1:ncol(b))
  
  
  
  #seleccionamos pares entre k,2*k-1
  flag <- 0
  for(i in 1:nrow(a)){
    
    if( ( (a[i,1] + k ) <= a[i,2]) && ( a[i,2] < (a[i,1] + 2*k)) ) {
      
      if(flag == 0){
        
        out <- rbind(a[i,])
        flag <- 1
        
        
      }else{
        
        out <- rbind(out, a[i,])
      }
      
      
    }
    
  }
  
  #eliminamos a, comentar si se quiere mantener la tabla intermedia  
  remove(a)
  
  #añadimos los arcos segun MHM
  for(i in 1:nrow(out)){
    
    b[out[i,1]+1,out[i,2]+1] <- obtener_sse(out[i,1], out[i,2], elementos_ordenados, datos_standard) # el más 1 R comienza las matrices en 1
    
  }  
  
  #convertimos la matriz b en un grafo pesado, calculamos el camino mas corto
  grafo <- graph_from_adjacency_matrix(b, weighted=TRUE)
  path <- get.shortest.paths(grafo,1,(numRecords + 1)) # el temita
  camino <- as.numeric(path$vpath[[1]])
  
  #R no tiene indice 0 en las columnas, realizamos arreglo
  camino <- camino - 1
  
  #Convertimos el camino de un dataframe de salida a una lista de elementos
  for(i in 1:(length(camino)-1)){
    if(i == 1){  
      salida <- c(camino[i],camino[i+1])
    }else{
      
      salida <- rbind(salida, c(camino[i],camino[i+1]))
      
    }  
  }
  
  salida <- as.matrix(salida)
  dimnames(salida) = list(1:nrow(salida),1:ncol(salida))
  grupos_salida <- obtener_grupos(salida, elementos_ordenados, k)
  remove(salida)
  
  #preparamos datos para calcular SST y SSE
  grupos <- unique(grupos_salida)
  dimnames(grupos) <- list(1:nrow(grupos),1:ncol(grupos))
  
  
  
  #preparamos datos para calcular SST y SSE
  grupos <- unique(grupos_salida)
  dimnames(grupos) <- list(1:nrow(grupos),1:ncol(grupos))
  
  #write.table(grupos, file = paste0(root_path, nombre_dataset,"_resultados/", nombre_dataset, "_gruposk_",k,".csv") , sep = ";", row.names = F, col.names = F)
  
  #Creamos matriz de centroides de los grupos
  
  datos <- as.matrix(datos_standard)
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
  
  dimnames(matrix_centroides) = list(1:nrow(matrix_centroides),1:ncol(matrix_centroides))
  
  #Calculamos Centroide del dataset
  datos <- as.matrix(datos_standard)
  centroide_dataset <- as.numeric(centroids(datos,rownames(datos), verbose=F)$means[,"(pooled)"])#calculamos el centrodie
  
  #Detenemos contador
  time1 <- Sys.time() 
  
  #Inicializamos SST a 0
  sst <- 0
  #datos <- as.matrix(datos_standard)
  #Calculamos SST
  for(i in 1:nrow(datos)){
    
    sst <-  sst + (datos[i,]-centroide_dataset)%*%(datos[i,]-centroide_dataset)
    
  } 
  
  
  #Inicializamos SSE
  sse <- 0	
  #datos <- as.matrix(datos_standard)
  
  for(i in 1:nrow(grupos)){
    
    for(j in 1:((2*k)-1)){
      
      if(!is.na(grupos[i,j])){
        
        sse <- sse + (datos[grupos[i,j],]- matrix_centroides[i,])%*%(datos[grupos[i,j],]- matrix_centroides[i,])
        
      }
    }	
  }
  
  
  
  #Calculamos Iloss
  iloss <- round((sse/sst)*100, digits = 4)
  sse <- round(sse, digits = 4)
  sst <- round(sst, digits = 4)
  
  #Compute time
  
  tiempo <- round(as.numeric(difftime(time1, time0, units = "secs")), digits = 2)
  print(paste0("MicroTSP --> K_value: ", k, " SSE: ", sse , " SST: ", sst , " Iloss: ", iloss, " Tiempo: ", tiempo))
  
  retornar <- rbind(c(k,sse,sst,iloss,tiempo))
  retornar <- as.matrix(retornar)
  dimnames(retornar) = list(1:nrow(retornar),1:ncol(retornar))
  
  return(retornar)

}# fin funcion

########################################################################

# microConcorde("/Users/samsa/URV/Databases/census.csv")  

#microConcordeAll("/Users/samsa/URV/Databases/File_2_ID_2015_Domains_of_deprivation.csv", "concorde", ";", F)

########################################################################


microConcordeAllPaths <- function(nombre_fichero_datos, metodo, separador, cabeceras){
  
  #Dataset de entrada
  separador_fichero_in = separador #fichero datos
  cabeceras_fichero_in = cabeceras #fichero datos
  
  #Leemos datos
  datos <- read.table(nombre_fichero_datos, header = cabeceras_fichero_in, sep = separador_fichero_in)
  datos_standard <- datos 
  
  #Iniciamos contador
  time0 <- Sys.time() 
  #Estandarizamos
  for (i in 1:ncol(datos)){
    media <- mean(datos[,i])
    desviacion_tipica <- sqrt(var(datos[,i]))
    datos_standard[,i] <- (datos[,i]-media)/desviacion_tipica 	
  }
  
  #remove(datos)
  dimnames(datos_standard) = list(1:nrow(datos_standard),1:ncol(datos_standard))
  
  distancias<-dist(datos_standard)
  
  #Calculo de los Hamiltonian path Concorde
  
  distancias <- as.matrix(distancias)

  flag_sse <- 0
  
  numRecords <- nrow(datos_standard)
  
  print(numRecords)
  
  for(i in 1:numRecords){#for tsp
    
    #print(i)
    atsp <- as.ATSP(distancias)  	
    atsp[,i] <- 0
    
    tsp <- reformulate_ATSP_as_TSP(atsp)
    
    tour <- solve_TSP(tsp, method = "concorde", control = list(precision = 6, verbose = FALSE))
    
    tour <- as.TOUR(tour[tour <= n_of_cities(atsp)])
    
    print(paste0("No cities: ", n_of_cities(atsp)))
    
    path_solution <- cut_tour(tour, i, exclude_cut = FALSE)
    
    elementos_ordenados <- as.numeric(labels(path_solution))
    
    
    

  
  ##########################################
  ################## MHM ###################
  ##########################################
  
  numRecords <- nrow(datos_standard)
  flag_retornar <- 0
  #elementos_ordenados <- as.numeric(labels(path_solution))
  for(k in 3:3){
    #Creamos combination      
    a <- t(combn(c(0:numRecords), 2))
    
    b <- matrix(0,(numRecords + 1),(numRecords + 1)) #en R las matrices empiezan en 1
    
    dimnames(b) = list(1:nrow(b), 1:ncol(b))
    
    
    
    #seleccionamos pares entre k,2*k-1
    flag <- 0
    for(i in 1:nrow(a)){
      
      if( ( (a[i,1] + k ) <= a[i,2]) && ( a[i,2] < (a[i,1] + 2*k)) ) {
        
        if(flag == 0){
          
          out <- rbind(a[i,])
          flag <- 1
          
          
        }else{
          
          out <- rbind(out, a[i,])
        }
        
        
      }
      
    }
    
    #eliminamos a, comentar si se quiere mantener la tabla intermedia  
    remove(a)
    
    #añadimos los arcos segun MHM
    for(i in 1:nrow(out)){
      
      b[out[i,1]+1,out[i,2]+1] <- obtener_sse(out[i,1], out[i,2], elementos_ordenados, datos_standard) # el más 1 R comienza las matrices en 1
      
    }  
    
    #convertimos la matriz b en un grafo pesado, calculamos el camino mas corto
    grafo <- graph_from_adjacency_matrix(b, weighted=TRUE)
    path <- get.shortest.paths(grafo,1,(numRecords + 1)) # el temita
    camino <- as.numeric(path$vpath[[1]])
    
    #R no tiene indice 0 en las columnas, realizamos arreglo
    camino <- camino - 1
    
    #Convertimos el camino de un dataframe de salida a una lista de elementos
    for(i in 1:(length(camino)-1)){
      if(i == 1){  
        salida <- c(camino[i],camino[i+1])
      }else{
        
        salida <- rbind(salida, c(camino[i],camino[i+1]))
        
      }  
    }
    
    salida <- as.matrix(salida)
    dimnames(salida) = list(1:nrow(salida),1:ncol(salida))
    grupos_salida <- obtener_grupos(salida, elementos_ordenados, k)
    remove(salida)
    
    #preparamos datos para calcular SST y SSE
    grupos <- unique(grupos_salida)
    dimnames(grupos) <- list(1:nrow(grupos),1:ncol(grupos))
    
    
    
    #preparamos datos para calcular SST y SSE
    grupos <- unique(grupos_salida)
    dimnames(grupos) <- list(1:nrow(grupos),1:ncol(grupos))
    
    #write.table(grupos, file = paste0(root_path, nombre_dataset,"_resultados/", nombre_dataset, "_gruposk_",k,".csv") , sep = ";", row.names = F, col.names = F)
    
    #Creamos matriz de centroides de los grupos
    
    datos <- as.matrix(datos_standard)
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
    
    dimnames(matrix_centroides) = list(1:nrow(matrix_centroides),1:ncol(matrix_centroides))
    
    #Calculamos Centroide del dataset
    datos <- as.matrix(datos_standard)
    centroide_dataset <- as.numeric(centroids(datos,rownames(datos), verbose=F)$means[,"(pooled)"])#calculamos el centrodie
    
    #Detenemos contador
    time1 <- Sys.time() 
    
    #Inicializamos SST a 0
    sst <- 0
    #datos <- as.matrix(datos_standard)
    #Calculamos SST
    for(i in 1:nrow(datos)){
      
      sst <-  sst + (datos[i,]-centroide_dataset)%*%(datos[i,]-centroide_dataset)
      
    } 
    
    
    #Inicializamos SSE
    sse <- 0	
    #datos <- as.matrix(datos_standard)
    
    for(i in 1:nrow(grupos)){
      
      for(j in 1:((2*k)-1)){
        
        if(!is.na(grupos[i,j])){
          
          sse <- sse + (datos[grupos[i,j],]- matrix_centroides[i,])%*%(datos[grupos[i,j],]- matrix_centroides[i,])
          
        }
      }	
    }
    
    
    
    #Calculamos Iloss
    iloss <- round((sse/sst)*100, digits = 4)
    sse <- round(sse, digits = 4)
    sst <- round(sst, digits = 4)
    
    #Compute time
    
    tiempo <- round(as.numeric(difftime(time1, time0, units = "secs")), digits = 2)
    
    #print(paste0("MicroTSP --> K_value: ", k, " SSE: ", sse , " SST: ", sst , " Iloss: ", iloss, " Tiempo: ", tiempo))
    
    
    if(flag_retornar == 0){
      
      retornar <- rbind(c(k,sse,sst,iloss,tiempo))
      retornar <- as.matrix(retornar)
      dimnames(retornar) = list(1:nrow(retornar),1:ncol(retornar))
      
      flag_retornar = 1  
      
    }else{
      
      retornar <- rbind(retornar, c(k,sse,sst,iloss,tiempo))
      retornar <- as.matrix(retornar)
      dimnames(retornar) = list(1:nrow(retornar),1:ncol(retornar))
    }
    
    
    
  }#k

  if(flag_sse == 0){
    
    sse_minimo = sse
    retornar_minimo <- retornar
    flag_sse <- 1
    
  }else{
  
      if(sse_minimo > sse){
        
        sse_minimo = sse
        retornar_minimo <- retornar
      }  
  }    
    
      
  
  }#todos los paths  
  
  #debe ser el menor valor
  
  return(retornar_minimo)
  
  
  
  
}# fin funcion

########################################################################
