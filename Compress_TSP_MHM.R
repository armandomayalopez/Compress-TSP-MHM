rm(list = ls())
library("sda")
library("igraph")
library("stringr")
library("sdcMicro")
source('iloss_function.R')
source('microTSP.R')

file_separator = ";"
header = F 
filename_dataset <- "datasetname"
path2dataset <- "/path/to/dataset"
root_path <- "/path/to/out"

nombre_fichero_datos <- paste0(path2dataset, filename_dataset, ".csv")

timeIni <- Sys.time()
write(paste0("filename;k1step;macro;k;sse;sst;iloss;tiempo"), 
      file = paste0(root_path,"micro2step_",filename_dataset,"_sse_salida_",timeIni,".csv") , sep = ";", append = TRUE)
print("+-------------------------------------------+") 
dataset_entrada <- read.table(nombre_fichero_datos, header = header, sep = file_separator)
datos_standard <- dataset_entrada 
datos_standard <- as.matrix(datos_standard)
for (i in 1:ncol(dataset_entrada)){
  media <- mean(dataset_entrada[,i])
  desviacion_tipica <- sqrt(var(dataset_entrada[,i]))
  datos_standard[,i] <- (dataset_entrada[,i]-media)/desviacion_tipica 	
}
dimnames(datos_standard) = list(1:nrow(datos_standard),1:ncol(datos_standard))
centroide_dataset <- as.numeric(centroids(datos_standard,rownames(datos_standard), verbose=F)$means[,"(pooled)"])
sst <- 0

for(i in 1:nrow(datos_standard)){
  sst <-  sst + (datos_standard[i,]-centroide_dataset)%*%(datos_standard[i,]-centroide_dataset)
} 


for(k1step in 2:5){

  timea <- Sys.time()
  #Init Compress
  m1step <- microaggregation(dataset_entrada,aggr = k1step,method = "mdav")
  y1step <- m1step$mx
  grupos1step <- grupos_microagregados(dataset_entrada,y1step,k1step)
  grupos1step <- as.matrix(grupos1step)
  dimnames(grupos1step) = list(1:nrow(grupos1step),1:ncol(grupos1step))
  numgrupos1step <- nrow(grupos1step)
  y1Centroides <- rbind(na.omit(y1step[grupos1step[1,1],]))
  for(i in 2:numgrupos1step){
    y1Centroides <- rbind(y1Centroides, na.omit(y1step[grupos1step[i,1],]))
  }
  dimnames(y1Centroides) = list(1:nrow(y1Centroides),1:ncol(y1Centroides))
  centorides_standard <- y1Centroides
  for (i in 1:ncol(y1Centroides)){
    cmedia <- mean(y1Centroides[,i])
    cdesviacion_tipica <- sqrt(var(y1Centroides[,i]))
    centorides_standard[,i] <- (y1Centroides[,i]-cmedia)/cdesviacion_tipica 	
  }
  dimnames(centorides_standard) = list(1:nrow(centorides_standard),1:ncol(centorides_standard))
  distancias_centorides <- dist(centorides_standard)
  tsp <- as.TSP(distancias_centorides)
  tsp_dummy <- insert_dummy(tsp, label = "dummy")
  tour_solution <- solve_TSP(tsp_dummy, method = "concorde", control = list(precision = 2, verbose = FALSE))
  timeb <- Sys.time()
  remove(tsp)
  path_solution <- cut_tour(tour_solution, cut = "dummy")
  centroides_ordenados <- as.numeric(labels(path_solution))
  grupos1step <- grupos1step[centroides_ordenados[1:length(centroides_ordenados)],]
  
  flag <- 0
  for(i in 1:nrow(grupos1step)){
    
    if(flag == 0){
      
      fila_sin_na <- grupos1step[i,1]
      for(j in 2:ncol(grupos1step)){
        if(!is.na(grupos1step[i,j])){
          fila_sin_na <- c(fila_sin_na,grupos1step[i,j])
        }
      }
      camino_ordenado <- fila_sin_na
      flag <- 1
    }else{
      fila_sin_na <- grupos1step[i,1]
      for(j in 2:ncol(grupos1step)){
        if(!is.na(grupos1step[i,j])){
          fila_sin_na <- c(fila_sin_na,grupos1step[i,j])
        }
      }
      camino_ordenado <- c(camino_ordenado,fila_sin_na)
    }
  }
  elementos_ordenados <- as.numeric(camino_ordenado)
  numRecords <- nrow(datos_standard)
  flag_retornar <- 0
  for(k in 3:10){
    time0 <- Sys.time() 
    a <- t(combn(c(0:numRecords), 2))
    b <- matrix(0,(numRecords + 1),(numRecords + 1))
    dimnames(b) = list(1:nrow(b), 1:ncol(b))
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
    remove(a)
    
    for(i in 1:nrow(out)){
      b[out[i,1]+1,out[i,2]+1] <- obtener_sse(out[i,1], out[i,2], elementos_ordenados, datos_standard)
    }  
    
    grafo <- graph_from_adjacency_matrix(b, weighted=TRUE)
    path <- get.shortest.paths(grafo,1,(numRecords + 1))
    camino <- as.numeric(path$vpath[[1]])
    camino <- camino - 1
    for(i in 1:(length(camino)-1)){
      if(i == 1){  
        salida <- c(camino[i],camino[i+1])
      }else{
        salida <- rbind(salida, c(camino[i],camino[i+1]))
      }  
    }
    
    salida <- as.matrix(salida)
    dimnames(salida) = list(1:nrow(salida),1:ncol(salida))
    
    time1 <- Sys.time()
    grupos_salida <- obtener_grupos(salida, elementos_ordenados, k)
    remove(salida)
    grupos <- unique(grupos_salida)
    dimnames(grupos) <- list(1:nrow(grupos),1:ncol(grupos))

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
    datos <- as.matrix(datos_standard)
    centroide_dataset <- as.numeric(centroids(datos,rownames(datos), verbose=F)$means[,"(pooled)"])
    sst <- 0
    for(i in 1:nrow(datos)){
      sst <-  sst + (datos[i,]-centroide_dataset)%*%(datos[i,]-centroide_dataset)
    }
    sse <- 0	
    for(i in 1:nrow(grupos)){
      for(j in 1:((2*k)-1)){
        if(!is.na(grupos[i,j])){
          sse <- sse + (datos[grupos[i,j],]- matrix_centroides[i,])%*%(datos[grupos[i,j],]- matrix_centroides[i,])
        }
      }	
    }
    iloss <- round((sse/sst)*100, digits = 4)
    sse <- round(sse, digits = 4)
    sst <- round(sst, digits = 4)
    
    tiempo0 <- round(as.numeric(difftime(timeb, timea, units = "secs")), digits = 2)
    tiempo <- round(as.numeric(difftime(time1, time0, units = "secs")), digits = 2)
    
    print(paste0("Compress(" , k1step , ") --> K_value: ", k, " SSE: ", sse , " SST: ", sst , " Iloss: ", iloss, " TiempoMHM: ", tiempo, " Tiempo_MDAV&TSP: ", tiempo0))
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
  print("+-------------------------------------------+")  
}






