#DA ESEGUIRE DIRETTAMENTE DOPO remeshing_CAD_mesh.R

# computazione grafo sulla superficie della mesh

# calcolo vertici appartenenti ad ogni triangolo 'padre' della triangolazione originale
  # nota: un vertice può appartenere a più triangoli
father_to_subverts_count = rep(0, nrow(triangles))
for(i in 1:nrow(verts_tmp)){
  father_to_subverts_count[verts_tmp[[i,4]]] = father_to_subverts_count[verts_tmp[[i,4]]] + 1
  
  if(sum(verts_tmp[i,5:7]) == 1){ #contatto con un edge
    if(verts_tmp[[i,5]]){ #edge 1
      tr_v1 = vert_triangles[triangles[verts_tmp[[i,4]],1],] #triangoli a contatto con vertice 1
      tr_v1 = tr_v1[!is.na(tr_v1) & tr_v1!=verts_tmp[[i,4]]]
      tr_v2 = vert_triangles[triangles[verts_tmp[[i,4]],2],] #triangoli a contatto con vertice 2
      tr_v2 = tr_v2[!is.na(tr_v2) & tr_v2!=verts_tmp[[i,4]]]
      
      other_tr = tr_v1[tr_v1%in%tr_v2] #triangolo a contatto con edge 1
    }else{
      if(verts_tmp[[i,6]]){ #edge 2
        tr_v2 = vert_triangles[triangles[verts_tmp[[i,4]],2],] #triangoli a contatto con vertice 2
        tr_v2 = tr_v2[!is.na(tr_v2) & tr_v2!=verts_tmp[[i,4]]]
        tr_v3 = vert_triangles[triangles[verts_tmp[[i,4]],3],] #triangoli a contatto con vertice 3
        tr_v3 = tr_v3[!is.na(tr_v3) & tr_v3!=verts_tmp[[i,4]]]
        
        other_tr = tr_v2[tr_v2%in%tr_v3] #triangolo a contatto con edge 2
      }else{
        if(verts_tmp[[i,7]]){ #edge 3
          tr_v1 = vert_triangles[triangles[verts_tmp[[i,4]],1],] #triangoli a contatto con vertice 1
          tr_v1 = tr_v1[!is.na(tr_v1) & tr_v1!=verts_tmp[[i,4]]]
          tr_v3 = vert_triangles[triangles[verts_tmp[[i,4]],3],] #triangoli a contatto con vertice 3
          tr_v3 = tr_v3[!is.na(tr_v3) & tr_v3!=verts_tmp[[i,4]]]
          
          other_tr = tr_v3[tr_v3%in%tr_v1] #triangolo a contatto con edge 3
        }
      }
    }
    
    for(oth in other_tr){
      father_to_subverts_count[oth] = father_to_subverts_count[oth] + 1
    }
  }
  
  if(sum(verts_tmp[i,5:7]) == 2){ #contatto tramite vertice
    if(verts_tmp[[i,5]] & verts_tmp[[i,7]]){ #vertice 1
      tr_v = vert_triangles[triangles[verts_tmp[[i,4]],1],] #triangoli a contatto con vertice 1
      tr_v = tr_v[!is.na(tr_v) & tr_v!=verts_tmp[[i,4]]]
    }else{
      if(verts_tmp[[i,5]] & verts_tmp[[i,6]]){ #vertice 2
        tr_v = vert_triangles[triangles[verts_tmp[[i,4]],2],] #triangoli a contatto con vertice 2
        tr_v = tr_v[!is.na(tr_v) & tr_v!=verts_tmp[[i,4]]]
      }else{
        if(verts_tmp[[i,6]] & verts_tmp[[i,7]]){ #vertice 3
          tr_v = vert_triangles[triangles[verts_tmp[[i,4]],3],] #triangoli a contatto con vertice 3
          tr_v = tr_v[!is.na(tr_v) & tr_v!=verts_tmp[[i,4]]]
        }
      }
    }
    
    for(oth in tr_v){
      father_to_subverts_count[oth] = father_to_subverts_count[oth] + 1
    }
  }
  
  if(i%%10000==0) print(i)
}
max(father_to_subverts_count)
father_to_subverts_list = list()
for(i in 1:nrow(triangles)){
  father_to_subverts_list[[i]] = rep(NA, father_to_subverts_count[i])
}
father_to_subverts_count = rep(0, nrow(triangles))
for(i in 1:nrow(verts_tmp)){
  father_to_subverts_count[verts_tmp[[i,4]]] = father_to_subverts_count[verts_tmp[[i,4]]] + 1
  father_to_subverts_list[[verts_tmp[[i,4]]]][father_to_subverts_count[verts_tmp[[i,4]]]] = i
  
  if(sum(verts_tmp[i,5:7]) == 1){ #contatto con un edge
    if(verts_tmp[[i,5]]){ #edge 1
      tr_v1 = vert_triangles[triangles[verts_tmp[[i,4]],1],] #triangoli a contatto con vertice 1
      tr_v1 = tr_v1[!is.na(tr_v1) & tr_v1!=verts_tmp[[i,4]]]
      tr_v2 = vert_triangles[triangles[verts_tmp[[i,4]],2],] #triangoli a contatto con vertice 2
      tr_v2 = tr_v2[!is.na(tr_v2) & tr_v2!=verts_tmp[[i,4]]]
      
      other_tr = tr_v1[tr_v1%in%tr_v2] #triangolo a contatto con edge 1
    }else{
      if(verts_tmp[[i,6]]){ #edge 2
        tr_v2 = vert_triangles[triangles[verts_tmp[[i,4]],2],] #triangoli a contatto con vertice 2
        tr_v2 = tr_v2[!is.na(tr_v2) & tr_v2!=verts_tmp[[i,4]]]
        tr_v3 = vert_triangles[triangles[verts_tmp[[i,4]],3],] #triangoli a contatto con vertice 3
        tr_v3 = tr_v3[!is.na(tr_v3) & tr_v3!=verts_tmp[[i,4]]]
        
        other_tr = tr_v2[tr_v2%in%tr_v3] #triangolo a contatto con edge 2
      }else{
        if(verts_tmp[[i,7]]){ #edge 3
          tr_v1 = vert_triangles[triangles[verts_tmp[[i,4]],1],] #triangoli a contatto con vertice 1
          tr_v1 = tr_v1[!is.na(tr_v1) & tr_v1!=verts_tmp[[i,4]]]
          tr_v3 = vert_triangles[triangles[verts_tmp[[i,4]],3],] #triangoli a contatto con vertice 3
          tr_v3 = tr_v3[!is.na(tr_v3) & tr_v3!=verts_tmp[[i,4]]]
          
          other_tr = tr_v3[tr_v3%in%tr_v1] #triangolo a contatto con edge 3
        }
      }
    }
    
    for(oth in other_tr){
      father_to_subverts_count[oth] = father_to_subverts_count[oth] + 1
      father_to_subverts_list[[oth]][father_to_subverts_count[oth]] = i
    }
  }
  
  if(sum(verts_tmp[i,5:7]) == 2){ #contatto tramite vertice
    if(verts_tmp[[i,5]] & verts_tmp[[i,7]]){ #vertice 1
      tr_v = vert_triangles[triangles[verts_tmp[[i,4]],1],] #triangoli a contatto con vertice 1
      tr_v = tr_v[!is.na(tr_v) & tr_v!=verts_tmp[[i,4]]]
    }else{
      if(verts_tmp[[i,5]] & verts_tmp[[i,6]]){ #vertice 2
        tr_v = vert_triangles[triangles[verts_tmp[[i,4]],2],] #triangoli a contatto con vertice 2
        tr_v = tr_v[!is.na(tr_v) & tr_v!=verts_tmp[[i,4]]]
      }else{
        if(verts_tmp[[i,6]] & verts_tmp[[i,7]]){ #vertice 3
          tr_v = vert_triangles[triangles[verts_tmp[[i,4]],3],] #triangoli a contatto con vertice 3
          tr_v = tr_v[!is.na(tr_v) & tr_v!=verts_tmp[[i,4]]]
        }
      }
    }
    
    for(oth in tr_v){
      father_to_subverts_count[oth] = father_to_subverts_count[oth] + 1
      father_to_subverts_list[[oth]][father_to_subverts_count[oth]] = i
    }
  }
  
  if(i%%10000==0) print(i)
}

# calcolo archi (ci saranno archi ridondanti, rimossi dopo)
BATCH_LENGTH = 10000

unique_con_list_mat_list = list()
for(b in 1:(floor(nrow(verts_tmp) / BATCH_LENGTH) + 1)){
  gc()
  
  inds = ((b-1)*BATCH_LENGTH+1):(b*BATCH_LENGTH)
  if(max(inds) > nrow(verts_tmp)){
    inds = min(inds):nrow(verts_tmp)
  }
  
  con_list = list()
  num_cons = 0
  for(i in inds){
    for(dest in father_to_subverts_list[[verts_tmp[[i,4]]]]){ # connetto con tutti i vertici del padre 
      #if(i < dest){
      num_cons = num_cons + 1
      dist = sqrt( sum((verts_tmp[i,1:3] - verts_tmp[dest,1:3])^2) )
      con_list[[num_cons]] = c(i, dest, dist) #origine, destinazione, distanza
      #}
    }
    
    if(sum(verts_tmp[i,5:7]) == 1){ #contatto con un edge
      if(verts_tmp[[i,5]]){ #edge 1
        tr_v1 = vert_triangles[triangles[verts_tmp[[i,4]],1],] #triangoli a contatto con vertice 1
        tr_v1 = tr_v1[!is.na(tr_v1) & tr_v1!=verts_tmp[[i,4]]]
        tr_v2 = vert_triangles[triangles[verts_tmp[[i,4]],2],] #triangoli a contatto con vertice 2
        tr_v2 = tr_v2[!is.na(tr_v2) & tr_v2!=verts_tmp[[i,4]]]
        
        other_tr = tr_v1[tr_v1%in%tr_v2] #triangolo a contatto con edge 1
      }else{
        if(verts_tmp[[i,6]]){ #edge 2
          tr_v2 = vert_triangles[triangles[verts_tmp[[i,4]],2],] #triangoli a contatto con vertice 2
          tr_v2 = tr_v2[!is.na(tr_v2) & tr_v2!=verts_tmp[[i,4]]]
          tr_v3 = vert_triangles[triangles[verts_tmp[[i,4]],3],] #triangoli a contatto con vertice 3
          tr_v3 = tr_v3[!is.na(tr_v3) & tr_v3!=verts_tmp[[i,4]]]
          
          other_tr = tr_v2[tr_v2%in%tr_v3] #triangolo a contatto con edge 2
        }else{
          if(verts_tmp[[i,7]]){ #edge 3
            tr_v1 = vert_triangles[triangles[verts_tmp[[i,4]],1],] #triangoli a contatto con vertice 1
            tr_v1 = tr_v1[!is.na(tr_v1) & tr_v1!=verts_tmp[[i,4]]]
            tr_v3 = vert_triangles[triangles[verts_tmp[[i,4]],3],] #triangoli a contatto con vertice 3
            tr_v3 = tr_v3[!is.na(tr_v3) & tr_v3!=verts_tmp[[i,4]]]
            
            other_tr = tr_v3[tr_v3%in%tr_v1] #triangolo a contatto con edge 3
          }
        }
      }
      
      for(dest in father_to_subverts_list[[other_tr]]){ # connetto con tutti i vertici del padre adiacente
        #if(i < dest){
        num_cons = num_cons + 1
        dist = sqrt( sum((verts_tmp[i,1:3] - verts_tmp[dest,1:3])^2) )
        con_list[[num_cons]] = c(i, dest, dist) #origine, destinazione, distanza
        #}
      }
    }
    
    if(sum(verts_tmp[i,5:7]) == 2){ #contatto tramite vertice
      if(verts_tmp[[i,5]] & verts_tmp[[i,7]]){ #vertice 1
        tr_v = vert_triangles[triangles[verts_tmp[[i,4]],1],] #triangoli a contatto con vertice 1
        tr_v = tr_v[!is.na(tr_v) & tr_v!=verts_tmp[[i,4]]]
      }else{
        if(verts_tmp[[i,5]] & verts_tmp[[i,6]]){ #vertice 2
          tr_v = vert_triangles[triangles[verts_tmp[[i,4]],2],] #triangoli a contatto con vertice 2
          tr_v = tr_v[!is.na(tr_v) & tr_v!=verts_tmp[[i,4]]]
        }else{
          if(verts_tmp[[i,6]] & verts_tmp[[i,7]]){ #vertice 3
            tr_v = vert_triangles[triangles[verts_tmp[[i,4]],3],] #triangoli a contatto con vertice 3
            tr_v = tr_v[!is.na(tr_v) & tr_v!=verts_tmp[[i,4]]]
          }
        }
      }
      
      for(other_tr in tr_v){
        for(dest in father_to_subverts_list[[other_tr]]){ # connetto con tutti i vertici del padre adiacente
          #if(i < dest){
          num_cons = num_cons + 1
          dist = sqrt( sum((verts_tmp[i,1:3] - verts_tmp[dest,1:3])^2) )
          con_list[[num_cons]] = c(i, dest, dist) #origine, destinazione, distanza
          #}
        }
      }
    }
    
    if(i%%500 == 0) print(i)
  }
  
  con_list_mat = list.rbind(con_list)
  print(nrow(con_list_mat))
  
  for(i in 1:nrow(con_list_mat)){
    if(con_list_mat[i,1] == con_list_mat[i,2]){
      con_list_mat[i,1] = -1
      con_list_mat[i,2] = -1
      con_list_mat[i,3] = -1
    }else{
      if(con_list_mat[i,1] > con_list_mat[i,2]){
        tmp = con_list_mat[i,1]
        con_list_mat[i,1] = con_list_mat[i,2]
        con_list_mat[i,2] = tmp
      }
    }
  }
  
  unique_con_list_mat = unique(con_list_mat)
  unique_con_list_mat = unique_con_list_mat[unique_con_list_mat[,1]!=(-1),]
  print(nrow(unique_con_list_mat))
  
  unique_con_list_mat_list[[b]] = unique_con_list_mat
  
  print(b)
}

for(i in 1:length(unique_con_list_mat_list)){
  unique_con_list_mat_list[[i]] = data.frame(v1 = unique_con_list_mat_list[[i]][,1],
                                             v2 = unique_con_list_mat_list[[i]][,2],
                                             v3 = unique_con_list_mat_list[[i]][,3])
}

con_mat_complete = bind_rows(unique_con_list_mat_list)

# rimozione archi duplicati
unique_con_mat_complete = unique(con_mat_complete)
nrow(unique_con_mat_complete) # numero archi

# salvo la lista degli archi
options(scipen = 999)
write.table(unique_con_mat_complete, "unique_con_mat_complete.csv", row.names=FALSE, col.names=FALSE, sep=",")
options(scipen = 5)

con = file(description = "neigh_verts_list_pert.txt", open = "r")
for(tr_id in 1:10){
  line = readLines(con = con, n = 1, ok = TRUE, warn = TRUE, encoding = "unknown", skipNul = FALSE)
  tmp = strsplit(line, " ")[[1]]
  
  
  if(length(tmp) > 0){
    neigh_verts = strtoi(tmp) + 1 #indici in un certo range dal triangolo 'tr_id' (+1 perchè viene da c++)
  }else{
    print("problema")
  }
}
close(con)

x <- as.matrix(verts_tmp %>% select(X, Y, Z))[neigh_verts,1]
y <- as.matrix(verts_tmp %>% select(X, Y, Z))[neigh_verts,2]
z <- as.matrix(verts_tmp %>% select(X, Y, Z))[neigh_verts,3]

k=1
col = as.numeric(neigh_verts %in% c(strtoi(c(unique_con_mat_complete[unique_con_mat_complete[,1] == father_to_subverts_list[[tr_id]][k],,drop=FALSE][,2], unique_con_mat_complete[unique_con_mat_complete[,2] == father_to_subverts_list[[tr_id]][k],,drop=FALSE][,1])))) + 1
col[neigh_verts == father_to_subverts_list[[tr_id]][k]] = 3

plot_ly(x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'markers',
        marker = list(size = 3, color = col)) # in rosso il vertice di partenza, in arancione i vertici raggiungibili con un arco
k=k+1


sqrt(qchisq(0.999999, df = 2)) * 0.2
cutoff_range = 1

# creazione grafo per calcolo degli intorni, tolgo gli archi più grandi del cutoff e metto self loops
unique_con_mat_complete_cutoff = unique_con_mat_complete[unique_con_mat_complete[,3] < cutoff_range + 1e-10,]
self_loops = cbind(1:nrow(verts_tmp), 1:nrow(verts_tmp), 0)
unique_con_mat_complete_cutoff = rbind(as.matrix(unique_con_mat_complete_cutoff), self_loops)
nrow(unique_con_mat_complete_cutoff)

con_mat = cbind(as.character(unique_con_mat_complete_cutoff[,1]), as.character(unique_con_mat_complete_cutoff[,2]))
graph = graph_from_edgelist(con_mat, directed = FALSE)
edge.attributes(graph)$weight = unique_con_mat_complete_cutoff[,3]


con = file(description = "neigh_verts_list_pert.txt", open = "r")
for(tr_id in 1:10){
  line = readLines(con = con, n = 1, ok = TRUE, warn = TRUE, encoding = "unknown", skipNul = FALSE)
  tmp = strsplit(line, " ")[[1]]
  
  if(length(tmp) > 0){
    neigh_verts = strtoi(tmp) + 1 #indici in un certo range dal triangolo 'tr_id' (+1 perchè viene da c++)
    #neigh_verts = unique(unlist(father_to_subverts_list[unique(verts_tmp[neigh_verts,] %>% pull(father_id))]))
    
    geod_dist_tmp = distances(subgraph(graph, v = as.character(neigh_verts)), v = as.character(father_to_subverts_list[[tr_id]]), to = as.character(neigh_verts))
  }else{
    print("problema")
  }
  
  if(tr_id%%10 == 0) print(tr_id)
}
close(con)

x <- as.matrix(verts_tmp %>% select(X, Y, Z))[neigh_verts,1]
y <- as.matrix(verts_tmp %>% select(X, Y, Z))[neigh_verts,2]
z <- as.matrix(verts_tmp %>% select(X, Y, Z))[neigh_verts,3]

k=1
colors = geod_dist_tmp[k,]
colors = colorRamp(c("green", "red"))((colors-min(colors))/(max(colors)-min(colors)))
colors <- rgb(colors/255)

plot_ly(x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'markers',
        marker = list(size = 4, color = colors))
k=k+1


# calcolo triangoli figli per ogni triangolo padre
father_to_subtriangles_count = as.vector((triangles_tmp %>% group_by(father_id) %>% count())[,2])$n
father_to_subtriangles_list = list()
cur_sum_tmp = 0
for(i in 1:nrow(triangles)){
  father_to_subtriangles_list[[i]] = (cur_sum_tmp + 1):(cur_sum_tmp + father_to_subtriangles_count[i])
  cur_sum_tmp = cur_sum_tmp + father_to_subtriangles_count[i]
}

# calcolo triangoli "figli" di cui ogni vertice fa parte
vert_subtriangles = matrix(NA, nrow=nrow(verts_tmp), ncol=12)
vert_subtriangles_count = rep(0, nrow(verts_tmp))
for(i in 1:nrow(triangles_tmp)){
  tmp = as.numeric(triangles_tmp[i,c(1,2,3)])
  vert_subtriangles_count[tmp] = vert_subtriangles_count[tmp] + 1
  vert_subtriangles[tmp[1], vert_subtriangles_count[tmp[1]]] = i
  vert_subtriangles[tmp[2], vert_subtriangles_count[tmp[2]]] = i
  vert_subtriangles[tmp[3], vert_subtriangles_count[tmp[3]]] = i
}
max(vert_subtriangles_count)
vert_subtriangles = vert_subtriangles[,1:max(vert_subtriangles_count)]

# calcolo dell'intorno di ogni triangolo (ovvero degli indici dei triangoli geodeticamente vicini ad ogni triangolo)
geod_neigh_list = list()
geo_neigh_count = rep(NA, nrow(triangles_tmp))
con = file(description = "neigh_verts_list_pert.txt", open = "r")
for(tr_id in 1:nrow(triangles)){
  line = readLines(con = con, n = 1, ok = TRUE, warn = TRUE, encoding = "unknown", skipNul = FALSE)
  tmp = strsplit(line, " ")[[1]]
  
  if(length(tmp) > 0){
    neigh_verts = strtoi(tmp) + 1 #indici in un certo range dal triangolo 'tr_id' (+1 perchè viene da c++)
    # ft_id_tmp = unique(verts_tmp[neigh_verts,] %>% pull(father_id))
    # neigh_verts = unique(unlist(father_to_subverts_list[ft_id_tmp]))
    
    geod_dist_tmp = distances(subgraph(graph, v = as.character(neigh_verts)), v = as.character(father_to_subverts_list[[tr_id]]), to = as.character(neigh_verts))
    
    subtr_id = father_to_subtriangles_list[[tr_id]]
    
    for(i in 1:length(subtr_id)){
      neigh_tmp = unique(as.vector(vert_subtriangles[neigh_verts[colSums(geod_dist_tmp[match(as.numeric(triangles_tmp[subtr_id[i],-4]), father_to_subverts_list[[tr_id]]),] < cutoff_range+1e-10) > 0],]))
      neigh_tmp = neigh_tmp[!is.na(neigh_tmp)]
      
      geod_neigh_list[[subtr_id[i]]] = neigh_tmp
      geo_neigh_count[subtr_id[i]] = length(neigh_tmp)
    }
  }else{
    print("problema")
  }
  
  if(tr_id%%10 == 0) print(tr_id)
}
close(con)


con = file(description = "neigh_verts_list_pert.txt", open = "r")
for(tr_id in 1:10){
  line = readLines(con = con, n = 1, ok = TRUE, warn = TRUE, encoding = "unknown", skipNul = FALSE)
  tmp = strsplit(line, " ")[[1]]
  
  if(length(tmp) > 0){
    neigh_verts = strtoi(tmp) + 1 #indici in un certo range dal triangolo 'tr_id' (+1 perchè viene da c++)
    # ft_id_tmp = unique(verts_tmp[neigh_verts,] %>% pull(father_id))
    # neigh_verts = unique(unlist(father_to_subverts_list[ft_id_tmp]))
    
    subtr_id = father_to_subtriangles_list[[tr_id]]
    subtr_id_to = unique(as.vector(vert_subtriangles[neigh_verts,]))
    subtr_id_to = subtr_id_to[!is.na(subtr_id_to)]
  }else{
    print("problema")
  }
  
  if(tr_id%%100 == 0) print(tr_id)
}
close(con)


options(scipen = 999)

write.table(matrix(geo_neigh_count, ncol=1), "CAD DATASETS/geod_neig_tr_count.txt", row.names=FALSE, col.names=FALSE, sep=",")
con = file(description = "CAD DATASETS/geod_neig_tr_list.txt", open = "w")
for(i in 1:length(geod_neigh_list)){
  writeLines(text = strcat(as.character(geod_neigh_list[[i]] - 1), collapse = ", "), con = con)

  if(i%%10000 == 0) print(i)
}
close(con)

write.table(matrix(father_to_subtriangles_count, ncol=1), "CAD DATASETS/child_tr_count.txt", row.names=FALSE, col.names=FALSE, sep=",")
con = file(description = "CAD DATASETS/child_tr_list.txt", open = "w")
for(i in 1:length(father_to_subtriangles_list)){
  writeLines(text = strcat(as.character(father_to_subtriangles_list[[i]] - 1), collapse = ", "), con = con)

  if(i%%10000 == 0) print(i)
}
close(con)

options(scipen = 5)

