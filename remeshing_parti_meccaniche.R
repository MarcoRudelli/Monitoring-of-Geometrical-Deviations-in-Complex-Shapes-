
# install.packages("sfdct")
library(sfdct)
library(dplyr)
library(sf)
library(readobj)
#install.packages("rlist")
library(rlist)
library(igraph)
library(pracma)
library(plotly)

uovo_cad = read.obj("Acceptable.obj")

verts = t(uovo_cad$shapes[[1]]$positions)
triangles = t(uovo_cad$shapes[[1]]$indices + 1)


tr_areas = rep(NA, nrow(triangles))
cross <- function(a, b) {
  return(c(a[2] * b[3] - a[3] * b[2],
           a[3] * b[1] - a[1] * b[3],
           a[1] * b[2] - a[2] * b[1]))
}

for(j in 1:nrow(triangles)){
  triangle = verts[triangles[j,], ]
  
  v1 = triangle[2,]-triangle[1,]
  v2 = triangle[3,]-triangle[1,]
  
  tr_areas[j] = sqrt(sum(cross(v1, v2)^2))/2
}

mean(tr_areas)

triangles_tibble = tibble(V1 = triangles[,1], V2 = triangles[,2], V3 = triangles[,3], id = 1:nrow(triangles))

vert_triangles = matrix(NA, nrow=nrow(verts), ncol=12)
vert_triangles_count = rep(0, nrow(verts))
for(i in 1:nrow(verts)){
  tmp = which(triangles[,1] == i | triangles[,2] == i | triangles[,3] == i)
  vert_triangles[i, 1:length(tmp)] = tmp
  vert_triangles_count[i] = length(tmp)
}
max(vert_triangles_count)


min_num_triangles = 0.5
res_list_triangles_tmp = list()
res_list_verts_tmp = list()
res_list_is_edge = list()

prec_num_verts = 0

for(i_tmp in 1:nrow(triangles)){
  triangle_pos = verts[triangles[i_tmp,],]
  triangle_pos = triangle_pos - t(t(rep(1,3))) %*% t(triangle_pos[1,])
  
  v = cross(triangle_pos[2,], triangle_pos[3,])
  
  theta = acos(dot(v, c(0, 0, 1))/(sqrt(sum(v^2))))
  
  axis = cross(v, c(0, 0, 1))
  tmp = sqrt(sum(axis^2))
  if(tmp>0){
    axis = axis/tmp
  }
  
  K = rbind(c(0, -axis[3], axis[2]),
            c(axis[3], 0, -axis[1]),
            c(-axis[2], axis[1], 0))
  R = diag(3) + sin(theta)*K + (1-cos(theta))*K%*%K
  
  tmp.coord = (t(R %*% t(triangle_pos)))[,-3]
  b0 = st_polygon(list(rbind(c(0,0), tmp.coord[-1,], c(0,0))))
  
  tr = ct_triangulate(b0, a=mean(tr_areas)/min_num_triangles, D = TRUE)
  #plot(tr, col = "grey")
  
  verts_tmp_all = tibble(X = rep(NA, 3*length(tr)), Y = rep(NA, 3*length(tr)))
  for(i in 1:length(tr)){
    verts_tmp_all[((i-1)*3+1):(i*3),] = st_coordinates(tr[[i]])[-4,1:2]
  }
  verts_tmp = distinct(verts_tmp_all)
  verts_tmp = mutate(verts_tmp, id = 1:nrow(verts_tmp))
  triangles_tmp = matrix(left_join(verts_tmp_all, verts_tmp, by = join_by(X, Y))$id, ncol=3, byrow=TRUE)
  
  res_list_triangles_tmp[[i_tmp]] = tibble(V1 = triangles_tmp[,1], V2 = triangles_tmp[,2], V3 = triangles_tmp[,3]) + prec_num_verts
  verts_tmp_res = t(solve(R) %*% t(cbind(as.matrix(verts_tmp[,-3]), 0))) + t(t(rep(1,nrow(verts_tmp)))) %*% t(verts[triangles[i_tmp,][1],])
  prec_num_verts = prec_num_verts + nrow(verts_tmp_res)
  res_list_verts_tmp[[i_tmp]] = tibble(X = verts_tmp_res[,1], Y = verts_tmp_res[,2], Z = verts_tmp_res[,3])
  
  #is_edge_1, is_edge_2, is_edge_3
  A = tmp.coord[1,2] - tmp.coord[2,2]
  B = tmp.coord[2,1] - tmp.coord[1,1]
  C = tmp.coord[1,1] * tmp.coord[2,2] - tmp.coord[2,1] * tmp.coord[1,2]
  
  is_edge_1 = abs(A*verts_tmp$X + B*verts_tmp$Y + C) < 1e-10
  
  A = tmp.coord[2,2] - tmp.coord[3,2]
  B = tmp.coord[3,1] - tmp.coord[2,1]
  C = tmp.coord[2,1] * tmp.coord[3,2] - tmp.coord[3,1] * tmp.coord[2,2]
  
  is_edge_2 = abs(A*verts_tmp$X + B*verts_tmp$Y + C) < 1e-10
  
  A = tmp.coord[1,2] - tmp.coord[3,2]
  B = tmp.coord[3,1] - tmp.coord[1,1]
  C = tmp.coord[1,1] * tmp.coord[3,2] - tmp.coord[3,1] * tmp.coord[1,2]
  
  is_edge_3 = abs(A*verts_tmp$X + B*verts_tmp$Y + C) < 1e-10
  
  res_list_is_edge[[i_tmp]] = tibble(E1 = is_edge_1, E2 = is_edge_2, E3 = is_edge_3)
  
  if(i_tmp%%100==0) print(i_tmp)
}

res_list_verts_tmp[[1]]
res_list_triangles_tmp[[1]]
res_list_is_edge[[1]]

verts_tmp_all = bind_rows(res_list_verts_tmp)
triangles_tmp_all = bind_rows(res_list_triangles_tmp)
is_edge_all = bind_rows(res_list_is_edge)

father_id = rep(NA, nrow(verts_tmp_all))
ind_tmp = 0
for(i in 1:length(res_list_verts_tmp)){
  father_id[(ind_tmp+1):(ind_tmp+nrow(res_list_verts_tmp[[i]]))] = i

  ind_tmp = ind_tmp + nrow(res_list_verts_tmp[[i]])
}

verts_tmp_all = mutate(verts_tmp_all, father_id = father_id, E1 = is_edge_all$E1, E2 = is_edge_all$E2, E3 = is_edge_all$E3, row_id = 1:nrow(verts_tmp_all))
sum(rowSums(verts_tmp_all[,5:7])==2)/3

verts_tmp_all_sort_x = verts_tmp_all %>% select(X, row_id)
verts_tmp_all_sort_y = verts_tmp_all %>% select(Y, row_id)
verts_tmp_all_sort_z = verts_tmp_all %>% select(Z, row_id)
verts_tmp_all_sort_x = arrange(verts_tmp_all_sort_x, X)
verts_tmp_all_sort_y = arrange(verts_tmp_all_sort_y, Y)
verts_tmp_all_sort_z = arrange(verts_tmp_all_sort_z, Z)

verts_tmp_all_sort_x = mutate(verts_tmp_all_sort_x, group = cumsum(!c(FALSE, verts_tmp_all_sort_x$X[-1]-verts_tmp_all_sort_x$X[-length(verts_tmp_all_sort_x$X)] < 1e-10)))
verts_tmp_all_sort_y = mutate(verts_tmp_all_sort_y, group = cumsum(!c(FALSE, verts_tmp_all_sort_y$Y[-1]-verts_tmp_all_sort_y$Y[-length(verts_tmp_all_sort_y$Y)] < 1e-10)))
verts_tmp_all_sort_z = mutate(verts_tmp_all_sort_z, group = cumsum(!c(FALSE, verts_tmp_all_sort_z$Z[-1]-verts_tmp_all_sort_z$Z[-length(verts_tmp_all_sort_z$Z)] < 1e-10)))

verts_tmp_all = left_join(verts_tmp_all, verts_tmp_all_sort_x %>% select(row_id, group) %>% rename(group_x = group), join_by(row_id))
verts_tmp_all = left_join(verts_tmp_all, verts_tmp_all_sort_y %>% select(row_id, group) %>% rename(group_y = group), join_by(row_id))
verts_tmp_all = left_join(verts_tmp_all, verts_tmp_all_sort_z %>% select(row_id, group) %>% rename(group_z = group), join_by(row_id))

verts_tmp_all_grp = group_by(verts_tmp_all, group_x, group_y, group_z) %>% mutate(grp_id = cur_group_id())

verts_tmp = verts_tmp_all_grp %>% dplyr::filter(row_number()==which.max(sum(c(E1, E2, E3)))) %>% ungroup() %>% select(X, Y, Z, father_id, E1, E2, E3, grp_id)
sum(rowSums(verts_tmp[,5:7])==2)

nrow(verts_tmp_all) - nrow(verts_tmp)
verts_tmp = mutate(verts_tmp, row_id_unique = 1:nrow(verts_tmp))

ind_map = left_join(verts_tmp_all_grp, verts_tmp, join_by(grp_id)) %>% ungroup() %>% select(row_id, row_id_unique)
triangles_tmp = left_join(left_join(left_join(triangles_tmp_all, ind_map, join_by(V1==row_id)), ind_map, join_by(V2==row_id)), ind_map, join_by(V3==row_id))[,4:6] %>% rename(V1 = row_id_unique.x, V2 = row_id_unique.y, V3 = row_id_unique)

father_id_triangles = rep(NA, nrow(triangles_tmp))
ind_tmp = 0
for(i in 1:length(res_list_triangles_tmp)){
  father_id_triangles[(ind_tmp+1):(ind_tmp+nrow(res_list_triangles_tmp[[i]]))] = i
  
  ind_tmp = ind_tmp + nrow(res_list_triangles_tmp[[i]])
}
triangles_tmp = mutate(triangles_tmp, father_id = father_id_triangles)

nrow(triangles_tmp)

verts_overs = as.matrix(verts_tmp %>% select(X, Y, Z))
triangles_overs = as.matrix(triangles_tmp %>% select(V1, V2, V3))

tr_areas_overs = rep(NA, nrow(triangles_overs))
tr_centroids_overs = matrix(NA, ncol=3, nrow=nrow(triangles_overs))

for(j in 1:nrow(triangles_overs)){
  triangle = verts_overs[triangles_overs[j,], ]
  
  v1 = triangle[2,]-triangle[1,]
  v2 = triangle[3,]-triangle[1,]
  
  tr_areas_overs[j] = sqrt(sum(cross(v1, v2)^2))/2
  tr_centroids_overs[j,] = colMeans(triangle)
}

boxplot(tr_areas_overs, tr_areas)
abline(h=mean(tr_areas/min_num_triangles), col="red")


x <- tr_centroids_overs[,1]
y <- tr_centroids_overs[,2]
z <- tr_centroids_overs[,3]


plot_ly(x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'markers',
        marker = list(size = 1))


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

# con = file(description = "neigh_verts_list.txt", open = "r")
# for(tr_id in 1:1604){
#   line = readLines(con = con, n = 1, ok = TRUE, warn = TRUE, encoding = "unknown", skipNul = FALSE)
#   tmp = strsplit(line, " ")[[1]]
#   
#   
#   if(length(tmp) > 0){
#     neigh_verts = strtoi(tmp) + 1 #indici in un certo range dal triangolo 'tr_id' (+1 perchè viene da c++)
#   }else{
#     print("problema")
#   }
# }
# close(con)

x <- as.matrix(verts_tmp %>% select(X, Y, Z))[,1]
y <- as.matrix(verts_tmp %>% select(X, Y, Z))[,2]
z <- as.matrix(verts_tmp %>% select(X, Y, Z))[,3]

col = as.numeric(1:nrow(verts_tmp) %in% father_to_subverts_list[[112]]) + 1

plot_ly(x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'markers',
        marker = list(size = 5, color = col))


BATCH_LENGTH = 10000
floor(nrow(verts_tmp) / BATCH_LENGTH) + 1

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

unique_con_mat_complete = unique(con_mat_complete)
nrow(unique_con_mat_complete)


con = file(description = "neigh_verts_list.txt", open = "r")
for(tr_id in 1:1405){
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
        marker = list(size = 3, color = col))
k=k+1

#sigma2_boot = 0.01
sqrt(qchisq(0.999999, df = 2)) * sqrt(0.04)
cutoff_range = 1

unique_con_mat_complete_cutoff = unique_con_mat_complete[unique_con_mat_complete[,3] < cutoff_range + 1e-10,]
self_loops = cbind(1:nrow(verts_tmp), 1:nrow(verts_tmp), 0)
unique_con_mat_complete_cutoff = rbind(as.matrix(unique_con_mat_complete_cutoff), self_loops)
nrow(unique_con_mat_complete_cutoff)

con_mat = cbind(as.character(unique_con_mat_complete_cutoff[,1]), as.character(unique_con_mat_complete_cutoff[,2]))
graph = graph_from_edgelist(con_mat, directed = FALSE)
edge.attributes(graph)$weight = unique_con_mat_complete_cutoff[,3]


con = file(description = "neigh_verts_list.txt", open = "r")
for(tr_id in 1:112){
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
colors = as.numeric(geod_dist_tmp[k,] < cutoff_range+1e-10)+1
colors[neigh_verts == father_to_subverts_list[[tr_id]][k]] = 3
# colors = geod_dist_tmp[k,]
# colors = colorRamp(c("green", "red"))((colors-min(colors))/(max(colors)-min(colors)))
# colors <- rgb(colors/255)

plot_ly(x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'markers',
        marker = list(size = 4, color = colors))
k=k+1

nrow(verts_tmp)
x <- as.matrix(verts_tmp %>% select(X, Y, Z))[father_to_subverts_list[[tr_id]],1]
y <- as.matrix(verts_tmp %>% select(X, Y, Z))[father_to_subverts_list[[tr_id]],2]
z <- as.matrix(verts_tmp %>% select(X, Y, Z))[father_to_subverts_list[[tr_id]],3]

plot_ly(x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'markers',
        marker = list(size = 3))



father_to_subtriangles_count = as.vector((triangles_tmp %>% group_by(father_id) %>% count())[,2])$n
father_to_subtriangles_list = list()
cur_sum_tmp = 0
for(i in 1:nrow(triangles)){
  father_to_subtriangles_list[[i]] = (cur_sum_tmp + 1):(cur_sum_tmp + father_to_subtriangles_count[i])
  cur_sum_tmp = cur_sum_tmp + father_to_subtriangles_count[i]
}
father_to_subtriangles_list


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

geod_neigh_list = list()
geo_neigh_count = rep(NA, nrow(triangles_tmp))
con = file(description = "neigh_verts_list.txt", open = "r")
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

object.size(geod_neigh_list)/1024/1024

con = file(description = "neigh_verts_list.txt", open = "r")

for(tr_id in 1:112){
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

x <- tr_centroids_overs[subtr_id_to,1]
y <- tr_centroids_overs[subtr_id_to,2]
z <- tr_centroids_overs[subtr_id_to,3]

k=1
subtr_id[k]
colors = as.numeric(subtr_id_to %in% geod_neigh_list[[subtr_id[k]]]) + 1
colors[subtr_id_to == subtr_id[k]] = 3

plot_ly(x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'markers',
        marker = list(size = 4, color = colors))
k=k+1




options(scipen = 999)
# write.table(verts, "CAD_verts.csv", row.names=FALSE, col.names=FALSE, sep=",")
# write.table(triangles - 1, "CAD_triangles.csv", row.names=FALSE, col.names=FALSE, sep=",")
# 
# write.table(as.matrix(verts_tmp %>% select(X, Y, Z)), "CAD_verts_overs.csv", row.names=FALSE, col.names=FALSE, sep=",")
# write.table(as.matrix(triangles_tmp %>% select(V1, V2, V3) - 1), "CAD_triangles_overs.csv", row.names=FALSE, col.names=FALSE, sep=",")
# write.table(unique_con_mat_complete, "unique_con_mat_complete.csv", row.names=FALSE, col.names=FALSE, sep=",")
# 
# write.table(matrix(geo_neigh_count, ncol=1), "geod_neig_tr_count.txt", row.names=FALSE, col.names=FALSE, sep=",")
# con = file(description = "geod_neig_tr_list.txt", open = "w")
# for(i in 1:length(geod_neigh_list)){
#   writeLines(text = strcat(as.character(geod_neigh_list[[i]] - 1), collapse = ", "), con = con)
# 
#   if(i%%10000 == 0) print(i)
# }
# close(con)
# 
# write.table(matrix(father_to_subtriangles_count, ncol=1), "child_tr_count.txt", row.names=FALSE, col.names=FALSE, sep=",")
# con = file(description = "child_tr_list.txt", open = "w")
# for(i in 1:length(father_to_subtriangles_list)){
#   writeLines(text = strcat(as.character(father_to_subtriangles_list[[i]] - 1), collapse = ", "), con = con)
# 
#   if(i%%10000 == 0) print(i)
# }
# close(con)
# 
# options(scipen = 5)

# unique_con_mat_complete = read.csv("unique_con_mat_complete.csv", header = FALSE)

nrow(triangles_tmp)
