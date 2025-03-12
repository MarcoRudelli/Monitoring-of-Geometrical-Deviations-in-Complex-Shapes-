
library(sfdct)
library(dplyr)
library(sf)
library(readobj)
library(rlist)
library(igraph)
library(pracma)
library(plotly)

# carico mesh CAD originale
mesh_cad = read.obj("Acceptable.obj")

# estraggo vertici e triangoli
verts = t(mesh_cad$shapes[[1]]$positions)
triangles = t(mesh_cad$shapes[[1]]$indices + 1)

# calcolo aree triangoli mesh CAD originale
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

# creo tibble triangoli
triangles_tibble = tibble(V1 = triangles[,1], V2 = triangles[,2], V3 = triangles[,3], id = 1:nrow(triangles))

# creo matrice che, per ogni vertice, fornisce i triangoli di cui tale vertice fa parte
vert_triangles = matrix(NA, nrow=nrow(verts), ncol=15)
vert_triangles_count = rep(0, nrow(verts))
for(i in 1:nrow(verts)){
  tmp = which(triangles[,1] == i | triangles[,2] == i | triangles[,3] == i)
  vert_triangles[i, 1:length(tmp)] = tmp
  vert_triangles_count[i] = length(tmp)
}
max(vert_triangles_count)
vert_triangles = vert_triangles[,1:max(vert_triangles_count)]


# fisso il valore massimo per l'area dei triangoli
max_area = mean(tr_areas) / 0.5

# computazione mesh CAD con triangolazione 'aumentata'
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
  
  tr = ct_triangulate(b0, a = max_area, D = TRUE)
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

# ci sono vertici contati più volte
verts_tmp_all = mutate(verts_tmp_all, father_id = father_id, E1 = is_edge_all$E1, E2 = is_edge_all$E2, E3 = is_edge_all$E3, row_id = 1:nrow(verts_tmp_all))

# elimino vertici ridondanti
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

# vertici triangolazione aumentata senza ridondanze 
verts_tmp = verts_tmp_all_grp %>% dplyr::filter(row_number()==which.max(sum(c(E1, E2, E3)))) %>% ungroup() %>% select(X, Y, Z, father_id, E1, E2, E3, grp_id)

nrow(verts_tmp_all) - nrow(verts_tmp)
verts_tmp = mutate(verts_tmp, row_id_unique = 1:nrow(verts_tmp))

# calcolo triangolazione aumentata con gli indici dei vertici nuovi
ind_map = left_join(verts_tmp_all_grp, verts_tmp, join_by(grp_id)) %>% ungroup() %>% select(row_id, row_id_unique)
triangles_tmp = left_join(left_join(left_join(triangles_tmp_all, ind_map, join_by(V1==row_id)), ind_map, join_by(V2==row_id)), ind_map, join_by(V3==row_id))[,4:6] %>% rename(V1 = row_id_unique.x, V2 = row_id_unique.y, V3 = row_id_unique)
nrow(triangles_tmp)

# calcolo indici "sotto-triangoli" per ogni triangolo "padre" della triangolazione originale
father_id_triangles = rep(NA, nrow(triangles_tmp))
ind_tmp = 0
for(i in 1:length(res_list_triangles_tmp)){
  father_id_triangles[(ind_tmp+1):(ind_tmp+nrow(res_list_triangles_tmp[[i]]))] = i
  
  ind_tmp = ind_tmp + nrow(res_list_triangles_tmp[[i]])
}
triangles_tmp = mutate(triangles_tmp, father_id = father_id_triangles)


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
abline(h=max_area, col="red")

x <- tr_centroids_overs[,1]
y <- tr_centroids_overs[,2]
z <- tr_centroids_overs[,3]
plot_ly(x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'markers',
        marker = list(size = 1))


options(scipen = 999)

dir.create("CAD DATASETS", showWarnings = FALSE)
write.table(verts, "CAD DATASETS/CAD_verts_old.csv", row.names=FALSE, col.names=FALSE, sep=",")
write.table(triangles - 1, "CAD DATASETS/CAD_triangles_old.csv", row.names=FALSE, col.names=FALSE, sep=",")
write.table(as.matrix(verts_tmp %>% select(X, Y, Z)), "CAD DATASETS/CAD_verts_overs.csv", row.names=FALSE, col.names=FALSE, sep=",")
write.table(as.matrix(triangles_tmp %>% select(V1, V2, V3) - 1), "CAD DATASETS/CAD_triangles_overs.csv", row.names=FALSE, col.names=FALSE, sep=",")

options(scipen = 5)



