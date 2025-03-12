
library(sfdct)
library(dplyr)
library(sf)
library(readobj)
library(rlist)
library(igraph)
library(pracma)
library(plotly)
library(timeSeries)
library(BBmisc)

NUM_SIMS = 10
NUM_PARTS = 150
m0 = 100

verts_father = as.matrix(read.csv("CAD DATASETS/CAD_verts_old.csv", header = FALSE))
triangles_father = as.matrix(read.csv("CAD DATASETS/CAD_triangles_old.csv", header = FALSE) + 1)

verts = as.matrix(read.csv("CAD DATASETS/CAD_verts_overs.csv", header = FALSE))
triangles = as.matrix(read.csv("CAD DATASETS/CAD_triangles_overs.csv", header = FALSE) + 1)

unique_con_mat_complete = read.csv("unique_con_mat_complete.csv", header = FALSE)
nrow(unique_con_mat_complete)

sqrt(qchisq(0.999, df = 2)) * 1.6
cutoff_smoothing_range = 6

unique_con_mat_complete_cutoff = unique_con_mat_complete[unique_con_mat_complete[,3] < cutoff_smoothing_range+1e-10,]
self_loops = cbind(1:nrow(verts), 1:nrow(verts), 0)
unique_con_mat_complete_cutoff = rbind(unique_con_mat_complete_cutoff, self_loops)
nrow(unique_con_mat_complete_cutoff)

con_mat = cbind(as.character(unique_con_mat_complete_cutoff[,1]), as.character(unique_con_mat_complete_cutoff[,2]))
graph = graph_from_edgelist(con_mat, directed = FALSE)
edge.attributes(graph)$weight = unique_con_mat_complete_cutoff[,3]

father_to_subtriangles_count = read.table("CAD DATASETS/child_tr_count.txt", header = FALSE)[,1]
father_to_subtriangles_list = list()
cur_sum_tmp = 0
for(i in 1:length(father_to_subtriangles_count)){
  father_to_subtriangles_list[[i]] = (cur_sum_tmp + 1):(cur_sum_tmp + father_to_subtriangles_count[i])
  cur_sum_tmp = cur_sum_tmp + father_to_subtriangles_count[i]
}

vert_subtriangles = matrix(NA, nrow=nrow(verts), ncol=12)
vert_subtriangles_count = rep(0, nrow(verts))
for(i in 1:nrow(triangles)){
  tmp = as.numeric(triangles[i,c(1,2,3)])
  vert_subtriangles_count[tmp] = vert_subtriangles_count[tmp] + 1
  vert_subtriangles[tmp[1], vert_subtriangles_count[tmp[1]]] = i
  vert_subtriangles[tmp[2], vert_subtriangles_count[tmp[2]]] = i
  vert_subtriangles[tmp[3], vert_subtriangles_count[tmp[3]]] = i
}
max(vert_subtriangles_count)
vert_subtriangles = vert_subtriangles[,1:max(vert_subtriangles_count)]

tr_centroids = matrix(NA, nrow=nrow(triangles), ncol=3)
for(i in 1:nrow(triangles)){
  tr_centroids[i,] = colMeans(verts[triangles[i,], ])
}

cross <- function(a, b) {
  return(c(a[2] * b[3] - a[3] * b[2],
           a[3] * b[1] - a[1] * b[3],
           a[1] * b[2] - a[2] * b[1]))
}

triangles_areasCAD = rep(NA, nrow(triangles))
for(j in 1:nrow(triangles)){
  triangle = verts[triangles[j,],]
  
  v1 = triangle[2,]-triangle[1,]
  v2 = triangle[3,]-triangle[1,]
  
  triangles_areasCAD[j] = sqrt(sum(cross(v1, v2)^2))/2
}

smoothing_list = list()

smooth_sd = 1.6
curve(dnorm(x, 0, smooth_sd), from = -2-cutoff_smoothing_range, to = 2+cutoff_smoothing_range, n=10000)
abline(v=c(-cutoff_smoothing_range, cutoff_smoothing_range))

# calcolo pesi per smoothing geodetico
con = file(description = "neigh_verts_list_smooth.txt", open = "r")
for(tr_id in 1:nrow(triangles_father)){
  line = readLines(con = con, n = 1, ok = TRUE, warn = TRUE, encoding = "unknown", skipNul = FALSE)
  tmp = strsplit(line, " ")[[1]]
  
  if(length(tmp) > 0){
    neigh_verts = strtoi(tmp) + 1 #indici in un certo range dal triangolo 'tr_id' (+1 perchè viene da c++)
    
    subtr_id = father_to_subtriangles_list[[tr_id]]
    vert_ids = unique(as.vector(triangles[subtr_id,]))
    geod_dist_tmp = distances(subgraph(graph, v = as.character(neigh_verts)), v = as.character(vert_ids), to = as.character(neigh_verts))
    
    for(i in 1:length(subtr_id)){
      cmins = colMins(geod_dist_tmp[match(as.numeric(triangles[subtr_id[i],, drop=FALSE]), vert_ids),])
      gaussian_w = dnorm(cmins, 0, smooth_sd)[cmins < cutoff_smoothing_range+1e-10]
      tbl_tmp = tibble(cur_tr_id = as.vector(vert_subtriangles[neigh_verts[cmins < cutoff_smoothing_range+1e-10],]),
                       cur_wh = rep(gaussian_w, ncol(vert_subtriangles))) %>% distinct() %>% group_by(cur_tr_id) %>% summarise(m_w = max(cur_wh))
      tbl_tmp = tbl_tmp[!is.na(tbl_tmp[,1]),]
      
      # neigh_tmp = unique(as.vector(vert_subtriangles[neigh_verts[colSums(geod_dist_tmp[match(as.numeric(triangles[subtr_id[i],]), vert_ids),] < cutoff_smoothing_range+1e-10) > 0],]))
      # neigh_tmp = neigh_tmp[!is.na(neigh_tmp)]
      neigh_tmp = tbl_tmp %>% pull(cur_tr_id)
      
      #wh.mat = matrix((tbl_tmp %>% pull(m_w)) / sum(tbl_tmp %>% pull(m_w)), ncol=nrow(tbl_tmp), nrow=nrow(area_warp_mat), byrow = TRUE)
      
      #area_warp_mat_smooth[,subtr_id[i]] = rowSums(area_warp_mat[,neigh_tmp] * wh.mat)
      #back_dist_new_2_mat_smooth[,subtr_id[i]] = rowSums(back_dist_new_2_mat[,neigh_tmp] * area_warp_mat[,neigh_tmp] * wh.mat, na.rm = TRUE) / rowSums(area_warp_mat[,neigh_tmp] * wh.mat)
      
      #fwd_count_smooth[,subtr_id[i]] = rowSums(fwd_count[,neigh_tmp] * wh.mat)
      #forward_dist_smooth[,subtr_id[i]] = rowSums(forward_dist[,neigh_tmp] * fwd_count[,neigh_tmp] * wh.mat, na.rm = TRUE) / rowSums(fwd_count[,neigh_tmp] * wh.mat)
      smoothing_list[[subtr_id[i]]] = rbind(neigh_tmp, tbl_tmp %>% pull(m_w) / dnorm(0, 0, smooth_sd)) #/ sum(tbl_tmp %>% pull(m_w))
    }
  }else{
    print("problema")
  }
  
  if(tr_id%%5000 == 0) print(tr_id)
}
close(con)

dir.create("area_warp_mats", showWarnings = FALSE)
dir.create("backward_dist_mats", showWarnings = FALSE)
dir.create("forward_dist_mats", showWarnings = FALSE)
for(jjj in 1:NUM_SIMS){
  s = ((jjj-1)*NUM_PARTS+1):(jjj*NUM_PARTS)
  s1 = s[1:m0]
  s2 = s[(m0+1):NUM_PARTS]
  
  back_dist_new_2_mat = matrix(NA, nrow=NUM_PARTS, ncol=nrow(triangles))
  area_warp_mat = matrix(NA, nrow=NUM_PARTS, ncol=nrow(triangles))
  
  forward_dist = matrix(NA, nrow=NUM_PARTS, ncol=nrow(triangles))
  fwd_count = matrix(NA, nrow=NUM_PARTS, ncol=nrow(triangles))
  
  j=0
  for(i in c(s1, s2)){
    j=j+1
    res = scan(strcat(strcat("RESULTS BACKWARD/weights_", itostr(i)), ".txt"), what = double()) # * tot_areas_scan[j]/min(tot_areas_scan)
    res2 = scan(strcat(strcat("RESULTS BACKWARD/weighted_dists_", itostr(i)), ".txt"), what = double())
    
    res_f = scan(strcat(strcat("RESULTS FORWARD/weights_", itostr(i)), ".txt"), what = double()) # * tot_areas_scan[j]/min(tot_areas_scan)
    res2_f = scan(strcat(strcat("RESULTS FORWARD/weighted_dists_", itostr(i)), ".txt"), what = double())
    
    back_dist_new_2_mat[j,] = res2/res
    area_warp_mat[j,] = res
    forward_dist[j,] = res2_f/res_f
    fwd_count[j,] = res_f
    print(j)
  }
  
  back_dist_new_2_mat_smooth = back_dist_new_2_mat
  area_warp_mat_smooth = area_warp_mat
  
  forward_dist_smooth = forward_dist
  fwd_count_smooth = fwd_count

  for(tr_id in 1:ncol(back_dist_new_2_mat_smooth)){
    wh.mat = matrix(smoothing_list[[tr_id]][2,], ncol=length(smoothing_list[[tr_id]][2,]), nrow=nrow(area_warp_mat), byrow = TRUE)
    neigh_tmp = smoothing_list[[tr_id]][1,]

    area_warp_mat_smooth[,tr_id] = rowSums(area_warp_mat[,neigh_tmp] * wh.mat)
    back_dist_new_2_mat_smooth[,tr_id] = rowSums(back_dist_new_2_mat[,neigh_tmp] * area_warp_mat[,neigh_tmp] * wh.mat, na.rm = TRUE) / rowSums(area_warp_mat[,neigh_tmp] * wh.mat)

    fwd_count_smooth[,tr_id] = rowSums(fwd_count[,neigh_tmp] * wh.mat)
    forward_dist_smooth[,tr_id] = rowSums(forward_dist[,neigh_tmp] * fwd_count[,neigh_tmp] * wh.mat, na.rm = TRUE) / rowSums(fwd_count[,neigh_tmp] * wh.mat)
  }
  
  # rimozione eventuali NA
  while(sum(is.na(back_dist_new_2_mat_smooth))>0){
    bd_smooth_tmp = back_dist_new_2_mat_smooth
    for(tr_id in which(colSums( is.na(bd_smooth_tmp) ) > 0)){
      to_set = which(is.na(bd_smooth_tmp[,tr_id]))
      wh.mat = matrix(smoothing_list[[tr_id]][2,], ncol=length(smoothing_list[[tr_id]][2,]), nrow=nrow(area_warp_mat), byrow = TRUE)
      neigh_tmp = smoothing_list[[tr_id]][1,]
      
      back_dist_new_2_mat_smooth[to_set,tr_id] = rowSums(bd_smooth_tmp[to_set,neigh_tmp,drop=FALSE] * wh.mat[to_set,,drop=FALSE], na.rm = TRUE) / rowSums((!is.na(bd_smooth_tmp[to_set,neigh_tmp,drop=FALSE]))*wh.mat[to_set,,drop=FALSE])
    }
  }
  
  # campione casuale del 25% delle colonne
  sss = sample(1:ncol(back_dist_new_2_mat_smooth), ncol(back_dist_new_2_mat_smooth)/4)

  write.table(area_warp_mat_smooth[,sss], file = strcat(strcat("area_warp_mats/area_warp_mat_smooth_", itostr(jjj)), ".csv"), col.names = FALSE, row.names = FALSE, sep = ",")
  write.table(back_dist_new_2_mat_smooth[,sss], file = strcat(strcat("backward_dist_mats/backward_dist_smooth_", itostr(jjj)), ".csv"), col.names = FALSE, row.names = FALSE, sep = ",")
  write.table(forward_dist_smooth[,sss], file = strcat(strcat("forward_dist_mats/forward_dist_smooth_", itostr(jjj)), ".csv"), col.names = FALSE, row.names = FALSE, sep = ",")

  print(jjj)
}
