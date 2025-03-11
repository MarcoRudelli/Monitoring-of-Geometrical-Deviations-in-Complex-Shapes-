
verts_father = as.matrix(read.csv("CAD DATASETS/CAD_verts_old.csv", header = FALSE))
triangles_father = as.matrix(read.csv("CAD DATASETS/CAD_triangles_old.csv", header = FALSE) + 1)

verts = as.matrix(read.csv("CAD DATASETS/CAD_verts_overs.csv", header = FALSE))
triangles = as.matrix(read.csv("CAD DATASETS/CAD_triangles_overs.csv", header = FALSE) + 1)

unique_con_mat_complete = read.csv("unique_con_mat_complete.csv", header = FALSE)
nrow(unique_con_mat_complete)

# #sigma2_boot = 0.01
# sqrt(qchisq(0.999999, df = 2)) * sqrt(0.04)
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
father_to_subtriangles_list

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

con = file(description = "neigh_verts_list_6.txt", open = "r")
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

# ccc = rep(0, ncol(back_dist_new_2_mat))
# ccc[smoothing_list[[1]][1,]] = smoothing_list[[1]][2,]

#samp = 1:100#c(sample(1:1500),sample(1501:2250))
#write.table(samp, file="samp3.txt", col.names = FALSE, row.names = FALSE)


39*150
for(jjj in 1:10){
  # s.tmp = ((jjj-1)*150+1):(jjj*150)
  #s1 = ((jjj-1)*100+1):(jjj*100) #s.tmp[1:100] #((jjj-1)*100+1):(jjj*100) # samp[ ((jjj-1)*100+1):(jjj*100) ]
  #s2 = 10000 + ((jjj-1)*50+1):(jjj*50) #s.tmp[101:150] #10000 + ((jjj-1)*50+1):(jjj*50) # samp[ 10000 + ((jjj-1)*50+1):(jjj*50) ]
  s = ((jjj-1)*150+1):(jjj*150)
  s1 = s[1:100]
  s2 = s[101:150]
  
  back_dist_new_2_mat = matrix(NA, nrow=150, ncol=nrow(triangles))
  area_warp_mat = matrix(NA, nrow=150, ncol=nrow(triangles))
  
  forward_dist = matrix(NA, nrow=150, ncol=nrow(triangles))
  fwd_count = matrix(NA, nrow=150, ncol=nrow(triangles))
  
  j=0
  for(i in c(s1, s2)){
    j=j+1
    res = scan(strcat(strcat("/Volumes/EXTERNAL_USB/sim_2.51/RESULTS SIMULAZIONE_SIM2.51/weights_", itostr(i)), ".txt"), what = double()) # * tot_areas_scan[j]/min(tot_areas_scan)
    res2 = scan(strcat(strcat("/Volumes/EXTERNAL_USB/sim_2.51/RESULTS SIMULAZIONE_SIM2.51/weighted_dists_", itostr(i)), ".txt"), what = double())
    
    res_f = scan(strcat(strcat("/Volumes/EXTERNAL_USB/sim_2.51/RESULTS FORWARD SIMULAZIONE_SIM2.51/weights_", itostr(i)), ".txt"), what = double()) # * tot_areas_scan[j]/min(tot_areas_scan)
    res2_f = scan(strcat(strcat("/Volumes/EXTERNAL_USB/sim_2.51/RESULTS FORWARD SIMULAZIONE_SIM2.51/weighted_dists_", itostr(i)), ".txt"), what = double())
    
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
  
  # prova oversmoothing
  #####
  # for(sm_i in 1:10){
  #   back_dist_new_2_mat = back_dist_new_2_mat_smooth
  #   area_warp_mat = area_warp_mat_smooth
  # 
  #   forward_dist = forward_dist_smooth
  #   fwd_count = fwd_count_smooth
  # 
  #   for(tr_id in 1:ncol(back_dist_new_2_mat_smooth)){
  #     wh.mat = matrix(smoothing_list[[tr_id]][2,], ncol=length(smoothing_list[[tr_id]][2,]), nrow=nrow(area_warp_mat), byrow = TRUE)
  #     neigh_tmp = smoothing_list[[tr_id]][1,]
  # 
  #     area_warp_mat_smooth[,tr_id] = rowSums(area_warp_mat[,neigh_tmp] * wh.mat)
  #     back_dist_new_2_mat_smooth[,tr_id] = rowSums(back_dist_new_2_mat[,neigh_tmp] * area_warp_mat[,neigh_tmp] * wh.mat, na.rm = TRUE) / rowSums(area_warp_mat[,neigh_tmp] * wh.mat)
  # 
  #     fwd_count_smooth[,tr_id] = rowSums(fwd_count[,neigh_tmp] * wh.mat)
  #     forward_dist_smooth[,tr_id] = rowSums(forward_dist[,neigh_tmp] * fwd_count[,neigh_tmp] * wh.mat, na.rm = TRUE) / rowSums(fwd_count[,neigh_tmp] * wh.mat)
  #   }
  # 
  #   print(sm_i)
  # }
  #####
  
  while(sum(is.na(back_dist_new_2_mat_smooth))>0){
    bd_smooth_tmp = back_dist_new_2_mat_smooth
    for(tr_id in which(colSums( is.na(bd_smooth_tmp) ) > 0)){
      to_set = which(is.na(bd_smooth_tmp[,tr_id]))
      wh.mat = matrix(smoothing_list[[tr_id]][2,], ncol=length(smoothing_list[[tr_id]][2,]), nrow=nrow(area_warp_mat), byrow = TRUE)
      neigh_tmp = smoothing_list[[tr_id]][1,]
      
      back_dist_new_2_mat_smooth[to_set,tr_id] = rowSums(bd_smooth_tmp[to_set,neigh_tmp,drop=FALSE] * wh.mat[to_set,,drop=FALSE], na.rm = TRUE) / rowSums((!is.na(bd_smooth_tmp[to_set,neigh_tmp,drop=FALSE]))*wh.mat[to_set,,drop=FALSE])
    }
  }
  
  #####
 #  # num_pc = 3 
 #  #####
 #  # pcs = matrix(NA, ncol=num_pc*3, nrow=150)
 #  # sc_err_back_dist_new_2_mat_smooth = back_dist_new_2_mat_smooth
 #  # sc_err_area_warp_mat_smooth = area_warp_mat_smooth
 #  # sc_err_forward_dist_smooth = forward_dist_smooth
 #  # 
 #  # win_ids = 1:100
 #  # tmp1 = scale(back_dist_new_2_mat_smooth[win_ids,],
 #  #                                                     center = apply(back_dist_new_2_mat_smooth[win_ids,], 2, median),
 #  #                                                     scale = apply(back_dist_new_2_mat_smooth[win_ids,], 2, mad))
 #  # tmp2 = scale(area_warp_mat_smooth[win_ids,],
 #  #                                                     center = apply(area_warp_mat_smooth[win_ids,], 2, median),
 #  #                                                     scale = apply(area_warp_mat_smooth[win_ids,], 2, mad))
 #  # tmp3 = scale(forward_dist_smooth[win_ids,],
 #  #                                                     center = apply(forward_dist_smooth[win_ids,], 2, median),
 #  #                                                     scale = apply(forward_dist_smooth[win_ids,], 2, mad))
 #  # tmpa = prcomp(cbind(tmp1), center = FALSE, scale. = FALSE)
 #  # tmpb = prcomp(cbind(tmp2), center = FALSE, scale. = FALSE)
 #  # tmpc = prcomp(cbind(tmp3), center = FALSE, scale. = FALSE)
 #  # pcs[win_ids,] = tmp$x[,1:num_pc]
 #  # err = cbind(tmp1, tmp2, tmp3) - pcs[win_ids,] %*% t(tmp$rotation[,1:num_pc])
 #  # 
 #  # sc_err_back_dist_new_2_mat_smooth[win_ids,] = err[,1:ncol(tmp1)]
 #  # sc_err_area_warp_mat_smooth[win_ids,] = err[,(ncol(tmp1)+1):(ncol(tmp1)+ncol(tmp2))]
 #  # sc_err_forward_dist_smooth[win_ids,] = err[,(ncol(tmp1)+ncol(tmp2)+1):ncol(err)]
 #  # 
 #  # for(i in 101:150){
 #  #   win_ids = sort(c(sample(1:min(100, i-10), size = 90), (i-9):i))
 #  #   
 #  #   tmp1 = scale(back_dist_new_2_mat_smooth[win_ids,],
 #  #                center = apply(back_dist_new_2_mat_smooth[win_ids,], 2, median),
 #  #                scale = apply(back_dist_new_2_mat_smooth[win_ids,], 2, mad))
 #  #   tmp2 = scale(area_warp_mat_smooth[win_ids,],
 #  #                center = apply(area_warp_mat_smooth[win_ids,], 2, median),
 #  #                scale = apply(area_warp_mat_smooth[win_ids,], 2, mad))
 #  #   tmp3 = scale(forward_dist_smooth[win_ids,],
 #  #                center = apply(forward_dist_smooth[win_ids,], 2, median),
 #  #                scale = apply(forward_dist_smooth[win_ids,], 2, mad))
 #  #   tmpa = prcomp(cbind(tmp1), center = FALSE, scale. = FALSE)
 #  #   tmpb = prcomp(cbind(tmp2), center = FALSE, scale. = FALSE)
 #  #   tmpc = prcomp(cbind(tmp3), center = FALSE, scale. = FALSE)
 #  #   pcs[i,] = tmp$x[100, 1:num_pc]
 #  #   err = cbind(tmp1[100,, drop=FALSE], tmp2[100,, drop=FALSE], tmp3[100,, drop=FALSE]) - pcs[i,, drop=FALSE] %*% t(tmp$rotation[,1:num_pc])
 #  #   
 #  #   sc_err_back_dist_new_2_mat_smooth[i,] = err[,1:ncol(tmp1)]
 #  #   sc_err_area_warp_mat_smooth[i,] = err[,(ncol(tmp1)+1):(ncol(tmp1)+ncol(tmp2))]
 #  #   sc_err_forward_dist_smooth[i,] = err[,(ncol(tmp1)+ncol(tmp2)+1):ncol(err)]
 #  #   
 #  #   print(i)
 #  # }
 #  #####
 #  
 #  # pcs = matrix(NA, ncol=num_pc, nrow=150)
 #  # sc_err_back_dist_new_2_mat_smooth = back_dist_new_2_mat_smooth
 #  # sc_err_area_warp_mat_smooth = area_warp_mat_smooth
 #  # sc_err_forward_dist_smooth = forward_dist_smooth
 #  # for(i in 1:25){ # (i in 1:10){
 #  #   ids_to_set = ((i-1)*4+1):(i*4) # ((i-1)*10+1):(i*10)
 #  #   ids_to_set = c(ids_to_set, 100+((i-1)*2+1):(i*2)) # , 100+((i-1)*5+1):(i*5))
 #  #   ids_to_drop = c(((i-1)*4+1):(i*4), 101:150) # c(((i-1)*10+1):(i*10), 101:150)
 #  #   
 #  #   sc_area_warp_phase1 = scale(area_warp_mat_smooth[-ids_to_drop,]) #[1:100,])
 #  #   sc_area_warp_phase2 = scale(area_warp_mat_smooth[ids_to_set,], center = attr(sc_area_warp_phase1,"scaled:center"), scale = attr(sc_area_warp_phase1,"scaled:scale")) #[1:150,]
 #  #   sc_back_dist_phase1 = scale(back_dist_new_2_mat_smooth[-ids_to_drop,]) #[1:100,])
 #  #   sc_back_dist_phase2 = scale(back_dist_new_2_mat_smooth[ids_to_set,], center = attr(sc_back_dist_phase1,"scaled:center"), scale = attr(sc_back_dist_phase1,"scaled:scale")) #[1:150,]
 #  #   sc_fwd_dist_phase1 = scale(forward_dist_smooth[-ids_to_drop,]) #[1:100,])
 #  #   sc_fwd_dist_phase2 = scale(forward_dist_smooth[ids_to_set,], center = attr(sc_fwd_dist_phase1,"scaled:center"), scale = attr(sc_fwd_dist_phase1,"scaled:scale")) #[1:150,]
 #  # 
 #  #   tmp = prcomp(cbind(sc_area_warp_phase1, sc_back_dist_phase1, sc_fwd_dist_phase1), scale. = FALSE, center = FALSE)
 #  # 
 #  #   pcs_phase2 = cbind(sc_area_warp_phase2, sc_back_dist_phase2, sc_fwd_dist_phase2) %*% tmp$rotation[,1:num_pc]
 #  #   err = cbind(sc_area_warp_phase2, sc_back_dist_phase2, sc_fwd_dist_phase2) - pcs_phase2 %*% t(tmp$rotation[,1:num_pc])
 #  #   
 #  #   sc_err_area_warp_mat_smooth[ids_to_set,] = err[,1:ncol(sc_area_warp_phase2)]
 #  #   sc_err_back_dist_new_2_mat_smooth[ids_to_set,] = err[,(ncol(sc_area_warp_phase2)+1):(ncol(sc_area_warp_phase2)+ncol(sc_back_dist_phase2))]
 #  #   sc_err_forward_dist_smooth[ids_to_set,] = err[,(ncol(sc_area_warp_phase2)+ncol(sc_back_dist_phase2)+1):ncol(err)]
 #  #   pcs[ids_to_set,] = pcs_phase2
 #  # 
 #  #   print(i)
 #  # }
 #  
 #  #####

  area_warp_mat_norm = area_warp_mat_smooth
  back_dist_norm = back_dist_new_2_mat_smooth
  forw_dist_norm = forward_dist_smooth
  # area_warp_mat_sc = area_warp_mat_smooth
  # back_dist_sc = back_dist_new_2_mat_smooth
  # forw_dist_sc = forward_dist_smooth

  for(i in 1:150){ # (i in 1:10){
    # ids_to_set = ((i-1)*4+1):(i*4) # ((i-1)*10+1):(i*10)
    # ids_to_set = c(ids_to_set, 100+((i-1)*2+1):(i*2)) # , 100+((i-1)*5+1):(i*5))
    # ids_to_drop = c(((i-1)*4+1):(i*4), 101:150) # c(((i-1)*10+1):(i*10), 101:150)
    
    if(i <= 100) tmp_s = sample((1:100)[-i], 80)
    else tmp_s = sample((1:100), 80)
    
    ids_to_set = i

    tmp_seq = seq(from=0, to=1, length.out = length(tmp_s) + 3)
    tmp_seq = qnorm(tmp_seq[-c(1,length(tmp_seq))], 0, 1)
    for(j2 in ids_to_set){
      tmp_rk_1 = rep(0, ncol(area_warp_mat_norm))
      tmp_rk_2 = rep(0, ncol(back_dist_norm))
      tmp_rk_3 = rep(0, ncol(forw_dist_norm))
      for(j1 in (1:150)[tmp_s]){
        tmp_rk_1 = tmp_rk_1 + (area_warp_mat_smooth[j2,] > area_warp_mat_smooth[j1,])
        tmp_rk_1 = tmp_rk_1 + (area_warp_mat_smooth[j2,] == area_warp_mat_smooth[j1,])/2
        tmp_rk_2 = tmp_rk_2 + (back_dist_new_2_mat_smooth[j2,] > back_dist_new_2_mat_smooth[j1,])
        tmp_rk_3 = tmp_rk_3 + (forward_dist_smooth[j2,] > forward_dist_smooth[j1,])
      }
      area_warp_mat_norm[j2,] = tmp_seq[tmp_rk_1 + 1]
      back_dist_norm[j2,] = tmp_seq[tmp_rk_2 + 1]
      forw_dist_norm[j2,] = tmp_seq[tmp_rk_3 + 1]
    }

    # means.tmp1 = colMeans(area_warp_mat_smooth[tmp_s,])
    # sd.tmp1 = apply(area_warp_mat_smooth[tmp_s,], 2, sd)
    # means.tmp2 = colMeans(back_dist_new_2_mat_smooth[tmp_s,])
    # sd.tmp2 = apply(back_dist_new_2_mat_smooth[tmp_s,], 2, sd)
    # means.tmp3 = colMeans(forward_dist_smooth[tmp_s,])
    # sd.tmp3 = apply(forward_dist_smooth[tmp_s,], 2, sd)
    #
    # area_warp_mat_sc[ids_to_set,] = scale(area_warp_mat_smooth[ids_to_set,, drop=FALSE], center = means.tmp1, scale = sd.tmp1)
    # back_dist_sc[ids_to_set,] = scale(back_dist_new_2_mat_smooth[ids_to_set,, drop=FALSE], center = means.tmp2, scale = sd.tmp2)
    # forw_dist_sc[ids_to_set,] = scale(forward_dist_smooth[ids_to_set,, drop=FALSE], center = means.tmp3, scale = sd.tmp3)

    print(i)
  }
  
  # dim(forward_dist_smooth)
  # boxplot(t(forward_dist_smooth))
  # plot(apply(forw_dist_norm, 1, quantile, probs=0.99))

  for(i in 1:nrow(back_dist_norm)){
    tmp_max = pmax(back_dist_norm[i,], forw_dist_norm[i,])
    back_dist_norm[i,] = tmp_max
  }
  
  qt1 = t(apply(area_warp_mat_norm, 1, function(x) { c(mean(x), var(x), skewness(x), kurtosis(x)) })) #, quantile, probs = c(0.1,0.5,0.9)))
  qt2 = t(apply(back_dist_norm, 1, function(x) { c(mean(x), var(x), skewness(x), kurtosis(x)) })) #, quantile, probs = c(0.1,0.5,0.9)))
  #qt3 = t(apply(forw_dist_norm, 1, function(x) { c(mean(x), var(x), skewness(x), kurtosis(x)) })) #, quantile, probs = c(0.1,0.5,0.9)))
  results = cbind(qt1, qt2)

  write.table(results, file = strcat(strcat("/Volumes/EXTERNAL_USB/sim_2.51/results_SIM2.51_CV_moments/res_SIM2.51_", itostr(jjj)), ".csv"), col.names = FALSE, row.names = FALSE, sep = ",")

  qt1 = t(apply(area_warp_mat_norm, 1, quantile, probs = c(0.01,0.1,0.5,0.9,0.99)))
  qt2 = t(apply(back_dist_norm, 1, quantile, probs = c(0.01,0.1,0.5,0.9,0.99)))
  #qt3 = t(apply(forw_dist_norm, 1, quantile, probs = c(0.01,0.1,0.5,0.9,0.99)))
  results = cbind(qt1, qt2)

  write.table(results, file = strcat(strcat("/Volumes/EXTERNAL_USB/sim_2.51/results_SIM2.51_CV_quantiles/res_SIM2.51_", itostr(jjj)), ".csv"), col.names = FALSE, row.names = FALSE, sep = ",")
  #####
  
  
  sss = sample(1:ncol(back_dist_new_2_mat_smooth), ncol(back_dist_new_2_mat_smooth)/4)

  write.table(area_warp_mat_smooth[,sss], file = strcat(strcat("/Volumes/EXTERNAL_USB/sim_2.51/area_warp_smooth_mats/area_warp_mat_smooth_", itostr(jjj)), ".csv"), col.names = FALSE, row.names = FALSE, sep = ",")
  write.table(back_dist_new_2_mat_smooth[,sss], file = strcat(strcat("/Volumes/EXTERNAL_USB/sim_2.51/backward_dist_mats/backward_dist_smooth_", itostr(jjj)), ".csv"), col.names = FALSE, row.names = FALSE, sep = ",")
  write.table(forward_dist_smooth[,sss], file = strcat(strcat("/Volumes/EXTERNAL_USB/sim_2.51/forward_dist_mats/forward_dist_smooth_", itostr(jjj)), ".csv"), col.names = FALSE, row.names = FALSE, sep = ",")

  print(jjj)
}
length(sss)
# a sto punto... potrei aggiungere un indicatore di outlyingness (tipo con stime kernel) ed avere media, varianza, asimmetria, kurtosi e indicatore di outlyingness...
  # con l'idea di avere uno strumento molto generico e distribution free (quello dell'outlyingness) e quattro strumenti che diano potenza specifica a variazioni interessanti e note (media, varianza, ecc.)

# 
# boxplot(t(tttmp))
# boxplot(t(sc_area_warp))
# 
# wm = matrix(triangles_areasCAD, ncol=ncol(area_warp_mat_smooth), nrow=150, byrow=TRUE)
# tttmp = area_warp_mat_smooth/sum(area_warp_mat_smooth) - wm/sum(wm)
# # prova = apply(area_warp_mat_smooth, 2, rank)
# 
# boxplot(t(prova-50))
# boxplot(prova[3,]-50, prova[110,]-50)
# abline(h=0)
# plot(apply(prova, 1, sd))
# 
# plot(apply(prova, 1, function(x){mean(x)}))
# plot(apply(sc_err_area_warp_mat_smooth, 1, quantile, prob=0.5))
# abline(v=c(3,110))
# 
# boxplot(prova[2,], prova[141,])
# var(prova[141,])
# #sc_area_warp = scale(area_warp_mat_smooth)
# colors2 = colMeans(prova[2,, drop=FALSE])
# mean(prova[2,, drop=FALSE])
# #colors2 = colMeans(area_warp_mat_smooth[104,, drop=FALSE])
# #boxplot(colors, colors2)
# # ca = colors-colors2
# # cb = colors-colors2
# # boxplot(cb, ca)
# colors = colorRamp(c("green", "red"))((colors2-min(c(1)))/(max(c(150))-min(c(1))))
# colors <- rgb(colors/255)
# x = tr_centroids[,1]
# y = tr_centroids[,2]
# z = tr_centroids[,3]
# c = colors
# plot_ly(x = ~x, 
#         y = ~y, 
#         z =~ z, type = 'scatter3d', mode = 'markers',
#         marker = list(size = 3, color=c))
# 
# 

kkk=1
#plot(results[,kkk], type="b")
# boxplot(pcs[1:100,kkk], pcs[101:150,kkk])
plot(results[,kkk], type="b")
rsp(results[,kkk])
#plot(results2[,kkk])
#t.test(results2[1:100,kkk], results2[101:150,kkk])
kkk=kkk+1
# library("dfphase1")


rsp(apply(area_warp_mat, 1, mean))


area_warp_mat_smooth_loglik = area_warp_mat_smooth
for(i in 1:nrow(triangles)){
  for(j in 1:10){
    ids_to_set = ((j-1)*10+1):(j*10)
    ids_to_set = c(ids_to_set, 100+((j-1)*5+1):(j*5))
    ids_to_drop = c(((j-1)*10+1):(j*10), 101:150)
    
    sd.tmp = bw.nrd(area_warp_mat_smooth[-ids_to_drop, i])
    
    for(j2 in 1:length(ids_to_set)){
      area_warp_mat_smooth_loglik[ids_to_set[j2], i] = sum( dnorm(area_warp_mat_smooth[ids_to_set[j2], i], area_warp_mat_smooth[-ids_to_drop, i], sd.tmp, log = TRUE) )
    }
  }
 print(i)
}


sc_area_warp = area_warp_mat_smooth
for(i in 1:100){
  if(i <= 50){
    ids_to_set = c(i, i+100)
  }
  else{
    ids_to_set = i
  }
  ids_to_drop = c(i, 101:150)
  
  means.tmp = colMeans(area_warp_mat_smooth[-ids_to_drop,])
  sd.tmp = apply(area_warp_mat_smooth[-ids_to_drop,], 2, sd)
  
  sc_area_warp[ids_to_set,] = scale(area_warp_mat_smooth[ids_to_set,, drop=FALSE], center = means.tmp, scale = sd.tmp)
  print(i)
}

ks_stat = rep(NA, 150)
rg = (1:999)/1000
for(i in 1:100){
  if(i <= 50){
    ids_to_set = c(i, i+100)
  }
  else{
    ids_to_set = i
  }
  ids_to_drop = c(i, 101:150)
  
  tmp = quantile(as.numeric(sc_area_warp[-ids_to_drop,]), probs=rg)
  #tmp = as.numeric(sc_area_warp[-ids_to_drop,])
  for(j in ids_to_set){
    # plot((1:999)/1000, quantile(tmp, probs = (1:999)/1000), type="l")
    # lines((1:999)/1000, quantile(sc_area_warp[90,], probs = (1:999)/1000), col="red")
    # abline(v=range(rg))

    ks_stat[j] = sqrt(sum( (tmp-quantile(sc_area_warp[j,], probs=rg))^2)) #max(abs(tmp-ecdf(sc_area_warp[j,])(rg)))
   
    # o_rk = rank(c(tmp, sc_area_warp[j,]))
    # rx = o_rk[1:length(tmp)]
    # ry = o_rk[(length(tmp)+1):length(o_rk)]
    #   
    # N = length(rx)
    # M = length(ry)
    # 
    # U = ( N * sum( (rx - 1:N)^2 ) + M * sum( (ry - 1:M)^2 ) ) / (N + M)
    #ks_stat[j] = mean(sc_area_warp[j,])
  }
  
  print(i)
}

plot(ks_stat, type="b")
abline(v=90)
rsp(ks_stat)
rsp(colSums(apply(sc_area_warp, 1, quantile, probs=rg)))
rsp(rowMeans(sc_area_warp))
abline(v=115)

plot(ks_stat)
plot(colSums(apply(sc_area_warp, 1, quantile, probs=rg)))


sd.tmp = 0
for(i in 1:nrow(triangles)){
  sd.tmp = sd.tmp + bw.SJ(sc_area_warp[, i]) / nrow(triangles)
}
sd.tmp


(sum(sc_area_warp[ids_to_set[j2], i] > sc_area_warp[-ids_to_drop, i]) + 0.5) / (length(sc_area_warp[-ids_to_drop, i] + 1))

area_warp_mat_smooth_loglik = sc_area_warp #apply(sc_area_warp, 2, function(x){ qnorm(p = (rank(x)-0.5)/150, 0, 1) })
for(i in 1:nrow(triangles)){
  for(j in 1:25){
    ids_to_set = ((j-1)*4+1):(j*4)
    ids_to_set = c(ids_to_set, 100+((j-1)*2+1):(j*2))
    ids_to_drop = c(((j-1)*2+1):(j*2), 101:150)
    
    sd.tmp = bw.SJ(sc_area_warp[-ids_to_drop, i])
    
    for(j2 in 1:length(ids_to_set)){
      area_warp_mat_smooth_loglik[ids_to_set[j2], i] = sum( dnorm(sc_area_warp[ids_to_set[j2], i], sc_area_warp[-ids_to_drop, i], sd.tmp, log = TRUE) )
    }
  }
  print(i)
}

plot(rowMeans(area_warp_mat_smooth_loglik), type="b")
rsp(rowMeans(area_warp_mat_smooth_loglik))
rsp(rowMeans(area_warp_mat_smooth))


prova = area_warp_mat_smooth
for(i in 1:100){
  if(i <= 50){
    ids_to_set = c(i, i+100)
  }
  else{
    ids_to_set = i
  }
  ids_to_drop = c(i, 101:150)
  
  for(j2 in ids_to_set){
    tmp_rk = rep(0, ncol(area_warp_mat_smooth))
    for(j1 in (1:150)[-ids_to_drop]){
      tmp_rk = tmp_rk + (area_warp_mat_smooth[j2,] > area_warp_mat_smooth[j1,])
    }
    prova[j2,] = tmp_rk
  }
  print(i)
}

boxplot(t(prova))
rsp(apply(sc_area_warp, 1, median))

# rl1 = read.table("/Users/marcorudelli/Desktop/tesi magistrale/utch_a_1772114_sm2376/Supplementary Code/Code/cylinder_simulated_SIM1/rl1_cp.csv")[,1]
# rl2 = read.table("/Users/marcorudelli/Desktop/tesi magistrale/utch_a_1772114_sm2376/Supplementary Code/Code/cylinder_simulated_SIM1/rl2_cp.csv")[,1]
# 
# sim_res = read.table("/Users/marcorudelli/Desktop/tesi magistrale/utch_a_1772114_sm2376/Supplementary Code/Code/cylinder_simulated_SIM1/simulation_res.csv", header=FALSE, sep=",")
# 
# boxplot(cbind(sim_res, rl1, rl2))
# 
# tmp = density(rl1)
# plot(tmp)
# lines(density(sim_res[,1], bw = tmp$bw), col="red")
# 
# tmp = density(sim_res[,4])
# plot(tmp)
# lines(density(sim_res[,2], bw = tmp$bw), col="red")
# 
# apply(sim_res, 2, summary)
# apply(sim_res, 2, sd)


sim_res = read.table("/Volumes/EXTERNAL_USB/sim_2.1/simulation2.1_res.csv", header=FALSE, sep=",")
boxplot(sim_res[,c(2,4,5,6)])

apply(sim_res[,c(2,4,5,6)], 2, summary)
apply(sim_res[,c(2,4,5,6)], 2, mean)
apply(sim_res[,c(2,4,5,6)], 2, sd)

sim_res = read.table("/Volumes/EXTERNAL_USB/sim_2.2/simulation2.2_res.csv", header=FALSE, sep=",")
boxplot(sim_res[,c(2,4,5,6)])

apply(sim_res[,c(2,4,5,6)], 2, summary)
apply(sim_res[,c(2,4,5,6)], 2, mean)
apply(sim_res[,c(2,4,5,6)], 2, sd)

sim_res = read.table("/Volumes/EXTERNAL_USB/sim_2.3/simulation2.3_res.csv", header=FALSE, sep=",")
boxplot(sim_res[,c(2,4,5,6)])

apply(sim_res[,c(2,4,5,6)], 2, summary)
apply(sim_res[,c(2,4,5,6)], 2, mean)
apply(sim_res[,c(2,4,5,6)], 2, sd)

sim_res = read.table("/Volumes/EXTERNAL_USB/sim_2.9/simulation2.9_res.csv", header=FALSE, sep=",")
boxplot(sim_res[,c(2,4,5,6)])

apply(sim_res[,c(2,4,5,6)], 2, summary)
apply(sim_res[,c(2,4,5,6)], 2, mean)
apply(sim_res[,c(2,4,5,6)], 2, sd)



sim_res = read.table("/Volumes/EXTERNAL_USB/parts_simulated_SIM2.5/simulation2.5_res.csv", header=FALSE, sep=",")

boxplot(cbind(sim_res))
apply(sim_res, 2, summary)
apply(sim_res, 2, sd)


sim_res = read.table("/Volumes/EXTERNAL_USB/parts_simulated_SIM2.6/simulation2.6_res.csv", header=FALSE, sep=",")

boxplot(cbind(sim_res))
apply(sim_res, 2, summary)
apply(sim_res, 2, sd)


tmp = density(sim_res[,4])
plot(tmp, xlim=c(0,50))
lines(density(sim_res[,6], bw = tmp$bw), col="red")


sim_res = read.table("/Volumes/EXTERNAL_USB/sim_2.9/simulation2.9_res.csv", header=FALSE, sep=",")

boxplot(cbind(sim_res))
apply(sim_res, 2, summary)
apply(sim_res, 2, sd)



# rsp(apply(sc_err_forward_dist_smooth, 1, quantile, probs=0.9999))
# rsp(apply(sc_err_forward_dist_smooth, 1, kurtosis))

rsp(area_warp_mat_smooth[,19000])
rsp(apply(area_warp_mat_smooth, 1, mean))
rsp(area_warp_mat[,19000])
rsp(apply(area_warp_mat, 1, mean))
rsp(apply(sc_err_area_warp_mat_smooth, 1, mean))
plot(apply(sc_err_area_warp_mat_smooth, 1, mean), apply(area_warp_mat, 1, mean))


www = 1/apply(forward_dist[1:100,], 2, sd)
www = www / max(www)
www2 = 1/apply(area_warp_mat[1:100,], 2, sd)
www2 = www2 / max(www2)
boxplot(www, www2)

mmm = apply(area_warp_mat_smooth[1:100,], 2, mean)
boxplot(t(back_dist_norm))
#sss = sample(1:ncol(back_dist_new_2_mat_smooth), ncol(back_dist_new_2_mat_smooth)/4)
kkk = 110

plot(apply(back_dist_norm, 1, kurtosis)) 
colors = back_dist_new_2_mat_smooth[kkk, ] #forw_dist_sc[kkk,] #area_warp_mat_norm2[kkk,] #- area_warp_mat_norm2[kkk,] # 
#ccc[sss]
colors = colorRamp(c("green", "red"))((colors-min(colors))/(max(colors)-min(colors)))
colors <- rgb(colors/255)
x = verts[,1]#tr_centroids[sss,1]
y = verts[,2]#tr_centroids[sss,2]
z = verts[,3]#tr_centroids[sss,3]
i = triangles[,1]-1
j = triangles[,2]-1
k = triangles[,3]-1
c = colors
plot_ly(x = ~x,
        y = ~y,
        z =~ z, 
        i = ~i,
        j = ~j,
        k =~ k, 
        type = 'mesh3d',
        facecolor = c)
kkk = kkk + 1



lbs_res = read.csv("/Users/marcorudelli/Desktop/cartella senza nome/cylinder_simulated_SIM2.10/cylinder_simulated_res_LBS_SIM2.10_150.csv", header = FALSE)
dim(lbs_res)
kkk=1
plot(lbs_res[1:150,kkk], type="b")
rsp(lbs_res[1:150,kkk], alpha = 0.2)
kkk=kkk+1





plot(rchisq(1000, 3)/10000)
abline(h=0, col="red")
abline(h=0.0005, col="red")


boxplot(t(area_warp_mat))
boxplot(t(area_warp_mat_smooth))
boxplot(t(area_warp_mat_sc))
boxplot(t(area_warp_mat_sc2))


plot(apply(area_warp_mat, 1, mean))
plot(apply(area_warp_mat_smooth, 1, mean))
plot(apply(area_warp_mat_norm, 1, sd))
#plot(apply(area_warp_mat_sc2, 1, sd))
t.test(apply(area_warp_mat, 1, sd)[1:100], apply(area_warp_mat, 1, sd)[101:150])

plot(apply(area_warp_mat_norm2, 1, mean), apply(area_warp_mat_norm, 1, mean))



boxplot(t(forward_dist))
plot(apply(forward_dist, 1, mean))
boxplot(t(forward_dist_smooth))
plot(apply(forward_dist_smooth, 1, mean))
boxplot(t(forw_dist_sc))
plot(apply(forw_dist_norm, 1, mean))
t.test(apply(forw_dist_sc, 1, kurtosis)[1:100], apply(forw_dist_sc, 1, kurtosis)[101:150])
t.test(apply(forw_dist_norm, 1, quantile, probs = 0.8)[1:100], apply(forw_dist_norm, 1, quantile, probs = 0.8)[101:150])
#boxplot(t(area_warp_mat_sc2)) 



#NEI CILINDRI, QUANDO OFF CONTROL E' VARIAZIONE DI DELTA:

# PERCHE' LO SMOOTHING PEGGIORA LA PERFORMANCE (RISPETTO A NON SMOOTHING)
  # potenzialmente complesso, sto trasformando la funzione iniziale sperando di cancellare rumore e tenere il segnale
    # ciò che rimane dopo lo smoothing suppongo sia il segnale specifico di quell'uovo...
  # nella metrica raw la media non cambia, la varianza ha poco valore perchè influenzata dalla media funzionale (ma potenzialmente fare smoothing abbassa i picchi... meno differenza)
  # NELLA METRICA SCALATA/NORMALIZZATA DEVO RICORDARE CHE:
    # fare smoothing influenza il calcolo della media finale SOLO attraverso le s.d. calcolate per scalare (penso alla metrica scalata)
      # rispetto alla media RAW sto solo cambiando i pesi della media
        # e (in questo scenario) più faccio smoothing più aumento il peso di zone che contengono solo noise
    # in generale sto sempre dando lo stesso peso ad ogni triangolo (dopo aver scalato)
      # ma la varianza della media dei triangoli nelle facce orizzontali è minore pre smoothing rispetto a una media delle medie
        # perchè: scalo a var unitaria tutti i tirangoli e poi faccio la media, rispetto a, faccio tante medie, le scalo a varianza unitaria, e poi faccio la media
          # la riscalatura è più grande nel secondo caso, aumento l'errore (è come prima avere tot triangoli "incorrelati" e poi averli correlati)
    # QUINDI (sempre per la metrica scalata o normalizzata):
      # sorveglianza in media: discorso appena fatto (in questo caso è controproducente, ma in generale è più robusto!)
        # dò più peso a parti che normalmente variano poco... cogliendo la possibilità di segnalare ovunque
          # più robusto a spese di un po' do potenza nelle parti che variano tanto
      # invece, per la sorveglianza della varianza: forse facendo smoothing liscio i picchi, che sono maggiormente presenti nei dati IC, rendendo meno evidente il calo della varianza quando vado OC

# PERCHE' NORMALIZZARE PEGGIORA LA PERFORMANCE (RISPETTO A SCALARE)
  # se il segnale che sposta la media sta nei picchi, normalizzando li sostituisco con valori meno "estremi" e quindi rendo la differenza in media meno evidente
  # INOLTRE, sto rendendo la distribuzione simmetrica, dando più peso alla coda inferiore, non è detto che questo sia di interesse per la rilevazione in questo caso... ma comunque rendo il metodo più robusto!
  # inoltre miglioro la possibilità di cogliere variazioni in skewness e curtosi (ha senso se voglio monitorare i primi 4 momenti)

# PERCHE' SCALARE (o normalizzare) PEGGIORA LA PERFORMANCE (RISPETTO A METRICA RAW)
  # il segnale sta in una zona a varianza elevata, sto dando più potenza alla zona priva di segnale, amplificando il rumore... (ma concettualmente ha senso, quindi lo tengo)
    # in pratica rendo il metodo più robusto a segnalare cambiamenti, ovunque si verifichino nella mesh
  # inoltre tolgo la media, quindi il monitoraggio della varianza marginale non risente più delle variazioni della media funzionale (che ha senso, perchè inquinavano la rilevazione...)



# SE FACCIO SMOOTHING ABBATTO LA VARIABILITA' DEL RUMORE AD ALTA FREQUENZA (isotropico) MA NON ABBATTO LA VARIABILITA' DEL RUMORE A BASSA FREQUENZA!!!
  # NEL CASO DEI CINILDRI, ABBATTO LA VARIABILITA' DEL RUMORE NELLE DUE FACCE CIRCOLARI MA NON QUELLA DELLA PARETE VERTICALE
    # POI CON LO SCALE RIPORTO IL RUMORE AD AVERE LA STESSA VARIANZA OVUNQUE
      # RAGIONARE SE È QUELLO CHE VOGLIO!!!
# LA QUESTIONE: che senso ha fare smoothing per abbattere il rumore, potenzialmente facendo più smoothing per abbattere di più il rumore
  # se poi lo riscalo in ogni caso a varianza unitaria, addirittura peggiorando la media rispetto a se non avessi fatto smoothing???
    # l'idea è che se c'è segnale, fare smoothing e cancellare il meglio possibile la varianza, ti porta a rivelare il segnale
      # ovviamente non posso sapere dove sta il segnale... allora sono costretto a rischiare di amplificare rumore per rendere il metodo robusto...
        # POTREI PENSARCI MEGLIO MA PER ORA LA CONCLUSIONE E' QUESTA...


# la differenza tra scalare e normalizzare è che quando scalo ottengo più potenza se le distribuzioni sono già di loro abbastanza normali
  # quindi va bene codificare i dati di fase 2 attraverso il numero di deviazioni standard
    # e non perdo gli outliers molto estremi
# se invece normalizzo è perchè essere a Y deviazioni standard dalla media non ha significato se le distribuzioni sono diverse 



# piccolo scenario per capire perchè fare smoothing potrebbe mascherare differenze in varianza
# tmpa = rep(NA, 100)
# tmpa2 = rep(NA, 100)
# tmp = rep(NA, 100)
# tmp2 = rep(NA, 100)
# for(i in 1:50){
#   y = rnorm(200,0,1)
#   tmpa[i] = mean(y)
#   tmpa2[i] = sd(y)
#   y[1:100] = mean(y[1:100])
#   y[101:200] = mean(y[101:200])
#   tmp[i] = mean(y)
#   tmp2[i] = sd(y)
# }
# for(i in 51:100){
#   y = c(rnorm(5,5,1), rnorm(190,0,1), rnorm(5,5,1))
#   tmpa[i] = mean(y)
#   tmpa2[i] = sd(y)
#   y[1:100] = mean(y[1:100])
#   y[101:200] = mean(y[101:200])
#   tmp[i] = mean(y)
#   tmp2[i] = sd(y)
# }
# plot(tmpa)
# plot(tmp)
# plot(tmpa2)
# plot(tmp2)


# dim(area_warp_mat_norm)

boxplot(t(area_warp_mat_norm))

#library(glmnet)
# tmp = cv.glmnet(area_warp_mat_norm, c(rep(0,100), rep(1,50)), family = "binomial", lambda.min.ratio = 1e-14, nfolds = 50)
# plot(tmp)
# 
# colors = as.numeric(coef(tmp, s="lambda.min"))!=0
# colors = colorRamp(c("green", "red"))((colors-min(colors))/(max(colors)-min(colors)))
# colors <- rgb(colors/255)
# x = tr_centroids[,1]
# y = tr_centroids[,2]
# z = tr_centroids[,3]
# c = colors
# plot_ly(x = ~x,
#         y = ~y,
#         z =~ z, type = 'scatter3d', mode = 'markers',
#         marker = list(size = 4, color=c))

# which(colors>0)
# 
# # install.packages("genlasso")
# # library(genlasso)
# tr_adj_list = list()
# 
# for(i in 1:length(smoothing_list)){
#   for(j in which(smoothing_list[[i]][2,] > 1-1e-10)){
#     if(smoothing_list[[i]][1,j] > i){
#       tr_adj_list[[length(tr_adj_list) + 1]] = c(i, smoothing_list[[i]][1,j])
#       names(tr_adj_list[[length(tr_adj_list)]]) = c("t1", "t2")
#     }
#   }
# }
# length(tr_adj_list)
# tr_adj_mat = bind_rows(tr_adj_list)
# 
# #library(Matrix)
# gr = graph_from_edgelist(as.matrix(tr_adj_mat), directed = FALSE)
# 
# A = rbind(cbind(seq(1,nrow(tr_adj_mat)),as.matrix(tr_adj_mat)[,1],-1),cbind(seq(1,nrow(tr_adj_mat)),as.matrix(tr_adj_mat)[,2],1))
# D = sparseMatrix(i=A[,1], j=A[,2], x=A[,3], dims=c(ecount(gr),vcount(gr)))
# 
# prova = fusedlasso(area_warp_mat_norm[1,], D = D, approx = TRUE, verbose = TRUE, )
# dim(D)



#library(data.table)

rl = rep(NA, ncol=100)
for(kkk in 1:100){
  prova = fread(strcat(strcat("/Volumes/EXTERNAL_USB/sim_2.1/area_warp_smooth_mats/area_warp_mat_smooth_", itostr(kkk)), ".csv"))
  #library(mvtnorm)
  # prova_IC = matrix(NA, ncol=4000, nrow=150)
  # for(i in 1:(4000/4)){
  #   prova_IC[, ((i-1)*4+1):(i*4)] = rmvnorm(150, c(1,2,0,4), matrix(c(1.0,0.4,0.2,0.0,
  #                                                                     0.4,2.0,0.2,0.1,
  #                                                                     0.2,0.2,1.0,0.3,
  #                                                                     0.0,0.1,0.3,0.8), byrow=TRUE, ncol=4))
  # }
  
  minwin = 5
  maxwin = 10
  win = 1:50
  win = pmax(minwin, win)
  win = pmin(win, maxwin)
  
  sampleNum = 101:150
  
  mu = (1-(1-lambda)^win)/lambda*(sampleNum + 1)/2
  
  lambda = 0.01
  lam = 1-lambda
  temp = lam^win
  
  sd = (1-temp^2)*(sampleNum+1)*(sampleNum-1)/((2*lambda-lambda^2)*12)
  sd = sd - ( (lam-temp)/lambda^2-(lam^2-temp^2)/(2*lambda^2-lambda^3) )*(sampleNum+1)/6;
  sd = sqrt(sd);
  
  working_weight = (maxwin-1):0
  working_weight = (1-lambda)^working_weight
  working_wm = matrix(working_weight, ncol=ncol(prova), nrow=length(working_weight), byrow=FALSE)
  
  tmp2 = matrix(NA, ncol=ncol(prova), nrow=50)
  t.stat = matrix(NA, ncol=50, nrow=5) #rep(NA, 50)
  lo_lim = matrix(NA, ncol=50, nrow=5) #rep(NA, 50)
  up_lim = matrix(NA, ncol=50, nrow=5) #rep(NA, 50)
  for(i in 101:150){
    win_ids = (1+i-win[i-100]):i
   
    tmp = apply(prova[1:i,, drop=FALSE], 2, rank)
    
    tmp2[i-100,] = (colSums( tmp[win_ids,, drop=FALSE]*working_wm[(maxwin-length(win_ids)+1):maxwin,] ) - mu[i-100] ) / sd[i-100]
    t.stat[,i-100] = quantile(tmp2[i-100,], probs=c(0.01, 0.1, 0.5, 0.9, 0.99)) # mean(tmp2[i-100,])
   
    count = 1
    t_perm = matrix(NA, ncol=10000, nrow=5)
    while(count <= 10000){
      sss = sample(1:i)
     
      perm_ok = TRUE
      if(i > 101){
        working_tmp = tmp[sss,]
        to_check = (i-1):max(i-win[i-100]+1, 101)
        
        min_id = 1+min(to_check)-win[min(to_check)-100]
        
        for(j in to_check){
          working_win_ids = (1+j-win[j-100]):j
          
          for(k in max(working_win_ids):min_id){
            sel = working_tmp[j+1,] < working_tmp[k,]
            sel2 = (working_tmp[j+1,] == working_tmp[k,]) / 2
            working_tmp[k,] = working_tmp[k,] - sel - sel2
          }
          
          t.stat.temp = quantile( (colSums( working_tmp[working_win_ids,, drop=FALSE]*working_wm[(maxwin-length(working_win_ids)+1):maxwin,] ) - mu[j-100] ) / sd[j-100], probs=c(0.01, 0.1, 0.5, 0.9, 0.99))
                        # mean((colSums( working_tmp[working_win_ids,, drop=FALSE]*working_wm[(maxwin-length(working_win_ids)+1):maxwin,] ) - mu[j-100] ) / sd[j-100])
          
          if(any(t.stat.temp > up_lim[,j-100]) || any(t.stat.temp < lo_lim[,j-100])){
            perm_ok = FALSE
            break
          }
        }
      }
      
      if(perm_ok){
        t_perm[,count] = quantile( (colSums( tmp[sss[win_ids],, drop=FALSE]*working_wm[(maxwin-length(win_ids)+1):maxwin,] ) - mu[i-100] ) / sd[i-100], probs=c(0.01, 0.1, 0.5, 0.9, 0.99))
                        # mean( (colSums( tmp[sss[win_ids],, drop=FALSE]*working_wm[(maxwin-length(win_ids)+1):maxwin,] ) - mu[i-100] ) / sd[i-100] )
        count = count + 1
        if(count%%100==0) print(count)
      }
    }
   
    num_ok = 10000
    alfa = 10000
    t_perm_sort = t(apply(t_perm, 1, sort))
    while(num_ok > 10000*(1-0.005)){
      alfa = alfa - 1
      lims_up = t_perm_sort[,alfa]
      lims_down = t_perm_sort[,10000-alfa]
      
      num_ok = sum(apply(t_perm, 2, function(x){sum(x > lims_up | x < lims_down)}) == 0)
    }
    
    lo_lim[,i-100] = lims_down # quantile(t_perm, 0.1/2)
    up_lim[,i-100] = lims_up # quantile(t_perm, 1-0.1/2)
    if( any(t.stat[,i-100] > up_lim[,i-100]) || any(t.stat[,i-100] < lo_lim[,i-100]) ){
      break
    }
    
    print(i)
  }
  
  rl[kkk] = i-100
  print(kkk)
  
  par(mfrow=c(2,3))
  boxplot(t(tmp2[1:(i-100),]))
  for(jjj in 1:5){
    plot(t.stat[jjj,1:(i-100)], type="b", ylim=range(c(lo_lim[jjj,1:(i-100)], up_lim[jjj,1:(i-100)], t.stat[jjj,1:(i-100)]), na.rm = TRUE))
    lines(lo_lim[jjj,1:(i-100)], type="b", col="red")
    lines(up_lim[jjj,1:(i-100)], type="b", col="red")
  }
  par(mfrow=c(1,1))
}

# mean(rl)
# sd(rl)
# curve(ecdf(rl)(x), from=0, to=50)
# curve(pgeom(x-1, 0.1), add=TRUE, col="red")

boxplot(rl)
summary(rl)



# library(BBmisc)
# library(pracma)

for(jjj in 1:300){
  area_warp_mat_smooth = matrix(rnorm(10*1000), ncol=10, nrow=1000)
  back_dist_new_2_mat_smooth = matrix(rnorm(10*1000), ncol=10, nrow=1000)
  forward_dist_smooth = matrix(rnorm(10*1000), ncol=10, nrow=1000)
  
  
  write.table(area_warp_mat_smooth, file = strcat(strcat("/Volumes/EXTERNAL_USB/sim_IC/area_warp_smooth_mats/area_warp_mat_smooth_", itostr(jjj)), ".csv"), col.names = FALSE, row.names = FALSE, sep = ",")
  write.table(back_dist_new_2_mat_smooth, file = strcat(strcat("/Volumes/EXTERNAL_USB/sim_IC/backward_dist_mats/backward_dist_smooth_", itostr(jjj)), ".csv"), col.names = FALSE, row.names = FALSE, sep = ",")
  write.table(forward_dist_smooth, file = strcat(strcat("/Volumes/EXTERNAL_USB/sim_IC/forward_dist_mats/forward_dist_smooth_", itostr(jjj)), ".csv"), col.names = FALSE, row.names = FALSE, sep = ",")
  
  print(jjj)
}

#library(data.table)
rl = rep(NA, 7)
for(kkk in 1:7){
  # lo_lim = as.matrix(read.table(strcat(strcat("/Volumes/EXTERNAL_USB/sim_IC/res_surv_new/lo_lim_", itostr(kkk)), ".txt"), header = FALSE))
  # up_lim = as.matrix(read.table(strcat(strcat("/Volumes/EXTERNAL_USB/sim_IC/res_surv_new/up_lim_", itostr(kkk)), ".txt"), header = FALSE))
  t.stat = as.matrix(read.table(strcat(strcat("/Volumes/EXTERNAL_USB/sim_IC/res_surv_new/t_stat_", itostr(kkk)), ".txt"), header = FALSE))
  # tmp2_1 = fread(strcat(strcat("/Volumes/EXTERNAL_USB/sim_IC/res_surv_new/tmp2_1_", itostr(kkk)), ".txt"), header = FALSE, sep = " ")
  # tmp2_2 = fread(strcat(strcat("/Volumes/EXTERNAL_USB/sim_IC/res_surv_new/tmp2_2_", itostr(kkk)), ".txt"), header = FALSE, sep = " ")
  i = max(which(t.stat[1,]!=0))
  
  # par(mfrow=c(3,5))
  # boxplot(t(tmp2_1[1:i,]))
  # for(jjj in 1:4){
  #   plot(t.stat[jjj,1:i], type="b", ylim=range(c(lo_lim[jjj,1:i], up_lim[jjj,1:i], t.stat[jjj,1:i]), na.rm = TRUE))
  #   lines(lo_lim[jjj,1:i], type="b", col="red")
  #   lines(up_lim[jjj,1:i], type="b", col="red")
  # }
  # boxplot(t(tmp2_2[1:i,]))
  # for(jjj in 5:8){
  #   plot(t.stat[jjj,1:i], type="b", ylim=range(c(lo_lim[jjj,1:i], up_lim[jjj,1:i], t.stat[jjj,1:i]), na.rm = TRUE))
  #   lines(lo_lim[jjj,1:i], type="b", col="red")
  #   lines(up_lim[jjj,1:i], type="b", col="red")
  # }
  # boxplot(t(tmp2_3[1:i,]))
  # for(jjj in 9:12){
  #   plot(t.stat[jjj,1:i], type="b", ylim=range(c(lo_lim[jjj,1:i], up_lim[jjj,1:i], t.stat[jjj,1:i]), na.rm = TRUE))
  #   lines(lo_lim[jjj,1:i], type="b", col="red")
  #   lines(up_lim[jjj,1:i], type="b", col="red")
  # }
  # par(mfrow=c(1,1))
  rl[kkk] = i
}


mean(rl==1, na.rm = TRUE)*100
mean(rl, na.rm = TRUE)
sd(rl, na.rm = TRUE)
mean(rl==900, na.rm = TRUE)*100

mean(rgeom(300, 0.005)+1)
mean(rgeom(300, 0.005)+1 == 1)*100
mean(rgeom(300, 0.005)+1 >= 900)*100


curve(pgeom(x-1,0.005), col="red", from=0, to=899, n=1000)
for(i in 1:100) curve(ecdf(rgeom(300, 0.005)+1)(x), from=0, to=899, add=TRUE, lwd=0.5, col="green", n=1000)
curve(ecdf(rl)(x), from=0, to=899, add=TRUE, n=1000)
curve(pgeom(x-1,0.005), col="red", from=0, to=899, add=TRUE, n=1000)


rl2 = rep(NA, 100)
for(kkk in 1:100){
  lo_lim = as.matrix(read.table(strcat(strcat("/Volumes/EXTERNAL_USB/sim_2.2/res_surv_new/lo_lim_", itostr(kkk)), ".txt"), header = FALSE))
  up_lim = as.matrix(read.table(strcat(strcat("/Volumes/EXTERNAL_USB/sim_2.2/res_surv_new/up_lim_", itostr(kkk)), ".txt"), header = FALSE))
  t.stat = as.matrix(read.table(strcat(strcat("/Volumes/EXTERNAL_USB/sim_2.2/res_surv_new/t_stat_", itostr(kkk)), ".txt"), header = FALSE))
  # tmp2_1 = fread(strcat(strcat("/Volumes/EXTERNAL_USB/sim_2.1/res_surv_new/tmp2_1_", itostr(kkk)), ".txt"), header = FALSE, sep = " ")
  # tmp2_2 = fread(strcat(strcat("/Volumes/EXTERNAL_USB/sim_2.1/res_surv_new/tmp2_2_", itostr(kkk)), ".txt"), header = FALSE, sep = " ")
  i = max(which(t.stat[1,]!=0))
  
  if(i==50){
    i = i + (up_lim[50] >= t.stat[50] & lo_lim[50] <= t.stat[50])
  }
  # kkk=10
  # par(mfrow=c(2,5))
  # boxplot(t(tmp2_1[1:i,]))
  # for(jjj in 1:4){
  #   plot(t.stat[jjj,1:i], type="b", ylim=range(c(lo_lim[jjj,1:i], up_lim[jjj,1:i], t.stat[jjj,1:i]), na.rm = TRUE))
  #   lines(lo_lim[jjj,1:i], type="b", col="red")
  #   lines(up_lim[jjj,1:i], type="b", col="red")
  # }
  # boxplot(t(tmp2_2[1:i,]))
  # for(jjj in 5:8){
  #   plot(t.stat[jjj,1:i], type="b", ylim=range(c(lo_lim[jjj,1:i], up_lim[jjj,1:i], t.stat[jjj,1:i]), na.rm = TRUE))
  #   lines(lo_lim[jjj,1:i], type="b", col="red")
  #   lines(up_lim[jjj,1:i], type="b", col="red")
  # }
  # par(mfrow=c(1,1))
  rl2[kkk] = i
  # 
  # kkk=kkk+1
}
# rl22 = rl2
# table(rl2)
# mean(rl2)
# mean(rl22)
# plot(rl22, rl2)
rl2
mean(rl2>50)
boxplot(rl2)

mean(rl2[1:max(which(!is.na(rl2)))])
sd(rl2[1:max(which(!is.na(rl2)))])
plot(rl2, ylim=c(0,50))
(sum(rl2==50)*25+sum(rl2))/100
mean(sim_res[1:max(which(!is.na(rl2))),4])

cbind(sim_res[1:max(which(!is.na(rl2))),4], rl2[1:max(which(!is.na(rl2)))])
plot(sim_res[1:max(which(!is.na(rl2))),4], rl2[1:max(which(!is.na(rl2)))])
abline(0,1)

boxplot(sim_res[1:max(which(!is.na(rl2))),4], rl2[1:max(which(!is.na(rl2)))])

boxplot(sim_res[1:max(which(!is.na(rl2))),6], rl2[1:max(which(!is.na(rl2)))])


rl = rep(NA, 100)
for(kkk in 1:100){
  lim = as.matrix(read.table(strcat(strcat("/Volumes/EXTERNAL_USB/sim_2.11/res_surv_old_lbs/lim_", itostr(kkk)), ".txt"), header = FALSE))
  t.stat = as.matrix(read.table(strcat(strcat("/Volumes/EXTERNAL_USB/sim_2.11/res_surv_old_lbs/t_stat_", itostr(kkk)), ".txt"), header = FALSE))
  # tmp2 = fread(strcat(strcat("/Volumes/EXTERNAL_USB/sim_2.1/res_surv_old_mm/tmp2_", itostr(kkk)), ".txt"), header = FALSE, sep = " ")
  i = max(which(t.stat!=0))
  
  if(i==50){
    i = i + (lim[50] >= t.stat[50])
  }
  # par(mfrow=c(1,2))
  # matplot(tmp2[1:i,], type="b")
  # 
  # plot(t.stat[1:i], type="b", ylim=range(c(0, lim[1:i], t.stat[1:i]), na.rm = TRUE))
  # lines(lim[1:i], type="b", col="red")
  # par(mfrow=c(1,1))
  
  rl[kkk] = i
  
  kkk=kkk+1
}

rl
boxplot(rl)
mean(rl>50)
mean(rl[1:max(which(!is.na(rl)))])
sd(rl[1:max(which(!is.na(rl)))])
mean(sim_res[1:max(which(!is.na(rl))),4])

cbind(sim_res[1:max(which(!is.na(rl))),4], rl[1:max(which(!is.na(rl)))])
plot(sim_res[1:max(which(!is.na(rl))),4], rl[1:max(which(!is.na(rl)))])
abline(0,1)

boxplot(sim_res[1:max(which(!is.na(rl))),4], rl[1:max(which(!is.na(rl)))])

boxplot(sim_res[1:max(which(!is.na(rl))),6], rl[1:max(which(!is.na(rl)))])
plot(sim_res[1:max(which(!is.na(rl))),6], rl[1:max(which(!is.na(rl)))])

table(rl)
boxplot(rl, rl2, ylim=c(0,50))
abline(h=0)
mean(rl)
mean(rl2)
plot(rl, rl2, xlim=c(0,50), ylim=c(0,50))
abline(0,1)
cbind(rl, rl2)
plot(sim_res[,6], rl)
# aaa_mat = matrix(NA, ncol=10000, nrow=4)
# for(i in 1:10000){
#   tmp = t(apply(matrix(rep(1:110, 10), ncol=110, byrow=TRUE), 1, sample))
#   aaa = (rowSums(tmp[,106:110])-(55.5*5)) / 69.8
#   aaa_mat[1,i] = mean(aaa)
#   aaa_mat[2,i] = var(aaa)
#   aaa_mat[3,i] = skewness(aaa)
#   aaa_mat[4,i] = kurtosis(aaa)
# }
# cor(t(aaa_mat))
# 
# 
# 
# tmp = apply(matrix(rnorm(10*110), ncol=110), 1, rank)
# sss = (colSums(tmp[101:110,]) - (10*55.5)) / (96.2)
# 
# sss = matrix(NA, nrow=4, ncol=25000)
# for(i in 1:25000){
#   tmp2 = (colSums(tmp[sample(1:110)[101:110],]) - (10*55.5)) / (96.2)
#   sss[1,i] = mean(tmp2)
#   sss[2,i] = var(tmp2)
#   sss[3,i] = skewness(tmp2)
#   sss[4,i] = kurtosis(tmp2)
# }
# 
# num_ok = 25000
# alfa = 25000
# num_ok_pr = num_ok
# alfa_pr = alfa
# lim_up_pr = 0
# lim_down_pr = 0
# lims_up = 0
# lims_down = 0
# t_perm_sort = t(apply(sss, 1, sort))
# while(num_ok > 25000*(1-0.005)){
#   alfa_pr = alfa
#   alfa = alfa - 1
#   lim_up_pr = lims_up
#   lim_down_pr = lims_down
#   lims_up = t_perm_sort[,alfa]
#   lims_down = t_perm_sort[,25000-alfa]
#   
#   num_ok_pr = num_ok
#   num_ok = sum(apply(sss, 2, function(x){sum(x > lims_up | x < lims_down)}) == 0)
# }
# num_ok_pr/25000
# alfa_pr/25000
# 
# 
# 0.99928^4
# 
# plot(sss[1,], sss[2,])
# plot(sss[2,], sss[3,])
# plot(sss[3,], sss[4,])

# norm_quantile_mat = matrix(0, nrow=50, ncol=150)
# for(j in 101:150) {
#   tmp_seq = seq(from=0, to=1, length.out = j + 2)
#   tmp_seq = tmp_seq[-c(1,length(tmp_seq))]
#   norm_quantile_mat[j-100, 1:j] = qnorm(tmp_seq, 0, 1)
# }
# write.table(norm_quantile_mat, file="/Volumes/EXTERNAL_USB/norm_quantile_mat.csv", col.names = FALSE, row.names = FALSE, sep = ",")
# 
# 
# 
# norm_quantile_mat = matrix(0, nrow=900, ncol=1000)
# for(j in 101:1000) {
#   tmp_seq = seq(from=0, to=1, length.out = j + 2)
#   tmp_seq = tmp_seq[-c(1,length(tmp_seq))]
#   norm_quantile_mat[j-100, 1:j] = qnorm(tmp_seq, 0, 1)
# }
# write.table(norm_quantile_mat, file="/Volumes/EXTERNAL_USB/norm_quantile_mat_IC.csv", col.names = FALSE, row.names = FALSE, sep = ",")




