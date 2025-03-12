NUM_SIMS = 10
cens_lim = 50

rl_new = rep(NA, NUM_SIMS)
for(kkk in 1:NUM_SIMS){
  lo_lim = as.matrix(read.table(strcat(strcat("res_surv_new/lo_lim_", itostr(kkk)), ".txt"), header = FALSE))
  up_lim = as.matrix(read.table(strcat(strcat("res_surv_new/up_lim_", itostr(kkk)), ".txt"), header = FALSE))
  t.stat = as.matrix(read.table(strcat(strcat("res_surv_new/t_stat_", itostr(kkk)), ".txt"), header = FALSE))
  # tmp2_1 = fread(strcat(strcat("/Volumes/EXTERNAL_USB/sim_2.1/res_surv_new/tmp2_1_", itostr(kkk)), ".txt"), header = FALSE, sep = " ")
  # tmp2_2 = fread(strcat(strcat("/Volumes/EXTERNAL_USB/sim_2.1/res_surv_new/tmp2_2_", itostr(kkk)), ".txt"), header = FALSE, sep = " ")
  i = max(which(t.stat[1,]!=0))
  
  if(i==cens_lim){
    i = i + (up_lim[cens_lim] >= t.stat[cens_lim] & lo_lim[cens_lim] <= t.stat[cens_lim])
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
  rl_new[kkk] = i
  # 
  # kkk=kkk+1
}

mean(rl_new>cens_lim)
boxplot(rl_new)

mean(rl_new)
sd(rl_new)

rl_old = rep(NA, NUM_SIMS)
for(kkk in 1:NUM_SIMS){
  lim = as.matrix(read.table(strcat(strcat("res_surv_old/lim_", itostr(kkk)), ".txt"), header = FALSE))
  t.stat = as.matrix(read.table(strcat(strcat("res_surv_old/t_stat_", itostr(kkk)), ".txt"), header = FALSE))
  # tmp2 = fread(strcat(strcat("/Volumes/EXTERNAL_USB/sim_2.1/res_surv_old_mm/tmp2_", itostr(kkk)), ".txt"), header = FALSE, sep = " ")
  i = max(which(t.stat!=0))
  
  if(i==cens_lim){
    i = i + (lim[cens_lim] >= t.stat[cens_lim])
  }
  # par(mfrow=c(1,2))
  # matplot(tmp2[1:i,], type="b")
  # 
  # plot(t.stat[1:i], type="b", ylim=range(c(0, lim[1:i], t.stat[1:i]), na.rm = TRUE))
  # lines(lim[1:i], type="b", col="red")
  # par(mfrow=c(1,1))
  
  rl_old[kkk] = i
  
  kkk=kkk+1
}

mean(rl_old>cens_lim)
boxplot(rl_old)

mean(rl_old)
sd(rl_old)


boxplot(rl_new, rl_old, ylim=c(0,cens_lim))
