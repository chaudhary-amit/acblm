

#' Prepare the finction for test run
#' intialize the model parameters (helper function)
#' @export
ModelInitializer <- function(){
  model = m2.mixt.new(nk=2,nf=4,nb=1)

  #model = m2.mixt.new(nk=2,nf=6,nb=1)

  model$A1[1,,] = seq(0.1,1,l=model$nf) %o% seq(0.1,1,l=model$nk)
  model$A2 = model$A1

  model$S1[] = 0.1
  model$S2[] = 0.1


  model$pk1[,,]=0.5
  model$pk0[,,]=0.5
  model$NNm[1,1,,] = 1000
  return(model)
}


#' Simulate the data consiting of three sided heterogeniety
#' @export
Simulate.data.threeSided <- function(model_test){

  test_data = m2.mixt.simulate.sim.clust(model_test,fsize=50,msize=100)
  # get the transformed 2-d model
  #model_ret <- m2.mixt.transform.model(model_test,model)
  # get the transformed data
  #model_trans_data <- m2.mixt.transform.data(model_test,model,test_data)
  return(test_data)

}

#' Set the solver controls
#' @export
set.solver.controls <- function(model_test=model_test, n_startValues=1, stayers_sample=0.1){
  model <- ModelInitializer()
  model_ret <- m2.mixt.transform.model(model_test,model)

  ctrl  = em.control( nplot=10000,         # how often to plot (either wages, or model versus model0)
                      ncat =2000,          # how often to show step increases (large to keep output small)
                      model0 = model_ret ,      # a model to compare estimates too (this is for testing)
                      fixb=TRUE,           # impose fixed interactions over time
                      tol=1e-7,            # tolerance to stop EM
                      est_rep=1,           # number of starting values to use (we usually use 50)
                      est_nbest= n_startValues,         # how many best liklihood to choose the best connected model from (we usually use 10)
                      sdata_subsample=stayers_sample, # whether to sub-sample the stayers in estimation
                      maxiter=1000)

  return(ctrl)
}

#' Clustering (Step 1 of the estimation: Seep paper for details)
#' @export
threeSided.Clustering <- function( test_data, # my data
                                   moments="ecdf", # moments type
                                   Nw=20, #support
                                   y_var="y1", # variable
                                   three_way="F"){

  flog.info("here i am ",name='logger.a')
  #clone data
  sdata_clone = test_data$sdata
  jdata_clone = test_data$jdata

  #paste manager to firm
  sdata_clone$f1=as.character(test_data$sdata$g1)
  sdata_clone$f2=as.character(test_data$sdata$g2)
  jdata_clone$f1=as.character(test_data$jdata$g1)
  jdata_clone$f2=as.character(test_data$jdata$g2)

  sdata_clone$j1true=as.character(test_data$sdata$g1true)

  ad_new_em_clone <- list(sdata=sdata_clone,jdata=jdata_clone)
  flog.info("generating measures for firms")
  ms_new_em = grouping.getMeasures.em(test_data,"ecdf",Nw=20,y_var = "y1",three_way = F)
  grps_new_em  = grouping.classify.once(ms_new_em,k = 2,nstart = 1000,iter.max = 200,step=250)
  ad_employee_em   = grouping.append(test_data,grps_new_em$best_cluster,drop=T)

  #use cloned data to add manager classes
  flog.info("generating measures for managers")
  ms_new_em_clone    = grouping.getMeasures.em(ad_new_em_clone,"ecdf",Nw=20,y_var = "y1",three_way = F)
  grps_new_em_clone  = grouping.classify.once(ms_new_em_clone,k = 2,nstart = 1000,iter.max = 200,step=250)
  ad_employee_em_clone    = grouping.append(ad_new_em_clone,grps_new_em_clone$best_cluster,drop=T)

  ad_employee_em$sdata$m1 = ad_employee_em_clone$sdata$j1[match(unlist(ad_employee_em$sdata$g1),ad_employee_em_clone$sdata$g1)]
  ad_employee_em$sdata$m2 = ad_employee_em_clone$sdata$j2[match(unlist(ad_employee_em$sdata$g2),ad_employee_em_clone$sdata$g2)]
  ad_employee_em$jdata$m1 = ad_employee_em_clone$jdata$j1[match(unlist(ad_employee_em$jdata$g1),ad_employee_em_clone$jdata$g1)]
  ad_employee_em$jdata$m2 = ad_employee_em_clone$jdata$j2[match(unlist(ad_employee_em$jdata$g2),ad_employee_em_clone$jdata$g2)]
  return(ad_employee_em)
}

#' Mixture model estimation of three sided model
#' @export
estimation.threeSided.model <- function(model_test,ad_employee_em,ctrl){
  # model is 2d initiliazer
  # model_test is 3d
  # model_ret is 2d transformed
  model <- ModelInitializer()
  model_ret <- m2.mixt.transform.model(model_test,model)
  model_trans_data_new <- m2.mixt.transform.data(model_test,model,ad_employee_em)
  my_model_test_cluster = m2.mixt.estimate.all(sim=model_trans_data_new,         # the data
                                               nk=model_ret$nk,    # number of worker types (we use here the same as in simulation)
                                               ctrl=ctrl,nbb=2)      # parameters to control the EM

  return(my_model_test_cluster)
}

#' Proportion plot
#' @export
threeSided.proportion.plot <- function(my_model_test_cluster, m){
  if(m==1) {
    m2.mixt.pplot(my_model_test_cluster$model$pk0[1,1:2,]) +
      xlab("police station class") + ylab("Proportions") +
      ggtitle("Proportions (estimated) of worker type \n (Manager class = 1)") + scale_fill_discrete(name="Worker\nType")+
      theme(plot.title = element_text(hjust = 0.5))
  } else if (m==2) {
    m2.mixt.pplot(my_model_test_cluster$model$pk0[1,3:4,]) +
      xlab("police station class") + ylab("Proportions") +
      ggtitle("Proportions (estimated) of worker type \n (Manager class = 2)") + scale_fill_discrete(name="Worker\nType")+
      theme(plot.title = element_text(hjust = 0.5))
  }

}

#' Estimated means plot
#' @export
threeSided.means.plot <- function(my_model_test_cluster, m){
  if(m==1){
    m2.mixt.wplot(my_model_test_cluster$model$A1[1,1:2,]) +
      geom_point() +
      ggtitle("Mean productivity (Manager class = 1)") +
      xlab("police station class") +
      theme_bw() +
      scale_y_continuous("log productivity")
  } else if(m==2){
    m2.mixt.wplot(my_model_test_cluster$model$A1[1,3:4,]) +
      geom_point() +
      ggtitle("Mean productivity (Manager class = 2)") +
      xlab("police station class") +
      theme_bw() +
      scale_y_continuous("log productivity")
  }
}

############################################################


######### Model dimension conversion ##################

##########################################################

# function takes 3d model object
# transforms it to 2d model object
# returns 2d model object
#' Model transformation for solver (3d array object to 2d array objects)
#' used in multicore effitiency
#' @export
m2.mixt.transform.model <- function(model_3d, model_2d_init, data=NA) {

  nb_3d  = model_3d$nb
  nf_3d  = model_3d$nf
  nk_3d  = model_3d$nk
  A1_3d  = model_3d$A1
  S1_3d  = model_3d$S1
  A2_3d  = model_3d$A1
  S2_3d  = model_3d$S1

  NNm_3d = model_3d$NNm
  pk1_3d = model_3d$pk1
  pk0_3d = model_3d$pk0


  nb_2d  = model_2d_init$nb
  nf_2d  = model_2d_init$nf
  nk_2d  = model_2d_init$nk

  A1_2d  = model_2d_init$A1
  S1_2d  = model_2d_init$S1
  A2_2d  = model_2d_init$A2
  S2_2d  = model_2d_init$S2
  NNm_2d = model_2d_init$NNm
  pk1_2d = model_2d_init$pk1
  pk0_2d = model_2d_init$pk0

  mapping = array(0,c(nb_3d*nf_3d,3))
  # create the mappings
  ci=1
  for (b1 in 1:nb_3d) for (l1 in 1:nf_3d){
    mapping[ci,1]= b1
    mapping[ci,2]= l1
    mapping[ci,3]= ci
    ci=ci+1
  }

  map2d = data.table(b_3d = mapping[,1],f_3d = mapping[,2],b_2d=1,f_2d=mapping[,3])
  #print(map2d)
  #transform model

  #tranform means and Sd
  for (b1 in 1:nb_3d) for (l1 in 1:nf_3d){
    b1_2d = map2d[b_3d == b1 & f_3d == l1, b_2d]
    l1_2d = map2d[b_3d == b1 & f_3d == l1, f_2d]
    A1_2d[b1_2d,l1_2d,] = A1_3d[b1,l1,]
    A2_2d[b1_2d,l1_2d,] = A2_3d[b1,l1,]
    S1_2d[b1_2d,l1_2d,] = S1_3d[b1,l1,]
    S2_2d[b1_2d,l1_2d,] = S2_3d[b1,l1,]

  }

  model_2d_init$A1 = A1_2d
  model_2d_init$A2 = A2_2d
  model_2d_init$S1 = S1_2d
  model_2d_init$S2 = S2_2d

  # transform num movers here
  for (b1 in 1:nb_3d) for (l1 in 1:nf_3d) for (b2 in 1:nb_3d) for (l2 in 1:nf_3d) {
    b1_2d = map2d[b_3d == b1 & f_3d == l1, b_2d]
    l1_2d = map2d[b_3d == b1 & f_3d == l1, f_2d]
    b2_2d = map2d[b_3d == b2 & f_3d == l2, b_2d]
    l2_2d = map2d[b_3d == b2 & f_3d == l2, f_2d]
    NNm_2d[b1_2d,b2_2d,l1_2d,l2_2d] = NNm_3d[b1,b2,l1,l2]

  }

  model_2d_init$NNm = NNm_2d


  # transform pk1
  for (b1 in 1:nb_3d) for (l1 in 1:nf_3d) for (b2 in 1:nb_3d) for (l2 in 1:nf_3d) {
    mm = b1 + nb_3d*(b2 -1)
    jj = l1 + nf_3d*(l2 -1)

    b1_2d = map2d[b_3d == b1 & f_3d == l1, b_2d]
    l1_2d = map2d[b_3d == b1 & f_3d == l1, f_2d]
    b2_2d = map2d[b_3d == b2 & f_3d == l2, b_2d]
    l2_2d = map2d[b_3d == b2 & f_3d == l2, f_2d]

    mm_2d = b1_2d + nb_2d*(b2_2d -1)
    jj_2d = l1_2d + nf_2d*(l2_2d -1)

    pk1_2d[mm_2d,jj_2d,] = pk1_3d[mm,jj,]
  }


  model_2d_init$pk1 = pk1_2d

  # transform pk0
  for (b1 in 1:nb_3d) for (l1 in 1:nf_3d){
    b1_2d = map2d[b_3d == b1 & f_3d == l1, b_2d]
    l1_2d = map2d[b_3d == b1 & f_3d == l1, f_2d]
    pk0_2d[b1_2d,l1_2d,] = pk0_3d[b1,l1,]
  }

  model_2d_init$pk0 = pk0_2d

  model_2d = model_2d_init
  return (model_2d)


}

#' Data tranformation ( 3d array to 2d array)
#' @export
m2.mixt.transform.data <- function(model_3d, model_2d_init, data=NA) {

  data_3d <- copy(data)
  sdata_3d <- data_3d$sdata
  jdata_3d <-  data_3d$jdata

  nb_3d  = model_3d$nb
  nf_3d  = model_3d$nf
  nk_3d  = model_3d$nk

  nb_2d  = model_2d_init$nb
  nf_2d  = model_2d_init$nf
  nk_2d  = model_2d_init$nk

  mapping = array(0,c(nb_3d*nf_3d,3))
  # create the mappings
  ci=1
  for (b1 in 1:nb_3d) for (l1 in 1:nf_3d){
    mapping[ci,1]= b1
    mapping[ci,2]= l1
    mapping[ci,3]= ci
    ci=ci+1
    #print(b1)
  }
  map2d = data.table(b_3d = mapping[,1],f_3d = mapping[,2],b_2d=1,f_2d=mapping[,3])
  #print(map2d)
  # transform num movers here
  for (b1 in 1:nb_3d) for (l1 in 1:nf_3d) for (b2 in 1:nb_3d) for (l2 in 1:nf_3d) {
    b1_2d = map2d[b_3d == b1 & f_3d == l1, b_2d]
    l1_2d = map2d[b_3d == b1 & f_3d == l1, f_2d]
    b2_2d = map2d[b_3d == b2 & f_3d == l2, b_2d]
    l2_2d = map2d[b_3d == b2 & f_3d == l2, f_2d]

    # code change: data.table lib style ref multi condition (non pythonic)
    jdata_3d[(m1==b1) & (j1==l1) & (m2==b2) & (j2==l2),m1_2d:=b1_2d]
    jdata_3d[(m1==b1) & (j1==l1) & (m2==b2) & (j2==l2),j1_2d:=l1_2d]
    jdata_3d[(m1==b1) & (j1==l1) & (m2==b2) & (j2==l2),m2_2d:=b2_2d]
    jdata_3d[(m1==b1) & (j1==l1) & (m2==b2) & (j2==l2),j2_2d:=l2_2d]
    #print(jdata_3d)
    #browser()
  }
  jdata_3d[,m1:=m1_2d]
  jdata_3d[,j1:=j1_2d]
  jdata_3d[,m2:=m2_2d]
  jdata_3d[,j2:=j2_2d]
  jdata_3d[,j1true:=j1]
  jdata_3d[,j2true:=j2]
  jdata_3d[,g1true:=m1]
  jdata_3d[,g2true:=m2]

  jdata_3d[,c("m1_2d","j1_2d","m2_2d","j2_2d"):=NULL]

  for (b1 in 1:nb_3d) for (l1 in 1:nf_3d) for (b2 in 1:nb_3d) for (l2 in 1:nf_3d) {
    b1_2d = map2d[b_3d == b1 & f_3d == l1, b_2d]
    l1_2d = map2d[b_3d == b1 & f_3d == l1, f_2d]
    b2_2d = map2d[b_3d == b2 & f_3d == l2, b_2d]
    l2_2d = map2d[b_3d == b2 & f_3d == l2, f_2d]

    # code change: data.table lib style ref multi condition (non pythonic)
    sdata_3d[(m1==b1) & (j1==l1) & (m2==b2) & (j2==l2),m1_2d:=b1_2d]
    sdata_3d[(m1==b1) & (j1==l1) & (m2==b2) & (j2==l2),j1_2d:=l1_2d]
    sdata_3d[(m1==b1) & (j1==l1) & (m2==b2) & (j2==l2),m2_2d:=b2_2d]
    sdata_3d[(m1==b1) & (j1==l1) & (m2==b2) & (j2==l2),j2_2d:=l2_2d]
    #print(jdata_3d)
    #browser()
  }
  sdata_3d[,m1:=m1_2d]
  sdata_3d[,j1:=j1_2d]
  sdata_3d[,m2:=m2_2d]
  sdata_3d[,j2:=j2_2d]
  sdata_3d[,j1true:=j1]
  sdata_3d[,g1true:=m1]

  sdata_3d[,c("m1_2d","j1_2d","m2_2d","j2_2d"):=NULL]

  data_2d = list(sdata= sdata_3d,jdata = jdata_3d)
  return(data_2d)

}

