# conduct permutation test to determine direction between continuous variable and discrete variable
permcontidANM<-function(contiX,Y,number_of_permutations=5000,level=0.05,cycX=1,cycY=1,nbins=10){
  require(infotheo)
  X=unlist(discretize(contiX,nbins=nbins))
  output=fit_both_dir_discrete(X,cycX=cycX,Y,cycY=cycY,level=level)
  P_X2Y=output$p_val_fw
  P_Y2X=output$p_val_bw
  diff_estimated=P_Y2X-P_X2Y
  abs_diff=abs(diff_estimated)
  
  # P value of causation test by permutation test based on difference
  perm_X2Y<-perm_Y2X<-c()
  sizeY=length(Y)
  for (i in 1:number_of_permutations){
    Y_perm=sample(Y,sizeY)
    perm_output=fit_both_dir_discrete(X,cycX=cycX,Y=Y_perm,cycY=cycY,level=level)
    perm_X2Y=c(perm_X2Y,perm_output$p_val_fw)
    perm_Y2X=c(perm_Y2X,perm_output$p_val_bw)
  }
  perm_diff<-perm_Y2X-perm_X2Y
  abs_perm_diff=abs(perm_diff)
  
  Top_Y2X=(sum(perm_diff>diff_estimated)+.5*sum(perm_diff==diff_estimated))/number_of_permutations
  Top_X2Y=1-Top_Y2X
  dir=ifelse(Top_Y2X<level,-1,ifelse(Top_X2Y<level,1,0))
  
  Pc=(sum(abs_perm_diff>abs_diff)+.5*sum(abs_perm_diff==abs_diff))/number_of_permutations
  perm_P_X2Y=(sum(perm_X2Y<P_X2Y)+.5*sum(perm_X2Y==P_X2Y))/number_of_permutations
  perm_P_Y2X=(sum(perm_Y2X<P_Y2X)+.5*sum(perm_Y2X==P_Y2X))/number_of_permutations
  
  list(dir=dir,P_no_causation=Pc,p_val_ind=output$p_val_ind,fct_fw=output$fct_fw,fct_bw=output$fct_bw,P_X2Y=P_X2Y,P_Y2X=P_Y2X,diff_percent_X2Y=Top_X2Y,diff_percent_Y2X=Top_Y2X,perm_P_X2Y=perm_P_X2Y,perm_P_Y2X=perm_P_Y2X)
}

########################### source script ##################################
fit_both_dir_discrete<-function(X,cycX,Y,cycY,level,display=c(FALSE,TRUE)){
  if (cycY==0){
    fit.output=fit_discrete(X,Y,level)
    fct_fw=fit.output$fct
    p_val_fw=fit.output$p_val
  }else if(cycY==1){
    fit.output=fit_discrete_cyclic(X,Y,level)
    fct_fw=fit.output$fct
    p_val_fw=fit.output$p_val
  }
  
  if (cycX==0){
    fit.output=fit_discrete(Y,X,level)
    fct_bw=fit.output$fct
    p_val_bw=fit.output$p_val
  }else if(cycX==1){
    fit.output=fit_discrete_cyclic(Y,X,level)
    fct_bw=fit.output$fct
    p_val_bw=fit.output$p_val
  }
  
  options(warn=-1)
  p_val_ind=ifelse((length(unique(Y))==1|length(unique(X))==1),1,chisq.test(Y,X,correct=FALSE)$p.value)
  if (display==TRUE){
    # p_val_fw
    if (p_val_fw>level){
      cat("fct_fw",fct_fw,"\n")
      cat("ANM could be fitted in the direction X->Y using fct_fw. \n")
    }
    # p_val_bw
    if (p_val_bw>level){
      cat("fct_bw",fct_bw,"\n")
      cat("ANM could be fitted in the direction Y->X using fct_bw. \n")
    }
    if (p_val_bw>level & p_val_fw<level){
      cat("Only one ANM could be fit. The method infers Y->X. \n")
    }else if(p_val_bw<level & p_val_fw>level){
      cat("Only one ANM could be fit. The method infers X->Y. \n")
    }else if(p_val_bw<level & p_val_fw<level){
      cat("No ANM could be fit. The method does not know the causal direction. \n")
    }else{
      cat("Both ANM could be fit. The method does not know the causal direction. \n")
    }
    # are X and Y independent?
    if (p_val_ind>level){
      cat("But note that X and Y are considered to be independent anyway. (Thus no causal relation) \n")
    }
  }
  options(warn=0)
  
  list(fct_fw=fct_fw,fct_bw=fct_bw,p_val_fw=p_val_fw,p_val_bw=p_val_bw,p_val_ind=p_val_ind)
}

#####################
# cyclic
fit_discrete_cyclic<-function(X,Y,level){
  options(warn=-1)
  require(descr)
  # parameter
  num_iter=10
  num_pos_fct=min(max(Y)-min(Y),10)
  
  # rescaling
  # X_new takes values from 1 ... X_new_max
  # Y_values are everything between Y_min and Y_max
  X_values=unique(X)
  Y_values=seq(min(Y),max(Y),by=1)
  
  if (length(X_values)==1|length(Y_values)==1){
    fct=rep(1,length(X_values))*Y_values[1]
    p_val=1
  }else{
    p<-CrossTable(c(X,rep(NA,length(Y_values))),c(Y,Y_values),prop.chisq = FALSE)$t
    
    fct=c()
    cand=list()
    for (i in 1:length(X_values)){
      b=order(p[i,])
      for (k in 1:ncol(p)){
        p[i,k]=ifelse(k==b[length(b)],p[i,k]+1,p[i,k]+1/(2*abs(k-b[length(b)])))
      }
      b=order(p[i,])
      cand[[i]]=b
      fct=c(fct,Y_values[b[length(b)]])
    }
    
    X_new=X
    for (i in 1:nrow(p)){
      X_new[X==rownames(p)[i]]=i
    }
    
    yhat=fct[X_new]
    eps=(Y-yhat)%%(max(Y)-min(Y)+1)
    p_val=ifelse((length(unique(eps))==1),1,chisq.test(eps,X,correct=FALSE)$p.value)
    # correct=TRUE as default; if correct=FALSE, completely consistant to original MATLAB scripts
    i=0
    while(p_val<level & i<num_iter){
      for (j_new in sample.int(length(X_values))){
        pos_fct=list()
        p_val_comp<-p_val_comp2<-c()
        for (j in 1:(num_pos_fct+1)){
          pos_fct[[j]]=fct
          pos_fct[[j]][j_new]=Y_values[cand[[j_new]][length(cand[[j_new]])-(j-1)]]
          yhat=pos_fct[[j]][X_new]
          eps=(Y-yhat)%%(max(Y)-min(Y)+1)
          if (length(unique(eps))==1){
            p_val_comp=c(p_val_comp,1)
            p_val_comp2=c(p_val_comp2,0)
          }else{
            chi_sq=chisq.test(eps,X,correct=FALSE)
            p_val_comp=c(p_val_comp,chi_sq$p.value)
            p_val_comp2=c(p_val_comp2,chi_sq$statistic)
          }
        }
        aa=max(p_val_comp)
        j_max=which(p_val_comp==aa)
        if (aa<exp(-3)){
          j_max=which(p_val_comp2==min(p_val_comp2))
        }
        fct=pos_fct[[min(j_max)]]
        yhat=fct[X_new]
        eps=(Y-yhat)%%(max(Y)-min(Y)+1)
        p_val=ifelse((length(unique(eps))==1),1,chisq.test(eps,X,correct=FALSE)$p.value)
      }
      i=i+1
    }

  }
  options(warn=0)
  list(fct=fct,p_val=p_val)
}

###################
# non_cyclic
fit_discrete<-function(X,Y,level){
  options(warn=-1)
  require(descr)
  # parameter
  num_iter=10
  num_pos_fct=min(max(Y)-min(Y),20)
  
  # rescaling
  # X_new takes values from 1 ... X_new_max
  # Y_values are everything between Y_min and Y_max
  X_values=unique(X)
  Y_values=seq(min(Y),max(Y),by=1)
  
  if (length(X_values)==1|length(Y_values)==1){
    fct=rep(1,length(X_values))*Y_values[1]
    p_val=1
  }else{
    p<-CrossTable(c(X,rep(NA,length(Y_values))),c(Y,Y_values),prop.chisq = FALSE)$t
    
    fct=c()
    cand=list()
    for (i in 1:length(X_values)){
      b=order(p[i,])
      for (k in 1:ncol(p)){
        p[i,k]=ifelse(k==b[length(b)],p[i,k]+1,p[i,k]+1/(2*abs(k-b[length(b)])))
      }
      b=order(p[i,])
      cand[[i]]=b
      fct=c(fct,Y_values[b[length(b)]])
    }
    # the following script more convenient compared to MATLAB
    X_new=X
    for (i in 1:nrow(p)){
      X_new[X==rownames(p)[i]]=i
    }
    
    yhat=fct[X_new]
    eps=Y-yhat
    if (length(unique(eps))==1){
      cat("Warning!!there is a deterministic relation between X and Y \n")
      p_val=1
    }else{
      p_val=chisq.test(eps,X,correct=FALSE)$p.value
    }

    i=0
    while(p_val<level & i<num_iter){
      for (j_new in sample.int(length(X_values))){
        pos_fct=list()
        p_val_comp<-p_val_comp2<-c()
        for (j in 1:(num_pos_fct+1)){
          pos_fct[[j]]=fct
          pos_fct[[j]][j_new]=Y_values[cand[[j_new]][length(cand[[j_new]])-(j-1)]]
          yhat=pos_fct[[j]][X_new]
          eps=Y-yhat
          
          if (length(unique(eps))==1){
            p_val_comp=c(p_val_comp,1)
            p_val_comp2=c(p_val_comp2,0)
          }else{
            chi_sq=chisq.test(eps,X,correct=FALSE)
            p_val_comp=c(p_val_comp,chi_sq$p.value)
            p_val_comp2=c(p_val_comp2,chi_sq$statistic)
          }
        }
        aa=max(p_val_comp)
        j_max=which(p_val_comp==aa)
        if (aa<exp(-3)){
          j_max=which(p_val_comp2==min(p_val_comp2))
        }
        fct=pos_fct[[min(j_max)]]
        yhat=fct[X_new]
        eps=Y-yhat
        p_val=ifelse((length(unique(eps))==1),1,chisq.test(eps,X,correct=FALSE)$p.value)
      }
      i=i+1
    }
    fct=fct+round(mean(eps))

  }
  options(warn=0)
  list(fct=fct,p_val=p_val)
}