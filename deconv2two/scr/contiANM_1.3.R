# conduct permutation test to determine direction
permANM<-function(X,Y,number_of_permutations=5000,level=0.05,dplot=FALSE,dword=FALSE,
                  fit.method=c("ssr","B.spline","loess")){
  P_X2Y=contiANM(X,Y,fit.method=fit.method)$pvalue
  P_Y2X=contiANM(Y,X,fit.method=fit.method)$pvalue
  diff_estimated=P_Y2X-P_X2Y
  abs_diff=abs(diff_estimated)
  
  # P value of causation test by permutation test based on difference
  perm_X2Y<-perm_Y2X<-c()
  for (i in 1:number_of_permutations){
    if (dword){
      if (i %% 10 ==0){
        cat (i,"\n")
      }
    }
    Y_perm=sample(Y,length(Y),FALSE)
    perm_X2Y=c(perm_X2Y,contiANM(X,Y_perm,fit.method=fit.method)$pvalue)
    perm_Y2X=c(perm_Y2X,contiANM(Y_perm,X,fit.method=fit.method)$pvalue)
  }
  
  perm_diff<-perm_Y2X-perm_X2Y
  
  Top_Y2X=(sum(perm_diff>diff_estimated)+.5*sum(perm_diff==diff_estimated))/number_of_permutations
  Top_X2Y=1-Top_Y2X
  
  # summation from right to left
  abs_perm_diff=abs(perm_diff)
  Pc=(sum(abs_perm_diff>abs_diff)+.5*sum(abs_perm_diff==abs_diff))/number_of_permutations
  perm_P_X2Y=(sum(perm_X2Y<P_X2Y)+.5*sum(perm_X2Y==P_X2Y))/number_of_permutations
  perm_P_Y2X=(sum(perm_Y2X<P_Y2X)+.5*sum(perm_Y2X==P_Y2X))/number_of_permutations
  
  if (dplot){
    plot(density(perm_diff),main="permutation of P_Y2X-P_X2Y")
    abline(v=diff_estimated,col="blue")
  }
  if (Top_Y2X<level){
    dir=-1
  }else if(Top_X2Y<level){
    dir=1
  }else(dir=0)
  list(dir=dir,P_X2Y=P_X2Y,P_Y2X=P_Y2X,diff_percent_X2Y=Top_X2Y,diff_percent_Y2X=Top_Y2X,P_no_causation=Pc,perm_P_X2Y=perm_P_X2Y,perm_P_Y2X=perm_P_Y2X)
}

# continuous ANM
contiANM<-function(x,y,m=c(2,3,4),
                      odd=TRUE,N=40,limnla=c(-10,3),fit.method=c("ssr","B.spline","loess")){
  # m: m=2, cubic spline
  #    m=3, quintic spline
  #    m=4, septic spline
  
  # limnla: a vector of length one or two, specifying a 
  #         search range for log10(n*lambda), where lambda is the 
  #         smoothing parameter and n is the sample size. If it is 
  #         a single value, the smoothing parameter will be fixed at 
  #         this value. This option is only applicable to spline 
  #         smoothing with a single smoothing parameter.
  
  require(assist)
  require(psych)
  require(dHSIC)
  require(fANCOVA)

  if(length(m)==3){m<-2}

  n.sub=length(x)
  if (n.sub!=length(y)){
    stop("lengths of x and y do not match")
  }else{
    n.train=floor(n.sub/5*4)
    if (fit.method=="ssr"){
      x=x-min(x)
      y=y-min(y)
    }
    x.train=x[1:n.train]
    y.train=y[1:n.train]
    x.test=x[(n.train+1):n.sub]
    y.test=y[(n.train+1):n.sub]
    m.test=length(x.test)
    
    ##### B spline #####
    if (fit.method=="B.spline"){
      BS=smooth.spline(x.train,y.train,nknots=N)
      new.x.train=x.train[order(x.train)]
      new.fit=predict(BS,new.x.train)
      if (odd==TRUE){
        #lines(new.fit$x, new.fit$y,lty="solid", col="green")
      }else{
        #lines(new.fit$y,new.fit$x, lty="solid", col="darkred")
      }
      fitted <- predict(BS,x.test)$y
      eps=y.test-fitted
    }
    
    ##### smoothing splines #####
    if (fit.method=="ssr"){
      if (m==2){
        B<-ssr(y.train~x.train,rk=cubic2(x.train),limnla=limnla)
      }else if(m==3){
        B<-ssr(y.train~x.train+I(x.train^2),rk=quintic2(x.train),limnla=limnla)
      }else if(m==4){
        B<-ssr(y.train~x.train+I(x.train^2)+I(x.train^3),rk=septic2(x.train),limnla=limnla)
      }
      
      fitted=predict(B,data.frame(x.train=x.test),pstd=FALSE)$fit
      eps=y.test-fitted
      
      # draw fitted curve
      new.x.train=x.train[order(x.train)]
      new.fit=B$fit[order(x.train)]
      if (odd==TRUE){
        #lines(new.x.train, new.fit,lty="solid", col="green")
      }else{
        #lines(new.fit,new.x.train, lty="solid", col="darkred")
      }
      
    }
    
    ##### LOESS #####
    if (fit.method=="loess"){
      new.x.train=x.train[order(x.train)]
      new.y.train=y.train[order(x.train)]
      # training span parameter based on AICC
      span=(loess.as(new.y.train,new.x.train)$pars)$span
      cars.lo <- loess(new.y.train ~ new.x.train,span=span)
      
      if (odd==TRUE){
        #lines(new.x.train, cars.lo$fitted,lty="solid", col="green")
      }else{
        #lines(cars.lo$fitted,new.x.train, lty="solid", col="darkred")
      }
      fitted <- predict(cars.lo,x.test)
      eps=y.test-fitted
    }
 
    ########################
    x.test <- x.test[!is.na(eps)]
    eps <- eps[!is.na(eps)]
    pvalue=dhsic.test(x.test,eps)$p.value
    #pvalue
    statistic=dhsic.test(x.test,eps)$statistic
    #statistic
    list(statistic=statistic,pvalue=pvalue)
  }
}
