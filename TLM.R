

# packages <- c("tensorflow","reshape2","telefit","reticulate","keras","matrixNormal","survival","rTensor","pracma","parallel", "foreach", "snow", "doSNOW","tictoc","TRES")
packages <- c("reshape2","telefit","matrixNormal","survival","rTensor","pracma","parallel", "foreach", "snow", "doSNOW","tictoc","TRES")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


# transform to length 1 unit vector
unit_vec <- function(vec){
  return(vec/(sqrt(sum(vec**2,na.rm=TRUE))))
}

# AR structure matrix
genAR<-function(dim,corr){
  temp.mat<-matrix(0,dim,dim)
  for (i in 1:dim){
    for (j in 1:dim){
      if (i==j) {temp.mat[i,j]=1
      }
      else {temp.mat[i,j]=corr^abs(i-j)
      }
    }
  }	
  return(temp.mat)
}

# Breslow function
breslow<- function(time, status, X, B,z,eta){
  # z: n x pz, eta: pz vector
  pz=dim(z)[2]
  data <-data.frame(time, status,z, X)
  data <-data[order(data$time), ]
  t<-unique(data$time)
  k<-length(t)
  h<- rep(0,k)
  for(i in 1:k) {
    lp <- (data.matrix(data[,c((2+1):(2+pz))])%*%eta+data.matrix(data[,-c(1:(2+pz))]) %*% B)[data$time>=t[i]]
    risk <- exp(lp)
    h[i] <- sum(data$status[data$time==t[i]]) / sum(risk)
  }
  res <- cumsum(h)
  output<-cbind(time,time,time)
  irep=1
  for (irep in 1:k){
    pos<-which(output[,1]==t[irep],arr.ind=T)
    output[pos,2]<-res[irep] # \Lambda_{0}
    output[pos,3]<-h[irep] # delta \Lambda_{0}
  }
  return(output)
}

TLM<-function(X,z,R,Id,time,e,max.iter=500,tol=1e-2){
  ##############
  # Start from here
  data.tensor<-as.tensor(X)
  n<-dim(X)[1]
  p1<-dim(X)[2]
  p2<-dim(X)[3]
  pz<-dim(z)[2]
  Beta.est = array(0,dim=c(p1,p2,R))
  Beta.est.r = array(0,dim=c(p1,p2,R))
  conv.mat<-matrix(0,1,(5+pz))
  colnames(conv.mat)<-c("Rank","niter","FullLiklihood","Lik.Diff.","Scalar.update",c(rep("eta",pz)))
  
  
  r=1
  Beta.update<-matrix(0,p1,p2)
  AIC<-rep(0,R)
  BIC<-rep(0,R)
  
  eta.init<-matrix(rep(0,pz),ncol=1)
  eta.update<-eta.init
  
  
  mB1.init=matrix(rnorm(p1,0,0.01),ncol=1)
  mB2.init=matrix(rnorm(p2,0,0.01),ncol=1)
  scalar.init<-1
  Beta.init<-scalar.init*mB1.init%*%t(mB2.init)
  Beta.update<-Beta.init
  r=1
  for (r in 1:R){
    
    # initialization the rank r coefficient for CP decomposition
    
    
    scalar.update<-scalar.init
    mB1.update<-mB1.init
    mB2.update<-mB2.init
    
    norm.diff=Inf;max.iter=max.iter;tol=tol
    niter=1
    while(norm.diff > tol && niter < max.iter){
      
      # transformed X1
      temp.x1<-matrix(0,n,p1)
      X.tensor=rTensor::as.tensor(X)
      temp1<-scalar.update*rTensor::ttm(X.tensor,t(mB2.update),m=3)
      temp.x1<-temp1[,,1]@data

      tt<-rTensor::k_unfold(X.tensor,m=1)
      temp.obs<-tt@data
      
      breslow.est<-breslow(time=time, status=e, X=temp.obs, B=matrix(Beta.update+scalar.update*(mB1.update%*%t(mB2.update)),ncol=1),z=z,eta=eta.update)
      lambda0<-breslow.est[,2]
      
      mu<-lambda0*exp(z%*%eta.update+rTensor::ttm(rTensor::k_unfold(X.tensor,m=1),matrix(Beta.update,1),m=2)@data+rTensor::ttm(rTensor::k_unfold(X.tensor,m=1),matrix((scalar.update*mB1.update%*%t(mB2.update)),1),m=2)@data)
      K=z%*%eta.update+rTensor::ttm(rTensor::k_unfold(X.tensor,m=1),matrix(Beta.update,1),m=2)@data+rTensor::ttm(rTensor::k_unfold(X.tensor,m=1),matrix((scalar.update*mB1.update%*%t(mB2.update)),1),m=2)@data
    
      mu<-matrix(mu,ncol=1)
      # mu[which(mu==0,arr.ind=T)[,1]]<-0.001
      mu[which(mu <1e-3,arr.ind=T)[,1]]<-0.01
      K<-matrix(K,ncol=1)
      V<-diag(as.vector(mu))
      V.inv<-diag(1/as.vector(mu))
      
      f=K+(V.inv%*%(as.matrix(e)-mu))
      # f[is.nan(f)]<-200
      f[which(f >1e+5,arr.ind=T)[,1]]<-100
      
      mB1.update2=MASS::ginv(t(temp.x1)%*%V%*%temp.x1)%*%t(temp.x1)%*%V%*%f
      
      scalar1.mat<-mB1.update2/unit_vec(mB1.update2)
      scalar1<-scalar1.mat[1,1]
      mB1.update2<-unit_vec(mB1.update2)
      
      temp.x2<-matrix(0,n,p2)
      temp2<-scalar.update*rTensor::ttm(X.tensor,t(mB1.update2),m=2)
      temp.x2<-temp2[,1,]@data
      
      Beta.update1<-Beta.update+(scalar.update*mB1.update2%*%t(mB2.update))
      breslow.est2<-breslow(time=time, status=e, X=temp.obs, B=matrix(Beta.update1,ncol=1),z=z,eta=eta.update)
      lambda02<-breslow.est2[,2]
      
      mu2<-lambda02*exp(z%*%eta.update+rTensor::ttm(rTensor::k_unfold(X.tensor,m=1),matrix(Beta.update1,1),m=2)@data)
      K2=z%*%eta.update+rTensor::ttm(rTensor::k_unfold(X.tensor,m=1),matrix(Beta.update1,1),m=2)@data
      
      mu2<-matrix(mu2,ncol=1)
      # mu2[which(mu2==0,arr.ind=T)[,1]]<-0.0001
      mu2[which(mu2 <1e-3,arr.ind=T)[,1]]<-0.01
      K2<-matrix(K2,ncol=1)
      V2<-diag(as.vector(mu2))
      V2.inv<-diag(1/as.vector(mu2))
      
      f2=K2+(V2.inv%*%(as.matrix(as.matrix(e))-mu2))
      # f2[is.nan(f2)]<-200
      f2[which(f2 >1e+5,arr.ind=T)[,1]]<-100
      
      mB2.update2=MASS::ginv(t(temp.x2)%*%V2%*%temp.x2)%*%t(temp.x2)%*%V2%*%f2
      
      scalar2.mat<-mB2.update2/unit_vec(mB2.update2)[1]
      scalar2<-scalar2.mat[1,1]
      mB2.update2<-unit_vec(mB2.update2)
      
      # update scalar
      temp.x3<-matrix(0,n,1)
      temp3<-rTensor::ttm(rTensor::ttm(X.tensor,t(mB1.update2),m=2),t(mB2.update2),m=3)
      temp.x3<-temp3[,1,1]@data
      
      
      Beta.update3<-Beta.update+(mB1.update2%*%t(mB2.update2))
      breslow.est3<-breslow(time=time, status=e, X=temp.obs, B=matrix(Beta.update3,ncol=1),z=z,eta=eta.update)
      lambda03<-breslow.est3[,2]
      
      mu3<-lambda03*exp(z%*%eta.update+rTensor::ttm(rTensor::k_unfold(X.tensor,m=1),matrix(Beta.update3,1),m=2)@data)
      K3<-z%*%eta.update+rTensor::ttm(rTensor::k_unfold(X.tensor,m=1),matrix(Beta.update3,1),m=2)@data
      
      mu3<-matrix(mu3,ncol=1)
      # mu3[which(mu3==0,arr.ind=T)[,1]]<-0.0001
      mu3[which(mu3 <1e-3,arr.ind=T)[,1]]<-0.01
      K3<-matrix(K3,ncol=1)
      V3<-diag(as.vector(mu3))
      V3.inv<-diag(1/as.vector(mu3))
      
      f3=K3+(V3.inv%*%(as.matrix(as.matrix(e))-mu3))
      
      # f3[is.nan(f3)]<-200
      f3[which(f3 >1e+5,arr.ind=T)[,1]]<-100
      
      scalar.update2=MASS::ginv(t(temp.x3)%*%V3%*%temp.x3)%*%t(temp.x3)%*%V3%*%f3
      scalar.update2<-as.numeric(scalar.update2)
      
      # eta update
      Beta.update4<-Beta.update+(scalar.update2*mB1.update2%*%t(mB2.update2))
      breslow.est4<-breslow(time=time, status=e, X=temp.obs, B=matrix(Beta.update4,ncol=1),z=z,eta=eta.update)
      lambda04<-breslow.est4[,2]
      
      mu4<-lambda04*exp(z%*%eta.update+rTensor::ttm(rTensor::k_unfold(X.tensor,m=1),matrix(Beta.update4,1),m=2)@data)
      
      mu4<-matrix(mu4,ncol=1)
      # mu4[which(mu4==0,arr.ind=T)[,1]]<-0.0001
      mu4[which(mu4 <1e-3,arr.ind=T)[,1]]<-0.01
      V4<-diag(as.vector(mu4))
      V4.inv<-diag(1/as.vector(mu4))
      
      
      f4=z%*%eta.update+(V4.inv%*%(as.matrix(as.matrix(e))-mu4))
      
      f4[which(f4 >1e+3,arr.ind=T)[,1]]<-100
      
      eta.update2=MASS::ginv(t(z)%*%V4%*%z)%*%t(z)%*%V4%*%f4
      
      # calculate convergence
      # Convergence of the Frobenius norm difference
      # Beta.update2<-Beta.update+(scalar1*scalar2)*(mB1.update2%*%t(mB2.update2))
      Beta.new<-Beta.update+(scalar.update2*mB1.update2%*%t(mB2.update2))
      Beta.old<-Beta.update+(scalar.update*mB1.update%*%t(mB2.update))
      
      # Convergence of the partial likelihood
      breslow.est.new<-breslow(time=time, status=e, X=temp.obs, B=matrix(Beta.new,ncol=1),z=z,eta=eta.update2)
      breslow.est.old<-breslow(time=time, status=e, X=temp.obs, B=matrix(Beta.old,ncol=1),z=z,eta=eta.update)
      
      
      fullik.new<-sum(log(((breslow.est.new[,3]*exp(z%*%eta.update2+temp.obs%*%as.vector(Beta.new)))^{e})*(exp(-breslow.est.new[,2]*exp(z%*%eta.update2+temp.obs%*%as.vector(Beta.new))))))
      fullik.old<-sum(log(((breslow.est.old[,3]*exp(z%*%eta.update+temp.obs%*%as.vector(Beta.old)))^{e})*(exp(-breslow.est.old[,2]*exp(z%*%eta.update+temp.obs%*%as.vector(Beta.old))))))
      fullik.diff=-(fullik.new-fullik.old)
      norm.diff<-abs(fullik.diff)
      
      conv.mat<-rbind(conv.mat,c(as.integer(r),niter,fullik.new,fullik.diff,scalar.update2,eta.update2))
      
      # iteration output 
      mB1.update<-mB1.update2
      mB2.update<-mB2.update2
      scalar.update<-scalar.update2
      eta.update<-eta.update2
      norm.diff<-round(norm.diff,4)
      
      niter=niter+1
    }
    
    
    #AIC
    # K=r*(p1+p2)-r^2
    K=r*(p1+p2)+pz
    AIC[r]=round(-2*fullik.new+2*K,4)
    #BIC
    BIC[r]=round(-2*fullik.new+K*log(sum(e)),4)
    
    Beta.est[,,r]=scalar.update*mB1.update%*%t(mB2.update)
    Beta.update<-Beta.update+Beta.est[,,r]
    eta.update<-eta.update
    r=r+1
  }
  Beta.hat<-matrix(0,p1,p2)
  for (rrep in 1:R){
    Beta.hat<-Beta.hat+Beta.est[,,rrep]
    Beta.est.r[,,rrep]<-Beta.hat
  }
  output<-list(R=R,eta.hat=eta.update,Beta.hat=Beta.hat,logLikelihood=fullik.new,Beta.hat.r=Beta.est.r,AIC=AIC,BIC=BIC)
  return(output)
}


# function for Cindex 
# input multivariate form of X,z
Cindex<-function(X,z,Beta,eta,time,status){
  library(survival)
  y <- Surv(time,status)
  p <- ncol(y)
  time <- y[, p - 1]
  status <- y[, p]
  x <- z%*%matrix(eta,ncol=1)+X%*%matrix(Beta,ncol=1)
  n <- length(time)
  ord <- order(time, -status)
  time <- time[ord]
  status <- status[ord]
  x <- x[ord]
  wh <- which(status == 1)
  total <- concordant <- 0
  for (i in wh) {
    for (j in ((i + 1):n)) {
      if (time[j] > time[i]) {
        total <- total + 1
        if (x[j] < x[i]) 
          concordant <- concordant + 1
        if (x[j] == x[i]) 
          concordant <- concordant + 0.5
      }
    }
  }
  return(list(concordant = concordant, total = total, cindex = concordant/total))
}

Cindex2<-function(X,Beta,time,status){
  library(survival)
  y <- Surv(time,status)
  p <- ncol(y)
  time <- y[, p - 1]
  status <- y[, p]
  x <- X%*%matrix(Beta,ncol=1)
  n <- length(time)
  ord <- order(time, -status)
  time <- time[ord]
  status <- status[ord]
  x <- x[ord]
  wh <- which(status == 1)
  total <- concordant <- 0
  for (i in wh) {
    for (j in ((i + 1):n)) {
      if (time[j] > time[i]) {
        total <- total + 1
        if (x[j] < x[i]) 
          concordant <- concordant + 1
        if (x[j] == x[i]) 
          concordant <- concordant + 0.5
      }
    }
  }
  return(list(concordant = concordant, total = total, cindex = concordant/total))
}

# function for the IPCW Brier score
IPCW_brier<-function(t_p,X,z,Beta,eta,time,status){
  n=dim(X)[1] # set sample size
  # Create the censoring survival object
  censoring_surv <- Surv(time, 1 - status)
  # Fit the Kaplan-Meier estimator for the censoring distribution
  km_censoring_fit <- survfit(censoring_surv ~ 1)
  # Function to extract survival probability at a specific time
  get_surv_prob <- function(time, km_fit) {
    surv_probs <- summary(km_fit, times = time)$surv
    if (length(surv_probs) == 0) return(NA) else return(surv_probs)
  }
  # Calculate IPCW weights
  ipcw_weight <- sapply(time, get_surv_prob, km_fit = km_censoring_fit)
  ipcw_weight <- 1 / ipcw_weight
  ipcw_weight[is.infinite(ipcw_weight)] <- NA 
  
  # if (length(which(ipcw_weight==Inf,arr.ind=T))!=0){
  #   ipcw_weight[which(ipcw_weight==Inf,arr.ind=T)]=sort(ipcw_weight)[n-length(which(ipcw_weight==Inf,arr.ind=T))]
  # }
  
  surv_prob<-exp(-breslow(time, status, X, Beta,z,eta)[,2]*exp(z%*%matrix(eta,ncol=1)+X%*%matrix(Beta,ncol=1)))
  event_indicators <- as.numeric(time <= t_p & status == 1)
  ipcw_brier<-mean(status*ipcw_weight * (surv_prob - event_indicators)^2, na.rm = TRUE)
  
  
  output=list(ipcw_brier=ipcw_brier)
  return(output)
}


IPCW_brier2<-function(t_p,X,z,Beta,eta,time,status){
  n=dim(X)[1] # set sample size
  # Create the censoring survival object
  censoring_surv <- Surv(time, 1 - status)
  # Fit the Kaplan-Meier estimator for the censoring distribution
  km_censoring_fit <- survfit(censoring_surv ~ 1)
  # Function to extract survival probability at a specific time
  get_surv_prob <- function(time, km_fit) {
    surv_probs <- summary(km_fit, times = time)$surv
    if (length(surv_probs) == 0) return(NA) else return(surv_probs)
  }
  
  
  ipcw_weight<-rep(0,n)
  for (nrep in 1:n){
    if (time[nrep] <= t_p){
      ipcw_weight[nrep]=status[nrep]/get_surv_prob(time[nrep], km_censoring_fit)
    } else if (time[nrep] > t_p){
      ipcw_weight[nrep]=1/get_surv_prob(t_p, km_censoring_fit)
    }
  }
  
  ipcw_weight[is.infinite(ipcw_weight)] <- NA 
  
  # if (length(which(ipcw_weight==Inf,arr.ind=T))!=0){
  #   ipcw_weight[which(ipcw_weight==Inf,arr.ind=T)]=sort(ipcw_weight)[n-length(which(ipcw_weight==Inf,arr.ind=T))]
  # }
  
  surv_prob<-exp(-breslow(t_p, status, X, Beta,z,matrix(eta,ncol=1))[,2]*exp(z%*%matrix(eta,ncol=1)+X%*%matrix(Beta,ncol=1)))
  event_indicators <- as.numeric(time > t_p)
  ipcw_brier<-mean(status*ipcw_weight * (surv_prob - event_indicators)^2, na.rm = TRUE)
  
  
  output=list(ipcw_brier=ipcw_brier)
  return(output)
}

# function for the Brier score
brier<-function(t_p,X,z,Beta,eta,time,status){
  n=dim(X)[1] # set sample size
  
  
  surv_prob<-exp(-breslow(time, status, X, Beta,z,eta)[,2]*exp(z%*%matrix(eta,ncol=1)+X%*%matrix(Beta,ncol=1)))
  
  brier<-mean((surv_prob-status)^2,na.rm=TRUE)
  
  
  output=list(brier=brier)
  return(output)
}

# function for indexing for a cross-validation
cv.index<-function(Data,nfolds){
  # stratified sampling 
  unique.e1<-Data$ID[which(Data$status==1,arr.ind=T)] # subject observed to fail during the study
  unique.e0<-Data$ID[which(Data$status==0,arr.ind=T)] # subject observed to fail during the study
  
  # indexing for k-fold cross-validation
  n.e1<-length(unique.e1)
  n.e0<-length(unique.e0)
  
  id.folds.e1 <- cut(seq(1:n.e1), breaks = nfolds, labels = 1:nfolds)
  id.folds.e1 <- sample(id.folds.e1, n.e1, replace = FALSE)
  id.folds.e1 <- as.numeric(id.folds.e1)
  fold.index.e1<-data.frame(uniq.id=unique.e1,id.folds=id.folds.e1)
  
  id.folds.e0 <- cut(seq(1:n.e0), breaks = nfolds, labels = 1:nfolds)
  id.folds.e0 <- sample(id.folds.e0, n.e0, replace = FALSE)
  id.folds.e0 <- as.numeric(id.folds.e0)
  fold.index.e0<-data.frame(uniq.id=unique.e0,id.folds=id.folds.e0)
  
  output<-list(fold.index.e1=fold.index.e1,fold.index.e0=fold.index.e0)
  return(output)
}