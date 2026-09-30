# stage2.R
# Cross-fitted debiased estimate of the average derivative of log gas
# with respect to log income (the income elasticity).

L=5

rrr<-function(Y,X,X.up,X.down,delta,p0,D_LB,D_add,max_iter,alpha_estimator,gamma_estimator,bias){

  n=nrow(X)
  folds <- split(sample(n, n,replace=FALSE), as.factor(1:L))

  Psi_tilde=numeric(0)

  for (l in 1:L){

    idx.l=folds[[l]]

    Y.l=Y[idx.l]
    Y.nl=Y[-idx.l]

    X.l=X[idx.l,,drop=FALSE]
    X.nl=X[-idx.l,,drop=FALSE]

    X.up.l=X.up[idx.l,,drop=FALSE]
    X.up.nl=X.up[-idx.l,,drop=FALSE]

    X.down.l=X.down[idx.l,,drop=FALSE]
    X.down.nl=X.down[-idx.l,,drop=FALSE]

    # drop columns that are constant in the training fold (keep the intercept, column 1)
    keep=which(apply(X.nl,2,sd)>0 | seq_len(ncol(X.nl))==1)
    X.l=X.l[,keep,drop=FALSE];           X.nl=X.nl[,keep,drop=FALSE]
    X.up.l=X.up.l[,keep,drop=FALSE];     X.up.nl=X.up.nl[,keep,drop=FALSE]
    X.down.l=X.down.l[,keep,drop=FALSE]; X.down.nl=X.down.nl[,keep,drop=FALSE]

    n.l=nrow(X.l)
    n.nl=nrow(X.nl)

    # get stage 1 (on nl)
    stage1_estimators<-get_stage1(Y.nl,X.nl,X.up.nl,X.down.nl,delta,p0,D_LB,D_add,max_iter,alpha_estimator,gamma_estimator)
    alpha_hat=stage1_estimators[[1]]
    gamma_hat=stage1_estimators[[2]]

    print(paste0('fold: ',l))

    #get stage 2 (on l)
    Psi_tilde.l=rep(0,n.l)
    for (i in 1:n.l){
      if(bias){ #plug-in
        Psi_tilde.l[i]=psi_tilde_bias(Y.l[i],X.l[i,],X.up.l[i,],X.down.l[i,],delta,m,alpha_hat,gamma_hat)
      }else{ #DML
        Psi_tilde.l[i]=psi_tilde(Y.l[i],X.l[i,],X.up.l[i,],X.down.l[i,],delta,m,alpha_hat,gamma_hat)
      }
    }

    Psi_tilde=c(Psi_tilde,Psi_tilde.l)

  }

  #point estimation: average derivative of log gas wrt log income (income elasticity)
  ate=mean(Psi_tilde)

  #influences
  Psi=Psi_tilde-ate

  var=mean(Psi^2)
  se=sqrt(var/n)

  out<-c(n,ate,se)

  return(out)
}

printer<-function(spec1){
  print(paste("   n:    ",spec1[1],"   AD:    ",round(spec1[2],3), "   SE:   ", round(spec1[3],3),
              "   95% CI: [", round(spec1[2]-1.96*spec1[3],3), ", ", round(spec1[2]+1.96*spec1[3],3), "]", sep=""))
}

for_tex<-function(spec1){
  print(paste(" & ",spec1[1]," & ",round(spec1[2],3), "   &   ", round(spec1[3],3), sep=""))
}
