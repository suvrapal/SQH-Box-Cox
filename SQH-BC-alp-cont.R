# generating data from Box-Cox for alpha in (0,1] using continuous covariate

library(numDeriv)

data_gen_cont_Wei=function(m,g1,g2,p01,p02,c,x_min,x_max,alpha){
  
  b1= (log((1/alpha)*(((1/p01)^alpha)-1))-log((1/alpha)*(((1/p02)^alpha)-1)))/(x_max-x_min)
  b0= log((1/alpha)*(((1/p01)^alpha)-1))-(x_max*b1)
  
  t=rep(NA,m)
  d=rep(NA,m)
  x=runif(m,min=0.1,max=20)
  
  phi=exp(b0+b1*x)/(1+(alpha*exp(b0+b1*x)))
  p0=(1-(alpha*phi))^(1/alpha)
  
  C=rexp(m,rate=c)
  
  countp=0
  countc=0
  
  for(i in 1:m){
    U=runif(1,min=0,max=1)
    if(U<=p0[i]){
      t[i]=C[i]
      d[i]=0
      countp=countp+1
    }else{
      z=qweibull((1-((p0[i]+((1-p0[i])*runif(1,min=0,max=1)))^alpha))/(alpha*phi[i]),shape=1/g1,scale=1/g2)
      t[i]=min(z,C[i])
      if(min(z,C[i])==z){
        d[i]=1
      }else{
        d[i]=0
        countc=countc+1
      }
    }
  }
  
  return(data.frame(t,d,x))
  
}

# log-likelihood function for SQH

log.lik = function(b10,b11,g1,g2,alp,data_obs,data_cens,xt,xc){ 
  
  Fc=1-exp(-((g2*data_cens)^(1/g1)))
  
  phic=exp(b10+b11*xc)/(1+(alp*exp(b10+b11*xc)))
  Ft=1-exp(-((g2*data_obs)^(1/g1)))
  phit=exp(b10+b11*xt)/(1+(alp*exp(b10+b11*xt)))
  
  ft=(1/(data_obs*g1))*((g2*data_obs)^(1/g1))*(1-Ft)
  
  Spc=(1-(alp*phic*Fc))^(1/alp)
  Spt=(1-(alp*phit*Ft))^(1/alp)
  fpt=Spt*phit*ft*((1-(alp*phit*Ft))^(-1))
  
  log.lik=sum(log(fpt))+sum(log(Spc))
  
  return(log.lik)
  
}


# Augmented log-likelihood w.r.t. b10

l.aug.b10 = function(par=c(bb10),b10,b11,g1,g2,alp,data_obs,data_cens,xt,xc,eps){ 
  
  Fc=1-exp(-((g2*data_cens)^(1/g1)))
  
  phic=exp(par[1]+b11*xc)/(1+(alp*exp(par[1]+b11*xc)))
  Ft=1-exp(-((g2*data_obs)^(1/g1)))
  phit=exp(par[1]+b11*xt)/(1+(alp*exp(par[1]+b11*xt)))
  
  ft=(1/(data_obs*g1))*((g2*data_obs)^(1/g1))*(1-Ft)
  
  Spc=(1-(alp*phic*Fc))^(1/alp)
  Spt=(1-(alp*phit*Ft))^(1/alp)
  fpt=Spt*phit*ft*((1-(alp*phit*Ft))^(-1))
  
  log.lik.aug.b10 = sum(log(fpt))+sum(log(Spc)) - (eps*(par[1]-b10)^2)
  
  return(-log.lik.aug.b10)
  
}


# Augmented log-likelihood w.r.t. b11


l.aug.b11 = function(par=c(bb11),b10,b11,g1,g2,alp,data_obs,data_cens,xt,xc,eps){ 
  
  Fc=1-exp(-((g2*data_cens)^(1/g1)))
  
  phic=exp(b10+par[1]*xc)/(1+(alp*exp(b10+par[1]*xc)))
  Ft=1-exp(-((g2*data_obs)^(1/g1)))
  phit=exp(b10+par[1]*xt)/(1+(alp*exp(b10+par[1]*xt)))
  
  ft=(1/(data_obs*g1))*((g2*data_obs)^(1/g1))*(1-Ft)
  
  Spc=(1-(alp*phic*Fc))^(1/alp)
  Spt=(1-(alp*phit*Ft))^(1/alp)
  fpt=Spt*phit*ft*((1-(alp*phit*Ft))^(-1))
  
  log.lik.aug.b11 = sum(log(fpt))+sum(log(Spc)) - (eps*(par[1]-b11)^2)
  
  return(-log.lik.aug.b11)
  
}

# Augmented log-likelihood w.r.t. g1

l.aug.g1 = function(par=c(gg1),b10,b11,g1,g2,alp,data_obs,data_cens,xt,xc,eps){ 
  
  Fc=1-exp(-((g2*data_cens)^(1/par[1])))
  
  phic=exp(b10+b11*xc)/(1+(alp*exp(b10+b11*xc)))
  Ft=1-exp(-((g2*data_obs)^(1/par[1])))
  phit=exp(b10+b11*xt)/(1+(alp*exp(b10+b11*xt)))
  
  ft=(1/(data_obs*par[1]))*((g2*data_obs)^(1/par[1]))*(1-Ft)
  
  Spc=(1-(alp*phic*Fc))^(1/alp)
  Spt=(1-(alp*phit*Ft))^(1/alp)
  fpt=Spt*phit*ft*((1-(alp*phit*Ft))^(-1))
  
  log.lik.aug.g1 = sum(log(fpt))+sum(log(Spc)) - (eps*(par[1]-g1)^2)
  
  return(-log.lik.aug.g1)
  
}


# Augmented log-likelihood w.r.t. g2


l.aug.g2 = function(par=c(gg2),b10,b11,g1,g2,alp,data_obs,data_cens,xt,xc,eps){ 
  
  Fc=1-exp(-((par[1]*data_cens)^(1/g1)))
  
  phic=exp(b10+b11*xc)/(1+(alp*exp(b10+b11*xc)))
  Ft=1-exp(-((par[1]*data_obs)^(1/g1)))
  phit=exp(b10+b11*xt)/(1+(alp*exp(b10+b11*xt)))
  
  ft=(1/(data_obs*g1))*((par[1]*data_obs)^(1/g1))*(1-Ft)
  
  Spc=(1-(alp*phic*Fc))^(1/alp)
  Spt=(1-(alp*phit*Ft))^(1/alp)
  fpt=Spt*phit*ft*((1-(alp*phit*Ft))^(-1))
  
  log.lik.aug.g2 = sum(log(fpt))+sum(log(Spc)) - (eps*(par[1]-g2)^2)
  
  return(-log.lik.aug.g2)
  
}

# Augmented log-likelihood w.r.t. alpha


l.aug.alp = function(par=c(aalp),b10,b11,g1,g2,alp,data_obs,data_cens,xt,xc,eps){ 
  
  Fc=1-exp(-((g2*data_cens)^(1/g1)))
  
  phic=exp(b10+b11*xc)/(1+(par[1]*exp(b10+b11*xc)))
  Ft=1-exp(-((g2*data_obs)^(1/g1)))
  phit=exp(b10+b11*xt)/(1+(par[1]*exp(b10+b11*xt)))
  
  ft=(1/(data_obs*g1))*((g2*data_obs)^(1/g1))*(1-Ft)
  
  Spc=(1-(par[1]*phic*Fc))^(1/par[1])
  Spt=(1-(par[1]*phit*Ft))^(1/par[1])
  fpt=Spt*phit*ft*((1-(par[1]*phit*Ft))^(-1))
  
  log.lik.aug.alp = sum(log(fpt))+sum(log(Spc)) - (eps*(par[1]-alp)^2)
  
  return(-log.lik.aug.alp)
  
}

# Function to perform one-dimensional maximization

max.fn = function(p.init,b10,b11,g1,g2,alp,data_obs,data_cens,xt,xc,eps){
  
  # Augmented log-likelihood w.r.t. b10
  l.aug.b10 = function(par=c(bb10)){ 
    
    Fc=1-exp(-((g2*data_cens)^(1/g1)))
    
    phic=exp(par[1]+b11*xc)/(1+(alp*exp(par[1]+b11*xc)))
    Ft=1-exp(-((g2*data_obs)^(1/g1)))
    phit=exp(par[1]+b11*xt)/(1+(alp*exp(par[1]+b11*xt)))
    
    ft=(1/(data_obs*g1))*((g2*data_obs)^(1/g1))*(1-Ft)
    
    Spc=(1-(alp*phic*Fc))^(1/alp)
    Spt=(1-(alp*phit*Ft))^(1/alp)
    fpt=Spt*phit*ft*((1-(alp*phit*Ft))^(-1))
    
    log.lik.aug.b10 = sum(log(fpt))+sum(log(Spc)) - (eps*(par[1]-b10)^2)
    
    return(-log.lik.aug.b10)
    
  }
  
  # Augmented log-likelihood w.r.t. b11
  l.aug.b11 = function(par=c(bb11)){ 
    
    Fc=1-exp(-((g2*data_cens)^(1/g1)))
    
    phic=exp(b10+par[1]*xc)/(1+(alp*exp(b10+par[1]*xc)))
    Ft=1-exp(-((g2*data_obs)^(1/g1)))
    phit=exp(b10+par[1]*xt)/(1+(alp*exp(b10+par[1]*xt)))
    
    ft=(1/(data_obs*g1))*((g2*data_obs)^(1/g1))*(1-Ft)
    
    Spc=(1-(alp*phic*Fc))^(1/alp)
    Spt=(1-(alp*phit*Ft))^(1/alp)
    fpt=Spt*phit*ft*((1-(alp*phit*Ft))^(-1))
    
    log.lik.aug.b11 = sum(log(fpt))+sum(log(Spc)) - (eps*(par[1]-b11)^2)
    
    return(-log.lik.aug.b11)
    
  }
  
  # Augmented log-likelihood w.r.t. g1
  l.aug.g1 = function(par=c(gg1)){ 
    
    Fc=1-exp(-((g2*data_cens)^(1/par[1])))
    
    phic=exp(b10+b11*xc)/(1+(alp*exp(b10+b11*xc)))
    Ft=1-exp(-((g2*data_obs)^(1/par[1])))
    phit=exp(b10+b11*xt)/(1+(alp*exp(b10+b11*xt)))
    
    ft=(1/(data_obs*par[1]))*((g2*data_obs)^(1/par[1]))*(1-Ft)
    
    Spc=(1-(alp*phic*Fc))^(1/alp)
    Spt=(1-(alp*phit*Ft))^(1/alp)
    fpt=Spt*phit*ft*((1-(alp*phit*Ft))^(-1))
    
    log.lik.aug.g1 = sum(log(fpt))+sum(log(Spc)) - (eps*(par[1]-g1)^2)
    
    return(-log.lik.aug.g1)
    
  }
  
  # Augmented log-likelihood w.r.t. g2
  l.aug.g2 = function(par=c(gg2)){ 
    
    Fc=1-exp(-((par[1]*data_cens)^(1/g1)))
    
    phic=exp(b10+b11*xc)/(1+(alp*exp(b10+b11*xc)))
    Ft=1-exp(-((par[1]*data_obs)^(1/g1)))
    phit=exp(b10+b11*xt)/(1+(alp*exp(b10+b11*xt)))
    
    ft=(1/(data_obs*g1))*((par[1]*data_obs)^(1/g1))*(1-Ft)
    
    Spc=(1-(alp*phic*Fc))^(1/alp)
    Spt=(1-(alp*phit*Ft))^(1/alp)
    fpt=Spt*phit*ft*((1-(alp*phit*Ft))^(-1))
    
    log.lik.aug.g2 = sum(log(fpt))+sum(log(Spc)) - (eps*(par[1]-g2)^2)
    
    return(-log.lik.aug.g2)
    
  }
  
  # Augmented log-likelihood w.r.t. alpha
  l.aug.alp = function(par=c(aalp)){ 
    
    Fc=1-exp(-((g2*data_cens)^(1/g1)))
    
    phic=exp(b10+b11*xc)/(1+(par[1]*exp(b10+b11*xc)))
    Ft=1-exp(-((g2*data_obs)^(1/g1)))
    phit=exp(b10+b11*xt)/(1+(par[1]*exp(b10+b11*xt)))
    
    ft=(1/(data_obs*g1))*((g2*data_obs)^(1/g1))*(1-Ft)
    
    Spc=(1-(par[1]*phic*Fc))^(1/par[1])
    Spt=(1-(par[1]*phit*Ft))^(1/par[1])
    fpt=Spt*phit*ft*((1-(par[1]*phit*Ft))^(-1))
    
    log.lik.aug.alp = sum(log(fpt))+sum(log(Spc)) - (eps*(par[1]-alp)^2)
    
    return(-log.lik.aug.alp)
    
  }
  
  # b10.new = tryCatch({optim(par=p.init[1],fn=l.aug.b10, method="Nelder-Mead")$par
  
  b10.new = tryCatch({nlm(p=p.init[1],f=l.aug.b10)$estimate
  },error=function(e){
    b10.new = c(0)
    return(b10.new)
  }
  )
  
  b11.new = tryCatch({nlm(p=p.init[2],f=l.aug.b11)$estimate
  },error=function(e){
    b11.new = c(0)
    return(b11.new)
  }
  )
  
  g1.new = tryCatch({nlm(p=p.init[3],f=l.aug.g1)$estimate
  },error=function(e){
    g1.new = c(0)
    return(g1.new)
  }
  )
  
  g2.new = tryCatch({nlm(p=p.init[4],f=l.aug.g2)$estimate
  },error=function(e){
    g2.new = c(0)
    return(g2.new)
  }
  )
  
  alp.new = tryCatch({nlm(p=p.init[5],f=l.aug.alp)$estimate
  },error=function(e){
    alp.new = c(0)
    return(alp.new)
  }
  )
  
  out = c(b10.new,b11.new,g1.new,g2.new,alp.new)
  return(out)
}


# Development of the SQH algorithm 

SQH_BC_Gen_Wei=function(data_obs,data_cens,xt,xc,tol,maxit,b10,b11,g1,g2,alp,epsilon1,lambda1,eta1,rho1){
  
  p.new=rep(0,5)
  p.old=rep(0,5)
  
  p.old = c(b10,b11,g1,g2,alp)
  eps = epsilon1
  
  continue = TRUE
  iter = 1
  
  #p.new = max.fn(p.init=p.old,b10=p.old[1],b11=p.old[2],g1=p.old[3],g2=p.old[4],alp=p.old[5],data_obs=data_obs,data_cens=data_cens,xt=xt,xc=xc,eps=eps)
  
  while(continue){
   #print("iter")
     #print(iter)
    
    p.new = max.fn(p.init=p.old,b10=p.old[1],b11=p.old[2],g1=p.old[3],g2=p.old[4],alp=p.old[5],data_obs=data_obs,data_cens=data_cens,xt=xt,xc=xc,eps=eps)
    
    if(0%in%p.new | p.new[3]<0 | p.new[4]<0 | p.new[5]<0 | p.new[5]>1){
      p.new = c(0,0,0,0,0)
      continue = FALSE
    }else{
      tau = sum((p.new-p.old)^2)
      if((log.lik(b10=p.new[1],b11=p.new[2],g1=p.new[3],g2=p.new[4],alp=p.new[5],data_obs=data_obs,data_cens=data_cens,xt=xt,xc=xc) -
             log.lik(b10=p.old[1],b11=p.old[2],g1=p.old[3],g2=p.old[4],alp=p.old[5],data_obs=data_obs,data_cens=data_cens,xt=xt,xc=xc)) < (rho1*tau)){
        eps = eps*lambda1
        #p.new = max.fn(p.init=p.old,b10=p.old[1],b11=p.old[2],g1=p.old[3],g2=p.old[4],alp=p.old[5],data_obs=data_obs,data_cens=data_cens,xt=xt,xc=xc,eps=eps)
      }else{
        eps = eps*eta1
        p.old = p.new
        continue = (tau>tol) & (iter<maxit)
        iter = iter+1
        if(iter==maxit){
          p.new=c(0,0,0,0,0)
        }
      }#end of else
    }#end of else
}#end of while

result = p.new

return(result)

}# end of the main function



################# NCG Functions ########################



log.lik.fun = function(par=c(b10,b11,g1,g2,alp),data_obs,data_cens,xt,xc){ 
  
  Fc=1-exp(-((par[4]*data_cens)^(1/par[3])))
  
  phic=exp(par[1]+par[2]*xc)/(1+(par[5]*exp(par[1]+par[2]*xc)))
  Ft=1-exp(-((par[4]*data_obs)^(1/par[3])))
  phit=exp(par[1]+par[2]*xt)/(1+(par[5]*exp(par[1]+par[2]*xt)))
  
  ft=(1/(data_obs*par[3]))*((par[4]*data_obs)^(1/par[3]))*(1-Ft)
  
  Spc=(1-(par[5]*phic*Fc))^(1/par[5])
  Spt=(1-(par[5]*phit*Ft))^(1/par[5])
  fpt=Spt*phit*ft*((1-(par[5]*phit*Ft))^(-1))
  
  log.lik=sum(log(fpt))+sum(log(Spc))
  
  return(-log.lik)
  
}


lambda = function(pars,d.k,g.k,data_obs,data_cens,xt,xc,del){
  
  k = 1
  cont = TRUE
  
  while(cont){
    lam.k = 1/(2^(k-1))
    pp = pars+(d.k*lam.k)
    pp.new = c(pp[1],pp[2],max(pp[3],0.01),max(pp[4],0.01),min(max(pp[5],0.01),1))
    cc1 = log.lik.fun(pp.new,data_obs=data_obs,data_cens=data_cens,xt=xt,xc=xc)
    cc2 = (log.lik.fun(pars,data_obs=data_obs,data_cens=data_cens,xt=xt,xc=xc)+(del*lam.k*sum(d.k*g.k))) 
    
    #cat("s = ", k ,",", "cc1 = ", cc1,",", "cc2 = ", cc2, "\n")
    
    if(  (cc1 != "NaN" & cc1 != Inf & cc2 != "NaN" & cc2 != "Inf") ){
      cont = (cc1 > cc2) & (k < 17) 
      k = k + 1
    }else{
      k = k + 1
      if(k>16){
        cont = FALSE
      }
    }
    
  }#end of while
  return(c(lam.k,k))
}

#NCG algorithm for a fixed value of alpha in (0,1]

NCG_BC_Gen_Wei=function(data_obs,data_cens,xt,xc,tol,maxit,b0,b1,gam1,gam2,alpha,del){
  
  p.new=rep(0,5)
  p.old=rep(0,5)
  d.old=rep(0,5)
  d.new=rep(0,5)
  
  p.old = c(b0,b1,gam1,gam2,alpha)
  d.old = -1*grad(log.lik.fun,p.old,method="Richardson",data_obs=data_obs,data_cens=data_cens,xt=xt,xc=xc)
  
  
  continue = TRUE
  iter=1
  
  while(continue){
    
    #print(iter)
    
    g.old = grad(log.lik.fun,p.old,method="Richardson",data_obs=data_obs,data_cens=data_cens,xt=xt,xc=xc)
    lam.vec =  lambda(p.old,d.old,g.old,data_obs=data_obs,data_cens=data_cens,xt=xt,xc=xc,del=del) # lam calculated from line search algorithm
    
    if(lam.vec[2]>16){
      p.new = p.old
      
      #cat("lam = ", lam.vec[1],",", "k = ", lam.vec[2], ",", "pnew = ", p.new, ",","norm = ",sqrt(sum(g.old*g.old)), "," , "value of fn = ", 
      #log.lik.fun(p.old,data_obs=data_obs,data_cens=data_cens,xt=xt,xc=xc), "\n")
      
      continue = FALSE
    }else{
      
      p.new = p.old + (lam.vec[1]*d.old)
      dum1.new = (-1*grad(log.lik.fun,p.new,method="Richardson",data_obs=data_obs,data_cens=data_cens,xt=xt,xc=xc))
      g.new = -dum1.new
      
      #cat("lam = ", lam.vec[1], ",", "k = ", lam.vec[2], ",", "pnew = ", p.new, ",", "norm = ", sqrt(sum(g.new*g.new)), ",", "value of fn = ", 
      #log.lik.fun(p.new,data_obs=data_obs,data_cens=data_cens,xt=xt,xc=xc), "\n")
      
      y =  g.new - g.old
      dum2 = y - (2*d.old*sum(y*y)/sum(d.old*y))
      HG = sum(dum2*g.new)/sum(d.old*y)
      d.new = dum1.new + (HG*d.old)
      
      iter = iter + 1
      
      continue=(abs((p.new[1]-p.old[1])/p.old[1])>tol | abs((p.new[2]-p.old[2])/p.old[2])>tol | abs((p.new[3]-p.old[3])/p.old[3])>tol | 
                  abs((p.new[4]-p.old[4])/p.old[4])>tol | abs((p.new[5]-p.old[5])/p.old[5])>tol) & (iter<maxit)
      
      p.old = p.new
      d.old = d.new
      #g.old = g.new
      
      if(iter==maxit){
        p.new=matrix(c(0,0,0,0,0))
      }
      
    }# end of else
    
  }#end of while
  
  
  result = p.new
  
  return(result)
  
}# end of the main function


# calling the functions

n = 500

m=300
p01_true=0.05
p02_true=0.65
g1_true=0.316
g2_true=0.179
alpha_true=0.5

tol=0.001
maxit=500
incr=0.20
c=0.10
del = 0.1
x_min=0.1
x_max=20

b1_true= (log((1/alpha_true)*(((1/p01_true)^alpha_true)-1))-log((1/alpha_true)*(((1/p02_true)^alpha_true)-1)))/(x_max-x_min)
b0_true= log((1/alpha_true)*(((1/p01_true)^alpha_true)-1))-(x_max*b1_true)


beta0_hat_sqh=rep(NA,n)
temp_b0_sqh=rep(NA,n)
beta1_hat_sqh=rep(NA,n)
temp_b1_sqh=rep(NA,n)
g1_hat_sqh=rep(NA,n)
temp_g1_sqh=rep(NA,n)
g2_hat_sqh=rep(NA,n)
temp_g2_sqh=rep(NA,n)
alpha_hat_sqh=rep(NA,n)
temp_alpha_sqh=rep(NA,n)


beta0_hat_ncg=rep(NA,n)
temp_b0_ncg=rep(NA,n)
beta1_hat_ncg=rep(NA,n)
temp_b1_ncg=rep(NA,n)
g1_hat_ncg=rep(NA,n)
temp_g1_ncg=rep(NA,n)
g2_hat_ncg=rep(NA,n)
temp_g2_ncg=rep(NA,n)
alpha_hat_ncg=rep(NA,n)
temp_alpha_ncg=rep(NA,n)


count_ncg=0
count_sqh=0

for(i in 1:n){
  print(i)
  set.seed(100+i)
  
  data=data_gen_cont_Wei(m,g1_true,g2_true,p01_true,p02_true,c,x_min,x_max,alpha_true) # m,g1,g2,p01,p02,c,x_min,x_max,alpha
  
  obs_data=data[data$d==1,]
  cens_data=data[data$d==0,]
  data_obs=obs_data$t
  data_cens=cens_data$t
  xt=obs_data$x
  xc=cens_data$x
  
  b0_init=sample(seq((b0_true-(incr*abs(b0_true))),(b0_true+(incr*abs(b0_true))),by=0.01),1)
  b1_init=sample(seq((b1_true-(incr*abs(b1_true))),(b1_true+(incr*abs(b1_true))),by=0.01),1)
  g1_init=sample(seq((g1_true-(incr*abs(g1_true))),(g1_true+(incr*abs(g1_true))),by=0.01),1)
  g2_init=sample(seq((g2_true-(incr*abs(g2_true))),(g2_true+(incr*abs(g2_true))),by=0.01),1)
  alpha_init=sample(seq((alpha_true-(incr*abs(alpha_true))),(alpha_true+(incr*abs(alpha_true))),by=0.01),1)
  
  nr_sqh = SQH_BC_Gen_Wei(data_obs=data_obs,data_cens=data_cens,xt=xt,xc=xc,tol=tol,maxit=1000,b10=b0_init,b11=b1_init,g1=g1_init,g2=g2_init,alp=alpha_init,
                          epsilon1=1000,lambda1=1000,eta1=0.5,rho1=1000)
  
  nr_ncg = NCG_BC_Gen_Wei(data_obs,data_cens,xt,xc,0.001,500,b0_init,b1_init,g1_init,g2_init,alpha_init,del)

  if(0%in%nr_sqh){
    
    count_sqh=count_sqh+1
    
    beta0_hat_sqh[i]=0
    temp_b0_sqh[i]=0
    beta1_hat_sqh[i]=0
    temp_b1_sqh[i]=0
    g1_hat_sqh[i]=0
    temp_g1_sqh[i]=0
    g2_hat_sqh[i]=0
    temp_g2_sqh[i]=0
    alpha_hat_sqh[i]=0
    temp_alpha_sqh[i]=0
    
  }else{
    
    beta0_hat_sqh[i]=nr_sqh[1]
    temp_b0_sqh[i]=abs(beta0_hat_sqh[i]-b0_true)
    beta1_hat_sqh[i]=nr_sqh[2]
    temp_b1_sqh[i]=abs(beta1_hat_sqh[i]-b1_true)
    g1_hat_sqh[i]=nr_sqh[3]
    temp_g1_sqh[i]=abs(g1_hat_sqh[i]-g1_true)
    g2_hat_sqh[i]=nr_sqh[4]
    temp_g2_sqh[i]=abs(g2_hat_sqh[i]-g2_true)
    alpha_hat_sqh[i]=nr_sqh[5]
    temp_alpha_sqh[i]=abs(alpha_hat_sqh[i]-alpha_true)
    
  }#end of else
  
  
  if(0%in%nr_ncg){
    
    count_ncg=count_ncg+1
    
    beta0_hat_ncg[i]=0
    temp_b0_ncg[i]=0
    beta1_hat_ncg[i]=0
    temp_b1_ncg[i]=0
    g1_hat_ncg[i]=0
    temp_g1_ncg[i]=0
    g2_hat_ncg[i]=0
    temp_g2_ncg[i]=0
    alpha_hat_ncg[i]=0
    temp_alpha_ncg[i]=0
    
  }else{
    
    beta0_hat_ncg[i]=nr_ncg[1]
    temp_b0_ncg[i]=abs(beta0_hat_ncg[i]-b0_true)
    beta1_hat_ncg[i]=nr_ncg[2]
    temp_b1_ncg[i]=abs(beta1_hat_ncg[i]-b1_true)
    g1_hat_ncg[i]=nr_ncg[3]
    temp_g1_ncg[i]=abs(g1_hat_ncg[i]-g1_true)
    g2_hat_ncg[i]=nr_ncg[4]
    temp_g2_ncg[i]=abs(g2_hat_ncg[i]-g2_true)
    alpha_hat_ncg[i]=nr_ncg[5]
    temp_alpha_ncg[i]=abs(alpha_hat_ncg[i]-alpha_true)
    
  }#end of else
  
  
}#end of for


avg_b0_sqh=sum(beta0_hat_sqh)/(n-count_sqh)
avg_b1_sqh=sum(beta1_hat_sqh)/(n-count_sqh)
avg_g1_sqh=sum(g1_hat_sqh)/(n-count_sqh)
avg_g2_sqh=sum(g2_hat_sqh)/(n-count_sqh)
avg_alpha_sqh = sum(alpha_hat_sqh)/(n-count_sqh)


bias_b0_sqh=sum(temp_b0_sqh)/(n-count_sqh)
bias_b1_sqh=sum(temp_b1_sqh)/(n-count_sqh)
bias_g1_sqh=sum(temp_g1_sqh)/(n-count_sqh)
bias_g2_sqh=sum(temp_g2_sqh)/(n-count_sqh)
bias_alpha_sqh=sum(temp_alpha_sqh)/(n-count_sqh)


mse_b0_sqh=sum(temp_b0_sqh^2)/(n-count_sqh-1)
mse_b1_sqh=sum(temp_b1_sqh^2)/(n-count_sqh-1)
mse_g1_sqh=sum(temp_g1_sqh^2)/(n-count_sqh-1)
mse_g2_sqh=sum(temp_g2_sqh^2)/(n-count_sqh-1)
mse_alpha_sqh=sum(temp_alpha_sqh^2)/(n-count_sqh-1)


rmse_b0_sqh=sqrt(mse_b0_sqh)
rmse_b1_sqh=sqrt(mse_b1_sqh)
rmse_g1_sqh=sqrt(mse_g1_sqh)
rmse_g2_sqh=sqrt(mse_g2_sqh)
rmse_alpha_sqh=sqrt(mse_alpha_sqh)



avg_b0_sqh
avg_b1_sqh
avg_g1_sqh
avg_g2_sqh
avg_alpha_sqh

bias_b0_sqh
bias_b1_sqh
bias_g1_sqh
bias_g2_sqh
bias_alpha_sqh

rmse_b0_sqh
rmse_b1_sqh
rmse_g1_sqh
rmse_g2_sqh
rmse_alpha_sqh

#### NCG #######

avg_b0_ncg=sum(beta0_hat_ncg)/(n-count_ncg)
avg_b1_ncg=sum(beta1_hat_ncg)/(n-count_ncg)
avg_g1_ncg=sum(g1_hat_ncg)/(n-count_ncg)
avg_g2_ncg=sum(g2_hat_ncg)/(n-count_ncg)
avg_alpha_ncg = sum(alpha_hat_ncg)/(n-count_ncg)


bias_b0_ncg=sum(temp_b0_ncg)/(n-count_ncg)
bias_b1_ncg=sum(temp_b1_ncg)/(n-count_ncg)
bias_g1_ncg=sum(temp_g1_ncg)/(n-count_ncg)
bias_g2_ncg=sum(temp_g2_ncg)/(n-count_ncg)
bias_alpha_ncg=sum(temp_alpha_ncg)/(n-count_ncg)


mse_b0_ncg=sum(temp_b0_ncg^2)/(n-count_ncg-1)
mse_b1_ncg=sum(temp_b1_ncg^2)/(n-count_ncg-1)
mse_g1_ncg=sum(temp_g1_ncg^2)/(n-count_ncg-1)
mse_g2_ncg=sum(temp_g2_ncg^2)/(n-count_ncg-1)
mse_alpha_ncg=sum(temp_alpha_ncg^2)/(n-count_ncg-1)


rmse_b0_ncg=sqrt(mse_b0_ncg)
rmse_b1_ncg=sqrt(mse_b1_ncg)
rmse_g1_ncg=sqrt(mse_g1_ncg)
rmse_g2_ncg=sqrt(mse_g2_ncg)
rmse_alpha_ncg=sqrt(mse_alpha_ncg)



avg_b0_ncg
avg_b1_ncg
avg_g1_ncg
avg_g2_ncg
avg_alpha_ncg

bias_b0_ncg
bias_b1_ncg
bias_g1_ncg
bias_g2_ncg
bias_alpha_ncg

rmse_b0_ncg
rmse_b1_ncg
rmse_g1_ncg
rmse_g2_ncg
rmse_alpha_ncg


bias_sqh = c(bias_b0_sqh,bias_b1_sqh,bias_g1_sqh,bias_g2_sqh,bias_alpha_sqh)
bias_ncg = c(bias_b0_ncg,bias_b1_ncg,bias_g1_ncg,bias_g2_ncg,bias_alpha_ncg)
rmse_sqh = c(rmse_b0_sqh,rmse_b1_sqh,rmse_g1_sqh,rmse_g2_sqh,rmse_alpha_sqh)
rmse_ncg = c(rmse_b0_ncg,rmse_b1_ncg,rmse_g1_ncg,rmse_g2_ncg,rmse_alpha_ncg)

round(data.frame(bias_sqh,bias_ncg,rmse_sqh,rmse_ncg),3)

print("count_sqh is:")
print(count_sqh)

print("count_ncg is:")
print(count_ncg)




