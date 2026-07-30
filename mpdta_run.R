###############################################################################
## mpdta_run.R -- lead empirical application: Callaway & Sant'Anna (2021)
## county minimum-wage / teen-employment data (mpdta, public in the `did` pkg).
##
## We estimate the K=7 post-treatment cohort-time effects ATT(g,t) with the
## Callaway-Sant'Anna estimator, take their EXACT joint covariance from the
## influence functions (cells share counties, hence are correlated), and apply
## the l0 (BIC) and DP estimators using that covariance -- the same exact-
## covariance treatment validated in the simulations.
###############################################################################
dir.create("figures", showWarnings = FALSE)
suppressMessages(library(did)); set.seed(7)
data(mpdta, package="did")

res <- att_gt(yname="lemp", tname="year", idname="countyreal",
              gname="first.treat", control_group="notyettreated",
              data=mpdta, bstrap=FALSE, cband=FALSE)

keep <- res$t >= res$group                    # post-treatment cells
g <- res$group[keep]; ti <- res$t[keep]
tau <- res$att[keep]; se <- res$se[keep]
lab_cell <- paste0(g, ":", ti)
K <- length(tau)
n  <- res$n                                   # number of units
Sigma <- (as.matrix(res$V_analytical) / n)[keep, keep, drop=FALSE]  # exact K x K cov
cat(sprintf("K=%d cells, n=%d units. SE check (max|IF - did|)= %.2e\n",
            K, n, max(abs(sqrt(diag(Sigma)) - se))))
Om <- solve(Sigma)                            # precision matrix
cvec <- Om %*% tau; yty <- as.numeric(t(tau) %*% Om %*% tau)
prec_d <- 1/se^2                              # diagonal precision (path moves)

memb <- function(lab){u<-sort(unique(lab));A<-matrix(0,length(lab),length(u))
  for(j in seq_along(u))A[lab==u[j],j]<-1;A}
## GLS grouped fit using exact Omega
fit <- function(lab){A<-memb(lab);GA<-t(A)%*%Om%*%A;cA<-t(A)%*%cvec
  GAi<-solve(GA);phi<-GAi%*%cA;tg<-as.vector(A%*%phi)
  list(tau=tg,dev=yty-sum(cA*phi),m=ncol(A),cov=A%*%GAi%*%t(A))}

## ---- l0: greedy path (diagonal moves) + BIC on exact GLS deviance ----
groups<-as.list(seq_len(K));gt<-tau;gn<-prec_d;LA<-vector("list",K);LA[[K]]<-seq_len(K)
for(s in seq_len(K-1)){ng<-length(groups);best<-Inf;bi<-bj<-NA
 for(a in 1:(ng-1))for(b in (a+1):ng){c0<-(gn[a]*gn[b]/(gn[a]+gn[b]))*(gt[a]-gt[b])^2
  if(c0<best){best<-c0;bi<-a;bj<-b}}
 nn<-gn[bi]+gn[bj];tt<-(gn[bi]*gt[bi]+gn[bj]*gt[bj])/nn
 groups[[bi]]<-c(groups[[bi]],groups[[bj]]);gt[bi]<-tt;gn[bi]<-nn
 groups[[bj]]<-NULL;gt<-gt[-bj];gn<-gn[-bj]
 lb<-integer(K);for(gi in seq_along(groups))lb[groups[[gi]]]<-gi;LA[[K-s]]<-lb}
bic<-sapply(1:K,function(m){fit(LA[[m]])$dev + m*log(K)})
ms<-which.min(bic); l0<-fit(LA[[ms]]); l0lab<-LA[[ms]]

## ---- DP: diagonal moves, EXACT-covariance recorded draws ----
mu0<-median(tau); s0sq<-(10*sd(tau))^2; alpha<-1; iters<-4000; burn<-1000
z<-rep(1L,K); b0<-1/s0sq; post<-numeric(K); coc<-matrix(0,K,K); kept<-0
md<-numeric(iters-burn); attdr<-numeric(iters-burn)
Ncoh<-c("2004"=100,"2006"=200,"2007"=655)[as.character(g)]      # cohort sizes
w<-Ncoh/sum(Ncoh)                                              # ATT weights
for(it in seq_len(iters)){
 for(k in seq_len(K)){z[k]<-NA;cl<-sort(unique(z[!is.na(z)]));lp<-numeric(length(cl)+1)
  for(j in seq_along(cl)){me<-which(z==cl[j]);Pc<-b0+sum(prec_d[me])
   mc<-(mu0*b0+sum(prec_d[me]*tau[me]))/Pc
   lp[j]<-log(length(me))+dnorm(tau[k],mc,sqrt(1/prec_d[k]+1/Pc),log=TRUE)}
  lp[length(cl)+1]<-log(alpha)+dnorm(tau[k],mu0,sqrt(1/prec_d[k]+s0sq),log=TRUE)
  p<-exp(lp-max(lp));p<-p/sum(p);pk<-sample.int(length(p),1,prob=p)
  z[k]<-if(pk<=length(cl))cl[pk] else max(cl,0)+1}
 z<-match(z,sort(unique(z)))
 if(it>burn){kept<-kept+1;A<-memb(z);m<-ncol(A)
  GA<-t(A)%*%Om%*%A;Vs<-solve(GA+diag(m)/s0sq);mu<-Vs%*%(t(A)%*%cvec+mu0/s0sq)
  phi<-as.vector(mu+t(chol(Vs))%*%rnorm(m));tg<-as.vector(A%*%phi)
  post<-post+tg;coc<-coc+outer(z,z,`==`);attdr[kept]<-sum(w*tg);md[kept]<-m}}
dp_tau<-post/kept;coc<-coc/kept

## ---- efficiency: per-cell variance, flexible vs l0-grouped ----
var_flex<-diag(Sigma); var_l0<-diag(l0$cov)
att_flex<-sum(w*tau); v_attflex<-as.numeric(t(w)%*%Sigma%*%w)
att_l0<-sum(w*l0$tau); v_attl0<-as.numeric(t(w)%*%l0$cov%*%w)

cat("\n================ mpdta: cohort-time effects & grouping ================\n")
tab<-data.frame(cell=lab_cell, att=round(tau,4), se=round(se,4),
                l0_group=l0lab, l0_effect=round(l0$tau,4),
                dp_effect=round(dp_tau,4),
                var_ratio=round(var_l0/var_flex,2))
print(tab, row.names=FALSE)
cat(sprintf("\nSignal-to-noise sd(att)/median(se) = %.2f\n", sd(tau)/median(se)))
cat(sprintf("l0 (BIC): %d groups ; DP (alpha=1): %.2f groups (posterior mean)\n",
            ms, mean(md)))
cat(sprintf("Per-cell variance reduction (l0 vs flexible), pooled cells: mean ratio %.2f\n",
            mean((var_l0/var_flex)[table(l0lab)[as.character(l0lab)]>1])))
cat(sprintf("Overall ATT: flexible %.4f (se %.4f) ; l0 %.4f (se %.4f) ; DP %.4f [%.4f, %.4f]\n",
            att_flex, sqrt(v_attflex), att_l0, sqrt(v_attl0),
            mean(attdr), quantile(attdr,.025), quantile(attdr,.975)))
cat("\nDP co-clustering probabilities (rows/cols = cells):\n")
dimnames(coc)<-list(lab_cell,lab_cell); print(round(coc,2))

## ---- figures ----
png("figures/mpdta_effects.png",1700,950,res=200)
o<-order(tau);cols<-c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e")
plot(seq_len(K),tau[o],pch=19,cex=1.3,col=cols[l0lab[o]],
     ylim=range(tau-1.96*se,tau+1.96*se),xaxt="n",
     xlab="cohort:year cell (sorted by effect)",ylab="ATT(g,t)  (log teen employment)",
     main=sprintf("mpdta: cohort-time effects grouped by l0 (%d groups)",ms))
axis(1,at=seq_len(K),labels=lab_cell[o],cex.axis=.8)
arrows(seq_len(K),(tau-1.96*se)[o],seq_len(K),(tau+1.96*se)[o],length=0,col="grey70")
segments(seq_len(K)-.3,l0$tau[o],seq_len(K)+.3,l0$tau[o],col=cols[l0lab[o]],lwd=3)
abline(h=0,lty=3)
legend("bottomleft",paste("group",sort(unique(l0lab))),col=cols,pch=19,bty="n")
dev.off()

png("figures/mpdta_coclust.png",1200,1100,res=200)
o<-order(tau)
image(1:K,1:K,coc[o,o],axes=FALSE,col=grey.colors(20,start=1,end=0),
      xlab="",ylab="",main="DP posterior co-clustering probability")
axis(1,1:K,lab_cell[o],las=2,cex.axis=.8);axis(2,1:K,lab_cell[o],las=1,cex.axis=.8)
dev.off()

## ---- dynamic (event-time) profile: explains the recovered structure ----
dyn <- aggte(res, type="dynamic", na.rm=TRUE)
cat("\nDynamic (event-time) aggregation:\n")
print(data.frame(e=dyn$egt, att=round(dyn$att.egt,4), se=round(dyn$se.egt,4)))
png("figures/mpdta_dynamic.png",1500,900,res=200)
e<-dyn$egt; a<-dyn$att.egt; s<-dyn$se.egt
plot(e,a,pch=19,type="b",ylim=range(a-1.96*s,a+1.96*s),
     xlab="event time (years since minimum-wage increase)",
     ylab="ATT (log teen employment)",main="mpdta: dynamic effect profile")
arrows(e,a-1.96*s,e,a+1.96*s,length=0,col="grey60")
abline(h=0,lty=3); abline(v=-0.5,lty=2,col="grey50")
dev.off()

write.csv(tab,"mpdta_grouped.csv",row.names=FALSE)
cat("\nWrote mpdta_grouped.csv, mpdta_effects.png, mpdta_coclust.png, mpdta_dynamic.png\n")
