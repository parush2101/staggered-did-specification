###############################################################################
## cengiz_run.R  -- apply l0 (BIC) and DP to the 138 event-level employment
## effects from Cengiz et al. (2019), extracted from the Harvard Dataverse
## replication package (DOI 10.7910/DVN/TJCTC7).
##
## Per-event effect  tau_g = mean(p1..p4)  (avg post-treatment effect on affected
## employment, event-study coefficients).  Per-event noise is gauged INTERNALLY
## from the pre-treatment coefficients m1..m4, whose true value is zero under
## no-anticipation, so their cross-event dispersion is pure sampling noise.
###############################################################################
dir.create("figures", showWarnings = FALSE)
suppressMessages(library(haven))
set.seed(42)
d <- as.data.frame(read_dta("cengiz_repl/event_fulltable.dta"))
pre <- c("m1","m2","m3","m4"); post <- c("p1","p2","p3","p4")

tau  <- rowMeans(d[, post])
se_c <- sd(rowMeans(d[, pre]))          # common per-event SE from pre-trends
E    <- length(tau)
prec <- rep(1/se_c^2, E)
wt   <- d$population / sum(d$population) # population weights for the overall ATT

att_pool <- sum(prec*tau)/sum(prec)     # single common effect (pooled)
att_flex <- sum(wt*tau)                 # weighted avg of flexible CATTs

## ---- helpers (summary-statistic form; events are independent blocks) -------
grp_fit <- function(lab){ phi <- tapply(prec*tau,lab,sum)/tapply(prec,lab,sum)
  tg <- phi[as.character(lab)]; list(tau=as.numeric(tg),
    rss=sum(prec*(tau-tg)^2), m=length(phi)) }

l0_run <- function(){
  groups <- as.list(seq_len(E)); gt <- tau; gn <- prec
  labs_at <- vector("list",E); labs_at[[E]] <- seq_len(E)
  for(s in seq_len(E-1)){ ng<-length(groups); best<-Inf; bi<-bj<-NA
    for(a in 1:(ng-1)) for(b in (a+1):ng){
      c0<-(gn[a]*gn[b]/(gn[a]+gn[b]))*(gt[a]-gt[b])^2
      if(c0<best){best<-c0;bi<-a;bj<-b} }
    nn<-gn[bi]+gn[bj]; tt<-(gn[bi]*gt[bi]+gn[bj]*gt[bj])/nn
    groups[[bi]]<-c(groups[[bi]],groups[[bj]]); gt[bi]<-tt; gn[bi]<-nn
    groups[[bj]]<-NULL; gt<-gt[-bj]; gn<-gn[-bj]
    lab<-integer(E); for(gi in seq_along(groups)) lab[groups[[gi]]]<-gi
    labs_at[[E-s]]<-lab }
  bic<-sapply(1:E,function(m){f<-grp_fit(labs_at[[m]]); f$rss+m*log(E)})
  ms<-which.min(bic); list(labels=labs_at[[ms]],m=ms,fit=grp_fit(labs_at[[ms]]),bic=bic)
}

dp_run <- function(alpha, mu0=median(tau), s0sq=(4*sd(tau))^2,
                   iters=3000, burn=1000){
  z<-rep(1L,E); b0<-1/s0sq; post_tau<-numeric(E); coc<-matrix(0,E,E); kept<-0
  att<-numeric(iters-burn); md<-numeric(iters-burn)
  for(it in seq_len(iters)){
    for(k in seq_len(E)){ z[k]<-NA; cl<-sort(unique(z[!is.na(z)])); lp<-numeric(length(cl)+1)
      for(j in seq_along(cl)){ me<-which(z==cl[j]); Pc<-b0+sum(prec[me])
        mc<-(mu0*b0+sum(prec[me]*tau[me]))/Pc
        lp[j]<-log(length(me))+dnorm(tau[k],mc,sqrt(1/prec[k]+1/Pc),log=TRUE) }
      lp[length(cl)+1]<-log(alpha)+dnorm(tau[k],mu0,sqrt(1/prec[k]+s0sq),log=TRUE)
      p<-exp(lp-max(lp)); p<-p/sum(p); pk<-sample.int(length(p),1,prob=p)
      z[k]<-if(pk<=length(cl)) cl[pk] else max(cl,0)+1 }
    z<-match(z,sort(unique(z)))
    if(it>burn){ kept<-kept+1; phi<-numeric(max(z))
      for(j in seq_len(max(z))){ me<-which(z==j); Pc<-b0+sum(prec[me])
        mc<-(mu0*b0+sum(prec[me]*tau[me]))/Pc; phi[j]<-rnorm(1,mc,sqrt(1/Pc)) }
      post_tau<-post_tau+phi[z]; coc<-coc+outer(z,z,`==`)
      att[kept]<-sum(wt*phi[z]); md[kept]<-max(z) } }
  list(tau=post_tau/kept, coclust=coc/kept, labels=z, m_mean=mean(md),
       att=mean(att), att_lo=quantile(att,.025), att_hi=quantile(att,.975))
}

l0 <- l0_run()
dp1 <- dp_run(alpha = 1)     # parsimonious prior
dp5 <- dp_run(alpha = 5)     # weaker prior (robustness)

cat("================ Cengiz et al. (2019): specification ================\n")
cat(sprintf("Events: %d\n", E))
cat(sprintf("Per-event noise (pre-trend sd): %.5f ; post-effect spread: %.5f ; S/N: %.2f\n",
            se_c, sd(tau), sd(tau)/se_c))
cat(sprintf("Pooled (common) effect          : %.5f\n", att_pool))
cat(sprintf("Flexible ATT (pop-wtd avg CATT) : %.5f\n", att_flex))
cat(sprintf("l0 (BIC)   : %d group(s); ATT = %.5f\n", l0$m, sum(wt*l0$fit$tau)))
cat(sprintf("DP (alpha=1): %.2f groups; ATT = %.5f  [%.5f, %.5f]\n",
            dp1$m_mean, dp1$att, dp1$att_lo, dp1$att_hi))
cat(sprintf("DP (alpha=5): %.2f groups; ATT = %.5f  [%.5f, %.5f]\n",
            dp5$m_mean, dp5$att, dp5$att_lo, dp5$att_hi))

## flexible-CATT variance vs l0/DP (efficiency): report the SE of the ATT
se_flex <- sqrt(sum(wt^2)) * se_c        # equal common noise
cat(sprintf("ATT SE: flexible ~ %.5f ; pooled/l0 ~ %.5f (1 group)\n",
            se_flex, se_c*sqrt(sum(wt^2 * (l0$m==1)))))

res <- data.frame(event=d$event_no, statenum=d$statenum, tau=tau, se=se_c,
                  l0_group=l0$labels, dp_effect=dp1$tau)
write.csv(res, "cengiz_events_grouped.csv", row.names=FALSE)

## figure: sorted event effects with common CI and the pooled line
png("figures/cengiz_effects.png", 1700, 900, res=200)
o <- order(tau)
plot(seq_len(E), tau[o], pch=19, col="grey40", ylim=range(tau-1.96*se_c, tau+1.96*se_c),
     xlab="event (sorted by effect)", ylab="avg post effect on affected employment",
     main=sprintf("Cengiz events: flexible CATTs vs pooled (l0 selects %d group)", l0$m))
arrows(seq_len(E),(tau-1.96*se_c)[o],seq_len(E),(tau+1.96*se_c)[o],length=0,col="grey80")
abline(h=att_pool, col="firebrick", lwd=2, lty=2)
legend("topleft", c("flexible CATT (±1.96 SE)","pooled effect"),
       pch=c(19,NA), lty=c(NA,2), col=c("grey40","firebrick"), bty="n")
dev.off()
cat("\nWrote cengiz_events_grouped.csv, cengiz_effects.png\n")
