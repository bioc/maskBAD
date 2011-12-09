##library(hgu95av2cdf)
mask <- function(affy,exprlist=NULL,useExpr=TRUE,ind,PM=FALSE,verbose=TRUE){
  v1 <- which(ind==1)
  v2 <- which(ind==2)
  .method1 <- function(x,v1,v2){
    a <- dim(x)
    change <- matrix(1,ncol=a[1],nrow=a[1])
    for(i in 1:a[1]){
      for(j in 1:a[1]){
        if(i!=j){
          lfit <- .regression(x[i,],x[j,],reg="RMA")
          if(lfit$a>0){
            interc <- lfit$b
          }
          if(lfit$a<=0){
            interc <- min(x[i,])
          }
          sum.c <- (x[i,v1]-interc)/x[j,v1]
          sum.h <- (x[i,v2]-interc)/x[j,v2]
          change[i,j] <- t.test(sum.c,sum.h)$p.value
        }
      }
    }
    dummy <- .sort.matrix(change)
    sort.p.value <- dummy$p.value
    names(sort.p.value) <- dummy$sum.order
    return(list(change=change,sort.p.value=sort.p.value))
  }
  .sort.matrix <- function(matrix,SUM=FALSE){
    m <- matrix
    a <- dim(matrix)
    sum.order <- rep(0,a[1])
    tt <- seq(1,a[1])
    p.value <- rep(0,a[1])
    if(SUM==FALSE){
      m <- log(m)
    }
    for(j in seq(1,a[1])){
      sum.sum <- colSums(m)+rowSums(m)
      if(j<a[1]){
        p.value[j] <- min(sum.sum)/(2*a[1]-2*j)
        max.p <- max(sum.sum)/(2*a[1]-2*j)
      }
      if(j==a[1]){
        p.value[j] <- max.p
      }
      order.sum <- order(sum.sum)
      vv <- tt[order.sum[1]]
      tt <- tt[tt!=tt[order.sum[1]]]
      uu <- seq(1,(a[1]-j+1))
      uu <- uu[uu!=order.sum[1]]
      m <- matrix(m[uu,uu],a[1]-j,a[1]-j)
      rownames(m) <- tt
      colnames(m) <- tt
      sum.order[j] <- vv
    }
    if(SUM==FALSE){
      p.value <- exp(p.value)
    }
    names(p.value) <- sum.order
    
    return(list(sum.order=sum.order,p.value=p.value))
  }
  .regression <- function(x,y,reg="RMA",siex=1,siey=1){
    mx <- mean(x)
    my <- mean(y)
    sxsq <- sum((x-mx)*(x-mx))
    sxy <- sum((x-mx)*(y-my))
    sysq <- sum((y-my)*(y-my))
    n <- length(x)
    if(reg=="RMA"){
      beta1 <- sxy/sxsq
      beta2 <- sysq/sxy
      if(sxy!=0){
        a <- sign(sxy)*sqrt(beta1*beta2)
        b <- my-a*mx
      }
      if(sxy==0){
        lm <- lm(y~x)
        a <- lm$coefficients[2]
        b <- lm$coefficients[1]
      }
    }
    return(list(a=a,b=b,n=n,sxsq=sxsq,sysq=sysq,sxy=sxy))
  }
  if(verbose) cat("extract expression information\n")
  g <- geneNames(affy)
  if(length(exprlist)>0){
    g <- g[g %in% exprlist]
  }
  if(useExpr){
    aa <- mas5calls(affy,verbose=FALSE)
    expr.probes <- exprs(aa)
    exprval <- 0.9
    .P <- function(x){
      pp1 <- sum(x[v1]=="P")
      pp2 <- sum(x[v2]=="P")
      PP <- min(pp1/length(v1),pp2/length(v2))
      return(PP)
    }
    PP <- apply(expr.probes,1,FUN=.P)
    list <- which(PP>=exprval)
    ## if(PM) mm <- mm[g %in% names(list)]
    ## l <- l[g %in% names(list)]
    g <- g[g %in% names(list)]
  }
  if(length(g)<50000){
    l <- lapply(g,function(gene) pm(affy,gene))
    if(PM) mm <- lapply(g,function(gene) mm(affy,gene))
  }
  if(length(g)>=50000){
    if(verbose) cat("huge amount of probesets - extraction will take some time...\n")
    for(subs in 1:(trunc(length(g)/5000+1))){
      cat(signif(subs/(trunc(length(g)/5000+1)),3)*100,"%\n")
      if(subs==1){
        l <- lapply(g[(5000*(subs-1)+1):(5000*subs)],function(gene) pm(affy,gene))
        if(PM) mm <- lapply(g[(5000*(subs-1)+1):(5000*subs)],function(gene) mm(affy,gene))
      }
      if(!(subs %in% c(1,(trunc(length(g)/5000+1))))){
        l <- c(l,lapply(g[(5000*(subs-1)+1):(5000*subs)],function(gene) pm(affy,gene)))
        if(PM) mm <- c(mm,lapply(g[(5000*(subs-1)+1):(5000*subs)],function(gene) mm(affy,gene)))
      }
      if(subs==(trunc(length(g)+1))){
        l <- c(l,lapply(g[(5000*(subs-1)+1):length(g)],function(gene) pm(affy,gene)))
        if(PM) mm <- c(mm,lapply(g[(5000*(subs-1)+1):length(g)],function(gene) mm(affy,gene)))
      }
    }
  }
  indices <- indexProbes(affy,"pm",unique(probeNames(affy)))
  locations <- indices2xy(unlist(indices),abatch=affy)
  rownames <- rownames(locations)
  PR <- 0
  PRxy <- matrix(0,ncol=2,nrow=1)
  GENE.id <- ""
  pm1 <- ""
  if(PM) not.used.probes <- ""
  lel <- length(l)
  if(verbose) cat("mask probes\n")
  for(k in 1:lel){
    if(length(l[[k]])!=0){
      if(verbose) cat(k,"/",lel,"\n")
      if(PM){
        diff <- l[[k]]-mm[[k]]
        pm.mm <- apply(diff,1,mean)
        pm.m <- apply(l[[k]],1,mean)
        mm.m <- apply(mm[[k]],1,mean)
        disc.pr <- names(pm.mm)[which(sign(pm.mm)==-1)]
        if(length(dim(l[[k]][-which(rownames(l[[k]]) %in% disc.pr),]))!=0){
          if(length(intersect(rownames(l[[k]]),disc.pr))!=0){
            l[[k]] <- l[[k]][-which(rownames(l[[k]]) %in% disc.pr),]
            mm[[k]] <- mm[[k]][-which(rownames(mm[[k]]) %in% disc.pr),]
          }
        }
      }
      if(nrow(l[[k]])==1){
        sort.p.value=c("1"=1)
      }
      else{
        sort.p.value <- .method1(l[[k]],v1,v2)$sort.p.value
        names(sort.p.value) <- as.numeric(gsub(g[k],"",rownames(l[[k]])[as.numeric(names(sort.p.value))]))
      }

      pr <- sort.p.value[order(as.numeric(names(sort.p.value)))]
      prxy <- matrix(0,ncol=2,nrow=length(sort.p.value))
      for(i in 1:length(sort.p.value)){
        if(nrow(l[[k]])>1) name <- paste(g[k],sort(as.numeric(names(sort.p.value)))[i],sep="") else name=g[k]
        prxy[i,] <- locations[which(rownames==name),,drop=FALSE]
      }
      gene.id <- rep(g[k],length(sort.p.value))
      PR <- c(PR,pr)
      PRxy <- rbind(PRxy,prxy)
      GENE.id <- c(GENE.id,gene.id)
      if(PM) pm1 <- c(pm1,disc.pr)
    }
  }
  probes <- data.frame(x=PRxy[-1,1],y=PRxy[-1,2],quality.score=PR[-1],probeset=GENE.id[-1])
  probes[,4] <- as.character(probes[,4])
  if(PM){
    return(list(probes=probes,notUsed=unlist(pm1[-1])))
  }
  if(!PM){
    return(list(probes=probes))
  }
}

overlapExprExtMasks <- function(probes,seqdata,cutoffs="none",wilcox.ks=FALSE,sample=10,plotCutoffs=TRUE,verbose=TRUE){
  .konfcp <- function(j,n,alpha){
    pd <- rep(0,length(j))
    pu <- rep(0,length(j))
    po <- rep(0,length(j))
    pun <- rep(0,length(j))
    pob <- rep(1,length(j))
    for(i in c(1:length(j))){
      pd[i] <- j[i]/n
      pu[i] <- 0
      if(j[i]>0){
        pu[i] <- qbeta((alpha/2),j[i],(n-j[i]+1))
      }
      po[i] <- 1
      if(j[i]<n){
        po[i] <- qbeta((1-(alpha/2)),(j[i]+1),(n-j[i]))
      }
      ##      if(n>100){
      ##        if(pd[i]>(qnorm((1-(alpha/2)))*sqrt(pd[i]*(1-pd[i])/n))){
      ##          pun[i] <- pd[i]-qnorm((1-(alpha/2)))*sqrt(pd[i]*(1-pd[i])/n)
      ##        }
      ##        if(pd[i]<(1-qnorm((1-(alpha/2)))*sqrt(pd[i]*(1-pd[i])/n))){
      ##          pob[i] <- pd[i]+qnorm((1-(alpha/2)))*sqrt(pd[i]*(1-pd[i])/n)
      ##    }
      ##  }
    }
    ##    return(list(pdach=pd,poben=po,punten=pu,pobenn=pob,puntenn=pun))
    return(c(pdach=pd,poben=po,punten=pu))
  }
  if(length(cutoffs)==1){
    cutoffs <- quantile(probes[,3],seq(0,1,0.01))
  }
  nr <- 1
  wiltest <- c()
  kstest <- c()
  wiltest1 <- c()
  kstest1 <- c()
  type1 <- rep(0,length(cutoffs))
  type2 <- rep(0,length(cutoffs))
  used.cutoff <- c()
  used.cutoff1 <- c()
  conf.t1 <- matrix(0,ncol=2,nrow=length(cutoffs))
  conf.t2 <- matrix(0,ncol=2,nrow=length(cutoffs))
  wilcoxo <- c()
  kso <- c()
  UC <- c()
  for(cut in cutoffs){
    if(verbose) cat(signif((which(cutoffs==cut)/length(cutoffs))*100,2),"% \r")
    if(wilcox.ks==TRUE){
      wilcox <- c()
      ks <- c()
      psd <- probes[paste(probes[,1],".",probes[,2],sep="") %in% paste(seqdata[seqdata[,3]==0,1],".",seqdata[seqdata[,3]==0,2],sep=""),3]
      pss <- probes[paste(probes[,1],".",probes[,2],sep="") %in% paste(seqdata[seqdata[,3]==1,1],".",seqdata[seqdata[,3]==1,2],sep=""),3]
      if(sum(pss>=cut,na.rm=TRUE)>0 & sum(psd>=cut,na.rm=TRUE)>0) wilcoxo[nr] <- suppressWarnings(wilcox.test(pss[pss>=cut],psd[psd>=cut])$p.value) else wilcoxo[nr] <- NA
      if(sum(pss>=cut,na.rm=TRUE)>0 & sum(psd>=cut,na.rm=TRUE)>0) kso[nr] <- suppressWarnings(ks.test(pss[pss>=cut],psd[psd>=cut])$p.value) else kso[nr] <- NA
      if(sum(pss>=cut,na.rm=TRUE)>0 & sum(psd>=cut,na.rm=TRUE)>0) UC[nr] <- cut else UC[nr] <- NA
      for(j in 1:sample){
        if(sum(pss>=cut & !is.na(pss))>0 & sum(psd>=cut & !is.na(psd))>0){
          wilcox[j] <- suppressWarnings(wilcox.test(sample(pss[pss>=cut],sum(pss>=cut,na.rm=TRUE),replace=TRUE),
                                   sample(psd[psd>=cut],sum(psd>=cut,na.rm=TRUE),replace=TRUE))$p.value)
          ks[j] <- suppressWarnings(ks.test(sample(pss[pss>=cut],sum(pss>=cut,na.rm=TRUE),replace=TRUE),
                                            sample(psd[psd>=cut],sum(psd>=cut,na.rm=TRUE),replace=TRUE))$p.value)
          used.cutoff <- c(used.cutoff,cut)
        }
      }
      wiltest <- c(wiltest,wilcox)
      kstest <- c(kstest,ks)
    }
    ## if(wilcox.ks==TRUE){
    ##   wilcox1 <- c()
    ##   ks1 <- c()
    ##   psd <- probes[paste(probes[,1],".",probes[,2],sep="") %in% paste(seqdata[seqdata[,3]==0,1],".",seqdata[seqdata[,3]==0,2],sep=""),3]
    ##   pss <- probes[paste(probes[,1],".",probes[,2],sep="") %in% paste(seqdata[seqdata[,3]==1,1],".",seqdata[seqdata[,3]==1,2],sep=""),3]
    ##   for(j in 1:sample){
    ##     if(sum(pss<=cut & !is.na(pss))>0 & sum(psd<=cut & !is.na(psd))>0){
    ##       wilcox1[j] <- wilcox.test(sample(pss[pss<=cut],100,replace=TRUE),
    ##                                sample(psd[psd<=cut],100,replace=TRUE))$p.value
    ##       ks1[j] <- ks.test(sample(pss[pss<=cut],100,replace=TRUE),sample(psd[psd<=cut],100,replace=TRUE))$p.value
    ##       used.cutoff1 <- c(used.cutoff1,cut)
    ##     }
    ##   }
    ##   wiltest1 <- c(wiltest1,wilcox1)
    ##   kstest1 <- c(kstest1,ks1)
    ## }
    probes <- probes[paste(probes[,1],".",probes[,2],sep="") %in% paste(seqdata[,1],".",seqdata[,2],sep=""),]
    mp <- order(probes[,1],probes[,2])
    probes <- probes[mp,]
    seqdata <- seqdata[paste(seqdata[,1],".",seqdata[,2],sep="")  %in% paste(probes[,1],".",probes[,2],sep=""),]
    ms <- order(seqdata[,1],seqdata[,2])
    seqdata <- seqdata[ms,]
    t11 <- length(probes[probes[,3]>cut & seqdata[,3]==1,3])
    t10 <- length(probes[probes[,3]>cut & seqdata[,3]==0,3])
    t01 <- length(probes[probes[,3]<cut & seqdata[,3]==1,3])
    t00 <- length(probes[probes[,3]<cut & seqdata[,3]==0,3])
    if(length(probes[probes[,3]<cut,3])==0){
      type2[nr] <- 0
      type1[nr] <- 1
      conf.t1[nr,] <- 1
      conf.t2[nr,] <- 0
    }
    if(length(probes[probes[,3]>cut,3])==0){
      type2[nr] <- 1
      type1[nr] <- 0
      conf.t1[nr,] <- 0
      conf.t2[nr,] <- 1
    }
    if(length(probes[probes[,3]<cut,3])!=0 & length(probes[probes[,3]>cut,3])!=0){
      type2[nr] <- t01/(t11+t01)
      type1[nr] <- t10/(t10+t00)
      conf.t1[nr,] <- .konfcp(t10,(t10+t00),0.05)[2:3]
      conf.t2[nr,] <- .konfcp(t01,(t11+t01),0.05)[2:3]
    }
    nr <- nr+1
  }
  plot(type1,type2,type="l",col="red",main="Overlap expression based mask - external mask",xlab="Type 1",ylab="Type 2")
  lines(conf.t1[,1],conf.t2[,1],lty=2)
  lines(conf.t1[,2],conf.t2[,2],lty=2)
  if(plotCutoffs) text(type1,type2,labels=signif(cutoffs,2))
  abline(1,-1,col="gray")
  uc=used.cutoff
  uuc=unique(used.cutoff)
  for(ucc in 1:length(uuc)){
    uc[uc==uuc[ucc]]=ucc
  }
  uc=list(uc,used.cutoff)
  ## if(wilcox.ks){
  ##   X11()
  ##   layout(matrix(1:2,ncol=1))
  ##   plot(uc[[1]],kstest,col="red",main="Kolmogorov-Smirnov Test",xlab="Quality score cutoff",
  ##        ylab="p value (Kolmogorov-Smirnov Test)",ylim=c(0,1),pch=16,xaxt="n")
  ##   axis(1,at=1:length(unique(uc[[2]])),labels=signif(unique(uc[[2]]),2),las=3)
  ##   lines(which(unique(uc[[2]]) %in% used.cutoff),kso[!is.na(kso)],type="p",pch=16,cex=0.8)
  ##   plot(uc[[1]],wiltest,col="green",main="Wilcoxon Rank Sum Test",xlab="Quality score cutoff",
  ##        ylab="p value (Wilcoxon Rank Sum Test)",ylim=c(0,1),pch=16,xaxt="n")
  ##   axis(1,at=1:length(unique(uc[[2]])),labels=signif(unique(uc[[2]]),2),las=3)
  ##   lines(which(unique(uc[[2]]) %in% used.cutoff),wilcoxo[!is.na(wilcoxo)],type="p",pch=16,cex=0.8)
  ##   ## X11()
  ##   ## layout(matrix(1:2,ncol=1))
  ##   ## plot(used.cutoff1,kstest1,col="red",main="Kolmogorov-Smirnov Test",xlab="Quality score cutoff",ylab="p value (Kolmogorov-Smirnov Test)",ylim=c(0,1))
  ##   ## plot(used.cutoff1,wiltest1,col="blue",main="Wilcoxon Rank Sum Test",xlab="Quality score cutoff",ylab="p value (Wilcoxon Rank Sum Test)",ylim=c(0,1))
  ## }
  if(wilcox.ks==TRUE){
    return(list(type1=type1,type2=type2,wilcoxonP=wilcoxo,ksP=kso,ksBoot=kstest,wilcoxonBoot=wiltest,cutoffs=cutoffs,testCutoff=uc))
  }
  else{
    return(list(type1=type1,type2=type2,confT1=conf.t1,confT2=conf.t2,cutoffs=cutoffs))
  }
}

prepareMaskedAffybatch <- function(affy,cdfTablePath,exmask="none",cdfName="new_cdf",exclude=NA,cutoff=0.2){
  grouping="probeset"
  .cdf.env.fromTable = function (
    mask.object="none",
    chip.size,
    cdfTablePath="",
    env = new.env(hash=TRUE,parent=baseenv()), 
    make.env = TRUE,
    EXCLUDE.SINGLE
    )  
    {
      if(length(mask.object)==1){
        Ftab = read.table(cdfTablePath, as.is=TRUE, header=TRUE, comment.char="")
        Ftab = data.frame(list(probeset=as.character(Ftab[,1]),x=as.numeric(Ftab[,2]),y=as.numeric(Ftab[,3])))
      }
      else {
        if(sum(paste(mask.object[mask.object[,3]>=cutoff,1],mask.object[mask.object[,3]>=cutoff,2],sep="_") %in% paste(mask.object[mask.object[,3]<cutoff,1],mask.object[mask.object[,3]<cutoff,2],sep="_"))>0){
          cat("Removing probes with contradictorily quality scores\n")
          mask.object=mask.object[!(paste(mask.object[,1],mask.object[,2],sep="_") %in% intersect(paste(mask.object[mask.object[,3]>=cutoff,1],mask.object[mask.object[,3]>=cutoff,2],sep="_"),paste(mask.object[mask.object[,3]<cutoff,1],mask.object[mask.object[,3]<cutoff,2],sep="_"))),]
        }
        Ftab = mask.object[mask.object[,3]>=cutoff,c(4,1,2)]
      }
      cnFtab = colnames(Ftab)

      if (ncol(Ftab) > 3) {	#if there are more than 3 columns, choose the appropriate ones given "probeset"
        Ftab = Ftab[, cnFtab %in% c("probeset", "x", "y", "px", "py", "mx", "my")]
        cnFtab = colnames(Ftab)
      }
      
      
      if (!is.na(EXCLUDE.SINGLE)) {
        gn = Ftab[,1]	#grouping names
        tgn = table(gn)
        include = gn %in% names(tgn)[tgn > EXCLUDE.SINGLE]	#groupings with > 1 probes
        if (sum(include) == 0) {  	#if the CDF table is based on probes 
          stop("Wrong choice: Change exclude.single to FALSE!")
        } else  cat(paste("Excluded ", length(gn) - sum(include), " groupings with ",EXCLUDE.SINGLE," probes\n", sep=""))
        Ftab = Ftab[include,] 
      }
      
      if (ncol(Ftab)==3) {
        cat("Inferring mismatch positions from perfect match\n")
        if (sum(cnFtab == c("probeset", "px", "py")) != ncol(Ftab))
          cat(paste("Renaming columns: ", paste(cnFtab, collapse=" "), "as", "probeset", "px py \n", collapse=" "))
        colnames(Ftab) = c("probeset", "px", "py")
        nn.pm=Ftab$"px"+ chip.size*Ftab$"py"+1
        nn.mm=Ftab$"px" + chip.size*(Ftab$"py"+1)+1
      } else if (ncol(Ftab)==5) {
        if (sum(colnames(Ftab) == c("probeset", "px", "py", "mx", "my")) != ncol(Ftab))
          cat(paste("Renaming columns: ", paste(colnames(Ftab), collapse=" "), "as", "probeset", "px py mx my\n", collapse=" ")) 
        colnames(Ftab) = c("probeset", "px", "py", "mx", "my")
        nn.pm=Ftab$"px"+ chip.size*Ftab$"py"+1
        nn.mm=Ftab$"mx" + chip.size*Ftab$"my"+1
      } else {cat("Wrong table format! Stopping! \n"); stop}
      
      pmm = cbind(nn.pm,nn.mm)
      colnames(pmm) = c("pm","mm")
      
      if (length(unique(Ftab[,1])) != length(Ftab[,1])) {
        cat(paste("Regular run: more than ",EXCLUDE.SINGLE," probe per grouping\n"),sep="")
        l=tapply(1:(dim(pmm)[1]),Ftab[,1],
          function(v) pmm[v,,drop=FALSE])
        cat("Made list.\n")
        if( make.env ) {
          cat("Making env\n")
          nn=names(l)
          for(i in(1:length(nn) ) ) {
            assign(nn[i],l[[i]],pos=env,inherit=FALSE)
          }
          cat("finished making CDF!\n")
          return(env)
        } else {
          cat("finished making CDF!\n")
          return(l)
        }
      } else {
        cat("One probe per grouping\n")
        if( make.env ) {
          cat("Making env\n")
          nn=as.character(Ftab[,1])
          for(i in 1:length(nn)) {
            assign(nn[i],pmm[i,,drop=FALSE],pos=env,inherit=FALSE)
          }
          cat("finished making CDF!\n")
          return(env)
        } else {
          nn=as.character(Ftab[,1])
          l=lapply(1:length(nn), function(i) pmm[i,,drop=FALSE])
          names(l) = nn
          return(l)
        }
        
      }
    }
  
  if(length(exmask)==1){
    CDFenv = .cdf.env.fromTable(chip.size=affy@nrow, cdfTablePath=cdfTablePath, EXCLUDE.SINGLE = exclude)
  }
  else{
    CDFenv = .cdf.env.fromTable( chip.size=affy@nrow, mask.object=exmask, EXCLUDE.SINGLE = exclude)
  }
  affy@cdfName = cdfName
  assign( cdfName(affy), CDFenv, pos=globalenv())
  newAffyBatch = affy
  return(list(newAffyBatch=newAffyBatch,newCdf=CDFenv))
}

plotProbe <- function(affy,probeset,probe=NA,probeXY=NA,scan=TRUE,ind,exmask="none",seqmask="none",names=FALSE){
  v1 <- which(ind==1)
  v2 <- which(ind==2)
  .method1 <- function(x,v1,v2){
    a <- dim(x)
    change <- matrix(1,ncol=a[1],nrow=a[1])
    for(i in 1:a[1]){
      for(j in 1:a[1]){
        if(i!=j){
          lfit <- .regression(x[i,],x[j,],reg="RMA")
          if(lfit$a>0){
            interc <- lfit$b
          }
          if(lfit$a<=0){
            interc <- min(x[i,])
          }
          sum.c <- (x[i,v1]-interc)/x[j,v1]
          sum.h <- (x[i,v2]-interc)/x[j,v2]
          change[i,j] <- t.test(sum.c,sum.h)$p.value
        }
      }
    }
    dummy <- .sort.matrix(change)
    sort.p.value <- dummy$p.value
    names(sort.p.value) <- dummy$sum.order
    return(list(change=change,sort.p.value=sort.p.value))
  }
  .sort.matrix <- function(matrix,SUM=FALSE){
    m <- matrix
    a <- dim(matrix)
    sum.order <- rep(0,a[1])
    tt <- seq(1,a[1])
    p.value <- rep(0,a[1])
    if(SUM==FALSE){
      m <- log(m)
    }
    for(j in seq(1,a[1])){
      sum.sum <- colSums(m)+rowSums(m)
      if(j<a[1]){
        p.value[j] <- min(sum.sum)/(2*a[1]-2*j)
        max.p <- max(sum.sum)/(2*a[1]-2*j)
      }
      if(j==a[1]){
        p.value[j] <- max.p
      }
      order.sum <- order(sum.sum)
      vv <- tt[order.sum[1]]
      tt <- tt[tt!=tt[order.sum[1]]]
      uu <- seq(1,(a[1]-j+1))
      uu <- uu[uu!=order.sum[1]]
      m <- matrix(m[uu,uu],a[1]-j,a[1]-j)
      rownames(m) <- tt
      colnames(m) <- tt
      sum.order[j] <- vv
    }
    if(SUM==FALSE){
      p.value <- exp(p.value)
    }
    names(p.value) <- sum.order
    
    return(list(sum.order=sum.order,p.value=p.value))
  }
  .regression <- function(x,y,reg="RMA",siex=1,siey=1){
    mx <- mean(x)
    my <- mean(y)
    sxsq <- sum((x-mx)*(x-mx))
    sxy <- sum((x-mx)*(y-my))
    sysq <- sum((y-my)*(y-my))
    n <- length(x)
    if(reg=="RMA"){
      beta1 <- sxy/sxsq
      beta2 <- sysq/sxy
      if(sxy!=0){
        a <- sign(sxy)*sqrt(beta1*beta2)
      b <- my-a*mx
      }
      if(sxy==0){
        lm <- lm(y~x)
      a <- lm$coefficients[2]
        b <- lm$coefficients[1]
      }
    }
    return(list(a=a,b=b,n=n,sxsq=sxsq,sysq=sysq,sxy=sxy))
  }
  .mod <- function(x,m){
    t1<-floor(x/m)
    return(x-t1*m)
  }   
  if(is.na(probe) & is.na(probeXY)) cat("Please give the main probe number or its x and y coordinates")
  if(is.na(probe) & !is.na(probeXY) & length(exmask)==1) cat("Please provide exmask object with x/y information")
  if(!is.na(probe) | (!is.na(probeXY) & length(exmask)>1)){
    pm.matrix <- pm(affy,probeset)
    change.ma <- .method1(pm.matrix,v1=v1,v2=v2)$change
    mask.q <- .sort.matrix(change.ma)
    
    nr.probes <- dim(pm.matrix)[1]
    if(length(exmask)!=1){
      probes <- exmask[exmask[,4]==probeset,]
      probes <- cbind(probes,NA)
    }
    if(length(exmask)==1){
      probes <- data.frame(x=NA,y=NA,quality.score=as.numeric(mask.q$p.value[as.character(1:nrow(pm.matrix))]),probeset=probeset,NA)
    }
    if(is.na(probe) & !is.na(probeXY) & !(probeXY %in% paste(probes[,1],".",probes[,2],sep=""))) cat("No probe with the given x/y coordinates in probeset")
    if(!is.na(probe)) random.probe <- probe else random.probe <- which(paste(probes[,1],".",probes[,2],sep="")==probeXY)
    if(is.list(seqmask)){
      seqdata <- seqmask$seqmask
      seqdata <- seqdata[seqdata[,3]==probeset,]
      probes[,5] <- apply(probes,1,function(x){if(sum(paste(seqdata[,1],".",seqdata[,2],sep="")==paste(x[1],".",x[2],sep=""))>0) {
        seqdata[paste(seqdata[,1],".",seqdata[,2],sep="")==paste(x[1],".",x[2],sep=""),4]}
      else NA})
    }
    
    ## layout(matrix(1:15,ncol=3))
    other.probes <- seq(1,nr.probes)[-random.probe]
    if(!scan){
      if(.mod(nr.probes-1,3)==0) layout(matrix(1:(nr.probes-1),ncol=3))
      else{
        message("Reducing number of probes to ",trunc(nr.probes/3)*3)
        other.probes <- sort(sample(other.probes,trunc(nr.probes/3)*3))
        layout(matrix(1:(trunc(nr.probes/3)*3),ncol=3))
        nr.probes <- (trunc(nr.probes/3)*3)+1
      }
    }
    for(i in other.probes){
      if(length(exmask)!=1){
        pxy=paste(" (x: ",probes[random.probe,1]," y: ",probes[random.probe,2],") ",sep="")
        pxy1=paste(" (x: ",probes[i,1]," y: ",probes[i,2],") ",sep="")
      }
      else{
        pxy=""
        pxy1=""
      }
      if(sum(!is.na(probes[c(random.probe,i),5]))>0){
        pxy=paste(pxy," sequence mask: ",probes[random.probe,5],sep="")
        pxy1=paste(pxy1,"sequence mask: ",probes[i,5],sep="")
      }
      if(!names){
         plot(pm.matrix[random.probe,v1],pm.matrix[i,v1],
             main=paste("Masking quality score for pairwise comparison: ",signif(change.ma[random.probe,i],3),"\n",probeset,sep=""),
             xlab=paste("Probe",random.probe,pxy," , expression mask: ",signif(probes[random.probe,3],3)),
             ylab=paste("Probe",i,pxy1," , expression mask: ",signif(probes[i,3],3)),
             xlim=c(min(pm.matrix[random.probe,]),max(pm.matrix[random.probe,])),
             ylim=c(min(pm.matrix[i,]),max(pm.matrix[i,])),col="blue")
        lines(pm.matrix[random.probe,v2],pm.matrix[i,v2],col="red",type="p")
        mtext(c("group 1                   ","                  group 2"),3,col=c("blue","red"))
      }
      else{
        plot(1,1,
             main=paste("Masking quality score for pairwise comparison: ",signif(change.ma[random.probe,i],3),"\n",probeset,sep=""),
             xlab=paste("Probe",random.probe,pxy," , expression mask: ",signif(probes[random.probe,3],3)),
             ylab=paste("Probe",i,pxy1," , expression mask: ",signif(probes[i,3],3)),
             xlim=c(min(pm.matrix[random.probe,])*0.95,max(pm.matrix[random.probe,])*1.05),
             ylim=c(min(pm.matrix[i,]),max(pm.matrix[i,])),col="blue")
        mtext(c("group 1                   ","                  group 2"),3,col=c("blue","red"))
        coltex <- rep("white",max(v1,v2))
        coltex[1:length(v1)] <- "blue"
        coltex[(length(v1)+1):(length(v1)+length(v2))] <- "red"
        text(pm.matrix[random.probe,c(v1,v2)],pm.matrix[i,c(v1,v2)],labels=rownames(affy@phenoData@data)[c(v1,v2)],col=coltex)
      }
      if(scan) scan()
    }
  }
}
