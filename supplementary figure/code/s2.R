library(data.table)
library(ggplot2)
library(patchwork)
library(eoffice)

gdir <- "D:/bald/data/gwas"
rfile <- "D:/bald/postgwas/ldsc/all.rg.res"
out <- "D:/bald/sup"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

labs <- c(bald.qt="Continuous", bald12.bt="M shape", bald13.bt="U shape", bald14.bt="O shape")
files <- file.path(gdir, paste0(names(labs), ".gz"))

readg <- function(f){
  d <- fread(f); setnames(d, names(d), toupper(names(d)))
  n <- names(d); p <- function(x) x[x %in% n][1]
  setnames(d, c(p(c("SNP","RSID","ID")), p(c("EA","A1")), p(c("NEA","A2")),
                p(c("BETA","EFFECT","LOGOR")), p(c("SE","STDERR")), p(c("P","PVAL"))),
           c("SNP","EA","NEA","BETA","SE","P"))
  d[, .(SNP=as.character(SNP), EA=toupper(EA), NEA=toupper(NEA),
        BETA=as.numeric(BETA), SE=as.numeric(SE), P=as.numeric(P))]
}

qq <- function(d, nm){
  p <- sort(d[is.finite(P)&P>0&P<=1, P])
  q <- data.table(e=-log10((seq_along(p)-.5)/length(p)), o=-log10(p))
  ggplot(q, aes(e,o)) + geom_point(size=.25, alpha=.5, color="#2C7FB8") +
    geom_abline(lty=2, linewidth=.3) + labs(x="Expected -log10(P)", y="Observed -log10(P)", title=labs[nm]) +
    theme_bw(11) + theme(panel.grid=element_blank(), plot.title=element_text(hjust=.5))
}

zz <- function(a,b,n1,n2){
  x <- merge(a[P<=5e-8, .(SNP,EA1=EA,NEA1=NEA,Z1=BETA/SE)],
             b[P<=5e-8, .(SNP,EA2=EA,NEA2=NEA,Z2=BETA/SE)], by="SNP")
  x[, f := fifelse(EA1==EA2 & NEA1==NEA2, 1, fifelse(EA1==NEA2 & NEA1==EA2, -1, NA_real_))]
  ggplot(x[!is.na(f)], aes(Z1, Z2*f)) + geom_point(size=1, alpha=.55, color="#3167CD") +
    geom_hline(yintercept=0,lty=2,color="grey60",linewidth=.3) + geom_vline(xintercept=0,lty=2,color="grey60",linewidth=.3) +
    geom_abline(lty=2,color="grey45",linewidth=.3) + labs(x=labs[n1], y=labs[n2], title=paste(labs[n1],"vs",labs[n2])) +
    theme_bw(11) + theme(panel.grid=element_blank(), plot.title=element_text(hjust=.5))
}

g <- setNames(lapply(files, readg), names(labs))

rg <- fread(rfile)[, .(p1=sub(".*/","",p1), p2=sub(".*/","",p2), rg)]
tr <- names(labs); m <- matrix(1,4,4,dimnames=list(tr,tr))
for(i in 1:nrow(rg)) if(rg$p1[i]%in%tr & rg$p2[i]%in%tr) m[rg$p1[i],rg$p2[i]] <- m[rg$p2[i],rg$p1[i]] <- rg$rg[i]
dt <- as.data.table(as.table(m)); setnames(dt,c("x","y","rg"))
dt[, `:=`(x=factor(labs[x], labs[tr]), y=factor(labs[y], rev(labs[tr])), lab=sprintf("%.2f",rg))]
p1 <- ggplot(dt,aes(x,y,fill=rg)) + geom_tile(color="white") + geom_text(aes(label=lab),size=4) +
  scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,name=expression(r[g])) +
  coord_fixed() + labs(x=NULL,y=NULL) + theme_classic(13) + theme(axis.text.x=element_text(angle=45,hjust=1), axis.ticks=element_blank(), axis.line=element_blank())
topptx(p1, file.path(out,"s2.ldsc.pptx"), width=6.5, height=5.5)

p2 <- wrap_plots(lapply(names(g), \(x) qq(g[[x]], x)), ncol=2)
topptx(p2, file.path(out,"s2.qq.pptx"), width=10, height=8)

cmb <- combn(names(g),2,simplify=FALSE)
p3 <- wrap_plots(lapply(cmb, \(x) zz(g[[x[1]]], g[[x[2]]], x[1], x[2])), ncol=3)
topptx(p3, file.path(out,"s2.zz.pptx"), width=13, height=8)
