#setwd("C:/Users/Royce/Documents/R_files/GUI")
#package.skeleton(name="CCTpack", code_files="cctgui.r")
#source("C:/Users/Royce/Documents/R_files/GUI/CCTpack/R/cctgui.r")
#Check effect of thinning on sample numbers  n.keep, n.samp
#DIC, the pD is different from JAGS but when I use their formula on their data, I get the same as mine.

#Example
#cctgui() 
#write.csv(x=testdat, file="testdat.csv",row.names=FALSE)

cctgui <- function(){
x <- y <- z <- dat <- clusters <- itemdiff <- whmodel <- Mode <- jagsfit.prealg <- jagsfit <- checksrunbefore <- NULL
alltraceplot <- dtraceplot <- NULL
rm(x,y,z,dat,clusters,itemdiff,whmodel,Mode,jagsfit.prealg,jagsfit,checksrunbefore)
rm(alltraceplot,dtraceplot)
instant_pkgs <- function(pkgs) { 
    pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
    if (length(pkgs_miss) > 0) {
	message("\n ...Some packages are missing, please follow the prompt to install them.\n")
        install.packages(pkgs_miss)
    }
    
	if (length(pkgs_miss) == 0) {
	#message("\n ...Packages were already installed!\n")
	}

    # install packages not already loaded:
    pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
    if (length(pkgs_miss) > 0) {
        install.packages(pkgs_miss)
    }
    
    # load packages not already loaded:
    attached <- search()
    attached_pkgs <- attached[grepl("package", attached)]
    need_to_attach <- pkgs[which(!pkgs %in% gsub("package:", "", attached_pkgs))]
    
    if (length(need_to_attach) > 0) {
      for (i in 1:length(need_to_attach)) require(need_to_attach[i], character.only = TRUE)
    }

	if (length(need_to_attach) == 0) {
	#message("\n ...Packages were already loaded!\n")
	message("\n ...Starting CCT Inference Software\n")
	}
}

instant_pkgs(c("tcltk", "psych","rjags","R2jags","mvtnorm","polycor"))


######################
#MC-CCT Application
#suppressPackageStartupMessages(library("tcltk"))
#suppressPackageStartupMessages(library("psych"))
#suppressPackageStartupMessages(library(rjags))
#suppressPackageStartupMessages(library(R2jags))


######################
#Variables
samplesvar <- tclVar("10000") #good default for CRM is 12000 samples, 3000 burn in
chainsvar <- tclVar("3")
burninvar <- tclVar("2000")
thinvar <- tclVar("1")
culturesvar <- tclVar("1") #cultures.entry
#cctpackdir <- "CCTpack/"

tt <- tktoplevel()
#datframe <- tkframe(tt, borderwidth = 0, width=800, height=100,relief="ridge")
datframe <- tkframe(tt, borderwidth = 0)
datframe2 <- tkframe(tt, borderwidth = 0)
datframe3 <- tkframe(tt, borderwidth = 0)
applyframe <- tkframe(tt, borderwidth = 0)
applyframe2 <- tkframe(tt, borderwidth = 0)
applyframe3 <- tkframe(tt, borderwidth = 0)
applyframe4 <- tkframe(tt, borderwidth = 0)
resultsframe <- tkframe(tt, borderwidth = 0)
settingsframe <- tkframe(tt, borderwidth = 0)


tkwm.title(tt,"CCT Model Application Software")
samples.entry <- tkentry(settingsframe, textvariable=samplesvar,width="6")
chains.entry <- tkentry(settingsframe, textvariable=chainsvar,width="2")
burnin.entry <- tkentry(settingsframe, textvariable=burninvar,width="6")
thin.entry <- tkentry(settingsframe, textvariable=thinvar,width="2")
cultures.entry <- tkentry(applyframe, textvariable=culturesvar,width="2")
#cultures.entry <- tktext(applyframe, textvariable=culturesvar,width="2")


#######################
#Buttons
loadfilefunc <- function(){
## <<- convention sets it globally
fileName <<- file.path(tclvalue(tkgetOpenFile(filetypes = "{{csv Files} {.csv .txt}}")))
if (!nchar(fileName)) {
#tkmessageBox(message = "No file was selected!")
return()
}
else {
dat <- as.matrix(read.csv(fileName,header=FALSE))
options(warn=-3)
if(!is.numeric(dat[1,1])){dat <- matrix(as.numeric(dat),dim(dat)[1],dim(dat)[2])
if(all(is.na(dat[,1]))){dat <- dat[,-1]}
if(all(is.na(dat[1,]))){dat <- dat[-1,]}
}
options(warn=0)
if(all(dat[1:dim(dat)[1],1] == 1:dim(dat)[1])){dat <- dat[,-1]}

setwd(file.path(dirname(fileName)))

tkconfigure(screeplot.but, state="normal") #"normal" / "disabled"
tkconfigure(applymodel.but, state="normal") #"normal" / "disabled"

#tkmessageBox(message = paste("The file selected was", fileName))
}
tkconfigure(datafiletxt, state="normal") #"disabled"
tkconfigure(resptxt, state="normal") #"disabled"
tkconfigure(itemtxt, state="normal") #"disabled"
tkconfigure(dattypetxt, state="normal") #"disabled"
tkconfigure(modeltxt, state="normal") #"disabled"

tkdelete(datafiletxt,"1.0","800.0")
tkdelete(resptxt,"1.0","800.0")
tkdelete(itemtxt,"1.0","800.0")
tkdelete(dattypetxt,"1.0","800.0")
tkdelete(modeltxt,"1.0","800.0")


tkinsert(datafiletxt,"end",basename(fileName)) #basename(), dirname()

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

if(!all(is.wholenumber(dat))){ 
tkinsert(modeltxt,"end","CRM") 
whmodel <<- "CRM"
if(min(dat) >= 0 && max(dat) <= 1){
tkinsert(dattypetxt,"end","Continuous in [0,1]")
#dat convert to logit
dat[dat==0] <- .01; dat[dat==1] <- .99
dat <- logit(dat) 
}else{tkinsert(dattypetxt,"end","Continuous")}     

}else{if(min(dat) >= 0 && max(dat) <= 1){  
tkinsert(modeltxt,"end","GCM") 
whmodel <<- "GCM"
tkinsert(dattypetxt,"end","Dichotomous 0 & 1")
}else{ 
tkinsert(modeltxt,"end","LTRM") 
whmodel <<- "LTRM"
tkinsert(dattypetxt,"end","Categorical 1, 2, ...")
}}

tkinsert(resptxt,"end",dim(dat)[1]) #basename(), dirname()
tkinsert(itemtxt,"end",dim(dat)[2]) #basename(), dirname()

tkconfigure(datafiletxt, state="disabled") #"disabled"
tkconfigure(resptxt, state="disabled") #"disabled"
tkconfigure(itemtxt, state="disabled") #"disabled"
tkconfigure(dattypetxt, state="disabled") #"disabled"
tkconfigure(modeltxt, state="disabled") #"disabled"

dat <<- dat

}
submitfunc <- function() {
samples <- as.numeric(tclvalue(samplesvar))
chains <- as.numeric(tclvalue(chainsvar))
burnin <- as.numeric(tclvalue(burninvar))
thin <- as.numeric(tclvalue(thinvar))
tkmessageBox(message=paste("x + y + z = ", x+y+z, ""))
}
resetfunc <- function() {
tclvalue(samplesvar)<-""
tclvalue(chainsvar)<-""
tclvalue(burninvar)<-""
}
simdatfunc <- function(V=3,n=10,m=30, id=0,stats=1) {
library(cluster); library(psych)
dat <- list(); dat$V <- V; dat$n <- V*n; if(length(n) > 1){dat$n <- sum(n)}; dat$m <- m

datagencrm <- function(h=1,n=10,m=20,id=1,opt=0,lam=0) {
dat <- list()
dat$Tmu <- runif(1,-3,3); dat$Ttau <- runif(1,.5,2)
dat$amu <- 1; dat$atau <- runif(1,.5,3)
dat$bmu <- 0; dat$btau <- runif(1,.5,3)
dat$Emu <- runif(1,1,3); dat$Etau <- runif(1,.5,2)
dat$T <- rnorm(m,dat$Tmu,dat$Ttau^-.5)
dat$E <- sort(rgamma(n,(dat$Emu^2)*dat$Etau,dat$Emu*dat$Etau))
dat$a <- rgamma(n,(dat$amu^2)*dat$atau,dat$amu*dat$atau)
dat$b <- rnorm(n,dat$bmu,dat$btau^-.5) 
if(id==1){
dat$lammu <- 1; dat$lamtau <- runif(1,.5,2)
dat$lam <- rgamma(m,(dat$lammu^2)*dat$lamtau,dat$lammu*dat$lamtau)
if(lam[1] > 0){dat$lam <- lam}
}
dat$Y <- matrix(-1,n,m); dat$X <- matrix(-1,n,m); dat$ilogitY <- matrix(-1,n,m); dat$ilogitX <- matrix(-1,n,m)
invlogit <- function(x){x <- 1 / (1 + exp(-x)); return(x)}

for(i in 1:n) {
    for(k in 1:m) {
 	dat$X[i,k] <- rnorm(1,dat$T[k],dat$E[i]^(-.5) )
	if(id==1){dat$X[i,k] <- rnorm(1,dat$T[k],(dat$E[i]/dat$lam[k])^(-.5) )}
        dat$Y[i,k] <- (dat$a[i]*dat$X[i,k])+dat$b[i]
	dat$ilogitY[i,k] <- invlogit(dat$Y[i,k]); dat$ilogitX[i,k] <- invlogit(dat$X[i,k]) 		 	 		 
       }}
return(dat)
}

lam <- 0
for(v in 1:V){
if(length(n) > 1){
sd1 <- datagencrm(h=1,n=n[v],m=m,id=id,lam) 
} else{sd1 <- datagencrm(h=1,n=n,m=m,id=id)}

# A single item difficulty for all cultures needs to be input, rather than separate

dat$T <- cbind(dat$T,sd1$T); #if(id==1){dat$lam <- c(dat$lam,sd1$lam); dat$lammu <- c(dat$lammu,sd1$lammu); dat$lamtau <- c(dat$lamtau,sd1$lamtau) }
if(id==1 && v==1){dat$lam <- sd1$lam; lam <- sd1$lam; dat$lammu <- sd1$lammu; dat$lamtau <- sd1$lamtau }
dat$E <- c(dat$E,sd1$E); dat$a <- c(dat$a,sd1$a); dat$b <- c(dat$b,sd1$b)
sd1$e <- runif(length(sd1$E),v,v); dat$e <- c(dat$e,sd1$e)
dat$Tmu <- c(dat$Tmu,sd1$Tmu); dat$Ttau <- c(dat$Ttau,sd1$Ttau)
dat$amu <- c(dat$amu,sd1$amu); dat$atau <- c(dat$atau,sd1$atau)
dat$bmu <- c(dat$bmu,sd1$bmu); dat$btau <- c(dat$btau,sd1$btau)
dat$Emu <- c(dat$Emu,sd1$Emu); dat$Etau <- c(dat$Etau,sd1$Etau)
dat$pi <- c(dat$pi,sd1$pi); dat$Y <- rbind(dat$Y,sd1$Y); dat$X <- rbind(dat$X,sd1$X)
dat$V <- sd1$V
}
rm(sd1)
if(stats==1){
par(mfrow=c(1,3))
#plot(fanny(cor(t(dat$Y)),T,maxit = 1000))
plot(fa(cor(t(dat$Y)))$values[1:8])}

return(dat)
}
recovplot <- function(jagsfit,gen,mfrow=c(2,3)) {
par(mfrow = mfrow,family="serif",mar=c(3,4,2.5,2),oma=c(0,0,0,0),mgp=c(1.75,.5,0)) 

hdi <- function(sampleVec, credMass=0.95){
sortedPts = sort(sampleVec)
ciIdxInc = floor(credMass * length( sortedPts) )
nCIs = length( sortedPts) - ciIdxInc
ciwidth = rep(0, nCIs)
for(i in 1:nCIs){
	ciwidth[i] = sortedPts[i + ciIdxInc] - sortedPts[i]}
HDImin = sortedPts[which.min(ciwidth)]
HDImax = sortedPts[which.min(ciwidth)+ciIdxInc]
HDIlim = c(HDImin,HDImax)
return(HDIlim) 
}

CIplot <- function(jagsfit, gen, param="T"){

if(length(dim(jagsfit$BUGSoutput$mean[[param]])) == 1){jagsfit$BUGSoutput$mean[[param]]<-array(jagsfit$BUGSoutput$mean[[param]],c(length(jagsfit$BUGSoutput$mean[[param]]),1))
jagsfit$BUGSoutput$sims.list[[param]] <- array(jagsfit$BUGSoutput$sims.list[[param]],c(dim(jagsfit$BUGSoutput$sims.list[[param]])[1],dim(jagsfit$BUGSoutput$sims.list[[param]])[2],1))
gen[[param]] <- array(gen[[param]],c(length(gen[[param]]),1))
} 
for(v in 1:dim(jagsfit$BUGSoutput$mean[[param]])[2]){
hditmp <- apply(jagsfit$BUGSoutput$sims.list[[param]][,,v],2,hdi)
plot(1:length(jagsfit$BUGSoutput$mean[[param]][,v]),jagsfit$BUGSoutput$mean[[param]][,v],ylim=c(min(hditmp[1,]),max(hditmp[2,])),main=paste(param,"[",v,"k] Parameter Node"),ylab="Posterior Mean",xlab="Node")
#points(1:length(jagsfit$BUGSoutput$mean[[param]][,v]),gen[[param]][,v],pch=22)
points(1:length(jagsfit$BUGSoutput$mean[[param]][,v]),gen[[param]][,which(cor(gen[[param]],jagsfit$BUGSoutput$mean[[param]][,v]) == max(cor(gen[[param]],jagsfit$BUGSoutput$mean[[param]][,v])))
],pch=22)


segments(1:jagsfit$m,hditmp[1,], 1:jagsfit$m, hditmp[2,])
arrows(1:jagsfit$m,hditmp[1,], 1:jagsfit$m, hditmp[2,],code=3,angle=90,length=.025)
} }

CIplot(jagsfit,gen,param="T")
CIplot(jagsfit,gen,param="lam")
CIplot(jagsfit,gen,param="E")
CIplot(jagsfit,gen,param="a")
CIplot(jagsfit,gen,param="b")
}
quitfunc <- function() {
#q(save = "no")  #quits R
tkdestroy(tt)
}
screeplotfunc <- function(saveplots=0,savedir=0) {
options(warn=-3)

if(whmodel != "LTRM"){
if(saveplots==0){message("\n ...Producing Scree Plot\n")}
#if(saveplots==1){jpeg(file.path(getwd(),"CCTpack","CCTpackscree.jpg"),width = 6, height = 6, units = "in", pointsize = 12,quality=100,res=400)}
#if(saveplots==2){postscript(file=file.path(getwd(),"CCTpack","CCTpackscree.eps"), onefile=FALSE, horizontal=FALSE, width = 6, height = 6, paper="special", family="Times")}
if(saveplots==1){jpeg(file.path(gsub(".Rdata","scree.jpg",savedir)),width = 6, height = 6, units = "in", pointsize = 12,quality=100,res=400)}
if(saveplots==2){postscript(file=file.path(gsub(".Rdata","scree.eps",savedir)), onefile=FALSE, horizontal=FALSE, width = 6, height = 6, paper="special", family="Times")}
par(oma=c(0,0,0,0),mar=c(4,4,3,1),mgp=c(2.25,.75,0),mfrow=c(1,1))
suppressMessages(plot(fa(cor(t(dat)))$values[1:8],las=1,type="b",bg="black",pch=21,xlab="Factor",ylab="Magnitude",main="Scree Plot of Data"))
if(saveplots==1 || saveplots ==2){dev.off()}
}else{
if(saveplots==0){message("\n ...Producing Scree Plot\n")}
if(saveplots==1){jpeg(file.path(gsub(".Rdata","scree.jpg",savedir)),width = 6, height = 6, units = "in", pointsize = 12,quality=100,res=400)}
if(saveplots==2){postscript(file=file.path(gsub(".Rdata","scree.eps",savedir)), onefile=FALSE, horizontal=FALSE, width = 6, height = 6, paper="special", family="Times")}
if(!exists('datfactors')){datfactors <- suppressMessages(fa(polychoric(t(dat),polycor=TRUE)$rho)$values[1:8])}
par(oma=c(0,0,0,0),mar=c(4,4,3,1),mgp=c(2.25,.75,0),mfrow=c(1,1))
plot(datfactors,las=1,type="b",bg="black",pch=21,xlab="Factor",ylab="Magnitude",main="Scree Plot of Data")
datfactors <<- datfactors
#suppressMessages(plot(fa(polychoric(t(dat),polycor=TRUE)$rho)$values[1:8],las=1,type="b",bg="black",pch=21,xlab="Factor",ylab="Magnitude",main="Scree Plot of Data"))

if(saveplots==1 || saveplots ==2){dev.off()}
}

options(warn=0)
}
applymodelfunc <- function() {
jags.iter <- as.numeric(tclvalue(samplesvar))
jags.chains <- as.numeric(tclvalue(chainsvar))
jags.burnin <- as.numeric(tclvalue(burninvar))
jags.thin <- as.numeric(tclvalue(thinvar))
clusters <<- as.numeric(tclvalue(culturesvar))
itemdiff <<- as.numeric(tclvalue(itemdiffvar))
if(exists("checksrunbefore") == TRUE){rm(list = ls(envir=globalenv())[
             grep("checksrunbefore", ls(envir=globalenv()))], envir = globalenv())}
#rm("checksrunbefore")}


if(whmodel == "GCM"){
Y <- dat; n <- dim(Y)[1]; m <- dim(Y)[2]; T <- clusters
jags.data <- list("Y","n","m","T")

if(itemdiff==0){
model.file <- mcgcm
jags.params <- c("z","D","g","p","dmu","dth","gmu","gth","e","pi")
if(clusters>1){
jags.inits <- function(){ list("z"=matrix(rbinom(m*T,1,.5),m,T),"D"= runif(n,.2,.8), "g"= runif(n,.2,.8),"e"= sample(1:T,n,replace=TRUE) )}
}else{
jags.inits <- function(){ list("z"=matrix(rbinom(m*T,1,.5),m,T),"D"= runif(n,.2,.8), "g"= runif(n,.2,.8) )}
}

}
if(itemdiff==1){
model.file <- mcgcmid
jags.params <- c("z","D","g","lam","p","dmu","dth","gmu","gth","lammu","lamth","e","pi")
if(clusters>1){
jags.inits <- function(){ list("z"=matrix(rbinom(m*T,1,.5),m,T),"D"= runif(n,.2,.8), "g"= runif(n,.2,.8), "lam"= runif(m,.2,.8), "e"= sample(1:T,n,replace=TRUE) )}
}else{
jags.inits <- function(){ list("z"=matrix(rbinom(m*T,1,.5),m,T),"D"= runif(n,.2,.8), "g"= runif(n,.2,.8), "lam"= runif(m,.2,.8) )}
}
}
if(clusters==1){
model.file <- gsub(pattern="pi\\[1\\:T\\] ~ ddirch\\(L\\)", "pi <- 1", model.file)
model.file <- gsub(pattern="e\\[i\\] ~ dcat\\(pi\\)", "e\\[i\\] <- 1", model.file)
}
}

if(whmodel == "LTRM"){
Y <- dat; n <- dim(Y)[1]; m <- dim(Y)[2]; V <- clusters; C <- max(Y)
jags.data <- list("Y","n","m","C","V")

if(itemdiff==0){
model.file <- mcltrm
jags.params <- c("T","gam","E","a","b","Tmu","Ttau","Emu","Etau","amu","atau","bmu","btau","e","pi")
if(clusters>1){
jags.inits <- function(){ list("T"=matrix(rnorm(m*V,0,1),m,V), "tgam"=matrix(rnorm((C-2)*V,0,1),(C-2),V), "E"= runif(n,.8,1.2), "a"= runif(n,.8,1.2), "b"= rnorm(n,0,1), "e"= sample(1:V,n,replace=TRUE)  )}
}else{
jags.inits <- function(){ list("T"=matrix(rnorm(m*V,0,1),m,V), "tgam"=matrix(rnorm((C-2)*V,0,1),(C-2),V), "E"= runif(n,.8,1.2), "a"= runif(n,.8,1.2), "b"= rnorm(n,0,1)  )}
}

if(C==2){
if(clusters>1){
jags.inits <- function(){ list("T"=matrix(rnorm(m*V,0,1),m,V), "E"= runif(n,.8,1.2), "b"= rnorm(n,0,1), "e"= sample(1:V,n,replace=TRUE)  )}#"tgam"=matrix(rnorm((C-2)*V,0,1),(C-2),V), "e"= (rbinom(n,(V-1),.5)+1))}
}else{
jags.inits <- function(){ list("T"=matrix(rnorm(m*V,0,1),m,V), "E"= runif(n,.8,1.2), "b"= rnorm(n,0,1)  )}
}
model.file <- gsub(pattern="a\\[i\\] ~ dgamma\\(atau\\[e\\[i\\]\\],atau\\[e\\[i\\]\\]\\)", "a\\[i\\] <- 1", model.file)
model.file <- gsub(pattern="atau\\[v\\] ~ dgamma\\(4,4\\)", "atau\\[v\\] <- 0", model.file)
model.file <- gsub(pattern="for \\(c in 2\\:\\(C-1\\)\\)", " #for \\(c in 2\\:\\(C-1\\)\\)", model.file)
model.file <- gsub(pattern="gam\\[1\\:\\(C-1\\),v\\] <- sort\\(tgam2\\[1\\:\\(C-1\\),v\\]\\)", "gam\\[1,v\\] <- 0 }", model.file)
model.file <- gsub(pattern="for \\(c in 1\\:\\(C-2\\)\\)", "#for \\(c in 1\\:\\(C-2\\)\\)", model.file)
model.file <- gsub(pattern="\\(1 - sum\\(pY\\[i,k,1\\:\\(C-1\\)\\]\\)\\)", "1 - pY\\[i,k,1\\]", model.file)
model.file <- gsub(pattern="tgam2\\[1", "#tgam2\\[1", model.file)
model.file <- gsub(pattern="tgam2\\[C", "#tgam2\\[C", model.file)
}
}

if(itemdiff==1){
model.file <- mcltrmid
jags.params <- c("T","lam","gam","E","a","b","Tmu","Ttau","Emu","Etau","amu","atau","bmu","btau","lammu","lamtau","e","pi")
if(clusters>1){
jags.inits <- function(){ list("T"=matrix(rnorm(m*V,0,1),m,V), "lam"= runif(m,.8,1.2), "tgam"=matrix(rnorm((C-2)*V,0,1),(C-2),V), "E"= runif(n,.8,1.2), "a"= runif(n,.8,1.2), "b"= rnorm(n,0,1), "e"= sample(1:V,n,replace=TRUE) )}#,"e"= (rbinom(n,(V-1),.5)+1))}
}else{
jags.inits <- function(){ list("T"=matrix(rnorm(m*V,0,1),m,V), "lam"= runif(m,.8,1.2), "tgam"=matrix(rnorm((C-2)*V,0,1),(C-2),V), "E"= runif(n,.8,1.2), "a"= runif(n,.8,1.2), "b"= rnorm(n,0,1) )}
}
if(C==2){
if(clusters>1){
jags.inits <- function(){ list("T"=matrix(rnorm(m*V,0,1),m,V), "lam"= runif(m,.8,1.2), "E"= runif(n,.8,1.2), "b"= rnorm(n,0,1), "e"= sample(1:V,n,replace=TRUE) )}#"tgam"=matrix(rnorm((C-1)*V,0,1),(C-1),V),"e"= (rbinom(n,(V-1),.5)+1))} 
}else{
jags.inits <- function(){ list("T"=matrix(rnorm(m*V,0,1),m,V), "lam"= runif(m,.8,1.2), "E"= runif(n,.8,1.2), "b"= rnorm(n,0,1) )}
}
model.file <- gsub(pattern="a\\[i\\] ~ dgamma\\(atau\\[e\\[i\\]\\],atau\\[e\\[i\\]\\]\\)", "a\\[i\\] <- 1", model.file)
model.file <- gsub(pattern="atau\\[v\\] ~ dgamma\\(4,4\\)", "atau\\[v\\] <- 0", model.file)
model.file <- gsub(pattern="for \\(c in 2\\:\\(C-1\\)\\)", " #for \\(c in 2\\:\\(C-1\\)\\)", model.file)
model.file <- gsub(pattern="gam\\[1\\:\\(C-1\\),v\\] <- sort\\(tgam2\\[1\\:\\(C-1\\),v\\]\\)", "gam\\[1,v\\] <- 0 }", model.file)
model.file <- gsub(pattern="for \\(c in 1\\:\\(C-2\\)\\)", "#for \\(c in 1\\:\\(C-2\\)\\)", model.file)
model.file <- gsub(pattern="\\(1 - sum\\(pY\\[i,k,1\\:\\(C-1\\)\\]\\)\\)", "1 - pY\\[i,k,1\\]", model.file)
model.file <- gsub(pattern="tgam2\\[1", "#tgam2\\[1", model.file)
model.file <- gsub(pattern="tgam2\\[C", "#tgam2\\[C", model.file)
}
}

if(clusters==1){ 
model.file <- gsub(pattern="pi\\[1\\:V\\] ~ ddirch\\(L\\)", "pi <- 1", model.file)
model.file <- gsub(pattern="e\\[i\\] ~ dcat\\(pi\\)", "e\\[i\\] <- 1", model.file)
}
}

if(whmodel == "CRM"){
Y <- dat; n <- dim(Y)[1]; m <- dim(Y)[2]; V <- clusters
jags.data <- list("Y","n","m","V")

if(itemdiff==0 ){
model.file <- mccrm
jags.params <- c("T","E","a","b","Tmu","Ttau","Emu","Etau","amu","atau","bmu","btau","e","pi")
if(clusters>1){
jags.inits <- function(){ list("T"=matrix(rnorm(m*V,0,1),m,V),"E"= runif(n,.8,1.2), "a"= runif(n,.8,1.2), "b"= rnorm(n,0,1), "e"= sample(1:V,n,replace=TRUE)   )}#, "e"= (rbinom(n,(V-1),.5)+1))}
}else{
jags.inits <- function(){ list("T"=matrix(rnorm(m*V,0,1),m,V),"E"= runif(n,.8,1.2), "a"= runif(n,.8,1.2), "b"= rnorm(n,0,1)  )}
}
}

if(itemdiff==1){
model.file <- mccrmid
jags.params <- c("T","E","a","b","lam","Tmu","Ttau","Emu","Etau","amu","atau","bmu","btau","lammu","lamtau","e","pi")
if(clusters>1){
jags.inits <- function(){ list("T"=matrix(rnorm(m*V,0,1),m,V),"E"= runif(n,.8,1.2), "a"= runif(n,.8,1.2), "b"= rnorm(n,0,1), "lam"= runif(m,.8,1.2), "e"= sample(1:V,n,replace=TRUE)   )}#,"e"= (rbinom(n,(V-1),.5)+1))}
}else{
jags.inits <- function(){ list("T"=matrix(rnorm(m*V,0,1),m,V),"E"= runif(n,.8,1.2), "a"= runif(n,.8,1.2), "b"= rnorm(n,0,1), "lam"= runif(m,.8,1.2)  )}
}
}
if(clusters==1){ 
model.file <- gsub(pattern="pi\\[1\\:V\\] ~ ddirch\\(L\\)", "pi <- 1", model.file)
model.file <- gsub(pattern="e\\[i\\] ~ dcat\\(pi\\)", "e\\[i\\] <- 1", model.file)
} 
}

jagsfit <- jags(data=jags.data, inits=jags.inits, parameters.to.save=jags.params,
n.chains=jags.chains, n.iter=(jags.iter+jags.burnin), n.burnin=jags.burnin, 
n.thin=jags.thin, model.file=textConnection(model.file))
jagsfit$data <- Y; jagsfit$n <- n; jagsfit$m <- m; 
if(whmodel=="GCM"){jagsfit$T <- T}
if(whmodel=="LTRM"){jagsfit$V <- V; jagsfit$C <- C}
if(whmodel=="CRM"){jagsfit$V <- V}

#jagsfit1 <<- jagsfit

Rhat1 <- function(mat) {
  m <- ncol(mat)
  n <- nrow(mat)
  b <- apply(mat,2,mean)
  B <- sum((b-mean(mat))^2)*n/(m-1)
  w <- apply(mat,2,var)
  W <- mean(w)
  s2hat <- (n-1)/n*W + B/n
  Vhat <- s2hat + B/m/n 
  covWB <- n /m * (cov(w,b^2)-2*mean(b)*cov(w,b))
  varV <- (n-1)^2 / n^2 * var(w)/m +
          (m+1)^2 / m^2 / n^2 * 2*B^2/(m-1) +
          2 * (m-1)*(n-1)/m/n^2 * covWB
  df <- 2 * Vhat^2 / varV
  R <- sqrt((df+3) * Vhat / (df+1) / W)
  return(R)
}

Rhat <- function(arr) {
  dm <- dim(arr)
  if (length(dm)==2) return(Rhat1(arr))
  if (dm[2]==1) return(NULL)
  if (dm[3]==1) return(Rhat1(arr[,,1]))
  return(apply(arr,3,Rhat1))
}

# Rhat1discrete <- function(mat) {
  # m <- ncol(mat)
  # n <- nrow(mat)
  # b <- apply(mat,2,mean)
  # B <- sum((b-mean(mat))^2)*n/(m-1)
  # w <- apply(mat,2,var)
  # W <- mean(w)
  # s2hat <- (n-1)/n*W + B/n
  # Vhat <- s2hat + B/m/n 
  # covWB <- n /m * (cov(w,b^2)-2*mean(b)*cov(w,b))
  # varV <- (n-1)^2 / n^2 * var(w)/m +
          # (m+1)^2 / m^2 / n^2 * 2*B^2/(m-1) +
          # 2 * (m-1)*(n-1)/m/n^2 * covWB
  # df <- 2 * Vhat^2 / varV
  # R <- sqrt((df+3) * Vhat / (df+1) / W)
  # #R <- sqrt(Vhat/W)
  # return(R)
# }

# Rhatdiscrete <- function(arr) {
  # dm <- dim(arr)
  # if (length(dm)==2) return(Rhat1discrete(arr))
  # if (dm[2]==1) return(NULL)
  # if (dm[3]==1) return(Rhat1discrete(arr[,,1]))
  # return(apply(arr,3,Rhat1discrete))
# }

acfsum <- function(jagsfit){ 
ACF <- matrix(-1,dim(jagsfit$BUGSoutput$sims.array)[3],5)
for(i in 1:dim(jagsfit$BUGSoutput$sims.array)[3]){ 
ACF[i,1:5] <- acf(jagsfit$BUGSoutput$sims.array[,,i],plot=FALSE,lag.max=6,na.action=na.pass)$acf[2:6]
} 
acmat <- array(-1,5,dimnames=list(c("% >.2","% >.15","% >.1","% >.1","% >.05")) ) 
acmat[1] <- 100*sum(abs(ACF[,1]) > .2,na.rm=TRUE)/ sum(ACF[,1] >0,na.rm=TRUE)
acmat[2] <- 100*sum(abs(ACF[,2]) > .15,na.rm=TRUE)/ sum(ACF[,2] >0,na.rm=TRUE)
acmat[3] <- 100*sum(abs(ACF[,3]) > .1,na.rm=TRUE)/ sum(ACF[,3] >0,na.rm=TRUE)
acmat[4] <- 100*sum(abs(ACF[,4]) > .1,na.rm=TRUE)/ sum(ACF[,4] >0,na.rm=TRUE)
acmat[5] <- 100*sum(abs(ACF[,5]) > .05,na.rm=TRUE)/ sum(ACF[,5] >0,na.rm=TRUE)
return(acmat) 
}

labelswitchalggcm <- function(jagsfit,chnind=0){

jagsfit2 <- jagsfit

nch <- jagsfit2$BUGSoutput$n.chains
nsamp <- jagsfit2$BUGSoutput$n.keep

if(nch != 1){
ntruths <- jagsfit2$T
truths <- array(NA,c(jagsfit2$m,nch,ntruths))
inds <- jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "z")]]
inds <- matrix(inds,jagsfit2$m,jagsfit2$T)
for(t in 1:jagsfit$T){truths[,,t] <- t(apply(jagsfit$BUGSoutput$sims.array[ ,,inds[,t]],c(2,3),mean))}

T <- jagsfit2$T

chstart <- 1
if(length(chnind)==1){ 
chstart <- 2
chnind <- array(NA,c(T,nch))
chnind[1:T,1] <- 1:T

for(t in 1:T){
for(ch in chstart:nch){
Tind <- c(1:T)[-chnind[,ch][!is.na(chnind[,ch])]]
if(length(Tind)==0){Tind <- c(1:T)}

chnind[t,ch] <- which(max(cor(truths[1:jagsfit2$m, 1, t],truths[1:jagsfit2$m, ch, Tind]))==cor(truths[1:jagsfit2$m, 1, t],truths[1:jagsfit2$m, ch, ]))
}}
}

nsamp <- jagsfit$BUGSoutput$n.keep

inds <- rbind(inds,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "dmu")]],
jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "dth")]],
jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "gmu")]],
jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "gth")]],
jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "p")]],
jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "pi")]])

einds <- jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]]
tmpeinds <- jagsfit$BUGSoutput$sims.array[ ,,einds]

for(t in 1:jagsfit$T){
for(ch in chstart:nch){
jagsfit2$BUGSoutput$sims.array[ ,ch,inds[,t]] <- jagsfit$BUGSoutput$sims.array[ ,ch,inds[,chnind[t,ch]]]  #this jagsfit2 vs. jagsfit difference is intentional
if(chstart==2){jagsfit2$BUGSoutput$sims.array[ ,ch,einds][jagsfit$BUGSoutput$sims.array[ ,ch,einds] == chnind[t,ch]] <- chnind[t,1]}
if(chstart==1){jagsfit2$BUGSoutput$sims.array[ ,ch,einds][jagsfit$BUGSoutput$sims.array[ ,ch,einds] == t] <- chnind[t,1]}
}}

}

#jagsfit2$BUGSoutput$sims.array[ ,,einds] <- tmpeinds
jagsfit2$BUGSoutput$sims.list[["e"]] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]]], c(nsamp*nch,jagsfit$n))
jagsfit2$BUGSoutput$mean$e <- apply(jagsfit2$BUGSoutput$sims.list[["e"]],2,mean)

jagsfit2$BUGSoutput$sims.matrix <- array(jagsfit2$BUGSoutput$sims.array,c(nsamp*nch,dim(jagsfit$BUGSoutput$sims.array)[3]))  #this jagsfit2 vs. jagsfit difference is intentional

jagsfit2$BUGSoutput$sims.list[["D"]] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "D")]]], c(nsamp*nch,jagsfit$n))
jagsfit2$BUGSoutput$sims.list[["g"]] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "g")]]], c(nsamp*nch,jagsfit$n))
if(itemdiff == 1){jagsfit2$BUGSoutput$sims.list[["lam"]] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "lam")]]], c(nsamp*nch,jagsfit$m))}

jagsfit2$BUGSoutput$sims.list[["z"]] <- array(NA, c(nsamp*nch,jagsfit$m,jagsfit$T))
jagsfit2$BUGSoutput$sims.list[["dmu"]] <- array(NA, c(nsamp*nch,jagsfit$T))
jagsfit2$BUGSoutput$sims.list[["dth"]] <- array(NA, c(nsamp*nch,jagsfit$T))
jagsfit2$BUGSoutput$sims.list[["gmu"]] <- array(NA, c(nsamp*nch,jagsfit$T))
jagsfit2$BUGSoutput$sims.list[["gth"]] <- array(NA, c(nsamp*nch,jagsfit$T))
jagsfit2$BUGSoutput$sims.list[["p"]] <- array(NA, c(nsamp*nch,jagsfit$T))
jagsfit2$BUGSoutput$sims.list[["pi"]] <- array(NA, c(nsamp*nch,jagsfit$T))

for(t in 1:T){
jagsfit2$BUGSoutput$sims.list[["z"]][,,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "z")]][1:jagsfit$m +((t-1)*jagsfit$m)]], c(nsamp*nch,jagsfit$m))
#jagsfit2$BUGSoutput$sims.list[["z"]][,,t] <- array(jagsfit2$BUGSoutput$sims.array[,,inds[1:jagsfit$m,t]], c(nsamp*nch,jagsfit$m))
jagsfit2$BUGSoutput$sims.list[["dmu"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "dmu")]][t]], c(nsamp*nch))
jagsfit2$BUGSoutput$sims.list[["dth"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "dth")]][t]], c(nsamp*nch))
jagsfit2$BUGSoutput$sims.list[["gmu"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "gmu")]][t]], c(nsamp*nch))
jagsfit2$BUGSoutput$sims.list[["gth"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "gth")]][t]], c(nsamp*nch))
jagsfit2$BUGSoutput$sims.list[["p"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "p")]][t]], c(nsamp*nch))
jagsfit2$BUGSoutput$sims.list[["pi"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "pi")]][t]], c(nsamp*nch))

jagsfit2$BUGSoutput$mean$z[,t] <- apply(jagsfit2$BUGSoutput$sims.list[["z"]][,,t],2,mean)
jagsfit2$BUGSoutput$mean$dmu[t] <- mean(jagsfit2$BUGSoutput$sims.list[["dmu"]][,t])
jagsfit2$BUGSoutput$mean$dth[t] <- mean(jagsfit2$BUGSoutput$sims.list[["dth"]][,t])
jagsfit2$BUGSoutput$mean$gmu[t] <- mean(jagsfit2$BUGSoutput$sims.list[["gmu"]][,t])
jagsfit2$BUGSoutput$mean$gth[t] <- mean(jagsfit2$BUGSoutput$sims.list[["gth"]][,t])
jagsfit2$BUGSoutput$mean$p[t] <- mean(jagsfit2$BUGSoutput$sims.list[["p"]][,t])
jagsfit2$BUGSoutput$mean$pi[t] <- mean(jagsfit2$BUGSoutput$sims.list[["pi"]][,t])
}

if(nch != 1){
jagsfit2$BUGSoutput$summary[,1] <- apply(apply(jagsfit2$BUGSoutput$sims.array,c(2,3),mean),2,mean)
jagsfit2$BUGSoutput$summary[,2] <- apply(apply(jagsfit2$BUGSoutput$sims.array,c(2,3),sd),2,sd)
jagsfit2$BUGSoutput$summary[,3:7] <- t(apply(jagsfit2$BUGSoutput$sims.matrix,2,function(x) quantile(x,probs=c(.025,.25,.50,.75,.975))))
jagsfit2$BUGSoutput$summary[,8] <- Rhat(jagsfit2$BUGSoutput$sims.array)
jagsfit2$BUGSoutput$summary[,8][is.nan(jagsfit2$BUGSoutput$summary[,8])] <- 1.000000

dimnames(jagsfit2$BUGSoutput$sims.array) <- list(NULL,NULL,rownames(jagsfit$BUGSoutput$summary))
dimnames(jagsfit2$BUGSoutput$sims.matrix) <- list(NULL,rownames(jagsfit$BUGSoutput$summary))
}else{

jagsfit2$BUGSoutput$summary[,1] <- apply(jagsfit2$BUGSoutput$sims.array,2,mean)
jagsfit2$BUGSoutput$summary[,2] <- apply(jagsfit2$BUGSoutput$sims.array,2,sd)
jagsfit2$BUGSoutput$summary[,3:7] <- t(apply(jagsfit2$BUGSoutput$sims.matrix,2,function(x) quantile(x,probs=c(.025,.25,.50,.75,.975))))
jagsfit2$BUGSoutput$summary <- jagsfit2$BUGSoutput$summary[,-c(8,9)]
dimnames(jagsfit2$BUGSoutput$sims.array) <- list(NULL,NULL,rownames(jagsfit$BUGSoutput$summary))
dimnames(jagsfit2$BUGSoutput$sims.matrix) <- list(NULL,rownames(jagsfit$BUGSoutput$summary))
}

#jagsfit.prealg <<- jagsfit
jagsfit <- jagsfit2; rm(jagsfit2)

return(jagsfit)
}

labelswitchalgltrm <- function(jagsfit,chnind=0){

jagsfit2 <- jagsfit

nch <- jagsfit2$BUGSoutput$n.chains
nsamp <- jagsfit2$BUGSoutput$n.keep

if(nch != 1){
ntruths <- jagsfit2$V
truths <- array(NA,c(jagsfit2$m,nch,ntruths))
inds <- jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "T")]]
inds <- matrix(inds,jagsfit2$m,jagsfit2$V)
for(t in 1:jagsfit$V){truths[,,t] <- t(apply(jagsfit$BUGSoutput$sims.array[ ,,inds[,t]],c(2,3),mean))}

T <- jagsfit2$V

chstart <- 1
if(length(chnind)==1){ 
chstart <- 2
chnind <- array(NA,c(T,nch))
chnind[1:T,1] <- 1:T

for(t in 1:T){
for(ch in chstart:nch){
Tind <- c(1:T)[-chnind[,ch][!is.na(chnind[,ch])]]
if(length(Tind)==0){Tind <- c(1:T)}

chnind[t,ch] <- which(max(cor(truths[1:jagsfit2$m, 1, t],truths[1:jagsfit2$m, ch, Tind]))==cor(truths[1:jagsfit2$m, 1, t],truths[1:jagsfit2$m, ch, ]))
#print(chnind); print(Tind)
}}
}

nsamp <- jagsfit$BUGSoutput$n.keep

inds <- rbind(inds,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "Tmu")]],
jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "Ttau")]],
matrix(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "gam")]],jagsfit2$C-1,jagsfit2$V),
jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "Emu")]],
jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "Etau")]],
jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "amu")]],
jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "atau")]],
jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "bmu")]],
jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "btau")]],
jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "pi")]])

einds <- jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]]
tmpeinds <- jagsfit$BUGSoutput$sims.array[ ,,einds]

for(t in 1:jagsfit$V){
for(ch in chstart:nch){
jagsfit2$BUGSoutput$sims.array[ ,ch,inds[,t]] <- jagsfit$BUGSoutput$sims.array[ ,ch,inds[,chnind[t,ch]]]  #this jagsfit2 vs. jagsfit difference is intentional
if(chstart==2){jagsfit2$BUGSoutput$sims.array[ ,ch,einds][jagsfit$BUGSoutput$sims.array[ ,ch,einds] == chnind[t,ch]] <- chnind[t,1]}
if(chstart==1){jagsfit2$BUGSoutput$sims.array[ ,ch,einds][jagsfit$BUGSoutput$sims.array[ ,ch,einds] == t] <- chnind[t,1]}
}}

}
#jagsfit2$BUGSoutput$sims.array[ ,,einds] <- tmpeinds
jagsfit2$BUGSoutput$sims.list[["e"]] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]]], c(nsamp*nch,jagsfit$n))
jagsfit2$BUGSoutput$mean$e <- apply(jagsfit2$BUGSoutput$sims.list[["e"]],2,mean)

jagsfit2$BUGSoutput$sims.matrix <- array(jagsfit2$BUGSoutput$sims.array,c(nsamp*nch,dim(jagsfit$BUGSoutput$sims.array)[3]))  #this jagsfit2 vs. jagsfit difference is intentional

jagsfit2$BUGSoutput$sims.list[["E"]] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "E")]]], c(nsamp*nch,jagsfit$n))
jagsfit2$BUGSoutput$sims.list[["a"]] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "a")]]], c(nsamp*nch,jagsfit$n))
jagsfit2$BUGSoutput$sims.list[["b"]] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "b")]]], c(nsamp*nch,jagsfit$n))
if(itemdiff == 1){jagsfit2$BUGSoutput$sims.list[["lam"]] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "lam")]]], c(nsamp*nch,jagsfit$m))}

jagsfit2$BUGSoutput$sims.list[["T"]] <- array(NA, c(nsamp*nch,jagsfit$m,jagsfit$V))
jagsfit2$BUGSoutput$sims.list[["Tmu"]] <- array(NA, c(nsamp*nch,jagsfit$V))
jagsfit2$BUGSoutput$sims.list[["Ttau"]] <- array(NA, c(nsamp*nch,jagsfit$V))
jagsfit2$BUGSoutput$sims.list[["gam"]] <- array(NA, c(nsamp*nch,jagsfit$C-1,jagsfit$V))
jagsfit2$BUGSoutput$sims.list[["Emu"]] <- array(NA, c(nsamp*nch,jagsfit$V))
jagsfit2$BUGSoutput$sims.list[["Etau"]] <- array(NA, c(nsamp*nch,jagsfit$V))
jagsfit2$BUGSoutput$sims.list[["amu"]] <- array(NA, c(nsamp*nch,jagsfit$V))
jagsfit2$BUGSoutput$sims.list[["atau"]] <- array(NA, c(nsamp*nch,jagsfit$V))
jagsfit2$BUGSoutput$sims.list[["bmu"]] <- array(NA, c(nsamp*nch,jagsfit$V))
jagsfit2$BUGSoutput$sims.list[["btau"]] <- array(NA, c(nsamp*nch,jagsfit$V))
jagsfit2$BUGSoutput$sims.list[["pi"]] <- array(NA, c(nsamp*nch,jagsfit$V))

if(jagsfit$C == 2 && length(dim(jagsfit2$BUGSoutput$sims.list[["gam"]])) < 3 ){
jagsfit2$BUGSoutput$sims.list[["gam"]] <- array(jagsfit2$BUGSoutput$sims.list[["gam"]], c(nsamp*nch,jagsfit$C-1,jagsfit$V))
}

#jagsfit.prealg <<- jagsfit2
for(t in 1:jagsfit2$V){
jagsfit2$BUGSoutput$sims.list[["T"]][,,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "T")]][1:jagsfit$m +((t-1)*jagsfit$m)]], c(nsamp*nch,jagsfit$m))
jagsfit2$BUGSoutput$sims.list[["Tmu"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "Tmu")]][t]], c(nsamp*nch))
jagsfit2$BUGSoutput$sims.list[["Ttau"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "Ttau")]][t]], c(nsamp*nch))
jagsfit2$BUGSoutput$sims.list[["gam"]][,,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "gam")]][1:(jagsfit$C-1) +((t-1)*(jagsfit$C-1))]], c(nsamp*nch,jagsfit$C-1))
jagsfit2$BUGSoutput$sims.list[["Emu"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "Emu")]][t]], c(nsamp*nch))
jagsfit2$BUGSoutput$sims.list[["Etau"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "Etau")]][t]], c(nsamp*nch))
jagsfit2$BUGSoutput$sims.list[["amu"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "amu")]][t]], c(nsamp*nch))
jagsfit2$BUGSoutput$sims.list[["atau"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "atau")]][t]], c(nsamp*nch))
jagsfit2$BUGSoutput$sims.list[["bmu"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "bmu")]][t]], c(nsamp*nch))
jagsfit2$BUGSoutput$sims.list[["btau"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "btau")]][t]], c(nsamp*nch))
jagsfit2$BUGSoutput$sims.list[["pi"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "pi")]][t]], c(nsamp*nch))

jagsfit2$BUGSoutput$mean$T[,t] <- apply(jagsfit2$BUGSoutput$sims.list[["T"]][,,t],2,mean)
jagsfit2$BUGSoutput$mean$Tmu[t] <- mean(jagsfit2$BUGSoutput$sims.list[["Tmu"]][,t])
jagsfit2$BUGSoutput$mean$Ttau[t] <- mean(jagsfit2$BUGSoutput$sims.list[["Ttau"]][,t])
if(jagsfit$C == 2){
jagsfit2$BUGSoutput$mean$gam[t] <- mean(jagsfit2$BUGSoutput$sims.list[["gam"]][,,t])
}else{jagsfit2$BUGSoutput$mean$gam[,t] <- apply(jagsfit2$BUGSoutput$sims.list[["gam"]][,,t],2,mean)}
jagsfit2$BUGSoutput$mean$Emu[t] <- mean(jagsfit2$BUGSoutput$sims.list[["Emu"]][,t])
jagsfit2$BUGSoutput$mean$Etau[t] <- mean(jagsfit2$BUGSoutput$sims.list[["Etau"]][,t])
jagsfit2$BUGSoutput$mean$amu[t] <- mean(jagsfit2$BUGSoutput$sims.list[["amu"]][,t])
jagsfit2$BUGSoutput$mean$atau[t] <- mean(jagsfit2$BUGSoutput$sims.list[["atau"]][,t])
jagsfit2$BUGSoutput$mean$bmu[t] <- mean(jagsfit2$BUGSoutput$sims.list[["bmu"]][,t])
jagsfit2$BUGSoutput$mean$btau[t] <- mean(jagsfit2$BUGSoutput$sims.list[["btau"]][,t])
jagsfit2$BUGSoutput$mean$pi[t] <- mean(jagsfit2$BUGSoutput$sims.list[["pi"]][,t])
}

if(nch != 1){
jagsfit2$BUGSoutput$summary[,1] <- apply(apply(jagsfit2$BUGSoutput$sims.array,c(2,3),mean),2,mean)
jagsfit2$BUGSoutput$summary[,2] <- apply(apply(jagsfit2$BUGSoutput$sims.array,c(2,3),sd),2,sd)
jagsfit2$BUGSoutput$summary[,3:7] <- t(apply(jagsfit2$BUGSoutput$sims.matrix,2,function(x) quantile(x,probs=c(.025,.25,.50,.75,.975))))
jagsfit2$BUGSoutput$summary[,8] <- Rhat(jagsfit2$BUGSoutput$sims.array)
jagsfit2$BUGSoutput$summary[,8][is.nan(jagsfit2$BUGSoutput$summary[,8])] <- 1.000000
dimnames(jagsfit2$BUGSoutput$sims.array) <- list(NULL,NULL,rownames(jagsfit$BUGSoutput$summary))
dimnames(jagsfit2$BUGSoutput$sims.matrix) <- list(NULL,rownames(jagsfit$BUGSoutput$summary))
}else{

jagsfit2$BUGSoutput$summary[,1] <- apply(jagsfit2$BUGSoutput$sims.array,2,mean)
jagsfit2$BUGSoutput$summary[,2] <- apply(jagsfit2$BUGSoutput$sims.array,2,sd)
jagsfit2$BUGSoutput$summary[,3:7] <- t(apply(jagsfit2$BUGSoutput$sims.matrix,2,function(x) quantile(x,probs=c(.025,.25,.50,.75,.975))))
jagsfit2$BUGSoutput$summary <- jagsfit2$BUGSoutput$summary[,-c(8,9)]
dimnames(jagsfit2$BUGSoutput$sims.array) <- list(NULL,NULL,rownames(jagsfit$BUGSoutput$summary))
dimnames(jagsfit2$BUGSoutput$sims.matrix) <- list(NULL,rownames(jagsfit$BUGSoutput$summary))
}

#jagsfit.prealg <<- jagsfit
jagsfit <- jagsfit2; rm(jagsfit2)

return(jagsfit)
}


labelswitchalgcrm <- function(jagsfit,chnind=0){

jagsfit2 <- jagsfit

nch <- jagsfit2$BUGSoutput$n.chains
nsamp <- jagsfit2$BUGSoutput$n.keep

if(nch != 1){
ntruths <- jagsfit2$V
truths <- array(NA,c(jagsfit2$m,nch,ntruths))
inds <- jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "T")]]
inds <- matrix(inds,jagsfit2$m,jagsfit2$V)
for(t in 1:jagsfit$V){truths[,,t] <- t(apply(jagsfit$BUGSoutput$sims.array[ ,,inds[,t]],c(2,3),mean))}

T <- jagsfit2$V

chstart <- 1
if(length(chnind)==1){ 
chstart <- 2
chnind <- array(NA,c(T,nch))
chnind[1:jagsfit2$V,1] <- 1:jagsfit2$V

for(t in 1:jagsfit2$V){
for(ch in chstart:nch){
Tind <- c(1:jagsfit2$V)[-chnind[,ch][!is.na(chnind[,ch])]]
if(length(Tind)==0){Tind <- c(1:jagsfit2$V)}

chnind[t,ch] <- which(max(cor(truths[1:jagsfit2$m, 1, t],truths[1:jagsfit2$m, ch, Tind]))==cor(truths[1:jagsfit2$m, 1, t],truths[1:jagsfit2$m, ch, ]))
#print(chnind); print(Tind)
}}
}

nsamp <- jagsfit$BUGSoutput$n.keep


inds <- rbind(inds,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "Tmu")]],
jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "Ttau")]],
jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "Emu")]],
jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "Etau")]],
jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "amu")]],
jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "atau")]],
jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "bmu")]],
jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "btau")]],
jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "pi")]])

einds <- jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]]
tmpeinds <- jagsfit$BUGSoutput$sims.array[ ,,einds]

for(t in 1:jagsfit$V){
for(ch in chstart:nch){
jagsfit2$BUGSoutput$sims.array[ ,ch,inds[,t]] <- jagsfit$BUGSoutput$sims.array[ ,ch,inds[,chnind[t,ch]]]  #this jagsfit2 vs. jagsfit difference is intentional

if(chstart==2){jagsfit2$BUGSoutput$sims.array[ ,ch,einds][jagsfit$BUGSoutput$sims.array[ ,ch,einds] == chnind[t,ch]] <- chnind[t,1]}
if(chstart==1){jagsfit2$BUGSoutput$sims.array[ ,ch,einds][jagsfit$BUGSoutput$sims.array[ ,ch,einds] == t] <- chnind[t,1]}
}}

}

#jagsfit2$BUGSoutput$sims.array[ ,,einds] <- tmpeinds
jagsfit2$BUGSoutput$sims.list[["e"]] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]]], c(nsamp*nch,jagsfit$n))
jagsfit2$BUGSoutput$mean$e <- apply(jagsfit2$BUGSoutput$sims.list[["e"]],2,mean)

jagsfit2$BUGSoutput$sims.matrix <- array(jagsfit2$BUGSoutput$sims.array,c(nsamp*nch,dim(jagsfit$BUGSoutput$sims.array)[3]))  #this jagsfit2 vs. jagsfit difference is intentional

jagsfit2$BUGSoutput$sims.list[["E"]] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "E")]]], c(nsamp*nch,jagsfit$n))
jagsfit2$BUGSoutput$sims.list[["a"]] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "a")]]], c(nsamp*nch,jagsfit$n))
jagsfit2$BUGSoutput$sims.list[["b"]] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "b")]]], c(nsamp*nch,jagsfit$n))
if(itemdiff == 1){jagsfit2$BUGSoutput$sims.list[["lam"]] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "lam")]]], c(nsamp*nch,jagsfit$m))}

jagsfit2$BUGSoutput$sims.list[["T"]] <- array(NA, c(nsamp*nch,jagsfit$m,jagsfit$V))
jagsfit2$BUGSoutput$sims.list[["Tmu"]] <- array(NA, c(nsamp*nch,jagsfit$V))
jagsfit2$BUGSoutput$sims.list[["Ttau"]] <- array(NA, c(nsamp*nch,jagsfit$V))
jagsfit2$BUGSoutput$sims.list[["Emu"]] <- array(NA, c(nsamp*nch,jagsfit$V))
jagsfit2$BUGSoutput$sims.list[["Etau"]] <- array(NA, c(nsamp*nch,jagsfit$V))
jagsfit2$BUGSoutput$sims.list[["amu"]] <- array(NA, c(nsamp*nch,jagsfit$V))
jagsfit2$BUGSoutput$sims.list[["atau"]] <- array(NA, c(nsamp*nch,jagsfit$V))
jagsfit2$BUGSoutput$sims.list[["bmu"]] <- array(NA, c(nsamp*nch,jagsfit$V))
jagsfit2$BUGSoutput$sims.list[["btau"]] <- array(NA, c(nsamp*nch,jagsfit$V))
jagsfit2$BUGSoutput$sims.list[["pi"]] <- array(NA, c(nsamp*nch,jagsfit$V))

for(t in 1:jagsfit2$V){
jagsfit2$BUGSoutput$sims.list[["T"]][,,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "T")]][1:jagsfit$m +((t-1)*jagsfit$m)]], c(nsamp*nch,jagsfit$m))
jagsfit2$BUGSoutput$sims.list[["Tmu"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "Tmu")]][t]], c(nsamp*nch))
jagsfit2$BUGSoutput$sims.list[["Ttau"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "Ttau")]][t]], c(nsamp*nch))
jagsfit2$BUGSoutput$sims.list[["Emu"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "Emu")]][t]], c(nsamp*nch))
jagsfit2$BUGSoutput$sims.list[["Etau"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "Etau")]][t]], c(nsamp*nch))
jagsfit2$BUGSoutput$sims.list[["amu"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "amu")]][t]], c(nsamp*nch))
jagsfit2$BUGSoutput$sims.list[["atau"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "atau")]][t]], c(nsamp*nch))
jagsfit2$BUGSoutput$sims.list[["bmu"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "bmu")]][t]], c(nsamp*nch))
jagsfit2$BUGSoutput$sims.list[["btau"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "btau")]][t]], c(nsamp*nch))
jagsfit2$BUGSoutput$sims.list[["pi"]][,t] <- array(jagsfit2$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "pi")]][t]], c(nsamp*nch))

jagsfit2$BUGSoutput$mean$T[,t] <- apply(jagsfit2$BUGSoutput$sims.list[["T"]][,,t],2,mean)
jagsfit2$BUGSoutput$mean$Tmu[t] <- mean(jagsfit2$BUGSoutput$sims.list[["Tmu"]][,t])
jagsfit2$BUGSoutput$mean$Ttau[t] <- mean(jagsfit2$BUGSoutput$sims.list[["Ttau"]][,t])
jagsfit2$BUGSoutput$mean$Emu[t] <- mean(jagsfit2$BUGSoutput$sims.list[["Emu"]][,t])
jagsfit2$BUGSoutput$mean$Etau[t] <- mean(jagsfit2$BUGSoutput$sims.list[["Etau"]][,t])
jagsfit2$BUGSoutput$mean$amu[t] <- mean(jagsfit2$BUGSoutput$sims.list[["amu"]][,t])
jagsfit2$BUGSoutput$mean$atau[t] <- mean(jagsfit2$BUGSoutput$sims.list[["atau"]][,t])
jagsfit2$BUGSoutput$mean$bmu[t] <- mean(jagsfit2$BUGSoutput$sims.list[["bmu"]][,t])
jagsfit2$BUGSoutput$mean$btau[t] <- mean(jagsfit2$BUGSoutput$sims.list[["btau"]][,t])
jagsfit2$BUGSoutput$mean$pi[t] <- mean(jagsfit2$BUGSoutput$sims.list[["pi"]][,t])
}

if(nch != 1){
jagsfit2$BUGSoutput$summary[,1] <- apply(apply(jagsfit2$BUGSoutput$sims.array,c(2,3),mean),2,mean)
jagsfit2$BUGSoutput$summary[,2] <- apply(apply(jagsfit2$BUGSoutput$sims.array,c(2,3),sd),2,sd)
jagsfit2$BUGSoutput$summary[,3:7] <- t(apply(jagsfit2$BUGSoutput$sims.matrix,2,function(x) quantile(x,probs=c(.025,.25,.50,.75,.975))))
jagsfit2$BUGSoutput$summary[,8] <- Rhat(jagsfit2$BUGSoutput$sims.array)
jagsfit2$BUGSoutput$summary[,8][is.nan(jagsfit2$BUGSoutput$summary[,8])] <- 1.000000
dimnames(jagsfit2$BUGSoutput$sims.array) <- list(NULL,NULL,rownames(jagsfit$BUGSoutput$summary))
dimnames(jagsfit2$BUGSoutput$sims.matrix) <- list(NULL,rownames(jagsfit$BUGSoutput$summary))
}else{

jagsfit2$BUGSoutput$summary[,1] <- apply(jagsfit2$BUGSoutput$sims.array,2,mean)
jagsfit2$BUGSoutput$summary[,2] <- apply(jagsfit2$BUGSoutput$sims.array,2,sd)
jagsfit2$BUGSoutput$summary[,3:7] <- t(apply(jagsfit2$BUGSoutput$sims.matrix,2,function(x) quantile(x,probs=c(.025,.25,.50,.75,.975))))
jagsfit2$BUGSoutput$summary <- jagsfit2$BUGSoutput$summary[,-c(8,9)]
dimnames(jagsfit2$BUGSoutput$sims.array) <- list(NULL,NULL,rownames(jagsfit$BUGSoutput$summary))
dimnames(jagsfit2$BUGSoutput$sims.matrix) <- list(NULL,rownames(jagsfit$BUGSoutput$summary))
}

#jagsfit.prealg <<- jagsfit
jagsfit <- jagsfit2; rm(jagsfit2)

return(jagsfit)
}

message("\n ...Inference complete, data is saved as 'jagsfit' \n")
message(" ...Performing final calculations")
#message(paste("Number of Rhats above 1.05 : ",sum(jagsfit$BUGSoutput$summary[,8]>1.05),"/",length(jagsfit$BUGSoutput$summary[,8]),"\nNumber of Rhats above 1.10 : ",sum(jagsfit$BUGSoutput$summary[,8]>1.1),"/",length(jagsfit$BUGSoutput$summary[,8]),sep=""))
if(jagsfit$BUGSoutput$n.chains > 1){
jagsfit$BUGSoutput$summary[,8] <- Rhat(jagsfit$BUGSoutput$sims.array)
jagsfit$BUGSoutput$summary[,8][is.nan(jagsfit$BUGSoutput$summary[,8])] <- 1.000000
#message(paste("Number of Rhats above 1.10 : ",sum(jagsfit$BUGSoutput$summary[,8]>1.10),"/",length(jagsfit$BUGSoutput$summary[,8]),"\nNumber of Rhats above 1.05 : ",sum(jagsfit$BUGSoutput$summary[,8]>1.05),"/",length(jagsfit$BUGSoutput$summary[,8]),sep=""))
if(whmodel== "GCM"){
message(paste("\nFor Continuous Parameters"))
message(paste("Number of Rhats above 1.10 : ",sum(jagsfit$BUGSoutput$summary[-c(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]],jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "z")]]),8]>1.10),"/",length(jagsfit$BUGSoutput$summary[-c(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]],jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "z")]]),8]),"\nNumber of Rhats above 1.05 : ",sum(jagsfit$BUGSoutput$summary[-c(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]],jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "z")]]),8]>1.05),"/",length(jagsfit$BUGSoutput$summary[-c(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]],jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "z")]]),8]),sep=""))
}else{
message(paste("Number of Rhats above 1.10 : ",sum(jagsfit$BUGSoutput$summary[-c(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]]),8]>1.10),"/",length(jagsfit$BUGSoutput$summary[-c(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]]),8]),"\nNumber of Rhats above 1.05 : ",sum(jagsfit$BUGSoutput$summary[-c(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]]),8]>1.10),"/",length(jagsfit$BUGSoutput$summary[-c(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]]),8]),sep=""))
}
}

if(clusters>1 && jagsfit$BUGSoutput$n.chains == 1){
#Some note whether the chain contains fewer mixtures than requested

Mode <<- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }
tmp <- unique(apply(jagsfit$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]]],c(2),Mode))
#print(tmp)
if(length(tmp) == 1){
message(paste("\n ...This chain has ",length(tmp)," culture rather than the ",clusters," cultures requested", sep=""))
message(paste("\n ...Try running the inference again",sep="" ))
}
if(length(tmp) != 1 && length(tmp) < clusters){
message(paste("\n ...This chain has ",tmp," cultures rather than the ",clusters," cultures requested", sep=""))
message(paste("\n ...Try running the inference again",sep="" ))
}
}

if(clusters>1 && jagsfit$BUGSoutput$n.chains > 1){
message("\n ...More than 1 culture applied with more than 1 chain")

einds <- jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]]

Mode <<- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

tmp <- apply(apply(jagsfit$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]]],c(2,3),Mode),1,unique)
#print(tmp)

chntoremove <- NULL
if(is.list(tmp)){
for(i in 1:jagsfit$BUGSoutput$n.chains){
if(length(tmp[[i]]) != clusters){chntoremove <- c(chntoremove,i)}
}
}

if(length(dim(tmp))==2){
for(i in 1:jagsfit$BUGSoutput$n.chains){
if(length(tmp[,i]) != clusters){chntoremove <- c(chntoremove,i)}
}
}

if(is.null(dim(tmp)) && !is.list(tmp)){chntoremove <- c(1:jagsfit$BUGSoutput$n.chains)}


if(length(chntoremove)>0){
if(length(chntoremove) < jagsfit$BUGSoutput$n.chains){

if(length(chntoremove)==1){
message(paste("\n ...",length(chntoremove)," chain out of ",jagsfit$BUGSoutput$n.chains," had fewer than ",clusters," cultures requested",sep=""))
message(paste("\n ...", "removing the ", length(chntoremove)," chain", sep=""))
}else{
message(paste("\n ...",length(chntoremove)," chains out of ",jagsfit$BUGSoutput$n.chains," had fewer than ",clusters," cultures requested",sep=""))
message(paste("\n ...", "removing these ", length(chntoremove)," chains", sep=""))
}

jagsfit$BUGSoutput$n.chains <- jagsfit$BUGSoutput$n.chains-length(chntoremove)
jagsfit$BUGSoutput$n.chain <- jagsfit$BUGSoutput$n.chains
jagsfit$BUGSoutput$n.sims <- jagsfit$BUGSoutput$n.chains*jagsfit$BUGSoutput$n.keep
if(jagsfit$BUGSoutput$n.chain == 1){
jagsfit$BUGSoutput$sims.array <- array(jagsfit$BUGSoutput$sims.array[,-chntoremove,], c(dim(jagsfit$BUGSoutput$sims.array)[1],1,dim(jagsfit$BUGSoutput$sims.array)[3]))
}else{jagsfit$BUGSoutput$sims.array <- jagsfit$BUGSoutput$sims.array[,-chntoremove,]}

}else{
message(paste("\n ...All chains out of ",jagsfit$BUGSoutput$n.chains," had fewer than ",clusters," cultures requested", sep=""))
message(paste("\n ...Try running the inference again",sep="" ))
}
}

message("\n ...Computing the most-consistent labeling across chains")


if(whmodel=="GCM"){

jagsfit <- labelswitchalggcm(jagsfit)

Mode <<- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }
jagsfit$respmem <- apply(jagsfit$BUGSoutput$sims.list$e[,],2,Mode)
tmeans <- jagsfit$BUGSoutput$mean$z
tmeans[tmeans<.5] <- tmeans[tmeans<.5]+1
ind <- rank(apply(abs(1-tmeans),2,mean))

# if(sum(c(1:length(ind))-ind==0) != length(ind)){
# chnind <- array(NA,c(jagsfit$T,jagsfit$BUGSoutput$n.chains))
# for(t in 1:jagsfit$T){chnind[t,] <- ind[t]}
# jagsfit <- labelswitchalggcm(jagsfit,chnind) }
} 


if(whmodel=="LTRM"){

jagsfit <- labelswitchalgltrm(jagsfit)

Mode <<- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
jagsfit$respmem <- apply(jagsfit$BUGSoutput$sims.list$e[,],2,Mode)
tmeans <- jagsfit$BUGSoutput$mean$T
ind <- rank(-apply(tmeans,2,sd))

# if(sum(c(1:length(ind))-ind==0) != length(ind)){
# chnind <- array(NA,c(jagsfit$V,jagsfit$BUGSoutput$n.chains))
# for(t in 1:jagsfit$V){chnind[t,] <- ind[t]}
# jagsfit <- labelswitchalgltrm(jagsfit,chnind) }
}

if(whmodel=="CRM"){

jagsfit <- labelswitchalgcrm(jagsfit)

Mode <<- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
jagsfit$respmem <- apply(jagsfit$BUGSoutput$sims.list$e[,],2,Mode)
tmeans <- jagsfit$BUGSoutput$mean$T
#tmeans[tmeans<.5] <- tmeans[tmeans<.5]+1
#ind <- rank(apply(abs(1-tmeans),2,mean))
ind <- rank(-apply(tmeans,2,sd))

# if(sum(c(1:length(ind))-ind==0) != length(ind)){
# chnind <- array(NA,c(jagsfit$V,jagsfit$BUGSoutput$n.chains))
# for(t in 1:jagsfit$V){chnind[t,] <- ind[t]}
# jagsfit <- labelswitchalgcrm(jagsfit,chnind) }

}

if(jagsfit$BUGSoutput$n.chains > 1){
#message(paste("Number of Rhats above 1.10 : ",sum(jagsfit$BUGSoutput$summary[,8]>1.10),"/",length(jagsfit$BUGSoutput$summary[,8]),"\nNumber of Rhats above 1.05 : ",sum(jagsfit$BUGSoutput$summary[,8]>1.05),"/",length(jagsfit$BUGSoutput$summary[,8]),sep=""))
message(paste("\nFor Continuous Parameters:"))
if(whmodel== "GCM"){
message(paste("Number of Rhats above 1.10 : ",sum(jagsfit$BUGSoutput$summary[-c(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]],jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "z")]]),8]>1.10),"/",length(jagsfit$BUGSoutput$summary[-c(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]],jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "z")]]),8]),"\nNumber of Rhats above 1.05 : ",sum(jagsfit$BUGSoutput$summary[-c(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]],jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "z")]]),8]>1.05),"/",length(jagsfit$BUGSoutput$summary[-c(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]],jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "z")]]),8]),sep=""))
}else{
message(paste("Number of Rhats above 1.10 : ",sum(jagsfit$BUGSoutput$summary[-c(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]]),8]>1.10),"/",length(jagsfit$BUGSoutput$summary[-c(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]]),8]),"\nNumber of Rhats above 1.05 : ",sum(jagsfit$BUGSoutput$summary[-c(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]]),8]>1.10),"/",length(jagsfit$BUGSoutput$summary[-c(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]]),8]),sep=""))
}
}
if(jagsfit$BUGSoutput$n.chains > 1){
message(paste("\nFor Discrete Parameters:"))
message(paste("Type 'dtraceplot()' to see their trace plots"))
alltraceplot <<- function(){traceplot(jagsfit,mfrow=c(4,4))}
dtraceplot <<- function(){
if(whmodel == "GCM"){
if(jagsfit$T == 1){
jagsfit2 <- jagsfit
jagsfit2$BUGSoutput$sims.array <- jagsfit2$BUGSoutput$sims.array[,,c(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "z")]])]
traceplot(jagsfit2,mfrow=c(4,4))
}else{
jagsfit2 <- jagsfit
jagsfit2$BUGSoutput$sims.array <- jagsfit2$BUGSoutput$sims.array[,,c(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "z")]],jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]])]
traceplot(jagsfit2,mfrow=c(4,4))
}
}else{
if(jagsfit$V == 1){
message("\n There are no discrete nodes in this inference")
return()
}else{
jagsfit2 <- jagsfit
jagsfit2$BUGSoutput$sims.array <- jagsfit2$BUGSoutput$sims.array[,,c(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]])]
traceplot(jagsfit2,mfrow=c(4,4))
}
}
}
}

}else{

if(whmodel == "GCM" && jagsfit$BUGSoutput$n.chains > 1){
message(paste("\nFor Discrete Parameters:"))
message(paste("Type 'dtraceplot()' to see their trace plots"))
alltraceplot <<- function(){traceplot(jagsfit,mfrow=c(4,4))}
dtraceplot <<- function(){
if(whmodel == "GCM"){
if(jagsfit$T == 1){
jagsfit2 <- jagsfit
jagsfit2$BUGSoutput$sims.array <- jagsfit2$BUGSoutput$sims.array[,,c(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "z")]])]
traceplot(jagsfit2,mfrow=c(4,4))
}else{
jagsfit2 <- jagsfit
jagsfit2$BUGSoutput$sims.array <- jagsfit2$BUGSoutput$sims.array[,,c(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "z")]],jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]])]
traceplot(jagsfit2,mfrow=c(4,4))
}
}else{
if(jagsfit$V == 1){
message("\n There are no discrete nodes in this inference")
return()
}else{
jagsfit2 <- jagsfit
jagsfit2$BUGSoutput$sims.array <- jagsfit2$BUGSoutput$sims.array[,,c(jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "e")]])]
traceplot(jagsfit2,mfrow=c(4,4))
}
}
}
}

}

message("\n ...Calculating DIC")
#DIC Calculations
jagsfit$BUGSoutput$pD <- var(jagsfit$BUGSoutput$sims.list$deviance[,1])/2
if(whmodel=="GCM"){
if(clusters == 1){
jagsfit$T <- 1;
jagsfit$BUGSoutput$sims.list$e <- array(1,c(jagsfit$BUGSoutput$n.sims,jagsfit$n));
if(length(dim(jagsfit$BUGSoutput$sims.list$z)) < 3){
jagsfit$BUGSoutput$sims.list$z <- array(jagsfit$BUGSoutput$sims.list[["z"]][,],c(jagsfit$BUGSoutput$n.sims,jagsfit$m,1))}
}
if(itemdiff == 0){jagsfit$BUGSoutput$sims.list[["lam"]] <- array(.5,c(jagsfit$BUGSoutput$n.sims,jagsfit$m))
}
jagsfit$BUGSoutput$sims.list$ID <- array(-1,c(jagsfit$BUGSoutput$n.sims,jagsfit$n,jagsfit$m))
jagsfit$Lik  <- array(NA, c(jagsfit$n,jagsfit$m,jagsfit$BUGSoutput$n.sims));
for(samp in 1:jagsfit$BUGSoutput$n.sims){
jagsfit$BUGSoutput$sims.list[["ID"]][samp,,] <- t( (1-jagsfit$BUGSoutput$sims.list[["lam"]][samp,])%o%(jagsfit$BUGSoutput$sims.list[["D"]][samp,]) / 
( (1-jagsfit$BUGSoutput$sims.list[["lam"]][samp,])%o%(jagsfit$BUGSoutput$sims.list[["D"]][samp,]) + 
(jagsfit$BUGSoutput$sims.list[["lam"]][samp,])%o%((1-jagsfit$BUGSoutput$sims.list[["D"]][samp,])) ) )
jagsfit$Lik[,,samp] <- t( ( (t(jagsfit$BUGSoutput$sims.list[["ID"]][samp,,])+(t(1-jagsfit$BUGSoutput$sims.list[["ID"]][samp,,])%*%diag(jagsfit$BUGSoutput$sims.list[["g"]][samp,])))^(t(jagsfit$data[,])*jagsfit$BUGSoutput$sims.list[["z"]][samp,,jagsfit$BUGSoutput$sims.list[["e"]][samp,]]))*(t(1-jagsfit$BUGSoutput$sims.list[["ID"]][samp,,])%*%diag(jagsfit$BUGSoutput$sims.list[["g"]][samp,]))^(t(jagsfit$data[,])*(1-jagsfit$BUGSoutput$sims.list[["z"]][samp,,jagsfit$BUGSoutput$sims.list[["e"]][samp,]]))*(t(1-jagsfit$BUGSoutput$sims.list[["ID"]][samp,,])%*%diag(1-jagsfit$BUGSoutput$sims.list[["g"]][samp,]))^((1-t(jagsfit$data[,]))*jagsfit$BUGSoutput$sims.list[["z"]][samp,,jagsfit$BUGSoutput$sims.list[["e"]][samp,]])*(t(jagsfit$BUGSoutput$sims.list[["ID"]][samp,,])+t(1-jagsfit$BUGSoutput$sims.list[["ID"]][samp,,])%*%diag(1-jagsfit$BUGSoutput$sims.list[["g"]][samp,]))^((1-t(jagsfit$data[,]))*(1-jagsfit$BUGSoutput$sims.list[["z"]][samp,,jagsfit$BUGSoutput$sims.list[["e"]][samp,]]))  )
}
 
jagsfit$BUGSoutput$sims.list$deviance <- array(-2*apply(log(jagsfit$Lik),3,sum),c(jagsfit$BUGSoutput$n.sims,1))
jagsfit$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "deviance")]]] <- array(jagsfit$BUGSoutput$sims.list$deviance,c(jagsfit$BUGSoutput$n.keep,jagsfit$BUGSoutput$n.chains))
jagsfit$BUGSoutput$sims.matrix[,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "deviance")]]] <- jagsfit$BUGSoutput$sims.list$deviance[,1]
 
jagsfit$BUGSoutput$mean$deviance <- mean(jagsfit$BUGSoutput$sims.list$deviance) #Dbar, also known as deviance
jagsfit$BUGSoutput$pD <- var(jagsfit$BUGSoutput$sims.list$deviance[,1])/2  #pD, variance of the deviance, divided by 2
jagsfit$BUGSoutput$DIC <- jagsfit$BUGSoutput$mean$deviance + jagsfit$BUGSoutput$pD
}

if(whmodel=="LTRM"){
if(clusters == 1){
jagsfit$V <- 1;
jagsfit$BUGSoutput$sims.list$e <- array(1,c(jagsfit$BUGSoutput$n.sims,jagsfit$n));
if(length(dim(jagsfit$BUGSoutput$sims.list$T)) < 3){
jagsfit$BUGSoutput$sims.list$T <- array(jagsfit$BUGSoutput$sims.list[["T"]][,],c(jagsfit$BUGSoutput$n.sims,jagsfit$m,1))}
if(length(dim(jagsfit$BUGSoutput$sims.list$gam))<3){jagsfit$BUGSoutput$sims.list$gam <- array(jagsfit$BUGSoutput$sims.list[["gam"]][,],c(jagsfit$BUGSoutput$n.sims,jagsfit$C-1,1))}
}
if(itemdiff == 0){jagsfit$BUGSoutput$sims.list[["lam"]] <- array(1,c(jagsfit$BUGSoutput$n.sims,jagsfit$m))
}
jagsfit$BUGSoutput$sims.list$ID <- array(-1,c(jagsfit$BUGSoutput$n.sims,jagsfit$n,jagsfit$m))
jagsfit$ppdelta <- array(NA, c(jagsfit$n,jagsfit$C-1,jagsfit$BUGSoutput$n.sims))
jagsfit$ppdeltafull <- array(NA, c(jagsfit$n,jagsfit$C+1,jagsfit$BUGSoutput$n.sims))
jagsfit$ppdeltafull[,1,] <- -1000000; jagsfit$ppdeltafull[,jagsfit$C+1,] <- 1000000
jagsfit$Lik  <- array(NA, c(jagsfit$n,jagsfit$m,jagsfit$BUGSoutput$n.sims));

for(samp in 1:jagsfit$BUGSoutput$n.sims){
jagsfit$BUGSoutput$sims.list[["ID"]][samp,,] <- jagsfit$BUGSoutput$sims.list[["E"]][samp,]%o%(1/jagsfit$BUGSoutput$sims.list[["lam"]][samp,])
jagsfit$ppdeltafull[,2:(jagsfit$C),samp] <- jagsfit$BUGSoutput$sims.list[["a"]][samp,]*t(jagsfit$BUGSoutput$sims.list[["gam"]][samp,,jagsfit$BUGSoutput$sims.list[["e"]][samp,]])+jagsfit$BUGSoutput$sims.list[["b"]][samp,]

for(i in 1:jagsfit$n){
jagsfit$Lik[i,,samp] <-  pnorm(jagsfit$ppdeltafull[i,jagsfit$data[i,]+1,samp] ,jagsfit$BUGSoutput$sims.list[["T"]][samp,,jagsfit$BUGSoutput$sims.list[["e"]][samp,i]],jagsfit$BUGSoutput$sims.list[["ID"]][samp,i,]^-.5)-pnorm(jagsfit$ppdeltafull[i,jagsfit$data[i,],samp] ,jagsfit$BUGSoutput$sims.list[["T"]][samp,,jagsfit$BUGSoutput$sims.list[["e"]][samp,i]],jagsfit$BUGSoutput$sims.list[["ID"]][samp,i,]^-.5)
}
}

if(jagsfit$C==2){jagsfit$Lik[jagsfit$Lik == 0] <- 0.001}

jagsfit$BUGSoutput$sims.list$deviance <- array(-2*apply(log(jagsfit$Lik),3,sum),c(jagsfit$BUGSoutput$n.sims,1))
jagsfit$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "deviance")]]] <- array(jagsfit$BUGSoutput$sims.list$deviance,c(jagsfit$BUGSoutput$n.keep,jagsfit$BUGSoutput$n.chains))
jagsfit$BUGSoutput$sims.matrix[,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "deviance")]]] <- jagsfit$BUGSoutput$sims.list$deviance[,1]
 
jagsfit$BUGSoutput$mean$deviance <- mean(jagsfit$BUGSoutput$sims.list$deviance) #Dbar, also known as deviance
jagsfit$BUGSoutput$pD <- var(jagsfit$BUGSoutput$sims.list$deviance[,1])/2  #pD, variance of the deviance, divided by 2
jagsfit$BUGSoutput$DIC <- jagsfit$BUGSoutput$mean$deviance + jagsfit$BUGSoutput$pD
}

if(whmodel=="CRM"){
if(clusters == 1){
jagsfit$V <- 1;
jagsfit$BUGSoutput$sims.list$e <- array(1,c(jagsfit$BUGSoutput$n.sims,jagsfit$n));
if(length(dim(jagsfit$BUGSoutput$sims.list$T)) < 3){
jagsfit$BUGSoutput$sims.list$T <- array(jagsfit$BUGSoutput$sims.list[["T"]][,],c(jagsfit$BUGSoutput$n.sims,jagsfit$m,1))}
}
if(itemdiff == 0){jagsfit$BUGSoutput$sims.list[["lam"]] <- array(1,c(jagsfit$BUGSoutput$n.sims,jagsfit$m))
}
jagsfit$BUGSoutput$sims.list$ID <- array(-1,c(jagsfit$BUGSoutput$n.sims,jagsfit$n,jagsfit$m))
jagsfit$Lik  <- array(NA, c(jagsfit$n,jagsfit$m,jagsfit$BUGSoutput$n.sims));

for(samp in 1:jagsfit$BUGSoutput$n.sims){
jagsfit$BUGSoutput$sims.list[["ID"]][samp,,] <- jagsfit$BUGSoutput$sims.list[["E"]][samp,]%o%(1/jagsfit$BUGSoutput$sims.list[["lam"]][samp,])
jagsfit$Lik[,,samp] <-  dnorm(jagsfit$data,mean=(jagsfit$BUGSoutput$sims.list[["a"]][samp,]*t(jagsfit$BUGSoutput$sims.list[["T"]][samp,,jagsfit$BUGSoutput$sims.list[["e"]][samp,]]))+array(jagsfit$BUGSoutput$sims.list[["b"]][samp,],c(jagsfit$n,jagsfit$m)),sd=jagsfit$BUGSoutput$sims.list[["a"]][samp,]*(jagsfit$BUGSoutput$sims.list[["ID"]][samp,,]^-.5))
}

jagsfit$BUGSoutput$sims.list$deviance <- array(-2*apply(log(jagsfit$Lik),3,sum),c(jagsfit$BUGSoutput$n.sims,1))
jagsfit$BUGSoutput$sims.array[,,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "deviance")]]] <- array(jagsfit$BUGSoutput$sims.list$deviance,c(jagsfit$BUGSoutput$n.keep,jagsfit$BUGSoutput$n.chains))
jagsfit$BUGSoutput$sims.matrix[,jagsfit$BUGSoutput$long.short[[which(jagsfit$BUGSoutput$root.short == "deviance")]]] <- jagsfit$BUGSoutput$sims.list$deviance[,1]
 
jagsfit$BUGSoutput$mean$deviance <- mean(jagsfit$BUGSoutput$sims.list$deviance) #Dbar, also known as deviance
jagsfit$BUGSoutput$pD <- var(jagsfit$BUGSoutput$sims.list$deviance[,1])/2  #pD, variance of the deviance, divided by 2
jagsfit$BUGSoutput$DIC <- jagsfit$BUGSoutput$mean$deviance + jagsfit$BUGSoutput$pD
}

message(paste("DIC : ",round(jagsfit$BUGSoutput$DIC,2),"   pD : ",round(jagsfit$BUGSoutput$pD,2),sep=""))
#print(var(sort(jagsfit$BUGSoutput$sims.list$deviance))/2); print(mean(jagsfit$BUGSoutput$sims.list$deviance))

jagsfit <<- jagsfit

tkconfigure(plotresults.but, state="normal") #"normal" / "disabled"
tkconfigure(doppc.but, state="normal") #"normal" / "disabled"
tkconfigure(exportresults.but, state="normal") #"normal" / "disabled"


}

plotresultsfunc <- function(saveplots=0,savedir=0) {

hdi <- function(sampleVec, credMass=0.95){
sortedPts = sort(sampleVec)
ciIdxInc = floor(credMass * length( sortedPts) )
nCIs = length( sortedPts) - ciIdxInc
ciwidth = rep(0, nCIs)
for(i in 1:nCIs){
	ciwidth[i] = sortedPts[i + ciIdxInc] - sortedPts[i]}
HDImin = sortedPts[which.min(ciwidth)]
HDImax = sortedPts[which.min(ciwidth)+ciIdxInc]
HDIlim = c(HDImin,HDImax)
return(HDIlim) 
}

Mode <<- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
jagsfit$respmem <- apply(jagsfit$BUGSoutput$sims.list$e[,],2,Mode)
#jagsfit$BUGSoutput$mode$e <- jagsfit$respmem 


#if(saveplots==1){jpeg(file.path(getwd(),"CCTpack","CCTpackplot.jpg"),width = 6, height = 6, units = "in", pointsize = 12,quality=100,res=400)}
#if(saveplots==2){postscript(file=file.path(getwd(),"CCTpack","CCTpackplot.eps"), onefile=FALSE, horizontal=FALSE, width = 6, height = 6, paper="special", family="Times")}
if(saveplots==1){jpeg(file.path(gsub(".Rdata","results.jpg",savedir)),width = 6, height = 6, units = "in", pointsize = 12,quality=100,res=400)}
if(saveplots==2){postscript(file=file.path(gsub(".Rdata","results.eps",savedir)), onefile=FALSE, horizontal=FALSE, width = 6, height = 6, paper="special", family="Times")}

if(whmodel=="GCM"){
par(oma=c(0,0,0,0),mar=c(4,4,3,1),mgp=c(2.25,.75,0),mfrow=c(2,2))

sym <- c(21,22, 23, 24, 25, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14) # possible symbols: 

if(clusters == 1){
if(length(dim(jagsfit$BUGSoutput$sims.list[["z"]])) < 3){
jagsfit$BUGSoutput$sims.list[["z"]] <- array(jagsfit$BUGSoutput$sims.list[["z"]],c(jagsfit$BUGSoutput$n.sims,jagsfit$m,1))}
plot(1:jagsfit$m,rep(NA,jagsfit$m),main=expression(paste("Item Truth (",z[tk],")")),xlab="Item",ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=21,bg="white")
for(i in 1:dim(jagsfit$BUGSoutput$sims.list[["z"]])[3]){
points(apply(jagsfit$BUGSoutput$sims.list[["z"]][,,i],c(2),mean),xlab=expression(paste("Item Truth (",z[tk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg="white")
}}else{
plot(1:jagsfit$m,rep(NA,jagsfit$m),main=expression(paste("Item Truth (",z[tk],") Per Culture")),xlab="Item",ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=21,bg="white")
for(i in 1:dim(jagsfit$BUGSoutput$sims.list[["z"]])[3]){
points(apply(jagsfit$BUGSoutput$sims.list[["z"]][,,i],c(2),mean),xlab=expression(paste("Item Truth (",z[tk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg="white")
}
}

if(itemdiff==1){plot(jagsfit$BUGSoutput$mean$lam,main=expression(paste("Item Difficulty (",lambda[k],")")),xlab="Item",ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[1],bg="white")
hditmp <- apply(jagsfit$BUGSoutput$sims.list$lam,2,hdi)
segments(1:jagsfit$m,hditmp[1,], 1:jagsfit$m, hditmp[2,])
arrows(1:jagsfit$m,hditmp[1,], 1:jagsfit$m, hditmp[2,],code=3,angle=90,length=.025)
}else{
plot(c(1:jagsfit$m),rep(.5,jagsfit$m),main=expression(paste("Item Difficulty (",lambda[k],")")),xlab="Item",ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[1],bg="white",col="grey")
points(c(par("usr")[1],par("usr")[2]),c(par("usr")[3],par("usr")[4]),type="l",col="grey")
points(c(par("usr")[1],par("usr")[2]),c(par("usr")[4],par("usr")[3]),type="l",col="grey")
}

plot(jagsfit$BUGSoutput$mean$D,main=expression(paste("Respondent Competency (",D[i],")")),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[jagsfit$respmem],bg="white")
hditmp <- apply(jagsfit$BUGSoutput$sims.list$D,2,hdi)
segments(1:jagsfit$n,hditmp[1,], 1:jagsfit$n, hditmp[2,])
arrows(1:jagsfit$n,hditmp[1,], 1:jagsfit$n, hditmp[2,],code=3,angle=90,length=.025)

plot(jagsfit$BUGSoutput$mean$g,main=expression(paste("Respondent Guessing Bias (",g[i],")")),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[jagsfit$respmem],bg="white")
hditmp <- apply(jagsfit$BUGSoutput$sims.list$g,2,hdi)
segments(1:jagsfit$n,hditmp[1,], 1:jagsfit$n, hditmp[2,])
arrows(1:jagsfit$n,hditmp[1,], 1:jagsfit$n, hditmp[2,],code=3,angle=90,length=.025)
}

if(whmodel=="LTRM"){
#par(oma=c(0,0,0,0),mar=c(4,4,3,1),mgp=c(2.25,.75,0),mfrow=c(2,2))
invlogit <- function(x){x <- 1 / (1 + exp(-x)); return(x)}

#par(oma=c(0,0,0,0),mar=c(4,4,3,1),mgp=c(2.25,.75,0))
#layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5), 2, 6, byrow = TRUE),heights=c(1,1,1,1,1))

par(oma=c(0,0,0,0),mar=c(4,4,3,1),mgp=c(2.25,.75,0),mfrow=c(2,2))

sym <- c(21,22, 23, 24, 25, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14) # possible symbols: 

if(jagsfit$C==2){useinvlogit <- 1}else{useinvlogit <- 0}

if(useinvlogit == 1){
tmp1 <- jagsfit$BUGSoutput$sims.list[["T"]]; jagsfit$BUGSoutput$sims.list[["T"]] <- invlogit(jagsfit$BUGSoutput$sims.list[["T"]])
tmp2 <- jagsfit$BUGSoutput$sims.list[["gam"]]; jagsfit$BUGSoutput$sims.list[["gam"]] <- invlogit(jagsfit$BUGSoutput$sims.list[["gam"]])
tmp3 <- jagsfit$BUGSoutput$sims.list[["b"]]; jagsfit$BUGSoutput$sims.list[["b"]] <- invlogit(jagsfit$BUGSoutput$sims.list[["b"]])
}

if(clusters == 1){
newm <- ceiling(jagsfit$m*1.092)
if(length(dim(jagsfit$BUGSoutput$sims.list[["T"]])) <3){
jagsfit$BUGSoutput$sims.list[["T"]] <- array(jagsfit$BUGSoutput$sims.list[["T"]],c(jagsfit$BUGSoutput$n.sims,jagsfit$m,1))
}
if(length(dim(jagsfit$BUGSoutput$sims.list[["gam"]])) <3){
jagsfit$BUGSoutput$sims.list[["gam"]] <- array(jagsfit$BUGSoutput$sims.list[["gam"]],c(jagsfit$BUGSoutput$n.sims,jagsfit$C-1,1))
}
hditmp <- apply(jagsfit$BUGSoutput$sims.list[["T"]][,,1],2,hdi)
if(useinvlogit == 1){
plot(1:newm,rep(NA,newm),main=expression(paste("Item Truth (",T[k],") and Thresholds (",gamma[c],")")),xlab="Item",ylim=c(0,1),ylab="Posterior Mean Value",las=1,pch=21,bg="white",axes=FALSE)
}else{plot(1:newm,rep(NA,newm),main=expression(paste("Item Truth (",T[k],") and Thresholds (",gamma[c],")")),xlab="Item",ylim=c(min(hditmp,min(apply(jagsfit$BUGSoutput$sims.list[["gam"]],c(2,3),mean))),max(hditmp,max(apply(jagsfit$BUGSoutput$sims.list[["gam"]],c(2,3),mean)))),ylab="Posterior Mean Value",las=1,pch=21,bg="white",axes=FALSE)
}
for(i in 1:dim(jagsfit$BUGSoutput$sims.list[["T"]])[3]){
points(apply(jagsfit$BUGSoutput$sims.list[["T"]][,,i],c(2),mean),xlab=expression(paste("Item Truth (",T[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg="white")

if(jagsfit$C == 2){
text(newm,mean(jagsfit$BUGSoutput$sims.list[["gam"]][,,i]),labels=sapply(c(1:(jagsfit$C-1)),function(x) as.expression(substitute(list(gamma[x]),list(x=x)))))
segments(1,mean(jagsfit$BUGSoutput$sims.list[["gam"]][,,i]),max(jagsfit$m+1,newm-2),mean(jagsfit$BUGSoutput$sims.list[["gam"]][,,i]),lty=2);
}
else{
text(newm,apply(jagsfit$BUGSoutput$sims.list[["gam"]][,,i],c(2),mean),labels=sapply(c(1:(jagsfit$C-1)),function(x) as.expression(substitute(list(gamma[x]),list(x=x)))))
segments(1,apply(jagsfit$BUGSoutput$sims.list[["gam"]][,,i],c(2),mean),max(jagsfit$m+1,newm-2),apply(jagsfit$BUGSoutput$sims.list[["gam"]][,,i],c(2),mean),lty=2);
}
segments(1:jagsfit$m,hditmp[1,], 1:jagsfit$m, hditmp[2,])
arrows(1:jagsfit$m,hditmp[1,], 1:jagsfit$m, hditmp[2,],code=3,angle=90,length=.025)
box()
axis(2, labels = TRUE,las=1)
axis(side = 1,labels=TRUE)
axis(side = 1, at = newm, labels = expression(gamma[c]) )
}}else{

newm <- ceiling(jagsfit$m*1.14)

if(jagsfit$C == 2){

if(useinvlogit == 1){
plot(1:newm,rep(NA,newm),main=expression(paste("Item Truth (",T[vk],") and Thresholds (",gamma[vc],") Per Culture")),xlab="Item",ylim=c(0,1),ylab="Posterior Mean Value",las=1,pch=21,bg="white",axes=FALSE)
}else{plot(1:newm,rep(NA,newm),main=expression(paste("Item Truth (",T[vk],") and Thresholds (",gamma[vc],") Per Culture")),xlab="Item",ylim=c(min(min(apply(jagsfit$BUGSoutput$sims.list[["T"]],c(2,3),mean))
,min(apply(jagsfit$BUGSoutput$sims.list[["gam"]],c(3),mean))),max(max(apply(jagsfit$BUGSoutput$sims.list[["T"]],c(2,3),mean))
,max(apply(jagsfit$BUGSoutput$sims.list[["gam"]],c(3),mean)))),ylab="Posterior Mean Value",las=1,pch=21,bg="white",axes=FALSE)
}

for(i in 1:dim(jagsfit$BUGSoutput$sims.list[["T"]])[3]){
points(apply(jagsfit$BUGSoutput$sims.list[["T"]][,,i],c(2),mean),xlab=expression(paste("Item Truth (",T[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg="white")
text(jagsfit$m+(((newm-jagsfit$m)/dim(jagsfit$BUGSoutput$sims.list[["T"]])[3])*i),mean(jagsfit$BUGSoutput$sims.list[["gam"]][,,i]),labels=sapply(c(1:(jagsfit$C-1)),function(x) as.expression(substitute(list(gamma[x]),list(x=x)))))
#segments(1,apply(jagsfit$BUGSoutput$sims.list[["gam"]][,,i],c(2),mean),(jagsfit$m+jagsfit$m+(((newm-jagsfit$m)/dim(jagsfit$BUGSoutput$sims.list[["T"]])[3])*1))/2,apply(jagsfit$BUGSoutput$sims.list[["gam"]][,,i],c(2),mean),lty=2);
}
}
else{

if(useinvlogit == 1){
plot(1:newm,rep(NA,newm),main=expression(paste("Item Truth (",T[vk],") and Thresholds (",gamma[vc],") Per Culture")),xlab="Item",ylim=c(0,1),ylab="Posterior Mean Value",las=1,pch=21,bg="white",axes=FALSE)
}else{
plot(1:newm,rep(NA,newm),main=expression(paste("Item Truth (",T[vk],") and Thresholds (",gamma[vc],") Per Culture")),xlab="Item",ylim=c(min(min(apply(jagsfit$BUGSoutput$sims.list[["T"]],c(2,3),mean))
,min(apply(jagsfit$BUGSoutput$sims.list[["gam"]][,,],c(2,3),mean))),max(max(apply(jagsfit$BUGSoutput$sims.list[["T"]],c(2,3),mean))
,max(apply(jagsfit$BUGSoutput$sims.list[["gam"]][,,],c(2,3),mean)))),ylab="Posterior Mean Value",las=1,pch=21,bg="white",axes=FALSE)
}

for(i in 1:dim(jagsfit$BUGSoutput$sims.list[["T"]])[3]){
points(apply(jagsfit$BUGSoutput$sims.list[["T"]][,,i],c(2),mean),xlab=expression(paste("Item Truth (",T[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg="white")
text(jagsfit$m+(((newm-jagsfit$m)/dim(jagsfit$BUGSoutput$sims.list[["T"]])[3])*i),apply(jagsfit$BUGSoutput$sims.list[["gam"]][,,i],c(2),mean),labels=sapply(c(1:(jagsfit$C-1)),function(x) as.expression(substitute(list(gamma[x]),list(x=x)))))
#segments(1,apply(jagsfit$BUGSoutput$sims.list[["gam"]][,,i],c(2),mean),(jagsfit$m+jagsfit$m+(((newm-jagsfit$m)/dim(jagsfit$BUGSoutput$sims.list[["T"]])[3])*1))/2,apply(jagsfit$BUGSoutput$sims.list[["gam"]][,,i],c(2),mean),lty=2);
}
}
box()
axis(2, labels = TRUE,las=1)
axis(side=1,at=axTicks(side=1)[axTicks(side=1)<= jagsfit$m])
axis(side = 1, at = jagsfit$m+(((newm-jagsfit$m)/dim(jagsfit$BUGSoutput$sims.list[["T"]])[3])*(1:2)), 
labels = sapply(1:2,function(x) as.expression(substitute(list(gamma[x*c]),list(x=x)))) )
}

if(itemdiff==1){
hditmp <- apply(jagsfit$BUGSoutput$sims.list$lam,2,hdi)
plot(jagsfit$BUGSoutput$mean$lam,main=expression(paste("Item Difficulty (",lambda[k],")")),xlab="Item",ylab="Posterior Mean Value",las=1,ylim=c(min(hditmp),max(hditmp)),pch=sym[1],bg="white")
segments(1:jagsfit$m,hditmp[1,], 1:jagsfit$m, hditmp[2,])
arrows(1:jagsfit$m,hditmp[1,], 1:jagsfit$m, hditmp[2,],code=3,angle=90,length=.025)
}else{
plot(c(1:jagsfit$m),rep(.5,jagsfit$m),main=expression(paste("Item Difficulty (",lambda[k],")")),xlab="Item",ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[1],bg="white",col="grey")
points(c(par("usr")[1],par("usr")[2]),c(par("usr")[3],par("usr")[4]),type="l",col="grey")
points(c(par("usr")[1],par("usr")[2]),c(par("usr")[4],par("usr")[3]),type="l",col="grey")
}

hditmp <- apply(jagsfit$BUGSoutput$sims.list$E,2,hdi)
plot(jagsfit$BUGSoutput$mean$E,main=expression(paste("Respondent Competency (",E[i],")")),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(0,max(hditmp)),pch=sym[jagsfit$respmem],bg="white")
segments(1:jagsfit$n,hditmp[1,], 1:jagsfit$n, hditmp[2,])
arrows(1:jagsfit$n,hditmp[1,], 1:jagsfit$n, hditmp[2,],code=3,angle=90,length=.025)




#hditmp <- apply(jagsfit$BUGSoutput$sims.list$b,2,hdi)

if(useinvlogit == 1){
plot(invlogit(jagsfit$BUGSoutput$mean$b)-.5,main=expression(paste("Respondent Shift Bias (",b[i],")")),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(-.5,.5),pch=sym[jagsfit$respmem],bg="white")
hditmp <- apply(jagsfit$BUGSoutput$sims.list$b-.5,2,hdi) # because this is invlogit transformed still
segments(1:jagsfit$n,hditmp[1,], 1:jagsfit$n, hditmp[2,])
arrows(1:jagsfit$n,hditmp[1,], 1:jagsfit$n, hditmp[2,],code=3,angle=90,length=.025)
}else{
plot(-5000, -5000, main=expression(paste("Category Usage Bias (",a[i]," and ",b[i],")")),xlim=c(min(-1.6,-abs(jagsfit$BUGSoutput$mean$b)),max(1.6,jagsfit$BUGSoutput$mean$b)), ylim=c(min(-.6,jagsfit$BUGSoutput$mean$a),min(3,max(2.4,jagsfit$BUGSoutput$mean$a))),xlab=expression(paste("Respondent Shift Bias (",b[i],")")),
ylab=expression(paste("Respondent Scale Bias (",a[i],")")),las=1,pch=21)
points(cbind(jagsfit$BUGSoutput$mean$b,jagsfit$BUGSoutput$mean$a),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(0,max(hditmp)),pch=sym[jagsfit$respmem],bg="white")
segments(-1,1,1,1)
segments(0,0,0,2)
text(0,2.25,"Middle Categories",cex=.8)
text(0,-.25,"Outer Categories",cex=.8)
text(-1.3,1,"Left \n Categories",cex=.8)
text(1.3,1,"Right \n Categories",cex=.8)
# text(0,2.4,"Middle",cex=.8)
# text(0,-.4,"Outer",cex=.8)
# text(-1.3,1,"Left",cex=.8)
# text(1.3,1,"Right",cex=.8)
}

# hditmp <- apply(jagsfit$BUGSoutput$sims.list$a,2,hdi)
# if(jagsfit$C == 2){
# plot(jagsfit$BUGSoutput$mean$a,main=expression(paste("Respondent Scale Bias (",a[i],")")),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(0,2),pch=sym[jagsfit$respmem],bg="white")
# }else{
# plot(jagsfit$BUGSoutput$mean$a,main=expression(paste("Respondent Scale Bias (",a[i],")")),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(0,max(hditmp)),pch=sym[jagsfit$respmem],bg="white")
# segments(1:jagsfit$n,hditmp[1,], 1:jagsfit$n, hditmp[2,])
# arrows(1:jagsfit$n,hditmp[1,], 1:jagsfit$n, hditmp[2,],code=3,angle=90,length=.025)
# }


# hditmp <- apply(jagsfit$BUGSoutput$sims.list$b,2,hdi)

# if(useinvlogit == 1){
# plot(invlogit(jagsfit$BUGSoutput$mean$b)-.5,main=expression(paste("Respondent Shift Bias (",b[i],")")),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(-.5,.5),pch=sym[jagsfit$respmem],bg="white")
# hditmp <- apply(jagsfit$BUGSoutput$sims.list$b-.5,2,hdi) # because this is invlogit transformed still
# }else{
# plot(jagsfit$BUGSoutput$mean$b,main=expression(paste("Respondent Shift Bias (",b[i],")")),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(min(hditmp),max(hditmp)),pch=sym[jagsfit$respmem],bg="white")
# }
# segments(1:jagsfit$n,hditmp[1,], 1:jagsfit$n, hditmp[2,])
# arrows(1:jagsfit$n,hditmp[1,], 1:jagsfit$n, hditmp[2,],code=3,angle=90,length=.025)

if(useinvlogit == 1){
jagsfit$BUGSoutput$sims.list[["T"]] <- tmp1 
jagsfit$BUGSoutput$sims.list[["gam"]] <- tmp2
jagsfit$BUGSoutput$sims.list[["b"]] <- tmp3
}

}

if(whmodel=="CRM"){

#par(oma=c(0,0,0,0),mar=c(4,4,3,1),mgp=c(2.25,.75,0))
#layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5), 2, 6, byrow = TRUE),heights=c(1,1,1,1,1))

par(oma=c(0,0,0,0),mar=c(4,4,3,1),mgp=c(2.25,.75,0),mfrow=c(2,2))

sym <- c(21,22, 23, 24, 25, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14) # possible symbols: 

if(clusters == 1){
if(length(dim(jagsfit$BUGSoutput$sims.list[["T"]])) < 3){
jagsfit$BUGSoutput$sims.list[["T"]] <- array(jagsfit$BUGSoutput$sims.list[["T"]],c(jagsfit$BUGSoutput$n.sims,jagsfit$m,1))}
hditmp <- apply(jagsfit$BUGSoutput$sims.list[["T"]][,,1],2,hdi)
plot(1:jagsfit$m,rep(NA,jagsfit$m),main=expression(paste("Item Truth (",T[vk],")")),xlab="Item",ylim=c(min(hditmp),max(hditmp)),ylab="Posterior Mean Value",las=1,pch=21,bg="white")
segments(1:jagsfit$m,hditmp[1,], 1:jagsfit$m, hditmp[2,])
arrows(1:jagsfit$m,hditmp[1,], 1:jagsfit$m, hditmp[2,],code=3,angle=90,length=.025)
}else{plot(1:jagsfit$m,rep(NA,jagsfit$m),main=expression(paste("Item Truth (",T[vk],") Per Culture")),xlab="Item",ylim=c(min(apply(jagsfit$BUGSoutput$sims.list[["T"]],c(2,3),mean)),max(apply(jagsfit$BUGSoutput$sims.list[["T"]],c(2,3),mean))),ylab="Posterior Mean Value",las=1,pch=21,bg="white")}
for(i in 1:dim(jagsfit$BUGSoutput$sims.list[["T"]])[3]){
points(apply(jagsfit$BUGSoutput$sims.list[["T"]][,,i],c(2),mean),xlab=expression(paste("Item Truth (",T[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg="white")
}

if(itemdiff==1){
hditmp <- apply(jagsfit$BUGSoutput$sims.list$lam,2,hdi)
plot(jagsfit$BUGSoutput$mean$lam,main=expression(paste("Item Difficulty (",lambda[k],")")),xlab="Item",ylab="Posterior Mean Value",las=1,ylim=c(min(hditmp),max(hditmp)),pch=sym[1],bg="white")
segments(1:jagsfit$m,hditmp[1,], 1:jagsfit$m, hditmp[2,])
arrows(1:jagsfit$m,hditmp[1,], 1:jagsfit$m, hditmp[2,],code=3,angle=90,length=.025)
}else{
plot(c(1:jagsfit$m),rep(.5,jagsfit$m),main=expression(paste("Item Difficulty (",lambda[k],")")),xlab="Item",ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[1],bg="white",col="grey")
points(c(par("usr")[1],par("usr")[2]),c(par("usr")[3],par("usr")[4]),type="l",col="grey")
points(c(par("usr")[1],par("usr")[2]),c(par("usr")[4],par("usr")[3]),type="l",col="grey")
}

hditmp <- apply(jagsfit$BUGSoutput$sims.list$E,2,hdi)
plot(jagsfit$BUGSoutput$mean$E,main=expression(paste("Respondent Competency (",E[i],")")),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(0,max(hditmp)),pch=sym[jagsfit$respmem],bg="white")
segments(1:jagsfit$n,hditmp[1,], 1:jagsfit$n, hditmp[2,])
arrows(1:jagsfit$n,hditmp[1,], 1:jagsfit$n, hditmp[2,],code=3,angle=90,length=.025)

# hditmp <- apply(jagsfit$BUGSoutput$sims.list$a,2,hdi)
# plot(jagsfit$BUGSoutput$mean$a,main=expression(paste("Respondent Scale Bias (",a[i],")")),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(0,max(hditmp)),pch=sym[jagsfit$respmem],bg="white")
# segments(1:jagsfit$n,hditmp[1,], 1:jagsfit$n, hditmp[2,])
# arrows(1:jagsfit$n,hditmp[1,], 1:jagsfit$n, hditmp[2,],code=3,angle=90,length=.025)

# hditmp <- apply(jagsfit$BUGSoutput$sims.list$b,2,hdi)
# plot(jagsfit$BUGSoutput$mean$b,main=expression(paste("Respondent Shift Bias (",b[i],")")),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(min(hditmp),max(hditmp)),pch=sym[jagsfit$respmem],bg="white")
# segments(1:jagsfit$n,hditmp[1,], 1:jagsfit$n, hditmp[2,])
# arrows(1:jagsfit$n,hditmp[1,], 1:jagsfit$n, hditmp[2,],code=3,angle=90,length=.025)

plot(-5000, -5000, main=expression(paste("Appraisal Bias (",a[i]," and ",b[i],")")),xlim=c(min(-2.8,-abs(jagsfit$BUGSoutput$mean$b)),max(2.8,jagsfit$BUGSoutput$mean$b)), ylim=c(min(-.6,jagsfit$BUGSoutput$mean$a),min(3,max(2.4,jagsfit$BUGSoutput$mean$a))),xlab=expression(paste("Respondent Shift Bias (",b[i],")")),
ylab=expression(paste("Respondent Scale Bias (",a[i],")")),las=1,pch=21)
points(cbind(jagsfit$BUGSoutput$mean$b,jagsfit$BUGSoutput$mean$a),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(0,max(hditmp)),pch=sym[jagsfit$respmem],bg="white")
segments(-2,1,2,1)
segments(0,0,0,2)
text(0,2.25,"Middle Values",cex=.8)
text(0,-.25,"Outer Values",cex=.8)
text(-2.5,1,"Left \n Values",cex=.8)
text(2.5,1,"Right \n Values",cex=.8)
# text(0,2.4,"Middle",cex=.8)
# text(0,-.4,"Outer",cex=.8)
# text(-1.3,1,"Left",cex=.8)
# text(1.3,1,"Right",cex=.8)



}

if(saveplots==1 || saveplots==2){dev.off()}

saveplots <- 0
}

ppcfunc <- function(saveplots=0,savedir=0) {
if(exists("checksrunbefore") == FALSE){
message("\n ...One moment, calculating posterior predictive checks\n")
if(whmodel=="GCM"){
if(clusters == 1){
jagsfit$T <- 1;
jagsfit$BUGSoutput$sims.list$e <- array(1,c(jagsfit$BUGSoutput$n.sims,jagsfit$n));
if(length(dim(jagsfit$BUGSoutput$sims.list$z)) < 3){
jagsfit$BUGSoutput$sims.list$z <- array(jagsfit$BUGSoutput$sims.list[["z"]][,],c(jagsfit$BUGSoutput$n.sims,jagsfit$m,1))}
}
if(itemdiff == 0){
jagsfit$BUGSoutput$sims.list[["lam"]] <- array(.5,c(jagsfit$BUGSoutput$n.sims,jagsfit$m))
}
jagsfit$ppY <- array(NA, c(jagsfit$n,jagsfit$m,jagsfit$BUGSoutput$n.sims));
jagsfit$BUGSoutput$sims.list$ID <- array(-1,c(jagsfit$BUGSoutput$n.sims,jagsfit$n,jagsfit$m))
jagsfit$Lik  <- array(NA, c(jagsfit$n,jagsfit$m,jagsfit$BUGSoutput$n.sims));

for(samp in 1:jagsfit$BUGSoutput$n.sims){
jagsfit$BUGSoutput$sims.list[["ID"]][samp,,] <- t( (1-jagsfit$BUGSoutput$sims.list[["lam"]][samp,])%o%(jagsfit$BUGSoutput$sims.list[["D"]][samp,]) / 
( (1-jagsfit$BUGSoutput$sims.list[["lam"]][samp,])%o%(jagsfit$BUGSoutput$sims.list[["D"]][samp,]) + 
(jagsfit$BUGSoutput$sims.list[["lam"]][samp,])%o%((1-jagsfit$BUGSoutput$sims.list[["D"]][samp,])) ) )

jagsfit$ppY[,,samp] <- matrix(rbinom((jagsfit$n*jagsfit$m),1,
(t(jagsfit$BUGSoutput$sims.list[["z"]][samp,,jagsfit$BUGSoutput$sims.list[["e"]][samp,]])*jagsfit$BUGSoutput$sims.list[["ID"]][samp,,] ) +
 t(t(1-jagsfit$BUGSoutput$sims.list[["ID"]][samp,,])%*%diag(jagsfit$BUGSoutput$sims.list[["g"]][samp,])) ), jagsfit$n,jagsfit$m)


 }

#sumloglik <- apply(log(jagsfit$Lik),3,sum) # -2*sumloglik = the deviance
#jagsfit$BUGSoutput$mean$deviance <- mean(-2*sumloglik) #Dbar, also known as deviance
#jagsfit$BUGSoutput$pD <- var(-2*sumloglik)/2  #pD, variance of the deviance, divided by 2
#pD <- var(sumloglik)/2 
#jagsfit$BUGSoutput$DIC <- jagsfit$BUGSoutput$mean$deviance + jagsfit$BUGSoutput$pD

jagsfit$mean$ppY <- rowMeans(jagsfit$ppY[,,],dims=2)

leigv <- 12
options(warn=-3)
ind <- sample(jagsfit$BUGSoutput$n.sims,min(500,jagsfit$BUGSoutput$n.sims)) #500 is the number of samples
eigv <- matrix(-1,length(ind),leigv)
tmp <- apply(jagsfit$ppY[,,ind],c(1,3),function(x) t(x))
tmp2 <- apply(tmp[,,],3,function(x) cor(x))
tmp3 <- array(tmp2,c(jagsfit$n,jagsfit$n,jagsfit$BUGSoutput$n.sims))
for(i in 1:length(ind)){
suppressMessages(try(eigv[i,] <- fa(tmp3[,,i])$values[1:leigv],silent=TRUE))}
wch <- -which(eigv[,1] == -1); 
if(length(wch)==0){jagsfit$ppeig <-  eigv}else{jagsfit$ppeig <-  eigv[-which(eigv[,1] == -1),]}
jagsfit$dateig <-  suppressMessages(fa(cor(t(jagsfit$data)))$values[1:leigv])
options(warn=0)

varvec <- matrix(-1,length(ind),dim(jagsfit$ppY)[2])
vdi <- matrix(-1,dim(jagsfit$ppY)[3],1)
#for(i in 1:dim(jagsfit$ppY)[3]){
for(i in 1:length(ind)){
varvec[i,] <- apply(jagsfit$ppY[,,i],2,var)
}
vdi <- apply(varvec,1,var)
jagsfit$ppVDI <- vdi
jagsfit$datVDI <- var(apply(jagsfit$data,2,var))

rm(vdi, varvec)
rm(eigv,leigv,tmp,tmp2,tmp3,ind);

checksrunbefore <<- 1
jagsfit <<- jagsfit
}

if(whmodel=="LTRM"){
if(clusters == 1){
jagsfit$V <- 1;
jagsfit$BUGSoutput$sims.list$e <- array(1,c(jagsfit$BUGSoutput$n.sims,jagsfit$n));
if(length(dim(jagsfit$BUGSoutput$sims.list[["T"]])) < 3){
jagsfit$BUGSoutput$sims.list$T <- array(jagsfit$BUGSoutput$sims.list[["T"]][,],c(jagsfit$BUGSoutput$n.sims,jagsfit$m,1))}
if(length(dim(jagsfit$BUGSoutput$sims.list[["gam"]])) < 3){
jagsfit$BUGSoutput$sims.list$gam <- array(jagsfit$BUGSoutput$sims.list[["gam"]][,],c(jagsfit$BUGSoutput$n.sims,jagsfit$C-1,1))}
}
if(itemdiff == 0){
jagsfit$BUGSoutput$sims.list[["lam"]] <- array(1,c(jagsfit$BUGSoutput$n.sims,jagsfit$m))
}
jagsfit$ppY <- array(NA, c(jagsfit$n,jagsfit$m,jagsfit$BUGSoutput$n.sims));
jagsfit$ppX <- array(NA, c(jagsfit$n,jagsfit$m,jagsfit$BUGSoutput$n.sims));
jagsfit$ppdelta <- array(NA, c(jagsfit$n,jagsfit$C-1,jagsfit$BUGSoutput$n.sims))

jagsfit$BUGSoutput$sims.list$ID <- array(-1,c(jagsfit$BUGSoutput$n.sims,jagsfit$n,jagsfit$m))

for(samp in 1:jagsfit$BUGSoutput$n.sims){

jagsfit$BUGSoutput$sims.list[["ID"]][samp,,] <- jagsfit$BUGSoutput$sims.list[["E"]][samp,]%o%(1/jagsfit$BUGSoutput$sims.list[["lam"]][samp,])

jagsfit$ppX[,,samp] <- matrix(
rnorm((jagsfit$n*jagsfit$m),t(jagsfit$BUGSoutput$sims.list[["T"]][samp,,jagsfit$BUGSoutput$sims.list[["e"]][samp,]]),(jagsfit$BUGSoutput$sims.list[["ID"]][samp,,]^(-.5))), 
jagsfit$n,jagsfit$m)

jagsfit$ppdelta[,,samp] <- jagsfit$BUGSoutput$sims.list[["a"]][samp,]*t(jagsfit$BUGSoutput$sims.list[["gam"]][samp,,jagsfit$BUGSoutput$sims.list[["e"]][samp,]])+jagsfit$BUGSoutput$sims.list[["b"]][samp,]

rc <- which(jagsfit$ppX[,,samp] < jagsfit$ppdelta[,1,samp],arr.ind=TRUE)
jagsfit$ppY[cbind(rc,array(samp,dim(rc)[1]))] <- 1 
for(c in 1:(jagsfit$C-1)){
rc <- which(jagsfit$ppX[,,samp] > jagsfit$ppdelta[,c,samp],arr.ind=TRUE)
jagsfit$ppY[cbind(rc,array(samp,dim(rc)[1]))] <- (c+1) }

}
jagsfit$mean$ppY <- rowMeans(jagsfit$ppY[,,],dims=2)

leigv <- 12
eigv <- matrix(-1,dim(jagsfit$ppY)[3],leigv)
options(warn=-3)
ind <- sample(jagsfit$BUGSoutput$n.sims,min(250,jagsfit$BUGSoutput$n.sims)) #250 is the number of samples
tmp <- apply(jagsfit$ppY[,,ind],c(1,3),function(x) t(x))
tmp2 <- apply(tmp[,,],3,function(x) polychoric(x,polycor=TRUE)$rho)
tmp3 <- array(tmp2,c(jagsfit$n,jagsfit$n,length(ind))) 
for(i in 1:length(ind)){
suppressMessages(try(eigv[i,] <- fa(tmp3[,,i])$values[1:leigv],silent=TRUE))}
wch <- -which(eigv[,1] == -1); 
if(length(wch)==0){jagsfit$ppeig <-  eigv}else{jagsfit$ppeig <-  eigv[-which(eigv[,1] == -1),]}
jagsfit$dateig <-  suppressMessages(fa(cor(t(jagsfit$data)))$values[1:leigv])
options(warn=0)
rm(eigv,leigv,tmp,tmp2,tmp3,ind);

varvec <- matrix(-1,dim(jagsfit$ppY)[3],dim(jagsfit$ppY)[2])
vdi <- matrix(-1,dim(jagsfit$ppY)[3],1)
for(i in 1:dim(jagsfit$ppY)[3]){
varvec[i,] <- apply(jagsfit$ppY[,,i],2,var)
}
vdi <- apply(varvec,1,var)
jagsfit$ppVDI <- vdi
jagsfit$datVDI <- var(apply(jagsfit$data,2,var))
rm(vdi, varvec)

checksrunbefore <<- 1
jagsfit <<- jagsfit
}

if(whmodel=="CRM"){
if(clusters == 1){
jagsfit$V <- 1;
jagsfit$BUGSoutput$sims.list$e <- array(1,c(jagsfit$BUGSoutput$n.sims,jagsfit$n));
if(length(dim(jagsfit$BUGSoutput$sims.list[["T"]])) < 3){
jagsfit$BUGSoutput$sims.list$T <- array(jagsfit$BUGSoutput$sims.list[["T"]][,],c(jagsfit$BUGSoutput$n.sims,jagsfit$m,1))}
}
if(itemdiff == 0){
jagsfit$BUGSoutput$sims.list[["lam"]] <- array(1,c(jagsfit$BUGSoutput$n.sims,jagsfit$m))
}
jagsfit$ppY <- array(NA, c(jagsfit$n,jagsfit$m,jagsfit$BUGSoutput$n.sims));
jagsfit$BUGSoutput$sims.list$ID <- array(-1,c(jagsfit$BUGSoutput$n.sims,jagsfit$n,jagsfit$m))

for(samp in 1:jagsfit$BUGSoutput$n.sims){

jagsfit$BUGSoutput$sims.list[["ID"]][samp,,] <- jagsfit$BUGSoutput$sims.list[["E"]][samp,]%o%(1/jagsfit$BUGSoutput$sims.list[["lam"]][samp,])

jagsfit$ppY[,,samp] <- matrix(
rnorm((jagsfit$n*jagsfit$m),(jagsfit$BUGSoutput$sims.list[["a"]][samp,]*t(jagsfit$BUGSoutput$sims.list[["T"]][samp,,jagsfit$BUGSoutput$sims.list[["e"]][samp,]]))+matrix(rep(jagsfit$BUGSoutput$sims.list[["b"]][samp,],jagsfit$m),jagsfit$n,jagsfit$m),
(jagsfit$BUGSoutput$sims.list[["a"]][samp,]*(jagsfit$BUGSoutput$sims.list[["ID"]][samp,,]^(-.5)))), 
jagsfit$n,jagsfit$m)
}
jagsfit$mean$ppY <- rowMeans(jagsfit$ppY[,,],dims=2)

leigv <- 12
options(warn=-3)
ind <- sample(jagsfit$BUGSoutput$n.sims,min(500,jagsfit$BUGSoutput$n.sims)) #500 is the number of samples
eigv <- matrix(-1,length(ind),leigv)
tmp <- apply(jagsfit$ppY[,,ind],c(1,3),function(x) t(x))
tmp2 <- apply(tmp[,,],3,function(x) cor(x))
tmp3 <- array(tmp2,c(jagsfit$n,jagsfit$n,jagsfit$BUGSoutput$n.sims))
for(i in 1:length(ind)){
suppressMessages(try(eigv[i,] <- fa(tmp3[,,i])$values[1:leigv],silent=TRUE))}
wch <- -which(eigv[,1] == -1); 
if(length(wch)==0){jagsfit$ppeig <-  eigv}else{jagsfit$ppeig <-  eigv[-which(eigv[,1] == -1),]}
jagsfit$dateig <-  suppressMessages(fa(cor(t(jagsfit$data)))$values[1:leigv])
options(warn=0)

varvec <- matrix(-1,length(ind),dim(jagsfit$ppY)[2])
vdi <- matrix(-1,dim(jagsfit$ppY)[3],1)
#for(i in 1:dim(jagsfit$ppY)[3]){
for(i in 1:length(ind)){
varvec[i,] <- apply(jagsfit$ppY[,,i],2,var)
}
vdi <- apply(varvec,1,var)
jagsfit$ppVDI <- vdi
jagsfit$datVDI <- var(apply(jagsfit$data,2,var))

rm(vdi, varvec)
rm(eigv,leigv,tmp,tmp2,tmp3,ind);

checksrunbefore <<- 1
jagsfit <<- jagsfit
}

message("\n ...Posterior predictive checks complete\n")
}


#if(saveplots==1){jpeg(file.path(getwd(),"CCTpack","CCTpackPPC.jpg"),width = 6, height = 3, units = "in", pointsize = 12,quality=100,res=400)}
#if(saveplots==2){postscript(file=file.path(getwd(),"CCTpack","CCTpackPPC.eps"), onefile=FALSE, horizontal=FALSE, width = 6, height = 6, paper="special", family="Times")}
if(saveplots==1){jpeg(file.path(gsub(".Rdata","ppc.jpg",savedir)),width = 6, height = 3, units = "in", pointsize = 12,quality=100,res=400)}
if(saveplots==2){postscript(file=file.path(gsub(".Rdata","ppc.eps",savedir)), onefile=FALSE, horizontal=FALSE, width = 6, height = 6, paper="special", family="Times")}

par(oma=c(0,0,0,0),mar=c(4,4,3,1),mgp=c(2.25,.75,0),mfrow=c(1,2))

plot(jagsfit$ppeig[1,], main="Culture Number Check",xlim=c(1,min(jagsfit$n,8)),ylim=c(min(jagsfit$dateig[min(jagsfit$n,8)],jagsfit$ppeig[,min(jagsfit$n,8)]),ceiling(max(jagsfit$ppeig[,1],jagsfit$dateig[1]))),xlab="Eigenvalue",ylab="Value",las=1,pch=21,type="l",col="black")
for(i in 2:dim(jagsfit$ppeig)[1]){
points(jagsfit$ppeig[i,],col="grey",type="l") }
points(jagsfit$dateig,col="black",type="l")

tmpdist <-ecdf(jagsfit$ppVDI)

color <- "black"; color2 <- "black"
if(clusters>1){color<-"grey"; color2 <- "white"}
plot(density(jagsfit$ppVDI),main="Item Difficulty Check", xlab= "VDI", xlim=c(min(jagsfit$datVDI,quantile(tmpdist, .005)),max(jagsfit$datVDI,quantile(tmpdist, .995))), ylim=c(0,max(density(jagsfit$ppVDI)$y)), ylab="Value",las=1,col=color2)
segments(jagsfit$datVDI,0,jagsfit$datVDI,par("usr")[4],col=color2)
text(.7*par("usr")[2],.7*par("usr")[4],labels=paste(round(100*tmpdist(jagsfit$datVDI),digits=2)," percentile"),col=color2); rm(tmpdist)

if(clusters>1){
points(c(par("usr")[1],par("usr")[2]),c(par("usr")[3],par("usr")[4]),type="l",col="grey")
points(c(par("usr")[1],par("usr")[2]),c(par("usr")[4],par("usr")[3]),type="l",col="grey")
#text(.5*(par("usr")[1]+par("usr")[2]),.5*(par("usr")[3]+par("usr")[4]),labels=paste("Compare Runs with DIC"),col="black")
}
#text(.8*par("usr")[2],.0275*par("usr")[4],labels=paste("DIC: ",round(jagsfit$BUGSoutput$DIC,digits=1)),col="black")
rm(color,color2)

if(saveplots==1 || saveplots==2){dev.off()}
saveplots <- 0

}

exportfunc <- function() {

#savedir <- tclvalue(tkgetSaveFile())
savedir <- tclvalue(tkgetSaveFile(initialfile = "CCTpackdata.Rdata",filetypes = "{{Rdata Files} {.Rdata}}"))
if (!nchar(savedir)) {
message("\n ...Export cancelled\n")
return()} 

if (savedir == 0 && !file.exists("CCTpack")){dir.create(file.path(getwd(), "CCTpack"));
savedir <- file.path(getwd(),"CCTpack","CCTpackdata.Rdata")
} 

if(substr(savedir, nchar(savedir)-6+1, nchar(savedir)) != ".Rdata"){
savedir <- paste(savedir,".Rdata",sep="")
}

message("\n ...Exporting results \n")
#save(jagsfit,file=file.path(getwd(),) )
write.csv(jagsfit$BUGSoutput$summary,file.path(gsub(".Rdata","posterior.csv",savedir)))
save(jagsfit,file=file.path(savedir))

if(exists("checksrunbefore") == TRUE){
ppcfunc(saveplots=1,savedir); ppcfunc(saveplots=2,savedir)}
plotresultsfunc(saveplots=1,savedir); plotresultsfunc(saveplots=2,savedir);
screeplotfunc(saveplots=1,savedir); screeplotfunc(saveplots=2,savedir);

message("\n ...Export complete\n")

#if (!file.exists("CCTpack")){dir.create(file.path(getwd(), "CCTpack"))} 
# if(exists("checksrunbefore") == TRUE){
# ppcfunc(saveplots=1); ppcfunc(saveplots=2)}
# plotresultsfunc(saveplots=1); plotresultsfunc(saveplots=2);
# screeplotfunc(saveplots=1); screeplotfunc(saveplots=2);
# save(jagsfit,file="CCTpack/CCTpackdata.Rdata")
# message("\n ...Results saved to 'CCTguidata.Rdata' , type 'jagsfit' to see the data \n")
}

loaddata.but <- tkbutton(datframe, text="Load Data", command=loadfilefunc)
screeplot.but <- tkbutton(datframe, text="Scree Plot", command=screeplotfunc)
applymodel.but <- tkbutton(applyframe, text="Apply CCT Model", command=applymodelfunc)
plotresults.but <- tkbutton(resultsframe, text = "Plot Results", command = plotresultsfunc)
doppc.but <- tkbutton(resultsframe, text = "Run Checks", command = ppcfunc)
exportresults.but <- tkbutton(resultsframe, text = "Export Results", command = exportfunc)

#submit.but <- tkbutton(tt, text="submit", command=submitfunc)
quit.but <- tkbutton(settingsframe, text = "Quit", command = quitfunc)
#reset.but <- tkbutton(tt, text="Reset", command=resetfunc)

sameid <- tkradiobutton(applyframe)
diffid <- tkradiobutton(applyframe)
#sameid <- tkradiobutton(tt)
#diffid <- tkradiobutton(tt)
itemdiffvar <- tclVar("0")
tkconfigure(sameid,variable=itemdiffvar, value="0")
tkconfigure(diffid,variable=itemdiffvar, value="1")

##########################
#Grid Setup


datafiletxt <- tktext(tt,bg="white",width=20,height=1) #font="courier"
resptxt <- tktext(tt,bg="white",width=4,height=1) #font="courier"
itemtxt <- tktext(tt,bg="white",width=4,height=1) #font="courier"
dattypetxt <- tktext(tt,bg="white",width=11,height=1) #font="courier"
modeltxt <- tktext(tt,bg="white",width=4,height=1) #font="courier"


tkgrid(datframe,columnspan=3,row=1,column=0)
tkgrid(datframe2,columnspan=4,row=2,column=0)
tkgrid(datframe3,columnspan=4,row=3,column=0)
tkgrid(applyframe,columnspan=10,row=4,column=0)
tkgrid(resultsframe,columnspan=3,row=5)
tkgrid(settingsframe,columnspan=8,row=6)


tkgrid(tklabel(datframe,text="    "))
tkgrid(tklabel(datframe,text="Data Input"),columnspan=3, pady = 5) 

tkgrid(loaddata.but,datafiletxt,screeplot.but,pady= 10, padx= 10)

tkgrid(tklabel(datframe2,text="Number of Respondents"),resptxt,tklabel(datframe2,text="Number of Items"),itemtxt, padx = 2, pady = 5) 

tkgrid(tklabel(datframe3,text="Data Type Detected"),dattypetxt,tklabel(datframe3,text="CCT Model"),modeltxt, padx = 2, pady = 5) 


#wm geometry . 346x334

tkgrid(tklabel(applyframe,text="Model Application"),columnspan=8, pady = 5) 
tkgrid(tklabel(applyframe,text="Number of Cultures to Assume:"), cultures.entry, tklabel(applyframe,text="Item Difficulty:"),tklabel(applyframe,text="No"),sameid,tklabel(applyframe,text="Yes"),diffid,applymodel.but,pady= 10, padx= 2)

tkconfigure(screeplot.but, state="disabled") #"disabled"

tkgrid(tklabel(resultsframe,text="Application Results"),columnspan=3, pady = 5) 
#tkgrid(plotresults.but,doppc.but,exportresults.but,pady= 10, padx= 10)
tkgrid(doppc.but,plotresults.but,exportresults.but,pady= 10, padx= 10)
tkconfigure(applymodel.but, state="disabled") #"disabled"
tkconfigure(plotresults.but, state="disabled") #"disabled"
tkconfigure(doppc.but, state="disabled") #"disabled"
tkconfigure(exportresults.but, state="disabled") #"disabled"

tkgrid(tklabel(settingsframe,text="Sampler Settings (Optional)"),columnspan=8, pady = 5) 
tkgrid(tklabel(settingsframe,text="Samples"), samples.entry, tklabel(settingsframe,text="Chains"), chains.entry, tklabel(settingsframe,text="Burn-in"), burnin.entry, tklabel(settingsframe,text="Thinning"), thin.entry, pady= 10, padx= 2)

fileName <- "(load data file)"
tkinsert(datafiletxt,"end",fileName)
tkconfigure(datafiletxt, state="disabled") #"disabled"
tkconfigure(resptxt, state="disabled") #"disabled"
tkconfigure(itemtxt, state="disabled") #"disabled"
tkconfigure(dattypetxt, state="disabled") #"disabled"
tkconfigure(modeltxt, state="disabled") #"disabled"

tkgrid(tklabel(tt,text="    "))
tkgrid.columnconfigure(tt,0,weight=1) 
tkgrid.rowconfigure(tt,0,weight=1) 

tkwm.resizable(tt,0,0)

##########################
#Model Code
mcgcm <-
"model{
for (i in 1:n){
 for (k in 1:m){
  tau[i,k] <- (D[i]*(1-lam[k])) / ((D[i]*(1-lam[k]))+(lam[k]*(1-D[i]))) 
  pY[i,k] <- (tau[i,k]*z[k,e[i]]) +((1-tau[i,k])*g[i])
  Y[i,k] ~ dbern(pY[i,k]) }} 

for (i in 1:n){
 e[i] ~ dcat(pi) 
 D[i] ~ dbeta(dmu[e[i]]*dth[e[i]],(1-dmu[e[i]])*dth[e[i]])
 g[i] ~ dbeta(gmu[e[i]]*gth[e[i]],(1-gmu[e[i]])*gth[e[i]]) }

 for (k in 1:m){
  lam[k] <- .5
 for (t in 1:T){
  z[k,t] ~ dbern(p[t]) }}

#Hyper Parameters
 alpha <- 2
 gsmu <- 10
 gssig <- 10
 dsmu <- 10
 dssig <- 10
 pi[1:T] ~ ddirch(L)

for (t in 1:T){
 L[t] <- 1
 gmu[t] <- .5
 gth[t] ~ dgamma(pow(gsmu,2)/pow(gssig,2),gsmu/pow(gssig,2))
 dmu[t] ~ dbeta(alpha,alpha)
 dth[t] ~ dgamma(pow(dsmu,2)/pow(dssig,2),dsmu/pow(dssig,2))
 p[t] ~ dunif(0,1) }}"

mcgcmid <-
"model{
for (i in 1:n){
 for (k in 1:m){
  tau[i,k] <- (D[i]*(1-lam[k])) / ((D[i]*(1-lam[k]))+(lam[k]*(1-D[i]))) 
  pY[i,k] <- (tau[i,k]*z[k,e[i]]) +((1-tau[i,k])*g[i])
  Y[i,k] ~ dbern(pY[i,k]) }} 

for (i in 1:n){
 e[i] ~ dcat(pi) 
 D[i] ~ dbeta(dmu[e[i]]*dth[e[i]],(1-dmu[e[i]])*dth[e[i]])
 g[i] ~ dbeta(gmu[e[i]]*gth[e[i]],(1-gmu[e[i]])*gth[e[i]]) }

 for (k in 1:m){
  lam[k] ~ dbeta(lammu*lamth,(1-lammu)*lamth)
 for (t in 1:T){
  z[k,t] ~ dbern(p[t]) }}

#Hyper Parameters
 alpha <- 2
 gsmu <- 10
 gssig <- 10
 dsmu <- 10
 dssig <- 10
 lamsmu <- 10
 lamssig <- 10
 lammu <- .5
 lamth ~ dgamma(pow(lamsmu,2)/pow(lamssig,2),lamsmu/pow(lamssig,2))
 pi[1:T] ~ ddirch(L)

for (t in 1:T){
 L[t] <- 1
 gmu[t] <- .5
 gth[t] ~ dgamma(pow(gsmu,2)/pow(gssig,2),gsmu/pow(gssig,2))
 dmu[t] ~ dbeta(alpha,alpha)
 dth[t] ~ dgamma(pow(dsmu,2)/pow(dssig,2),dsmu/pow(dssig,2))
 p[t] ~ dunif(0,1) }}"

mcltrm <-
"model{
   for (i in 1:n){
      for (k in 1:m){  
	tau[i,k] <- E[i]

	pY[i,k,1] <- pnorm((a[i]*gam[1,e[i]]) + b[i],T[k,e[i]],tau[i,k])
	for (c in 2:(C-1)){pY[i,k,c] <- pnorm((a[i]*gam[c,e[i]]) + b[i],T[k,e[i]],tau[i,k]) - sum(pY[i,k,1:(c-1)])}
	pY[i,k,C] <- (1 - sum(pY[i,k,1:(C-1)]))

	Y[i,k] ~ dcat(pY[i,k,1:C]) }}

#Parameters
   for (i in 1:n){
      e[i] ~ dcat(pi) 
      E[i] ~ dgamma(pow(Emu[e[i]],2)*Etau[e[i]],Emu[e[i]]*Etau[e[i]])
      a[i] ~ dgamma(atau[e[i]],atau[e[i]])
      b[i] ~ dnorm(bmu[e[i]],btau[e[i]]) }

   for (k in 1:m){
    for (v in 1:V){
     T[k,v] ~ dnorm(Tmu[v],Ttau[v]) }}

   for (v in 1:V){
     gam[1:(C-1),v] <- sort(tgam2[1:(C-1),v])
     for (c in 1:(C-2)){tgam[c,v] ~ dnorm(0,.1)}
     tgam2[1:(C-2),v] <- tgam[1:(C-2),v]
     tgam2[C-1,v] <- -sum(tgam[1:(C-2),v]) }

     pi[1:V] ~ ddirch(L)

#Hyperparameters	
 for (v in 1:V){
  L[v] <- 1 
 #Tmu[v] <- 0
  Tmu[v] ~ dnorm(0,.001)
  Ttau[v] ~ dgamma(1,.1)
  Emu[v] ~ dgamma(1,1)
  Etau[v] ~ dgamma(1,1)
  amu[v] <- 1
  atau[v] ~ dgamma(4,4)
  bmu[v] <- 0
  btau[v] ~ dgamma(4,4)
}}"

mcltrmid <-
"model{
   for (i in 1:n){
      for (k in 1:m){  
	tau[i,k] <- E[i]/lam[k]

	pY[i,k,1] <- pnorm((a[i]*gam[1,e[i]]) + b[i],T[k,e[i]],tau[i,k])
	for (c in 2:(C-1)){pY[i,k,c] <- pnorm((a[i]*gam[c,e[i]]) + b[i],T[k,e[i]],tau[i,k]) - sum(pY[i,k,1:(c-1)])}
	pY[i,k,C] <- (1 - sum(pY[i,k,1:(C-1)]))

	Y[i,k] ~ dcat(pY[i,k,1:C]) }}

#Parameters
   for (i in 1:n){
      e[i] ~ dcat(pi) 
      E[i] ~ dgamma(pow(Emu[e[i]],2)*Etau[e[i]],Emu[e[i]]*Etau[e[i]])
      a[i] ~ dgamma(atau[e[i]],atau[e[i]])
      b[i] ~ dnorm(bmu[e[i]],btau[e[i]]) }

   for (k in 1:m){
     lam[k] ~ dgamma(lamtau,lamtau)  
    for (v in 1:V){
     T[k,v] ~ dnorm(Tmu[v],Ttau[v]) }}

   for (v in 1:V){
     gam[1:(C-1),v] <- sort(tgam2[1:(C-1),v])
     for (c in 1:(C-2)){tgam[c,v] ~ dnorm(0,.1)}
     tgam2[1:(C-2),v] <- tgam[1:(C-2),v]
     tgam2[C-1,v] <- -sum(tgam[1:(C-2),v]) }

     pi[1:V] ~ ddirch(L)

#Hyperparameters	
 for (v in 1:V){
  L[v] <- 1 
 #Tmu[v] <- 0
  Tmu[v] ~ dnorm(0,.001)
  Ttau[v] ~ dgamma(1,.1)
  Emu[v] ~ dgamma(1,1)
  Etau[v] ~ dgamma(1,1)
  amu[v] <- 1
  atau[v] ~ dgamma(4,4)
  bmu[v] <- 0
  btau[v] ~ dgamma(4,4)
 }
 lammu <- 1
 lamtau ~ dgamma(4,4) 
}"

mccrm <-
"model{
   for (i in 1:n){
      for (k in 1:m){  
	tau[i,k] <- pow( a[i]*pow(E[i],-.5),-2)
	Y[i,k] ~ dnorm((a[i]*T[k,e[i]])+b[i],tau[i,k]) }}

#Parameters
   for (i in 1:n){
      e[i] ~ dcat(pi) 
      E[i] ~ dgamma(pow(Emu[e[i]],2)*Etau[e[i]],Emu[e[i]]*Etau[e[i]])
      a[i] ~ dgamma(atau[e[i]],atau[e[i]])
      b[i] ~ dnorm(bmu[e[i]],btau[e[i]]) }

   for (k in 1:m){
    for (v in 1:V){
     T[k,v] ~ dnorm(Tmu[v],Ttau[v]) }}

     pi[1:V] ~ ddirch(L)

#Hyperparameters	
 for (v in 1:V){
  L[v] <- 1 
 #Tmu[v] <- 0
 #Tmu[v] ~ dnorm(0,.001)
  Tmu[v] ~ dnorm(0,.1)
 #Ttau[v] ~ dgamma(0.01,0.01)
 Ttau[v] ~ dgamma(1,.1)
 Emu[v] ~ dgamma(1,1)
 Etau[v] ~ dgamma(1,1)
  #Emu[v] ~ dgamma(0.01,0.01)
  #Etau[v] ~ dgamma(0.01,0.01)
  amu[v] <- 1
 atau[v] ~ dgamma(4,4)
 # atau[v] ~ dgamma(0.01,0.01)
  bmu[v] <- 0
 btau[v] ~ dgamma(4,4)
 # btau[v] ~ dgamma(0.01,0.01)
}}"
mccrmid <-
"model{
   for (i in 1:n){
      for (k in 1:m){  
	tau[i,k] <- pow( a[i]*pow(E[i]/lam[k],-.5),-2)
	Y[i,k] ~ dnorm((a[i]*T[k,e[i]])+b[i],tau[i,k]) }}

#Parameters
   for (i in 1:n){
      e[i] ~ dcat(pi) 
      E[i] ~ dgamma(pow(Emu[e[i]],2)*Etau[e[i]],Emu[e[i]]*Etau[e[i]])
      a[i] ~ dgamma(atau[e[i]],atau[e[i]])
      b[i] ~ dnorm(bmu[e[i]],btau[e[i]]) }

   for (k in 1:m){
     lam[k] ~ dgamma(lamtau,lamtau)  
    for (v in 1:V){
     T[k,v] ~ dnorm(Tmu[v],Ttau[v]) }}

     pi[1:V] ~ ddirch(L)

#Hyperparameters	
 for (v in 1:V){
  L[v] <- 1 
 #Tmu[v] <- 0
  Tmu[v] ~ dnorm(0,.001)
  Ttau[v] ~ dgamma(0.01,0.01)
 #Ttau[v] ~ dgamma(1,.1)
 #Emu[v] ~ dgamma(1,1)
 #Etau[v] ~ dgamma(1,1)
  Emu[v] ~ dgamma(0.01,0.01)
  Etau[v] ~ dgamma(0.01,0.01)
  amu[v] <- 1
  #atau[v] ~ dgamma(4,4)
  atau[v] ~ dgamma(0.01,0.01)
  bmu[v] <- 0
  #btau[v] ~ dgamma(4,4)
  btau[v] ~ dgamma(0.01,0.01)
 }
 lammu <- 1
 #lamtau ~ dgamma(4,4)
 lamtau ~ dgamma(0.01,0.01) 
}"
} 
#cctgui()
####

