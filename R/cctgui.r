######################
#CCTpack GUI Application, begin with invoking the command cctgui()
######################

cctgui <- function(){
options(warn=-3)
suppressMessages(try(memory.limit(10000),silent=TRUE))
suppressMessages(try(memory.limit(20000),silent=TRUE))
options(warn=0)
Mode <- cctfit <- datob <- alltraceplot <- dtraceplot <- polywind <- cctscree <- cctapply <- cctresults <- cctppc <- cctexport <- NULL
rm(Mode,cctfit, datob, alltraceplot,dtraceplot,cctscree,cctapply,cctresults,cctppc,cctexport)
Mode <<- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }
polywind <- 0


######################
#Loads Packages, installs packages automatically if missing
######################
loadpkgs <- function(pkgs) { 
    pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
    if (length(pkgs_miss) > 0) {
	message("\n ...Some packages are missing, please follow the prompt to install them.\n")
    install.packages(pkgs_miss) }

    # install packages not already loaded:
    pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
    if (length(pkgs_miss) > 0) {install.packages(pkgs_miss) }
    
    # load packages not already loaded:
    attached <- search()
    attached_pkgs <- attached[grepl("package", attached)]
    need_to_attach <- pkgs[which(!pkgs %in% gsub("package:", "", attached_pkgs))]
    if (length(need_to_attach) > 0) {for (i in 1:length(need_to_attach)) require(need_to_attach[i], character.only = TRUE) }

	if (length(need_to_attach) == 0) {message("\n ...Starting CCT Inference Software\n")}
}

loadpkgs(c("tcltk", "psych","rjags","R2jags","mvtnorm","polycor"))

######################
#Sets Default GUI Variables and Frames
######################

samplesvar <- tclVar("10000") 
chainsvar <- tclVar("3")
burninvar <- tclVar("2000")
thinvar <- tclVar("1")
culturesvar <- tclVar("1") 

tt <- tktoplevel()
datframe <- tkframe(tt, borderwidth = 0)
datframe2 <- tkframe(tt, borderwidth = 0)
datframe3 <- tkframe(tt, borderwidth = 0)
datframe4 <- tkframe(tt, borderwidth = 0)
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

#######################
#Functions
######################



######################
#Function for the 'Load Data' Button
#- Detects the type of data
#- Detects the number of respondents and number of items
#- Selects the appropriate model for the data
#- Detects the number of missing data points (NA values)
#- Estimates initial estimates for missing data (if any) for the initial scree plot
#- Reports this information back to the GUI or R console
#- Enables the Apply Model Button on the GUI
#- Disables buttons GUI buttons: 'Run Checks', 'Plot Results', 'Export Results' 
#    if they are active from a previous run
######################
loadfilefuncbutton <- function(){
loadfilefunc(gui=1)
}

loadfilefunc <- function(data=0,gui=0,polych=0){
datob <- list()
datob$polych <- polych; datob$datfactorsp <- 1; datob$datfactorsc

if(gui==1){
datob$fileName <- file.path(tclvalue(tkgetOpenFile(filetypes = "{{csv Files} {.csv .txt}}")))
if (!nchar(datob$fileName)) {
return()
}else{
datob$dat <- as.matrix(read.csv(datob$fileName,header=FALSE))
options(warn=-3)
if(!is.numeric(datob$dat[1,1])){datob$dat <- matrix(as.numeric(datob$dat),dim(datob$dat)[1],dim(datob$dat)[2])
if(all(is.na(datob$dat[,1]))){datob$dat <- datob$dat[,-1];
}
if(all(is.na(datob$dat[1,]))){datob$dat <- datob$dat[-1,];
}
}

if(sum(is.na(datob$dat))>0){
datob$thena <- which(is.na(datob$dat),arr.ind=TRUE)
datob$thenalist <- which(is.na(datob$dat))
datob$dat[datob$thena] <- datob$dat[max(which(!is.na(datob$dat)),arr.ind=TRUE)]
datob$mval <- 1
}else{datob$mval <- 0}

options(warn=0)
if(all(datob$dat[1:dim(datob$dat)[1],1] == 1:dim(datob$dat)[1])){datob$dat <- datob$dat[,-1]; if(datob$mval==1){datob$thena[,2] <- datob$thena[,2] - 1}}

setwd(file.path(dirname(datob$fileName)))

tkconfigure(screeplot.but, state="normal") 
tkconfigure(applymodel.but, state="normal") 
}
tkconfigure(datafiletxt, state="normal") 
tkconfigure(resptxt, state="normal")
tkconfigure(itemtxt, state="normal") 
tkconfigure(dattypetxt, state="normal") 
tkconfigure(modeltxt, state="normal") 

tkdelete(datafiletxt,"1.0","800.0")
tkdelete(resptxt,"1.0","800.0")
tkdelete(itemtxt,"1.0","800.0")
tkdelete(dattypetxt,"1.0","800.0")
tkdelete(modeltxt,"1.0","800.0")

tkinsert(datafiletxt,"end",basename(datob$fileName))
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

if(!all(is.wholenumber(datob$dat))){ 
tkinsert(modeltxt,"end","CRM")
datob$whmodel <- "CRM"
if(min(datob$dat) >= 0 && max(datob$dat) <= 1){
message("\n ...Continuous data detected")
tkinsert(dattypetxt,"end","Continuous in [0,1]")
if(datob$mval==1){
datob$dat[datob$thena] <- NA
datob$thenalist <- which(is.na(datob$dat)) 
datob$dat[datob$thena] <- colMeans(datob$dat[,datob$thena[,2]],na.rm=TRUE)
}
datob$dat[datob$dat==0] <- .01; datob$dat[datob$dat==1] <- .99
datob$dat <- logit(datob$dat) 
}else{tkinsert(dattypetxt,"end","Continuous");
if(datob$mval==1){
datob$dat[datob$thena] <- NA
datob$thenalist <- which(is.na(datob$dat)) 
datob$dat[datob$thena] <- colMeans(datob$dat[,datob$thena[,2]],na.rm=TRUE)
}
message("\n ...Continuous data detected")}     
}else{if(min(datob$dat) >= 0 && max(datob$dat) <= 1){  
tkinsert(modeltxt,"end","GCM") 
datob$whmodel <- "GCM"
tkinsert(dattypetxt,"end","Binary")
message("\n ...Binary (dichotomous) data detected")
if(datob$mval==1){
datob$dat[datob$thena] <- NA
datob$thenalist <- which(is.na(datob$dat)) 
for(i in unique(datob$thena[,2])){
datob$dat[,i][is.na(datob$dat[,i])] <- Mode(datob$dat[,i][which(!is.na(datob$dat[,i]))]) 
}
}

}else{ 
tkinsert(modeltxt,"end","LTRM") 
datob$whmodel <- "LTRM"
tkinsert(dattypetxt,"end","Ordinal")
message("\n ...Ordinal (categorical) data detected")
if(polywind==0){
tkgrid(tklabel(datframe4,text="Use Polychoric Correlations:"),tklabel(datframe4,text="Yes"),polyyes,tklabel(datframe4,text="No"),polyno,pady= 10, padx= 2)
polywind <<- 1
}


if(datob$mval==1){
datob$dat[datob$thena] <- NA
datob$thenalist <- which(is.na(datob$dat)) 
for(i in unique(datob$thena[,2])){
datob$dat[,i][is.na(datob$dat[,i])] <- Mode(datob$dat[,i][which(!is.na(datob$dat[,i]))]) 
}
}
}

}

tkinsert(resptxt,"end",dim(datob$dat)[1])
tkinsert(itemtxt,"end",dim(datob$dat)[2])

tkconfigure(datafiletxt, state="disabled") 
tkconfigure(resptxt, state="disabled") 
tkconfigure(itemtxt, state="disabled") 
tkconfigure(dattypetxt, state="disabled") 
tkconfigure(modeltxt, state="disabled") 

datob$datind <- cbind(expand.grid(t(row(datob$dat))),expand.grid(t(col(datob$dat))),expand.grid(t(datob$dat)))
if(datob$mval==1){message("\n ...Data has ",dim(datob$thena)[1]," missing values out of ",length(datob$dat))
datob$datind <- datob$datind[-datob$thenalist,]
datob$datna <- datob$dat 
datob$datna[datob$thena] <- NA
}

message("\n ...Data loaded")

if(datob$whmodel=="LTRM"){	 
tkconfigure(polyyes, state="normal") 
tkconfigure(polyno, state="normal") 
}else{
if(polywind==1){
tkconfigure(polyyes, state="disabled") 
tkconfigure(polyno, state="disabled") 
} }			 

datob <<- datob
}
if(gui==0){

datob$dat <- data
options(warn=-3)
if(!is.numeric(datob$dat[1,1])){datob$dat <- matrix(as.numeric(datob$dat),dim(datob$dat)[1],dim(datob$dat)[2])
if(all(is.na(datob$dat[,1]))){datob$dat <- datob$dat[,-1];
}
if(all(is.na(datob$dat[1,]))){datob$dat <- datob$dat[-1,];
}
}

if(sum(is.na(datob$dat))>0){
datob$thena <- which(is.na(datob$dat),arr.ind=TRUE)
datob$thenalist <- which(is.na(datob$dat))
datob$dat[datob$thena] <- datob$dat[max(which(!is.na(datob$dat)),arr.ind=TRUE)]
datob$mval <- 1
}else{datob$mval <- 0}

options(warn=0)
if(all(datob$dat[1:dim(datob$dat)[1],1] == 1:dim(datob$dat)[1])){datob$dat <- datob$dat[,-1]; if(datob$mval==1){datob$thena[,2] <- datob$thena[,2] - 1}}
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

if(!all(is.wholenumber(datob$dat))){ 
datob$whmodel <- "CRM"
if(min(datob$dat) >= 0 && max(datob$dat) <= 1){
message("\n ...Continuous data detected")
if(datob$mval==1){
datob$dat[datob$thena] <- NA
datob$thenalist <- which(is.na(datob$dat)) 
datob$dat[datob$thena] <- colMeans(datob$dat[,datob$thena[,2]],na.rm=TRUE)
}
datob$dat[datob$dat==0] <- .01; datob$dat[datob$dat==1] <- .99
datob$dat <- logit(datob$dat) 
}else{
if(datob$mval==1){
datob$dat[datob$thena] <- NA
datob$thenalist <- which(is.na(datob$dat)) 
datob$dat[datob$thena] <- colMeans(datob$dat[,datob$thena[,2]],na.rm=TRUE)
}
message("\n ...Continuous data detected")}     
}else{if(min(datob$dat) >= 0 && max(datob$dat) <= 1){  
datob$whmodel <- "GCM"
message("\n ...Binary (dichotomous) data detected")
if(datob$mval==1){
datob$dat[datob$thena] <- NA
datob$thenalist <- which(is.na(datob$dat)) 
for(i in unique(datob$thena[,2])){
datob$dat[,i][is.na(datob$dat[,i])] <- Mode(datob$dat[,i][which(!is.na(datob$dat[,i]))]) 
}
}
}else{ 
datob$whmodel <- "LTRM"
message("\n ...Ordinal (categorical) data detected")

if(datob$mval==1){
datob$dat[datob$thena] <- NA
datob$thenalist <- which(is.na(datob$dat)) 
for(i in unique(datob$thena[,2])){
datob$dat[,i][is.na(datob$dat[,i])] <- Mode(datob$dat[,i][which(!is.na(datob$dat[,i]))]) 
}
}
}

}

message("... ",dim(datob$dat)[1]," respondents and ",dim(datob$dat)[2]," items")

datob$datind <- cbind(expand.grid(t(row(datob$dat))),expand.grid(t(col(datob$dat))),expand.grid(t(datob$dat)))
if(datob$mval==1){message("\n ...Data has ",dim(datob$thena)[1]," missing values out of ",length(datob$dat))
datob$datind <- datob$datind[-datob$thenalist,]
datob$datna <- datob$dat 
datob$datna[datob$thena] <- NA
}

return(datob)
}
			 
}

######################
#Function for the 'Scree Plot' Button
#- Uses the data object from loaddatafunc
#- Performs factor analysis on the Pearson correlations of the data 
#    using the fa() function from psych package
#- If polychoric correlations are picked for the LTRM, uses fa() on the polychoric 
#    correlations instead using polychoric() from the polycor package     
#- Creates a plot of 8 eigenvalues with appropriate titles and labels
#- When exporting the results, this function is run and produces .eps and .jpeg's of the plot
######################
screeplotfuncbutton <- function() {
screeplotfunc(datob=datob,saveplots=0,savedir=0,polych=as.numeric(tclvalue(polyvar)))
}
screeplotfunc <- function(datob,saveplots=0,savedir=0,gui=0,polych=0) {
options(warn=-3)

if(datob$whmodel != "LTRM"){
if(saveplots==0){
tmp <- ""
if(datob$mval==1 && datob$whmodel == "GCM"){tmp <- ", missing data handled by mode of respective columns"}
if(datob$mval==1 && datob$whmodel == "CRM"){tmp <- ", missing data handled by mean of respective columns"}
message("\n ...Producing Scree Plot",tmp)
}
if(saveplots==1){jpeg(file.path(gsub(".Rdata","scree.jpg",savedir)),width = 6, height = 6, units = "in", pointsize = 12,quality=100,res=400)}
if(saveplots==2){postscript(file=file.path(gsub(".Rdata","scree.eps",savedir)), onefile=FALSE, horizontal=FALSE, width = 6, height = 6, paper="special", family="Times")}
par(oma=c(0,0,0,0),mar=c(4,4,3,1),mgp=c(2.25,.75,0),mfrow=c(1,1))
suppressMessages(plot(fa(cor(t(datob$dat)))$values[1:8],las=1,type="b",bg="black",pch=21,xlab="Factor",ylab="Magnitude",main="Scree Plot of Data"))

if(saveplots==1 || saveplots ==2){dev.off()}
}else{
if(saveplots==0){
tmp <- ""
if(datob$mval==1 && datob$whmodel == "LTRM"){tmp <- ", missing data handled by mode of respective columns"}
message("\n ...Producing Scree Plot",tmp)
}
if(saveplots==1){jpeg(file.path(gsub(".Rdata","scree.jpg",savedir)),width = 6, height = 6, units = "in", pointsize = 12,quality=100,res=400)}
if(saveplots==2){postscript(file=file.path(gsub(".Rdata","scree.eps",savedir)), onefile=FALSE, horizontal=FALSE, width = 6, height = 6, paper="special", family="Times")}
if(polych == 1){
if(length(datob$datfactorsp)==1){
datob$datfactorsp <- suppressMessages(fa(polychoric(t(datob$dat),polycor=TRUE)$rho)$values[1:8])
}
datob$datfactors <- datob$datfactorsp
}else{
if(length(datob$datfactorsp)==1){
datob$datfactorsc <- suppressMessages(fa(cor(t(datob$dat)))$values[1:8])
}
datob$datfactors <- datob$datfactorsc
}
par(oma=c(0,0,0,0),mar=c(4,4,3,1),mgp=c(2.25,.75,0),mfrow=c(1,1))
plot(datob$datfactors,las=1,type="b",bg="black",pch=21,xlab="Factor",ylab="Magnitude",main="Scree Plot of Data")

if(gui==1){datob <<- datob}
if(saveplots==1 || saveplots ==2){dev.off()}
}

if(gui==0){return(datob)}
options(warn=0)
}

######################
#Function for the 'Apply CCT Model' Button
#- Uses the data object from loaddatafunc
#- Applies hierarchical Bayesian inference for the model using JAGS by packages rjags and R2jags
#- Model code for each of the three models are at the bottom of this rcode file
#- Reads in user preferences for model specifications from the GUI
#- Uses these specifications to apply different form(s) of the model
#- During mixture model cases, applies an algorithm that corrects for label-switching/mixing issues
#- Recalculates the statistics for all parameters after correcting for label-switching, as well as
#    the Rhat statistics, and DIC
#- Enables the user to look at traceplots for discrete nodes via dtraceplots() 
#- Provides the model-based clustering by cctfit$respmem   (respondent membership)
#- Provides cctfit$Lik, which is the likelihood of the model evaluated at each sample
#- Enables GUI buttons: 'Run Checks', 'Plot Results', 'Export Results'
######################
applymodelfuncbutton <- function() {
applymodelfunc(datob=datob,clusters=as.numeric(tclvalue(culturesvar)),itemdiff=as.numeric(tclvalue(itemdiffvar)),
jags.iter=as.numeric(tclvalue(samplesvar)),jags.chains= as.numeric(tclvalue(chainsvar)),jags.burnin=as.numeric(tclvalue(burninvar)),
jags.thin=as.numeric(tclvalue(thinvar)),gui=1)
}

applymodelfunc <- function(datob,clusters=1,itemdiff=0,jags.iter=10000,jags.chains=3,jags.burnin=2000,jags.thin=1,gui=0) {

######################
#Sets up Parameters, Variables, and Data for JAGS for the model selected
######################
if(datob$whmodel == "GCM"){
Y <- datob$datind; n <- dim(datob$dat)[1]; m <- dim(datob$dat)[2]; V <- clusters; nobs <- dim(datob$datind)[1]
jags.data <- list("Y","n","m","V","nobs")

if(itemdiff==0){
model.file <- mcgcm
jags.params <- c("Z","th","g","p","thmu","thtau","gmu","gtau","e","pi")
if(clusters>1){
jags.inits <- function(){ list("Z"=matrix(rbinom(m*V,1,.5),m,V),"th"= runif(n,.2,.8), "g"= runif(n,.2,.8),"e"= sample(1:V,n,replace=TRUE) )}
}else{
jags.inits <- function(){ list("Z"=matrix(rbinom(m*V,1,.5),m,V),"th"= runif(n,.2,.8), "g"= runif(n,.2,.8) )}
}

}
if(itemdiff==1){
model.file <- mcgcmid
jags.params <- c("Z","th","g","lam","p","thmu","thtau","gmu","gtau","lammu","lamtau","e","pi")
if(clusters>1){
jags.inits <- function(){ list("Z"=matrix(rbinom(m*V,1,.5),m,V),"th"= runif(n,.2,.8), "g"= runif(n,.2,.8), "lam"= matrix(runif(m*V,.2,.8),m,V), "e"= sample(1:V,n,replace=TRUE) )}
}else{
jags.inits <- function(){ list("Z"=matrix(rbinom(m*V,1,.5),m,V),"th"= runif(n,.2,.8), "g"= runif(n,.2,.8), "lam"= matrix(runif(m*V,.2,.8),m,V) )}
}
}
if(clusters==1){
model.file <- gsub(pattern="pi\\[1\\:V\\] ~ ddirch\\(L\\)", "pi <- 1", model.file)
model.file <- gsub(pattern="e\\[i\\] ~ dcat\\(pi\\)", "e\\[i\\] <- 1", model.file)
}
}

if(datob$whmodel == "LTRM"){
Y <- datob$datind; n <- dim(datob$dat)[1]; m <- dim(datob$dat)[2]; V <- clusters; C <- max(datob$datind[,3]); nobs <- dim(datob$datind)[1]
jags.data <- list("Y","n","m","C","V","nobs")

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
jags.inits <- function(){ list("T"=matrix(rnorm(m*V,0,1),m,V), "E"= runif(n,.8,1.2), "b"= rnorm(n,0,1), "e"= sample(1:V,n,replace=TRUE)  )} 
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
jags.inits <- function(){ list("T"=matrix(rnorm(m*V,0,1),m,V), "lam"= matrix(runif(m*V,.8,1.2),m,V), "tgam"=matrix(rnorm((C-2)*V,0,1),(C-2),V), "E"= runif(n,.8,1.2), "a"= runif(n,.8,1.2), "b"= rnorm(n,0,1), "e"= sample(1:V,n,replace=TRUE) )}
}else{
jags.inits <- function(){ list("T"=matrix(rnorm(m*V,0,1),m,V), "lam"= matrix(runif(m*V,.8,1.2),m,V), "tgam"=matrix(rnorm((C-2)*V,0,1),(C-2),V), "E"= runif(n,.8,1.2), "a"= runif(n,.8,1.2), "b"= rnorm(n,0,1) )}
}
if(C==2){
if(clusters>1){
jags.inits <- function(){ list("T"=matrix(rnorm(m*V,0,1),m,V), "lam"= matrix(runif(m*V,.8,1.2),m,V), "E"= runif(n,.8,1.2), "b"= rnorm(n,0,1), "e"= sample(1:V,n,replace=TRUE) )}
}else{
jags.inits <- function(){ list("T"=matrix(rnorm(m*V,0,1),m,V), "lam"= matrix(runif(m*V,.8,1.2),m,V), "E"= runif(n,.8,1.2), "b"= rnorm(n,0,1) )}
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

if(datob$whmodel == "CRM"){
Y <- datob$datind; n <- dim(datob$dat)[1]; m <- dim(datob$dat)[2]; V <- clusters; nobs <- dim(datob$datind)[1]
jags.data <- list("Y","n","m","V","nobs")

if(itemdiff==0 ){
model.file <- mccrm
jags.params <- c("T","E","a","b","Tmu","Ttau","Emu","Etau","amu","atau","bmu","btau","e","pi")
if(clusters>1){
jags.inits <- function(){ list("T"=matrix(rnorm(m*V,0,1),m,V),"E"= runif(n,.8,1.2), "a"= runif(n,.8,1.2), "b"= rnorm(n,0,1), "e"= sample(1:V,n,replace=TRUE)   )}
}else{
jags.inits <- function(){ list("T"=matrix(rnorm(m*V,0,1),m,V),"E"= runif(n,.8,1.2), "a"= runif(n,.8,1.2), "b"= rnorm(n,0,1)  )}
}
}

if(itemdiff==1){
model.file <- mccrmid
jags.params <- c("T","E","a","b","lam","Tmu","Ttau","Emu","Etau","amu","atau","bmu","btau","lammu","lamtau","e","pi")
if(clusters>1){
jags.inits <- function(){ list("T"=matrix(rnorm(m*V,0,1),m,V),"E"= runif(n,.8,1.2), "a"= runif(n,.8,1.2), "b"= rnorm(n,0,1), "lam"= matrix(runif(m*V,.8,1.2),m,V), "e"= sample(1:V,n,replace=TRUE)   )}
}else{
jags.inits <- function(){ list("T"=matrix(rnorm(m*V,0,1),m,V),"E"= runif(n,.8,1.2), "a"= runif(n,.8,1.2), "b"= rnorm(n,0,1), "lam"= matrix(runif(m*V,.8,1.2),m,V)  )}
}
}
if(clusters==1){ 
model.file <- gsub(pattern="pi\\[1\\:V\\] ~ ddirch\\(L\\)", "pi <- 1", model.file)
model.file <- gsub(pattern="e\\[i\\] ~ dcat\\(pi\\)", "e\\[i\\] <- 1", model.file)
} 
}

######################
#Runs the Model in JAGS
#- saves the data object from loadfilefunc within the jags object
#- saves other useful details used later into the jags object
######################
cctfit <- jags(data=jags.data, inits=jags.inits, parameters.to.save=jags.params,
n.chains=jags.chains, n.iter=(jags.iter+jags.burnin), n.burnin=jags.burnin, 
n.thin=jags.thin, model.file=textConnection(model.file))
cctfit$dataind <- datob$datind; cctfit$data <- datob$dat; cctfit$n <- n; cctfit$m <- m; cctfit$V <- V
cctfit$mval <- datob$mval; cctfit$itemdiff <- itemdiff; cctfit$checksrun <- 0; cctfit$whmodel <- datob$whmodel; cctfit$datob <- datob
if(cctfit$mval==1){cctfit$datamiss <- datob$datna}
if(cctfit$whmodel=="LTRM"){cctfit$C <- C}

######################
#Function used to calculate the Rhats
######################
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

######################
#Algorithm that corrects for label switching for the GCM
######################
labelswitchalggcm <- function(cctfit,chnind=0){

cctfit2 <- cctfit

nch <- cctfit2$BUGSoutput$n.chains
nsamp <- cctfit2$BUGSoutput$n.keep

if(nch != 1){
ntruths <- cctfit2$V
truths <- array(NA,c(cctfit2$m,nch,ntruths))
inds <- cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Z")]]
inds <- matrix(inds,cctfit2$m,cctfit2$V)
for(v in 1:cctfit$V){truths[,,v] <- t(apply(cctfit$BUGSoutput$sims.array[ ,,inds[,v]],c(2,3),mean))}

V <- cctfit2$V

chstart <- 1
if(length(chnind)==1){ 
chstart <- 2
chnind <- array(NA,c(V,nch))
chnind[1:V,1] <- 1:V

for(v in 1:V){
for(ch in chstart:nch){
Tind <- c(1:V)[-chnind[,ch][!is.na(chnind[,ch])]]
if(length(Tind)==0){Tind <- c(1:V)}

chnind[v,ch] <- which(max(cor(truths[1:cctfit2$m, 1, v],truths[1:cctfit2$m, ch, Tind]))==cor(truths[1:cctfit2$m, 1, v],truths[1:cctfit2$m, ch, ]))
}}
}

nsamp <- cctfit$BUGSoutput$n.keep

if(cctfit$itemdiff == 1){
inds2 <- cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "lam")]]
inds2 <- matrix(inds2,cctfit2$m,cctfit2$V)

inds <- rbind(inds,
inds2,
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "lammu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "lamtau")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "thmu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "thtau")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "gmu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "gtau")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "p")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "pi")]])
}else{
inds <- rbind(inds,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "thmu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "thtau")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "gmu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "gtau")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "p")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "pi")]])
}

einds <- cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]]
tmpeinds <- cctfit$BUGSoutput$sims.array[ ,,einds]

for(v in 1:cctfit$V){
for(ch in chstart:nch){
cctfit2$BUGSoutput$sims.array[ ,ch,inds[,v]] <- cctfit$BUGSoutput$sims.array[ ,ch,inds[,chnind[v,ch]]]  #this cctfit2 vs. cctfit difference is intentional
if(chstart==2){cctfit2$BUGSoutput$sims.array[ ,ch,einds][cctfit$BUGSoutput$sims.array[ ,ch,einds] == chnind[v,ch]] <- chnind[v,1]}
if(chstart==1){cctfit2$BUGSoutput$sims.array[ ,ch,einds][cctfit$BUGSoutput$sims.array[ ,ch,einds] == v] <- chnind[v,1]}
}}

}

cctfit2$BUGSoutput$sims.list[["e"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]]], c(nsamp*nch,cctfit$n))
cctfit2$BUGSoutput$mean$e <- apply(cctfit2$BUGSoutput$sims.list[["e"]],2,mean)

cctfit2$BUGSoutput$sims.matrix <- array(cctfit2$BUGSoutput$sims.array,c(nsamp*nch,dim(cctfit$BUGSoutput$sims.array)[3]))  #this cctfit2 vs. cctfit difference is intentional

cctfit2$BUGSoutput$sims.list[["th"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "th")]]], c(nsamp*nch,cctfit$n))
cctfit2$BUGSoutput$sims.list[["g"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "g")]]], c(nsamp*nch,cctfit$n))
if(cctfit$itemdiff == 1){cctfit2$BUGSoutput$sims.list[["lam"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "lam")]]], c(nsamp*nch,cctfit$m))}

cctfit2$BUGSoutput$sims.list[["Z"]] <- array(NA, c(nsamp*nch,cctfit$m,cctfit$V))
if(cctfit$itemdiff == 1){
cctfit2$BUGSoutput$sims.list[["lam"]] <- array(NA, c(nsamp*nch,cctfit$m,cctfit$V))
cctfit2$BUGSoutput$sims.list[["lammu"]] <- array(NA, c(nsamp*nch,cctfit$V))
cctfit2$BUGSoutput$sims.list[["lamtau"]] <- array(NA, c(nsamp*nch,cctfit$V))
}
cctfit2$BUGSoutput$sims.list[["thmu"]] <- array(NA, c(nsamp*nch,cctfit$V))
cctfit2$BUGSoutput$sims.list[["thtau"]] <- array(NA, c(nsamp*nch,cctfit$V))
cctfit2$BUGSoutput$sims.list[["gmu"]] <- array(NA, c(nsamp*nch,cctfit$V))
cctfit2$BUGSoutput$sims.list[["gtau"]] <- array(NA, c(nsamp*nch,cctfit$V))
cctfit2$BUGSoutput$sims.list[["p"]] <- array(NA, c(nsamp*nch,cctfit$V))
cctfit2$BUGSoutput$sims.list[["pi"]] <- array(NA, c(nsamp*nch,cctfit$V))

for(v in 1:V){
cctfit2$BUGSoutput$sims.list[["Z"]][,,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Z")]][1:cctfit$m +((v-1)*cctfit$m)]], c(nsamp*nch,cctfit$m))
if(cctfit$itemdiff == 1){
cctfit2$BUGSoutput$sims.list[["lam"]][,,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "lam")]][1:cctfit$m +((v-1)*cctfit$m)]], c(nsamp*nch,cctfit$m))
cctfit2$BUGSoutput$sims.list[["lammu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "lammu")]][v]], c(nsamp*nch))
cctfit2$BUGSoutput$sims.list[["lamtau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "lamtau")]][v]], c(nsamp*nch))
}
cctfit2$BUGSoutput$sims.list[["thmu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "thmu")]][v]], c(nsamp*nch))
cctfit2$BUGSoutput$sims.list[["thtau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "thtau")]][v]], c(nsamp*nch))
cctfit2$BUGSoutput$sims.list[["gmu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "gmu")]][v]], c(nsamp*nch))
cctfit2$BUGSoutput$sims.list[["gtau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "gtau")]][v]], c(nsamp*nch))
cctfit2$BUGSoutput$sims.list[["p"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "p")]][v]], c(nsamp*nch))
cctfit2$BUGSoutput$sims.list[["pi"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "pi")]][v]], c(nsamp*nch))

cctfit2$BUGSoutput$mean$Z[,v] <- apply(cctfit2$BUGSoutput$sims.list[["Z"]][,,v],2,mean)
if(cctfit$itemdiff == 1){
cctfit2$BUGSoutput$mean$lam[,v] <- apply(cctfit2$BUGSoutput$sims.list[["lam"]][,,v],2,mean)
cctfit2$BUGSoutput$mean$lammu[v] <- mean(cctfit2$BUGSoutput$sims.list[["lammu"]][,v])
cctfit2$BUGSoutput$mean$lamtau[v] <- mean(cctfit2$BUGSoutput$sims.list[["lamtau"]][,v])
}
cctfit2$BUGSoutput$mean$thmu[v] <- mean(cctfit2$BUGSoutput$sims.list[["thmu"]][,v])
cctfit2$BUGSoutput$mean$thtau[v] <- mean(cctfit2$BUGSoutput$sims.list[["thtau"]][,v])
cctfit2$BUGSoutput$mean$gmu[v] <- mean(cctfit2$BUGSoutput$sims.list[["gmu"]][,v])
cctfit2$BUGSoutput$mean$gtau[v] <- mean(cctfit2$BUGSoutput$sims.list[["gtau"]][,v])
cctfit2$BUGSoutput$mean$p[v] <- mean(cctfit2$BUGSoutput$sims.list[["p"]][,v])
cctfit2$BUGSoutput$mean$pi[v] <- mean(cctfit2$BUGSoutput$sims.list[["pi"]][,v])
}

if(nch != 1){
cctfit2$BUGSoutput$summary[,1] <- apply(apply(cctfit2$BUGSoutput$sims.array,c(2,3),mean),2,mean)
cctfit2$BUGSoutput$summary[,2] <- apply(apply(cctfit2$BUGSoutput$sims.array,c(2,3),sd),2,sd)
cctfit2$BUGSoutput$summary[,3:7] <- t(apply(cctfit2$BUGSoutput$sims.matrix,2,function(x) quantile(x,probs=c(.025,.25,.50,.75,.975))))
cctfit2$BUGSoutput$summary[,8] <- Rhat(cctfit2$BUGSoutput$sims.array)
cctfit2$BUGSoutput$summary[,8][is.nan(cctfit2$BUGSoutput$summary[,8])] <- 1.000000

dimnames(cctfit2$BUGSoutput$sims.array) <- list(NULL,NULL,rownames(cctfit$BUGSoutput$summary))
dimnames(cctfit2$BUGSoutput$sims.matrix) <- list(NULL,rownames(cctfit$BUGSoutput$summary))
}else{

cctfit2$BUGSoutput$summary[,1] <- apply(cctfit2$BUGSoutput$sims.array,2,mean)
cctfit2$BUGSoutput$summary[,2] <- apply(cctfit2$BUGSoutput$sims.array,2,sd)
cctfit2$BUGSoutput$summary[,3:7] <- t(apply(cctfit2$BUGSoutput$sims.matrix,2,function(x) quantile(x,probs=c(.025,.25,.50,.75,.975))))
cctfit2$BUGSoutput$summary <- cctfit2$BUGSoutput$summary[,-c(8,9)]
dimnames(cctfit2$BUGSoutput$sims.array) <- list(NULL,NULL,rownames(cctfit$BUGSoutput$summary))
dimnames(cctfit2$BUGSoutput$sims.matrix) <- list(NULL,rownames(cctfit$BUGSoutput$summary))
}

cctfit <- cctfit2; rm(cctfit2)

return(cctfit)
}

######################
#Algorithm that corrects for label switching for the LTRM
######################
labelswitchalgltrm <- function(cctfit,chnind=0){

cctfit2 <- cctfit

nch <- cctfit2$BUGSoutput$n.chains
nsamp <- cctfit2$BUGSoutput$n.keep

if(nch != 1){
ntruths <- cctfit2$V
truths <- array(NA,c(cctfit2$m,nch,ntruths))
inds <- cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "T")]]
inds <- matrix(inds,cctfit2$m,cctfit2$V)
for(v in 1:cctfit$V){truths[,,v] <- t(apply(cctfit$BUGSoutput$sims.array[ ,,inds[,v]],c(2,3),mean))}

V <- cctfit2$V

chstart <- 1
if(length(chnind)==1){ 
chstart <- 2
chnind <- array(NA,c(V,nch))
chnind[1:V,1] <- 1:V

for(v in 1:V){
for(ch in chstart:nch){
Tind <- c(1:V)[-chnind[,ch][!is.na(chnind[,ch])]]
if(length(Tind)==0){Tind <- c(1:V)}

chnind[v,ch] <- which(max(cor(truths[1:cctfit2$m, 1, v],truths[1:cctfit2$m, ch, Tind]))==cor(truths[1:cctfit2$m, 1, v],truths[1:cctfit2$m, ch, ]))
}}
}

nsamp <- cctfit$BUGSoutput$n.keep

if(cctfit$itemdiff == 1){
inds2 <- cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "lam")]]
inds2 <- matrix(inds2,cctfit2$m,cctfit2$V)

inds <- rbind(inds,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Tmu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Ttau")]],
inds2,
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "lammu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "lamtau")]],
matrix(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "gam")]],cctfit2$C-1,cctfit2$V),
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Emu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Etau")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "amu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "atau")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "bmu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "btau")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "pi")]])

rm(inds2)
}else{
inds <- rbind(inds,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Tmu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Ttau")]],
matrix(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "gam")]],cctfit2$C-1,cctfit2$V),
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Emu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Etau")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "amu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "atau")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "bmu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "btau")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "pi")]])
}

einds <- cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]]
tmpeinds <- cctfit$BUGSoutput$sims.array[ ,,einds]

for(v in 1:cctfit$V){
for(ch in chstart:nch){
cctfit2$BUGSoutput$sims.array[ ,ch,inds[,v]] <- cctfit$BUGSoutput$sims.array[ ,ch,inds[,chnind[v,ch]]]  #this cctfit2 vs. cctfit difference is intentional
if(chstart==2){cctfit2$BUGSoutput$sims.array[ ,ch,einds][cctfit$BUGSoutput$sims.array[ ,ch,einds] == chnind[v,ch]] <- chnind[v,1]}
if(chstart==1){cctfit2$BUGSoutput$sims.array[ ,ch,einds][cctfit$BUGSoutput$sims.array[ ,ch,einds] == v] <- chnind[v,1]}
}}

}
cctfit2$BUGSoutput$sims.list[["e"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]]], c(nsamp*nch,cctfit$n))
cctfit2$BUGSoutput$mean$e <- apply(cctfit2$BUGSoutput$sims.list[["e"]],2,mean)

cctfit2$BUGSoutput$sims.matrix <- array(cctfit2$BUGSoutput$sims.array,c(nsamp*nch,dim(cctfit$BUGSoutput$sims.array)[3]))  #this cctfit2 vs. cctfit difference is intentional

cctfit2$BUGSoutput$sims.list[["E"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "E")]]], c(nsamp*nch,cctfit$n))
cctfit2$BUGSoutput$sims.list[["a"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "a")]]], c(nsamp*nch,cctfit$n))
cctfit2$BUGSoutput$sims.list[["b"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "b")]]], c(nsamp*nch,cctfit$n))
if(cctfit$itemdiff == 1){cctfit2$BUGSoutput$sims.list[["lam"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "lam")]]], c(nsamp*nch,cctfit$m))}

cctfit2$BUGSoutput$sims.list[["T"]] <- array(NA, c(nsamp*nch,cctfit$m,cctfit$V))
cctfit2$BUGSoutput$sims.list[["Tmu"]] <- array(NA, c(nsamp*nch,cctfit$V))
cctfit2$BUGSoutput$sims.list[["Ttau"]] <- array(NA, c(nsamp*nch,cctfit$V))
if(cctfit$itemdiff == 1){
cctfit2$BUGSoutput$sims.list[["lam"]] <- array(NA, c(nsamp*nch,cctfit$m,cctfit$V))
cctfit2$BUGSoutput$sims.list[["lammu"]] <- array(NA, c(nsamp*nch,cctfit$V))
cctfit2$BUGSoutput$sims.list[["lamtau"]] <- array(NA, c(nsamp*nch,cctfit$V))
}
cctfit2$BUGSoutput$sims.list[["gam"]] <- array(NA, c(nsamp*nch,cctfit$C-1,cctfit$V))
cctfit2$BUGSoutput$sims.list[["Emu"]] <- array(NA, c(nsamp*nch,cctfit$V))
cctfit2$BUGSoutput$sims.list[["Etau"]] <- array(NA, c(nsamp*nch,cctfit$V))
cctfit2$BUGSoutput$sims.list[["amu"]] <- array(NA, c(nsamp*nch,cctfit$V))
cctfit2$BUGSoutput$sims.list[["atau"]] <- array(NA, c(nsamp*nch,cctfit$V))
cctfit2$BUGSoutput$sims.list[["bmu"]] <- array(NA, c(nsamp*nch,cctfit$V))
cctfit2$BUGSoutput$sims.list[["btau"]] <- array(NA, c(nsamp*nch,cctfit$V))
cctfit2$BUGSoutput$sims.list[["pi"]] <- array(NA, c(nsamp*nch,cctfit$V))

if(cctfit$C == 2 && length(dim(cctfit2$BUGSoutput$sims.list[["gam"]])) < 3 ){
cctfit2$BUGSoutput$sims.list[["gam"]] <- array(cctfit2$BUGSoutput$sims.list[["gam"]], c(nsamp*nch,cctfit$C-1,cctfit$V))
}

for(v in 1:cctfit2$V){
cctfit2$BUGSoutput$sims.list[["T"]][,,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "T")]][1:cctfit$m +((v-1)*cctfit$m)]], c(nsamp*nch,cctfit$m))
cctfit2$BUGSoutput$sims.list[["Tmu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Tmu")]][v]], c(nsamp*nch))
cctfit2$BUGSoutput$sims.list[["Ttau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Ttau")]][v]], c(nsamp*nch))
if(cctfit$itemdiff == 1){
cctfit2$BUGSoutput$sims.list[["lam"]][,,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "lam")]][1:cctfit$m +((v-1)*cctfit$m)]], c(nsamp*nch,cctfit$m))
cctfit2$BUGSoutput$sims.list[["lammu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "lammu")]][v]], c(nsamp*nch))
cctfit2$BUGSoutput$sims.list[["lamtau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "lamtau")]][v]], c(nsamp*nch))
}
cctfit2$BUGSoutput$sims.list[["gam"]][,,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "gam")]][1:(cctfit$C-1) +((v-1)*(cctfit$C-1))]], c(nsamp*nch,cctfit$C-1))
cctfit2$BUGSoutput$sims.list[["Emu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Emu")]][v]], c(nsamp*nch))
cctfit2$BUGSoutput$sims.list[["Etau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Etau")]][v]], c(nsamp*nch))
cctfit2$BUGSoutput$sims.list[["amu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "amu")]][v]], c(nsamp*nch))
cctfit2$BUGSoutput$sims.list[["atau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "atau")]][v]], c(nsamp*nch))
cctfit2$BUGSoutput$sims.list[["bmu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "bmu")]][v]], c(nsamp*nch))
cctfit2$BUGSoutput$sims.list[["btau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "btau")]][v]], c(nsamp*nch))
cctfit2$BUGSoutput$sims.list[["pi"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "pi")]][v]], c(nsamp*nch))

cctfit2$BUGSoutput$mean$T[,v] <- apply(cctfit2$BUGSoutput$sims.list[["T"]][,,v],2,mean)
cctfit2$BUGSoutput$mean$Tmu[v] <- mean(cctfit2$BUGSoutput$sims.list[["Tmu"]][,v])
cctfit2$BUGSoutput$mean$Ttau[v] <- mean(cctfit2$BUGSoutput$sims.list[["Ttau"]][,v])
if(cctfit$itemdiff == 1){
cctfit2$BUGSoutput$mean$lam[,v] <- apply(cctfit2$BUGSoutput$sims.list[["lam"]][,,v],2,mean)
cctfit2$BUGSoutput$mean$lammu[v] <- mean(cctfit2$BUGSoutput$sims.list[["lammu"]][,v])
cctfit2$BUGSoutput$mean$lamtau[v] <- mean(cctfit2$BUGSoutput$sims.list[["lamtau"]][,v])
}
if(cctfit$C == 2){
cctfit2$BUGSoutput$mean$gam[v] <- mean(cctfit2$BUGSoutput$sims.list[["gam"]][,,v])
}else{cctfit2$BUGSoutput$mean$gam[,v] <- apply(cctfit2$BUGSoutput$sims.list[["gam"]][,,v],2,mean)}
cctfit2$BUGSoutput$mean$Emu[v] <- mean(cctfit2$BUGSoutput$sims.list[["Emu"]][,v])
cctfit2$BUGSoutput$mean$Etau[v] <- mean(cctfit2$BUGSoutput$sims.list[["Etau"]][,v])
cctfit2$BUGSoutput$mean$amu[v] <- mean(cctfit2$BUGSoutput$sims.list[["amu"]][,v])
cctfit2$BUGSoutput$mean$atau[v] <- mean(cctfit2$BUGSoutput$sims.list[["atau"]][,v])
cctfit2$BUGSoutput$mean$bmu[v] <- mean(cctfit2$BUGSoutput$sims.list[["bmu"]][,v])
cctfit2$BUGSoutput$mean$btau[v] <- mean(cctfit2$BUGSoutput$sims.list[["btau"]][,v])
cctfit2$BUGSoutput$mean$pi[v] <- mean(cctfit2$BUGSoutput$sims.list[["pi"]][,v])
}

if(nch != 1){
cctfit2$BUGSoutput$summary[,1] <- apply(apply(cctfit2$BUGSoutput$sims.array,c(2,3),mean),2,mean)
cctfit2$BUGSoutput$summary[,2] <- apply(apply(cctfit2$BUGSoutput$sims.array,c(2,3),sd),2,sd)
cctfit2$BUGSoutput$summary[,3:7] <- t(apply(cctfit2$BUGSoutput$sims.matrix,2,function(x) quantile(x,probs=c(.025,.25,.50,.75,.975))))
cctfit2$BUGSoutput$summary[,8] <- Rhat(cctfit2$BUGSoutput$sims.array)
cctfit2$BUGSoutput$summary[,8][is.nan(cctfit2$BUGSoutput$summary[,8])] <- 1.000000
dimnames(cctfit2$BUGSoutput$sims.array) <- list(NULL,NULL,rownames(cctfit$BUGSoutput$summary))
dimnames(cctfit2$BUGSoutput$sims.matrix) <- list(NULL,rownames(cctfit$BUGSoutput$summary))
}else{

cctfit2$BUGSoutput$summary[,1] <- apply(cctfit2$BUGSoutput$sims.array,2,mean)
cctfit2$BUGSoutput$summary[,2] <- apply(cctfit2$BUGSoutput$sims.array,2,sd)
cctfit2$BUGSoutput$summary[,3:7] <- t(apply(cctfit2$BUGSoutput$sims.matrix,2,function(x) quantile(x,probs=c(.025,.25,.50,.75,.975))))
cctfit2$BUGSoutput$summary <- cctfit2$BUGSoutput$summary[,-c(8,9)]
dimnames(cctfit2$BUGSoutput$sims.array) <- list(NULL,NULL,rownames(cctfit$BUGSoutput$summary))
dimnames(cctfit2$BUGSoutput$sims.matrix) <- list(NULL,rownames(cctfit$BUGSoutput$summary))
}

cctfit <- cctfit2; rm(cctfit2)

return(cctfit)
}

######################
#Algorithm that corrects for label switching for the CRM
######################
labelswitchalgcrm <- function(cctfit,chnind=0){

cctfit2 <- cctfit

nch <- cctfit2$BUGSoutput$n.chains
nsamp <- cctfit2$BUGSoutput$n.keep

if(nch != 1){
ntruths <- cctfit2$V
truths <- array(NA,c(cctfit2$m,nch,ntruths))
inds <- cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "T")]]
inds <- matrix(inds,cctfit2$m,cctfit2$V)
for(v in 1:cctfit$V){truths[,,v] <- t(apply(cctfit$BUGSoutput$sims.array[ ,,inds[,v]],c(2,3),mean))}

V <- cctfit2$V

chstart <- 1
if(length(chnind)==1){ 
chstart <- 2
chnind <- array(NA,c(V,nch))
chnind[1:cctfit2$V,1] <- 1:cctfit2$V

for(v in 1:cctfit2$V){
for(ch in chstart:nch){
Tind <- c(1:cctfit2$V)[-chnind[,ch][!is.na(chnind[,ch])]]
if(length(Tind)==0){Tind <- c(1:cctfit2$V)}

chnind[v,ch] <- which(max(cor(truths[1:cctfit2$m, 1, v],truths[1:cctfit2$m, ch, Tind]))==cor(truths[1:cctfit2$m, 1, v],truths[1:cctfit2$m, ch, ]))
}}
}

nsamp <- cctfit$BUGSoutput$n.keep

if(cctfit$itemdiff == 1){
inds2 <- cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "lam")]]
inds2 <- matrix(inds2,cctfit2$m,cctfit2$V)

inds <- rbind(inds,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Tmu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Ttau")]],
inds2,
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "lammu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "lamtau")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Emu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Etau")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "amu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "atau")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "bmu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "btau")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "pi")]])

rm(inds2)
}else{
inds <- rbind(inds,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Tmu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Ttau")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Emu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Etau")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "amu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "atau")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "bmu")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "btau")]],
cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "pi")]])
}
einds <- cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]]
tmpeinds <- cctfit$BUGSoutput$sims.array[ ,,einds]

for(v in 1:cctfit$V){
for(ch in chstart:nch){
cctfit2$BUGSoutput$sims.array[ ,ch,inds[,v]] <- cctfit$BUGSoutput$sims.array[ ,ch,inds[,chnind[v,ch]]]  #this cctfit2 vs. cctfit difference is intentional

if(chstart==2){cctfit2$BUGSoutput$sims.array[ ,ch,einds][cctfit$BUGSoutput$sims.array[ ,ch,einds] == chnind[v,ch]] <- chnind[v,1]}
if(chstart==1){cctfit2$BUGSoutput$sims.array[ ,ch,einds][cctfit$BUGSoutput$sims.array[ ,ch,einds] == v] <- chnind[v,1]}
}}

}

cctfit2$BUGSoutput$sims.list[["e"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]]], c(nsamp*nch,cctfit$n))
cctfit2$BUGSoutput$mean$e <- apply(cctfit2$BUGSoutput$sims.list[["e"]],2,mean)

cctfit2$BUGSoutput$sims.matrix <- array(cctfit2$BUGSoutput$sims.array,c(nsamp*nch,dim(cctfit$BUGSoutput$sims.array)[3]))  #this cctfit2 vs. cctfit difference is intentional

cctfit2$BUGSoutput$sims.list[["E"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "E")]]], c(nsamp*nch,cctfit$n))
cctfit2$BUGSoutput$sims.list[["a"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "a")]]], c(nsamp*nch,cctfit$n))
cctfit2$BUGSoutput$sims.list[["b"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "b")]]], c(nsamp*nch,cctfit$n))
if(cctfit$itemdiff == 1){cctfit2$BUGSoutput$sims.list[["lam"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "lam")]]], c(nsamp*nch,cctfit$m))}

cctfit2$BUGSoutput$sims.list[["T"]] <- array(NA, c(nsamp*nch,cctfit$m,cctfit$V))
cctfit2$BUGSoutput$sims.list[["Tmu"]] <- array(NA, c(nsamp*nch,cctfit$V))
cctfit2$BUGSoutput$sims.list[["Ttau"]] <- array(NA, c(nsamp*nch,cctfit$V))
if(cctfit$itemdiff == 1){
cctfit2$BUGSoutput$sims.list[["lam"]] <- array(NA, c(nsamp*nch,cctfit$m,cctfit$V))
cctfit2$BUGSoutput$sims.list[["lammu"]] <- array(NA, c(nsamp*nch,cctfit$V))
cctfit2$BUGSoutput$sims.list[["lamtau"]] <- array(NA, c(nsamp*nch,cctfit$V))
}
cctfit2$BUGSoutput$sims.list[["Emu"]] <- array(NA, c(nsamp*nch,cctfit$V))
cctfit2$BUGSoutput$sims.list[["Etau"]] <- array(NA, c(nsamp*nch,cctfit$V))
cctfit2$BUGSoutput$sims.list[["amu"]] <- array(NA, c(nsamp*nch,cctfit$V))
cctfit2$BUGSoutput$sims.list[["atau"]] <- array(NA, c(nsamp*nch,cctfit$V))
cctfit2$BUGSoutput$sims.list[["bmu"]] <- array(NA, c(nsamp*nch,cctfit$V))
cctfit2$BUGSoutput$sims.list[["btau"]] <- array(NA, c(nsamp*nch,cctfit$V))
cctfit2$BUGSoutput$sims.list[["pi"]] <- array(NA, c(nsamp*nch,cctfit$V))

for(v in 1:cctfit2$V){
cctfit2$BUGSoutput$sims.list[["T"]][,,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "T")]][1:cctfit$m +((v-1)*cctfit$m)]], c(nsamp*nch,cctfit$m))
cctfit2$BUGSoutput$sims.list[["Tmu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Tmu")]][v]], c(nsamp*nch))
cctfit2$BUGSoutput$sims.list[["Ttau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Ttau")]][v]], c(nsamp*nch))
if(cctfit$itemdiff == 1){
cctfit2$BUGSoutput$sims.list[["lam"]][,,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "lam")]][1:cctfit$m +((v-1)*cctfit$m)]], c(nsamp*nch,cctfit$m))
cctfit2$BUGSoutput$sims.list[["lammu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "lammu")]][v]], c(nsamp*nch))
cctfit2$BUGSoutput$sims.list[["lamtau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "lamtau")]][v]], c(nsamp*nch))
}
cctfit2$BUGSoutput$sims.list[["Emu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Emu")]][v]], c(nsamp*nch))
cctfit2$BUGSoutput$sims.list[["Etau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Etau")]][v]], c(nsamp*nch))
cctfit2$BUGSoutput$sims.list[["amu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "amu")]][v]], c(nsamp*nch))
cctfit2$BUGSoutput$sims.list[["atau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "atau")]][v]], c(nsamp*nch))
cctfit2$BUGSoutput$sims.list[["bmu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "bmu")]][v]], c(nsamp*nch))
cctfit2$BUGSoutput$sims.list[["btau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "btau")]][v]], c(nsamp*nch))
cctfit2$BUGSoutput$sims.list[["pi"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "pi")]][v]], c(nsamp*nch))

cctfit2$BUGSoutput$mean$T[,v] <- apply(cctfit2$BUGSoutput$sims.list[["T"]][,,v],2,mean)
cctfit2$BUGSoutput$mean$Tmu[v] <- mean(cctfit2$BUGSoutput$sims.list[["Tmu"]][,v])
cctfit2$BUGSoutput$mean$Ttau[v] <- mean(cctfit2$BUGSoutput$sims.list[["Ttau"]][,v])
if(cctfit$itemdiff == 1){
cctfit2$BUGSoutput$mean$lam[,v] <- apply(cctfit2$BUGSoutput$sims.list[["lam"]][,,v],2,mean)
cctfit2$BUGSoutput$mean$lammu[v] <- mean(cctfit2$BUGSoutput$sims.list[["lammu"]][,v])
cctfit2$BUGSoutput$mean$lamtau[v] <- mean(cctfit2$BUGSoutput$sims.list[["lamtau"]][,v])
}
cctfit2$BUGSoutput$mean$Emu[v] <- mean(cctfit2$BUGSoutput$sims.list[["Emu"]][,v])
cctfit2$BUGSoutput$mean$Etau[v] <- mean(cctfit2$BUGSoutput$sims.list[["Etau"]][,v])
cctfit2$BUGSoutput$mean$amu[v] <- mean(cctfit2$BUGSoutput$sims.list[["amu"]][,v])
cctfit2$BUGSoutput$mean$atau[v] <- mean(cctfit2$BUGSoutput$sims.list[["atau"]][,v])
cctfit2$BUGSoutput$mean$bmu[v] <- mean(cctfit2$BUGSoutput$sims.list[["bmu"]][,v])
cctfit2$BUGSoutput$mean$btau[v] <- mean(cctfit2$BUGSoutput$sims.list[["btau"]][,v])
cctfit2$BUGSoutput$mean$pi[v] <- mean(cctfit2$BUGSoutput$sims.list[["pi"]][,v])
}

if(nch != 1){
cctfit2$BUGSoutput$summary[,1] <- apply(apply(cctfit2$BUGSoutput$sims.array,c(2,3),mean),2,mean)
cctfit2$BUGSoutput$summary[,2] <- apply(apply(cctfit2$BUGSoutput$sims.array,c(2,3),sd),2,sd)
cctfit2$BUGSoutput$summary[,3:7] <- t(apply(cctfit2$BUGSoutput$sims.matrix,2,function(x) quantile(x,probs=c(.025,.25,.50,.75,.975))))
cctfit2$BUGSoutput$summary[,8] <- Rhat(cctfit2$BUGSoutput$sims.array)
cctfit2$BUGSoutput$summary[,8][is.nan(cctfit2$BUGSoutput$summary[,8])] <- 1.000000
dimnames(cctfit2$BUGSoutput$sims.array) <- list(NULL,NULL,rownames(cctfit$BUGSoutput$summary))
dimnames(cctfit2$BUGSoutput$sims.matrix) <- list(NULL,rownames(cctfit$BUGSoutput$summary))
}else{

cctfit2$BUGSoutput$summary[,1] <- apply(cctfit2$BUGSoutput$sims.array,2,mean)
cctfit2$BUGSoutput$summary[,2] <- apply(cctfit2$BUGSoutput$sims.array,2,sd)
cctfit2$BUGSoutput$summary[,3:7] <- t(apply(cctfit2$BUGSoutput$sims.matrix,2,function(x) quantile(x,probs=c(.025,.25,.50,.75,.975))))
cctfit2$BUGSoutput$summary <- cctfit2$BUGSoutput$summary[,-c(8,9)]
dimnames(cctfit2$BUGSoutput$sims.array) <- list(NULL,NULL,rownames(cctfit$BUGSoutput$summary))
dimnames(cctfit2$BUGSoutput$sims.matrix) <- list(NULL,rownames(cctfit$BUGSoutput$summary))
}

cctfit <- cctfit2; rm(cctfit2)

return(cctfit)
}


######################
#This is run after the jags inference, we are still in the 'Apply CCT Model' function
#- Reports the number of Rhats above 1.05 and 1.10   (before applying the algorithm that corrects for label-switching)
#- Detects if a mixture model was applied
#- If so, the appropriate label-correcting algorithm is applied
#- Recalculates the Rhats, DIC, and statistics for the parameters
#- Reports the number of Rhats above 1.05 and 1.10 
#- Outputs the DIC that is calculated after the label-correcting algorithm (if applicable)
######################
message("\n ...Inference complete, data is saved as 'cctfit'")
if(clusters > 1){
message("\n    'cctfit$respmem' provides the respondent clustering")
}

message("\n ...Performing final calculations")
if(cctfit$BUGSoutput$n.chains > 1){
cctfit$BUGSoutput$summary[,8] <- Rhat(cctfit$BUGSoutput$sims.array)
cctfit$BUGSoutput$summary[,8][is.nan(cctfit$BUGSoutput$summary[,8])] <- 1.000000
if(cctfit$whmodel== "GCM"){
message(paste("\nFor Continuous Parameters"))
message(paste("Number of Rhats above 1.10 : ",sum(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]],cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Z")]]),8]>1.10),"/",length(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]],cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Z")]]),8]),"\nNumber of Rhats above 1.05 : ",sum(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]],cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Z")]]),8]>1.050),"/",length(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]],cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Z")]]),8]),sep=""))
}else{
message(paste("Number of Rhats above 1.10 : ",sum(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]]),8]>1.10),"/",length(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]]),8]),"\nNumber of Rhats above 1.05 : ",sum(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]]),8]>1.050),"/",length(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]]),8]),sep=""))
}
}

if(cctfit$V>1 && cctfit$BUGSoutput$n.chains == 1){

Mode <<- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }
tmp <- unique(apply(cctfit$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]]],c(2),Mode))
if(length(tmp) == 1){
message(paste("\n ...This chain has ",length(tmp)," culture rather than the ",cctfit$V," cultures requested", sep=""))
message(paste("\n ...Try running the inference again",sep="" ))
}
if(length(tmp) != 1 && length(tmp) < cctfit$V){
message(paste("\n ...This chain has ",tmp," cultures rather than the ",cctfit$V," cultures requested", sep=""))
message(paste("\n ...Try running the inference again",sep="" ))
}
}

if(cctfit$V>1 && cctfit$BUGSoutput$n.chains > 1){
message("\n ...More than 1 culture applied with more than 1 chain")

einds <- cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]]

Mode <<- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

tmp <- apply(apply(cctfit$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]]],c(2,3),Mode),1,unique)

chntoremove <- NULL
if(is.list(tmp)){
for(i in 1:cctfit$BUGSoutput$n.chains){
if(length(tmp[[i]]) != cctfit$V){chntoremove <- c(chntoremove,i)}
}
}

if(length(dim(tmp))==2){
for(i in 1:cctfit$BUGSoutput$n.chains){
if(length(tmp[,i]) != cctfit$V){chntoremove <- c(chntoremove,i)}
}
}

if(is.null(dim(tmp)) && !is.list(tmp)){chntoremove <- c(1:cctfit$BUGSoutput$n.chains)}


if(length(chntoremove)>0){
if(length(chntoremove) < cctfit$BUGSoutput$n.chains){

if(length(chntoremove)==1){
message(paste("\n ...",length(chntoremove)," chain out of ",cctfit$BUGSoutput$n.chains," had fewer than ",cctfit$V," cultures requested",sep=""))
message(paste("\n ...", "removing the ", length(chntoremove)," chain", sep=""))
}else{
message(paste("\n ...",length(chntoremove)," chains out of ",cctfit$BUGSoutput$n.chains," had fewer than ",cctfit$V," cultures requested",sep=""))
message(paste("\n ...", "removing these ", length(chntoremove)," chains", sep=""))
}

cctfit$BUGSoutput$n.chains <- cctfit$BUGSoutput$n.chains-length(chntoremove)
cctfit$BUGSoutput$n.chain <- cctfit$BUGSoutput$n.chains
cctfit$BUGSoutput$n.sims <- cctfit$BUGSoutput$n.chains*cctfit$BUGSoutput$n.keep
if(cctfit$BUGSoutput$n.chain == 1){
cctfit$BUGSoutput$sims.array <- array(cctfit$BUGSoutput$sims.array[,-chntoremove,], c(dim(cctfit$BUGSoutput$sims.array)[1],1,dim(cctfit$BUGSoutput$sims.array)[3]))
}else{cctfit$BUGSoutput$sims.array <- cctfit$BUGSoutput$sims.array[,-chntoremove,]}

}else{
message(paste("\n ...All chains out of ",cctfit$BUGSoutput$n.chains," had fewer than ",cctfit$V," cultures requested", sep=""))
message(paste("\n ...Try running the inference again",sep="" ))
}
}

message("\n ...Computing the most-consistent labeling across chains")

if(cctfit$whmodel=="GCM"){

cctfit <- labelswitchalggcm(cctfit)

Mode <<- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }
cctfit$respmem <- apply(cctfit$BUGSoutput$sims.list$e[,],2,Mode)
tmeans <- cctfit$BUGSoutput$mean$Z
tmeans[tmeans<.5] <- tmeans[tmeans<.5]+1
ind <- rank(apply(abs(1-tmeans),2,mean))
} 


if(cctfit$whmodel=="LTRM"){

cctfit <- labelswitchalgltrm(cctfit)

Mode <<- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
cctfit$respmem <- apply(cctfit$BUGSoutput$sims.list$e[,],2,Mode)
tmeans <- cctfit$BUGSoutput$mean$T
ind <- rank(-apply(tmeans,2,sd))

}

if(cctfit$whmodel=="CRM"){

cctfit <- labelswitchalgcrm(cctfit)

Mode <<- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
cctfit$respmem <- apply(cctfit$BUGSoutput$sims.list$e[,],2,Mode)
tmeans <- cctfit$BUGSoutput$mean$T
ind <- rank(-apply(tmeans,2,sd))
}

if(cctfit$BUGSoutput$n.chains > 1){
message(paste("\nFor Continuous Parameters:"))
if(cctfit$whmodel== "GCM"){
message(paste("Number of Rhats above 1.10 : ",sum(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]],cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Z")]]),8]>1.10),"/",length(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]],cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Z")]]),8]),"\nNumber of Rhats above 1.05 : ",sum(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]],cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Z")]]),8]>1.050),"/",length(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]],cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Z")]]),8]),sep=""))
}else{
message(paste("Number of Rhats above 1.10 : ",sum(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]]),8]>1.10),"/",length(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]]),8]),"\nNumber of Rhats above 1.05 : ",sum(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]]),8]>1.050),"/",length(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]]),8]),sep=""))
}
}
if(cctfit$BUGSoutput$n.chains > 1){
if(gui==1){message(paste("\nFor Discrete Parameters:"))
message(paste("Type 'dtraceplot()' to see their trace plots")) }
alltraceplot <<- function(){traceplot(cctfit,mfrow=c(4,4))}
dtraceplot <<- function(){
if(cctfit$whmodel == "GCM"){
if(cctfit$V == 1){
cctfit2 <- cctfit
cctfit2$BUGSoutput$sims.array <- cctfit2$BUGSoutput$sims.array[,,c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Z")]])]
traceplot(cctfit2,mfrow=c(4,4))
}else{
cctfit2 <- cctfit
cctfit2$BUGSoutput$sims.array <- cctfit2$BUGSoutput$sims.array[,,c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Z")]],cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]])]
traceplot(cctfit2,mfrow=c(4,4))
}
}else{
if(cctfit$V == 1){
if(gui==1){message("\n There are no discrete nodes in this inference")}
return()
}else{
cctfit2 <- cctfit
cctfit2$BUGSoutput$sims.array <- cctfit2$BUGSoutput$sims.array[,,c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]])]
traceplot(cctfit2,mfrow=c(4,4))
}
}
}
}

}else{

if(cctfit$whmodel == "GCM" && cctfit$BUGSoutput$n.chains > 1){
message(paste("\nFor Discrete Parameters:"))
message(paste("Type 'dtraceplot()' to see their trace plots"))
alltraceplot <<- function(){traceplot(cctfit,mfrow=c(4,4))}
dtraceplot <<- function(){
if(cctfit$whmodel == "GCM"){
if(cctfit$V == 1){
cctfit2 <- cctfit
cctfit2$BUGSoutput$sims.array <- cctfit2$BUGSoutput$sims.array[,,c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Z")]])]
traceplot(cctfit2,mfrow=c(4,4))
}else{
cctfit2 <- cctfit
cctfit2$BUGSoutput$sims.array <- cctfit2$BUGSoutput$sims.array[,,c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "Z")]],cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]])]
traceplot(cctfit2,mfrow=c(4,4))
}
}else{
if(cctfit$V == 1){
message("\n There are no discrete nodes in this inference")
return()
}else{
cctfit2 <- cctfit
cctfit2$BUGSoutput$sims.array <- cctfit2$BUGSoutput$sims.array[,,c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "e")]])]
traceplot(cctfit2,mfrow=c(4,4))
}
}
}
}

}

######################
#DIC is calculated for the model that was applied
#The likelihood of the model is evaluated at each node of each sample and saved to cctfit$Lik
######################
message("\n ...Calculating DIC")
cctfit$BUGSoutput$pD <- var(cctfit$BUGSoutput$sims.list$deviance[,1])/2
nsimstouse <- min(cctfit$BUGSoutput$n.sims,1000)
ind <- sample(1:cctfit$BUGSoutput$n.sims,nsimstouse)

storelik <- 0

if(cctfit$whmodel=="GCM"){

if(cctfit$V == 1){
cctfit$V <- 1;
cctfit$BUGSoutput$sims.list$e <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$n));
if(length(dim(cctfit$BUGSoutput$sims.list$Z)) < 3){
cctfit$BUGSoutput$sims.list$Z <- array(cctfit$BUGSoutput$sims.list[["Z"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
if(cctfit$itemdiff == 1){
if(length(dim(cctfit$BUGSoutput$sims.list$lam)) < 3){
cctfit$BUGSoutput$sims.list$lam <- array(cctfit$BUGSoutput$sims.list[["lam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
}
}
if(cctfit$itemdiff == 0){
if(cctfit$V == 1){
cctfit$BUGSoutput$sims.list$lam <- array(.5,c(cctfit$BUGSoutput$n.sims,cctfit$m,1))
}
if(cctfit$V > 1){
cctfit$BUGSoutput$sims.list$lam <- array(.5,c(cctfit$BUGSoutput$n.sims,cctfit$m,cctfit$V))
}
}

cctfit$BUGSoutput$sims.list$D <- array(-1,c(cctfit$BUGSoutput$n.sims,cctfit$n,cctfit$m))

storefulllik <- 0
if(storefulllik == 1){
cctfit$Lik  <- array(NA, c(cctfit$n,cctfit$m,cctfit$BUGSoutput$n.sims));
for(samp in 1:cctfit$BUGSoutput$n.sims){
cctfit$BUGSoutput$sims.list[["D"]][samp,,] <- (cctfit$BUGSoutput$sims.list[["th"]][samp,]*t(1-cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])) / 
( (cctfit$BUGSoutput$sims.list[["th"]][samp,]*t(1-cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])) + 
((1-cctfit$BUGSoutput$sims.list[["th"]][samp,])*t(cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])) ) 

cctfit$Lik[,,samp] <- t( ( (t(cctfit$BUGSoutput$sims.list[["D"]][samp,,])+(t(1-cctfit$BUGSoutput$sims.list[["D"]][samp,,])%*%diag(cctfit$BUGSoutput$sims.list[["g"]][samp,])))^(t(cctfit$data[,])*cctfit$BUGSoutput$sims.list[["Z"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]]))*(t(1-cctfit$BUGSoutput$sims.list[["D"]][samp,,])%*%diag(cctfit$BUGSoutput$sims.list[["g"]][samp,]))^(t(cctfit$data[,])*(1-cctfit$BUGSoutput$sims.list[["Z"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]]))*(t(1-cctfit$BUGSoutput$sims.list[["D"]][samp,,])%*%diag(1-cctfit$BUGSoutput$sims.list[["g"]][samp,]))^((1-t(cctfit$data[,]))*cctfit$BUGSoutput$sims.list[["Z"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])*(t(cctfit$BUGSoutput$sims.list[["D"]][samp,,])+t(1-cctfit$BUGSoutput$sims.list[["D"]][samp,,])%*%diag(1-cctfit$BUGSoutput$sims.list[["g"]][samp,]))^((1-t(cctfit$data[,]))*(1-cctfit$BUGSoutput$sims.list[["Z"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]]))  )
if(cctfit$mval==1){cctfit$Lik[cbind(datob$thena,samp)] <- 1}
}
 
cctfit$BUGSoutput$sims.list$deviance <- array(-2*apply(log(cctfit$Lik),3,sum),c(cctfit$BUGSoutput$n.sims,1))
}else{
cctfit$LogLik  <- array(NA, c(cctfit$BUGSoutput$n.sims));
for(samp in 1:cctfit$BUGSoutput$n.sims){
Dtmp <- (cctfit$BUGSoutput$sims.list[["th"]][samp,]*t(1-cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])) / 
( (cctfit$BUGSoutput$sims.list[["th"]][samp,]*t(1-cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])) + 
((1-cctfit$BUGSoutput$sims.list[["th"]][samp,])*t(cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])) ) 

Liktmp <- t( ( (t(Dtmp)+(t(1-Dtmp)%*%diag(cctfit$BUGSoutput$sims.list[["g"]][samp,])))^(t(cctfit$data[,])*cctfit$BUGSoutput$sims.list[["Z"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]]))*(t(1-Dtmp)%*%diag(cctfit$BUGSoutput$sims.list[["g"]][samp,]))^(t(cctfit$data[,])*(1-cctfit$BUGSoutput$sims.list[["Z"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]]))*(t(1-Dtmp)%*%diag(1-cctfit$BUGSoutput$sims.list[["g"]][samp,]))^((1-t(cctfit$data[,]))*cctfit$BUGSoutput$sims.list[["Z"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])*(t(Dtmp)+t(1-Dtmp)%*%diag(1-cctfit$BUGSoutput$sims.list[["g"]][samp,]))^((1-t(cctfit$data[,]))*(1-cctfit$BUGSoutput$sims.list[["Z"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]]))  )
if(cctfit$mval==1){Liktmp[datob$thena] <- 1}

cctfit$LogLik[samp] <- sum(log(Liktmp))
}

cctfit$BUGSoutput$sims.list$deviance <- array(-2*(cctfit$LogLik),c(cctfit$BUGSoutput$n.sims,1))
}

cctfit$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "deviance")]]] <- array(cctfit$BUGSoutput$sims.list$deviance,c(cctfit$BUGSoutput$n.keep,cctfit$BUGSoutput$n.chains))
cctfit$BUGSoutput$sims.matrix[,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "deviance")]]] <- cctfit$BUGSoutput$sims.list$deviance[,1]
 
if(sum(cctfit$BUGSoutput$sims.list$deviance == Inf) >0){ 
cctfit$BUGSoutput$mean$deviance <- mean(cctfit$BUGSoutput$sims.list$deviance[cctfit$BUGSoutput$sims.list$deviance != Inf,1],na.rm=TRUE) #Dbar, also known as deviance
cctfit$BUGSoutput$pD <- var(cctfit$BUGSoutput$sims.list$deviance[cctfit$BUGSoutput$sims.list$deviance != Inf,1],na.rm=TRUE)/2  #pD, variance of the deviance, divided by 2
cctfit$BUGSoutput$DIC <- cctfit$BUGSoutput$mean$deviance + cctfit$BUGSoutput$pD
}else{
cctfit$BUGSoutput$mean$deviance <- mean(cctfit$BUGSoutput$sims.list$deviance) #Dbar, also known as deviance
cctfit$BUGSoutput$pD <- var(cctfit$BUGSoutput$sims.list$deviance[,1])/2  #pD, variance of the deviance, divided by 2
cctfit$BUGSoutput$DIC <- cctfit$BUGSoutput$mean$deviance + cctfit$BUGSoutput$pD
}

}

if(cctfit$whmodel=="LTRM"){
if(cctfit$V == 1){
cctfit$V <- 1;
cctfit$BUGSoutput$sims.list$e <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$n));
if(length(dim(cctfit$BUGSoutput$sims.list$T)) < 3){
cctfit$BUGSoutput$sims.list$T <- array(cctfit$BUGSoutput$sims.list[["T"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
if(length(dim(cctfit$BUGSoutput$sims.list$gam))<3){cctfit$BUGSoutput$sims.list$gam <- array(cctfit$BUGSoutput$sims.list[["gam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$C-1,1))}
if(cctfit$itemdiff == 1){
if(length(dim(cctfit$BUGSoutput$sims.list$lam)) < 3){
cctfit$BUGSoutput$sims.list$lam <- array(cctfit$BUGSoutput$sims.list[["lam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
}
}
if(cctfit$itemdiff == 0){
if(cctfit$V == 1){
cctfit$BUGSoutput$sims.list$lam <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$m,1))
}
if(cctfit$V > 1){
cctfit$BUGSoutput$sims.list$lam <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$m,cctfit$V))
}
}


storefulllik <- 0
if(storefulllik == 1){

cctfit$BUGSoutput$sims.list$tau <- array(-1,c(cctfit$BUGSoutput$n.sims,cctfit$n,cctfit$m))
cctfit$ppdelta <- array(NA, c(cctfit$n,cctfit$C-1,cctfit$BUGSoutput$n.sims))
cctfit$ppdeltafull <- array(NA, c(cctfit$n,cctfit$C+1,cctfit$BUGSoutput$n.sims))
cctfit$ppdeltafull[,1,] <- -1000000; cctfit$ppdeltafull[,cctfit$C+1,] <- 1000000
cctfit$Lik  <- array(NA, c(cctfit$n,cctfit$m,cctfit$BUGSoutput$n.sims));

for(samp in 1:cctfit$BUGSoutput$n.sims){
cctfit$BUGSoutput$sims.list[["tau"]][samp,,] <- cctfit$BUGSoutput$sims.list[["E"]][samp,]*t(1/cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])
cctfit$ppdeltafull[,2:(cctfit$C),samp] <- cctfit$BUGSoutput$sims.list[["a"]][samp,]*t(cctfit$BUGSoutput$sims.list[["gam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])+cctfit$BUGSoutput$sims.list[["b"]][samp,]

for(i in 1:cctfit$n){
cctfit$Lik[i,,samp] <-  pnorm(cctfit$ppdeltafull[i,cctfit$data[i,]+1,samp] ,cctfit$BUGSoutput$sims.list[["T"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,i]],cctfit$BUGSoutput$sims.list[["tau"]][samp,i,]^-.5)-pnorm(cctfit$ppdeltafull[i,cctfit$data[i,],samp] ,cctfit$BUGSoutput$sims.list[["T"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,i]],cctfit$BUGSoutput$sims.list[["tau"]][samp,i,]^-.5)
if(cctfit$mval==1){cctfit$Lik[cbind(datob$thena,samp)] <- 1}
}
}

if(cctfit$C==2){cctfit$Lik[cctfit$Lik == 0] <- 0.001}

cctfit$BUGSoutput$sims.list$deviance <- array(-2*apply(log(cctfit$Lik),3,sum),c(cctfit$BUGSoutput$n.sims,1))
}else{
cctfit$LogLik  <- array(NA, c(cctfit$BUGSoutput$n.sims));
Liktmp <- array(NA, c(cctfit$n,cctfit$m))
deltatmp <- array(NA, c(cctfit$n,cctfit$C+1))
deltatmp[,1] <- -100000; deltatmp[,cctfit$C+1] <- 1000000

for(samp in 1:cctfit$BUGSoutput$n.sims){
tautmp <- cctfit$BUGSoutput$sims.list[["E"]][samp,]*t(1/cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])
deltatmp[,2:(cctfit$C)] <- cctfit$BUGSoutput$sims.list[["a"]][samp,]*t(cctfit$BUGSoutput$sims.list[["gam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])+cctfit$BUGSoutput$sims.list[["b"]][samp,]

for(i in 1:cctfit$n){
Liktmp[i,] <-  pnorm(deltatmp[i,cctfit$data[i,]+1] ,cctfit$BUGSoutput$sims.list[["T"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,i]],tautmp[i,]^-.5)-pnorm(deltatmp[i,cctfit$data[i,]] ,cctfit$BUGSoutput$sims.list[["T"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,i]],tautmp[i,]^-.5)
if(cctfit$mval==1){cctfit$Lik[datob$thena] <- 1}
}
if(cctfit$C==2){Liktmp[Liktmp == 0] <- 0.001}

cctfit$LogLik[samp] <- sum(log(Liktmp))
}

cctfit$BUGSoutput$sims.list$deviance <- array(-2*(cctfit$LogLik),c(cctfit$BUGSoutput$n.sims,1))
}

cctfit$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "deviance")]]] <- array(cctfit$BUGSoutput$sims.list$deviance,c(cctfit$BUGSoutput$n.keep,cctfit$BUGSoutput$n.chains))
cctfit$BUGSoutput$sims.matrix[,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "deviance")]]] <- cctfit$BUGSoutput$sims.list$deviance[,1]
 
 
if(sum(cctfit$BUGSoutput$sims.list$deviance == Inf) >0){ 
cctfit$BUGSoutput$mean$deviance <- mean(cctfit$BUGSoutput$sims.list$deviance[cctfit$BUGSoutput$sims.list$deviance != Inf,1],na.rm=TRUE) #Dbar, also known as deviance
cctfit$BUGSoutput$pD <- var(cctfit$BUGSoutput$sims.list$deviance[cctfit$BUGSoutput$sims.list$deviance != Inf,1],na.rm=TRUE)/2  #pD, variance of the deviance, divided by 2
cctfit$BUGSoutput$DIC <- cctfit$BUGSoutput$mean$deviance + cctfit$BUGSoutput$pD
}else{
cctfit$BUGSoutput$mean$deviance <- mean(cctfit$BUGSoutput$sims.list$deviance) #Dbar, also known as deviance
cctfit$BUGSoutput$pD <- var(cctfit$BUGSoutput$sims.list$deviance[,1])/2  #pD, variance of the deviance, divided by 2
cctfit$BUGSoutput$DIC <- cctfit$BUGSoutput$mean$deviance + cctfit$BUGSoutput$pD
}

}

if(cctfit$whmodel=="CRM"){
if(cctfit$V == 1){
cctfit$V <- 1;
cctfit$BUGSoutput$sims.list$e <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$n));
if(length(dim(cctfit$BUGSoutput$sims.list$T)) < 3){
cctfit$BUGSoutput$sims.list$T <- array(cctfit$BUGSoutput$sims.list[["T"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
if(cctfit$itemdiff == 1){
if(length(dim(cctfit$BUGSoutput$sims.list$lam)) < 3){
cctfit$BUGSoutput$sims.list$lam <- array(cctfit$BUGSoutput$sims.list[["lam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
}
}
if(cctfit$itemdiff == 0){
if(cctfit$V == 1){
cctfit$BUGSoutput$sims.list$lam <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$m,1))
}
if(cctfit$V > 1){
cctfit$BUGSoutput$sims.list$lam <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$m,cctfit$V))
}
}

storefulllik <- 0
if(storefulllik == 1){
cctfit$BUGSoutput$sims.list$tau <- array(-1,c(cctfit$BUGSoutput$n.sims,cctfit$n,cctfit$m))
cctfit$Lik  <- array(NA, c(cctfit$n,cctfit$m,cctfit$BUGSoutput$n.sims));

for(samp in 1:cctfit$BUGSoutput$n.sims){
cctfit$BUGSoutput$sims.list[["tau"]][samp,,] <- cctfit$BUGSoutput$sims.list[["E"]][samp,]*t(1/cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])
cctfit$Lik[,,samp] <-  dnorm(cctfit$data,mean=(cctfit$BUGSoutput$sims.list[["a"]][samp,]*t(cctfit$BUGSoutput$sims.list[["T"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]]))+array(cctfit$BUGSoutput$sims.list[["b"]][samp,],c(cctfit$n,cctfit$m)),sd=cctfit$BUGSoutput$sims.list[["a"]][samp,]*(cctfit$BUGSoutput$sims.list[["tau"]][samp,,]^-.5))
if(cctfit$mval==1){cctfit$Lik[cbind(datob$thena,samp)] <- 1}
}
cctfit$BUGSoutput$sims.list$deviance <- array(-2*apply(log(cctfit$Lik),3,sum),c(cctfit$BUGSoutput$n.sims,1))
}else{
cctfit$LogLik  <- array(NA, c(cctfit$BUGSoutput$n.sims));

for(samp in 1:cctfit$BUGSoutput$n.sims){
tautmp <- cctfit$BUGSoutput$sims.list[["E"]][samp,]*t(1/cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])
Liktmp <-  dnorm(cctfit$data,mean=(cctfit$BUGSoutput$sims.list[["a"]][samp,]*t(cctfit$BUGSoutput$sims.list[["T"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]]))+array(cctfit$BUGSoutput$sims.list[["b"]][samp,],c(cctfit$n,cctfit$m)),sd=cctfit$BUGSoutput$sims.list[["a"]][samp,]*(tautmp^-.5))
if(cctfit$mval==1){Liktmp[datob$thena] <- 1}
cctfit$LogLik[samp] <- sum(log(Liktmp))
}
cctfit$BUGSoutput$sims.list$deviance <- array(-2*(cctfit$LogLik),c(cctfit$BUGSoutput$n.sims,1))
}

cctfit$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "deviance")]]] <- array(cctfit$BUGSoutput$sims.list$deviance,c(cctfit$BUGSoutput$n.keep,cctfit$BUGSoutput$n.chains))
cctfit$BUGSoutput$sims.matrix[,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short == "deviance")]]] <- cctfit$BUGSoutput$sims.list$deviance[,1]
 
if(sum(cctfit$BUGSoutput$sims.list$deviance == Inf) >0){ 
cctfit$BUGSoutput$mean$deviance <- mean(cctfit$BUGSoutput$sims.list$deviance[cctfit$BUGSoutput$sims.list$deviance != Inf,1],na.rm=TRUE) #Dbar, also known as deviance
cctfit$BUGSoutput$pD <- var(cctfit$BUGSoutput$sims.list$deviance[cctfit$BUGSoutput$sims.list$deviance != Inf,1],na.rm=TRUE)/2  #pD, variance of the deviance, divided by 2
cctfit$BUGSoutput$DIC <- cctfit$BUGSoutput$mean$deviance + cctfit$BUGSoutput$pD
}else{
cctfit$BUGSoutput$mean$deviance <- mean(cctfit$BUGSoutput$sims.list$deviance) #Dbar, also known as deviance
cctfit$BUGSoutput$pD <- var(cctfit$BUGSoutput$sims.list$deviance[,1])/2  #pD, variance of the deviance, divided by 2
cctfit$BUGSoutput$DIC <- cctfit$BUGSoutput$mean$deviance + cctfit$BUGSoutput$pD
}

}

message(paste("DIC : ",round(cctfit$BUGSoutput$DIC,2),"   pD : ",round(cctfit$BUGSoutput$pD,2),sep=""))

if(gui==1){
cctfit <<- cctfit
tkconfigure(plotresults.but, state="normal") 
tkconfigure(doppc.but, state="normal") 
tkconfigure(exportresults.but, state="normal") 
}else{
cctfit <<- cctfit
return(cctfit)
}

}


######################
#Function for the 'Plot Results' Button
#- Creates a sophisticated plot of the posterior results, which includes
#    the posterior means of each parameters and their highest density intervals (HDIs, see Kruschke 2011)
#- Denotes a different symbol for the parameters pertaining to each culture
#- This function is also called during the file export, and .eps and .jpeg's of the plot are saved
######################
plotresultsfuncbutton <- function() {
plotresultsfunc(cctfit=cctfit,gui=1)
}

plotresultsfunc <- function(cctfit,saveplots=0,savedir=0,gui=0) {

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
cctfit$respmem <- apply(cctfit$BUGSoutput$sims.list$e[,],2,Mode)

if(saveplots==1){jpeg(file.path(gsub(".Rdata","results.jpg",savedir)),width = 6, height = 6, units = "in", pointsize = 12,quality=100,res=400)}
if(saveplots==2){postscript(file=file.path(gsub(".Rdata","results.eps",savedir)), onefile=FALSE, horizontal=FALSE, width = 6, height = 6, paper="special", family="Times")}

if(cctfit$whmodel=="GCM"){
par(oma=c(0,0,0,0),mar=c(4,4,3,1),mgp=c(2.25,.75,0),mfrow=c(2,2))

sym <- c(21,22, 23, 24, 25, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14) # possible symbols: 

if(cctfit$V == 1){
if(length(dim(cctfit$BUGSoutput$sims.list[["Z"]])) < 3){
cctfit$BUGSoutput$sims.list[["Z"]] <- array(cctfit$BUGSoutput$sims.list[["Z"]],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
plot(1:cctfit$m,rep(NA,cctfit$m),main=expression(paste("Item Truth (",Z[vk],")")),xlab="Item",ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=21,bg="white")
for(i in 1:dim(cctfit$BUGSoutput$sims.list[["Z"]])[3]){
points(apply(cctfit$BUGSoutput$sims.list[["Z"]][,,i],c(2),mean),xlab=expression(paste("Item Truth (",Z[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg="white")
}}else{
plot(1:cctfit$m,rep(NA,cctfit$m),main=expression(paste("Item Truth (",Z[vk],") Per Culture")),xlab="Item",ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=21,bg="white")
for(i in 1:dim(cctfit$BUGSoutput$sims.list[["Z"]])[3]){
points(apply(cctfit$BUGSoutput$sims.list[["Z"]][,,i],c(2),mean),xlab=expression(paste("Item Truth (",Z[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg="white")
}
}

if(cctfit$itemdiff==1){
if(cctfit$V == 1){
if(length(dim(cctfit$BUGSoutput$sims.list[["lam"]])) < 3){
cctfit$BUGSoutput$sims.list[["lam"]] <- array(cctfit$BUGSoutput$sims.list[["lam"]],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
hditmp <- apply(cctfit$BUGSoutput$sims.list[["lam"]][,,1],2,hdi)
plot(1:cctfit$m,rep(NA,cctfit$m),main=expression(paste("Item Difficulty (",lambda[vk],")")),xlab="Item",ylim=c(0,1),ylab="Posterior Mean Value",las=1,pch=21,bg="white")
for(i in 1:dim(cctfit$BUGSoutput$sims.list[["lam"]])[3]){
points(apply(cctfit$BUGSoutput$sims.list[["lam"]][,,i],c(2),mean),xlab=expression(paste("Item Difficulty (",lambda[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg="white")
}
segments(1:cctfit$m,hditmp[1,], 1:cctfit$m, hditmp[2,])
arrows(1:cctfit$m,hditmp[1,], 1:cctfit$m, hditmp[2,],code=3,angle=90,length=.025)
}else{plot(1:cctfit$m,rep(NA,cctfit$m),main=expression(paste("Item Difficulty (",lambda[vk],") Per Culture")),xlab="Item",ylim=c(0,1),ylab="Posterior Mean Value",las=1,pch=21,bg="white")
for(i in 1:dim(cctfit$BUGSoutput$sims.list[["lam"]])[3]){
points(apply(cctfit$BUGSoutput$sims.list[["lam"]][,,i],c(2),mean),xlab=expression(paste("Item Difficulty (",lambda[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg="white")
}
}
}else{
plot(c(1:cctfit$m),rep(.5,cctfit$m),main=expression(paste("Item Difficulty (",lambda[vk],")")),xlab="Item",ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[1],bg="white",col="grey")
points(c(par("usr")[1],par("usr")[2]),c(par("usr")[3],par("usr")[4]),type="l",col="grey")
points(c(par("usr")[1],par("usr")[2]),c(par("usr")[4],par("usr")[3]),type="l",col="grey")
}

plot(cctfit$BUGSoutput$mean$th,main=expression(paste("Respondent Competency (",theta[i],")")),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[cctfit$respmem],bg="white")
hditmp <- apply(cctfit$BUGSoutput$sims.list$th,2,hdi)
segments(1:cctfit$n,hditmp[1,], 1:cctfit$n, hditmp[2,])
arrows(1:cctfit$n,hditmp[1,], 1:cctfit$n, hditmp[2,],code=3,angle=90,length=.025)

plot(cctfit$BUGSoutput$mean$g,main=expression(paste("Respondent Guessing Bias (",g[i],")")),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[cctfit$respmem],bg="white")
hditmp <- apply(cctfit$BUGSoutput$sims.list$g,2,hdi)
segments(1:cctfit$n,hditmp[1,], 1:cctfit$n, hditmp[2,])
arrows(1:cctfit$n,hditmp[1,], 1:cctfit$n, hditmp[2,],code=3,angle=90,length=.025)
}

if(cctfit$whmodel=="LTRM"){
invlogit <- function(x){x <- 1 / (1 + exp(-x)); return(x)}

par(oma=c(0,0,0,0),mar=c(4,4,3,1),mgp=c(2.25,.75,0),mfrow=c(2,2))

sym <- c(21,22, 23, 24, 25, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14) # possible symbols: 

if(cctfit$C==2){useinvlogit <- 1}else{useinvlogit <- 0}

if(useinvlogit == 1){
tmp1 <- cctfit$BUGSoutput$sims.list[["T"]]; cctfit$BUGSoutput$sims.list[["T"]] <- invlogit(cctfit$BUGSoutput$sims.list[["T"]])
tmp2 <- cctfit$BUGSoutput$sims.list[["gam"]]; cctfit$BUGSoutput$sims.list[["gam"]] <- invlogit(cctfit$BUGSoutput$sims.list[["gam"]])
tmp3 <- cctfit$BUGSoutput$sims.list[["b"]]; cctfit$BUGSoutput$sims.list[["b"]] <- invlogit(cctfit$BUGSoutput$sims.list[["b"]])
}

if(cctfit$V == 1){
newm <- ceiling(cctfit$m*1.092)
if(length(dim(cctfit$BUGSoutput$sims.list[["T"]])) <3){
cctfit$BUGSoutput$sims.list[["T"]] <- array(cctfit$BUGSoutput$sims.list[["T"]],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))
}
if(length(dim(cctfit$BUGSoutput$sims.list[["gam"]])) <3){
cctfit$BUGSoutput$sims.list[["gam"]] <- array(cctfit$BUGSoutput$sims.list[["gam"]],c(cctfit$BUGSoutput$n.sims,cctfit$C-1,1))
}
hditmp <- apply(cctfit$BUGSoutput$sims.list[["T"]][,,1],2,hdi)
if(useinvlogit == 1){
plot(1:newm,rep(NA,newm),main=expression(paste("Item Truth (",T[vk],") and Thresholds (",gamma[vc],")")),xlab="Item",ylim=c(0,1),ylab="Posterior Mean Value",las=1,pch=21,bg="white",axes=FALSE)
}else{plot(1:newm,rep(NA,newm),main=expression(paste("Item Truth (",T[vk],") and Thresholds (",gamma[vc],")")),xlab="Item",ylim=c(min(hditmp,min(apply(cctfit$BUGSoutput$sims.list[["gam"]],c(2,3),mean))),max(hditmp,max(apply(cctfit$BUGSoutput$sims.list[["gam"]],c(2,3),mean)))),ylab="Posterior Mean Value",las=1,pch=21,bg="white",axes=FALSE)
}
for(i in 1:dim(cctfit$BUGSoutput$sims.list[["T"]])[3]){
points(apply(cctfit$BUGSoutput$sims.list[["T"]][,,i],c(2),mean),xlab=expression(paste("Item Truth (",T[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg="white")

if(cctfit$C == 2){
text(newm,mean(cctfit$BUGSoutput$sims.list[["gam"]][,,i]),labels=sapply(c(1:(cctfit$C-1)),function(x) as.expression(substitute(list(gamma[x]),list(x=x)))))
segments(1,mean(cctfit$BUGSoutput$sims.list[["gam"]][,,i]),max(cctfit$m+1,newm-2),mean(cctfit$BUGSoutput$sims.list[["gam"]][,,i]),lty=2);
}
else{
text(newm,apply(cctfit$BUGSoutput$sims.list[["gam"]][,,i],c(2),mean),labels=sapply(c(1:(cctfit$C-1)),function(x) as.expression(substitute(list(gamma[x]),list(x=x)))))
segments(1,apply(cctfit$BUGSoutput$sims.list[["gam"]][,,i],c(2),mean),max(cctfit$m+1,newm-2),apply(cctfit$BUGSoutput$sims.list[["gam"]][,,i],c(2),mean),lty=2);
}
segments(1:cctfit$m,hditmp[1,], 1:cctfit$m, hditmp[2,])
arrows(1:cctfit$m,hditmp[1,], 1:cctfit$m, hditmp[2,],code=3,angle=90,length=.025)
box()
axis(2, labels = TRUE,las=1)
axis(side = 1,labels=TRUE)
axis(side = 1, at = newm, labels = expression(gamma[c]) )
}}else{

newm <- ceiling(cctfit$m*1.14)

if(cctfit$C == 2){

if(useinvlogit == 1){
plot(1:newm,rep(NA,newm),main=expression(paste("Item Truth (",T[vk],") and Thresholds (",gamma[vc],") Per Culture")),xlab="Item",ylim=c(0,1),ylab="Posterior Mean Value",las=1,pch=21,bg="white",axes=FALSE)
}else{plot(1:newm,rep(NA,newm),main=expression(paste("Item Truth (",T[vk],") and Thresholds (",gamma[vc],") Per Culture")),xlab="Item",ylim=c(min(min(apply(cctfit$BUGSoutput$sims.list[["T"]],c(2,3),mean))
,min(apply(cctfit$BUGSoutput$sims.list[["gam"]],c(3),mean))),max(max(apply(cctfit$BUGSoutput$sims.list[["T"]],c(2,3),mean))
,max(apply(cctfit$BUGSoutput$sims.list[["gam"]],c(3),mean)))),ylab="Posterior Mean Value",las=1,pch=21,bg="white",axes=FALSE)
}

for(i in 1:dim(cctfit$BUGSoutput$sims.list[["T"]])[3]){
points(apply(cctfit$BUGSoutput$sims.list[["T"]][,,i],c(2),mean),xlab=expression(paste("Item Truth (",T[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg="white")
text(cctfit$m+(((newm-cctfit$m)/dim(cctfit$BUGSoutput$sims.list[["T"]])[3])*i),mean(cctfit$BUGSoutput$sims.list[["gam"]][,,i]),labels=sapply(c(1:(cctfit$C-1)),function(x) as.expression(substitute(list(gamma[x]),list(x=x)))))
}
}
else{

if(useinvlogit == 1){
plot(1:newm,rep(NA,newm),main=expression(paste("Item Truth (",T[vk],") and Thresholds (",gamma[vc],") Per Culture")),xlab="Item",ylim=c(0,1),ylab="Posterior Mean Value",las=1,pch=21,bg="white",axes=FALSE)
}else{
plot(1:newm,rep(NA,newm),main=expression(paste("Item Truth (",T[vk],") and Thresholds (",gamma[vc],") Per Culture")),xlab="Item",ylim=c(min(min(apply(cctfit$BUGSoutput$sims.list[["T"]],c(2,3),mean))
,min(apply(cctfit$BUGSoutput$sims.list[["gam"]][,,],c(2,3),mean))),max(max(apply(cctfit$BUGSoutput$sims.list[["T"]],c(2,3),mean))
,max(apply(cctfit$BUGSoutput$sims.list[["gam"]][,,],c(2,3),mean)))),ylab="Posterior Mean Value",las=1,pch=21,bg="white",axes=FALSE)
}

for(i in 1:dim(cctfit$BUGSoutput$sims.list[["T"]])[3]){
points(apply(cctfit$BUGSoutput$sims.list[["T"]][,,i],c(2),mean),xlab=expression(paste("Item Truth (",T[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg="white")
text(cctfit$m+(((newm-cctfit$m)/dim(cctfit$BUGSoutput$sims.list[["T"]])[3])*i),apply(cctfit$BUGSoutput$sims.list[["gam"]][,,i],c(2),mean),labels=sapply(c(1:(cctfit$C-1)),function(x) as.expression(substitute(list(gamma[x]),list(x=x)))))
}
}
box()
axis(2, labels = TRUE,las=1)
axis(side=1,at=axTicks(side=1)[axTicks(side=1)<= cctfit$m])
axis(side = 1, at = cctfit$m+(((newm-cctfit$m)/dim(cctfit$BUGSoutput$sims.list[["T"]])[3])*(1:2)), 
labels = sapply(1:2,function(x) as.expression(substitute(list(gamma[x*c]),list(x=x)))) )
}

if(cctfit$itemdiff==1){
if(cctfit$V == 1){
if(length(dim(cctfit$BUGSoutput$sims.list[["lam"]])) < 3){
cctfit$BUGSoutput$sims.list[["lam"]] <- array(cctfit$BUGSoutput$sims.list[["lam"]],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
hditmp <- apply(cctfit$BUGSoutput$sims.list[["lam"]][,,1],2,hdi)
plot(1:cctfit$m,rep(NA,cctfit$m),main=expression(paste("Item Difficulty (",lambda[vk],")")),xlab="Item",ylim=c(min(hditmp),max(hditmp)),ylab="Posterior Mean Value",las=1,pch=21,bg="white")
for(i in 1:dim(cctfit$BUGSoutput$sims.list[["lam"]])[3]){
points(apply(cctfit$BUGSoutput$sims.list[["lam"]][,,i],c(2),mean),xlab=expression(paste("Item Difficulty (",lambda[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg="white")
}
segments(1:cctfit$m,hditmp[1,], 1:cctfit$m, hditmp[2,])
arrows(1:cctfit$m,hditmp[1,], 1:cctfit$m, hditmp[2,],code=3,angle=90,length=.025)
}else{plot(1:cctfit$m,rep(NA,cctfit$m),main=expression(paste("Item Difficulty (",lambda[vk],") Per Culture")),xlab="Item",ylim=c(min(apply(cctfit$BUGSoutput$sims.list[["lam"]],c(2,3),mean)),max(apply(cctfit$BUGSoutput$sims.list[["lam"]],c(2,3),mean))),ylab="Posterior Mean Value",las=1,pch=21,bg="white")
for(i in 1:dim(cctfit$BUGSoutput$sims.list[["lam"]])[3]){
points(apply(cctfit$BUGSoutput$sims.list[["lam"]][,,i],c(2),mean),xlab=expression(paste("Item Difficulty (",lambda[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg="white")
}
}

}else{
plot(c(1:cctfit$m),rep(.5,cctfit$m),main=expression(paste("Item Difficulty (",lambda[vk],")")),xlab="Item",ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[1],bg="white",col="grey")
points(c(par("usr")[1],par("usr")[2]),c(par("usr")[3],par("usr")[4]),type="l",col="grey")
points(c(par("usr")[1],par("usr")[2]),c(par("usr")[4],par("usr")[3]),type="l",col="grey")
}

hditmp <- apply(cctfit$BUGSoutput$sims.list$E,2,hdi)
plot(cctfit$BUGSoutput$mean$E,main=expression(paste("Respondent Competency (",E[i],")")),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(0,max(hditmp)),pch=sym[cctfit$respmem],bg="white")
segments(1:cctfit$n,hditmp[1,], 1:cctfit$n, hditmp[2,])
arrows(1:cctfit$n,hditmp[1,], 1:cctfit$n, hditmp[2,],code=3,angle=90,length=.025)

if(useinvlogit == 1){
plot(invlogit(cctfit$BUGSoutput$mean$b)-.5,main=expression(paste("Respondent Shift Bias (",b[i],")")),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(-.5,.5),pch=sym[cctfit$respmem],bg="white")
hditmp <- apply(cctfit$BUGSoutput$sims.list$b-.5,2,hdi) # because this is invlogit transformed still
segments(1:cctfit$n,hditmp[1,], 1:cctfit$n, hditmp[2,])
arrows(1:cctfit$n,hditmp[1,], 1:cctfit$n, hditmp[2,],code=3,angle=90,length=.025)
}else{
plot(-5000, -5000, main=expression(paste("Category Usage Bias (",a[i]," and ",b[i],")")),xlim=c(min(-1.6,-abs(cctfit$BUGSoutput$mean$b)),max(1.6,cctfit$BUGSoutput$mean$b)), ylim=c(min(-.6,cctfit$BUGSoutput$mean$a),min(3,max(2.4,cctfit$BUGSoutput$mean$a))),xlab=expression(paste("Respondent Shift Bias (",b[i],")")),
ylab=expression(paste("Respondent Scale Bias (",a[i],")")),las=1,pch=21)
points(cbind(cctfit$BUGSoutput$mean$b,cctfit$BUGSoutput$mean$a),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(0,max(hditmp)),pch=sym[cctfit$respmem],bg="white")
segments(-1,1,1,1)
segments(0,0,0,2)
text(0,2.25,"Middle Categories",cex=.8)
text(0,-.25,"Outer Categories",cex=.8)
text(-1.3,1,"Left \n Categories",cex=.8)
text(1.3,1,"Right \n Categories",cex=.8)
}

if(useinvlogit == 1){
cctfit$BUGSoutput$sims.list[["T"]] <- tmp1 
cctfit$BUGSoutput$sims.list[["gam"]] <- tmp2
cctfit$BUGSoutput$sims.list[["b"]] <- tmp3
}

}

if(cctfit$whmodel=="CRM"){

par(oma=c(0,0,0,0),mar=c(4,4,3,1),mgp=c(2.25,.75,0),mfrow=c(2,2))

sym <- c(21,22, 23, 24, 25, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14) # possible symbols: 

if(cctfit$V == 1){
if(length(dim(cctfit$BUGSoutput$sims.list[["T"]])) < 3){
cctfit$BUGSoutput$sims.list[["T"]] <- array(cctfit$BUGSoutput$sims.list[["T"]],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
hditmp <- apply(cctfit$BUGSoutput$sims.list[["T"]][,,1],2,hdi)
plot(1:cctfit$m,rep(NA,cctfit$m),main=expression(paste("Item Truth (",T[vk],")")),xlab="Item",ylim=c(min(hditmp),max(hditmp)),ylab="Posterior Mean Value",las=1,pch=21,bg="white")
}else{plot(1:cctfit$m,rep(NA,cctfit$m),main=expression(paste("Item Truth (",T[vk],") Per Culture")),xlab="Item",ylim=c(min(apply(cctfit$BUGSoutput$sims.list[["T"]],c(2,3),mean)),max(apply(cctfit$BUGSoutput$sims.list[["T"]],c(2,3),mean))),ylab="Posterior Mean Value",las=1,pch=21,bg="white")}
for(i in 1:dim(cctfit$BUGSoutput$sims.list[["T"]])[3]){
points(apply(cctfit$BUGSoutput$sims.list[["T"]][,,i],c(2),mean),xlab=expression(paste("Item Truth (",T[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg="white")
}
if(cctfit$V == 1){
segments(1:cctfit$m,hditmp[1,], 1:cctfit$m, hditmp[2,])
arrows(1:cctfit$m,hditmp[1,], 1:cctfit$m, hditmp[2,],code=3,angle=90,length=.025)
}

if(cctfit$itemdiff==1){
if(cctfit$V == 1){
if(length(dim(cctfit$BUGSoutput$sims.list[["lam"]])) < 3){
cctfit$BUGSoutput$sims.list[["lam"]] <- array(cctfit$BUGSoutput$sims.list[["lam"]],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
hditmp <- apply(cctfit$BUGSoutput$sims.list[["lam"]][,,1],2,hdi)
plot(1:cctfit$m,rep(NA,cctfit$m),main=expression(paste("Item Difficulty (",lambda[vk],")")),xlab="Item",ylim=c(min(hditmp),max(hditmp)),ylab="Posterior Mean Value",las=1,pch=21,bg="white")
for(i in 1:dim(cctfit$BUGSoutput$sims.list[["lam"]])[3]){
points(apply(cctfit$BUGSoutput$sims.list[["lam"]][,,i],c(2),mean),xlab=expression(paste("Item Difficulty (",lambda[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg="white")
}
segments(1:cctfit$m,hditmp[1,], 1:cctfit$m, hditmp[2,])
arrows(1:cctfit$m,hditmp[1,], 1:cctfit$m, hditmp[2,],code=3,angle=90,length=.025)
}else{plot(1:cctfit$m,rep(NA,cctfit$m),main=expression(paste("Item Difficulty (",lambda[vk],") Per Culture")),xlab="Item",ylim=c(min(apply(cctfit$BUGSoutput$sims.list[["lam"]],c(2,3),mean)),max(apply(cctfit$BUGSoutput$sims.list[["lam"]],c(2,3),mean))),ylab="Posterior Mean Value",las=1,pch=21,bg="white")
for(i in 1:dim(cctfit$BUGSoutput$sims.list[["lam"]])[3]){
points(apply(cctfit$BUGSoutput$sims.list[["lam"]][,,i],c(2),mean),xlab=expression(paste("Item Difficulty (",lambda[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg="white")
}
}
}else{
plot(c(1:cctfit$m),rep(.5,cctfit$m),main=expression(paste("Item Difficulty (",lambda[vk],")")),xlab="Item",ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[1],bg="white",col="grey")
points(c(par("usr")[1],par("usr")[2]),c(par("usr")[3],par("usr")[4]),type="l",col="grey")
points(c(par("usr")[1],par("usr")[2]),c(par("usr")[4],par("usr")[3]),type="l",col="grey")
}

hditmp <- apply(cctfit$BUGSoutput$sims.list$E,2,hdi)
plot(cctfit$BUGSoutput$mean$E,main=expression(paste("Respondent Competency (",E[i],")")),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(0,max(hditmp)),pch=sym[cctfit$respmem],bg="white")
segments(1:cctfit$n,hditmp[1,], 1:cctfit$n, hditmp[2,])
arrows(1:cctfit$n,hditmp[1,], 1:cctfit$n, hditmp[2,],code=3,angle=90,length=.025)

plot(-5000, -5000, main=expression(paste("Appraisal Bias (",a[i]," and ",b[i],")")),xlim=c(min(-2.8,-abs(cctfit$BUGSoutput$mean$b)),max(2.8,cctfit$BUGSoutput$mean$b)), ylim=c(min(-.6,cctfit$BUGSoutput$mean$a),min(3,max(2.4,cctfit$BUGSoutput$mean$a))),xlab=expression(paste("Respondent Shift Bias (",b[i],")")),
ylab=expression(paste("Respondent Scale Bias (",a[i],")")),las=1,pch=21)
points(cbind(cctfit$BUGSoutput$mean$b,cctfit$BUGSoutput$mean$a),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(0,max(hditmp)),pch=sym[cctfit$respmem],bg="white")
segments(-2,1,2,1)
segments(0,0,0,2)
text(0,2.25,"Middle Values",cex=.8)
text(0,-.25,"Outer Values",cex=.8)
text(-2.5,1,"Left \n Values",cex=.8)
text(2.5,1,"Right \n Values",cex=.8)
}

if(saveplots==1 || saveplots==2){dev.off()}

saveplots <- 0
}

######################
#Function for the 'Run Checks' Button
#- Calculates the posterior predictive data for the model that was applied
#- Calculates 2 important posterior predictive checks and plots them
#- The VDI check percentiles are reported in the R console, and saved to cctfit$VDIperc
#- Performs factor analysis using fa() and polychoric() (if selected)
#- This function is also called in the file export and saves .eps and .jpegs of the plot
######################
ppcfuncbutton <- function(){

if(cctfit$whmodel=="LTRM" && cctfit$checksrun == 1){ 
if(cctfit$polycor != as.numeric(tclvalue(polyvar))){
cctfit$checksrun <- 0; cctfit <<- cctfit
}
}
ppcfunc(cctfit=cctfit,gui=1,polych=as.numeric(tclvalue(polyvar)))
}

ppcfunc <- function(cctfit,saveplots=0,savedir=0,gui=0,polych=0) {

if(cctfit$checksrun == 0){
message("\n ...One moment, calculating posterior predictive checks")

usesubset <- 1

if(usesubset == 1){
subsetsize <- min(500,cctfit$BUGSoutput$n.sims)
indices <- sample(cctfit$BUGSoutput$n.sims,min(subsetsize,cctfit$BUGSoutput$n.sims))

if(cctfit$whmodel=="GCM"){
if(cctfit$V == 1){
cctfit$V <- 1;
cctfit$BUGSoutput$sims.list$e <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$n));
if(length(dim(cctfit$BUGSoutput$sims.list$Z)) < 3){
cctfit$BUGSoutput$sims.list$Z <- array(cctfit$BUGSoutput$sims.list[["Z"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
if(cctfit$itemdiff == 1){
if(length(dim(cctfit$BUGSoutput$sims.list$lam)) < 3){
cctfit$BUGSoutput$sims.list$lam <- array(cctfit$BUGSoutput$sims.list[["lam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
}
}
if(cctfit$itemdiff == 0){
if(cctfit$V == 1){
cctfit$BUGSoutput$sims.list$lam <- array(.5,c(cctfit$BUGSoutput$n.sims,cctfit$m,1))
}
if(cctfit$V > 1){
cctfit$BUGSoutput$sims.list$lam <- array(.5,c(cctfit$BUGSoutput$n.sims,cctfit$m,cctfit$V))
} }

cctfit$ppY <- array(NA, c(cctfit$n,cctfit$m,subsetsize));

for(samp in indices){
tautmp <- (cctfit$BUGSoutput$sims.list[["th"]][samp,]*t(1-cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])) / 
( (cctfit$BUGSoutput$sims.list[["th"]][samp,]*t(1-cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])) + 
((1-cctfit$BUGSoutput$sims.list[["th"]][samp,])*t(cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])) ) 

cctfit$ppY[,,which(indices == samp)] <- matrix(rbinom((cctfit$n*cctfit$m),1,
(t(cctfit$BUGSoutput$sims.list[["Z"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])*tautmp ) +
 t(t(1-tautmp)%*%diag(cctfit$BUGSoutput$sims.list[["g"]][samp,])) ), cctfit$n,cctfit$m)
}

if(cctfit$mval==1){
for(i in 1:dim(cctfit$datob$thena)[1]){cctfit$data[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2]] <- Mode(cctfit$ppY[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2],])}
cctfit$MVest <- cbind(cctfit$datob$thena,cctfit$data[cctfit$datob$thena])
colnames(cctfit$MVest) <- c("Pers","Item","Resp")
}

cctfit$mean$ppY <- rowMeans(cctfit$ppY[,,],dims=2)

leigv <- 12
options(warn=-3)
ind <- indices 
eigv <- matrix(-1,length(ind),leigv)
tmp <- apply(cctfit$ppY,c(1,3),function(x) t(x))
tmp2 <- apply(tmp[,,],3,function(x) cor(x))
tmp3 <- array(tmp2,c(cctfit$n,cctfit$n,subsetsize))
for(i in 1:length(ind)){
suppressMessages(try(eigv[i,] <- fa(tmp3[,,i])$values[1:leigv],silent=TRUE))}
wch <- -which(eigv[,1] == -1); 
if(length(wch)==0){cctfit$ppeig <-  eigv}else{cctfit$ppeig <-  eigv[-which(eigv[,1] == -1),]}
cctfit$dateig <-  suppressMessages(fa(cor(t(cctfit$data)))$values[1:leigv])
options(warn=0)

if(cctfit$V==1){
varvec <- matrix(-1,length(ind),dim(cctfit$ppY)[2])
vdi <- matrix(-1,dim(cctfit$ppY)[3],1)
for(i in 1:length(ind)){
varvec[i,] <- apply(cctfit$ppY[,,i],2,var)
}
vdi <- apply(varvec,1,var)
cctfit$ppVDI <- vdi
cctfit$datVDI <- var(apply(cctfit$data,2,var))
}

options(warn=-3)
if(cctfit$V>1){
varvec <- array(NA,c(length(ind),dim(cctfit$ppY)[2],cctfit$V))
vdi <- matrix(NA,dim(cctfit$ppY)[3],cctfit$V)
cctfit$datVDI  <- array(-1,c(1,cctfit$V))
for(v in 1:cctfit$V){
for(i in 1:length(ind)){
if(sum(cctfit$BUGSoutput$sims.list$e[ind[i],] == v)>1){
suppressMessages(try(varvec[i,,v] <- apply(cctfit$ppY[cctfit$BUGSoutput$sims.list$e[ind[i],] == v,,i],2,var),silent=TRUE))
}else{
suppressMessages(try(varvec[i,,v] <- var(cctfit$ppY[cctfit$BUGSoutput$sims.list$e[ind[i],] == v,,i]),silent=TRUE))
}
}
if(sum(cctfit$respmem == v)>1){
cctfit$datVDI[v] <- var(apply(cctfit$data[cctfit$respmem == v,],2,var))
}else{
cctfit$datVDI[v] <- var(cctfit$data[cctfit$respmem == v,])
}
}
vdi <- apply(varvec,c(1,3),var)
cctfit$ppVDI <- vdi

options(warn=0)
if(sum(apply(cctfit$ppVDI,2,function(x) all(x == 0,na.rm=TRUE)))>=1){
for(i in which(apply(cctfit$ppVDI,2,function(x) all(x == 0)))){
cctfit$ppVDI[,i] <- NA
}
}
for(v in 1:dim(cctfit$ppVDI)[2]){
if(identical(unique(cctfit$ppVDI[,v]),c(NA,0)) || identical(unique(cctfit$ppVDI[,v]),c(0,NA))){cctfit$ppVDI[,v] <- NA}
}
}
rm(vdi, varvec)
rm(eigv,leigv,tmp,tmp2,tmp3,ind);

cctfit$checksrun <- 1
}

if(cctfit$whmodel=="LTRM"){
if(cctfit$V == 1){
cctfit$V <- 1;
cctfit$BUGSoutput$sims.list$e <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$n));
if(length(dim(cctfit$BUGSoutput$sims.list[["T"]])) < 3){
cctfit$BUGSoutput$sims.list$T <- array(cctfit$BUGSoutput$sims.list[["T"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
if(length(dim(cctfit$BUGSoutput$sims.list[["gam"]])) < 3){
cctfit$BUGSoutput$sims.list$gam <- array(cctfit$BUGSoutput$sims.list[["gam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$C-1,1))}
if(cctfit$itemdiff == 1){
if(length(dim(cctfit$BUGSoutput$sims.list$lam)) < 3){
cctfit$BUGSoutput$sims.list$lam <- array(cctfit$BUGSoutput$sims.list[["lam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
}
}
if(cctfit$itemdiff == 0){
if(cctfit$V == 1){
cctfit$BUGSoutput$sims.list$lam <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$m,1))
}
if(cctfit$V > 1){
cctfit$BUGSoutput$sims.list$lam <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$m,cctfit$V))
}
}
cctfit$ppY <- array(NA, c(cctfit$n,cctfit$m,subsetsize));
cctfit$ppX <- array(NA, c(cctfit$n,cctfit$m,subsetsize));

for(samp in indices){
tautmp <- cctfit$BUGSoutput$sims.list[["E"]][samp,]*t(1/cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])

cctfit$ppX[,,which(indices == samp)] <- matrix(
rnorm((cctfit$n*cctfit$m),t(cctfit$BUGSoutput$sims.list[["T"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]]),(tautmp^(-.5))), 
cctfit$n,cctfit$m)

deltatmp <- cctfit$BUGSoutput$sims.list[["a"]][samp,]*t(cctfit$BUGSoutput$sims.list[["gam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])+cctfit$BUGSoutput$sims.list[["b"]][samp,]

rc <- which(cctfit$ppX[,,which(indices == samp)] < deltatmp[,1],arr.ind=TRUE)
cctfit$ppY[cbind(rc,array(which(indices == samp),dim(rc)[1]))] <- 1 
for(c in 1:(cctfit$C-1)){
rc <- which(cctfit$ppX[,,which(indices == samp)] > deltatmp[,c],arr.ind=TRUE)
cctfit$ppY[cbind(rc,array(which(indices == samp),dim(rc)[1]))] <- (c+1) }

}
cctfit$mean$ppY <- rowMeans(cctfit$ppY[,,],dims=2)

if(cctfit$mval==1){
for(i in 1:dim(cctfit$datob$thena)[1]){cctfit$data[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2]] <- Mode(cctfit$ppY[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2],])}
cctfit$MVest <- cbind(cctfit$datob$thena,cctfit$data[cctfit$datob$thena])
colnames(cctfit$MVest) <- c("Pers","Item","Resp")
}

leigv <- 12
eigv <- matrix(-1,dim(cctfit$ppY)[3],leigv)
options(warn=-3)
ind <- sample(indices,min(cctfit$BUGSoutput$n.sims,250,subsetsize)) #250 is the number of samples
tmp <- apply(cctfit$ppY,c(1,3),function(x) t(x))
if(polych != 0){tmp2 <- apply(tmp[,,],3,function(x) polychoric(x,polycor=TRUE)$rho)}else{
tmp2 <- apply(tmp[,,],3,function(x) cor(x))
}
tmp3 <- array(tmp2,c(cctfit$n,cctfit$n,length(ind))) 
for(i in 1:length(ind)){
suppressMessages(try(eigv[i,] <- fa(tmp3[,,i])$values[1:leigv],silent=TRUE))}
wch <- -which(eigv[,1] == -1); 
if(length(wch)==0){cctfit$ppeig <-  eigv}else{cctfit$ppeig <-  eigv[-which(eigv[,1] == -1),]}
if(polych != 0){
cctfit$dateig <-  suppressMessages(fa(polychoric(t(cctfit$data),polycor=TRUE)$rho)$values[1:leigv])}else{
cctfit$dateig <-  suppressMessages(fa(cor(t(cctfit$data)))$values[1:leigv])
}

options(warn=0)

if(cctfit$V==1){
varvec <- matrix(NA,length(ind),dim(cctfit$ppY)[2])
vdi <- matrix(NA,dim(cctfit$ppY)[3],1)
for(i in 1:length(ind)){
varvec[i,] <- apply(cctfit$ppY[,,i],2,var)
}
vdi <- apply(varvec,1,var)
cctfit$ppVDI <- vdi
cctfit$datVDI <- var(apply(cctfit$data,2,var))
}
options(warn=-3)
if(cctfit$V>1){
varvec <- array(NA,c(length(ind),dim(cctfit$ppY)[2],cctfit$V))
vdi <- matrix(NA,dim(cctfit$ppY)[3],cctfit$V)
cctfit$datVDI  <- array(-1,c(1,cctfit$V))
for(v in 1:cctfit$V){
for(i in 1:length(ind)){
if(sum(cctfit$BUGSoutput$sims.list$e[ind[i],] == v)>1){
suppressMessages(try(varvec[i,,v] <- apply(cctfit$ppY[cctfit$BUGSoutput$sims.list$e[ind[i],] == v,,i],2,var),silent=TRUE))
}else{
suppressMessages(try(varvec[i,,v] <- var(cctfit$ppY[cctfit$BUGSoutput$sims.list$e[ind[i],] == v,,i]),silent=TRUE))
}
}
if(sum(cctfit$respmem == v)>1){
cctfit$datVDI[v] <- var(apply(cctfit$data[cctfit$respmem == v,],2,var))
}else{
cctfit$datVDI[v] <- var(cctfit$data[cctfit$respmem == v,])
}
}
vdi <- apply(varvec,c(1,3),var)
cctfit$ppVDI <- vdi

options(warn=0)
if(sum(apply(cctfit$ppVDI,2,function(x) all(x == 0,na.rm=TRUE)))>=1){
for(i in which(apply(cctfit$ppVDI,2,function(x) all(x == 0)))){
cctfit$ppVDI[,i] <- NA
}
}
for(v in 1:dim(cctfit$ppVDI)[2]){
if(identical(unique(cctfit$ppVDI[,v]),c(NA,0)) || identical(unique(cctfit$ppVDI[,v]),c(0,NA))){cctfit$ppVDI[,v] <- NA}
}
}
cctfit$polycor <- as.numeric(tclvalue(polyvar))

rm(vdi, varvec)
rm(eigv,leigv,tmp,tmp2,tmp3,ind);

cctfit$checksrun <- 1
}

if(cctfit$whmodel=="CRM"){
invlogit <- function(x){x <- 1 / (1 + exp(-x)); return(x)}
logit <- function(x){x <- log(x/(1-x)); return(x)}

if(cctfit$V == 1){
cctfit$V <- 1;
cctfit$BUGSoutput$sims.list$e <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$n));
if(length(dim(cctfit$BUGSoutput$sims.list[["T"]])) < 3){
cctfit$BUGSoutput$sims.list$T <- array(cctfit$BUGSoutput$sims.list[["T"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
if(cctfit$itemdiff == 1){
if(length(dim(cctfit$BUGSoutput$sims.list$lam)) < 3){
cctfit$BUGSoutput$sims.list$lam <- array(cctfit$BUGSoutput$sims.list[["lam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
}
}
if(cctfit$itemdiff == 0){
if(cctfit$V == 1){
cctfit$BUGSoutput$sims.list$lam <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$m,1))
}
if(cctfit$V > 1){
cctfit$BUGSoutput$sims.list$lam <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$m,cctfit$V))
}
}

cctfit$ppY <- array(NA, c(cctfit$n,cctfit$m,subsetsize));
cctfit$ppX <- array(NA, c(cctfit$n,cctfit$m,subsetsize));

for(samp in indices){
tautmp <- cctfit$BUGSoutput$sims.list[["E"]][samp,]*t(1/cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])

cctfit$ppY[,,which(indices == samp)] <- matrix(
rnorm((cctfit$n*cctfit$m),(cctfit$BUGSoutput$sims.list[["a"]][samp,]*t(cctfit$BUGSoutput$sims.list[["T"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]]))+matrix(rep(cctfit$BUGSoutput$sims.list[["b"]][samp,],cctfit$m),cctfit$n,cctfit$m),
(cctfit$BUGSoutput$sims.list[["a"]][samp,]*(tautmp^(-.5)))), 
cctfit$n,cctfit$m)

cctfit$ppX[,,which(indices == samp)] <- (array(rep(1/cctfit$BUGSoutput$sims.list[["a"]][samp,],times=cctfit$m),c(cctfit$n,cctfit$m))*(cctfit$ppY[,,which(indices == samp)]))-matrix(rep(cctfit$BUGSoutput$sims.list[["b"]][samp,],times=cctfit$m),c(cctfit$n,cctfit$m))
}
cctfit$mean$ppY <- rowMeans(cctfit$ppY[,,],dims=2)

if(cctfit$mval==1){
for(i in 1:dim(cctfit$datob$thena)[1]){cctfit$data[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2]] <- mean(cctfit$ppY[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2],])}
cctfit$MVest <- cbind(cctfit$datob$thena,cctfit$data[cctfit$datob$thena])
colnames(cctfit$MVest) <- c("Pers","Item","Resp")
}

leigv <- 12
options(warn=-3)
ind <- indices #500 is the number of samples
eigv <- matrix(-1,length(ind),leigv)
tmp <- apply(cctfit$ppY,c(1,3),function(x) t(x))
tmp2 <- apply(tmp[,,],3,function(x) cor(x))
tmp3 <- array(tmp2,c(cctfit$n,cctfit$n,subsetsize))
for(i in 1:length(ind)){
suppressMessages(try(eigv[i,] <- fa(tmp3[,,i])$values[1:leigv],silent=TRUE))}
wch <- -which(eigv[,1] == -1); 
if(length(wch)==0){cctfit$ppeig <-  eigv}else{cctfit$ppeig <-  eigv[-which(eigv[,1] == -1),]}
cctfit$dateig <-  suppressMessages(fa(cor(t(cctfit$data)))$values[1:leigv])
options(warn=0)

if(cctfit$V==1){
vdi <- matrix(-1,dim(cctfit$ppY)[3],1)
varvec <- array(-1,c(length(ind),dim(cctfit$ppY)[2]))

for(i in 1:length(ind)){
varvec[i,] <- apply(invlogit(cctfit$ppY[,,i]),2,var)  
}
vdi <- apply(varvec,1,var)
cctfit$ppVDI <- apply(varvec,1,var) 

cctfit$datVDI <- var(apply(invlogit(cctfit$data),2,var))  

}
if(cctfit$V>1){
varvec <- array(NA,c(length(ind),dim(cctfit$ppY)[2],cctfit$V))
cctfit$datVDI  <- array(NA,c(1,cctfit$V))
options(warn=-3)
for(v in 1:cctfit$V){
for(i in 1:length(ind)){
if(sum(cctfit$BUGSoutput$sims.list$e[ind[i],] == v)>1){
suppressMessages(try(varvec[i,,v] <- apply(invlogit(cctfit$ppY[cctfit$BUGSoutput$sims.list$e[ind[i],] == v,,i]),2,var),silent=TRUE))
}else{
suppressMessages(try(varvec[i,,v] <- var(invlogit(cctfit$ppY[cctfit$BUGSoutput$sims.list$e[ind[i],] == v,,i])),silent=TRUE)) 
}
}
if(sum(cctfit$respmem == v)>1){
cctfit$datVDI[v] <- var(apply(invlogit(cctfit$data[cctfit$respmem == v,]),2,var)) 
}else{
cctfit$datVDI[v] <- var(invlogit(cctfit$data[cctfit$respmem == v,])) 
}
}
options(warn=0)
vdi <- apply(varvec,c(1,3),function(x) var(x,na.rm=TRUE))
cctfit$ppVDI <- vdi

if(sum(apply(cctfit$ppVDI,2,function(x) all(x == 0,na.rm=TRUE)))>=1){
for(i in which(apply(cctfit$ppVDI,2,function(x) all(x == 0)))){
cctfit$ppVDI[,i] <- NA
}
}
for(v in 1:dim(cctfit$ppVDI)[2]){
if(identical(unique(cctfit$ppVDI[,v]),c(NA,0)) || identical(unique(cctfit$ppVDI[,v]),c(0,NA))){cctfit$ppVDI[,v] <- NA}
}
}

rm(vdi, varvec)
rm(eigv,leigv,tmp,tmp2,tmp3,ind);

cctfit$checksrun <- 1
}


}else{
#########################################
##### Calculate Posterior Predictive Data From All Samples (memory intensive)
#########################################
if(cctfit$whmodel=="GCM"){
if(cctfit$V == 1){
cctfit$V <- 1;
cctfit$BUGSoutput$sims.list$e <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$n));
if(length(dim(cctfit$BUGSoutput$sims.list$Z)) < 3){
cctfit$BUGSoutput$sims.list$Z <- array(cctfit$BUGSoutput$sims.list[["Z"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
if(cctfit$itemdiff == 1){
if(length(dim(cctfit$BUGSoutput$sims.list$lam)) < 3){
cctfit$BUGSoutput$sims.list$lam <- array(cctfit$BUGSoutput$sims.list[["lam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
}
}
if(cctfit$itemdiff == 0){
if(cctfit$V == 1){
cctfit$BUGSoutput$sims.list$lam <- array(.5,c(cctfit$BUGSoutput$n.sims,cctfit$m,1))
}
if(cctfit$V > 1){
cctfit$BUGSoutput$sims.list$lam <- array(.5,c(cctfit$BUGSoutput$n.sims,cctfit$m,cctfit$V))
} }

cctfit$ppY <- array(NA, c(cctfit$n,cctfit$m,cctfit$BUGSoutput$n.sims));
cctfit$BUGSoutput$sims.list$tau <- array(-1,c(cctfit$BUGSoutput$n.sims,cctfit$n,cctfit$m))

for(samp in 1:cctfit$BUGSoutput$n.sims){
cctfit$BUGSoutput$sims.list[["tau"]][samp,,] <- (cctfit$BUGSoutput$sims.list[["th"]][samp,]*t(1-cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])) / 
( (cctfit$BUGSoutput$sims.list[["th"]][samp,]*t(1-cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])) + 
((1-cctfit$BUGSoutput$sims.list[["th"]][samp,])*t(cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])) ) 

cctfit$ppY[,,samp] <- matrix(rbinom((cctfit$n*cctfit$m),1,
(t(cctfit$BUGSoutput$sims.list[["Z"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])*cctfit$BUGSoutput$sims.list[["tau"]][samp,,] ) +
 t(t(1-cctfit$BUGSoutput$sims.list[["tau"]][samp,,])%*%diag(cctfit$BUGSoutput$sims.list[["g"]][samp,])) ), cctfit$n,cctfit$m)
}

if(cctfit$mval==1){
for(i in 1:dim(cctfit$datob$thena)[1]){cctfit$data[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2]] <- Mode(cctfit$ppY[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2],])}
cctfit$MVest <- cbind(cctfit$datob$thena,cctfit$data[cctfit$datob$thena])
colnames(cctfit$MVest) <- c("Pers","Item","Resp")
}

cctfit$mean$ppY <- rowMeans(cctfit$ppY[,,],dims=2)

leigv <- 12
options(warn=-3)
ind <- sample(cctfit$BUGSoutput$n.sims,min(500,cctfit$BUGSoutput$n.sims)) #500 is the number of samples
eigv <- matrix(-1,length(ind),leigv)
tmp <- apply(cctfit$ppY[,,ind],c(1,3),function(x) t(x))
tmp2 <- apply(tmp[,,],3,function(x) cor(x))
tmp3 <- array(tmp2,c(cctfit$n,cctfit$n,cctfit$BUGSoutput$n.sims))
for(i in 1:length(ind)){
suppressMessages(try(eigv[i,] <- fa(tmp3[,,i])$values[1:leigv],silent=TRUE))}
wch <- -which(eigv[,1] == -1); 
if(length(wch)==0){cctfit$ppeig <-  eigv}else{cctfit$ppeig <-  eigv[-which(eigv[,1] == -1),]}
cctfit$dateig <-  suppressMessages(fa(cor(t(cctfit$data)))$values[1:leigv])
options(warn=0)

if(cctfit$V==1){
varvec <- matrix(-1,length(ind),dim(cctfit$ppY)[2])
vdi <- matrix(-1,dim(cctfit$ppY)[3],1)
for(i in 1:length(ind)){
varvec[i,] <- apply(cctfit$ppY[,,i],2,var)
}
vdi <- apply(varvec,1,var)
cctfit$ppVDI <- vdi
cctfit$datVDI <- var(apply(cctfit$data,2,var))
}

options(warn=-3)
if(cctfit$V>1){
varvec <- array(NA,c(length(ind),dim(cctfit$ppY)[2],cctfit$V))
vdi <- matrix(NA,dim(cctfit$ppY)[3],cctfit$V)
cctfit$datVDI  <- array(-1,c(1,cctfit$V))
for(v in 1:cctfit$V){
for(i in 1:length(ind)){
if(sum(cctfit$BUGSoutput$sims.list$e[ind[i],] == v)>1){
suppressMessages(try(varvec[i,,v] <- apply(cctfit$ppY[cctfit$BUGSoutput$sims.list$e[ind[i],] == v,,i],2,var),silent=TRUE))
}else{
suppressMessages(try(varvec[i,,v] <- var(cctfit$ppY[cctfit$BUGSoutput$sims.list$e[ind[i],] == v,,i]),silent=TRUE))
}
}
if(sum(cctfit$respmem == v)>1){
cctfit$datVDI[v] <- var(apply(cctfit$data[cctfit$respmem == v,],2,var))
}else{
cctfit$datVDI[v] <- var(cctfit$data[cctfit$respmem == v,])
}
}
vdi <- apply(varvec,c(1,3),var)
cctfit$ppVDI <- vdi

options(warn=0)
if(sum(apply(cctfit$ppVDI,2,function(x) all(x == 0,na.rm=TRUE)))>=1){
for(i in which(apply(cctfit$ppVDI,2,function(x) all(x == 0)))){
cctfit$ppVDI[,i] <- NA
}
}
for(v in 1:dim(cctfit$ppVDI)[2]){
if(identical(unique(cctfit$ppVDI[,v]),c(NA,0)) || identical(unique(cctfit$ppVDI[,v]),c(0,NA))){cctfit$ppVDI[,v] <- NA}
}
}
rm(vdi, varvec)
rm(eigv,leigv,tmp,tmp2,tmp3,ind);

cctfit$checksrun <- 1
}

if(cctfit$whmodel=="LTRM"){
if(cctfit$V == 1){
cctfit$V <- 1;
cctfit$BUGSoutput$sims.list$e <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$n));
if(length(dim(cctfit$BUGSoutput$sims.list[["T"]])) < 3){
cctfit$BUGSoutput$sims.list$T <- array(cctfit$BUGSoutput$sims.list[["T"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
if(length(dim(cctfit$BUGSoutput$sims.list[["gam"]])) < 3){
cctfit$BUGSoutput$sims.list$gam <- array(cctfit$BUGSoutput$sims.list[["gam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$C-1,1))}
if(cctfit$itemdiff == 1){
if(length(dim(cctfit$BUGSoutput$sims.list$lam)) < 3){
cctfit$BUGSoutput$sims.list$lam <- array(cctfit$BUGSoutput$sims.list[["lam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
}
}
if(cctfit$itemdiff == 0){
if(cctfit$V == 1){
cctfit$BUGSoutput$sims.list$lam <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$m,1))
}
if(cctfit$V > 1){
cctfit$BUGSoutput$sims.list$lam <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$m,cctfit$V))
}
}
cctfit$ppY <- array(NA, c(cctfit$n,cctfit$m,cctfit$BUGSoutput$n.sims));
cctfit$ppX <- array(NA, c(cctfit$n,cctfit$m,cctfit$BUGSoutput$n.sims));
cctfit$ppdelta <- array(NA, c(cctfit$n,cctfit$C-1,cctfit$BUGSoutput$n.sims))

cctfit$BUGSoutput$sims.list$tau <- array(-1,c(cctfit$BUGSoutput$n.sims,cctfit$n,cctfit$m))

for(samp in 1:cctfit$BUGSoutput$n.sims){
cctfit$BUGSoutput$sims.list[["tau"]][samp,,] <- cctfit$BUGSoutput$sims.list[["E"]][samp,]*t(1/cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])

cctfit$ppX[,,samp] <- matrix(
rnorm((cctfit$n*cctfit$m),t(cctfit$BUGSoutput$sims.list[["T"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]]),(cctfit$BUGSoutput$sims.list[["tau"]][samp,,]^(-.5))), 
cctfit$n,cctfit$m)

cctfit$ppdelta[,,samp] <- cctfit$BUGSoutput$sims.list[["a"]][samp,]*t(cctfit$BUGSoutput$sims.list[["gam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])+cctfit$BUGSoutput$sims.list[["b"]][samp,]

rc <- which(cctfit$ppX[,,samp] < cctfit$ppdelta[,1,samp],arr.ind=TRUE)
cctfit$ppY[cbind(rc,array(samp,dim(rc)[1]))] <- 1 
for(c in 1:(cctfit$C-1)){
rc <- which(cctfit$ppX[,,samp] > cctfit$ppdelta[,c,samp],arr.ind=TRUE)
cctfit$ppY[cbind(rc,array(samp,dim(rc)[1]))] <- (c+1) }

}
cctfit$mean$ppY <- rowMeans(cctfit$ppY[,,],dims=2)

if(cctfit$mval==1){
for(i in 1:dim(cctfit$datob$thena)[1]){cctfit$data[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2]] <- Mode(cctfit$ppY[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2],])}
cctfit$MVest <- cbind(cctfit$datob$thena,cctfit$data[cctfit$datob$thena])
colnames(cctfit$MVest) <- c("Pers","Item","Resp")
}

leigv <- 12
eigv <- matrix(-1,dim(cctfit$ppY)[3],leigv)
options(warn=-3)
ind <- sample(cctfit$BUGSoutput$n.sims,min(250,cctfit$BUGSoutput$n.sims)) #250 is the number of samples
tmp <- apply(cctfit$ppY[,,ind],c(1,3),function(x) t(x))
if(polych != 0){tmp2 <- apply(tmp[,,],3,function(x) polychoric(x,polycor=TRUE)$rho)}else{
tmp2 <- apply(tmp[,,],3,function(x) cor(x))
}
tmp3 <- array(tmp2,c(cctfit$n,cctfit$n,length(ind))) 
for(i in 1:length(ind)){
suppressMessages(try(eigv[i,] <- fa(tmp3[,,i])$values[1:leigv],silent=TRUE))}
wch <- -which(eigv[,1] == -1); 
if(length(wch)==0){cctfit$ppeig <-  eigv}else{cctfit$ppeig <-  eigv[-which(eigv[,1] == -1),]}
if(polych != 0){
cctfit$dateig <-  suppressMessages(fa(polychoric(t(cctfit$data),polycor=TRUE)$rho)$values[1:leigv])}else{
cctfit$dateig <-  suppressMessages(fa(cor(t(cctfit$data)))$values[1:leigv])
}

options(warn=0)

if(cctfit$V==1){
varvec <- matrix(NA,length(ind),dim(cctfit$ppY)[2])
vdi <- matrix(NA,dim(cctfit$ppY)[3],1)
for(i in 1:length(ind)){
varvec[i,] <- apply(cctfit$ppY[,,i],2,var)
}
vdi <- apply(varvec,1,var)
cctfit$ppVDI <- vdi
cctfit$datVDI <- var(apply(cctfit$data,2,var))
}
options(warn=-3)
if(cctfit$V>1){
varvec <- array(NA,c(length(ind),dim(cctfit$ppY)[2],cctfit$V))
vdi <- matrix(NA,dim(cctfit$ppY)[3],cctfit$V)
cctfit$datVDI  <- array(-1,c(1,cctfit$V))
for(v in 1:cctfit$V){
for(i in 1:length(ind)){
if(sum(cctfit$BUGSoutput$sims.list$e[ind[i],] == v)>1){
suppressMessages(try(varvec[i,,v] <- apply(cctfit$ppY[cctfit$BUGSoutput$sims.list$e[ind[i],] == v,,i],2,var),silent=TRUE))
}else{
suppressMessages(try(varvec[i,,v] <- var(cctfit$ppY[cctfit$BUGSoutput$sims.list$e[ind[i],] == v,,i]),silent=TRUE))
}
}
if(sum(cctfit$respmem == v)>1){
cctfit$datVDI[v] <- var(apply(cctfit$data[cctfit$respmem == v,],2,var))
}else{
cctfit$datVDI[v] <- var(cctfit$data[cctfit$respmem == v,])
}
}
vdi <- apply(varvec,c(1,3),var)
cctfit$ppVDI <- vdi

options(warn=0)
if(sum(apply(cctfit$ppVDI,2,function(x) all(x == 0,na.rm=TRUE)))>=1){
for(i in which(apply(cctfit$ppVDI,2,function(x) all(x == 0)))){
cctfit$ppVDI[,i] <- NA
}
}
for(v in 1:dim(cctfit$ppVDI)[2]){
if(identical(unique(cctfit$ppVDI[,v]),c(NA,0)) || identical(unique(cctfit$ppVDI[,v]),c(0,NA))){cctfit$ppVDI[,v] <- NA}
}
}
cctfit$polycor <- as.numeric(tclvalue(polyvar))

rm(vdi, varvec)
rm(eigv,leigv,tmp,tmp2,tmp3,ind);

cctfit$checksrun <- 1
}

if(cctfit$whmodel=="CRM"){
invlogit <- function(x){x <- 1 / (1 + exp(-x)); return(x)}
logit <- function(x){x <- log(x/(1-x)); return(x)}

if(cctfit$V == 1){
cctfit$V <- 1;
cctfit$BUGSoutput$sims.list$e <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$n));
if(length(dim(cctfit$BUGSoutput$sims.list[["T"]])) < 3){
cctfit$BUGSoutput$sims.list$T <- array(cctfit$BUGSoutput$sims.list[["T"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
if(cctfit$itemdiff == 1){
if(length(dim(cctfit$BUGSoutput$sims.list$lam)) < 3){
cctfit$BUGSoutput$sims.list$lam <- array(cctfit$BUGSoutput$sims.list[["lam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
}
}
if(cctfit$itemdiff == 0){
if(cctfit$V == 1){
cctfit$BUGSoutput$sims.list$lam <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$m,1))
}
if(cctfit$V > 1){
cctfit$BUGSoutput$sims.list$lam <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$m,cctfit$V))
}
}

cctfit$ppY <- array(NA, c(cctfit$n,cctfit$m,cctfit$BUGSoutput$n.sims));
cctfit$ppX <- array(NA, c(cctfit$n,cctfit$m,cctfit$BUGSoutput$n.sims));

cctfit$BUGSoutput$sims.list$tau <- array(-1,c(cctfit$BUGSoutput$n.sims,cctfit$n,cctfit$m))

for(samp in 1:cctfit$BUGSoutput$n.sims){
cctfit$BUGSoutput$sims.list[["tau"]][samp,,] <- cctfit$BUGSoutput$sims.list[["E"]][samp,]*t(1/cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]])

cctfit$ppY[,,samp] <- matrix(
rnorm((cctfit$n*cctfit$m),(cctfit$BUGSoutput$sims.list[["a"]][samp,]*t(cctfit$BUGSoutput$sims.list[["T"]][samp,,cctfit$BUGSoutput$sims.list[["e"]][samp,]]))+matrix(rep(cctfit$BUGSoutput$sims.list[["b"]][samp,],cctfit$m),cctfit$n,cctfit$m),
(cctfit$BUGSoutput$sims.list[["a"]][samp,]*(cctfit$BUGSoutput$sims.list[["tau"]][samp,,]^(-.5)))), 
cctfit$n,cctfit$m)

cctfit$ppX[,,samp] <- (array(rep(1/cctfit$BUGSoutput$sims.list[["a"]][samp,],times=cctfit$m),c(cctfit$n,cctfit$m))*(cctfit$ppY[,,samp]))-matrix(rep(cctfit$BUGSoutput$sims.list[["b"]][samp,],times=cctfit$m),c(cctfit$n,cctfit$m))
}
cctfit$mean$ppY <- rowMeans(cctfit$ppY[,,],dims=2)

if(cctfit$mval==1){
for(i in 1:dim(cctfit$datob$thena)[1]){cctfit$data[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2]] <- mean(cctfit$ppY[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2],])}
cctfit$MVest <- cbind(cctfit$datob$thena,cctfit$data[cctfit$datob$thena])
colnames(cctfit$MVest) <- c("Pers","Item","Resp")
}

leigv <- 12
options(warn=-3)
ind <- sample(cctfit$BUGSoutput$n.sims,min(500,cctfit$BUGSoutput$n.sims)) #500 is the number of samples
eigv <- matrix(-1,length(ind),leigv)
tmp <- apply(cctfit$ppY[,,ind],c(1,3),function(x) t(x))
tmp2 <- apply(tmp[,,],3,function(x) cor(x))
tmp3 <- array(tmp2,c(cctfit$n,cctfit$n,cctfit$BUGSoutput$n.sims))
for(i in 1:length(ind)){
suppressMessages(try(eigv[i,] <- fa(tmp3[,,i])$values[1:leigv],silent=TRUE))}
wch <- -which(eigv[,1] == -1); 
if(length(wch)==0){cctfit$ppeig <-  eigv}else{cctfit$ppeig <-  eigv[-which(eigv[,1] == -1),]}
cctfit$dateig <-  suppressMessages(fa(cor(t(cctfit$data)))$values[1:leigv])
options(warn=0)

if(cctfit$V==1){
vdi <- matrix(-1,dim(cctfit$ppY)[3],1)
varvec <- array(-1,c(length(ind),dim(cctfit$ppY)[2]))

for(i in 1:length(ind)){
varvec[i,] <- apply(invlogit(cctfit$ppY[,,i]),2,var)  
}
vdi <- apply(varvec,1,var)
cctfit$ppVDI <- apply(varvec,1,var) 

cctfit$datVDI <- var(apply(invlogit(cctfit$data),2,var))  

}
if(cctfit$V>1){
varvec <- array(NA,c(length(ind),dim(cctfit$ppY)[2],cctfit$V))
cctfit$datVDI  <- array(NA,c(1,cctfit$V))
options(warn=-3)
for(v in 1:cctfit$V){
for(i in 1:length(ind)){
if(sum(cctfit$BUGSoutput$sims.list$e[ind[i],] == v)>1){
suppressMessages(try(varvec[i,,v] <- apply(invlogit(cctfit$ppY[cctfit$BUGSoutput$sims.list$e[ind[i],] == v,,i]),2,var),silent=TRUE))
}else{
suppressMessages(try(varvec[i,,v] <- var(invlogit(cctfit$ppY[cctfit$BUGSoutput$sims.list$e[ind[i],] == v,,i])),silent=TRUE)) 
}
}
if(sum(cctfit$respmem == v)>1){
cctfit$datVDI[v] <- var(apply(invlogit(cctfit$data[cctfit$respmem == v,]),2,var)) 
}else{
cctfit$datVDI[v] <- var(invlogit(cctfit$data[cctfit$respmem == v,])) 
}
}
options(warn=0)
vdi <- apply(varvec,c(1,3),function(x) var(x,na.rm=TRUE))
cctfit$ppVDI <- vdi

if(sum(apply(cctfit$ppVDI,2,function(x) all(x == 0,na.rm=TRUE)))>=1){
for(i in which(apply(cctfit$ppVDI,2,function(x) all(x == 0)))){
cctfit$ppVDI[,i] <- NA
}
}
for(v in 1:dim(cctfit$ppVDI)[2]){
if(identical(unique(cctfit$ppVDI[,v]),c(NA,0)) || identical(unique(cctfit$ppVDI[,v]),c(0,NA))){cctfit$ppVDI[,v] <- NA}
}
}

rm(vdi, varvec)
rm(eigv,leigv,tmp,tmp2,tmp3,ind);

cctfit$checksrun <- 1
}

}

if(cctfit$mval==1 && gui==1){
message("\n ...Type 'cctfit$MVest' to display the posterior predictive estimates for missing data cells \n ...     'cctfit$data' to view the full data matrix \n ...     'cctfit$datamiss' to view the matrix with missing values")
}
message("\n ...Posterior predictive checks complete")
}

if(gui==1){cctfit <<- cctfit}

if(saveplots==1){jpeg(file.path(gsub(".Rdata","ppc.jpg",savedir)),width = 6, height = 3, units = "in", pointsize = 12,quality=100,res=400)}
if(saveplots==2){postscript(file=file.path(gsub(".Rdata","ppc.eps",savedir)), onefile=FALSE, horizontal=FALSE, width = 6, height = 6, paper="special", family="Times")}

par(oma=c(0,0,0,0),mar=c(4,4,3,1),mgp=c(2.25,.75,0),mfrow=c(1,2))

plot(NA, main="Culture Number Check",xlim=c(1,min(cctfit$n,8)),ylim=c(min(cctfit$dateig[min(cctfit$n,8)],cctfit$ppeig[,min(cctfit$n,8)]),ceiling(max(cctfit$ppeig[,1],cctfit$dateig[1]))),xlab="Eigenvalue",ylab="Value",las=1,pch=21,type="l",col="black")
for(i in 1:dim(cctfit$ppeig)[1]){
points(cctfit$ppeig[i,],col="grey",type="l") }
points(cctfit$dateig,col="black",type="l")

if(cctfit$V>1){
cctfit$VDIperc <- array(NA,c(cctfit$V,1))
color <- "black"; color2 <- "black"
color <- c("black","dark grey","light grey",rainbow(cctfit$V))
minmax <- array(NA,c(cctfit$V,2))
minmax2 <- array(NA,c(cctfit$V,1))

for(v in which(sort(which(!apply(cctfit$ppVDI,2,function(x) all(is.na(x))))) %in% sort(unique(cctfit$respmem)))){
tmpdist <-ecdf(cctfit$ppVDI[,v])
minmax[v,] <- c(quantile(tmpdist, .005),quantile(tmpdist, .995))
minmax2[v] <- max(density(cctfit$ppVDI[,v],na.rm=TRUE)$y)
}
plot(NA,main="Item Difficulty Check Resp", xlab= "VDI", xlim=c(min(cctfit$datVDI[!is.na(cctfit$datVDI)],minmax,na.rm=TRUE),max(cctfit$datVDI[!is.na(cctfit$datVDI)],minmax,na.rm=TRUE)), ylim=c(0,max(minmax2,na.rm=TRUE)), ylab="Value",las=1,col=color[1],type="l")
for(v in which(sort(which(!apply(cctfit$ppVDI,2,function(x) all(is.na(x))))) %in% sort(unique(cctfit$respmem)))){
tmpdist <-ecdf(cctfit$ppVDI[,v])
points(density(cctfit$ppVDI[,v],na.rm=TRUE),lwd=1.5,main="Item Difficulty Check Resp", xlab= "VDI", xlim=c(min(cctfit$datVDI[!is.na(cctfit$datVDI)],quantile(tmpdist, .005),na.rm=TRUE),max(cctfit$datVDI[!is.na(cctfit$datVDI)],quantile(tmpdist, .995),na.rm=TRUE)), ylim=c(0,max(density(cctfit$ppVDI,na.rm=TRUE)$y,na.rm=TRUE)), ylab="Value",las=1,col=color[v],type="l")
if(saveplots != 1 && saveplots != 2){print(paste("VDI Culture ",v," : ",round(100*tmpdist(cctfit$datVDI[,v]),digits=2)," percentile"),col=color2) 
cctfit$VDIperc[v,1] <- round(100*tmpdist(cctfit$datVDI[,v]),digits=2)
}
}
rm(tmpdist)
segments(cctfit$datVDI,0,cctfit$datVDI,par("usr")[4],col=color,lwd=1.5)
}else{
color <- "black"; color2 <- "black"
tmpdist <-ecdf(cctfit$ppVDI)
plot(density(cctfit$ppVDI,na.rm=TRUE),main="Item Difficulty Check Resp", xlab= "VDI", xlim=c(min(cctfit$datVDI,quantile(tmpdist, .005)),max(cctfit$datVDI,quantile(tmpdist, .995))), ylim=c(0,max(density(cctfit$ppVDI,na.rm=TRUE)$y,na.rm=TRUE)), ylab="Value",las=1,col=color2)
segments(cctfit$datVDI,0,cctfit$datVDI,par("usr")[4],col=color2)
text(.7*par("usr")[2],.7*par("usr")[4],labels=paste(round(100*tmpdist(cctfit$datVDI),digits=2)," percentile"),col=color2) 
if(saveplots != 1 && saveplots != 2){cctfit$VDIperc <- round(100*tmpdist(cctfit$datVDI),digits=2)
print(paste("VDI Culture 1 : ",round(100*tmpdist(cctfit$datVDI),digits=2)," percentile"),col=color2)
}; rm(tmpdist)
}

rm(color,color2)
if(saveplots==1 || saveplots==2){dev.off()}
saveplots <- 0

if(gui==1){cctfit <<- cctfit}
if(gui==0){return(cctfit)}

}

######################
#Function for the 'Export Results' Button
#- Prompts user where to save, and the filename to use
#- Saves the inference results (cctfit) to an .Rdata file of that filename
#- Saves .eps and .jpeg plots of the scree plot, plot results button, and run checks button
#    if the checks were calculated (the checks button was pressed)
######################
exportfuncbutton <- function() {
exportfunc(cctfit=cctfit,gui=1)
}

exportfunc <- function(cctfit,filename=paste(getwd(),"/CCTpackdata.Rdata",sep=""),gui=0) {

if(gui==1){
savedir <- tclvalue(tkgetSaveFile(initialfile = "CCTpackdata.Rdata",filetypes = "{{Rdata Files} {.Rdata}}"))
if (!nchar(savedir)) {
message("\n ...Export cancelled\n")
return()}
}else{
savedir <- filename
}

if (savedir == 0 && !file.exists("CCTpack")){dir.create(file.path(getwd(), "CCTpack"));
savedir <- file.path(getwd(),"CCTpack","CCTpackdata.Rdata")
} 

if(substr(savedir, nchar(savedir)-6+1, nchar(savedir)) != ".Rdata"){
savedir <- paste(savedir,".Rdata",sep="")
}

message("\n ...Exporting results")
write.csv(cctfit$BUGSoutput$summary,file.path(gsub(".Rdata","posterior.csv",savedir)))
save(cctfit,file=file.path(savedir))

if(cctfit$checksrun == 1){
ppcfunc(cctfit=cctfit,saveplots=1,savedir); ppcfunc(cctfit=cctfit,saveplots=2,savedir)}
plotresultsfunc(cctfit=cctfit,saveplots=1,savedir); plotresultsfunc(cctfit=cctfit,saveplots=2,savedir);
screeplotfunc(datob = cctfit$datob,saveplots=1,savedir); screeplotfunc(datob = cctfit$datob,saveplots=2,savedir)

message("\n ...Export complete\n")
}


######################
#Manual function to do all of these things if the gui is not compatible for the user, or not preferred
#- Takes in a data as a respondent by item array or matrix
#- 'clusters' defines # of clusters, and 'itemdiff' if item difficulty is wanted
#- 'jags' defines for the JAGS MCMC, the number of samples, chains, burn-in, thinning
#- 'runchecks' if one wants the posterior predictive checks calculated after inference
#- 'exportfilename' set a name different from "" if one wants to automatically export the results 
#    to the working directory
######################
cctapply <<- function(data,clusters=1,itemdiff=0,jags.iter=10000,jags.chains=3,jags.burnin=2000,jags.thin=1,runchecks=1,exportfilename="",polych=0){
datob <- loadfilefunc(data)
datob <- screeplotfunc(datob)
cctfit <- applymodelfunc(datob,clusters=clusters,itemdiff=itemdiff,jags.iter=jags.iter,jags.chains=jags.chains,jags.burnin=jags.burnin,jags.thin=jags.thin)
plotresultsfunc(cctfit)
if(runchecks==1){
cctfit <- ppcfunc(cctfit,polych=polych)
}
if(exportfilename!=""){
exportfunc(cctfit,filename=exportfilename)
}
return(cctfit)
}

######################
#Manual function to do the scree plot
#- Takes in a data as a respondent by item array or matrix
######################
cctscree <<- function(data,polych=0){
datob <- loadfilefunc(data)
datob <- screeplotfunc(datob,polych)
}

######################
#Manual function to plot the results, equivalent to the 'Plot Results Button'
#- Takes the cctfit object from cctapply() or the 'Apply CCT Model' button
######################
cctresults <<- function(cctfit){
plotresultsfunc(cctfit)
}

######################
#Manual function to plot the posterior predictive check plots, equivalent to the 'Run Checks Button'
#- Takes the cctfit object from cctapply() or the 'Apply CCT Model' button
#- Plots the posterior predictive checks, or calculates them and then plots them (if they have not been calculated yet)
######################
cctppc <<- function(cctfit,polych=0){
cctfit <- ppcfunc(cctfit,polych=polych)
return(cctfit)
}

######################
#Manual function to export the results, equivalent to the 'Export Results Button'
#- Takes the cctfit object from cctapply() or the 'Apply CCT Model' button
#- Exports the cctfit object as a .Rdata file, and plots of the scree plot, 
#    posterior results, and posterior predictive checks
######################
cctexport <<- function(cctfit,filename="CCTpackdata.Rdata"){
cctfit <- exportfunc(cctfit,filename=filename)
}

######################
#The GUI Grid Setup
######################
loaddata.but <- tkbutton(datframe, text="Load Data", command=loadfilefuncbutton)
screeplot.but <- tkbutton(datframe, text="Scree Plot", command=screeplotfuncbutton)
applymodel.but <- tkbutton(applyframe, text="Apply CCT Model", command=applymodelfuncbutton)
plotresults.but <- tkbutton(resultsframe, text = "Plot Results", command = plotresultsfuncbutton)
doppc.but <- tkbutton(resultsframe, text = "Run Checks", command = ppcfuncbutton)
exportresults.but <- tkbutton(resultsframe, text = "Export Results", command = exportfuncbutton)

sameid <- tkradiobutton(applyframe)
diffid <- tkradiobutton(applyframe)
itemdiffvar <- tclVar("0")
tkconfigure(sameid,variable=itemdiffvar, value="0")
tkconfigure(diffid,variable=itemdiffvar, value="1")
polyyes <- tkradiobutton(datframe4)
polyno <- tkradiobutton(datframe4)
polyvar <- tclVar("0")
tkconfigure(polyyes,variable=polyvar, value="1")
tkconfigure(polyno,variable=polyvar, value="0")

datafiletxt <- tktext(tt,bg="white",width=20,height=1)
resptxt <- tktext(tt,bg="white",width=4,height=1)
itemtxt <- tktext(tt,bg="white",width=4,height=1) 
dattypetxt <- tktext(tt,bg="white",width=11,height=1) 
modeltxt <- tktext(tt,bg="white",width=4,height=1)

tkgrid(datframe,columnspan=3,row=1,column=0)
tkgrid(datframe2,columnspan=4,row=2,column=0)
tkgrid(datframe3,columnspan=4,row=3,column=0)
tkgrid(datframe4,columnspan=5,row=4,column=0)
tkgrid(applyframe,columnspan=10,row=5,column=0)
tkgrid(resultsframe,columnspan=3,row=6)
tkgrid(settingsframe,columnspan=8,row=7)

tkgrid(tklabel(datframe,text="    "))
tkgrid(tklabel(datframe,text="Data Input"),columnspan=3, pady = 5) 

tkgrid(loaddata.but,datafiletxt,screeplot.but,pady= 10, padx= 10)
tkgrid(tklabel(datframe2,text="Number of Respondents"),resptxt,tklabel(datframe2,text="Number of Items"),itemtxt, padx = 2, pady = 5) 
tkgrid(tklabel(datframe3,text="Data Type Detected"),dattypetxt,tklabel(datframe3,text="CCT Model"),modeltxt, padx = 2, pady = 5) 
tkgrid(tklabel(applyframe,text="Model Application"),columnspan=8, pady = 5) 
tkgrid(tklabel(applyframe,text="Number of Cultures to Assume:"), cultures.entry, tklabel(applyframe,text="Item Difficulty:"),tklabel(applyframe,text="No"),sameid,tklabel(applyframe,text="Yes"),diffid,applymodel.but,pady= 10, padx= 2)

tkconfigure(screeplot.but, state="disabled") 
tkgrid(tklabel(resultsframe,text="Application Results"),columnspan=3, pady = 5) 
tkgrid(doppc.but,plotresults.but,exportresults.but,pady= 10, padx= 10)
tkconfigure(applymodel.but, state="disabled")
tkconfigure(plotresults.but, state="disabled")
tkconfigure(doppc.but, state="disabled")
tkconfigure(exportresults.but, state="disabled")
tkgrid(tklabel(settingsframe,text="Sampler Settings (Optional)"),columnspan=8, pady = 5) 
tkgrid(tklabel(settingsframe,text="Samples"), samples.entry, tklabel(settingsframe,text="Chains"), chains.entry, tklabel(settingsframe,text="Burn-in"), burnin.entry, tklabel(settingsframe,text="Thinning"), thin.entry, pady= 10, padx= 2)

tkinsert(datafiletxt,"end","(load data file)")
tkconfigure(datafiletxt, state="disabled")
tkconfigure(resptxt, state="disabled")
tkconfigure(itemtxt, state="disabled")
tkconfigure(dattypetxt, state="disabled")
tkconfigure(modeltxt, state="disabled")

tkgrid(tklabel(tt,text="    "))
tkgrid.columnconfigure(tt,0,weight=1) 
tkgrid.rowconfigure(tt,0,weight=1) 
tkwm.resizable(tt,0,0)

######################
#Model Code
#- Used by the 'Apply CCT Model' button/function
#- Each model has 2 code versions, 1 without item difficulty, 1 with item difficulty
######################
mcgcm <-
"model{
for (l in 1:nobs){
  D[Y[l,1],Y[l,2]] <- (th[Y[l,1]]*(1-lam[Y[l,2]])) / ((th[Y[l,1]]*(1-lam[Y[l,2]]))+(lam[Y[l,2]]*(1-th[Y[l,1]]))) 
  pY[Y[l,1],Y[l,2]] <- (D[Y[l,1],Y[l,2]]*Z[Y[l,2],e[Y[l,1]]]) +((1-D[Y[l,1],Y[l,2]])*g[Y[l,1]])
  Y[l,3] ~ dbern(pY[Y[l,1],Y[l,2]]) }
  
for (i in 1:n){
 e[i] ~ dcat(pi) 
 th[i] ~ dbeta(thmu[e[i]]*thtau[e[i]],(1-thmu[e[i]])*thtau[e[i]])
 g[i] ~ dbeta(gmu[e[i]]*gtau[e[i]],(1-gmu[e[i]])*gtau[e[i]]) }

 for (k in 1:m){
  lam[k] <- .5
 for (v in 1:V){
  Z[k,v] ~ dbern(p[v]) }}

#Hyper Parameters
 gsmu <- 10
 gssig <- 10
 dsmu <- 10
 dssig <- 10
 pi[1:V] ~ ddirch(L)
 alpha <- 2
 
for (v in 1:V){
 p[v] ~ dunif(0,1)
 gmu[v] <- .5
 gtau[v] ~ dgamma(pow(gsmu,2)/pow(gssig,2),gsmu/pow(gssig,2))
 thmu[v] ~ dbeta(alpha,alpha)
 thtau[v] ~ dgamma(pow(dsmu,2)/pow(dssig,2),dsmu/pow(dssig,2))
 L[v] <- 1 }}"

mcgcmid <-
"model{
for (l in 1:nobs){
  D[Y[l,1],Y[l,2]] <- (th[Y[l,1]]*(1-lam[Y[l,2],e[Y[l,1]]])) / ((th[Y[l,1]]*(1-lam[Y[l,2],e[Y[l,1]]]))+(lam[Y[l,2],e[Y[l,1]]]*(1-th[Y[l,1]])))   
  pY[Y[l,1],Y[l,2]] <- (D[Y[l,1],Y[l,2]]*Z[Y[l,2],e[Y[l,1]]]) +((1-D[Y[l,1],Y[l,2]])*g[Y[l,1]])
  Y[l,3] ~ dbern(pY[Y[l,1],Y[l,2]]) }  

for (i in 1:n){
 e[i] ~ dcat(pi) 
 th[i] ~ dbeta(thmu[e[i]]*thtau[e[i]],(1-thmu[e[i]])*thtau[e[i]])
 g[i] ~ dbeta(gmu[e[i]]*gtau[e[i]],(1-gmu[e[i]])*gtau[e[i]]) }

 for (k in 1:m){
 for (v in 1:V){
  Z[k,v] ~ dbern(p[v])
  lam[k,v] ~ dbeta(lammu[v]*lamtau[v],(1-lammu[v])*lamtau[v])
  }}

#Hyper Parameters
 lamsmu <- 10
 lamssig <- 10
 gsmu <- 10
 gssig <- 10
 dsmu <- 10
 dssig <- 10
 alpha <- 2
 pi[1:V] ~ ddirch(L)
 
for (v in 1:V){
 p[v] ~ dunif(0,1)
 lammu[v] <- .5
 lamtau[v] ~ dgamma(pow(lamsmu,2)/pow(lamssig,2),lamsmu/pow(lamssig,2))
 gmu[v] <- .5
 gtau[v] ~ dgamma(pow(gsmu,2)/pow(gssig,2),gsmu/pow(gssig,2))
 thmu[v] ~ dbeta(alpha,alpha)
 thtau[v] ~ dgamma(pow(dsmu,2)/pow(dssig,2),dsmu/pow(dssig,2))
 L[v] <- 1
 }}"

mcltrm <-
"model{
for (l in 1:nobs){
  tau[Y[l,1],Y[l,2]] <- E[Y[l,1]]
 	pY[Y[l,1],Y[l,2],1] <- pnorm((a[Y[l,1]]*gam[1,e[Y[l,1]]]) + b[Y[l,1]],T[Y[l,2],e[Y[l,1]]],tau[Y[l,1],Y[l,2]])
	for (c in 2:(C-1)){pY[Y[l,1],Y[l,2],c] <- pnorm((a[Y[l,1]]*gam[c,e[Y[l,1]]]) + b[Y[l,1]],T[Y[l,2],e[Y[l,1]]],tau[Y[l,1],Y[l,2]]) - sum(pY[Y[l,1],Y[l,2],1:(c-1)])}
	pY[Y[l,1],Y[l,2],C] <- (1 - sum(pY[Y[l,1],Y[l,2],1:(C-1)]))
	Y[l,3] ~ dcat(pY[Y[l,1],Y[l,2],1:C])
  }

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
  Tmu[v] ~ dnorm(0,.01)
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
for (l in 1:nobs){
  tau[Y[l,1],Y[l,2]] <- E[Y[l,1]]/lam[Y[l,2],e[Y[l,1]]]
 	pY[Y[l,1],Y[l,2],1] <- pnorm((a[Y[l,1]]*gam[1,e[Y[l,1]]]) + b[Y[l,1]],T[Y[l,2],e[Y[l,1]]],tau[Y[l,1],Y[l,2]])
	for (c in 2:(C-1)){pY[Y[l,1],Y[l,2],c] <- pnorm((a[Y[l,1]]*gam[c,e[Y[l,1]]]) + b[Y[l,1]],T[Y[l,2],e[Y[l,1]]],tau[Y[l,1],Y[l,2]]) - sum(pY[Y[l,1],Y[l,2],1:(c-1)])}
	pY[Y[l,1],Y[l,2],C] <- (1 - sum(pY[Y[l,1],Y[l,2],1:(C-1)]))
	Y[l,3] ~ dcat(pY[Y[l,1],Y[l,2],1:C])
  }

#Parameters
   for (i in 1:n){
      e[i] ~ dcat(pi) 
      E[i] ~ dgamma(pow(Emu[e[i]],2)*Etau[e[i]],Emu[e[i]]*Etau[e[i]])
      a[i] ~ dgamma(atau[e[i]],atau[e[i]])
      b[i] ~ dnorm(bmu[e[i]],btau[e[i]]) }

   for (k in 1:m){       
    for (v in 1:V){
     T[k,v] ~ dnorm(Tmu[v],Ttau[v])
     lam[k,v] ~ dgamma(lamtau[v],lamtau[v])	 }}

   for (v in 1:V){
     gam[1:(C-1),v] <- sort(tgam2[1:(C-1),v])
     for (c in 1:(C-2)){tgam[c,v] ~ dnorm(0,.1)}
     tgam2[1:(C-2),v] <- tgam[1:(C-2),v]
     tgam2[C-1,v] <- -sum(tgam[1:(C-2),v]) }

     pi[1:V] ~ ddirch(L)

#Hyperparameters	
 for (v in 1:V){
  L[v] <- 1 
  Tmu[v] ~ dnorm(0,.01)
  Ttau[v] ~ dgamma(1,.1)
  lammu[v] <- 1
  lamtau[v] ~ dgamma(1,0.1) 
  Emu[v] ~ dgamma(1,1)
  Etau[v] ~ dgamma(1,1)
  amu[v] <- 1
  atau[v] ~ dgamma(4,4)
  bmu[v] <- 0
  btau[v] ~ dgamma(4,4)
 }
}"

mccrm <-
"model{
   for (l in 1:nobs){ 
	Y[l,3] ~ dnorm((a[Y[l,1]]*T[Y[l,2],e[Y[l,1]]])+b[Y[l,1]],pow(a[Y[l,1]],-2)*E[Y[l,1]]) 
	}
	
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
  Tmu[v] ~ dnorm(0,1)
  Ttau[v] ~ dgamma(1,1)
  Emu[v] ~ dgamma(1,1)
  Etau[v] ~ dgamma(0.01,0.01)
  amu[v] <- 1
  atau[v] ~ dgamma(0.01,0.01)
  bmu[v] <- 0
  btau[v] ~ dgamma(0.01,0.01)
}}"

mccrmidlog <-
"model{
   for (l in 1:nobs){ 
	Y[l,3] ~ dnorm((a[Y[l,1]]*T[Y[l,2],e[Y[l,1]]])+b[Y[l,1]],pow(a[Y[l,1]],-2)*E[Y[l,1]]/lam[Y[l,2],e[Y[l,1]]])
	}	  

#Parameters
   for (i in 1:n){
      e[i] ~ dcat(pi) 
      Elog[i] ~ dnorm(Emu[e[i]],Etau[e[i]])
	  E[i] <- exp(Elog[i])
	
      alog[i] ~ dnorm(amu[e[i]],atau[e[i]])
	  a[i] <- exp(alog[i])
	 
      b[i] ~ dnorm(bmu[e[i]],btau[e[i]]) }

   for (k in 1:m){ 
    for (v in 1:V){
     T[k,v] ~ dnorm(Tmu[v],Ttau[v]) 
	 lamlog[k,v] ~ dnorm(lammu[v],lamtau[v]) 
	 lam[k,v] <- exp(lamlog[k,v])
	 }}

     pi[1:V] ~ ddirch(L)

#Hyperparameters	
 for (v in 1:V){
  L[v] <- 1 
  Tmu[v] ~ dnorm(0,1)
  Ttau[v] ~ dgamma(1,1)
  Emu[v] ~ dnorm(0,.01)
  Etau[v] ~ dgamma(.01,.01)
  lammu[v] <- 0
  lamtau[v] ~ dgamma(.01,.01)
  amu[v] <- 0
  atau[v] ~ dgamma(.01,.01)
  bmu[v] <- 0
  btau[v] ~ dgamma(.01,.01)
 }
}" 
	
mccrmid <-
"model{
   for (l in 1:nobs){ 
	Y[l,3] ~ dnorm((a[Y[l,1]]*T[Y[l,2],e[Y[l,1]]])+b[Y[l,1]],pow(a[Y[l,1]],-2)*E[Y[l,1]]/lam[Y[l,2],e[Y[l,1]]])
	}	  

#Parameters
   for (i in 1:n){
      e[i] ~ dcat(pi) 
      E[i] ~ dgamma(pow(Emu[e[i]],2)*Etau[e[i]],Emu[e[i]]*Etau[e[i]])
      a[i] ~ dgamma(atau[e[i]],atau[e[i]])
      b[i] ~ dnorm(bmu[e[i]],btau[e[i]]) }

   for (k in 1:m){ 
    for (v in 1:V){
     T[k,v] ~ dnorm(Tmu[v],Ttau[v]) 
	 lam[k,v] ~ dgamma(lamtau[v],lamtau[v]) 
	 }}

     pi[1:V] ~ ddirch(L)

#Hyperparameters	
 for (v in 1:V){
  L[v] <- 1 
  Tmu[v] ~ dnorm(0,1)
  Ttau[v] ~ dgamma(1,1)
  lammu[v] <- 1
  lamtau[v] ~ dgamma(0.01,0.01)
  Emu[v] ~ dgamma(1,1)
  Etau[v] ~ dgamma(0.01,0.01)
  amu[v] <- 1
  atau[v] ~ dgamma(0.01,0.01)
  bmu[v] <- 0
  btau[v] ~ dgamma(0.01,0.01)
  }
}"
class(cctscree) <<- "occult"; class(cctapply) <<- "occult"; class(cctresults) <<- "occult"; class(cctppc) <<- "occult"; class(cctexport) <<- "occult"
}
class(cctgui) <- "occult"; 
print.occult <- function(x, ...){
	print(args(x))
	invisible(x)
}
######################
#Occult Class hides the very long output that would be generated from typing in the cctgui function accidentally
#The arguments are still displayed for the functions however, for the user to see what can be input
######################