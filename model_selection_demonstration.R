##################################################################
####                                                          ####
####             DEMONSTRATING MODEL SELECTION                ####
####                                                          ####
##################################################################

library(lme4)


## 1. LOAD IN THE NECESSARY DATA AND HAVE A LOOK
load("reduced_colonisation_dataset.RData")
head(col_reduced)


## 2. SET UP ALL THE COMBINATIONS OF VARIABLES
trait_variables  <- c("lprev", "log10(ABW)", "log10(ALD)", "factor(Habitat)")

a <- combn(trait_variables,1)  #4
b <- combn(trait_variables,2)  #6
c <- combn(trait_variables,3)  #4
d <- combn(trait_variables,4)  #1

a
b
c
d


## 3. RUN THE MODEL SELECTION FOR ALL THE COMBINATIONS- FOR LOOPS

## 3a. 1 Variable
nvars <- matrix(0,4,1)
LL <- matrix(0,4,1)
mAIC <- matrix(0,4,1)
mBIC <- matrix(0,4,1)
mVARS <- vector(mode = "list", length = 4)

for (x in 1:4) {
  
  fmla <- paste("transition ~ lambda_pt15 +", (paste(a[,x],collapse="+")),
                " + (1|species_code) + (1|brc.name)")
  
  
  lme_mod <-glmer(fmla, family="binomial", data=col_reduced, 
                  control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 50000)))
  
  
  nvars <- 1
  LL[x] <- logLik(lme_mod)
  mAIC[x] <- AIC(lme_mod)
  mBIC[x] <- BIC(lme_mod)
  mVARS <- a[,x]
  rm(lme_mod, fmla)
}

col_1variable <- data.frame(nvars=nvars, LL=LL, mAIC=mAIC,  mBIC= mBIC, model= t(a))
rm(nvars,LL,mAIC,mVARS,mBIC)




## 3b. 2 Variables
nvars <- matrix(0,6,1)
LL <- matrix(0,6,1)
mAIC <- matrix(0,6,1)
mBIC <- matrix(0,6,1)
mVARS <- vector(mode = "list", length = 6)

for (x in 1:6) {
  
  fmla <- paste("transition ~ lambda_pt15 +", (paste(b[1:2,x],collapse="+")),
                " + (1|species_code) + (1|brc.name)")
  
  
  lme_mod <-glmer(fmla, family="binomial", data=col_reduced, 
                  control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 50000)))
  
  
  nvars <- 2
  LL[x] <- logLik(lme_mod)
  mAIC[x] <- AIC(lme_mod)
  mBIC[x] <- BIC(lme_mod)
  mVARS <- b[1:2,x]
  rm(lme_mod, fmla)
}

col_2variables <- data.frame(nvars=nvars, LL=LL, mAIC=mAIC,  mBIC= mBIC, model= t(b))
rm(nvars,LL,mAIC,mVARS,mBIC)





## 3c. 3 Variables
nvars <- matrix(0,4,1)
LL <- matrix(0,4,1)
mAIC <- matrix(0,4,1)
mBIC <- matrix(0,4,1)
mVARS <- vector(mode = "list", length = 4)

for (x in 1:4) {
  
  fmla <- paste("transition ~ lambda_pt15 +", (paste(c[1:3,x],collapse="+")),
                " + (1|species_code) + (1|brc.name)")
  
  
  lme_mod <-glmer(fmla, family="binomial", data=col_reduced, 
                  control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 50000)))
  
  
  nvars <- 1
  LL[x] <- logLik(lme_mod)
  mAIC[x] <- AIC(lme_mod)
  mBIC[x] <- BIC(lme_mod)
  mVARS <- c[1:3,x]
  rm(lme_mod, fmla)
}

col_3variables <- data.frame(nvars=nvars, LL=LL, mAIC=mAIC,  mBIC= mBIC, model= t(c))
rm(nvars,LL,mAIC,mVARS,mBIC)



## 3d. 4 Variables

fmla <- paste("transition ~lambda_pt15 +", (paste(d,collapse="+")),
              " + (1|species_code) + (1|brc.name)")


lme_mod <-glmer(fmla, family="binomial", data=col_reduced, 
                control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 50000)))

nvars<- 4
LL <- logLik(lme_mod)
mAIC <- AIC(lme_mod)
mBIC <- BIC(lme_mod)
mVARS <- d

col_4variables <- data.frame(nvars=nvars, LL=LL, mAIC=mAIC,  mBIC= mBIC, model= t(d))
rm(LL,mAIC,mVARS, mBIC)


## 3e. BINDING THE RESULTS TOGETHER

ncols7t <- ncol(col_4variables)#final no. collums is that of the 7 variable model
nms7t <- names(col_4variables) # final col names are same as those for the 7 variable model
col_trait_models_all<-do.call(rbind, lapply(list(col_1variable,col_2variables,col_3variables,
                                                 col_4variables), function(x){
                    if((coldifft<-ncols7t-ncol(x))!=0){
                    adding_colst<-cbind(x,
                        data.frame(matrix(rep(NA, times=coldifft*nrow(x)),ncol=coldifft)))
                        } else {adding_colst <- x}
                        names(adding_colst) <- nms7t
                        return(adding_colst)
                                  }))
#function loop looks at all the dataframes in the do.call part that is rbinding all the necessary things. it asks what are the difference in the number of columns in comparison to the 7 variable one(max number of collums) and then creates extra collums (see adding_cols)


col_trait_models_mBIC<- col_trait_models_all[order(col_trait_models_all$mBIC),]




## 4. USING AN LAPPLY

traitlist <- list(a,b,c,d)

mymodel<- function(comb, x, data){
  
  fmla <- paste("transition ~ lambda_pt15 +", 
                (paste(comb[,x],collapse="+")),
                " + (1|species_code) +(1|brc.name)")
  
  lme_mod <-glmer(fmla, family="binomial", data=data,
                  control=glmerControl(optimizer="bobyqa",
                          optCtrl = list(maxfun = 50000)))
  
  return(data.frame(nvars = ncol(comb),
                    LL = logLik(lme_mod),
                    mAIC = AIC(lme_mod),
                    mBIC = BIC(lme_mod),
                    mVARS = paste(comb[,x], collapse="+")))
  
  
}




COL_TRAIT<- do.call(rbind, lapply(traitlist, function(COMB){
    do.call(rbind, lapply(1:ncol(COMB),function(X) {
      return(mymodel(comb=COMB,x=X,data=col_reduced))
    }))
  }))





















