setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(11)

# package
list.of.packages <- c("easypackages", "ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

library("easypackages")
libraries(list.of.packages)

source("utils.r")

num_to_return <- 1
upper_num_limit <- 1000
exe_num <- sample(1:upper_num_limit, num_to_return)

alpha <- 1

# Matthews correlation coefficient
mcc <- function(thisTP, thisTN, thisFP, thisFN)
{
  # Compute the Matthews correlation coefficient (MCC) score
  # Jeff Hebert 9/1/2016
  # Geoffrey Anderson 10/14/2016 
  # Added zero denominator handling.
  # Avoided overflow error on large-ish products in denominator.
  #
  # actual = vector of true outcomes, 1 = Positive, 0 = Negative
  # predicted = vector of predicted outcomes, 1 = Positive, 0 = Negative
  # function returns MCC

  
  #TP;TN;FP;FN # for debugging
  sum1 <- thisTP+thisFP; 
  sum2 <-thisTP+thisFN ; 
  sum3 <-thisTN+thisFP ; 
  sum4 <- thisTN+thisFN;
  denom <- as.double(sum1)*sum2*sum3*sum4 # as.double to avoid overflow error on large products
  if (any(sum1==0, sum2==0, sum3==0, sum4==0)) {
    denom <- 1
  }
  mcc <- ((thisTP*thisTN)-(thisFP*thisFN)) / sqrt(denom)
  
  # cat("\nMCC = ", (mcc), "\n\n", sep="")
  
  return(mcc)
}

# Normalized Matthews correlation coefficient
normMCC <- function(thisTP, thisTN, thisFP, thisFN)
{
        thisMCC <- mcc(thisTP, thisTN, thisFP, thisFN)
        thisNormMCC <- (thisMCC + 1) /2
        # cat("normMCC = ", dec_three(thisNormMCC), "\n", sep="")
        # cat("MCC = ", dec_three(thisMCC), "\n\n", sep="")
        
        return(thisNormMCC)
}

# prevalence_threshold
prevalence_threshold <-  function(thisTP, thisTN, thisFP, thisFN)
{

        this_PT <- NULL
        thisTPR <- thisTP / (thisTP + thisFN)
        thisTNR <- thisTN / (thisTN + thisFP)

        this_PT_numerator <- (sqrt(thisTPR * (-thisTNR+1))) + thisTNR - 1
        this_PT_demonimator <- thisTNR + thisTPR -1
        
        if(this_PT_demonimator != 0)  this_PT <-  this_PT_numerator / this_PT_demonimator

        return(this_PT)
}


# complPT
complPT <-  function(thisTP, thisTN, thisFP, thisFN)
{
        this_PT <- prevalence_threshold(thisTP, thisTN, thisFP, thisFN)
        
        if(is.null(this_PT)) return(NULL)
        else {
                this_complPT <- 1 - this_PT
                return(this_complPT)
        }
}


# compute DOR
normDOR <- function(thisTP, thisTN, thisFP, thisFN) { 

#         cat("\nTP = ", thisTP, "\n", sep="")
#         cat("TN = ", thisTN, "\n", sep="")
#         cat("FP = ", thisFP, "\n", sep="")
#         cat("FN = ", thisFN, "\n", sep="")
        
        ground_truth_positives <- TP + FN
        ground_truth_negatives <- TN + FP
        # cat("ground_truth_positives = ", ground_truth_positives, "\n", sep="")
        # cat("ground_truth_negatives = ", ground_truth_negatives, "\n", sep="")
                
        thisM <- thisTP + thisTN + thisFP + thisFN
         #cat("M = ", thisM, "\n", sep="")
        
#         # max TP and TP for M, with FP and FN = 1
#         quota <- (thisM - 2)/2
#         maxTP <- quota
#         maxTN <- quota
#         
#         cat("maxTP = ", (maxTP), "\n", sep="")
#         cat("maxTN = ", (maxTN), "\n", sep="")

        thisDOR <- (thisTP * thisTN) / (thisFP * thisFN)        
        
#          cat("DOR = ", dec_three(thisDOR), "\n", sep="")
        
         minDORhere <- alpha / ((ground_truth_positives-1) * (ground_truth_negatives-1))
        # cat("minDORhere = ", dec_three(minDORhere), "\n", sep="")
        
        maxDORhere <- (ground_truth_positives-1) * (ground_truth_negatives-1) / alpha
        # cat("maxDORhere = ", maxDORhere, "\n", sep="")
        
        # thisNormDOR <- ( thisDOR  / (((thisM/2)-1)^2) ) 
#         # thisNormDOR <- (thisDOR + maxDORhere ) / ((maxDORhere - minDORhere)*2)
#         thisNormDOR <- (thisDOR + maxDORhere) / (maxDORhere * maxDORhere)
#         thisNormDOR <- (thisDOR - minDORhere ) / (maxDORhere - minDORhere)
                
        thisNormDOR = thisDOR / ( thisDOR + alpha)
                
#          cat("normDOR = ", dec_three(thisNormDOR), "\n", sep="")
        
        return(thisNormDOR)
}

# # TN = 1; TP = 1; FP = 10; FN = 10;
# # normDOR(TP, TN, FP, FN)
# 
# TN = 10; TP = 10; FP = 1; FN = 1;
# normDORhere <- normDOR(TP, TN, FP, FN)
# normMCChere <- normMCC(TP, TN, FP, FN)
# 
# TN = 0; TP = 10; FP = 5; FN = 5;
# normDORhere <- normDOR(TP, TN, FP, FN)
# normMCChere <- normMCC(TP, TN, FP, FN)
# 
# TN = 5; TP = 10; FP = 5; FN = 0;
# normDORhere <- normDOR(TP, TN, FP, FN)
# normMCChere <- normMCC(TP, TN, FP, FN)
# 
# TP=2; TN=4; FP=8; FN=12  
# normDORhere <- normDOR(TP, TN, FP, FN)
# normMCChere <- normMCC(TP, TN, FP, FN)
# 
# TP=2; TN=2; FP=2; FN=24
# normDORhere <- normDOR(TP, TN, FP, FN)
# normMCChere <- normMCC(TP, TN, FP, FN)
# 
# TN = 999; TP = 999; FP =1 ; FN = 2;
# normDORhere <- normDOR(TP, TN, FP, FN)
# normMCChere <- normMCC(TP, TN, FP, FN)
# 
# TN = 99999; TP = 99999; FP = 2; FN = 2;
# normDORhere <- normDOR(TP, TN, FP, FN)
# normMCChere <- normMCC(TP, TN, FP, FN)
# 
# TN = 99999; TP = 99999; FP = 1; FN = 1;
# normDORhere <- normDOR(TP, TN, FP, FN)
# normMCChere <- normMCC(TP, TN, FP, FN)
# 
# TN = 40; TP = 10; FP = 5; FN = 45;
# normDORhere <- normDOR(TP, TN, FP, FN)
# normMCChere <- normMCC(TP, TN, FP, FN)
# 
# TN = 40; TP = 40; FP = 40; FN = 40;
# normDORhere <- normDOR(TP, TN, FP, FN)
# normMCChere <- normMCC(TP, TN, FP, FN)
# 
# 
# TP = 1000000; FP = 1000; FN = 1; TN = 1; 
# normDORhere <- normDOR(TP, TN, FP, FN)
# normMCChere <- normMCC(TP, TN, FP, FN)
# 
# TP = 100000; FP = 1000000; FN = 10; TN = 1; 
# normDORhere <- normDOR(TP, TN, FP, FN)
# normMCChere <- normMCC(TP, TN, FP, FN)
# 
# TP = 0; FP = 1; FN = 2; TN = 1000000; 
# normDORhere <- normDOR(TP, TN, FP, FN)
# normMCChere <- normMCC(TP, TN, FP, FN)
# 
# TP = 0; FP = 2; FN = 2; TN = 1000; 
# normDORhere <- normDOR(TP, TN, FP, FN)
# normMCChere <- normMCC(TP, TN, FP, FN)
# 
# 
# TP = 1; FP = 1000; FN = 555; TN = 10; 
# normDORhere <- normDOR(TP, TN, FP, FN)
# normMCChere <- normMCC(TP, TN, FP, FN)
# 
# TP = 1; FP = 1000; FN = 1; TN = 10000; 
# normDORhere <- normDOR(TP, TN, FP, FN)
# normMCChere <- normMCC(TP, TN, FP, FN)

STEP <- 1000
START <- 1
END <- 100000

this_threshold <- 0.4

cat("TP TN FP FN TPR TNR PPV NPV PT MCC complPT normMCC diff\n")

for(TP in seq(from = START, to = END, by = STEP)) {
    for(TN in seq(from = START, to = END, by = STEP)) {
        for(FP in seq(from = START, to = END, by = STEP)) {
            for(FN in seq(from = START, to = END, by = STEP)) {
                
                PThere <- prevalence_threshold(TP, TN, FP, FN)
                MCChere <- mcc(TP, TN, FP, FN)
                
                complPThere <- complPT(TP, TN, FP, FN)
                normMCChere <- normMCC(TP, TN, FP, FN)
                
                if(is.null(complPThere) == FALSE) {
                
                    diff <- abs(normMCChere - complPThere)
                
                    if(diff >= this_threshold) {
                        cat(TP, " ",sep="")
                        cat(TN, " ", sep="")
                        cat(FP, " ", sep="")
                        cat(FN, " ", sep="")
                        
                        innerTPR <- TP / (TP + FN)
                        innerTNR <- TN / (TN + FP)
                        innerPPV <- TP / (TP + FP)
                        innerNPV <- TN / (TN + FN)
                        
                        cat(dec_three(innerTPR)," ", sep="")
                        cat(dec_three(innerTNR)," ", sep="")
                        cat(dec_three(innerPPV)," ", sep="")
                        cat(dec_three(innerNPV)," ", sep="")
                        
                        cat(dec_three(PThere)," ", sep="")
                        cat(dec_three(MCChere)," ", sep="")
                    
                        cat(dec_three(complPThere)," ", sep="")
                        cat(dec_three(normMCChere)," ", sep="")
                    
                    
                        cat(diff, "\n", sep="")
                    }
                }
            }
        }
    }
}

computeExecutionTime()
