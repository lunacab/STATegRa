
####### DATA PREPROCESSING #######

# weightBlocks ------------------------------------------------------------
## Function to weigth the different blocks
## INPUT:
##    xx: (list) List of omics data sets (matrices)
## OUTPUT:
##    (list) List of weigthed omics data sets
setGeneric(
    name="weightBlocks",
    def=function(xx){standardGeneric("weightBlocks")}
)

setMethod(
    f="weightBlocks",
    signature=signature(xx="list"),
    definition=function(xx){
        Data <- do.call("rbind",xx)
        #dataWeight <- lapply(xx,function(X){X*(sqrt(as.numeric(nrow(Data))*as.numeric(nrow(Data))/ssq(X)))})
        dataWeight <- lapply(xx,function(X){X/sqrt(ssq(X))})
        names(dataWeight) <- names(Data)
        return(dataWeight)
    }
)

####### DISCO-SCA APPROACH #######

# scaClass ----------------------------------------------------------------
## Class for Simultaneous Components Analysis Results
## SLOTS:
##     U: (matrix) a matrix whose columns contain the left singular vectors of concatenated data
##     S: (matrix) diagonal matrix  containing the singular values of concatenated data
##     V: (matrix) a matrix whose columns contain the right singular vectors of concatenated data
##     scores: (matrix) scores for SCA decomposition
##     loadings: (matrix) loadings for SCA decomposition
##     ssq: (vector) sum-of-squares of each block
##     VAF: (matrix) VAF for each component in each block and in the contatenated data
setClass(
    Class="scaClass",
    slots=list(U="matrix",
               S="matrix",
               V="matrix",
               scores="matrix",
               loadings="matrix",
               ssqs="numeric",
               VAF="matrix"),
    prototype=list(U=NULL,
                   S=NULL,
                   V=NULL,
                   scores=NULL,
                   loadings=NULL,
                   ssqs=NULL,
                   VAF=NULL)
)

# SCA ------------------------------------------------------------
## Function to calculate the Simultaneous Component Analysis (SCA) and VAF of joined data composed by two blocks
## INPUT:
##     block1: (matrix) data matrix for block1
##     block2: (matrix) data matrix for block2
## OUTPUT:
##     (scaClass) object with the following information
##         U: (matrix) a matrix whose columns contain the left singular vectors of concatenated data
##         S: (matrix) diagonal matrix  containing the singular values of concatenated data
##         V: (matrix) a matrix whose columns contain the right singular vectors of concatenated data
##         scores: (matrix) scores for SCA decomposition
##         loadings: (matrix) loadings for SCA decomposition
##         ssq: (vector) sum-of-squares of each block
##         VAF: (matrix) VAF for each component in each block and in the contatenated data
setGeneric(
    name="SCA",
    def=function(block1,block2){standardGeneric("SCA")}
)

setMethod(
    f="SCA",
    signature=signature(
        block1="matrix",
        block2="matrix"),
    definition=function(block1,block2){
        mat <- rbind(block1,block2)
        svds <- svd(mat)
        scores <- svds$v
        loadings <- svds$u%*%diag(svds$d)
        scares <- list(U=svds$u,S=diag(svds$d),V=svds$v,scores=scores,loadings=loadings)
        ## Sum of squares of each block
        ssq1 <- ssq(block1)
        ssq2 <- ssq(block2)
        ## Calculating VAF
        vaf_block1 <- NULL
        for(i in 1:(ncol(block1)-1)){
            datahat <- scares$S[i,i]*scares$U[1:nrow(block1),i]%*%t(scares$V[,i])
            vaf_block1[i] <- ssq(datahat)/ssq1
        }
        vaf_block2 <- NULL
        for(i in 1:(ncol(block2)-1)){
            datahat <- scares$S[i,i]*scares$U[(nrow(block1)+1):(nrow(block1)+nrow(block2)),i]%*%t(scares$V[,i])
            vaf_block2[i] <- ssq(datahat)/ssq2
        }
        vaf_ow <- (vaf_block1+vaf_block2)/2
        ssqs <- c(ssq1,ssq2)
        names(ssqs) <- c("Block1","Block2")
        VAF <- rbind(vaf_block1,vaf_block2,vaf_ow)
        rownames(VAF) <- c("Block1","Block2","Block1+Block2")
        colnames(VAF) <- paste("C",1:ncol(VAF),sep="")
        ## Assigning names to loadings and scores
        rownames(scores) <- colnames(block1)
        rownames(loadings) <- c(rownames(block1),rownames(block2))
        colnames(scores) <- colnames(loadings) <- paste("Comp.",1:ncol(scores),sep="")
        ## Put info in scaClass
        res <- new("scaClass",U=svds$u,S=diag(svds$d),V=svds$v,scores=scores,loadings=loadings,ssqs=ssqs,VAF=t(VAF))
        return(res)
    }
)

# ssq ---------------------------------------------------------------------
## Function to calculate the sum-of-squares of a given matrix
## INPUT:
##    A: (matrix) input matrix to calculate ssq
## OUTPUT:
##    (numeric) Sum of squares of A matrix
setGeneric(
    name="ssq",
    def=function(A){standardGeneric("ssq")}
)

setMethod(
    f="ssq",
    signature=signature(
        A="matrix"),
    definition=function(A){
        res <- sum(sum(A^2))
        return(res)
    }
)

# PSTRLOSS ----------------------------------------------------------------
# Function to calculate the objective function associated to the pstr2 function
## INPUTS:
##    B: (matrix) rotation matrix
##    mat: (matrix) matrix to rotate
##    Target: (matrix) target to reach
##    W: (matrix) weight matrix (0/1 for resp. unspecified and specified elements of the target)
## OUTPUTS:
##    (numeric) Value of objetive function to minimize [2]
setGeneric(
    name="PSTRLOSS",
    def=function(B,mat,Target,W){standardGeneric("PSTRLOSS")}
)

setMethod(
    f="PSTRLOSS",
    signature=signature(
        B="matrix",
        mat="matrix",
        Target="matrix",
        W="matrix"),
    definition=function(B,mat,Target,W){
        A <- W*(mat%*%B)
        Loss <- ssq(A)
        return(Loss)
    }
)

# PSTR --------------------------------------------------------------------
## Function to calculate an orthogonal rotation matrix (from [2])
## INPUTS:
##    mat: (matrix) matrix to rotate
##    Target: (matrix) target to reach
##    W: (matrix) weight matrix (0/1 for resp. unspecified and specified elements of the target)
##    maxiter: (numeric) maximum number of iterations
##    convergence: (numeric) minimal loss between current and previous iteration
## OUTPUTS:
##    B: (matrix) Orthogonal rotation matrix
##    Loss: (numeric) Value of objetive function to minimize [2]
##    iter: (numeric) Number of necessary iterations for convergence
##    Diff: (numeric) Value of the improvement in the last iteration
setGeneric(
    name="PSTR",
    def=function(mat,Target,W,maxiter,convergence){standardGeneric("PSTR")}
)

setMethod(
    f="PSTR",
    signature=signature(
        mat="matrix",
        Target="matrix",
        W="matrix",
        maxiter="numeric",
        convergence="numeric"),
    definition=function(mat,Target,W,maxiter,convergence){
        ## Random initialization of B
        aux <- svd(matrix(rnorm(ncol(mat)*ncol(mat)),ncol=ncol(mat)))
        Lossc <- PSTRLOSS(aux$u,mat,Target,W)
        Bcurrent <- aux$u
        iter <- 1
        STOP <- 0
        B1 <- t(mat)%*%mat
        alpha <- max(eigen(B1)$values)
        while(STOP==0){
            TargetSter <- mat%*%Bcurrent+W*(Target-mat%*%Bcurrent)
            aux2 <- svd(t(TargetSter)%*%mat)
            B <- aux2$v%*%t(aux2$u)
            if(iter==maxiter){
                STOP <- 1
            }
            Loss <- PSTRLOSS(B,mat,Target,W)
            Diff <- Lossc-Loss
            if (abs(Diff)<convergence){
                STOP <- 1
            }
            iter <- iter+1
            Lossc <- Loss
            Bcurrent <- B
        }
        res <- list(B=B,Loss=Loss,iter=iter,Diff=Diff)
        return(res)
    }
)

# DISCO_SCA ---------------------------------------------------------------
## Function to calculate rotation matrix given the number of exclusive and common components by using the pstr2 function of the official DISCO-SCA page
## INPUT:
##    d1: (numeric) number of Block1 exclusive components
##    d2: (numeric) number of Block2 exclusive components
##    common: (numeric) number of common components between both blocks
##    block1: (matrix) First block of data
##    block2: (matrix) Second block of data
##    maxiter: (numeric) maximum number of iterations (stop criterium)
##    convergence: (numeric) maximum value for convergence (stop criterium)
## OUTPUT:
##    rotMatrix: (matrix) Rotation matrix according the given model
##    loadings1: (matrix) First block loadings
##    loadings2: (matrix) Second block loadings
##    maxdev: (vector) Maximal deviations for each component
##    VAF: (matrix) Proportion of VAF for each component in each block and in the concatenation
##    sca: Results on initial SCA
setGeneric(
    name="DISCO_SCA",
    def=function(d1,d2,common,block1,block2,maxiter=600,convergence=0.0000001){standardGeneric("DISCO_SCA")}
)

setMethod(
    f="DISCO_SCA",
    signature=signature(
        d1="numeric",
        d2="numeric",
        common="numeric",
        block1="matrix",
        block2="matrix"),
    definition=function(d1,d2,common,block1,block2,maxiter,convergence){
        R <- d1+d2+common
        sca_total <- SCA(block1,block2)
        ssq1 <- sca_total@ssqs[1]
        ssq2 <- sca_total@ssqs[2]
        W=rbind(cbind(matrix(0,nrow=nrow(block1),ncol=d1),
                      matrix(1,nrow=nrow(block1),ncol=d2),
                      matrix(0,nrow=nrow(block1),ncol=R-d1-d2)),
                cbind(matrix(1,nrow=nrow(block2),ncol=d1),
                      matrix(0,nrow=nrow(block2),ncol=d2),
                      matrix(0,nrow=nrow(block2),ncol=R-d1-d2)))
        TARGET <- 1-W
        LOSS <- NULL
        BMAT <- NULL
        for(i in 1:2){
            B <- PSTR(mat=sca_total@loadings[,1:R],Target=TARGET,W,maxiter,convergence)$B
            loss <- PSTRLOSS(B,mat=sca_total@loadings[,1:R],Target=TARGET,W)
            LOSS <- cbind(LOSS,loss)
            BMAT <- cbind(BMAT,B)
        }
        k <- which(LOSS==min(LOSS))
        B <- BMAT[,((k[1]-1)*R+1):(k[1]*R)] ## B <- BMAT[((k[1]-1)*R+1):(k[1]*R),]
        Trot <- sca_total@scores[,1:R]%*%B
        Prot <- sca_total@loadings[,1:R]%*%B
        B <- diag(R)
        if(d1>0){
            res1 <- svd(as.matrix(Prot[,1:d1]))
            B[1:d1,1:d1] <- res1$v
        }
        if (d2>0){
            res2 <- svd(as.matrix(Prot[,(d1+1):(d1+d2)]))
            B[(d1+1):(d1+d2),(d1+1):(d1+d2)] <- res2$v
        }
        if (d1+d2<R){
            res3 <- svd(as.matrix(Prot[,(d1+d2+1):R]))
            B[(d1+d2+1):R,(d1+d2+1):R] <- res3$v
        }
        Trot <- Trot[,1:R]%*%B
        Prot <- Prot[,1:R]%*%B
        P1rot <- Prot[1:nrow(block1),]
        P2rot <- Prot[(nrow(block1)+1):(nrow(block1)+nrow(block2)),]
        P1 <- sca_total@loadings[1:nrow(block1),]
        P2 <- sca_total@loadings[(nrow(block1)+1):(nrow(block1)+nrow(block2)),]
        P <- sca_total@loadings
        Ts <- sca_total@scores
        ## Calculate the maximal deviation of distintive and common components
        d1exp <- 0
        d2exp <- 0
        common <- 0
        for(i in 1:R){
            dhat1unr <- sca_total@scores[,i]%*%t(sca_total@loadings[1:nrow(block1),i])
            dhat2unr <- sca_total@scores[,i]%*%t(sca_total@loadings[(nrow(block1)+1):(nrow(block1)+nrow(block2)),i])
            if(i<=R){
                dhat1 <- Trot[,i]%*%t(P1rot[,i])
                dhat2 <- Trot[,i]%*%t(P2rot[,i])
                if(d1>=i){
                    d1exp <- max(d1exp,ssq(dhat1)/ssq1)
                }else if(d1+d2>=i){
                    d2exp <- max(d2exp,ssq(dhat2)/ssq2)
                }
                if(i>d1+d2){
                    common <- max(common,abs((ssq(dhat1)/ssq1)-(ssq(dhat2)/ssq2)))
                }
            }
        }
        maxdev <- c(d1exp,d2exp,common)
        names(maxdev) <- c("Block1","Block2","Common")
        ## Calculate VAF of each component
        if(d1==0 & d2==0){
            if (R==(d1+d2+1)){
                Tc <- as.matrix(Trot[,(d1+d2+1):R]) ## Common scores
                colnames(Tc) <- 1
                P1c <- as.matrix(P1rot[,(d1+d2+1):R]) ## Common loadings Block1
                colnames(P1c) <- 1
                P2c <- as.matrix(P2rot[,(d1+d2+1):R]) ## Common loadings Block2
                colnames(P2c) <- 1
            }else{
                Tc <- Trot[,(d1+d2+1):R] ## Common scores
                colnames(Tc) <- 1:ncol(Tc)
                P1c <- P1rot[,(d1+d2+1):R] ## Common loadings Block1
                colnames(P1c) <- 1:ncol(P1c)
                P2c <- P2rot[,(d1+d2+1):R] ## Common loadings Block2
                colnames(P2c) <- 1:ncol(P2c)
            }
            T1a <- T2a <- P1a <- P2a <- NULL
            ## Common VAF
            ssq_common1 <- ssq_common2 <- NULL
            for(i in 1:(R-d1-d2)){
                ssq_common1 <- c(ssq_common1,ssq(P1c[,i]%*%t(Tc[,i])))
                ssq_common2 <- c(ssq_common2,ssq(P2c[,i]%*%t(Tc[,i])))
            }
            ssq_common <- rbind(ssq_common1,ssq_common2)
            rownames(ssq_common) <- c("Block1","Block2")
            colnames(ssq_common) <- paste("Comp.",1:(R-d1-d2),sep="")
            ssq_block <- NULL
            ssq_cross <- NULL
        } else{
            if (R==(d1+d2+1)){
                Tc <- as.matrix(Trot[,(d1+d2+1):R]) ## Common scores
                colnames(Tc) <- 1
                P1c <- as.matrix(P1rot[,(d1+d2+1):R]) ## Common loadings Block1
                colnames(P1c) <- 1
                P2c <- as.matrix(P2rot[,(d1+d2+1):R]) ## Common loadings Block2
                colnames(P2c) <- 1
            }else{
                Tc <- Trot[,(d1+d2+1):R] ## Common scores
                colnames(Tc) <- 1:ncol(Tc)
                P1c <- P1rot[,(d1+d2+1):R] ## Common loadings Block1
                colnames(P1c) <- 1:ncol(P1c)
                P2c <- P2rot[,(d1+d2+1):R] ## Common loadings Block2
                colnames(P2c) <- 1:ncol(P2c)
            }
            if (d1==1){
                T1a <- as.matrix(Trot[,1:d1]) ## Individual scores Block1
                colnames(T1a) <- 1
                P1a <- as.matrix(P1rot[,1:d1]) ## Individual loadings block1
                colnames(P1a) <- 1
                P2a1 <- as.matrix(P2rot[,1:d1]) ## Cross loadings Block2
                colnames(P2a1) <- 1
            }else{
                T1a <- Trot[,1:d1] ## Individual scores Block1
                colnames(T1a) <- 1:ncol(T1a)
                P1a <- P1rot[,1:d1] ## Individual loadings block1
                colnames(P1a) <- 1:ncol(P1a)
                P2a1 <- P2rot[,1:d1] ## Cross loadings Block2
                colnames(P2a1) <- 1:ncol(P2a1)
            }
            if (d2==1){
                T2a <- as.matrix(Trot[,(d1+1):(d1+d2)]) ## Individual scores Block2
                colnames(T2a) <- 1
                P1a2 <- as.matrix(P1rot[,(d1+1):(d1+d2)]) ## Cross loadings Block1
                colnames(P1a2) <- 1
                P2a <- as.matrix(P2rot[,(d1+1):(d1+d2)]) ## Individual loadings Block2
                colnames(P2a) <- 1
            }else{
                T2a <- Trot[,(d1+1):(d1+d2)] ## Individual scores Block2
                colnames(T2a) <- 1:ncol(T2a)
                P1a2 <- P1rot[,(d1+1):(d1+d2)] ## Cross loadings Block1
                colnames(P1a2) <- 1:ncol(P1a2)
                P2a <- P2rot[,(d1+1):(d1+d2)] ## Individual loadings Block2
                colnames(P2a) <- 1:ncol(P2a)
            }
            ## Common VAF
            ssq_common1 <- ssq_common2 <- NULL
            for(i in 1:(R-d1-d2)){
                ssq_common1 <- c(ssq_common1,ssq(P1c[,i]%*%t(Tc[,i])))
                ssq_common2 <- c(ssq_common2,ssq(P2c[,i]%*%t(Tc[,i])))
            }
            ssq_common <- rbind(ssq_common1,ssq_common2)
            rownames(ssq_common) <- c("Block1","Block2")
            colnames(ssq_common) <- paste("Comp.",1:(R-d1-d2),sep="")
            ## Distinctive VAF
            ssq_block1 <- ssq_block2 <- NULL
            if(d1>0){
                for(i in 1:d1){
                    ssq_block1 <- c(ssq_block1,ssq(P1a[,i]%*%t(T1a[,i])))
                }
            }
            if(d2>0){
                for(i in 1:d2){
                    ssq_block2 <- c(ssq_block2,ssq(P2a[,i]%*%t(T2a[,i])))
                }
            }
            if(d1>d2){
                ssq_block2 <- c(ssq_block2,rep(0,d1-d2))
            } else if (d1<d2){
                ssq_block1 <- c(ssq_block1,rep(0,d2-d1))
            }
            ssq_block <- rbind(ssq_block1,ssq_block2)
            rownames(ssq_block) <- c("Block1","Block2")
            colnames(ssq_block) <- paste("Comp.",1:max(d1,d2),sep="")
            ## Cross VAF
            ssq_cross1 <- ssq_cross2 <- NULL
            if(d2>0){
                for(i in 1:d2){
                    ssq_cross1 <- c(ssq_cross1,ssq(P1a2[,i]%*%t(T2a[,i])))
                }
            }
            if(d1>0){
                for(i in 1:d1){
                    ssq_cross2 <- c(ssq_cross2,ssq(P2a1[,i]%*%t(T1a[,i])))
                }
            }
            if(d1>d2){
                ssq_cross1 <- c(ssq_cross1,rep(0,d1-d2))
            } else if (d1<d2){
                ssq_cross2 <- c(ssq_cross2,rep(0,d2-d1))
            }
            ssq_cross <- rbind(ssq_cross1,ssq_cross2)
            rownames(ssq_cross) <- c("Block1","Block2")
            colnames(ssq_cross) <- paste("Comp.",1:max(d1,d2),sep="")
        }
        ## Creating output
        res <- list(scores=list(common=Tc,dist=list(Block1=T1a,Block2=T2a)),
                    loadings=list(common=list(Block1=P1c,Block2=P2c),dist=list(Block1=P1a,Block2=P2a)),
                    VAF=list(common=ssq_common,dist=list(Block1=ssq_block["Block1",],Block2=ssq_block["Block2",],cross=ssq_cross)),
                    rotMatrix=B,
                    maxdev=maxdev,
                    sca=sca_total)
        return(res)
    }
)

####### JIVE APPROACH #######

# PCA.genes ---------------------------------------------------------------
## INPUT:
##    Data: (matrix) matrix to perform PCA
## OUTPUT:
##    eigen: Eigen values
##    var.exp: Explained variance
##    scores: Matrix of scores
##    loadings: Matrix of loadings
setGeneric(
    name="PCA.genes",
    def=function(Data){standardGeneric("PCA.genes")}
)

setMethod(
    f="PCA.genes",
    signature=signature(
        Data="matrix"),
    definition=function(Data){
        X <- t(Data)
        n<-ncol(X)
        p<-nrow(X)
        offset<-apply(X,2,mean)
        Xoff<-X-(cbind(matrix(1,p,1))%*%rbind(offset))
        #eigen command sorts the eigenvalues in decreasing orden.
        eigen<-eigen(Xoff%*%t(Xoff)/(p-1))
        var<-cbind(eigen$values/sum(eigen$values),cumsum(eigen$values/sum(eigen$values)))
        colnames(var) <- c("%var","%accumvar")
        loadings2<-eigen$vectors
        scores2<-t(Xoff)%*%loadings2
        normas2<-sqrt(apply(scores2^2,2,sum))
        scores1<-loadings2%*%diag(normas2)
        loadings1<-scores2%*%diag(1/normas2)
        rownames(var) <- colnames(loadings1) <- colnames(scores1) <- paste("Comp.",1:ncol(scores1),sep="")
        output<-list(eigen,var,scores1,loadings1)
        names(output)<-c("eigen","var.exp","scores","loadings")
        return(output)
    }
)

# JIVE --------------------------------------------------------------------
## Function to do JIVE analysis of any number of blocks of data
## INPUT:
##    Data: (List) List of data matrices (one for each block and sharing number of columns)
##    r: (numeric) number of common components
##    rIndiv: (vector) Vector with number of individual components of each block (same length than length of Data)
##    ConvergenceThresh: (numeric) Threshold for convergence (accepted LOSS value)
##    maxIter: (numeric) Maximum number of iterations
## OUTPUT:
##    J: (matrix) Joint structure, of dimension (d1+d2+...+dk, n)
##    A: (matrix) Individual structure
setGeneric(
    name="JIVE",
    def=function(Data,r,rIndiv,ConvergenceThresh=10^(-10),maxIter=3000){standardGeneric("JIVE")}
)

setMethod(
    f="JIVE",
    signature=signature(
        Data="list",
        r="numeric",
        rIndiv="numeric"),
    definition=function(Data,r,rIndiv,ConvergenceThresh=10^(-10),maxIter=3000){
        ## Number of blocks (datasets)
        nData <- length(Data)
        ## Getting dimensions of matrices
        dataSizesOr <- as.data.frame(do.call("rbind",lapply(Data,dim)))
        colnames(dataSizesOr) <- c("rows","columns")
        ## Checking if all datasets has the same number o columns (samples)
        n <- unique(dataSizesOr$columns)
        if(length(n)>1){
            stop("Data must have same number of columns(samples)")
        }
        if (length(Data) != length(rIndiv)){
            stop("Number of individual components for each block must be specify")
        }
        ## Dimension reducing transformation for high-dimensional data
        origU <- list()
        dataSizes <- dataSizesOr
        for(j in 1:nData){
            if(dataSizes$rows[j]>n){
                aux <- svd(Data[[j]],nu=ncol(Data[[j]])-1,nv=ncol(Data[[j]])-1)
                Data[[j]] <- diag(aux$d[1:(ncol(Data[[j]])-1)])%*%t(aux$v)
                dataSizes$rows[j] <- nrow(Data[[j]])
                origU[[j]] <- aux$u
            }
        }
        Tot <- do.call("rbind",Data)
        A <- matrix(0,nrow=sum(dataSizes$rows),ncol=n)
        X_est <- 0
        for(iter in 1:maxIter){
            if(r>0){
                aux2 <- svd(Tot,nu=r,nv=r)
                J <- aux2$u%*%diag(aux2$d[1:r],nrow=r)%*%t(aux2$v)
            } else{
                aux2 <- svd(Tot,nu=r,nv=r)
                J <- Tot
            }
            for(k in 1:nData){
                if(rIndiv[k]>0){
                    indexes <- max(1,dataSizes$rows[k-1]+1):sum(dataSizes$rows[1:k])
                    U <- Data[[k]]-J[indexes,]
                    aux3 <- svd(U-U%*%aux2$v%*%t(aux2$v),nu=rIndiv[k],nv=rIndiv[k])
                    A[indexes,] <- aux3$u%*%diag(aux3$d[1:rIndiv[k]],nrow=rIndiv[k])%*%t(aux3$v)
                    Tot[indexes,] <- Data[[k]] - A[indexes,]
                }
            }
            if(norm(X_est-J-A,type="F")^2<ConvergenceThresh){
                break;
            }
            X_est <- J+A
            if(iter==maxIter){
                warning("Maximum number of iterations reached before convergence")
            }
        }
        ## Transform back to original spapce
        origJ <- origA <- NULL
        for(m in 1:nData){
            indexes <- max(1,dataSizes$rows[m-1]+1):sum(dataSizes$rows[1:m])
            Joint <- J[indexes,]
            Indiv <- A[indexes,]
            if(length(origU)!=0){ ## A dimension reducing is applied
                Joint <- origU[[m]]%*%Joint
                Indiv <- origU[[m]]%*%Indiv
            }
            origJ <- rbind(origJ,Joint)
            origA <- rbind(origA,Indiv)
        }
        ## PCA of each matrix to extract the components (scores & loadings)
        pcaJ <- PCA.genes(origJ)
        commonVAF <- as.matrix(pcaJ$var.exp[1:r,"%var"])
        rownames(commonVAF) <- 1:nrow(commonVAF)
        colnames(commonVAF) <- "common"
        A1 <- origA[1:dataSizesOr$rows[1],]
        A2 <- origA[(dataSizesOr$rows[1]+1):sum(dataSizesOr$rows),]
        pcaA1 <- PCA.genes(A1)
        pcaA2 <- PCA.genes(A2)
        if(!all(A==0)){
            distVAF1 <- as.matrix(pcaA1$var.exp[1:rIndiv[1],"%var"])
            distVAF2 <- as.matrix(pcaA2$var.exp[1:rIndiv[2],"%var"])
            rownames(distVAF1) <- 1:nrow(distVAF1)
            rownames(distVAF2) <- 1:nrow(distVAF2)
            colnames(distVAF1) <- colnames(distVAF2) <-"distintive"
        }else{
            distVAF <- NULL
        }
        res <- list(J=J,
                    A=A,
                    commonComps=r,
                    distComps=rIndiv,
                    scores=list(common=pcaJ$scores[,1:r],
                                dist=list(Block1=as.matrix(pcaA1$scores[,1:rIndiv[1]]),Block2=as.matrix(pcaA2$scores[,1:rIndiv[2]]))),
                    loadings=list(common=list(Block1=pcaJ$loadings[1:dataSizes$rows[1],1:r],
                                              Block2=pcaJ$loadings[dataSizes$rows[1]:nrow(pcaJ$loadings),1:r]),
                                  dist=list(Block1=pcaA1$loadings[,1:rIndiv[1]],Block2=pcaA2$loadings[,1:rIndiv[2]])),
                    VAF=list(common=commonVAF,dist=list(Block1=distVAF1,Block2=distVAF2)))
        return(res)
    }
)

####### O2PLS APPROACH #######

# PSVD --------------------------------------------------------------------
## Calculates rank R approximation of svd of AB' -- SVD of covariance matrix X_1'X_2
# Calculates a base for the predictive space
## INPUT:
##    A: (matrix) First block of data
##    B: (matrix) Second block of data
##    R: (numeric) rank.
## OUTPUT:
##    A: (matrix) Left singular vectors of the decomposition
##    S: (matrix) Diagonal matrix of decomposition eigen values
##    B: (matrix) Rigth singular vector of the decomposition
setGeneric(
    name="PSVD",
    def=function(A,B,R){standardGeneric("PSVD")}
)

setMethod(
    f="PSVD",
    signature=signature(
        A="matrix",
        B="matrix",
        R="numeric"),
    definition=function(A,B,R){
        qrA <- qr(A,0)
        qrB <- qr(B,0)
        svd1 <- svd(qr.R(qrA)%*%t(qr.R(qrB)))
        res <- list(A=qr.Q(qrA)%*%svd1$u,S=diag(svd1$d),B=qr.Q(qrB)%*%svd1$v)
        return(res)
    }
)

# O2PLS -------------------------------------------------------------------
## Function to perform O2PLS analysis
## INPUT:
##    X1: (matrix) First block of data
##    X2: (matrix) Second block of data
##    Rcommon: (numeric) Number of common components
##    Rspecific: (numeric) Number of specific components
## OUTPUT:
##    weights: (list) Weigth matrices
##
##    scores: (matrix) Scores matrix
##    loadings: (matrix) LOadings matrix
##    vaf: (matrix)
setGeneric(
    name="O2PLS",
    def=function(X1,X2,Rcommon,Rspec){standardGeneric("O2PLS")}
)

setMethod(
    f="O2PLS",
    signature=signature(
        X1="matrix",
        X2="matrix",
        Rcommon="numeric",
        Rspec="numeric"),
    definition=function(X1,X2,Rcommon,Rspec){
        ## The code is performed for matrices with samples in rows
        X <- t(X1)
        Y <- t(X2)
        LV <- Rcommon
        LVX <- Rspec[1]
        LVY <- Rspec[2]
        if(LVX==0){
            T_Yosc <- matrix(0,nrow=nrow(Y),ncol=LVX)
            P_Yosc <- matrix(0,nrow=ncol(X),ncol=LVX)
            W_Yosc <- matrix(0,nrow=ncol(X),ncol=LVX)
        }
        if(LVY==0){
            T_Yosc <- matrix(0,nrow=nrow(X),ncol=LVY)
            P_Yosc <- matrix(0,nrow=ncol(Y),ncol=LVY)
            W_Yosc <- matrix(0,nrow=ncol(Y),ncol=LVY)
        }
        X_true <- X
        Y_true <- Y
        svd1 <- svd(t(Y)%*%X,nu=LV,nv=LV)
        T_Yosc <- NULL
        P_Yosc <- NULL
        W_Yosc <- NULL
        ## Remove orthogonal components from X sequentially
        for(lv in 1:LVX){
            TT <- X%*%svd1$v
            E_XY <- X-TT%*%t(svd1$v)
            svd2 <- svd(t(E_XY)%*%TT,nu=1,nv=1)
            t_Yosc <- X%*%svd2$u
            p_Yosc <- (t(X)%*%t_Yosc)/as.numeric((t(t_Yosc)%*%t_Yosc))
            X <- X - t_Yosc%*%t(p_Yosc)
            ## Collect the orthogonal components
            T_Yosc <- cbind(T_Yosc,t_Yosc)
            P_Yosc <- cbind(P_Yosc,p_Yosc)
            W_Yosc <- cbind(W_Yosc,svd2$u)
        }
        ## Update T again (since X is changed)
        TT <- X%*%svd1$v
        U_Xosc <- NULL
        P_Xosc <- NULL
        C_Xosc <- NULL
        ## calculate LVY orthogonal Y components
        for(lv in 1:LVY){
            UU <- Y%*%svd1$u
            F_XY <- Y-UU%*%t(svd1$u)
            svd3 <- svd(t(F_XY)%*%UU,nu=1,nv=1)
            u_Xosc <- Y%*%svd3$u
            p_Xosc <- (t(Y)%*%u_Xosc)/as.numeric((t(u_Xosc)%*%u_Xosc))
            Y <- Y-u_Xosc%*%t(p_Xosc)
            ## Collect the orthogonal components
            U_Xosc <- cbind(U_Xosc,u_Xosc)
            P_Xosc <- cbind(P_Xosc,p_Xosc)
            C_Xosc <- cbind(C_Xosc,svd3$u)
        }
        ## Update U again (since Y is changed)
        UU <- Y%*%svd1$u
        ## Repeat steps 1,2,4 before step 6
        svd4 <- svd(t(Y)%*%X,nu=LV,nv=LV)
        TT <- X%*%svd4$v
        UU <- Y%*%svd4$u
        B_U <- solve(t(UU)%*%UU)%*%t(UU)%*%TT
        B_T <- solve(t(TT)%*%TT)%*%t(TT)%*%UU
        ## Other model components
        EE <- X_true-TT%*%t(svd4$v)-T_Yosc%*%t(P_Yosc)
        FF <- Y_true-UU%*%t(svd4$u)-U_Xosc%*%t(P_Xosc)
        H_TU <- TT-UU%*%B_U
        H_UT <- UU-TT%*%B_T
        Y_hat <- TT%*%B_T%*%t(svd4$u)
        X_hat <- UU%*%B_U%*%t(svd4$v)
        ## Statistics
        ssqX_true <- ssq(X_true)
        ssqY_true <- ssq(Y_true)
        R2X <- 1-(ssq(EE)/ssqX_true)
        R2Y <- 1-(ssq(FF)/ssqY_true)
        R2Xcorr <- ssq(TT%*%t(svd4$v))/ssqX_true
        R2Ycorr <- ssq(UU%*%t(svd4$u))/ssqY_true
        R2X_Y0 <- ssq(T_Yosc%*%t(P_Yosc))/ssqX_true
        R2Y_X0 <- ssq(U_Xosc%*%t(P_Xosc))/ssqY_true
        R2Xhat <- 1-(ssq(UU%*%B_U%*%t(svd4$v)-X_true)/ssqX_true)
        R2Yhat <- 1-(ssq(TT%*%B_T%*%t(svd4$u)-Y_true)/ssqY_true)
        res <- list(scores=list(common=list(Block1=TT,Block2=UU),
                                dist=list(Block1=T_Yosc,Block2=U_Xosc)),
                    loadings=list(common=list(Block1=W_Yosc,Block2=C_Xosc),
                                  dist=list(Block1=P_Yosc,Block2=P_Xosc)))
        return(res)
    }
)

# omicsCompAnalysis -------------------------------------------------------
## INPUT:
##    ######...: Set of expression sets one for each omic dataset ## De momento este no que no se como declararlo
##    Input: (list) List of expression sets one for each omic dataset
##    Names: (vector) Names assigned to each omic dataset
##    method=c("DISCOSCA","JIVE","O2PLS"): (character) Method used for the analysis
##    Rcommon: (numeric) NUmber of common components
##    Rspecific: (vector) Number of specific components, one for each omic dataset
##    convThres: (numeric) Convergence threshold (only used with DISCOSCA and O2PLS)
##    maxIter: (numeric) Maximum number of iterations (only used with DISCOSCA and O2PLS)
##    center=c("PERBLOCKS","ALLBLOCKS"): Has center be applied before analysis?
##    scale=c("PERBLOCKS","ALLBLOCKS"): Has scale be applied before analysis?
##    weight: (logical): Have blocks be weighted?
## OUTPUT:
##    An object of class 'caClass'

#' @export
#' @title Components analysis for multiple objects
#' @aliases omicsCompAnalysis,list,character,character,numeric,numeric-method
#' @description
#' This function performs a components analysis of object wise omics data to understand the mechanisms that underlay all the data blocks under study (common mechanisms) and the mechanisms underlying each of the data block independently (distinctive mechanisms). This analysis include both, the preprocessing of data and the components analysis by using three different methodologies.
#' @usage omicsCompAnalysis(Input, Names, method, Rcommon, Rspecific, 
#'                          convThres=1e-10, maxIter=600, center=FALSE, 
#'                          scale=FALSE, weight=FALSE)
#' 
#' @param Input List of \code{ExpressionSet} objects, one for each block of data.
#' @param Names Character vector giving names for each Input object.
#' @param method Method to use for analysis (either "DISCOSCA", "JIVE", or "O2PLS").
#' @param Rcommon Number of common components between all blocks
#' @param Rspecific Vector giving number of unique components for each input block
#' @param convThres Stop criteria for convergence
#' @param maxIter Maximum number of iterations
#' @param center Character (or FALSE) specifying which (if any) centering will be applied before analysis. Choices are "PERBLOCKS" (each block separately) or "ALLBLOCKS" (all data together).
#' @param scale Character (or FALSE) specifying which (if any) scaling will be applied before analysis. Choices are "PERBLOCKS" (each block separately) or "ALLBLOCKS" (all data together).
#' @param weight Logical indicating whether weighting is to be done.
#' 
#' @return An object of class \code{\link{caClass-class}}.
#' 
#' @author Patricia Sebastian Leon
#' 
#' @examples
#' data("STATegRa_S3")
#' B1 <- createOmicsExpressionSet(Data=Block1.PCA,pData=ed.PCA,
#' pDataDescr=c("classname"))
#' B2 <- createOmicsExpressionSet(Data=Block2.PCA,
#'                                pData=ed.PCA,pDataDescr=c("classname"))
#' # Omics components analysis
#' discoRes <- omicsCompAnalysis(Input=list(B1,B2),Names=c("expr","mirna"),
#'                               method="DISCOSCA",Rcommon=2,Rspecific=c(2,2),
#'                               center=TRUE,scale=TRUE,weight=TRUE)
#' jiveRes <- omicsCompAnalysis(Input=list(B1,B2),Names=c("expr","mirna"),
#'                              method="JIVE",Rcommon=2,Rspecific=c(2,2),
#'                              center=TRUE,scale=TRUE,weight=TRUE)
#' o2plsRes <- omicsCompAnalysis(Input=list(B1,B2),Names=c("expr","mirna"),
#'                               method="O2PLS",Rcommon=2,Rspecific=c(2,2),
#'                               center=TRUE,scale=TRUE,weight=TRUE)
setGeneric(
    name="omicsCompAnalysis",
    def=function(Input,Names,method,Rcommon,Rspecific,convThres=1e-10,maxIter=600,center=FALSE,scale=FALSE,weight=FALSE){
        standardGeneric("omicsCompAnalysis")}
)

setMethod(
    f="omicsCompAnalysis",
    signature=signature(
        Input="list",
        Names="character",
        method="character",
        Rcommon="numeric",
        Rspecific="numeric"),
    definition=function(Input,Names,method,Rcommon,Rspecific,convThres,maxIter,center,scale,weight){
        ## STEP1: Reading and checking provided data
        ## Expression sets
        orData <- Input
        if(!all(do.call("c",lapply(orData,FUN=class))=="ExpressionSet")){
            warning("Please provide a list of expression sets")
            return();
        }
        ## Converting expression sets list in matrices list
        Data <- lapply(orData,exprs)
        ## Number of Expression sets
        nData <- length(Data)
        ## Expression sets names
        if (is.null(Names)){
            Names <- paste("Block",1:nData)
        }
        ## Check is number of columns/samples is the same in all datasets
        if(length(unique(do.call("c",lapply(Data,FUN=ncol))))!=1){
            warning("Number of samples is not the same for all data sets")
            return();
        }
        prepros <- NULL
        ## STEP2: Removing constant rows if any
        Data <- lapply(Data,function(xx){
            index <- which(apply(xx,1,sd)==0)
            if(length(index)>0){
                print(paste("The following features are eliminated because are constant:",length(index),"\n",paste(names(index),collapse=", ")))
                xx <- xx[-index,]
            }
            return(xx)
        })
        ## STEP3: Center & scale datasets if selected
        if(center & scale){
            prepros <- c(prepros,"centered & scaled")
            Data <- lapply(Data,function(x){t(scale(t(x),center=TRUE,scale=TRUE))})
        }else{
            if(center){
                prepros <- c(prepros,"centered")
                Data <- lapply(Data,function(x){t(scale(t(x),center=TRUE,scale=FALSE))})
            }
            ## STEP3: Scale datasets if selected
            if(scale){
                prepros <- c(prepros,"scaled")
                Data <- lapply(Data,function(x){t(scale(t(x),center=FALSE,scale=TRUE))})
            }
        }
        ## STEP4: Weight datasets if selected
        if(weight){
            prepros <- c(prepros,"weighted")
            Data <- weightBlocks(Data)
        }
        if(is.null(prepros)){
            prepros <- "none"
        }
        ## STEP5: Apply selected method to data
        method <- match.arg(method,choices=c("DISCOSCA","JIVE","O2PLS"))
        if(method=="DISCOSCA"){
            auxres <- DISCO_SCA(d1=Rspecific[1],d2=Rspecific[2],common=Rcommon,block1=Data[[1]],block2=Data[[2]],maxiter=maxIter,
                                convergence=convThres)
            res <- new("caClass",
                       InitialData=orData,
                       Names=Names,
                       preprocessing=prepros,
                       preproData=Data,
                       caMethod="DISCO-SCA",
                       commonComps=Rcommon,
                       distComps=Rspecific,
                       scores=auxres$scores,
                       loadings=auxres$loadings,
                       VAF=auxres$VAF,
                       others=list(rotMatrix=auxres$rotMatrix,maxdev=auxres$maxdev,sca=auxres$sca))
        } else if (method=="JIVE"){
            auxres <- JIVE(Data=Data,r=Rcommon,rIndiv=Rspecific,ConvergenceThresh=convThres,maxIter=maxIter)
            res <- new("caClass",
                       InitialData=orData,
                       Names=Names,
                       preprocessing=prepros,
                       preproData=Data,
                       caMethod="JIVE",
                       commonComps=Rcommon,
                       distComps=Rspecific,
                       scores=auxres$scores,
                       loadings=auxres$loadings,
                       VAF=auxres$VAF,
                       others=list(J=auxres$J,A=auxres$A))
        } else if (method=="O2PLS"){
            auxres <- O2PLS(X1=Data[[1]],X2=Data[[2]],Rcommon=Rcommon,Rspec=Rspecific)
            res <- new("caClass",
                       InitialData=orData,
                       Names=Names,
                       preprocessing=prepros,
                       preproData=Data,
                       caMethod="O2PLS",
                       commonComps=Rcommon,
                       distComps=Rspecific,
                       scores=auxres$scores,
                       loadings=auxres$loadings,
                       VAF=list(),
                       others=list())
        }
        return(res)
    }
)

####### MODEL SELECTION #######

# selectCommonComps ------------------------------------------------------------
## INPUT:
## X: (matrix) First block of data
## Y: (matrix) Second block of data
## Rmax: (numeric) Maximum number of common components
## OUTPUT:
## (numeric) Optimal number of common components
#' @import ggplot2
#' @import MASS
#' @export
#' @title Select common components in two data blocks
#' @aliases selectCommonComps,matrix,matrix,numeric-method
#' @description 
#' This function applies a Simultaneous Component Analysis (SCA). The idea is that the scores for both blocks should have a similar behaviour if the components are in the common mode. Evaluation is by the ratios between the explained variances (SSQ) of each block and the estimator. The highest component count with 0.8 < ratio < 1.5 is selected.
#' @usage selectCommonComps(X, Y, Rmax)
#' @param X Matrix of omics data
#' @param Y Matrix of omics data
#' @param Rmax Maximum number of common components to find
#' @return A list with components:
#' \describe{
#'      \item{common}{Optimal number of common components}
#'      \item{ssqs}{Matrix of SSQ for each block and estimator}
#'      \item{pssq}{\code{\link{ggplot}} object showing SSQ for each block and estimator}
#'      \item{pratios}{\code{\link{ggplot}} object showing SSQ ratios between each block and estimator}
#' }
#' @author Patricia Sebastian-Leon
#' @examples
#' data(STATegRa_S3)
#' cc <- selectCommonComps(X=Block1.PCA, Y=Block2.PCA, Rmax=3)
#' cc$common
#' cc$pssq
#' cc$pratios
setGeneric(
    name="selectCommonComps",
    def=function(X,Y,Rmax){standardGeneric("selectCommonComps")}
)

setMethod(
    f="selectCommonComps",
    signature=signature(
        X="matrix",
        Y="matrix",
        Rmax="numeric"),
    definition=function(X,Y,Rmax){
        ## Center and scale matrices
        X <- scale(t(X),center=TRUE,scale=FALSE)
        X <- t(X/sqrt(ssq(X)))
        Y <- scale(t(Y),center=TRUE,scale=FALSE)
        Y <- t(Y/sqrt(ssq(Y)))
        cb <- rbind(X,Y)
        svd1 <- svd(cb,nu=Rmax,nv=Rmax)
        ssqs <- NULL
        for(i in 1:Rmax){
            if(i ==1){
                loadings <- as.matrix(svd1$u[,1:i]*diag(svd1$d)[1:i,1:i])
            }else{
                loadings <- svd1$u[,1:i]%*%diag(svd1$d)[1:i,1:i]
            }
            l1 <- as.matrix(loadings[1:nrow(X),]) ## Block1 loadings
            l2 <- as.matrix(loadings[(nrow(X)+1):(nrow(X)+nrow(Y)),]) ## Block2 loadings
            ## Calculate scores per block from the loadings
            scoresX <- ginv(t(l1)%*%l1)%*%t(l1)%*%X
            scoresY <- ginv(t(l2)%*%l2)%*%t(l2)%*%Y
            ## Estimate common blocks using their scores
            eX <- l1%*%scoresX
            eY <- l2%*%scoresY
            ## Estimate common blocks using their counterpart loadings
            rX <- l1%*%scoresY
            rY <- l2%*%scoresX
            ssqs <- rbind(ssqs,c(ssq(eX),ssq(eY),ssq(rX),ssq(rY)))
        }
        rownames(ssqs) <- 1:Rmax
        colnames(ssqs) <- c("estX","estY","estXY","estYX")
        ## Plot ssqs for estimated blocks
        dfssq <- data.frame(ssq=c(ssqs),comps=rep(rownames(ssqs),ncol(ssqs)),block=rep(colnames(ssqs),each=nrow(ssqs)))
        p <- ggplot(dfssq,aes(x=comps,y=ssq,fill=block))+
            geom_bar(stat="identity",position="dodge")+
            xlab("Number of common components")+
            ylab("ssq of estimated block")
        ## Determine acceptable ratios for allowed differences between 1st and 2nd block predictions
        r13 <- ssqs[,1]/ssqs[,3]
        r24 <- ssqs[,2]/ssqs[,4]
        ssqratio <- data.frame(ratio=c(r13,r24),comp=c(names(r13),names(r24)),block=c(rep("ratioX",length(r13)),rep("ratioY",length(r24))))
        p2 <- ggplot(ssqratio,aes(x=comp,y=ratio,fill=block))+
            geom_bar(stat="identity",position="dodge")+
            xlab("Number of common components")+
            ylab("Ratio between SSQ/estSSQ")+
            geom_hline(aes(yintercept=0.8),colour="black",linetype="dashed")+
            geom_hline(aes(yintercept=1.2),colour="black",linetype="dashed")
        xx13 <- r13>0.8 & r13<1.25
        xx24 <- r24>0.8 & r24<1.25
        indexes <- which((r13>0.8 & r13<1.25) & (r24>0.8 & r24<1.25))
        if (length(indexes)==0){
            common <- 0
        }else {
            common <- max(indexes)
        }
        res <- list(common=common,ssqs=ssqs,pssq=p,pratios=p2)
        return(res)
    }
)

# PCA.selection -----------------------------------------------------------
## INPUT:
## Data: (matrix) matrix with the omic data
## fac.sel: (character) Criterium for selecting number of components c("%accum", "single%", "rel.abs", "fixed.num")
## varthreshold: (numeric) Threshold for the selection of components in "%accum", "single%" criteriums
## nvar: (numeric) Threshold applied when the option "abs.val" is selected
## PCnum: (numeric) Fixed number of components to select with the option"fixed.num"
## OUTPUT: (list)
## - PCAres: (list) Results of the PCA of the data
## - numComps: (numeric) Number of selected components applying the selected criterium
#' @export
#' @title Select an optimal number of components using PCA
#' @aliases PCA.selection,matrix,character-method
#' @description 
#' Selects the optimal number of components from data using PCA. There are four different criteria available: accumulated variance explained, individual explained variance of each component, absolute value of variability or fixed number of components.
#' @usage PCA.selection(Data, fac.sel, varthreshold=NULL, nvar=NULL, PCnum=NULL)
#' @param Data Data matrix (with samples in columns and features in rows)
#' @param fac.sel Selection criteria ("\%accum", "single\%", "rel.abs", "fixed.num")
#' @param varthreshold Threshold for "\%accum" or "single\%" criteria
#' @param nvar Threshold for "rel.abs"
#' @param PCnum Fixed number of components for "fixed.num"
#' @return List containing:
#' \describe{
#'      \item{PCAres}{List containing results of PCA, with fields "eigen", "var.exp", "scores" and "loadings"}
#'      \item{numComps}{Number of components selected}
#' }
#' @author Patricia Sebastian Leon
#' @examples 
#' data(STATegRa_S3)
#' ps <- PCA.selection(Data=Block2.PCA, fac.sel="single%", varthreshold=0.03)
#' ps$numComps
setGeneric(
    name="PCA.selection",
    def=function(Data,fac.sel,varthreshold=NULL,nvar=NULL,PCnum=NULL){standardGeneric("PCA.selection")}
)

setMethod(
    f="PCA.selection",
    signature=signature(
        Data="matrix",
        fac.sel="character"),
    definition=function(Data,fac.sel,varthreshold,nvar,PCnum){
        fac.sel <- match.arg(fac.sel, c("%accum", "single%", "rel.abs", "fixed.num"))
        pca.sel <- PCA.genes(Data)
        eigen <- pca.sel$eigen$values
        tot.var <- sum(eigen)
        rank <- length(which(eigen > 1e-16))
        if (fac.sel == "%accum") {
            fac <- max(length(which(pca.sel$var.exp[,2] <= varthreshold)),1)
        } else if (fac.sel == "single%"){
            fac <- length(which(pca.sel$var.exp[,1] >= (varthreshold)))
        } else if (fac.sel == "rel.abs"){
            mean.expl.var <- tot.var/ nrow(Data)
            fac <- length(which(eigen >= (mean.expl.var*nvar)))
        } else if (fac.sel == "fixed.num"){
            fac <- PCnum
        }
        res <- list(PCAres=pca.sel,numComps=fac)
        return(res)
    }
)

# modelSelection ----------------------------------------------------------
## INPUT:
## Input: (list) List of ExpressionSets
## Rmax: (numeric) Number of maximum common components
## (*) fac.sel: (character) Criterium for selecting number of components c("%accum", "single%", "rel.abs", "fixed.num")
## (*) varthreshold: (numeric) Threshold for the selection of components in "%accum", "single%" criteriums
## (*) nvar: (numeric) Threshold applied when the option "rel.abs" is selected
## (*) PCnum: (numeric) Fixed number of components to select with the option"fixed.num"
## OUTPUT: (list)
## - common: Number of optimal common components
## - dist: Number of optimal distictive components for each block
#' @export
#' @title Find optimal common and distinctive components
#' @aliases modelSelection,list,numeric,character-method
#' @description 
#' Uses \code{\link{selectCommonComps}} and \code{\link{PCA.selection}} to estimate the optimal number of common and distinctive components according to given selection criteria.
#' @usage modelSelection(Input, Rmax, fac.sel, varthreshold=NULL, nvar=NULL, PCnum=NULL)
#' @param Input List of two ExpressionSet objects
#' @param Rmax Maximum common components (see \code{\link{selectCommonComps}})
#' @param fac.sel PCA criteria (see \code{\link{PCA.selection}})
#' @param varthreshold Cumulative variance criteria (see \code{\link{PCA.selection}})
#' @param nvar Relative variance criteria (see \code{\link{PCA.selection}})
#' @param PCnum Fixed component number (see \code{\link{PCA.selection}})
#' @return List containing:
#' \describe{
#'      \item{common}{Number of common components}
#'      \item{dist}{Number of distinct components per input block}
#' }
#' @author Patricia Sebastian-Leon
#' @seealso \code{\link{selectCommonComps}},\code{\link{PCA.selection}},\code{\link{omicsCompAnalysis}}
#' @examples
#' data(STATegRa_S3)
#' B1 <- createOmicsExpressionSet(Data=Block1.PCA,pData=ed.PCA,pDataDescr=c("classname"))
#' B2 <- createOmicsExpressionSet(Data=Block2.PCA,pData=ed.PCA,pDataDescr=c("classname"))
#' ms <- modelSelection(Input=list(B1, B2), Rmax=4, fac.sel="single\%", varthreshold=0.03)
#' ms
setGeneric(
    name="modelSelection",
    def=function(Input,Rmax,fac.sel,varthreshold=NULL,nvar=NULL,PCnum=NULL){standardGeneric("modelSelection")}
)

setMethod(
    f="modelSelection",
    signature=signature(
        Input="list",
        Rmax="numeric",
        fac.sel="character"),
    definition=function(Input,Rmax,fac.sel,varthreshold,nvar,PCnum){
        X <- exprs(Input[[1]])
        Y <- exprs(Input[[2]])
        ## For the moment we did not cross validate de number of individual components!!
        n1 <- PCA.selection(X,fac.sel=fac.sel,varthreshold=varthreshold,nvar=nvar,PCnum=PCnum)$numComps
        n2 <- PCA.selection(Y,fac.sel=fac.sel,varthreshold=varthreshold,nvar=nvar,PCnum=PCnum)$numComps
        if (Rmax > min(n1,n2)){
            Rmax <- min(n1,n2)
            warning(paste("Rmax cannot be higher than the minimum of components selected for each block. Rmax fixed to:",Rmax))
        }
        common <- selectCommonComps(X,Y,Rmax)$common
        ## Calculate number of optimal components
        res <- list(common=common,dist=c(n1-common,n2-common))
        return(res)
    }
)


