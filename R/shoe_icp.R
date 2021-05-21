
icp <- function(mesh1, mesh2, iterations=3,lm1=NULL, lm2=NULL, uprange=1, maxdist=NULL, minclost=50, distinc=0.5, rhotol=pi, k=50, reflection=FALSE, silent=FALSE,subsample=NULL,subsampletype=c("km","pd"),type=c("rigid","similarity","affine"),getTransform=FALSE,pcAlign=FALSE,pcOptim=TRUE,threads=0,weights=NULL) {

  meshorig <- mesh1 <- vcgUpdateNormals(mesh1,silent=silent)
  mesh2 <- vcgUpdateNormals(mesh2)
  if (pcAlign) {
    mesh1 <- pcAlign(mesh1,mesh2,optim=pcOptim)
    if (!is.null(lm1))
      lm1 <- applyTransform(lm1,computeTransform(mesh2,mesh1))
  }
  mysample <- NULL

  KDtree <- vcgCreateKDtreeFromBarycenters(mesh2)
  starticks <- 10

  type <- match.arg(type,c("rigid","similarity","affine"))
  if (!is.null(lm1) && !pcAlign){## perform initial rough registration
    trafo <- computeTransform(lm2,lm1,type=type,reflection=reflection,weights=weights)
    mesh1 <- applyTransform(mesh1,trafo)
  }
  ## create subsample to speed up registration
  if (!is.null(subsample)) {
    subsampletype <- match.arg(subsampletype[1],c("pd","km"))
    if (subsampletype == "pd")
      mysample <- Rvcg::vcgSample(mesh1,type=subsampletype,SampleNum=subsample,MCsamp = 20)
    else
      mysample <- fastKmeans(mesh1,k=subsample,threads=threads)$centers
    mysample <- vcgClostKD(mysample,mesh1,threads=threads)
  }
  origsample <- mysample
  if (!silent)
    cat(paste0("\n performing ",type," registration\n\n") )
  count <- 0
  while (count < iterations) {
    if (!silent) {
      if ((count %% starticks)  == 0 && count != 0)
        cat(paste0(" ",count," "))
      if ((count %% 50)  == 0 && count != 0)
        cat("\n")
      cat("*")
    }
    copymesh <- mesh1
    if (!is.null(subsample) ) {
      minclost <- min(minclost,subsample)
      copymesh <- mysample
    }

    proMesh <- vcgClostOnKDtreeFromBarycenters(KDtree,copymesh,sign=F,k=k,threads=threads) ## project mesh1 onto mesh2
    x1 <- vert2points(copymesh)
    x2 <- vert2points(proMesh)
    dists <- abs(proMesh$quality)
    good <- 1:nrow(x1)

    ## check if normals angles are below rhotol
    if (rhotol < pi) {
      normchk <- normcheck(copymesh,proMesh,threads)
      goodnorm <- which(normchk < rhotol)
      x1 <- x1[goodnorm,]
      x2 <- x2[goodnorm,]
      dists <- dists[goodnorm]
      good <- 1:nrow(x1)
    }
    if (!is.null(maxdist) || uprange < 1) {
      ## check distances of remaining points and select points
      if (is.null(maxdist)) {
        qud <- quantile(dists,probs=uprange)
        good <- which(dists <= qud)
      } else {
        qud <- maxdist
        good <- which(dists <= qud)
        increase <- distinc
        while (length(good) < minclost) {
          good <- which(dists <= (qud+increase))
          if (!silent)
            cat(paste("distance increased to",qud+increase,"\n"))
          increase <- increase+distinc
        }
      }
    }

    #if reverse is true apply the switch of the transform to get the reversal
    ## get transform for current iteration
    trafo <- computeTransform(x1[good,],x2[good,],type=type) #applying mesh 2 to mesh 1 (the fixed mesh object) and giving the transformation:

    ## apply transformation to mesh1 if no subsampling
    if (is.null(subsample)) {
      mesh2 <- applyTransform(mesh2,trafo)
    }
    if (!is.null(subsample)) { #right now we are not working with this, so it can be changed later
      ## hack until changes from Morpho::applyTransform are published
      ntrafo <- trafo
      ntrafo[1:3,4] <- 0
      orignorms <- mysample$normals
      orignorms[1:3,] <- t(applyTransform(t(orignorms[1:3,]),ntrafo))
      mysample <- applyTransform(mysample,trafo)
      mysample$normals <- orignorms
    }
    count <- count+1
  }
  if (!is.null(subsample)) {
    trafo <- computeTransform(mysample,origsample)
    mesh2 <- applyTransform(mesh2,trafo)
  }
  if (!silent) {
    if ((count %% 50)  == 0 && count != 0)
      cat(paste0(" ",count," \n"))
    cat("\n")
  }
  if (getTransform) {
    trafo <- computeTransform(vert2points(mesh1),vert2points(meshorig),type=type)
    if (!is.null(lm1))
      lm1 <- applyTransform(lm1,trafo)
    return(list(mesh=mesh1,transform=trafo,landmarks=lm1, distance=distinc)) #also add get distance to return from the last iteration
  } else {
    return(list(mesh1,dists))
  }
}

icp_rev <- function(mesh1, mesh2, iterations=3,lm1=NULL, lm2=NULL, uprange=1, maxdist=NULL, minclost=50, distinc=0.5, rhotol=pi, k=50, reflection=FALSE, silent=FALSE,subsample=NULL,subsampletype=c("km","pd"),type=c("rigid","similarity","affine"),getTransform=FALSE,pcAlign=FALSE,pcOptim=TRUE,threads=0,weights=NULL) {

  mesh1 <- vcgUpdateNormals(mesh1,silent=silent)
  meshorig <- mesh2 <- vcgUpdateNormals(mesh2)
  mysample <- NULL

  KDtree <- vcgCreateKDtreeFromBarycenters(mesh1)
  starticks <- 10

  type <- match.arg(type,c("rigid","similarity","affine"))

  origsample <- mysample
  if (!silent)
    cat(paste0("\n performing ",type," registration\n\n") )
  distances <- rep(NA, iterations)
  count <- 0
  while (count < iterations) {
    if (!silent) {
      if ((count %% starticks)  == 0 && count != 0)
        cat(paste0(" ",count," "))
      if ((count %% 50)  == 0 && count != 0)
        cat("\n")
      cat("*")
    }
    copymesh <- mesh2

    proMesh <- vcgClostOnKDtreeFromBarycenters(KDtree,copymesh,sign=F,k=k,threads=threads) ## project mesh1 onto mesh2
    x1 <- vert2points(copymesh)
    x2 <- vert2points(proMesh)
    dists <- abs(proMesh$quality)
    distances[count+1] <- mean(dists)

    good <- 1:nrow(x1)

    ## check if normals angles are below rhotol
    if (rhotol < pi) {
      normchk <- normcheck(copymesh,proMesh,threads)
      goodnorm <- which(normchk < rhotol)
      x1 <- x1[goodnorm,]
      x2 <- x2[goodnorm,]
      dists <- dists[goodnorm]
      good <- 1:nrow(x1)
    }
    if (!is.null(maxdist) || uprange < 1) {
      ## check distances of remaining points and select points
      if (is.null(maxdist)) {
        qud <- quantile(dists,probs=uprange)
        good <- which(dists <= qud)
      } else {
        qud <- maxdist
        good <- which(dists <= qud)
        increase <- distinc
        while (length(good) < minclost) {
          good <- which(dists <= (qud+increase))
          if (!silent)
            cat(paste("distance increased to",qud+increase,"\n"))
          increase <- increase+distinc
        }
      }
    }
    ## get transform for current iteration
    trafo <- computeTransform(x2[good,],x1[good,],type=type)

    ## apply transformation to mesh1 if no subsampling
    if (is.null(subsample)) {
      mesh2 <- applyTransform(mesh2,trafo)
    }

    count <- count+1
  }

  if (!silent) {
    if ((count %% 50)  == 0 && count != 0)
      cat(paste0(" ",count," \n"))
    cat("\n")
  }
  if (getTransform) {
    trafo <- computeTransform(vert2points(mesh1),vert2points(meshorig),type=type)
    if (!is.null(lm1))
      lm1 <- applyTransform(lm1,trafo)
    return(list(mesh=mesh1,transform=trafo,landmarks=lm1))
  } else {
    return(list(mesh=mesh2, dists = distances))
  }
}
