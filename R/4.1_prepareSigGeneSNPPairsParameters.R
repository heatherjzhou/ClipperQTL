


prepareSigGeneSNPPairsParameters<-function(approach,sigGeneSNPPairMethod,percent){
  # approach<-"standard"
  # sigGeneSNPPairMethod<-NULL
  # percent<-0.1
  # #To comment out.

  if(approach=="standard"){
    if(is.null(sigGeneSNPPairMethod)){
      sigGeneSNPPairMethod<-"FastQTL"
    }
  }else if(approach=="Clipper"){
    if(is.null(sigGeneSNPPairMethod)){
      sigGeneSNPPairMethod<-"topPercent"
    }else if(sigGeneSNPPairMethod=="FastQTL"){
      stop('If the Clipper approach is used for eGene identification, then sigGeneSNPPairMethod can only be "topPercent".\n')
    }
  }

  if(!sigGeneSNPPairMethod%in%c("FastQTL","topPercent")){
    stop('sigGeneSNPPairMethod must be "FastQTL" or "topPercent".\n')
  }

  toReturn<-list(sigGeneSNPPairMethod=sigGeneSNPPairMethod,percent=percent)
  return(toReturn)
}

# prepareSigGeneSNPPairsParameters(approach="standard",sigGeneSNPPairMethod=NULL,percent=0.1) #Test the function.
# prepareSigGeneSNPPairsParameters(approach="Clipper",sigGeneSNPPairMethod=NULL,percent=0.1) #Test the function.
# prepareSigGeneSNPPairsParameters(approach="Clipper",sigGeneSNPPairMethod="FastQTL",percent=0.1) #Test the function.
# prepareSigGeneSNPPairsParameters(approach="Clipper",sigGeneSNPPairMethod="test",percent=0.1) #Test the function.


