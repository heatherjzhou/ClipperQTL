


prepareMethodParameters<-function(approach,B,sampleSize){
  # approach<-NULL
  # B<-NULL
  # sampleSize<-515
  # #To comment out.

  #Prepare approach.
  if(is.null(approach)){
    if(sampleSize<450){
      approach<-"standard"
    }else{
      approach<-"Clipper"
    }
  }else if(approach=="standard"){
    #Do nothing.
  }else if(approach=="Clipper"){
    if(sampleSize<450){
      message('\nWarning: since the sample size is under 450, we strongly recommend using approach="standard".\n') #Include a new line at the beginning of this print message to make it stand out from the rest.
      message('\nContinuing with approach="Clipper"...\n') #Include a new line at the beginning of this print message to make it stand out from the rest.
    }
  }else{
    stop('approach must be NULL, "standard", or "Clipper".\n')
  }

  #Prepare B.
  if(is.null(B)){
    if(approach=="standard"){
      B<-1000
    }else if(approach=="Clipper"){
      B<-20
    }
  }

  toReturn<-list(approach=approach,B=B)
  return(toReturn)
}
