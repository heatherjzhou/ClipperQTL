


prepareMethodParameters<-function(approach,B,sampleSize){
  # approach<-"standard"
  # B<-1000
  # sampleSize<-515
  # #To comment out.

  if(approach=="standard"){ #Default is approach="standard" and B=1000.
    if(B<1000){
      message('\nWarning: if approach="standard", we strongly recommend using B>=1000.\n') #Include a new line at the beginning of this print message to make it stand out from the rest.
      message('Continuing with the user-specified B...\n')
    }
  }else if(approach=="Clipper"){
    if(sampleSize<450){
      message('\nWarning: since the sample size is under 450, we strongly recommend using approach="standard".\n') #Include a new line at the beginning of this print message to make it stand out from the rest.
      message('Continuing with approach="Clipper"...\n')
    }

    if(is.null(B)){
      B<-20
    }else if(B<20||B>100){
      message('\nWarning: if approach="Clipper", we strongly recommend using B between 20 and 100.\n') #Include a new line at the beginning of this print message to make it stand out from the rest.
      message('Continuing with the user-specified B...\n')
    }
  }else{
    stop('Approach must be "standard" or "Clipper".\n')
  }

  toReturn<-list(approach=approach,B=B)
  return(toReturn)
}


# prepareMethodParameters(approach="standard",B=500,sampleSize=200) #Test print messages.
# prepareMethodParameters(approach="Clipper",B=20,sampleSize=200) #Test print messages.
# prepareMethodParameters(approach="Clipper",B=10,sampleSize=500) #Test print messages.
# prepareMethodParameters(approach="Clipper",B=200,sampleSize=500) #Test print messages.
# prepareMethodParameters(approach="test",B=1000,sampleSize=500) #Test print messages.


