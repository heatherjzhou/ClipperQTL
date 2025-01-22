


prepareMethodParameters<-function(approach,B,sampleSize){
  # approach<-"standard"
  # B<-1000
  # sampleSize<-515
  # #To comment out.

  if(approach=="standard"){
    if(B<1000){
      cat('\n\n\n\n\n!!!Warning: if approach="standard", we strongly recommend using B>=1000.\n') #Include 5 empty lines before and after the warning to emphasize it.
      cat('\n\n\n\n\nContinuing with the user-specified B...\n') #Include 5 empty lines before and after the warning to emphasize it.
    }
  }else if(approach=="Clipper"){
    if(sampleSize<450){
      cat('\n\n\n\n\n!!!Warning: since the sample size is under 450, we strongly recommend using approach="standard".\n')
      cat('\n\n\n\n\nContinuing with approach="Clipper"...\n')
    }
    if(!B%in%c(1,20:100)){
      cat('\n\n\n\n\n!!!Warning: if approach="Clipper", we strongly recommend using B=1 or B between 20 and 100.\n')
      cat('\n\n\n\n\nContinuing with the user-specified B...\n')
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
# prepareMethodParameters(approach="test",B=1000,sampleSize=500) #Test print messages.


