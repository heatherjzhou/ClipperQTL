



Rcpp::sourceCpp("~/2022.03.14_ClipperQTL/ClipperQTL/test/2023.02.13_testNABehavior/23.02.18_func1.1_testFunctions.cpp")

#max(,1) cannot handle NAs correctly when the first entry of a row is NA.
if(TRUE){
  X<-matrix(data=sample(1:9,9),nrow=3,byrow=TRUE)
  X
  testFunction1(X)

  X[1,1]<-NA
  X[2,2]<-NA
  X[3,3]<-NA
  testFunction1(X)
}

#max() can handle NAs correctly when the argument is a vector.
if(TRUE){
  x<-sample(1:9,9)
  x
  testFunction2(x)

  x[1]<-NA
  testFunction2(x)
}





