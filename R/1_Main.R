#' Test function
#'
#' @export

testFunction<-function(){
  X<-matrix(c(1,3,1,4),nrow=2,byrow=TRUE)
  y<-cbind(c(1,2))
  myFastLM(X,y)
}
