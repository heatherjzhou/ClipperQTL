#include <RcppArmadillo.h>
using namespace arma;


//Given X, the n*K predictor matrix (without a column of ones),
//get the design matrix (with a column of ones).
mat getDesignMatrix(const mat X){

  int n=X.n_rows; //Sample size.
  vec vecOfOnes(n,fill::ones); //n*1.
  mat D=join_rows(vecOfOnes,X); //n*(K+1).

  return D;
}


//Given Y, the gene by sample expression matrix (not residualized, no informational columns),
//XRaw, the SNP by sample genotype matrix (not residualized, no informational columns),
//dataCovariates, the sample by covariate covariate matrix (without a column of ones),
//geneTSSs, the vector of gene TSSs, corresponding to the rows of Y,
//SNPPositions, the vector of SNP positions, corresponding to the rows of X,
//cisDistance,
//approach,
//and B,
//get tableMaxAbsCors.
//[[Rcpp::export()]]
arma::mat getTableMaxAbsCorsCpp(const arma::mat Y,
                                const arma::mat XRaw, //To be residualized.
                                const arma::mat dataCovariates,
                                const arma::vec geneTSSs,
                                const arma::vec SNPPositions,
                                const int cisDistance,
                                const std::string approach,
                                const int B){

  //Create tableMaxAbsCors. To be filled and returned.
  mat tableMaxAbsCors(Y.n_rows,1+B); //257*21. Each row corresponds to a gene. The columns are: exp, bg1, ..., bg20.

  //Get the design matrix.
  mat D=getDesignMatrix(dataCovariates); //515*53. Design matrix.
  // mat H=D*(D.t()*D).i()*D.t(); //515*515. Projection matrix. H=D*(D^T*D)^(-1)*D^T. Symmetric matrix. Not calculating the projection matrix in advance cuts the run time in about half.

  //Create YCube.
  //If approach is "standard", then the first slice is YResid, and each of the remaining slices is a YPermResid.
  //If approach is "Clipper", then the first slice is YResid, and each of the remaining slices is a YResidPerm.
  cube YCube(Y.n_rows,Y.n_cols,1+B); //257*515*21.
  if(approach=="standard"){ //Permute first and then residualize.
    for(int indexOfExperAndBg=0;indexOfExperAndBg<(1+B);indexOfExperAndBg++){
      if(indexOfExperAndBg==0){
        // mat YResid=Y-Y*H; //257*515.
        mat YResid=Y-Y*D*(D.t()*D).i()*D.t(); //257*515.
        YCube.slice(indexOfExperAndBg)=YResid;
      }else{
        mat YPerm=shuffle(Y,1); //1 means to shuffle each row.
        // mat YPermResid=YPerm-YPerm*H; //257*515.
        mat YPermResid=YPerm-YPerm*D*(D.t()*D).i()*D.t(); //257*515.
        YCube.slice(indexOfExperAndBg)=YPermResid;
      }
    }
  }else if(approach=="Clipper"){ //Residualize first and then permute.
    // mat YResid=Y-Y*H; //257*515.
    mat YResid=Y-Y*D*(D.t()*D).i()*D.t(); //257*515.
    for(int indexOfExperAndBg=0;indexOfExperAndBg<(1+B);indexOfExperAndBg++){
      if(indexOfExperAndBg==0){
        YCube.slice(indexOfExperAndBg)=YResid;
      }else{
        YCube.slice(indexOfExperAndBg)=shuffle(YResid,1); //1 means to shuffle each row.
      }
    }
  }

  //Residualize XRaw against the covariates. Use the following sequence of steps for the fastest speed. Inspired by residualizeCpp().
  //Using two steps instead of only one step reduces the run time by about 4 seconds.
  //Ensuring that the smaller matrices are multiplied together first in the calculation of coef_T reduces the run time slightly.
  //The same strategy does not make residualizing Y faster (in fact, it makes it slightly slower).
  mat coef_T=XRaw*(D*(D.t()*D).i());
  mat X=XRaw-coef_T*D.t(); //131,682*515.

  //Given YCube and X, fill tableMaxAbsCors.
  for(int indexOfGene=0;indexOfGene<Y.n_rows;indexOfGene++){
    //Get YCubeSub.
    mat YCubeSub=YCube.row(indexOfGene); //515*21. Sample by experAndBg.

    //Get XSub.
    int geneTSS=geneTSSs(indexOfGene); //75,786,699.
    uvec indicesOfSNPs=find(abs(SNPPositions-geneTSS)<=cisDistance); //"u" stands for unsigned integer.
    mat XSub=X.rows(indicesOfSNPs); //7261*515.

    //Calculate absCorMatrix.
    mat absCorMatrix=abs(cor(YCubeSub,XSub.t())); //21*7261.

    //Fill tableMaxAbsCors.
    // tableMaxAbsCors.row(indexOfGene)=max(absCorMatrix,1).t(); //1 means to return the statistic for each row. Chose not to use this because max(,1) cannot handle NAs correctly when the first entry of a row is NA.
    for(int indexOfExperAndBg=0;indexOfExperAndBg<(1+B);indexOfExperAndBg++){
      tableMaxAbsCors(indexOfGene,indexOfExperAndBg)=max(absCorMatrix.row(indexOfExperAndBg)); //max() can handle NAs correctly when the argument is a vector.
    }
  }

  return tableMaxAbsCors;
}








