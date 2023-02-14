#include <RcppArmadillo.h>
#include "3.4_residualizeCpp.h"
using namespace arma;


//Given Y, gene by sample expression matrix (not residualized, no informational columns),
//X, SNP by sample genotype matrix (residualized, no informational columns),
//dataCovariates, sample by covariate covariate matrix,
//geneTSSs, vector of gene TSSs, corresponding to the rows of Y,
//SNPPositions, vector of SNP positions, corresponding to the rows of X,
//and B,
//get tableMaxAbsCors.
//[[Rcpp::export()]]
mat getTableMaxAbsCorsCpp(const mat Y,
                          const mat X,
                          const mat dataCovariates,
                          const vec geneTSSs,
                          const vec SNPPositions,
                          const int B){

  // auto start=std::chrono::high_resolution_clock::now();

  //Create tableMaxAbsCors. To be filled and returned.
  mat tableMaxAbsCors(Y.n_rows,1+B); //257*21. Each row corresponds to a gene. The columns are: exp, bg1, ..., bg20.

  // //Create YCube, where the first slice is Y, and each of the remaining slices is a YPerm.
  // cube YCube(Y.n_rows,Y.n_cols,1+B); //257*515*21.
  // for(int indexOfExperAndBg=0;indexOfExperAndBg<(1+B);indexOfExperAndBg++){
  //   if(indexOfExperAndBg==0){
  //     YCube.slice(indexOfExperAndBg)=Y;
  //   }else{
  //     YCube.slice(indexOfExperAndBg)=shuffle(Y,1); //1 means to shuffle each row.
  //   }
  // }

  //Create YCube, where the first slice is YResid, and each of the remaining slices is a YPermResid.
  cube YCube(Y.n_rows,Y.n_cols,1+B); //257*515*21.
  for(int indexOfExperAndBg=0;indexOfExperAndBg<(1+B);indexOfExperAndBg++){
    if(indexOfExperAndBg==0){
      mat YResid=residualizeCpp(Y.t(),dataCovariates); //257*515.
      YCube.slice(indexOfExperAndBg)=YResid;
    }else{
      mat YPermResid=residualizeCpp(shuffle(Y,1).t(),dataCovariates); //257*515. 1 means to shuffle each row.
      YCube.slice(indexOfExperAndBg)=YPermResid;
    }
  }

  //Fill tableMaxAbsCors.
  for(int indexOfGene=0;indexOfGene<Y.n_rows;indexOfGene++){
    // cout<<"indexOfGene="<<indexOfGene<<"...\n";

    //Get YCubeSub.
    mat YCubeSub=YCube.row(indexOfGene); //515*21. Sample by experAndBg.

    //Get XSub.
    int geneTSS=geneTSSs(indexOfGene); //75,786,699.
    uvec indicesOfSNPs=find(abs(SNPPositions-geneTSS)<=1e6); //"u" stands for unsigned integer.
    mat XSub=X.rows(indicesOfSNPs); //7261*515.

    //Calculate absCorMatrix.
    mat absCorMatrix=abs(cor(YCubeSub,XSub.t())); //21*7261.

    //Fill tableMaxAbsCors.
    // tableMaxAbsCors.row(indexOfGene)=max(absCorMatrix,1).t(); //1 means to return the statistic for each row. Chose not to use this because max(,1) cannot handle NAs correctly when the first entry of a row is NA.
    for(int indexOfExperAndBg=0;indexOfExperAndBg<(1+B);indexOfExperAndBg++){
      tableMaxAbsCors(indexOfGene,indexOfExperAndBg)=max(absCorMatrix.row(indexOfExperAndBg)); //max() can handle NAs correctly when the argument is a vector.
    }
  }

  // auto stop=std::chrono::high_resolution_clock::now();
  // auto duration=std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
  // cout<<"The C++ code took "<<duration.count()<<" milliseconds.\n";

  return tableMaxAbsCors;
}








