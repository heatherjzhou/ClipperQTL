#include <RcppArmadillo.h>
using namespace arma;


//Given X, the n*K predictor matrix (without a column of ones),
//get the design matrix (with a column of ones).
mat getDesignMatrix2(const mat X){ //Name as getDesignMatrix2() to avoid duplicate function names, which would cause package compilation to fail.

  int n=X.n_rows; //Sample size.
  vec vecOfOnes(n,fill::ones); //n*1.
  mat D=join_rows(vecOfOnes,X); //n*(K+1).

  return D;
}


//Given Y, the gene by sample expression matrix (not residualized, no informational columns),
//X, the SNP by sample genotype matrix (not residualized, no informational columns),
//dataCovariates, the sample by covariate covariate matrix (without a column of ones),
//geneTSSs, the vector of gene TSSs, corresponding to the rows of Y,
//SNPPositions, the vector of SNP positions, corresponding to the rows of X,
//cisDistance, e.g., 1e6,
//sigGeneSNPPairMethod,
//absCorThresholds,
//percent,
//indexOfChunk,
//outputDirChunk,
//save temporary result files, one for each gene.
//[[Rcpp::export()]]
int saveTemporaryResultFilesCpp(const arma::mat Y,
                                const arma::mat X,
                                const arma::mat dataCovariates,
                                const arma::vec geneTSSs,
                                const arma::vec SNPPositions,
                                const int cisDistance, //For filtering SNPs for each gene.
                                const std::string sigGeneSNPPairMethod,
                                const arma::vec absCorThresholds,
                                const double percent,
                                const int indexOfChunk, //For naming temporary result files, e.g., Chunk5_gene1.csv.
                                const std::string outputDirChunk){

  // auto start=std::chrono::high_resolution_clock::now();

  //Get the design matrix.
  mat D=getDesignMatrix2(dataCovariates); //515*53. Design matrix.
  // mat H=D*(D.t()*D).i()*D.t(); //515*515. Projection matrix. H=D*(D^T*D)^(-1)*D^T. Symmetric matrix.

  //Get YResid (see getTableMaxAbsCorsCpp()).
  mat YResid=Y-Y*D*(D.t()*D).i()*D.t(); //127*515.

  //Get XResid (see getTableMaxAbsCorsCpp()).
  mat coef_T=X*(D*(D.t()*D).i());
  mat XResid=X-coef_T*D.t(); //131,649*515.

  // auto stop=std::chrono::high_resolution_clock::now();
  // auto duration=std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
  // cout<<"Residualizing Y and X took "<<duration.count()<<" milliseconds.\n";

  // start=std::chrono::high_resolution_clock::now();

  //For each gene, save a temporary result file.
  for(int indexOfGene=0;indexOfGene<Y.n_rows;indexOfGene++){
    //Get yResid.
    arma::rowvec yResid=YResid.row(indexOfGene); //1*515.

    //Get XResidSub.
    int geneTSS=geneTSSs(indexOfGene); //75,796,882.
    uvec indicesOfSNPs=find(abs(SNPPositions-geneTSS)<=cisDistance); //"u" stands for unsigned integer.
    mat XResidSub=XResid.rows(indicesOfSNPs); //7278*515.

    //Calculate cors and absCors.
    mat corMatrix=cor(yResid.t(),XResidSub.t()); //1*7278.
    vec cors=corMatrix.t(); //7278*1.
    vec absCors=abs(cors); //7278*1.

    uvec indicesOfSNPs2; //Indices of significant SNPs among the local SNPs of the gene.
    if(sigGeneSNPPairMethod=="FastQTL"){
      indicesOfSNPs2=find(absCors>=absCorThresholds(indexOfGene));
    }else{
      vec temp={1-percent/100}; //Convert double to vec because quantile() cannot take a double. Use 1-percent/100 because larger absolute correlations are more significant.
      vec thresholdVec=quantile(absCors,temp); //Save the result in a vector first because it cannot be saved in a double directly.
      indicesOfSNPs2=find(absCors>=thresholdVec(0));
    }

    uvec indicesFinal=indicesOfSNPs(indicesOfSNPs2); //Indices of significant SNPs among all SNPs in X.
    vec corsFinal=cors(indicesOfSNPs2);

    //Create toSave, which contains two columns.
    //The first column is the indices of significant SNPs among all SNPs in X.
    //The second column is the partial correlations.
    mat toSave=join_rows(conv_to<vec>::from(indicesFinal)+1,corsFinal); //Need to convert uvec to vec before using join_rows(). Add 1 for use in R.

    std::string filename=outputDirChunk+"Chunk"+std::to_string(indexOfChunk)+"_Gene"+std::to_string(indexOfGene+1)+".csv"; //Add 1 for use in R. The file name is Chunk5_gene1.csv, for example.
    toSave.save(filename,csv_ascii); //Save in CSV format without a header.
  }

  // stop=std::chrono::high_resolution_clock::now();
  // duration=std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
  // cout<<"Looping through the genes took "<<duration.count()<<" milliseconds.\n";

  return 0;
}








