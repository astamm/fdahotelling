#include "random.h"

arma::mat
rmvnorm(const unsigned int n,
        const arma::colvec& mean,
        const arma::mat& squareRootSigma)
{
    unsigned int numComponents = mean.size();
    arma::mat randomMatrix = arma::randn<arma::mat>(n,numComponents);
    
    arma::mat resVal(n,numComponents);
    for (unsigned int i = 0;i < n;++i)
    {
        for (unsigned int j = 0;j < numComponents;++j)
        {
            double tmpVal = 0;
            for (unsigned int k = 0;k < numComponents;++k)
                tmpVal += randomMatrix(i,k) * squareRootSigma(k,j);
            
            resVal(i,j) = mean[j] + tmpVal;
        }
    }
    
    return resVal;
}
