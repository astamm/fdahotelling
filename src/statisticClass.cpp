#include "statisticClass.h"

namespace fdahotelling
{
    
    void Statistic::SetInput1(const MatrixType& x)
    {
        m_DataMatrix1 = x;
        m_ModifiedData = true;
    }
    
    void Statistic::SetInput2(const MatrixType& x)
    {
        m_DataMatrix2 = x;
        m_ModifiedData = true;
    }

    void Statistic::SetDelta0(const VectorType& x)
    {
        m_Delta0 = x;
        m_ModifiedData = true;
    }

    void Statistic::Update()
    {
        if (!m_ModifiedData)
            return;

        m_NumberOfVariables = m_DataMatrix1.n_cols;

        // Check number of variables
        // If no variable, no possible inference
        if (m_NumberOfVariables < 1)
            Rcpp::stop("Data must have at least one variable.");

        m_SampleSize1 = m_DataMatrix1.n_rows;
        m_SampleSize2 = m_DataMatrix2.n_rows;
        m_TwoSampleTest = (m_SampleSize2 != 0);

        // Check sample sizes:
        // - for two-sample tests, inference requires
        // Na>= 1, Nb>=1 and Na+Nb>=3
        // - for one-sample tests, inference requires N>=2
        if (m_TwoSampleTest)
        {
            if (m_DataMatrix2.n_cols != m_NumberOfVariables)
                Rcpp::stop("The second dataframe needs to have the same number of variables (columns).");

            if (m_SampleSize1 < 1 || m_SampleSize2 < 1 || m_SampleSize1 + m_SampleSize2 < 3)
                Rcpp::stop("At least one of the two groups does not have enough samples for making inference (Na>=1, Nb>=1, Na+Nb>=3).");
        }
        else
        {
            if (m_SampleSize1 < 2)
                Rcpp::stop("Not enough samples for making inference (N>=2).");
        }
        
        if (m_UsePairedTest && m_TwoSampleTest && m_SampleSize1 != m_SampleSize2)
            Rcpp::stop("The two samples should have the same number of observations for a paired test.");
        
        if (m_UsePairedTest && m_TwoSampleTest)
        {
            m_DataMatrix1 -= m_DataMatrix2;
            m_DataMatrix2.clear();
            m_TwoSampleTest = false;
        }
        
        // TO DO: Set delta vector, is it still needed?
        if (m_Delta0.size() == 1)
        {
            double delta0 = arma::as_scalar(m_Delta0);
            m_Delta0 = arma::ones<VectorType>(m_NumberOfVariables);
            m_Delta0 *= delta0;
        }
        
        // Compute (difference between) sample mean(s)
        m_Delta = arma::mean(m_DataMatrix1, 0).t();
        if (m_TwoSampleTest)
            m_Delta -= arma::mean(m_DataMatrix2, 0).t();
        
        // Compute (pooled) sample covariance matrix
        m_CovarianceMatrix = arma::cov(m_DataMatrix1);

        if (m_TwoSampleTest)
        {
            m_CovarianceMatrix *= (m_SampleSize1 - 1.0);
            m_CovarianceMatrix += (m_SampleSize2 - 1.0) * arma::cov(m_DataMatrix2);
            m_CovarianceMatrix /= (m_SampleSize1 + m_SampleSize2 - 2.0);
        }
        
        // Account for basis functions
        if (m_StepSize != 0.0)
        {
            double halfStep = std::sqrt(m_StepSize / 2.0);
            double wholeStep = std::sqrt(m_StepSize);
            for (unsigned int i = 0;i < m_NumberOfVariables;++i)
            {
                double iMultVal = 1.0;
                if (i == 0 || i == m_NumberOfVariables-1)
                {
                    m_Delta[i] *= halfStep;
                    m_Delta0[i] *= halfStep;
                    iMultVal = halfStep;
                }
                else
                {
                    m_Delta[i] *= wholeStep;
                    m_Delta0[i] *= wholeStep;
                    iMultVal = wholeStep;
                }
                
                for (unsigned int j = i;j < m_NumberOfVariables;++j)
                {
                    double jMultVal = 1.0;
                    if (j == 0 || j == m_NumberOfVariables-1)
                        jMultVal = halfStep;
                    else
                        jMultVal = wholeStep;
                    
                    m_CovarianceMatrix(i,j) *= iMultVal * jMultVal;
                    
                    if (i != j)
                        m_CovarianceMatrix(j,i) = m_CovarianceMatrix(i,j);
                }
            }
        }

        m_ModifiedData = false;
    }

    void Statistic::SetPseudoInverse()
    {
        // Compute pseudo-inverse of the covariance matrix
        // using Divide and Conquer method
        MatrixType eigvec;
        arma::eig_sym(m_EigenValues, eigvec, m_CovarianceMatrix, "dc");

        // Find acceptable range of eigenvalues to keep
        unsigned int sampleSize2 = (m_TwoSampleTest) ? m_SampleSize2 : 1;
        int startPos = m_NumberOfVariables - (m_SampleSize1 + sampleSize2 - 2);
        bool traditionalForm = (startPos < 0);
        if (traditionalForm)
            startPos = 0;
        unsigned int endPos = m_NumberOfVariables - 1;

        double trSigma = 0;
        for (unsigned int i = startPos;i <= endPos;++i)
            trSigma += m_EigenValues[i];

        while (m_EigenValues[startPos] < m_Tolerance * trSigma && startPos < (int)endPos - 1)
            ++startPos;

        m_EigenValues = m_EigenValues.subvec(startPos,endPos);
        eigvec = eigvec.cols(startPos,endPos);

        m_PseudoInverse.set_size(m_NumberOfVariables,m_NumberOfVariables);
        m_PseudoInverse.fill(0.0);

        for (unsigned int i = 0;i < m_NumberOfVariables;++i)
        {
            for (unsigned int j = i;j < m_NumberOfVariables;++j)
            {
                for (unsigned int k = 0;k < endPos-startPos+1;++k)
                    m_PseudoInverse(i,j) += eigvec(i,k) * eigvec(j,k) / m_EigenValues[k];

                if (i != j)
                    m_PseudoInverse(j,i) = m_PseudoInverse(i,j);
            }
        }
    }

    void Statistic::SetTraces()
    {
        if (m_NumberOfVariables == 1)
        {
            m_TraceOfSigma = m_EigenValues[0];
            m_TraceOfSigmaSquare = m_EigenValues[0] * m_EigenValues[0];
            return;
        }
        
        // Taken from Srivastava, 2005
        unsigned int rank = m_EigenValues.size();

        m_TraceOfSigma = 0.0;
        m_TraceOfSigmaSquare = 0.0;
        for (unsigned int i = 0;i < rank;++i)
        {
            double eigenValue = m_EigenValues[i];
            m_TraceOfSigma += eigenValue;
            m_TraceOfSigmaSquare += eigenValue * eigenValue;
        }

        m_TraceOfSigmaSquare -= m_TraceOfSigma * m_TraceOfSigma / rank;
        m_TraceOfSigmaSquare *= (rank * rank) / ((rank - 1.0) * (rank + 2.0));
    }

    double Statistic::GetHotellingStatistic() const
    {
        double constant = 1.0;
        if (m_UseCorrection)
        {
            unsigned int rank = m_EigenValues.size();
            unsigned int df1 = std::min(rank, m_NumberOfVariables);
            unsigned int df2 = (m_NumberOfVariables > rank) ? m_NumberOfVariables - rank : rank - m_NumberOfVariables;
            df2++;
            constant = m_TraceOfSigma * m_TraceOfSigma * df2 / (m_TraceOfSigmaSquare * rank * m_NumberOfVariables * df1);
        }

        double F = 0;
        for (unsigned int i = 0;i < m_NumberOfVariables;++i)
        {
            double idiff = m_Delta[i] - m_Delta0[i];
            F += m_PseudoInverse(i,i) * idiff * idiff;
            for (unsigned int j = i+1;j < m_NumberOfVariables;++j)
                F += 2.0 * m_PseudoInverse(i,j) * idiff * (m_Delta[j] - m_Delta0[j]);
        }
        
        if (m_TwoSampleTest)
            F *= m_SampleSize1 * m_SampleSize2 / (m_SampleSize1 + m_SampleSize2);
        else
            F *= m_SampleSize1;

        return constant * F;
    }

    double Statistic::GetL1Statistic() const
    {
        return arma::as_scalar(arma::sum(arma::abs(m_Delta - m_Delta0)));
    }

    double Statistic::GetL2Statistic() const
    {
        VectorType deltaSq = arma::square(m_Delta - m_Delta0);
        return arma::as_scalar(arma::sum(deltaSq));
    }

    double Statistic::GetLinfStatistic() const
    {
        return arma::as_scalar(arma::max(arma::abs(m_Delta - m_Delta0)));
    }

    double Statistic::GetStandardizedL1Statistic() const
    {
        VectorType deltaSq = arma::square(m_Delta - m_Delta0);
        VectorType variances = arma::diagvec(m_CovarianceMatrix);
        return arma::as_scalar(arma::sum(arma::sqrt(deltaSq / variances)));
    }

    double Statistic::GetStandardizedL2Statistic() const
    {
        VectorType deltaSq = arma::square(m_Delta - m_Delta0);
        VectorType variances = arma::diagvec(m_CovarianceMatrix);
        return arma::as_scalar(arma::sum(deltaSq / variances));
    }

    double Statistic::GetStandardizedLinfStatistic() const
    {
        VectorType deltaSq = arma::square(m_Delta - m_Delta0);
        VectorType variances = arma::diagvec(m_CovarianceMatrix);
        return arma::as_scalar(arma::max(arma::sqrt(deltaSq / variances)));
    }

    Rcpp::NumericVector Statistic::GetValue(const statistic_choice statChoice)
    {
        this->Update();

        Rcpp::NumericVector resVal;

        switch (statChoice)
        {
            case Hotelling:
            {
                this->SetPseudoInverse();

                if (m_UseCorrection)
                    this->SetTraces();

                resVal = Rcpp::NumericVector::create(Rcpp::Named("Hotelling") = this->GetHotellingStatistic());
                break;
            }

            case L1:
            {
                resVal = Rcpp::NumericVector::create(Rcpp::Named("L1") = this->GetL1Statistic());
                break;
            }

            case L2:
            {
                resVal = Rcpp::NumericVector::create(Rcpp::Named("L2") = this->GetL2Statistic());
                break;
            }

            case Linf:
            {
                resVal = Rcpp::NumericVector::create(Rcpp::Named("Linf") = this->GetLinfStatistic());
                break;
            }

            case StdL1:
            {
                resVal = Rcpp::NumericVector::create(Rcpp::Named("StdL1") = this->GetStandardizedL1Statistic());
                break;
            }

            case StdL2:
            {
                resVal = Rcpp::NumericVector::create(Rcpp::Named("StdL2") = this->GetStandardizedL2Statistic());
                break;
            }

            case StdLinf:
            {
                resVal = Rcpp::NumericVector::create(Rcpp::Named("StdLinf") = this->GetStandardizedLinfStatistic());
                break;
            }

            case All:
            default:
            {
                this->SetPseudoInverse();

                if (m_UseCorrection)
                    this->SetTraces();

                resVal = Rcpp::NumericVector::create(Rcpp::Named("Hotelling") = this->GetHotellingStatistic(), Rcpp::Named("L1") = this->GetL1Statistic(), Rcpp::Named("L2") = this->GetL2Statistic(), Rcpp::Named("Linf") = this->GetLinfStatistic(), Rcpp::Named("StdL1") = this->GetStandardizedL1Statistic(), Rcpp::Named("StdL2") = this->GetStandardizedL2Statistic(), Rcpp::Named("StdLinf") = this->GetStandardizedLinfStatistic());
                break;
            }
        }

        return resVal;
    }

    Statistic::MatrixType Statistic::GetPseudoInverse()
    {
        this->Update();
        this->SetPseudoInverse();
        return m_PseudoInverse;
    }

    Statistic::VectorType Statistic::GetEigenValues()
    {
        this->Update();
        this->SetPseudoInverse();
        return m_EigenValues;
    }

} // end namespace fdahotelling
