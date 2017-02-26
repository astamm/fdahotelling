#include <RcppArmadillo.h>

namespace fdahotelling
{
    class Statistic {
    public:
        typedef arma::mat MatrixType;
        typedef arma::colvec VectorType;
        
        typedef enum
        {
            Hotelling = 0,
            L1,
            L2,
            Linf,
            StdL1,
            StdL2,
            StdLinf,
            All
        } statistic_choice;
        
        Statistic()
        {
            m_DataMatrix1.clear();
            m_DataMatrix2.clear();
            m_CovarianceMatrix.clear();
            m_PseudoInverse.clear();
            m_Delta.clear();
            m_Delta0.clear();
            m_EigenValues.clear();
            
            m_NumberOfVariables = 1;
            m_SampleSize1 = 1;
            m_SampleSize2 = 0;
            
            m_TwoSampleTest = false;
            m_UseCorrection = false;
            m_ModifiedData = false;
            
            m_Tolerance = 0;
            m_TraceOfSigma = 0;
            m_TraceOfSigmaSquare = 0;
            m_StepSize = 0.0; // vector case, not functional by default
        }
        
        ~Statistic() {}
        
        void SetInput1(const MatrixType& x);
        void SetInput2(const MatrixType& x);
        void SetDelta0(const VectorType& x);
        void SetUsePairedTest(const bool x) {m_UsePairedTest = x;};
        void SetStepSize(const double x) {m_StepSize = x;};
        void SetUseCorrection(const bool x) {m_UseCorrection = x;};
        void SetTolerance(const double x) {m_Tolerance = x;};
        void Update();
        
        Rcpp::NumericVector GetValue(const statistic_choice statChoice);
        MatrixType GetPseudoInverse();
        VectorType GetEigenValues();
        
    protected:
        double GetHotellingStatistic() const;
        double GetL1Statistic() const;
        double GetL2Statistic() const;
        double GetLinfStatistic() const;
        double GetStandardizedL1Statistic() const;
        double GetStandardizedL2Statistic() const;
        double GetStandardizedLinfStatistic() const;
        void SetPseudoInverse();
        void SetTraces();
        
    private:
        MatrixType m_DataMatrix1, m_DataMatrix2;
        MatrixType m_CovarianceMatrix, m_PseudoInverse;
        VectorType m_Delta, m_Delta0, m_EigenValues;
        unsigned int m_NumberOfVariables, m_SampleSize1, m_SampleSize2;
        bool m_TwoSampleTest, m_UseCorrection, m_ModifiedData, m_UsePairedTest;
        double m_Tolerance, m_TraceOfSigma, m_TraceOfSigmaSquare, m_StepSize;
    };
    
} // end namespace fdahotelling
