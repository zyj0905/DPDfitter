//Minimize class
#ifndef MINIMIZE_H 
#define MINIMIZE_H

#include "TObject.h"
#include "NLL_estimator.h"
#include "TMinuit.h"
#include "TMatrixD.h"

class Minimize: public TObject{
    public:
        Minimize(){
            num_NLL = 0;
            npar_tot = 0;
            cov_matrix = nullptr;
            par_fitted = nullptr;
            par_error = nullptr;
        }

        void Add_NLL(NLL_estimator* m_NLL_estimator){
            NLL_estimator_array[num_NLL] = m_NLL_estimator;
            num_NLL++;
            npar_tot = npar_tot +  m_NLL_estimator->N_totpar;
        }

        void Initial_optimize();
        void Do_SCAN();
        void Do_MIGRAD();
        void Do_HESSE();
        void SavePars(string filename, int my_idx_nll);
        void SaveCovMatrix();
        Json::Value WritePar(Para* mypar);
        void UpdateMypar();
        void GenerateAlterPars(int rndseed, double* genpars);
        void Cal_FitFraction(int my_idx_nll, double** Fit_fraction, double** Fit_fraction_std, TString save_file);
        void Print_FitFraction(int my_idx_nll, TString save_file);

        NLL_estimator* NLL_estimator_array[10];
        int num_NLL;
        TMinuit *gMinuit;
        double* cov_matrix;
        int npar_tot;
        double* par_fitted;
        double *par_error;
        
};
#endif //MINIMIZE_H