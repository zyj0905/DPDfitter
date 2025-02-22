//Minimize class
#ifndef MINIMIZE_H 
#define MINIMIZE_H

#include "TObject.h"
#include "NLL_estimator.h"
#include "TMinuit.h"

class Minimize: public TObject{
    public:
        Minimize(){
            num_NLL = 0;
        }

        void Add_NLL(NLL_estimator* m_NLL_estimator){
            NLL_estimator_array[num_NLL] = m_NLL_estimator;
            num_NLL++;
        }

        void Initial_optimize();
        void Do_SCAN();
        void Do_MIGRAD();
        void Do_HESSE();
        void SavePars(string filename, int my_idx_nll);
        Json::Value WritePar(Para* mypar);
        void UpdateMypar();

        NLL_estimator* NLL_estimator_array[10];
        int num_NLL;
        TMinuit *gMinuit;
        
};
#endif //MINIMIZE_H