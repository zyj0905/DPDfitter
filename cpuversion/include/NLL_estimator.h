#ifndef NLL_ESTIMATOR_H 
#define NLL_ESTIMATOR_H

//====ROOT LIB======
#include "TLorentzVector.h"
#include "TChain.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TObject.h"
//====ROOT LIB======

#include <fstream>
#include <iostream>
#include <string>
#include "json.h"
#include "Amplitude.h"
#include "Parameter.h"
using namespace std;

class NLL_estimator: public TObject
{
    public:
        NLL_estimator(Amplitude* m_amp_obj, LSCoeff_par* m_LS1Coeffpar_list, LSCoeff_par* m_LS2Coeffpar_list, Reson_par* m_Respar_list){
            amp_obj = m_amp_obj; Respar_list = m_Respar_list; N_Respar = m_amp_obj->Get_total_dynamic_par();
            LS1Coeffpar_list = m_LS1Coeffpar_list; LS2Coeffpar_list = m_LS2Coeffpar_list;
            N1_LSCoeff = m_amp_obj->Get_total_LS1coeff_par(); N2_LSCoeff = m_amp_obj->Get_total_LS2coeff_par();
            N_totpar = 2*(N1_LSCoeff + N2_LSCoeff) + N_Respar;
            
            evt_dt = 0; evt_bg = 0; evt_mc = 0;
            array_evt_dt = nullptr; array_evt_bg = nullptr; array_evt_mc = nullptr;
            with_sec = false; p4_dau1_name_sec = "void"; p4_dau2_name_sec = "void";
        }
        void set_sec_decay(TString p4_dau1_name, TString p4_dau2_name){
            with_sec = true;
            p4_dau1_name_sec = p4_dau1_name;
            p4_dau2_name_sec = p4_dau2_name;
        }

        void Load_file(int file_type, TString file_name, TString chain_name, TString p4_1_name, TString p4_2_name, TString p4_3_name);
        void CalPDF(Event* evt_arr, int num_evt, double* amp2);
        void Update_Paras();
        void Clear_LScoeff();
        
        LSCoeff_par* Search_LSCoeff_par(int type_ls, int m_idx_chain, int m_idx_LS, LSCoeff_par* m_LSCoeffpar_list);
        LSCoeff_par* Search_LSCoeff_par(int type_ls, int m_idx_chain, int m_idx_LS){return Search_LSCoeff_par(type_ls,m_idx_chain,m_idx_LS,LS1Coeffpar_list);};

        Reson_par* Search_Reson_par(int m_idx_chain, int m_idx_dynamic_par, Reson_par* m_Respar_list);
        Reson_par* Search_Reson_par(int m_idx_chain, int m_idx_dynamic_par){return Search_Reson_par(m_idx_chain,m_idx_dynamic_par,Respar_list);};
        void save_root(int file_type, TString file_name, int save_component);
        Para* Search_idx_minuit(int m_idx_minuit, bool& isfind);
        Para* Search_par_name(string m_name, bool& isfind);
        

        Amplitude* amp_obj;
        LSCoeff_par* LS1Coeffpar_list; int N1_LSCoeff;
        LSCoeff_par* LS2Coeffpar_list; int N2_LSCoeff;
        Reson_par* Respar_list; int N_Respar;
        int N_totpar;
        
        int evt_dt;
        int evt_mc;
        int evt_mcT;
        int evt_bg;

        Event* array_evt_dt;
        Event* array_evt_mc;
        Event* array_evt_mcT;
        Event* array_evt_bg;

        bool with_sec;
        TString p4_dau1_name_sec;
        TString p4_dau2_name_sec;
        
        
        
};

#endif //NLL_ESTIMATOR_H
