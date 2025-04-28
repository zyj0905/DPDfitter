#ifndef NLL_ESTIMATOR_H 
#define NLL_ESTIMATOR_H

#define CHECK_ERROR(call)                                                               \
    do                                                                                  \
    {                                                                                   \
        int _error = (call);                                                            \
        if (_error)                                                                     \
        {                                                                               \
            printf("*** Error *** at [%s:%d] error=%d \n", __FILE__, __LINE__, _error); \
        }                                                                               \
    } while (0)

#define CUDA_CHECK_ERROR(call)                                              \
    do                                                                      \
    {                                                                       \
        cudaError_t _error = (cudaError_t)(call);                           \
        if (_error != cudaSuccess)                                          \
        {                                                                   \
            printf("*** CUDA Error *** at [%s:%d] error=%d, reason:%s \n",  \
                   __FILE__, __LINE__, _error, cudaGetErrorString(_error)); \
        }                                                                   \
    } while (0)

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
            array_device_evt_dt = nullptr; array_device_evt_bg = nullptr; array_device_evt_mc = nullptr;
            with_sec = false; p4_dau1_name_sec = "void"; p4_dau2_name_sec = "void";
            amp_obj_device = nullptr;

            fit_type = 0; bg_ratio = 0.0;Normalization_factor_bg=0.0;
        }
        void set_sec_decay(TString p4_dau1_name, TString p4_dau2_name){
            with_sec = true;
            p4_dau1_name_sec = p4_dau1_name;
            p4_dau2_name_sec = p4_dau2_name;
        }
        void set_fit_type(int m_fit_type, TString m_PDF_bg_name, double m_bg_ratio){
            fit_type = m_fit_type;
            PDF_bg_name = m_PDF_bg_name;
            bg_ratio = m_bg_ratio;
        }
        void Update_ampobj_device(){
            CUDA_CHECK_ERROR(cudaFree(amp_obj_device));
            amp_obj_device = nullptr;
            CUDA_CHECK_ERROR(cudaMalloc( &amp_obj_device, sizeof(Amplitude)));
            CUDA_CHECK_ERROR(cudaMemcpy( amp_obj_device, amp_obj, sizeof(Amplitude), cudaMemcpyHostToDevice));
        }
        void InitDevice(){
            Update_ampobj_device();
            CUDA_CHECK_ERROR(cudaFree(array_device_evt_dt));
            CUDA_CHECK_ERROR(cudaFree(array_device_evt_bg));
            CUDA_CHECK_ERROR(cudaFree(array_device_evt_mc));
            array_device_evt_dt = nullptr; array_device_evt_bg = nullptr; array_device_evt_mc = nullptr;
            CUDA_CHECK_ERROR(cudaMalloc( &array_device_evt_dt, evt_dt*sizeof(Event)));
            CUDA_CHECK_ERROR(cudaMemcpy( array_device_evt_dt, array_evt_dt, evt_dt*sizeof(Event), cudaMemcpyHostToDevice));
            CUDA_CHECK_ERROR(cudaMalloc( &array_device_evt_bg, evt_bg*sizeof(Event)));
            CUDA_CHECK_ERROR(cudaMemcpy( array_device_evt_bg, array_evt_bg, evt_bg*sizeof(Event), cudaMemcpyHostToDevice));
            CUDA_CHECK_ERROR(cudaMalloc( &array_device_evt_mc, evt_mc*sizeof(Event)));
            CUDA_CHECK_ERROR(cudaMemcpy( array_device_evt_mc, array_evt_mc, evt_mc*sizeof(Event), cudaMemcpyHostToDevice));
        }
        void Load_file(int file_type, TString file_name, TString chain_name, TString p4_1_name, TString p4_2_name, TString p4_3_name);
        void CalPDF(Event* evt_arr, int num_evt, double* amp2);
        void CalPDFComponent(int idx_ch1, int idx_ch2, double *PDF);
        void GetFitFraction(double** Fit_Fraction);
        void Update_Paras();
        void Clear_LScoeff();
        double Cal_log_likelihood_nFit();
        double Cal_log_likelihood_cFit();
        double Cal_log_likelihood();
        
        LSCoeff_par* Search_LSCoeff_par(int type_ls, int m_idx_chain, int m_idx_LS, LSCoeff_par* m_LSCoeffpar_list);
        LSCoeff_par* Search_LSCoeff_par(int type_ls, int m_idx_chain, int m_idx_LS){return Search_LSCoeff_par(type_ls,m_idx_chain,m_idx_LS,LS1Coeffpar_list);};

        Reson_par* Search_Reson_par(int m_idx_chain, int m_idx_dynamic_par, Reson_par* m_Respar_list);
        Reson_par* Search_Reson_par(int m_idx_chain, int m_idx_dynamic_par){return Search_Reson_par(m_idx_chain,m_idx_dynamic_par,Respar_list);};
        void save_root(int file_type, TString file_name_out, TString file_name_in, TString chain_name, TString p4_1_name, TString p4_2_name, TString p4_3_name);
        Para* Search_idx_minuit(int m_idx_minuit, bool& isfind);
        Para* Search_par_name(string m_name, bool& isfind);
        

        Amplitude* amp_obj;
        Amplitude* amp_obj_device;
        LSCoeff_par* LS1Coeffpar_list; int N1_LSCoeff;
        LSCoeff_par* LS2Coeffpar_list; int N2_LSCoeff;
        Reson_par* Respar_list; int N_Respar;
        int N_totpar;
        
        int evt_dt;
        int evt_mc;
        int evt_mcT;
        int evt_bg;

        int fit_type;
        TString PDF_bg_name;
        double bg_ratio;
        double Normalization_factor_bg;

        Event* array_evt_dt; Event* array_device_evt_dt;
        Event* array_evt_mc; Event* array_device_evt_mc;
        Event* array_evt_bg; Event* array_device_evt_bg;

        bool with_sec;
        TString p4_dau1_name_sec;
        TString p4_dau2_name_sec;
        
        
        
};

#endif //NLL_ESTIMATOR_H
