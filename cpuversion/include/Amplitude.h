//Amplitude class
#ifndef AMPLITUDE_H 
#define AMPLITUDE_H

#include "DecayChain.h"
#include "DPD.h"

class Amplitude{
    public:
        double spin_mom; double spin_dau[3];
        double Spin_Density_Matrix[100];
        int nchain;
        DecayChain array_chain[20];
        //Secondary decays
        int idx_sec; int type_sec;
        
        Amplitude(particle particle_list[4]){
            spin_mom = particle_list[0].spin; spin_dau[0] = particle_list[1].spin; spin_dau[1] = particle_list[2].spin; spin_dau[2] = particle_list[3].spin;
            nchain = 0;
            idx_sec = -1; type_sec = -1;
        }
        void AddDecayChain(DecayChain chain){array_chain[nchain] = chain;nchain++;}
        void SetSDM(double* SDM){for(int i=0;i<pow(int(spin_mom*2+1),2);i++){ Spin_Density_Matrix[i] = SDM[i];}}
        void Add_Secondary_Decay(int m_type_sec, int m_idx_sec){type_sec = m_type_sec;idx_sec = m_idx_sec;}

        int Get_total_LS1coeff_par(){
            int sum(0);
            for(int idx_chain=0;idx_chain<nchain;idx_chain++){sum = sum + array_chain[idx_chain].N1_LS;}
            return sum;
        }

        int Get_total_LS2coeff_par(){
            int sum(0);
            for(int idx_chain=0;idx_chain<nchain;idx_chain++){sum = sum + array_chain[idx_chain].N2_LS;}
            return sum;
        }

        int Get_total_dynamic_par(){
            int sum(0);
            for(int idx_chain=0;idx_chain<nchain;idx_chain++){
                sum = sum + array_chain[idx_chain].intermediate.dynamic.num_par;
            }
            return sum;
        }
        
        DeviceComplex Amp_Secondary_Decay(Event* evt, double Lamb, int type);
        DeviceComplex SumDecayChain(Event* evt, double nu, double lam[3]);
        DeviceComplex SumOverNu(Event* evt, double Lamb, double lam[3]);
        DeviceComplex SumOverlam(Event* evt, double Lamb, double Lambp);
        double SumOverLam(Event* evt);

};
#endif //AMPLITUDE_H