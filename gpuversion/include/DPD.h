//DPD class
#ifndef DPD_H 
#define DPD_H

#include "Event.h"
#include "DecayChain.h"
#include "DeviceComplex.h"

class DPD{
    public:
        __device__ DPD();
        __device__ double factorial(int x);
        __device__ double power_double_int(double a, int n){
            double result(1.0);
            for(int i=0;i<n;i++){result*=a;}
            return result;
        }
        __device__ double Wigner_smallD(double J, double _M, double _Mp, double beta);
        __device__ double Wigner_smallD_Jle4(int J, int _M, int _Mp, double beta);
        __device__ DeviceComplex Wigner_bigD(double j, double mp, double m, double alpha, double beta, double gamma);
        __device__ double CG_coeff_le224(double j1, double m1, double j2, double m2, double j3, double m3);
        __device__ double CG_coeff(double j1, double m1, double j2, double m2, double j3, double m3);
        __device__ double Blatt_Weisskopf_factor(double Q, int L);
        __device__ double Helicity_LScoupling(double p, int L);
        __device__ DeviceComplex Helicity_HCcoupling(double J, double j1, double lam1, double j2, double lam2, double mom, DecayChain* mychain, int type_ls);
        __device__ DeviceComplex Dalita_plot_function(double J, double s, double nu, double lam[3], double spin[3], Event* evt, int k, DecayChain* mychain);
};
#endif //DPD_H
