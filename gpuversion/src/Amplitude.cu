//The overall amplitude based on DPD and sum over the decay chains
#include "../include/Amplitude.h"

__device__ DeviceComplex Amplitude::SumDecayChain(Event* evt, double nu, double lam[3]){
    DPD DPD_obj(myQ0);
    DeviceComplex DPD_sum(0.0,0.0);
    for(int idx_chain=0;idx_chain<nchain;idx_chain++){
        double spin_isobar = array_chain[idx_chain].intermediate.spin;
        int idx_isobar = array_chain[idx_chain].idx_isobar;
        double sqrt_s(0.0);
        if(idx_isobar==1){sqrt_s = sqrt(evt->_sigma1);}
        if(idx_isobar==2){sqrt_s = sqrt(evt->_sigma2);}
        if(idx_isobar==3){sqrt_s = sqrt(evt->sigma3_func());}

        DPD_sum = DPD_sum + DPD_obj.Dalita_plot_function(spin_mom,spin_isobar,nu,lam,spin_dau,evt,idx_isobar,&(array_chain[idx_chain]))
        *(array_chain[idx_chain].intermediate.dynamic.eval(sqrt_s,evt,idx_isobar));
    }

    return DPD_sum;

}

__device__ DeviceComplex Amplitude::SumOverNu(Event* evt, double Lamb, double lam[3]){
    DPD DPD_obj(myQ0);
    DeviceComplex Nu_sum(0.0);
    for(double nu=-spin_mom;nu<=spin_mom;nu++){
        DeviceComplex wignerD = (DPD_obj.Wigner_bigD(spin_mom,Lamb,nu,evt->_alpha,evt->_beta,evt->_gamma)).conjugate();
        Nu_sum = Nu_sum + wignerD*SumDecayChain(evt,nu,lam);
    }
    return Nu_sum;
}

__device__ DeviceComplex Amplitude::Amp_Secondary_Decay(Event* evt, double lam, int type){
    DPD DPD_obj(myQ0);
    //Vpp vertex
    if(type==0){
        return (DPD_obj.Wigner_bigD(1,lam,0,evt->_second_phi,evt->_second_theta,0)).conjugate();
    }
    return 1;
}

__device__ DeviceComplex Amplitude::SumSecondary(Event* evt, double Lamb, double lam[3]){

    if(idx_sec==-1){return SumOverNu(evt,Lamb,lam);}

    DeviceComplex lam_sec_sum(0.0);
    for(double lam_sec=-spin_dau[idx_sec-1];lam_sec<=spin_dau[idx_sec-1];lam_sec++){
        lam[idx_sec-1] = lam_sec;
        DeviceComplex sec_decay(1,0);
        sec_decay = Amp_Secondary_Decay(evt,lam_sec,type_sec);
        lam_sec_sum = lam_sec_sum + SumOverNu(evt,Lamb,lam)*sec_decay;
    }
    return lam_sec_sum;
}

__device__ DeviceComplex Amplitude::SumOverlam(Event* evt, double Lamb, double Lambp){
    DeviceComplex lam_sum(0.0);
    double lam[3] = {0.0};

    for(double lam1=-spin_dau[0];lam1<=spin_dau[0];lam1++){
        for(double lam2=-spin_dau[1];lam2<=spin_dau[1];lam2++){
            for(double lam3=-spin_dau[2];lam3<=spin_dau[2];lam3++){
                lam[0] = lam1; lam[1] = lam2; lam[2] = lam3;
                if(idx_sec==1 && lam1!=spin_dau[0]){continue;}
                if(idx_sec==2 && lam2!=spin_dau[1]){continue;}
                if(idx_sec==3 && lam3!=spin_dau[2]){continue;}

                if(Lambp==Lamb){lam_sum = lam_sum + ((SumSecondary(evt,Lamb,lam)).rho2());}
                else{lam_sum = lam_sum + (SumSecondary(evt,Lamb,lam))*((SumSecondary(evt,Lambp,lam)).conjugate());}
            }
        }
    }
    return lam_sum;
}


__device__ double Amplitude::SumOverLam(Event* evt){
    DeviceComplex Lam_sum(0.0);
    for(double Lam=-spin_mom;Lam<=spin_mom;Lam++){
        for(double Lamp=-spin_mom;Lamp<=spin_mom;Lamp++){
            int idx_Lam = int(Lam+spin_mom);
            int idx_Lamp = int(Lamp+spin_mom);
            int idx_matrix = int(idx_Lam*(2*spin_mom+1)+idx_Lamp);
            if(Spin_Density_Matrix[idx_matrix]==0){continue;}
            Lam_sum = Lam_sum + SumOverlam(evt,Lam,Lamp)*Spin_Density_Matrix[idx_matrix];
        }
    }
    return Lam_sum.real;
}
