#include "../include/Event.h"

/////////////////
//DEVICE PART
////////////////

double Event::sigma3_func(){
    double _mass2_sum = _mass2_mom0 + _mass2_dau1 + _mass2_dau2 + _mass2_dau3;
    return _mass2_sum - _sigma1 - _sigma2;
}

double Event::sqrt_lamb_func(double x, double y, double z){
    return sqrt((x*x+y*y+z*z)-2.0*(x*y+x*z+y*z));
}

double Event::break_mom1_func(int idx){
    if(idx==1){
        return sqrt_lamb_func(_mass2_mom0,_mass2_dau1,_sigma1)/2.0/sqrt(_mass2_mom0);
    }
    if(idx==2){
        return sqrt_lamb_func(_mass2_mom0,_mass2_dau2,_sigma2)/2.0/sqrt(_mass2_mom0);
    }
    if(idx==3){
        double _sigma3 = sigma3_func();
        return sqrt_lamb_func(_mass2_mom0,_mass2_dau3,_sigma3)/2.0/sqrt(_mass2_mom0);
    }
    printf("break_mom1_func: OUT OF INDEX!");
    return 0.0;
    
}

double Event::break_mom2_func(int idx){
    if(idx==1){
        return sqrt_lamb_func(_sigma1,_mass2_dau2,_mass2_dau3)/2/sqrt(_sigma1);
    }
    if(idx==2){
        return sqrt_lamb_func(_sigma2,_mass2_dau1,_mass2_dau3)/2/sqrt(_sigma2);
    }
    if(idx==3){
        double _sigma3 = sigma3_func();
        return sqrt_lamb_func(_sigma3,_mass2_dau1,_mass2_dau2)/2/sqrt(_sigma3);
    }
    printf("break_mom2_func: OUT OF INDEX!");
    return 0.0;
    
}

double Event::alignment_angle_func(int idx){
    if(idx==1){
        return 0.0;
    }
    if(idx==2){
        double _sigma3 = sigma3_func();
        double dom = sqrt_lamb_func(_mass2_mom0,_mass2_dau2,_sigma2)*sqrt_lamb_func(_mass2_mom0,_sigma1,_mass2_dau1);
        double num = (_mass2_mom0+_mass2_dau1-_sigma1)*(_mass2_mom0+_mass2_dau2-_sigma2)-2.0*_mass2_mom0*(_sigma3-_mass2_dau1-_mass2_dau2);
        return -acos(1.0*(num/dom));
    }
    if(idx==3){
        double _sigma3 = sigma3_func();
        double dom = sqrt_lamb_func(_mass2_mom0,_mass2_dau1,_sigma1)*sqrt_lamb_func(_mass2_mom0,_sigma3,_mass2_dau3);
        double num = (_mass2_mom0+_mass2_dau3-_sigma3)*(_mass2_mom0+_mass2_dau1-_sigma1)-2.0*_mass2_mom0*(_sigma2-_mass2_dau3-_mass2_dau1);
        return acos((num/dom));
    }
    printf("alignment_angle_func: OUT OF INDEX!");
    return 0.0;
}

double Event::scatter_angle_func(int idx){
    double _sigma3 = sigma3_func();
    if(idx==1){
        double dom = sqrt_lamb_func(_mass2_mom0,_mass2_dau1,_sigma1)*sqrt_lamb_func(_sigma1,_mass2_dau2,_mass2_dau3);
        double num = 2.0*_sigma1*(_sigma3-_mass2_dau1-_mass2_dau2)-(_sigma1+_mass2_dau2-_mass2_dau3)*(_mass2_mom0-_sigma1-_mass2_dau1);
        return acos(num/dom);
    }
    if(idx==2){
        double dom = sqrt_lamb_func(_mass2_mom0,_mass2_dau2,_sigma2)*sqrt_lamb_func(_sigma2,_mass2_dau3,_mass2_dau1);
        double num = 2.0*_sigma2*(_sigma1-_mass2_dau2-_mass2_dau3)-(_sigma2+_mass2_dau3-_mass2_dau1)*(_mass2_mom0-_sigma2-_mass2_dau2);
        return acos(num/dom);
    }
    if(idx==3){
        double dom = sqrt_lamb_func(_mass2_mom0,_mass2_dau3,_sigma3)*sqrt_lamb_func(_sigma3,_mass2_dau1,_mass2_dau2);
        double num = 2.0*_sigma3*(_sigma2-_mass2_dau3-_mass2_dau1)-(_sigma3+_mass2_dau1-_mass2_dau2)*(_mass2_mom0-_sigma3-_mass2_dau3);
        return acos(num/dom);
    }
    printf("scatter_angle_func: OUT OF INDEX!");
    return 0.0;
}

double Event::wrotation_angle_func(int idx1, int idx2){
    double _sigma3 = sigma3_func();
    if(idx1==1){
        if(idx2==1){
            return 0.0;
        }
        if(idx2==2){
            double dom = sqrt_lamb_func(_mass2_mom0,_mass2_dau1,_sigma1)*sqrt_lamb_func(_sigma2,_mass2_dau1,_mass2_dau3);
            double num = 2*_mass2_dau1*(_sigma3-_mass2_mom0-_mass2_dau3) + (_mass2_mom0+_mass2_dau1-_sigma1)*(_sigma2-_mass2_dau1-_mass2_dau3);
            return acos(num/dom);
        }
        if(idx2==3){
            double dom = sqrt_lamb_func(_mass2_mom0,_mass2_dau1,_sigma1)*sqrt_lamb_func(_sigma3,_mass2_dau1,_mass2_dau2);
            double num = 2*_mass2_dau1*(_sigma2-_mass2_mom0-_mass2_dau2) + (_mass2_mom0+_mass2_dau1-_sigma1)*(_sigma3-_mass2_dau1-_mass2_dau2);
            return -acos(1.0*(num/dom));
        }
    }
    if(idx1==2){
        if(idx2==1){
            double dom = sqrt_lamb_func(_mass2_mom0,_mass2_dau2,_sigma2)*sqrt_lamb_func(_sigma1,_mass2_dau2,_mass2_dau3);
            double num = 2*_mass2_dau2*(_sigma3-_mass2_mom0-_mass2_dau3) + (_mass2_mom0+_mass2_dau2-_sigma2)*(_sigma1-_mass2_dau2-_mass2_dau3);
            return -acos(1.0*(num/dom));
        }
        if(idx2==2){
            return 0.0;
        }
        if(idx2==3){
            double dom = sqrt_lamb_func(_mass2_mom0,_mass2_dau2,_sigma2)*sqrt_lamb_func(_sigma3,_mass2_dau2,_mass2_dau1);
            double num = 2*_mass2_dau2*(_sigma1-_mass2_mom0-_mass2_dau1) + (_mass2_mom0+_mass2_dau2-_sigma2)*(_sigma3-_mass2_dau2-_mass2_dau1);
            return acos(num/dom);
        }
    }
    if(idx1==3){
        if(idx2==1){
            double dom = sqrt_lamb_func(_mass2_mom0,_mass2_dau3,_sigma3)*sqrt_lamb_func(_sigma1,_mass2_dau3,_mass2_dau2);
            double num = 2*_mass2_dau3*(_sigma2-_mass2_mom0-_mass2_dau2) + (_mass2_mom0+_mass2_dau3-_sigma3)*(_sigma1-_mass2_dau3-_mass2_dau2);
            return acos(num/dom);
        }
        if(idx2==2){
            double dom = sqrt_lamb_func(_mass2_mom0,_mass2_dau3,_sigma3)*sqrt_lamb_func(_sigma2,_mass2_dau3,_mass2_dau1);
            double num = 2*_mass2_dau3*(_sigma1-_mass2_mom0-_mass2_dau1) + (_mass2_mom0+_mass2_dau3-_sigma3)*(_sigma2-_mass2_dau3-_mass2_dau1);
            return -acos(1.0*(num/dom));
            
        }
        if(idx2==3){
            return 0.0;
        }
    }

    printf("wrotation_angle_func: OUT OF INDEX!");
    return 0.0;
}