//Parameter class
#ifndef PARAMETER_H 
#define PARAMETER_H
#include <string>
#include <iostream>
using namespace std;

class Para{
    public:
        double Val;
        double Step;
        double up_bound;
        double low_bound;
        string para_name;
        int idx_minuit;

        bool is_fixed_to;
        Para* fixed_to_this;
        double fix_ratio;

    Para(){
        Val = 0; up_bound = 0; low_bound = 0; Step = 0; idx_minuit = -1;
        para_name = "void"; is_fixed_to = false; fixed_to_this = nullptr; fix_ratio = 0.0;
    }

    Para(string m_name, double m_Val, double m_Step, double m_low_bound, double m_up_bound){
        para_name = m_name;
        Val = m_Val;
        Step = m_Step;
        low_bound = m_low_bound;
        up_bound = m_up_bound;
        is_fixed_to = false; fixed_to_this = nullptr; idx_minuit = -1;
        fix_ratio = 0.0;
    }

    Para &operator=(const Para& other){
        if(this == &other) return *this;
        para_name = other.para_name;
        Val = other.Val;
        Step = other.Step;
        idx_minuit = other.idx_minuit;
        low_bound = other.low_bound;
        up_bound = other.up_bound;
        is_fixed_to = other.is_fixed_to;
        fixed_to_this = other.fixed_to_this;
        fix_ratio = other.fix_ratio;
        
        return *this;
    }

    void fixed_to(Para* other_cpar, double ratio){
        fixed_to_this = other_cpar;
        is_fixed_to = true;
        Step = 0;
        fix_ratio = ratio;
    }

    double Get_Val(){
        if(is_fixed_to==false){return Val;}
        else{return fix_ratio*fixed_to_this->Get_Val();}
    }

    void Save_Fix(){if(is_fixed_to==true){Val = Get_Val();}}

};


class Reson_par{
    public:
        Para par;
        int idx_chain;
        int idx_dynamic_par;

        Reson_par(){idx_chain = -1; idx_dynamic_par = -1;}
        Reson_par(Para m_par, int m_idx_chain, int m_idx_dynamic_par){
            par = m_par;
            idx_chain = m_idx_chain;
            idx_dynamic_par = m_idx_dynamic_par;
        }

        Reson_par &operator=(const Reson_par& other){
            if(this == &other) return *this;
            par = other.par;
            idx_chain = other.idx_chain;
            idx_dynamic_par = other.idx_dynamic_par;
            return *this;
        }
};

class LSCoeff_par{
    public:
        Para par_rho;
        Para par_phi;
        int idx_chain;
        int idx_LS;

        LSCoeff_par(){idx_chain = -1; idx_LS = -1;}
        LSCoeff_par(Para m_par_rho, Para m_par_phi, int m_idx_chain, int m_idx_LS){
            par_rho = m_par_rho; par_phi = m_par_phi;
            idx_chain = m_idx_chain;
            idx_LS = m_idx_LS;
        }

        LSCoeff_par &operator=(const LSCoeff_par& other){
            if(this == &other) return *this;
            par_rho = other.par_rho;
            par_phi = other.par_phi;
            idx_chain = other.idx_chain;
            idx_LS = other.idx_LS;
            return *this;
        }
};

#endif //PARAMETER_H