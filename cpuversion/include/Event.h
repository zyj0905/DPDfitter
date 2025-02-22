//Event class
#ifndef EVENT_H 
#define EVENT_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
using namespace std;


class Event{
    public:
        void Init(){
            _alpha = 0; _beta = 0; _gamma = 0;
            _sigma1 = 0; _sigma2 = 0;
            _mass2_mom0 = 0; _mass2_dau1 = 0; _mass2_dau2 = 0; _mass2_dau3 = 0;
            _second_theta = 0; _second_phi = 0;
        }
        Event(){
            Init();
        }
        Event(const Event &other){
            Init();
            _mass2_mom0 = other._mass2_mom0;
            _mass2_dau1 = other._mass2_dau1; _mass2_dau2 = other._mass2_dau2; _mass2_dau3 = other._mass2_dau3;
            _alpha = other._alpha; _beta = other._beta; _gamma = other._gamma;
            _second_theta = other._second_theta; _second_phi = other._second_phi;
            _sigma1 = other._sigma1; _sigma2 = other._sigma2;
        }
        Event& operator=(const Event &other)
        {
            if(this == &other)
                return *this;
            Init();
            _mass2_mom0 = other._mass2_mom0;
            _mass2_dau1 = other._mass2_dau1; _mass2_dau2 = other._mass2_dau2; _mass2_dau3 = other._mass2_dau3;
            _alpha = other._alpha; _beta = other._beta; _gamma = other._gamma;
            _second_theta = other._second_theta; _second_phi = other._second_phi;
            _sigma1 = other._sigma1; _sigma2 = other._sigma2;
            return *this;
        }
        
        double _alpha; double _beta; double _gamma;
        double _sigma1; double _sigma2;
        double _mass2_mom0; double _mass2_dau1; double _mass2_dau2; double _mass2_dau3;  
        double _second_theta; double _second_phi;

        double sigma3_func();
        double sqrt_lamb_func(double x, double y, double z);
        double alignment_angle_func(int idx);
        double scatter_angle_func(int idx);
        double wrotation_angle_func(int idx1, int idx2);
        double break_mom1_func(int idx);
        double break_mom2_func(int idx);
};

#endif //EVENT_H
