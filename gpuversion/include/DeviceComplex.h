//======================================================================
//The complex number class used in GPU
//======================================================================

#ifndef DEVICECOMPLEX_H 
#define DEVICECOMPLEX_H

#include <cmath>
#include <stdio.h>
#include <iostream>

class DeviceComplex {
public:
    double real;
    double imag;
     __device__ DeviceComplex(double r = 0, double i = 0) : real(r), imag(i) {}
    // Simple operators complex - complex
     __device__ DeviceComplex operator+(const DeviceComplex& other) const {
        return DeviceComplex(real + other.real, imag + other.imag);
    }
     __device__ DeviceComplex operator-(const DeviceComplex& other) const {
        return DeviceComplex(real - other.real, imag - other.imag);
    }
     __device__ DeviceComplex operator*(const DeviceComplex& other) const {
        return DeviceComplex(real * other.real - imag * other.imag, real * other.imag + imag * other.real);
    }
     __device__ DeviceComplex operator/(const DeviceComplex& other) const {
        double denom = other.real * other.real + other.imag * other.imag;
        return DeviceComplex((real * other.real + imag * other.imag) / denom, (imag * other.real - real * other.imag) / denom);
    }

    // Simple operators complex - double
     __device__ DeviceComplex operator *(double c) const {
        return DeviceComplex(real * c, imag * c);
    }
     __device__ DeviceComplex operator +(double c) const {
        return DeviceComplex(real + c, imag);
    }
     __device__ DeviceComplex operator /(double c) const {
        return DeviceComplex(real / c, imag / c);
    }
     __device__ DeviceComplex operator -(double c) const {
        return DeviceComplex(real - c, imag);
    }

    // Simple operators double - complex
     __device__ friend DeviceComplex operator *(double c, const DeviceComplex & other) {
        return DeviceComplex(c * other.real, c * other.imag);
    }
     __device__ friend DeviceComplex operator +(double c, const DeviceComplex & other) {
        return DeviceComplex(c + other.real, other.imag);
    }
     __device__ friend DeviceComplex operator /(double c, const DeviceComplex & other) {
        double denom = other.real * other.real + other.imag * other.imag;
        return DeviceComplex(c * other.real, -c * other.imag) / denom;
    }
     __device__ friend DeviceComplex operator -(double c, const DeviceComplex & other) {
        return DeviceComplex(c - other.real, -other.imag);
    }

     __device__ double rho() const {
        return std::sqrt(real * real + imag * imag);
    }
     __device__ double rho2() const {
        return real * real + imag * imag;
    }
     __device__ double phi() const {
        double theta = acos(real/std::sqrt(real * real + imag * imag));
        if(imag>=0){return theta;}
        else{return -1.0*theta;}
     }

     __device__ DeviceComplex reciprocal() const {
        double denom = real * real + imag * imag;
        return DeviceComplex(real / denom, -imag / denom);
    }
     __device__ DeviceComplex conjugate() const {
        return DeviceComplex(real, -imag);
    }

    __device__ DeviceComplex sqrt() const {
        double myrho = std::sqrt(rho());
        double myphi = phi() * 0.5;
        return DeviceComplex(myrho*cos(myphi),myrho*sin(myphi));
    }

    __device__ DeviceComplex ln() const {
        double myrho = rho();
        double myphi = phi();
        return log(myrho) + DeviceComplex(0,1) * myphi;
    }

    __device__ DeviceComplex atan() const {
        return 1.0/(DeviceComplex(0,2)) * ((DeviceComplex(0,1)-DeviceComplex(real,imag))/(DeviceComplex(0,1)+DeviceComplex(real,imag))).ln();
    }

     __device__ void print() const {
        printf("%f + %fi\n", real, imag);
    }
     __device__ DeviceComplex& operator=(const DeviceComplex& other) {
        if (this != &other) {
            real = other.real;
            imag = other.imag;
        }
        return *this;
    }
    
};

#endif //DEVICECOMPLEX_H