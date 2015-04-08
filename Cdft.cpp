/* 
 * File:   CFourier.cpp
 * Author: lixun
 * 
 * Created on 2015年4月1日, 上午12:14
 */
#include "Config.h"

#include "CFastFourier.h"
#include "CFastBluestein.h"

#include "Cdft.h"







bool Cdft::IsRadix2(int Length) {
    if((Length & (Length-1)) == 0) {
        return true;
    } else {
        return false;
    }
}





void Cdft::dft_f(int Points, float * In, float * Out) {
    
    if(IsRadix2(Points)) {
        CFastFourier dft;
        dft.fft_f(Points, In, Out);
    } else {
        CFastBluestein dft;
        dft.chirpz_f(Points, In, Out);
    }
}

void Cdft::idft_f(int Points, float * In, float * Out) {
    
    if(IsRadix2(Points)) {
        CFastFourier dft;
        dft.ifft_f(Points, In, Out);
    } else {
        CFastBluestein dft;
        dft.ichirpz_f(Points, In, Out);
    }    
}

void Cdft::dft_d(int Points, double * In, double * Out) {

    if(IsRadix2(Points)) {
        CFastFourier dft;
        dft.fft_d(Points, In, Out);
    } else {
        CFastBluestein dft;
        dft.chirpz_d(Points, In, Out);
    }    
}

void Cdft::idft_d(int Points, double * In, double * Out) {
    
    if(IsRadix2(Points)) {
        CFastFourier dft;
        dft.ifft_d(Points, In, Out);
    } else {
        CFastBluestein dft;
        dft.ichirpz_d(Points, In, Out);
    }    
}





Cdft::Cdft() {
}

Cdft::Cdft(const Cdft& orig) {
}

Cdft::~Cdft() {
}

