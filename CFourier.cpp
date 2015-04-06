/* 
 * File:   CFourier.cpp
 * Author: lixun
 * 
 * Created on 2015年4月1日, 上午12:14
 */

#include "CFourier.h"
#include "CFastFourier.h"
#include "Config.h"





bool CFourier::IsRadix2(int Length) {
    if((Length & (Length-1)) == 0) {
        return true;
    } else {
        return false;
    }
}













void CFourier::dft_f(int Points, float * In, float * Out) {
    
    if(IsRadix2(Points)) {
        CFastFourier dft;
        dft.fft_f(Points, In, Out);
    } else {
        //Chirp-Z conversion...
    }
}

void CFourier::idft_f(int Points, float * In, float * Out) {
    
    if(IsRadix2(Points)) {
        CFastFourier dft;
        dft.ifft_f(Points, In, Out);
    } else {
        //Chirp-Z conversion...
    }    
}

void CFourier::dft_d(int Points, double * In, double * Out) {

    if(IsRadix2(Points)) {
        CFastFourier dft;
        dft.fft_d(Points, In, Out);
    } else {
        //Chirp-Z conversion...
    }    
}

void CFourier::idft_d(int Points, double * In, double * Out) {
    
    if(IsRadix2(Points)) {
        CFastFourier dft;
        dft.ifft_d(Points, In, Out);
    } else {
        //Chirp-Z conversion...
    }    
}





CFourier::CFourier() {
}

CFourier::CFourier(const CFourier& orig) {
}

CFourier::~CFourier() {
}

