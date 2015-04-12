/* 
 * File:   CChirpZd.cpp
 * Author: lixun
 * 
 * Created on 2015年4月5日, 下午6:40
 */

#include <cstring>
#include <memory>
#include <cmath>
#include <new>
#include <cstdlib>
#include <mm_malloc.h>
#include <iostream>

#include "instrset.h"
#include "vectorclass.h"
#include "complexvec.h"

#include "CChirpZd.h"
#include "Config.h"
#include "CFastFourier.h"

#include "fftw3.h"



CChirpZd::CChirpZd(int Points) {
    
    mPoints = Points;
    mPointsChirpZ = 1;
    while(mPointsChirpZ < (2 * mPoints - 1)) {
        mPointsChirpZ <<= 1;
    }
    
    C = NULL;
    X = NULL;
    Cs = NULL;
    Xs = NULL;
    
    void * PlaceC = NULL;
    void * PlaceX = NULL;
    void * PlaceCs = NULL;
    void * PlaceXs = NULL;
    
    try {
        int RetC = posix_memalign(&PlaceC, CacheAlign, sizeof(double) * mPointsChirpZ * Cx2Re);
        int RetX = posix_memalign(&PlaceX, CacheAlign, sizeof(double) * mPointsChirpZ * Cx2Re);        
        int RetCs = posix_memalign(&PlaceCs, CacheAlign, sizeof(double) * mPointsChirpZ * Cx2Re);
        int RetXs = posix_memalign(&PlaceXs, CacheAlign, sizeof(double) * mPointsChirpZ * Cx2Re);        
        
        if(RetC || RetX || RetCs || RetXs) {
            throw RetC + RetX + RetCs + RetXs ;
        }
        
        C = new(PlaceC) double[mPointsChirpZ * Cx2Re] {};
        X = new(PlaceX) double[mPointsChirpZ * Cx2Re] {};
        Cs = new(PlaceCs) double[mPointsChirpZ * Cx2Re] {};
        Xs = new(PlaceXs) double[mPointsChirpZ * Cx2Re] {};        
        
    } catch(int& RetComb) {
        free(PlaceC);
        free(PlaceX);
        free(PlaceCs);
        free(PlaceXs);        
        std::cout << "ChirpZ _memalign failed: " << RetComb << std::endl;      
    }
}








CChirpZd::CChirpZd(const CChirpZd& orig) {
}


CChirpZd::~CChirpZd() {
    free(C);
    free(X);
    free(Cs);
    free(Xs);      
}


double * CChirpZd::GetInputDouble() {
    return X;
}

double * CChirpZd::GetOutputDouble() {
    return X;
}






//      -i 2pi/N
// W = e
//
//                2             2
//           1/2 n     -i pi/N n
// C[n] =   W       = e
//      
// n = 0, 1, ..., N-1 where N is the number of points.
//
// modulate X[n] and then zero padding to M, which is radix 2 and larger than or equal to 2N-1
// Y[n] = X[n] * C[n], then Y[n] is written back to X[n].  

void CChirpZd::ConvertDouble(long double MathPiWithSign) {
    
    // assume input data have been stored in X from 0 to mPoints-1.
    // pad the rest with zeros
    memset(X + (mPoints * Cx2Re), 0, sizeof(double) * Cx2Re * (mPointsChirpZ - mPoints));

    
    // generate C[n]
    long long int n  = 0;
    long long int n2 = 0;
    long long int mPoints2 = 2 * mPoints;
    long double Index = 0.0L;
    
    for(int i = 0; i < mPoints; i++) {
        
        n2 = n % mPoints2;
        Index = (long double)n2 / (long double)mPoints * MathPiWithSign;
        n += (2*i+1);
        
        //std::cout << ":" << n << std::endl;        
        C[re(i)] = 0.0;
        C[im(i)] = (double)Index;
    }

    Complex4d ax;
    Complex4d bx;
    Complex4d cx;
    Complex4d dx;
    for(int i = 0; i < mPoints*Cx2Re; i+=4) {
        cx.load(X+i);
        ax.load(C+i);
        bx = cexp(-ax);  // for modulation of X
        dx = bx * cx;
        dx.store(X+i);
        
        bx = cexp(ax);   // for de-modulation of Y and convolution
        bx.store(C+i);
    }
    
    // so far X has been modulated and zero padded
    // coefficients for de-modulation of Y is ready
    // start to generate h[n] for convolution based on its even symmetry.
    for(int i = 1; i < mPoints; i++) {
        C[re(mPointsChirpZ-i)] = C[re(i)];
        C[im(mPointsChirpZ-i)] = C[im(i)];        
    }

    std::unique_ptr<CFastFourier> pFFTL(new CFastFourier);
    pFFTL->fft_d(mPointsChirpZ, X, Xs);
    pFFTL->fft_d(mPointsChirpZ, C, Cs);

    
    for(int i = 0; i < mPointsChirpZ*Cx2Re; i+=4) {
        ax.load(Xs+i);
        bx.load(Cs+i);
        cx = ax * bx;
        cx.store(Xs+i);
    }  
    
    pFFTL->ifft_d(mPointsChirpZ, Xs, X);
    
    Complex4d GainIFFTx((double)mPointsChirpZ);
    for(int i = 0; i < mPoints*Cx2Re; i+=4) {
        ax.load(X+i);
        bx.load(C+i);
        cx = ax / bx;
        dx = cx / GainIFFTx;
        dx.store(X+i);
    }    
}










