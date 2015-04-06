/* 
 * File:   CCoefRadix4.cpp
 * Author: lixun
 * 
 * Created on 2015年3月10日, 下午7:53
 */



//  
//  radix 4 twiddle factor lookup table generation for qudrant 1, 2 and 3. 
//  qudrant 0 has coefficient to be all ones.
//
//       Wn(0x1) Wn(1x1) ... Wn((N/4-1)x1)
//       Wn(0x2) Wn(1x2) ... Wn((N/4-1)x2)
//       Wn(0x3) Wn(1x3) ... Wn((N/4-1)x3)
//
//  calculate twiddle factor lookup table for each
//  butter-fly pipeline of depth log2N
//
//        -i 2Pi k/N
//       e              k=[0..N/4-1]
//
//  data structure would be heap for each quadrant, for instance N=64
//  and quadrant 1. note that each cell is a complex number.
//
//  +===+===+===+===+===+===+===+===+===+===+===+===+
//  | 0 | 4 | 8 | 12| 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |--------o
//  +===+===+===+===+===+===+===+===+===+===+===+===+        |
//                                                           |
//  o--------------------------------------------------------o
//  |
//  |   +===+===+===+===+===+===+===+===+
//  o---| 8 | 9 | 10| 11| 12| 13| 14| 15|
//      +===+===+===+===+===+===+===+===+
//
//
//
//#define _USE_MATH_DEFINES
#include <cmath>
#include <new>
#include <stdlib.h>
#include <mm_malloc.h>
#include <iostream>

#include "instrset.h"
#include "vectorclass.h"
#include "complexvec.h"

#include "CCoef4.h"
#include "Config.h"




//#define DEBUG_CCOEFRADIX4




CCoef4::CCoef4(int Length) {
    
    mN = Length;

    // allocate cache-aligned memory for factor coefficient lookup table
    Wn = NULL;

    CoWnHeap1 = NULL;
    CoWnHeap2 = NULL;
    CoWnHeap3 = NULL;    

    CoWnHeap1I = NULL;
    CoWnHeap2I = NULL;
    CoWnHeap3I = NULL;    
    
    void * PlaceWn = NULL;

    void * PlaceCoWnHeap1 = NULL;
    void * PlaceCoWnHeap2 = NULL;
    void * PlaceCoWnHeap3 = NULL;    

    void * PlaceCoWnHeap1I = NULL;
    void * PlaceCoWnHeap2I = NULL;
    void * PlaceCoWnHeap3I = NULL;    
   
    //find out quadrant and heap length in term of complex elements
    mQuadPts = mN / 4;
    mHeapPts = 0;
    for(int i = mQuadPts; i >= 4; i = i / 4) {
        mHeapPts += i;
    }
    
    try{
        int RetWn  = posix_memalign(&PlaceWn,        CacheAlign, sizeof(double) * mQuadPts * Cx2Re);
        
        int RetCo1 = posix_memalign(&PlaceCoWnHeap1, CacheAlign, sizeof(float) * mHeapPts * Cx2Re);
        int RetCo2 = posix_memalign(&PlaceCoWnHeap2, CacheAlign, sizeof(float) * mHeapPts * Cx2Re);
        int RetCo3 = posix_memalign(&PlaceCoWnHeap3, CacheAlign, sizeof(float) * mHeapPts * Cx2Re);        
        
        int RetCo1I = posix_memalign(&PlaceCoWnHeap1I,CacheAlign, sizeof(float) * mHeapPts * Cx2Re);
        int RetCo2I = posix_memalign(&PlaceCoWnHeap2I,CacheAlign, sizeof(float) * mHeapPts * Cx2Re);
        int RetCo3I = posix_memalign(&PlaceCoWnHeap3I,CacheAlign, sizeof(float) * mHeapPts * Cx2Re);        
        
        if(RetWn || RetCo1 || RetCo2 || RetCo3 || RetCo1I || RetCo2I || RetCo3I) {
            throw RetWn + RetCo1 + RetCo2 + RetCo3 + RetCo1I + RetCo2I + RetCo3I;
        }

        // allocate the memory for twiddle factor and distribute them into
        // heap structure for pipeline-stage calculations. Note that there
        // is no need for pipeline zero stage in the FFT. See the comments
        // at the top of this file.

        Wn          = new(PlaceWn)          double[mQuadPts * Cx2Re] {};

        CoWnHeap1   = new(PlaceCoWnHeap1)   float[mHeapPts * Cx2Re] {};
        CoWnHeap2   = new(PlaceCoWnHeap2)   float[mHeapPts * Cx2Re] {};
        CoWnHeap3   = new(PlaceCoWnHeap3)   float[mHeapPts * Cx2Re] {};
        
        CoWnHeap1I  = new(PlaceCoWnHeap1I) float[mHeapPts * Cx2Re] {};                
        CoWnHeap2I  = new(PlaceCoWnHeap2I) float[mHeapPts * Cx2Re] {};                
        CoWnHeap3I  = new(PlaceCoWnHeap3I) float[mHeapPts * Cx2Re] {};                        
    }
    catch(int& RetComb) {
        free(PlaceWn);
        free(PlaceCoWnHeap1);
        free(PlaceCoWnHeap2);
        free(PlaceCoWnHeap3);        
        free(PlaceCoWnHeap1I);
        free(PlaceCoWnHeap2I);
        free(PlaceCoWnHeap3I);        
        std::cout << "_memalign failed: " << RetComb << std::endl;      
    }
        
}






void CCoef4::CoefGen(int Quadrant, float * CoWnHeap, double Math2) {
    
    // calculate the complex exponential function for Wn,
    // doing nothing if allocation fails.

    if(Wn != NULL) {
        for(int i = 0; i < mQuadPts * Cx2Re; i += 2) {
            Wn[i]   = 0.0;
            Wn[i+1] = (double)(i/2) * (double)Quadrant;
        }
        
        Vec4d coeff(Math2 * M_PI / (double)mN);
        Vec4d a;
        Vec4d b;
        for(int i = 0; i < mQuadPts * Cx2Re; i += 4) {
            a.load(Wn+i);
            b = a * coeff;
            b.store(Wn+i);
        }
        
        Complex4d ax;
        Complex4d bx;
        for(int i = 0; i < mQuadPts * Cx2Re; i += 4) {
            ax.load(Wn+i);
            bx = cexp(ax);
            bx.store(Wn+i);
        }
    }

    // map Wn to the heap, which could be all zeros if allocation fails.

#ifdef DEBUG_CCOEFRADIX4    
    for(int i = 0; i < mQuadPts * Cx2Re; i += 1) {
        Wn[i] = i;
    }
#endif
    
    if(CoWnHeap != NULL && Wn != NULL) {
        int WnStride   = mQuadPts/4;  //use N/4 instead of N/2 because the optimization:
        int WnElements = 4;           //the heap starts from the second conversion stage  
        int WnCnt = 0;                            
        
        while (WnStride > 0) {
            for(int i = 0; i < WnElements; i++) {
                CoWnHeap[WnCnt++] = (float)Wn[re(WnStride * i)];
                CoWnHeap[WnCnt++] = (float)Wn[im(WnStride * i)];            
            }
            WnElements = WnElements * 4;
            WnStride = WnStride / 4;
        }
    }
#ifdef DEBUG_CCOEFRADIX4
    for(int i = 0; i < mHeapPts * Cx2Re; i += 1) {
        std::cout << CoWnHeap[i] << std::endl;
    }
#endif
    //end of twiddle factor calculation and mapping  
}




void CCoef4::CoefComp(void) {
  
    CoefGen(1, CoWnHeap1, -2.0);
    CoefGen(2, CoWnHeap2, -2.0);
    CoefGen(3, CoWnHeap3, -2.0);
    
}


void CCoef4::CoefCompI(void) {
  
    CoefGen(1, CoWnHeap1I, 2.0);
    CoefGen(2, CoWnHeap2I, 2.0);
    CoefGen(3, CoWnHeap3I, 2.0);
}









// getters and getters


const float * CCoef4::GetCoefComp(int Quadrant) const {
    switch(Quadrant) {
        case 1: return (const float *) CoWnHeap1; 
        case 2: return (const float *) CoWnHeap2;
        case 3: return (const float *) CoWnHeap3;
        default: return NULL;        
    }

}

const float * CCoef4::GetCoefCompI(int Quadrant) const {
    switch(Quadrant) {
        case 1: return (const float *) CoWnHeap1I; 
        case 2: return (const float *) CoWnHeap2I;
        case 3: return (const float *) CoWnHeap3I;
        default: return NULL;                
    }
}


int CCoef4::GetPoints(void) const {
    return mN;
}









CCoef4::CCoef4(const CCoef4& orig) {
}


CCoef4::~CCoef4() {
    free(CoWnHeap3);
    free(CoWnHeap2);
    free(CoWnHeap1);
    free(Wn);
    free(CoWnHeap3I);
    free(CoWnHeap2I);
    free(CoWnHeap1I);
}

