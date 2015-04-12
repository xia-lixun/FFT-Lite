/* 
 * File:   CCoef4d.cpp
 * Author: lixun
 * 
 * Created on 2015年3月16日, 上午12:21
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
#include "CCoef4d.h"
#include "Config.h"






CCoef4d::CCoef4d(){
}




CCoef4d::CCoef4d(int Length) : CCoef4() {
    
    mN = Length;

    // allocate cache-aligned memory for factor coefficient lookup table
    CoWnHeap1 = NULL;
    CoWnHeap2 = NULL;
    CoWnHeap3 = NULL;    

    CoWnHeap1I = NULL;
    CoWnHeap2I = NULL;
    CoWnHeap3I = NULL;    

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
        int RetCo1 = posix_memalign(&PlaceCoWnHeap1, CacheAlign, sizeof(double) * mHeapPts * Cx2Re);
        int RetCo2 = posix_memalign(&PlaceCoWnHeap2, CacheAlign, sizeof(double) * mHeapPts * Cx2Re);
        int RetCo3 = posix_memalign(&PlaceCoWnHeap3, CacheAlign, sizeof(double) * mHeapPts * Cx2Re);        
        
        int RetCo1I = posix_memalign(&PlaceCoWnHeap1I,CacheAlign, sizeof(double) * mHeapPts * Cx2Re);
        int RetCo2I = posix_memalign(&PlaceCoWnHeap2I,CacheAlign, sizeof(double) * mHeapPts * Cx2Re);
        int RetCo3I = posix_memalign(&PlaceCoWnHeap3I,CacheAlign, sizeof(double) * mHeapPts * Cx2Re);        
        
        if(RetCo1 || RetCo2 || RetCo3 || RetCo1I || RetCo2I || RetCo3I) {
            throw RetCo1 + RetCo2 + RetCo3 + RetCo1I + RetCo2I + RetCo3I;
        }

        // allocate the memory for twiddle factor and distribute them into
        // heap structure for pipeline-stage calculations. Note that there
        // is no need for pipeline zero stage in the FFT. See the comments
        // at the top of this file.

        CoWnHeap1   = new(PlaceCoWnHeap1)   double[mHeapPts * Cx2Re] {};
        CoWnHeap2   = new(PlaceCoWnHeap2)   double[mHeapPts * Cx2Re] {};
        CoWnHeap3   = new(PlaceCoWnHeap3)   double[mHeapPts * Cx2Re] {};
        
        CoWnHeap1I  = new(PlaceCoWnHeap1I) double[mHeapPts * Cx2Re] {};                
        CoWnHeap2I  = new(PlaceCoWnHeap2I) double[mHeapPts * Cx2Re] {};                
        CoWnHeap3I  = new(PlaceCoWnHeap3I) double[mHeapPts * Cx2Re] {};                        
    }
    catch(int& RetComb) {
       
        free(PlaceCoWnHeap1);
        free(PlaceCoWnHeap2);
        free(PlaceCoWnHeap3);        
        free(PlaceCoWnHeap1I);
        free(PlaceCoWnHeap2I);
        free(PlaceCoWnHeap3I);        
        
        std::cout << "CCoef4d:: _memalign failed: " << RetComb << std::endl;      
    }
            
}







#if(0)
void CCoef4d::CoefGenWn(int Quadrant, double Math2) {
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
}
#endif







void CCoef4d::CoefGenHeap(double * CoWnHeap, const long double * Wn) {
        
    if(CoWnHeap != NULL && Wn != NULL) {
        int WnStride   = mQuadPts/4;  //use N/4 instead of N/2 because the optimization:
        int WnElements = 4;           //the heap starts from the second conversion stage  
        int WnCnt = 0;                            
        
        while (WnStride > 0) {
            for(int i = 0; i < WnElements; i++) {
                CoWnHeap[WnCnt++] = (double)Wn[re(WnStride * i)];
                CoWnHeap[WnCnt++] = (double)Wn[im(WnStride * i)];            
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









void CCoef4d::CoefComp(void) {
  
    long double * WnQ1Array = NULL;
    long double * WnQ2Array = NULL;
    long double * WnQ3Array = NULL;
    
    WnQ1Array = new long double[mQuadPts * Cx2Re] {};
    WnQ2Array = new long double[mQuadPts * Cx2Re] {};
    WnQ3Array = new long double[mQuadPts * Cx2Re] {};

    CoefGenWnQ1(WnQ1Array);
    CoefGenHeap(CoWnHeap1, WnQ1Array);
    
    CoefGenWnQ2(WnQ2Array, WnQ1Array);    
    CoefGenHeap(CoWnHeap2, WnQ2Array);
    
    CoefGenWnQ3(WnQ3Array, WnQ1Array);
    CoefGenHeap(CoWnHeap3, WnQ3Array);

    delete[] WnQ1Array;    
    delete[] WnQ2Array;    
    delete[] WnQ3Array; 
}


void CCoef4d::CoefCompI(void) {

    long double * WnQ1Array = NULL;
    long double * WnQ2Array = NULL;
    long double * WnQ3Array = NULL;
    
    WnQ1Array = new long double[mQuadPts * Cx2Re] {};
    WnQ2Array = new long double[mQuadPts * Cx2Re] {};
    WnQ3Array = new long double[mQuadPts * Cx2Re] {};

    CoefGenWnQ1I(WnQ1Array);
    CoefGenHeap(CoWnHeap1I, WnQ1Array);
    
    CoefGenWnQ2I(WnQ2Array, WnQ1Array);    
    CoefGenHeap(CoWnHeap2I, WnQ2Array);
    
    CoefGenWnQ3I(WnQ3Array, WnQ1Array);
    CoefGenHeap(CoWnHeap3I, WnQ3Array);

    delete[] WnQ1Array;    
    delete[] WnQ2Array;    
    delete[] WnQ3Array; 
}





// getters and getters

const double * CCoef4d::GetCoefComp(int Quadrant) const {
    switch(Quadrant) {
        case 1: return (const double *) CoWnHeap1; 
        case 2: return (const double *) CoWnHeap2;
        case 3: return (const double *) CoWnHeap3;
        default: return NULL;        
    }

}

const double * CCoef4d::GetCoefCompI(int Quadrant) const {
    switch(Quadrant) {
        case 1: return (const double *) CoWnHeap1I; 
        case 2: return (const double *) CoWnHeap2I;
        case 3: return (const double *) CoWnHeap3I;
        default: return NULL;                
    }
}











CCoef4d::CCoef4d(const CCoef4d& orig) {
}

CCoef4d::~CCoef4d() {
    
    free(CoWnHeap3);
    free(CoWnHeap2);
    free(CoWnHeap1);

    free(CoWnHeap3I);
    free(CoWnHeap2I);
    free(CoWnHeap1I);    
}

