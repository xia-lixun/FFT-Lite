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




//default constructor
CCoef4::CCoef4() {

    CoWnHeap1 = NULL;
    CoWnHeap2 = NULL;
    CoWnHeap3 = NULL;    

    CoWnHeap1I = NULL;
    CoWnHeap2I = NULL;
    CoWnHeap3I = NULL;     
}





CCoef4::CCoef4(int Length) {
    
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
    mN = Length;
    mQuadPts = mN / 4;
    mHeapPts = 0;
    for(int i = mQuadPts; i >= 4; i = i / 4) {
        mHeapPts += i;
    }
    
    try{        
        int RetCo1 = posix_memalign(&PlaceCoWnHeap1, CacheAlign, sizeof(float) * mHeapPts * Cx2Re);
        int RetCo2 = posix_memalign(&PlaceCoWnHeap2, CacheAlign, sizeof(float) * mHeapPts * Cx2Re);
        int RetCo3 = posix_memalign(&PlaceCoWnHeap3, CacheAlign, sizeof(float) * mHeapPts * Cx2Re);        
        
        int RetCo1I = posix_memalign(&PlaceCoWnHeap1I,CacheAlign, sizeof(float) * mHeapPts * Cx2Re);
        int RetCo2I = posix_memalign(&PlaceCoWnHeap2I,CacheAlign, sizeof(float) * mHeapPts * Cx2Re);
        int RetCo3I = posix_memalign(&PlaceCoWnHeap3I,CacheAlign, sizeof(float) * mHeapPts * Cx2Re);        
        
        if(RetCo1 || RetCo2 || RetCo3 || RetCo1I || RetCo2I || RetCo3I) {
            throw RetCo1 + RetCo2 + RetCo3 + RetCo1I + RetCo2I + RetCo3I;
        }

        // allocate the memory for twiddle factor and distribute them into
        // heap structure for pipeline-stage calculations. Note that there
        // is no need for pipeline zero stage in the FFT. See the comments
        // at the top of this file.

        CoWnHeap1   = new(PlaceCoWnHeap1)   float[mHeapPts * Cx2Re] {};
        CoWnHeap2   = new(PlaceCoWnHeap2)   float[mHeapPts * Cx2Re] {};
        CoWnHeap3   = new(PlaceCoWnHeap3)   float[mHeapPts * Cx2Re] {};
        
        CoWnHeap1I  = new(PlaceCoWnHeap1I) float[mHeapPts * Cx2Re] {};                
        CoWnHeap2I  = new(PlaceCoWnHeap2I) float[mHeapPts * Cx2Re] {};                
        CoWnHeap3I  = new(PlaceCoWnHeap3I) float[mHeapPts * Cx2Re] {};                        
    }
    catch(int& RetComb) {
        
        free(PlaceCoWnHeap1);
        free(PlaceCoWnHeap2);
        free(PlaceCoWnHeap3);        
        
        free(PlaceCoWnHeap1I);
        free(PlaceCoWnHeap2I);
        free(PlaceCoWnHeap3I);        
        std::cout << "CCoef4:: _memalign failed: " << RetComb << std::endl;      
    }
        
}






#if(0)
void CCoef4::CoefGenWn(int Quadrant, double Math2) {

    // calculate the complex exponential function for Wn,
    // doing nothing if allocation fails.
    
    if(WnQ1 != NULL) {
        for(int i = 0; i < mQuadPts * Cx2Re; i += 2) {
            WnQ1[i]   = 0.0;
            WnQ1[i+1] = (double)(i/2) * (double)Quadrant;
        }
        
        Vec4d coeff(Math2 * M_PI / (double)mN);
        Vec4d a;
        Vec4d b;
        for(int i = 0; i < mQuadPts * Cx2Re; i += 4) {
            a.load(WnQ1+i);
            b = a * coeff;
            b.store(WnQ1+i);
        }
        
        Complex4d ax;
        Complex4d bx;
        for(int i = 0; i < mQuadPts * Cx2Re; i += 4) {
            ax.load(WnQ1+i);
            bx = cexp(ax);
            bx.store(WnQ1+i);
        }
    }

    // map Wn to the heap, which could be all zeros if allocation fails.

#ifdef DEBUG_CCOEFRADIX4    
    for(int i = 0; i < mQuadPts * Cx2Re; i += 1) {
        WnQ1[i] = i;
    }
#endif    
}
#endif




// { 0..N/4-1 } / N * 2 * Pi

void CCoef4::CoefGenWnQ1(long double * WnQ1) {

    long double PhaseCosine;
    long double PhaseSine;
    
    // calculate the complex exponential function for Wn,
    // doing nothing if allocation fails.
    
    if(WnQ1 != NULL) {        
        WnQ1[0] = 1.0L;  //when phase = 0
        WnQ1[1] = 0.0L;
        
        for(int i = 2; i < mQuadPts * Cx2Re; i += 2) {  
            PhaseCosine = (long double)i / (long double)mN * M_PIl;     
            PhaseSine = (long double)(mN/2 + i) / (long double)mN * M_PIl;
            WnQ1[i] = cos(PhaseCosine);
            WnQ1[i+1] = cos(PhaseSine);
        }
    }    
}

void CCoef4::CoefGenWnQ1I(long double * WnQ1) {

    long double PhaseCosine;
    long double PhaseSine;
    
    // calculate the complex exponential function for Wn,
    // doing nothing if allocation fails.
    
    if(WnQ1 != NULL) {        
        WnQ1[0] = 1.0L;  //when phase = 0
        WnQ1[1] = 0.0L;
        
        for(int i = 2; i < mQuadPts * Cx2Re; i += 2) {  
            PhaseCosine = (long double)i / (long double)mN * M_PIl;     
            PhaseSine = (long double)(i - mN/2) / (long double)mN * M_PIl;
            WnQ1[i] = cos(PhaseCosine);
            WnQ1[i+1] = cos(PhaseSine);
        }
    }    
}







// { 0..N/4-1 } * 2 / N * 2 * Pi
//
//    cos(2a) =  cos(a)*cos(a) - sin(a)*sin(a)
//   -sin(2a) = -2sin(a)cos(a)
//
// Calculation is based on the results of CCoef4::CoefGenWnQ1()

void CCoef4::CoefGenWnQ2(long double * WnQ2, const long double * WnQ1) {

    long double mCosine;
    long double mSineMinus;
    
    if(WnQ1 != NULL && WnQ2 != NULL) {        
        WnQ2[0] = 1.0L;  //when phase = 0
        WnQ2[1] = 0.0L;
        
        for(int i = 2; i < mQuadPts * Cx2Re; i += 2) {  
            mCosine = WnQ1[i];
            mSineMinus = WnQ1[i+1];
            WnQ2[i] = (mCosine * mCosine) - (mSineMinus * mSineMinus);
            WnQ2[i+1] = 2.0L * mSineMinus * mCosine;
        }
    }    
}

void CCoef4::CoefGenWnQ2I(long double * WnQ2, const long double * WnQ1) {

    long double mCosine;
    long double mSine;
    
    if(WnQ1 != NULL && WnQ2 != NULL) {        
        WnQ2[0] = 1.0L;  //when phase = 0
        WnQ2[1] = 0.0L;
        
        for(int i = 2; i < mQuadPts * Cx2Re; i += 2) {  
            mCosine = WnQ1[i];
            mSine = WnQ1[i+1];
            WnQ2[i] = (mCosine * mCosine) - (mSine * mSine);
            WnQ2[i+1] = 2.0L * mSine * mCosine;
        }
    }    
}







// { 0..N/4-1 } * 3 / N * 2 * Pi
//
//    cos(3a) = cos(a)*cos(a)*cos(a) - 3*cos(a)*sin(a)*sin(a)
//   -sin(3a) = sin(a)*sin(a)*sin(a) - 3*sin(a)*cos(a)*cos(a)
//
// Calculation is based on the results of CCoef4::CoefGenWnQ1()

void CCoef4::CoefGenWnQ3(long double * WnQ3, const long double * WnQ1) {

    long double mCosine;
    long double mSineMinus;
    
    if(WnQ1 != NULL && WnQ3 != NULL) {        
        WnQ3[0] = 1.0L;  //when phase = 0
        WnQ3[1] = 0.0L;
        
        for(int i = 2; i < mQuadPts * Cx2Re; i += 2) {  
            mCosine = WnQ1[i];
            mSineMinus = WnQ1[i+1];            
            WnQ3[i] = (mCosine * mCosine * mCosine) - (3.0L * mCosine * mSineMinus * mSineMinus);
            WnQ3[i+1] = (3.0L * mSineMinus * mCosine * mCosine) - (mSineMinus * mSineMinus * mSineMinus);
        }
    }    
}

void CCoef4::CoefGenWnQ3I(long double * WnQ3, const long double * WnQ1) {

    long double mCosine;
    long double mSine;
    
    if(WnQ1 != NULL && WnQ3 != NULL) {        
        WnQ3[0] = 1.0L;  //when phase = 0
        WnQ3[1] = 0.0L;
        
        for(int i = 2; i < mQuadPts * Cx2Re; i += 2) {  
            mCosine = WnQ1[i];
            mSine = WnQ1[i+1];            
            WnQ3[i] = (mCosine * mCosine * mCosine) - (3.0L * mCosine * mSine * mSine);
            WnQ3[i+1] = (3.0L * mSine * mCosine * mCosine) - (mSine * mSine * mSine);
        }
    }    
}









void CCoef4::CoefGenHeap(float * CoWnHeap, const long double * WnQ) {
        
    if(CoWnHeap != NULL && WnQ != NULL) {
        int WnStride   = mQuadPts/4;  //use N/4 instead of N/2 because the optimization:
        int WnElements = 4;           //the heap starts from the second conversion stage  
        int WnCnt = 0;                            
        
        while (WnStride > 0) {
            for(int i = 0; i < WnElements; i++) {
                CoWnHeap[WnCnt++] = (float)WnQ[re(WnStride * i)];
                CoWnHeap[WnCnt++] = (float)WnQ[im(WnStride * i)];            
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


void CCoef4::CoefCompI(void) {
    
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

    free(CoWnHeap3I);
    free(CoWnHeap2I);
    free(CoWnHeap1I);
}

