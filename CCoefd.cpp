/* 
 * File:   CTwiddleFactorDP.cpp
 * Author: lixun
 * 
 * Created on January 30, 2015, 7:57 PM
 */

//
//  twiddle factor lookup table generation: Wn(0) Wn(1) ... Wn(N/2-1)
//
//  calculate twiddle factor lookup table for each
//  butter-fly pipeline of depth log2N
//
//        -i 2Pi k/N
//       e              k=[0..N/2-1]
//
//  data structure would be heap, for instance N=16
//  note that each cell is a complex number.
//
//  +===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+
//  | 0 | 1 | 1 | 2 | 2 | 2 | 2 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
//  +===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+
//
//  For the sake of better performance, there is no need to store
//  the first element because of the following reasons:
//    (1) the value of the coefficient is always one
//    (2) cache-line alignment is easier to be achieved
//  Hence the resultant heap would look like:
//
//  +===+===+===+===+===+===+===+===+===+===+===+===+===+===+
//  | 1 | 1 | 2 | 2 | 2 | 2 | 4 | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
//  +===+===+===+===+===+===+===+===+===+===+===+===+===+===+    
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

#include "CCoefd.h"
#include "Config.h"



CCoefd::CCoefd(int Length) {
    mN = Length;
    //
    // allocate cache-aligned memory for factor coefficient lookup table
    //
    Wn = NULL;
    CoWnHeap = NULL;
    CoWnHeapI = NULL;
    
    void * PlaceWn        = NULL;
    void * PlaceCoWnHeap  = NULL;
    void * PlaceCoWnHeapI = NULL;
    
    try{
        int RetWn        = posix_memalign(&PlaceWn,        CacheAlign, sizeof(double)*mN);
        int RetCoWnHeap  = posix_memalign(&PlaceCoWnHeap,  CacheAlign, sizeof(double)*(mN-2)*Cx2Re);
        int RetCoWnHeapI = posix_memalign(&PlaceCoWnHeapI, CacheAlign, sizeof(double)*(mN-2)*Cx2Re);
        
        if(RetWn || RetCoWnHeap || RetCoWnHeapI) {
            throw RetWn + RetCoWnHeap + RetCoWnHeapI;
        }
        //
        // allocate the memory for twiddle factor and distribute them into
        // heap structure for pipeline-stage calculations. Note that there
        // is no need for pipeline zero stage in the FFT. See the comments
        // at the top of this file.
        //
        Wn        = new(PlaceWn)        double[mN]           {};
        CoWnHeap  = new(PlaceCoWnHeap)  double[(mN-2)*Cx2Re] {};        
        CoWnHeapI = new(PlaceCoWnHeapI) double[(mN-2)*Cx2Re] {};                
    }
    catch(int& RetComb) {
        free(PlaceWn);
        free(PlaceCoWnHeap);
        free(PlaceCoWnHeapI);        
        std::cout << "_memalign failed: " << RetComb << std::endl;      
    }
}





void CCoefd::CoefGen(double * CoWnHp, double Math2) {
    
    //
    // calculate the complex exponential function for Wn,
    // doing nothing if allocation fails.
    //
    if(Wn != NULL) {
        for(int i = 0; i < mN; i+=2) {
            Wn[i]   = 0.0;
            Wn[i+1] = (double)(i/2);
        }
        
        Vec4d coeff(Math2 * M_PI / (double)mN);
        Vec4d a;
        Vec4d b;
        for(int i = 0; i < mN; i+=4) {
            a.load(Wn+i);
            b = a * coeff;
            b.store(Wn+i);
        }
        
        Complex4d ax;
        Complex4d bx;
        for(int i = 0; i < mN; i+=4) {
            ax.load(Wn+i);
            bx = cexp(ax);
            bx.store(Wn+i);
        }
    }
    //
    // map Wn to the heap, which could be all zeros if allocation fails.
    //
    if(CoWnHp != NULL && Wn != NULL) {
        int WnStride  = mN/4;  //use N/4 instead of N/2 because the optimization:
        int WnElements = 2;   //the heap starts from the second conversion stage  
        int WnCnt = 0;                            
        
        while (WnStride > 0) {
            for(int i = 0; i < WnElements; i++) {
                CoWnHp[WnCnt++] = Wn[re(WnStride * i)];
                CoWnHp[WnCnt++] = Wn[im(WnStride * i)];            
            }
            WnElements = WnElements * 2;
            WnStride = WnStride / 2;
        }
    }
    //
    //end of twiddle factor calculation and mapping    
        
}





void CCoefd::CoefComp(void) {
    CoefGen(CoWnHeap, -2.0);
}


void CCoefd::CoefCompI(void) {
    CoefGen(CoWnHeapI, 2.0);
}






//
// getters and getters
//
const double * CCoefd::GetCoefComp(void) const {
    return (const double *) CoWnHeap;
}

const double * CCoefd::GetCoefCompI(void) const {
    return (const double *) CoWnHeapI;
}

int CCoefd::GetPoints(void) const {
    return mN;
}


CCoefd::CCoefd(const CCoefd& orig) {
}

CCoefd::~CCoefd() {
    free(CoWnHeap);
    free(Wn);
    free(CoWnHeapI);
}
