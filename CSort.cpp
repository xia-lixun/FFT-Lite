/* 
 * File:   CBitReverse.cpp
 * Author: lixun
 * 
 * Created on January 24, 2015, 11:11 AM
 */

#include <stdlib.h>
#include <mm_malloc.h>
#include <new>
#include <iostream>

#include "CSort.h"
#include "Config.h"







CSort::CSort(int Length) {
    
    // length for FFT sort
    mN = Length;

    // allocate cache-aligned memory for bit-reverse lookup table
    mTreeEven = NULL;
    mLutComp = NULL;
    mLutReal = NULL;
    
    void * PlaceTreeEven = NULL;
    void * PlaceLutComp = NULL;
    void * PlaceLutReal = NULL;
    try{
        int RetTreeEven = posix_memalign(&PlaceTreeEven, CacheAlign, sizeof(int) * mN/4);
        int RetLutComp  = posix_memalign(&PlaceLutComp,  CacheAlign, sizeof(int) * mN * Cx2Re);
        int RetLutReal  = posix_memalign(&PlaceLutReal,  CacheAlign, sizeof(int) * mN);
        
        if(RetTreeEven || RetLutComp || RetLutReal) {
            throw RetTreeEven + RetLutComp + RetLutReal;
        }

        // allocate memories for the bit-reverse mapping scheme
        // both complex and real cases are considered.

        mTreeEven = new(PlaceTreeEven) int[mN/4] {};
        mLutComp  = new(PlaceLutComp) int[mN*2] {};
        mLutReal  = new(PlaceLutReal) int[mN] {};
    }
    catch(int& RetCombined) {
        free(PlaceTreeEven);
        free(PlaceLutComp);
        free(PlaceLutReal);
        std::cout << "_memalign failed: " << RetCombined << std::endl;
    }
}





void CSort::TreeEven(void) {
    //
    // even number index generated with growing tree method:
    // N=4    {0}
    // N=8    {0} {N/4}
    // N=16   {0   N/4} {0+N/8   N/4+N/8}
    // N=32   {0   N/4   0+N/8   N/4+N/8} {0+N/16   N/4+N/16   0+N/8+N/16   N/4+N/8+N/16}
    // N=...   ...
    //
    if(mTreeEven != NULL) {
        int Nd = mN/4;
        mTreeEven[0] = 0;
        for(int i = 1; i <= mN/8; i <<= 1) {
            //i is the same as the segment length
            for(int k = 0; k < i; k++) {
                mTreeEven[i + k] = mTreeEven[k] + Nd;
            }
            Nd = Nd/2;
        }
    }// doing nothing if memory is not ready 
}// end of GenerateTreeEven()






void CSort::LutComp(void) {

    TreeEven();
    //
    //generate global bit reverse map
    //
    if(mLutComp != NULL && mTreeEven != NULL) {        
        for(int i = 0; i < mN/4; i++) {
            mLutComp[(i*2) * 2 + 0] = mTreeEven[i] * 2 + 0;
            mLutComp[(i*2) * 2 + 1] = mTreeEven[i] * 2 + 1;
            
            mLutComp[(i*2+1) * 2 + 0] = (mTreeEven[i] + mN/2) * 2 + 0;
            mLutComp[(i*2+1) * 2 + 1] = (mTreeEven[i] + mN/2) * 2 + 1;
        }
        for(int i = 0; i < mN; i++) {
            mLutComp[i+mN] = mLutComp[i] + 2;
        }    
    }// do nothing if memory is not ready
}// end of GenerateComplex()






void CSort::LutReal(void) {
    
    TreeEven();

    //generate global bit reverse map

    if(mLutReal != NULL && mTreeEven != NULL) {        
        for(int i = 0; i < mN/4; i++) {
            mLutReal[i*2] = mTreeEven[i];
            mLutReal[i*2+1] = mTreeEven[i] + mN/2;
        }
        for(int i = 0; i < mN/2; i++) {
            mLutReal[i+mN/2] = mLutReal[i] + 1;
        }    
    }// do nothing is memory is not ready    
}// end of GenerateReal()








// getters and setters

const int * CSort::GetLutComp(void) const {
    return (const int *) mLutComp;
}

const int * CSort::GetLutReal(void) const {
    return (const int *) mLutReal;
}

int CSort::GetPoints(void) const {
    return mN;
}








CSort::CSort(const CSort& orig) {
}

CSort::~CSort() {
    free(mLutReal);
    free(mLutComp);
    free(mTreeEven);
}



    
    
