/* 
 * File:   CSort4.cpp
 * Author: lixun
 * 
 * Created on 2015年3月14日, 下午11:25
 */
#include <stdlib.h>
#include <mm_malloc.h>
#include <new>
#include <iostream>
#include <bitset>

#include "CSort4.h"
#include "Config.h"




//#define DEBUG_CSORT_RADIX_4




CSort4::CSort4(int Length) {

    // length for FFT sort
    mN = Length;

    // allocate cache-aligned memory for bit-reverse lookup table
    mTreeQuad = NULL;
    mLutComp = NULL;
    mLutReal = NULL;
    
    void * PlaceTreeQuad = NULL;
    void * PlaceLutComp = NULL;
    void * PlaceLutReal = NULL;
    try{
        int RetTreeQuad = posix_memalign(&PlaceTreeQuad, CacheAlign, sizeof(int) * mN / 4);
        int RetLutComp  = posix_memalign(&PlaceLutComp,  CacheAlign, sizeof(int) * mN * Cx2Re);
        int RetLutReal  = posix_memalign(&PlaceLutReal,  CacheAlign, sizeof(int) * mN);
        
        if(RetTreeQuad || RetLutComp || RetLutReal) {
            throw RetTreeQuad + RetLutComp + RetLutReal;
        }

        // allocate memories for the bit-reverse mapping scheme
        // both complex and real cases are considered.

        mTreeQuad = new(PlaceTreeQuad) int[mN/4] {};
        mLutComp  = new(PlaceLutComp)  int[mN*2] {};
        mLutReal  = new(PlaceLutReal)  int[mN] {};
    }
    catch(int& RetComb) {
        free(PlaceTreeQuad);
        free(PlaceLutComp);
        free(PlaceLutReal);
        std::cout << "_memalign failed: " << RetComb << std::endl;
    }    
}



void CSort4::TreeQuad() {
    
    // quad number index generated with growing tree method:
    // [0 1/4 2/4 3/4] = [A]
    // [A][A+1/16][A+2/16][A+3/16] = [B]
    // [B][B+1/64][B+2/64][B+3/64] = [C]
    // ...

    if(mTreeQuad != NULL) {

        int Nd = mN / 16;
        
        mTreeQuad[0] = 0;
        mTreeQuad[1] = mN / 4;
        mTreeQuad[2] = mN / 2;
        mTreeQuad[3] = mN / 4 * 3;
        
        for(int i = 4; i < mN/4; i <<= 2) {
            
            for(int k = 0; k < i; k++) {
                mTreeQuad[1*i + k] = mTreeQuad[k] + Nd * 1;
                mTreeQuad[2*i + k] = mTreeQuad[k] + Nd * 2;
                mTreeQuad[3*i + k] = mTreeQuad[k] + Nd * 3;
            }
            Nd = Nd / 4;
        }
    }// doing nothing if memory is not ready 
    
    
#ifdef DEBUG_CSORT_RADIX_4
    for(int i = 0; i < mN / 4; i++) {
        //std::bitset<6> x(mTreeQuad[i]); 
        //std::cout << x << std::endl;
        std::cout << mTreeQuad[i] << std::endl;
    }
#endif
}




void CSort4::LutReal(void) {
    
    TreeQuad();
    
    //generate global bit reverse map

    if(mLutReal != NULL && mTreeQuad != NULL) {        
        for(int i = 0; i < mN/4; i++) {
            mLutReal[i] = mTreeQuad[i];
            mLutReal[1 * mN/4 + i] = mTreeQuad[i] + 1;
            mLutReal[2 * mN/4 + i] = mTreeQuad[i] + 2;            
            mLutReal[3 * mN/4 + i] = mTreeQuad[i] + 3;               
        }
    }// do nothing is memory is not ready        
}





void CSort4::LutComp(void) {

    TreeQuad();

    //generate global bit reverse map

    if(mLutComp != NULL && mTreeQuad != NULL) {    
        
        for(int i = 0; i < mN/4; i++) {

            mLutComp[re(i)] = re(mTreeQuad[i]);
            mLutComp[im(i)] = im(mTreeQuad[i]);
            
            mLutComp[re(1 * mN/4 + i)] = re(mTreeQuad[i] + 1);
            mLutComp[im(1 * mN/4 + i)] = im(mTreeQuad[i] + 1);
            
            mLutComp[re(2 * mN/4 + i)] = re(mTreeQuad[i] + 2);
            mLutComp[im(2 * mN/4 + i)] = im(mTreeQuad[i] + 2);

            mLutComp[re(3 * mN/4 + i)] = re(mTreeQuad[i] + 3);
            mLutComp[im(3 * mN/4 + i)] = im(mTreeQuad[i] + 3);
                       
        }
    }// do nothing if memory is not ready
    
#ifdef DEBUG_CSORT_RADIX_4
    for(int i = 0; i < mN*Cx2Re; i++) {
        //std::bitset<6> x(mTreeQuad[i]); 
        //std::cout << x << std::endl;
        std::cout << mLutComp[i] << std::endl;
    }    
#endif
}// end of GenerateComplex()






// getters and setters

const int * CSort4::GetLutComp(void) const {
    return (const int *) mLutComp;
}

const int * CSort4::GetLutReal(void) const {
    return (const int *) mLutReal;
}

int CSort4::GetPoints(void) const {
    return mN;
}






CSort4::CSort4(const CSort4& orig) {
}

CSort4::~CSort4() {
    free(mLutReal);
    free(mLutComp);
    free(mTreeQuad);    
}

