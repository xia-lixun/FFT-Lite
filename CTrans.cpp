/* 
 * File:   CFastFourierTransform.cpp
 * Author: lixun
 * 
 * Created on January 21, 2015, 11:51 PM
 */



//#define _USE_MATH_DEFINES
#include <cmath>
#include <new>
#include <cstdlib>
#include <mm_malloc.h>
#include <iostream>

#include "CTrans.h"
#include "Config.h"






//constructor for float type
CTrans::CTrans(const CSort& Map, const CCoef& Wn) {
    
    xSIMDd[0] = NULL;
    xSIMDd[1] = NULL;

    // Wn and bit-reverse map shall be handled by other classes

    SortLutComp = Map.GetLutComp();
    CoWnHeap    = Wn.GetCoefComp(); 
    CoWnHeapI   = Wn.GetCoefCompI();

    // allocate butter-fly pipelines for out-of-place conversions
    // note that they are address aligned to 16 byte boundaries

    xSIMD[0] = NULL;
    xSIMD[1] = NULL;
    
    void * PlacePin = NULL;
    void * PlacePon = NULL;
    
    if( Map.GetPoints() == Wn.GetPoints() ) {
        mN = Map.GetPoints();
        try{        
            int RetPin = posix_memalign(&PlacePin, CacheAlign, sizeof(float) * mN * Cx2Re);
            int RetPon = posix_memalign(&PlacePon, CacheAlign, sizeof(float) * mN * Cx2Re);    
            if(RetPin || RetPon) {
                throw RetPin + RetPon;
            }
            
            xSIMD[0] = new(PlacePin) float[mN * Cx2Re] {};
            xSIMD[1] = new(PlacePon) float[mN * Cx2Re] {};    
        }
        catch(int& RetComb) {
            free(PlacePin);  //free(NULL) does nothing.
            free(PlacePon);
            std::cout << "_memalign failed: " << RetComb << std::endl;              
        }
    } else {
        std::cout << "CTrans:: _map_N / _coef_N mismatch: no mem allocated" << Map.GetPoints() << "/" << Wn.GetPoints() << std::endl;              
    }

    // simple route info for bit-reversing operations of the first
    // and second stage butter-fly pipelines.
#if(0)
    RouteLUT[0] = 0;
    RouteLUT[1] = 1;
    RouteLUT[2] = 4;
    RouteLUT[3] = 5;
                                
    RouteLUT[4] = 2;
    RouteLUT[5] = 3;
    RouteLUT[6] = 6;
    RouteLUT[7] = 7;

    XeRoute.load(RouteLUT);   
    XoRoute.load(RouteLUT+4); 
#endif
    // load twiddle factor of Wn(0) and Wn(N/4) from the heap for the second 
    // stage pipeline.
}




//constructor for double type
CTrans::CTrans(const CSort& Map, const CCoefd& Wn) {
    
    xSIMD[0] = NULL;
    xSIMD[1] = NULL;  

    // Wn and bit-reverse map shall be handled by other classes

    SortLutComp = Map.GetLutComp();
    CoWnHeapd   = Wn.GetCoefComp(); 
    CoWnHeapId  = Wn.GetCoefCompI();

    // allocate butter-fly pipelines for out-of-place conversions
    // note that they are address aligned to 16 byte boundaries

    xSIMDd[0] = NULL;
    xSIMDd[1] = NULL;
    
    void * PlacePin = NULL;
    void * PlacePon = NULL;
    if( Map.GetPoints() == Wn.GetPoints() ) {
        mN = Map.GetPoints();
        try{
            int RetPin = posix_memalign(&PlacePin, CacheAlign, sizeof(double) * mN * Cx2Re);
            int RetPon = posix_memalign(&PlacePon, CacheAlign, sizeof(double) * mN * Cx2Re);    
            if(RetPin || RetPon) {
                throw RetPin + RetPon;
            }
            
            xSIMDd[0] = new(PlacePin) double[mN * Cx2Re] {};
            xSIMDd[1] = new(PlacePon) double[mN * Cx2Re] {};    
        }
        catch(int& RetComb) {
            free(PlacePin);  //free(NULL) does nothing.
            free(PlacePon);
            std::cout << "_memalign failed: " << RetComb << std::endl;              
        }
    } else {
        std::cout << "CTrans:: _map_N / _coef_N mismatch: no mem allocated" << Map.GetPoints() << "/" << Wn.GetPoints() << std::endl;              
    }

    // simple route info for bit-reversing operations of the first
    // and second stage butter-fly pipelines.

#if(0)
    RouteLUT[0] = 0;
    RouteLUT[1] = 1;
    RouteLUT[2] = 4;
    RouteLUT[3] = 5;
                                
    RouteLUT[4] = 2;
    RouteLUT[5] = 3;
    RouteLUT[6] = 6;
    RouteLUT[7] = 7;

    XeRoute.load(RouteLUT);   
    XoRoute.load(RouteLUT+4); 
#endif
    // load twiddle factor of Wn(0) and Wn(N/4) from the heap for the second 
    // stage pipeline.    
}








void CTrans::convert(const double * CoWnHpd) {

    // execution state: decimation in time.
    // step 1: remap the input array in bit-reverse pattern
    //         xSIMDd[0] ==/bit reverse indexed/==> xSIMDd[1]

    for(int i = 0; i < mN*Cx2Re; i++) {
        xSIMDd[1][i] = xSIMDd[0][SortLutComp[i]];
    }    

    // step2: iterate over all butter-fly patterns
    //
    //        2.0 calculate the first conversion pipeline
    //            xSIMDd[1] ==//==> xSIMDd[0]

    for(int i = 0; i < mN*2; i += 4) {
        // d-cache friendly
        XeValueInC2d.load(xSIMDd[1] + i);
        XoValueInC2d.load(xSIMDd[1] + i + 2);
        
        XeValueUtC2d = XeValueInC2d + XoValueInC2d;
        XoValueUtC2d = XeValueInC2d - XoValueInC2d;
        
        XeValueUtC2d.store(xSIMDd[0] + i);
        XoValueUtC2d.store(xSIMDd[0] + i + 2);
    }    

    //        2.1 calculate the second conversion pipeline
    //            xSIMDd[0] ==//==> xSIMDd[1]    

    CoWn4d.load(CoWnHpd);
    for(int i = 0; i < mN*2; i += 8) {  
        //        
          XoComplexd.load(xSIMDd[0] + i + 4);        
        XeComplexInd.load(xSIMDd[0] + i);
        XoComplexInd = XoComplexd * CoWn4d;
        
        XeComplexUtd = XeComplexInd + XoComplexInd;
        XoComplexUtd = XeComplexInd - XoComplexInd;
        
        XeComplexUtd.store(xSIMDd[1] + i);
        XoComplexUtd.store(xSIMDd[1] + i + 4);
    }
 
    //        2.2 iterate over the rest of the stages using Complex8f
    //            input: xSIMDd[1]
    //              end: xSIMDd[Result]
 
    pind = xSIMDd[1];
    pond = xSIMDd[0];
    pinpond = NULL;
    Result = 0x1;
    
    CoWnSet = 4;      // initial offset in the Wn heap
    BiN = mN/8;        // initial number of butter-fly groups
 
    // loop for butter-fly pipelines
 
    for(int BuLn = 4; BuLn < mN; BuLn = BuLn << 1) {   
 
        // loop for butter-fly SIMD groups
 
        XoeSet  = 0;
        for(int BiP = 0; BiP < BuLn/2; BiP++) {
 
            // loop for all butter-fly groups
 
            CoWn4d.load(CoWnHpd + CoWnSet);
            for(int Bi = 0; Bi < BiN; Bi++) {        
                
                  XoComplexd.load(pind + im(Bi)*BuLn*Cx2Re + XoeSet);
                XeComplexInd.load(pind + re(Bi)*BuLn*Cx2Re + XoeSet);
                XoComplexInd = XoComplexd * CoWn4d;
                
                XeComplexUtd = XeComplexInd + XoComplexInd;
                XoComplexUtd = XeComplexInd - XoComplexInd;
                
                XeComplexUtd.store(pond + re(Bi)*BuLn*Cx2Re + XoeSet);
                XoComplexUtd.store(pond + im(Bi)*BuLn*Cx2Re + XoeSet);
            }
            CoWnSet += (2*Cx2Re);
            XoeSet  += (2*Cx2Re);
        }
        pinpond = pind;
        pind = pond;
        pond = pinpond;
        BiN = BiN >> 1;
        Result = (Result + 1) & 0x1;        
    }    
}





void CTrans::ForwardDouble(void) {
    convert(CoWnHeapd);
}


void CTrans::InverseDouble(void) {
    convert(CoWnHeapId);
}










void CTrans::convert(const float * CoWnHp) {
    
    // execution state: decimation in time.
    // step 1: remap the input array in bit-reverse pattern
    //         xSIMD[0] ==/bit reverse indexed/==> xSIMD[1]
    //
    // I don't see obvious loop-carried dependency chain here
    // so no duplications will be made here. Am I right?
    //
    // I can see there is a lot of d-cache misses here, as
    // extracting non-contiguous samples from xSIMD[0] is the  
    // cause. However, writing data back to xSIMD[1] must
    // be smooth for cache-blocking.
    //
    //    for(int i = 0; i < N*2; i += 8) {
    //        ReMapRoute8.load(BitReverseMapComplex + i);
    //        ReMapValue8 = lookup<N*2>(ReMapRoute8, xSIMD[0]);
    //        ReMapValue8.store(xSIMD[1] + i);
    //    }    
    //
    // alternative re-mapping method:

    for(int i = 0; i < mN * Cx2Re; i++) {
        xSIMD[1][i] = xSIMD[0][SortLutComp[i]];
    }    

    // step2: iterate over all butter-fly patterns
    //
    //     2.0 calculate the first conversion pipeline
    //         xSIMD[1] ==//==> xSIMD[0]

    for(int i = 0; i < mN * Cx2Re; i += 4) {
        // d-cache friendly
        XeValueInC2f.load(xSIMD[1] + i);
        XoValueInC2f.load(xSIMD[1] + i + 2);
        
        XeValueUtC2f = XeValueInC2f + XoValueInC2f;
        XoValueUtC2f = XeValueInC2f - XoValueInC2f;
        
        XeValueUtC2f.store(xSIMD[0] + i);
        XoValueUtC2f.store(xSIMD[0] + i + 2);
    }    

    //     2.1 calculate the second conversion pipeline
    //         xSIMD[0] ==//==> xSIMD[1]    

    CoWn.load(CoWnHp);
    for(int i = 0; i < mN * Cx2Re; i += 8) {  

          XoComplex.load(xSIMD[0] + i + 4);        
        XeComplexIn.load(xSIMD[0] + i);
        XoComplexIn = XoComplex * CoWn;
        
        XeComplexUt = XeComplexIn + XoComplexIn;
        XoComplexUt = XeComplexIn - XoComplexIn;
        
        XeComplexUt.store(xSIMD[1] + i);
        XoComplexUt.store(xSIMD[1] + i + 4);
    }

    //        2.2 iterate over the rest of the stages using Complex8f
    //            input: xSIMD[1]
    //              end: xSIMD[Result]

    pin = xSIMD[1];
    pon = xSIMD[0];
    pinpon = NULL;
    Result = 0x1;
    
    CoWnSet = 4;     // initial offset in the Wn heap
    BiN = mN/8;      // initial number of butter-fly groups

    // loop for butter-fly pipelines

    for(int BuLn = 4; BuLn < mN; BuLn = BuLn << 1) {   

        // loop for butter-fly SIMD groups

        XoeSet  = 0;
        for(int BiP = 0; BiP < BuLn/4; BiP++) {

            // loop for all butter-fly groups

            CoWn8.load(CoWnHp + CoWnSet);
            for(int Bi = 0; Bi < BiN; Bi++) {        
                
                  XoComplex8.load(pin + im(Bi)*BuLn*Cx2Re + XoeSet);
                XeComplexIn8.load(pin + re(Bi)*BuLn*Cx2Re + XoeSet);
                XoComplexIn8 = XoComplex8 * CoWn8;
                
                XeComplexUt8 = XeComplexIn8 + XoComplexIn8;
                XoComplexUt8 = XeComplexIn8 - XoComplexIn8;
                
                XeComplexUt8.store(pon + re(Bi)*BuLn*Cx2Re + XoeSet);
                XoComplexUt8.store(pon + im(Bi)*BuLn*Cx2Re + XoeSet);
            }
            CoWnSet += (4*Cx2Re);
            XoeSet  += (4*Cx2Re);
        }
        pinpon = pin;
        pin = pon;
        pon = pinpon;
        BiN = BiN >> 1;
        Result = (Result + 1) & 0x1;        
    }
}






void CTrans::ForwardFloat(void) {
    convert(CoWnHeap);
}


void CTrans::InverseFloat(void) {
    convert(CoWnHeapI);
}












#if(0)
void CTrans::Forward4F(void) {
    //
    // execution stage: decimation in time.
    // step 1: re-map the input array in bit-reverse pattern
    //         xSIMD[0] ==/bit reverse indexed/==> xSIMD[1]
    //
    //    Vec4i ReMapRoute;
    //    Vec4f ReMapValue;
    //
    for(int i = 0; i < mN*2; i += 4) {
        ReMapRoute.load(SortLutComp + i);
        ReMapValue = lookup<mN*2>(ReMapRoute, xSIMD[0]);
        ReMapValue.store(xSIMD[1] + i);
    }
    //
    //
    //step 2: iterate over butter-fly patterns
    //
    //     2.0: calculate the first conversion
    //          xSIMD[1] ==//==> xSIMD[0]
    //
    //    Vec4f XeValueIn;
    //    Vec4f XoValueIn;
    //    Vec4i XeRoute(0, 1, 4, 5);
    //    Vec4i XoRoute(2, 3, 6, 7);
    //    Vec4f XeValueUt;
    //    Vec4f XoValueUt;
    //    Vec4f XeValue;
    //    Vec4f XoValue; 
    //
    for(int i = 0; i < mN*2; i += 8) {
        XeValueIn = lookup<8>(XeRoute, xSIMD[1] + i);
        XoValueIn = lookup<8>(XoRoute, xSIMD[1] + i);
        
        XeValueUt = XeValueIn + XoValueIn;
        XoValueUt = XeValueIn - XoValueIn;
        
        XeValue = blend4f<0, 1, 4, 5>(XeValueUt, XoValueUt);
        XoValue = blend4f<2, 3, 6, 7>(XeValueUt, XoValueUt);
        
        XeValue.store(xSIMD[0] + i);
        XoValue.store(xSIMD[0] + i + 4);    
    }
    //
    //step 2.1: calculate the second conversion
    //          
    //
    //          apply Wn coefficients to the input
    //
    //    Complex4f CoWn;
    //    Complex4f XoComplexIn;
    //    Complex4f XeComplexIn;
    //    Complex4f XoComplexUt;
    //    Complex4f XeComplexUt;
    //    Complex4f XoComplex;
    //    Complex4f XeComplex;    
    
    //    int XoeSet;
    //    int CoWnSet;
    //    int BiN;
    //    float * pin = xSIMD[0];
    //    float * pon = xSIMD[1];
    //    float * pinpon = NULL;
    //
    //
    //       Butter-fly group: repeating butter-fly structure.
    //Butter-fly group length: number of different twiddle coefficients
    //                         at every conversion stage.
       
    //step 2.2: calculate the second and the rest of the conversions
    //          
    //
    pin = xSIMD[0];
    pon = xSIMD[1];
    pinpon = NULL;
    Result = 0x0;
    
    CoWnSet = 0;
    BiN = mN/4;
    for(int BuLn = 2; BuLn < mN; BuLn = BuLn << 1) {
        
        XoeSet  = 0;
        for(int BiP = 0; BiP < BuLn/2; BiP++) {

            CoWn.load(CoWnHeap + CoWnSet);
            for(int Bi = 0; Bi < BiN; Bi++) {        
                
                  XoComplex.load(pin + im(Bi)*BuLn*Cx2Re + XoeSet);
                XeComplexIn.load(pin + re(Bi)*BuLn*Cx2Re + XoeSet);
                XoComplexIn = XoComplex * CoWn;
                
                XeComplexUt = XeComplexIn + XoComplexIn;
                XoComplexUt = XeComplexIn - XoComplexIn;
                
                XeComplexUt.store(pon + re(Bi)*BuLn*Cx2Re + XoeSet);
                XoComplexUt.store(pon + im(Bi)*BuLn*Cx2Re + XoeSet);
            }
            CoWnSet += (2*Cx2Re);
            XoeSet  += (2*Cx2Re);
        }
        pinpon = pin;
        pin = pon;
        pon = pinpon;
        BiN = BiN >> 1;
        Result = (Result + 1) & 0x1;        
    }
}
#endif




//
// setters and getters
// Note: input vector could be destroyed!
//       resultant vector can be modified directly!
//
float * CTrans::GetInputFloat(void) {
    return xSIMD[0];
}


float * CTrans::GetOutputFloat(void) {
    return xSIMD[Result];
}


double * CTrans::GetInputDouble(void) {
    return xSIMDd[0];
}


double * CTrans::GetOutputDouble(void) {
    return xSIMDd[Result];
}

int CTrans::GetPoints(void) {
    return mN;
}







CTrans::CTrans(const CTrans& orig) {
}

CTrans::~CTrans() {
    free(xSIMD[0]);
    free(xSIMD[1]);
    free(xSIMDd[0]);
    free(xSIMDd[1]);    
}

