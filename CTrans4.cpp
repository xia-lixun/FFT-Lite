/* 
 * File:   CTransCompRadix4.cpp
 * Author: lixun
 * 
 * Created on 2015年3月12日, 上午12:34
 */



//#define _USE_MATH_DEFINES
#include <cmath>
#include <new>
#include <cstdlib>
#include <mm_malloc.h>
#include <iostream>

#include "Config.h"
#include "CTrans4.h"



//#define DEBUG_TRANS_RADIX_4



#define FFTL_CONVERT_INIT_FORWARD_FLOAT    (0)
#define FFTL_CONVERT_INIT_BACKWARD_FLOAT   (1)
#define FFTL_CONVERT_INIT_FORWARD_DOUBLE   (2)
#define FFTL_CONVERT_INIT_BACKWARD_DOUBLE  (3)





CTrans4::CTrans4(const CSort4& Map, const CCoef4& Wn) {
    
    xSIMDd[0] = NULL;
    xSIMDd[1] = NULL;
    
        
    // Wn and bit-reverse map shall be handled by other classes
    
    SortLutComp = Map.GetLutComp();
    
    CoWnHeapQ1  = Wn.GetCoefComp(1);
    CoWnHeapQ2  = Wn.GetCoefComp(2);
    CoWnHeapQ3  = Wn.GetCoefComp(3);
    
    CoWnHeapIQ1 = Wn.GetCoefCompI(1);
    CoWnHeapIQ2 = Wn.GetCoefCompI(2);
    CoWnHeapIQ3 = Wn.GetCoefCompI(3);    
 
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
        std::cout << "CTrans4:: _map_N / _coef_N mismatch: no mem allocated" << Map.GetPoints() << "/" << Wn.GetPoints() << std::endl;              
    }
    
    ComplexJ = Complex8f(0.0f, 1.0f);
    ComplexJ_ = Complex8f(0.0f, -1.0f);
}







CTrans4::CTrans4(const CSort4& Map, const CCoef4d& Wn) {
    
    xSIMD[0] = NULL;
    xSIMD[1] = NULL;    
        
    // Wn and bit-reverse map shall be handled by other classes
    
    SortLutComp = Map.GetLutComp();
    
    CoWnHeapQ1d  = Wn.GetCoefComp(1);
    CoWnHeapQ2d  = Wn.GetCoefComp(2);
    CoWnHeapQ3d  = Wn.GetCoefComp(3);
    
    CoWnHeapIQ1d = Wn.GetCoefCompI(1);
    CoWnHeapIQ2d = Wn.GetCoefCompI(2);
    CoWnHeapIQ3d = Wn.GetCoefCompI(3);    
 
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
        std::cout << "CTrans4:: _map_N / _coef_N mismatch: no mem allocated" << Map.GetPoints() << "/" << Wn.GetPoints() << std::endl;              
    }
    
    ComplexJd = Complex4d(0.0, 1.0);
    ComplexJd_ = Complex4d(0.0, -1.0);    
}












void CTrans4::convert_init(int FFTL_CONVERT_INIT_CODE) {

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
    
    if(FFTL_CONVERT_INIT_CODE == FFTL_CONVERT_INIT_FORWARD_FLOAT || FFTL_CONVERT_INIT_CODE == FFTL_CONVERT_INIT_BACKWARD_FLOAT) {
        for(int i = 0; i < mN * Cx2Re; i++) {
            xSIMD[1][i] = xSIMD[0][SortLutComp[i]];
#ifdef DEBUG_TRANS_RADIX_4        
            std::cout << xSIMD[1][i] << std::endl;
#endif        
        }   
    } else if (FFTL_CONVERT_INIT_CODE == FFTL_CONVERT_INIT_FORWARD_DOUBLE || FFTL_CONVERT_INIT_CODE == FFTL_CONVERT_INIT_BACKWARD_DOUBLE) {
        for(int i = 0; i < mN * Cx2Re; i++) {
            xSIMDd[1][i] = xSIMDd[0][SortLutComp[i]];
#ifdef DEBUG_TRANS_RADIX_4        
            std::cout << xSIMD[1][i] << std::endl;
#endif        
        }           
    }
    
    
    
    // step2: iterate over all butter-fly patterns
    //
    //     2.0 calculate the first conversion pipeline
    //         xSIMD[1] ==//==> xSIMD[0]

    switch(FFTL_CONVERT_INIT_CODE) {
        
        case FFTL_CONVERT_INIT_FORWARD_FLOAT:
            
            for(int i = 0; i < mN * Cx2Re; i += 8) {
                // d-cache friendly
                // (1, 1, 1, 1)
                xSIMD[0][i+0] = xSIMD[1][i+0] + xSIMD[1][i+2] + xSIMD[1][i+4] + xSIMD[1][i+6]; 
                xSIMD[0][i+1] = xSIMD[1][i+1] + xSIMD[1][i+3] + xSIMD[1][i+5] + xSIMD[1][i+7];
                // (1, -j, -1, j)
                xSIMD[0][i+2] = xSIMD[1][i+0] + xSIMD[1][i+3] - xSIMD[1][i+4] - xSIMD[1][i+7];
                xSIMD[0][i+3] = xSIMD[1][i+1] - xSIMD[1][i+2] - xSIMD[1][i+5] + xSIMD[1][i+6];
                // (1, -1, 1, -1)
                xSIMD[0][i+4] = xSIMD[1][i+0] - xSIMD[1][i+2] + xSIMD[1][i+4] - xSIMD[1][i+6];
                xSIMD[0][i+5] = xSIMD[1][i+1] - xSIMD[1][i+3] + xSIMD[1][i+5] - xSIMD[1][i+7];
                // (1, j, -1, -j)
                xSIMD[0][i+6] = xSIMD[1][i+0] - xSIMD[1][i+3] - xSIMD[1][i+4] + xSIMD[1][i+7];
                xSIMD[0][i+7] = xSIMD[1][i+1] + xSIMD[1][i+2] - xSIMD[1][i+5] - xSIMD[1][i+6];        
            }    
#ifdef DEBUG_TRANS_RADIX_4
            std::cout << "=========" << std::endl;
            for(int i = 0; i < mN * Cx2Re; i++) {
                std::cout << xSIMD[0][i] << std::endl;
            }   
            std::cout << "=========" << std::endl;    
#endif               
            break;
            
        case FFTL_CONVERT_INIT_BACKWARD_FLOAT:
            
            for(int i = 0; i < mN * Cx2Re; i += 8) {
                // d-cache friendly
                // (1, 1, 1, 1)
                xSIMD[0][i+0] = xSIMD[1][i+0] + xSIMD[1][i+2] + xSIMD[1][i+4] + xSIMD[1][i+6]; 
                xSIMD[0][i+1] = xSIMD[1][i+1] + xSIMD[1][i+3] + xSIMD[1][i+5] + xSIMD[1][i+7];
                // (1, j, -1, -j)
                xSIMD[0][i+2] = xSIMD[1][i+0] - xSIMD[1][i+3] - xSIMD[1][i+4] + xSIMD[1][i+7];
                xSIMD[0][i+3] = xSIMD[1][i+1] + xSIMD[1][i+2] - xSIMD[1][i+5] - xSIMD[1][i+6];
                // (1, -1, 1, -1)
                xSIMD[0][i+4] = xSIMD[1][i+0] - xSIMD[1][i+2] + xSIMD[1][i+4] - xSIMD[1][i+6];
                xSIMD[0][i+5] = xSIMD[1][i+1] - xSIMD[1][i+3] + xSIMD[1][i+5] - xSIMD[1][i+7];
                // (1, -j, -1, j)
                xSIMD[0][i+6] = xSIMD[1][i+0] + xSIMD[1][i+3] - xSIMD[1][i+4] - xSIMD[1][i+7];
                xSIMD[0][i+7] = xSIMD[1][i+1] - xSIMD[1][i+2] - xSIMD[1][i+5] + xSIMD[1][i+6];        
            }              
            break;          
            
        case FFTL_CONVERT_INIT_FORWARD_DOUBLE:
            
            for(int i = 0; i < mN * Cx2Re; i += 8) {
                // d-cache friendly
                // (1, 1, 1, 1)
                xSIMDd[0][i+0] = xSIMDd[1][i+0] + xSIMDd[1][i+2] + xSIMDd[1][i+4] + xSIMDd[1][i+6]; 
                xSIMDd[0][i+1] = xSIMDd[1][i+1] + xSIMDd[1][i+3] + xSIMDd[1][i+5] + xSIMDd[1][i+7];
                // (1, -j, -1, j)
                xSIMDd[0][i+2] = xSIMDd[1][i+0] + xSIMDd[1][i+3] - xSIMDd[1][i+4] - xSIMDd[1][i+7];
                xSIMDd[0][i+3] = xSIMDd[1][i+1] - xSIMDd[1][i+2] - xSIMDd[1][i+5] + xSIMDd[1][i+6];
                // (1, -1, 1, -1)
                xSIMDd[0][i+4] = xSIMDd[1][i+0] - xSIMDd[1][i+2] + xSIMDd[1][i+4] - xSIMDd[1][i+6];
                xSIMDd[0][i+5] = xSIMDd[1][i+1] - xSIMDd[1][i+3] + xSIMDd[1][i+5] - xSIMDd[1][i+7];
                // (1, j, -1, -j)
                xSIMDd[0][i+6] = xSIMDd[1][i+0] - xSIMDd[1][i+3] - xSIMDd[1][i+4] + xSIMDd[1][i+7];
                xSIMDd[0][i+7] = xSIMDd[1][i+1] + xSIMDd[1][i+2] - xSIMDd[1][i+5] - xSIMDd[1][i+6];        
            }    
#ifdef DEBUG_TRANS_RADIX_4
            std::cout << "=========" << std::endl;
            for(int i = 0; i < mN * Cx2Re; i++) {
                std::cout << xSIMD[0][i] << std::endl;
            }   
            std::cout << "=========" << std::endl;    
#endif               
            break;
            
        case FFTL_CONVERT_INIT_BACKWARD_DOUBLE:
            
            for(int i = 0; i < mN * Cx2Re; i += 8) {
                // d-cache friendly
                // (1, 1, 1, 1)
                xSIMDd[0][i+0] = xSIMDd[1][i+0] + xSIMDd[1][i+2] + xSIMDd[1][i+4] + xSIMDd[1][i+6]; 
                xSIMDd[0][i+1] = xSIMDd[1][i+1] + xSIMDd[1][i+3] + xSIMDd[1][i+5] + xSIMDd[1][i+7];
                // (1, -j, -1, j)
                xSIMDd[0][i+2] = xSIMDd[1][i+0] - xSIMDd[1][i+3] - xSIMDd[1][i+4] + xSIMDd[1][i+7];
                xSIMDd[0][i+3] = xSIMDd[1][i+1] + xSIMDd[1][i+2] - xSIMDd[1][i+5] - xSIMDd[1][i+6];
                // (1, -1, 1, -1)
                xSIMDd[0][i+4] = xSIMDd[1][i+0] - xSIMDd[1][i+2] + xSIMDd[1][i+4] - xSIMDd[1][i+6];
                xSIMDd[0][i+5] = xSIMDd[1][i+1] - xSIMDd[1][i+3] + xSIMDd[1][i+5] - xSIMDd[1][i+7];
                // (1, j, -1, -j)
                xSIMDd[0][i+6] = xSIMDd[1][i+0] + xSIMDd[1][i+3] - xSIMDd[1][i+4] - xSIMDd[1][i+7];
                xSIMDd[0][i+7] = xSIMDd[1][i+1] - xSIMDd[1][i+2] - xSIMDd[1][i+5] + xSIMDd[1][i+6];        
            }             
            break;
            
        default:
            break;
    }
}







void CTrans4::convert(const float * WnHpQ1, const float * WnHpQ2, const float * WnHpQ3, Complex8f SignI) {
        
#if(0)
    //     2.1 calculate the second conversion pipeline
    //         xSIMD[0] ==//==> xSIMD[1]    

    CoWnQ1.load(WnHpQ1);
    CoWnQ2.load(WnHpQ2);
    CoWnQ3.load(WnHpQ3);    
    
    for(int i = 0; i < mN * Cx2Re; i += 32) {  

        X1Com8.load(xSIMD[0]+i+8);        
        X2Com8.load(xSIMD[0]+i+16);
        X3Com8.load(xSIMD[0]+i+24);        
        
        X0ComIn8.load(xSIMD[0]+i);
        X1ComIn8 = X1Com8 * CoWnQ1;
        X2ComIn8 = X2Com8 * CoWnQ2;
        X3ComIn8 = X3Com8 * CoWnQ3;        
        
        //X1VecIn8  = X1ComIn8.to_vector();
        //X1VecIn8J = permute8f<1,0,3,2,5,4,7,6>(X1VecIn8);
        //X1VecIn8  = change_sign<1,0,1,0,1,0,1,0>(X1VecIn8J);
        //X1ComIn8J = X1VecIn8.to_complex();        

        X1ComIn8J = X1ComIn8 * ComplexJ;
        X3ComIn8J = X3ComIn8 * ComplexJ;
                        
        X0ComUt8 = X0ComIn8 + X1ComIn8  + X2ComIn8 + X3ComIn8;
        X1ComUt8 = X0ComIn8 - X1ComIn8J - X2ComIn8 + X3ComIn8J;
        X2ComUt8 = X0ComIn8 - X1ComIn8  + X2ComIn8 - X3ComIn8;
        X3ComUt8 = X0ComIn8 + X1ComIn8J - X2ComIn8 - X3ComIn8J;
        
        X0ComUt8.store(xSIMD[1]+i);
        X1ComUt8.store(xSIMD[1]+i+8);        
        X2ComUt8.store(xSIMD[1]+i+16);
        X3ComUt8.store(xSIMD[1]+i+24); 
        
#ifdef DEBUG_TRANS_RADIX_4
    std::cout << "====%====" << std::endl;
    for(int i = 0; i < mN * Cx2Re; i++) {
        std::cout << xSIMD[1][i] << std::endl;
    }   
    std::cout << "=========" << std::endl;    
#endif          
    }
#endif
    
    
    
    //     2.2 iterate over the rest of the stages using Complex8f
    //         input: xSIMD[1]
    //           end: xSIMD[Result]

    pin = xSIMD[0];
    pon = xSIMD[1];
    pinpon = NULL;
    Result = 0x0;
    
    CoWnSet = 0;     // initial offset in the Wn heap
    BiN = mN / 16;   // initial number of butter-fly groups

    // loop for butter-fly pipelines

    for(int BuLn = 4; BuLn < mN; BuLn = BuLn << 2) {   

        // loop for butter-fly SIMD groups

        XoeSet  = 0;
        for(int BiP = 0; BiP < BuLn/4; BiP++) {

            // loop for all butter-fly groups

            CoWnQ1.load(WnHpQ1 + CoWnSet);
            CoWnQ2.load(WnHpQ2 + CoWnSet);
            CoWnQ3.load(WnHpQ3 + CoWnSet);  

            for(int Bi = 0; Bi < BiN; Bi++) {        
                
                X1Com8.load(pin + (Bi * 4 + 1) * BuLn * Cx2Re + XoeSet);        
                X2Com8.load(pin + (Bi * 4 + 2) * BuLn * Cx2Re + XoeSet);
                X3Com8.load(pin + (Bi * 4 + 3) * BuLn * Cx2Re + XoeSet);        
                
                X0ComIn8.load(pin + (Bi * 4 + 0) * BuLn * Cx2Re + XoeSet);
                X1ComIn8 = X1Com8 * CoWnQ1;
                X2ComIn8 = X2Com8 * CoWnQ2;
                X3ComIn8 = X3Com8 * CoWnQ3;   
                
                //XoComplex8.load(pin + im(Bi)*BuLn*Cx2Re + XoeSet);
                //XeComplexIn8.load(pin + re(Bi)*BuLn*Cx2Re + XoeSet);
                //XoComplexIn8 = XoComplex8 * CoWn8;
                
                //XeComplexUt8 = XeComplexIn8 + XoComplexIn8;
                //XoComplexUt8 = XeComplexIn8 - XoComplexIn8;
                
                X1ComIn8J = X1ComIn8 * SignI;
                X3ComIn8J = X3ComIn8 * SignI;
                
                X0ComUt8 = X0ComIn8 + X1ComIn8  + X2ComIn8 + X3ComIn8;
                X1ComUt8 = X0ComIn8 - X1ComIn8J - X2ComIn8 + X3ComIn8J;
                X2ComUt8 = X0ComIn8 - X1ComIn8  + X2ComIn8 - X3ComIn8;
                X3ComUt8 = X0ComIn8 + X1ComIn8J - X2ComIn8 - X3ComIn8J;
                
                //XeComplexUt8.store(pon + re(Bi)*BuLn*Cx2Re + XoeSet);
                //XoComplexUt8.store(pon + im(Bi)*BuLn*Cx2Re + XoeSet);

                X0ComUt8.store(pon + (Bi * 4 + 0) * BuLn * Cx2Re + XoeSet);
                X1ComUt8.store(pon + (Bi * 4 + 1) * BuLn * Cx2Re + XoeSet);        
                X2ComUt8.store(pon + (Bi * 4 + 2) * BuLn * Cx2Re + XoeSet);
                X3ComUt8.store(pon + (Bi * 4 + 3) * BuLn * Cx2Re + XoeSet);                  
            }
            CoWnSet += (4*Cx2Re);
            XoeSet  += (4*Cx2Re);
        }
        pinpon = pin;
        pin = pon;
        pon = pinpon;
        BiN = BiN >> 2;
        Result = (Result + 1) & 0x1;        
    }
    //end of conversion    
}













void CTrans4::convert(const double * WnHpQ1d, const double * WnHpQ2d, const double * WnHpQ3d, Complex4d SignI) {

     
    //     2.2 iterate over the rest of the stages using Complex8f
    //         input: xSIMD[0]
    //           end: xSIMD[Result]

    pind = xSIMDd[0];
    pond = xSIMDd[1];
    pinpond = NULL;
    Result = 0x0;
    
    CoWnSet = 0;     // initial offset in the Wn heap
    BiN = mN / 16;   // initial number of butter-fly groups

    // loop for butter-fly pipelines

    for(int BuLn = 4; BuLn < mN; BuLn = BuLn << 2) {   

        // loop for butter-fly SIMD groups

        XoeSet  = 0;
        for(int BiP = 0; BiP < BuLn/2; BiP++) {

            // loop for all butter-fly groups

            CoWnQ1d.load(WnHpQ1d + CoWnSet);
            CoWnQ2d.load(WnHpQ2d + CoWnSet);
            CoWnQ3d.load(WnHpQ3d + CoWnSet);  

            for(int Bi = 0; Bi < BiN; Bi++) {        
                
                X1Com8d.load(pind + (Bi * 4 + 1) * BuLn * Cx2Re + XoeSet);        
                X2Com8d.load(pind + (Bi * 4 + 2) * BuLn * Cx2Re + XoeSet);
                X3Com8d.load(pind + (Bi * 4 + 3) * BuLn * Cx2Re + XoeSet);        
                
                X0ComIn8d.load(pind + (Bi * 4 + 0) * BuLn * Cx2Re + XoeSet);
                X1ComIn8d = X1Com8d * CoWnQ1d;
                X2ComIn8d = X2Com8d * CoWnQ2d;
                X3ComIn8d = X3Com8d * CoWnQ3d;   
        
                
                X1ComIn8Jd = X1ComIn8d * SignI;
                X3ComIn8Jd = X3ComIn8d * SignI;
                
                X0ComUt8d = X0ComIn8d + X1ComIn8d  + X2ComIn8d + X3ComIn8d;
                X1ComUt8d = X0ComIn8d - X1ComIn8Jd - X2ComIn8d + X3ComIn8Jd;
                X2ComUt8d = X0ComIn8d - X1ComIn8d  + X2ComIn8d - X3ComIn8d;
                X3ComUt8d = X0ComIn8d + X1ComIn8Jd - X2ComIn8d - X3ComIn8Jd;
                

                X0ComUt8d.store(pond + (Bi * 4 + 0) * BuLn * Cx2Re + XoeSet);
                X1ComUt8d.store(pond + (Bi * 4 + 1) * BuLn * Cx2Re + XoeSet);        
                X2ComUt8d.store(pond + (Bi * 4 + 2) * BuLn * Cx2Re + XoeSet);
                X3ComUt8d.store(pond + (Bi * 4 + 3) * BuLn * Cx2Re + XoeSet);                  
            }
            CoWnSet += (2*Cx2Re);
            XoeSet  += (2*Cx2Re);
        }
        pinpond = pind;
        pind = pond;
        pond = pinpond;
        BiN = BiN >> 2;
        Result = (Result + 1) & 0x1;        
    }
    //end of conversion    
}









void CTrans4::ForwardFloat(void) {
    convert_init(FFTL_CONVERT_INIT_FORWARD_FLOAT);
    convert(CoWnHeapQ1, CoWnHeapQ2, CoWnHeapQ3, ComplexJ);
}


void CTrans4::InverseFloat(void) {
    convert_init(FFTL_CONVERT_INIT_BACKWARD_FLOAT);
    convert(CoWnHeapIQ1, CoWnHeapIQ2, CoWnHeapIQ3, ComplexJ_);
}


void CTrans4::ForwardDouble(void) {
    convert_init(FFTL_CONVERT_INIT_FORWARD_DOUBLE);
    convert(CoWnHeapQ1d, CoWnHeapQ2d, CoWnHeapQ3d, ComplexJd);    
}

void CTrans4::InverseDouble(void) {
    convert_init(FFTL_CONVERT_INIT_BACKWARD_DOUBLE);
    convert(CoWnHeapIQ1d, CoWnHeapIQ2d, CoWnHeapIQ3d, ComplexJd_);    
}   








//
// setters and getters
// Note: input vector could be destroyed!
//       resultant vector can be modified directly!
//
float * CTrans4::GetInputFloat(void) {
    return xSIMD[0];
}


float * CTrans4::GetOutputFloat(void) {
    return xSIMD[Result];
}



double * CTrans4::GetInputDouble(void) {
    return xSIMDd[0];
}


double * CTrans4::GetOutputDouble(void) {
    return xSIMDd[Result];
}


int CTrans4::GetPoints(void) {
    return mN;
}










CTrans4::CTrans4(const CTrans4& orig) {
}


CTrans4::~CTrans4() {
    
    free(xSIMD[0]);
    free(xSIMD[1]);
    free(xSIMDd[0]);
    free(xSIMDd[1]);    
}

