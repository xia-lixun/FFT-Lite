/* 
 * File:   main.cpp
 * Author: root
 *
 * Created on January 10, 2015, 10:19 PM
 */

#include "instrset.h"

//#define _USE_MATH_DEFINES
#include <new>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
#include "vectorclass.h"
#include "complexvec.h"



#define  N          (8192)                 //N must be greater than or equal to 4
#define  Cx2Re      (2)                 //Conversion of Index from real to complex domain
#define  re(x)     ((x) * Cx2Re + 0)    //real index of the complex number indexed by x
#define  im(x)     ((x) * Cx2Re + 1)    //real index of the complex number indexed by x




using namespace std;



//unsigned int iLog2(unsigned int value)
//{
//    unsigned int l = 0;
//    while( (value >> l) > 1 ) ++l;
//    return l;
//}




/*  
 *  Complex Number FFT DIT
 *  
 * 
 *           N-1
 *         +----+ 
 *          \         -i 2Pi k n/N
 *  X(k) =   +  x(n) e   
 *          /    
 *         +----+
 *           n=0
 * 
 * 
 *  Decimate x(n) with even and odd numbers:
 * 
 *  ======================
 *  = For k = [0..N/2-1] =
 *  ======================
 * 
 *   N/2-1                                    N/2-1
 *  +----+                                   +----+
 *   \         -i 2Pi k p/(N/2)   -i 2Pi k/N  \           -i 2Pi k p/(N/2)
 *    +  x(2p)e                + e             +  x(2p+1)e                
 *   /                                        /
 *  +----+                                   +----+
 *    p=0                                      p=0
 * 
 * 
 * 
 *  ==========================
 *  = For k = [0..N/2-1]+N/2 =
 *  ==========================
 * 
 *   N/2-1                                    N/2-1
 *  +----+                                   +----+
 *   \         -i 2Pi k p/(N/2)   -i 2Pi k/N  \           -i 2Pi k p/(N/2)
 *    +  x(2p)e                - e             +  x(2p+1)e                
 *   /                                        /
 *  +----+                                   +----+
 *    p=0                                      p=0
 */

int main(int argc, char** argv) {
    
    float * xSIMD[2];
    xSIMD[0] = new float[N*2];
    xSIMD[1] = new float[N*2]; 

    for(int i = 0; i < N*2; i++) {
        xSIMD[0][i] = (float)i;        
    }    

    
    
    //  preparation phase:
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
    //generate Wn(0) Wn(1) ... Wn(N/2-1)
    float * Wn = new float[N];
    for(int i = 0; i < N; i+=2) {
        Wn[i]   = 0.0f;
        Wn[i+1] = (float)(i/2);
    }
    
    Vec4f coeff(-2.0f * M_PI / (float)N);
    Vec4f a;
    Vec4f b;
    for(int i = 0; i < N; i+=4) {
        a.load(Wn+i);
        b = a * coeff;
        b.store(Wn+i);
    }
    
    Complex4f ax;
    Complex4f bx;
    for(int i = 0; i < N; i+=4) {
        ax.load(Wn+i);
        bx = cexp(ax);
        bx.store(Wn+i);
    }
    
    //map Wn to the heap
    float * CoWnHeap = new float[(N-2)*Cx2Re];
    int WnStride  = N/4;                      //use N/4 instead of N/2 because the optimization of Wn heap
    int WnElements = 2;                       //the heap starts from the second conversion stage  
    int WnCnt = 0;                            
    
    while (WnStride > 0) {
        for(int i = 0; i < WnElements; i++) {
            CoWnHeap[WnCnt++] = Wn[re(WnStride * i)];
            CoWnHeap[WnCnt++] = Wn[im(WnStride * i)];            
        }
        WnElements = WnElements * 2;
        WnStride = WnStride / 2;
    }
    //
    //
    //end of twiddle factor calculation
    
    
   
    
    //generate the bit-reverse map for input complex array
    int * BitReverseTreeEven   = new int[N/4];
    int * BitReverseMapComplex = new int[N*2];
    
    
    //even number index generated with growing tree method:
    //N=4    {0}
    //N=8    {0} {N/4}
    //N=16   {0   N/4} {0+N/8   N/4+N/8}
    //N=32   {0   N/4   0+N/8   N/4+N/8} {0+N/16   N/4+N/16   0+N/8+N/16   N/4+N/8+N/16}
    //N=...   ...
    int Nd = N/4;
    BitReverseTreeEven[0] = 0;
    for(int i = 1; i <= N/8; i <<= 1) {
        //i is the same as the segment length
        for(int k = 0; k < i; k++) {
            BitReverseTreeEven[i + k] = BitReverseTreeEven[k] + Nd;
        }
        Nd = Nd/2;
    }
    //generate global bit reverse map
    for(int i = 0; i < N/4; i++) {
        BitReverseMapComplex[(i*2) * 2 + 0] = BitReverseTreeEven[i] * 2 + 0;
        BitReverseMapComplex[(i*2) * 2 + 1] = BitReverseTreeEven[i] * 2 + 1;
        
        BitReverseMapComplex[(i*2+1) * 2 + 0] = (BitReverseTreeEven[i] + N/2) * 2 + 0;
        BitReverseMapComplex[(i*2+1) * 2 + 1] = (BitReverseTreeEven[i] + N/2) * 2 + 1;
    }
    for(int i = 0; i < N; i++) {
        BitReverseMapComplex[i+N] = BitReverseMapComplex[i] + 2;
    }
    
    
    
    timespec TimeZero, TimeOne;
    long Elapsed;
    
    clock_gettime(CLOCK_MONOTONIC_RAW, &TimeZero);
    
    //execution stage.
    //step 1: re-map the input array with bit-reverse pattern
    //        xSIMD[0] ==/bit reverse indexed/==> xSIMD[1]
    Vec4i ReMapRoute;
    Vec4f ReMapValue;
    
    for(int i = 0; i < N*2; i += 4) {
        ReMapRoute.load(BitReverseMapComplex + i);
        ReMapValue = lookup<N*2>(ReMapRoute, xSIMD[0]);
        ReMapValue.store(xSIMD[1] + i);
    }

    //step 2: iterate over butter-fly patterns
    
    //step 2.0: calculate the first conversion
    //          xSIMD[1] ==//==> xSIMD[0]
    Vec4f XeValueIn;
    Vec4f XoValueIn;
    Vec4i XeRoute(0, 1, 4, 5);
    Vec4i XoRoute(2, 3, 6, 7);
    Vec4f XeValueUt;
    Vec4f XoValueUt;
    Vec4f XeValue;
    Vec4f XoValue;    
    
    for(int i = 0; i < N*2; i += 8) {
        XeValueIn = lookup<8>(XeRoute, xSIMD[1] + i);
        XoValueIn = lookup<8>(XoRoute, xSIMD[1] + i);
        
        XeValueUt = XeValueIn + XoValueIn;
        XoValueUt = XeValueIn - XoValueIn;
        
        XeValue = blend4f<0, 1, 4, 5>(XeValueUt, XoValueUt);
        XoValue = blend4f<2, 3, 6, 7>(XeValueUt, XoValueUt);
        
        XeValue.store(xSIMD[0] + i);
        XoValue.store(xSIMD[0] + i + 4);    
    }
    
    //step 2.1: calculate the second conversion
    //          xSIMD[0] ==//==> xSIMD[1]
    //
    //apply Wn coefficients to the input
    Complex4f CoWn;
    Complex4f XoComplexIn;
    Complex4f XeComplexIn;
    Complex4f XoComplexUt;
    Complex4f XeComplexUt;
    Complex4f XoComplex;
    Complex4f XeComplex;    
    
    
    int XoeSet;
    int CoWnSet;
    int BiN;
    float * pin = xSIMD[0];
    float * pon = xSIMD[1];
    float * pinpon = NULL;
    unsigned int Result = 0x0;
    //       Butter-fly group: repeating butter-fly structure.
    //Butter-fly group length: number of different twiddle coefficients
    //                         at every conversion stage.
    
    
//int BuLn = 2;                                     
//    CoWn.load(CoWnHeap);
//    for(int Bi = 0; Bi < N/4; Bi++) {            //butter-fly group index
//        
//          XoComplex.load(xSIMD[0] + im(Bi)*BuLn*Cx2Re);        
//        XeComplexIn.load(xSIMD[0] + re(Bi)*BuLn*Cx2Re);
//        XoComplexIn = XoComplex * CoWn;
//        
//        XeComplexUt = XeComplexIn + XoComplexIn;
//        XoComplexUt = XeComplexIn - XoComplexIn;
//        
//        XeComplexUt.store(xSIMD[1] + re(Bi)*BuLn*Cx2Re);
//        XoComplexUt.store(xSIMD[1] + im(Bi)*BuLn*Cx2Re);
//    }
    
//step 2.2: calculate the third conversion
//          xSIMD[1] ==//==> xSIMD[0]
//
//BuLn    = (4);
    CoWnSet = 0;
    BiN = N/4;
    for(int BuLn = 2; BuLn < N; BuLn = BuLn << 1) {

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
    
//    //portion II
//    CoWn.load(CoWnHeap + (2)*Cx2Re + (2*Cx2Re));
//    for(int Bi = 0; Bi < N/8; Bi++) {
//        
//        XoComplex.load(xSIMD[1] + (Bi*2+1)*BuLn*Cx2Re + (2*Cx2Re));   
//        XeComplexIn.load(xSIMD[1] + (Bi*2+0)*BuLn*Cx2Re + (2*Cx2Re));
//        XoComplexIn = XoComplex * CoWn;
//        
//        XeComplexUt = XeComplexIn + XoComplexIn;
//        XoComplexUt = XeComplexIn - XoComplexIn;
//        
//        XeComplexUt.store(xSIMD[0] + (Bi*2+0)*BuLn*Cx2Re + (2*Cx2Re));
//        XoComplexUt.store(xSIMD[0] + (Bi*2+1)*BuLn*Cx2Re + (2*Cx2Re));                
//    }

    
    clock_gettime(CLOCK_MONOTONIC_RAW, &TimeOne);
    Elapsed = (TimeOne.tv_sec - TimeZero.tv_sec) * 1000000000 + (TimeOne.tv_nsec - TimeZero.tv_nsec);
    
    for(int i = 0; i < N; i++) {
        cout << xSIMD[Result][i*2+0] << " " << xSIMD[Result][i*2+1] << endl;        
    }
    cout << "Time Elapsed:" << Elapsed << " ns" << endl; 
    
    
    delete BitReverseMapComplex;
    delete BitReverseTreeEven;
    delete Wn;
    delete xSIMD[1];
    delete xSIMD[0];    
    return 0;
}



