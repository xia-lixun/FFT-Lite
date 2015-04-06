/* 
 * File:   CFastFourierTransform.h
 * Author: lixun
 *
 * Created on January 21, 2015, 11:51 PM
 */

#ifndef CFOURIERTRANSCOMPLEX_H
#define	CFOURIERTRANSCOMPLEX_H



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
#include "instrset.h"
#include "vectorclass.h"
#include "complexvec.h"

#include "CSort.h"
#include "CCoef.h"
#include "CCoefd.h"




class CTrans {
    
    
public:
    CTrans(const CSort& Map, const CCoef& Wn);
    CTrans(const CSort& Map, const CCoefd& Wn);
    virtual ~CTrans();
    
    //void Forward4F(void);  //deprecated
    void ForwardFloat(void);
    void InverseFloat(void);
    
    void ForwardDouble(void);
    void InverseDouble(void);
    
    float * GetInputFloat(void);
    float * GetOutputFloat(void);
    
    double * GetInputDouble(void);
    double * GetOutputDouble(void);    

    // Some resources cannot or should not be copied, such as file handles or 
    // mutexes. In that case, simply declare the copy constructor and copy 
    // assignment operator as private without giving a definition:

    int GetPoints(void);
    
    
private:
    CTrans(const CTrans& orig);
    CTrans& operator=(const CTrans& assign);
    
    void convert(const float * CoWnHp);
    void convert(const double * CoWnHpd);    
    
protected:
    const int * SortLutComp;
    
    /////////////////////////////////////////////////////////////////
    //                       IEEE754 FLOAT                         //
    /////////////////////////////////////////////////////////////////
    float * xSIMD[2];
    const float * CoWnHeap;
    const float * CoWnHeapI;

    /////////////////////////////////////////////////////////////////
    //                       IEEE754 DOUBLE                        //
    /////////////////////////////////////////////////////////////////
    double * xSIMDd[2];
    const double * CoWnHeapd;
    const double * CoWnHeapId;    

    // working SIMD vectors

#if(0)
    // 128-bit category    
    int RouteLUT[8];
    
    Vec4i ReMapRoute;        
    Vec4f ReMapValue;        
    
    Vec4f XeValueIn;         
    Vec4f XoValueIn;         
    Vec4i XeRoute;     
    Vec4i XoRoute;     
    Vec4f XeValueUt;         
    Vec4f XoValueUt;         
    Vec4f XeValue;     
    Vec4f XoValue;   
#endif
    
    Complex4f CoWn;          
    Complex4f XoComplexIn;   
    Complex4f XeComplexIn;   
    Complex4f XoComplexUt;   
    Complex4f XeComplexUt;   
    Complex4f XoComplex;     
    Complex4f XeComplex;     

  
    //256-bit category
    
    //Vec8i ReMapRoute8;
    //Vec8f ReMapValue8;   
    
    Complex2f XeValueInC2f;    Complex2d XeValueInC2d;    
    Complex2f XoValueInC2f;    Complex2d XoValueInC2d;    
    Complex2f XeValueUtC2f;    Complex2d XeValueUtC2d;
    Complex2f XoValueUtC2f;    Complex2d XoValueUtC2d;
    
    Complex8f CoWn8;           Complex4d CoWn4d;    
    Complex8f XoComplexIn8;    Complex4d XoComplexInd;
    Complex8f XeComplexIn8;    Complex4d XeComplexInd;
    Complex8f XoComplexUt8;    Complex4d XoComplexUtd;
    Complex8f XeComplexUt8;    Complex4d XeComplexUtd;
    Complex8f XoComplex8;      Complex4d XoComplexd;
    Complex8f XeComplex8;      Complex4d XeComplexd;
    
    int XoeSet;
    int CoWnSet;
    int BiN;
    float * pin;                                        double * pind;
    float * pon;                                        double * pond;
    float * pinpon;                                     double * pinpond;
    unsigned int Result;
    int mN;
};

#endif	/* CFASTFOURIERTRANSFORM_H */

