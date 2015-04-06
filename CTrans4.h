/* 
 * File:   CTransCompRadix4.h
 * Author: lixun
 *
 * Created on 2015年3月12日, 上午12:34
 */

#ifndef CTRANSCOMPRADIX4_H
#define	CTRANSCOMPRADIX4_H



#include "instrset.h"
#include "vectorclass.h"
#include "complexvec.h"

#include "CSort4.h"
#include "CCoef4.h"
#include "CCoef4d.h"





class CTrans4 {
    
    
public:
    CTrans4(const CSort4& Map, const CCoef4& Wn);
    CTrans4(const CSort4& Map, const CCoef4d& Wn);
    virtual ~CTrans4();
    
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
    CTrans4(const CTrans4& orig);
    CTrans4& operator=(const CTrans4& assign);
    
    
    void convert_init(int FFTL_CONVERT_INIT_CODE);
    void convert(const float * WnHpQ1, const float * WnHpQ2, const float * WnHpQ3, Complex8f SignI);
    void convert(const double * WnHpQ1d, const double * WnHpQ2d, const double * WnHpQ3d, Complex4d SignI);    
    
    
protected: 
    const int * SortLutComp;
    
    /////////////////////////////////////////////////////////////////
    //                       IEEE754 FLOAT                         //
    /////////////////////////////////////////////////////////////////
    float * xSIMD[2];
    double * xSIMDd[2];
    
    const float * CoWnHeapQ1;
    const float * CoWnHeapQ2;
    const float * CoWnHeapQ3;
    
    const float * CoWnHeapIQ1;
    const float * CoWnHeapIQ2;
    const float * CoWnHeapIQ3;


    const double * CoWnHeapQ1d;
    const double * CoWnHeapQ2d;
    const double * CoWnHeapQ3d;
    
    const double * CoWnHeapIQ1d;
    const double * CoWnHeapIQ2d;
    const double * CoWnHeapIQ3d;
    
    
    // working SIMD vectors
        
    Complex8f CoWnQ1;          Complex4d CoWnQ1d;
    Complex8f CoWnQ2;          Complex4d CoWnQ2d;
    Complex8f CoWnQ3;          Complex4d CoWnQ3d;    
    
    Complex8f X0ComIn8;        Complex4d X0ComIn8d;
    Complex8f X1ComIn8;        Complex4d X1ComIn8d;
    Complex8f X2ComIn8;        Complex4d X2ComIn8d;
    Complex8f X3ComIn8;        Complex4d X3ComIn8d;
    
    Complex8f X0ComUt8;        Complex4d X0ComUt8d;
    Complex8f X1ComUt8;        Complex4d X1ComUt8d;
    Complex8f X2ComUt8;        Complex4d X2ComUt8d;
    Complex8f X3ComUt8;        Complex4d X3ComUt8d;
    
    Complex8f X1Com8;          Complex4d X1Com8d;
    Complex8f X2Com8;          Complex4d X2Com8d;
    Complex8f X3Com8;          Complex4d X3Com8d;
    
    //Vec8f     X1VecIn8;
    //Vec8f     X1VecIn8J;
    
    Complex8f ComplexJ;        Complex4d ComplexJd;
    Complex8f ComplexJ_;       Complex4d ComplexJd_;
    Complex8f X1ComIn8J;       Complex4d X1ComIn8Jd;  
    Complex8f X3ComIn8J;       Complex4d X3ComIn8Jd;
   
    int XoeSet;
    int CoWnSet;
    int BiN;
    float * pin;               double * pind;                         
    float * pon;               double * pond;           
    float * pinpon;            double * pinpond;  
    unsigned int Result;
    int mN;    
    
};

#endif	/* CTRANSCOMPRADIX4_H */

