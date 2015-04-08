/* 
 * File:   CFourier.h
 * Author: lixun
 *
 * if length is 2 to the power of N, invoke CFastFourier directly
 * if not (arbitrary length), use CFastFourier based on Bluestein algorithm.
 * 
 * Created on 2015年4月1日, 上午12:14
 */

#ifndef CFOURIER_H
#define	CFOURIER_H



class Cdft {
    
public:
    Cdft();
    virtual ~Cdft();
    
    void dft_f(int Points, float * In, float * Out);
    void idft_f(int Points, float * In, float * Out);
    void dft_d(int Points, double * In, double * Out);
    void idft_d(int Points, double * In, double * Out);     
    
private:
    Cdft(const Cdft& orig);
    
protected:
    bool IsRadix2(int Length);
    
};

#endif	/* CFOURIER_H */

