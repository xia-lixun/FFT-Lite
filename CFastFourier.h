/* 
 * File:   CFastFourier.h
 * Author: lixun
 *
 * Fast Fourier transform with length 2 to the power of N.
 * 
 * Created on 2015年4月1日, 上午12:09
 */

#ifndef CFASTFOURIER_H
#define	CFASTFOURIER_H










class CFastFourier {
    
    
public:
    CFastFourier();
    virtual ~CFastFourier();
    
    void fft_f(int Points, float * In, float * Out);
    void ifft_f(int Points, float * In, float * Out);
    void fft_d(int Points, double * In, double * Out);
    void ifft_d(int Points, double * In, double * Out);    

private:
    CFastFourier(const CFastFourier& orig);
    bool IsRadix4(int Length);
    
    
protected:

    
};

#endif	/* CFASTFOURIER_H */

