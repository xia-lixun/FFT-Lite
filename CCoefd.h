/* 
 * File:   CTwiddleFactorDP.h
 * Author: lixun
 *
 * Created on January 30, 2015, 7:57 PM
 */

#ifndef CTWIDDLEFACTORDP_H
#define	CTWIDDLEFACTORDP_H

#include "CCoef.h"



class CCoefd : public CCoef {
    
public:
    CCoefd();
    CCoefd(int Length);
    virtual ~CCoefd();
    
    void CoefComp(void);
    void CoefCompI(void);    
    
    const double * GetCoefComp(void) const;
    const double * GetCoefCompI(void) const;

    
private:

    double * CoWnHeap;    
    double * CoWnHeapI;    
    
    CCoefd(const CCoefd& orig);
    void CoefGenHeap(double * CoWnHp, const long double * Wn);


};

#endif	/* CTWIDDLEFACTORDP_H */

