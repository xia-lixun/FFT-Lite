/* 
 * File:   CCoef4d.h
 * Author: lixun
 *
 * Created on 2015年3月16日, 上午12:21
 */

#ifndef CCOEF4D_H
#define	CCOEF4D_H

#include "CCoef4.h"




class CCoef4d : public CCoef4 {

    
public:
    
    CCoef4d();
    CCoef4d(int Length);
    virtual ~CCoef4d();
    
    void CoefComp(void);
    void CoefCompI(void);

    const double * GetCoefComp(int Quadrant) const;
    const double * GetCoefCompI(int Quadrant) const;
      
    
private:
    
    double * CoWnHeap1;
    double * CoWnHeap2;
    double * CoWnHeap3;    

    double * CoWnHeap1I;
    double * CoWnHeap2I;
    double * CoWnHeap3I; 
    
    CCoef4d(const CCoef4d& orig);
    void CoefGenHeap(double * CoWnHeap, const long double * Wn);      
    
};

#endif	/* CCOEF4D_H */

