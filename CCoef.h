/* 
 * File:   CTwiddleFactor.h
 * Author: lixun
 *
 * Created on January 24, 2015, 11:10 AM
 */

#ifndef CTWIDDLEFACTOR_H
#define	CTWIDDLEFACTOR_H

class CCoef {
    
public:
    CCoef(int Length);
    virtual ~CCoef();
    
    void CoefComp(void);
    void CoefCompI(void);
    
    const float * GetCoefComp(void) const;
    const float * GetCoefCompI(void) const;
    int           GetPoints(void) const;
    
private:
    CCoef(const CCoef& orig);
    void CoefGen(float * CoWnHp, double Math2);
    
protected:
    double * Wn;
    float  * CoWnHeap;    
    float  * CoWnHeapI;
    int      mN;
};

#endif	/* CTWIDDLEFACTOR_H */

