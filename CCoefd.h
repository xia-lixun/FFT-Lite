/* 
 * File:   CTwiddleFactorDP.h
 * Author: lixun
 *
 * Created on January 30, 2015, 7:57 PM
 */

#ifndef CTWIDDLEFACTORDP_H
#define	CTWIDDLEFACTORDP_H


class CCoefd {
    
public:
    CCoefd(int Length);
    virtual ~CCoefd();
    void CoefComp(void);
    void CoefCompI(void);    
    const double * GetCoefComp(void) const;
    const double * GetCoefCompI(void) const;
    int            GetPoints(void) const;
    
private:
    CCoefd(const CCoefd& orig);
    void CoefGen(double * CoWnHp, double Math2);
    
protected:
    double * Wn;
    double * CoWnHeap;    
    double * CoWnHeapI;    
    int      mN;
};

#endif	/* CTWIDDLEFACTORDP_H */

