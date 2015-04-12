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
    
    CCoef();
    CCoef(int Length);
    virtual ~CCoef();
    
    virtual void CoefComp(void);
    virtual void CoefCompI(void);
    
    const float * GetCoefComp(void) const;
    const float * GetCoefCompI(void) const;
    
    int GetPoints(void) const;
    
    
private:

    float * CoWnHeap;    
    float * CoWnHeapI;
    
    CCoef(const CCoef& orig);
    void CoefGenHeap(float * CoWnHp, const long double * Wn);
    
    
protected:
    
    int mN;    

    void CoefGenWnForward(long double * Wn);
    void CoefGenWnBackward(long double * Wn);
    
};

#endif	/* CTWIDDLEFACTOR_H */

