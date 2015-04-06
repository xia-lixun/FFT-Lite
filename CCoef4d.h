/* 
 * File:   CCoef4d.h
 * Author: lixun
 *
 * Created on 2015年3月16日, 上午12:21
 */

#ifndef CCOEF4D_H
#define	CCOEF4D_H





class CCoef4d {
    
public:
    CCoef4d(int Length);
    virtual ~CCoef4d();
    
    void CoefComp(void);
    void CoefCompI(void);

    const double * GetCoefComp(int Quadrant) const;
    const double * GetCoefCompI(int Quadrant) const;
    int            GetPoints(void) const;        
    
private:
    CCoef4d(const CCoef4d& orig);
    void CoefGen(int Quadrant, double * CoWnHeap, double Math2);    
    
protected:
    double * Wn;
    
    double * CoWnHeap1;
    double * CoWnHeap2;
    double * CoWnHeap3;    

    double * CoWnHeap1I;
    double * CoWnHeap2I;
    double * CoWnHeap3I;    
    
    int     mN;
    int     mQuadPts;  //number of points per quadrant. mN/4
    int     mHeapPts;  //number of points for heaps of all quadrants        
    
};

#endif	/* CCOEF4D_H */

