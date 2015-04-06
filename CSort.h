/* 
 * File:   CBitReverse.h
 * Author: lixun
 *
 * Created on January 24, 2015, 11:11 AM
 */

#ifndef CBITREVERSE_H
#define	CBITREVERSE_H

class CSort {
    
public:
    CSort(int Length);
    virtual ~CSort();
    
    void LutComp(void);
    void LutReal(void);
    const int * GetLutComp(void) const;
    const int * GetLutReal(void) const;
          int   GetPoints(void) const;
    
private:
    CSort(const CSort& orig);
    
protected:
    void TreeEven(void);
    
    int * mTreeEven;
    int * mLutComp;
    int * mLutReal;
    int   mN;
};

#endif	/* CBITREVERSE_H */

