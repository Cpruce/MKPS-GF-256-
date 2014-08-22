#ifndef SYM_KEY_H
#define SYM_KEY_H

#include "poly_aux.h"

class SymKey{
    
private:
    vector<unsigned char> *expns;
    vector<unsigned char> *coeffs;
    int d;		// see if can be changed to a char
    int m;
public:
    SymKey(){};
    SymKey(int dim, int mNum){
        d = dim;
        m = mNum;
        expns = new vector<unsigned char>();
        coeffs = new vector<unsigned char>();
    };
    vector<unsigned char> *getExpns() {   
        return expns;
    }
    vector<unsigned char> *getCoeffs() {   
        return coeffs;
    }
    int getSize() {   
        return d;
    }
    int getM(){
        return m;
    }
    ~SymKey(){ 
        delete coeffs;
        delete expns;
    };
    
};

#endif
