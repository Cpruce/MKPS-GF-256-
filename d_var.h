#ifndef D_VAR_H
#define D_VAR_H

#include "poly_aux.h"

class D_var{
    
private:
    vector < vector <unsigned char > > *expns; 
    vector <unsigned char > *coeffs;
    unsigned char d;
    unsigned char m;
    unsigned char L;
    unsigned char t;
public:
    D_var(){};
    D_var(unsigned char dim, unsigned char mNum){
        d = dim;
        m = mNum;
        t = calcLT(d);
        expns = new vector<vector<unsigned char> > ();
        coeffs = new vector<unsigned char> ();
    }; 
    vector<unsigned char> *getCoeffs() {   
        return coeffs;
    }
    vector< vector<unsigned char> > *getExpns() {   
        return expns;
    }
    unsigned char getDim() {   
        return d;
    }
    unsigned char getLM(){
        return m;
    }
    unsigned char getLT(){
        return t;
    }
	// make sure elem was malloc'd
    void setDPolyEx(vector< vector<unsigned char> > *elem){
        delete expns;
		expns = elem;
    }
    void setDPolyCo(vector<unsigned char> *elem){
        delete coeffs;
		coeffs = elem;
    }
    ~D_var(){ 
        delete coeffs;
        delete expns;
    };
    
};

#endif
