#ifndef UNI_VAR_H
#define UNI_VAR_H

//using namespace std;

/* Constants */
//#define M 1000     // i.e.

//#define n 10000    // i.e.

//#define two_exp8 256 // 2^8
//#define two_exp8_m_one 255 // 2^8 - 1

//int d -> dimension of hypercube

//int n -> # of nodes in field

//int M -> maximum storage space of each sensor node

#include "poly_aux.h"

class Uni_var{
    
private:
    vector<unsigned char> *coeffs;
    vector<unsigned char> *expns;
    unsigned char dmo;
    unsigned char m;
public:
    Uni_var(){};
    Uni_var(unsigned char d, unsigned char mNum){
        dmo = d-1;
        m = mNum;
        coeffs = new vector<unsigned char>();
        expns = new vector<unsigned char>();
    }; 
    vector<unsigned char> *getCoeffs() {   
        return coeffs;
    }
    vector<unsigned char> *getExpns() {   
        return expns;
    }
    unsigned char getSize() {   
        return dmo;
    }
    unsigned char getM(){
        return m;
    }
    void setUniPolyEx(unsigned char elem){
        expns->push_back(elem);
    }
    void setUniPolyCo(unsigned char elem){
        coeffs->push_back(elem);
    }
    ~Uni_var(){ 
        delete coeffs;
        delete expns;
    };
    
};

#endif
