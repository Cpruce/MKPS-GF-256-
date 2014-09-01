#ifndef D_VAR_H
#define D_VAR_H

#include "poly_aux.h"

class D_var{
    
private:
    unsigned char d;
    unsigned char m;
    unsigned char L;
    unsigned char t;
public:
    vector < vector <unsigned char > > expns; 
    vector <unsigned char > coeffs;
	D_var(){};
    D_var(unsigned char dim, unsigned char mNum){
        d = dim;
        m = mNum;
        t = calcLT(d);		
        //expns = new vector<vector<unsigned char> > ();
        //coeffs = new vector<unsigned char> (2); //alpha beta
    }; 
    vector<unsigned char> getCoeffs() {   
        return coeffs;
    }
    vector< vector<unsigned char> > getExpns() {   
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
    void setDPolyEx(vector< vector<unsigned char> > elem){
		//expns = elem;
    	int elem_len = elem.size();
		vector<unsigned char> temp;
		int in_arr_len;
	
		expns.clear();
		expns.resize(elem_len);

		for(int i = 0; i < elem_len; i++){
			
			temp = elem.at(i);

			in_arr_len = temp.size();
			
			expns.at(i).resize(in_arr_len);

			for(int j = 0; j < in_arr_len; j++){
				
				expns.at(i).push_back(temp.at(j));

			}
		}
	}
    void setDPolyCo(vector<unsigned char> elem){
		//coeffs = elem;
		
		int elem_len = elem.size();

		coeffs.clear();
		coeffs.resize(elem_len);

		for(int i = 0; i < elem_len; i++){
			coeffs.push_back(elem.at(i));
		}
    }
    ~D_var(){ 
        delete &coeffs;
        delete &expns;
    };
    
};

#endif
