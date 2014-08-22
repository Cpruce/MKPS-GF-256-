#ifndef HYPERCUBE_H
#define HYPERCUBE_H

#include "poly_aux.h"

template<size_t d>
class hypercube{
	private:
    	vector<hypercube<d-1> > *hQ;
    	vector<ID> *ids;
	public: 
    	hypercube(){
        	if(d != 1){
            	hQ = new vector<hypercube<d-1> >(d);
            	int i = 0;
            	int m = calcLM(d);
            	double md = pow((double)m, (double)d);
            	ids = new vector<ID>(md);
	    		while(i < md){
                
                	//*(&ids + i) = (*hQ).push_back(hypercube<d-1>());
                	ids[i] = hQ->push_back(hypercube<d-1>());
					i++;
            	}
        	}
    	}
    	~hypercube(){
			delete hQ;
			delete ids;	
		}
};

template<size_t d>
class hypercubeAux{
	private:
    	vector<hypercubeAux<d-1> > *hQ;
	public: 
    	hypercubeAux(){
        	if(d != 1){
            	hQ = new vector<hypercubeAux<d-1> >(d);
            	int i = 0;
            	while(i < TWO_EXP8){
                
                	hQ->push_back(hypercubeAux<d-1>());
                	i++;
            	}
        	}
    	}
    	~hypercubeAux(){delete hQ;}
};

#endif
