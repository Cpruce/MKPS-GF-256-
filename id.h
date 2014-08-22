#ifndef ID_H
#define ID_H

class ID{
    
private:
    unsigned char * ls;
    int d;
    
public:
    ID(){};
    ID(int dim){
        d = dim;
        ls = new unsigned char[d];
    };
    
    unsigned char * get() {   
        return ls;
    }
    int getSize() {   
        return d;
    }
    
    ~ID(){ 
        delete[] ls;
    };
    
};


#endif
