//
//  MKPS.cpp
//  GFChar
//
//  Created by Cory Pruce on 7/12/13.
//  Copyright (c) 2013 NYIT REU. All rights reserved.
//

#include <cstdlib>                      // rand function
#include <vector>                       // vector class
#include <cmath>                        // pow function
#include "GFChar.cpp"                   // GF(256) operations

//int d -> dimension of hypercube
    
//int n -> # of nodes in field

//int M -> maximum storage space of each sensor node

const int M = 1000;     // i.e.

const int n = 10000;    // i.e.

unsigned char calcLT(unsigned char M, unsigned char d);

unsigned char calcLM(unsigned char d, unsigned char n);

int factorial(int ent);


/* Setup */

class ID{
    
    private:
        unsigned char * ls;
        unsigned char d;
        
    public:
        ID(){};
        ID(unsigned char d){
            (*this).d = d;
            ls = new unsigned char[d]; // d
        };
        
        unsigned char * get() {   
            return (*this).ls;
        }
        unsigned char getSize() {   
            return (*this).d;
        }
        
        ~ID(){ delete this;};
        
};

/*class Var{
    private:
        unsigned char d;
    public:
        Var(){}
        Var(unsigned char d){ (*this).d = d;}
        unsigned char get(){
            return (*this).d;
        }
        ~Var(){delete this;}
};*/


class Uni_var{
    
    private:
        unsigned char * coeffs;
        unsigned char * expns;
        unsigned char dmo;
        unsigned char m;
    public:
        Uni_var(){};
        Uni_var(unsigned char d, unsigned char m, unsigned char L){
            (*this).dmo = d-1;
            coeffs = new unsigned char[L];
            expns = new unsigned char[L];
            (*this).m = m;
        }; 
        Uni_var * jRemove(ID id, unsigned char j, unsigned char m); // remove jth element
        unsigned char * getCoeffs() {   
            return (*this).coeffs;
        }
        unsigned char * getExpns() {   
            return (*this).expns;
        }
        unsigned char getSize() {   
            return (*this).dmo;
        }
        unsigned char getM(){
            return (*this).m;
        }
        void setUniPoly(unsigned char i, unsigned char elem){
            *((*this).expns + i) = elem; 
        }
        ~Uni_var(){ delete this;};
    
};

class D_var{
    
private:
    vector < vector <unsigned char > > expns; 
    vector <unsigned char > coeffs;
    unsigned char d;
    unsigned char m;
    unsigned char L;
    unsigned char t;
public:
    D_var(){};
    D_var(unsigned char d, unsigned char m){
        (*this).d = d;
        (*this).m = m;
        (*this).t = calcLT(M, d);
    }; 
    vector<unsigned char> getCoeffs() {   
        return (*this).coeffs;
    }
    vector< vector<unsigned char> > getExpns() {   
        return (*this).expns;
    }
    unsigned char getDim() {   
        return (*this).d;
    }
    unsigned char getM(){
        return (*this).m;
    }
    unsigned char getLT(){
        return (*this).t;
    }
    ~D_var(){ delete this;};
    
};

class SymKey{
    
    private:
        unsigned char * ls;
        unsigned char d;
        unsigned char m;
    public:
        SymKey(){};
        SymKey(unsigned char d, unsigned char m){
            (*this).d = d;
            ls = new unsigned char[d]; // d
            (*this).m = m;
        };
        SymKey * jRemove(Uni_var id, unsigned char j, unsigned char m);
        
        unsigned char * get() {   
            return (*this).ls;
        }
        unsigned char getSize() {   
            return (*this).d;
        }
        unsigned char getM(){
            return (*this).m;
        }
        void setSymKey(unsigned char i, unsigned char elem){
            *((*this).ls + i) = elem; 
        }
        ~SymKey(){ delete this;};
    
};


/*template<size_t d, typename T>
class hypercube{
    private:
        vector<hypercube<d-1,T> > * hQ;
    public: 
        hypercube(){
            hQ = new vector<hypercube<d-1,T> >();
            int i = 0;
            while(i < 256){
                (*hQ).push_back(hypercube<d-1, T>());
                i++;
            }
        }
};

template<typename T>
class hypercube<1, T>{
    private:
        vector<T> * hQ; 
    public:
        hypercube<0, T>(){
            hQ = new vector<T>();
            int i = 0;
            while(i < 256){
                (*hQ).push_back(static_cast<T>(i));
                i++;
            }
        }
};*/

template<size_t d>
class hypercube{
private:
    vector<hypercube<d-1> > hQ;
    vector<ID> ids;
public: 
    hypercube(){
        if(d != 1){
        hQ = new vector<hypercube<d-1> >(d);
        int i = 0;
        int m = calcLM(d, n);
        double md = pow((double)m, (double)d);
        vector<ID> ids = *new vector<ID>(md);
        while(i < md){

            *(&ids + i) = (*hQ).push_back(hypercube<d-1>()); //*(id + i) = hypercubeAux
            i++;
            }
        }
    }
};

template<size_t d>
class hypercubeAux{
private:
    vector<hypercubeAux<d-1> > hQ;
public: 
    hypercubeAux(){
        if(d != 1){
            hQ = new vector<hypercubeAux<d-1> >(d);
            int i = 0;
            while(i < 256){
                
                (*hQ).push_back(hypercubeAux<d-1>()); //*(id + i) = hypercubeAux
                i++;
            }
        }
    }
};

/*D_var * createD_varSymPolys(unsigned char d, unsigned char m){
    D_var * D_vars = new D_var[d*m];
    unsigned char * cfs;
    int L = factorial(2*d);
    int c = 0;
    unsigned char * eps = new unsigned char[d];
    D_var * temp = new D_var(d, m, L);
    unsigned char t = (*temp).getLT();
    int r = rand() % (t - 1) + 1;                   // 1 <= (exponents = psuedo-random number) < t
    
    for(int v = 0; v < d*m; v++){
    
        for(int i = 0; i < d-1; i++){
            *(eps + i) = r;
            r = rand() % (t - 1) + 1;
        }
    
        *(eps + (d-1)) = t;
    
        do {
            for(int i = 0; i < L; i++){
                
                for(int j = 0; j < d; j++){
                    *((*temp).getExpns()->at(i) + j) = *(eps + j);
                }
            }
            
        
        } while (next_permutation(eps, eps + d - 1));
        
        *(D_vars + c) = *(temp);
        c++;
    }
    
    return D_vars;
}*/

D_var * createD_varSymPolys(unsigned char d, unsigned char m){
    
    D_var * D_vars = new D_var[d*m];
    vector< vector < vector <unsigned char > > > cfs = *new vector< vector < vector <unsigned char> > > ();
    vector< vector< vector < vector <unsigned char > > > > eps = *new vector< vector< vector < vector <unsigned char> > > >();
    D_var temp = *new D_var(d, m);
    unsigned char t = temp.getLT();
    unsigned char alpha = rand() % 255;          // coefficients from GF(256)
    unsigned char beta = rand() % 255;
    int a1 = rand() % (t - 1);                   // 0 <= (exponents = psuedo-random number) < t
    int a2 = rand() % (t - 1);
    int b0 = rand() % (t - 1);
    int b1 = rand() % (t - 1);
    int b2 = rand() % (t - 1);
    
    for(int i = 0; i < d; i++){
        
        eps.push_back(*new vector< vector < vector <unsigned char > > >());
        cfs.push_back(*new vector < vector <unsigned char > >());
        
        for(int j = 0; j < m; j++){
            
            eps.at(i).push_back(*new vector < vector <unsigned char > >());
            cfs.at(i).push_back(*new vector<unsigned char>());
            cfs.at(i).at(j).push_back(alpha);
            cfs.at(i).at(j).push_back(beta);
            eps.at(i).at(j).push_back(*new vector<unsigned char>());
            eps.at(i).at(j).push_back(*new vector<unsigned char>());
            eps.at(i).at(j).at(0).push_back(t);
            eps.at(i).at(j).at(0).push_back(a1);
            eps.at(i).at(j).at(0).push_back(a2);
            eps.at(i).at(j).at(1).push_back(b0);
            eps.at(i).at(j).at(1).push_back(b1);
            eps.at(i).at(j).at(1).push_back(b2);
            alpha = rand() % 255;
            beta = rand() % 255;
            a1 = rand() % (t - 1);                   
            a2 = rand() % (t - 1);
            b0 = rand() % (t - 1);
            b1 = rand() % (t - 1);
            b2 = rand() % (t - 1);
            temp.getCoeffs() = cfs[i][j];
            temp.getExpns() = eps[i][j];
            *(D_vars + i*d + j) = temp;
        }
    }
    return D_vars;
}

/*template<size_t d, typename T>
void buildRecAux(ID id, hypercube<d, T> hq){
    for(int i = 0; i < d; i++){
        for(int j = 0; j < 256; j++){
            buildRecAux(id, count, hq, index);
            
        }
        hq = static_cast<hypercube<d-1, T> >((*hq).hQ);
        
    }
    
}

template<size_t d, typename T>
vector<ID> buildRec(vector<ID> ids, int count, hypercube<d, T> hq, int index){
    for(int i = 0; i < d; i++){
        ID temp;
        for(int j = 0; j < 256; j++){
            buildRecAux(ids, count, hq, index);
            count++;
            index++;
        }
        hq = static_cast<hypercube<d-1, T> >((*hq).hQ);
        
    }

    
    
    return buildRec(ids, count, hq, index);
}

template<size_t d, typename T>
vector<ID> build(vector<ID> ids, hypercube<d, T> hq){
    int m = calcLM(d, n);
    int md = pow(m, d);
    ids = *new vector<ID>(md);
        
    buildRec(ids, 0, hq, 0);
}*/


/* Link-Key Establishment */

Uni_var * jRemove(ID id, unsigned char j, unsigned char m){                          // remove jth element
    unsigned char d = id.getSize();
    int L = 2*factorial(d);
    Uni_var * jR = new Uni_var(d, m, L);                // d - 1
    unsigned char jFlag = 0;
    
    for(unsigned char i = 0; i < d; i++){
        
        if(i != j && jFlag == 0){
            (*jR).setUniPoly(i, *(id.get() + i));     // before jth element is passed
        }
        else if (i != j) {
            (*jR).setUniPoly(i-1, *(id.get() + i));   // after jth element is passed
        }
        else {
            jFlag = 1;                                  // jth element is passed
            j = -1;                                     
        }
        
    }
    
    return jR;
}

SymKey * jRemove(Uni_var id, unsigned char j, unsigned char m){                          // remove jth element
    unsigned char d = id.getSize();
    SymKey * jR = new SymKey(d, m);                // d - 1
    unsigned char jFlag = 0;
    
    for(unsigned char i = 0; i < d; i++){
        if(i != j && jFlag == 0){
            (*jR).setSymKey(i, *(id.getExpns() + i));     // before jth element is passed   OR GETCOEFFSS!!!!
        }
        else if (i != j) {
            (*jR).setSymKey(i-1, *(id.getExpns() + i));   // after jth element is passed
        }
        else {
            jFlag = 1;                                  // jth element is passed
            j = -1;                                     
        }
    }
    
    return jR;
}

D_var * simplify(D_var * lst){
    
    
    
    return lst;
}

/*Uni_var * createKeyRing(ID* id){                        // constructs d-univariate keys for a given ID
    unsigned char d = (*id).getSize();
    Uni_var * ring = new Uni_var[d];
    Uni_var * kk = new Uni_var(d, *(*id).get());
    unsigned char k;

    for(unsigned char i = 0; i < d; i++){
        
        k = *((*id).get() + i);
        
        *kk = *jRemove(id, i, k);   
        
        *(ring + i) = *kk;
        
    }
    
    return ring;
}*/

Uni_var * createKeyRing(ID id, D_var dv[], unsigned char m){                        // constructs d-univariate keys for a given ID
    unsigned char d = id.getSize();
    int L = 2*factorial(d);
    Uni_var * ring = new Uni_var[d];
    Uni_var temp = *new Uni_var(d, *id.get(), L);
    
    
    for(int x = 0; x < d; x++){
        
        for (int y = 0; y < m; y++) {

            int count = 0;
            for(int i = 0; i < 2; i++){
        
                for(int j = 0; j < d-1; j++){
                    temp.getExpns()[count] = dv[(x*d)+y].getExpns()[i][j];
                    for(int k = 1; k < d-1; k++){
                         
                    }
            
                    count++;
            
                }
                
            }
        }
    }
    
    return ring;
}

int hasHamOne(Uni_var A, Uni_var B){                                // returns -1 if doesn't have a Hamming distance of one, else returns
    unsigned char d = A.getSize();                                // the jth element that is the discrepency between the two ID's
    int count = 0;
    int j = 0;
    
    if(d != B.getSize()){
        return -1;
    }
    
    for(int i = 0; i < d; i++){
        if(*(A.getExpns() + i) == *(B.getExpns() + i)){ //or getCoeffs
            j = i;
            //cout << "j = " << static_cast<int>(j) << endl;
            count++;
        }
        if(count > 1){
            return -1;
        }
    }
    if(count == 0){
        return -1;
    }                                                        // count will never be 0 since every ID is unique
    
    return j;
}

SymKey * establishLinkKey(ID A, ID B, D_var * dv, unsigned char m){
    
    unsigned char d = A.getSize();
    
    SymKey * symKeys = new SymKey[d-1];
    
    Uni_var * aRing = createKeyRing(A, dv, m);
    Uni_var * bRing = createKeyRing(B, dv, m);
    
    int common = -1;
    unsigned char sInd = 0;
    
    for(unsigned char i = 0; i < d; i++){
        
        if((common = hasHamOne(aRing[i], bRing[i])) != -1){   // create key
            
            m = (*(aRing + i)).getM();
            
            *(symKeys + sInd) = *new SymKey(d, m);
            
            int j = 0;
            int ind = 0;
            
            while(j < d){
                
                if(j != common){
                    
                    *((symKeys + sInd)->get() + j) = *((*(aRing + i)).get() + ind);
                    *((symKeys + sInd)->get() + j+1) = *((*(bRing + i)).get() + ind);
                    j++;
                    
                }
                else {
            
                    *((symKeys + sInd)->get() + j) = *((*(aRing + i)).get() + ind);
                    
                }
                ind++;
                j++;
            }
        
        }
        
        else {
            sInd--;
        }
        common = -1;
        sInd++;
    }
    
    unsigned char sum;
    for(int i = 0; i < d-1; i++){
        cout << "SymKey[" << i << "]'s m = " << static_cast<int>((*(symKeys + i)).getM()) << endl;
        for(int j = 0; j < d; j++){
            cout << "SymKey[" << i << "][" << j << "] = " << static_cast<int>(*((*(symKeys + i)).get() + j)) << endl;
            sum^=*((*(symKeys + i)).get() + j);
        }
    }
    //cout << static_cast<int>(sum) << endl;
    
    return symKeys;
}

unsigned char calcLM(unsigned char d, unsigned char n){
    return static_cast<int>(ceil(pow(n, 1.0/d)));                           //double m = ceiling of the d route of n
}

unsigned char calcLT(unsigned char M, unsigned char d){                                                   
    return M/(d - 1);                                                       //int t =  floor of M/(d - 1)
}


/* Network Connectivity and Resilience */

// double R = min Communication radius ~= sqrt((ln(n) + c)/(pi * n * probabilityLinkKeyEstablishment)) where c > 0 is a constant

int factorial(int ent){
    int acc = 1;
    
    while(ent > 0){
        acc *= ent;
        ent--;
    }
    
    return acc;
}

double binomialCoefficient(int n, int j){               // n choose j
    return ((double)factorial(n))/(factorial(n - j)*factorial(j));
}

double probabilityLinkKeyEstablishment(double d, double m, double v){
    double theta = 1 - ((double)(1/m));
    return ((d*(m-2)+v)*(pow(m-1, d)))/(n*(n-1)) + (theta*pow(m, d)*(m*(d-1) + v*pow(theta, v-1) + m*(1-d)*pow(theta, v))/(n*(n-1)));
}

/*int sizeOfLargeComponent(){
    int size = 0;
    
    return size;
}*/

double probabilityLinkKeyCompromise(int t, int m){
    double acc = 0.0;
    double pnc = 0.0;
    
    for(int i = 0; i < t; i++){
        acc += binomialCoefficient(m, i)*(pow(pnc, i))*(pow(1 - pnc, m - i)); 
    }
    
    return 1 - acc;
}


//print function

template<class T>
void print(T obj){
    int d = (*obj).getSize();
    unsigned char * l = (*obj).get();
    //cout << "d is " << d << endl;
    for(int i = 0; i < d; i++){
        cout << "The " << i << "th element is " << static_cast<int>(*(l + i)) << endl;
    }
    cout << endl;
}


int main(){
    
    ID * id = new ID(3);
    ID * lp = new ID(3);
    
    /*for(int h = 0; h < (*(*lp).get()).size(); h++){
        (*(*lp).get()).at(h) = h;
    }*/
    
    //*((*lp).get() + 1) = 4;
    //(*(*lp).get()).at(2) = 9;
    
    *((*id).get()) = 3;
    *((*id).get() + 1) = 2;
    *((*id).get() + 2) = 1;
    
    *((*lp).get()) = 3;
    *((*lp).get() + 1) = 4;
    *((*lp).get() + 2) = 1;
    
    print(id);
    
    print(lp);
    
    /*Uni_var * a = new Uni_var(3, 0);
    Uni_var * b = new Uni_var(3, 0);
    
    *((*a).get()) = 3;
    *((*a).get() + 1) = 2;
    
    *((*b).get()) = 3;
    *((*b).get() + 1) = 4;
    
    

    
    cout  << "hasHamOne(a, b) is " << static_cast<int>(hasHamOne(a, b)) << endl;*/
    
    //D_var * jK = D_varove(id, 8);
    
    //print(id);
    
    //print(jK);
    
    //Uni_var * jk = createKeyRing(lp);
    
    //print(jk);
    //cout << static_cast<int>((*jk).getSize()) << endl;
    
    D_var * ds = createD_varSymPolys(3, 3);
    
    for(int i = 0; i < 9; i++){
     
     for(int j = 0; j < 12; j++){
         cout << "ds[" << i << "][" << j << "] = " << static_cast<int>((*(ds + i)).getCoeffs().at(j)) << endl;
         //cout << "ds[" << i << "][" << j << "] = " << static_cast<int>((*(ds + i)).getExpns().at(j)) << endl;
        }
     }
    
    //Uni_var * jp = createKeyRing(lp);
    
    //print(jp);
    
    //establishLinkKey(id, lp);
    
    //cout << "binomialCoefficient(6,2) = " << binomialCoefficient(6,2) << endl;
    
    //cout << endl;
    
    //hypercube<3, unsigned char> * hT = new hypercube<3, unsigned char>();
    
    //matrix<3, vector<unsigned char> > * hQ = new matrix<3, vector<unsigned char> >();
    
    //cout << "sizeof(unsigned char) is " << sizeof(unsigned char) << endl;
    
    //cout << "size of everything is " << sizeof((*id).get()) + (*(*id).get()).capacity() << endl;
    
    return 0;
}
