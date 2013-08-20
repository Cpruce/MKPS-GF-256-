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
    ID(unsigned char dim){
        d = dim;
        ls = new unsigned char[d];
    };
    
    unsigned char * get() {   
        return ls;
    }
    unsigned char getSize() {   
        return d;
    }
    
    ~ID(){ 
        delete[] ls;
    };
    
};

class Uni_var{
    
private:
    vector<unsigned char> coeffs;
    vector<unsigned char> expns;
    unsigned char dmo;
    unsigned char m;
public:
    Uni_var(){};
    Uni_var(unsigned char d, unsigned char mNum){
        dmo = d-1;
        m = mNum;
        coeffs = *new vector<unsigned char>();
        expns = *new vector<unsigned char>();
    }; 
    vector<unsigned char> getCoeffs() {   
        return coeffs;
    }
    vector<unsigned char> getExpns() {   
        return expns;
    }
    unsigned char getSize() {   
        return dmo;
    }
    unsigned char getM(){
        return m;
    }
    void setUniPolyEx(unsigned char elem){
        expns.push_back(elem);
    }
    void setUniPolyCo(unsigned char elem){
        coeffs.push_back(elem);
    }
    ~Uni_var(){ 
        delete &coeffs;
        delete &expns;
    };
    
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
    D_var(unsigned char dim, unsigned char mNum){
        d = dim;
        m = mNum;
        t = calcLT(M, d);
        expns = *new vector<vector<unsigned char> > ();
        coeffs = *new vector<unsigned char> ();
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
    unsigned char getM(){
        return m;
    }
    unsigned char getLT(){
        return t;
    }
    void setDPolyEx(vector< vector<unsigned char> > elem){
        expns = elem;
    }
    void setDPolyCo(vector<unsigned char> elem){
        coeffs = elem;
    }
    ~D_var(){ 
        delete &coeffs;
        delete &expns;
        };
    
};

class SymKey{
    
private:
    vector<unsigned char> expns;
    vector<unsigned char> coeffs;
    unsigned char d;
    unsigned char m;
public:
    SymKey(){};
    SymKey(unsigned char dim, unsigned char mNum){
        d = dim;
        m = mNum;
        expns = *new vector<unsigned char>();
        coeffs = *new vector<unsigned char>();
    };
    vector<unsigned char> getExpns() {   
        return expns;
    }
    vector<unsigned char> getCoeffs() {   
        return coeffs;
    }
    unsigned char getSize() {   
        return d;
    }
    unsigned char getM(){
        return m;
    }
    ~SymKey(){ 
        delete &coeffs;
        delete &expns;
        };
    
};

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
                
                *(&ids + i) = (*hQ).push_back(hypercube<d-1>());
                i++;
            }
        }
    }
    ~hypercube(){delete &hQ;}
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
                
                (*hQ).push_back(hypercubeAux<d-1>());
                i++;
            }
        }
    }
    ~hypercubeAux(){delete &hQ;}
};


D_var * createD_varSymPolys(unsigned char d, unsigned char m){
    
    D_var * D_vars = new D_var[d*m];
    vector< vector < vector <unsigned char > > > cfs = *new vector< vector < vector <unsigned char> > > ();
    vector< vector< vector < vector <unsigned char > > > > eps = *new vector< vector< vector < vector <unsigned char> > > >();
    D_var * temp = new D_var(d, m);
    unsigned char t = (*temp).getLT();
    
    unsigned char alpha = rand() % 255;                    // coefficients from GF(256)
    unsigned char beta = rand() % 255;
    unsigned char a1 = rand() % (t - 1);                   // 0 <= (exponents = psuedo-random number) < t
    unsigned char a2 = rand() % (t - 1);
    unsigned char b0 = rand() % (t - 1);
    unsigned char b1 = rand() % (t - 1);
    unsigned char b2 = rand() % (t - 1);
    
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
            
            (*temp).setDPolyCo(cfs[i][j]);
            (*temp).setDPolyEx(eps[i][j]);
            (*(D_vars + i*d + j)).setDPolyCo((*temp).getCoeffs());
            (*(D_vars + i*d + j)).setDPolyEx((*temp).getExpns());
            
        }
        
        
        
        
        
        
    }

    return D_vars; 
}

/* Link-Key Establishment */

Uni_var * simplify(Uni_var * lst){
    
    int d = (*lst).getSize();
    int L = 2*factorial(d);
    for(int i = 0; i < d; i++){   
        for (int j = 0; j < L; j++) {
            for(int l = j+1; l < L; l++){
                if((*(lst + i)).getExpns().at(j) == (*(lst + i)).getExpns().at(l)){
                    (*(lst + i)).getCoeffs().at(j) ^= (*(lst + i)).getCoeffs().at(l);
                    (*(lst + i)).getCoeffs().at(l) = (*(lst + i)).getCoeffs().at(j);
                }
            }
            
        }
        
    }
     
    
    return lst;
}

vector<unsigned char> remove (vector<unsigned char> lst, int j){
    vector<unsigned char> n = *new vector<unsigned char>();
    for(int i = 0; i < lst.size(); i++){
        if(i != j){
            n.push_back(lst.at(i));
        }
    }
    return n;
}

Uni_var * createKeyRing(ID id, D_var * dv, int m){                        // constructs d-univariate keys for a given ID
    int d = static_cast<int>(id.getSize());
    //int L = 2*factorial(d);
    Uni_var * ring = new Uni_var[d];
    vector<unsigned char> ex = *new vector<unsigned char>();
    unsigned char val;

    for(int x = 0; x < d; x++){
        
        for (int y = 0; y < m; y++) {
            
            for(int i = 0; i < 2; i++){
                
                for(int j = 0; j < d; j++){  
                    
                    if(j != x){
                        
                        ex = remove(dv[(x*d)+y].getExpns()[i], j);
                        
                        for(int z = 0; z < factorial(d-1)-1; z++){
                                
                                for(int q = 0; q < 2; q++){
                                    next_permutation(ex.begin(), ex.end());
                                    val += Product(dv[(x*d)+y].getCoeffs()[q], power(*(id.get() + q), ex.at(q)));
                                    
                                }
                            
                        }
                        
                        (*(ring + x)).setUniPolyCo(val);
                        
                        (*(ring + x)).setUniPolyEx(dv[(x*d)+y].getExpns()[i][j]);
                    
                        }
                    
                    }

                        
                }
            
            }
        
        ex.clear();  
    
    }
    
    //ring = simplify(ring);
                       
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
        if(A.getExpns().at(i) == B.getExpns().at(i)){ //or getCoeffs
            j = i;
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
                       
SymKey * establishLinkKey(ID A, ID B, D_var * dv, int m){
    
        unsigned char d = A.getSize();
                           
        SymKey * symKeys = new SymKey[d-1];
    
        Uni_var * aRing = createKeyRing(A, dv, m);
        Uni_var * bRing = createKeyRing(B, dv, m);
            
        int common = -1;
        unsigned char sInd = 0;
                           
        for(unsigned char i = 0; i < d; i++){
                               cout << "here" << endl; 
            if((common = hasHamOne(aRing[i], bRing[i])) != -1){   // create key
                                  
                m = (*(aRing + i)).getM();
                                   
                *(symKeys + sInd) = *new SymKey(d, m);
                                   
                int j = 0;
                int ind = 0;
                                   
                while(j < d){
                                       
                    if(j != common){
                                           
                        (symKeys + sInd)->getExpns().at(j) = (*(aRing + i)).getExpns().at(ind);
                        (symKeys + sInd)->getCoeffs().at(j) = (*(aRing + i)).getCoeffs().at(ind);
                        (symKeys + sInd)->getExpns().at(j+1) = (*(bRing + i)).getExpns().at(ind);
                        (symKeys + sInd)->getCoeffs().at(j+1) = (*(bRing + i)).getCoeffs().at(ind);
                        j++;
                                           
                    }
                    else {
                                    
                        (symKeys + sInd)->getExpns().at(j) = (*(aRing + i)).getExpns().at(ind);
                        (symKeys + sInd)->getCoeffs().at(j) = (*(aRing + i)).getCoeffs().at(ind);
                                           
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
                           
        /*unsigned char sum;
        for(int i = 0; i < d-1; i++){
            cout << "SymKey[" << i << "]'s m = " << static_cast<int>((*(symKeys + i)).getM()) << endl;
            for(int j = 0; j < d; j++){
                cout << "SymKey[" << i << "][" << j << "] = " << static_cast<int>((*(symKeys + i)).getCoeffs().at(j)) << endl;
                sum^=(*(symKeys + i)).getCoeffs().at(j);
            }
        }*/
                           
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

double probabilityLinkKeyCompromise(int t, int m){
    double acc = 0.0;
    double pnc = 0.0;
                           
    for(int i = 0; i < t; i++){
        acc += binomialCoefficient(m, i)*(pow(pnc, i))*(pow(1 - pnc, m - i)); 
    }
                           
    return 1 - acc;
}
                       
//print function
                       
/*template<class T>
void print(T obj){
    int d = (*obj).getSize();
    vector<unsigned char> l = (*obj).get();
    //cout << "d is " << d << endl;
    for(int i = 0; i < d; i++){
        cout << "The " << i << "th element is " << static_cast<int>(l.at(i)) << endl;
    }
    cout << endl;
}*/
                       
                       
int main(){
                           FillLogArrays();
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
                           
                           //print(id);
                           
                           //print(lp);
                           
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
                           
                           D_var * dvs = createD_varSymPolys(3, 3);
                           
                           //Uni_var * kks = createKeyRing(*id, ds, 3);
                           
                           /*for(int i = 0; i < 9; i++){
                               
                               for(int j = 0; j < 2; j++){
                                   for(int s = 0; s < 3; s++){
                                   cout << "ds[" << i << "][" << j << "][" << s << "] = " << static_cast<int>((*(dvs + i)).getExpns().at(j).at(s)) << endl;
                                   //cout << "ds[" << i << "][" << j << "] = " << static_cast<int>((*(ds + i)).getExpns().at(j)) << endl;
                                   }
                               }
                           }*/
    /*for(int i = 0; i < 9; i++){
        
        for(int j = 0; j < 2; j++){
            //for(int s = 0; s < 3; s++){
                cout << "ds[" << i << "][" << j << "] = " << static_cast<int>((*(dvs + i)).getCoeffs().at(j)) << endl;
                //cout << "ds[" << i << "][" << j << "] = " << static_cast<int>((*(ds + i)).getExpns().at(j)) << endl;
            //}
        }
    }
    */
    
                           //cout <<"id size = " << (*id).getSize() << endl;
                           
                           //Uni_var * rg = createKeyRing(*id, dvs, 3);
    
    //cout << static_cast<int>((rg + 0)->getCoeffs().at(0)) << endl;
    
    /*for(int i = 0; i < 3; i++){
        for(int j = 0; j < rg->getExpns().size(); j++){
            cout << "rg[" << i << "].expnAt(" << j << ") is " << static_cast<int>((rg + i)->getExpns().at(j)) << endl; 
        }
        for (int k = 0; k < rg->getCoeffs().size(); k++) {
            cout << "rg[" << i << "].coefAt(" << k << ") is " << static_cast<int>((rg + i)->getCoeffs().at(k)) << endl;
        }
    }*/
                           
                           //print(jp);
                           
                           SymKey * sks = establishLinkKey(*id, *lp, dvs, 3);
                            
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < sks->getExpns().size(); j++){
            cout << "sks[" << i << "].expnAt(" << j << ") is " << static_cast<int>((sks + i)->getExpns().at(j)) << endl; 
        }
        for (int k = 0; k < sks->getCoeffs().size(); k++) {
            cout << "sks[" << i << "].coefAt(" << k << ") is " << static_cast<int>((sks + i)->getCoeffs().at(k)) << endl;
        }
    }
                
                           
                           //cout << "binomialCoefficient(6,2) = " << binomialCoefficient(6,2) << endl;
                           
                           //cout << endl;
                           
                           //hypercube<3, unsigned char> * hT = new hypercube<3, unsigned char>();
                           
                           //matrix<3, vector<unsigned char> > * hQ = new matrix<3, vector<unsigned char> >();
                           
                           //cout << "sizeof(unsigned char) is " << sizeof(unsigned char) << endl;
                           
                           //cout << "size of everything is " << sizeof((*id).get()) + (*(*id).get()).capacity() << endl;
                           
                           //return 0;
                       }
