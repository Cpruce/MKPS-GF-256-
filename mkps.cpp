//
//  MKPS.cpp
//  GFChar
//
//  Created by Cory Pruce on 7/12/13.
//  Copyright (c) 2013 NYIT REU. All rights reserved.
//

#include "mkps.h"

D_var ** createD_varSymPolys(unsigned char d, unsigned char m){
		D_var ** D_vars = new D_var*[d*m];

		vector<unsigned char> cfs;
		vector<unsigned char> eps;

		vector<unsigned char> a_shuffle;
		vector<unsigned char> b_shuffle;

		D_vars[0] = new D_var(d, m);
		unsigned char t = D_vars[0]->getLT();
		delete D_vars[0];
		unsigned char t_m_one = t - 1; // NEED CHECK FOR t = 0
		unsigned char alpha = rand() % TWO_EXP8_M_ONE;                    // coefficients from GF(256)
		unsigned char beta = rand() % TWO_EXP8_M_ONE;
		unsigned char a1 = rand() % t_m_one;                   // 0 <= (exponents = psuedo-random number) < t
		unsigned char a2 = rand() % t_m_one;
		unsigned char b0 = rand() % t_m_one;
		unsigned char b1 = rand() % t_m_one;
		unsigned char b2 = rand() % t_m_one;
		int maj_ind;
		int i;
		int j;

		for(i = 0; i < d; i++){

				maj_ind = i * d; 

				for(j = 0; j < m; j++){ // d*m d-variate polynomials
						
						D_vars[maj_ind + j] = new D_var(d, m);
						D_vars[maj_ind + j]->coeffs->push_back(alpha);
						D_vars[maj_ind + j]->coeffs->push_back(beta);

						// update here
						std::exit(1);
						//for(int k = 0; k < d-1; k++){
								D_vars[maj_ind + j]->expns->push_back(t);
								D_vars[maj_ind + j]->expns->push_back(a1);
								D_vars[maj_ind + j]->expns->push_back(a2);	
						//}
					
						D_vars[maj_ind + j]->expns->push_back(b0);
						D_vars[maj_ind + j]->expns->push_back(b1);	
						D_vars[maj_ind + j]->expns->push_back(b2);	

						alpha = rand() % TWO_EXP8_M_ONE;
						beta = rand() % TWO_EXP8_M_ONE;
						a1 = rand() % t_m_one;                   
						a2 = rand() % t_m_one;
						b0 = rand() % t_m_one;
						b1 = rand() % t_m_one;
						b2 = rand() % t_m_one;



						eps.clear();
						cfs.clear();

				}

		}
		return D_vars; 
}

/* Link-Key Establishment */

Uni_var * simplify(Uni_var * lst){

		int d = lst->getSize();
		int L = 2*factorial(d);
		for(int i = 0; i < d; i++){   
				for (int j = 0; j < L; j++) {
						for(int l = j+1; l < L; l++){
								if(lst[i].getExpns()->at(j) == lst[i].getExpns()->at(l)){
										lst[i].getCoeffs()->at(j) ^= lst[i].getCoeffs()->at(l);
										lst[i].getCoeffs()->at(l) = lst[i].getCoeffs()->at(j);
								}
						}

				}

		}


		return lst;
}

/*void remove (vectover<unsigned char> * ls, vector<unsigned char> lst, int j){
//vector<unsigned char> * ls = new vector<unsigned char>();
for(int i = 0; i < lst.size(); i++){
if(i != j){
ls->push_back(lst.at(i));
}
}

}*/

Uni_var ** createKeyRing(ID id, D_var ** dv, int m){                        // constructs d-univariate keys for a given ID
		cout << "createKeyRing\n" << endl;	
		int d = static_cast<int>(id.getSize());
		Uni_var ** ring = new Uni_var*[d];
		vector<unsigned char> * ex = new vector<unsigned char>();
		unsigned char val;
		int maj_ind;
		
		for(int x = 0; x < d; x++){

				maj_ind = x * d;  // major index
				ring[x] = new Uni_var(d, m);

				for (int y = 0; y < m; y++) { // iterate through the d*m symmetric d-variate polynomials

						for(int i = 0; i < 2; i++){ // alpha beta

								for(int j = 0; j < d; j++){  // d variables


										if(j != x){


												*ex = *dv[maj_ind+y]->getExpns();//.at(i);
												//cout << "size is " << dv[maj_ind+y]->getExpns()->size() <<  endl;
												ex->erase(ex->begin() + j);

												for(int z = 0; z < factorial(d-1)-1; z++){

														for(int q = 0; q < 2; q++){

																next_permutation(ex->begin(), ex->end());
																//cout << "size is " << ex->size() << endl;//dv[maj_ind+y]->getCoeffs().size() << endl;
																val += Product(dv[maj_ind+y]->getCoeffs().at(q), power(*(id.get() + q), ex->at(q)));
														}

												}
												//cout << "ex is " << dv[maj_ind+y] << endl;
												//cout << "size is " << ex->size() << endl;

												ring[x]->setUniPolyCo(val);
												cout << "dv size is " <<dv[maj_ind+y]->getExpns()->size() << endl; 
												//ring[x]->getExpns()->push_back(dv[maj_ind+y]->getExpns()->at(i));
												ring[x]->getExpns()->push_back(ex->at(i));
												
												cout << "j is " << j  << " and ring[" << x << "]->getExpns()->size() is " << ring[x]->getExpns()->size()<< endl;
										}

								}



						}

				}

				ex->clear();  

		}

		//ring = simplify(ring);
		delete ex;
		return ring;
}

int hasHamOne(Uni_var * A, Uni_var * B){                                // returns -1 if doesn't have a Hamming distance of one, else returns
		unsigned char dmo = A->getSize();                                // the jth element that is the discrepency between the two ID's
		int count = 0;
		int j = 0;

		if(dmo != B->getSize()){
				return -1;
		}
		cout << "expns is " << static_cast<int>(B->getExpns()->size()) << " and dmo is " << static_cast<int>(dmo) << endl;

		for(int i = 0; i < std::min((int)dmo, (int)B->getExpns()->size()); i++){
				cout << "A->getExpns()->at(i) is " << static_cast<int>(A->getExpns()->at(i)) << " and for B: " << static_cast<int>(B->getExpns()->at(i)) << endl; 
				if(A->getExpns()->at(i) == B->getExpns()->at(i)){ //or getCoeffs
						j = i;
						count++;
				}
				if(count > 1){
				cout << "wat" << endl;
						return -1;
				}
		}               
		if(count == 0){
				return -1;
		}                                                        // count will never be 0 since every ID is unique

		return j;
}

SymKey ** establishLinkKey(ID A, ID B, D_var ** dv, int m){
		/*		for(int i = 0; i < 9; i++){

				for(int j = 0; j < 2; j++){
				for(int s = 0; s < 3; s++){
		//std::cout << "ds[" << i << "][" << j << "][" << s << "] = " << static_cast<int>(dvs[i].getExpns().at(j).at(s)) << endl;
		std::cout << "ds[" << i << "][" << j << "] = " << static_cast<int>((*(dv + i))->getExpns()->at(j)) << endl;
		}
		}
		} */

		unsigned char d = A.getSize();
		SymKey ** symKeys = new SymKey*[d-1];

		Uni_var ** aRing = createKeyRing(A, dv, m);
		Uni_var ** bRing = createKeyRing(B, dv, m);
		

		
		int common = -1;
		int sInd = 0;
		int j;
		int ind;
		cout << static_cast<int>(d) << endl;
		for(int i = 0; i < d; i++){
				std::cout << "bRing[i]->getM()" << static_cast<int>(bRing[i]->getM()) << endl; 
				if((common = hasHamOne(aRing[i], bRing[i])) != -1){   // create key of d-1 symmetric keys
						cout << "return ? " << endl;
						
						m = aRing[i]->getM();

						symKeys[sInd] = new SymKey(d, m);

						j = 0;
						ind = 0;

						while(j < d){

								if(j != common){
										cout << "ind is " << ind << " and aRing[i]->getSize() is " << static_cast<int>(aRing[i]->getSize()) << endl;
										cout << "symKeys[sInd] size is " << symKeys[sInd]->getExpns()->size() << " and aRing[i]->getCoeffs()->size() is " << aRing[i]->getCoeffs()->size()<< endl;
										symKeys[sInd]->getExpns()->push_back(aRing[i]->getExpns()->at(ind));
										symKeys[sInd]->getCoeffs()->push_back(aRing[i]->getCoeffs()->at(ind));
										symKeys[sInd]->getExpns()->push_back(bRing[i]->getExpns()->at(ind));
										symKeys[sInd]->getCoeffs()->push_back(bRing[i]->getCoeffs()->at(ind));


										j++;

								}
								else {
										cout << "symKeys[sInd] is " << symKeys[sInd]->getExpns()->size() << endl;

										symKeys[sInd]->getExpns()->push_back(aRing[i]->getExpns()->at(ind));
										symKeys[sInd]->getCoeffs()->push_back(aRing[i]->getCoeffs()->at(ind));	
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

		//delete[] aRing;
		//delete[] bRing;		

		//assert that there really are d-1 common keys!

		/*unsigned char sum;
		  for(int i = 0; i < d-1; i++){
		  std::cout << "SymKey[" << i << "]'s m = " << static_cast<int>((*(symKeys + i)).getM()) << endl;
		  for(int j = 0; j < d; j++){
		  std::cout << "SymKey[" << i << "][" << j << "] = " << static_cast<int>((*(symKeys + i)).getCoeffs().at(j)) << endl;
		  sum^=(*(symKeys + i)).getCoeffs().at(j);
		  }
		  }*/
		return symKeys;
}

unsigned char calcLM(unsigned char d){ // took out unsigned char n, check if ok
		return static_cast<int>(ceil(pow(N, 1.0/d)));                           //double m = ceiling of the d route of n
}

unsigned char calcLT(unsigned char d){                                                   
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
		int n_m_one = N - 1;
		return ((d*(m-2)+v)*(pow(m-1, d)))/(N*n_m_one) + (theta*pow(m, d)*(m*(d-1) + v*pow(theta, v-1) + m*(1-d)*pow(theta, v))/(N*n_m_one));
}

double probabilityLinkKeyCompromise(int t, int m){
		double acc = 0.0;
		double pnc = 0.0;

		for(int i = 0; i < t; i++){
				acc += binomialCoefficient(m, i)*(pow(pnc, i))*(pow(1 - pnc, m - i)); 
		}

		return 1 - acc;
}

//print functions

void printID(ID * id){
		int d = id->getSize();
		unsigned char * l = id->get();
		//std::cout << "d is " << d << endl;
		for(int i = 0; i < d; i++){
				std::cout << "The " << i << "th element is " << static_cast<int>(l[i]) << endl;
		}
		std::cout << endl;
}

void printD_var(D_var * dv){}

void printUni_var(Uni_var * uv){}

void printSymKey(SymKey * sk){


}

int main(){
		FillLogArrays();
		ID * id = new ID(3);
		ID * lp = new ID(3);

		/*for(int h = 0; h < (*(*lp).get()).size(); h++){
		 * (*(*lp).get()).at(h) = h;
		 * }*/

		//*((*lp).get() + 1) = 4;
		//(*(*lp).get()).at(2) = 9;

		*(id->get()) = 3;
		*(id->get() + 1) = 2;
		*(id->get() + 2) = 1;

		*(lp->get()) = 3;
		*(lp->get() + 1) = 4;
		*(lp->get() + 2) = 1;

		//print(id);

		//print(lp);
		////
		/*Uni_var * a = new Uni_var(3, 0);
		  Uni_var * b = new Uni_var(3, 0);

		 *((*a).get()) = 3;
		 *((*a).get() + 1) = 2;

		 *((*b).get()) = 3;
		 *((*b).get() + 1) = 4;




		 std::cout  << "hasHamOne(a, b) is " << static_cast<int>(hasHamOne(a, b)) << endl;*/

		//D_var * jK = D_varove(id, 8);

		//print(id);

		//print(jK);

		//Uni_var * jk = createKeyRing(lp);

		//print(jk);
		//std::cout << static_cast<int>((*jk).getSize()) << endl;
		printID(id);
		printID(lp);		
		D_var ** dvs = createD_varSymPolys(3, 3);
		//Uni_var * kks = createKeyRing(*id, ds, 3);

		/*		for(int i = 0; i < 9; i++){

				for(int j = 0; j < 2; j++){
				for(int s = 0; s < 3; s++){
		//std::cout << "ds[" << i << "][" << j << "][" << s << "] = " << static_cast<int>(dvs[i].getExpns().at(j).at(s)) << endl;
		std::cout << "ds[" << i << "][" << j << "] = " << static_cast<int>((*(dvs + i))->getExpns()->at(j)) << endl;
		}
		}
		}*/


		/*for(int i = 0; i < 9; i++){

		  for(int j = 0; j < 2; j++){
		//for(int s = 0; s < 3; s++){
		std::cout << "ds[" << i << "][" << j << "] = " << static_cast<int>((*(dvs + i)).getCoeffs().at(j)) << endl;
		//std::cout << "ds[" << i << "][" << j << "] = " << static_cast<int>((*(ds + i)).getExpns().at(j)) << endl;
		//}
		}
		}
		*/

		//std::cout <<"id size = " << (*id).getSize() << endl;

		//Uni_var * rg = createKeyRing(*id, dvs, 3);

		//std::cout << static_cast<int>((rg + 0)->getCoeffs().at(0)) << endl;

		/*for(int i = 0; i < 3; i++){
		  for(int j = 0; j < rg->getExpns().size(); j++){
		  std::cout << "rg[" << i << "].expnAt(" << j << ") is " << static_cast<int>((rg + i)->getExpns().at(j)) << endl; 
		  }
		  for (int k = 0; k < rg->getCoeffs().size(); k++) {
		  std::cout << "rg[" << i << "].coefAt(" << k << ") is " << static_cast<int>((rg + i)->getCoeffs().at(k)) << endl;
		  }
		  }*/

		//print(jp);

		SymKey ** sks = establishLinkKey(*id, *lp, dvs, 3);
		// expns seem correct with 3 (d) * 2 (d-1) * 2 (alpha beta)


		for(int i = 0; i < 2; i++){
						cout << "huh" << endl;
						cout << "size is " <<sks[i]->getExpns()->size() << endl; 
				for(int j = 0; j < sks[i]->getExpns()->size(); j++){
						std::cout << "sks[" << i << "].expnAt(" << j << ") is " << static_cast<int>(sks[i]->getExpns()->at(j)) << endl; 
				}
				for (int k = 0; k < sks[i]->getCoeffs()->size(); k++) {
						std::cout << "sks[" << i << "].coefAt(" << k << ") is " << static_cast<int>(sks[i]->getCoeffs()->at(k)) << endl;
				}
		}


		//std::cout << "binomialCoefficient(6,2) = " << binomialCoefficient(6,2) << endl;

		//std::cout << endl;

		//hypercube<3, unsigned char> * hT = new hypercube<3, unsigned char>();

		//matrix<3, vector<unsigned char> > * hQ = new matrix<3, vector<unsigned char> >();

		//std::cout << "sizeof(unsigned char) is " << sizeof(unsigned char) << endl;

		//std::cout << "size of everything is " << sizeof((*id).get()) + (*(*id).get()).capacity() << endl;

		//return
}
