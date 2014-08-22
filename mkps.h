//	mkps.h
//	GfMKPS


#ifndef MKPS_H
#define MKPS_H

#include <cstdlib>                      // rand function
#include <vector>                       // vector class
#include <cmath>                        // pow function
#include <algorithm>
#include <iostream>
#include "gfchar.h"                   // GF(256) operations
#include "id.h"
#include "uni_var.h"
#include "d_var.h"
#include "sym_key.h"
#include "hypercube.h"

D_var * createD_varSymPolys(unsigned char d, unsigned char m);  
Uni_var * simplify(Uni_var * lst);
vector<unsigned char> * remove (vector<unsigned char> lst, int j);
Uni_var * createKeyRing(ID id, D_var * dv, int m);
int hasHamOne(Uni_var A, Uni_var B);
SymKey * establishLinkKey(ID A, ID B, D_var * dv, int m);
int factorial(int ent);
double binomialCoefficient(int n, int j);
double probabilityLinkKeyEstablishment(double d, double m, double v);
double probabilityLinkKeyCompromise(int t, int m);

#endif 
