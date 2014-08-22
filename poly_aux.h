#ifndef POLY_AUX_H
#define POLY_AUX_H

using namespace std;

/* Constants */
#define M 1000     // i.e.

#define N 10000    // i.e.

#define TWO_EXP8 256 // 2^8
#define TWO_EXP8_M_ONE 255 // 2^8 - 1

//int d -> dimension of hypercube

//int n -> # of nodes in field

//int M -> maximum storage space of each sensor node

unsigned char calcLT(unsigned char d);
unsigned char calcLM(unsigned char d);

#endif
