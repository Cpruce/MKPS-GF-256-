//
//  GFChar.h
//  GFChar
//
#ifndef GFCHAR_H
#define GFCHAR_H

#define GF 256 // define the Size & Prime Polynomial of this Galois field 
#define PP 285 

void FillLogArrays ();
unsigned char Product (unsigned char A, unsigned char B); 
unsigned char Quotient (unsigned char A, unsigned char B);
unsigned char inverse(unsigned char y);
unsigned char power (unsigned char a, unsigned char b);

#endif
