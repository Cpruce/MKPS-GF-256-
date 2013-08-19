//
//  GFChar.cpp
//  GFChar
//
//  Created by Cory Pruce on 7/10/13.
//  Copyright (c) 2013 NYIT REU. All rights reserved.
//  Code modified from "Arithmetic Operations in a Power-of-Two Galois Field" by AIM Inc.

#include <iostream>

using namespace std;

#define GF 256 // define the Size & Prime Polynomial of this Galois field 
#define PP 285 

unsigned char Log[GF], ALog[GF]; // establish global Log and Antilog arrays 
        // fill the Log[] and ALog[] arrays with appropriate integer values 


void FillLogArrays () { 
    int i; 
    Log[0] = GF - 1;                        //initialization
    ALog[0] = 1; 
    for (i=1; i<GF; i++) { 
        if (ALog[i-1] >= 128) { 
            
            ALog[i] = (ALog[i-1] * 2) ^ PP; // (2^i) mod PP, where (2^i) > 255
            
        }
        else {
            ALog[i] = ALog[i-1] * 2;        // (2^i) mod PP
        }
        
        Log[ALog[i]] = i;                   // corresponding exponent
    } 
} 

unsigned char Product (unsigned char A, unsigned char B) { 
    if ((A == 0) || (B == 0)) return (0); 
    else return (ALog[(Log[A] + Log[B]) % (GF-1)]); 
}

unsigned char Quotient (unsigned char A, unsigned char B) {  
    if (B == 0) return (1-GF); 
    else if (A == 0) return (0); 
    else return (ALog[(Log[A] - Log[B] + (GF-1)) % (GF-1)]); 
} 

unsigned char inverse(unsigned char y)
{
    if (y == 0) return -1;
    return Quotient(1, y);
}
unsigned char power (unsigned char a, unsigned char b){
    if ((a == 1) || (b == 0)) return (1);
    return ALog[Product(b, Log[a])];
}
/*int main(){
    FillLogArrays();
    //cout << static_cast<int>(Product(Product(8, 8), 8)) << endl;
    cout << static_cast<int>(Product(Product(32, 32), 32)) << endl;
    cout << static_cast<int>(power(32, 3)) << endl;
    
}*/
