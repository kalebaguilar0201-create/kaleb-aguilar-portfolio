#ifndef INTEGRACIONNUMERICA_H
#define INTEGRACIONNUMERICA_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"

// C_ij
static const long double C[][5] = {
		{2.0, 0, 0, 0, 0},
		{1.0, 1.0,0,0,0},
		{0.5555555556, 0.8888888889, 0.5555555556, 0, 0},
		{0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451, 0},
		{0.2369268850, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268850}
	};
	
	
// R_ij Raices
static const long double R[][5] = { 
		{0, 0, 0, 0, 0},
		{0.5773502692, -0.5773502692,0,0,0},
		{0.7745966692, 0, -0.7745966692,0,0},
		{0.8611363116, 0.3399810436, -0.3399810436, -0.8611363116, 0},
		{0.9061798459, 0.5384693101, 0.0, -0.5384693101, -0.9061798459}
	};


long double Newton_Cotes_Abierto(int n,Vector *X,long double (*f) (long double),long double a, long double b){
	long double out;
	long double h = (b - a) / (n+2);
	switch (n){
		case 0:
			out = 2*h*f(X->v[0]);
			break;
		case 1:
			out = 3*h/2*(f(X->v[0]) + f(X->v[1]));
			break;
		case 2:
			out = 4*h/3*(2*f(X->v[0]) - f(X->v[1]) + 2*f(X->v[2]));
			break;
		case 3:
			out = 5*h/24*(11*f(X->v[0]) + f(X->v[1]) + f(X->v[2]) + 11*f(X->v[3]));
			break;
		default:
			out = 0;
	}
	return out;
}
long double Newton_Cotes_Cerrado(int n,Vector *X,long double (*f) (long double),long double a, long double b){
	long double out;
	long double h = (b - a) / n;
	switch (n){
		case 1:
			out = h/2*(f(X->v[0]) + f(X->v[1]));
			break;
		case 2:
			out = h/3 * (f(X->v[0]) + 4*f(X->v[1]) + f(X->v[2]));
			break;
		case 3:
			out = 3*h/8*(f(X->v[0]) + 3*f(X->v[1]) + 3*f(X->v[2]) + f(X->v[3]));
			break;
		case 4:
			out = 2*h/45*(7*f(X->v[0]) + 32*f(X->v[1]) + 12*f(X->v[2]) + 32*f(X->v[3]) + 5*f(X->v[4]));
			break;
		default:
			out = 0;
	}
	return out;
}

// \sum c_i*f(x_i)
long double Cuadrática_Gaussiana(int n,Vector *X,long double (*f) (long double), long double a, long double b){
	long double out = 0;
	for(int i = 0; i < n ; i++){
		out += C[n-1][i] * f(0.5*(b - a) * R[n-1][i] + 0.5 * (b + a));
	}
	out *= (b-a);
	out /= 2;
	return out;
}


#endif
