//  Created by Kaleb Aguilar on 13/08/25.
//
#ifndef SOLVE_BY_H
#define SOLVE_BY_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Matrix.h"

int Factorize_by_Crouth(Matrix *A, Vector *b,Matrix *L, Matrix *U){
    if (A->X != A->Y || A->X != L->X || A->Y != L->Y || A->X != U->X || A->Y != U->Y) {
        printf("Dimensiones no válidas\n");
        return -1;
    }
    double sum,sum2;
    
    //PASO 1
    L->m[0][0] = A->m[0][0];
    U->m[0][0] = 1;
    if (L->m[0][0]==0){
		//printf("Factorización Imposible A[0][0] = 0\n");
		return -1;
	} 
	
	//PASO 2
    for (int j = 1; j < A->Y ; j++){ // PRIMER RENGLÓN Y COLUMNA
		U->m[0][j] = A->m[0][j]/L->m[0][0];
		L->m[j][0] = A->m[j][0]/U->m[0][0];
	}
	//PASO 3
    for (int i = 1; i < A->Y - 1; i++){
		//PASO 4
        U->m[i][i] = 1;
        sum = 0;
        for(int k = 0 ; k < i ; k++){
			sum += L->m[i][k]*U->m[k][i];
		}
		L->m[i][i] = A->m[i][i] - sum;
		if (L->m[i][i] == 0){
			//printf("Factorización Imposible L[%d][%d] = 0\n",i,i);
			return -1;
		}
		//PASO 4
		for(int j = i + 1 ; j < A->Y ; j++){
			sum = 0;
			sum2 = 0;
			for(int k = 0 ; k < i ; k++){
				sum += L->m[i][k]*U->m[k][j];
				sum2 += L->m[j][k]*U->m[k][i];
			}
			U->m[i][j] = (A->m[i][j] - sum) / L->m[i][i];
			L->m[j][i] = (A->m[j][i] - sum2) / U->m[i][i];
			
		}
    }
    //PASO 5
    U->m[A->Y-1][A->Y-1] = 1;
	sum = 0;
	for(int k = 0 ; k < A->Y - 1; k++){
		sum += L->m[A->Y-1][k]*U->m[k][A->Y-1];
	}
	L->m[A->Y-1][A->Y-1] = A->m[A->Y-1][A->Y-1] - sum;
    return 0;
}
Vector solve_by_Crouth(Matrix *A, Vector *b){
	Matrix L = initialize_Matrix(A->Y,A->X);
	Matrix U = initialize_Matrix(A->Y,A->X);
	if(Factorize_by_Crouth(A, b, &L, &U) != 0){
        free_Matrix(&L);
        free_Matrix(&U);
        return initialize_Vector(0); 
    }
	
	Vector y = solve_L(&L,b);
	Vector x = solve_U(&U,&y);
	//print_Vector(y);
	//print_Vector(x);
	free_Matrix(&L);
	free_Matrix(&U);
	free_Vector(&y);
	return x;
}

int Factorize_by_Cholensky(Matrix *A, Vector *b,Matrix *L, Matrix *LT){
	
	int n = A->Y;
	double sum;
	
	//PASO1
	if (A->m[0][0] < 0){
		//printf("Factorización Imposible A[0][0] < 0\n");
		return -1;
	}
	
	L->m[0][0] = sqrt(A->m[0][0]);
	LT->m[0][0] = L->m[0][0];
	
	//PASO2
	for(int j = 1 ; j < n ; j++){
		L->m[j][0] = A->m[j][0] / L->m[0][0];
		LT->m[0][j] = L->m[j][0];
	}
	
	//PASO3
	for(int i = 1 ; i < n - 1 ; i++){
		//PASO4
		sum = 0;
		for(int k = 0; k < i ; k++){
			sum += L->m[i][k]*L->m[i][k];
		}
		sum = A->m[i][i] - sum;
		if(sum < 0){
			//printf("Factorización Imposible A[%d][%d] < 0\n",i,i);
			return -1;
		}
		L->m[i][i] = sqrt(sum);
		LT->m[i][i] = L->m[i][i];
		//PASO5
		for(int j = i + 1 ; j < n ; j++){
			sum = 0;
			for(int k = 0; k < i ; k++){
				sum += L->m[j][k]*L->m[i][k];
			}
			L->m[j][i] = (A->m[j][i] - sum) / L->m[i][i];
			LT->m[i][j] = L->m[j][i];	
		}	
		
		//PASO6
		sum = 0;
		for(int k = 0 ; k < n - 1; k++){
			sum += L->m[n-1][k]*L->m[n-1][k] ; 
		}
		sum = A->m[n-1][n-1] - sum;
		if (sum < 0){
			//printf("Factorización Imposible A[%d][%d] < 0\n",n-1,n-1);
			return-1;
		}
		L->m[n-1][n-1] = sqrt(sum);
		LT->m[n-1][n-1] = L->m[n-1][n-1];
	}
	return 0;
}
Vector solve_by_Cholensky(Matrix *A, Vector *b){
	Matrix L = initialize_Matrix(A->Y,A->X);
	Matrix LT = initialize_Matrix(A->Y,A->X);
	Factorize_by_Cholensky(A,b ,&L, &LT);if(Factorize_by_Cholensky(A, b, &L, &LT) != 0){
        free_Matrix(&L);
        free_Matrix(&LT);
        return initialize_Vector(0); 
    }
	
	Vector y = solve_L(&L,b);
	Vector x = solve_U(&LT,&y);
	//print_Vector(y);
	//print_Vector(x);
	free_Matrix(&L);
	free_Matrix(&LT);
	free_Vector(&y);
	return x;
}

Vector solve_by_Jacobi(Matrix *A ,Vector *b){
	int n = b->n;
	double sum, error,norm_x;
	double norma = 1;
	double toler = 0.00000001;
	
	Vector x = initialize_Vector(n);
	Vector x_old = initialize_Vector(n);
	
	//PIVOTEO EN CASO DE A_ii = 0
	for(int i = 0; i < n ; i++){
		if(A->m[i][i] == 0){
			pivote(A,b,i);
		}
	}
	
	// X_old_i = bi/ Aii
	for(int i = 0; i < n ; i++){
		x_old.v[i] = b->v[i] / A->m[i][i];
	}
	
	//RECURSION
	while(norma > toler){
		for(int i = 0; i < n ; i++){
			sum = 0;
			for(int j = 0 ; j < n; j++){
				if(i != j){
					sum += A->m[i][j]*x_old.v[j];
				}
			}
			x.v[i] = (b->v[i] - sum) / A->m[i][i] ;
		}
		error = 0;
		norm_x = 0;
		for(int i = 0; i < n ; i++){
			error += (x.v[i] - x_old.v[i])*(x.v[i] - x_old.v[i]);
			norm_x += x.v[i]*x.v[i];
		}
		norma = sqrt(error/norm_x);
		for(int i = 0; i < n ; i++){
			x_old.v[i] = x.v[i];
		}
	}
	free_Vector(&x_old);
	return x;
}

Vector solve_by_Gauss_Seidel(Matrix *A ,Vector *b){
	int n = b->n;
	double sum, error,norm_x;
	double norma = 1;
	double toler = 0.00000001;
	
	//PIVOTEO EN CASO DE A_ii = 0 
	Vector x = initialize_Vector(n);
	Vector x_old = initialize_Vector(n);
	for(int i = 0; i < n ; i++){
		if(A->m[i][i] == 0){
			pivote(A,b,i);
		}
	}
	// X_old_i = bi/ Aii
	for(int i = 0; i < n ; i++){
		x_old.v[i] = b->v[i] / A->m[i][i];
	}
	//RECURSION
	while(norma > toler){
		for(int i = 0; i < n ; i++){
			sum = 0;
			for(int j = 0 ; j < i; j++){
				sum += A->m[i][j]*x.v[j];
			}
			for(int j = i + 1 ; j < n; j++){
				sum += A->m[i][j]*x_old.v[j];
			}
			x.v[i] = (b->v[i] - sum) / A->m[i][i] ;
		}
		error = 0;
		norm_x = 0;
		for(int i = 0; i < n ; i++){
			error += (x.v[i] - x_old.v[i])*(x.v[i] - x_old.v[i]);
			norm_x += x.v[i]*x.v[i];
		}
		norma = sqrt(error/norm_x);
		for(int i = 0; i < n ; i++){
			x_old.v[i] = x.v[i];
		}
	}
	free_Vector(&x_old);
	return x;
}

int solve_by(Matrix *A, Vector *b, Vector*out) {
    Vector x;
    
    x = solve_by_Cholensky(A, b);
    if (x.n > 0){
		copy_vector(&x,out);
		return 0;
	}

    x = solve_by_Crouth(A, b);
    if (x.n > 0){
		copy_vector(&x,out);
		return 0;
	}

    x = solve_by_Gauss_Seidel(A, b);
    if (x.n > 0){
		copy_vector(&x,out);
		return 0;
	}

    x = solve_by_Jacobi(A, b);
    if (x.n > 0){
		copy_vector(&x,out);
		return 0;
	}

    printf("Error: no se pudo resolver el sistema\n");
    return -1; 
}

#endif
