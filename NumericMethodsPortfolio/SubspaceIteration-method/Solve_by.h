//  Created by Kaleb Aguilar on 13/08/25.
//
#ifndef SOLVE_BY_H
#define SOLVE_BY_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "Matrix.h"

/// SOLUCIÓN DE SISTEMAS 
Vector solve_D(Matrix *m, Vector *b) {
    Vector x = initialize_Vector(b->n);
    for (int i = 0; i < b->n; i++) {
        x.v[i] =  b->v[i] / m->m[i][i];
    }
    return x;
}

Vector solve_U(Matrix *m, Vector *b) {
    Vector x = initialize_Vector(b->n);
    for (int i = b->n - 1; i >= 0; i--) {
        x.v[i] = b->v[i];
        for (int j = b->n - 1; j > i; j--) {
            x.v[i] -= m->m[i][j] * x.v[j];
        }
        x.v[i] /= m->m[i][i];
    }
    return x;
}

Vector solve_L(Matrix *m, Vector *b) {
    Vector x = initialize_Vector(b->n);
    for (int i = 0; i < b->n; i++) {
        x.v[i] = b->v[i];
        for (int j = 0; j < i; j++) {
            x.v[i] -= m->m[i][j] * x.v[j];
        }
        x.v[i] /= m->m[i][i];
    }
    return x;
}

int reduce_col_from_row(Matrix *m, Vector *b, int row_1, int row_2) {
    // ELIMINATE COL FROM ROW_2 USING ROW_1.
    int col = row_2;
    if (row_1 >= m->Y || row_2 >= m->Y) return -1;
    if (m->m[row_1][col] == 0) return -1;

    double factor = m->m[row_2][col] / m->m[row_1][col];
    for (int i = 0; i < m->X; i++) {
        m->m[row_2][i] -= factor * m->m[row_1][i];
    }
    b->v[row_2] -= factor * b->v[row_1];
    return 0;
}

void swap_rows(Matrix *m, Vector *b, int row_1, int row_2) {
    double aux;
    for (int i = 0; i < m->X; i++) {
        aux = m->m[row_1][i];
        m->m[row_1][i] = m->m[row_2][i];
        m->m[row_2][i] = aux;
    }
    aux = b->v[row_1];
    b->v[row_1] = b->v[row_2];
    b->v[row_2] = aux;
}

void pivote(Matrix *A, Vector *b, int col){
    if(A->m[col][col] != 0) return;

    for(int i = col + 1; i < A->Y; i++){
        if(A->m[i][col] != 0){
            // SWAP DE FILAS i <--> col
            for(int j = 0; j < A->X; j++){
                double tmp = A->m[col][j];
                A->m[col][j] = A->m[i][j];
                A->m[i][j] = tmp;
            }
            // SWAP EN VECTOR B
            if(b){
                double tmpb = b->v[col];
                b->v[col] = b->v[i];
                b->v[i] = tmpb;
            }
            return;
        }
    }
}

Vector solve_by_Gauss(Matrix *m, Vector *b) {
    for (int i = 0; i < m->Y; i++) {
        for (int j = i + 1; j < m->Y; j++) {
            if (i == j) continue;
            reduce_col_from_row(m, b, i, j);
        }
    }
    Vector x = solve_U(m, b);
    return x;
}

void print_solution(Matrix m, Vector b, Vector x) {
    for (int i = 0; i < m.Y; i++) {
        for (int j = 0; j < m.X; j++) {
            printf("%10.Lf ", m.m[i][j]);
        }
        printf("   |x%2d: %10.10Lf |b%2d: %10.10Lf", i, x.v[i], i, b.v[i]);
        printf("\n");
    }
    printf("\n");
}



///SOLVER SIMPLE
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
			printf("Factorización Imposible A[%d][%d] < 0\n",n-1,n-1);
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
    
    

    x = solve_by_Crouth(A, b);
    if (x.n > 0){
		printf("Solved by Crouth\n");
		copy_vector(&x,out);
		return 0;
	}

    x = solve_by_Gauss_Seidel(A, b);
    if (x.n > 0){
		printf("Solved by Gauss-Seidel\n");
		copy_vector(&x,out);
		return 0;
	}

    x = solve_by_Jacobi(A, b);
    if (x.n > 0){
		printf("Solved by Jacobi\n");
		copy_vector(&x,out);
		return 0;
	}
	x = solve_by_Cholensky(A, b);
    if (x.n > 0){
		printf("Solved by Cholensky\n");
		copy_vector(&x,out);
		return 0;
	}

    printf("Error: no se pudo resolver el sistema\n");
    return -1; 
}


int hide_solve_by(Matrix *A, Vector *b, Vector*out) {
    Vector x;
    
    

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
	x = solve_by_Cholensky(A, b);
    if (x.n > 0){
		copy_vector(&x,out);
		return 0;
	}

    printf("Error: no se pudo resolver el sistema\n");
    return -1; 
}

void Normalize(Vector *v){
	double normal = 0;
	for(int i = 0; i < v->n ; i++){
		normal += v->v[i]*v->v[i];
	}
	normal = sqrt(normal);
	for(int i = 0; i < v->n ; i++){
		v->v[i]/= normal;
	}
} 
void Metodo_de_Potencia_seed(Matrix *A, int m,Matrix *v_seed ,Vector *eigvals, Vector **eigvecs){ 
	int n = A->Y; 
	Vector v1 = initialize_Vector(n); 
	for (int j = 0; j < m; j++) { 
		Vector v0 = initialize_Vector(n); 
		for (int i = 0; i < n; i++) v0.v[i] = v_seed->m[i][j]; 
		double lambda = 0.0, lambda_old = 0.0; int iter = 0, max_iter = 1000000; 
		do { 
			lambda_old = lambda; 
			Matrix_times_Vector(A, &v0, &v1); 
					// V1 -= <V1|V_prev> 
			for (int prev = 0; prev < j; prev++) { 
				double prod_inter = 0.0; 
				//<V1|V_prev> 
				for (int i = 0; i < n; i++) prod_inter += v1.v[i] * eigvecs[prev]->v[i]; 
				
				// V1 -= <V1|V_prev> 
				for (int i = 0; i < n; i++) v1.v[i] -= prod_inter * eigvecs[prev]->v[i]; 
				
			} 
			
			double num = 0.0, den = 0.0; 
				
			for (int i = 0; i < n; i++) { 
				num += v1.v[i] * v1.v[i]; 
				den += v1.v[i] * v0.v[i]; 
			} 
			
			lambda = num / den; Normalize(&v1); 
			iter++; 
			copy_vector(&v1, &v0); 
		} while (fabs(lambda - lambda_old) > 1/INT_MAX && iter < max_iter); 
		
		if (iter >= max_iter) printf("Límite de iteraciones excedido para eigenvector %d\n", j+1); 
		eigvals->v[j] = lambda; 
		copy_vector(&v0, eigvecs[j]); 
		free_Vector(&v0); 
		
	} 
	free_Vector(&v1); 
}



#endif
