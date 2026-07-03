//  Created by Kaleb Aguilar on 13/08/25.
//
#ifndef SOLVE_BY_H
#define SOLVE_BY_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "Matrix.h"

typedef struct {
	int i;
	int j;
	double value;
}max_value;

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
double error_norm(Matrix *A, Vector *x, Vector *b) {
    int n = b->n;
    double norm = 0.0;

    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += A->m[i][j] * x->v[j];
        }
        double r_i = sum - b->v[i];  // residuo en componente i
        norm += r_i * r_i;
    }

    return sqrt(norm);
}



/// SOLVERS DE FACTORIZACIÓN 
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
int Factorize_by_Cholensky(Matrix *A, Matrix *L, Matrix *LT) {
    int n = A->Y;
    double sum;

    // inicializar a cero
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            L->m[i][j] = 0.0;
            LT->m[i][j] = 0.0;
        }
    }

    // Paso 1
    if (A->m[0][0] <= 0) {
        return -1; // no es SPD
    }

    L->m[0][0] = sqrt(A->m[0][0]);
    LT->m[0][0] = L->m[0][0];

    // Paso 2: 
    for (int j = 1; j < n; j++) {
        L->m[j][0] = A->m[j][0] / L->m[0][0];
        LT->m[0][j] = L->m[j][0];
    }

    // Paso 3: 
    for (int i = 1; i < n; i++) {
        // diagonal
        sum = 0.0;
        for (int k = 0; k < i; k++) {
            sum += L->m[i][k] * L->m[i][k];
        }
        double diag = A->m[i][i] - sum;
        if (diag <= 0) {
            return -1;
        }
        L->m[i][i] = sqrt(diag);
        LT->m[i][i] = L->m[i][i];

        for (int j = i + 1; j < n; j++) {
            sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += L->m[j][k] * L->m[i][k];
            }
            L->m[j][i] = (A->m[j][i] - sum) / L->m[i][i];
            LT->m[i][j] = L->m[j][i];
        }
    }

    return 0; // ok
}

Vector solve_by_Cholensky(Matrix *A, Vector *b){
	Matrix L = initialize_Matrix(A->Y,A->X);
	Matrix LT = initialize_Matrix(A->Y,A->X);
	Factorize_by_Cholensky(A ,&L, &LT);
	if(Factorize_by_Cholensky(A, &L, &LT) != 0){
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
void solve_by_QR_inplace(Matrix *A, Vector *b, Vector *out) {
    int n = A->Y;
    // Q: n x m, R: m x m
    Matrix Q = initialize_Matrix(n, n);
    Matrix R = initialize_Matrix(n, n);
    Vector a_j = initialize_Vector(n);
    Vector q_i = initialize_Vector(n);
    Vector r_ij_q_i = initialize_Vector(n);
    Vector aux_row = initialize_Vector(n);

    for (int j = 0; j < n; j++) {
        Vector_col(A, j, &a_j);
        // aux_row = a_j
        for (int k = 0; k < n; k++) aux_row.v[k] = a_j.v[k];

        for (int i = 0; i < j; i++) {
            Vector_col(&Q, i, &q_i);                    // q_i (long n)
            long double R_ij = (long double) Vector_times_Vector(&q_i, &aux_row);
            R.m[i][j] = R_ij;
            scalar_times_Vector((double)R_ij, &q_i, &r_ij_q_i);
            for (int k = 0; k < n; k++) aux_row.v[k] -= r_ij_q_i.v[k];
        }

        long double R_jj = (long double) Vector_norm(&aux_row);
        R.m[j][j] = R_jj;
        if (R_jj < 1e-18L) {
            for (int k = 0; k < n; k++) Q.m[k][j] = 0.0;
        } else {
            for (int k = 0; k < n; k++) Q.m[k][j] = aux_row.v[k] / R_jj;
        }
    }

    Vector y = initialize_Vector(n);
    for (int j = 0; j < n ; j++) {
        Vector_col(&Q, j, &q_i);
        y.v[j] = (long double) Vector_times_Vector(&q_i, b);
    }

    Vector x = solve_U(&R, &y);

    for (int i = 0; i < n; i++) out->v[i] = x.v[i];

    free_Matrix(&Q);
    free_Matrix(&R);
    free_Vector(&a_j);
    free_Vector(&q_i);
    free_Vector(&r_ij_q_i);
    free_Vector(&aux_row);
    free_Vector(&y);
    free_Vector(&x);

    return;
}

//SOLVERS ITERATIVOS
Vector solve_by_Jacobi(Matrix *A ,Vector *b){
	int n = b->n;
	double sum, error,norm_x;
	double norma = 1;
	double toler = 0.000001;
	
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
	double toler = 0.000001;
	
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
void solve_by_Gradiente_Conjugado(Matrix *A, Vector *b, Vector *x) {
	
	int k = 0;
	long double norma_r = 10000;
	long double alpha,beta;
	
	Vector aux = initialize_Vector(b->n);
	Vector x_k = initialize_Vector(b->n);
	Vector r_k = initialize_Vector(b->n);  
	Vector P_k = initialize_Vector(b->n);
	Vector w_k = initialize_Vector(b->n);
	
	copy_Vector(b,&r_k);
	copy_Vector(b,&P_k);
	
	Vector x_next = initialize_Vector(b->n);
	Vector r_next = initialize_Vector(b->n); 
	Vector P_next = initialize_Vector(b->n); 
	
	while (norma_r > pow(10,-12) && k < 100000){
		// w_k = A P_k
		Matrix_times_Vector(A,&P_k,&w_k);
		
		//alpha = P_k*r_k / P_k*w_k
		long double num = Vector_times_Vector(&P_k,&r_k);
		long double den = Vector_times_Vector(&P_k,&w_k);
		alpha = num / den;
		
		// X_k+1 = X_k + alpha*P_k
		scalar_times_Vector(alpha,&P_k,&aux);
		sum_Vectors(&x_k,&aux,&x_next);
		
		//X_k+1 = X_k - alpha*w_k
		alpha *= -1;
		scalar_times_Vector(alpha,&w_k,&aux);
		sum_Vectors(&r_k,&aux,&r_next);
		
		//beta = P_k*r_next / P_k*P_k
		copy_Vector(&P_k,&aux);
		num = Vector_times_Vector(&P_k,&r_next);	
		den = Vector_times_Vector(&P_k,&aux);
		beta = num / den;
		
		// P_next = r_next + beta*P_k;
		scalar_times_Vector(beta,&P_k,&aux);
		sum_Vectors(&r_next,&aux,&P_next);
		
		norma_r = Vector_norm(&r_next);
		
		copy_Vector(&x_next,&x_k);
		copy_Vector(&r_next,&r_k);
		copy_Vector(&P_next,&P_k);
		
		k++;
	}
	copy_Vector(&x_next,x);
	
}


//SOLVER_SIMPLE
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
    x = solve_by_Cholensky(A, b);
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
    
    

    x = solve_by_Crouth(A, b);
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




// VECTORES Y VALORES PROPIOS
void verify_eigen(Matrix *A, Vector *v, double lambda) {
	double tol = 0.01;
    int n = v->n;
    Vector Av = initialize_Vector(n);
    Matrix_times_Vector(A, v, &Av);
    for (int i = 0; i < n; i++) {
        Av.v[i] -= lambda * v->v[i];
    }
    double res_norm = Vector_norm(&Av);

    if (res_norm < tol) 
        printf("Eigenvector verificado correctamente con lambda = %lf y |Av - λv| = %e \n", lambda,res_norm);
    else 
        printf("ERROR: No es un eigenvector on lambda = %lf y |Av - λv| = %e\n", lambda, res_norm);

    free_Vector(&Av);
}

// POTENCIA
void Potencia_K(Matrix *A, int k, Vector *eigvals, Vector **eigvecs){ 
	int n = A->Y; 
	Vector v1 = initialize_Vector(n); 
	for (int j = 0; j < k; j++) { 
		Vector v0 = initialize_Vector(n); 
		for (int i = 0; i < n; i++) v0.v[i] = 1.0/sqrt(n); 
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

//POTENCIA INVERSA
void Potencia_Inversa_K(Matrix *A, int k, Vector *eigvals, Vector **eigvecs) {
    int n = A->Y;
	Vector v1 = initialize_Vector(n);
	Vector Av = initialize_Vector(n);
    for (int j = 0; j < k; j++) {
        Vector v0 = initialize_Vector(n);
		for (int i = 0; i < n; i++) v0.v[i] =1.0/sqrt(n);
		
        double lambda = 0.0, lambda_old = 0.0;
        int iter = 0, max_iter = 1000000;
		
        do {
            lambda_old = lambda;
			hide_solve_by(A, &v0, &v1);
			
			
			//print_Vector(v0);
            /*if(hide_solve_by(A, &v0, &v1) != 0){
				printf("No se pudo resolver el sistema en la iteración de potencia inversa.\n");
				break;
			}*/ 
			
			// V1 -= <V1|V_prev> 
            for (int prev = 0; prev < j; prev++) {
                double prod_inter = 0.0;
                //<V1|V_prev> 
                for (int i = 0; i < n; i++)
                    prod_inter += v1.v[i] * eigvecs[prev]->v[i];
                // V1 -= <V1|V_prev> 
                for (int i = 0; i < n; i++)
                    v1.v[i] -= prod_inter * eigvecs[prev]->v[i];
            }

			//Matrix_times_Vector(A,&v1,&Av);
            double num = 0.0, den = 0.0;
            for (int i = 0; i < n; i++) {
				
                num += v1.v[i] * v1.v[i]; 
				den += v1.v[i] * v0.v[i]; 
            }
            lambda = num / den;

            Normalize(&v1);
            iter++;
			printf("LAMBDA:%lf\n",lambda);
            copy_vector(&v1, &v0);

        } while (fabs(lambda - lambda_old) > 1e-12 && iter < max_iter);
	
        if (iter >= max_iter)
            printf("Límite de iteraciones excedido para eigenvector %d\n", j+1);
		
        eigvals->v[j] = 1.0 / lambda; 
        copy_vector(&v0, eigvecs[j]);
        free_Vector(&v0);
       
    }
    free_Vector(&Av);
    free_Vector(&v1);
}


//ITERACION DE SUBESPACIO
max_value max_without_diagonal(Matrix *m){
	max_value out;
	out.i = -1;
	out.j = -1;
	out.value = 0;
	for(int i = 0; i < m->Y ; i++){
		for(int j = 0; j < m->X; j++){
			if (i == j) continue;
			if ( fabsl(out.value) < fabsl(m->m[i][j])){ 
				out.value = m->m[i][j];
				out.i = i;
				out.j = j;
			}
		}
	}
	return out;
}
void Metodo_Jacobi(Matrix *A, Matrix *I){
    int n = A->Y;
    double tol = 1e-6; // tolerancia estricta
    max_value max = max_without_diagonal(A);
    int iter = 0;

    while(fabsl(max.value) > tol && iter < 1e6){
        int p = max.i;
        int q = max.j;

        if(p == -1 || q == -1) break;

        double app = A->m[p][p];
        double aqq = A->m[q][q];
        double apq = A->m[p][q];

        // Ángulo de rotación
        double theta = 0.5 * atan2(2*apq, aqq - app);
        double c = cos(theta);
        double s = sin(theta);
		double t = s/c;
        // ROTACVION DE A
        A->m[p][p] = app - t*apq; 
        A->m[q][q] = aqq + t*apq; 
        A->m[p][q] = 0.0;
        A->m[q][p] = 0.0;

        for(int k=0;k<n;k++){
            if(k != p && k != q){
                double aik = A->m[p][k];
                double aqi = A->m[q][k];
                A->m[p][k] = c*aik - s*aqi;
                A->m[k][p] = A->m[p][k]; 
                A->m[q][k] = s*aik + c*aqi;
                A->m[k][q] = A->m[q][k]; 
            }
        }

        //ROTACION DE I 
        for(int k=0;k<n;k++){
            double vip = I->m[k][p];
            double viq = I->m[k][q];
            I->m[k][p] = c*vip - s*viq;
            I->m[k][q] = s*vip + c*viq;
        }

        iter++;
        max = max_without_diagonal(A);
    }
}
void normalize_column(Matrix *M, int col) {
    double norm = 0.0;
    for (int i = 0; i < M->Y; i++) {
        norm += M->m[i][col] * M->m[i][col];
    }
    norm = sqrt(norm);

    if (norm < 1e-14) return; 
    for (int i = 0; i < M->Y; i++) {
        M->m[i][col] /= norm;
    }
}
void gram_schmidt_modified(Matrix *M) {
    int n = M->Y;
    int m = M->X;

    for (int k = 0; k < m; k++) {
        // Normalizar la columna k antes de proyectar sobre siguientes
        normalize_column(M, k);

        for (int j = k+1; j < m; j++) {
            double dot = 0.0;
            for (int i = 0; i < n; i++)
                dot += M->m[i][j] * M->m[i][k];

            for (int i = 0; i < n; i++)
                M->m[i][j] -= dot * M->m[i][k];
        }
    }
}
void Iteracion_de_subespacio(Matrix *A, int m, Vector *eigval, Vector **eigvecs) {
    int n = A->Y; 
    double tol = 1e-4;
	
	// Inicializamos matrices auxiliares
	Matrix I_0 = initialize_Matrix(n, m);
	Matrix I_1 = initialize_Matrix(n, m);
	Matrix I_1_T = initialize_Matrix(m, n);
	Matrix B = initialize_Matrix(m, n);
	Matrix A_I_0 = initialize_Matrix(n,m);
	Matrix C = initialize_Matrix(m, m);
	Matrix V = initialize_Matrix(n,m);
	Vector v_i = initialize_Vector(n);
	
	Vector tmp = initialize_Vector(n);
	
	
	// Inicializamos I0 como n x m con valores aleatorios
	
	for (int i = 0; i < n; i++)
	for (int j = 0; j < m; j++)
		I_0.m[i][j] = 2.0*((double)rand() / RAND_MAX) - 1.0; // valores en [-1,1]

	gram_schmidt_modified(&I_0); // normaliza y ortogonaliza
			
	double max_err;
	int iter = 0;
	do {
		
		//SHIFT I_0 -> I_1 
		Matrix_times_Matrix(A, &I_0, &I_1);
		
		for (int j=0; j<m; j++) normalize_column(&I_1, j);
		gram_schmidt_modified(&I_1);
		// C = I_1_T * A * I1
		transpose(&I_1, &I_1_T);
		Matrix_times_Matrix(&I_1_T, A, &B); // B = I_1^T * A   (m x n)
		Matrix_times_Matrix(&B, &I_1, &C);  // C = I_1^T * A * I_1  (m x m)

		Matrix S = initialize_I(m, m); // EIGENVALORES DE JACOBI
		
		// Jacobi :=  C (I_1^T * A * I_1) ,S (identidad) -> C (Eigenvalores)   S (eigenvectores),
		Metodo_Jacobi(&C, &S);

		//ORDENAR S Y C
		for (int i = 0; i < m-1; i++) {
			int idx_max = i;
			for (int j = i+1; j < m; j++) {
				if (C.m[j][j] > C.m[idx_max][idx_max])
					idx_max = j;
			}
			if (idx_max != i) {
				//SWAP EN C
				double tmp_val = C.m[i][i];
				C.m[i][i] = C.m[idx_max][idx_max];
				C.m[idx_max][idx_max] = tmp_val;

				//SWAP EN S
				for (int k = 0; k < m; k++) {
					double tmp_s = S.m[k][i];
					S.m[k][i] = S.m[k][idx_max];
					S.m[k][idx_max] = tmp_s;
				}
			}
		}

		Matrix_times_Matrix(&I_1,&S,&V);
		
		max_err = 0.0;
		
		for (int i=0; i<m; i++){
			double lambda_i = C.m[i][i]; 
			// v = V[:,i]
			Vector_col(&V, i, &v_i);
			// tmp = A * v_i
			Matrix_times_Vector(A, &v_i, &tmp);  
			for (int k=0;k<n;k++) tmp.v[k] -= lambda_i * v_i.v[k];
			double err_i = Vector_norm(&tmp);
			if (err_i > max_err) max_err = err_i;
		}
		
		if (iter % 50 == 0){
			printf("iter: %d error: %f \n",iter,max_err);
		}
		
		// Preparar siguiente iteración
		copy_Matrix(&V,&I_0);
		for(int j=0;j<m;j++) normalize_column(&I_0, j);
		gram_schmidt_modified(&I_0);

		iter++;
		
		free_Matrix(&S);
	} while (max_err > tol && iter < 1000);
	
	
	//OUTPUT
	for (int i = 0; i < m; i++) { 
		eigval->v[i] = C.m[i][i]; // eigenvalue 
		for (int j = 0; j < n; j++) { 
			eigvecs[i]->v[j] = V.m[j][i]; // eigenvector i, componente j 
		} 
	}
	
	//FREE
	free_Vector(&tmp);	
	free_Vector(&v_i);
	free_Matrix(&I_0); 
	free_Matrix(&A_I_0);
	free_Matrix(&I_1); 
	free_Matrix(&I_1_T);
	free_Matrix(&B);
	free_Matrix(&C);
}


//ITERACION DE SUBESPACIO INVERSO
void solve_Matrix_system(Matrix *A, Matrix *B, Matrix *X) {
    int n = A->Y;
    int m = B->X;

    for (int j = 0; j < m; j++) {
        Vector b = initialize_Vector(n);
        Vector_col(B, j,&b);
        Vector x = initialize_Vector(n);

		solve_by_QR_inplace(A,&b ,&x);
		
        // Guardar en X
        for (int i = 0; i < n; i++) {
            X->m[i][j] = x.v[i];
        }

        free_Vector(&b);
        free_Vector(&x);
    }

}
void Iteracion_de_subespacio_inverso(Matrix *A, int m, Vector *eigval, Vector **eigvecs) {
	Potencia_Inversa_K(A,1, eigval, eigvecs);
	long double sigma = 16;
	
    int n = A->Y; 
    int iter = 0;
    
    // Auxiliares
    Matrix I_0 = initialize_Matrix(n, m);
    Matrix I_1 = initialize_Matrix(n, m);
    Matrix I_1_T = initialize_Matrix(m, n);
    Matrix B = initialize_Matrix(m, n);
    Matrix S = initialize_I(m, m);
    Matrix C = initialize_Matrix(m, m);
    Matrix V = initialize_Matrix(n,m);
    Matrix AV = initialize_Matrix(n,m);
    Matrix VL = initialize_Matrix(n,m);
    Matrix R  = initialize_Matrix(n,m);

    Matrix A_modified = initialize_Matrix(n,n);
    copy_Matrix(A,&A_modified);


    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            I_0.m[i][j] = 2.0*((double)rand() / RAND_MAX) - 1.0;

    gram_schmidt_modified(&I_0);

    for (int i = 0; i < n; i++) {
        A_modified.m[i][i] = A->m[i][i] - sigma;
    }

    do {
        solve_Matrix_system(&A_modified, &I_0, &I_1);
        gram_schmidt_modified(&I_1);

        copy_Matrix(&I_1,&I_0);
        iter++;
    } while (iter < 100);

    // REDUCCIÓN DE ESPACIO
    transpose(&I_1, &I_1_T);
    Matrix_times_Matrix(&I_1_T, A, &B); 
    Matrix_times_Matrix(&B, &I_1, &C);  

    Metodo_Jacobi(&C, &S);
    Matrix_times_Matrix(&I_1,&S,&V);

    // Eigenvalores y eigenvectores
    for (int i = 0; i < m; i++) { 
        eigval->v[i] = C.m[i][i]; 
        for (int j = 0; j < n; j++) { 
            eigvecs[i]->v[j] = V.m[j][i]; 
        } 
    }
    // FREE
    free_Matrix(&I_0); free_Matrix(&I_1); free_Matrix(&I_1_T);
    free_Matrix(&B); free_Matrix(&C); free_Matrix(&S);
    free_Matrix(&V); free_Matrix(&AV); free_Matrix(&VL); free_Matrix(&R);
    free_Matrix(&A_modified);
}



#endif
