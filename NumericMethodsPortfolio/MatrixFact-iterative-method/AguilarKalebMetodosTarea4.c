
//AguilarKalebMetodosTarea4.c
#include <stdio.h>
#include <math.h>
#include "Matrix.h"

//FUNCIONES PROBLEMA 1
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
		printf("Factorización Imposible A[0][0] = 0\n");
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
			printf("Factorización Imposible L[%d][%d] = 0\n",i,i);
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
	Factorize_by_Crouth(A,b ,&L, &U);
	
	Vector y = solve_L(&L,b);
	Vector x = solve_U(&U,&y);
	//print_Vector(y);
	//print_Vector(x);
	free_Matrix(&L);
	free_Matrix(&U);
	free_Vector(&y);
	return x;
}

//FUNCIONES PROBLEMA 3
int Factorize_by_Cholensky(Matrix *A, Vector *b,Matrix *L, Matrix *LT){
	
	int n = A->Y;
	double sum;
	
	//PASO1
	if (A->m[0][0] < 0){
		printf("Factorización Imposible A[0][0] < 0\n");
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
			printf("Factorización Imposible A[%d][%d] < 0\n",i,i);
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
	Factorize_by_Cholensky(A,b ,&L, &LT);
	
	Vector y = solve_L(&L,b);
	Vector x = solve_U(&LT,&y);
	//print_Vector(y);
	//print_Vector(x);
	free_Matrix(&L);
	free_Matrix(&LT);
	free_Vector(&y);
	return x;
}

//FUNCIONES PROBLEMA 4
double BD_matrix(int i, int j , Vector *Band, int rows){
	if (i < 0 || i >= rows || j < 0 || j >= rows){
		printf("Indices de BD_matrix no válidos i = %d, j = %d",i,j);
		return -3000;
	}
	int diff = abs(i - j);
	if ( (diff < Band->n) && (i-j >= 0) ){
		return Band->v[diff];
	}
	return 0;
}
int Factorize_Band_Diagonal_by_Cholensky(Vector *Band,int rows, Vector *b,Matrix *L, Matrix *LT){
	int n = rows;
	double sum;
	
	//PASO1
	if (BD_matrix(0,0 ,Band,rows) < 0){
		printf("Factorización Imposible A[0][0] < 0\n");
		return -1;
	}
	
	L->m[0][0] = sqrt(BD_matrix(0,0 ,Band,rows));
	LT->m[0][0] = L->m[0][0];
	
	//PASO2
	for(int j = 1 ; j < n ; j++){
		L->m[j][0] = BD_matrix(j,0 ,Band,rows) / L->m[0][0];
		LT->m[0][j] = L->m[j][0];
	}
	
	//PASO3
	for(int i = 1 ; i < n - 1 ; i++){
		//PASO4
		sum = 0;
		for(int k = 0; k < i ; k++){
			sum += L->m[i][k]*L->m[i][k];
		}
		sum = BD_matrix(i,i,Band,rows) - sum;
		if(sum < 0){
			printf("Factorización Imposible A[%d][%d] < 0\n",i,i);
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
			L->m[j][i] = (BD_matrix(j,i ,Band,rows) - sum) / L->m[i][i];
			LT->m[i][j] = L->m[j][i];	
		}	
		
		//PASO6
		sum = 0;
		for(int k = 0 ; k < n - 1; k++){
			sum += L->m[n-1][k]*L->m[n-1][k] ; 
		}
		sum = BD_matrix(n-1,n-1 ,Band,rows) - sum;
		if (sum < 0){
			printf("Factorización Imposible A[%d][%d] < 0\n",n-1,n-1);
			return-1;
		}
		L->m[n-1][n-1] = sqrt(sum);
		LT->m[n-1][n-1] = L->m[n-1][n-1];
	}
	return 0;
}
Vector eq_calor(double Q, double K , double L, double Phi_0, double Phi_n){
	
	Vector Band = initialize_Vector(2);
	Band.v[0] = 2;
	Band.v[1] = -1; // Band = {2,-1}
	
	int n = 5;
	double E_dPhi = 1;
	double E_norm = fabs(Phi_0 - Phi_n);
	Vector Phi;
	while (E_dPhi > 0.001){
		double dx = L / (double)n;
		Matrix L = initialize_Matrix(n,n);
		Matrix LT = initialize_Matrix(n,n);
		Vector b = initialize_Vector(n);
		for(int i = 0 ; i < b.n ; i++){
			b.v[i] = -Q*dx*dx/K;
		}
		b.v[0] += Phi_0;
		b.v[n-1] += Phi_n; 
		
		Factorize_Band_Diagonal_by_Cholensky(&Band,n,&b ,&L, &LT);
		
		Vector y = solve_L(&L,&b);
		Phi = solve_U(&LT,&y);
		double count_dPhi = 0;
		for(int i = 0 ; i < Phi.n-1; i++){
			count_dPhi += fabs(Phi.v[i] - Phi.v[i+1]);
		}
		E_dPhi = (count_dPhi / (n-1)) / E_norm;
		
		free_Matrix(&L);
		free_Matrix(&LT);
		free_Vector(&y);
		printf("n: %d , EdPhi: %f\n",n,E_dPhi);
		n = 2*n +1;
	}
	return Phi;
}

//FUNCIONES PROBLEMA 5
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

//FUNCIONES PROBLEMA 6
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


//IMPRESIÓN DE RESULTADOS
void Problema1(){
	printf("\n\n\n----------Problema 1: Método de Crouth ------------\n");

	FILE *file_A = fopen("A.txt","r");
	FILE *file_b = fopen("b.txt","r");
	Matrix A = create_Matrix_from_FILE(file_A,' ');
	Vector b = create_Vector_from_FILE(file_b,' ');
	Vector x = solve_by_Crouth(&A,&b);
	printf("Dimensiones de A: %d x %d\n", A.Y, A.X);
	printf("Dimensiones de b: %d\n", b.n);
	printf("Dimensiones de x: %d\n", x.n);
	//print_Vector(x);
	// error = ||Ax - b||
    Vector Ax = Matrix_times_Vector(&A, &x);       
    Vector error = diff_Vector(&Ax, &b);           
    double diff = Vector_norm(&error);
    printf("||Ax - b|| = %lf\n", diff);
	//print_solution(A,b,x);

    free_Vector(&Ax);
    free_Vector(&error);
    free_Vector(&x);
    free_Vector(&b);
    free_Matrix(&A);
}
void Problema3(){
	printf("\n\n\n----------Problema 3: Método de Cholensky ------------\n");

	FILE *file_A = fopen("SPD.txt","r");
	FILE *file_b = fopen("b_spd.txt","r");
	Matrix A = create_Matrix_from_FILE(file_A,' ');
	Vector b = create_Vector_from_FILE(file_b,' ');
	Vector x = solve_by_Cholensky(&A,&b);
	// Imprimir dimensiones
	printf("Dimensiones de A: %d x %d\n", A.Y, A.X);
	printf("Dimensiones de b: %d\n", b.n);
	printf("Dimensiones de x: %d\n", x.n);
	//print_Vector(b);
	// error = ||Ax - b||
    Vector Ax = Matrix_times_Vector(&A, &x);       
    Vector error = diff_Vector(&Ax, &b);           
    double diff = Vector_norm(&error);
    printf("||Ax - b|| = %lf\n",diff);
	//print_solution(A,b,x);

    free_Vector(&Ax);
    free_Vector(&error);
    free_Vector(&x);
    free_Vector(&b);
    free_Matrix(&A);
}
void Problema4(){
	printf("\n\n\n----------Problema 4: Ecuación de Calor ------------\n");
	double Q = 5;
	double K = 1;
	double L = 1;
	double Phi_0 = 10;
	double Phi_n = 100;
	Vector Phi = eq_calor(Q,K,L,Phi_0,Phi_n);
	save_Vector_to_csv(&Phi, "Phi.csv");
	free_Vector(&Phi);
}
void Problema5(){
    printf("\n\n\n----------Problema 5: Método de Jacobi ------------\n");

    // A_small
    FILE *file_A_small = fopen("A_small.txt","r");
    FILE *file_b_small = fopen("b_small.txt","r");

    Matrix A_small = create_Matrix_from_FILE(file_A_small,' ');
    Vector b_small = create_Vector_from_FILE(file_b_small,' ');
    Vector x_small = solve_by_Jacobi(&A_small,&b_small);

    // DIMENSIONES
    printf("Dimensiones de A_small: %d x %d\n", A_small.Y, A_small.X);
    printf("Dimensiones de b_small: %d\n", b_small.n);
    printf("Dimensiones de x_small: %d\n", x_small.n);

    //REPORTE DE ERROR
    Vector Ax_small = Matrix_times_Vector(&A_small, &x_small);       
    Vector error_small = diff_Vector(&Ax_small, &b_small);           
    double diff_small = Vector_norm(&error_small);
    printf("||A_small * x_small - b_small|| = %lf\n", diff_small);
	print_solution(A_small,b_small,x_small);
    //FREE
    free_Vector(&Ax_small);
    free_Vector(&error_small);
    free_Vector(&x_small);
    free_Vector(&b_small);
    free_Matrix(&A_small);

    //A_big
    FILE *file_A_big = fopen("A_big.txt","r");
    FILE *file_b_big = fopen("b_big.txt","r");

    Matrix A_big = create_Matrix_from_FILE(file_A_big,' ');
    Vector b_big = create_Vector_from_FILE(file_b_big,' ');
    Vector x_big = solve_by_Cholensky(&A_big,&b_big);

    //DIMENSIONES
    printf("\n\nDimensiones de A_big: %d x %d\n", A_big.Y, A_big.X);
    printf("Dimensiones de b_big: %d\n", b_big.n);
    printf("Dimensiones de x_big: %d\n", x_big.n);

    //ERROR
    Vector Ax_big = Matrix_times_Vector(&A_big, &x_big);       
    Vector error_big = diff_Vector(&Ax_big, &b_big);           
    double diff_big = Vector_norm(&error_big);
    printf("||A_big * x_big - b_big|| = %lf\n", diff_big);

    //FREE
    free_Vector(&Ax_big);
    free_Vector(&error_big);
    free_Vector(&x_big);
    free_Vector(&b_big);
    free_Matrix(&A_big);

}
void Problema6(){
    printf("\n\n\n----------Problema 6: Método de Gauss Seidel ------------\n");

    // A_small
    FILE *file_A_small = fopen("A_small.txt","r");
    FILE *file_b_small = fopen("b_small.txt","r");

    Matrix A_small = create_Matrix_from_FILE(file_A_small,' ');
    Vector b_small = create_Vector_from_FILE(file_b_small,' ');
    Vector x_small = solve_by_Gauss_Seidel(&A_small,&b_small);

    // DIMENSIONES
    printf("Dimensiones de A_small: %d x %d\n", A_small.Y, A_small.X);
    printf("Dimensiones de b_small: %d\n", b_small.n);
    printf("Dimensiones de x_small: %d\n", x_small.n);

    //REPORTE DE ERROR
    Vector Ax_small = Matrix_times_Vector(&A_small, &x_small);       
    Vector error_small = diff_Vector(&Ax_small, &b_small);           
    double diff_small = Vector_norm(&error_small);
    printf("||A_small * x_small - b_small|| = %lf\n", diff_small);
	print_solution(A_small,b_small,x_small);
    //FREE
    free_Vector(&Ax_small);
    free_Vector(&error_small);
    free_Vector(&x_small);
    free_Vector(&b_small);
    free_Matrix(&A_small);

    //A_big
    FILE *file_A_big = fopen("A_big.txt","r");
    FILE *file_b_big = fopen("b_big.txt","r");

    Matrix A_big = create_Matrix_from_FILE(file_A_big,' ');
    Vector b_big = create_Vector_from_FILE(file_b_big,' ');
    Vector x_big = solve_by_Gauss_Seidel(&A_big,&b_big);

    //DIMENSIONES
    printf("\n\nDimensiones de A_big: %d x %d\n", A_big.Y, A_big.X);
    printf("Dimensiones de b_big: %d\n", b_big.n);
    printf("Dimensiones de x_big: %d\n", x_big.n);

    //ERROR
    Vector Ax_big = Matrix_times_Vector(&A_big, &x_big);       
    Vector error_big = diff_Vector(&Ax_big, &b_big);           
    double diff_big = Vector_norm(&error_big);
    printf("||A_big * x_big - b_big|| = %lf\n", diff_big);

    //FREE
    free_Vector(&Ax_big);
    free_Vector(&error_big);
    free_Vector(&x_big);
    free_Vector(&b_big);
    free_Matrix(&A_big);
}


int main(int argc, char **argv)
{
	Problema1();
	Problema3();
	Problema4();
	Problema5();
	Problema6();
	return 0;
}

