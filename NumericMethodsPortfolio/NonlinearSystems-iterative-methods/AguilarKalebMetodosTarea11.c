#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"
#include "Solve_by.h"

#define _USE_MATH_DEFINES

//PROBLEMA 1
int NoLinealSolvebyPuntoFijo(Vector *X0, long double (**G)(Vector*),Vector *Xout,long double tol){
	Vector error = initialize_Vector(X0->n);
	int iter = 0;
	do{
		for(int i = 0; i < X0->n;i++){
			Xout->v[i] = G[i](X0);
		}
		diff_Vector_inplace(X0,Xout,&error);
		for(int i = 0; i < X0->n;i++){
			X0->v[i] = Xout->v[i];
		}
		iter++;
	}while(Vector_norm(&error) > tol && iter < 1e7);
	
	free_Vector(&error);
	return iter;
}

//PROBLEMA 2
void Jacobian(
	Vector *X,
	long double (**F) (Vector *),
	long double dx,
	Matrix *out
	){
	for(int i = 0 ; i < X->n ; i++){
		
		long double f_x0 = F[i](X);
		 
		for(int j = 0 ; j < X->n ; j++){
			long double aux = X->v[j] ;
			X->v[j] += dx;
			long double f_x1 = F[i](X); 
			X->v[j]  =  aux;

			out->m[i][j] = (f_x1 - f_x0)/dx;
			
		}
	}
}

int NoLinealSolvebyBroyden(Vector *X0, long double (**F)(Vector*),Vector *Xout,long double tol, long double dx){
	
	Matrix A0 = initialize_Matrix(X0->n,X0->n);
	Matrix A1 = initialize_Matrix(X0->n,X0->n);
	Vector A_S = initialize_Vector(X0->n);
	Vector Y = initialize_Vector(X0->n); // Y = F(X^1) - F(X^0) 
	Vector F_X0 = initialize_Vector(X0->n);
	Vector F_X1 = initialize_Vector(X0->n);
	Vector S = initialize_Vector(X0->n); // S = X^1 - X^0 
	Vector X1 = initialize_Vector(X0->n);
	
	long double norm_S = 1e8;
	
	//==== A0,X0,X1====
	//-F(X0)
	for(int i = 0 ; i < X0->n ; i++){
		F_X0.v[i] = -F[i](X0);
	}
	//A0 = J(X^{k-1})
	Jacobian(X0,F,dx,&A0);
		
	int iter = 0;
	do{
		// NUEVO S
		solve_by_QR_inplace(&A0,&F_X0,&S);
		
		//X1= X0 + S
		sum_Vectors(X0,&S,&X1);
		//F1
		for(int i = 0 ; i < X1.n ; i++){
			F_X1.v[i] = F[i](&X1);
		}
		// Y = F1 + F0
		sum_Vectors(&F_X1,&F_X0,&Y);
		
		//A0*S
		Matrix_times_Vector(&A0,&S,&A_S);
		
		// norm = ||S||^2
		long double norm = Vector_norm(&S);
		norm_S = norm*norm; // |S|^2
		
		if (sqrt(norm_S) <  tol) break;
		
		// A1
		for(int i = 0; i < X0->n ; i++){
			for(int j = 0; j < X0->n ; j++){
				A1.m[i][j] = A0.m[i][j] + (Y.v[i] - A_S.v[i] ) *S.v[j]/norm;
			}
		}
		
		//X1 -> X0
		//F_X0 -> F_X1
		for(int i = 0 ; i < X1.n ; i++){
			X0->v[i] = X1.v[i];
			F_X0.v[i] = -F_X1.v[i];
			for (int j = 0 ; j < X1.n ; j++){
				A0.m[i][j] = A1.m[i][j];
			}
		}
		
		iter++;
	}while(iter < 1e7);
	
	for(int i = 0 ; i < X1.n ; i++){
		Xout->v[i] = X1.v[i];
	}
	
	free_Matrix(&A0);
	free_Matrix(&A1);
	free_Vector(&A_S);
	free_Vector(&X1);
	free_Vector(&F_X0);
	free_Vector(&F_X1);
	free_Vector(&S);
	free_Vector(&Y);
	return iter;
}

int NoLinealSolvebyNewton(Vector *X0, long double (**F)(Vector*),Vector *Xout,long double tol, long double dx){
	
	Vector F_X = initialize_Vector(X0->n);
	Matrix J = initialize_Matrix(X0->n,X0->n);
	Vector Y = initialize_Vector(X0->n);
	
	int iter = 0;
	do{
		//-F(X^{k-1})
		for(int i = 0 ; i < X0->n ; i++){
			F_X.v[i] = -1*F[i](X0);
		}

		//J(X^{k-1})
		Jacobian(X0,F,dx,&J);
		
		// J(X) * Y = -F(X)
		solve_by_QR_inplace(&J,&F_X,&Y);
		
		//X_k = X^{k-1} + Y
		sum_Vectors(X0,&Y,Xout);
		for(int i = 0; i < X0->n;i++){
			X0->v[i] = Xout->v[i];
		}
		iter++;
	}while(Vector_norm(&Y) > tol && iter < 1e7);
	free_Matrix(&J);
	free_Vector(&F_X);
	free_Vector(&Y);
	return iter;
}

//PROBLEMA 3
long double compute_error(Vector *X, long double (**F)(Vector*), int n){
    long double sum = 0;
    for(int i = 0; i < n; i++){
        long double val = F[i](X);
        sum += val*val;
    }
    return sqrtl(sum);
}

long double fa1(Vector *X){
	return X->v[0] + X->v[1] - 3;
}
long double fa2(Vector *X){
	return X->v[0]*X->v[0] + X->v[1]*X->v[1] - 9;
}
long double ga1(Vector *X){
	return 3 - X->v[1];;
}
long double ga2(Vector *X){
	long double val = 9 - X->v[0]*X->v[0];
    if(val < 0) val = 0; // evitar sqrt de negativo
    return sqrt(val);
}
void Problema3a(){
	int n = 2;
	Vector X0 = initialize_Vector(n);
	Vector XoutPuntoFijo = initialize_Vector(n);
	Vector XoutBroyden = initialize_Vector(n);
	Vector XoutNewton = initialize_Vector(n);
	Vector lambda = initialize_Vector(n);
	
	long double (*F[n])(Vector*);
	F[0] = fa1,F[1] = fa2;
	
	long double (*G[2])(Vector*);
	G[0] = ga1,G[1] = ga2;
	
	X0.v[0] = 1,X0.v[1] = 5;
	int iter_PuntoFijo = NoLinealSolvebyPuntoFijo(&X0,G,&XoutPuntoFijo,1e-6);
	long double error_PF = compute_error(&XoutPuntoFijo, F, n);
	
	X0.v[0] = 1,X0.v[1] = 5;
	int iter_Broyden = NoLinealSolvebyBroyden(&X0,F,&XoutBroyden,1e-6,1e-8);
	long double error_B = compute_error(&XoutBroyden, F, n);
	
	X0.v[0] = 1,X0.v[1] = 5;
	int iter_Newton = NoLinealSolvebyNewton(&X0,F,&XoutNewton,1e-6,1e-8);
	long double error_N = compute_error(&XoutNewton, F, n);
	printf("\n---- PROBLEMA 3a ----\n");
	printf("\n%-18s %-32s %-12s %-12s\n","Metodo", "Resultado (x1,x2,x3)", "Iteraciones", "   Error");
    printf("--------------------------------------------------------------------------------\n");
    printf("%-18s (%.8Lf, %.8Lf) 		 %-12d %.3Le\n","Punto Fijo",XoutPuntoFijo.v[0], XoutPuntoFijo.v[1],iter_PuntoFijo, error_PF);
    printf("%-18s (%.8Lf, %.8Lf) 		 %-12d %.3Le\n","Broyden",XoutBroyden.v[0], XoutBroyden.v[1],iter_Broyden, error_B);
    printf("%-18s (%.8Lf, %.8Lf) 		 %-12d %.3Le\n", "Newton",XoutNewton.v[0], XoutNewton.v[1],iter_Newton, error_N);

	
	free_Vector(&X0);
	free_Vector(&lambda);
	free_Vector(&XoutBroyden);
	free_Vector(&XoutNewton );
	
}



long double fb1(Vector *X){
	return 3*X->v[0]-cos(X->v[0]*X->v[2]) - 0.5;
}
long double fb2(Vector *X){
	return X->v[0]*X->v[0] - 81*(X->v[1] + 0.1)*(X->v[1] + 0.1) + sin(X->v[2]) + 1.06;
}
long double fb3(Vector *X){
	return exp(-1*X->v[0]*X->v[1]) + 20*X->v[2] + (10*M_PI-3)/3;
}
long double gb1(Vector *X){
	long double out = (cos(X->v[0]*X->v[2]) + 0.5)/ 3;
	//printf("GB1: %Lf\n",out);
	return out;
}
long double gb2(Vector *X){
	long double out = (sqrt( X->v[0]*X->v[0] + sin(X->v[2]) + 1.06) / 9) - 0.1;
	//printf("GB2: %Lf\n",out);
	return out;
}
long double gb3(Vector *X){
	long double out = (-1*exp(-1*X->v[0]*X->v[1]) - ((10*M_PI - 3) / 3) ) / 20;
	//printf("GB3: %Lf\n",out);
	return out;
}
void Problema3b(){
	int n = 3;
	Vector X0 = initialize_Vector(n);
	Vector XoutPuntoFijo = initialize_Vector(n);
	Vector XoutBroyden = initialize_Vector(n);
	Vector XoutNewton = initialize_Vector(n);
	Vector lambda = initialize_Vector(n);
	
	long double (*F[n])(Vector*);
	F[0] = fb1,F[1] = fb2, F[2] = fb3;
	
	long double (*G[n])(Vector*);
	G[0] = gb1,G[1] = gb2 , G[2] = gb3; 
	
	X0.v[0] = 0.1,X0.v[1] = 0.1, X0.v[2] = -0.1;
	int iter_PuntoFijo = NoLinealSolvebyPuntoFijo(&X0,G,&XoutPuntoFijo,1e-6);
	long double error_PF = compute_error(&XoutPuntoFijo, F, n);
	
	X0.v[0] = 0.1,X0.v[1] = 0.1, X0.v[2] = -0.1;
	int iter_Broyden = NoLinealSolvebyBroyden(&X0,F,&XoutBroyden,1e-6,1e-8);
	long double error_B = compute_error(&XoutBroyden, F, n);
	
	X0.v[0] = 0.1,X0.v[1] = 0.1, X0.v[2] = -0.1;
	int iter_Newton = NoLinealSolvebyNewton(&X0,F,&XoutNewton,1e-6,1e-8);
	long double error_N = compute_error(&XoutNewton, F, n);
	
	printf("\n---- PROBLEMA 3b ----\n");
	printf("\n%-11s %-32s %-12s %-12s\n","Metodo", "Resultado (x1,x2,x3)", "     Iteraciones", "	  Error");
    printf("--------------------------------------------------------------------------------\n");
    printf("%-11s (%.8Lf, %.8Lf, %.8Lf)   %-12d %.3Le\n","Punto Fijo",XoutPuntoFijo.v[0], XoutPuntoFijo.v[1], XoutPuntoFijo.v[2],iter_PuntoFijo, error_PF);
    printf("%-11s (%.8Lf, %.8Lf, %.8Lf)   %-12d %.3Le\n","Broyden",XoutBroyden.v[0], XoutBroyden.v[1], XoutBroyden.v[2],iter_Broyden, error_B);
    printf("%-11s (%.8Lf, %.8Lf, %.8Lf)   %-12d %.3Le\n", "Newton",XoutNewton.v[0], XoutNewton.v[1], XoutNewton.v[2],iter_Newton, error_N);

	
	
	free_Vector(&X0);
	free_Vector(&lambda);
	free_Vector(&XoutBroyden);
	free_Vector(&XoutNewton );
	
}


int main(int argc, char **argv)
{
    if (argc == 1) {
		Problema3a();
		Problema3b();
		return 0;
    }
    //CASO INDIVIDUAL
    if (argc == 2) {
        if (strcmp(argv[1], "P3a") == 0) {
            Problema3a();
            return 0;
        }
        if (strcmp(argv[1], "P3b") == 0) {
            Problema3b();
            return 0;
        }

        printf("Uso:\n");
        printf("  %s       -> Ejecuta Problema3a y Problema3b\n", argv[0]);
        printf("  %s P3a   -> Ejecuta solo Problema3a\n", argv[0]);
        printf("  %s P3b   -> Ejecuta solo Problema3b\n", argv[0]);
        return 1;
    }
    printf("Uso:\n");
	printf("  %s        -> Ejecuta Problema3a y Problema3b\n", argv[0]);
	printf("  %s P3a    -> Ejecuta solo Problema3a\n", argv[0]);
	printf("  %s P3b    -> Ejecuta solo Problema3b\n", argv[0]);
	return 1;
}

