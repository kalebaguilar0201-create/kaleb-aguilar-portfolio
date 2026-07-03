
//AguilarKalebMetodosTarea4.c
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include "Matrix.h"
#include "Solve_by.h"


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

void solve_by_Potencia_K(Matrix *A, int k, Vector *eigvals, Vector **eigvecs) {
    int n = A->Y;
	Vector v1 = initialize_Vector(n);
	
    for (int j = 0; j < k; j++) {
        Vector v0 = initialize_Vector(n);
        for (int i = 0; i < n; i++) v0.v[i] = 1.0/sqrt(n);
        double lambda = 1000.0, lambda_old = 0.0;
        int iter = 0, max_iter = 1000000;

        do {
            lambda_old = lambda;

            
            Matrix_times_Vector(A, &v0, &v1);
			
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

            double num = 0.0, den = 0.0;
            for (int i = 0; i < n; i++) {
                num += v1.v[i] * v1.v[i];
                den += v1.v[i] * v0.v[i];
            }
            lambda = num / den;

            Normalize(&v1);
            iter++;

            copy_vector(&v1, &v0);

        } while (fabs(lambda - lambda_old) > 1e-18 && iter < max_iter);
		
        if (iter >= max_iter)
            printf("Límite de iteraciones excedido para eigenvector %d\n", j+1);

        eigvals->v[j] = lambda;
        copy_vector(&v0, eigvecs[j]);
        free_Vector(&v0);
        free_Vector(&v1);
    }
}

void solve_by_Potencia_Inversa_K(Matrix *A, int k, Vector *eigvals, Vector **eigvecs) {
    int n = A->Y;
	Vector v1 = initialize_Vector(n);

    for (int j = 0; j < k; j++) {
        Vector v0 = initialize_Vector(n);
        for (int i = 0; i < n; i++) v0.v[i] = 1.0/sqrt(n);

        double lambda = 0.0, lambda_old = 0.0;
        int iter = 0, max_iter = 1000000;

        do {
            lambda_old = lambda;

            if(solve_by(A, &v0, &v1) != 0){
				printf("No se pudo resolver el sistema en la iteración de potencia inversa.\n");
				break;
			} 
			
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

            double num = 0.0, den = 0.0;
            for (int i = 0; i < n; i++) {
                num += v1.v[i] * v1.v[i];
                den += v1.v[i] * v0.v[i];
            }
            lambda = num / den;

            Normalize(&v1);
            iter++;

            copy_vector(&v1, &v0);

        } while (fabs(lambda - lambda_old) > 1e-18 && iter < max_iter);

        if (iter >= max_iter)
            printf("Límite de iteraciones excedido para eigenvector %d\n", j+1);

        eigvals->v[j] = 1.0 / lambda; 
        copy_vector(&v0, eigvecs[j]);
        free_Vector(&v0);
        free_Vector(&v1);
    }
}

void Problema1() {
	printf("\n-------------PROBLEMA 1: MÉTODO DE POTENCIA------------\n");
    const char *archivos[] =  {"Eigen_3x3.txt", "Eigen_5x5.txt","Eigen_50x50.txt", "Eigen_125x125.txt"};
   
	int k ;
    for(int idx = 0; idx < 4; idx++){
		if (idx < 2){
			k = 2;
		}else{
			k = 7;
		}
        FILE *file = fopen(archivos[idx], "r");
        if (!file) {
            printf("Error al abrir el archivo %s\n", archivos[idx]);
            continue;
        }

        Matrix A = create_Matrix_from_FILE(file, ' ');
        fclose(file);

        int n = A.Y;

        Vector eigvals = initialize_Vector(k);
        Vector *eigvecs[k];
        for (int i = 0; i < k; i++) {
            eigvecs[i] = (Vector *)malloc(sizeof(Vector));
            *eigvecs[i] = initialize_Vector(n);
        }

        printf("\n\n\nArchivo: %s\nDimensiones de A: %d x %d\n", archivos[idx], A.Y, A.X);
        
        solve_by_Potencia_K(&A, k, &eigvals, eigvecs);

        printf("\nEigenvalores:\n");
        print_Vector(eigvals);

        for (int i = 0; i < k; i++) {
            printf("\nEigenvector %d:\n", i+1);
            print_Vector(*eigvecs[i]);
        }

        // VERIFICACION
        for (int i = 0; i < k; i++) {
            printf("\nVerificando eigenvector %d:\n", i+1);
            verify_eigen(&A, eigvecs[i], eigvals.v[i]);
        }

        // FREE
        free_Vector(&eigvals);
        for (int i = 0; i < k; i++) {
            free_Vector(eigvecs[i]);
            free(eigvecs[i]);
        }
        free_Matrix(&A);
    }
}

void Problema2() {
	printf("\n\n-------------PROBLEMA 2; MÉTODO DE POTENCIA INVERSA------------\n");
    const char *archivos[] = {"Eigen_3x3.txt", "Eigen_5x5.txt","Eigen_50x50.txt", "Eigen_125x125.txt"};
    int k ;
    for(int idx = 0; idx < 4; idx++){
		if (idx < 2){
			k = 2;
		}else{
			k = 7;
		}
        FILE *file = fopen(archivos[idx], "r");
        if (!file) {
            printf("Error al abrir el archivo %s\n", archivos[idx]);
            continue;
        }

        Matrix A = create_Matrix_from_FILE(file, ' ');
        fclose(file);

        int n = A.Y;

        Vector eigvals = initialize_Vector(k);
        Vector *eigvecs[k];
        for (int i = 0; i < k; i++) {
            eigvecs[i] = (Vector *)malloc(sizeof(Vector));
            *eigvecs[i] = initialize_Vector(n);
        }

        printf("\n\nArchivo: %s\nDimensiones de A: %d x %d\n", archivos[idx], A.Y, A.X);
        solve_by_Potencia_Inversa_K(&A, k, &eigvals, eigvecs);

        printf("\nEigenvalores:\n");
        print_Vector(eigvals);

        for (int i = 0; i < k; i++) {
            printf("\nEigenvector %d:\n", i+1);
            print_Vector(*eigvecs[i]);
        }

        // VERIFICACION
        for (int i = 0; i < k; i++) {
            printf("\nVerificando eigenvector %d:\n", i+1);
            verify_eigen(&A, eigvecs[i], eigvals.v[i]);
        }

        // FREE
        free_Vector(&eigvals);
        for (int i = 0; i < k; i++) {
            free_Vector(eigvecs[i]);
            free(eigvecs[i]);
        }
        free_Matrix(&A);
    }
}

int main(int argc, char **argv)
{
	Problema1();
	Problema2();
    return 0;
	return 0;
}

