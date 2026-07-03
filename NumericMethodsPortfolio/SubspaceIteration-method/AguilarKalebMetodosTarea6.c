
//AguilarKalebMetodosTarea5.c
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include "Matrix.h"
#include "Solve_by.h"

typedef struct {
	int i;
	int j;
	double value;
}max_value;


//VERIFICACION DE IGENVECTORES
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


//PROBLEMA 1: POTENCIA
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


//PROBLEMA 2: POTENCIA INVERSA
void Potencia_Inversa_K(Matrix *A, int k, Vector *eigvals, Vector **eigvecs) {
    int n = A->Y;
	Vector v1 = initialize_Vector(n);

    for (int j = 0; j < k; j++) {
        Vector v0 = initialize_Vector(n);
        for (int i = 0; i < n; i++) v0.v[i] = 1.0/sqrt(n);

        double lambda = 0.0, lambda_old = 0.0;
        int iter = 0, max_iter = 1000000;

        do {
            lambda_old = lambda;

            if(hide_solve_by(A, &v0, &v1) != 0){
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
       
    }
     free_Vector(&v1);
}


//PROBLEMA 3: ITERACION DE SUBESPACIO
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


//PROBLEMA 4: ITERACION DE SUBESPACIO INVERSO
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
	long double sigma = eigval->v[0]*(0.99);
	
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


//PROBLEMA 5: GRADIENTE CONJUGADO
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


// INPUT-OUTPUT
void Problema1(char *filename) {
    printf("\n-------------PROBLEMA 1: MÉTODO DE POTENCIA------------\n");

    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Error al abrir el archivo %s\n", filename);
        return;
    }

    Matrix A = create_Matrix_from_FILE(file);
    fclose(file);
	
	int n = A.Y;
    int k = (n < 6) ? 2 : 7; 
    
    Vector eigvals = initialize_Vector(k);
    Vector *eigvecs[k];
    for (int i = 0; i < k; i++) {
        eigvecs[i] = (Vector *)malloc(sizeof(Vector));
        *eigvecs[i] = initialize_Vector(n);
    }
	
    printf("\n\nArchivo: %s\nDimensiones de A: %d x %d\n", filename, A.Y, A.X);
    Potencia_K(&A, k, &eigvals, eigvecs);
	
    // VERIFICACIÓN
    for (int i = 0; i < k; i++) {
        printf("\nVerificando eigenvector %d:\n", i + 1);
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
void Problema2(char *filename) {
    printf("\n-------------PROBLEMA 2: MÉTODO DE POTENCIA INVERSA------------\n");

    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Error al abrir el archivo %s\n", filename);
        return;
    }

    Matrix A = create_Matrix_from_FILE(file);
    fclose(file);

    int n = A.Y;
    int k = (n < 6) ? 2 : 7; 
    
    Vector eigvals = initialize_Vector(k);
    Vector *eigvecs[k];
    for (int i = 0; i < k; i++) {
        eigvecs[i] = (Vector *)malloc(sizeof(Vector));
        *eigvecs[i] = initialize_Vector(n);
    }

    printf("\n\nArchivo: %s\nDimensiones de A: %d x %d\n", filename, A.Y, A.X);
    Potencia_Inversa_K(&A, k, &eigvals, eigvecs);
	

    // VERIFICACIÓN
    for (int i = 0; i < k; i++) {
        printf("\nVerificando eigenvector %d:\n", i + 1);
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
void Problema3(char *filename) {
    printf("\n-------------PROBLEMA 3: ITERACIÓN DE SUBESPACIOS------------\n");

    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Error al abrir el archivo %s\n", filename);
        return;
    }

    Matrix A = create_Matrix_from_FILE(file);
    fclose(file);

    int n = A.Y;
    int m = (n < 6) ? 2 : 7; 
    Vector eigvals = initialize_Vector(m);
    Vector *eigvecs[m];
    for (int i = 0; i < m; i++) {
        eigvecs[i] = (Vector *)malloc(sizeof(Vector));
        *eigvecs[i] = initialize_Vector(n);
    }
	
	
    Iteracion_de_subespacio(&A,m, &eigvals, eigvecs);
	
    // VERIFICACIÓN
    for (int i = 0; i < m; i++) {
        printf("\nVerificando eigenvector %d:\n", i + 1);
        verify_eigen(&A, eigvecs[i], eigvals.v[i]);
    }

    // FREE
    free_Vector(&eigvals);
    for (int i = 0; i < m; i++) {
        free_Vector(eigvecs[i]);
        free(eigvecs[i]);
    }
    free_Matrix(&A);
}
void Problema4(char *filename) {
    printf("\n-------------PROBLEMA 4: ITERACIÓN DE SUBESPACIOS MÍNIMOS------------\n");

    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Error al abrir el archivo %s\n", filename);
        return;
    }

    Matrix A = create_Matrix_from_FILE(file);
    fclose(file);

    int n = A.Y;
    int m = (n < 6) ? 2 : 7;

    Vector eigvals = initialize_Vector(m);
    Vector *eigvecs[m];
    for (int i = 0; i < m; i++) {
        eigvecs[i] = (Vector *)malloc(sizeof(Vector));
        *eigvecs[i] = initialize_Vector(n);
    }
    
    // Ejecutamos iteración de subespacios mínimos
    Iteracion_de_subespacio_inverso(&A,m,&eigvals, eigvecs);
	
    // VERIFICACIÓN
    for (int i = 0; i < m; i++) {
        printf("\nVerificando eigenvector %d:\n", i + 1);
        verify_eigen(&A, eigvecs[i], eigvals.v[i]);
    }

    // FREE
    free_Vector(&eigvals);
    for (int i = 0; i < m; i++) {
        free_Vector(eigvecs[i]);
        free(eigvecs[i]);
    }
    free_Matrix(&A);
}
void Problema5(char *filename_A, char *filename_B) {
    printf("\n------------ PROBLEMA 5: MÉTODO DE GRADIENTE CONJUGADO ------------\n");

    // Abrir archivo de A
    FILE *fileA = fopen(filename_A, "r");
    if (!fileA) {
        printf("Error al abrir el archivo %s\n", filename_A);
        return;
    }
    Matrix A = create_Matrix_from_FILE(fileA);
    fclose(fileA);

    // Abrir archivo de b
    FILE *fileB = fopen(filename_B, "r");
    if (!fileB) {
        printf("Error al abrir el archivo %s\n", filename_B);
        free_Matrix(&A);
        return;
    }
    Vector b = create_Vector_from_FILE(fileB);  
    fclose(fileB);

    int n = A.Y;

    Vector x = initialize_Vector(n);
	
    solve_by_Gradiente_Conjugado(&A, &b, &x);
	
	if (A.Y < 10) print_solution(A,b,x);
	
	Vector Ax = initialize_Vector(A.Y);
	Matrix_times_Vector(&A, &x,&Ax);
	Vector error_vector = diff_Vector(&Ax, &b);           
    long double error = Vector_norm(&error_vector);
    printf("||A_small * x_small - b_small|| = %Lf\n", error);
	
	free_Vector(&Ax);
	free_Vector(&error_vector);

    // FREE
    free_Vector(&x);
    free_Vector(&b);
    free_Matrix(&A);
}

//MAIN
int main(int argc, char **argv) {
    if (argc < 2) {
        printf("Ejecutando todos los problemas con archivos por defecto...\n");
        
        // P1-P4 con Eigen_3x3.txt
        
        Problema1("Eigen_3x3.txt");
        Problema2("Eigen_3x3.txt");
        
        Problema3("Eigen_3x3.txt");
        Problema4("Eigen_3x3.txt");

        // P5 con A1.txt y b1.txt
        Problema5("A1.txt","b1.txt");
        
        return 0;
    }

    char *problema = argv[1];

    if (strcmp(problema, "P1") == 0) {
        char *archivo = (argc >= 3) ? argv[2] : "Eigen_3x3.txt";
        Problema1(archivo);
    } 
    else if (strcmp(problema, "P2") == 0) {
        char *archivo = (argc >= 3) ? argv[2] : "Eigen_3x3.txt";
        Problema2(archivo);
    } 
    else if (strcmp(problema, "P3") == 0) {
        char *archivo = (argc >= 3) ? argv[2] : "Eigen_3x3.txt";
        Problema3(archivo);
    } 
    else if (strcmp(problema, "P4") == 0) {
        char *archivo = (argc >= 3) ? argv[2] : "Eigen_3x3.txt";
        Problema4(archivo);
    } 
    else if (strcmp(problema, "P5") == 0) {
        char *archivoA = (argc >= 3) ? argv[2] : "A1.txt";
        char *archivoB = (argc >= 4) ? argv[3] : "b1.txt";
        Problema5(archivoA, archivoB);
    } 
    else {
        printf("Problema no reconocido: %s\n", problema);
        printf("Opciones válidas: P1, P2, P3, P4, P5\n");
        return 1;
    }

    return 0;
}


