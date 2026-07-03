
		||----        README - Aguilar_Kaleb_Metodos_Tarea6      ----||

En este archivo se describe la estructura del documento Aguilar_Kaleb_MetodosTarea6.c
y las librerías Matrix.h y Solve_by.h. En este código se implementan métodos para 
cálculo de eigenvalores/eigenvectores de matrices y resolución de sistemas lineales 
usando técnicas iterativas y factorizaciones directas.



COMPILACIÓN Y EJECUCIÓN:

	Para compilar el programa principal usar:

		gcc AguilarKalebMetodosTarea6.c -o Metodos_Tarea6 -lm

	El ejecutable generado será:   ./Metodos_Tarea6

EJECUCIÓN POR DEFECTO (sin argumentos):

		./Metodos_Tarea6
		
		Esto ejecuta todos los problemas con los archivos de ejemplo:
		   • Eigen_3x3.txt (problemas 1 a 4)
		   • A1.txt y b1.txt (problema 5)

EJECUCIÓN DE PROBLEMAS INDIVIDUALES:

	Formato:
		./Metodos_Tarea6 <Problema> <Archivo(s)>

	Donde <Problema> puede ser:
		P1  -> Método de Potencia
		P2  -> Método de Potencia Inversa
		P3  -> Iteración de Subespacios
		P4  -> Iteración de Subespacios Inverso
		P5  -> Gradiente Conjugado

	Ejemplos:

		./Metodos_Tarea6 P1 Eigen_5x5.txt

		./Metodos_Tarea6 P2 Eigen_50x50.txt

		./Metodos_Tarea6 P3 Eigen_125x125.txt

		./Metodos_Tarea6 P4 Eigen_3x3.txt

		./Metodos_Tarea6 P5 A1.txt b1.txt


ARCHIVOS ADJUNTOS:
	-AguilarKalebMetodosTarea6.c : Archivo principal con implementación de los métodos.
	-Matrix.h : Librería para manejo de estructuras Matrix y Vector, y operaciones básicas.
	-Solve_by.h : Librería con funciones para resolver sistemas lineales.
	-Eigen_3x3.txt, Eigen_5x5.txt, Eigen_50x50.txt, Eigen_125x125.txt : Archivos de entrada para pruebas de eigenvalores.
	-A1.txt, b1.txt : Archivos de entrada para pruebas del método de Gradiente Conjugado.


LIBRERÍA Matrix.h:

	MATRIX:
		typedef struct {
    			long double **m; 
    			int X; // columnas
    			int Y; // filas
		}Matrix;

		Matrix initialize_Matrix(int Y, int X)
		Matrix initialize_I(int X,int Y)
		Matrix create_Matrix_from_FILE(FILE *f)
		Matrix multiply_Matrix(Matrix *A, Matrix *B)
		Matrix add_Matrix(Matrix *A, Matrix *B)
		Matrix diff_Matrix(Matrix *A, Matrix *B)
	
		void diff_Matrix_inplace(Matrix *A, Matrix *B, Matrix *out)
		void free_Matrix(Matrix *m)
		void print_Matrix(Matrix m)
		void copy_Matrix(Matrix *from, Matrix *to)
		void Matrix_times_Vector(Matrix *A, Vector *x, Vector *result)
		void Matrix_times_Matrix(Matrix *A, Matrix *B, Matrix *result)		
		void transpose(Matrix *M, Matrix *T)					
		double Matrix_Norm_Frobenius(Matrix *Diff)				

	VECTOR:
		typedef struct {
    			int n;        
    			long double *v;
		}Vector;

		Vector initialize_Vector(int n)
		Vector create_Vector_from_FILE(FILE *f)
		Vector diff_Vector(Vector *a, Vector *b)
		
		double Vector_norm(Vector *v)
		double Vector_times_Vector(Vector *u, Vector *v)

		void free_Vector(Vector *a)
		void print_Vector(Vector a)
		void save_Vector_to_csv(Vector *v, const char *filename)
		void copy_Vector(Vector *from, Vector *to)
		
		void diff_Vector_inplace(Vector *a, Vector *b, Vector *out)
		void scalar_times_Vector(double alpha, Vector *v, Vector *out)
		void sum_Vectors(Vector *a, Vector *b, Vector *out)
		
		void Vector_row(Matrix *A, int row , Vector *out)
		void Vector_col(Matrix *A, int col, Vector *out)

LIBRERÍA Solve_by.h:

	Funciones principales de resolución de sistemas lineales:

	◦	solve_D(Matrix *m, Vector *b)        					-> Sustitución directa en matriz diagonal.
	◦	solve_L(Matrix *m, Vector *b)        					-> Sustitución progresiva en triangular inferior.
	◦	solve_U(Matrix *m, Vector *b)        					-> Sustitución regresiva en triangular superior.
	◦	reduce_col_from_row(Matrix *m, Vector *b, int row_1, int row_2) 	-> Eliminación Gaussiana.
	◦	swap_rows(Matrix *m, Vector *b, int row_1, int row_2) 			-> Intercambio de filas.
	◦	pivote(Matrix *A, Vector *b, int col) 					-> Pivoteo parcial.
	◦	solve_by_Gauss(Matrix *m, Vector *b) 					-> Método de eliminación Gaussiana.
	◦	Factorize_by_Crouth / solve_by_Crouth 					-> Factorización LU y resolución.
	◦	Factorize_by_Cholensky / solve_by_Cholensky 				-> Factorización de Cholesky.
	◦	solve_by_Jacobi(Matrix *A, Vector *b) 					-> Método iterativo de Jacobi.
	◦	solve_by_Gauss_Seidel(Matrix *A, Vector *b) 				-> Método iterativo de Gauss-Seidel.
	◦	solve_by(Matrix *A, Vector *b, Vector *out)				-> Intenta resolver por Crouth, GS, Jacobi, Cholensky.
	◦	hide_solve_by(Matrix *A, Vector *b, Vector *out) 			-> Versión silenciosa de solve_by.
	◦	solve_by_Gradiente_Conjugado(Matrix *A, Vector *b, Vector *x) 		-> Método de Gradiente Conjugado.
	◦	Normalize(Vector *v) -> Normaliza un vector.

PROBLEMAS IMPLEMENTADOS:

	PROBLEMA 1: MÉTODO DE POTENCIA

		Función principal:

			Potencia_K(Matrix *A, int k, Vector *eigvals, Vector **eigvecs)

		Funciones auxiliares:

			verify_eigen(Matrix *A, Vector *v, double lambda)
				-> Verifica que v sea eigenvector con |Av - λv| < tol.
			Normalize(Vector *v)
				-> Normaliza vector a norma unitaria.

	PROBLEMA 2: MÉTODO DE POTENCIA INVERSA

		Función principal:

			Potencia_Inversa_K(Matrix *A, int k, Vector *eigvals, Vector **eigvecs)

		Funciones auxiliares:

			hide_solve_by(Matrix *A, Vector *b, Vector *out)	-> Resuelve Ax=b con los métodos disponibles (sin imprimir en terminal).
			verify_eigen(Matrix *A, Vector *v, double lambda)	-> Verifica eigenvector.
			Normalize(Vector *v)					-> Normalización.

	PROBLEMA 3: ITERACIÓN DE SUBESPACIOS

		Función principal:

			Iteracion_de_subespacio(Matrix *A, int m, Vector *eigval, Vector **eigvecs)

		Funciones auxiliares:

			max_without_diagonal(Matrix *m)				-> Devuelve el valor máximo fuera de la diagonal (p,q).
			Metodo_Jacobi(Matrix *A, Matrix *I)			-> Método de Jacobi para diagonalización de matrices simétricas.
			normalize_column(Matrix *M, int col)			-> Normaliza columna de una matriz.
			gram_schmidt_modified(Matrix *M)			-> Ortogonalización de columnas por Gram-Schmidt modificado.

	PROBLEMA 4: ITERACIÓN DE SUBESPACIOS INVERSO

		Función principal:

			Iteracion_de_subespacio_inverso(Matrix *A, int m, Vector *eigval, Vector **eigvecs)

		Funciones auxiliares:

			solve_by_QR_inplace(Matrix *A, Vector *b, Vector *out)	-> Resuelve Ax=b mediante factorización QR (in-place).
			solve_Matrix_system(Matrix *A, Matrix *B, Matrix *X)	-> Resuelve múltiples sistemas AX=B por QR.
			Metodo_Jacobi(Matrix *A, Matrix *I)			-> Para reducción de C en subespacios.
			gram_schmidt_modified(Matrix *M)			-> Ortogonalización de columnas.
			normalize_column(Matrix *M, int col)			-> Normaliza columnas.

	PROBLEMA 5: MÉTODO DE GRADIENTE CONJUGADO
		Función principal:
			solve_by_Gradiente_Conjugado(Matrix *A, Vector *b, Vector *x)



IMPRESIÓN EN TERMINAL:
	◦	Se imprimen dimensiones de matrices de entrada y eigenvalores/eigenvectores.
	◦	Se verifica cada eigenvector con la función verify_eigen().
	◦	Se imprime error ||Ax - b|| para gradiente conjugado.
	◦	Se reporta si se exceden iteraciones máximas en métodos iterativos.

	||---------------------------------------------------------------------||
