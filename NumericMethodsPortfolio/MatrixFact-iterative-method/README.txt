


		||----        README - Aguilar_Kaleb_Metodos_Tarea4      ----||


En este archivo se describe la estructura del documento  Aguilar_Kaleb_Metodos_Tarea4.c.
En este código  se implementan diferentes métodos numéricos para la resolución de sistemas de ecuaciones lineales y la ecuación de calor 1D.
 
COMPILACIÓN:
		gcc Aguilar_Kaleb_Metodos_Tarea4.c -o Metodos_Tarea4 -lm


ARCHIVOS ADJUNTOS:
	•	AguilarKalebMetodosTarea4.c : Archivo principal con la implementación de todos los métodos.
	•	Matrix.h : Librería para manejo de las estructuras Matrix y Vector.
	•	A.txt, b.txt (Entrada de Problema 1)
	•	SPD.txt, b_spd.txt (Entrada de Problema 3)
	•	A_small.txt, b_small.txt, A_big.txt, b_big.txt (Entrada de Problemas 5 y 6)


LIBRERÍA Matrix.h:

	MATRIX:
		typedef struct {
    			double **m; // Matriz 
    			int X; // columnas
    			int Y; // filas
		}Matrix;
	

		Matrix initialize_Matrix(int Y, int X) -> crea una matriz Y×X inicializada en 0.

		Matrix initialize_I(int X) -> crea matriz identidad X×X.

		Matrix create_Matrix_from_FILE(FILE *f, char c) -> carga matriz desde archivo con delimitador c.

		Matrix add_Matrix(Matrix *A, Matrix *B) -> suma de matrices.

		Matrix diff_Matrix(Matrix *A, Matrix *B) -> resta de matrices.

		Matrix multiply_Matrix(Matrix *A, Matrix *B) -> multiplicación de matrices.

		void free_Matrix(Matrix *m) -> libera la memoria de la matriz.

		void print_Matrix(Matrix m) -> imprime la matriz en consola.



	VECTOR:
		typedef struct {
    			int n;  // Número de elemento
    			double *v;  // Vector
		}Vector;
 
		
		Vector initialize_Vector(int n) -> crea un vector de tamaño n.

		Vector create_Vector_from_FILE(FILE *f, char c) -> carga vector desde archivo con delimitador del caracter c.

		double Vector_norm(Vector *v) -> norma del vector.

		Vector diff_Vector(Vector *a, Vector *b) -> resta de vectores.

		void print_Vector(Vector a) -> imprime el vector.

		void save_Vector_to_csv(Vector *v, const char *filename) -> guarda vector en CSV.

		void free_Vector(Vector *v) -> libera memoria del vector.


PROBLEMA 1: Factorización por Crouth

	◦	Funciones: Factorize_by_Crouth, solve_by_Crouth
	◦	Entrada: Matriz A y vector b.
	◦	Salida: Vector solución x y norma del error ||Ax-b||.


PROBLEMA 3: Factorización de Cholensky 

	◦	Funciones: Factorize_by_Cholensky, solve_by_Cholensky
	◦	Entrada: Matriz A SPD y vector b.
	◦	Salida: Vector solución x y norma del error  ||Ax-b||..


PROBLEMA 4: Ecuación de calor 1D

	◦	Funciones: BD_matrix, Factorize_Band_Diagonal_by_Cholensky, eq_calor, Problema4
	◦	Método iterativo con refinamiento de la malla hasta convergencia E_dPhi < 0.001 (La definición de E_dPhi se declara en el reporte).
	◦	Salida: Vector Phi guardado en Phi.csv.


PROBLEMA 5: Método de Jacobi

	◦	Función: solve_by_Jacobi
	◦	Entrada: Matriz A, vector b.
	◦	Salida: Vector solución x y reporte de error ||Ax-b||.

PROBLEMA 6: Método de Gauss-Seidel

	◦	Función: solve_by_Gauss_Seidel
	◦	Entrada: Matriz A, vector b.
	◦	Salida: Vector solución x y reporte de error ||Ax-b||.
	
IMPRESIÓN EN TERMINAL:

	◦	Para problemas 1,3,5,6 se imprimen en consola dimensiones y errores.
	◦	Para problema 4, se imprime Phi en consola y se guarda en Phi.csv.


	||---------------------------------------------------------------------||
