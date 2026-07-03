||----        README - AguilarKalebMetodosTarea11     ----||

En este archivo se describe la estructura del código usado para resolver 
los Problemas 3a y 3b mediante Punto Fijo, Newton y Broyden.
El documento contiene funciones para manejo de vectores, matrices,
sistemas no lineales y métodos iterativos.

NO EXISTEN ARCHIVOS ADJUNTOS. TODO ESTÁ CONTENIDO EN UN SOLO CÓDIGO C.


COMPILACIÓN:
		gcc AguilarKalebMetodosTarea11.c -o Programa -lm

USO:
	./Programa             (Ejecuta modo default: Problema 3a y 3b)
	./Programa P3a         (Ejecuta únicamente Problema 3a)
	./Programa P3b         (Ejecuta únicamente Problema 3b)



ARCHIVOS ADJUNTOS:
	•	Aguilar_Kaleb_Metodos_Tarea11.c : Archivo principal con la implementación de todos los métodos.
	•	Matrix.h : Librería para manejo de las estructuras Matrix y Vector.
	•	Solve_by.h : Librería con métodos de resolución de sistemas lineales.



LIBRERÍA Matrix.h:

	MATRIX:
		typedef struct {
    			long double **m; // Matriz 
    			int X; // columnas
    			int Y; // filas
		}Matrix;
	

		Matrix initialize_Matrix(int Y, int X) -> crea una matriz Y×X inicializada en 0.

		Matrix initialize_I(int X, int Y) -> crea matriz identidad X×Y.

		Matrix create_Matrix_from_FILE(FILE *f) -> carga matriz desde archivo.

		Matrix add_Matrix(Matrix *A, Matrix *B) -> suma de matrices.

		Matrix diff_Matrix(Matrix *A, Matrix *B) -> resta de matrices.

		Matrix multiply_Matrix(Matrix *A, Matrix *B) -> multiplicación de matrices.

		void free_Matrix(Matrix *m) -> libera la memoria de la matriz.

		void print_Matrix(Matrix m) -> imprime la matriz en consola.

		void copy_Matrix(Matrix *from, Matrix *to) -> copia una matriz.

		void transpose(Matrix *M, Matrix *T) -> calcula la transpuesta.

		double Matrix_Norm_Frobenius(Matrix *Diff) -> norma de Frobenius.


	VECTOR:
		typedef struct {
    			int n;  // Número de elementos
    			long double *v;  // Vector
		}Vector;
 
		
		Vector initialize_Vector(int n) -> crea un vector de tamaño n.

		Vector create_Vector_from_FILE(FILE *f) -> carga vector desde archivo.

		double Vector_norm(Vector *v) -> norma del vector.

		Vector diff_Vector(Vector *a, Vector *b) -> resta de vectores.

		void diff_Vector_inplace(Vector *a, Vector *b, Vector *out) -> resta in-place.

		void print_Vector(Vector a) -> imprime el vector.

		void save_Vector_to_csv(Vector *v, const char *filename) -> guarda vector en CSV.

		void free_Vector(Vector *v) -> libera memoria del vector.

		void copy_Vector(Vector *from, Vector *to) -> copia un vector.

		void copy_vector(Vector *from, Vector *to) -> copia un vector (alternativa).


	OPERACIONES MATRIZ-VECTOR:
		void Matrix_times_Vector(Matrix *A, Vector *x, Vector *result) -> producto matriz-vector.
		void Vector_row(Matrix *A, int row, Vector *out) -> extrae fila como vector.
		void Vector_col(Matrix *A, int col, Vector *out) -> extrae columna como vector.
		void scalar_times_Vector(double alpha, Vector *v, Vector *out) -> producto por escalar.
		void sum_Vectors(Vector *a, Vector *b, Vector *out) -> suma de vectores.
		double Vector_times_Vector(Vector *u, Vector *v) -> producto punto.


MÉTODOS ITERATIVOS NO LINEALES:

	• Punto Fijo:
		void NoLinealSolvebyPuntoFijo(Vector *X0,void (*G[])(Vector*,Vector*), Vector *Xout)

	• Newton:
		void NoLinealSolvebyNewton(Vector *X0,void (*F)(Vector*,Vector*),void (*J)(Vector*,long double**), Vector *Xout)

	• Broyden:
		void NoLinealSolvebyBroyden(Vector *X0,void (*F)(Vector*,Vector*), void (*J)(Vector*,long double**), Vector *Xout)

FORMATO DE SALIDA EN TERMINAL:

	• Problema 3a:
		- Imprime iteraciones, error, y solución final para cada método.

	• Problema 3b:
		- Imprime iteraciones, error, y solución final para cada método.


	||---------------------------------------------------------------------||

