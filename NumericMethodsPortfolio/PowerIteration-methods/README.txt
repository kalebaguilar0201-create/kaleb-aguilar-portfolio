

		||----        README - Aguilar_Kaleb_Metodos_Tarea5      ----||


En este archivo se describe la estrutura del documento  Aguilar_Kaleb_MetodosTarea5.c.
En este código se implementan métodos para cálculo de eigenvalores y eigenvectores de matrices usando los métodos de potencia y potencia inversa.

COMPILACIÓN:
		gcc AguilarKalebMetodosTarea5.c -o Metodos_Tarea5 -lm


ARCHIVOS ADJUNTOS:
	•	AguilarKalebMetodosTarea5.c : Archivo principal con la implementación de todos los métodos de eigenvalores y eigenvectores.
	•	Matrix.h : Librería para manejo de las estructuras Matrix y Vector.
	•	Solve_by.h : Librería con funciones de resolución de sistemas lineales: Cholesky, Crouth, Jacobi, Gauss-Seidel.
	•	Eigen_3x3.txt, Eigen_5x5.txt, Eigen_50x50.txt, Eigen_125x125.txt : Archivos de entrada para pruebas de eigenvalores.


LIBRERÍA Matrix.h:

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

LIBRERÍA Solve_by.h:

	Funciones principales de resolución de sistemas lineales:

	◦	Factorize_by_Crouth(Matrix *A, Vector *b, Matrix *L, Matrix *U)

		- Entrada: Matriz A, vector b, matrices L y U para factorización.
		- Salida: 0 si factoriza correctamente, -1 si falla.
		- Lógica: Factorización LU mediante Crouth, con verificación de ceros en pivote.


	◦	solve_by_Crouth(Matrix *A, Vector *b)

		- Entrada: Matriz A y vector b.
		- Salida: Vector x solución de Ax=b.
		- Lógica: Usa Factorize_by_Crouth y resuelve L y U mediante sustitución hacia adelante y hacia atrás.


	◦	Factorize_by_Cholensky(Matrix *A, Vector *b, Matrix *L, Matrix *LT)

		- Entrada: Matriz A SPD, vector b, matrices L y LT.
		- Salida: 0 si factoriza correctamente, -1 si falla.
		- Lógica: Factorización de Cholesky con verificación de positividad.


	◦	solve_by_Cholensky(Matrix *A, Vector *b)

		- Entrada: Matriz A SPD y vector b.
		- Salida: Vector solución x.
		- Lógica: Factoriza con Cholesky y resuelve L y LT.


	◦	solve_by_Jacobi(Matrix *A, Vector *b)

		- Entrada: Matriz A y vector b.
		- Salida: Vector solución x.
		- Lógica: Método iterativo de Jacobi con pivoteo si A[ii]=0, hasta tolerancia 1e-8.


	◦	solve_by_Gauss_Seidel(Matrix *A, Vector *b)

		- Entrada: Matriz A y vector b.
		- Salida: Vector solución x.
		- Lógica: Método iterativo de Gauss-Seidel con pivoteo si A[ii]=0, hasta tolerancia 1e-8.


	◦	solve_by(Matrix *A, Vector *b, Vector *out)

		- Entrada: Matriz A y vector b.
		- Salida: Vector out con solución y código de retorno 0 si tiene éxito, -1 si falla.
		- Lógica: Intenta resolver Ax=b primero por Cholesky, luego Crouth, luego Gauss-Seidel, y finalmente Jacobi. Imprime mensaje de error si ninguna funciona.


PROBLEMA 1: Método de Potencia

	◦	Función: solve_by_Potencia_K(Matrix *A, int k, Vector *eigvals, Vector **eigvecs)
	◦	Entrada: Matriz A y número de eigenvalores k a calcular.
	◦	Salida: Vector eigvals con eigenvalores, array eigvecs con eigenvectores.
	◦	Lógica: Método de potencia clásico para calcular los k mayores eigenvalores y sus eigenvectores. 

PROBLEMA 2: Método de Potencia Inversa

	◦	Función: solve_by_Potencia_Inversa_K(Matrix *A, int k, Vector *eigvals, Vector **eigvecs)
	◦	Entrada: Matriz A y número de eigenvalores k a calcular.
	◦	Salida: Vector eigvals con eigenvalores inversos, array eigvecs con eigenvectores.
	◦	Lógica: Método de potencia inversa para aproximar los k menores eigenvalores de A. 

FUNCIONES AUXILIARES:

	◦	Normalize(Vector *v)
		- Entrada: Vector v.
		- Salida: Vector v normalizado en norma 2.
		- Lógica: Calcula la norma Euclidiana del vector y divide cada componente entre dicha norma, dejando el vector con magnitud igual a 1.

	◦	verify_eigen(Matrix *A, Vector *v, double lambda)
		- Entrada: Matriz A, vector v y valor lambda.
		- Salida: Impresión en consola indicando si v es un eigenvector de A con eigenvalor lambda, evaluando |Av - λv| < tol.
		- Lógica: Multiplica la matriz por el vector, resta lambda*v, calcula la norma de la diferencia y verifica si está por debajo de la tolerancia establecida.


IMPRESIÓN EN TERMINAL:

	◦	Se imprimen dimensiones de las matrices de entrada, eigenvalores y eigenvectores.
	◦	Se imprime la verificación de cada eigenvector con |Av - λv|.

	||---------------------------------------------------------------------||

