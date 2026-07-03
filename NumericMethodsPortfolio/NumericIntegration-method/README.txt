		||----        README - Aguilar_Kaleb_Metodos_Tarea10      ----||

En este archivo se describe la estructura del documento Aguilar_Kaleb_Metodos_Tarea10.c.
En este código se implementan métodos de interpolación por mínimos cuadrados y métodos de integración numérica.

ARCHIVOS ADJUNTOS:
	•	Aguilar_Kaleb_Metodos_Tarea10.c : Archivo principal con la implementación de todos los métodos.
	•	Matrix.h : Librería para manejo de las estructuras Matrix y Vector.
	•	IntegracionNumerica.h : Librería con métodos de integración numérica.
	•	Solve_by.h : Librería con métodos de resolución de sistemas lineales.
	•	ValoresSutherland1.txt, ValoresSutherland2.txt (Archivos de entrada para Problema 2)

COMPILACIÓN:
		gcc Aguilar_Kaleb_Metodos_Tarea10.c -o Metodos_Tarea10 -lm

USO:
	./Metodos_Tarea10                    (Ejecuta modo default - Problemas 2 y 3)
	./Metodos_Tarea10 P2 [Cos|Polinomial|RBF] [lambda] [FILENAME] [OUTNAME] [m] [a] [b] [dx]
	./Metodos_Tarea10 P3 [f1|f2|f3] [NCAbierto|NCCerrado|Gaussiana|ALL] [n]

	Ejemplos:
		./Metodos_Tarea10                    (modo default completo)
		./Metodos_Tarea10 P2 Polinomial 0 ValoresSutherland1.txt Salida_Polinomial_L0_Sutherland1.txt 3 0.273 1.6 0.05
		./Metodos_Tarea10 P3 f1 ALL
		./Metodos_Tarea10 P3 f2 Gaussiana 3

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


LIBRERÍA IntegracionNumerica.h:

	CONSTANTES:
		• C[][5]: Coeficientes para cuadratura gaussiana (n=1 a 5)
		• R[][5]: Raices para cuadratura gaussiana (n=1 a 5)

	FUNCIONES DE INTEGRACIÓN:

	•	long double Newton_Cotes_Abierto(int n, Vector *X, long double (*f)(long double), long double a, long double b)
		- Implementa fórmulas de Newton-Cotes abiertas para n=0,1,2,3

	•	long double Newton_Cotes_Cerrado(int n, Vector *X, long double (*f)(long double), long double a, long double b)
		- Implementa fórmulas de Newton-Cotes cerradas para n=1,2,3,4

	•	long double Cuadrática_Gaussiana(int n, Vector *X, long double (*f)(long double), long double a, long double b)
		- Implementa cuadratura gaussiana para n=1 a 5


PROBLEMA 2: Interpolación por Mínimos Cuadrados

	FUNCIONES BASE (\phi):
		• long double phi_Polinomial(long double x, Vector *X, int j)
			- \phi_j(x) = x^j (base polinomial)
		
		• long double phi_Cos(long double x, Vector *X, int j)
			- \phi_j(x) = cos(j*π*x/6) (base coseno)
		
		• long double phi_RBF(long double x, Vector *X, int j)
			- \phi_j(x) = exp(-(x - X_j)^2) (funciones de base radial)

	FUNCIONES PRINCIPALES:
		• void Create_Phi(long double (*phi)(long double, Vector*, int), int m, Vector *X, Matrix *out)
			- Construye la matriz \Phi de funciones base
		
		• void Create_W(int p, int m, long double (*phi)(long double, Vector*, int), Matrix *Lambda, Vector *X, Vector *Y, Vector *W)
			- Calcula los coeficientes W 
		
		• long double Interpolacion_Min_Cuadrados(long double x, long double (*phi)(long double, Vector*, int), Vector *X, Vector *W)
			- Evalúa el polinomio de interpolación en un punto x
		
		• void GraficarInter_Min_Cuadrados(int m, const char *filename, const char *outname, Matrix *Lambda, long double (*phi)(long double, Vector*, int), long double a, long double b, long double dx, char *encabezado)
			- Genera archivo con evaluación del polinomio en un intervalo

	ARCHIVOS DE SALIDA (Formato):
		Salida_{Tipo}_L{Lambda}_Sutherland{Num}.txt
		Ejemplo: Salida_Polinomial_L0_Sutherland1.txt


PROBLEMA 3: Integración Numérica

	FUNCIONES DE PRUEBA:
		• long double f1(long double x) -> sin(x)
		• long double f2(long double x) -> x^2 * ln(x)
		• long double f3(long double x) -> x^2 * exp(-x)

	INTERVALOS DE INTEGRACIÓN:
		• f1: [0, pi/4]
		• f2: [1, 1.5]
		• f3: [0, 1]

	MÉTODOS IMPLEMENTADOS:
		• Newton-Cotes Abierto (n = 0,1,2,3)
		• Newton-Cotes Cerrado (n = 1,2,3,4)
		• Cuadratura Gaussiana (n = 1,2,3,4,5)


FUNCIONES AUXILIARES:

	•	main(int argc, char **argv) -> función principal del programa
		- Modo default: Ejecuta Problemas 2 y 3 completos
		- Modo P2: Interpolación por mínimos cuadrados con parámetros personalizados
		- Modo P3: Integración numérica con parámetros personalizados


SALIDAS EN TERMINAL:

	◦	Problema 2:
		- Evaluación en x=1.2 para diferentes φ y λ
		- Generación de archivos de salida con la interpolación

	◦	Problema 3:
		- Tabla comparativa de métodos de integración
		- Aproximaciones para diferentes valores de n


ARCHIVOS DE ENTRADA POR DEFECTO:
	•	ValoresSutherland1.txt
	•	ValoresSutherland2.txt

ARCHIVOS DE SALIDA POR DEFECTO:
	•	Salida_{Pol|Cos|RBF}_L{0|1e-5|1e-7}_Sutherland{1|2}.txt


	||---------------------------------------------------------------------||