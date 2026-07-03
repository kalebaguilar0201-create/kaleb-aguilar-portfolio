		||----        README - AguilarKalebMetodosTarea9     ----||

En este archivo se describe la estructura y funcionamiento del programa "AguilarKalebMetodosTarea9.c", 
el cual implementa distintos métodos de Interpolación Numérica en C:

- Interpolación de Hermite  
- Interpolación Spline Cúbico Natural  
- Interpolación Spline Cúbico con Derivadas Fijas  

------------------------------------------------------------
ARCHIVOS ADJUNTOS:
------------------------------------------------------------
	•	Aguilar_Kaleb_Metodos_Tarea9_Interpolacion.c : Archivo principal con la implementación completa.  
	•	P1a.txt : Archivo de entrada para el Problema 1 (Interpolación de Hermite - caso 1).  
	•	P1b.txt : Archivo de entrada para el Problema 1 (Interpolación de Hermite - caso 2).  
	•	P1propuesto.txt : Archivo de entrada para el Problema 1 (Interpolación de Hermite - caso propuesto).  
	•	Pato.txt : Archivo de entrada para el Problema 2 (Spline Cúbico Natural).  
	•	Snoopy.txt : Archivo de entrada para el Problema 3 (Spline Cúbico Fijo).  

------------------------------------------------------------
COMPILACIÓN:
------------------------------------------------------------
		gcc AguilarKalebMetodosTarea9.c -o Interpolacion -lm

------------------------------------------------------------
EJECUCIÓN:
------------------------------------------------------------
	./Interpolacion [Pn] [archivo(s)] [parámetros]

	Si no se proporcionan argumentos, el programa ejecuta los problemas (P1, P2 y P3)
	con los archivos por defecto: P1a.txt, P1b.txt, P1propuesto.txt, Pato.txt y Snoopy.txt.

	Ejemplos:

		./Interpolacion 							(Ejecuta P1, P2 y P3 con los archivos por defecto)
		./Interpolacion P1 P1a.txt salidaP1a.txt 0.3 0.355 0.001
		./Interpolacion P2 Pato.txt salidaPato.txt 1.0 13.4 0.1
		./Interpolacion P3 Snoopy.txt salidaSnoopy.txt 0.1


------------------------------------------------------------
PROBLEMAS:
------------------------------------------------------------

	PROBLEMA 1: Interpolación de Hermite  
	
	-	Función principal: GraficarInterHermite()  
	
	-	Entrada: archivo con n puntos (x, f(x), f’(x)).  
	
	-	Procesa los datos duplicando nodos y genera una tabla de diferencias divididas modificadas.  
	
	-	Salida: archivo con los valores interpolados P(x) en el rango especificado.  


	PROBLEMA 2: Interpolación Spline Cúbico Natural  
	
	-	Función principal: GraficarInterSpliceCubico()  
	
	-	Entrada: archivo con n puntos (x, y).  
	
	-	Asume condiciones naturales en los extremos (segunda derivada nula).  
	
	-	Salida: archivo con los valores de S(x) evaluados en [a, b] con paso dx.  


	PROBLEMA 3: Interpolación Spline Cúbico con Derivadas Fijas  

	-	Función principal: GraficarInterSpliceCubicoFijo()  
	
	-	Entrada: archivo con m curvas, cada una definida por:
			n FP0 FPN
			x_0 f(x_0)
			x_1 f(x_1)
			...
			x_n f(x_n)
	
	-	Evalúa cada curva usando condiciones de frontera dadas (f’(x_0)=FP0, f’(x_n)=FPN).  
	
	-	Salida: archivo con valores de P(x) concatenados para cada curva.  


------------------------------------------------------------
FUNCIONES PRINCIPALES:
------------------------------------------------------------

	-	param_Hermite P_Hermite(int n, long double *X, long double *f, long double *df)
	-	long double Interpolacion_Hermite(long double x, param_Hermite H)
	-	void GraficarInterHermite(char *filename, char *outname, long double a, long double b, long double dx)

	-	param_SpliceCubico P_SpliceCubico(int n, long double *X, long double *f)
	-	long double Interpolacion_SpliceCubico(long double x, long double *X, param_SpliceCubico S)
	-	void GraficarInterSpliceCubico(char *filename, char *outname, long double a, long double b, long double dx)

	-	param_SpliceCubico P_SpliceCubicoFijo(int n, long double *X, long double *f, long double FP0, long double FPN)
	-	long double Interpolacion_SpliceCubicoFijo(param_SpliceCubico S, long double x, long double *X)
	-	void GraficarInterSpliceCubicoFijo(char *filename, char *outname, long double dx)

------------------------------------------------------------
Formato de Archivos:
------------------------------------------------------------

	-	P1a.txt / P1b.txt / P1propuesto.txt:
			n
			x_0	f(x_0)	f’(x_0)
			x_1	f(x_1)	f’(x_1)
			...
			x_n	f(x_n)	f’(x_n)

	-	Pato.txt:
			n
			x_0	f(x_0)
			x_1	f(x_1)
			...
			x_n	f(x_n)

	-	Snoopy.txt:
			m
			n_1 FP0_1 FPN_1
			x_0 f(x_0)
			x_1 f(x_1)
			...
			x_n_1 f(x_n_1)
			n_2 FP0_2 FPN_2
			x_0 f(x_0)
			x_1 f(x_1)
			...
			x_n_2 f(x_n_2)

			(repite para la m curvas)


------------------------------------------------------------
SALIDAS:
------------------------------------------------------------
	
	Los resultados se guardan en archivos .txt con el formato:

		# Método de interpolación
		# n, rango, dx
		x   P(x)
		x_0 P(x_0)
		x_1 P(x_1)
		...
		x_n P(x_n)

------------------------------------------------------------
IMPRESIÓN EN TERMINAL:
------------------------------------------------------------

	-	Para el Problema 1: Imprime los parámetros de interpolación y confirma la creación del archivo.  
	-	Para el Problema 2: Imprime los parámetros del spline natural y confirma la salida generada.  
	-	Para el Problema 3: Imprime la información de cada curva evaluada y el nombre del archivo final.  

------------------------------------------------------------
NOTAS:
------------------------------------------------------------
	- Todos los cálculos se realizan con tipo **long double** para mayor precisión.  
	- Es necesario incluir la librería **math.h** y compilar con **-lm**.  
	- El rango [a,b] y dx determinan el dominio de evaluación.  
	- El programa soporta ejecución automática o personalizada mediante argumentos.  

	||---------------------------------------------------------------------||
