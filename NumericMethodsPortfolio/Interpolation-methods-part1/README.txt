		||----        README - AguilarKalebMetodosTarea8     ----||

En este archivo se describe la estructura y funcionamiento del programa "AguilarKalebMetodosTarea8.c", 
el cual implementa distintos métodos de Interpolación Numérica en C:

- Interpolación de Taylor  
- Interpolación de Lagrange  
- Interpolación de Neville  
- Interpolación de Newton  

------------------------------------------------------------
ARCHIVOS ADJUNTOS:
------------------------------------------------------------
	•	Aguilar_Kaleb_Metodos_Tarea8_Interpolacion.c : Archivo principal con la implementación completa.  
	•	Taylor.txt : Archivo de entrada para el Problema 1 (coeficientes de derivadas).  
	•	ValoresExp.txt : Archivo de entrada para el Problema 5 (puntos para evaluar e^x).  
	•	IV.txt : Archivo de entrada con valores proporcionados para el Problema 6.  

------------------------------------------------------------
COMPILACIÓN:
------------------------------------------------------------
		gcc AguilarKalebMetodosTarea8.c -o Interpolacion -lm

------------------------------------------------------------
EJECUCIÓN:
------------------------------------------------------------
	./Interpolacion [Pn] [archivo(s)] [parámetros]

	Si no se proporcionan argumentos, el programa ejecuta los problemas (P1a,P1b, P5 y P6)
	con los archivos por defecto: Taylor.txt (P1a y P1b), ValoresExp.txt (P5) e IV.txt (P6).

	Ejemplos:

		./Interpolacion 						(Ejecuta P1a, P1b , P5 y P6 con los archivos por defecto)
		./Interpolacion P1a Taylor.txt
		./Interpolacion P1b Taylor.txt Salida 5 0 2 0.001
		./Interpolacion P5 ValoresExp.txt
		./Interpolacion P6 IV.txt SalidaIV 0 0.025 0.0005


------------------------------------------------------------
PROBLEMAS:
------------------------------------------------------------

	PROBLEMA 1: Interpolación de Taylor  
	
	-	Función principal: Problema1a()  // interactivo
	
	-	Función principal: Problema1b()  // evaluación completa de [0,2] para n = {1,3,5,10}
	
	-	Entrada: archivo con derivadas f[0..n], punto x_0, orden n, rango [a,b], y paso dx.  
	
	-	Salida: archivo con los valores de P(x) en el rango especificado.   


	PROBLEMA 5: Comparación de métodos  
		
	-	Función principal: Problema5()  
	
	-	Entrada: archivo con n puntos (x, y).  
	
	-	Procesa f(x)=eˣ en los puntos y compara los resultados obtenidos por:
		- Lagrange  
		- Neville  
		- Newton  
	
	-	Salida: imprime en terminal una tabla con errores absolutos.  


	PROBLEMA 6: Gráficas de interpolación  

	-	Función principal: Problema6() 
 
	-	Ejecuta simultáneamente los tres métodos de interpolación (Lagrange, Neville y Newton).  
		
	-	Entrada: archivo con puntos (x,y), rango [a,b] y paso dx.  
	
	-	Salida: tres archivos:
			- SalidaIV_Lagrange.txt  
			- SalidaIV_Neville.txt  
			- SalidaIV_Newton.txt  

------------------------------------------------------------
FUNCIONES PRINCIPALES:
------------------------------------------------------------

	-	long double Interpolacion_Taylor(int n, long double x₀, long double x, long double *f)

	-	long double Interpolacion_Lagrange(int n, long double x, long double *X, long double *f)

	-	long double Interpolacion_Neville(int i, int j, long double x, long double *X, long double *f)

	-	long double Interpolacion_Newton(int n, long double x, long double *X, long double *f)

	-	Problema1(), Problema5(), Problema6()  

------------------------------------------------------------
Formato de Archivos:
------------------------------------------------------------

	-	Taylor.txt:
			n
			x_0
			f'[x_0]
			f''[x_0]
			...
			f^n[x_0]

	-	ValoresExp.txt:
			n
			x_0	y_0
			x_1	y_1
			...
			x_n	y_n

	-	IV.txt:
			n
			x_0	y_0
			x_1	y_1
			...
			x_n	y_n

------------------------------------------------------------
SALIDAS:
------------------------------------------------------------
	
	Los resultados se guardan en archivos .txt con el formato:

		# Método de interpolación
		# Archivo base: <archivo>
		# n, rango, dx
		# x			P(x)
		x_0	y_0
		x_1	y_1
		...
		x_n	y_n

------------------------------------------------------------
IMPRESIÓN EN TERMINAL:
------------------------------------------------------------

	-	Para el Problema 1a: Imprime en terminal los valores interpolados según se soliciten.  
	-	Para el Problema 1b: Imprime parámetros de interpolación y el nombre de los archivos generados.  
	-	Para el Problema 5: Tabla con errores comparativos entre métodos.  
	-	Para el Problema 6: Imprime parámetros de interpolación y el nombre de los archivos generados.  

------------------------------------------------------------
NOTAS:
------------------------------------------------------------
	- Todos los cálculos se realizan con tipo **long double** para mayor precisión.  
	- Es necesario incluir la librería **math.h** y compilar con **-lm**.  
	- El rango [a,b] y dx determinan el dominio de evaluación.  

	||---------------------------------------------------------------------||
