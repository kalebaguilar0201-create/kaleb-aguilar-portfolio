||----        README - AguilarKalebMetodosTarea12     ----||

En este archivo se describe la estructura del código usado para resolver 
los Problemas P1,P2 y P3 mediante los métodos de Euler,Heun,Taylor y RK4 para solucionar 
de Ecuaciones Diferenciales Ordinarías de primer orden.

NO EXISTEN ARCHIVOS ADJUNTOS. TODO ESTÁ CONTENIDO EN UN SOLO CÓDIGO C y sus librerías anexas


COMPILACIÓN:
		gcc AguilarKalebMetodosTarea12.c -o Programa -lm

USO:
	./Programa             (Ejecuta modo default: Problema P1,P2,P3)
	./Programa P1         (Ejecuta únicamente Problema P1)
	./Programa P2         (Ejecuta únicamente Problema P2)
	./Programa P3         (Ejecuta únicamente Problema P3)



ARCHIVOS ADJUNTOS:
	+	Aguilar_Kaleb_Metodos_Tarea11.	c : Archivo principal con la implementación de todos los métodos.
	+	Matrix.h : Librería para manejo de las estructuras Matrix y Vector.
	+	Solve_by.h : Librería con métodos de resolución de sistemas lineales.
	+	SolveEDO.h : Librería con métodos de resolución Ecuaciones Diferenciales.


MAIN:
	+	En el archivo 	AguilarKalebMetodosTarea12.c se muestras las funciones:
		
	P1: 	*f1(long double x, long double y) , *df1dx(long double x, long double y), *df1dy(long double x, long double y)
		
	P2: 	*f2(long double t, long double y) , *f2_gauss(long double t) ,*df2dt(long double t, long double y), *df2dy(long double t, long double y)
		
	P3: 	*f3a(long double x, long double y) , *df3adx(long double x, long double y), *df3ady(long double x, long double y)
		
		*f3b(long double x, long double y) , *df3bdx(long double x, long double y), *df3bdy(long double x, long double y)
		
	RESULTADOS:
			Problema1(), Problema2(),Problema3()



SolveEDO.h
	
	========== EDO ================

		Solve_EDO — Función general que itera en el intervalo [a,b] con N pasos y aplica un método dado.

		EDO_Euler_Method — Método de Euler explícito.

		EDO_Houden_Method — Método de Heun (predictor-corrector de orden 2).

		EDO_Taylor_2ndOrder_Method — Taylor de segundo orden (requiere derivadas parciales precalculadas: f, df/dx, df/dy).

		EDO_RK4_Method — Método de Runge–Kutta de orden 4.
	
		Graficar_Solved_EDO — Imprime solución de EDO .

	
	======= SISTEMA EDO 2X2 =======

		Definición de tipo step2x2_fn para pasos de sistemas 2x2.

		EDO_System2x2_euler, EDO_System2x2_heun, EDO_System2x2_taylor2, EDO_System2x2_RK4 — implementaciones por método.
		
		Solve_EDO_System2x2 — Calculo de solucionar a partir de algún método dado .

		Graficar_Solved_EDO_System2x2 —  Imprime solución del sistema.


Archivos de salida


	PROBLEMA 1:  problema1_euler.csv, problema1_heun.csv,  problema1_taylor2.csv, problema1_rk4.csv

	PROBLEMA 2: problema2_euler.csv, problema2_heun.csv, problema2_taylor2.csv, problema2_rk4.csv, problema2_CuadraturaGaussiana.csv

	PROBLEMA 3: problema3_euler.csv, problema3_heun.csv, problema3_taylor2.csv, problema3_RK4.csv

	
	||---------------------------------------------------------------------||

