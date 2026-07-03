#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Matrix.h"
#include "IntegracionNumerica.h"
#include "SolveEDO.h"

#define _USE_MATH_DEFINES

// PROBLEMA 1

long double f1(long double x, long double y){
	return y;
}
long double df1dx(long double x, long double y){
	return 0;
}
long double df1dy(long double x, long double y){
	return 1;
}
void Problema1(){
    printf("\n=== PROBLEMA 1 ===\n\n");

    int N = 50;
    long double a = 0;
    long double b = 4;
    long double alpha = 1;

    Vector x = initialize_Vector(N);
    Vector y = initialize_Vector(N);

    long double (*f[3]) (long double, long double);
    f[0] = f1;  
    f[1] = df1dx; 
    f[2] = df1dy;

	//EULER
    Solve_EDO(1, EDO_Euler_Method, a, b, N, alpha, f, &x, &y);
    Graficar_Solved_EDO("problema1_euler.csv","Método de Euler",a, b,&x, &y);

	//HOUDEN
    Solve_EDO(1, EDO_Houden_Method, a, b, N, alpha, f, &x, &y);
    Graficar_Solved_EDO("problema1_heun.csv","Método de Heun/Houden",a, b,&x, &y);
    
	//TAYLOR2
    Solve_EDO(3, EDO_Taylor_2ndOrder_Method, a, b, N, alpha, f, &x, &y);
    Graficar_Solved_EDO("problema1_taylor2.csv","Método de Taylor 2° Orden", a, b,&x, &y);

	//RK4
    Solve_EDO(1, EDO_RK4_Method, a, b, N, alpha, f, &x, &y);
    Graficar_Solved_EDO("problema1_rk4.csv","Método de Runge-Kutta 4",a, b, &x, &y);

	//FREE
    free_Vector(&y);
    free_Vector(&x);
}


// PROBLEMA 2
long double f2(long double t, long double y){
	return sqrt(1 + t*t*t);
}
long double f2_gauss(long double t){
	return sqrt(1 + t*t*t);
}
long double df2dt(long double t, long double y){
	return 3*t*t/(2*f2(t,y));
}
long double df2dy(long double t, long double y){
	return 1;
}
void Problema2(){
	    printf("\n=== PROBLEMA 2 ===\n\n");

    int N = 11;
    long double a = 0;
    long double b = 2.2;
    long double alpha = 0;
	long double h = (b-a)/N;
	
    Vector t = initialize_Vector(N);
    Vector y = initialize_Vector(N);
    
    long double (*f[3]) (long double, long double);
    f[0] = f2;  
    f[1] = df2dt; 
    f[2] = df2dy;
	
	//EULER
    Solve_EDO(1, EDO_Euler_Method, a, b, N, alpha, f, &t, &y);
    Graficar_Solved_EDO("problema2_euler.csv","Método de Euler",a, b,&t, &y);
	
	
	//HOUDEN
    Solve_EDO(1, EDO_Houden_Method, a, b, N, alpha, f, &t, &y);
    Graficar_Solved_EDO("problema2_heun.csv","Método de Heun/Houden x->t",a, b,&t, &y);
    
	//TAYLOR2
    Solve_EDO(3, EDO_Taylor_2ndOrder_Method, a, b, N, alpha, f, &t, &y);
    Graficar_Solved_EDO("problema2_taylor2.csv","Método de Taylor 2° Orden", a, b,&t, &y);

	//RK4
    Solve_EDO(1, EDO_RK4_Method, a, b, N, alpha, f, &t, &y);
    Graficar_Solved_EDO("problema2_rk4.csv","Método de Runge-Kutta 4",a, b, &t, &y);

	//Integracion por Cuadratura Gaussiana
	
	Vector y_gauss = initialize_Vector(N);
	for(int i = 0; i < N; i++){
		long double xi = a + i*h;
		y_gauss.v[i] = Cuadrática_Gaussiana(4,NULL,f2_gauss,a,xi);
	}
	
	Graficar_Solved_EDO("problema2_CuadraturaGaussiana.csv","Integración numérica",a, b,&t, &y_gauss);

	
	//FREE
	
	free_Vector(&y_gauss);
    free_Vector(&y);
    free_Vector(&t);
}

// PROBLEMA 3
long double f3a(long double x, long double y){
	return x*(0.4 - 0.018*y);
}
long double df3adx(long double x, long double y){
	return 0.4 - 0.018*y;
}
long double df3ady(long double x, long double y){
	return -0.018*x;
}
long double f3b(long double x, long double y){
	return (-0.8+ 0.023*x)*y;
}
long double df3bdx(long double x, long double y){
	return 0.023*y;
}
long double df3bdy(long double x, long double y){
	return -0.8 + 0.023*x;
}
void Problema3(){
	printf("\n=== PROBLEMA 3 ===\n\n");
    int N = 200;
    long double a = 0;
    long double b = 25;

    long double alpha_x = 30;
    long double alpha_y = 4;

    Vector t = initialize_Vector(N);
    Vector X = initialize_Vector(N);
    Vector Y = initialize_Vector(N);
	
    long double (*Fx[3])(long double, long double);
	long double (*Fy[3])(long double, long double);
    // Para f3a (ecuación de X)
    Fx[0] = f3a;
    Fx[1] = df3adx;
    Fx[2] = df3ady;
	
	Fy[0] = f3b;
    Fy[1] = df3bdx;
    Fy[2] = df3bdy;
    
    
    //EULER
    Solve_EDO_System2x2(EDO_System2x2_euler, a, b, N,alpha_x,alpha_y, Fx,Fy, &t,&X,&Y);
    Graficar_Solved_EDO_System2x2("problema3_euler.csv","P3 : Método de Euler", a, b, &t, &X,&Y);
    
    //HEUN
    Solve_EDO_System2x2(EDO_System2x2_heun, a, b, N,alpha_x,alpha_y, Fx,Fy, &t,&X,&Y);
    Graficar_Solved_EDO_System2x2("problema3_heun.csv","P3 : Método de Heun", a, b, &t, &X,&Y);
    
    //TAYLOR
    Solve_EDO_System2x2(EDO_System2x2_taylor2, a, b, N,alpha_x,alpha_y, Fx,Fy, &t,&X,&Y);
    Graficar_Solved_EDO_System2x2("problema3_taylor2.csv","P3 : Método de Taylo 2do Orden", a, b, &t, &X,&Y);
    
    //RK4
    Solve_EDO_System2x2(EDO_System2x2_RK4, a, b, N,alpha_x,alpha_y, Fx,Fy, &t,&X,&Y);
    Graficar_Solved_EDO_System2x2("problema3_RK4.csv","P3 : Método de RK4", a, b, &t, &X,&Y);

	//FREE
    free_Vector(&t);
    free_Vector(&X);
    free_Vector(&Y);
}


long double f4a(Vector *X){
	return X->v[0]*(0.4 - 0.018*X->v[1]);
}
long double f4b(Vector *X){
	return (-0.8+ 0.023*X->v[0])*X->v[1];
}


int main(int argc, char **argv)
{
    if (argc == 1) {
		Problema1();
		Problema2();
		Problema3();
		return 0;
    }
    //CASO INDIVIDUAL
    if (argc == 2) {
        if (strcmp(argv[1], "P1") == 0) {
			Problema1();
            return 0;
        }
        if (strcmp(argv[1], "P2") == 0) {
			Problema2();
            return 0;
        }
		if (strcmp(argv[1], "P3") == 0) {
			Problema3();
            return 0;
        }
    }
    printf("Uso:\n");
	printf("  %s        -> Ejecuta P1,P2,P3\n", argv[0]);
	printf("  %s P1    -> Ejecuta solo Problema1\n", argv[0]);
	printf("  %s P2    -> Ejecuta solo Problema2\n", argv[0]); 
	printf("  %s P3    -> Ejecuta solo Problema3\n", argv[0]);
	return 1;
}

