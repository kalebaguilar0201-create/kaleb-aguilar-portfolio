#ifndef SOLVEEDO_H
#define SOLVEEDO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"

//EDO 1 DIM
void Solve_EDO(int n_f, long double (*method)(int, long double,long double, long double,  long double (**f) (long double , long double)) ,long double a, long double b, int N, long double alpha, long double (**f) (long double , long double), Vector *x, Vector *y){
	long double h = (b - a)/N;
	x->v[0] = a;
	y->v[0] = alpha;
	for(int i = 1; i < y->n; i++){
		y->v[i] = method(n_f,h,x->v[i-1],y->v[i-1],f);
		x->v[i] = a + i*h;
	}
}
long double EDO_Euler_Method(int n_f, long double h,  long double last_x ,long double last_y, long double (**f) (long double , long double)){
	return last_y + h*f[0](last_x,last_y); //[i][orden]
}
long double EDO_Houden_Method(int n_f, long double h,  long double last_x ,long double last_y, long double (**f) (long double , long double)){
	long double aprox_y = last_y + h*f[0](last_x,last_y);
	return last_y + (h/2)*(f[0](last_x,last_y) + f[0](last_x,aprox_y)) ; //[i][orden]
}
long double EDO_Taylor_2ndOrder_Method(int n_f, long double h,  long double last_x ,long double last_y, long double (**f) (long double , long double)){
	//f[0] := f(x,y)  f[1] := df(x,y)/dx  f[2] := df(x,y)/dy
	if (n_f < 3){
		printf("ERROR en método de Taylor df/dx o df/dy no existentes");
		return nan("");
	}
	return last_y + h*f[0](last_x,last_y) + (h*h/2)*(f[1](last_x,last_y) + f[2](last_x,last_y)*f[0](last_x,last_y));
}
long double EDO_RK4_Method(int n_f, long double h,  long double last_x ,long double last_y, long double (**f) (long double , long double)){
	long double k1 = f[0](last_x,last_y);
	long double k2 = f[0](last_x + (h/2),last_y + (h/2)*k1);
	long double k3 = f[0](last_x + h/2, last_y + h*k2/2);
	long double k4 = f[0](last_x + h,   last_y + h*k3);
	return last_y + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
}

void Graficar_Solved_EDO(const char *filename,const char *method_header,long double a,long double b,Vector *X,Vector *Y){
    FILE *fp = fopen(filename, "w");
    if (!fp){
        printf("ERROR: No se pudo crear el archivo %s\n", filename);
        return;
    }

    fprintf(fp, "# Metodo: %s\n", method_header);
    fprintf(fp, "# Intervalo: a = %.15Lf, b = %.15Lf\n", a, b);
    fprintf(fp, "# N = %d\n", X->n);
    fprintf(fp, "x,y\n");

    for (int i = 0; i < X->n; i++){
        fprintf(fp, "%.15Lf,%.15Lf\n", X->v[i], Y->v[i]);
    }

    fclose(fp);
    printf("Archivo CSV generado: %s\n", filename);
}



//SISTEMA EDO 2X2
typedef void (*step2x2_fn)(
    long double h,
    long double x, long double y,
    long double *x_new, long double *y_new,
    long double (**fx)(long double,long double),
    long double (**fy)(long double,long double)
);
void EDO_System2x2_euler(long double h,long double x, long double y,long double *x_new, long double *y_new,long double (**fx)(long double,long double),long double (**fy)(long double,long double)){
    long double fx0 = fx[0](x,y);
    long double fy0 = fy[0](x,y);
    *x_new = x + h * fx0;
    *y_new = y + h * fy0;
}
void EDO_System2x2_heun(long double h,long double x, long double y,long double *x_new, long double *y_new,long double (**fx)(long double,long double),long double (**fy)(long double,long double)){
    
    long double k1x = fx[0](x,y);
    long double k1y = fy[0](x,y);

    long double xp = x + h * k1x;
    long double yp = y + h * k1y;

    long double k2x = fx[0](xp,yp);
    long double k2y = fy[0](xp,yp);

    *x_new = x + (h/2.0L) * (k1x + k2x);
    *y_new = y + (h/2.0L) * (k1y + k2y);
}
void EDO_System2x2_taylor2(long double h,long double x, long double y,long double *x_new, long double *y_new,long double (**fx)(long double,long double),long double (**fy)(long double,long double)){
    long double fx0 = fx[0](x,y);
    long double fx_x = fx[1](x,y);
    long double fx_y = fx[2](x,y);

    long double fy0 = fy[0](x,y);
    long double fy_x = fy[1](x,y);
    long double fy_y = fy[2](x,y);

    *x_new = x + h * fx0 + (h*h/2) * (fx_x + fx_y * fx0);
    *y_new = y + h * fy0 + (h*h/2) * (fy_x + fy_y * fy0);
}
void EDO_System2x2_RK4(long double h,long double x, long double y,long double *x_new, long double *y_new,long double (**fx)(long double,long double),long double (**fy)(long double,long double)){
    long double k1x = fx[0](x,y);
    long double k1y = fy[0](x,y);

    long double x2 = x + h*k1x/2;
    long double y2 = y + h*k1y/2;
    long double k2x = fx[0](x2,y2);
    long double k2y = fy[0](x2,y2);

    long double x3 = x + h*k2x/2;
    long double y3 = y + h*k2y/2;
    long double k3x = fx[0](x3,y3);
    long double k3y = fy[0](x3,y3);

    long double x4 = x + h*k3x;
    long double y4 = y + h*k3y;
    long double k4x = fx[0](x4,y4);
    long double k4y = fy[0](x4,y4);

    *x_new = x + (h/6) * (k1x + 2*k2x + 2*k3x + k4x);
    *y_new = y + (h/6) * (k1y + 2*k2y + 2*k3y + k4y);
}
void Solve_EDO_System2x2(step2x2_fn step, long double a, long double b, int N, long double alpha_x, long double alpha_y, long double (**fx)(long double,long double),long double (**fy)(long double,long double),Vector *t,Vector *X, Vector *Y){
    long double h = (b - a) / N;
    X->v[0] = alpha_x;
    Y->v[0] = alpha_y;
	t->v[0] = a;
    for(int i = 1; i < N; i++){
        long double x_old = X->v[i-1];
        long double y_old = Y->v[i-1];
        long double x_new, y_new;

        step(h, x_old, y_old, &x_new, &y_new, fx, fy);

        X->v[i] = x_new;
        Y->v[i] = y_new;
        t->v[i] = a + i*h;
    }
}
void Graficar_Solved_EDO_System2x2(
    const char *filename,
    const char *method_header,
    long double a,
    long double b,
    Vector *t,
    Vector *X,
    Vector *Y
){
    FILE *fp = fopen(filename, "w");
    if (!fp){
        printf("ERROR: No se pudo crear el archivo %s\n", filename);
        return;
    }

    fprintf(fp, "# Metodo: %s\n", method_header);
    fprintf(fp, "# Intervalo: a = %.15Lf, b = %.15Lf\n", a, b);
    fprintf(fp, "# N = %d\n", X->n);
    fprintf(fp, "t,X,Y\n");

    for (int i = 0; i < X->n; i++){
        fprintf(fp, "%.15Lf,%.15Lf,%.15Lf\n",t->v[i], X->v[i], Y->v[i]);
    }

    fclose(fp);
    printf("Archivo CSV generado: %s\n", filename);
}


#endif
