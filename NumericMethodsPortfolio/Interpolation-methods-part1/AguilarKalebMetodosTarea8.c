#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


//PROBLEMA 1
long double Interpolacion_Taylor(int n,long double x_0,long double x ,long double *f){
	long double out = f[0];
	long double aux_den = 1;
	long double diff = x-x_0;
	for(int i = 1; i < n ; i++){
		aux_den *= i;
		out += f[i]*pow(diff,i)/aux_den;
	}
	return out;
}

void Problema1a(char *filename) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Error: no se pudo abrir el archivo '%s'.\n", filename);
        return;
    }

    int n;
    long double x_0;

    // Leer n
    if (fscanf(fp, "%d", &n) != 1) {
        printf("Error al leer n.\n");
        fclose(fp);
        return;
    }

    // Leer x_0
    if (fscanf(fp, "%Lf", &x_0) != 1) {
        printf("Error al leer x_0.\n");
        fclose(fp);
        return;
    }
	
	int n_aux;
	while(1){
		printf("Ingrese el grado n que dese utilizar (menor que %d): ",n);
		scanf("%d", &n_aux);
		if (n_aux < n) break;
		else printf("Valor no válido \n");
    }
	n = n_aux;
    // Reservar memoria para f(x_0), f'(x_0), ..., f^n(x_0)
    long double *f = (long double *)malloc((n + 1) * sizeof(long double));
    if (!f) {
        printf("Error de memoria.\n");
        fclose(fp);
        return;
    }

    // Leer los valores de las derivadas
    for (int i = 0; i < n; i++) {
        if (fscanf(fp, "%Lf", &f[i]) != 1) {
            printf("Error al leer f[%d].\n", i);
            free(f);
            fclose(fp);
            return;
        }
    }

    fclose(fp);

    printf("\nDatos cargados correctamente desde '%s'\n", filename);
    printf("n = %d\nx₀ = %.10Lf\n", n, x_0);
    printf("f(x_0) = %.10Lf, f'(x_0) = %.10Lf, ... f^n(x_0) = %.10Lf\n", f[0], f[1], f[n]);
    printf("\nIngrese puntos x para interpolar (escriba 'q' para salir):\n");

    while (1) {
        char entrada[64];
        printf("\nx = ");
        if (scanf("%s", entrada) != 1) break;

        if (entrada[0] == 'q' || entrada[0] == 'Q') {
            printf("Saliendo...\n");
            break;
        }

        long double x;
        if (sscanf(entrada, "%Lf", &x) != 1) {
            printf("Entrada inválida.\n");
            continue;
        }

        long double resultado = Interpolacion_Taylor(n, x_0, x, f);
        printf(">> P(%.6Lf) = %.10Lf\n", x, resultado);
    }

    free(f);
}

void Problema1b(char *filename, char *outname, int n_aux, long double a, long double b, long double dx) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Error: no se pudo abrir el archivo '%s'.\n", filename);
        return;
    }

    int n;
    long double x_0;

    if (fscanf(fp, "%d", &n) != 1) {
        printf("Error al leer n.\n");
        fclose(fp);
        return;
    }

    if (fscanf(fp, "%Lf", &x_0) != 1) {
        printf("Error al leer x_0.\n");
        fclose(fp);
        return;
    }

    long double *f = (long double *)malloc((n + 1) * sizeof(long double));
    if (!f) {
        printf("Error de memoria.\n");
        fclose(fp);
        return;
    }

    for (int i = 0; i < n; i++) {
        if (fscanf(fp, "%Lf", &f[i]) != 1) {
            printf("Error al leer f[%d].\n", i);
            free(f);
            fclose(fp);
            return;
        }
    }

    fclose(fp);
	
    FILE *out = fopen(outname, "w");
    if (!out) {
        printf("Error al crear archivo de salida '%s'\n", outname);
        free(f);
        return;
    }

    printf("\nDatos cargados correctamente desde '%s'\n", filename);
    
    printf("n = %d\nx_0 = %.10Lf\n", n_aux, x_0);
    printf("Evaluando P(x) con dx = %.6Lf en el rango [%.6Lf, %.6Lf]\n", dx, a, b);
    printf("Guardando resultados en: '%s'\n\n", outname);

    fprintf(out, "# Interpolación de Taylor\n");
    fprintf(out, "# Archivo base: %s\n", filename);
    fprintf(out, "# x0 = %.10Lf, n = %d, dx = %.10Lf, rango = [%.10Lf, %.10Lf]\n", x_0, n_aux, dx, a, b);
    fprintf(out, "# x\t\tP(x)\n");
	
    for (long double x = a; x <= b + 1e-12; x += dx) {
        long double Px = Interpolacion_Taylor(n_aux, x_0, x, f);
        fprintf(out, "%.10Lf\t%.10Lf\n", x, Px);
    }

    fclose(out);
    free(f);

    printf("Evaluación completada y guardada correctamente.\n");
}


//PROBLEMA 2
long double Interpolacion_Lagrange(int n,long double x,long double *X ,long double *f){
	long double out = 0;
	for(int i = 0; i < n; i++){
		long double  num = 1;
		long double den = 1;
		for(int j = 0; j < n; j++){
			if (i == j) continue;
			num *= (x - X[j]);
			den *= (X[i] - X[j]);
			
		}
		out += f[i] * (num/den) ;
	}
	
	return out;
}

//PROBLEMA 3
long double Interpolacion_Neville(int i ,int j ,long double x,long double *X ,long double *f){
	long double out;
	if (i == j) return f[i];
	long double Pinext_j = Interpolacion_Neville(i+1, j ,x,X ,f);
	long double Pi_jprev = Interpolacion_Neville(i, j-1 ,x,X ,f);
	out = (Pinext_j * (x-X[i]) - Pi_jprev * (x-X[j])) / (X[j] - X[i]);
	return out;
}

//PROBLEMA 4
long double Interpolacion_Newton(int n, long double x,long double *X ,long double *f){
	
	long double **m = (long double **)malloc(n * sizeof(long double *));
	for(int i = 0; i < n;i++ ) m[i] = (long double *)malloc(n * sizeof(long double));
	
	for(int i = 0; i < n;i++ ) m[i][0] = f[i];
	
	for(int j = 1; j < n; j++){// columna
		for(int i = 0; i < n - j; i++){//fila
			long double num = m[i+1][j-1] - m[i][j-1];
			long double den = X[i+j] - X[i];
			m[i][j] = num/den;
		}
	}
	
	long double out = 0;
	long double prod_aux = 1;
	for(int i = 0; i < n; i++){
		out += m[0][i]*prod_aux;
		prod_aux *= (x-X[i]);
	}
	for (int i = 0; i < n; i++)
		free(m[i]);
	free(m);

	return out;
}

//IMPRESION DE ARCHIVOS
void  GraficarInterNeville(char *filename, char *outname, long double a, long double b, long double dx) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Error: no se pudo abrir el archivo '%s'.\n", filename);
        return;
    }

    int n;
    if (fscanf(fp, "%d", &n) != 1) {
        printf("Error al leer n.\n");
        fclose(fp);
        return;
    }

    long double *X = (long double *)malloc(n * sizeof(long double));
    long double *f = (long double *)malloc(n * sizeof(long double));
    if (!X || !f) {
        printf("Error de memoria.\n");
        fclose(fp);
        return;
    }

    for (int i = 0; i < n; i++) {
        if (fscanf(fp, "%Lf %Lf", &X[i], &f[i]) != 2) {
            printf("Error al leer par (x%d, y%d).\n", i, i);
            free(X);
            free(f);
            fclose(fp);
            return;
        }
    }

    fclose(fp);

    FILE *out = fopen(outname, "w");
    if (!out) {
        printf("Error al crear archivo de salida '%s'\n", outname);
        free(X);
        free(f);
        return;
    }

    printf("\nDatos cargados correctamente desde '%s'\n", filename);
    printf("n = %d\n", n);
    printf("Evaluando P(x) con dx = %.6Lf en el rango [%.6Lf, %.6Lf]\n", dx, a, b);
    printf("Guardando resultados en: '%s'\n\n", outname);

    fprintf(out, "# Interpolación de Neville\n");
    fprintf(out, "# Archivo base: %s\n", filename);
    fprintf(out, "# n = %d, dx = %.10Lf, rango = [%.10Lf, %.10Lf]\n", n, dx, a, b);
    fprintf(out, "# x\t\tP(x)\n");

    for (long double x = a; x <= b + 1e-12; x += dx) {
        long double Px = Interpolacion_Neville(0,n - 1, x, X, f);
        fprintf(out, "%.10Lf\t%.10Lf\n", x, Px);
    }

    fclose(out);
    free(X);
    free(f);

    printf("Evaluación completada y guardada correctamente.\n");
}

void GraficarInterLagrange(char *filename, char *outname, long double a, long double b, long double dx) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Error: no se pudo abrir el archivo '%s'.\n", filename);
        return;
    }

    int n;
    if (fscanf(fp, "%d", &n) != 1) {
        printf("Error al leer n.\n");
        fclose(fp);
        return;
    }

    long double *X = (long double *)malloc(n * sizeof(long double));
    long double *f = (long double *)malloc(n * sizeof(long double));
    if (!X || !f) {
        printf("Error de memoria.\n");
        fclose(fp);
        return;
    }

    for (int i = 0; i < n; i++) {
        if (fscanf(fp, "%Lf %Lf", &X[i], &f[i]) != 2) {
            printf("Error al leer par (x%d, y%d).\n", i, i);
            free(X);
            free(f);
            fclose(fp);
            return;
        }
    }

    fclose(fp);

    FILE *out = fopen(outname, "w");
    if (!out) {
        printf("Error al crear archivo de salida '%s'\n", outname);
        free(X);
        free(f);
        return;
    }

    printf("\nDatos cargados correctamente desde '%s'\n", filename);
    printf("n = %d\n", n);
    printf("Evaluando P(x) con dx = %.6Lf en el rango [%.6Lf, %.6Lf]\n", dx, a, b);
    printf("Guardando resultados en: '%s'\n\n", outname);

    fprintf(out, "# Interpolación de Lagrange\n");
    fprintf(out, "# Archivo base: %s\n", filename);
    fprintf(out, "# n = %d, dx = %.10Lf, rango = [%.10Lf, %.10Lf]\n", n, dx, a, b);
    fprintf(out, "# x\t\tP(x)\n");

    for (long double x = a; x <= b + 1e-12; x += dx) {
        long double Px = Interpolacion_Lagrange(n, x, X, f);
        fprintf(out, "%.10Lf\t%.10Lf\n", x, Px);
    }

    fclose(out);
    free(X);
    free(f);

    printf("Evaluación completada y guardada correctamente.\n");
}

void GraficarInterNewton(char *filename, char *outname, long double a, long double b, long double dx) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Error: no se pudo abrir el archivo '%s'.\n", filename);
        return;
    }

    int n;
    if (fscanf(fp, "%d", &n) != 1) {
        printf("Error al leer n.\n");
        fclose(fp);
        return;
    }

    long double *X = (long double *)malloc(n * sizeof(long double));
    long double *f = (long double *)malloc(n * sizeof(long double));
    if (!X || !f) {
        printf("Error de memoria.\n");
        fclose(fp);
        return;
    }

    for (int i = 0; i < n; i++) {
        if (fscanf(fp, "%Lf %Lf", &X[i], &f[i]) != 2) {
            printf("Error al leer par (x%d, y%d).\n", i, i);
            free(X);
            free(f);
            fclose(fp);
            return;
        }
    }

    fclose(fp);

    FILE *out = fopen(outname, "w");
    if (!out) {
        printf("Error al crear archivo de salida '%s'\n", outname);
        free(X);
        free(f);
        return;
    }

    printf("\nDatos cargados correctamente desde '%s'\n", filename);
    printf("n = %d\n", n);
    printf("Evaluando P(x) con dx = %.6Lf en el rango [%.6Lf, %.6Lf]\n", dx, a, b);
    printf("Guardando resultados en: '%s'\n\n", outname);

    fprintf(out, "# Interpolación de Newton\n");
    fprintf(out, "# Archivo base: %s\n", filename);
    fprintf(out, "# n = %d, dx = %.10Lf, rango = [%.10Lf, %.10Lf]\n", n, dx, a, b);
    fprintf(out, "# x\t\tP(x)\n");

    for (long double x = a; x <= b + 1e-12; x += dx) {
        long double Px = Interpolacion_Neville(0,n - 1, x, X, f);
        fprintf(out, "%.10Lf\t%.10Lf\n", x, Px);
    }

    fclose(out);
    free(X);
    free(f);

    printf("Evaluación completada y guardada correctamente.\n");
}

//PROBLEMA 5
void Problema5(char *filename) {
	FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Error: no se pudo abrir el archivo '%s'.\n", filename);
        return;
    }

    int n;
    if (fscanf(fp, "%d", &n) != 1) {
        printf("Error al leer n.\n");
        fclose(fp);
        return;
    }

    long double *X = (long double *)malloc(n * sizeof(long double));
    long double *f = (long double *)malloc(n * sizeof(long double));
    if (!X || !f) {
        printf("Error de memoria.\n");
        fclose(fp);
        return;
    }

    for (int i = 0; i < n; i++) {
        if (fscanf(fp, "%Lf %Lf", &X[i], &f[i]) != 2) {
            printf("Error al leer par (x%d, y%d).\n", i, i);
            free(X);
            free(f);
            fclose(fp);
            return;
        }
    }

    fclose(fp);

    for (int i = 0; i < n; i++)
        f[i] = expl(X[i]); // e^x

    //printf("x, f(x)=e^x , Lagrange, Error_L, Neville, Error_N, Newton, Error_New");
    printf("\n%-10s %-15s %-15s %-15s %-15s %-15s %-15s %-15s\n","x", "f(x)=e^x", "Lagrange", "Error_L", "Neville", "Error_N", "Newton", "Error_New");
    printf("----------------------------------------------------------------------------------------------\n");
	
	long double z_values[5] = {0.4,0.8,1.2,1.6,1.9};
    for (int i = 0; i < 5; i++) {
        long double fx = expl(z_values[i]);
        long double PL = Interpolacion_Lagrange(n, z_values[i], X, f);
        long double PN = Interpolacion_Neville(0, n - 1, z_values[i], X, f);
        long double PNew = Interpolacion_Newton(n, z_values[i], X, f);

        long double errL = fabsl(fx - PL);
        long double errN = fabsl(fx - PN);
        long double errNew = fabsl(fx - PNew);

        printf("%-10.2Lf %-15.10Lf %-15.10Lf %-15.10Lf %-15.10Lf %-15.10Lf %-15.10Lf %-15.10Lf\n", z_values[i], fx, PL, errL, PN, errN, PNew,errNew);
    }
}

//PROBLEMA 6
void Problema6(char *filename , char *outname, long double a, long double b, long double dx){
	
	char outLagrange[128], outNeville[128], outNewton[128];
	
	snprintf(outLagrange, sizeof(outLagrange), "%s_Lagrange.txt", outname);
    snprintf(outNeville,  sizeof(outNeville),  "%s_Neville.txt",  outname);
    snprintf(outNewton,   sizeof(outNewton),   "%s_Newton.txt",   outname);
    
    GraficarInterLagrange(filename,outLagrange,a,b,dx);
	GraficarInterNeville(filename,outNeville,a,b,dx);
	GraficarInterNewton(filename,outNewton,a,b,dx);
	
}


int main(int argc, char **argv)
{
    if (argc == 1)
    {
		char filename1[100] = "Taylor.txt";

        char outname1_3[100] = "InterTaylor3";
        char outname1_5[100] = "InterTaylor5";
        char outname1_10[100] = "InterTaylor10";
        
        char filename5[100] = "ValoresExp.txt";
        
        char filename6[100] = "IV.txt";
        char outname6[100] = "SalidaIV";


        printf("Ejecución por defecto:\n");
        
        printf("\n\n---------    PROBLEMA 1a -----------\n\n");
        Problema1a(filename1);
        
        printf("\n\n---------    PROBLEMA 1b -----------\n\n");
        Problema1b(filename1, outname1_3,3, 0, 2, 0.001);
		Problema1b(filename1, outname1_5,5, 0, 2, 0.001);
		Problema1b(filename1, outname1_10,10, 0, 2, 0.001);
		
		
        printf("\n\n---------    PROBLEMA 5 -----------\n\n");
        Problema5(filename5);

        printf("\n\n---------    PROBLEMA 6 -----------\n\n");
        Problema6(filename6, outname6, 0.0, 0.026, 0.0005);

        return 0;
    }

    if (argc >= 3)
    {
        char *problema = argv[1];

        if (strcmp(problema, "P1b") == 0 && argc == 7)
        {
            char *filename = argv[2];
            char *outname = argv[3];
            int n = atof(argv[4]);
            long double a = atof(argv[5]);
            long double b = atof(argv[6]);
            long double dx = atof(argv[7]);
            Problema1b(filename, outname,n, a, b, dx);
        }else if(strcmp(problema, "P1a") == 0 && argc == 3){
			char *filename = argv[2];
            Problema1a(filename);
		}
        else if (strcmp(problema, "P5") == 0 && argc == 3)
        {
            char *filename = argv[2];
            Problema5(filename);
        }
        else if (strcmp(problema, "P6") == 0 && argc == 7)
        {
            char *filename = argv[2];
            char *outname = argv[3];
            long double a = atof(argv[4]);
            long double b = atof(argv[5]);
            long double dx = atof(argv[6]);
            Problema6(filename, outname, a, b, dx);
        }
        else
        {
            printf("Uso incorrecto.\n");
            printf("Formas válidas:\n");
            printf("  ./Programa                  -> ejecuta P1, P5 y P6 con los archivos por defecto\n");
            printf("  ./Programa P1 filename outname n a b dx\n");
            printf("  ./Programa P5 filename\n");
            printf("  ./Programa P6 filename outname a b dx\n");
            return 1;
        }
    }
    else
    {
        printf("Uso incorrecto.\n");
        printf("Formas válidas:\n");
        printf("  ./Programa                  -> ejecuta P1, P5 y P6 por defecto\n");
        printf("  ./Programa P1 filename outname n a b dx\n");
        printf("  ./Programa P5 filename\n");
        printf("  ./Programa P6 filename outname a b dx\n");
        return 1;
    }

    return 0;
}
