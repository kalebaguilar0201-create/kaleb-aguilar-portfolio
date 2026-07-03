#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//PROBLEMA 1 HERMITE
typedef struct param_Hermite{
	int n;			//número de puntos
	long double *q; // VECTOR QUE ES LA DIAGONAL Qii;
	long double *z; // VECTOR zi;
}param_Hermite;

param_Hermite P_Hermite(int n,long double *X ,long double *f,long double *df){
	
	param_Hermite out;
	out.n = n;
	long double **Q = (long double **)malloc((2*n+2) * sizeof(long double *));
	for(int i = 0; i <= 2*n+2;i++ ) Q[i] = (long double *)malloc((2*n+2) * sizeof(long double));
	
	out.z = (long double *)malloc((2*n+2) * sizeof(long double));
	out.q = (long double *)malloc((2*n+2) * sizeof(long double));
	
	//STEP1
	for(int i = 0; i < n; i++){
		
		//STEP2
		out.z[2*i] = X[i];
		out.z[2*i+1] = X[i];
		Q[2*i][0] = f[i];
		Q[2*i+1][0] = f[i];
		
		//STEP3
		if(i != 0)
			Q[2*i+1][1] = df[i];
	}
	

	for(int i = 2; i <= 2*n+1; i++){ // columna
		for(int j = 2; j <= i; j++){ //fila
			long double num = Q[i][j-1] -  Q[i-1][j-1];
			long double den = out.z[i] -  out.z[i-j];
			Q[i][j] = num/den;
		}
	}
	

	for(int i = 0; i <= 2*n+1 ; i++){
		out.q[i] = Q[i][i];
		free(Q[i]);
	}
	
	free(Q);
	return out;
}

long double Interpolacion_Hermite(long double x,param_Hermite H){
	long double out = H.q[0];
	long double aux = 1;
	for(int k = 1; k <= 2*H.n +1 ; k++){
		aux *= x - H.z[k];
		out += H.q[k]*aux;
	}	
	return out;
}

void GraficarInterHermite(char *filename, char *outname, long double a, long double b, long double dx) {
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
    long double *df =(long double *) malloc(n * sizeof(long double));

    for (int i = 0; i < n; i++) {
        if (fscanf(fp, "%Lf %Lf %Lf", &X[i], &f[i], &df[i]) != 3) {
            printf("Error al leer línea %d.\n", i + 1);
            free(X); free(f); free(df);
            fclose(fp);
            return;
        }
    }
    fclose(fp);

    param_Hermite H = P_Hermite(n, X, f, df);
	
    FILE *out = fopen(outname, "w");
    if (!out) {
        printf("Error al crear archivo '%s'\n", outname);
        free(X); free(f); free(df);
        free(H.z); free(H.q);
        return;
    }

    fprintf(out, "# Interpolación de Hermite\n");
    fprintf(out, "# n = %d, dx = %.10Lf, rango = [%.10Lf, %.10Lf]\n", n, dx, a, b);
    fprintf(out, "# x\t\tP(x)\n");

    for (long double x = a; x <= b ; x += dx) {
        long double Px = Interpolacion_Hermite(x, H);
        fprintf(out, "%.10Lf\t%.10Lf\n", x, Px);
    }

    fclose(out);

    printf("\n Evaluación completada y guardada correctamente en '%s'\n", outname);

    //FREE
    free(X); 
    free(f); 
    free(df);
    free(H.z); 
    free(H.q);
}

//PROBLEMA 2 SPLICE CÚBICO NATURAL
typedef struct param_SpliceCubico{
	int n;			//número de puntos
	long double *a; 
	long double *b; 
	long double *c; 
	long double *d; 
}param_SpliceCubico;

param_SpliceCubico P_SpliceCubico(int n, long double *X, long double *f){
	
	//MEMORIA
	long double *h =  (long double *)malloc((n+1) * sizeof(long double ));
	long double *alpha =  (long double *)malloc((n+1) * sizeof(long double ));
	long double *l =  (long double *)malloc((n+1) * sizeof(long double ));
	long double *m = (long double *)malloc((n+1) * sizeof(long double ));
	long double *z = (long double *)malloc((n+1) * sizeof(long double ));

	param_SpliceCubico out;
	
	out.n = n;
	out.a = (long double *)malloc((n+1) * sizeof(long double ));
	out.b = (long double *)malloc((n+1) * sizeof(long double ));
	out.c = (long double *)malloc((n+1) * sizeof(long double ));
	out.d = (long double *)malloc((n+1) * sizeof(long double ));
	
	//CONIDCIONES INICIALES
	out.a[0] = f[0];
	out.a[1] = f[1];
	
	for(int i = 0; i <= n; i++){
		out.a[i] = f[i];
		
		if(i < n){
			h[i] = X[i+1] - X[i];
		}
		
		if(i > 0){
			alpha[i] = 3.0L*(f[i+1] - f[i])/h[i] - 3*(f[i] - f[i-1])/h[i-1];
		}
	}
	
	//CONDICIONES DE FRONTERA
	l[0] = 1, m[0] = 0, z[0] = 0;
	l[n] = 1, m[n] = 0, out.c[n] = 0;
	
	//L,M,Z
	for(int i =1; i < n; i++){
		l[i] = 2*(X[i+1] - X[i-1]) - h[i-1]*m[i-1];
		m[i] = h[i]/l[i];
		z[i] = (alpha[i] - h[i-1]*z[i-1])/l[i];
	}
	
	for(int j = n- 1 ; j >= 0; j--){
		out.c[j] = z[j] - m[j]*out.c[j+1]; 
		out.b[j] = (out.a[j+1] - out.a[j])/h[j] - h[j]*(out.c[j+1] + 2*out.c[j])/3;
		out.d[j] = (out.c[j+1] - out.c[j])/(3*h[j]);
	}

	//FREE
	free(h);
	free(alpha);
	free(l);
	free(m);
	free(z);
	return out;
}

long double Interpolacion_SpliceCubico(long double x,long double *X,param_SpliceCubico S){
	for(int i = 0; i < S.n; i++){
		if(x < X[i+1]){
			return S.a[i] + S.b[i]*(x-X[i]) + S.c[i]*pow(x-X[i],2) + S.d[i]*pow(x-X[i],3);
		}
	}
	return S.a[S.n-1] + S.b[S.n-1]*(x-X[S.n-1]) + S.c[S.n-1]*pow(x-X[S.n-1],2) + S.d[S.n-1]*pow(x-X[S.n-1],3);;
}

void GraficarInterSpliceCubico(char *filename, char *outname, long double a, long double b, long double dx) {
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

    long double *X = (long double *)malloc((n+1) * sizeof(long double));
    long double *f = (long double *)malloc((n+1) * sizeof(long double));
    if (!X || !f) {
        printf("Error de memoria.\n");
        fclose(fp);
        return;
    }

    for (int i = 0; i < n; i++) {
        if (fscanf(fp, "%Lf %Lf", &X[i], &f[i]) != 2) {
            printf("Error al leer línea %d.\n", i + 1);
            free(X);
            free(f);
            fclose(fp);
            return;
        }
    }
    fclose(fp);

    // Calcular parámetros del spline cúbico
    param_SpliceCubico S = P_SpliceCubico(n, X, f);

    FILE *out = fopen(outname, "w");
    if (!out) {
        printf("Error al crear archivo '%s'\n", outname);
        free(X); free(f);
        free(S.a); free(S.b); free(S.c); free(S.d);
        return;
    }

    fprintf(out, "# Interpolación Spline Cúbico Natural\n");
    fprintf(out, "# n = %d, dx = %.10Lf, rango = [%.10Lf, %.10Lf]\n", n, dx, a, b);
    fprintf(out, "# x\t\tS(x)\n");

    for (long double x = a; x <= b; x += dx) {
        long double Sx = Interpolacion_SpliceCubico(x, X, S);
        fprintf(out, "%.10Lf\t%.10Lf\n", x, Sx);
    }

    fclose(out);

    printf("\n Evaluación spline completada y guardada en '%s'\n", outname);

    // Liberar memoria
    free(X); free(f);
    free(S.a); free(S.b); free(S.c); free(S.d);
}

//PROBLEMA3 SPLICE CÚBICO FIJO
param_SpliceCubico P_SpliceCubicoFijo(int n, long double *X, long double *f, long double FP0, long double FPN){
	
	//MEMORIA
	long double *h =  (long double *)malloc((n+1) * sizeof(long double ));
	long double *alpha =  (long double *)malloc((n+1) * sizeof(long double ));
	long double *l =  (long double *)malloc((n+1) * sizeof(long double ));
	long double *m = (long double *)malloc((n+1) * sizeof(long double ));
	long double *z = (long double *)malloc((n+1) * sizeof(long double ));

	param_SpliceCubico out;
	
	out.n = n;
	out.a = (long double *)malloc((n+1) * sizeof(long double ));
	out.b = (long double *)malloc((n+1) * sizeof(long double ));
	out.c = (long double *)malloc((n+1) * sizeof(long double ));
	out.d = (long double *)malloc((n+1) * sizeof(long double ));
	
    for (int i = 0; i <= n; i++)
        out.a[i] = f[i];

    for (int i = 0; i < n; i++)
        h[i] = X[i+1] - X[i];

    alpha[0] = 3 * (f[1] - f[0] ) / h[0] - 3 * FP0;
	alpha[n] = 3* FPN - 3 * (f[n] - f[n-1]) / h[n-1];
	
    for (int i = 1; i < n; i++)
        alpha[i] = 3.0L * ((out.a[i+1] - out.a[i])/h[i] - (out.a[i] - out.a[i-1])/h[i-1]);

    l[0] = 2.0L * h[0];
    m[0] = 0.5L;
    z[0] = alpha[0] / l[0];

    for (int i = 1; i < n; i++) {
        l[i] = 2.0L *(X[i+1] - X[i-1]) - h[i-1]* m[i-1];
        m[i] = h[i] /l[i];
        z[i] = (alpha[i] - h[i-1]* z[i-1])/l[i];
    }

    l[n] = h[n-1] * (2.0L - m[n-1]);
    z[n] = (alpha[n] - h[n-1] *z[n-1])/l[n];
    out.c[n] = z[n];

    for (int j = n - 1; j >= 0; j--) {
        out.c[j] = z[j] - m[j] * out.c[j + 1];
        out.b[j] = (out.a[j+1] - out.a[j])/h[j] - h[j] * (out.c[j+1] + 2.0L * out.c[j])/3.0L;
        out.d[j] = (out.c[j+1] - out.c[j])/(3.0L * h[j]);
    }

	//FREE
	free(h);
	free(alpha);
	free(l);
	free(m);
	free(z);
	return out;
}


long double Interpolacion_SpliceCubicoFijo(param_SpliceCubico S, long double x, long double *X) {
    for (int i = 0; i < S.n; i++) {
        if (x <= X[i+1]) {
            return S.a[i] + S.b[i]*(x - X[i])
                 + S.c[i]*pow(x - X[i], 2)
                 + S.d[i]*pow(x - X[i], 3);
        }
    }
    int j = S.n - 1;
    return S.a[j] + S.b[j]*(x - X[j])
         + S.c[j]*pow(x - X[j], 2)
         + S.d[j]*pow(x - X[j], 3);
}


void GraficarInterSpliceCubicoFijo(char *filename, char *outname, long double dx) {
    FILE *fp = fopen(filename, "r");
    if (!fp) { printf("Error: no se pudo abrir '%s'\n", filename); return; }
	
	//M
    int m;
    if (fscanf(fp, "%d", &m) != 1) {
        printf("Error de lectura\n");
        fclose(fp);
        return;
    }
	
	// X , f, FP0, FPN, N_i 
    long double **X = (long double **)malloc(m * sizeof(long double *));
    long double **f = (long double **)malloc(m * sizeof(long double *));
    long double *FPa = (long double *)malloc(m * sizeof(long double));
    long double *FPb = (long double *)malloc(m * sizeof(long double));
    int *Npts = (int *)malloc(m * sizeof(int));
	
    for (int k = 0; k < m; k++) {
        int N;
        long double fp0, fpn;

		//N
        if (fscanf(fp, "%d %Lf %Lf", &N, &fp0, &fpn) != 3) {
            printf("Error de lectura %d\n", k+1);
            fclose(fp);
            return;
        }

        Npts[k] = N;
        FPa[k] = fp0;
        FPb[k] = fpn;

        X[k] = (long double *)malloc(N * sizeof(long double));
        f[k] = (long double *)malloc(N * sizeof(long double));

        for (int i = 0; i < N; i++) {
            if (fscanf(fp, "%Lf %Lf", &X[k][i], &f[k][i]) != 2) {
                printf("Error  de lectura \n");
                fclose(fp);
                return;
            }
        }
    }
    fclose(fp);

    //INTERPOLACION Y OUT
    FILE *out = fopen(outname, "w");
    if (!out) { printf("Error creando '%s'\n", outname); return; }

	fprintf(out, "x P(x)\n");
    for (int k = 0; k < m; k++) {
        int n = Npts[k] - 1;
        param_SpliceCubico S = P_SpliceCubicoFijo(n, X[k], f[k], FPa[k], FPb[k]);

        for (long double xx = X[k][0]; xx <= X[k][n]; xx += dx) {
            long double Sx = Interpolacion_SpliceCubicoFijo(S, xx, X[k]);
            fprintf(out, "%.12Lf %.12Lf \n", xx, Sx);
        }
        free(S.a); 
        free(S.b); 
        free(S.c); 
        free(S.d);
    }

    fclose(out);
    printf("Evaluación  spline fija  completada en '%s'\n", outname);

    for (int k = 0; k < m; k++) { 
		free(X[k]); 
		free(f[k]); 
	}
    
    free(X); 
    free(f); 
    free(FPa); 
    free(FPb); 
    free(Npts);
}




int main(int argc, char **argv)
{
	//DEFAULT
	if (argc == 1){
		printf("Ejecución por defecto:\n");

		printf("\n\n---------    PROBLEMA 1 -----------\n\n");
		char filename1b[100] = "P1b.txt";
		char filename1c[100] = "P1c.txt";
		char filename1p[100] = "P1propuesto.txt";
		char outname1b[100] = "salidaP1b.txt";
		char outname1c[100] = "salidaP1c.txt";
		char outname1p[100] = "salidaP1propuesto.txt";

		GraficarInterHermite(filename1b, outname1b, 0.3, 0.355, 0.001);
		GraficarInterHermite(filename1c, outname1c, 0.3, 0.355, 0.001);
		GraficarInterHermite(filename1p, outname1p, 0.3, 0.355, 0.001);

		printf("\n\n---------    PROBLEMA 2 -----------\n\n");
		char filename2[100] = "Pato.txt";
		char outname2[100] = "salidaPato.txt";
		GraficarInterSpliceCubico(filename2, outname2, 1.0, 13.4, 0.1);

		printf("\n\n---------    PROBLEMA 3 -----------\n\n");
		char filename3[100] = "Snoopy.txt";
		char outname3[100] = "salidaSnoopy.txt";
		GraficarInterSpliceCubicoFijo(filename3, outname3, 0.1);

		return 0;
	}

    if (argc >= 3)
    {
        char *problema = argv[1];

        // P1 [filename] [outname] [a] [b] [dx]
        if (strcmp(problema, "P1") == 0 && argc == 7)
        {
            char *filename = argv[2];
            char *outname = argv[3];
            long double a = atof(argv[4]);
            long double b = atof(argv[5]);
            long double dx = atof(argv[6]);
            GraficarInterHermite(filename, outname, a, b, dx);
            return 0;
        }
        // P2 [filename] [outname] [a] [b] [dx]
        else if (strcmp(problema, "P2") == 0 && argc == 7)
        {
            char *filename = argv[2];
            char *outname = argv[3];
            long double a = atof(argv[4]);
            long double b = atof(argv[5]);
            long double dx = atof(argv[6]);
            GraficarInterSpliceCubico(filename, outname, a, b, dx);
            return 0;
        }
        // P3 [filename] [outname] [dx]
        else if (strcmp(problema, "P3") == 0 && argc == 5)
        {
            char *filename = argv[2];
            char *outname = argv[3];
            long double dx = atof(argv[4]);
            GraficarInterSpliceCubicoFijo(filename, outname, dx);
            return 0;
        }
        else{
			printf("Uso incorrecto.\n");
			printf("Formas válidas:\n");
			printf("  ./Programa (default)\n");
			printf("  ./Programa P1 [filename] [outname] [a] [b] [dx]\n");
			printf("  ./Programa P2 [filename] [outname] [a] [b] [dx]\n");
			printf("  ./Programa P3 [filename] [outname] [dx]\n");
			return 1;
		}
    }
    else
	{
		printf("Uso incorrecto.\n");
		printf("Formas válidas:\n");
		printf("  ./Programa (default)\n");
		printf("  ./Programa P1 [filename] [outname] [a] [b] [dx]\n");
		printf("  ./Programa P2 [filename] [outname] [a] [b] [dx]\n");
		printf("  ./Programa P3 [filename] [outname] [dx]\n");
		return 1;
	}
    return 1;
}
