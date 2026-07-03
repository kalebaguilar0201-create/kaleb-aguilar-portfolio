#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"
#include "Solve_by.h"
#include "IntegracionNumerica.h"

#define _USE_MATH_DEFINES

//PROBLEMA 1 Y 2
long double phi_Polinomial(long double x,Vector *X,int j){
	return pow(x,j);
}
long double phi_Cos(long double x, Vector *X,int j){
	return cos(j*M_PI*x/6);
}
long double phi_RBF(long double x, Vector *X,int j){
	return exp(-1 * pow( (x - X->v[j]) , 2) );
}

void Create_Phi(
	long double (*phi) (long double ,Vector *X , int), 
	int m,
	Vector *X,
	Matrix *out
	){
	int p = X->n;
	for(int i = 0; i < p; i++){
		for(int j = 0; j < m; j++){
			out->m[i][j] = phi(X->v[i],X,j);
		}
	}
}	
void Create_W(
	int p, int m, 
	long double (*phi) (long double ,Vector *,int), 
	Matrix *Lambda,Vector *X  ,Vector *Y, Vector *W
	){
	//Phi
	Matrix Phi = initialize_Matrix(p,m);
	Matrix Phi_T = initialize_Matrix(m,p);
	Create_Phi(phi,m,X,&Phi);
	transpose(&Phi,&Phi_T);

	//A = (Phi_T*Phi + Lambda)
	Matrix A = initialize_Matrix(m,m);
	Matrix Phi_T_Phi = initialize_Matrix(m,m);
	Matrix_times_Matrix(&Phi_T,&Phi,&Phi_T_Phi);
	sum_Matrix(&Phi_T_Phi,Lambda,&A);
	
	//F = Phi_T*Y
	Vector F = initialize_Vector(m);
	Matrix_times_Vector(&Phi_T,Y,&F);
	
	// solve A*W = F
	solve_by_QR_inplace(&A,&F,W);
	
	free_Matrix(&Phi);
	free_Matrix(&Phi_T);
	free_Matrix(&A);
	free_Matrix(&Phi_T_Phi);
	
	free_Vector(&F);
}	

long double Interpolacion_Min_Cuadrados(long double x, long double (*phi) (long double ,Vector *,int), Vector *X,Vector *W){
	long double out = 0;
	for(int j = 0; j < W->n; j++){
		out += W->v[j]*phi(x,X,j);
	}
	return out;
}

void GraficarInter_Min_Cuadrados(
	int m,
    const char *filename,
    const char *outname,    
    Matrix *Lambda,    
    long double (*phi)(long double ,Vector *,int),              
    long double a, long double b, long double dx,
    char *encabezado 
){
	
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Error: no se pudo abrir el archivo '%s'.\n", filename);
        return;
    }

    int p;
    if (fscanf(fp, "%d", &p) != 1) {
        printf("Error al leer p.\n");
        fclose(fp);
        return;
    }

    Vector X = initialize_Vector(p);
    Vector Y = initialize_Vector(p);

    for (int i = 0; i < p; i++) {
        if (fscanf(fp, "%Lf %Lf", &X.v[i], &Y.v[i]) != 2) {
            printf("Error al leer línea %d.\n", i + 1);
            free_Vector(&X);
            free_Vector(&Y);
            fclose(fp);
            return;
        }
    }
    fclose(fp);
    	 
    Vector W = initialize_Vector(m);
    Create_W(p,m,phi,Lambda,&X,&Y,&W);
	
    FILE *out = fopen(outname, "w");
    if (!out) {
        printf("Error al crear archivo '%s'\n", outname);
        free_Vector(&X);
        free_Vector(&Y);
        free_Vector(&W);
        return;
    }

    fprintf(out, "# Interpolación por Mínimos Cuadrados\n");
    fprintf(out,"%s\n#",encabezado);
    
    for(int i = 0; i < p; i++){
		 fprintf(out,"a%d: %Lf , ",i,W.v[i]);
	}
    fprintf(out,"\n\n");
    fprintf(out, "# m = %d  ,p = %d,lambda = %Lf, dx = %.10Lf, rango = [%.10Lf, %.10Lf]\n",m, p,Lambda->m[0][0], dx, a, b);
    fprintf(out, "#\n\n");
    fprintf(out, "# x\t\tP(x)\n");
	
    for (long double x = a; x <= b; x += dx) {
        long double fx = Interpolacion_Min_Cuadrados(x,phi,&X,&W);
        fprintf(out, "%.10Lf\t%.10Lf\n", x, fx);
    }

    fclose(out);

    printf("\nEvaluación completada y guardada en '%s'\n", outname);

    free_Vector(&X);
    free_Vector(&Y);
    free_Vector(&W);
}

//PROBLEMA 3
long double f1(long double x){
	return sin(x);
}
long double f2(long double x){
	return pow(x,2)*log(x);
}
long double f3(long double x){
	return pow(x,2)*exp(-1*x);
}

int main(int argc, char **argv)
{
    if (argc == 1) {
        // ============================================
        // MODO DEFAULT
        // ============================================
        
        // === EJERCICIO 2 b) ===
		printf("\n\n=== EJERCICIO 2 b) Evaluación en x = 1.2 ===\n\n");

		const char *archivos_b[] = {"ValoresSutherland1.txt", "ValoresSutherland2.txt"};
		long double lambdas_b[] = {0.0L, 1e-5L, 1e-7L};
		long double x_eval_b = 1.2L;
		int m_b = 3;

		long double (*phis_b[])(long double, Vector *, int) = {phi_Polinomial, phi_Cos, phi_RBF};
		const char *nombres_phi_b[] = {"Pol", "Cos", "RBF"};

		for (int f_b = 0; f_b < 2; f_b++) {
			// Cargar datos desde archivo
			FILE *fp_b = fopen(archivos_b[f_b], "r");
			if (!fp_b) {
				printf("Error al abrir archivo '%s'\n", archivos_b[f_b]);
				continue;
			}

			int p_b;
			fscanf(fp_b, "%d", &p_b);
			Vector X_b = initialize_Vector(p_b);
			Vector Y_b = initialize_Vector(p_b);
			for (int i_b = 0; i_b < p_b; i_b++)
				fscanf(fp_b, "%Lf %Lf", &X_b.v[i_b], &Y_b.v[i_b]);
			fclose(fp_b);

			for (int p_phi_b = 0; p_phi_b < 3; p_phi_b++) {
				for (int l_b = 0; l_b < 3; l_b++) {
					Matrix Lambda_b = initialize_Matrix(m_b, m_b);
					for (int i_b = 0; i_b < m_b; i_b++)
						Lambda_b.m[i_b][i_b] = lambdas_b[l_b];

					Vector W_b = initialize_Vector(m_b);
					Create_W(p_b, m_b, phis_b[p_phi_b], &Lambda_b, &X_b, &Y_b, &W_b);

					long double fx_b = Interpolacion_Min_Cuadrados(x_eval_b, phis_b[p_phi_b], &X_b, &W_b);
					printf("Archivo %d | Phi %s | Lambda = %.0Le -> P(%.2Lf) = %.10Lf\n",
						   f_b + 1, nombres_phi_b[p_phi_b], lambdas_b[l_b], x_eval_b, fx_b);

					free_Vector(&W_b);
					free_Matrix(&Lambda_b);
				}
			}

			free_Vector(&X_b);
			free_Vector(&Y_b);
		}

        
		// === EJERCICIO 2 c) ===
		printf("\n\n=== EJERCICIO 2 c) ===\n\n");

        const char *archivos[] = {"ValoresSutherland1.txt", "ValoresSutherland2.txt"};
        int n_archivos = 2;

        long double lambdas[] = {0.0L, 1e-5L, 1e-7L};
        int n_lambdas = 3;

        //Funciones phi
        long double (*phis[])(long double, Vector *, int) = {
            phi_Polinomial, phi_Cos, phi_RBF
        };
        const char *nombres_phi[] = {"Polinomial", "Cos", "RBF"};
        int n_phi = 3;

        long double a = 0.273L, b = 1.6L, dx = 0.05L;
        int m = 3;

        for (int f = 0; f < n_archivos; f++) {
            for (int p = 0; p < n_phi; p++) {
                for (int l = 0; l < n_lambdas; l++) {
                    Matrix Lambda = initialize_Matrix(m, m);
                    for (int i = 0; i < m; i++)
                        Lambda.m[i][i] = lambdas[l];

                    char outname[200];
                    snprintf(outname, sizeof(outname),"Salida_%s_L%.0Le_Sutherland%d.txt",nombres_phi[p], lambdas[l], f + 1);

                    char encabezado[200];
                    snprintf(encabezado, sizeof(encabezado),"Interpolación %s con Lambda = %.0Le (archivo %d)",nombres_phi[p], lambdas[l], f + 1);

                    GraficarInter_Min_Cuadrados(
                        m,
                        archivos[f],
                        outname,
                        &Lambda,
                        phis[p],
                        a, b, dx,
                        encabezado
                    );
					
                    free_Matrix(&Lambda);
                }
            }
        }

        // === EJERCICIO 3 ===
		printf("\n\n=== EJERCICIO 3 ===\n\n");

		long double a_vals[] = {0.0L, 1.0L, 0.0L};
		long double b_vals[] = {M_PI / 4.0L, 1.5L, 1.0L};

		long double (*funcs[])(long double) = {f1, f2, f3};
		const char *names[] = {
			"f1(x) = sin(x)",
			"f2(x) = x^2 * ln(x)",
			"f3(x) = x^2 * e^{-x}"
		};
		int n_funcs = 3;

		for (int f = 0; f < n_funcs; f++) {
			long double a = a_vals[f];
			long double b = b_vals[f];

			printf("=============================================================\n");
			printf(">>> Integral %d: %s en [%.3Lf, %.3Lf]\n", f + 1, names[f], a, b);
			printf("-------------------------------------------------------------\n");
			printf("%-10s %-10s %-25s\n", "Método", "n", "Aproximación");
			printf("-------------------------------------------------------------\n");

			for (int n = 0; n <= 3; n++) {
				Vector X_abierto = initialize_Vector(n + 2);
				long double h_ab = (b - a) / (n + 2);
				for (int i = 0; i < n + 2; i++)
					X_abierto.v[i] = a + h_ab * (i + 1); 
				long double res = Newton_Cotes_Abierto(n, &X_abierto, funcs[f], a, b);
				printf("%-10s %-10d %-25.10Lf\n", "NC Abierto", n, res);
				free_Vector(&X_abierto);
			}

			for (int n = 1; n <= 4; n++) {
				Vector X_cerrado = initialize_Vector(n + 1);
				long double h_c = (b - a) / n;
				for (int i = 0; i < n + 1; i++)
					X_cerrado.v[i] = a + h_c * i;
				long double res = Newton_Cotes_Cerrado(n, &X_cerrado, funcs[f], a, b);
				printf("%-10s %-10d %-25.10Lf\n", "NC Cerrado", n, res);
				free_Vector(&X_cerrado);
			}

			for (int n = 1; n <= 5; n++) {
				long double res = Cuadrática_Gaussiana(n, NULL, funcs[f], a, b);
				printf("%-10s %-10d %-25.10Lf\n", "Gaussiana", n, res);
			}

			printf("=============================================================\n\n");
		}

        return 0;
    }

    // ============================================
    // MODO MANUAL P2
    // ============================================
    else if (strcmp(argv[1], "P2") == 0) {
        if (argc !=10) {
            fprintf(stderr, "Uso: ./Programa P2 [Cos|Polinomial|RBF] [lambda] [FILENAME] [OUTNAME] [m] [a] [b] [dx]\n");
            return 1;
        }

        const char *phi_name = argv[2];
        long double lambda = strtold(argv[3], NULL);
        const char *filename = argv[4];
        const char *outname = argv[5];
        int m = atoi(argv[6]);           
		long double a = strtold(argv[7], NULL);  
		long double b = strtold(argv[8], NULL);
		long double dx = strtold(argv[9], NULL);
		
        long double (*phi)(long double, Vector *, int) = NULL;
        if (strcmp(phi_name, "Cos") == 0)
            phi = phi_Cos;
        else if (strcmp(phi_name, "Polinomial") == 0)
            phi = phi_Polinomial;
        else if (strcmp(phi_name, "RBF") == 0)
            phi = phi_RBF;
        else {
            fprintf(stderr, "Error: función φ desconocida. Use 'Cos', 'Polinomial' o 'RBF'.\n");
            return 1;
        }

        Matrix Lambda = initialize_Matrix(m, m);
        for (int i = 0; i < m; i++)
            Lambda.m[i][i] = lambda;

        char encabezado[200];
        snprintf(encabezado, sizeof(encabezado),
                 "Interpolación %s con Lambda = %.0Le (manual)", phi_name, lambda);

        GraficarInter_Min_Cuadrados(
            m, filename, outname, &Lambda, phi, a, b, dx, encabezado);

        free_Matrix(&Lambda);
        return 0;
    }
	// ============================================
    // MODO MANUAL P3
    // ============================================
    else if (strcmp(argv[1], "P3") == 0) {
        if (argc < 3) {
            fprintf(stderr, "Uso: ./Programa P3 [f1|f2|f3] [NCAbierto|NCCerrado|Gaussiana|ALL] [n]\n");
            return 1;
        }

        const char *f_name = argv[2];            
        const char *method = (argc >= 4) ? argv[3] : "ALL"; 
        int n_arg = (argc >= 5) ? atoi(argv[4]) : -1;     
        long double a_vals[] = {0.0L, 1.0L, 0.0L};
        long double b_vals[] = {M_PI / 4.0L, 1.5L, 1.0L};
        long double (*funcs[])(long double) = {f1, f2, f3};
        int func_idx = -1;
        if (strcmp(f_name, "f1") == 0) func_idx = 0;
        else if (strcmp(f_name, "f2") == 0) func_idx = 1;
        else if (strcmp(f_name, "f3") == 0) func_idx = 2;
        else {
            fprintf(stderr, "Función desconocida: use f1, f2 o f3.\n");
            return 1;
        }

        long double a = a_vals[func_idx];
        long double b = b_vals[func_idx];

        printf("=============================================================\n");
        printf(">>> Integral: %s en [%.3Lf, %.3Lf]\n", f_name, a, b);
        printf("-------------------------------------------------------------\n");
        printf("%-10s %-10s %-25s\n", "Método", "n", "Aproximación");
        printf("-------------------------------------------------------------\n");

        // Newton-Codes Abierto 
        if (strcmp(method, "NCAbierto") == 0 || strcmp(method, "ALL") == 0) {
            if (strcmp(method, "NCAbierto") == 0 && (n_arg < 0 || n_arg > 3)) {
                fprintf(stderr, "NC Abierto: n debe estar en [0..3]\n");
                return 1;
            }
            int start = (strcmp(method, "ALL") == 0) ? 0 : n_arg;
            int end   = (strcmp(method, "ALL") == 0) ? 3 : n_arg;
            for (int n = start; n <= end; n++) {
                Vector X_abierto = initialize_Vector(n + 2);
                long double h_ab = (b - a) / (n + 2);
                for (int i = 0; i < n + 2; i++)
                    X_abierto.v[i] = a + h_ab * (i + 1);
                long double res = Newton_Cotes_Abierto(n, &X_abierto, funcs[func_idx], a, b);
                printf("%-10s %-10d %-25.10Lf\n", "NC Abierto", n, res);
                free_Vector(&X_abierto);
            }
        }

        // Newton-Codes Cerrado 
        if (strcmp(method, "NCCerrado") == 0 || strcmp(method, "ALL") == 0) {
            if (strcmp(method, "NCCerrado") == 0 && (n_arg < 1 || n_arg > 4)) {
                fprintf(stderr, "NC Cerrado: n debe estar en [1..4]\n");
                return 1;
            }
            int start = (strcmp(method, "ALL") == 0) ? 1 : n_arg;
            int end   = (strcmp(method, "ALL") == 0) ? 4 : n_arg;
            for (int n = start; n <= end; n++) {
                Vector X_cerrado = initialize_Vector(n + 1);
                long double h_c = (b - a) / n;
                for (int i = 0; i < n + 1; i++)
                    X_cerrado.v[i] = a + h_c * i;
                long double res = Newton_Cotes_Cerrado(n, &X_cerrado, funcs[func_idx], a, b);
                printf("%-10s %-10d %-25.10Lf\n", "NC Cerrado", n, res);
                free_Vector(&X_cerrado);
            }
        }

        // Cuadratura Gaussiana
        if (strcmp(method, "Gaussiana") == 0 || strcmp(method, "ALL") == 0) {
            if (strcmp(method, "Gaussiana") == 0 && (n_arg < 1 || n_arg > 5)) {
                fprintf(stderr, "Gaussiana: n debe estar en [1..5]\n");
                return 1;
            }
            int start = (strcmp(method, "ALL") == 0) ? 1 : n_arg;
            int end   = (strcmp(method, "ALL") == 0) ? 5 : n_arg;
            for (int n = start; n <= end; n++) {
                long double res = Cuadrática_Gaussiana(n, NULL, funcs[func_idx], a, b);
                printf("%-10s %-10d %-25.10Lf\n", "Gaussiana", n, res);
            }
        }

        printf("=============================================================\n\n");
        return 0;
    }
    else {
        printf("Uso incorrecto.\n");
        printf("Formas válidas:\n");
        printf("  ./Programa (sin argumentos) -> 18 corridas automáticas\n");
        printf("  ./Programa P2 [Cos|Polinomial|RBF] [lambda] [FILENAME] [OUTNAME]\n");
        printf("  ./Programa P3 [f1|f2|f3] [NCAbierto|NCCerrado|Gaussiana|ALL] [n]\n");
        return 1;
    }
}

