//  Created by Kaleb Aguilar on 13/08/25.
//
#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef struct {
    int n;        
    long double *v;
} Vector;

typedef struct{ 
    long double **m;
    int X;
    int Y;
} Matrix;


//MATRIX
Matrix initialize_Matrix(int Y,int X){
    Matrix m;
    m.X = X;
    m.Y = Y;
    m.m = (long double **)calloc(Y, sizeof(long double *));
    for(int i = 0; i < Y; i++){
        m.m[i] = (long double *)calloc(X, sizeof(long double));
    }
    return m;
}

Matrix create_Matrix_from_FILE(FILE *f) {
    char *line = NULL;
    size_t len = 0;

    // LEER DIMENSIONES
    if (getline(&line, &len, f) == -1) {
        printf("Error: archivo vacío\n");
        free(line);
        return initialize_Matrix(0,0);
    }

    int rows = 0, cols = 0;
    char *token = strtok(line, " \t\r\n");
    if (token) rows = atoi(token);
    token = strtok(NULL, " \t\r\n");
    if (token) cols = atoi(token);

    // INICIALIZAR MATRIZ (de long double)
    Matrix m = initialize_Matrix(rows, cols);

    // LLENAR MATRIZ
    for (int i = 0; i < rows; i++) {
        if (getline(&line, &len, f) == -1) {
            printf("Error: filas insuficientes en el archivo\n");
            break;
        }
        int j = 0;
        token = strtok(line, " \t\r\n");
        while (token && j < cols) {
            m.m[i][j] = strtold(token, NULL); // <-- long double
            token = strtok(NULL, " \t\r\n");
            j++;
        }
    }

    free(line);
    return m;
}

Matrix initialize_I(int X,int Y){
    Matrix m;
    m.X = X;
    m.Y = Y;
    m.m = (long double **)malloc(sizeof(long double *)*Y);
    for(int i = 0; i < m.Y; i++){
        m.m[i] = (long double *)malloc(sizeof(long double)*X);
        for(int j = 0; j < m.X; j++){
            if(i==j){
                m.m[i][j] = 1.0;
            }else{
                m.m[i][j] = 0.0;
            }
        }
    }
    return m;
}

void free_Matrix(Matrix *m){
    if(m->m != NULL){
        for(int i = 0; i < m->Y; i++){
            free(m->m[i]);
        }
        free(m->m);
    }
}

void print_Matrix(Matrix m){
    printf("Matrix (%d x %d):\n", m.Y, m.X);
    for(int i = 0; i < m.Y; i++){
        for(int j = 0; j < m.X; j++){
            printf("%10.10Lf ", m.m[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

Matrix multiply_Matrix(Matrix *A, Matrix *B) {
    if (A->X != B->Y) {
        printf("Error: dimensiones incompatibles para la multiplicación (%d x %d) * (%d x %d)\n",
               A->Y, A->X, B->Y, B->X);
        exit(1);
    }
    Matrix C = initialize_Matrix(A->Y, B->X);
    for (int i = 0; i < A->Y; i++) {
        for (int j = 0; j < B->X; j++) {
            C.m[i][j] = 0;
            for (int k = 0; k < A->X; k++) {
                C.m[i][j] += A->m[i][k] * B->m[k][j];
            }
        }
    }
    return C;
}

Matrix add_Matrix(Matrix *A, Matrix *B) {
    if (A->X != B->X || A->Y != B->Y) {
        printf("Error: dimensiones incompatibles para la suma (%d x %d) + (%d x %d)\n",
               A->Y, A->X, B->Y, B->X);
        exit(1);
    }
    Matrix C = initialize_Matrix(A->Y, A->X);
    for (int i = 0; i < A->Y; i++) {
        for (int j = 0; j < A->X; j++) {
            C.m[i][j] = A->m[i][j] + B->m[i][j];
        }
    }
    return C;
}
Matrix diff_Matrix(Matrix *A, Matrix *B) {
    if (A->X != B->X || A->Y != B->Y) {
        printf("Error: dimensiones incompatibles para la resta (%d x %d) - (%d x %d)\n",
               A->Y, A->X, B->Y, B->X);
        exit(1);
    }
    Matrix C = initialize_Matrix(A->Y, A->X);
    for (int i = 0; i < A->Y; i++) {
        for (int j = 0; j < A->X; j++) {
            C.m[i][j] = A->m[i][j] - B->m[i][j];
        }
    }
    return C;
}



//VECTOR
Vector initialize_Vector(int n) {
    Vector vec;
    vec.n = n;
    vec.v = (long double *)calloc(n, sizeof(long double));
    return vec;
}

void print_Vector(Vector a) {
    printf("[ ");
    for (int i = 0; i < a.n; i++) {
        printf("%6.10Lf ", a.v[i]);
    }
    printf("]\n");
}

Vector create_Vector_from_FILE(FILE *f) {
    char *line = NULL;
    size_t len = 0;

    // LEER DIMENSIÓN
    if (getline(&line, &len, f) == -1) {
        printf("Error: archivo vacío\n");
        free(line);
        return initialize_Vector(0);
    }

    int n = atoi(line);  // dimensión del vector

    // INICIALIZAR VECTOR
    Vector v = initialize_Vector(n);

    // LEER ELEMENTOS
    for (int i = 0; i < n; i++) {
        if (getline(&line, &len, f) == -1) {
            printf("Error: elementos insuficientes en el archivo\n");
            break;
        }
        v.v[i] = strtold(line, NULL); // convierte cada línea a long double
    }

    free(line);
    return v;
}



double Vector_norm(Vector *v) {
    double sum = 0.0;
    for (int i = 0; i < v->n; i++) {
        sum += v->v[i] * v->v[i];
    }
    return sqrt(sum);
}

Vector diff_Vector(Vector *a, Vector *b){
    if(a->n != b->n){
        printf("Error: dimensiones distintas en diff_Vector\n");
        Vector empty = {0, NULL};
        return empty;
    }

    Vector r = initialize_Vector(a->n);
    for(int i = 0; i < a->n; i++){
        r.v[i] = a->v[i] - b->v[i];
    }
    return r;
}
void diff_Vector_inplace(Vector *a, Vector *b,Vector *out){
    if(a->n != b->n){
        printf("Error: dimensiones distintas en diff_Vector\n");
        return;
    }
    for(int i = 0; i < a->n; i++){
        out->v[i] = a->v[i] - b->v[i];
    }
}
void save_Vector_to_csv(Vector *v, const char *filename) {
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        printf("Error al abrir el archivo %s\n", filename);
        return;
    }
    fprintf(f, "Index,Value\n");
    for (int i = 0; i < v->n; i++) {
        fprintf(f, "%d,%Lf\n", i, v->v[i]);
    }
    fclose(f);
    printf("Vector guardado en %s\n", filename);
}

void copy_vector(Vector *from, Vector *to){
	to->n = from->n;
	if (to->v != NULL) {
            free(to->v);
        }
    to->v = (long double *)malloc(from->n * sizeof(long double));
    for (int i = 0; i < from->n; i++) {
        to->v[i] = from->v[i];
    }
}

void free_Vector(Vector *a) {
    if (a->v != NULL) {
        free(a->v);
        a->v = NULL;
    }
    a->n = 0;
}




//OPERACIONES ENTRE VECTORES Y MATRIX
void scalar_times_Vector(double alpha, Vector *v, Vector *out) {
    if (v->n != out->n) {
        printf("Error: dimensiones incompatibles en scalar_times_Vector\n");
        return;
    }

    for (int i = 0; i < v->n; i++) {
        out->v[i] = alpha * v->v[i];
    }
}

void sum_Vectors(Vector *a, Vector *b, Vector *out) {
    if (a->n != b->n || a->n != out->n) {
        printf("Error: dimensiones incompatibles en sum_Vector\n");
        return;
    }

    for (int i = 0; i < a->n; i++) {
        out->v[i] = a->v[i] + b->v[i];
    }
}

void sum_Matrix(Matrix *A, Matrix *B, Matrix *out) {
    if (A->X != B->X || A->Y != B->Y || A->X != out->X || A->Y != out->Y) {
        printf("Error: dimensiones incompatibles en sum_Matrix\n");
        return;
    }

    for (int i = 0; i < A->Y; i++) {
        for (int j = 0; j < A->X; j++) {
            out->m[i][j] = A->m[i][j] + B->m[i][j];
        }
    }
}
void diff_Matrix_inplace(Matrix *A, Matrix *B, Matrix *out) {
    if (A->X != B->X || A->Y != B->Y || A->X != out->X || A->Y != out->Y) {
        printf("Error: dimensiones incompatibles en sum_Matrix\n");
        return;
    }

    for (int i = 0; i < A->Y; i++) {
        for (int j = 0; j < A->X; j++) {
            out->m[i][j] = A->m[i][j] - B->m[i][j];
        }
    }
}
double Vector_times_Vector(Vector *u, Vector *v) {
    if (u->n != v->n) {
        printf("Error: dimensiones incompatibles para producto de vectores (%d vs %d)\n", u->n, v->n);
        return 0.0;
    }

    double result = 0.0;
    for (int i = 0; i < u->n; i++) {
        result += u->v[i] * v->v[i];
    }
    return result;
}

void Matrix_times_Vector(Matrix *A, Vector *x, Vector *result){
    if (A->X != x->n) {
        printf("Error: dimensiones incompatibles para multiplicación.\n");
        return;
    }
    for (int i = 0; i < A->Y; i++) {
        result->v[i] = 0.0;
        for (int j = 0; j < A->X; j++) {
            result->v[i] += A->m[i][j] * x->v[j];
        }
    }
}

void Matrix_times_Matrix(Matrix *A, Matrix *B, Matrix *result) {
    if (A->X != B->Y) {
        printf("Error: dimensiones incompatibles para multiplicación (%d x %d) * (%d x %d)\n", A->Y, A->X, B->Y, B->X);
        return;
    }
    for (int i = 0; i < result->Y; i++) {
        for (int j = 0; j < result->X; j++) {
            result->m[i][j] = 0.0;
            for (int k = 0; k < A->X; k++) {
                result->m[i][j] += A->m[i][k] * B->m[k][j];
            }
        }
    }
}

void transpose(Matrix *M, Matrix *T) {
    for (int i = 0; i < M->Y; i++)
        for (int j = 0; j < M->X; j++)
            T->m[j][i] = M->m[i][j];
}

double Matrix_Norm_Frobenius(Matrix *Diff) {
    double sum = 0.0;
    for (int i = 0; i < Diff->Y; i++) {
        for (int j = 0; j < Diff->X; j++) {
            sum += Diff->m[i][j] * Diff->m[i][j];
        }
    }
    return sqrt(sum);
}

// COPY VECTOR Y MATRIX
void copy_Vector(Vector *from, Vector *to) {
    if (from->n != to->n) {
        printf("Error: dimensiones incompatibles en copy_Vector (%d vs %d)\n", from->n, to->n);
        return;
    }

    for (int i = 0; i < from->n; i++) {
        to->v[i] = from->v[i];
    }
}

void copy_Matrix(Matrix *from, Matrix *to) {
    if (from->X != to->X || from->Y != to->Y) {
        printf("Error: dimensiones incompatibles en copy_Matrix (%dx%d vs %dx%d)\n",
               from->X, from->Y, to->X, to->Y);
        return;
    }

    for (int i = 0; i < from->Y; i++) {
        for (int j = 0; j < from->X; j++) {
            to->m[i][j] = from->m[i][j];
        }
    }
}

void Vector_row(Matrix *A, int row , Vector *out){
	if (row >= A->Y) return;
	for(int i = 0; i < A->X ;i++){
		out->v[i] = A->m[row][i];
	}

}
void Vector_col(Matrix *A, int col, Vector *out){
	if (col >= A->X) return;
	for(int i = 0; i < A->Y ;i++){
		out->v[i] = A->m[i][col];
	}
}
#endif
