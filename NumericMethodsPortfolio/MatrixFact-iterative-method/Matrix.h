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
    double *v;
} Vector;



typedef struct{ 
    double **m;
    int X;
    int Y;
} Matrix;

Matrix initialize_Matrix(int Y,int X){
    Matrix m;
    m.X = X;
    m.Y = Y;
    m.m = (double **)calloc(Y, sizeof(double *));
    for(int i = 0; i < Y; i++){
        m.m[i] = (double *)calloc(X, sizeof(double));
    }
    return m;
}

Matrix create_Matrix_from_FILE(FILE *f , char c){
    char buffer[1024];
    int rows = 0, cols = 0;

    // CALCULAR DIMENSIONES
    while (fgets(buffer, sizeof(buffer), f)) {
        rows++;

        if (cols == 0) {
            // eliminar salto de línea
            buffer[strcspn(buffer, "\n")] = 0;

            // contar columnas con strtok
            cols = 0;
            char delim[2] = {c, '\0'};
            char *token = strtok(buffer, delim);
            while (token) {
                cols++;
                token = strtok(NULL, delim);
            }
        }
    }

    // volver al inicio del archivo
    rewind(f);

    // INICIALIZAR MATRIX
    Matrix m = initialize_Matrix(rows, cols);

    // LLENAR MATRIX
    int i = 0;
    while (fgets(buffer, sizeof(buffer), f)) {
        buffer[strcspn(buffer, "\n")] = 0;
        char delim[2] = {c, '\0'};
        int j = 0;
        char *token = strtok(buffer, delim);
        while (token && j < cols) {
            m.m[i][j] = strtod(token, NULL);
            token = strtok(NULL, delim);
            j++;
        }
        i++;
    }

    return m;
}


Matrix initialize_I(int X){
    Matrix m;
    m.X = X;
    m.Y = X;
    m.m = (double **)malloc(sizeof(double *)*X);
    for(int i = 0; i < X; i++){
        m.m[i] = (double *)malloc(sizeof(double)*X);
        for(int j = 0; j < X; j++){
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
            printf("%10.1f ", m.m[i][j]); // más espacio y 6 decimales
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
    vec.v = (double *)malloc(sizeof(double) * n);
    return vec;
}

// Imprimir vector
void print_Vector(Vector a) {
    printf("[ ");
    for (int i = 0; i < a.n; i++) {
        printf("%6.6lf ", a.v[i]);
    }
    printf("]\n");
}

Vector create_Vector_from_FILE(FILE *f, char c) {
    char buffer[1024];
    int count = 0;

    while (fgets(buffer, sizeof(buffer), f)) {
        char delim[2] = {c, '\0'};
        char *token = strtok(buffer, delim);
        while (token) {
            count++;
            token = strtok(NULL, delim);
        }
    }

    rewind(f);

    Vector vec;
    vec.n = count;
    vec.v = (double *)malloc(sizeof(double) * count);

    int idx = 0;
    while (fgets(buffer, sizeof(buffer), f)) {
        char delim[2] = {c, '\0'};
        char *token = strtok(buffer, delim);
        while (token && idx < count) {
            vec.v[idx++] = strtod(token, NULL);
            token = strtok(NULL, delim);
        }
    }
    return vec;
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

void save_Vector_to_csv(Vector *v, const char *filename) {
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        printf("Error al abrir el archivo %s\n", filename);
        return;
    }
    fprintf(f, "Index,Value\n");
    for (int i = 0; i < v->n; i++) {
        fprintf(f, "%d,%lf\n", i, v->v[i]);
    }
    fclose(f);
    printf("Vector guardado en %s\n", filename);
}

void free_Vector(Vector *a) {
    free(a->v);
}


/// SOLUCIÓN DE SISTEMAS 

Vector solve_D(Matrix *m, Vector *b) {
    Vector x = initialize_Vector(b->n);
    for (int i = 0; i < b->n; i++) {
        x.v[i] =  b->v[i] / m->m[i][i];
    }
    return x;
}

Vector solve_U(Matrix *m, Vector *b) {
    Vector x = initialize_Vector(b->n);
    for (int i = b->n - 1; i >= 0; i--) {
        x.v[i] = b->v[i];
        for (int j = b->n - 1; j > i; j--) {
            x.v[i] -= m->m[i][j] * x.v[j];
        }
        x.v[i] /= m->m[i][i];
    }
    return x;
}

Vector solve_L(Matrix *m, Vector *b) {
    Vector x = initialize_Vector(b->n);
    for (int i = 0; i < b->n; i++) {
        x.v[i] = b->v[i];
        for (int j = 0; j < i; j++) {
            x.v[i] -= m->m[i][j] * x.v[j];
        }
        x.v[i] /= m->m[i][i];
    }
    return x;
}


int reduce_col_from_row(Matrix *m, Vector *b, int row_1, int row_2) {
    // ELIMINATE COL FROM ROW_2 USING ROW_1.
    int col = row_2;
    if (row_1 >= m->Y || row_2 >= m->Y) return -1;
    if (m->m[row_1][col] == 0) return -1;

    double factor = m->m[row_2][col] / m->m[row_1][col];
    for (int i = 0; i < m->X; i++) {
        m->m[row_2][i] -= factor * m->m[row_1][i];
    }
    b->v[row_2] -= factor * b->v[row_1];
    return 0;
}

void swap_rows(Matrix *m, Vector *b, int row_1, int row_2) {
    double aux;
    for (int i = 0; i < m->X; i++) {
        aux = m->m[row_1][i];
        m->m[row_1][i] = m->m[row_2][i];
        m->m[row_2][i] = aux;
    }
    aux = b->v[row_1];
    b->v[row_1] = b->v[row_2];
    b->v[row_2] = aux;
}
void pivote(Matrix *A, Vector *b, int col){
    if(A->m[col][col] != 0) return;

    for(int i = col + 1; i < A->Y; i++){
        if(A->m[i][col] != 0){
            // SWAP DE FILAS i <--> col
            for(int j = 0; j < A->X; j++){
                double tmp = A->m[col][j];
                A->m[col][j] = A->m[i][j];
                A->m[i][j] = tmp;
            }
            // SWAP EN VECTOR B
            if(b){
                double tmpb = b->v[col];
                b->v[col] = b->v[i];
                b->v[i] = tmpb;
            }
            return;
        }
    }
}

Vector solve_by_Gauss(Matrix *m, Vector *b) {
    for (int i = 0; i < m->Y; i++) {
        for (int j = i + 1; j < m->Y; j++) {
            if (i == j) continue;
            reduce_col_from_row(m, b, i, j);
        }
    }
    Vector x = solve_U(m, b);
    return x;
}

void print_solution(Matrix m, Vector b, Vector x) {
    for (int i = 0; i < m.Y; i++) {
        for (int j = 0; j < m.X; j++) {
            printf("%10.6lf ", m.m[i][j]);
        }
        printf("   |x%2d: %10.6f |b%2d: %10.6lf", i, x.v[i], i, b.v[i]);
        printf("\n");
    }
    printf("\n");
}


Vector Matrix_times_Vector(Matrix *A, Vector *x) {
    if (A->X != x->n) {
        printf("Error: dimensiones incompatibles para multiplicación.\n");
        Vector empty = {0, NULL};
        return empty;
    }
    Vector result;
    result.n = A->Y;
    result.v = (double *)malloc(sizeof(double) * A->Y);

    for (int i = 0; i < A->Y; i++) {
        result.v[i] = 0.0;
        for (int j = 0; j < A->X; j++) {
            result.v[i] += A->m[i][j] * x->v[j];
        }
    }
    return result;
}


#endif
