#pragma once
#include<stdio.h>

int read_matrix(double* A, int n, FILE* input);
int fill_matrix(double* A, int n, int s);
void print_matrix(double* matr, int l, int n, int r);
int Check_Sym(double* a, int n);
void Rotate(double* Matrix, int SizeMatrix, double* Sin, double* Cos, double Norm);
int main_eigenvalues(double* a, double* value, int n, double eps, double left, double right, double Norm);
int Eigenvalues(double* a, double* value, int n, double eps, double left, double right, int &Iterations, int &k, double Norm);
double norma(double* a, int n);
double residual1(double* a, int n, double* value);
double residual2(double* a, int n, double* value);
