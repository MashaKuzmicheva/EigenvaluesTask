#include <stdio.h>
#include <time.h>
#include "Header.h"
int main(int argc, char* argv[]) {
	int n = 0, r = 0, s = 0;
	double r1 = -1, r2 = -1, eps = 0;
	char* name = 0;
	double* A = nullptr;
	int its = 0, count_value = 0, i = 0;
	double right = 0, left = 0;
	double* value = nullptr;
	int flag = 0;
	FILE* input;
	double t1 = 0, t2 = 0;
	if (!((argc == 5 || argc == 6)
		&& sscanf(argv[1], "%d", &n) == 1
		&& sscanf(argv[2], "%d", &r) == 1
		&& sscanf(argv[3], "%lf", &eps) == 1
		&& sscanf(argv[4], "%d", &s) == 1)) {
		printf("Usage %s n m r s [file]\n", argv[0]);
		return 1;
	}

	//создадим матрицы
	A = new double[n * n];
	if (!A) {
		delete[] A;
		printf("Matrix A creation error\n");
		return -1;
	}

	if (s == 0 && argc == 6) {
		name = argv[5];

		//откроем файл для чтения
		input = fopen(name, "r");
		if (!input) {
			delete[]A;
			printf("File opening error\n");
			return -1;
		}

		//прочитаем матрицу А
		if (read_matrix(A, n, input) == -1)
		{
			printf("Matrix reading error\n");
			fclose(input);
			delete[]A;
			return -1;
		}
		fclose(input);
	}

	if (s > 0 && s <= 4) {
		if (fill_matrix(A, n, s) == -1) {
			printf("Matrix filling error\n");
			delete[]A;
			return -1;
		}
	}

	flag = Check_Sym(A,n);
	if (!flag)
	{
	   printf("Matrix is asymmetrical\n");
	   printf("%s : Residual1 = %e Residual2 = %e Iterations = %d Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n", argv[0], r1, r2, its, its / n, t1, t2);
	   delete[]A;
	   return -2;
	}
	

	printf("Matrix A: \n");
	print_matrix(A, n, n, r);
	printf("\n");

	value = new double[n];
	if (!value) {
		delete[] A;
		delete[] value;
		printf("Array value creation error\n");
		return -1;
	}

	double norm = norma(A, n);


	double* sin = new double[n];
	if (sin == nullptr)
	{
		delete[] value;
		delete[] A;
		printf("memory was not allocated\n");
		return -1;
	}
	double* cos = new double[n];
	if (cos == nullptr)
	{
		delete[] value;
		delete[] sin;
		delete[] A;
		printf("memory was not allocated\n");
		return -1;
	}

	//Приводим к трехдиагональному виду
	t1 = clock();
	Rotate(A, n, sin, cos, norm);
	//tridiagonal_matrix(n, A);
	t1 = (clock() - t1)/1e+6;

	//print_matrix(A, n, n, r);

	//printf("left = %f right = %f\n", left, right);

	//print_matrix(A, n, n, r);
	right = norm * (1 + eps);
	left = -1 * right;


	//Ищем собственные значения

	t2 = clock();
	its = main_eigenvalues(A, value, n, eps, left, right, norm);
	t2 = (clock() - t2)/1e+6;
	printf("Eigenvalues:\n");
	count_value = n;
	if(count_value > r){
		count_value = r;
	}

	if (count_value == 0) {
		for (i = 0; i < r; ++i) {
			printf("0 ");
		}
		printf("\n");
	}
	else {
		for (i = 0; i < count_value; ++i) { //r
			printf("%10.3e ", value[i]);
			//printf("%f ", value[i]);
		}
		printf("\n");
	}
	printf("\n");


	if (s == 0 && argc == 6) {
		name = argv[5];

		//откроем файл для чтения
		input = fopen(name, "r");
		if (!input) {
			delete[]A;
			delete[] value;
			delete[] A;
			delete[] sin;
			delete[] cos;
			printf("File opening error\n");
			return -1;
		}

		//прочитаем матрицу А
		if (read_matrix(A, n, input) == -1)
		{
			printf("Matrix reading error\n");
			fclose(input);
			delete[] value;
			delete[] A;
			delete[] sin;
			delete[] cos;
			return -1;
		}
		fclose(input);
	}

	if (s > 0 && s <= 4) {
		if (fill_matrix(A, n, s) == -1) {
			printf("Matrix filling error\n");
			delete[] value;
			delete[] A;
			delete[] sin;
			delete[] cos;
			return -1;
		}
	}

	r1 = residual1(A, n, value);
	r2 = residual2(A, n, value);

	printf("%s : Residual1 = %e Residual2 = %e Iterations = %d Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n", argv[0], r1, r2, its, its / n, t1, t2);

	delete[] value;
	delete[] A;
	delete[] sin;
	delete[] cos;
	return 0;
}
