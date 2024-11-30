#include"Header.h"
#include <iostream>
#include <cmath>

int read_matrix(double* A, int n, FILE* input) {
	int i, j, c;
	double elem;
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			c = fscanf(input, "%lf", &elem);
			if (c == EOF || !c) {
				return -1;
			}
			A[i*n + j] = elem;
		}
	}
	return 0;
}

double formula(int s, int n, int i, int j) {
	if (s == 1) {
		return (i+1) > (j+1) ? (n - (i+1) + 1) : (n - (j+1) + 1);
	}
	else if (s == 2) {
        if(i == j){
            return 2;
        }
        int tmp = i - j;
        tmp = tmp > 0 ? tmp : -1* tmp; 
        if(tmp == 1){
            return -1;
        }
        else{
            return 0;
        }
	}
	else if (s == 3) {
        if(i == j && j + 1 < n){
            return 1;
        }
        if(j + 1 == n){
            return i + 1;
        }
        if(i + 1 == n){
            return j + 1;
        }
		return 0;
	}
	else if (s == 4) {
		return (1 / double(i+1 + j+1 - 1));
	}
	return 0;
}

int fill_matrix(double* A, int n, int s) {
	int i, j;
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			A[i * n + j] = formula(s, n, i, j);
		}
	}
	return 0;
}

int Check_Sym(double* a, int n){
    int i = 0, j = 0;
    double eps = 1e-16;
    for(i = 0; i < n; ++i){
        for(j = 0; j < n; ++j){
            if(fabs(a[i*n + j] - a[j*n + i]) > eps){
                return 0;
            }
        }
    }
    return 1;
}
void print_matrix(double* matr, int l, int n, int r) {
	int i = 0, j = 0;
	int m = (l > r ? r : l); //string
	int k = (n > r ? r : n); // column
	for (i = 0; i < m; ++i) {
		for (j = 0; j < k; ++j) {
			printf(" %10.3e", matr[i*l+j]);
		//	printf(" %f", matr[i*(l)+j]);
		}
		printf("\n");
	}
}

void Rotate(double* Matrix, int SizeMatrix, double* Sin, double* Cos, double Norm) // func to rotate the matrix
{
   int i, j, k;
   double sin, cos, d, a, b;
   for (i=0;i<SizeMatrix;i++)
    {
        for (j=i+2;j<SizeMatrix;j++)
        {
            if (fabs(Matrix[j*SizeMatrix+i])>1e-16*Norm)
            {
                d = sqrt(Matrix[(i+1)*SizeMatrix+i]*Matrix[(i+1)*SizeMatrix+i] + Matrix[j*SizeMatrix+i]*Matrix[j*SizeMatrix+i]);
                sin = Matrix[j*SizeMatrix+i] / d;
                cos = Matrix[(i+1)*SizeMatrix+i] / d;
                Sin[j] = sin;
                Cos[j] = cos;
         
                for (k=i;k<SizeMatrix;k++)
                {
                    a = Matrix[(i+1)*SizeMatrix+k];
                    b = Matrix[j*SizeMatrix+k];
                    Matrix[(i+1)*SizeMatrix+k] = cos*a+sin*b;
                    if (k!=i){
                        Matrix[j*SizeMatrix+k] = cos*b-sin*a;
                    }
                    else{
                        Matrix[j*SizeMatrix+k] = 0.;
                    }
                }
            }
            else 
            {
                Sin[j] = 2;
                Cos[j] = 2;
            }
        }
     
        for (j=i+2;j<SizeMatrix;j++)
        {
            if (!(Sin[j] > 1 || Cos[j] > 1))
            {
                sin = Sin[j];
                cos = Cos[j];
         
                for (k=i;k<SizeMatrix;k++)
                {
                    a = Matrix[k*SizeMatrix+i+1];
                    b = Matrix[k*SizeMatrix+j];
                    Matrix[k*SizeMatrix+i+1] = cos*a+sin*b;
                    if (k!=i){
                        Matrix[k*SizeMatrix+j] = cos*b-sin*a;
                    }
                    else{
                        Matrix[k*SizeMatrix+j] = 0.;
                    }
                }
            }
        }
    }
}


int n_lambda(double* a, int n, double lambda){
    double x = 0, y = 0, a_k = 0, b_k1 = 0, tmp = 0, tmp1 = 0, gamma = 0, u = 0, v = 0;
    int res = 0, k = 0;
    x = a[0 * n + 0] - lambda;
    y = 1;
    //printf("x = %f\n", x);
    //printf("lambda = %f\n", lambda);
    if(x  < 0){
        res = 1;
    }
    else{
        res = 0;
    }

    for(k = 1; k < n; ++k){
        a_k = a[k * n + k] - lambda;
        b_k1 = a[k * n + k - 1];
        tmp = b_k1 * b_k1 * y;
        tmp = tmp > 0? tmp : -1*tmp;
        tmp1 = x > 0? x: -1*x;
        if(tmp < tmp1){
            tmp = tmp1;
        }
        if (tmp < 1e-50){
			tmp = 1e-15;
        }
        gamma = 1e15/ tmp;
        u = gamma*(a_k*x - b_k1*b_k1*y);
        v = gamma * x;
        if((u*x) < 0){
            ++res;
        }
        x = u;
        y = v;
    }
    //printf("res = %d\n", res);
    return res;
}



//void Eigenvalues(double* a, int n, double eps, int* count_value, double* value, double left, double right, int* it){
int main_eigenvalues(double* a, double* value, int n, double eps, double left, double right, double Norm){
    int Iterations = 0;
    int k = 0;
    Eigenvalues(a, value, n, eps, left, right, Iterations, k, Norm);
   
    if (k!=n) return -1;
    return Iterations;
}

int Eigenvalues(double* a, double* value, int n, double eps, double left, double right, int &Iterations, int &k, double Norm) // recursive func to find eigen values in [a;b]
{
    int n_a = n_lambda(a, n, left);
    int n_b = n_lambda(a, n, right);
   
    if ((right-left>eps*Norm) && (n_b-n_a!=0)){
        Iterations++;
        Eigenvalues(a, value, n, eps, left, (left+right)/2, Iterations, k, Norm);
        Eigenvalues(a, value, n, eps, (left+right)/2, right, Iterations, k, Norm);
    }
    else{
        if (n_b-n_a!=0){
            //printf("value = %e\n", (a+b)/2); 
            for (int i=k; i<k+n_b-n_a; i++) value[i] = (left+right)/2;
            //printf("n_b-n_a=%d\n", n_b-n_a);
            k += n_b-n_a;
        }
        else{
            return 0;
        }
    }
    return 0;
}


double norma(double* a, int n){
    double res = 0, res1 = 0;
    int i = 0, j = 0;
    for(j = 0; j < n; ++j){
        res += fabs(a[j]);
    }
    //printf("res = %f\n", res);
    for(i = 1; i < n; ++i){
        res1 = 0;
        for(j = 0; j < n; ++j){
            res1 += fabs(a[i*n + j]);
        }
        if(res1 > res){
            res = res1;
        }
    }
    return res;
}

double norma_stol(double* a, int n){
    double res = 0, res1 = 0;
    int i = 0, j = 0;
    for(i = 0; i < n; ++i){
        res += fabs(a[i*n]);
    }
    //printf("res = %f\n", res);
    for(j = 1; j < n; ++j){
        res1 = 0;
        for(i = 0; i < n; ++i){
            res1 += fabs(a[i*n + j]);
        }
        if(res1 > res){
            res = res1;
        }
    }
    return res;
}

double residual1(double* a, int n, double* value){
    double res = 0, tr = 0, lambda = 0;
    double eps = 1e-25;
    int i = 0;
    res = norma_stol(a, n);
    if(res < eps){
        return -1;
    }
    for(i = 0; i < n; ++i){
        tr += a[i*n + i];
    }
    for(i = 0; i < n; ++i){
        lambda += value[i];
    }
    return fabs(tr - lambda)/res;
}

double residual2(double* a, int n, double* value){
    double res = 0, tr = 0, lambda = 0;
    double eps = 1e-25;
    int i = 0, j = 0;
    res = norma_stol(a, n);
    if(res < eps){
        return -1;
    }
    for(i = 0; i < n; ++i){
        for(j = 0; j < n; ++j)
        {
            tr += (a[i*n + j] * a[j*n + i]);
        }
    }
    for(i = 0; i < n; ++i){
        lambda += (value[i]*value[i]);
    }
    //printf("tr = %f, %f\n", tr, sqrt(tr));
    //printf("lambda = %f, %f\n", lambda, sqrt(lambda));
    return fabs(sqrt(tr) - sqrt(lambda))/res;
}
