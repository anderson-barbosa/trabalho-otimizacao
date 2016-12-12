#include <cmath>
#include <vector>
#include <iostream>
using namespace std;

#define Y 0.8
#define N 0.3


long double f1(long double x1, long double x2) {
	return pow(x1,2)+pow((exp(x1)-x2),2);
}

long double f2(long double x1, long double x2) {
	return log(1+pow(x1,2)+pow((exp(x1)-x2),2));
}

vector<long double> grad_f1(long double x1, long double x2) {
	vector<long double> ret(2);
	ret[0] = 2*x1+2*exp(x1)*(exp(x1)-x2);
	ret[1] = -2*(exp(x1)-x2);
	return ret;
}

vector<long double> grad_f2(long double x1, long double x2) {
	vector<long double> ret(2);
	ret[0] = (2*x1+2*exp(x1)*(exp(x1)-x2))/(1+pow(x1,2)+pow((exp(x1)-x2),2));
	ret[1] = (-2*(exp(x1)-x2))/(pow(x1,2)+pow((exp(x1)-x2),2)+1);
	return ret;
}

typedef long double (*Funcao)(long double, long double);
typedef vector<long double> (*Gradiente)(long double, long double);
typedef vector<vector<long double> > (*Hessiana)(long double, long double);

long double armijo(vector<long double> xbarra, vector<long double> d, long double y, long double n, Funcao func, Gradiente grad) {
	long double t = 1;
	while(func(xbarra[0]+t*d[0], xbarra[1]+t*d[1])>func(xbarra[0], xbarra[1])+n*t*(grad(xbarra[0], xbarra[1])[0]*d[0]+grad(xbarra[0], xbarra[1])[1]*d[1])) {
		t = y*t;
	}
	return t;
}

vector<long double> gradiente(vector<long double> x0, Funcao func, Gradiente grad) {
	int k = 0;
	vector<long double> xk=x0;
	while ((grad(xk[0], xk[1])[0]!=0 || grad(xk[0], xk[1])[1]!=0) && k<1000) {
		vector<long double> dk(2);
		dk[0] = -grad(xk[0], xk[1])[0];
		dk[1] = -grad(xk[0], xk[1])[1];
		long double tk = armijo(xk, dk, Y, N, func, grad);
		xk[0]+=tk*dk[0];
		xk[1]+=tk*dk[1];
		k+=1;
	}
	return xk;
}

vector<vector<long double> > inversa(vector<vector<long double> > m) {
	vector<vector<long double> > ret(2);
	ret[0].resize(2);
	ret[1].resize(2);
	long double d = ret[0][0]*ret[1][1]-ret[0][1]*ret[1][0];
	ret[0][0] = m[1][1]/d;
	ret[1][1] = m[0][0]/d;
	ret[0][1] = -m[0][1]/d;
	ret[1][0] = -m[1][0]/d;
	return ret;
}

vector<vector<long double> > h_f1(long double x1, long double x2) {
	vector<vector<long double> > ret(2);
	ret[0].resize(2);
	ret[1].resize(2);
	ret[0][0] = 2+2*(2*exp(2*x1)-exp(x1)*x2);
	ret[0][1] = -2*exp(x1);
	ret[1][0] = -2*exp(x1);
	ret[1][1] = 2;
	return ret;
}

vector<long double> mult(vector<vector<long double> > a, vector<long double> b) {
	vector<long double> ret(2);
	ret[0] = a[0][0]*b[0]+a[0][1]*b[1];
	ret[1] = a[1][0]*b[0]+a[1][1]*b[1];
	return ret;
}

vector<long double> newton(vector<long double> x0, Funcao func, Gradiente grad, Hessiana h) {
	int k = 0;
	vector<long double> xk = x0;
	while ((grad(xk[0], xk[1])[0]!=0 || grad(xk[0], xk[1])[1]!=0) && k<1000) {
		vector<long double> dk(2);
		dk[0] = -mult(inversa(h(xk[0], xk[1])), grad(xk[0], xk[1]))[0];
		dk[1] = -mult(inversa(h(xk[0], xk[1])), grad(xk[0], xk[1]))[1];
		long double tk = armijo(xk, dk, Y, N, func, grad);
		xk[0]+=tk*dk[0];
		xk[1]+=tk*dk[1];
		k+=1;
	}
	return xk;
}

vector<long double> quaseNewton(vector<long double> x0, vector<vector<long double> > h0, Funcao func, Gradiente grad) {
	int k = 0;
	vector<long double> xk = x0;
	vector<vector<long double> > hk=h0;
	while ((grad(xk[0], xk[1])[0]!=0 || grad(xk[0], xk[1])[1]!=0) && k<1000) {
		vector<long double> dk(2);
		dk[0] = -mult(hk, grad(xk[0], xk[1]))[0];
		dk[1] = -mult(hk, grad(xk[0], xk[1]))[1];
		long double tk = armijo(xk, dk, Y, N, func, grad);
		vector<long double> xl=xk;
		xk[0]+=tk*dk[0];
		xk[1]+=tk*dk[1];
		pk = sub(xk,xl); // a fazer
		qk = sub(grad(xk[0], xk[1]), grad(xl[0], xl[1])); //a fazer
		hk = calc_hk(pk, qk);  //a fazer
		k+=1;
	}
	return xk;
}



int main(){

	cout << f1(1,1) << endl;
	cout << f2(1,1) << endl;
	vector<long double> xbarra(2);
	xbarra[0]=2;
	xbarra[1]=2;
	vector<long double> x = newton(xbarra, f1, grad_f1, h_f1);
	cout << x[0] << " " << x[1] << endl;
	cout << "Erro = " << abs(f2(x[0], x[1])) << endl;
}
