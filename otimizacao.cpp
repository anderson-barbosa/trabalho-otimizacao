#include <cmath>
#include <vector>
#include <iostream>
using namespace std;

#define Y 0.8
#define N 0.3

typedef long double Tipo;

Tipo f1(Tipo x1, Tipo x2) {
	return pow(x1,2)+pow((exp(x1)-x2),2);
}

Tipo f2(Tipo x1, Tipo x2) {
	return log(1+pow(x1,2)+pow((exp(x1)-x2),2));
}

vector<Tipo> grad_f1(Tipo x1, Tipo x2) {
	vector<Tipo> ret(2);
	ret[0] = 2*x1+2*exp(x1)*(exp(x1)-x2);
	ret[1] = 2*(exp(x1)-x2);
	return ret;
}

vector<Tipo> grad_f2(Tipo x1, Tipo x2) {
	vector<Tipo> ret(2);
	ret[0] = (2*x1+2*exp(x1)*(exp(x1)-x2))/(1+pow(x1,2)+pow((exp(x1)-x2),2));
	ret[1] = (-2*(exp(x1)-x2))/(pow(x1,2)+pow((exp(x1)-x2),2)+1);
	return ret;
}

typedef Tipo (*Funcao)(Tipo, Tipo);
typedef vector<Tipo> (*Gradiente)(Tipo, Tipo);

Tipo armijo(vector<Tipo> xbarra, vector<Tipo> d, Tipo y, Tipo n, Funcao func, Gradiente grad) {
	Tipo t = 1;
	while(func(xbarra[0]+t*d[0], xbarra[1]+t*d[1])>func(xbarra[0], xbarra[1])+n*t*(grad(xbarra[0], xbarra[1])[0]*d[0]+grad(xbarra[0], xbarra[1])[1]*d[1])) {
		t = y*t;
	}
	return t;
}

vector<Tipo> gradiente(vector<Tipo> x0, Funcao func, Gradiente grad) {
	Tipo k = 0;
	vector<Tipo> xk=x0;
	while ((grad(xk[0], xk[1])[0]!=0 || grad(xk[0], xk[1])[1]!=0) && k<10000) {
		vector<Tipo> dk(2);
		dk[0] = -grad(xk[0], xk[1])[0];
		dk[1] = -grad(xk[0], xk[1])[1];
		Tipo tk = armijo(xk, dk, Y, N, func, grad);
		xk[0]+=tk*dk[0];
		xk[1]+=tk*dk[1];
		k+=1;
	}
	return xk;
}

int main(){

	cout << f1(1,1) << endl;
	cout << f2(1,1) << endl;
	vector<Tipo> xbarra(2);
	xbarra[0]=2;
	xbarra[1]=2;
	vector<Tipo> x = gradiente(xbarra, f1, grad_f1);
	cout << x[0] << " " << x[1] << endl;
	cout << "Erro = " << abs(f1(x[0], x[1])) << endl;
}