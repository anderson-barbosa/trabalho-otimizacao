#include <cmath>
#include <vector>
#include <iostream>
using namespace std;

#define Y 0.8
#define N 0.3

typedef long double (*Funcao)(long double, long double);
typedef vector<long double> (*Gradiente)(long double, long double);
typedef vector<vector<long double> > (*Hessiana)(long double, long double);

//expressão da primeira função f1
long double f1(long double x1, long double x2) {
	return pow(x1,2)+pow((exp(x1)-x2),2);
}

//expressão da segunda função f2
long double f2(long double x1, long double x2) {
	return log(1+pow(x1,2)+pow((exp(x1)-x2),2));
}

//calcula a primeira derivada da função f1, em relação a x1 e x2 respectivamente
vector<long double> grad_f1(long double x1, long double x2) {
	vector<long double> ret(2);
	ret[0] = 2*x1+2*exp(x1)*(exp(x1)-x2);
	ret[1] = -2*(exp(x1)-x2);
	return ret;
}

//calcula a primeira derivada da função f2, em relação a x1 e x2 respectivamente
vector<long double> grad_f2(long double x1, long double x2) {
	vector<long double> ret(2);
	ret[0] = (2*x1+2*exp(x1)*(exp(x1)-x2))/(1+pow(x1,2)+pow((exp(x1)-x2),2));
	ret[1] = (-2*(exp(x1)-x2))/(pow(x1,2)+pow((exp(x1)-x2),2)+1);
	//cout << ret[0]<< " " << ret[1] << " " << x1 << " " << x2 << endl;
	return ret;
}

//Algoritmo da Busca de Armijo
long double armijo(vector<long double> xbarra, vector<long double> d, long double y, long double n, Funcao func, Gradiente grad) {
	long double t = 1;
	while(func(xbarra[0]+t*d[0], xbarra[1]+t*d[1])>func(xbarra[0], xbarra[1])+n*t*(grad(xbarra[0], xbarra[1])[0]*d[0]+grad(xbarra[0], xbarra[1])[1]*d[1])) {
		t = y*t;
	}
	return t;
}

//Algoritmo do Método do Gradiente
vector<long double> gradiente(vector<long double> x0, Funcao func, Gradiente grad, long double tol) {
	int k = 0;
	vector<long double> xk=x0;
	while ((abs(grad(xk[0], xk[1])[0])>tol || abs(grad(xk[0], xk[1])[1])>tol) && k<1000) { //número de iterações limitada em 1000
		vector<long double> dk(2);
		dk[0] = -grad(xk[0], xk[1])[0];
		dk[1] = -grad(xk[0], xk[1])[1];
		long double tk = armijo(xk, dk, Y, N, func, grad);
		xk[0]+=tk*dk[0];
		xk[1]+=tk*dk[1];
		k+=1;
	}
	cout << "Número de iterações: " << k << endl;
	return xk;
}

//calcula a matriz inversa de uma matriz m 2x2
vector<vector<long double> > inversa(vector<vector<long double> > m) {
	vector<vector<long double> > ret(2);
	ret[0].resize(2);
	ret[1].resize(2);
	long double d = m[0][0]*m[1][1]-m[0][1]*m[1][0]; // determinante da matriz
	// cout << d << endl;
	ret[0][0] = m[1][1]/d;
	ret[1][1] = m[0][0]/d;
	ret[0][1] = -m[0][1]/d;
	ret[1][0] = -m[1][0]/d;
	return ret;
}

//Calcula a hessiana da função f1
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

//calcula a hessiana da função f2
vector<vector<long double> > h_f2(long double x1, long double x2) {
	vector<vector<long double> > ret(2);
	ret[0].resize(2);
	ret[1].resize(2);
	ret[0][0] = (-2*exp(x1)*pow(x2,3)-2*exp(x1)*pow(x1,2)*x2-6*exp(x1)*x2+8*exp(x1)*x1*x2-2*exp(3*x1)*x2+4*exp(2*x1)*pow(x2,2)+6*exp(2*x1)+4*exp(2*x1)*pow(x1,2)-8*exp(2*x1)*x1-2*pow(x1,2)+2*pow(x2,2)+2)/pow((1+pow(x1,2)+pow(exp(x1)-x2,2)),2);
	ret[0][1] = (2*exp(x1)*pow(x2,2)-2*exp(x1)*pow(x1,2)-2*exp(x1)+4*exp(x1)*x1-4*exp(2*x1)*x2+2*exp(3*x1)-4*x1*x2)/pow((1+pow(x1,2)+pow(exp(x1)-x2,2)),2);
	ret[1][0] = ret[0][1];
	ret[1][1] = -2*(exp(2*x1)-2*exp(x1)*x2+pow(x2,2)-pow(x1,2)-1)/pow((1+pow(x1,2)+pow(exp(x1)-x2,2)),2);
	return ret;
}

//calcula a multiplicação entre 2 matrizes, uma 2x2, a, e um vetor de 2 linhas e 1 coluna,b.
vector<long double> mult(vector<vector<long double> > a, vector<long double> b) {
	vector<long double> ret(2);
	ret[0] = a[0][0]*b[0]+a[0][1]*b[1];
	ret[1] = a[1][0]*b[0]+a[1][1]*b[1];
	// cout << a[0][0] << endl;
	return ret;
}

//Algoritmo do Método de Newton
vector<long double> newton(vector<long double> x0, Funcao func, Gradiente grad, Hessiana h, long double tol) {
	int k = 0;
	vector<long double> xk = x0;
	while ((abs(grad(xk[0], xk[1])[0])>tol || abs(grad(xk[0], xk[1])[1])>tol) && k<100) {
		vector<long double> dk(2);
		dk[0] = -mult(inversa(h(xk[0], xk[1])), grad(xk[0], xk[1]))[0];
		dk[1] = -mult(inversa(h(xk[0], xk[1])), grad(xk[0], xk[1]))[1];
		long double tk = armijo(xk, dk, Y, N, func, grad);
		//cout << dk[0] << endl;
		xk[0]+=tk*dk[0];
		xk[1]+=tk*dk[1];
		k+=1;
		//cout << tk << endl;
	}
	cout << "Número de iterações: " << k << endl;
	return xk;
}

// calcula a subtração de 2 vetores, ou seja, 2 matrizes de 2 linhas e 1 coluna, a e b.
vector<long double> sub(vector<long double> a, vector<long double> b) {
	vector<long double> ret(2);
	ret[0]=a[0]-b[0];
	ret[1]=a[1]-b[1];
	return ret;

}

//Método BFGS, que atualiza Hk de modo que, ao longo das iterações, a matriz se aproxime da inversa da Hessiana.
vector<vector<long double> > bfgs(vector<long double> pk, vector<long double> qk, vector<vector<long double> > hk) {
	long double a = pk[0]*qk[0]+pk[1]*qk[1];
	//cout << pk[0]<< " "<< pk[1] << " "<< qk[0] << " "<< qk[1] << endl;
	vector<long double> b(2);
	b[0] = hk[0][0]*qk[0]+hk[0][1]*qk[1];
	b[1] = hk[1][0]*qk[0]+hk[1][1]*qk[1];
	long double c = 1 + (qk[0]*b[0]+qk[1]*b[1])/a;
	vector<vector<long double> > d(2);
	d[0].resize(2);
	d[1].resize(2);
	d[0][0] = c*(pk[0]*pk[0])/a;
	d[0][1] = c*(pk[0]*pk[1])/a;
	d[1][0] = c*(pk[1]*pk[0])/a;
	d[1][1] = c*(pk[1]*pk[1])/a;
	
	vector<vector<long double> > e(2);
	e[0].resize(2);
	e[1].resize(2);
	e[0][0] = (pk[0]*qk[0]*hk[0][0]+pk[0]*qk[1]*hk[1][0])/a;
	e[0][1] = (pk[0]*qk[0]*hk[0][1]+pk[0]*qk[1]*hk[1][1])/a;
	e[1][0] = (pk[1]*qk[0]*hk[0][0]+pk[1]*qk[1]*hk[1][0])/a;
	e[1][1] = (pk[1]*qk[0]*hk[0][1]+pk[1]*qk[1]*hk[1][1])/a;
	vector<vector<long double> > f(2);
	f[0].resize(2);
	f[1].resize(2);
	f[0][0] = (qk[0]*pk[0]*hk[0][0]+qk[1]*pk[0]*hk[0][1])/a;
	f[0][1] = (qk[0]*pk[1]*hk[0][0]+qk[1]*pk[1]*hk[0][1])/a;
	f[1][0] = (qk[0]*pk[0]*hk[1][0]+qk[1]*pk[0]*hk[1][1])/a;
	f[1][1] = (qk[0]*pk[1]*hk[1][0]+qk[1]*pk[1]*hk[1][1])/a;
	vector<vector<long double> > g(2);
	g[0].resize(2);
	g[1].resize(2);
	g[0][0] = e[0][0]+f[0][0];
	g[0][1] = e[0][1]+f[0][1];
	g[1][0] = e[1][0]+f[1][0];
	g[1][1] = e[1][1]+f[1][1];
	vector<vector<long double> > ret(2);
	ret[0].resize(2);
	ret[1].resize(2);
	ret[0][0] = hk[0][0]+d[0][0]-g[0][0];
	ret[0][1] = hk[0][1]+d[0][1]-g[0][1];
	ret[1][0] = hk[1][0]+d[1][0]-g[1][0];
	ret[1][1] = hk[1][1]+d[1][1]-g[1][1];
	return ret;
}

//Algoritmo do Método Quase-Newton
vector<long double> quaseNewton(vector<long double> x0, vector<vector<long double> > h0, Funcao func, Gradiente grad, long double tol) {
	int k = 0;
	vector<long double> xk = x0;
	vector<vector<long double> > hk=h0;
	while ((abs(grad(xk[0], xk[1])[0])>tol || abs(grad(xk[0], xk[1])[1])>tol) && k<100) {
		vector<long double> dk(2);
		dk[0] = -mult(hk, grad(xk[0], xk[1]))[0];
		dk[1] = -mult(hk, grad(xk[0], xk[1]))[1];
		long double tk = armijo(xk, dk, Y, N, func, grad);
		vector<long double> xl=xk;
		//cout << xk[0] << " " << xk[1] << endl;
		xk[0]+=tk*dk[0];
		xk[1]+=tk*dk[1];
		//cout << dk[0] << " "<< dk[1] << endl;
		if (abs(xk[0]-xl[0])<1e-8 && abs(xk[0]-xl[0])<1e-8) {break;}
		vector<long double> pk = sub(xk,xl); 
		vector<long double>  qk = sub(grad(xk[0], xk[1]), grad(xl[0], xl[1]));
		hk = bfgs(pk, qk, hk); 
		k+=1;
	}
	cout << "Número de iterações: " << k << endl;
	return xk;
}

int main(){

	//cout << f1(1,1) << endl;
	//cout << f2(1,1) << endl;
	vector<long double> xbarra(2);
	//xbarra[0]=-0.9;
	//xbarra[1]=0;

	vector<vector<long double> > h0(2);
	h0[0].resize(2);
	h0[1].resize(2);
	h0[0][0] = 1;
	h0[0][1] = 0;
	h0[1][0] = 0;
	h0[1][1] = 1;
	vector<long double> x;
	long double pontos[][2] = {{0.5,1.5},{0.1,1.1},{0.9,1.7},{-0.5,-1.5},{0.35,1.05}};
	cout << "Função 1 Newton" << endl;
	for (int i=0; i<5; i++) {
		xbarra[0]=pontos[i][0];
		xbarra[1]=pontos[i][1];
		vector<long double> x = newton(xbarra, f1, grad_f1, h_f1,  0);
		cout << "Ponto inicial: " << "(" << xbarra[0] << ", " << xbarra[1] << "); ";
		cout << "Ponto ótimo: " << "(" << x[0] << ", " << x[1] << "); ";
		cout << "Valor ótimo: " << f1(x[0], x[1]) << "; ";
		cout << "Erro: " << abs(f1(x[0], x[1])-f1(0,1)) << endl;	
	}
	
	//vector<long double> x = gradiente(xbarra, f2, grad_f2);
	//vector<long double> x = newton(xbarra, f1, grad_f1, h_f1);
	//vector<long double> x = quaseNewton(xbarra, h0, f2, grad_f2);
	//cout << x[0] << " " << x[1] << endl;
	//cout << "Erro = " << abs(f2(x[0], x[1])-f2(0,1)) << endl;
}
