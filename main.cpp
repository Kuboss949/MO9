#include <iostream>
#include <functional>
#include <fstream>
#include <cmath>
#include "thomas.h"
using namespace std;

enum schemeType{
    CONVENTIONAL,
    NUMEROV
};

double Ufunc(double x){
    return (exp(2.0-2.0*x)-4.0*exp(4.0-2.0*x)+4.0*exp(2.0*x)-exp(2.0+2.0*x)-x+x*exp(4.0))/(4.0-4.0*exp(4.0));
}

void conventional(double h, double a, int n, double * l, double * d, double * u, double *result,
                  const function<double(double)>& p,
                  const function<double(double)>& q,
                  const function<double(double)>& r,
                  const function<double(double)>& s){
    double xi = a + h;
    for(int i=1; i<n-1; i++){
        l[i-1] = p(xi)/(h*h) - q(xi)/(2.0*h);
        d[i] = r(xi) - 2.0*p(xi)/(h*h);
        u[i] = p(xi)/(h*h) + q(xi)/(2.0*h);
        result[i] = -s(xi);
        xi+=h;
    }
}

void numerov(double h, double a, int n, double *l, double *d, double *u, double *result,
             const std::function<double(double)>& p,
             const std::function<double(double)>& q,
             const std::function<double(double)>& r,
             const std::function<double(double)>& s) {
    double xi = a + h;
    for (int i = 1; i < n - 1; i++) {
        l[i - 1] = p(xi) / (h*h) + r(xi-h)/12.;
        d[i] = (-2.*p(xi)) / (h*h) + r(xi) * 10./12.;
        u[i] = p(xi) / (h*h) + r(xi+h)/12.;
        result[i] = -(s(xi - h) + 10.*s(xi) + s(xi + h))/12.;

        xi += h;
    }
}

void diffScheme_error(double h, double a, double b, double alpha, double beta, double gamma, double phi, double psi, double theta,
                  const function<double(double)>& p,
                  const function<double(double)>& q,
                  const function<double(double)>& r,
                  const function<double(double)>& s, schemeType scheme, ofstream &plik){

    int n = static_cast<int>((b - a) / h)+1;
    auto *u = new double[n-1];
    auto *d = new double[n];
    auto *l = new double[n-1];
    auto *result = new double[n];
    u[0] = alpha/h;
    d[0] = beta - alpha/h;
    l[n-2] = -phi/h;
    d[n-1] = phi/h + psi;
    result[0] = -gamma;
    result[n-1] = -theta;
    if(scheme == CONVENTIONAL)
        conventional(h, a, n, l, d, u, result, p, q, r, s);
    else if (scheme == NUMEROV)
        numerov(h, a, n, l, d, u, result, p, q, r, s);
    matrixProcedure(l, d, u, n);
    vectorProcedure(result, l, d, u, n);

    double xi = a;
    double maxBlad=0;
    for(int i=0; i<n; i++){
        double blad = fabs(Ufunc(xi)-result[i]);
        if(blad>maxBlad)
            maxBlad = blad;

        xi+=h;
    }
    plik << log10(h) <<" "<< log10(maxBlad)<< '\n';
    delete[] u;
    delete[] d;
    delete[] l;
    delete[] result;
}

void diffScheme(double h, double a, double b, double alpha, double beta, double gamma, double phi, double psi, double theta,
                const function<double(double)>& p,
                const function<double(double)>& q,
                const function<double(double)>& r,
                const function<double(double)>& s, schemeType scheme, ofstream &plik){

    int n = static_cast<int>((b - a) / h)+1;
    auto *u = new double[n-1];
    auto *d = new double[n];
    auto *l = new double[n-1];
    auto *result = new double[n];
    u[0] = alpha/h;
    d[0] = beta - alpha/h;
    l[n-2] = -phi/h;
    d[n-1] = phi/h + psi;
    result[0] = -gamma;
    result[n-1] = -theta;
    if(scheme == CONVENTIONAL)
        conventional(h, a, n, l, d, u, result, p, q, r, s);
    else if (scheme == NUMEROV)
        numerov(h, a, n, l, d, u, result, p, q, r, s);
    matrixProcedure(l, d, u, n);
    vectorProcedure(result, l, d, u, n);
    double xi = a;
    for(int i=0; i<n; i++){
        plik << xi << " " << result[i] << '\n';
        xi+=h;
    }
    delete[] u;
    delete[] d;
    delete[] l;
    delete[] result;
}




int main() {
    double a = 0.0, b = 1.0, alpha = 0.0, beta = 1.0, gamma = -1.0, phi = 0.0, psi = 1.0, theta = 0.0;
    ofstream plik;

    plik.open("Konw.txt");
    diffScheme(0.001, a, b, alpha, beta, gamma, phi, psi, theta,
                     [](double x)->double{return 1.0;},
                     [](double x)->double{return 0.0;},
                     [](double x)->double{return -4.0;},
                     [](double x)->double{return -x;},
                     CONVENTIONAL, plik);
    plik.close();
    plik.open("Numerow.txt");
    diffScheme(0.001, a, b, alpha, beta, gamma, phi, psi, theta,
                     [](double x)->double{return 1.0;},
                     [](double x)->double{return 0.0;},
                     [](double x)->double{return -4.0;},
                     [](double x)->double{return -x;},
                     NUMEROV, plik);
    plik.close();

    plik.open("Analityczne.txt");
    double xi = a;
    for(int i=0; i<10001; i++) {
        plik << xi << " " << Ufunc(xi) << '\n';
        xi+=0.0001;
    }
    plik.close();
    plik.open("wynikKonw.txt");
    for (double h = 1; h>10e-9;h/=1.5)
        diffScheme_error(h, a, b, alpha, beta, gamma, phi, psi, theta,
                         [](double x)->double{return 1.0;},
                         [](double x)->double{return 0.0;},
                         [](double x)->double{return -4.0;},
                         [](double x)->double{return -x;},
                         CONVENTIONAL, plik);
    plik.close();

    plik.open("wynikNumerov.txt");
    for (double h = 1; h>10e-9;h/=1.5)
        diffScheme_error(h, a, b, alpha, beta, gamma, phi, psi, theta,
                         [](double x)->double{return 1.0;},
                         [](double x)->double{return 0.0;},
                         [](double x)->double{return -4.0;},
                         [](double x)->double{return -x;},
                         NUMEROV, plik);
    plik.close();
    return 0;
}
