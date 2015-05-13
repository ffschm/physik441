/* gcc -o uebung2 uebung2.c plot.c -lm
 * ./uebung2
 * Marvin Schmitz, Fabian Schmidt
 *
 * Dieses Programm implementiert die Besselfunktion 0.Ordnung
 * und deren Ableitung in Real- und Imaginärteil.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>

#include "plot.h"

/* Realteil der Besselfunktion 0-ter Ordnung 
 * Berechnung durch eine Potenzreihenentwicklung */
double ber0_ps (double x, double prec) {
	int k;
	double t,ber_neu,ber_alt;
	t = 1;
	ber_neu = t;
	k = 0;
	do {
		t = -pow(x/2,4)/pow(4*k*k+6*k+2,2)*t;
		ber_alt = ber_neu;
		ber_neu += t;
		k++;
//		printf("%i ber_neu = %.10f\n",k,ber_neu);
	}while (fabs(ber_neu/ber_alt-1) > prec);

	return ber_neu;
}


/* Imaginärteil der Besselfunktion 0-ter Ordnung
 * Berechnung durch eine Potenzreihenentwicklung */
double bei0_ps (double x, double prec) {
	int k;
	double t,bei_neu,bei_alt;
	t = pow(x/2,2);
	bei_neu = t;
	k = 0;
	do {
		t = -pow(x/2,4)/pow(4*k*k+10*k+6,2)*t;
		bei_alt = bei_neu;
		bei_neu += t;
		k++;
//		printf("%i ber_neu = %.10f\n",k,ber_neu);
	}while (fabs(bei_neu/bei_alt-1) > prec);

	return bei_neu;
}


/* Besselelfunktion 0-ter Ordnung
 * Berechnung durch asymptotische Näherung
 * J0(z) = sqrt(2/(pi*z))*cos(z-pi/4) */
double complex J0_an (double complex z) {
	return M_SQRT2*csqrt(M_1_PI/z)*ccos(z-M_PI_4);
}


/* Definiere den Realteil der Besselfunktion 0-ter Ordnung.
 * Nutze dazu eine Fallunterscheidung. */
double ber0(double x,double prec) {
	if (x<30) {
		return ber0_ps(x,prec);
	} else {
		return creal(J0_an(x*csqrt(-I)));
	}
}


/* Definiere den Imaginärteil der Besselfunktion 0-ter Ordnung.
 * Nutze dazu eine Fallunterscheidung. */
double bei0(double x,double prec) {
	if (x<30) {
		return bei0_ps(x,prec);
	} else {
		return cimag(J0_an(x*csqrt(-I)));
	}
}


/* Ableitung des Realteils der Besselfunktion 0-ter Ordnung
 * Berechnung durch h-Methode
 * unter Verwendung der Potenzreihendarstellung */
double derive_ber0(double x, double prec, double h) {
	return (ber0(x+h,prec)-ber0(x,prec))/h;
}


/* Ableitung des Imaginärteils der Besselfunktion 0-ter Ordnung
 * Berechnung durch h-Methode
 * unter Verwendung der Potenzreihendarstellung */
double derive_bei0(double x, double prec, double h) {
	return (bei0(x+h,prec)-bei0(x,prec))/h;
}


/* * * * * * * * Test-Funktionen * * * * * * */

double ber0_small_test(double x) { return ber0(x,1E-6); }
double bei0_small_test(double x) { return bei0(x,1E-6); }
double derive_ber0_small_test(double x) { return derive_ber0(x,1E-5,1E-8); }
double derive_ber0_big_test(double x) { return derive_ber0(x,1E-5,1E-8)/exp(x*M_SQRT1_2); }
double ber0_big_test(double x) { return ber0(x,1E-6)/exp(x*M_SQRT1_2); }
double bei0_big_test(double x) { return bei0(x,1E-6)/exp(x*M_SQRT1_2); }


void plotIt() {
	plot(&ber0_small_test,0,10,1000); /* Plottet ber(x) */
//	plot(&bei0_small_test,0,10,1000); /* Plottet ber(x) */
//	plot(&derive_ber0_small_test,0,10,1000); /* Plottet ber'(x) */
//	plot(&derive_ber0_big_test,0,100,1000); /* Plottet ber'(x) */
	plot(&ber0_big_test,0,100,1000);/* Plottet J0(x)/exp(x/sqrt(2)) zwischen 0 und 100 */
//	plot(&bei0_big_test,0,100,1000);/* Plottet J0(x)/exp(x/sqrt(2)) zwischen 0 und 100 */
}


int main(){
	plotIt();
	printf("Error: Not Implemented.\n");

	return 0;
}
