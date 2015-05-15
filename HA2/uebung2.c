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

const double SI_C = 299792458; /* m/s */
const double h = 1E-11; /* infinitissimale Länge der Ableitung */
const double prec = 1E-6; /* Abbruchbedingung der Reihendarstellung */


/* Realteil der Besselfunktion 0-ter Ordnung 
 * Berechnung durch eine Potenzreihenentwicklung */
double ber_ps (double x, double prec) {
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
/*		printf("%i ber_neu = %.10f\n",k,ber_neu); */
	}while (fabs(ber_neu/ber_alt-1) > prec);

	return ber_neu;
}


/* Imaginärteil der Besselfunktion 0-ter Ordnung
 * Berechnung durch eine Potenzreihenentwicklung */
double bei_ps (double x, double prec) {
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
/*		printf("%i ber_neu = %.10f\n",k,ber_neu); */
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
double ber(double x) {
	if (x<10) {
		return ber_ps(x,prec);
	} else {
		return creal(J0_an(x*csqrt(-I)));
	}
}


/* Definiere den Imaginärteil der Besselfunktion 0-ter Ordnung.
 * Nutze dazu eine Fallunterscheidung. */
double bei(double x) {
	if (x<10) {
		return bei_ps(x,prec);
	} else {
		return cimag(J0_an(x*csqrt(-I)));
	}
}


/* Ableitung des Realteils der Besselfunktion 0-ter Ordnung
 * Berechnung durch 3-Punkt-Methode
 * unter Verwendung der Potenzreihendarstellung */
double derive_ber(double x) {
	return (ber(x+h)-ber(x))/h;
}


/* Ableitung des Imaginärteils der Besselfunktion 0-ter Ordnung
 * Berechnung durch 3-Punkt-Methode
 * unter Verwendung der Potenzreihendarstellung */
double derive_bei(double x) {
	return (bei(x+h)-bei(x))/h;
}


/* * * * * * * * * * Physik  * * * * * * * * */


/* Berechne den konstanten Vorfaktor der Stromdichte. 
double prefactor(double I0, double kappa, double rho0) {
	return I0*kappa/2*M_1_PI/rho0;
} */


/* Berechne die Stromdichte eines Leiters ohne den konstanten Vorfaktor. */
double j_1(double kappa, double rho, double rho0) {
	return (ber(kappa*rho)+I*bei(kappa*rho))/(derive_ber(kappa*rho0)-I*derive_bei(kappa*rho0));
}


void show_skin_effect() {
	const double omega = 2*M_PI*5*1E10;  /* 1/s - Kreisfrequenz bei Frequenz f=50MHz. */
	const double sigma = 58.0*1E6;       /* 1/ohm/m - elektrische Leitfähigkeit von Kupfer */
	const double mu = 1-6.4*1E-6;        /* relative Ladungsträgerbeweglichkeit*/
	const double rho0 = 1*1E-2;          /* m - Leiterradius */
	const double I0 = 1;                 /* A - maximale Stromstärke */
	const double kappa = 2*sqrt(M_PI*mu*sigma*omega)/SI_C;
	const double prefactor = I0*kappa/2*M_1_PI/rho0;

	const double prec = 1E-8;
	const double xmin = 0;
	const double xmax = rho0;
	const int n = 1000;

	double xi[n+1];
	double fxi[n+1];

	int i;
	for(i=0;i<=n; i++){
		/* Berechne die n-te Stützstelle */
    		xi[i]=xmin+i*(xmax-xmin)/((double) n);
		/* Berechne den Funktionswert an der Stützstelle */
		fxi[i]=prefactor*j_1(kappa,xi[i],rho0);
    		
		printf("%f %f \n",xi[i],fxi[i]);
	}
}


/* * * * * * * * Test-Funktionen * * * * * * */

double ber_small_test(double x) { return ber(x); }
double bei_small_test(double x) { return bei(x); }
double derive_ber_small_test(double x) { return derive_ber(x); }
double derive_ber_big_test(double x) { return derive_ber(x)/exp(x*M_SQRT1_2); }
double derive_bei_small_test(double x) { return derive_bei(x); }
double derive_bei_big_test(double x) { return derive_bei(x)/exp(x*M_SQRT1_2); }
double ber_big_test(double x) { return ber(x)/exp(x*M_SQRT1_2); }
double bei_big_test(double x) { return bei(x)/exp(x*M_SQRT1_2); }


void plotIt() {
//	plot(&ber_small_test,0,10,1000); /* Plottet ber(x) */
//	plot(&bei_small_test,0,10,1000); /* Plottet ber(x) */
	plot(&derive_ber_small_test,0,10,1000); /* Plottet ber'(x) */
	plot(&derive_ber_big_test,0,100,1000); /* Plottet ber'(x) */
	plot(&derive_bei_small_test,0,10,1000); /* Plottet ber'(x) */
	plot(&derive_bei_big_test,0,100,1000); /* Plottet ber'(x) */
//	plot(&ber_big_test,0,100,1000);/* Plottet J0(x)/exp(x/sqrt(2)) zwischen 0 und 100 */
//	plot(&bei_big_test,0,100,1000);/* Plottet J0(x)/exp(x/sqrt(2)) zwischen 0 und 100 */
}




int main(){
/*	plotIt(); */
/*	show_skin_effect(); */
/*	printf("Error: Not Implemented.\n"); */
/*	test_continuation(); */

	int i;
	for (i = 5;i <15;i++) {
		printf("%.2f  %.2f %.2f\n",(double) i,(double) derive_ber(i),(double) derive_bei(i));
	}

	return 0;
}
