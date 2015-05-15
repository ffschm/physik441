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
/*		printf("%i ber_neu = %.10f\n",k,ber_neu); */
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
double ber0(double x,double prec) {
	if (x<10) {
		return ber0_ps(x,prec);
	} else {
		return creal(J0_an(x*csqrt(-I)));
	}
}


/* Definiere den Imaginärteil der Besselfunktion 0-ter Ordnung.
 * Nutze dazu eine Fallunterscheidung. */
double bei0(double x,double prec) {
	if (x<10) {
		return bei0_ps(x,prec);
	} else {
		return cimag(J0_an(x*csqrt(-I)));
	}
}


/* Ableitung des Realteils der Besselfunktion 0-ter Ordnung
 * Berechnung durch 3-Punkt-Methode
 * unter Verwendung der Potenzreihendarstellung */
double derive_ber0(double x, double prec) {
	const double h = prec;
	return (ber0(x+h,prec)-ber0(x-h,prec))/h/2;
}


/* Ableitung des Imaginärteils der Besselfunktion 0-ter Ordnung
 * Berechnung durch 3-Punkt-Methode
 * unter Verwendung der Potenzreihendarstellung */
double derive_bei0(double x, double prec) {
	const double h = prec;
	return (bei0(x+h,prec)-bei0(x-h,prec))/h/2;
}


/* * * * * * * * * * Physik  * * * * * * * * */


/* Berechne den konstanten Vorfaktor der Stromdichte. 
double prefactor(double I0, double kappa, double rho0) {
	return I0*kappa/2*M_1_PI/rho0;
} */


/* Berechne die Stromdichte eines Leiters ohne den konstanten Vorfaktor. */
double j_1(double kappa, double rho, double rho0, double prec) {
	return (ber0(kappa*rho,prec)+I*bei0(kappa*rho,prec))/(derive_ber0(kappa*rho0,prec)-I*derive_bei0(kappa*rho0,prec));
}


/* Berechne die Stromdichte eines Leiters.
 *
 * I0 - Stromamplitude
 * rho0 - Radius des Leiters
 * sigma - Elektrische Leitfähigkeit
 * mu - Ladungsträgerbeweglichkeit
 * omega - Kreisfrequenz des Stromes
 * rho - Abstand von der Drahtachse
 * prec - Genauigkeit
 */
double j(double I0, double rho0, double sigma, double mu, double omega, double rho, double prec) {
	const double kappa = 2*sqrt(M_PI*mu*sigma*omega)/SI_C;
	const double prefactor = I0*kappa/2*M_1_PI/rho0;

	return prefactor*j_1(kappa, rho, rho0, prec);
}


/* Berechne für die Kreisfrequenz omega die Verteilung der Stromdichte
 * über den Leiterquerschnitt. */
void calculate_current_density_test() {
	printf("Error: Not implemented.\n");
}


void show_skin_effect() {
	const double omega = 2*M_PI*5*1E10;    /* 1/s - Kreisfrequenz bei Frequenz f=50MHz. */
	const double mu = 58.0*1E6;    /* 1/ohm/m - elektrische Leitfähigkeit von Kupfer */
	const double sigma = 5*1E-3;   /* m*m/V/s - Ladungsträgerbeweglichkeit*/
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
		fxi[i]=prefactor*j_1(kappa,xi[i],rho0,prec);
    		
		printf("%f %f \n",xi[i],fxi[i]);
	}
}


/* * * * * * * * Test-Funktionen * * * * * * */

double ber0_small_test(double x) { return ber0(x,1E-6); }
double bei0_small_test(double x) { return bei0(x,1E-6); }
double derive_ber0_small_test(double x) { return derive_ber0(x,1E-5); }
double derive_ber0_big_test(double x) { return derive_ber0(x,1E-5)/exp(x*M_SQRT1_2); }
double ber0_big_test(double x) { return ber0(x,1E-6)/exp(x*M_SQRT1_2); }
double bei0_big_test(double x) { return bei0(x,1E-6)/exp(x*M_SQRT1_2); }


void plotIt() {
	plot(&ber0_small_test,0,10,1000); /* Plottet ber(x) */
	plot(&bei0_small_test,0,10,1000); /* Plottet ber(x) */
	plot(&derive_ber0_small_test,0,10,1000); /* Plottet ber'(x) */
	plot(&derive_ber0_big_test,0,100,1000); /* Plottet ber'(x) */
	plot(&ber0_big_test,0,100,1000);/* Plottet J0(x)/exp(x/sqrt(2)) zwischen 0 und 100 */
	plot(&bei0_big_test,0,100,1000);/* Plottet J0(x)/exp(x/sqrt(2)) zwischen 0 und 100 */
}


void test_continuation () {
	int i;
	const double prec = 1E-6;

	printf("ps      as\n");
	for (i=19;i<22;i++) {
		printf("%.2f  %.5f  %.5f\n",(double) 0.5*i, (double) ber0_ps(0.5*i,prec),(double) J0_an(0.5*i*csqrt(-I)));
	}
}


int main(){
/*	plotIt(); */
	show_skin_effect();
/*	printf("Error: Not Implemented.\n"); */
/*	test_continuation(); */

	return 0;
}
