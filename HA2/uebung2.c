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
double ber0 (double x, double prec) {
	int k;
	double t,ber_neu,ber_alt;
	t = 1;
	ber_neu = 0;
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

/* Besselelfunktion 0-ter Ordnung
 * Berechnung durch asymptotische Näherung */
double complex J0_an (double complex z) {
	return M_SQRT2*csqrt(M_1_PI/z)*ccos(z-M_PI_4);
}


void test_ber() {
	printf("Berechne ber0(0.1)...\n");
	printf("bar0_plot (0.1) = %.2f \n" , ber0(0.1,1E-6));

	int i;
	for (i = 0;i == 0; i++){
		printf("ber(%.3f) = %.5f\n",i,ber0(i,1E-10));
	}
}


double ber0_test(double x) {
	return ber0(x,1E-6);
}


double J0_test(double x) {
	return J0_an(x);
}


void plotIt() {
	/* Plottet ber(x) zwischen 0 und 10 */
//	plot(&ber0_test,0,10,1000);

	/* Plottet J0(x) zwischen 5 und 20 */
	plot(&J0_test,5,20,1000);
}


int main(){
	plotIt();
	printf("Error: Not Implemented.\n");

	return 0;
}
