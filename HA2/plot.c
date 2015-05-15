#include <stdlib.h>
#include <stdio.h>

/*
 * Plottet die übergebene Funktion.
 *
 * funcPtr - die betrachtete Funktion
 * xmin    - untere Grenze
 * xmax    - obere Grenze
 * n       - Anzahl der Stützstellen
 */
void plot(double (*funcPtr)(double),double xmin, double xmax, int n) {
	/* char *gnuplotCommands[] = {"set title \"Title01\"","plot 'data.dump'"};*/

	int i;
	double xi[n+1];
	double fxi[n+1];
	FILE *gnuplotPipe;

	for(i=0;i<=n; i++){
		/* Berechne die n-te Stützstelle */
    		xi[i]=xmin+i*(xmax-xmin)/((double) n);
		/* Berechne den Funktionswert an der Stützstelle */
		fxi[i]=(*funcPtr)(xi[i]);
	}

	
	gnuplotPipe = popen("gnuplot -persistent","w");

	fprintf(gnuplotPipe, "plot '-' \n");

	for (i=0; i<=n; i++) {
    		fprintf(gnuplotPipe, "%f %f \n",xi[i],fxi[i]);
	}

	fprintf(gnuplotPipe, "e");
}
