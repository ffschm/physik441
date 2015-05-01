#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*kompilieren mit: gcc -Wall -o test test.c -lm */

/* Ziel ist es die Energie eines Ions im Na-Cl Gitter zu bestimmen diese
 * approximieren wir mit der Evjen Methode wie in der Problembeschreibung
 * vermerkt */

/* definiere Funktion, die mir die Madlung Konstante errechnet im 2dim Fall */
double madlung_2dim (int n) {
	/* definiere Laufvariablen */
	int i;
	int j;
	double sum = 0.0;
	/* summiere ueber wuerfel */
	for (i=-n; i<=n; i++) {
		for (j=-n; j<=n; j++) {
			/* Sonderbehandlung Nullpunkt*/
			if (i==0 && j==0) {
				sum = sum;
			}
			else {	
				/*auf der Kante*/		
				if (abs(i)==n && abs(j)==n) {
					sum += pow(-1, i+j)/(4*sqrt(i*i+j*j));
				}
				/*auf der Seite*/
				else if (abs(i)==n || abs(j)==n) {
					sum += pow(-1, i+j)/(2*sqrt(i*i+j*j));
				}
				/*im Wuerfel*/
				else {
					sum += pow(-1, i+j)/(sqrt(i*i+j*j));
				}
			}
		}
	}
	return fabs(sum);
}

/*--------------------------------------------------------------------*/

/*analog zu madlung2_dim definiere ich fkt. für dreidimensionalen Wuerfel*/
double madlung_3dim (int n) {
	/* definiere Laufvariablen */
	int i;
	int j;
	int k;
	double sum = 0.0;
	/* summiere ueber wuerfel */
	for (i=-n; i<=n; i++) {
		for (j=-n; j<=n; j++) {
			for (k=-n; k<=n; k++) {
				if (i==0 && j==0 && k==0) {
					sum = sum;
				}
				else {			
					if (abs(i)==n && abs(j)==n && abs(k)==n) {
						sum += pow(-1, i+j+k)/(8*sqrt(i*i+j*j+k*k));
					}
					else if (abs(i)==n || abs(j)==n || abs(k)==n) {
						sum += pow(-1, i+j+k)/(2*sqrt(i*i+j*j+k*k));
					}
					else if ((abs(i)==n && abs(j)==n) || (abs(i)==n && 
							 abs(k)==n) || (abs(k)==n && abs(j)==n)) {
						sum += pow(-1, i+j+k)/(4*sqrt(i*i+j*j+k*k));
					}
					else {
						sum += pow(-1, i+j+k)/(sqrt(i*i+j*j+k*k));
					}
				}
			}
		}
	}
	return fabs(sum);
}

/*--------------------------------------------------------------------*/
			
/* Definiere Funktion, die Energie berechnet*/
double bindungsenergie (int dim, double d, int var) {
	/* d Gitterkonstante (Abstand) */
	/* var ist die Genauigkeit */
	int n = 1;
	double E_0;
	double E_neu;
	double e = 1.602176565*pow(10,-19);
	double eps_0 = 8.854187817*pow(10,-12);
	/* definiere Konstante k */
	double k = 1/(4*M_PI*eps_0)*pow(e,2);
	/* Fallunterscheidung fuer 3dim, 2dim */
	if (dim == 2) {
		E_0 = madlung_2dim(1)*k/d;
		E_neu = madlung_2dim(2)*k/d;
		/* benutze var um Iteration ueber n zu bestimmen (Groesse des Wuerfels) */
		while (fabs((E_neu)/(E_0)-1)>var) {
			n++;
			E_0 = E_neu;
			E_neu = madlung_2dim(n)*k/d;					
		}
	}
	else {
		E_0 = madlung_3dim(1)*k/d;
		E_neu = madlung_3dim(2)*k/d;
		/* benutze var um Iteration ueber n zu bestimmen (Groesse des Wuerfels) */
		while (fabs((E_neu)/(E_0)-1)>var) {
			n++;
			E_0 = E_neu;
			E_neu = madlung_2dim(n)*k/d;
		}	
	}
	return E_neu;
}

/*--------------------------------------------------------------------*/

int main () {
	/* Genauigkeit */
	double var = pow(10,-5);
	long double Gitterkonstante = 5.62*pow(10, -10);
	/*bindungsenergie hat 3 parameter: 1. dimension, 2. Gitterkonstante, 3. Genauigkeit*/
	double Energie = bindungsenergie(3, Gitterkonstante, var);  
	printf("%E\n", Energie);
	return 0;	
}
