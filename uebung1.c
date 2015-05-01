/* gcc -o uebung1 uebung1.c -lm
 * ./uebung1
 * Marvin Schmitz, Fabian Schmidt
 *
 * Dieses Programm bestimmt die Energie eines Ions
 * in einem NaCl-Gitter mithilfe der Evjen-Methode.
 */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>

const double e = 1.602176565E-19;     /* Elementarladung */
const double eps_0 = 8.854187817E-12; /* Permittivität */


/*
 * Berechnet (-1)^i.
 */
double powm1(int i){
	int ungerade = i & 1;
	return (1^(-ungerade))+ungerade;
}


/* definiere Funktion, die mir die Madlung Konstante errechnet im 2dim Fall */
double madelung_2dim (int n) {
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
					sum += powm1(i+j)/(4*sqrt(i*i+j*j));
				}
				/*auf der Seite*/
				else if (abs(i)==n || abs(j)==n) {
					sum += powm1(i+j)/(2*sqrt(i*i+j*j));
				}
				/*im Wuerfel*/
				else {
					sum += powm1(i+j)/(sqrt(i*i+j*j));
				}
			}
		}
	}
	return fabs(sum);
}


/*
 * Berechnet die Madelung-Konstante eines 3D-Quaders mit Kantenlänge n.
 */
double madelung3d (int n) {
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
				} else {
					if (abs(i)==n && abs(j)==n && abs(k)==n) {
						sum += powm1(i+j+k)/sqrt(i*i+j*j+k*k)/8;
					}
					else if (abs(i)==n || abs(j)==n || abs(k)==n) {
						sum += powm1(i+j+k)/sqrt(i*i+j*j+k*k)/2;
					}
					else if ((abs(i)==n && abs(j)==n) || (abs(i)==n && 
							 abs(k)==n) || (abs(k)==n && abs(j)==n)) {
						sum += powm1(i+j+k)/sqrt(i*i+j*j+k*k)/4;
					}
					else {
						sum += powm1(i+j+k)/sqrt(i*i+j*j+k*k);
					}
				}
			}
		}
	}
	return fabs(sum);
}


/* Berechnet die Bindungsenergie eines Ions in einem NaCl-Kristall
 * dim - 2: 2D-Kristall, 3: 3D-Kristall
 * d   - Gitterkonstante
 * var - gewünschte Genauigkeit
 */
double bindungsenergie (int dim, double d, double var) {
	int n = 1;
	double E_alt;
	double E_neu;

	/* definiere Konstante k */
	double k = 1/(4*M_PI*eps_0)*pow(e,2);
	/* Fallunterscheidung fuer 3dim, 2dim */
	if (dim == 2) {
		E_alt = madelung_2dim(1)*k/d;
		E_neu = madelung_2dim(2)*k/d;
		/* benutze var um Iteration ueber n zu bestimmen (Groesse des Wuerfels) */
		while (fabs((E_neu)/(E_alt)-1)>var) {
			n++;
			E_alt = E_neu;
			E_neu = madelung_2dim(n)*k/d;
		}
	} else {
		/* benutze var um Iteration ueber n zu bestimmen (Groesse des Wuerfels) */
		do {
			E_alt = E_neu;
			E_neu = madelung3d(n)*k/d;
			printf("E(n = %i) = %.20f\n",n,E_neu);
			n++;
		} while (fabs(E_neu/E_alt-1) > var);
	}
	return E_neu;
}


int main () {
	/* Genauigkeit */
	double var = pow(10,-5);
	long double Gitterkonstante = 5.62*pow(10, -10);
	int dim = 3;
	double Energie = bindungsenergie(dim, Gitterkonstante, var);
	/* TODO: Einheit!!*/
	printf("Bindungsenergie in %iD mit einer Genauigkeit von %E beträgt:\n %E Einheit??\n",dim,var,Energie);
	return 0;	
}
