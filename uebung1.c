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
double madelung2d (int n) {
	/* definiere Laufvariablen */
	int i;
	int j;
	double sum = 0.0;
	/* summiere ueber wuerfel */
	for (i=-n; i<=n; i++) {
		for (j=-n; j<=n; j++) {
			/* Sonderbehandlung Nullpunkt*/
			if (i==0 && j==0) {			/* mittleres Ion */
			} else if (abs(i)==n && abs(j)==n) {	/* auf einer Kante */
				sum += powm1(i+j)/sqrt(i*i+j*j)/4;
			} else if (abs(i)==n || abs(j)==n) {	/* auf einer Seite */
				sum += powm1(i+j)/sqrt(i*i+j*j)/2;
			} else {				/* im Wuerfel */
				sum += powm1(i+j)/sqrt(i*i+j*j);
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
				} else if (abs(i)==n && abs(j)==n && abs(k)==n) {
					sum += powm1(i+j+k)/sqrt(i*i+j*j+k*k)/8;

				} else if (abs(i)==n || abs(j)==n || abs(k)==n) {
					sum += powm1(i+j+k)/sqrt(i*i+j*j+k*k)/2;

				} else if ((abs(i)==n && abs(j)==n) || (abs(i)==n && 
					   abs(k)==n) || (abs(k)==n && abs(j)==n)) {
					sum += powm1(i+j+k)/sqrt(i*i+j*j+k*k)/4;

				} else {
					sum += powm1(i+j+k)/sqrt(i*i+j*j+k*k);
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
	int n=1;
	double m_alt;
	double m_neu;

	/* Berechne Vorfaktor k */
	double k = 1/(4*M_PI*eps_0)*e*e;

	/* Fallunterscheidung fuer 3dim, 2dim */
	if (dim == 2) {
		/* siehe unten */
		do {
			m_alt = m_neu;
			m_neu = madelung2d(n);
			printf("E(n = %i) = %.20f\n",n,m_neu);
			n++;
		} while (fabs(m_neu/m_alt-1) > var);
	} else {
		/* Iteriere über Kantenlänge n des Würfels,
		 * bis die Abweihung zwischen E(n) und E(n-1) kleiner als
		 * die gewünschte Genauigkeit ist.
		 */
		do {
			m_alt = m_neu;
			m_neu = madelung3d(n);
			printf("E(n = %i) = %.20f\n",n,m_neu);
			n++;
		} while (fabs(m_neu/m_alt-1) > var);
	}
	return m_neu*k/d;
}


int main () {
	/* Genauigkeit */
	double var = 1E-5;
	long double Gitterkonstante = 5.62E-10;
	int dim = 2;
	double Energie = bindungsenergie(dim, Gitterkonstante, var);
	/* TODO: Einheit!!*/
	printf("Bindungsenergie in %iD mit einer Genauigkeit von %E beträgt:\n %E Einheit??\n",dim,var,Energie);
	return 0;	
}
