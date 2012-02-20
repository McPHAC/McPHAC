#include "includes.h"

double GauntFactor(double localnu, double T) {
	// Calculates the Gaunt factor using bilinear interpolation from table

	extern double **gauntff, **gauntfactor, *ugrid, *g2grid;
	extern int NGaunt, NGauntu, NGauntg2;
	double  g2, u, utmp, g2tmp, gff;
	char ustr[50], g2str[50], gffstr[50], filename[500];
	FILE *fp;
	int i, n, m;
	double ua, ub, uc, ud, g2a, g2b, g2c, g2d, gffa, gffb, gffc, gffd, t, s;

	if (gauntff == NULL) { // Gaunt table is not populated, so read from file
		sprintf(filename,"%s", GFF_FILE);
		if ((fp = fopen(filename, "r")) == NULL) {
			printf("Can't open %s to retrieve gaunt factors\n", filename);
			exit(1);
		}

		// Allocate table
		n = 0;
		m = 3;
		while (fscanf(fp, "%s %s %s", g2str, ustr, gffstr) == 3) {
			n++;
		}
		gauntff = ddvector(n);
		for (i=0; i<n; i++) {
			gauntff[i] = dvector(m);
		}
		NGaunt = n;

		// Populate table
		fseek(fp, 0, SEEK_SET);
		i = 0;
		while (fscanf(fp, "%s %s %s", g2str, ustr, gffstr)==3 && i<NGaunt) {
			gauntff[i][0] = atof(ustr);
			gauntff[i][1] = atof(g2str);
			gauntff[i][2] = atof(gffstr);
			i++;
		}
		fclose(fp);

		// Allocate gauntfactor grid
		n = 0;
		m = 0;
		utmp = 0;
		g2tmp = 0;
		for (i=0; i<NGaunt; i++) {
			if (gauntff[i][0]>utmp) {
				utmp = gauntff[i][0];
				n++;
			}
			if (gauntff[i][1]>g2tmp) {
				g2tmp = gauntff[i][1];
				m++;
			}
		}
		gauntfactor = ddvector(n);
		for (i=0; i<n; i++) {
			gauntfactor[i] = dvector(m);
		}
		ugrid = dvector(n);
		g2grid = dvector(m);
		NGauntu = n;
		NGauntg2 = m;

		// Populate grid
		ugrid[0] = gauntff[0][0];
		g2grid[0] = gauntff[0][1];
		n = 0;
		m = 0;
		for (i=0; i<NGaunt; i++) {
			if (gauntff[i][0]>ugrid[n]) {
				n++;
				ugrid[n] = gauntff[i][0];
			} else if (gauntff[i][0]<ugrid[n]) {
				n=0;
			}
			if (gauntff[i][1]>g2grid[m]) {
				m++;
				g2grid[m] = gauntff[i][1];
			}
			gauntfactor[n][m] = gauntff[i][2];
		}
		fprintf(stderr, "Created gauntfactor array with %d points in u and %d points in g2\n", NGauntu, NGauntg2);
	}

	// Find requested Gaunt factor using bilinear interpolation
	u = localnu*HPLANCK/T/KBOLTZMAN;
	g2 = RY*ERGPEREV/T/KBOLTZMAN;
	if (u>ugrid[0] && g2>g2grid[0] && u<ugrid[NGauntu-1] && g2<g2grid[NGauntg2-1]) {
		n = 0;
		m = 0;
		while (ugrid[n]<u) n++;
		while (g2grid[m]<g2) m++;
		ua = ugrid[n-1];
		ub = ugrid[n];
		uc = ugrid[n];
		ud = ugrid[n-1];
		g2a = g2grid[m-1];
		g2b = g2grid[m-1];
		g2c = g2grid[m];
		g2d = g2grid[m];
		gffa = gauntfactor[n-1][m-1];
		gffb = gauntfactor[n][m-1];
		gffc = gauntfactor[n][m];
		gffd = gauntfactor[n-1][m];
		t = (u-ua)/(ub-ua);
		s = (g2-g2b)/(g2c-g2b);
	} else if (g2>g2grid[0] && u<ugrid[NGauntu-1] && g2<g2grid[NGauntg2-1]) {
		m = 0;
		while (g2grid[m]<g2) m++;
		g2a = g2grid[m-1];
		g2b = g2grid[m-1];
		g2c = g2grid[m];
		g2d = g2grid[m];
		gffa = gauntfactor[0][m-1];
		gffb = gauntfactor[0][m-1];
		gffc = gauntfactor[0][m];
		gffd = gauntfactor[0][m];
		t = 0;
		s = (g2-g2b)/(g2c-g2b);
	} else if (u>ugrid[0] && u<ugrid[NGauntu-1] && g2<g2grid[NGauntg2-1]) {
		n = 0;
		while (ugrid[n]<u) n++;
		ua = ugrid[n-1];
		ub = ugrid[n];
		uc = ugrid[n];
		ud = ugrid[n-1];
		gffa = gauntfactor[n-1][0];
		gffb = gauntfactor[n][0];
		gffc = gauntfactor[n][0];
		gffd = gauntfactor[n-1][0];
		t = (u-ua)/(ub-ua);
		s = 0.0;
	} else if (u>ugrid[0] && g2>g2grid[0] && g2<g2grid[NGauntg2-1]) {
		m = 0;
		while (g2grid[m]<g2) m++;
		g2a = g2grid[m-1];
		g2b = g2grid[m-1];
		g2c = g2grid[m];
		g2d = g2grid[m];
		gffa = gauntfactor[NGauntu-1][m-1];
		gffb = gauntfactor[NGauntu-1][m-1];
		gffc = gauntfactor[NGauntu-1][m];
		gffd = gauntfactor[NGauntu-1][m];
		t = 0;
		s = (g2-g2b)/(g2c-g2b);
	} else if (u>ugrid[0] && g2>g2grid[0] && u<ugrid[NGauntu-1]) {
		n = 0;
		while (ugrid[n]<u) n++;
		ua = ugrid[n-1];
		ub = ugrid[n];
		uc = ugrid[n];
		ud = ugrid[n-1];
		gffa = gauntfactor[n-1][NGauntg2-1];
		gffb = gauntfactor[n][NGauntg2-1];
		gffc = gauntfactor[n][NGauntg2-1];
		gffd = gauntfactor[n-1][NGauntg2-1];
		t = (u-ua)/(ub-ua);
		s = 0.0;
	}
	gff = (1.0-t)*(1.0-s)*gffa + t*(1.0-s)*gffb + t*s*gffc + (1.0-t)*s*gffd;
	return(gff);
}

