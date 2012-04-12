#include <R.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double sum(double *x, int dim) {
	double s;
	int i;
	s = 0.0;
	for (i=0; i<dim; ++i) s+=x[i];
	return(s);
}

void LifeTableC(int sex, double *mxm,
				double *LLm) {
	double am[21], lm[28], qmx[19];
	int i;
	
	if(sex > 1) {/* female*/
		if (mxm[0] < 0.107) {
			am[0] = 0.053 + 2.8 * mxm[0];
			am[1] = 1.522 - 1.518 * mxm[0];
		} else {
			am[0] = 0.35;
			am[1] = 1.361;
		}
	} else { /* male */
		if (mxm[0] < 0.107) {
			am[0] = 0.045 + 2.684 * mxm[0];
			am[1] = 1.651 - 2.816 * mxm[0];
		} else {
			am[0] = 0.33;
			am[1] = 1.352;
		}
	}
	/*Rprintf("\nam0=%f, am1=%f", am[0], am[1]);*/
	lm[0] = 1;
	lm[1] = 1 - mxm[0] / (1 + (1 - am[0]) * mxm[0]);
	lm[2] = (1 - 4 * mxm[1] / (1 + (4 - am[1]) * mxm[1]));
	lm[2] = lm[1] * lm[2];
	LLm[0] = (lm[1] + am[0] * (1 - lm[1])) + (4 * lm[2] + am[1] * (lm[1] - lm[2])); /* Life table male 4L0*/
    
	/* Age 5-9, .... 95-99 
	 Greville formula used in Mortpak and UN MLT (1982)*/
	for(i = 2; i < 21; ++i) {
		am[i] = 2.5 - (25 / 12.0) * (mxm[i] - 0.1 * log(fmax(mxm[i+1] / fmax(mxm[i-1], DBL_MIN), DBL_MIN)));
		/*Rprintf("am%i=%f, mxm%i=%f", i, am[i], i-1, mxm[i-1]);*/
		qmx[i-2] = 5 * mxm[i] / (1 + (5 - am[i]) * mxm[i]);
	}
    
    for(i = 1; i<20; ++i) {
		lm[i+2] = lm[i+1] * (1-qmx[i-1]);
		LLm[i] = 5 * lm[i+2] + am[i+1] * (lm[i+1] - lm[i+2]);
	}
	
	for(i = 20; i<26; ++i) {
		lm[i+2] = lm[i+1] * (1 - (1 - exp(-5 * mxm[i+1]))); /* Starts 104 */
		LLm[i] = (lm[i+1] - lm[i+2]) / fmax(mxm[i+1], DBL_MIN);
	}
	
	/* Age 130+ */
	LLm[26] = lm[27] / fmax(mxm[27], DBL_MIN); /* Assuming Mx levels off at age 130 */
	
}

void LCEoKtC(int sex, double *ax, double *bx, 
			 double eop, double kl, double ku, double *LLm) {
	double LTl[27], LTu[27], mxm[28], LTeo, k2;
	int i, dim;
	dim = 27;
	
	/* check if the eop lies outside of the bounds */
	for (i=0; i < 28; ++i) {
		mxm[i] = exp(ax[i] + bx[i]*kl);
	}
	LifeTableC(sex, mxm, LTl);
	
	if(eop < sum(LTl, dim)) {
		LLm = LTl;
		return;
	}
	for (i=0; i < 28; ++i) mxm[i] = exp(ax[i] + bx[i]*ku);
	LifeTableC(sex, mxm, LTu);

	if(eop > sum(LTu, dim)) {
		LLm = LTu;
		return;
	}
	/* Bi-section method */
	k2 = 0.5 * (kl + ku);
	for (i=0; i < 28; ++i) mxm[i] = exp(ax[i] + bx[i]*k2);
	LifeTableC(sex, mxm, LLm);
	/*Rprintf("\nLLm[2]=%lf, k2=%f, kl=%f, ku=%f, mxm0-2=%f %f %f", LLm[2], k2, kl, ku, mxm[0], mxm[1], mxm[2]);*/
	LTeo = sum(LLm, dim);
	while(abs(LTeo - eop) > 0.01) {
		if(LTeo < eop) kl = k2;
		else ku = k2;
		k2 = 0.5 * (kl + ku);
		for (i=0; i < 28; ++i) mxm[i] = exp(ax[i] + bx[i]*k2);
		LifeTableC(sex, mxm, LLm);
		/*Rprintf("\nLTeo=%lf, LLm[2]=%lf, k2=%f, kl=%f, ku=%f, mxm0-2=%f %f %f", LTeo, LLm[2], k2, kl, ku, mxm[0], mxm[1], mxm[2]);*/
		LTeo = sum(LLm, dim);
	}
	/*Rprintf("\nk2=%d, eop=%lf, LTeo=%lf", k2, eop, LTeo);*/
	
}

void get_sx(double *LLm, double *sx) {
	int i;
	/* Survival Ratios */
    sx[0] = LLm[0] / 5.0;
	for(i=1; i < 26; ++i) {
		if(LLm[i-1] == 0) sx[i] = 0;
		else sx[i] = LLm[i]/LLm[i-1];
	}
	/* Last age group */
	if(LLm[25] == 0 ||  (LLm[25]+LLm[26]) == 0) sx[26] = 0;
	else sx[26] = LLm[26] / (LLm[25] + LLm[26]);
	if(sx[26] > sx[25]) sx[26] = sx[25];
}

void LC(int *Npred, int *Sex, double *ax, double *bx, 
		double *Eop, double *Kl, double *Ku, double *LLm, double *Sr) {
	double eop, kl, ku, sx[27];
	int i, sex, npred, pred;
	
	npred = *Npred;
	sex=*Sex;
	ku=*Ku;
	kl=*Kl;
	
	for (pred=0; pred < npred; ++pred) {
		eop = Eop[pred];
		/*Rprintf("\n%i: eop=%lf", pred, eop);*/
		LCEoKtC(sex, ax, bx, eop, kl, ku, LLm);
		get_sx(LLm, sx);

		for (i=0; i < 27; ++i) {
			Sr[i + pred*27] = sx[i];
			/*Rprintf("\nLLm=%lf, Sr=%lf", LLm[i], Sr[i + pred*27]);*/
		}
	}
}


void TotalPopProj(int *npred, double *MIGm, double *MIGf, int *migr, int *migc,
				  int *MIGtype, double *srm, double *srf, double *asfr, double *srb, 
				  double *popm, double *popf, double *totp 
					) {
	double migm[*migr+6][*migc], migf[*migr+6][*migc], totmigm[*migr+6][*migc], totmigf[*migr+6][*migc];
	double b, bt[7], bm, bf, mmult;
	int i,j, adim, nrow, ncol, n;
	nrow = *migr;
	ncol = *migc;
	n = *npred;
	adim=27;
	for(j=0; j<ncol; ++j) {
		for(i=0; i<nrow; ++i) {		
			migm[i][j] = MIGm[i + j*nrow];
			migf[i][j] = MIGf[i + j*nrow];
		}
		for(i=nrow; i<nrow+6; ++i) {
			migm[i][j] = 0;
			migf[i][j] = 0;
		}
		switch (*MIGtype) {
			case 0: /* migration evenly distributed over each interval (MigCode=0) */
				for(i=1; i<nrow+6; ++i) {
					totmigm[i][j] = 0.5*(migm[i][j] + migm[i-1][j]*srm[i + j*adim]);
					totmigf[i][j] = 0.5*(migf[i][j] + migf[i-1][j]*srf[i + j*adim]);
				}
				mmult = 0.5;
				break;
			default: /* migration at the end of each interval (MigCode=9)*/
				for(i=0; i<nrow+6; ++i) {
					totmigm[i][j] = migm[i][j];
					totmigf[i][j] = migf[i][j];
				}
				mmult = 1;
				break;
		}
	}
	/* Population projection for one trajectory */
	for(j=1; j<(n+1); ++j) {
		/* Compute ages >=5 */
		for(i=1; i<(adim-1); ++i) {	
			popm[i + j*adim] = popm[i-1 + (j-1)*adim] * srm[i + (j-1)*adim] + totmigm[i][j-1];
			popf[i + j*adim] = popf[i-1 + (j-1)*adim] * srf[i + (j-1)*adim] + totmigf[i][j-1];
		}
		
		/* Age 100+ */
		popm[26 + j*adim] = (popm[26 + (j-1)*adim] + popm[25 + (j-1)*adim]) * srm[26 + (j-1)*adim] + migm[26][j-1];
		popf[26 + j*adim] = (popf[26 + (j-1)*adim] + popf[25 + (j-1)*adim]) * srf[26 + (j-1)*adim] + migf[26][j-1];
		
		/* birth in 5-yrs */
		for(i=3; i<10; ++i) {
			bt[i-3] = (popf[i + (j-1)*adim] + popf[i + j*adim]) * asfr[i-3 + (j-1)*7] * 0.5;
		}
		b = sum(bt, 7);
		bm = b * srb[j - 1] / (1 + srb[j - 1]);
		bf = b / (1 + srb[j - 1]);
		
		popm[j*adim] = bm * srm[(j-1)*adim] + mmult * migm[0][j-1];
		popf[j*adim] = bf * srf[(j-1)*adim] + mmult * migf[0][j-1];
		for(i=0; i<27; ++i) {
			totp[j] += popm[i + j*adim]+popf[i + j*adim];
		}
		/*Rprintf("\ntotp%d = %lf", j, totp[j]);*/
	}
}	

