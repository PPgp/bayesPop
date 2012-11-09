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
				double *LLm, double *L10) {
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
	lm[1] = 1 - mxm[0] / (1 + (1 - am[0]) * mxm[0]); /*add multiplication by lm[0] */
	lm[2] = (1 - 4 * mxm[1] / (1 + (4 - am[1]) * mxm[1]));
	lm[2] = lm[1] * lm[2];
	L10[0] =  lm[1] + am[0] * (lm[0] - lm[1]); /* 0L1 */
	LLm[0] = L10[0] + (4 * lm[2] + am[1] * (lm[1] - lm[2])); /* Life table 4L0*/
    
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

double get_constrained_mortality(double a, double b, double k, double constraint) {
	double mx;
	
	mx = exp(a + b*k);
	if(constraint > 0 && mx < constraint) mx = constraint;
	return mx;
}

void LCEoKtC(int sex, double *ax, double *bx, 
			 double eop, double kl, double ku, double *constraints, double *LLm, double *L10, double *Mx) {
	double LTl[27], LTu[27], mxm[28], LTeo, k2;
	int i, dim, debug;
	dim = 27;
	debug=0;
	/* check if the eop lies outside of the bounds */
	for (i=0; i < 28; ++i) {
		mxm[i] = get_constrained_mortality(ax[i], bx[i], kl, constraints[i]);
	}
	LifeTableC(sex, mxm, LTl, L10);
	
	if(eop < sum(LTl, dim)) {
		LLm = LTl;
		return;
	}
	for (i=0; i < 28; ++i) {
		mxm[i] = get_constrained_mortality(ax[i], bx[i], ku, constraints[i]);
	}
	LifeTableC(sex, mxm, LTu, L10);

	if(eop > sum(LTu, dim)) {
		LLm = LTu;
		if(debug==1) Rprintf("\nBreturn %f", sum(LTu, dim));
		return;
	}
	/* Bi-section method */
	k2 = 0.5 * (kl + ku);
	for (i=0; i < 28; ++i) {
		mxm[i] = get_constrained_mortality(ax[i], bx[i], k2, constraints[i]);
	}
	LifeTableC(sex, mxm, LLm, L10);
	if(debug==1) Rprintf("\nLLm[2]=%lf, k2=%f, kl=%f, ku=%f, mxm24-27=%f %f %f %f", LLm[2], k2, kl, ku, mxm[24], mxm[25], mxm[26], mxm[27]);
	LTeo = sum(LLm, dim);
	while(fabs(LTeo - eop) > 0.01) {
		if(LTeo < eop) kl = k2;
		else ku = k2;
		k2 = 0.5 * (kl + ku);
		for (i=0; i < 28; ++i) {
			mxm[i] = get_constrained_mortality(ax[i], bx[i], k2, constraints[i]);
		}
		LifeTableC(sex, mxm, LLm, L10);
		LTeo = sum(LLm, dim);
		if(debug==1) Rprintf("\nLTeo=%lf, dif=%lf, LLm0-2=%lf %lf %lf, k2=%f, kl=%f, ku=%f, mxm24-27=%f %f %f %f", LTeo, fabs(LTeo - eop), LLm[0], LLm[1], LLm[2], k2, kl, ku, mxm[24], mxm[25], mxm[26], mxm[27]);
	}
	if(debug==1) Rprintf("\nk2=%f, eop=%lf, LTeo=%lf, adif=%lf, LLm[0]=%lf", k2, eop, LTeo, fabs(LTeo - eop), LLm[0]);
	for (i=0; i < 28; ++i) Mx[i] = mxm[i];
}

void get_sx(double *LLm, double *sx, int n) {
	int i, l, oei;
	double sumLL;
	l = 27;
	oei=n-1;
	/* Survival Ratios */
    sx[0] = LLm[0] / 5.0;
	for(i=1; i < oei; ++i) {
		if(LLm[i-1] == 0) sx[i] = 0;
		else sx[i] = LLm[i]/LLm[i-1];
	}
	/* Last age group */
	sumLL = 0;
	for(i=oei-1; i < (l-1); ++i) {
		sumLL += LLm[i];
	}
	if((sumLL + LLm[l-1]) == 0 ||  sumLL == 0) sx[oei] = 0;
	else sx[oei] = sumLL/(sumLL+LLm[l-1]);
	if(sx[oei] > sx[oei-1]) sx[oei] = sx[oei-1];
}

void get_sx27(double *LLm, double *sx) {
	get_sx(LLm, sx, 27);
}

void get_sx21(double *LLm, double *sx) {
	get_sx(LLm, sx, 21);
}


void LC(int *Npred, int *Sex, double *ax, double *bx, 
		double *Eop, double *Kl, double *Ku, int *constrain, double *FMx, double *FEop, double *LLm, double *Sr, 
		double *L10, double *Mx) {
	double eop, kl, ku, sx[27], LL10[1], Lm[26], mxm[28], fmx[28];
	int i, sex, npred, pred;
	
	npred = *Npred;
	sex=*Sex;
	ku=*Ku;
	kl=*Kl;
	for (i=0; i < 28; ++i) fmx[i] = -1;
	for (pred=0; pred < npred; ++pred) {
		eop = Eop[pred];
		if(*constrain>0) {		
			if(FEop[pred] > eop) {
				for (i=22; i < 28; ++i) {fmx[i] = FMx[i + pred*28];}
			} else {
				for (i=22; i < 28; ++i) {fmx[i] = -1;}
			}
		}
		/*Rprintf("\n%i: eop=%lf", pred, eop);*/
		LCEoKtC(sex, ax, bx, eop, kl, ku, fmx, Lm, LL10, mxm);
		get_sx27(Lm, sx);

		for (i=0; i < 27; ++i) {
			Sr[i + pred*27] = sx[i];
			/*Rprintf("\nLLm=%lf, Sr=%lf", LLm[i], Sr[i + pred*27]);*/
			Mx[i + pred*28] = mxm[i];
		}
		Mx[27 + pred*28] = mxm[27];
		for (i=0; i < 26; ++i) {
			LLm[i + pred*26] = Lm[i];
		}
		L10[pred] = LL10[0];
	}
}

void get_sr_from_LT(int *N, int *Sex, double *Mx, double *Sr, double *L10) {
	double LLm[26], sx[21], mxm[28], LL10[1];
	int i, j;
	for (j=0; j < *N; ++j) {
		/*Rprintf("\nj=%i: mxm = ", j);*/
		for (i=0; i < 28; ++i) {
			mxm[i] = Mx[i + j*28];
			/*Rprintf(" %lf", mxm[i]);*/
		}
		LifeTableC(*Sex, mxm, LLm, LL10);
		get_sx21(LLm, sx);
		/*Rprintf("\n   sx = ");*/
		for (i=0; i < 21; ++i) {
			Sr[i + j*21] = sx[i];
			/*Rprintf(" %lf", sx[i]);*/
		}
		L10[j] = LL10[0];
	}
}

void get_sr_from_N(int *N, double *Pop, double *MIG, int *MIGtype, double *Births, double *Sr, double *Deaths) {
	double mig[21][*N], totmig[21][*N], mmult;
	int i, j, nrow, n;
	nrow = 21;
	n = *N;
	mmult = 1;
	if(*MIGtype == 0) mmult=0.5;
	for(j=0; j<n; ++j) {
		/* age < 5 */
		Sr[j*nrow] = (Pop[(j+1)*nrow] - mmult * MIG[j*nrow])/Births[j];
		Deaths[j*nrow] = Pop[(j+1)*nrow]*(1-Sr[j*nrow]);
		for(i=1; i<(nrow-1); ++i) {
			switch (*MIGtype) {
				case 0: /* migration evenly distributed over each interval (MigCode=0) */
					/* age >= 5 */
					Sr[i+j*nrow] = (Pop[i+(j+1)*nrow] - 0.5*MIG[i+j*nrow])/(Pop[i-1+j*nrow] + 0.5*MIG[i-1+j*nrow]);
					Deaths[i+j*nrow] = Pop[i-1+j*nrow]*(1-Sr[i+j*nrow]) + 0.5*MIG[i-1+j*nrow]*(1-Sr[i+j*nrow]);
					/*Rprintf(" %i %i: %lf %lf %lf\n", j, i, Pop[i-1+j*nrow]*(1-Sr[i+j*nrow]), 0.5*MIG[i-1+j*nrow]*(1-Sr[i+j*nrow]), Deaths[i+j*nrow]);*/
					/*totmig[i][j] = 0.5*(MIG[i + j*nrow] + MIG[i-1 + j*nrow]*Sr[i-1 + j*nrow]);*/
					break;
				default: /* migration at the end of each interval (MigCode=9)*/
					/* age >= 5 */
					Sr[i+j*nrow] = (Pop[i+(j+1)*nrow] - MIG[i+j*nrow])/Pop[i-1+j*nrow];
					Deaths[i+j*nrow] = Pop[i-1+j*nrow]*(1-Sr[i+j*nrow]);
					/*totmig[i][j] = MIG[i + j*nrow];*/
					break;
			}
		}
		/* Last open-ended age category */
		Sr[nrow-1+j*nrow] = (Pop[nrow-1+(j+1)*nrow] - MIG[nrow-1+j*nrow])/(Pop[nrow-2+j*nrow]+Pop[nrow-1+j*nrow]);
		Deaths[nrow-1+j*nrow] = Pop[nrow-2+j*nrow]*(1-Sr[i+j*nrow]);
	}
}

void TotalPopProj(int *npred, double *MIGm, double *MIGf, int *migr, int *migc,
				  int *MIGtype, double *srm, double *srf, double *asfr, double *srb, 
				  double *popm, double *popf, double *totp, 
				  double *btagem, double *btagef, double *deathsm, double *deathsf
					) {
	double migm[*migr+6][*migc], migf[*migr+6][*migc], totmigm[*migr+6][*migc], totmigf[*migr+6][*migc];
	double b, bt[7], bm, bf, mmult, srb_ratio, tmp;
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
			/* Time index (j) of survival ratio, migration and vital events is shifted by one in comparison to population,
			   i.e. pop[0] is the current period, whereas sr[0], mig[0] etc. is the first projection period.*/
			popm[i + j*adim] = popm[i-1 + (j-1)*adim] * srm[i + (j-1)*adim] + totmigm[i][j-1];
			popf[i + j*adim] = popf[i-1 + (j-1)*adim] * srf[i + (j-1)*adim] + totmigf[i][j-1];
		}
		
		
		/* Age 130+ */
		popm[26 + j*adim] = (popm[26 + (j-1)*adim] + popm[25 + (j-1)*adim]) * srm[26 + (j-1)*adim] + migm[26][j-1];
		popf[26 + j*adim] = (popf[26 + (j-1)*adim] + popf[25 + (j-1)*adim]) * srf[26 + (j-1)*adim] + migf[26][j-1];
		
		/* birth in 5-yrs */
		srb_ratio = srb[j - 1] / (1 + srb[j - 1]);
		for(i=3; i<10; ++i) {
			bt[i-3] = (popf[i + (j-1)*adim] + popf[i + j*adim]) * asfr[i-3 + (j-1)*7] * 0.5;
			btagem[i-3+(j-1)*7] = bt[i-3] * srb_ratio;
			btagef[i-3+(j-1)*7] = bt[i-3] - btagem[i-3+(j-1)*7];
		}
		b = sum(bt, 7);
		bm = b * srb_ratio;
		bf = b / (1 + srb[j - 1]);
		/* age 0-4 */	
		popm[j*adim] = bm * srm[(j-1)*adim] + mmult * migm[0][j-1];
		popf[j*adim] = bf * srf[(j-1)*adim] + mmult * migf[0][j-1];
		/* get total for all ages */
		for(i=0; i<adim; ++i) {
			totp[j] += popm[i + j*adim]+popf[i + j*adim];
			/*deathsm[i + (j-1)*adim] = popm[i + j*adim] * (1-srm[i + (j-1)*adim]);
			deathsf[i + (j-1)*adim] = popf[i + j*adim] * (1-srf[i + (j-1)*adim]);*/
		}
		deathsm[(j-1)*adim] = bm * (1-srm[(j-1)*adim]);
		deathsf[(j-1)*adim] = bf * (1-srf[(j-1)*adim]);
		for(i=1; i<(adim-1); ++i) {
			deathsm[i + (j-1)*adim] = popm[i-1 + (j-1)*adim]*(1-srm[i + (j-1)*adim]);
			deathsf[i + (j-1)*adim] = popf[i-1 + (j-1)*adim]*(1-srf[i + (j-1)*adim]);
		}
		i = 26;
		deathsm[i + (j-1)*adim] = (popm[i + (j-1)*adim]+popm[i-1 + (j-1)*adim])*(1-srm[i + (j-1)*adim]);
		deathsf[i + (j-1)*adim] = (popf[i + (j-1)*adim]+popf[i-1 + (j-1)*adim])*(1-srf[i + (j-1)*adim]);

		/*Rprintf("\ntotp%d = %lf", j, totp[j]);*/
	}
	/*Rprintf("\n\n\n");
	for(j=0; j<n; ++j) {
		Rprintf("\nj=%i: ",j);
		tmp = 0;
		for(i=0; i<adim; ++i) {
			Rprintf("%lf\t", deathsm[i + j*adim]+deathsf[i + j*adim]);
			tmp += deathsm[i + j*adim]+deathsf[i + j*adim];
		}
		Rprintf("\nj=%i sum deaths=%lf: ",j, tmp);
	}*/
}	

