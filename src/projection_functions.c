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

/*****************************************************************************
 * Function returns the ax for ages 0 and 1-4 (abridged life table)
 * based on mx(0) and sex
 * Formulas re-estimates from the Coale/Demeny separation factors
 * on mx(0) insted of qx(0)
 * Source from Preston et al. 2001, p.48
 *****************************************************************************/
double * get_a05(double mx0, int sex) {
	static double ax[2];
	if(sex > 1) {/* female*/
		if (mx0 < 0.107) {
			ax[0] = 0.053 + 2.8 * mx0;      /*1a0*/
			ax[1] = 1.522 - 1.518 * mx0;    /*4a1*/
		} else {
			ax[0] = 0.35;
			ax[1] = 1.361;
		}
	} else { /* male */
		if (mx0 < 0.107) {
			ax[0] = 0.045 + 2.684 * mx0;
			ax[1] = 1.651 - 2.816 * mx0;
		} else {
			ax[0] = 0.33;
			ax[1] = 1.352;
		}
	}
	return(ax);
}
/*****************************************************************************
 * Function calculates an abridged life table from age-specific mortality 
 * rates
 * Input: 
 * * mx   age specific mortality rates 
 * * sex  sex
 * * nage number of age groups
 * Output: 
 * qx Probabilities of dying
 * lx Survivors to exact age
 * Lx Person years lived
 * ax proportion of years lived by those who died
 *****************************************************************************/
void doLifeTable(int sex, int nage, double *mx, 
				double *Lx, double *lx, double *qx, double *ax) {
	
	int i;
  double k;     /* correcting factor in Greville approximation */
	double *tmpa; /* pointer to estimated ax[0] and ax[1] values */
	int nage1;
	nage1 = nage -1;
	
	tmpa = get_a05(mx[0], sex);
	ax[0] = tmpa[0];
	ax[1] = tmpa[1];
	qx[0] = mx[0] / (1 + (1 - ax[0]) * mx[0]);                    /* 1q0 */	
	qx[1] = 4 * mx[1] / (1 + (4 - ax[1]) * mx[1]);				        /* 4q1 */	
	lx[0] = 1;                                                    /* l0 */
	lx[1] = lx[0] * (1 - qx[0]);                                  /* l1 */
	lx[2] = lx[1] * (1 - qx[1]);                                  /* l5 = l1 * (1-4q1) */
	Lx[0] =  lx[1] + ax[0] * (lx[0] - lx[1]);                     /* 1L0 */
	Lx[1] =  4 * lx[2] + ax[1] * (lx[1] - lx[2]);                 /* 4L1 */
		
  /*Rprintf("\nnage=%i, L0=%f, ax0-1=%f %f, l1-2=%f %f, mx0-1=%f %f", nage, Lx[0], ax[0], ax[1], lx[1], lx[2], mx[0], mx[1]);*/
  /* Age 5-9, .... 125-129 
  Greville formula used in Mortpak and UN MLT (1982)*/
  /* TB: corrected for age group 3 (5-9), tentativley set to 2.5 or n/2 */
  /*  ax[2] = 2.5; */
  
  k     = 0.1 * log(fmax(mx[4] / fmax(mx[2], DBL_MIN), DBL_MIN));
  ax[2] = 2.5 - (25 / 12.0) * (mx[2] - k);
  
  for(i = 3; i < nage1; ++i) {
    k     = 0.1 * log(fmax(mx[i+1] / fmax(mx[i-1], DBL_MIN), DBL_MIN));
    ax[i] = 2.5 - (25 / 12.0) * (mx[i] - k);
    }
    /* penultimate ax calculated with k from previous age group */
    ax[nage1] = 2.5 - (25 / 12.0) * (mx[i] - k);
  
  /* correcting out-of (reasonable) bounds ax for older ages             */ 
  /* 0.97=1-5*exp(-5)/(1-exp(-5)), for constant mu=1, Kannisto assumption*/
  for(i = 10; i < nage; i++) {
    if(ax[i] < 0.97) {
       ax[i] = 0.97;
      }
    }
  
  /* caculate life table variables from mx and ax */
	for(i = 2; i < nage; ++i) {		
		/*Rprintf("ax%i=%f, mx%i=%f", i, ax[i], i-1, mx[i-1]);*/
		qx[i] = 5 * mx[i] / (1 + (5 - ax[i]) * mx[i]);
		lx[i+1] = lx[i] * (1-qx[i]);
		Lx[i] = 5 * lx[i+1] + ax[i] * (lx[i] - lx[i+1]);
	}
	
	/* Open ended age interval */
	Lx[nage] = lx[nage] / fmax(mx[nage], DBL_MIN); /* Assuming Mx levels off at age 130 */
	qx[nage] = 1.0;
	
	/* TB: added missing ax for last, open-ended age group */ 
  /*  worked after declaration in predict.pop was changed */
	ax[nage] = Lx[nage];
	
	/*Rprintf("\nLTend\n");*/
}

/* Function returns collapsed Lx and lx columns of life table*/
void LifeTableC(int sex, int nage, double *mxm, 
				double *LLm, double *lm) {
	double ax[27], qx[28], L[28];
	int i;
	doLifeTable(sex, nage, mxm, L, lm, qx, ax);
	/* collapse 1L0 and 4L1 into 5L0 */
	LLm[0] = L[0] + L[1];
	for(i = 1; i < nage; ++i) {
		LLm[i] = L[i+1];
	}			
}



double get_constrained_mortality(double a, double b, double k, double constraint) {
	double mx;
	
	mx = exp(a + b*k);
	if(constraint > 0 && mx < constraint) mx = constraint;
	return mx;
}

void LCEoKtC(int sex, double *ax, double *bx, 
			 double eop, double kl, double ku, double *constraints, double *LLm, double *lm, double *Mx) {
	double LTl[27], LTu[27], mxm[28], LTeo, k2;
	int i, dim, debug;
	dim = 27;
	debug=0;
	/* check if the eop lies outside of the bounds */
	for (i=0; i < 28; ++i) {
		mxm[i] = get_constrained_mortality(ax[i], bx[i], kl, constraints[i]);
	}
	LifeTableC(sex, 27, mxm, LTl, lm);
	
	if(eop < sum(LTl, dim)) {
		for (i=0; i < dim; ++i) LLm[i]=LTl[i];
		for (i=0; i < 28; ++i) Mx[i] = mxm[i];
		return;
	}
	for (i=0; i < 28; ++i) {
		mxm[i] = get_constrained_mortality(ax[i], bx[i], ku, constraints[i]);
	}
	LifeTableC(sex, 27, mxm, LTu, lm);

	if(eop > sum(LTu, dim)) {
		for (i=0; i < dim; ++i) LLm[i]=LTu[i]; 
		for (i=0; i < 28; ++i) Mx[i] = mxm[i];
		if(debug==1) Rprintf("\nBreturn %f", sum(LTu, dim));
		return;
	}
	/* Bi-section method */
	k2 = 0.5 * (kl + ku);
	for (i=0; i < 28; ++i) {
		mxm[i] = get_constrained_mortality(ax[i], bx[i], k2, constraints[i]);
	}
	LifeTableC(sex, 27, mxm, LLm, lm);
	if(debug==1) Rprintf("\nLLm[2]=%lf, k2=%f, kl=%f, ku=%f, mxm24-27=%f %f %f %f", LLm[2], k2, kl, ku, mxm[24], mxm[25], mxm[26], mxm[27]);
	LTeo = sum(LLm, dim);
	while(fabs(LTeo - eop) > 0.01) {
		if(LTeo < eop) kl = k2;
		else ku = k2;
		k2 = 0.5 * (kl + ku);
		for (i=0; i < 28; ++i) {
			mxm[i] = get_constrained_mortality(ax[i], bx[i], k2, constraints[i]);
		}
		LifeTableC(sex, 27, mxm, LLm, lm);
		LTeo = sum(LLm, dim);
		if(debug==1) Rprintf("\nLTeo=%lf, dif=%lf, LLm0-2=%lf %lf %lf, k2=%f, kl=%f, ku=%f, mxm24-27=%f %f %f %f", LTeo, fabs(LTeo - eop), LLm[0], LLm[1], LLm[2], k2, kl, ku, mxm[24], mxm[25], mxm[26], mxm[27]);
	}
	if(debug==1) Rprintf("\nk2=%f, eop=%lf, LTeo=%lf, adif=%lf, LLm[0]=%lf", k2, eop, LTeo, fabs(LTeo - eop), LLm[0]);
	for (i=0; i < 28; ++i) Mx[i] = mxm[i];
}

void get_sx(double *LLm, double *sx, int n, int Ldim) {
	/* compute survival ratios from Lx where the first age group is 0-5 (also for Lx)*/
	int i, oei;
	double sumLL;
	oei=n-1;
	/* Survival Ratios, radix of life table assumed to be 1.0 */
	sx[0] = LLm[0] / 5.0;
	for(i=1; i < oei; ++i) {
		if(LLm[i-1] == 0) sx[i] = exp(-5);
		else sx[i] = LLm[i]/LLm[i-1];
	}
	/* Last age group */
	sumLL = 0;
	for(i=oei; i < Ldim; ++i) {
		sumLL += LLm[i];
	}
	if((sumLL + LLm[oei-1]) == 0 ||  sumLL == 0) sx[oei] = exp(-5);
	else sx[oei] = sumLL/(sumLL+LLm[oei-1]);
	if(sx[oei] > sx[oei-1]) sx[oei] = sx[oei-1];
}

void get_sx27(double *LLm, double *sx) {
	get_sx(LLm, sx, 27, 27);
}

void get_sx21(double *LLm, double *sx) {
	get_sx(LLm, sx, 21, 27);
}

void get_sx21_21(double *LLm, double *sx) {
	get_sx(LLm, sx, 21, 21);
}

/*****************************************************************************
 * Wrapper for abridged life table function, 
 * called from R in predict.pop
 * TB: adding missing life table variables  Tx, sx and dx
 *     and correcting error in sx (larger than 1, shifted index)
 *****************************************************************************/

void LifeTable(int *sex, int *nage, double *mx, 
				double *Lx, double *lx, double *qx, double *ax, double *Tx, double *sx, double *dx) {
	int i;
  
  /* original implementation lacking Tx, sx, dx */
	doLifeTable(*sex, *nage, mx, Lx, lx, qx, ax);
	
	/* calculating additional life table columns dx, Tx, sx */
	
	/* TB: dx */
	for(i = 0; i < *nage; ++i) {
	  dx[i] = lx[i] - lx[i+1];
	}
	dx[*nage] = lx[*nage];
	
	/*TB: Tx */
	Tx[*nage] = Lx[*nage];
	for (i = *nage-1; i >= 0; i--) {
	  Tx[i] = Tx[i+1] + Lx[i];
	}	
	/* TB: sx */
	/* first age group is survival from births to age 0-5 */
	sx[0] = (Lx[0] + Lx[1]) / 5*lx[0];
	/* second age group is survival age 0-5 to age 5 - 10 */
	sx[1] = Lx[2] / (Lx[0] + Lx[1]);
	/* middle age groups  */
	for(i = 2; i < *nage-1; ++i) {
	  sx[i] = Lx[i+1] / Lx[i];
	}
	/* last but one age group */
	sx[*nage-1] = Lx[*nage] / (Lx[*nage-1]+Lx[*nage]);
	sx[*nage]= 0.0;
	
	/* obsolete
	Tx[*nage] = Lx[*nage];
	dx[*nage] = lx[*nage];
	for (i = *nage-1; i >= 0; i--) {
		Tx[i] = Tx[i+1] + Lx[i];
		dx[i] = lx[i] - lx[i+1];
	}	
	get_sx(Lx, sx, *nage, *nage);
	*/
	/* the above call assumed the first age category is 0-5 but here it is 0-1 */
	/*
	sx[0] = sx[0]*5;
	*/
}


void LC(int *Npred, int *Sex, double *ax, double *bx, 
		double *Eop, double *Kl, double *Ku, int *constrain, double *FMx, double *FEop, double *LLm, double *Sr, 
		double *lx, double *Mx) {
	double eop, sx[27], Lm[27], mxm[28], fmx[28], lm[28], locbx[28], locax[28];
	int i, sex, npred, pred;
	
	npred = *Npred;
	sex=*Sex;
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
		for (i=0; i < 28; ++i) {
			locbx[i] = bx[i + pred*28];
			locax[i] = ax[i + pred*28];
		}
		/*Rprintf("\n%i: eop=%lf", pred, eop);*/
		LCEoKtC(sex, locax, locbx, eop, Kl[pred], Ku[pred], fmx, Lm, lm, mxm);		
		get_sx27(Lm, sx);
		
		for (i=0; i < 27; ++i) {
			Sr[i + pred*27] = sx[i];
			/*Rprintf("\nLLm=%lf, Sr=%lf", LLm[i], Sr[i + pred*27]);*/
			Mx[i + pred*28] = mxm[i];
			lx[i + pred*28] = lm[i];
		}
		Mx[27 + pred*28] = mxm[27];
		lx[27 + pred*28] = lm[27];
		for (i=0; i < 26; ++i) {
			LLm[i + pred*27] = Lm[i];
		}
	}
}

void get_deaths_from_sr(double *Sr, int *N, double *Pop, double *MIG, int *MIGtype, double *Births, double *Deaths) {
	int i, j, nrow, n;
	nrow = 21;
	n = *N;
	for(j=0; j<n; ++j) {
		/* age < 5 */
		Deaths[j*nrow] = Pop[(j+1)*nrow]*(1-Sr[j*nrow]); /* this includes migration */
		/*Deaths[j*nrow] = Births[j*nrow]*(1-Sr[j*nrow]);*/
		for(i=1; i<(nrow-1); ++i) {
			switch (*MIGtype) {
				case 0: /* migration evenly distributed over each interval (MigCode=0) */
					/* age >= 5 */
					Deaths[i+j*nrow] = Pop[i-1+j*nrow]*(1-Sr[i+j*nrow]) + 0.5*MIG[i-1+j*nrow]*(1-Sr[i+j*nrow]);
					break;
				default: /* migration at the end of each interval (MigCode=9)*/
					/* age >= 5 */
					Deaths[i+j*nrow] = Pop[i-1+j*nrow]*(1-Sr[i+j*nrow]);
					break;
			}
		}
		/* Last open-ended age category */
		Deaths[nrow-1+j*nrow] = Pop[nrow-2+j*nrow]*(1-Sr[nrow-1+j*nrow]);
	}
}

void get_sr_from_N(int *N, double *Pop, double *MIG, int *MIGtype, double *Births, double *Sr, double *Deaths) {
	/* function not used */
	double mmult;
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
					break;
				default: /* migration at the end of each interval (MigCode=9)*/
					/* age >= 5 */
					Sr[i+j*nrow] = (Pop[i+(j+1)*nrow] - MIG[i+j*nrow])/Pop[i-1+j*nrow];
					Deaths[i+j*nrow] = Pop[i-1+j*nrow]*(1-Sr[i+j*nrow]);
					break;
			}
		}
		/* Last open-ended age category */
		Sr[nrow-1+j*nrow] = (Pop[nrow-1+(j+1)*nrow] - MIG[nrow-1+j*nrow])/(Pop[nrow-2+j*nrow]+Pop[nrow-1+j*nrow]);
		Deaths[nrow-1+j*nrow] = Pop[nrow-2+j*nrow]*(1-Sr[nrow-1+j*nrow]);
	}
}
/*****************************************************************************
 * Core population projection function TotalPopProjection
 * Called from function StoPopProj in predict.pop.R
 * Parameter
 * int *npred                 number of prediction intervals (? time points)
 * double *MIGm, *MIGf        male, female migration by age
 * int *migr, *migc           rows, columns of migration array
 * int *MIGtype               type of migration adjustement
 * double *srm, *srf          male, female survivor ration by age
 * double *asfr               age-specific fertility rates 
 * double *srb                sex ratio at birth
 * double *popm, *popf        male, female population by age
 * double *totp               total population by age 
 * double *btagem, *btagef    male, female births by age of mother
 * double *deathsm, *deathsf  male, female deaths by age
******************************************************************************/
void TotalPopProj(int *npred, double *MIGm, double *MIGf, int *migr, int *migc,
				  int *MIGtype, double *srm, double *srf, double *asfr, double *srb, 
				  double *popm, double *popf, double *totp, 
				  double *btagem, double *btagef, double *deathsm, double *deathsf
					) {
	double migm[*migr+6][*migc], migf[*migr+6][*migc], totmigm[*migr+6][*migc], totmigf[*migr+6][*migc];
	double b, bt[7], bm, bf, mmult, srb_ratio;
	int i,j, jve, adim, nrow, ncol, n;
	nrow = *migr;
	ncol = *migc;
	n = *npred;
	
	/* adim is an offset representing time, necessary to access and store two-dimensional data 
	organised by age and time/period into a one-dimensional vector requiered by accessing c from R
	the expression 
	 popm[i + j*adim] = popm[i-1 + (j-1)*adim] * srm[i + jve*adim];
	 reads forj = 1; j = 1; adim = 27
	 popm[1 + 1*adim] = popm[0 + 0*adim] * srm[1 + 1*adim];
	 popm[1 + 27] = popm[0 + 0] * srm[1 + 27];
	 the sam in multidimensional fashion:
 	 popm[1][t+1] = popm[0][t] * srm[1][t+1]
	 TB: By combining the population at time t=0 with survivor ratios from time t = 1, there could be 
	     an unwanted shift in projecting; or the survivor ratios for time t=0 must be empty/not defined 
	*/ 
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
		jve = j-1;
		/* Time index (j) of survival ratio, migration and vital events is shifted by one in comparison to population,
			   i.e. pop[0] is the current period, whereas sr[0], mig[0] etc. is the first projection period.*/
		/* Compute ages >=5 */
    
    /*TB: replace (j-1) with jve where appropriate */	
    /*    further simplification would mean to replace jve*adim 
          with another variable toffset (= jve*adim)*/
		for(i=1; i<(adim-1); ++i) {		
			popm[i + j*adim] = popm[i-1 + (j-1)*adim] * srm[i + jve*adim];
			popf[i + j*adim] = popf[i-1 + (j-1)*adim] * srf[i + jve*adim];
			totmigm[i][jve] = fmax(totmigm[i][jve], -1*popm[i + j*adim]); /* assures population is not negative */
			popm[i + j*adim] = popm[i + j*adim] + totmigm[i][jve];
			totmigf[i][jve] = fmax(totmigf[i][jve], -1*popf[i + j*adim]);
			popf[i + j*adim] = popf[i + j*adim] + totmigf[i][jve];
		}
		/* Age 130+ */
		popm[26 + j*adim] = (popm[26 + (j-1)*adim] + popm[25 + (j-1)*adim]) * srm[26 + jve*adim];
		popf[26 + j*adim] = (popf[26 + (j-1)*adim] + popf[25 + (j-1)*adim]) * srf[26 + jve*adim];
		totmigm[26][jve] = fmax(migm[26][jve], -1*popm[26 + j*adim]);
		popm[26 + j*adim] = popm[26 + j*adim] + totmigm[26][jve];
		totmigf[26][jve] = fmax(migf[26][jve], -1*popf[26 + j*adim]);
		popf[26 + j*adim] = popf[26 + j*adim] + totmigf[26][jve];
		/* birth in 5-yrs */
		srb_ratio = srb[jve] / (1 + srb[jve]);
		for(i=3; i<10; ++i) {
			bt[i-3] = (popf[i + (j-1)*adim] + popf[i + j*adim]) * asfr[i-3 + jve*7] * 0.5;
			btagem[i-3+jve*7] = bt[i-3] * srb_ratio;
			btagef[i-3+jve*7] = bt[i-3] - btagem[i-3+jve*7];
		}
		b = sum(bt, 7);
		bm = b * srb_ratio;
		bf = b / (1 + srb[jve]);
		/* age 0-4 */	
		popm[j*adim] = bm * srm[jve*adim];
		popf[j*adim] = bf * srf[jve*adim];
		totmigm[0][jve] = fmax(mmult * migm[0][jve], -1*popm[j*adim]);
		popm[j*adim] = popm[j*adim] + totmigm[0][jve];
		totmigf[0][jve] = fmax(mmult * migf[0][jve], -1*popf[j*adim]);
		popf[j*adim] = popf[j*adim] + totmigf[0][jve];
		
		/* get total for all ages */
		for(i=0; i<adim; ++i) {
			totp[j] += popm[i + j*adim]+popf[i + j*adim];
		}
		deathsm[jve*adim] = bm * (1-srm[jve*adim]);
		deathsf[jve*adim] = bf * (1-srf[jve*adim]);                
		for(i=1; i<(adim-1); ++i) {
			deathsm[i + jve*adim] = popm[i-1 + (j-1)*adim]*(1-srm[i + jve*adim]);
			deathsf[i + jve*adim] = popf[i-1 + (j-1)*adim]*(1-srf[i + jve*adim]);
		}
		i = 26;
		deathsm[i + jve*adim] = (popm[i + (j-1)*adim]+popm[i-1 + (j-1)*adim])*(1-srm[i + jve*adim]);
		deathsf[i + jve*adim] = (popf[i + (j-1)*adim]+popf[i-1 + (j-1)*adim])*(1-srf[i + jve*adim]);
	}	
}	

