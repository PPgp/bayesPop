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
 temporary function: prints content of an array to console
 *****************************************************************************/
void printArray(double *a, int count) {
	int i;
	for (i = 0; i < count; i++) {
		Rprintf("\n[%i] = %18.15f", i, a[i]);
	}
	Rprintf("\n");
}

/*****************************************************************************
 * Function returns the years lived by those who died in the age interval (ax)
 * for ages 0 and 1-4 (abridged life table)
 * based on mx(0,1) and sex
 * Formulas re-estimates from the Coale/Demeny separation factors
 * on mx(0,1) insted of qx(0,1)
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
	qx[1] = 4 * mx[1] / (1 + (4 - ax[1]) * mx[1]);				  /* 4q1 */	
	lx[0] = 1;                                                    /* l0 */
	lx[1] = lx[0] * (1 - qx[0]);                                  /* l1 */
	lx[2] = lx[1] * (1 - qx[1]);                                  /* l5 = l1 * (1-4q1) */
	Lx[0] =  lx[1] + ax[0] * (lx[0] - lx[1]);                     /* 1L0 */
	Lx[1] =  4 * lx[2] + ax[1] * (lx[1] - lx[2]);                 /* 4L1 */
		
    /*Rprintf("\nnage=%i, L0=%f, ax0-1=%f %f, l1-2=%f %f, mx0-1=%f %f", nage, Lx[0], ax[0], ax[1], lx[1], lx[2], mx[0], mx[1]);*/
    /* Age 5-9, .... 125-129 
	 Greville formula used in Mortpak and UN MLT (1982)*/
	/* TB: corrected for age group 3 (5-9), tentatively set to 2.5 or n/2 */
	/*  ax[2] = 2.5; */
																			  	/* test Mortpak/Abacus compatibility rule for age 5-9 and 10-14 with fixed ax
	k     = 0.1 * log(fmax(mx[4] / fmax(mx[2], DBL_MIN), DBL_MIN));
	ax[2] = 2.5 - (25 / 12.0) * (mx[2] - k);
	*/
	ax[2] = 2.5;
	ax[3] = 2.5;

	/*  test Mortpak/Abacus compatibility rule for age 5-9 and 10-14 with fixed ax
	for(i = 3; i < nage1; ++i) {
	*/
	for(i = 4; i < nage1; ++i) {
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



/* Function returns collapsed Lx and lx columns of life table */
/* function calls doLifeTable first, then collapsesLx and lx  */
/* TB: added missing collapsing of lx column                  */
/*     candidate for renaming parameters for consistency      */
/*     also check dimensioning of parameter                  */
void LifeTableC(int sex, int nage, double *mx,
                double *Lxx, double *lxx) {
  /* life table variables returned from doLifeTable */
  /* need to be declared with nage+1 elements       */
  /* Lx[nage+1]
   * lx[nage+1]
   * qx[nage+1]
   * ax[nage+1];
   * input todoLifeTableC mx is declared outside
   * output from LifeTableC Lxx, lxx is declared 
   * with nage-1 elements outside 
   */
  
  double Lx[nage+1], lx[nage+1], qx[nage+1], ax[nage+1];
  int i;
  
  /* do life table called with nage as last index */
  doLifeTable(sex, nage, mx, Lx, lx, qx, ax);
  /* collapse 1L0 and 4L1 into 5L0 */
  Lxx[0] = Lx[0] + Lx[1];
  lxx[0] = lx[0];
  for(i = 1; i < nage; ++i) {
    Lxx[i] = Lx[i+1];
    lxx[i] = lx[i+1];
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
	/* Survival Ratios, radix of life table assumed to be 1.0  */
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
 * Lee Carter model
 * Produces a projection of age -specific mortality rates
 * (more)
 * 
 *****************************************************************************/
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

void get_deaths_from_sr(int *Sex, double *Sr, int *N, double *Pop, double *MIG, int *MIGtype, 
                        double *Births, double *Deaths, double *Mx) {

	double Lx[27], lx[27];
	double cdeaths[27];
	double mxt[28];

	/* forward and backward estimated deaths */
	double dfw, dbw;        
	/* cohort separation factor males, females*/
	double csf[27]; 

	int i, j, nrow, n, sex;
	nrow = 21; /* age */
	n = *N; /* periods */
	sex=*Sex;
	
	/* loop by time period */	
	for(j=0; j<n; ++j) {
		/**************************************************************************/
		/* cohort deaths                                                          */
		/* Calculated by a combination of forward-backward estimation (UN ABACUS) */
		/* dfw deaths by forward projection    ABACUS: PART1                      */
		/* dbw deaths by backward projection   ABACUS: C                          */
		/* 1.  Deaths accuring to births                                          */
		/* 2.  Deaths occuring to closed age groups (middle age groups)           */
		/* 3.  Deaths accuring to last-opended age group                          */
		/**************************************************************************/
		/* age < 5 */
		/* cohort deaths forward no migration */
		dfw = Births[j] * (1-Sr[j*nrow]);
		/* cohort deaths backward with migration */		
		dbw = Pop[(j+1)*nrow] * ((1-Sr[j*nrow]) / Sr[j*nrow]);
		/* average of both types of deaths */
		/* cdeaths[j*nrow] = 0.5* (dfw + dbw); */
		cdeaths[0] = 0.5* (dfw + dbw);
		for(i=1; i<(nrow); ++i) {
			switch (*MIGtype) {
				case 0: /* migration evenly distributed over each interval (MigCode=0) */
					/* age >= 5 */
					/* Deaths[i+j*nrow] = Pop[i-1+j*nrow]*(1-Sr[i+j*nrow]) + 0.5*MIG[i-1+j*nrow]*(1-Sr[i+j*nrow]); */					
					dfw = Pop[i-1 + j*nrow] * (1-Sr[i + j*nrow]);
					dbw = Pop[i + (j+1)*nrow] * ((1-Sr[i + j*nrow]) / Sr[i + j*nrow]);
					/* cdeaths[i+j*nrow] = 0.5 * (dfw + dbw); */
					cdeaths[i] = 0.5 * (dfw + dbw);
					break;

				default: /* migration at the end of each interval (MigCode=9)*/
					/* age >= 5 */
					/* Deaths[i+j*nrow] = Pop[i-1+j*nrow]*(1-Sr[i+j*nrow]); */
					dfw = Pop[i-1 + j*nrow] * (1-Sr[i + j*nrow]);
					dbw = Pop[i + (j+1)*nrow] * ((1-Sr[i + j*nrow]) / Sr[i + j*nrow]);
					/* cdeaths[i+j*nrow] = 0.5 * (dfw + dbw); */
					cdeaths[i] = 0.5 * (dfw + dbw);
					break;
			}
		}
		/* Last open-ended age category */
		/* Deaths[nrow-1+j*nrow] = Pop[nrow-2+j*nrow]*(1-Sr[nrow-1+j*nrow]); */
		dfw = (Pop[nrow + j*nrow] + Pop[nrow-1 + j*nrow]) * (1-Sr[nrow + j*nrow]);		
		dbw = Pop[nrow + (j+1)*nrow] * ((1-Sr[nrow + j*nrow]) / Sr[nrow + j*nrow]);
		/* cdeaths[nrow-1+j*nrow] = 0.5 * (dfw + dbw); */
		cdeaths[nrow] = 0.5 * (dfw + dbw);


		/**************************************************************************/
	    /* period deaths                                                          */
	    /* 1. calculated cohort separation factors                                */
	    /* 2. split cohort deaths into Lexis triangles                            */
	    /* 3. rearrange lexis triangles into period deaths (period-age format)    */
	    /**************************************************************************/
	    
	    for(i=0; i<(nrow+1); ++i) {
	      mxt[i]=Mx[i+j*(nrow+1)];
	    }
	    /* Create an abridged life table from age-specific mortality rates
	     for columns Lx and lx alone; note the age format for abridged life table*/		 
	    LifeTableC(sex, nrow, mxt, Lx, lx);
	
	    /* cohort-period separation factors */
	    int ii;
		for(i=0; i<(nrow-1); ++i) {
			ii = i + 1;
			csf[i] = (Lx[i]-5.0*lx[ii])/(Lx[i] - Lx[ii]);	      			
		}

	    /* last age groups  */ 
	    csf[nrow-1] = 1.0; 
	    csf[nrow] = 0.0; 
	
		/***************************************************************************/
    	/* period deaths first age group*/
		i=0;
    	Deaths[j*nrow] = cdeaths[0] + cdeaths[1]*csf[0];

    	/* period deaths middle age groups and last, open ended age group */
		for(i=1; i<(nrow-1); ++i) {
			Deaths[i + j*nrow] = cdeaths[i] * (1-csf[i-1]) + cdeaths[i+1] * csf[i];
    	}
		Deaths[nrow-1 + j*nrow] = cdeaths[nrow-1] * (1-csf[nrow-2]);
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
 * Core population projection function TotalPopProj
 * Called from function StoPopProj in predict.pop.R
 * Parameter
 * int    *npred              number of prediction intervals (? time points)
 * double *MIGm, *MIGf        male, female migration by age
 * int    *migr, *migc        rows, columns of migration array
 * int    *MIGtype            type of migration adjustement
 * double *srm, *srf          male, female survivor ratios by age
 * double *asfr               age-specific fertility rates 
 * double *srb                sex ratio at birth
 * double *popm, *popf        male, female population by age
 * double *totp               total population by age 
 * double *btagem, *btagef    male, female births by age of mother
 * double *deathsm, *deathsf  male, female deaths by age
 ------------------------------------------------------------------------------
 TB: replaced expressions (j-1) with jve (jve <- j-1 )where appropriate        
     Note that more simplification is possible by replacing jve*adim 
     with another variable t_offset (<- jve*adim)
******************************************************************************/
void TotalPopProj(int *npred, double *MIGm, double *MIGf, int *migr, int *migc,
				  int *MIGtype, double *srm, double *srf, double *asfr, double *srb, 
				  double *mxm, double *mxf,
				  double *popm, double *popf, double *totp, 
				  double *btagem, double *btagef, double *deathsm, double *deathsf
					) {
	double migm[*migr+6][*migc], migf[*migr+6][*migc], totmigm[*migr+6][*migc], totmigf[*migr+6][*migc];
	double b, bt[7], bm, bf, mmult, srb_ratio;
	double Lxm[27], Lxf[27], lxm[27], lxf[27];
	double cdeathsm[27], cdeathsf[27];
	double mxtm[28], mxtf[28];
	/* forward and backward estimated deaths */
	double dfw, dbw;        
	/* cohort separation factor males, females*/
	double csfm[27],csff[27]; 
	int i, j, jve, adim, adim1, adimmx, nrow, ncol, n, t, t1, t_offset;
	int const male = 1;
	int const female = 2;


	int debug = 0;/* testing*/
	
	nrow = *migr;
	ncol = *migc;
	n = *npred;
	
	/* adim is an offset (representing time), which necessary to access and store two-dimensional data
	 organised by age and time/period into a one-dimensional vector requiered by accessing c from R
	 the expression
	*/

	adim = 27;         /* number of age groups of 5 year width up to age 130+ */
	adim1 = adim - 1;  /* Last index, number of age groups minus one           */
	adimmx = adim + 1; /* number of age groups of abridged life table input   */

	for(j=0; j<ncol; ++j) {
		for(i=0; i<nrow; ++i) {		
			migm[i][j] = MIGm[i + j*nrow];
			migf[i][j] = MIGf[i + j*nrow];
		}
		for(i=nrow; i<nrow+6; ++i) {
			migm[i][j] = 0;
			migf[i][j] = 0;
		}
		mmult = 1; /* warning if not initalized */
	    t = j*adim;
		switch (*MIGtype) {
			case 0: /* migration evenly distributed over each interval (MigCode=0) */
				for(i=1; i<nrow+6; ++i) {
					totmigm[i][j] = 0.5*(migm[i][j] + migm[i-1][j]*srm[i + t]);
					totmigf[i][j] = 0.5*(migf[i][j] + migf[i-1][j]*srf[i + t]);
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
		/* simplify notation for j and j-1 */
		t = j*adim;
		t1 = (j-1)*adim;
		t_offset = jve*adim;
		 
		/* Hana S.*/
		/* Time index (j) of survival ratio, migration and vital events is shifted by one in comparison to population,
			   i.e. pop[0] is the current period, whereas sr[0], mig[0] etc. is the first projection period.*/
		
		/* Compute ages >=5 */
		for(i=1; i<adim1; ++i) {

			popm[i + t] = popm[i-1 + t_offset] * srm[i + t_offset];
			popf[i + t] = popf[i-1 + t_offset] * srf[i + t_offset];
			totmigm[i][jve] = fmax(totmigm[i][jve], -1*popm[i + t]); /* assures population is not negative */
			popm[i + t] = popm[i + t] + totmigm[i][jve];
			totmigf[i][jve] = fmax(totmigf[i][jve], -1*popf[i + t]);
			popf[i + t] = popf[i + t] + totmigf[i][jve];
		}
		/* i = adim1 */
		popm[26 + t] = (popm[26 + t_offset] + popm[25 + t_offset]) * srm[26 + t_offset];
		popf[26 + t] = (popf[26 + t_offset] + popf[25 + t_offset]) * srf[26 + t_offset];
		totmigm[26][jve] = fmax(migm[26][jve], -1*popm[26 + t]);
		popm[26 + t] = popm[26 + t] + totmigm[26][jve];
		totmigf[26][jve] = fmax(migf[26][jve], -1*popf[26 + t]);
		popf[26 + t] = popf[26 + t] + totmigf[26][jve];
		
		/* calculating births, total, male, female */
		/* birth during 5-yrs */
		srb_ratio = srb[jve] / (1 + srb[jve]);
		for(i=3; i<10; ++i) {
			bt[i-3] = (popf[i + t_offset] + popf[i + t]) * asfr[i-3 + jve*7] * 0.5;
			btagem[i-3+jve*7] = bt[i-3] * srb_ratio;
			btagef[i-3+jve*7] = bt[i-3] - btagem[i-3+jve*7];
		}
		b = sum(bt, 7);
		bm = b * srb_ratio;
		bf = b - bm; /* avoids rounding errors, replaces bf = b / (1 + srb[jve]); */
		
		/* births surviving to age 0-4 */
		popm[t] = bm * srm[t_offset];
		popf[t] = bf * srf[t_offset];
		totmigm[0][jve] = fmax(mmult * migm[0][jve], -1*popm[t]);
		popm[t] = popm[t] + totmigm[0][jve];
		totmigf[0][jve] = fmax(mmult * migf[0][jve], -1*popf[t]);
		popf[t] = popf[t] + totmigf[0][jve];
		
		/* get total for all ages */
		for(i=0; i<adim; ++i) {
			totp[j] += popm[i + t]+popf[i + t];
		}
		
	    /**************************************************************************/
	    /* cohort deaths                                                          */
	    /* Calculated by a combination of forward-backward estimation (UN ABACUS) */
	    /* dfw deaths by forward projection    ABACUS: PART1                      */
	    /* dbw deaths by backward projection   ABACUS: C                          */
	    /* 1.  Deaths accuring to births                                           */
	    /* 2.  Deaths occuring to closed age groups (middle age groups)            */
	    /* 3.  Deaths accuring to last-opended age group                           */
	    /**************************************************************************/
	
	    /* deaths occuring to births */
	    i = 0;
	    /* males*/
	    dfw = bm * (1-srm[t_offset]);
	    dbw = popm[i + t] *((1-srm[i + t_offset])/srm[i + t_offset]);
	    cdeathsm[i] = 0.5* (dfw + dbw);
	
	    /* females*/
	    dfw = bf * (1-srf[t_offset]);
	    dbw = popf[i + t] *((1-srf[i + t_offset])/srf[i + t_offset]);
	    cdeathsf[i] = 0.5* (dfw + dbw);
	
	    /* more compact, but less readable */
	    /* cdeathsm[0] = 0.5* (bm * (1-srm[t_offset]) + (popm[1 + t] *(1-srm[t_offset])/srm[t_offset]));*/
	    /* cdeathsm[0] = 0.5* (bf * (1-srf[t_offset]) + (popf[1 + t] *(1-srf[t_offset])/srf[t_offset]));*/
	
	    /* closed age groups, index 1 to 25     */
	    for(i=1; i<adim1; ++i) {
	      /* males */
	      dfw = popm[i-1 + t_offset] * (1-srm[i + t_offset]);
	      dbw = popm[i + t] * (1-srm[i + t_offset])/srm[i + t_offset];
	      cdeathsm[i] = 0.5 * (dfw + dbw);
	
	      /* females */
	      dfw = popf[i-1 + t_offset] * (1-srf[i + t_offset]);
	      dbw = popf[i + t] * (1-srf[i + t_offset])/srf[i + t_offset];
	      cdeathsf[i] = 0.5 * (dfw + dbw);
	      
	    }
	    
	    /* last age group, index 26,  i = adim1;*/
	    /* males*/
	    dfw = (popm[i + t_offset]+popm[i-1 + t_offset])*(1-srm[i + t_offset]);
	    dbw = popm[i + t] * (1-srm[i + t_offset])/srm[i + t_offset];
	    cdeathsm[i] = 0.5 * (dfw + dbw);
	    
	    /* females*/
	    dfw = (popf[i + t_offset]+popf[i-1 + t_offset])*(1-srf[i + t_offset]);
	    dbw = popf[i + t] * (1-srf[i + t_offset])/srf[i + t_offset];
	    cdeathsf[i] = 0.5 * (dfw + dbw);
	
	    /**************************************************************************/
	    /* period deaths                                                          */
	    /* 1. calculated cohort separation factors                                */
	    /* 2. split cohort deaths into Lexis triangles                            */
	    /* 3. rearrange lexis triangles into period deaths (period-age format)    */
	    /**************************************************************************/
	    
	    /* Create an abridged life table from age-specific mortality rates
	     for columns Lx and lx alone; note the age format for abridged life table*/
	    for(i=0; i<adimmx; ++i) {
	      mxtm[i]=mxm[i + jve*adimmx];
	      mxtf[i]=mxf[i + jve*adimmx];
	    }
	    /* worked when calling with adim = 27 => last index of abridged life table */
	    /* Note that results have last index at adim1 = 26! */
	    LifeTableC(male, adim, mxtm, Lxm, lxm);
	    LifeTableC(female, adim, mxtf, Lxf, lxf);
	
	    /* cohort-period separation factors */
	    int ii;
	    for(i=0; i<(adim1-1); ++i) {
	      ii = i + 1;
	      csfm[i] = (Lxm[i]-5.0*lxm[ii])/(Lxm[i] - Lxm[ii]);
	      csff[i] = (Lxf[i]-5.0*lxf[ii])/(Lxf[i] - Lxf[ii]);
	    }
	    /* last age groups  */ 
	    csfm[adim1-1] = 1.0; 
	    csff[adim1-1] = 1.0;
	    csfm[adim1] = 0.0; 
	    csff[adim1] = 0.0;
	    
	    if((debug==1) && (j==1)){
	      Rprintf("\n csfm,adim1= %i", adim);
	      printArray(csfm,adim);
	    }

		/***************************************************************************/
    	/* period deaths first age group*/
    	deathsm[t_offset] = cdeathsm[0] + cdeathsm[1]*csfm[0];
    	deathsf[t_offset] = cdeathsf[0] + cdeathsf[1]*csff[0];
    	i=0;  

    	/* period deaths middle age groups and last, open ended age group */
    	for(i=1; i<adim1; ++i) {
      		deathsm[i + t_offset] = cdeathsm[i]*(1-csfm[i-1]) + cdeathsm[i+1] * csfm[i];
      		deathsf[i + t_offset] = cdeathsf[i]*(1-csff[i-1]) + cdeathsf[i+1] * csff[i];
    	}
	}	
}	
