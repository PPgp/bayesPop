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
double * get_a05_cd(double mx0, int sex) {
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

double get_a0_ak_female(double mx0){
    if (mx0 < 0.01724) {return(0.14903 - 2.05527 * mx0);}
    if (mx0 < 0.06891) {return(0.04667 + 3.88089 * mx0);}
    return(0.31411);
}

double get_a0_ak_male(double mx0){
    if (mx0 < 0.0230) {return(0.14929 - 1.99545 * mx0);}
    if (mx0 < 0.08307) {return(0.02832 + 3.26021 * mx0);}
    return(0.29915);
}

/*****************************************************************************
 * Function returns the years lived by those who died in the age interval (ax)
 * for ages 0 and 1-4 (abridged life table)
 * based on mx(0,1) and sex
 * using the Andreev-Kingkade method (Demographic Research 2015)
 *****************************************************************************/
double * get_a05_ak(double mx0, int sex) {
    static double ax[2];
    double *ax_cd, sr, af, am;
    
    ax_cd = get_a05_cd(mx0, sex); /* for ages 1-4 */ 
    if(sex == 2) {/* female*/
        ax[0] = get_a0_ak_female(mx0);
    } else { 
        if(sex == 1) { /* male */
            ax[0] = get_a0_ak_male(mx0);
        } else { /* total */
            af = get_a0_ak_female(mx0);
            am = get_a0_ak_male(mx0);
            sr =   1.05 / (1 + 1.05);
            ax[0] = sr * am + (1-sr) * af;
        }
    }
    ax[1] = ax_cd[1];    /*4a1*/
    return(ax);
}

/*****************************************************************************
 * Function returns the years lived by those who died in the age interval (ax)
 * for ages 0 and 1-4 (abridged life table)
 * based on mx(0,1) and sex using 
 * either the Andreev-Kingkade (a0rule = 1) or Coale-Demeny (a0rule = 2) method
 *****************************************************************************/
double * get_a05(int a0rule, double mx0, int sex) {
    if(a0rule == 1) {
        return(get_a05_ak(mx0, sex));
    } else {
        return(get_a05_cd(mx0, sex));
    }
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
void doLifeTable(int sex, int nage, double *mx, int a0rule,
				double *Lx, double *lx, double *qx, double *ax) {
	
	int i;
	double k;     /* correcting factor in Greville approximation */
	double *tmpa; /* pointer to estimated ax[0] and ax[1] values */
	int nage1;
	nage1 = nage -1;

	tmpa = get_a05(a0rule, mx[0], sex);
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
void LifeTableC(int sex, int nage, double *mx, int a0rule,
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
  doLifeTable(sex, nage, mx, a0rule, Lx, lx, qx, ax);
  /* collapse 1L0 and 4L1 into 5L0 */
  Lxx[0] = Lx[0] + Lx[1];
  lxx[0] = lx[0];
  for(i = 1; i < nage; ++i) {
    Lxx[i] = Lx[i+1];
    lxx[i] = lx[i+1];
  }
}


/*****************************************************************************
 * Function calculates a life table for one-year age groups 
 * from age-specific mortality 
 *****************************************************************************/
void doLifeTable1y(int sex, int nage, double *mx, int a0rule, 
                   double *Lx, double *lx, double *qx, double *ax) {
    
    int i;
    //    double k;     /* correcting factor in Greville approximation */
    double *tmpa; /* pointer to estimated ax[0] and ax[1] values */
    int nage1;
    nage1 = nage -1;
    tmpa = get_a05(a0rule, mx[0], sex); 
    ax[0] = tmpa[0]; /* use only ax[0] */
    for(i = 1; i < nage; ++i) {
        /* k = 0.5 * log(fmax(mx[i+1] / fmax(mx[i-1], DBL_MIN), DBL_MIN));*/
        /*ax[i] = 0.5 - (1 / 12.0) * (mx[i] - k);*/
        ax[i] = 0.5;
    }
    /* penultimate ax calculated with k from previous age group */
    /*ax[nage1] = 0.5 - (1 / 12.0) * (mx[nage1] - k);*/

    /* correcting out-of (reasonable) bounds ax for older ages             */ 
    /* 0.42=1-exp(-1)/(1-exp(-1)), for constant mu=1, Kannisto assumption*/
    /*for(i = 10; i < nage; i++) {
        if(ax[i] < 0.97) {
        ax[i] = 0.97;
        }
    }*/
    lx[0] = 1;       /* l0 */
    /* calculate life table variables from mx and ax */
    for(i = 0; i < nage; ++i) {
        qx[i] = mx[i] / (1 + (1 - ax[i]) * mx[i]);
        lx[i+1] = fmax(lx[i] * (1-qx[i]), DBL_MIN);
        Lx[i] = lx[i+1] + ax[i] * (lx[i] - lx[i+1]);
    }
    /* Open ended age interval */
    Lx[nage] = lx[nage] / fmax(mx[nage], DBL_MIN); 
    qx[nage] = 1.0;
    ax[nage] = Lx[nage];
}

/* Wrapper around doLifeTable1y 
 * Used when qx and ax is not needed.
 */
void LifeTable1yC(int sex, int nage, double *mx, int a0rule,
                  double *Lx, double *lx) {
    double qx[nage], ax[nage];
    doLifeTable1y(sex, nage-1, mx, a0rule, Lx, lx, qx, ax);
}

/* Compute period deaths from observed population and survivor data */
void get_deaths_from_sr_1x1(double *Sr, int *N, double *Pop, 
                            int *MIGtype, double *Mig, double *Births, double *Lx, double *lx, 
                            double *Deaths) {
    int i, j, nrow, n, t;

    n = *N; /* periods */
	nrow = 101;
			
	double cdeaths[nrow];
	double migend[nrow][n], popadj[nrow][n];
			
	/* cohort separation factor males, females*/
	double csf[nrow]; 
			
	for(j=0; j<n; ++j) {
	    switch (*MIGtype) {
		case 0: /* migration evenly distributed over each interval (MigCode=0) */
		    for(i=0; i<nrow; ++i) {
			    migend[i][j] = 0.5 * Mig[i + j*nrow];
			}
			break;
		default: /* migration at the end of each interval (MigCode=9)*/
			for(i=0; i<nrow; ++i) {
			    migend[i][j] = 0;
			 }
		    break;
	    }
	}
	/* loop by time period */	
	for(j=0; j < n; ++j) {
	    /**************************************************************************/
	    /* cohort deaths                                                          */
	    /* 1.  Deaths accuring to births                                          */
	    /* 2.  Deaths occuring to closed age groups (middle age groups)           */
	    /* 3.  Deaths accuring to last-opended age group                          */
	    /**************************************************************************/
        t = j*nrow;
	    /* Adjust population by migrants */
	    for(i=0; i < nrow; ++i) {
	        popadj[i][j] = Pop[i + t] + migend[i][j];
	    }
	    cdeaths[0] = Births[j] * (1-Sr[t]);
	    for(i=1; i<(nrow-1); ++i) {
	        cdeaths[i] = popadj[i-1][j] * (1 - Sr[i + t]);
	    }
	    /* Last open-ended age category */
	    i = nrow - 1;
	    cdeaths[i] = (popadj[i-1][j] + popadj[i][j]) * (1 - Sr[i + t]);

	    /**************************************************************************/
	    /* period deaths                                                          */
		/* 1. calculated cohort separation factors                                */
		/* 2. split cohort deaths into Lexis triangles                            */
	    /* 3. rearrange lexis triangles into period deaths (period-age format)    */
		/**************************************************************************/
			    
	    /* cohort-period separation factors */
	    int ii;
		int fct = 1; /* for abridged LT it should be 5 */
		for(i=0; i<(nrow - 2); ++i) {
		    ii = i + 1;
			csf[i] = (Lx[i + t]-fct*lx[ii +  t])/(Lx[i + t] - Lx[ii + t]);
		}
		/* last two age groups  */ 
	    i = nrow - 2;
		csf[i] = (Lx[i + t]-fct*lx[i + 1 +  t])/Lx[i + t]; 
		csf[nrow-1] = 1.0; 
		/*Rprintf("\nj = %i", j);*/
        /* period deaths first age group*/
	    Deaths[t] = cdeaths[0] + cdeaths[1]*csf[0];

		/* period deaths middle age groups and last, open ended age group */
		for(i=1; i<nrow-1; ++i) {
		    Deaths[i + t] = cdeaths[i]*(1-csf[i-1]) + cdeaths[i+1] * csf[i];
		    /*Rprintf("\n i= %i, Deaths=%f, csf=%f, cdeaths=%f", i,Deaths[i + t], csf[i], cdeaths[i]);*/
		}
		Deaths[nrow-1 + t] = cdeaths[nrow-1] * (1-csf[nrow-2]);
		
	}
}

void get_deaths_from_sr_1x1_touse(int *Sex, double *Sr, int *N, double *Pop, 
                            double *Births, double *Deaths, double *Mx) {
    int i, j, nrow, n, sex;
    int a0rule = 1;
    
    n = *N; /* periods */
		sex=*Sex;
		nrow = 101;
		
		double Lx[nrow], lx[nrow], cdeaths[nrow], mxt[nrow];
		
		/* forward and backward estimated deaths */
		double dfw, dbw;        
		/* cohort separation factor males, females*/
		double csf[nrow]; 
		
		/* loop by time period */	
		for(j=0; j < n; ++j) {
		    /**************************************************************************/
		    /* cohort deaths                                                          */
		    /* 1.  Deaths accuring to births                                          */
		    /* 2.  Deaths occuring to closed age groups (middle age groups)           */
		    /* 3.  Deaths accuring to last-opended age group                          */
		    /**************************************************************************/
		    cdeaths[0] = Births[j] * (1-Sr[j*nrow]);
		    for(i=1; i<(nrow-1); ++i) {
		        dfw = Pop[i-1 + j*nrow] * (1-Sr[i + j*nrow]);
		        dbw = Pop[i + (j+1)*nrow] * ((1-Sr[i + j*nrow]) / Sr[i + j*nrow]);
		        cdeaths[i] = 0.5 * (dfw + dbw);
		    }
		    /* Last open-ended age category */
		    i = nrow - 1;
		    dfw = Pop[i + j*nrow] + Pop[i - 1 + j*nrow] * (1 - Sr[i + j*nrow]);
		    dbw = Pop[i + (j+1)*nrow] * ((1-Sr[i + j*nrow]) / Sr[i  + j*nrow]);
		    cdeaths[i] = 0.5 * (dfw + dbw);
		    /*Rprintf("\nj = %i: Pop[nrow-1]=%f Pop[nrow-1, j+1]=%f Sr[nrow-1]=%f dfw=%f dbw=%f ", j, 
		     Pop[nrow -1 + j*nrow], Pop[nrow -1 + (j+1)*nrow], Sr[nrow - 1 + j*nrow], dfw, dbw);*/
		    
		    
		    /**************************************************************************/
		    /* period deaths                                                          */
		    /* 1. calculated cohort separation factors                                */
		    /* 2. split cohort deaths into Lexis triangles                            */
		    /* 3. rearrange lexis triangles into period deaths (period-age format)    */
		    /**************************************************************************/
		    
		    for(i=0; i<nrow; ++i) {
		        mxt[i]=Mx[i+j*(nrow)];
		    }
		    /* Create an abridged life table from age-specific mortality rates
		     for columns Lx and lx alone; note the age format for abridged life table*/
		    
		    LifeTable1yC(sex, nrow, mxt, a0rule, Lx, lx);
		    
		    /* cohort-period separation factors */
		    int ii;
		    for(i=0; i<(nrow-2); ++i) {
		        ii = i + 1;
		        csf[i] = (Lx[i] - lx[ii])/(Lx[i] - Lx[ii]);  			
		    }
		    /* last age groups  */ 
		    csf[nrow-2] = 1.0; 
		    csf[nrow-1] = 0.0; 
		    
		    /***************************************************************************/
		    /* period deaths first age group*/
		    i=0;
		    Deaths[j*nrow] = cdeaths[0] + cdeaths[1]*csf[0];
		    /*Deaths[j*nrow] = cdeaths[0];*/
		    
		    /* period deaths middle age groups and last, open ended age group */
		    for(i=1; i<(nrow-1); ++i) {
		        Deaths[i + j*nrow] = cdeaths[i] * (1-csf[i-1]) + cdeaths[i+1] * csf[i];
		        /*Deaths[i + j*nrow] = cdeaths[i];*/
		        
		        /*Rprintf("\n i= %i, Deaths=%f, csf=%f, cdeaths=%f", i,Deaths[i + j*nrow], csf[i], cdeaths[i]);*/
		    }
		    Deaths[nrow-1 + j*nrow] = cdeaths[nrow-1] * (1-csf[nrow-1]);
		    /*Deaths[nrow-1 + j*nrow] = cdeaths[nrow-1];*/
		    /*Rprintf("\n i= %i, Deaths=%f, csf=%f, cdeaths=%f", i,Deaths[i + j*nrow], csf[i], cdeaths[i]);*/
		}
}


void get_deaths_from_sr_abridged(int *Sex, double *Sr, int *N, double *Pop, 
                        double *Births, double *Deaths, double *Mx) {
    double Lx[21], lx[21];
    double cdeaths[21];
    double mxt[22];
    int a0rule = 1;
    
    /* forward and backward estimated deaths */
    double dfw, dbw;        
    /* cohort separation factor males, females*/
    double csf[21]; 
    
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
		cdeaths[0] = 0.5* (dfw + dbw);
		for(i=1; i< nrow - 1; ++i) {
            /* age >= 5 */
			dfw = Pop[i-1 + j*nrow] * (1-Sr[i + j*nrow]);
			dbw = Pop[i + (j+1)*nrow] * ((1-Sr[i + j*nrow]) / Sr[i + j*nrow]);
			cdeaths[i] = 0.5 * (dfw + dbw);
		}
		/* Last open-ended age category */
        i = nrow - 1;
		dfw = Pop[i + j*nrow] + Pop[i - 1 + j*nrow] * (1-Sr[i + j*nrow]);	
		dbw = Pop[i + (j+1)*nrow] * ((1-Sr[i + j*nrow]) / Sr[i + j*nrow]);
		cdeaths[nrow - 1] = 0.5 * (dfw + dbw);
		/*Rprintf("\nPop[nrow]=%i Sr[nrow]=%f", Pop[nrow + j*nrow], Sr[nrow + j*nrow]);*/

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
	    LifeTableC(sex, nrow, mxt, a0rule, Lx, lx);
	
	    /* cohort-period separation factors */
	    int ii;
		for(i=0; i<(nrow-1); ++i) {
			ii = i + 1;
			csf[i] = (Lx[i]-5.0*lx[ii])/(Lx[i] - Lx[ii]);	      			
		}

	    /* last age groups  */ 
	    csf[nrow-2] = 1.0; 
	    csf[nrow-1] = 0.0; 
	
		/***************************************************************************/
    	/* period deaths first age group*/
		i=0;
    	Deaths[j*nrow] = cdeaths[0] + cdeaths[1]*csf[0];

    	/* period deaths middle age groups and last, open ended age group */
		for(i=1; i<(nrow-1); ++i) {
			Deaths[i + j*nrow] = cdeaths[i] * (1-csf[i-1]) + cdeaths[i+1] * csf[i];
    	}
		/*Deaths[nrow-1 + j*nrow] = cdeaths[nrow-1] * (1-csf[nrow-2]);*/
		Deaths[nrow-1 + j*nrow] = cdeaths[nrow-1];
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
	int a0rule = 1;
	
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
	    LifeTableC(male, adim, mxtm, a0rule, Lxm, lxm);
	    LifeTableC(female, adim, mxtf, a0rule, Lxf, lxf);
	
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


/*****************************************************************************
 * Core population projection function TotalPopProj1x1 
 * for 1-year age groups and 1-year time periods
 * Called from function StoPopProj in predict.pop.R
 * This function should match the UN implementation in the ccmppWPP package,
 * function project_ccmpp_z_by_z
 * 
 * Parameters:
 * int    *npred              number of prediction time points (excluding present time)
 * double *MIGm, *MIGf        male, female migration by age
 * int    *migr, *migc        rows, columns of migration array
 * int    *MIGtype            type of migration adjustment (0: evenly distributed, 9: at the end of time interval)
 * double *srm, *srf          male, female survivor ratios by age
 * double *asfr               age-specific fertility rates 
 * double *srb                sex ratio at birth
 * double *Lm, *Lf, *lxm, *lxf  male, female life table quantities Lx and lx
 * double *popm, *popf        male, female population by age:
 *                              - values for the first time point are passed as inputs, the rest is filled in by this function
 * double *totp               total population (first value is passed as input, the rest is filled in here)
 * double *btagem, *btagef    male, female births by age of mother (output)
 * double *deathsm, *deathsf  male, female period deaths by age (output)
 * ****************************************************************************
 */
 
void TotalPopProj1x1(int *npred, double *MIGm, double *MIGf, int *migr, int *migc,
                     int *MIGtype, double *srm, double *srf, double *asfr, double *srb, 
                     double *Lm, double *Lf, double *lxm, double *lxf,
                     double *popm, double *popf, double *totp, 
                     double *btagem, double *btagef, double *deathsm, double *deathsf
) {
    
    double b, bm, bf, srb_ratio;
    int i, j, jve, adim, adim1, nrow, ncol, n, t, t1, t_offset;

    int debug = 0; /* for testing*/
    
    adim = 131;         /* number of age groups up to age 130+ */
    adim1 = adim - 1;   /* Number of age groups minus one      */
    int adimdif = adim - *migr; /* Difference between number of ages in migration data and age 130 */
    int adimfert = 45;     /* Number of reproductive age groups */
    int adimfert_start = 10; /* Start index of reproductive age */
    
    /*double Lxm[adim], Lxf[adim], lxm[adim], lxf[adim], */
    double cdeathsm[adim], cdeathsf[adim], mxtm[adim], mxtf[adim], bt[adimfert];
    double migm[adim][*migc], migf[adim][*migc];
    double popadjm[adim][*migc], popadjf[adim][*migc];
    double migendm[adim][*migc], migendf[adim][*migc];
    
    /* cohort separation factor males, females*/
    double csfm[adim],csff[adim]; 
    
    nrow = *migr;
    ncol = *migc;
    n = *npred;
    
    for(j=0; j<ncol; ++j) {
        /* Fill-in migration up to 130 */
        for(i=0; i<nrow; ++i) {		
            migm[i][j] = MIGm[i + j*nrow];
            migf[i][j] = MIGf[i + j*nrow];
        }
        for(i=nrow; i<nrow+adimdif; ++i) {
            migm[i][j] = 0;
            migf[i][j] = 0;
        }
        switch (*MIGtype) {
        case 0: /* migration evenly distributed over each interval (MigCode=0) */
            for(i=0; i<nrow+adimdif; ++i) {
                migendm[i][j] = 0.5 * migm[i][j];
                migendf[i][j] = 0.5 * migf[i][j];
            }
            break;
        default: /* migration at the end of each interval (MigCode=9)*/
            for(i=0; i<nrow+adimdif; ++i) {
                migendm[i][j] = 0;
                migendf[i][j] = 0;
            }
            break;
        }
    }
    
    /* Population projection for one trajectory */
    for(j=1; j<(n+1); ++j) {
        jve = j-1;
        t = j*adim;
        t1 = (j-1)*adim; /* used to access the previous time period */
        t_offset = jve*adim; /* although the same as t1, using a different name for clarity, see comment below */
        
        /* Note that time index (j) of survival ratio, migration and vital events/rates is shifted by one in comparison to population,
         * i.e. pop[0] is the current period, whereas sr[0], mig[0] etc. is the first projection period.
         * Thus, j and t for pop vs. and jve and t_offset for vital rates/events refer to the same time period.
         */
        
        /* Adjust population by migrants */
        for(i=0; i < adim; ++i) {
            popadjm[i][j] = popm[i + t1] + migendm[i][jve];
            popadjf[i][j] = popf[i + t1] + migendf[i][jve];
        }
        
        /* Compute cohort deaths and population for closed ages >= 1 */
        for(i=1; i < adim1; ++i) {
            cdeathsm[i] = popadjm[i-1][j] * (1 - srm[i + t_offset]);
            cdeathsf[i] = popadjf[i-1][j] * (1 - srf[i + t_offset]);
            popm[i + t] = popadjm[i-1][j] - cdeathsm[i] + migendm[i][jve];
            popf[i + t] = popadjf[i-1][j] - cdeathsf[i] + migendf[i][jve];
            if((debug==2) && (j==1)){
                Rprintf("\ni = %i", i);
                Rprintf("\npopadjm[i-1][j]= %f, srm[i+t_offset]= %f, cdeathsm[i]= %f, migendm[i][jve]= %f, popm[i+t]= %f",
                        popadjm[i-1][j], srm[i + t_offset], cdeathsm[i], migendm[i][jve], popm[i + t]);
            }
        }
        /* open age group */
        i = adim1;
        cdeathsm[i] = (popadjm[i-1][j] + popadjm[i][j]) * (1 - srm[i + t_offset]);
        cdeathsf[i] = (popadjf[i-1][j] + popadjf[i][j]) * (1 - srf[i + t_offset]);
        popm[i + t] = popadjm[i-1][j] + popadjm[i][j] - cdeathsm[i] + migendm[i][jve];
        popf[i + t] = popadjf[i-1][j] + popadjf[i][j] - cdeathsf[i] + migendf[i][jve];
        
        /* calculating births (total & sex-specific)*/
        srb_ratio = srb[jve] / (1 + srb[jve]);
        for(i=adimfert_start; i<adimfert+adimfert_start; ++i) {
            bt[i-adimfert_start] = (popadjf[i][j] + popf[i + t]) * asfr[i-adimfert_start + jve*adimfert] * 0.5;
            btagem[i-adimfert_start+jve*adimfert] = bt[i-adimfert_start] * srb_ratio;
            btagef[i-adimfert_start+jve*adimfert] = bt[i-adimfert_start] - btagem[i-adimfert_start+jve*adimfert];
        }
        b = sum(bt, adimfert); /* total over ages */
        bm = b * srb_ratio; /* male total */
        bf = b - bm; /* female total; avoids rounding errors, replaces bf = b / (1 + srb[jve]); */
        
        /* cohort deaths for the first age group */
        cdeathsm[0] = bm * (1-srm[t_offset]);
        cdeathsf[0] = bf * (1-srf[t_offset]);
        
        /* births surviving to age 1 */
        popm[t] = bm - cdeathsm[0] + migendm[0][j];
        popf[t] = bf - cdeathsf[0] + migendf[0][j];

        /* add migration if needed, adjust negative population and compute total pop */
        totp[j] = 0;
        for(i=0; i < adim; ++i) {
            if(*MIGtype != 0) {
                popm[i+t] = popm[i + t] + migm[i][jve];
                popf[i+t] = popf[i + t] + migf[i][jve];
            }
            popm[i+t] = fmax(popm[i + t], 0.0005);
            popf[i+t] = fmax(popf[i + t], 0.0005);
            totp[j] += popm[i + t]+popf[i + t]; 
        }
        
        /**************************************************************************/
        /* period deaths                                                          */
        /* 1. calculated cohort separation factors                                */
        /* 2. split cohort deaths into Lexis triangles                            */
        /* 3. rearrange lexis triangles into period deaths (period-age format)    */
        /**************************************************************************/
        
        /* cohort-period separation factors */
        int ii;
        int fct = 1; /* for abridged LT it should be 5 */
        for(i=0; i<(adim1-1); ++i) {
            ii = i + 1;
            csfm[i] = (Lm[i + t_offset]-fct*lxm[ii +  t_offset])/(Lm[i + t_offset] - Lm[ii + t_offset]);
            csff[i] = (Lf[i + t_offset]-fct*lxf[ii + t_offset])/(Lf[i + t_offset] - Lf[ii + t_offset]);
        }
        /* last two age groups  */ 
        i = adim1-1;
        csfm[i] = (Lm[i + t_offset]-fct*lxm[i + 1 +  t_offset])/Lm[i + t_offset]; 
        csff[adim1-1] = (Lf[i + t_offset]-fct*lxf[i + 1 +  t_offset])/Lf[i + t_offset];
        csfm[adim1] = 1.0; 
        csff[adim1] = 1.0;
        
        if((debug==1) && (j==1)){
            Rprintf("\n csfm,adim1= %i", adim);
            printArray(csfm,adim);
        }
        
        /***************************************************************************/
        /* period deaths first age group*/
        deathsm[t_offset] = cdeathsm[0] + cdeathsm[1]*csfm[0];
        deathsf[t_offset] = cdeathsf[0] + cdeathsf[1]*csff[0];
        
        /* period deaths middle age groups and last, open ended age group */
        for(i=1; i<adim1; ++i) {
            deathsm[i + t_offset] = cdeathsm[i]*(1-csfm[i-1]) + cdeathsm[i+1] * csfm[i];
            deathsf[i + t_offset] = cdeathsf[i]*(1-csff[i-1]) + cdeathsf[i+1] * csff[i];
        }
        deathsm[adim1 + t_offset] = cdeathsm[adim1] * (1-csfm[adim1-1]); 
        deathsf[adim1 + t_offset] = cdeathsf[adim1] * (1-csff[adim1-1]);
    }
}
