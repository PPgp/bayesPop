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
 * Core population projection function CCM (cohort-component method) 
 * Called from function StoPopProj in predict.pop.R
 * This function should match the UN implementation in the ccmppWPP package,
 * function project_ccmpp_z_by_z
 * 
 * Parameters:
 * int    *nobserved          0 if the CCM is used for projections, or number of time points for which
 *                            population is observed and only vital events are to be computed.
 * int    *abridged           0 if this is a single year & single age CCM, 1 if it is 5x5 CCM
 * int    *npred              number of time points the CCM should run for (excluding present time)
 * double *MIGm, *MIGf        male, female migration by age for npred time points
 * int    *migr, *migc        rows, columns of migration array
 * int    *MIGtype            type of migration adjustment (0: evenly distributed, 9: at the end of time interval)
 * double *srm, *srf          male, female survivor ratios by age for npred time points
 * double *asfr               age-specific fertility rates for npred time points
 * double *srb                sex ratio at birth for npred time points
 * double *Lm, *Lf, *lxm, *lxf  male, female life table quantities Lx and lx for npred time points
 * int *nages                 number of population age groups
 * int *nfages                number of child-bearing age groups
 * int *fstart                index of the first child-bearing age group (R-like, i.e. assuming index starts with 1)
 * double *popm, *popf        male, female population by age for npred+1 time points
 *                              If *nobserved is 0, values for the first time point are passed as inputs, 
 *                              the rest is filled in by this function. 
 *                              Otherwise *nobserved values are used as input.
 * double *totp               total population
 *                              If *nobserved is 0, first value is passed as input, the rest is filled in here
 *                              If *nobserved > 0, it is not used.
 * double *btagem, *btagef    male, female births by age of mother (output)
 * double *deathsm, *deathsf  male, female period deaths by age (output)
 * ****************************************************************************
 */
 
void CCM(int *nobserved, int *abridged, int *npred, double *MIGm, double *MIGf, int *migr, int *migc,
                     int *MIGtype, double *MIGratem, double *MIGratef, int *MIGratecode, 
                     double *srm, double *srf, double *asfr, double *srb, 
                     double *Lm, double *Lf, double *lxm, double *lxf,
                     int *nages, int *nfages, int *fstart, 
                     double *popm, double *popf, double *totp, 
                     double *btagem, double *btagef, double *deathsm, double *deathsf,
                     double *finmigm, double *finmigf
    ) {
     
    double b, bm, bf, srb_ratio;
    int i, j, jve, adim, adim1, nrow, ncol, n, t, t1, t_offset;

    int debug = 0; /* for testing*/
    
    adim = *nages;         /* number of age groups */
    adim1 = adim - 1;      /* Number of age groups minus one */
    int adimdif = adim - *migr; /* Difference between number of ages in migration data and last age group */
    int adimfert = *nfages;     /* Number of reproductive age groups */
    int adimfert_start = *fstart - 1; /* Start index of reproductive age (from R to C) */
    
    double cdeathsm[adim], cdeathsf[adim], bt[adimfert];
    double migm[adim][*migc], migf[adim][*migc]; /*migration counts */
    double popadjm[adim][*npred+1], popadjf[adim][*npred+1], current_popf;
    double migendm[adim][*migc], migendf[adim][*migc], migmidm[adim][*migc], migmidf[adim][*migc];
    double migstartm[adim][*migc], migstartf[adim][*migc];
    double totmigcount, tpopm, tpopf, trmig, trmigpos, trmigneg, tmp;
    
    double csfm[adim],csff[adim]; /* cohort separation factor males, females*/
    
    const double minpop = 0.0005; /* minimum accepted population */
    const double max_out_rate = -0.8; /* the maximum portion of population that can leave when using migration rates */ 
    
    nrow = *migr;
    ncol = *migc;
    n = *npred;
    
    if(debug>=2)
        Rprintf("\nadim = %i adimfert = %i adimfert_start = %i migc = %i npred = %i adimdif = %i", 
                adim, adimfert, adimfert_start, *migc, n, adimdif);
    
    /* If migration is given as rates in MIGratem and MIGratef, then migm and migf
     * are assumed to contain the age-schedules.
     * 
     * MIGratecode:
     *      0: migration given as counts in MIGm and MIGf
     *      > 0: migration given as total rates in MIGratem and MIGratef,
     *              while MIGm and MIGf contains the age schedules.
     *      1: the age schedules are proportions
     *      2: the age schedules are totals (used for schedules that have positive as well as negative parts)
     *      3: the age schedules are the actual final age-specific rates (i.e. rate * schedule)
     */
    
    for(j=0; j<ncol; ++j) {
        /* If migration has less age groups then population, fill-in the remaining age groups with 0 */
        for(i=0; i<nrow; ++i) {		
            migm[i][j] = MIGm[i + j*nrow];
            migf[i][j] = MIGf[i + j*nrow];
        }
        for(i=nrow; i<nrow+adimdif; ++i) {
            migm[i][j] = 0;
            migf[i][j] = 0;
        }
        /* Migration is going to be added to population at the start (for computing deaths), 
         * middle (for computing births) and at the end of the interval.
         * Migration type (MIGtype) determines what is added when.
         * If migration is given as rates, counts are computed and added at the end of the interval.
         */
        if(MIGratecode[j] == 0){ /* migration is given as counts */
            switch (*MIGtype) { 
            case 0: /* migration evenly distributed over each interval (MigCode=0) */
                for(i=0; i<nrow+adimdif; ++i) { /* half at the start, half in the middle, nothing at the end*/
                    migstartm[i][j] = 0.5 * migm[i][j];
                    migstartf[i][j] = 0.5 * migf[i][j];
                    migmidm[i][j] = migstartm[i][j];
                    migmidf[i][j] = migstartf[i][j];
                    migendm[i][j] = 0;
                    migendf[i][j] = 0;
                }
                break;
            default: /* migration at the end of each interval (MigCode=9)*/
                for(i=0; i<nrow+adimdif; ++i) { /* everything is added at the end */
                    migstartm[i][j] = 0;
                    migstartf[i][j] = 0;
                    migmidm[i][j] = 0;
                    migmidf[i][j] = 0;
                    migendm[i][j] = migm[i][j];
                    migendf[i][j] = migf[i][j];
                }
                break;
            }
        } else { /* migration is given as rates */
            for(i=0; i<nrow+adimdif; ++i) { 
                migstartm[i][j] = 0;
                migstartf[i][j] = 0;
                migmidm[i][j] = 0;
                migmidf[i][j] = 0;
                migendm[i][j] = migm[i][j]; /* these are now age schedules, depending on MIGratecode */
                migendf[i][j] = migf[i][j];
            }
        }
    }

 

    /* Population projection for one trajectory */
    for(j=1; j<(n+1); ++j) {
        jve = j-1;
        if(debug==1) Rprintf("\nj = %i, jve = %i", j, jve);
        t = j*adim;
        t1 = (j-1)*adim; /* used to access the previous time period */
        t_offset = jve*adim; /* although the same as t1, using a different name for clarity, see comment below */
        
        /* Note that time index (j) of survival ratio, migration and vital events/rates is shifted by one in comparison to population,
         * i.e. pop[0] is the current period, whereas sr[0], mig[0] etc. is the first projection period.
         * Thus, j and t for pop vs. and jve and t_offset for vital rates/events refer to the same time period.
         */
        
        /* Adjust population by migrants at the start of time interval*/
        for(i=0; i < adim; ++i) {
            popadjm[i][j] = popm[i + t1] + migstartm[i][jve];
            popadjf[i][j] = popf[i + t1] + migstartf[i][jve];
        }
        
        /* Compute cohort deaths and population for closed ages >= 1 */
        for(i=1; i < adim1; ++i) {
            cdeathsm[i] = popadjm[i-1][j] * (1 - srm[i + t_offset]);
            cdeathsf[i] = popadjf[i-1][j] * (1 - srf[i + t_offset]);
            if(*nobserved <= j){ /* add migration at the middle of the interval */ 
                popm[i + t] = popadjm[i-1][j] - cdeathsm[i] + migmidm[i][jve];
                popf[i + t] = popadjf[i-1][j] - cdeathsf[i] + migmidf[i][jve];
            }
            if((debug>=2) && (j==1)){
                Rprintf("\ni = %i", i);
                Rprintf("\npopadjm[i-1][j]= %f, srm[i+t_offset]= %f, cdeathsm[i]= %f, migendm[i][jve]= %f, popm[i+t]= %f",
                        popadjm[i-1][j], srm[i + t_offset], cdeathsm[i], migendm[i][jve], popm[i + t]);
            }
        }
        /* open age group */
        i = adim1;
        cdeathsm[i] = (popadjm[i-1][j] + popadjm[i][j]) * (1 - srm[i + t_offset]);
        cdeathsf[i] = (popadjf[i-1][j] + popadjf[i][j]) * (1 - srf[i + t_offset]);

        if(*nobserved <= j){
            popm[i + t] = popadjm[i-1][j] + popadjm[i][j] - cdeathsm[i] + migmidm[i][jve];
            popf[i + t] = popadjf[i-1][j] + popadjf[i][j] - cdeathsf[i] + migmidf[i][jve];
        }
        /* calculating births (total & sex-specific)*/
        srb_ratio = srb[jve] / (1 + srb[jve]);
        for(i=adimfert_start; i<adimfert+adimfert_start; ++i) {
            current_popf = popf[i + t];
            if(*nobserved > j) current_popf = popf[i + t] - migendf[i][jve]; /* migration at the end of the interval should not count into births */
            bt[i-adimfert_start] = (popadjf[i][j] + current_popf) * asfr[i-adimfert_start + jve*adimfert] * 0.5;
            btagem[i-adimfert_start+jve*adimfert] = bt[i-adimfert_start] * srb_ratio;
            btagef[i-adimfert_start+jve*adimfert] = bt[i-adimfert_start] - btagem[i-adimfert_start+jve*adimfert];
            if((debug>=2) && (*nobserved > 1)){
                Rprintf("\ni = %i", i);
                Rprintf("\npopadjm[i-1][j]= %f, srm[i+t_offset]= %f, cdeathsm[i]= %f, migendm[i][jve]= %f, popm[i+t]= %f",
                        popadjm[i-1][j], srm[i + t_offset], cdeathsm[i], migendm[i][jve], popm[i + t]);
            }
        }
        b = sum(bt, adimfert); /* total over ages */
        bm = b * srb_ratio; /* male total */
        bf = b - bm; /* female total; avoids rounding errors, replaces bf = b / (1 + srb[jve]); */
        
        /* cohort deaths for the first age group */
        cdeathsm[0] = bm * (1-srm[t_offset]);
        cdeathsf[0] = bf * (1-srf[t_offset]);
        
        if(*nobserved <= j){
            /* births surviving to age 1 */
            popm[t] = bm - cdeathsm[0] + migmidm[0][jve];
            popf[t] = bf - cdeathsf[0] + migmidf[0][jve];

            /* If migration is a rate, compute the migration counts */
            if(MIGratecode[jve] > 0){
                trmig = 0;
                trmigpos = 0;
                trmigneg = 0;
                if(MIGratecode[jve] == 2){ /* for schedules that have positive and negative parts get the sum of these */
                    for(i=0; i < adim; ++i) {  
                        trmig = trmig + migendm[i][jve] + migendf[i][jve];
                        trmigpos = trmigpos + fmax(0, migendm[i][jve]) + fmax(0, migendf[i][jve]);
                        trmigneg = trmigneg + fmin(0, migendm[i][jve]) + fmin(0, migendf[i][jve]);
                    }
                }
                
                if(MIGratecode[jve] < 3){ /* need to disaggregate into ages */
                    /* sum population together */
                    tpopm = 0; tpopf = 0;
                    for(i=0; i < adim; ++i) {
                        tpopm = tpopm + popm[i + t];
                        tpopf = tpopf + popf[i + t];
                    }
                    /* total migration count; not more than 80% of total population can leave */ 
                    totmigcount = fmax(MIGratem[jve], max_out_rate) * tpopm + fmax(MIGratef[jve], max_out_rate) * tpopf; 
                    tmp = 0;
                    /* distribute into ages */
                    for(i=0; i < adim; ++i) { 
                        if(MIGratecode[jve] == 1){ /* multiply the total rate with age schedule */
                            migendm[i][jve] = fmax(fmax(0, popm[i+t])*max_out_rate, totmigcount*migendm[i][jve]); /* assures there is no depopulation */
                            migendf[i][jve] = fmax(fmax(0, popf[i+t])*max_out_rate, totmigcount*migendf[i][jve]);
                            tmp = tmp + migendm[i][jve];
                        }
                        if(MIGratecode[jve] == 2){ /* these schedules are actual counts, so we shift them up and down, relative to the sum */
                            if(totmigcount > 0){
                                if(migendm[i][jve] > 0) /* distribute the difference across the positive part of the schedule */
                                    migendm[i][jve] = migendm[i][jve] + (totmigcount - trmig)*migendm[i][jve]/trmigpos;
                                if(migendf[i][jve] > 0) 
                                    migendf[i][jve] = migendf[i][jve] + (totmigcount - trmig)*migendf[i][jve]/trmigpos;
                            } else { /* rate is negative */
                                if(migendm[i][jve] < 0) /* distribute the difference across the negative part of the schedule */
                                    migendm[i][jve] = migendm[i][jve] + (totmigcount - trmig)*migendm[i][jve]/trmigneg;
                                if(migendf[i][jve] < 0) 
                                    migendf[i][jve] = migendf[i][jve] + (totmigcount - trmig)*migendf[i][jve]/trmigneg;   
                            }
                        }
                    }
                    if(debug>=1 && jve > 9 && totmigcount > 0){
                        Rprintf("\ntotmigcount = %f, tpopm = %f, tpopf = %f, sum(migendm) = %f, migendm[5] = %f, migendm[10] = %f, migendf[10] = %f",  
                                totmigcount, tpopm, tpopf, tmp, migendm[5][jve], migendm[10][jve], migendf[10][jve]);
                    }
                } else { /* MIGratecode[jve] == 3; the schedules are the final age-specific rates */
                    for(i=0; i < adim; ++i) { 
                        migendm[i][jve] = fmax(fmax(0,popm[i+t])*max_out_rate, migendm[i][jve]);
                        migendf[i][jve] = fmax(fmax(0,popf[i+t])*max_out_rate, migendf[i][jve]);
                    }
                }
            }
            
            /* add migration at the end of the interval, adjust negative population and compute total pop */
            totp[j] = 0;
            for(i=0; i < adim; ++i) {
                migendm[i][jve] = fmax(migendm[i][jve], -popm[i+t] + minpop); /* adjust migration if it would yield negative population */
                popm[i+t] = popm[i + t] + migendm[i][jve];
                migendf[i][jve] = fmax(migendf[i][jve], -popf[i+t] + minpop);
                popf[i+t] = popf[i + t] + migendf[i][jve];
                totp[j] += popm[i + t] + popf[i + t]; 
            }
        }
        
        for(i=0; i < nrow; ++i) {
            finmigm[i + jve*nrow] = migstartm[i][jve] + migmidm[i][jve] + migendm[i][jve];
            finmigf[i + jve*nrow] = migstartf[i][jve] + migmidf[i][jve] + migendf[i][jve];
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
        if(*abridged == 1) fct = 5;
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
        
        if((debug > 2) && (j==1)){
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
