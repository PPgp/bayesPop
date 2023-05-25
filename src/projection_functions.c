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
                     int *MIGtype, int *MIGratecode, double *srm, double *srf, double *asfr, double *srb, 
                     double *Lm, double *Lf, double *lxm, double *lxf,
                     int *nages, int *nfages, int *fstart, 
                     double *popm, double *popf, double *totp, 
                     double *btagem, double *btagef, double *deathsm, double *deathsf,
                     double *finmigm, double *finmigf
    ) {
     
    double b, bm, bf, srb_ratio;
    int i, j, jve, adim, adim1, nrow, ncol, n, t, t1, t_offset;

    int debug = 1; /* for testing*/
    
    adim = *nages;         /* number of age groups */
    adim1 = adim - 1;      /* Number of age groups minus one */
    int adimdif = adim - *migr; /* Difference between number of ages in migration data and last age group */
    int adimfert = *nfages;     /* Number of reproductive age groups */
    int adimfert_start = *fstart - 1; /* Start index of reproductive age (from R to C) */
    
    double cdeathsm[adim], cdeathsf[adim], bt[adimfert];
    double migm[adim][*migc], migf[adim][*migc];
    double popadjm[adim][*npred+1], popadjf[adim][*npred+1], current_popf;
    double migendm[adim][*migc], migendf[adim][*migc], migmidm[adim][*migc], migmidf[adim][*migc];
    double migstartm[adim][*migc], migstartf[adim][*migc];
    double trmigstart, trmigmid, trmigend, tmigstart, tmigmid, tmigend, tpop, trmig;
    double mmigage_schedule[101], fmigage_schedule[101];
    
    double csfm[adim],csff[adim]; /* cohort separation factor males, females*/
    
    const double minpop = 0.0005; /* minimum accepted population */
    /* Rogers-Castro schedule in case migration is given as total rate */
    const double migrc_schedule[27] = {0.06133, 0.02667, 0.02067, 0.10467, 0.188, 
                                        0.18067, 0.13733, 0.09533, 0.064, 0.04267, 
                                        0.028, 0.01867, 0.012, 0.008, 0.00533, 
                                        0.00333, 0.00333, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    const double migrc_schedule1y[101] = {0.01604, 0.01362, 0.01157, 0.00983, 0.00837, 0.00712, 0.00606, 0.00517, 0.00442, 0.00382, 
                                          0.00346, 0.00349, 0.00418, 0.00582, 0.00856, 0.01236, 0.01695, 0.02192, 0.02681, 0.03121, 
                                          0.03485, 0.03757, 0.03933, 0.0402, 0.04027, 0.03968, 0.03858, 0.03709, 0.03534, 0.03341, 
                                          0.0314, 0.02936, 0.02734, 0.02538, 0.02349, 0.02169, 0.01999, 0.01841, 0.01692, 0.01555, 
                                          0.01427, 0.01309, 0.012, 0.01101, 0.01008, 0.00924, 0.00847, 0.00776, 0.00711, 0.00652, 
                                          0.00598, 0.00548, 0.00502, 0.00461, 0.00423, 0.00388, 0.00356, 0.00327, 0.00301, 0.00277, 
                                          0.00255, 0.00235, 0.00216, 0.00199, 0.00184, 0.0017, 0.00157, 0.00145, 0.00135, 0.00125, 
                                          0.00116, 0.00108, 0.001, 0.00093, 0.00087, 0.00082, 0.00076, 0.00071, 0.00067, 0.00063, 
                                          6e-04, 0.00056, 0.00053, 0.00051, 0.00048, 0.00046, 0.00044, 0.00041, 4e-04, 0.00038, 
                                          0.00037, 3e-05, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    nrow = *migr;
    ncol = *migc;
    n = *npred;
    
    if(debug>=2)
        Rprintf("\nadim = %i adimfert = %i adimfert_start = %i migc = %i npred = %i adimdif = %i", 
                adim, adimfert, adimfert_start, *migc, n, adimdif);
    
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
         * Migration type determines what is added when.
         */
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
    }
    if(*MIGratecode > 0 && *MIGratecode < 4){
        if(*abridged == 1) { /* 5-year Rogers Castro */ 
            for(i=0; i < adim; ++i) {
                mmigage_schedule[i] = 0.5 * migrc_schedule[i];
                fmigage_schedule[i] = 0.5 * migrc_schedule[i];
            }
        } else {              /* annual Rogers Castro */ 
            for(i=0; i < adim; ++i) {
                mmigage_schedule[i] = 0.5 * migrc_schedule1y[i];
                fmigage_schedule[i] = 0.5 * migrc_schedule1y[i];
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
        
        /* If migration is a rate, compute the migration counts */
        if(*MIGratecode > 0){
            if(*MIGratecode >= 3) { /* migration rates are disaggregated totals */
                trmigstart = 0;
                trmigmid = 0;
                trmigend = 0;
                trmig = 0;
                tpop = 0;
                for(i=0; i < adim; ++i) { /* first sum together */
                    trmigstart = trmigstart + migstartm[i][jve] + migstartf[i][jve];
                    trmigmid = trmigmid + migmidm[i][jve] + migmidf[i][jve];
                    trmigend = trmigend + migendm[i][jve] + migendf[i][jve];
                    trmig = trmig + migm[i][jve] + migf[i][jve];
                    tpop = tpop + popm[i + t1] + popf[i + t1];
                }
                tmigstart = trmigstart*tpop;
                tmigmid = trmigmid*tpop;
                tmigend = trmigend*tpop;
                if(*MIGratecode == 4) { /* use migration schedule passed in the mig objects */
                    for(i=0; i < adim; ++i) {
                        mmigage_schedule[i] = migm[i][jve]/trmig;
                        fmigage_schedule[i] = migf[i][jve]/trmig;
                    }
                }
                for(i=0; i < adim; ++i) { /* distribute into ages using either Rogers Castro or passed schedule */
                    migstartm[i][jve] = tmigstart*mmigage_schedule[i];
                    migstartf[i][jve] = tmigstart*fmigage_schedule[i];
                    migmidm[i][jve] = tmigmid*mmigage_schedule[i];
                    migmidf[i][jve] = tmigmid*fmigage_schedule[i];
                    migendm[i][jve] = tmigend*mmigage_schedule[i];
                    migendf[i][jve] = tmigend*fmigage_schedule[i];
                }
                if(debug>=1){
                    Rprintf("\ntrmig = %f, trmigend = %f",  trmig, trmigend);
                    Rprintf("\ntmigend = %f, tpop = %f, migendm[5] = %f, migendm[10] = %f",  tmigend, tpop, migendm[5][jve], migendm[10][jve]);
                }
            }
        }
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
