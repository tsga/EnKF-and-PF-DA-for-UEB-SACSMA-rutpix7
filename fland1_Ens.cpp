/* fland1.f -- translated by f2c (version 20100827).

//4.15.18 revised to avoid having to link with libf2c

You must link the resulting object file with libf2c:
on Microsoft Windows system, link with libf2c.lib;
on Linux or Unix systems, link with .../path/to/libf2c.a -lm
or, if you install libf2c.a in a standard place, with -lf2c -lm
-- in that order, at the end of the command line, as in
cc *.o -lf2c -lm
Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

http://www.netlib.org/f2c/libf2c.zip
*/

//#include "f2c.h"
#include "uebDAsrPFuebpgdecls.h"
#include<cmath>

__host__ __device__
void fland1Device(float pxv, float dt, float aesc,
	float sacst[6], float sacpar[17],
	float &edmnd, float &surf, float &grnd, float &tet,
	float &sif, float &bfp, float &bfs, float &ssur, float &sdro, float &sperc,
	float hrapx, float hrapy, int &error)
{
	
/* System generated locals */
	int i__1;
	float r__1;
	double d__1, d__2;        //==double d__1, d__2;

	/* Local variables */
	int i__;

    float e1, e2, e3, e4, e5, bf, xx, xx1, del, red, sbf, tbf, tci, sfh, hpl, duz, sur, twx, bfcc, dinc, defr;
    int ninc;
    float pinc, spbf, perc;

    float dlzp, dlzs, bfncc, adimc, check, parea, addro, adimp, fracp, percf, eused, percm, percp, lzfpc, adsur, ratio,
	     lzfsc, lzfsh, lzfph, lzfsm, lzfpm, uzfwc, ratlz, roimp, perct, lztwc, uzfwh, lztwh, uzfwm, uztwc, lztwm, 
	    uzrat, ratlp, ratls, uztwh, percs, uztwm, excess, smcdry, duztwc, ratlzt, simpvt;

	/*static float e1, e2, e3, e4, e5, bf, xx, xx1, del, red, sbf, tbf, tci, sfh, hpl, duz, sur, twx, bfcc, dinc, defr;
	static int ninc;
	static float pinc, spbf, perc;
	static float dlzp, dlzs, bfncc, adimc, check, parea, addro, adimp, fracp, percf, eused, percm, percp, lzfpc, adsur, ratio,
		lzfsc, lzfsh, lzfph, lzfsm, lzfpm, uzfwc, ratlz, roimp, perct, lztwc, uzfwh, lztwh, uzfwm, uztwc, lztwm,
		uzrat, ratlp, ratls, uztwh, percs, uztwm, excess, smcdry, duztwc, ratlzt, simpvt;*/

/* vk  1/2008  Introduced option to generate normalized SM */
/* vk  1/2008  To select this option, change the next parameter from 0 to 1 */
/* vk  1/2008 */
/* vk Normalized soil moisture content */
/* c      PARAMETER (NORMALIZE = 1) */
/* vk Not normalized soil moisture content */
/* c      PARAMETER (NORMALIZE = 0) */
/*  DT is in days here */
/*  DTFRZ IN SEC., IDTFRZ IS # FRZ_STEPS */
/* ....................................... */
/*     THIS SUBROUTINE EXECUTES THE 'SAC-SMA ' OPERATION FOR ONE TIME */
/*         PERIOD. */
/* ....................................... */
/*     SUBROUTINE INITIALLY WRITTEN BY. . . */
/*            ERIC ANDERSON - HRL     APRIL 1979     VERSION 1 */
/* ....................................... */
/* VK  FROZEN GROUND CHANGES */
/* VK  UZTWC,UZFWC,LZTWC,LZFSC,LZFPC ARE TOTAL WATER STORAGES */
/* VK  UZTWH,UZFWH,LZTWH,LZFSH,LZFPH ARE UNFROZEN WATER STORAGES */
/* SACPAR() is array of original SAC parameters, and FRZPAR() is array */
/* of frozen ground parameters and calculated constants */
/* SACST() and FRZST() same for states */
/*  delited float FGCO(6),ZSOIL(6),TSOIL(8),FGPM(11) */
/* VK---------------------------------------------------------------- */
/* VK_02  NEW COMMON STATEMENT FOR DESIRED SOIL LAYERS */
/* VK     THIS VERSION HAS HARD CODED OUTPUT SOIL LAYERS */
/* VK     LATER ON IT SHOULD BE CHANGED TO MAKE THEM VARIABLE */
/*      INTEGER NINT/5/,NINTW/5/ */
/* k      REAL DSINT(10)/0.075,0.15,0.35,0.75,1.5,0.0,0.0,0.,0.,0./ */
/*      REAL DSINT(10)/0.10,0.40,0.6,0.75,1.5,0.0,0.0,0.,0.,0./ */
/* k      REAL DSINTW(10)/0.075,0.15,0.35,0.75,1.5,0.0,0.0,0.,0.,0./ */
/*      REAL DSINTW(10)/0.10,0.40,0.6,0.75,1.5,0.0,0.0,0.,0.,0./ */
/*      SAVE DSINT,DSINTW,NINT,NINTW */
/* VK---------------------------------------------------------------- */

/* The following variables need to be added to the output */
/* variables */
/*   Variable name:             Name in source code:    Units: */
/*  --------------------------------------------------------------------- */
/*   Interflow                        SIF                 mm/dt */
/*   PrimaryFlow                      BFP                 mm/dt */
/*   SupplementalFlow                 BFS                 mm/dt */
/*   ExcessFlow                       SSUR                mm/dt */
/*   DirectFlow                       SDRO                mm/dt */
/*   Percolation                      SPERC               mm/dt */

/* c      DIMENSION EPDIST(24) */
/*     COMMON BLOCKS */
/* c      COMMON/FSMPM1/UZTWM,UZFWM,UZK,PCTIM,ADIMP,RIVA,ZPERC,REXP,LZTWM, */
/* c     1              LZFSM,LZFPM,LZSK,LZPK,PFREE,SIDE,SAVED,PAREA */
/* VK */
/* VK      COMMON/FPMFG1/FGPM(10) */
/* --      COMMON/FPMFG1/itta,FGPM(15),ivers,ifrze */
/* VK_02  NEW COMMON BLOCK FOR INTERPOLATED SOIL TEMP AND SOIL PARAMETERS */
/* --      COMMON/TSLINT/TSINT(10),NINT,SWINT(10),SWHINT(10),NINTW */
/* --      COMMON/FRDSTFG/SMAX,PSISAT,BRT,SWLT,QUARTZ,STYPE,NUPL,NSAC, */
/* --     +               RTUZ,RTLZ,DZUP,DZLOW */
/* c      COMMON/FRZCNST/ FRST_FACT,CKSOIL,ZBOT */
/* --      COMMON/FRZCNST/ FRST_FACT,ZBOT */
/* VK */
/* VK  NEW FG VERSION PARAMETERS & SOIL LAYER DEFINITION: */
/* VK          FGPM(1) - SOIL TEXTURE CLASS */
/* VK          FGPM(2) - OPTIONAL, SOIL TEMPERATURE AT THE 3M DEPTH */
/* VK          FGPM(3) - OPTIONAL, POROSITY OF RESIDUE LAYER */
/* VK          PAR(18) [if no calb=FGPM(4)] - RUNOFF REDUCTION PARAMETER 1 */
/* VK          PAR(19) [if no calb=FGPM(5)] - RUNOFF REDUCTION PARAMETER 2 */
/* VK          PAR(20) [if no calb=FGPM(6)] - RUNOFF REDUCTION PARAMETER 3 */
/* VK          FGPM(6) - RUNOFF REDUCTION PARAMETER 3 (FOR ERIC'S VERSION ONLY) */
/* VK          FGPM(7) - NUMBER OF SOIL LAYERS */
/* VK          FGPM(8)-FGPM(15) - DEPTHS OF SOIL LAYERS (M), NEGATIVE. */
/* VK                             FIRST LAYER (RESIDUE) DEPTH=-0.03M */
/* --      COMMON/FSMCO1/UZTWC,UZFWC,LZTWC,LZFSC,LZFPC,ADIMC,FGCO(6),RSUM(7), */
/* --     1   PPE,PSC,PTA,PWE,PSH,TSOIL(8) */
/* --      COMMON/FSUMS1/SROT,SIMPVT,SRODT,SROST,SINTFT,SGWFP,SGWFS,SRECHT, */
/* --     1              SETT,SE1,SE3,SE4,SE5 */

/*    ================================= RCS keyword statements ========== */
/* --      CHARACTER*68     RCSKW1,RCSKW2 */
/* --      DATA             RCSKW1,RCSKW2 /                                 ' */
/* --     .$Source: /fs/hsmb5/hydro/CVS_root/gorms/sac/fland1.f,v $ */
/* --     . $',                                                             ' */
/* --     .$Id: fland1.f,v 3.1 2012-02-02 21:32:13 zcui Exp $ */
/* --     . $' / */
/*    =================================================================== */
/* VK  ADDITION FOR FROZEN DEPTH ESTIMATION */
/* --      SAVE IIPREV,FRSTPREV,FRZPREV */
/* --      IF(FGCO(1) .EQ. 0.) THEN */
/* --       IIPREV=0 */
/* --       FRSTPREV=0. */
/* --       FRZPREV=0. */
/* --      ENDIF */
/* VK-------------------------------------- */
/*       write(*,*) ' FRZPAR:',(frzpar(ii),ii=1,6) */
/* define major parameters from the array */
/*      WRITE(*,*) 'FLAND1: SACPAR: ', (SACPAR(i), i=1, 17 ) */
    /* Parameter adjustments */

    //--sacpar;
    //--sacst;

	error = 0;

	/* Function Body */
	if (pxv < 0)
	{
		printf(" -ve precipitation at x = %f, y = %f ", hrapx, hrapy);
		error = 1;
		return;
	}

	for (i__ = 0; i__ < 6; ++i__)
	{
		if (sacst[i__] < -1.f)
		{
			printf(" -ve SAC State i= %d,  %f  at fland1 start\n", i__, sacst[i__]);
			error = 1;
			//return;
		}
	}
	uztwm = sacpar[0];
	uzfwm = sacpar[1];
	adimp = sacpar[4];
	lztwm = sacpar[8];
	lzfsm = sacpar[9];
	lzfpm = sacpar[10];
	parea = 1.f - sacpar[3] - adimp;
	
	smcdry = 0.f;
	/* define states from the array */
	uztwc = sacst[0];
	uzfwc = sacst[1];
	lztwc = sacst[2];
	lzfsc = sacst[3];
	lzfpc = sacst[4];
	adimc = sacst[5];
/*      WRITE(*,*) 'FLAND1: UZFWC= ', UZFWC */
	
		/* VK  OLD FROZEN GROUND VERSION: KEEP UNFROZEN WATER = TOTAL */
		/* --       RUZICE=UZK */
		/* --       RLZICE=LZSK */
		/* --       RUZPERC=1.0 */
	uztwh = uztwc;
	uzfwh = uzfwc;
	lztwh = lztwc;
	lzfsh = lzfsc;
	lzfph = lzfpc;
	
/*      if(uztwc .ne. uztwh) WRITE(*,*) 'ST1=',uztwc,uztwh */
/*      if(uzfwc .ne. uzfwh) WRITE(*,*) 'ST2=',uzfwc,uzfwh */
/*      if(lztwc .ne. lztwh) WRITE(*,*) 'ST3=',lztwc,lztwh */
/*      if(lzfsc .ne. lzfsh) WRITE(*,*) 'ST4=',lzfsc,lzfsh */
/*      if(lzfpc .ne. lzfph) WRITE(*,*) 'ST5=',lzfpc,lzfph */
/* ....................................... */
/*     COMPUTE EVAPOTRANSPIRATION LOSS FOR THE TIME INTERVAL. */
/*        EDMND IS THE ET-DEMAND FOR THE TIME INTERVAL */
/* c      EDMND=EP*EPDIST(KINT) */
/* VK ADJUST EDMND FOR EFFECT OF SNOW & FOREST COVER. */
/* from EX1 OFS subroutine */
	edmnd = (1.f - (1.f - sacpar[16]) * aesc) * edmnd;

	/*     COMPUTE ET FROM UPPER ZONE. */
	/* VK      E1=EDMND*(UZTWC/UZTWM) */
	/* VK  ONLY UNFROZEN WATER CAN BE EVAPORATED */
	e1 = edmnd * (uztwh / uztwm);
	red = edmnd - e1;
	/*     RED IS RESIDUAL EVAP DEMAND */
	/* VK      UZTWC=UZTWC-E1 */
	uztwh -= e1;
	e2 = 0.f;
	/* V.K      IF(UZTWC.GE.0.) THEN */
	if (uztwh >= 0.f) {
		/* V.K    SUBTRACT ET FROM TOTAL WATER STORAGE */
		uztwc -= e1;
		goto L220;
	}
	/*     E1 CAN NOT EXCEED UZTWC */
	/* V.K      E1=E1+UZTWC */
	/* V.K      UZTWC=0.0 */
	e1 += uztwh;
	uztwh = 0.f;
	/* V.K   REDUCE TOTAL TENSION WATER BY ACTUAL E1 */
	uztwc -= e1;
	if (uztwc < 0.f) {
		uztwc = 0.f;
	}
	red = edmnd - e1;
	/* V.K      IF(UZFWC.GE.RED) GO TO 221 */
	if (uzfwh >= red) {
		goto L221;
	}
	/*     E2 IS EVAP FROM UZFWC. */
	/* V.K      E2=UZFWC */
	/* V.K      UZFWC=0.0 */
	e2 = uzfwh;
	uzfwh = 0.f;
	/* V.K   REDUCE TOTAL FREE WATER BY ACTUAL E2 */
	uzfwc -= e2;
	if (uzfwc < 0.f) {
		uzfwc = 0.f;
	}
	red -= e2;
	goto L225;
L221:
	e2 = red;
	/* VK   SUBTRACT E2 FROM TOTAL & UNFROZEN FREE WATER STORAGES */
	uzfwc -= e2;
	uzfwh -= e2;
	red = 0.f;
L220:
	if (uztwc / uztwm >= uzfwc / uzfwm) {
		goto L225;
	}
	/*     UPPER ZONE FREE WATER RATIO EXCEEDS UPPER ZONE */
	/*     TENSION WATER RATIO, THUS TRANSFER FREE WATER TO TENSION */
	uzrat = (uztwc + uzfwc) / (uztwm + uzfwm);
	/* V.K  ACCOUNT FOR RATIO OF UNFROZEN WATER ONLY */
	/* V.K  AND ADJUST FOUR SOIL STATES */
	/* V.K      UZTWC=UZTWM*UZRAT */
	/* V.K      UZFWC=UZFWM*UZRAT */
	duztwc = uztwm * uzrat - uztwc;
	if (duztwc > uzfwh) {
		duztwc = uzfwh;
	}
	/* V.K  TRANSFERED WATER CAN NOT EXCEED UNFROZEN FREE WATER */
	uztwc += duztwc;
	uztwh += duztwc;
	uzfwc -= duztwc;
	uzfwh -= duztwc;
	/* V.K  CHECK UNFROZEN WATER STORAGES TOO */
L225:
	if (uztwc < 1e-5f) {
		uztwc = 0.f;
		uztwh = 0.f;
	}
	if (uzfwc < 1e-5f) {
		uzfwc = 0.f;
		uzfwh = 0.f;
	}

	/*     COMPUTE ET FROM THE LOWER ZONE. */
	/*     COMPUTE ET FROM LZTWC (E3) */
	/* V.K      E3=RED*(LZTWC/(UZTWM+LZTWM)) */
	/* V.K      LZTWC=LZTWC-E3 */
	/* V.K      IF(LZTWC.GE.0.0) THEN */
	/* V.K  ONLY UNFROZEN WATER CAN BE EVAPORATED */
	e3 = red * (lztwh / (uztwm + lztwm));
	lztwh -= e3;
	if (lztwh >= 0.f) {
		lztwc -= e3;
		goto L226;
	}
	/*     E3 CAN NOT EXCEED LZTWC */
	/* V.K      E3=E3+LZTWC */
	/* V.K      LZTWC=0.0 */
	e3 += lztwh;
	lztwh = 0.f;
	/* V.K   REDUCE TOTAL TENSION WATER BY E3 */
	lztwc -= e3;
L226:
	ratlzt = lztwc / lztwm;
	ratlz = (lztwc + lzfpc + lzfsc - sacpar[15]) / (lztwm + lzfpm + lzfsm - sacpar[15]);
	if (ratlzt >= ratlz) {
		goto L230;
	}
	/*     RESUPPLY LOWER ZONE TENSION WATER FROM LOWER */
	/*     ZONE FREE WATER IF MORE WATER AVAILABLE THERE. */
	del = (ratlz - ratlzt) * lztwm;
	/* V.K  ONLY UNFROZEN WATER CAN BE TRANSFERED */
	/*       if(lzfsc .ne. lzfsh) write(*,*) 'BST4=',lzfsc,lzfsh */
	sfh = lzfsh + lzfph;
	if (del > sfh) {
		del = sfh;
	}
	lzfsh -= del;
	if (lzfsh >= 0.f) {
		/*     TRANSFER FROM LZFSC TO LZTWC. */
		lzfsc -= del;
		/*         if(lzfsc .lt. lzfsh) then */
		/*          write(*,*) ' lzfsc1: ',lzfsc,lzfsh,del */
		/*          stop */
		/*         endif */
	}
	else {
		/*     IF TRANSFER EXCEEDS LZFSC THEN REMAINDER COMES FROM LZFPC */
		lzfpc += lzfsh;
		lzfph += lzfsh;
		xx = lzfsh + del;
		lzfsc -= xx;
		/*         if(lzfsc .lt. lzfsh) then */
		/*          write(*,*) ' lzfsc2: ',lzfsc,lzfsh,del,xx */
		/*          stop */
		/*         endif */
		lzfsh = 0.f;
	}
	lztwc += del;
	lztwh += del;
	/* V.K      LZTWC=LZTWC+DEL */
	/* V.K      LZFSC=LZFSC-DEL */
	/* V.K      IF(LZFSC.GE.0.0) GO TO 230 */
	/* V.K      LZFPC=LZFPC+LZFSC */
	/* V.K      LZFSC=0.0 */
	/* V.K  CHECK UNFROZEN WATER STORAGE */
L230:
	if (lztwc < 1e-5f) {
		lztwc = 0.f;
		lztwh = 0.f;
	}

	/*     COMPUTE ET FROM ADIMP AREA.-E5 */
	e5 = e1 + (red + e2) * ((adimc - e1 - uztwc) / (uztwm + lztwm));
	/*      ADJUST ADIMC,ADDITIONAL IMPERVIOUS AREA STORAGE, FOR EVAPORATION. */
	adimc -= e5;
	if (adimc >= 0.f) {
		goto L231;
	}
	/*     E5 CAN NOT EXCEED ADIMC. */
	e5 += adimc;
	adimc = 0.f;
L231:
	e5 *= adimp;
	/*     E5 IS ET FROM THE AREA ADIMP. */
	/* ....................................... */
	/*     COMPUTE PERCOLATION AND RUNOFF AMOUNTS. */
	twx = pxv + uztwc - uztwm;
	/*     TWX IS THE TIME INTERVAL AVAILABLE MOISTURE IN EXCESS */
	/*     OF UZTW REQUIREMENTS. */
	if (twx >= 0.f) {
		goto L232;
	}
	/*     ALL MOISTURE HELD IN UZTW--NO EXCESS. */
	uztwc += pxv;
	/* V.K  ADJUST UNFROZEN TENSION WATER */
	uztwh += pxv;
	twx = 0.f;
	goto L233;
	/*      MOISTURE AVAILABLE IN EXCESS OF UZTWC STORAGE. */
	/* V.K  232 UZTWC=UZTWM */
L232:
	uztwh += uztwm - uztwc;
	uztwc = uztwm;
L233:
	adimc = adimc + pxv - twx;

	/*     COMPUTE IMPERVIOUS AREA RUNOFF. */
	roimp = pxv * sacpar[3];
	/*      ROIMP IS RUNOFF FROM THE MINIMUM IMPERVIOUS AREA. */
	simpvt += roimp;

	/*     INITIALIZE TIME INTERVAL SUMS. */
	sbf = 0.f;
	ssur = 0.f;
	sif = 0.f;
	sperc = 0.f;
	sdro = 0.f;
	spbf = 0.f;

	/*     DETERMINE COMPUTATIONAL TIME INCREMENTS FOR THE BASIC TIME */
	/*     INTERVAL */
	/* V.K      NINC=1.0+0.2*(UZFWC+TWX) */
	/* V.K  PERCOLATE UNFROZEN WATER ONLY */
	ninc = (uzfwh + twx) * .2f + 1.f;
	/*     NINC=NUMBER OF TIME INCREMENTS THAT THE TIME INTERVAL */
	/*     IS DIVIDED INTO FOR FURTHER */
	/*     SOIL-MOISTURE ACCOUNTING.  NO ONE INCREMENT */
	/*     WILL EXCEED 5.0 MILLIMETERS OF UZFWC+PAV */
	dinc = 1.f / ninc * dt;
	/*     DINC=LENGTH OF EACH INCREMENT IN DAYS. */
	pinc = twx / ninc;
	/*     PINC=AMOUNT OF AVAILABLE MOISTURE FOR EACH INCREMENT. */
/*      COMPUTE FREE WATER DEPLETION FRACTIONS FOR */
/*     THE TIME INCREMENT BEING USED-BASIC DEPLETIONS */
/*      ARE FOR ONE DAY */
/* VK INTRODUCED REDUCTION (RUZICE & RLZICE) DUE FROZEN GROUND */
/* VK HOWEVER, PRIMARY RUNOFF IS UNCHANGED */
/* VK      DUZ=1.0-((1.0-UZK)**DINC) */
/* VK      DLZS=1.0-((1.0-LZSK)**DINC) */
/* VK  Linear transformation for frozen ground */
/* c      DUZ=1.0-((1.0-UZK*RUZICE)**DINC) */
/* c      DLZS=1.0-((1.0-LZSK*RLZICE)**DINC) */
/* VK  Non-linear (correct) transformation for frozen ground */
	
		d__1 = (double)(1.f - sacpar[2]);
		d__2 = (double)dinc;
		duz = 1.f - pow(d__1, d__2);
		d__1 = (double)(1.f - sacpar[11]);
		d__2 = (double)dinc;
		dlzs = 1.f - pow(d__1, d__2);
	
	d__1 = (double)(1.f - sacpar[12]);
	d__2 = (double)dinc;
	dlzp = 1.f - pow(d__1, d__2);

	/* CVK  ADJUSTMENT TO DEPLETIONS DUE TO FROZEN WATER */
	/* ....................................... */
	/*     START INCREMENTAL DO LOOP FOR THE TIME INTERVAL. */
	i__1 = ninc;
	for (i__ = 0; i__ < i__1; ++i__) 
	{
		adsur = 0.f;
		/*     COMPUTE DIRECT RUNOFF (FROM ADIMP AREA). */
		ratio = (adimc - uztwc) / lztwm;
		if (ratio < 0.f) {
			ratio = 0.f;
		}
		/* Computing 2nd power */
		r__1 = ratio;
		addro = pinc * (r__1 * r__1);
		/*     ADDRO IS THE AMOUNT OF DIRECT RUNOFF FROM THE AREA ADIMP. */

		/*     COMPUTE BASEFLOW AND KEEP TRACK OF TIME INTERVAL SUM. */
		/* V.K      BF=LZFPC*DLZP */
		/* V.K      LZFPC=LZFPC-BF */
		/* V.K      IF (LZFPC.GT.0.0001) GO TO 234 */
		/* V.K      BF=BF+LZFPC */
		/* V.K      LZFPC=0.0 */
		/* V.K  BASEFLOW FROM UNFROZEN WATER ONLY */
		bf = lzfph * dlzp;
		lzfph -= bf;
		if (lzfph > 1e-4f) {
			lzfpc -= bf;
			goto L234;
		}
		bf += lzfph;
		lzfph = 0.f;
		lzfpc -= bf;
		if (lzfpc <= 1e-4f) {
			lzfpc = 0.f;
		}
		/* V.K------------------------------------- */

	L234:
		sbf += bf;
		spbf += bf;
		/* V.K  SUPPLAMENTAL FLOW FROM UNFROZEN WATER ONLY (NOTE, DLZS */
		/* V.K  NOTE, DLZS IS REDUCED DUE FROZEN GROUND */
		/* V.K      BF=LZFSC*DLZS */
		/* V.K      LZFSC=LZFSC-BF */
		/* V.K      IF(LZFSC.GT.0.0001) GO TO 235 */
		/* V.K      BF=BF+LZFSC */
		/* V.K      LZFSC=0.0 */
		bf = lzfsh * dlzs;
		lzfsh -= bf;
		if (lzfsh > 1e-4f) {
			/* c?      IF(LZFSH.GT.0.0) THEN */
			lzfsc -= bf;
			/*         if(abs(lzfsc-lzfsh) .gt. 0.000001) then */
			/*         if(abs(lzfsc-lzfsh) .gt. 0.000001) then */
			/*          write(*,*) ' lzfsc3: ',lzfsc,lzfsh,bf */
			/*         endif */
			goto L235;
		}
		bf += lzfsh;
		lzfsh = 0.f;
		lzfsc -= bf;
		if (lzfsc <= 1e-4f) {
			lzfsc = 0.f;
		}
		/* V.K-------------------------------------------- */

	L235:
		sbf += bf;

		/*      COMPUTE PERCOLATION-IF NO WATER AVAILABLE THEN SKIP */
		/* cvk      IF((PINC+UZFWC).GT.0.01) GO TO 251 */
		xx1 = pinc + uzfwh;
		if (xx1 > .01f) {
			goto L251;
		}
		uzfwc += pinc;
		/* V.K  ADD TO UNFROZEN WATER ALSO */
		uzfwh += pinc;
		goto L249;
	L251:
		percm = lzfpm * dlzp + lzfsm * dlzs;
		/* VK      PERC=PERCM*(UZFWC/UZFWM) */
		/* V.K  USE ONLY UNFROZEN WATER RATIOS */
		/* cvk  new change: PERCOLATION REDUCED BY RUZPERC */
		/* C       PERC=PERCM*(UZFWH/UZFWM)*RUZICE */
		perc = percm * (uzfwh / uzfwm);
		
		/* --      PERC=PERCM*(UZFWH/UZFWM)*RUZPERC */
		/* V.K      DEFR=1.0-((LZTWC+LZFPC+LZFSC)/(LZTWM+LZFPM+LZFSM)) */
		/* vk 6/22/00      DEFR=1.0-((LZTWH+LZFPH+LZFSH)/(LZTWM+LZFPM+LZFSM)) */
		/* vk  better to keep original definition of DEFR using total water */
		defr = 1.f - (lztwc + lzfpc + lzfsc) / (lztwm + lzfpm + lzfsm);
		/*     DEFR IS THE LOWER ZONE MOISTURE DEFICIENCY RATIO */
		/* --      FR=1.0 */
		/*     FR IS THE CHANGE IN PERCOLATION WITHDRAWAL DUE TO FROZEN GROUND. */
		/* --      FI=1.0 */
		/*     FI IS THE CHANGE IN INTERFLOW WITHDRAWAL DUE TO FROZEN GROUND. */
		/* --      IF (IFRZE.EQ.0) GO TO 239 */
		/* --       UZDEFR=1.0-((UZTWC+UZFWC)/(UZTWM+UZFWM)) */
		/* VK */
		/* VK     CALL FGFR1(DEFR,FR,UZDEFR,FI) */
		/* VK      IF( IVERS .EQ. 1) THEN */
		/* VK  IF IVERS=1, OLD VERSION; IF IVERS=2, NEW VERS. FROST INDEX, */
		/* VK  BUT OLD VERS. OF PERCOLAT. AND INTERFLOW REDUCTION */
		/* --      IF( IVERS .LE. 2) CALL FGFR1(DEFR,FR,UZDEFR,FI) */
		/* --      IF(IVERS .EQ. 3 .AND. FGPM(5) .GT. 0.) THEN */
		/* VK  OPTIONAL VERSION TO ACCOUNT FOR ADDITIONAL IMPERVIOUS */
		/* VK  AREAS EFFECTS DUE FROZEN GROUND */
		/* --       FR=1-SURFRZ1(FGCO(1),FGPM(6),FGPM(5)) */
		/* --       FI=FR */
		/* --      ENDIF */
		/* --  239 PERC=PERC*(1.0+ZPERC*(DEFR**REXP))*FR */
		/* L239: */
		d__1 = (double)defr;
		d__2 = (double)sacpar[7];
		perc *= sacpar[6] * pow(d__1, d__2) + 1.f;
		/*     NOTE...PERCOLATION OCCURS FROM UZFWC BEFORE PAV IS ADDED. */
		/* V.K      IF(PERC.LT.UZFWC) GO TO 241 */
		if (perc < uzfwh) {
			goto L241;
		}
		/*      PERCOLATION RATE EXCEEDS UZFWH. */
		/* V.K      PERC=UZFWC */
		perc = uzfwh;
		/*     PERCOLATION RATE IS LESS THAN UZFWH. */
	L241:
		uzfwc -= perc;
		/* V.K  ADJUST UNFROZEN STORAGE ALSO */
		uzfwh -= perc;
		/*     CHECK TO SEE IF PERCOLATION EXCEEDS LOWER ZONE DEFICIENCY. */
		check = lztwc + lzfpc + lzfsc + perc - lztwm - lzfpm - lzfsm;
		if (check <= 0.f) {
			goto L242;
		}
		perc -= check;
		uzfwc += check;
		/* V.K  ADJUST UNFROZEN STARAGE ALSO */
		uzfwh += check;
	L242:
		sperc += perc;
		/*     SPERC IS THE TIME INTERVAL SUMMATION OF PERC */

		/*     COMPUTE INTERFLOW AND KEEP TRACK OF TIME INTERVAL SUM. */
		/*     NOTE...PINC HAS NOT YET BEEN ADDED */
		/* V.K      DEL=UZFWC*DUZ*FI */
		/* VK  INTERFLOW ALSO REDUCED DUE FROFEN GROUND (DUZ REDUCED BY RUZICE) */
		/* VK  ADDITIONAL REDUCTION DUE IMPERVIOUS FROZEN AREAS (FI) IS OPTIONAL */
		/* VK  IN THE NEW VERSION. BASIC OPTION IS FI=1 */
		/* --      DEL=UZFWH*DUZ*FI */
		del = uzfwh * duz;
		sif += del;
		uzfwc -= del;
		/* V.K  ADJUST UNFROZEN STORAGE ALSO */
		uzfwh -= del;
		/*     DISTRIBE PERCOLATED WATER INTO THE LOWER ZONES */
		/*     TENSION WATER MUST BE FILLED FIRST EXCEPT FOR THE PFREE AREA. */
		/*     PERCT IS PERCOLATION TO TENSION WATER AND PERCF IS PERCOLATION */
		/*         GOING TO FREE WATER. */
		perct = perc * (1.f - sacpar[13]);
		xx1 = perct + lztwc;
		if (xx1 > lztwm) {
			goto L243;
		}
		lztwc += perct;
		/* V.K  ADJUST UNFROZEN STORAGE ALSO */
		lztwh += perct;
		percf = 0.f;
		goto L244;
	L243:
		percf = perct + lztwc - lztwm;
		/* V.K  CHANGE UNFROZEN WATER STORAGE */
		lztwh = lztwh + lztwm - lztwc;
		lztwc = lztwm;

		/*      DISTRIBUTE PERCOLATION IN EXCESS OF TENSION */
		/*      REQUIREMENTS AMONG THE FREE WATER STORAGES. */
	L244:
		percf += perc * sacpar[13];
		if (percf == 0.f) {
			goto L245;
		}
		hpl = lzfpm / (lzfpm + lzfsm);
		/*     HPL IS THE RELATIVE SIZE OF THE PRIMARY STORAGE */
		/*     AS COMPARED WITH TOTAL LOWER ZONE FREE WATER STORAGE. */
		/* VK changed to account for ZERO MAX storage */
		if (lzfpm != 0.f) {
			ratlp = lzfpc / lzfpm;
		}
		else {
			ratlp = 1.f;
		}
		if (lzfsm != 0.f) {
			ratls = lzfsc / lzfsm;
		}
		else {
			ratls = 1.f;
		}
		/*     RATLP AND RATLS ARE CONTENT TO CAPACITY RATIOS, OR */
		/*     IN OTHER WORDS, THE RELATIVE FULLNESS OF EACH STORAGE */
		fracp = hpl * 2.f * (1.f - ratlp) / (1.f - ratlp + (1.f - ratls));
		/*     FRACP IS THE FRACTION GOING TO PRIMARY. */
		if (fracp > 1.f) {
			fracp = 1.f;
		}
		percp = percf * fracp;
		percs = percf - percp;
		/*     PERCP AND PERCS ARE THE AMOUNT OF THE EXCESS */
		/*     PERCOLATION GOING TO PRIMARY AND SUPPLEMENTAL */
		/*      STORGES,RESPECTIVELY. */
		lzfsc += percs;
		/* V.K      IF(LZFSC.LE.LZFSM) GO TO 246 */
		if (lzfsc <= lzfsm) {
			lzfsh += percs;
			goto L246;
		}
		percs = percs - lzfsc + lzfsm;
		/* V.K  ADJUST UNFROZEN STORAGE ALSO */
		lzfsh += percs;
		lzfsc = lzfsm;
	L246:
		lzfpc += percf - percs;
		/*     CHECK TO MAKE SURE LZFPC DOES NOT EXCEED LZFPM. */
		/* V.K      IF (LZFPC.LE.LZFPM) GO TO 245 */
		if (lzfpc <= lzfpm) {
			lzfph += percf - percs;
			goto L245;
		}
		excess = lzfpc - lzfpm;
		lztwc += excess;
		/* V.K  ADJUST UNFROZEN STORAGES ALSO */
		lztwh += excess;
		lzfph = lzfph + (percf - percs) - excess;
		lzfpc = lzfpm;

		/*     DISTRIBUTE PINC BETWEEN UZFWC AND SURFACE RUNOFF. */
	L245:
		if (pinc == 0.f) {
			goto L249;
		}
		/*     CHECK IF PINC EXCEEDS UZFWM */
		xx1 = pinc + uzfwc;
		if (xx1 > uzfwm) {
			goto L248;
		}
		/*     NO SURFACE RUNOFF */
		uzfwc += pinc;
		/* V.K  ADJUST UNFROZEN STORAGE ALSO */
		uzfwh += pinc;
		goto L249;

		/*     COMPUTE SURFACE RUNOFF (SUR) AND KEEP TRACK OF TIME INTERVAL SUM. */
	L248:
		sur = pinc + uzfwc - uzfwm;
		uzfwc = uzfwm;
		/* V.K  ADJUST UNFROZEN STORAGE ALSO */
		uzfwh = uzfwh + pinc - sur;
		ssur += sur * parea;
		adsur = sur * (1.f - addro / pinc);
		/*     ADSUR IS THE AMOUNT OF SURFACE RUNOFF WHICH COMES */
		/*     FROM THAT PORTION OF ADIMP WHICH IS NOT */
		/*     CURRENTLY GENERATING DIRECT RUNOFF.  ADDRO/PINC */
		/*     IS THE FRACTION OF ADIMP CURRENTLY GENERATING */
		/*     DIRECT RUNOFF. */
		ssur += adsur * adimp;

		/*     ADIMP AREA WATER BALANCE -- SDRO IS THE 6 HR SUM OF */
		/*          DIRECT RUNOFF. */
	L249:
		adimc = adimc + pinc - addro - adsur;
		xx1 = uztwm + lztwm;
		if (adimc <= xx1) {
			goto L247;
		}
		addro = addro + adimc - xx1;
		adimc = xx1;
	L247:
		sdro += addro * adimp;
		if (adimc < 1e-5f) {
			adimc = 0.f;
		}
		/* L240: */
	}
	/* ....................................... */
	/*     END OF INCREMENTAL DO LOOP. */
	/* ....................................... */
	/*     COMPUTE SUMS AND ADJUST RUNOFF AMOUNTS BY THE AREA OVER */
	/*     WHICH THEY ARE GENERATED. */
	eused = e1 + e2 + e3;
	/*     EUSED IS THE ET FROM PAREA WHICH IS 1.0-ADIMP-PCTIM */
	sif *= parea;

	/*     SEPARATE CHANNEL COMPONENT OF BASEFLOW */
	/*     FROM THE NON-CHANNEL COMPONENT */
	tbf = sbf * parea;
	/*     TBF IS TOTAL BASEFLOW */
	bfcc = tbf * (1.f / (sacpar[14] + 1.f));
	/*     BFCC IS BASEFLOW, CHANNEL COMPONENT */
	bfp = spbf * parea / (sacpar[14] + 1.f);
	bfs = bfcc - bfp;
	if (bfs < 0.f) {
		bfs = 0.f;
	}
	bfncc = tbf - bfcc;
	/*     BFNCC IS BASEFLOW,NON-CHANNEL COMPONENT */

/*     ADD TO MONTHLY SUMS. */
/* --      SINTFT=SINTFT+SIF */
/* --      SGWFP=SGWFP+BFP */
/* --      SGWFS=SGWFS+BFS */
/* --      SRECHT=SRECHT+BFNCC */
/* --      SROST=SROST+SSUR */
/* --      SRODT=SRODT+SDRO */

	/*     COMPUTE TOTAL CHANNEL INFLOW FOR THE TIME INTERVAL. */
	tci = roimp + sdro + ssur + sif + bfcc;
	grnd = sif + bfcc;
	/* C	GRND = BFCC         ! interflow is part of surface flow */
	/* interflow is part of ground flow */
	surf = tci - grnd;

	/*     COMPUTE E4-ET FROM RIPARIAN VEGETATION. */
	e4 = (edmnd - eused) * sacpar[5];

	/*     SUBTRACT E4 FROM CHANNEL INFLOW */
	tci -= e4;
	if (tci >= 0.f) {
		goto L250;
	}
	e4 += tci;
	tci = 0.f;
	/* c  250 SROT=SROT+TCI */
L250:
	grnd -= e4;
	if (grnd < 0.f) {
		surf += grnd;
		grnd = 0.f;
		if (surf < 0.f) {
			surf = 0.f;
		}
	}

	/*     COMPUTE TOTAL EVAPOTRANSPIRATION-TET */
	eused *= parea;
	tet = eused + e5 + e4;
/* --      SETT=SETT+TET */
/* --      SE1=SE1+E1*PAREA */
/* --      SE3=SE3+E3*PAREA */
/* --      SE4=SE4+E4 */
/* --      SE5=SE5+E5 */
	/*     CHECK THAT ADIMC.GE.UZTWC */
	if (adimc < uztwc) {
		adimc = uztwc;
	}

	/*  Return back SAC states */
	sacst[0] = uztwc;
	sacst[1] = uzfwc;
	sacst[2] = lztwc;
	sacst[3] = lzfsc;
	sacst[4] = lzfpc;
	sacst[5] = adimc;
	/* new change: check negative states */
	for (i__ = 0; i__ < 6; ++i__)
	{
		if (sacst[i__] < -1.f)
		{
			printf(" -ve SAC State i = %d,  %f at fland1 end\n", i__, sacst[i__]);
			error = 1;
			return;
		}
		if (sacst[i__] < 0.f)
		{
			sacst[i__] = 0.f;
		}
	}
	if (uztwh < 0.f) {
		uztwh = 0.f;
	}
	if (uzfwh < 0.f) {
		uzfwh = 0.f;
	}
	if (lztwh < 0.f) {
		lztwh = 0.f;
	}
	if (lzfsh < 0.f) {
		lzfsh = 0.f;
	}
	if (lzfph < 0.f) {
		lzfph = 0.f;
	}
	if (sacst[0] < uztwh) {
		uztwh = sacst[0];
	}
	if (sacst[1] < uzfwh) {
		uzfwh = sacst[1];
	}
	if (sacst[2] < lztwh) {
		lztwh = sacst[2];
	}
	if (sacst[3] < lzfsh) {
		lzfsh = sacst[3];
	}
	if (sacst[4] < lzfph) {
		lzfph = sacst[4];
	}

	return;
} /* fland1_ */