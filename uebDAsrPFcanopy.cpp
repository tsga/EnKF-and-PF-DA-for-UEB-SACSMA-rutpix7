
//copied from Canopy.f90
//6.6.13
#include "uebDAsrPFuebpgdecls.h"
//*******************************************************************************************************************
//*		               CANOPY RADIATION TRANSMISSION PROCESSES
//*******************************************************************************************************************
//
//  Copyright (C) 2012  David Tarboton, Utah State University, dtarb@usu.edu.  http://www.engineering.usu.edu/dtarb
//
//  This file is part of UEB.
//
//    UEB is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    UEB is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    A copy of the GNU General Public License is included in the file gpl.txt.
//    This is also available at: http://www.gnu.org/licenses/.
//
//  if you wish to use or incorporate this program (or parts of it) into 
//  other software that does not meet the GNU General Public License 
//  conditions contact the author to request permission.
//  David G. Tarboton  
//  Utah State University 
//  8200 Old Main Hill 
//  Logan, UT 84322-8200 
//  USA 
//  http://www.engineering.usu.edu/dtarb/ 
//  email:  dtarb@usu.edu 
//
//********** PARTITION OF INCOMING SOLAR RADIATION INTO DIRECT AND DIFFUSE BEFORE IT HITS THE CANOPY *************       
//     Partition the incoming solar radiation into direct and diffuse components
 __host__ __device__ 
 void uebCellDA::PSOLRAD( float Qsi, float atff, float *param, float cf,     
                float &Taufb, float &Taufd, float &Qsib, float &Qsid)     // Output variables:
{

    float  Taufc, as,bs,Lambda;
		
    as     = param[27];			// Fraction of extraterrestaial radiation on cloudy day,Shuttleworth (1993)  
	bs     = param[28];		// (as+bs):Fraction of extraterrestaial radiation on clear day, Shuttleworth (1993) 
	Lambda = param[29];		    // Ratio of diffuse atm radiation to direct for cler sky,from Dingman (1993)

   Taufc = findMax (atff, as+bs);     // Mac atmosphereci transmittance for clear sky  
   if (cf == 0)
   { 
		Taufb = Lambda*Taufc;				 // Direct solar radiation atm.transmittance
		Taufd = (1.0 - Lambda)*Taufc;			 // Diffuse solar radiation atm.transmittance
   }
   else if (cf == 1)
         {
		    Taufb = 0.0;						     
		    Taufd = atff;
   }
   else 
   {
		Taufb = Lambda*(1.0 - cf)*(as+bs);
		Taufd = atff- Taufb;
   }    
//     Direct canopy solar radiation incident at the canopy 
   	Qsib  = Taufb*Qsi/(Taufb+Taufd);
//	Diffuse canopy solar radiation incident at the canopy     
	Qsid  = Taufd*Qsi/(Taufb+Taufd);    

	if (snowdgtvariteflag2 == 1)
	{
		std::cout << "Outputs from PSOLRad" << std::endl;
		std::cout << std::setprecision(15) << Taufb << " " << Taufd << " "<<Qsi<<" "<< Qsib << " " << Qsid << std::endl;
	}

	return; 
}

//*********** TRANSMISSION OF DIRECT AND DIFFUSE SOLAR RADIATION THROUGH THE CANOPY ***************
//     Estimates the direct and diffuse solar radiation fractions transmitted through the canopy      
 __host__ __device__   
 void uebCellDA::TRANSRADCAN (float COSZEN, float *sitev, float *param,       
     	            float &Betab, float &Betad, float &Taub, float &Taud)                     // Output variables: Transmission and Reflection fractions
{    
  float  EXPI, kk, Taubh,Taudh, Beta,Rho;

	float Cc    = sitev[4],     // Leaf Area Index 
	Hcan  = sitev[5],     // Canopy height  
	LAI   = sitev[6],     // Leaf Area Index  
	alpha = param[23],    // 
      G     = param[25];    // 

//      For open area where no vegetaion canopy is present                                       
	if (LAI == 0 )
	{
	  Taub = 1.0;
	  Taud = 1.0;
      Betab = 0.0;
      Betad = 0.0;
	}
	else
	{
//   Transmission without scattering 
      LAI   = Cc*LAI;                      // Effective leaf area index
      Rho   = LAI/Hcan;                    // leaf density per unit height of the canopy
//     Multiple scattering light penetration
	 kk    = sqrt(1 - alpha);                             // k-prime in the paper
	 EXPI  = EXPINT(kk*G*LAI);                           // Exponential integral function	
//      Deep Canopy Solution
//   Avoid divide by 0 error when COSZEN is 0
       if(COSZEN <= 0)
           Taubh=0.0;
	   else
           Taubh =  expf(-1*kk*G*Rho*Hcan/COSZEN);                 // Transmission function for deep canopy : Direct

       Taudh =  (1 - kk*G*Rho*Hcan) * expf(-1*kk*G*Rho*Hcan) + powf((kk*G*Rho*Hcan),2.0)*EXPI;                     // Transmission function for deep canopy : Diffuse
       Beta   = (1-kk)/(1+kk);                                // BetaPrime : reflection coefficient for infinitely deep canopy
//     Finite Canopy Solution
       Taub   = Taubh*(1 - powf(Beta, 2.0)) / (1 - powf(Beta, 2.0)*powf(Taubh, 2.0));     // Direct fraction of radiation transmitted down through the canopy 
       Taud   = Taudh*(1-powf(Beta, 2.0))/(1-powf(Beta, 2.0)*powf(Taudh, 2.0));         // Diffuse fraction of radiation transmitted down through the canopy 
       Betab  = Beta*(1-powf(Taubh, 2.0)) /(1-powf(Beta, 2.0)*powf(Taubh, 2.0));    // Direct fraction of radiation scattered up from the canopy 
       Betad  = Beta*(1-powf(Taudh, 2.0)) /(1-powf(Beta, 2.0)*powf(Taudh, 2.0));       // Diffuse fraction of radiation scattered up from the canopy
	}
	return;
}

//****************** NET CANOPY & SUB CANOPY SOLAR RADIATION *********************************
//      Computes the net canopy and sub-canopy solar radiation
//      considering the multiple scatterings of radiation by the canopy and multile reflections 
//      between the canopy and the snow surface
 __host__ __device__   
 void  uebCellDA::NETSOLRAD(float Ta,float A,float Betab,float Betad,float Wc,float Taub,float Taud, 
			     float Qsib,float Qsid,float *param, float	Fs,float &Qsns,float &Qsnc) //  Output: Qsns,Qsnc (Net subcanopy, canopy solar radiation) 
{      
	float  f1b, f2b, f3b, f1d, f2d, f3d;
//	Net radiation at snow surface below canopy
	f1b = (1-A)*Taub/(1-Betad*(1-Taud)*A) ;
    f1d = (1-A)*Taud/(1-Betad*(1-Taud)*A);              
    Qsns =  f1b*Qsib + f1d*Qsid;

	if (Taub == 1)              // (For LAI=0, Qsnc gives the numerical errors, so)
		Qsnc = 0.0;
	else
	{      
//     Net canopy solar radiation
		f3b = ((1-A)*Betab + A*Taub*Taud)/(1-Betad*(1-Taud)*A);
		f3d = ((1-A)*Betad + A*Taud*Taud)/(1-Betad*(1-Taud)*A); 
		f2b = 1- f1b-f3b;  
		f2d = 1- f1d-f3d;
		Qsnc = f2b*Qsib + f2d*Qsid;
	}
	if (snowdgtvariteflag2 == 1)
	{
		std::cout << "Outputs from NetSOLRad" << std::endl;
		std::cout << std::setprecision(15) << Taub << " " << Taud << " " << Qsib << " " << Qsid << " " << Qsns << " " << Qsnc << std::endl;
	}
	return;
}

//****************** NET CANOPY & SUB CANOPY LONGWAVE RADIATION ******************************
//     Computes the net canopy and beneath canopy longwave radiation
 __host__ __device__   
 void  uebCellDA::NETLONGRAD(float RH,float Ta,float Tss,float Tc,float Tk,float Fs,float EmC,float EmS,float SBC,float cf,float *sitev,float Qli,float *param,
				    float &Qlis, float &Qlns,float &Qlnc )                     //  Output: Qsns,Qsnc 
{
 
	float  EXPI, kk, Beta, Betad, Tak, Tssk, EmCS, Tck, Qle, Qlcd, Qlcu, Rho, Taudh, Taud;

	float Cc    = sitev[4],     // Leaf Area Index 
	Hcan  = sitev[5],     // Canopy height  
	LAI   = sitev[6],     // Leaf Area Index  
	alphaL= param[24],  
    G     = param[25]; 
//	ea   = SVPW(TA)*RH          	// Air vapor pressure 
    Tak  = Ta + Tk;
	Tssk = Tss +Tk;  
    EmCS = EmC*(1-Fs)+EmS*Fs;
	Tck  = Tc+Tk;
    Qle  = EmS*SBC*powf(Tssk,4.0);
	if (LAI == 0)          //  For Open area where no vegetation caonpy is present
	{
	   Qlcd = 0.0;
	   Qlcu = 0.0;
       Qlns = Qli*EmS -Qle;                    // To avoid numerical errors
       Qlnc = 0.0;
	}
	else
	{
//      Transmission without scattering      
	  LAI      = Cc*LAI;
      Rho      = LAI/Hcan; 
//    Downward scattered 
	  kk    = sqrt(1-alphaL);
//    Upward scattered 
      Beta = (1-kk)/(1+kk);
	  EXPI  = EXPINT(kk*G*LAI);
//    Deep canopy solution
	  Taudh =  (1-kk*G*Rho*Hcan)*expf(-kk*G*Rho*Hcan)+ powf((kk*G*Rho*Hcan),2.0)*EXPI;                               // Transmission function for deep canopy : Diffuse
//	 Finite canopy solution
      Taud   = (Taudh-powf(Beta,2.0)*Taudh) /(1-powf(Beta,2.0)*powf(Taudh,2.0));         // Transmission function for finite canopy depth
      Betad  = (Beta-Taudh*Beta*Taudh) /(1-powf(Beta,2.0)*powf(Taudh,2.0));       //  Scatterd up  from top of the canopy at 0 depth  
//      Canopy emitted longwave radiation 
	  Qlcd = (1-Taud)*EmC*SBC*powf(Tck,4.0);
	  Qlcu = (1-Taud)*EmCS*SBC*powf(Tck,4.0);
// Net sub canopy longwave radiation
       Qlns = Taud*Qli*EmS + EmS*Qlcd-Qle + (1.0-Taud)*(1.0 - 1.0)*Qle;  //*****#$%^^%#_ Term on the right 1-1=0? 6.6.13
// Net canopy longwave radiation
      Qlnc = (((1-Taud)*1)+((1-EmS)*Taud*1)) * Qli + (1-Taud)*Qle*1        // Putting 1 for Emc
                  + (1-Taud)*(1-EmS)*1*Qlcu - Qlcu - Qlcd;  
	}
	if (snowdgtvariteflag2 == 1)
	{
		std::cout << "Outputs from NetLongRad" << std::endl;
		std::cout << std::setprecision(15) << Taud << " " << Taudh << " " << Betad << " " << Qlcd << " " << Qlcu << " ";
		std::cout << std::setprecision(15) << Qli << " " << Qle << " " << Qlns << " " << Qlnc << std::endl;
	}
    
	return;
}

//*		               	TURBULENT TRANSFER PROCESSES
//****************  BULK AERODYNAMIC AND CANOPY BOUNDARY LAYER RESISTANCES  *******************
//      Calculates resistances for the above and beneath canopy turbulent heat fluxes 
 __host__ __device__  
 void  uebCellDA::AeroRes(float P, float Wc, float V, float Ta, float Tss, float Tc, float Fs, float *param, float *sitev, float Tk,
			    float &d, float &Z0c, float &Vz, float &Rc, float &Ra, float &Rbc, float &Rl, float &RKINc, float &RKINa, float &RKINbc, float &RKINl)         // Output variables
{
	float  Zm, Kh, Lbmean, Rs, Ri , ndecay, Dc, Rastar, Rbcstar, Rcstar;
	float Vstar, Vh, Vc;

   float z = param[4],       // Nominal meas. height for air temp. and humidity [2m],
      zo    = param[5],       // Surface aerodynamic roughness [m],
	Rimax = param[30],      // Maximum value of Richardsion number for stability corretion
	fstab = param[18],      // Stability correction control parameter 0 = no corrections, 1 = full corrections
      Wcoeff= param[31],      // 

	Cc    = sitev[4],       // Canopy Cover 
	Hcan  = sitev[5],       // Canopy height  
	LAI   = sitev[6];       // Leaf Area Index 

	Zm    = Hcan + z;         // observed wind speed above the canopy
    LAI = Cc*LAI;              
    if(V <= 0) 
	{
        Ra      = 0.0;
        Rc      = 0.0;
        Rs      = 0.0;
        Rbc     = 0.0;
        Rl      = 0.0;
        RKINc   = 0.0;                       
        RKINa   = 0.0;  
		RKINbc  = 0.0;                       
        RKINl   = 0.0; 
		Vz      = V; 
		//  if no wind, all resistances are zero; giving fluxes zero
	}
	else
	{
        if(LAI == 0)      // For open area
		{
            Vz = V;
            Rbc = (1.0/(powf(K_vc,2.0)*Vz)*powf((logf(z/zo)),2.0) )*1/Hs_f;            // sub layer snow resistance; Hs_f: to convert resistance second to hour
		    Ri = Gra_v*(Ta-Tss)*z /(powf(Vz,2.0)*(0.5*(Ta+Tss)+273.15)); 
			if (Ri > Rimax) 
				Ri= Rimax;   
			if (Ri > 0) 
				Rbcstar = Rbc/powf((1-5*Ri),2.0);     // Bulk
			else
      			Rbcstar = Rbc/(powf((1-5*Ri),0.75));		

            RKINbc  = 1.0/Rbcstar;                     // Bulk conductance
		}
		else
		{
             WINDTRANab(V, Zm, param, sitev, Vstar, Vh, d, Z0c, Vz, Vc);               // Output wind speed within canopy at height 2 m
//      Aerodynamic Resistances for Canopy
             ndecay =   Wcoeff*LAI;              //(Wcoeff=0.5 for TWDEF, 0.6 for Niwot) // wind speed decay coefficient inside canopy 
			 Kh = powf(K_vc,2.0)*V*(Hcan-d)/(logf((Zm-d)/Z0c));           
//        Above canopy aerodynamic resistance 
			Ra = (1.0* 1.0/(powf(K_vc,2.0)*V)*logf((Zm-d)/(Z0c))*logf((Zm-d) /(Hcan-d))
				              + Hcan/(Kh*ndecay) * (expf(ndecay-ndecay*(d+Z0c)/Hcan)-1.0) )*1/Hs_f; // Increased ten times (Calder 1990, Lundberg and Halldin 1994)
//        Below canopy aerodynamic resistance 
			Rc =(Hcan*expf(ndecay)/(Kh*ndecay)*(expf(-ndecay*z/Hcan)-expf(-ndecay*(d+Z0c)/Hcan))) *1/Hs_f;
			Rs = (1.0/(powf(K_vc,2.0)*Vz)*powf((logf(z/zo)),2.0))*1/Hs_f;	               
            Rc = Rc+ Rs;                                                     // Total below canopy resistance
//        Bulk aerodynamic resistance (used when there is no snow in the canopy)
		     Rbc = Rc+Ra;                                               // Total aerodynamic resistance from ground to above canopy used when there is no snow on canopy                 
//     Boundary layer resistance of canopy(Rl)
            Dc      = .04;    	                                             // Deckinson et.al (1986);Bonan(1991): Characterstics dimension (m) of leaves (average width) // put in parameters later
            Lbmean  = 0.02/ndecay *sqrt(Vc/Dc)* (1-expf(-ndecay/2));               
	        Rl      = 1/(LAI*Lbmean)*1/Hs_f;
//        Leaf boundary layer conductance: reciprocal of resistance
	        RKINl  = 1.0/Rl;
//      Correction for Ra, Rc and Rac for stable and unstable condition
			Ri = Gra_v*(Ta-Tss)*z /(powf(Vz,2.0)*(0.5*(Ta+Tss)+273.15));    // wind speed unit should not be changed.                                                    
	    	if (Ri > Rimax)
				Ri= Rimax;	    	             
            if (Ri > 0) 
			{
			   Rastar  = Ra;  ///((1-5*Ri)**2)       // Above canopy
			   Rcstar  = Rc /(powf((1-5*Ri),2.0));      // Below canpy
			   Rbcstar = Rbc; //  /((1-5*Ri)**2)     // Bulk
			}                                                       // Rlstar  = Rl/((1-5*Ri)**2)         // leaf boundary layer
			else
			{
				Rastar  = Ra;  ///(1-5*Ri)**0.75 
				Rcstar  = Rc /powf((1-5*Ri),0.75);              
      			Rbcstar = Rbc;  // /(1-5*Ri)**0.75
			}	                                                //	Rlstar  = Rl/(1-5*Ri)**0.75 
//     Turbulent (heat and vapor) transfer coefficient (Conductance) // 
            RKINa   = 1.0/Rastar;              // Above canopy
	 		RKINc   = 1.0/Rcstar;              // Below canpy                  
            RKINbc  = 1.0/Rbcstar;             // Bulk
//	         RKINl   = 1.0/Rlstar              // leaf boundary layer   
		} 
		if(snowdgtvariteflag2 == 1)
	  {
		 std::cout<<"In aeroRes"<<std::endl;
		 std::cout<<std::setprecision(15)<<ndecay<<" "<<Kh<<" "<<Ra<<" "<<Rc<<" "<<Rs<<" "<<Rbc<<" "<<Lbmean<<" "<< Rl<<" "<< RKINc<<" "<< RKINa<<" "<< RKINbc<<" "<< RKINl<<std::endl;
	}

	}           
	return;
}

//********************* WIND PROFILES ****************************************************************
//     Calculates wind at the 2 m above the surface and at the sink using the input of measured wind at 2 m above the canopy
 __host__ __device__   
void  uebCellDA::WINDTRANab(float V, float Zm, float *param, float *sitev,float &Vstar,float &Vh,float &d,float &Z0c,float &Vz,float &Vc)               // Out put wind speed within canopy at height 2 m
{
	float Z = param[4],       //  Nominal meas. height for air temp. and humidity [2m],
      Wcoeff= param[31],       // 

	Cc     = sitev[4],       // Canopy Cover 
	Hcan   = sitev[5],       // Canopy height  
 	LAI    = sitev[6],       // Leaf Area Index 
	Ycage  = sitev[8];     // Parameters for wind speed transformation

    LAI = Cc*LAI;      
//     Roughness length (z0c) and zero plane displacement(d) for canopy [Periera, 1982]
    if(V  <= 0)
	{
	   Vz= 0.0;
	   Vc= 0.0;
	}
	else
	{
		d = Hcan*(0.05 + (powf(LAI,0.20))/2.0 + (Ycage-1)/20.0);                         
		if (LAI < 1) 
			Z0c= 0.1*Hcan;
		else
			Z0c = Hcan*(0.23 - (powf(LAI,0.25))/10.0 - (Ycage-1)/67.0);
		Vstar  = K_vc*V/logf((Zm-d)/Z0c);       //Friction velocity        
	//	Wind speed at height,Z from the ground/snow surface below canopy (Exponential profile from canopy height to below canopy)
		Vh  = 1.0/K_vc*Vstar*logf((Hcan-d)/Z0c);  //Wind speed at the height of canopy (Logarithmic profile above canopy to canopy)
	    Vz   = 1.0* Vh*expf(-Wcoeff*LAI*(1-Z/Hcan));
		Vc   = 1.0* Vh*expf(-Wcoeff*LAI*(1-(d+Z0c)/Hcan));     //To estimate canopy boundary layer conductance
	//	   Vc   = 1.0* Vh*expf(-Wcoeff*LAI*(1-6./Hcan)**Beta)  
	}	                          
	return; 
}

//************************ TURBULENT FLUXES (ABOVE AND BELOW THE CANOPY *******************************
//     Calculates the turbulent heat fluxes (sensible and latent heat fluxes) and condensation/sublimation.
 __host__ __device__   
void  uebCellDA::TURBFLUX (float Ws, float Wc, float A, float Tk,float Tc,float Ta,float Tss, float RH, float V,float Ea,float P,float *param,float *sitev,
                  float &d, float &Z0c, float &Vz, float &Rkinc, float &Rkinbc, float &Tac, float &Fs, float &Ess, float &Esc,  // Output variables				 
                  float &QHc, float &QEc, float &Ec, float &QHs, float &QEs, float &Es, float &QH, float &QE, float &E )   
{
                                                      
    float  RHOA, Eac,  Tak, Tck;  
	float Rc, Ra, Rbc, Rl, Rkina, Rkinl, Eak, Tack, RHOAc;

	float radFracdenom;
    
    float z     =   param[4],      //  Nominal meas. height for air temp. and humidity [2m],
    zo    =   param[5],      //  Surface aerodynamic roughness [m],
    fstab =   param[18],     //  Stability correction control parameter 0 = no corrections, 1 = full corrections

    Apr   =   sitev[1],      //  Atmospheric Pressure [Pa],
    Cc    =   sitev[4],      // Canopy Coverage
	LAI   =   sitev[6];      // Leaf Area Index     
 
    Tak   = Ta + Tk;
    Tck   = Tc + Tk;
	RHOA   = Apr/(Ra_g*Tak);	  // RHOA in kg/m3,  APr: Atm Press
	Ea     = svpw(Ta) * RH;          // Actual vapor pressure sorrounding canopy
	if(snowdgtvariteflag2 == 1)
	{
		  std::cout<<std::setprecision(15)<<Ea<<std::endl;
	}

//     Wind less coefficient:
    float EHoA  = 0.0,       // for sensibleheat flux
	EHoC  = 0.0,
	EHoS  = 0.0,
	EEoA  = 0.0,       // for latent heat flux
    EEoC  = 0.0,
    EEoS  = 0.0;

   AeroRes (P, Wc,V, Ta, Tss, Tc, Fs, param, sitev,Tk, 
            d, Z0c, Vz, Rc, Ra, Rbc, Rl, Rkinc, Rkina, Rkinbc, Rkinl);         // Output variables
    if(snowdgtvariteflag2 == 1)
	  {
		 std::cout<<"Outputs from aeroRes"<<std::endl;
		 std::cout<<std::setprecision(15)<<d<<" "<< Z0c<<" "<< Vz<<" "<< Rc<<" "<< Ra<<" "<< Rbc<<" "<< Rl<<" "<< Rkinc<<" "<< Rkina<<" "<< Rkinbc<<" "<< Rkinl<<std::endl;
	}

	if(snowdgtvariteflag2 == 1)
		std::cout<<"LAI: "<<LAI<<std::endl;
   
  if (V <= 0) 
  {
		QHc  = 0.0;
		QEc  = 0.0;
		Ec   = 0.0;
		QHs  = 0.0;
		QEs  = 0.0;
		Es   = 0.0;
        Tac  = Ta;        // Though We don't need Tac when V=0.
  }
  else
  {
		if (LAI == 0)
		{						        
	//     Flux estimation from the open areas when no vanopy is present
			QHs = (RHOA* C_p * Rkinbc + EHoS) * (Ta-Tss);                        //	0.0 added for zero wind condition  
			QEs = (0.622 * Hne_u/(Ra_g* Tak) * Rkinbc + EEoS) * (Ea - Ess); 
			Es  = -QEs/(Rho_w * Hne_u);                                       // Es in  m/hr
			QHc = 0.0;
			QEc = 0.0;
			Ec  = 0.0;

			if(snowdgtvariteflag2 == 1)
				std::cout<<std::setprecision(15)<<"QHs,QEs,Es,QHc,QEc,Ec"<<" "<<QHs<<" "<<QEs<<" "<<Es<<" "<<QHc<<" "<<QEc<<" "<<Ec<<std::endl;

		}
		else
		{ 
			
			radFracdenom = 1*Rkina + 1*Rkinl + 1*Rkinc;

//###### 9.21.13 to avoid div by 0 when Rkina + 1*Rkinl + 1*Rkinc evaluate to zero 0.01 added
//TBC_9.21.13
			if (radFracdenom == 0)
			{
				radFracdenom += 0.01;
				std::cout<<"Warning! a zero denominator! the sum of turbulent conductances Rkina + 1*Rkinl + 1*Rkinc evaluates to zero"<<std::endl; 
				std::cout<<"added 0.01 to avoid numerical error; Need checking results"<<std::endl;
				//getchar();				
			}
			Tac = (Tc*Rkinl + Tss*Rkinc + Ta*Rkina) / radFracdenom; //(1*Rkina + 1*Rkinl + 1*Rkinc);
			Eac = (Esc*Rkinl + Ess*Rkinc + Ea*Rkina) / radFracdenom; //(1*Rkina + 1*Rkinl + 1*Rkinc);
			Tack    = Tac + Tk;   

//###### 9.22.13 observed Tac values close to -273 which may lead to Tack = 0 
//TBC_9.22.13
			if (Tack == 0)
			{
				Tack += 1.0;
				std::cout<<"Error! Surface temp in function Turbflux() Tack = 0"<<std::endl; 
				std::cout<<"added 1 to avoid numerical error; Need checking results"<<std::endl;
				//getchar();				
			}

    		RHOAc   =  Apr/(Ra_g* Tack);
			//	 Flux from canopy
			QHc  = (RHOAc* C_p *Rkinl + EHoC)*(Tac - Tc); 	

			if(snowdgtvariteflag2 == 1)
				std::cout<<std::setprecision(15)<<"Tc, Tac,Eac,Tack: "<<Tc<<" "<<Tac<<" "<<Eac<<" "<<Tack<<std::endl;
	       
			if ((Wc == 0) && (P == 0))
			{
					QEc = 0.0;                         // No latent heat flux when there is no snow in the caanopy
					Ec  = 0.0;
			}
			else
			{
				QEc  = (0.622*Hne_u/(Ra_g*Tack)*Rkinl + EEoC)*(Eac-Esc)*Fs;
				Ec   = -QEc/(Rho_w*Hne_u);   
			} 			
			//     Flux from snow at ground
		    QHs  =(RHOAc*C_p*Rkinc + EHoS)* (Tac-Tss);                  //	EHoS added for zero wind condition [Jordan et al.(1999); Koivusalo (2002); Herrero et al.(2009) ] 
			QEs  =(0.622*Hne_u/(Ra_g*Tack)*Rkinc+ EEoS)*(Eac-Ess); 
			Es   = -QEs/(Rho_w*Hne_u);
			if(snowdgtvariteflag2 == 1)
				std::cout<<std::setprecision(15)<<"QHs,QEs,Es,QHc,QEc,Ec"<<" "<<QHs<<" "<<QEs<<" "<<Es<<" "<<QHc<<" "<<QEc<<" "<<Ec<<std::endl;
			
		}
	//	Total flux *** moved back 6.25.13
		QH  =  QHs + QHc; 
		QE  =  QEs + QEc; 
		E  =  Es  + Ec;  
  } 
//	Total flux*******moved here 6.6.13=> moved back up 6.25.13
	/*QH  =  QHs + QHc; 
	QE  =  QEs + QEc; 
	E  =  Es  + Ec;  */
    return;
}

//************************** RATE OF INTERCEPTION AND MAX INTERCEPTION ***********************
//     Calculates amount of snow intercepted by the canopy
 __host__ __device__   
void  uebCellDA::INTERCEPT (float Ta, float LAI,float P, float Wc, float dt, float Inmax, float Uc, float Cc,
     	         // Output variables
				 float &ieff, float &Ur, float &intc)
{

	float  Intma;
	//Uc and Cc directly called instead of the param/sitev arrays 6.6.13 
    //Uc   = param(27)		    // Unloading rate coefficient (Per hour) (Hedstrom and pomeroy, 1998) 
    //Cc   = sitev(5)            // Canopy Coverage

	if (LAI == 0)
		Cc = 0.0;
//	Interception rate (di/dt)
    if (P == 0)
	{
	    intc  = 0.0;
        ieff = 0.0;
	}
	else 
	{
        Intma = 0.0;
	    if (Wc < Inmax)
	    {
		  intc  = Cc*(1.0 - Wc *Cc/(Inmax)) * P;
	      ieff = Cc*(1.0 - Wc* Cc/(Inmax));          // Interception efficiency (di/dp) 
	      Intma = Inmax - Wc;
		  if (intc >= Intma)
			intc = Intma;
		}
		else
	    {
	      intc  = 0.0;
          ieff = 0.0;
	    }  
	} 
//	Unloading rate	
	Ur = Uc * (Wc + (intc*dt)/2);       // Melt and Evap. are subtracted 
                                          // half of the current interception is also considered for unloading
	return; 	
}
//******************************   ENERGY ADVECTED BY RAIN   **********************************
//     Calculates the heat advected to the snowpack due to rain
 __host__ __device__   
float uebCellDA::QPF(float PR, float TA, float TO, float PS, float RHOW, float HF, float CW, float CS)
{
     float TRAIN, TSNOW, QPF_v;
	 if(TA > TO)
	  {
         TRAIN = TA;
         TSNOW = TO;	
	  }// TO : Freezing Temp (int.e.0C)
	  else
	  {
         TRAIN = TO;
         TSNOW = TA;
	  }
      
      QPF_v = PR * RHOW *(HF + CW *(TRAIN-TO)) + PS*RHOW*CS*(TSNOW-TO);
      return QPF_v;
}
//****************************** PREHELP () ***************************************************
//      Routine to correct energy and mass fluxes when 
//      numerical overshoots dictate that W was changed in 
//      the calling routine - either because W went negative
//      or due to the liquid fraction being held constant.
 __host__ __device__   
void  uebCellDA::PREHELP( float W1,float W,float DT,float &FM,float FM1,float fac,float PS,float PRAIN,float &E,float RHOW,float HF,float &Q,float &QM,float &MR, float &QE, float HSF)
{
       float QOTHER;
	   FM = (W1-W)/DT*fac-FM1; 
       MR = findMax( 0.0 , (PS + PRAIN - FM - E));
       E = PS + PRAIN - FM - MR;
//    Because melt rate changes the advected energy also changes.  Here
//     advected and melt energy are separated,
       QOTHER = Q + QM - QE;
//     then the part due to melt recalculated
       QM = MR*RHOW*HF;
//     then the part due to evaporation recalculated
       QE = -E*RHOW*HSF;
//     Energy recombined
       Q = QOTHER - QM + QE;
       return;
}

//****************** EXPONENTIAL INTEGRAL FUNCTION *****************************************
//     Computes the expfonential integral function for the given value      
 __host__ __device__  
float uebCellDA::EXPINT(float LAI)
{
	float  a0,a1,a2,a3,a4,a5,b1,b2,b3,b4, EXPINT_v;
	if(LAI == 0 )
		EXPINT_v = 1.0;                   
	else if(LAI <= 1)
	{
		a0= -0.57721566;
		a1= 0.99999193;
		a2= -0.24991055;
		a3= 0.05519968;
		a4= -0.00976004;
		a5= 0.00107857;
		EXPINT_v = a0 + a1*LAI + a2*powf(LAI,2.0) + a3*powf(LAI,3.0) + a4*powf(LAI,4.0) + a5*powf(LAI,5.0) - logf(LAI);
	}
	else
	{
		a1= 8.5733287401;
		a2= 18.0590169730;
		a3= 8.6347637343;
		a4= 0.2677737343;
		b1= 9.5733223454;
		b2= 25.6329561486;
		b3= 21.0996530827;
		b4= 3.9584969228;
		EXPINT_v = (powf(LAI,4.0) + a1*powf(LAI,3.0) + a2*powf(LAI,2.0) + a3*LAI + a4)/
			           	((powf(LAI,4.0) + b1*powf(LAI,3.0) + b2*powf(LAI,2.0) + b3*LAI + b4)*LAI*expf(LAI));
	}
	return EXPINT_v;
}

//      Reference Book:Handbook of mathematical functions with Formulas, Graphs, and Mathematical Tables
//      Edited by Milton Abramowitz and Irene A. Stegun, Volume 55, Issue 1972, page no. 231
//      Dover Publications, Inc, Mineola, New York.
//*************************** SVP () ***********************************************************
//     Calculates the vapour pressure at a specified temperature over water or ice
//     depending upon temperature.  Temperature is celsius here.
 __host__ __device__  
float uebCellDA::svp(float T)
{
    float SVPv;
	if(T >= 0)
        SVPv = svpw(T);
	else
      SVPv = svpi(T);
	return SVPv;
}
//     Calculates the vapour pressure at a specified temperature over water
//     using polynomial from Lowe (1977).
 __host__ __device__  
float uebCellDA::svpw(float T)
{
    float SVPWv;
	 SVPWv = 6.107799961 + T * (0.4436518521 + T * (0.01428945805 + 
		            T * (0.0002650648471 + T * (3.031240936*powf(10.0,-6.0) +
					T * (2.034080948*powf(10.0,-8.0) + T * 6.136820929*powf(10.0,-11.0))))));
	 SVPWv = 100.0*SVPWv;										// convert from mb to Pa
	return SVPWv;
}
//     Calculates the vapour pressure at a specified temperature over ice.
//     using polynomial from Lowe (1977).
 __host__ __device__  
float uebCellDA::svpi(float T)
{
      float SVPIv;
      SVPIv = 6.109177956 + T * (0.503469897 + T * (0.01886013408 +  
		        T * (0.0004176223716 + T * (5.82472028*powf(10.0,-6.0) + 
				T * (4.838803174*powf(10.0,-8.0) + T * 1.838826904*powf(10.0,-10.0))))));
      SVPIv = SVPIv * 100;											// convert from mb to Pa
      return SVPIv;
}
//     Estimates reflection and scattering coefficient
 __host__ __device__  
float uebCellDA::Tau1(float Rho, float G, float h, float COSZEN, float kk)
{
    float Tau1v;
	if(h == 0) 
		Tau1v =1.0;         // To avoid computational error when coszen =0
	else if (COSZEN == 0) 
	     Tau1v = 0;
    else Tau1v  = expf(- kk*G*Rho*h/COSZEN); 
	return Tau1v;
}

 __host__ __device__  
float uebCellDA::Tau2(float Rho, float G, float h, float kk, float EXPI)
{
	float tau2v = (1-kk*G*Rho*h)*expf(-kk*G*Rho*h) + powf((kk*G*Rho*h), 2.0)*EXPI;     // is the penetration function for deep canopy: diffuse
	return tau2v;
}
	 