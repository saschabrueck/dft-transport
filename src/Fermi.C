#include "Fermi.H"

Fermi::Fermi()
{
    Eps0 = 8.854e-12;
    e    = 1.6022e-19;
    kB   = 1.38e-23;
}

/************************************************************************************************/

Fermi::~Fermi()
{
}

/************************************************************************************************/

double Fermi::find_fermi(double Ndop,double sign_rho,double *Ek,double dk,int Nk,double *dky,\
			 int Nky,double *dkz,int Nkz,int n_of_modes,double spin_factor,\
			 double Temp,double A)
{
    int IE;
    int NE      = 125;
    double Emin = get_edge(Ek,sign_rho,Nky*Nkz,Nk,n_of_modes)-sign_rho*1.0;
    double Emax = Emin+sign_rho*2.0;
    double Ef   = Emax;
    double *rho = new double[NE];
  
    density(rho,sign_rho,Ek,dk,Nk,dky,Nky,dkz,Nkz,n_of_modes,Emin,Emax,NE,spin_factor,Temp,A);
    
    for(IE=0;IE<NE;IE++){rho[IE] = log(fabs(rho[IE]-Ndop));}

    for(IE=1;IE<NE-1;IE++){
        if((rho[IE]<rho[IE-1])&&(rho[IE]<rho[IE+1])){
            Ef = Emin+IE*(Emax-Emin)/(NE-1);
        }
    }
 
    Ef = Newton(Ef,Ndop,sign_rho,Ek,dk,Nk,dky,Nky,dkz,Nkz,n_of_modes,spin_factor,Temp,A);
    
    delete[] rho;
    
    return Ef;
}

/************************************************************************************************/

double Fermi::Newton(double Ef0,double Ndop,double sign_rho,double *Ek,double dk,int Nk,\
                     double *dky,int Nky,double *dkz,int Nkz,int n_of_modes,double spin_factor,\
		     double Temp,double A)
{
    int IC           = 0;
    int max_iter     = 20;
    double condition = INF;
    double criterion = 1.0e-4;
    double Ef        = Ef0;
    double n,dn_dEf;

    density(&n,sign_rho,Ek,dk,Nk,dky,Nky,dkz,Nkz,n_of_modes,Ef,Ef,1,spin_factor,Temp,A);

    while((condition>criterion)&&(IC<max_iter)){

        derivate(&dn_dEf,sign_rho,Ek,dk,Nk,dky,Nky,dkz,Nkz,n_of_modes,Ef,Ef,1,spin_factor,\
		 Temp,A);

        Ef = Ef-(n-Ndop)/(dn_dEf);

        density(&n,sign_rho,Ek,dk,Nk,dky,Nky,dkz,Nkz,n_of_modes,Ef,Ef,1,spin_factor,Temp,A);

        condition = fabs(n-Ndop)/Ndop;
        
        IC++;
    }

    return Ef;
    
}

/************************************************************************************************/

double Fermi::get_edge(double *Ek,double sign,int Nkyz,int Nk,int NM)
{
    int I;
    double Emin = INF;
    double result;

    for(I=0;I<Nkyz*Nk*NM;I++){
        if(sign*Ek[I]<Emin){
            Emin   = sign*Ek[I];
            result = Ek[I];
        }
    }

    return result;
}

/************************************************************************************************/

void Fermi::density(double *rho,double sign_rho,double *Ek,double dk,int Nk,double *dky,int Nky,\
		    double *dkz,int Nkz,int n_of_modes,double Efmin,double Efmax,int NE,\
		    double spin_factor,double Temp,double A)
{
    int IE,IM,IK,IKY,IKZ,ind,indk;
    double Ef;

    for(IE=0;IE<NE;IE++){
        
        Ef      = Efmin + IE*(Efmax-Efmin)/(max(NE-1,1));
        rho[IE] = 0.0;
        
	for(IKY=0;IKY<Nky;IKY++){

	    for(IKZ=0;IKZ<Nkz;IKZ++){
            
	        for(IM=0;IM<n_of_modes;IM++){
                
		    for(IK=0;IK<Nk;IK++){
                    
		        ind  = IKY*n_of_modes*Nk*Nkz+IKZ*n_of_modes*Nk+IK*n_of_modes+IM;
			indk = IKY*Nkz+IKZ;
                    
			if(abs(Ek[ind]-Ef)<4.0){

			    if((IK==0)||(IK==Nk-1)){
                        
			        rho[IE] = rho[IE] + (dky[indk]*dkz[indk]*dk/2.0)/\
				  (1.0+exp(e*sign_rho*(Ek[ind]-Ef)/(kB*Temp)));
                        
			    }else{
                        
			        rho[IE] = rho[IE] + dky[indk]*dkz[indk]*dk/ \
				  (1.0+exp(e*sign_rho*(Ek[ind]-Ef)/(kB*Temp)));
                        
			    }
			}
		    }
		}
	    }
	}
        /*
        //2 because 0<k<kmax only instead of -kmax<k<kmax
        rho[IE] = 2.0*spin_factor/(2.0*PI*A*1e-18)*rho[IE];
	*/
	//Now -kmax<k<kmax
	rho[IE] = spin_factor/(2.0*PI*A*1e-18)*rho[IE];
    }
}

/************************************************************************************************/

void Fermi::derivate(double *drho_dEf,double sign_rho,double *Ek,double dk,int Nk,double *dky,\
		     int Nky,double *dkz,int Nkz,int n_of_modes,double Efmin,double Efmax,int NE,\
		     double spin_factor,double Temp,double A)
{
    int IE,IM,IK,IKY,IKZ,ind,indk;
    double Ef;

    for(IE=0;IE<NE;IE++){
        
        Ef      = Efmin + IE*(Efmax-Efmin)/(max(NE-1,1));
        drho_dEf[IE] = 0.0;

	for(IKY=0;IKY<Nky;IKY++){

	    for(IKZ=0;IKZ<Nkz;IKZ++){
            
	        for(IM=0;IM<n_of_modes;IM++){
                
		    for(IK=0;IK<Nk;IK++){

		        ind  = IKY*n_of_modes*Nk*Nkz+IKZ*n_of_modes*Nk+IK*n_of_modes+IM;
			indk = IKY*Nkz+IKZ;

			if(abs(Ek[ind]-Ef)<4.0){

			    if((IK==0)||(IK==Nk-1)){
                        
			        drho_dEf[IE] = drho_dEf[IE] + (dky[indk]*dkz[indk]*dk/2.0)* \
				  exp(e*sign_rho*(Ek[ind]-Ef)/(kB*Temp))/ \
				  pow(1.0+exp(e*sign_rho*(Ek[ind]-Ef)/(kB*Temp)),2.0);
                        
			    }else{
                        
			        drho_dEf[IE] = drho_dEf[IE] + dky[indk]*dkz[indk]*dk* \
				  exp(e*sign_rho*(Ek[ind]-Ef)/(kB*Temp))/ \
				  pow(1.0+exp(e*sign_rho*(Ek[ind]-Ef)/(kB*Temp)),2.0);
                        
			    }
			}
		    }
		}
	    }
	}
        /*
        //2 because 0<k<kmax only instead of -kmax<k<kmax
        drho_dEf[IE] = 2.0*e*sign_rho*spin_factor/(2.0*PI*A*1e-18*kB*Temp)*drho_dEf[IE];
	*/
        //Now -kmax<k<kmax
	drho_dEf[IE] = e*sign_rho*spin_factor/(2.0*PI*A*1e-18*kB*Temp)*drho_dEf[IE];
    }
}

/************************************************************************************************/
