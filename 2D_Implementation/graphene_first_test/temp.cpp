
void GrapheneFluid2D::BoundaryCond(int type){
	//       X        Y 
	//1:       DS + open
	//2:       DS + periodic
	//3:     open + open
	//4:     open + periodic
	//5: periodic + open
	//6: periodic + periodic
	switch(type){
		case 1 : 		
			for(int j=0;j<Ny;j++){	
				den[0+j*Nx]=1.0;               			//constant density at x=0
				den[Nx-1+j*Nx]=den[Nx-2+j*Nx]; 			//free density at x=L
				flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
				flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
				flxX[0+j*Nx] = flxX[1+j*Nx]*pow(den[1+j*Nx],-1.5);			//free flux at x=0
				flxX[Nx-1+j*Nx] = sqrt(den[Nx-1+j*Nx]);	//constant current at x=L (flux equals mass)
			}
			for (int i=0; i<Nx; i++){
				den[i+0*Nx] = den[i+1*Nx];
				flxX[i+0*Nx] = flxX[i+1*Nx];
				flxY[i+0*Nx] = flxY[i+1*Nx];
				den[i+(Ny-1)*Nx] = den[i+(Ny-2)*Nx];
				flxX[i+(Ny-1)*Nx] = flxX[i+(Ny-2)*Nx];
				flxY[i+(Ny-1)*Nx] = flxY[i+(Ny-2)*Nx];
			}	 	
			break;
		case 2 :
			for(int j=0;j<Ny;j++){	
				den[0+j*Nx]=1.0;               			//constant density at x=0
				den[Nx-1+j*Nx]=den[Nx-2+j*Nx]; 			//free density at x=L
				flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
				flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
				flxX[0+j*Nx] = flxX[1+j*Nx]*pow(den[1+j*Nx],-1.5);			//free flux at x=0
				flxX[Nx-1+j*Nx] = sqrt(den[Nx-1+j*Nx]);	//constant current at x=L (flux equals mass)
			}
			for (int i=0; i<Nx; i++){
				den[i+0*Nx] = den[i+(Ny-2)*Nx];
				flxX[i+0*Nx] = flxX[i+(Ny-2)*Nx];
				flxY[i+0*Nx] = flxY[i+(Ny-2)*Nx];
				den[i+(Ny-1)*Nx] = den[i+1*Nx];
				flxX[i+(Ny-1)*Nx] = flxX[i+1*Nx];
				flxY[i+(Ny-1)*Nx] = flxY[i+1*Nx];
			}
			break;
		case 3 :
			for(int j=0;j<Ny;j++){	
				den[0+j*Nx]=den[1+j*Nx];                	
				den[Nx-1+j*Nx]=den[Nx-2+j*Nx]; 			//free density at x=L
				flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
				flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
				flxX[0+j*Nx] = flxX[1+j*Nx]*pow(den[1+j*Nx],-1.5);			//free flux at x=0
				flxX[Nx-1+j*Nx] =  flxX[Nx-2+j*Nx];
			}		
			for (int i=0; i<Nx; i++){
				den[i+0*Nx] = den[i+1*Nx];
				flxX[i+0*Nx] = flxX[i+1*Nx];
				flxY[i+0*Nx] = flxY[i+1*Nx];
				den[i+(Ny-1)*Nx] = den[i+(Ny-2)*Nx];
				flxX[i+(Ny-1)*Nx] = flxX[i+(Ny-2)*Nx];
				flxY[i+(Ny-1)*Nx] = flxY[i+(Ny-2)*Nx];
			}	 			
			break;
		case 4 :
			for(int j=0;j<Ny;j++){	
				den[0+j*Nx]=den[1+j*Nx];                	
				den[Nx-1+j*Nx]=den[Nx-2+j*Nx]; 			//free density at x=L
				flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
				flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
				flxX[0+j*Nx] = flxX[1+j*Nx]*pow(den[1+j*Nx],-1.5);			//free flux at x=0
				flxX[Nx-1+j*Nx] =  flxX[Nx-2+j*Nx];
			}		
			for (int i=0; i<Nx; i++){
				den[i+0*Nx] = den[i+(Ny-2)*Nx];
				flxX[i+0*Nx] = flxX[i+(Ny-2)*Nx];
				flxY[i+0*Nx] = flxY[i+(Ny-2)*Nx];
				den[i+(Ny-1)*Nx] = den[i+1*Nx];
				flxX[i+(Ny-1)*Nx] = flxX[i+1*Nx];
				flxY[i+(Ny-1)*Nx] = flxY[i+1*Nx];
			}		
			break;
		case 5 :
			for(int j=0;j<Ny;j++){	
				den[0+j*Nx]=den[(Nx-2)+j*Nx];                	
				den[Nx-1+j*Nx]=den[1+j*Nx]; 			
				flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
				flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
				flxX[0+j*Nx] = flxX[(Ny-2)+j*Nx];
				flxX[Nx-1+j*Nx] =  flxX[1+j*Nx];
			}		
			for (int i=0; i<Nx; i++){
				den[i+0*Nx] = den[i+1*Nx];
				flxX[i+0*Nx] = flxX[i+1*Nx];
				flxY[i+0*Nx] = flxY[i+1*Nx];
				den[i+(Ny-1)*Nx] = den[i+(Ny-2)*Nx];
				flxX[i+(Ny-1)*Nx] = flxX[i+(Ny-2)*Nx];
				flxY[i+(Ny-1)*Nx] = flxY[i+(Ny-2)*Nx];
			}	 	
			break;
		case 6 :
			for(int j=0;j<Ny;j++){	
				den[0+j*Nx]=den[(Nx-2)+j*Nx];                	
				den[Nx-1+j*Nx]=den[1+j*Nx]; 			
				flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
				flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
				flxX[0+j*Nx] = flxX[(Ny-2)+j*Nx];
				flxX[Nx-1+j*Nx] =  flxX[1+j*Nx];
			}		
			for (int i=0; i<Nx; i++){
				den[i+0*Nx] = den[i+(Ny-2)*Nx];
				flxX[i+0*Nx] = flxX[i+(Ny-2)*Nx];
				flxY[i+0*Nx] = flxY[i+(Ny-2)*Nx];
				den[i+(Ny-1)*Nx] = den[i+1*Nx];
				flxX[i+(Ny-1)*Nx] = flxX[i+1*Nx];
				flxY[i+(Ny-1)*Nx] = flxY[i+1*Nx];
			}
			break;
		default : 			
			for(int j=0;j<Ny;j++){	
				den[0+j*Nx]=1.0;               			//constant density at x=0
				den[Nx-1+j*Nx]=den[Nx-2+j*Nx]; 			//free density at x=L
				flxY[0+j*Nx] = 0.0; 					//flux only on x at x=0
				flxY[Nx-1+j*Nx] = 0.0 ;					//idem at x=L
				flxX[0+j*Nx] = flxX[1+j*Nx]*pow(den[1+j*Nx],-1.5);			//free flux at x=0
				flxX[Nx-1+j*Nx] = sqrt(den[Nx-1+j*Nx]);	//constant current at x=L (flux equals mass)
			}
			for (int i=0; i<Nx; i++){
				den[i+0*Nx] = den[i+1*Nx];
				flxX[i+0*Nx] = flxX[i+1*Nx];
				flxY[i+0*Nx] = flxY[i+1*Nx];
				den[i+(Ny-1)*Nx] = den[i+(Ny-2)*Nx];
				flxX[i+(Ny-1)*Nx] = flxX[i+(Ny-2)*Nx];
				flxY[i+(Ny-1)*Nx] = flxY[i+(Ny-2)*Nx];
			}	 	
	}
}
