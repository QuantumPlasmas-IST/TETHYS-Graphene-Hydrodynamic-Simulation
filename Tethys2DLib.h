#ifndef TETHYSLIB_2D_H
#define TETHYSLIB_2D_H

#include <H5Cpp.h>
#include "TethysLib.h"
using namespace H5;

class Fluid2D : public TETHYSBase
{
	protected:	
		float dx=1.0;
		float dy=1.0;		
		float dt=1.0;		
		const float lengX=1.0;
		const float lengY=1.0;
		float vel_snd =50.0;
		float kin_vis =0.0;
		float * vel_snd_arr;	
		float * den_mid ; // 1st Aux. Grid (Nx-1)*(Ny-1)
		float * flxX_mid ;
		float * flxY_mid ;

        float * lap_flxX ; //new grids for the laplacians
        float * lap_flxY ;

		std::ofstream data_preview;
		virtual void SetFileName();

	public :	
		float * den ;
		float * velX ;
		float * velY ;
		float * flxX ;
		float * flxY ;
		float * curX ;
		float * curY ;
		explicit Fluid2D(int sizeNx, int sizeNy,float VELSND, float VISCO);
		~Fluid2D();
		void SetVelSnd(float x);
		void SetSound();	
		float GetVelSnd();
		void SetKinVis(float x);
		float GetKinVis();
		float GetDx();
		void SetDx(float x);
		float GetDy();
		void SetDy(float x);
		float GetDt();
		void SetDt(float x);

    virtual void SetSimulationTime();

		void InitialCondRand();
		void InitialCondTEST();
		void Richtmyer();	
		//void DimensionalSplittingMethod();
		virtual void CFLCondition();

		virtual float DensityFluxX(float n, float velX, float velY,float mass, float S);
		virtual float DensityFluxY(float n, float velX, float velY,float mass, float S);
		virtual float MassFluxXFluxX(float n, float flxX, float flxY,float mass, float S);
		virtual float MassFluxXFluxY(float n, float flxX, float flxY,float mass, float S);
		virtual float MassFluxYFluxX(float n, float flxX, float flxY,float mass, float S);
		virtual float MassFluxYFluxY(float n, float flxX, float flxY,float mass, float S);
		virtual float DensitySource(float n, float flxX, float flxY, float S);
		virtual float MassFluxXSource(float n, float flxX, float flxY, float S);
		virtual float MassFluxYSource(float n, float flxX, float flxY, float S);
		virtual void MassFluxToVelocity();
		
		void CreateFluidFile();
		void WriteFluidFile(float t) ;
};

class GrapheneFluid2D : public Fluid2D{
	protected : 
		float vel_fer =10.0;							
		float col_freq =0.0; 						
	public : 
		using Fluid2D::Fluid2D;
		
		GrapheneFluid2D(int sizeNx,int sizeNy,float VELSND, float FERMI,float VISCO,float COL);

		
		void SetVelFer(float x);
		float GetVelFer();
		void SetColFreq(float x);
		float GetColFreq();
		void CFLCondition() override;
	    
	    void SetSimulationTime() override ;
	    
	    void MassFluxToVelocity() override;	    
		float DensityFluxX(float n, float flxX, float flxY,float mass, float S) override;
		float DensityFluxY(float n, float flxX, float flxY,float mass, float S) override;
		float MassFluxXFluxX(float n, float flxX, float flxY,float mass, float S) override;
		float MassFluxXFluxY(float n, float flxX, float flxY,float mass, float S) override;
		float MassFluxYFluxX(float n, float flxX, float flxY,float mass, float S) override;
		float MassFluxYFluxY(float n, float flxX, float flxY,float mass, float S) override;
		float DensitySource(float n, float flxX, float flxY, float S) override;
		float MassFluxXSource(float n, float flxX, float flxY, float S)  override;
		float MassFluxYSource(float n, float flxX, float flxY, float S) override;
		
		void MagneticSource();
        void SourceFTCS();
        void ViscosityFTCS();
		void WriteAtributes();
};
#endif
