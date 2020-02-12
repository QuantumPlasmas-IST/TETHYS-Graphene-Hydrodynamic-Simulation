#include "LaxWen2D.h"

LaxWen2D::LaxWen2D()
{
	dx = .01;
	Nx = 101;
	dy = .01;
	Ny = 101;
	dt = .001;
	bv = new bool[12];
	bv_mat = new float*[12];
	for(int i=0; i<6; i++)
	{
		bv[i]=false;
		bv_mat[i] = new float[Ny];
		for(int j=0; j<Ny; j++)
			bv_mat[i][j] = 1;
	}
	for(int i=6; i<12; i++)
	{
		bv[i] = false;
		bv_mat[i] = new float[Nx];
		for(int j=0; j<Nx; j++)
			bv_mat[i][j] = 1;
	}
	mcx = false;
	mcy = false; 
}


LaxWen2D::LaxWen2D(float dx_, int Nx_, float dy_, int Ny_, float dt_, bool* bv_)
{
	dx = dx_;
	Nx = Nx_;
	dy = dy_;
	Ny = Ny_;
	dt = dt_;
	bv = new bool[12];
	bv_mat = new float*[12];
	if(bv_ == NULL)
	{ 
		for(int i=0; i<6; i++)
		{
			bv[i] = false;
			bv_mat[i] = new float[Ny];
			for(int j=0; j<Ny; j++)
				bv_mat[i][j] = 1;
		}
		for(int i=6; i<12; i++)
		{
			bv[i] = false;
			bv_mat[i] = new float[Nx];
			for(int j=0; j<Nx; j++)
				bv_mat[i][j] = 1;
		}
	}
	else
	{
		for(int i=0; i<6; i++)
		{
			bv[i] = bv_[i];
			bv_mat[i] = new float[Ny];
			for(int j=0; j<Ny; j++)
				bv_mat[i][j] = 1;
		}
		for(int i=6; i<12; i++)
		{
			bv[i] = bv_[i];
			bv_mat[i] = new float[Nx];
			for(int j=0; j<Nx; j++)
				bv_mat[i][j] = 1;
		}
	}
	mcx = false;
	mcy = false;
}

LaxWen2D::~LaxWen2D()
{

	for(int i=0; i<12; i++)
		delete[] bv_mat[i];
	delete[] bv_mat;
	delete[] bv;
}

void LaxWen2D::SetGrid(float dx_, float Nx_, float dy_, float Ny_)
{
	dx = dx_;
	Nx = Nx_;
	dy = dy_;
	Ny = Ny_;
}

void LaxWen2D::SetTimeStep(float dt_)
{
	dt = dt_;
}

void LaxWen2D::SetBV(bool* bv_)
{
	for(int i=0; i<12; i++)
		bv[i] = bv_[i];
}

void LaxWen2D::SetBV(short unsigned cases, float* var_bv)
{
	switch(cases)
	{
		case 0:
			for (int i=0; i<Ny; i++)
				bv_mat[0][i] = var_bv[i];
			break;
		case 1:
			for (int i=0; i<Ny; i++)
				bv_mat[1][i] = var_bv[i];
			break;
		case 2:
			for (int i=0; i<Ny; i++)
				bv_mat[2][i] = var_bv[i];
			break;
		case 10:
			for (int i=0; i<Ny; i++)
				bv_mat[3][i] = var_bv[i];
			break;
		case 11:
			for (int i=0; i<Ny; i++)
				bv_mat[4][i] = var_bv[i];
			break;
		case 12:
			for (int i=0; i<Ny; i++)
				bv_mat[5][i] = var_bv[i];
			break;
		case 20:
			for (int i=0; i<Nx; i++)
				bv_mat[6][i] = var_bv[i];
			break;
		case 21:
			for (int i=0; i<Nx; i++)
				bv_mat[7][i] = var_bv[i];
			break;
		case 22:
			for (int i=0; i<Nx; i++)
				bv_mat[8][i] = var_bv[i];
			break;
		case 30:
			for (int i=0; i<Nx; i++)
				bv_mat[9][i] = var_bv[i];
			break;
		case 31:
			for (int i=0; i<Nx; i++)
				bv_mat[10][i] = var_bv[i];
			break;
		case 32:
			for (int i=0; i<Nx; i++)
				bv_mat[11][i] = var_bv[i];
			break;
		default:
			cout << "\nValue of int out of range, please consult .h file and try again \n";
	}
}

void LaxWen2D::SetMCbool(bool mcx_, bool mcy_)
{
	mcx = mcx_;
	mcy = mcy_;
}

void LaxWen2D::Rich2(float*** var_cur, void FluxX(float* var_in, float* var_out), void FluxY(float* var_in, float* var_out))
{
	// Intermediate variables required
	//cout << "oi\n";
	float*** var_med; // aux grid (Nx-1 * Ny-1)
	var_med = new float**[Nx-1];
	for(int i=0; i<Nx-1; i++)
	{
		var_med[i] = new float*[Ny-1];
		for (int j=0; j<Ny-1; j++)
			var_med[i][j] = new float[3];
	}

	float*** var_med_x; // 2nd aux grid for x
	float*** flx_temp_y; // saves FluxX values 
	var_med_x = new float**[Nx-1];
	flx_temp_y = new float**[Nx-1];
	for (int i=0; i<Nx-1; i++)
	{
		var_med_x[i] = new float*[Ny];
		flx_temp_y[i] = new float*[Ny];
		for(int j=0; j<Ny; j++)
		{
			var_med_x[i][j] = new float[3];
			flx_temp_y[i][j] = new float[3];
			for(int k=0; k<3; k++)
				var_med_x[i][j][k] = 0.5*(var_cur[i][j][k] + var_cur[i+1][j][k]);
			FluxY(var_med_x[i][j], flx_temp_y[i][j]);
		}
	}

	float*** var_med_y; // 2nd aux grid for y
	float*** flx_temp_x; // saves FluxY values
	var_med_y = new float**[Nx];
	flx_temp_x = new float**[Nx];
	for (int i=0; i<Nx; i++)
	{
		var_med_y[i] = new float*[Ny-1];
		flx_temp_x[i] = new float*[Ny-1];
		for(int j=0; j<Ny-1; j++)
		{
			var_med_y[i][j] = new float[3];
			flx_temp_x[i][j] = new float[3];
			for(int k=0; k<3; k++)
				var_med_y[i][j][k] = 0.5*(var_cur[i][j][k] + var_cur[i][j+1][k]);
			FluxX(var_med_y[i][j], flx_temp_x[i][j]);
		}
	}
	// Step 1 implementation ( evaluate var_med )
	for(int i=0; i<Nx-1; i++)
	{
		for(int j=0; j<Ny-1; j++)
		{

			for(int k=0; k<3; k++)
			{
				var_med[i][j][k] = .25*(var_cur[i][j][k] + var_cur[i+1][j][k] + var_cur[i][j+1][k] + var_cur[i+1][j+1][k]);
				var_med[i][j][k] -= 0.5*dt*(flx_temp_x[i+1][j][k]-flx_temp_x[i][j][k])/dx;
				var_med[i][j][k] -= 0.5*dt*(flx_temp_y[i][j+1][k]-flx_temp_y[i][j][k])/dy;

			}
		}
		//cout << var_med[i][100][1]<< " ";
	}
	//cout << endl;

	// Step 2 implementation 

	// Change flx_temp_x
	for(int i=0; i<Nx-2; i++)
	{
		for(int j=0; j<Ny-1; j++)
		{
			for (int k=0; k<3; k++)
				var_med_x[i][j][k] = 0.5*(var_med[i][j][k] + var_med[i+1][j][k]);
			FluxY(var_med_x[i][j], flx_temp_y[i][j]);
		}
		//cout << flx_temp_y[i][100][1] << " ";
	}
	//cout << endl << "flx_temp_x:" << endl;
	// Change flx_temp_y
	for(int i=0; i<Nx-1; i++)
	{
		for(int j=0; j<Ny-2; j++)
		{
			for (int k=0; k<3; k++)
				var_med_y[i][j][k] = 0.5*(var_med[i][j][k] + var_med[i][j+1][k]);
			FluxX(var_med_y[i][j], flx_temp_x[i][j]);
		}
		//cout << flx_temp_x[i][100][1]<< " ";
	}
	//cout << endl;
	// Updating var_cur
	for(int i=0; i<Nx-2; i++)
	{
		for (int j = 0; j<Ny-2; j++)
		{
			for(int k=0; k<3; k++)
			{
				var_cur[i+1][j+1][k] -= dt*(flx_temp_x[i+1][j][k] - flx_temp_x[i][j][k])/dx;
				var_cur[i+1][j+1][k] -= dt*(flx_temp_y[i][j+1][k] - flx_temp_y[i][j][k])/dy;
			}
		}
	}

	// implementing boundary conditions  (Nx-2),(Ny -2) points
	for(int i=0; i<3; i++)
	{
		if(bv[i])
		{
			for(int j=1; j<Ny-1; j++)
				var_cur[0][j][i%3] = var_cur[1][j][i%3];
		}
		else
		{
			for(int j=1; j<Ny-1; j++)
				var_cur[0][j][i%3] = bv_mat[i][j];
		}

	}
	for(int i=3; i<6; i++)
	{
		if(bv[i])
		{
			for(int j=1; j<Ny-1; j++)
				var_cur[Nx-1][j][i%3] = var_cur[Nx-2][j][i%3];
		}
		else
		{
			for(int j=1; j<Ny-1; j++)
				var_cur[Nx-1][j][i%3] = bv_mat[i][j];
		}

	}
	for(int i=6; i<9; i++)
	{
		if(bv[i])
		{
			for(int j=1; j<Nx-1; j++)
				var_cur[j][0][i%3] = var_cur[j][1][i%3];
		}
		else
		{
			for(int j=1; j<Nx-1; j++)
				var_cur[j][0][i%3] = bv_mat[i][j];
		}

	}
	for(int i=9; i<12; i++)
	{
		if(bv[i])
		{
			for(int j=1; j<Nx-1; j++)
				var_cur[j][Ny-1][i%3] = var_cur[j][Ny-2][i%3];
		}
		else
		{
			for(int j=1; j<Nx-1; j++)
				var_cur[j][Ny-1][i%3] = bv_mat[i][j];
		}

	}

	// Boundary Conditions - Corners
	// (0,0)
	for(int k=0; k<3; k++)
	{
		if (bv[k] && bv[k+6])
			var_cur[0][0][k] = var_cur[1][1][k];
		else if (!bv[k] && bv[k+6])
			var_cur[0][0][k] = bv_mat[k][0];
		else if (bv[k] && !bv[k+6])
			var_cur[0][0][k] = bv_mat[k+6][0];
		else if (bv_mat[k][0] * bv_mat[k+6][0] < .000001)
			var_cur[0][0][k] = 0.;
		else
			var_cur[0][0][k] = 0.5*(bv_mat[k][0] + bv_mat[k+6][0]);
	}
	// (Nx,0)
	for(int k=0; k<3; k++)
	{
		if (bv[k+3] && bv[k+6])
			var_cur[Nx-1][0][k] = var_cur[Nx-2][1][k];
		else if (!bv[k+3] && bv[k+6])
			var_cur[Nx-1][0][k] = bv_mat[k+3][0];
		else if (bv[k+3] && !bv[k+6])
			var_cur[Nx-1][0][k] = bv_mat[k+6][Nx-1];
		else if (bv_mat[k+3][0] * bv_mat[k+6][Nx-1] < .000001)
			var_cur[Nx-1][0][k] = 0.;
		else
			var_cur[Nx-1][0][k] = 0.5*(bv_mat[k+3][0] + bv_mat[k+6][Nx-1]);
	}
	// (0,Ny)
	for(int k=0; k<3; k++)
	{
		if (bv[k] && bv[k+9])
			var_cur[0][Ny-1][k] = var_cur[1][Ny-2][k];
		else if (!bv[k] && bv[k+9])
			var_cur[0][Ny-1][k] = bv_mat[k][Ny-1];
		else if (bv[k] && !bv[k+9])
			var_cur[0][Ny-1][k] = bv_mat[k+9][0];
		else if (bv_mat[k][Ny-1] * bv_mat[k+9][0] < .000001)
			var_cur[0][Ny-1][k] = 0.;
		else
			var_cur[0][Ny-1][k] = 0.5*(bv_mat[k][Ny-1] + bv_mat[k+9][0]);
	}
	// (Nx,Ny)
	for(int k=0; k<3; k++)
	{
		if (bv[k+3] && bv[k+9])
			var_cur[Nx-1][Ny-1][k] = var_cur[Nx-2][1][k];
		else if (!bv[k+3] && bv[k+9])
			var_cur[Nx-1][Ny-1][k] = bv_mat[k+3][Ny-1];
		else if (bv[k+3] && !bv[k+9])
			var_cur[Nx-1][Ny-1][k] = bv_mat[k+9][Nx-1];
		else if (bv_mat[k+3][0] * bv_mat[k+9][Nx-1] < .000001)
			var_cur[Nx-1][Ny-1][k] = 0.;
		else
			var_cur[Nx-1][Ny-1][k] = 0.5*(bv_mat[k+3][Ny-1] + bv_mat[k+9][Nx-1]);
	}

	// Deletes
	for(int i=0; i<Nx-1; i++)
	{
		for(int j=0; j<Ny-1; j++)
		{
			delete[]var_med[i][j];
		}
		delete[]var_med[i];
	}
	delete[]var_med;

	for(int i=0; i<Nx-1; i++)
	{
		for(int j=0; j<Ny; j++)
		{
			delete[]var_med_x[i][j];
			delete[]flx_temp_y[i][j];
		}
		delete[]var_med_x[i];
		delete[]flx_temp_y[i];
	}
	delete[]var_med_x;
	delete[]flx_temp_y;

	for(int i=0; i<Nx; i++)
	{
		for(int j=0; j<Ny-1; j++)
		{
			delete[]var_med_y[i][j];
			delete[]flx_temp_x[i][j];
		}
		delete[]var_med_y[i];
		delete[]flx_temp_x[i];
	}
	delete[]var_med_y;
	delete[]flx_temp_x;
}

void LaxWen2D::DumpVariables()
{
	cout 
	<< "\n\nCalled function DumpVariables().\nDumping...\n"
	<< "dx = " << dx << endl
	<< "Nx = " << Nx << endl
	<< "dy = " << dy << endl
	<< "Ny = " << Ny << endl
	<< "dt = " << dt << endl;
	for (int i = 0; i < 12; ++i)
	{
		cout << "bv[" << i << "] = ";
		if(bv[i])
			cout << "true -> Neumann" << endl;
		else 
			cout << "false -> Dirichelet" << endl;
	}
}

void LaxWen2D::DumpBV(unsigned short cases)
{
	switch(cases)
	{
		case 0:
			for (int i=0; i<Ny; i++)
				cout << endl << bv_mat[0][i];
			break;
		case 1:
			for (int i=0; i<Ny; i++)
				cout << endl << bv_mat[1][i];
			break;
		case 2:
			for (int i=0; i<Ny; i++)
				cout << endl << bv_mat[2][i];
			break;
		case 10:
			for (int i=0; i<Ny; i++)
				cout << endl << bv_mat[3][i];
			break;
		case 11:
			for (int i=0; i<Ny; i++)
				cout << endl << bv_mat[4][i];
			break;
		case 12:
			for (int i=0; i<Ny; i++)
				cout << endl << bv_mat[5][i];
			break;
		case 20:
			for (int i=0; i<Nx; i++)
				cout << endl << bv_mat[6][i];
			break;
		case 21:
			for (int i=0; i<Nx; i++)
				cout << endl << bv_mat[7][i];
			break;
		case 22:
			for (int i=0; i<Nx; i++)
				cout << endl << bv_mat[8][i];
			break;
		case 30:
			for (int i=0; i<Nx; i++)
				cout << endl << bv_mat[9][i];
			break;
		case 31:
			for (int i=0; i<Nx; i++)
				cout << endl << bv_mat[10][i];
			break;
		case 32:
			for (int i=0; i<Nx; i++)
				cout << endl << bv_mat[11][i];
			break;
		default:
			cout << "\nValue of int out of range, please consult .h file and try again";
	}
	cout << endl;
}
/*void LaxWen2D::MacCor(float*** var_cur, float* FluxX(float* var), float* FluxY(float* var), bool fix_dir)
{

}*/
