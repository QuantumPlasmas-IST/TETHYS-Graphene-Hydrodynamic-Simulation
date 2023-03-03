/************************************************************************************************\
* 2020 Pedro Cosme , João Santos, Ivan Figueiredom, João Rebelo, Diogo Simões                    *
* DOI: 10.5281/zenodo.4319281																	 *
* Distributed under the MIT License (license terms are at http://opensource.org/licenses/MIT).   *
\************************************************************************************************/

#include "includes/BoundaryLib.h"
#include "includes/RobinBoundaryLib.h"

void RobinBoundaryCondition::SlipLength(Fluid2D &fluid_class, float slip_length) {
	RobinBoundaryCondition::SlipLengthBottom(fluid_class, slip_length);
	RobinBoundaryCondition::SlipLengthBottom(fluid_class, slip_length);
}

void RobinBoundaryCondition::SlipLengthTop(Fluid2D &fluid_class, float slip_length) {
	int nx=fluid_class.SizeX();
	int ny=fluid_class.SizeY();
	float dy=fluid_class.GetDy();
	for (int i=0; i < nx; i++){
		float aux_top, dn_top,l_top;
		int top= i + (ny - 1) * nx;
/*		dn_top = (-3.0f*fluid_class.Den[i+(ny - 1)*nx]+4.0f*fluid_class.Den[i+(ny - 2)*nx]-1.0f*fluid_class.Den[i+(ny - 3)*nx])/(2.0f*dy);
		dn_top = dn_top/sqrt(fluid_class.Den[top]);
		l_top = slip_length/(1.0f+slip_length*dn_top);
		aux_top = l_top/(2.0f*dy+3.0f*l_top);
		fluid_class.Den[top] = fluid_class.Den[top - nx];
		fluid_class.FlxY[top] =  0.0f;
		fluid_class.FlxX[top] = aux_top*(4.0f*fluid_class.FlxX[i+(ny-2)*nx]-1.0f*fluid_class.FlxX[i+(ny-3)*nx]);
*/
		dn_top = (-3.0f*fluid_class.Umain[top].n()+4.0f*fluid_class.Umain[top-nx].n()-1.0f*fluid_class.Umain[top-2*nx].n())/(2.0f*dy);
		dn_top = dn_top/sqrt(fluid_class.Umain[top].n());
		l_top = slip_length/(1.0f+slip_length*dn_top);
		aux_top = l_top/(2.0f*dy+3.0f*l_top);
		fluid_class.Umain[top].n() = fluid_class.Den[top - nx];
		fluid_class.Umain[top].py() =  0.0f;
		fluid_class.Umain[top].px() = aux_top*(4.0f*fluid_class.Umain[top-nx].px()-1.0f*fluid_class.Umain[top-2*nx].px());
	}
}

void RobinBoundaryCondition::SlipLengthBottom(Fluid2D &fluid_class, float slip_length) {
	int nx=fluid_class.SizeX();
	float dy=fluid_class.GetDy();
	for (int i=0; i < nx; i++){
		float aux_bot, dn_bot, l_bot;
		dn_bot = (-3.0f*fluid_class.Umain[i+0*nx].n()+4.0f*fluid_class.Umain[i+1*nx].n()-1.0f*fluid_class.Umain[i+2*nx].n())/(2.0f*dy);
		dn_bot = dn_bot/sqrt(fluid_class.Umain[i].n());
		l_bot = slip_length/(1.0f-slip_length*dn_bot);
		aux_bot = l_bot/(2.0f*dy+3.0f*l_bot);
		fluid_class.Umain[i].n() = fluid_class.Umain[i+nx].n();
		fluid_class.Umain[i].py() = 0.0f;
		fluid_class.Umain[i].px() = aux_bot*(4.0f*fluid_class.Umain[i+1*nx].px()-1.0f*fluid_class.Umain[i+2*nx].px());
	}
}