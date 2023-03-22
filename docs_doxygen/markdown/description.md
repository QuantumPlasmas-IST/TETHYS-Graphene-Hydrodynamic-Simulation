@page  description Software Description
@tableofcontents


@section purpose Purpose and scope

The main objective of this code is to simulate and analyse a hydrodynamic description of electronic conduction on (mono-layer) graphene. At the present time limiting the simulations (we refer to one spatial dimension simulation when the geometry and symmetry of the problem let one foreseen that the dynamics of the system will develop solely along the direction of the main electronic current i.e. along the GFET channel and that the quantities and derivatives along the transverse direction of the plane can be discarded accordingly. If, on the contrary, the fluid is expected to display a non-parallel flow or transverse variations a complete 2D simulation is required.), both 1D+1 and 2D+1, to the case of gated graphene system. The physical description of such problem can be found at \cite Cosme2020 and \cite Cosme2021 . In the bidimensional simulation is possible to add the effects of shear viscosity and/or uniform and static magnetic field along the transverse direction, we include also, the implementation of odd viscosity \cite Avron1998OddViscosity ,  commonly dubbed as *Hall viscosity*. 


The equations governing this problem can be nondimensionalised  which is done resorting to the channel length @f$L@f$, equilibrium number density @f$n_0@f$, velocity @f$v_0@f$ and the Fermi temperature @f$T_F=E_F/k_B=\hbar v_F\sqrt{\pi n_0}/k_B@f$, with @f$v_F=1\times10^6\,{\rm ms^{-1}}@f$ the Fermi velocity. Thus, writing

@f[ 
x^\ast\equiv x/L\quad t^\ast\equiv tv_0/L\quad \omega^\ast\equiv \omega L/v_0
@f]
@f[
    v^\ast\equiv v/v_0\quad n^\ast\equiv n/n_0 \quad T^\ast \equiv T/T_F
@f]
all quantities are henceforth are implied to be dimensionless, even though the asterisk sign  will be dropped for simplicity. 

The complete model that the code aims to simulate is given by the following equations, already written in the form of conservation laws;

@f[
  \frac{\partial n}{\partial t}+\bm{\nabla} \cdot \frac{\vec{p}}{\sqrt{n}}=0,
@f]
@f[
     \frac{\partial \vec{p}}{\partial t}+\bm{\nabla} \cdot \left( \frac{\vec{p}\otimes \vec{p}}{n^{3/2}}+\frac{v_{F}^{2}}{3}n^{3/2}+\frac{S^{2}}{2}n^{2}\right)
     -\nu_{s}\nabla ^{2} \frac{\vec{p}}{n^{3/2}}-\nu _{o}\nabla ^{2}\frac{\vec{p}^{\dagger}}{n^{3/2}}=-\omega _{c}\frac{\vec{p}^\dagger}{\sqrt{n}}-\frac{\vec{p}}{\tau}
 @f] and 
 @f[
    \frac{\partial T}{\partial t}+\bm{\nabla} \cdot \left[\left( \frac{T}{n^{\sfrac{3}{2}}}+\frac{3}{4\pi ^{2}}\right) \vec{p}\right]-\alpha \nabla ^{2}T=\frac{S^{2}}{v_{F}^{2}}\frac{\vec{p}}{\sqrt{n}}\cdot \bm{\nabla} n+\Pi,%
@f]%
for the number density @f$n@f$, momentum density @f$\vec{p}=m^\star n \vec{v}@f$ (with @f$m^\star@f$ the effective mass of carriers and @f$\vec{v}@f$ the fluid velocity) and temperature @f$T@f$, and where @f$\nu_s@f$ is the kinematic shear viscosity, @f$\nu_o@f$ the odd viscosity \cite Avron1998OddViscosity \cite Narozhny2019MagnetohydrodynamicsViscosities \cite Pellegrino2017NonlocalLiquids , @f$\vec{p}^\dagger=\left[-p_y,p_x\right]^\top@f$,  @f$\omega_c@f$ the cyclotron frequency, @f$\tau@f$ the typical time for inelastic collisions, @f$\alpha@f$ the thermal diffusivity and @f$\Pi@f$ encompasses any possible source terms for the heat flux. The parameter @f$S@f$ can be interpreted as the sound velocity for the plasmons \cite Cosme2021 \cite Cosme2020  and for a GFET with a capacitance per area @f$C_g=\varepsilon/d@f$ it is given by @f$S^2=e^2dv_F\sqrt{ n_0}/ \varepsilon\hbar\sqrt{\pi}@f$, with @f$d@f$ the distance between the gate and the graphene layer, @f$\varepsilon@f$ the medium permittivity, and @f$n_0@f$ the steady state or equilibrium number density.  Finally, @f$m^\star = \hbar \sqrt{\pi n}/v_F@f$ is the effective mass of the carriers, obtained resorting to a Drude model \cite Chaves2017 \cite Cosme2020. The system of equations presented is implemented by separating the hyperbolic, parabolic and source terms in the equation and applying appropriate numerical methods for each one. 

The employed numerical methods are conditionally stable,  the simulation routines are robust for a wide range of the input parameters, provided that the validity of model assumptions still hold true. That is, highly degenerate electronic fluid where @f$T_F\gg T@f$; being, consequently,  away from the charge neutrality point so that the carriers are only of one kind, either electrons or holes. Moreover, in the case of simulations with an applied magnetic field, it must not be high enough to induce Landau levels in the system. To ease the of compliance with such requirements, we suggest in the tables below the range of the physical parameters of the system and their typical dimensionless values for the simulations. 


### Input parameters typical ranges

| Parameter   |      Designation    |   Dimensionless input | 
| :---:   |      :---      |   :---: |
|@f$S@f$ | sound velocity | @f$> 10@f$|
|@f$v_F@f$ | Fermi velocity|@f$> 10@f$|
|@f$\nu_s@f$ | shear viscosity| @f$0.1\sim 1@f$|
|@f$\nu_o@f$ | odd viscosity | @f$<0.5@f$|
|@f$\omega_c@f$ | cyclotron frequency  | @f$< 10@f$|
|@f$1/\tau@f$ | collision frequency& | @f$\ll 1@f$ |
|@f$\alpha@f$ | thermal diffusivity | @f$\sim 1@f$ |


### Physical parameters

| Parameter   |      Description    |   Relations  | Requirements | Typical values | 
| :---:   |      :---      |  --- | :---:| ---  |
| @f$ v_F @f$   |  Fermi velocity @f${}^{(a)}@f$   |  |   |  @f$1\times 10^6 ms^{-1}@f$  |
| @f$E_F@f$   |  Fermi level    | @f$= \hbar v_F\sqrt{\pi n}@f$  |  @f$ k_BT\ll E_F \ll t \simeq 3 eV@f$ |  @f$0.0658 \sqrt{n} eV \dagger@f$ |
| @f$\mu@f$  |   Carrier mobility  |   |   | @f$2000~5000 cm^2/Vs@f$   |
| @f$\tau@f$  |   Relaxation time  | @f$= \lambda_i/v_F@f$  |   | @f$\sim 5ps   (0.1 \sim 10ps)@f$ |
| @f$1/\tau@f$  |   Collision frequency@f${}^{(a)}@f$  |   |   | @f$\sim0.2THz @f$ |
| @f$\lambda_i@f$  |  Mean free path (impurities)     | |  |  @f$\sim5\mu m@f$ |
| @f$\lambda_{ee}@f$  | Mean free path (electrons)       ||   | @f$\sim10nm@f$  |
| @f$ v_{sat} @f$  |   Saturation velocity  |     @f$\sim0.5v_F@f$| | @f$0.5\times10^6 ms^{-1}@f$   |
| @f$\nu_s@f$  |   Shear kinematic viscosity@f${}^{(a)}@f$     ||   |  @f$\sim 0.2  m^2s^{-1}@f$  |
| @f$\nu_o@f$  |   Odd kinematic viscosity     || broken time symmetry (@f$B\neq0@f$)  | @f$ \lesssim 0.1  m^2s^{-1}@f$  |
| @f$C_g@f$ | Capacitance per area   |  @f$= \varepsilon/d_0@f$|  | @f$0.1~0.6 \mu F cm^{-2}@f$  |
| @f$L@f$ |   Channel length@f${}^{(b)}@f$  |  | @f$\lambda_{ee} < L <  \lambda_i @f$   | |
| @f$W@f$ |  Channel width   |  | @f$\lambda_{ee} < W < \lambda_i  @f$  | |
| @f$x:y@f$ |  Aspect ratio@f${}^{(a)}@f$   | = @f$L/W@f$ |    | |
| @f$v_0@f$ |   Mean drift velocity@f${}^{(b)}@f$  |  |  @f$\ll v_F and < v_sat@f$| @f$~1x10^5 ms^{-1}@f$|
| @f$n_0@f$ |   Mean carrier density@f${}^{(b)}@f$  |  | | @f$~ 1\times10^{12} cm^{-2}@f$ |
| @f$m^\star@f$ | Effective mass  | @f$= \hbar\sqrt{\pi n}/v_F@f$ | | @f$\sim0.02m_e@f$ |
| @f$S@f$ |  Plasmons velocity@f${}^{(a)}@f$   | @f$= ev_F n_0^{1/4}/\hbar C_g^{1/2}\pi^{1/2}@f$    |  @f$> v_F@f$ |@f$1.17\times10^6 n_0^{1/4}C_g^{-1/2}  ms^{-1}\dagger @f$|
| @f$Re@f$  |    Reynolds number    | @f$= v_0L/\nu_s@f$   |   | |
| @f$t_0@f$ |  Time scale@f${}^{(b)}@f$  | = @f$L/v_0@f$   |   | |
| @f$\omega_0@f$ |  Frequency scale@f${}^{(b)}@f$ |@f$ = v_0/L@f$   |   |  |
| @f$\omega_c@f$ |  Cyclotron frequency@f${}^{(a)}@f$  |@f$ =eB/m^\star@f$   |  @f$B\ll1 T @f$| @f$\ll 9 THz@f$ |
| @f$\omega_{DS}@f$ |   Dyakonov-Shur frequency  |  @f$\sim \pi S/2L@f$  |   |  |
| @f$\gamma_{DS}@f$ |  Dyakonov-Shur growth rate   | @f$\sim0.75 v_0/L@f$ | > @f$1/\tau@f$  |  |
| @f$U_0@f$ |   Gate mean voltage  | @f$= en_0/C_g@f$  | | @f$160 n_0/C_g mV\dagger@f$ |
| @f$I_0@f$ |   Typical channel current  | @f$= en_0v_0W@f$ | @f$<1.3 mA/cm for L=25 \mu m@f$  | @f$160 n_0v_0W \mu A \ddagger @f$ |

<p>@f${}^{(a)}@f$ User defined parameter</p>
<p>@f${}^{(b)}@f$ Normalisation parameter</p>
<p>@f$\dagger@f$ If the density  @f$n_0@f$ is given in @f$10^12 cm^{-2}@f$  and the capacitance @f$C_g@f$ in @f$\mu Fcm^{-2}@f$ </p>
<p>@f$\ddagger@f$ If the density  @f$n_0 is@f$ given in @f$10^12 cm^{-2}@f$, the drift velocity @f$v_0@f$ in @f$10^5 ms^{-1}@f$ and the channel width W in @f$\mu m@f$ </p>



