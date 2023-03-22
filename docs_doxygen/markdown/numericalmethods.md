@page numericalmethods Numerical methods
@tableofcontents

The simulations of the model described before rely on the division of the problem into two separate routines -- one for the hyperbolic part (convective and conserved terms) and a second one for the parabolic (viscous and source terms). Viewing the model as having the form @f$\partial_t\vec{u}=\mathscr{L}_\text{conservation}(\vec{u})+\mathscr{L}_\text{dissipation}(\vec{u})@f$, so that an operator splitting method could be used, uncoupling the two terms, for instance, for the mass flux equation one has:
@f[
\frac{\partial \vec{p}}{\partial t}+\overbrace{\bm{\nabla}\!\cdot\!\left( \frac{ \vec{p}\otimes \vec{p}}{{n}^{3/2}} + v_F^2\frac{n^{3/2}}{3}\mathds{1}+S^2\frac{{n}^2}{2}\mathds{1}\right)}^{\mathscr{L}_\text{conservation}(\vec{p})}+\\
-\underbrace{\nu_s\nabla^2\frac{\vec{p}}{n^{3/2}} -\nu_o\nabla^2\frac{\vec{p}^\dagger}{n^{3/2}}
+\omega_c\frac{\vec{p}\times\vec{\hat{z}}}{\sqrt{n}}}_{\mathscr{L}_\text{dissipation}(\vec{p})}=0.
@f]

Then, if the solution, after a time @f$t@f$, for the conservation and dissipation operators is well approximated by the numerical evolution operators @f$\Lambda^t_\text{con}@f$ and @f$\Lambda^t_\text{dis}@f$, respectively, the solution of the complete system is given by the Godunov splitting. This method is generally first order in time. One can, of course, invert the order of the operators and in fact the splitting error depends on the commutator @f$[\mathscr{L}_\text{dis},\mathscr{L}_\text{con}]@f$. For second order accuracy one can use Strang splitting where @f$\vec{u}^{t_0+t}\approx\Lambda^{t/2}_\text{dis}\circ\Lambda^t_\text{con}\circ\Lambda^{t/2}_\text{dis}(\vec{u}^{t_0})@f$. However, the performed test so far have not find significant accuracy or stability improvement when using the more computation time costly Strang splitting; so, we recommend applying the simple Godunov algorithm:
@f[
\vec{u}^{t_0+t}\approx\Lambda^t_\text{dis}\circ\Lambda^t_\text{con}(\vec{u}^{t_0}).
@f]

The hyperbolic operator was solved with a finite volume Lax-Wendroff type method \cite Hirsch2007 \cite LeVeque1992  --  the two-step Richtmyer scheme for nonlinear systems \cite LeVeque1992 .  Such scheme was then elected to implement, given it is computationally light, as it is not required to explicitly evaluate the system Jacobian at each step, yet capable of returning accurate simulations. It consists in two distinct steps, a predictor where the mid nodes are calculated and then a corrector step that updates the simulated quantities at the central nodes.

Whereas the dissipation and source terms were then computed with a forward time centred space (FTCS) scheme of finite differences. This choice of a lightweight parabolic solver has impact on the viable set of parameters for simulation, nonetheless, it allows to keep the maximum computing time for a complete 2D simulation bellow 150 min.

@section onedrich One-dimensional Richtmyer method

In the onedimensional case, the hyperbolic system of equations for density and velocity fields is written in a conservation form:
@f[
\frac{\partial}{\partial t}\begin{bmatrix}n\\v\end{bmatrix}+\frac{\partial}{\partial x}
\begin{bmatrix}      nv\\ \frac{v^2}{4}+\frac{v_F^2}{2}\log n + 2S^2\sqrt{n}\\\end{bmatrix} = 0 \iff \frac{\partial}{\partial t}\vec{u}+\frac{\partial}{\partial x}\mathbf{F(u)}=0, \label{eq:conserv}
@f]

and was implemented for @f$n@f$ and @f$v@f$ according to

@f{align}
\vec{u}_{i+1/2}^{k+1/2} &= \frac{1}{2}\left(\vec{u}_{i+1}^k + \vec{u}_{i}^k\right) - \frac{\Delta t}{2\Delta x}\left( \mathbf{F}_{i+1}^k - \mathbf{F}_{i}^k\right)\label{eq:richt_1}\\
\vec{u}_i^{k+1} &= \vec{u}_i^k - \frac{\Delta t}{\Delta x} \left( \mathbf{F}_{i+1/2}^{k+1/2} - \mathbf{F}_{i-1/2}^{k+1/2}\right)\label{eq:richt_2}
@f}
using @f$\vec{u}@f$ and @f$\mathbf{F}(\vec{u})@f$ as defined by equation and space (indicated at subscript indices) and time (indicated at superscript indices) discretisation where @f$t=k\Delta t@f$ and @f$t=i\Delta x@f$. Since such scheme is second order both in time and space it will not introduce spurious diffusion in the solution as a first order scheme would. It will, however, introduce artificial oscillations, contiguous to discontinuities, which can be corrected by the application of a moving average smoothing filter.

To comply with Courant--Friedrichs--Lewy condition, the time step is chosen to be

@f[
\Delta t= \frac{\Delta x}{\lambda}  \quad\text{with}\quad
\lambda=\begin{cases}
1.2\,v_F&\text{if }S<0.36\,v_F\\
1.97\,S+v_F/2&\text{otherwise }
\end{cases}\label{eq:cfl}
@f]

such expression for the information speed @f$\lambda@f$ was empirically obtained for a wide variety of space discretizations and parameters.

@section twodrich Two-dimensional Richtmyer method


Regarding the two-dimensional case, one needs to consider the flux function along @f$y@f$
@f[
\frac{\partial}{\partial t}\vec{u}+ \frac{\partial}{\partial x}\vec{F}(\vec{u}) +\frac{\partial}{\partial y}\vec{G}(\vec{u})=0   
@f]
where now the state vector is @f$\vec{u}=[n,p_x,p_y]^\mathsf{T}@f$ and the fluxes functions are given by
@f[
\vec{F}(\vec{u})=\begin{bmatrix}      
p_xn^{-1/2}\\ \frac{p_x^2}{n^{3/2}}+\frac{v_F^2}{3}n^{3/2} + \frac{S^2}{2}n^2\\
p_xp_yn^{-3/2}
\end{bmatrix}
@f]
and
@f[
\vec{G}(\vec{u})=\begin{bmatrix}      
p_yn^{-1/2}\\
p_xp_yn^{-3/2} \\
\frac{p_y^2}{n^{3/2}}+\frac{v_F^2}{3}n^{3/2} + \frac{S^2}{2}n^2
\end{bmatrix}.
@f]

The generalisation of leads to a more complex method having, for the first step, the form
@f[
\vec{u}^{k+1/2}_{i+1/2,j+1/2}=\frac{1}{4}\left(\vec{u}^{k}_{i,j}+\vec{u}^{k}_{i+1,j}+\vec{u}^{k}_{i,j+1}}+\vec{u}^{k}_{i+1,j+1}\right)
-\frac{\Delta t}{2\Delta x}\left(\vec{F}^{k}_{i+1,j+1/2}-\vec{F}^{k}_{i,j+1/2} \right) -\frac{\Delta t}{2\Delta y}\left( \vec{G}^{k}_{i+1/2,j+1}-\vec{G}^{k}_{i+1/2, j}\right)
@f]
and then, the corrector is obtained by
@f[
\vec{u}^{k+1}_{i,j}=\vec{u}^{k}_{i,j}-\frac{\Delta t}{\Delta x}\left(\vec{F}^{k+1/2}_{i+1/2,j}-\vec{F}^{k+1/2}_{i-1/2,j}\right)-\frac{\Delta t}{\Delta y}\left(\vec{G}^{k+1/2}_{i,j+1/2}}-\vec{G}^{k+1/2}_{i,j-1/2}\right)
@f]
In the matter of this scheme stability, for equal spacing @f$\Delta y=\Delta x@f$, the previous condition, defined for the one-dimensional case, has hitherto been sufficient to guarantee stability.

@section twodftcs19 Two-dimensional weighted FTCS method


Here we briefly present the weighted FTCS method \cite Hayman1988, which is used to solve the parabolic part of the differential equations. Hence, we aim to solve an equation of the type
@f[
\frac{\partial \vec{u}}{\partial t}=\vec{\alpha}\cdot \nabla^2\vec{u},
@f]
with @f$vec{\alpha}= \left[0, \nu_s, \nu_s, \alpha\right]^\top@f$ a vector of diffusion coefficients. It must be noted that, for this method, the mass flux is converted to velocity before applying the algorithm, so @f$\vec{u}=\left[n, v_x, v_y, T\right]^\top@f$.
Instead of simply applying the central difference operator at each grid point, a stencil of 9 points (the central point and its 1<sup>st</sup> and 2<sup>nd</sup> neighbours) is used to calculate the Laplacian. The contribution of each term is taken into account using a weight parameter @f$\theta@f$ such that
@f[
\frac{\partial^2 u}{\partial x^2}\approx (1-2\theta)\delta^2_xu_{i,j}+\theta\left(\delta^2_xu_{i,j-1}+\delta^2_xu_{i,j+1}\right)
@f]
and
@f[
\frac{\partial^2 u}{\partial y^2}\approx (1-2\theta)\delta^2_yu_{i,j}+\theta\left(\delta^2_yu_{i-1,j}+\delta^2_yu_{i+1,j}\right),
@f]
where @f$\delta^2_x@f$ and @f$\delta^2_y@f$ are the central difference operators along @f$x@f$ and @f$y@f$, respectively. It is seen that the weights that minimize the second-order errors, following the approach of _modified equivalent partial differential equation_ \cite Warming1974, lead to the following expression for the time evolution \cite Hayman1988 :
@f[
\vec{u}^{k+1}_{i,j}=\Delta_x\Delta_y\left(\vec{u}^{k}_{i-1,j-1}+\vec{u}^{k}_{i-1,j+1}+\vec{u}^{k}_{i+1,j-1}+\vec{u}^{k}_{i+1,j+1}\right)+\Delta_y(1-2\Delta_x)\left(\vec{u}^{k}_{i,j-1}+\vec{u}^{k}_{i,j+1}\right)
+\Delta_x(1-2\Delta_y)\left(\vec{u}^{k}_{i-1,j} + \vec{u}^{k}_{i+1,j}\right)+(1-2\Delta_x)(1-2\Delta_y)\vec{u}^{k}_{i,j}
@f]
with @f$\Delta_x=\frac{\alpha \Delta t}{\Delta x^2}@f$ and  @f$\Delta_y=\frac{\alpha \Delta t}{\Delta y^2}@f$. 