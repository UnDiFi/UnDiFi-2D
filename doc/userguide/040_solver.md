\hypertarget{solver}{}

# Gasdynamic Solvers \label{chap:solver}
In this chapter, the gasdynamic solvers are coupled
with UnDiFi, are detailed.

Both solvers, EulFS and NEO, relies on *Fluctuation Splitting* (or
*Residual Distribution*) schemes, meaning that the way of discretizing
the equations is based on the same approach, though there are few
features that differentiate the two codes that will be addressed in the
following sub-sections.

Fluctuation Splitting, first introduced in the late 80s by Roe
[@roe1987linear] to study scalar convection problems, have emerged
as an alternative to Finite Volumes and Finite Elements methods. The
introduction of wave modeling, in the same years [@roe1986discrete],
allowed to adapt the approach to hyperbolic systems, and in particular
the compressible Euler equations. Indeed, it consists in decomposing
the system into the superposition of contributions from simple waves.
Although originally conceived for solving the compressible Euler
equations, this class of schemes is applicable to those systems whose
inviscid terms are of hyperbolic nature.

In general, the computational domain $\Omega\in \Re^d$ is tessellated
into triangles in the 2D space and tetrahedra in 3D. A dual tesselation
is also defined, which consists in the medial dual cells, obtained
by joining the centroids of gravity of all cells surrounding
a given grid-point. Both tesselation are shown in Fig. [RD]:
$C_i$ is the median dual cell centered about grid-point $i$ and
$\Omega_e$ corresponds to triangle $e$. To simplify the process,
the discretization of the inviscid fluxes will be described for the
following scalar conservation law:

\begin{equation}\label{solver.eq.1}
\int_{C_i} \frac{\partial u}{\partial t}dV = \int_{\partial C_i} \vec{F} \cdot d\vec{n} \text{    with    } \vec{F} = \vec{a}u
\end{equation}

The surface integral in Eq. \ref{solver.eq.1} is written as a weighted sum of the flux integrals of the elements sharing the node:

\begin{equation}\label{solver.eq.2}
\int_{\partial C_i} \vec{F} \cdot d\vec{n} = \sum_{e\ni i} \beta_i^T \int_{\partial e} \vec{F} \cdot d\vec{n}
\end{equation}

Introducing the inflow parameters:

\begin{equation}\label{solver.eq.3}
k_i^e = \frac{1}{d}(\vec{a}\vec{n_i^e})
\end{equation}

where $\vec{n_i^e}$ denotes the normal of the face opposite to vertex $i$ of element $e$, scaled by its measure.
Therefore, the fluctuation $\phi^e$ reads:

\begin{equation}\label{solver.eq.4}
\phi^e = \int_{\partial e} \vec{F} \cdot d\vec{n} = \int_{\partial e} \vec{a}u\cdot d\vec{n} = - \sum_{j\in e}k_j^e u_j
\end{equation}

where the summation in Eq.\ref{solver.eq.4} ranges over the $d+1$ vertices of element $e$. Replacing Eq. \ref{solver.eq.4}
in Eq. \ref{solver.eq.2}, the right hand side of Eq. \ref{solver.eq.1} can be written in the follwing form:

\begin{equation}\label{solver.eq.5}
\int_{\partial C_i} \vec{F} \cdot d\vec{n} = \sum_{e\ni i}\beta_i^e \phi^e
\end{equation}

Thus, it is the choice of the distribution coefficients $\beta_i^e$
that determines the properties of the discrete solution. Several
criteria have been proposed and used in the construction of the
discrete schemes such as:

\begin{itemize}
\item \textit{Upwinding} An upwind scheme distributes fractions of the cell fluctuation only among its downstream vertices;
\item \textit{Positivity} The positivity criterion ensures that the maximum principle holds also at the discrete level;
\item \textit{Linearity preservation} Linearity preservation is an accuracy requirement and refers to the ability of the discrete
scheme to reproduce exactly a linear polynomial, steady solution of Eq. \ref{solver.eq.1}.
\end{itemize}

![Residual Distribution concept](./images/RD.png){ width=65% }

Both gasdynamic solvers, EULFS and NEO have the purpose of updating
the solution at time level $t+\Delta t$. In particular, they use the
computational mesh generated as input that is cut in non-communicating
parts by the shock which is treated by the code as if it were an
internal boundary.

## EulFS
The EULFS code is an in-house, unstructured CFD solver that has
been developed over the last 20 years by Aldo Bonfiglioli (see
[@bonfiglioli2000fluctuation] for a detailed description). This
tool is able to work in both two and three space dimensions
and stores the solution at the vertices of triangles, in
2-D, and tetrahedra, in 3-D. In both cases, the solution
is assumed to vary linearly and continuously in space. The
inviscid cell fluctuation $\phi^e$ is evaluated over each
triangular/tetrahedral element $e$ by means of a conservative
linearization [Deconinck,Roe,Struijs,1993] based on the parameter
vector $Z = (\sqrt{\rho},\sqrt{\rho}H,\sqrt{\rho}u,\sqrt{\rho}v)^T$
and scattered to the element vertices using signals $\phi_i^e$ (see
Fig. [RD].a). Within a cell $e$, the signals have to sum up to the net
flux for conservation, $\sum_{i\in e} \phi_i^e=\phi^e$. The different
Fluctuation Splitting schemes proposed in the literature differ by the
way cell residuals are split into signals. The schemes that may used
are several and based on the different features that they present.
Starting from a monotonicity preserving but first-order-accurate
scheme, named N scheme, and a second-order accurate, which may
lead to unphysical oscillations in the neighbourhood of a captured
discontinuity, called LDA scheme; EULFS is also provided with a
non-linear scheme, which captures discontinuities monotically and
preserves second order of accuracy in smooth regions, that blends the
linear N and LDA schemes in such a way that the former is activated
only where discontinuities occur. The blend is based on a smoothness
sensor that makes the new scheme non-linear.

## NEO
The NEO code has been developed by Mario Ricchiuto
[@ricchiuto2010explicit] and has been mainly used to study
time-dependent problems. It is based on a different formulation
aimed at designing explicit Runge-Kutta residual distribution
schemes exhaustively described in [@ricchiuto2010explicit]. This
explicit approach is based on three main ingredients: first recast
the RD discretization as a stabilized Galerkin scheme, then use
a shifted time discretization in the stabilization operator, and
lastly apply high order mass lumping on the Galerkin component of
the discretization. In particular, this approach turned out to be
very useful in simulating unsteady flows by coupling NEO with the
presented shock-fitting technique. Some of the results obtained from
this work have been published in the chapter [@campoli2017shock],
whose contributions showed very promising results in the development
of the unsteady shock-fitting version. The computations shown in the
chapter show the possibility of using the aforementioned linear first-
(monotone) and second-order N and LDA schemes, which are based on
a multidimensional upwind distribution of the cell residual, their
non-linear blend (B scheme), and two non-upwind methods. In particular,
the explicit predictor-corrector formulation of the second order linear
Streamline-Upwind (SU) method proposed in [@ricchiuto2010explicit]
and the nonlinear blended central (Bc) discretization obtained when
blending the SU method with a limited Lax-Friedrich's distribution.
