---
title: UnDiFi-2D
logo: ./logo.png
subtitle: User Guide
version: Version 1.4
author:
  - Dipartimento di Ingegneria Meccanica e Aerospaziale, (DIMA), Sapienza University, Rome, Italy
  - Scuola di Ingegneria, Università degli Studi della Basilicata, Potenza, Italy
  - Team CARDAMOM, Inria Bordeaux Sud-Ouest, Talence, France
  - Saint-Petersburg State University, 7/9 Universitetskaya nab., St. Petersburg, Russia
institute:
date: \today
documentclass: scrreprt
lang: en-US
papersize: a4
fontsize: 11pt
geometry: "left=2.0cm,right=2.0cm,top=3.5cm,bottom=2.5cm"
colorlinks: yes
toc: yes
header-includes:
  - \input{header}
bibliography: ./refs.bib
csl: ./ieee.csl
link-citations: true
---

\hypertarget{introduction}{}

# Introduction

 [**UnDiFi-2D**](https://github.com/UnDiFi-2D/UnDiFi-2D) is an open source
 Unstructured Discontinuity Fitting algorithm code.
 It is written in in standard (compliant) Fortran 77/95 with highly
 modularity as design target The aim of UnDiFi-2D is to explicitely
 manage discontinuities in the flow field.

 In our unstructured shock-fitting approach the shock front is
 described using a double-sided, polygonal curve. Two sets of flow
 states, corresponding to the upstream and downstream sides of the
 discontinuity, are assigned to the grid-points located on either side
 of the shock front. This is allowed to move, while obeying to the
 Rankine-Hugoniot jump relations, throughout a background triangular
 mesh that covers the entire computational domain. At each time
 step, a local, constrained Delaunay triangulation is applied in the
 neighbourhood of the shock front to ensure that the edges that make up
 the shock front are part of the overall triangular grid. The fitted
 shock acts as an interior boundary for the shock-capturing solver that
 is used to solve the discretised governing equations in the smooth
 regions of the flow-field.

 A typical shock-fitting calculation starts with the initial shock
 position guessed using a previous shock-capturing calculation The
 simulation is then advanced in time until the shock front has reached
 its steady state position, corresponding to vanishing shock speed.

 Three different *in-house* gas-dynamic solvers, have been used so far:

 * **EulFS**, developed by Aldo Bonfiglioli [@bonfiglioli2000fluctuation;@bonfiglioli2013mass]
 * **NEO**, developed by Mario Ricchiuto at INRIA [@ricchiuto2010explicit;@ricchiuto2015explicit;@arpaia2015ale]
 * **COOLFluiD**, developed by Andrea Lani at VKI [@lani2005coolfluid;@lani2016]

 All three codes share the same vertex-centered, Fluctuation Splitting
 [@abgrall2006residual;@deconinck1993multidimensional] discretization
 of the governing PDEs on linear triangles and tetrahedra. However,
 virtually any vertex-centered FE or FV solver featuring a linear
 representation of the dependent variables can be used as gas-dynamic
 solver within the shock-fitting procedure. Only EulFS and NEO are
 included in the present repository and documentation.

 Three key software components (or modules) can be identified in the
 unstructured discontinuity-fitting approach:

 * Algorithm (described in the paper) <!--(Chapter \ref{chap:algorithm})-->
 * Solver (Chapter \ref{chap:solver})
 * Meshing (Chapter \ref{chap:mesh})

 Modularity stems from the fact that these three different components
 communicate through ad-hoc interfaces. This programming approach may
 not be the most efficient from the standpoint of computational speed,
 because, for instance, one has to switch, at each time step, among
 the different data-structures used by the three different modules.
 However, the approach is very convenient, since it allows us to use
 "off the shelf" gas-dynamics solvers and mesh generation tools that
 are treated as black boxes and can be replaced by similar ones only
 by changing the interfaces, with a modest coding effort. In the
 following, we shall describe each module.

## How this documentation is organized

This user guide is organized to both guide the first steps as well as
provide a complete overview of the simulation code's features from a
user and a developer point of view.

* Chapter \ref{chap:installation} contains step by step instructions
from obtaining the source code, installation, up to running a first
simulation and visualizing the simulation results. In addition,
it provides an overview of the whole simulation framework and the
currently implemented features.
* Chapter \ref{chap:workflow} outlines the workflow starting with mesh
generation and concluding with the visualization of results produced
with UnDiFi-2D.
 <!--* Chapter \ref{chap:algorithm} describes the discontinuity fitting
algorithm implemented in*UnDiFi-2D.-->
* Chapter \ref{chap:solver} presents main features of the gasdynamic
solvers currently coupled with UnDiFi-2D.
* Chapter \ref{chap:mesh} gives and overview of the meshing software
packages used in UnDiFi-2D.
<!--* Chapter \ref{chap:features_models} shall serve as a reference for the models and features implemented in **UnDiFi-2D**.-->
<!--* Chapter \ref{chap:visu_output} presents the options and parameters for the output data, field and flow variables.-->
<!-- * Chapter \ref{chap:tools} lists tools within the **UnDiFi-2D** repository, including the post-processing tools.-->
* Simulation tutorials are contained in Chapter \ref{chap:tutorials}.
<!-- * Cluster-specific user guidelines are given in Chapter \ref{chap:cluster_guide}. -->
<!-- * A complete list of all parameters is given in Chapter \ref{chap:parameterfile}. -->
<!-- * The unit test system used to test key routines with CTest is described in Chapter \ref{chap:unittest}. -->
