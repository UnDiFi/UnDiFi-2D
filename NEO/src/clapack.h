/***************************************************************************
                                  clapack.h
                                  ---------
        Here the functions of CLAPACK used in the code are reported
                             -------------------
    begin                : Wed May 15 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
 ***************************************************************************/

extern int dgetrf_(integer *, integer *, doublereal *, integer *, integer *, integer * ) ;
extern int dgetri_(integer *, doublereal *, integer *, integer *, doublereal *, integer *, integer * ) ;

/***************************************************************************
 *                                                                         *
 *   This program has been developed at the von Karman Institute for Fluid *
 *   Dynamics by Mario Ricchiuto. Any change or update must be done in a   *
 *   clean way and it must be clearly commented.                           *
 *   Mario Ricchiuto                                                       *
 *                                              23/04/2002                 *
 ***************************************************************************/
