/*
 * SolveurExpl.cpp : Application schemas semi-implicite sur les objets.
 * Copyright (C) 2016 Florence Zara, LIRIS
 *               florence.zara@liris.univ-lyon1.fr
 *               http://liris.cnrs.fr/florence.zara/
 *
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/** \file Calculs.cpp
 Fonctions de calculs communes aux objets simules.
 \brief Fonctions de calculs communes aux objets simules.
 */

#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>

#include "vec.h"
#include "ObjetSimule.h"
#include "Viewer.h"
#include "SolveurExpl.h"


using namespace std;



/**
 * Calcul de l acceleration des particules
 * avec ajout de la gravite aux forces des particules
 * et ajout de la force due au vent sur une des particules du maillage
 * et reinitialisation des forces.
 */
void SolveurExpl::CalculAccel_ForceGravite(Vector g,
                                           int nb_som,
                                           std::vector<Vector> &A,
                                           std::vector<Vector> &Force,
                                           std::vector<float> &M)
{
    for(int i=0; i<nb_som; i++)
    {
        if(M[i] == 0.0)
            M[i] = 1;
        Force[i] = Force[i] + M[i]*g;
        A[i] = Force[i]/M[i];
    }

}//void


/*! Calcul des vitesses et positions :
 *  Formule d Euler semi-implicite :
 *  x'(t+dt) = x'(t) + dt x"(t)
 *  x(t+dt) = x(t) + dt x'(t+dt)
 */
void SolveurExpl::Solve(float visco,
                        int nb_som,
                        int Tps,
                        std::vector<Vector> &A,
                        std::vector<Vector> &V,
                        std::vector<Vector> &P)
{

     for(int i=0; i<nb_som; i++)
    {
        V[i] = V[i] + A[i] * _delta_t;
        P[i] = P[i] + V[i] * _delta_t;
    }

}//void

