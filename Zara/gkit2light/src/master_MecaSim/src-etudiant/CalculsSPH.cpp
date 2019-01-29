/*
 * CalculsSPH.cpp :
 * Copyright (C) 2017 Florence Zara, LIRIS
 *               florence.zara@liris.univ-lyon1.fr
 *               http://liris.cnrs.fr/florence.zara/
 *
 * Utilisation du code :
 * https://www.cs.cornell.edu/~bindel/class/cs5220-f11/code/sph.pdf
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


/** \file CalculsSPH.cpp
Programme calculant pour  un fluide son etat au pas de temps suivant (methode d 'Euler semi-implicite) en utilisant la methode SPH (Smoothed Particles Hydrodynamics):
principales fonctions de calculs.
\brief Fonctions de calculs pour un fluide avec methode SPH.
*/

#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>

#include "vec.h"
#include "ObjetSimule.h"
#include "ObjetSimuleSPH.h"
#include "Viewer.h"

using namespace std;



/**
 * Calcul des densites des particules.
 * Formule :
 *  \rho_i = \frac{4m}{\pi h^8} \sum_{j \in N_i} (h^2 - r^2)^3.
 */
void ObjetSimuleSPH::CalculDensite()
{

    for(int i=0; i< _Nb_Sommets; i++)
    {
        float d = 0;
        for(int j=0; j< _Nb_Sommets; j++)
        {

            float val = h*h - length2(P[i] - P[j]);
            if( val > 0)
                d += pow(val,3);
        }
        rho[i] = ( (4*M[i]) / (M_PI*pow(h,8)) ) * d;

        std::cout<<" rho[" << i << "]= " << rho[i] <<std::endl;
    }


    /*
    for(int i= 0 ; i < _Nb_Sommets-1; i++)
    {
        float constant = (4*M[i]) / ( M_PI*pow(h,8)) ;

        for(int j= i+1 ; j < _Nb_Sommets; j++)
        {
            float d = h*h - length2( P[i] - P[j]) ;
            if( d > 0)
            {
                float pij = pow(d,3);

                rho[i] += ( constant ) * pij ;
                rho[j] += ( constant ) * pij ;
            }
            else
                continue;
        }
    }
    */

    for(int i=0; i< _Nb_Sommets; i++)
    {
        std::cout<<" rho[" << i << "]= " << rho[i] <<std::endl;
    }
}//void


/**
 * Calcul des forces d interaction entre particules.
  * Attention - Calcul direct de fij / rho_i
 */
void ObjetSimuleSPH::CalculInteraction(float visco)
{
     for(int i= 0 ; i < _Nb_Sommets-1; i++)
    {
        float constant = (4*M[i]) / ( M_PI*pow(h,8)) ;

        for(int j= i+1 ; j < _Nb_Sommets; j++)
        {
            float d = h*h - length2( P[i] - P[j]) ;
            if( d > 0)
            {
                Vector rij = P[i] - P[j] ;
                Vector vij = V[i] - V[j] ;

                float qij = length(rij) / h ;

                float left = M[j]*(1-qij) / (M_PI*pow(h,4)*rho[j]);

                Vector right = 15*bulk*(rho[i] + rho[j] - 2*rho0) * ( (1- qij) / qij) * rij -40*visco*vij;

               Force[i] = Force[i] + left*right;
               Force[j] = Force[j] + left*right;

               //std::cout<<" Force[" << i << "]= " << left*right <<std::endl;
            }
            else
                continue;
        }
        Force[i] = Force[i] / rho[i];
    }
    Force[_Nb_Sommets-1] = Force[_Nb_Sommets-1] / rho[_Nb_Sommets-1];
}//void


/**
 * Gestion des collisions.
 * Notre condition aux limites correspond à une frontière inélastique
 * avec un coefficient de restitution spécifié < 1. Quand
 * une particule frappe une barrière verticale ([[which = 0]]) ou horizontale
 * barrier ([[which = 1]]), nous le traitons avec [[damp_reflect]].
 * Cela réduit la distance totale parcourue en fonction du temps écoulé depuis
 * la collision reflète, amortit les vitesses, et reflète
 * quels que soient les composants de la solution à refléter.
 */

void ObjetSimuleSPH::damp_reflect(int frontiere, float barrier, int indice_part)
{
    /// frontiere : indique quelle frontiere (x, y, z) du domaine est concernee


}


/**
 * Gestion des collisions.
 * Pour chacune des particules nous verifions la reflection
 * avec chacun des 4 murs du domaine.
 */
void ObjetSimuleSPH::CollisionPlan()
{
    /*
    for(int i=0; i< _Nb_Sommets; i++)
    {
        if(P[i].y < -2)
        {
            P[i].y = -2;
            V[i] = -0.5*V[i];
        }
    }
    */

}

