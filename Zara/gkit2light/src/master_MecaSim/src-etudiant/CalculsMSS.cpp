/*
 * CalculsMSS.cpp :
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

/** \file CalculsMSS.cpp
Programme calculant pour chaque particule i d un MSS son etat au pas de temps suivant
 (methode d 'Euler semi-implicite) : principales fonctions de calculs.
\brief Fonctions de calculs de la methode semi-implicite sur un systeme masses-ressorts.
*/

#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>

#include "vec.h"
#include "ObjetSimule.h"
#include "ObjetSimuleMSS.h"
#include "Viewer.h"

using namespace std;





/**
* Calcul des forces appliquees sur les particules du systeme masses-ressorts.
 */
void ObjetSimuleMSS::CalculForceSpring()
{
	/// f = somme_i (ki * (l(i,j)-l_0(i,j)) * uij ) + (nuij * (vi - vj) * uij) + (m*g) + force_ext

	/// Rq : Les forces dues a la gravite et au vent sont ajoutees lors du calcul de l acceleration
     for(int i=0; i< _Nb_Sommets; i++)
    {
        Vector fe = Vector(0,0,0); Vector fv = Vector(0,0,0);
        for( auto ressort : _SytemeMasseRessort->GetParticule(i)->GetRessortList())
        {
            Particule *A = ressort->GetParticuleA();
            Particule *B = ressort->GetParticuleB();
            if(B == _SytemeMasseRessort->GetParticule(i)) std::swap(A,B);

            Vector AB = B->GetPosition() - A->GetPosition();
            float norm = length(AB);

            Vector u(0,0,0);
            if(norm != 0)
                u = AB/norm;
            else
                u = Vector(0.0,0.0,0.0);

            fe  = fe + ressort->GetRaideur() * (norm - ressort->GetLrepos() )*u;
            fv  = fv +  dot(ressort->GetAmortissement() *(V[B->GetId()]- V[A->GetId()]), u)*u;


        }

        Force[i] = fe + fv   ;

    }

}//void


/**
 * Gestion des collisions avec le sol.
 */
void ObjetSimuleMSS::CollisionPlan()
{
    /// Arret de la vitesse quand touche le plan
   for(int i=0; i< P.size(); i++)
    {
        if(P[i].y < -10.0)
        {
            P[i].y = -10;
            V[i] = -V[i]*0.1;
        }
    }

}// void


void ObjetSimuleMSS::Wind()
{
    for(int i=0; i< _Nb_Sommets; i++)
    {
        if(P[i].y > -1 && P[i].y < 0)
        {
            A[i] = Vector(0,0,0);
            Force[i] = Vector(0,0,0);
        }
    }
}
