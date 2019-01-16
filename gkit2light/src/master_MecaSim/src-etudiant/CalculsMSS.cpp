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
    std::cout<<"Size: "<<P.size()<<std::endl;

    for(int i=0; i< P.size()-1; i++)
    {
        int nb_voisin = _SytemeMasseRessort->_ParticuleList[i]->GetNbVoisins();
        for(int j=0; j<nb_voisin; j++)
        {

        }
    }
	/// Rq : Les forces dues a la gravite et au vent sont ajoutees lors du calcul de l acceleration


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
         V[i] = -V[i]*0.9;
       }
    }

}// void

