/*
 * CalculsRigidBody.cpp :
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
#include "ObjetSimuleRigidBody.h"
#include "Viewer.h"

using namespace std;




/*
 * Calcul de la masse de l objet rigide
 */
void ObjetSimuleRigidBody::CalculMasse()
{

    for(int i=0 ; i< _Nb_Sommets; i++)
    {

        if(M[i] == 0)
            M[i] = 1.0;
        _Mass += M[i];

        _Position = _Position + P[i] * M[i];
    }
    _Position = _Position/_Mass;

    std::cout<<"Sum of Mass: "<< _Mass << std::endl;

}


/*
 * Calcul du tenseur d inertie de l objet rigide - - partie constante Ibody
 * et remplissage tableau roi (position particules dans repere objet)
 */
void ObjetSimuleRigidBody::CalculIBody()
{

    Matrix I = Matrix::NullMatrix();
    std::cout<<I<<std::endl;
    std::cout<<Matrix::UnitMatrix()<<std::endl;
    for(int i=0; i< _Nb_Sommets; i++)
    {
        _ROi.push_back(P[i]);

        Vector ri = _ROi[i] - _Position;
        Matrix r = StarMatrix(ri);
        Matrix rt = StarMatrix(ri).TransposeConst();

        float rtr =  ri.x*ri.x + ri.y*ri.y + ri.z*ri.z;
        Matrix p1 = Matrix::UnitMatrix();

        p1 *= rtr;
        I += M[i]*( p1 - r*rt);


    }
    std::cout<<"-----------"<<std::endl;
    _Ibody = I;
    _IbodyInv = I.InverseConst();
    std::cout<<_Ibody<<std::endl;
    std::cout<<I.InverseConst()<<std::endl;
}


/*
 * Calcul de l etat de l objet rigide.
 */
void ObjetSimuleRigidBody::CalculStateX()
{
    _InertieTenseurInv = _Rotation*_IbodyInv*_Rotation.InverseConst();
    _VitesseAngulaire = _InertieTenseurInv*_MomentCinetique;
}



/*
 * Calcul de la derivee de l etat : d/dt X(t).
 */
void ObjetSimuleRigidBody::CalculDeriveeStateX(Vector gravite)
{
    _Vitesse = _QuantiteMouvement/_Mass;
    _InertieTenseurInv = _Rotation*_IbodyInv*_Rotation.TransposeConst();
    _VitesseAngulaire = _InertieTenseurInv*_MomentCinetique;
    _RotationDerivee = _VitesseAngulaire*_Rotation;

    _Force = _Mass*gravite;

    for(int i=0 ; i< _Nb_Sommets; i++)
    {
        _Torque = _Torque +  cross(( P[i] - _Position ), Force[i]);
    }
}


/**
 * Schema integration pour obbtenir X(t+dt) en fct de X(t) et d/dt X(t)
 */
void ObjetSimuleRigidBody::Solve(float visco)
{
    _Position = _Position + _Vitesse* _delta_t;
    _Rotation = _Rotation + _RotationDerivee *_delta_t;
    _QuantiteMouvement = _QuantiteMouvement + _Force*_delta_t;
    _MomentCinetique = _MomentCinetique + _Torque*_delta_t;

}//void



/**
 * Gestion des collisions avec le sol.
 */
void ObjetSimuleRigidBody::CollisionPlan()
{



}// void

