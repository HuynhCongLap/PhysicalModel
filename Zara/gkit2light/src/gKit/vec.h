
#ifndef _VEC_H
#define _VEC_H


//! \addtogroup math
///@{

//! \file
//! operations sur points et vecteurs

//! declarations anticipees.
struct vec2;
struct vec3;
struct vec4;
struct Vector;

//! representation d'un point 3d.
struct Point
{
    //! constructeur par defaut.
    Point( const float _x= 0, const float _y= 0, const float _z= 0 ) : x(_x), y(_y), z(_z) {}

    //! cree un point a partir des coordonnees du vecteur generique (v.x, v.y, v.z).
    explicit Point( const vec3& v );   // l'implementation se trouve en fin de fichier, la structure vec3 n'est pas encore connue.
    //! cree un point a partir des coordonnes du vecteur (v.x, v.y, v.z).
    explicit Point( const Vector& v );   // l'implementation se trouve en fin de fichier, la structure vector n'est pas encore connue.

    float x, y, z;
};

//! renvoie la distance etre 2 points.
float distance( const Point& a, const Point& b );
//! renvoie le carre de la distance etre 2 points.
float distance2( const Point& a, const Point& b );

//! renvoie le milieu du segment ab.
Point center( const Point& a, const Point& b );


//! representation d'un vecteur 3d.
struct Vector
{
    //! constructeur par defaut.
    Vector( const float _x= 0, const float _y= 0, const float _z= 0) : x(_x), y(_y), z(_z) {}
    //! cree le vecteur ab.
    Vector( const Point& a, const Point& b ) : x(b.x - a.x), y(b.y - a.y), z(b.z - a.z) {}

    //! cree un vecteur a partir des coordonnees du vecteur generique (v.x, v.y, v.z).
    explicit Vector( const vec3& v );   // l'implementation se trouve en fin de fichier, la structure vec3 n'est pas encore connue.
    //! cree un vecteur a partir des coordonnes du vecteur (v.x, v.y, v.z).
    explicit Vector( const Point& a );   // l'implementation se trouve en fin de fichier.

    float x, y, z;
};

//! renvoie un vecteur unitaire / longueur == 1.
Vector normalize( const Vector& v );
//! renvoie le produit vectoriel de 2 vecteurs.
Vector cross( const Vector& u, const Vector& v );
//! renvoie le produit scalaire de 2 vecteurs.
float dot( const Vector& u, const Vector& v );
//! renvoie la longueur d'un vecteur.
float length( const Vector& v );
//! renvoie la carre de la longueur d'un vecteur.
float length2( const Vector& v );

//! renvoie le vecteur a - b.
Vector operator- ( const Point& a, const Point& b );
//! renvoie le vecteur -v.
Vector operator- ( const Vector& v );

//! renvoie le point a+v.
Point operator+ ( const Point& a, const Vector& v );
//! renvoie le point a+v.
Point operator+ ( const Vector& v, const Point& a );
//! renvoie le point a-v.
Point operator- ( const Vector& v, const Point& a );
//! renvoie le point a-v.
Point operator- ( const Point& a, const Vector& v );
//! renvoie le vecteur u+v.
Vector operator+ ( const Vector& u, const Vector& v );
//! renvoie le vecteur u-v.
Vector operator- ( const Vector& u, const Vector& v );
//! renvoie le vecteur k*u;
Vector operator* ( const float k, const Vector& v );
//! renvoie le vecteur k*v;
Vector operator* ( const Vector& v, const float k );
//! renvoie le vecteur v/k;
Vector operator/ ( const Vector& v, const float k );


//! vecteur generique, utilitaire.
struct vec2
{
    //! constructeur par defaut.
    vec2( const float _x= 0, const float _y= 0 ) : x(_x), y(_y) {}

    float x, y;
};


//! vecteur generique, utilitaire.
struct vec3
{
    //! constructeur par defaut.
    vec3( const float _x= 0, const float _y= 0, const float _z= 0 ) : x(_x), y(_y), z(_z) {}
    //! constructeur par defaut.
    vec3( const vec2& a, const float _z ) : x(a.x), y(a.y), z(_z) {}

    //! cree un vecteur generique a partir des coordonnees du point a.
    explicit vec3( const Point& a );    // l'implementation se trouve en fin de fichier.
    //! cree un vecteur generique a partir des coordonnees du vecteur v.
    explicit vec3( const Vector& v );    // l'implementation se trouve en fin de fichier.

    float x, y, z;
};


//! vecteur generique 4d, ou 3d homogene, utilitaire.
struct vec4
{
    //! constructeur par defaut.
    vec4( const float _x= 0, const float _y= 0, const float _z= 0, const float _w= 0 ) : x(_x), y(_y), z(_z), w(_w) {}
    //! constructeur par defaut.
    vec4( const vec2& v, const float _z= 0, const float _w= 0 ) : x(v.x), y(v.y), z(_z), w(_w) {}
    //! constructeur par defaut.
    vec4( const vec3& v, const float _w= 0 ) : x(v.x), y(v.y), z(v.z), w(_w) {}

    //! cree un vecteur generique a partir des coordonnees du point a, (a.x, a.y, a.z, 1).
    explicit vec4( const Point& a );    // l'implementation se trouve en fin de fichier.
    //! cree un vecteur generique a partir des coordonnees du vecteur v, (v.x, v.y, v.z, 0).
    explicit vec4( const Vector& v );    // l'implementation se trouve en fin de fichier.

    float x, y, z, w;
};


// implementation des constructeurs explicites.
inline Point::Point( const vec3& v ) : x(v.x), y(v.y), z(v.z) {}
inline Point::Point( const Vector& v ) : x(v.x), y(v.y), z(v.z) {}

inline Vector::Vector( const vec3& v ) : x(v.x), y(v.y), z(v.z) {}
inline Vector::Vector( const Point& a ) : x(a.x), y(a.y), z(a.z) {}

inline vec3::vec3( const Point& a ) : x(a.x), y(a.y), z(a.z) {}
inline vec3::vec3( const Vector& v ) : x(v.x), y(v.y), z(v.z) {}

inline vec4::vec4( const Point& a ) : x(a.x), y(a.y), z(a.z), w(1.f) {}
inline vec4::vec4( const Vector& v ) : x(v.x), y(v.y), z(v.z), w(0.f) {}

#include <iostream>

inline std::ostream& operator<<(std::ostream& o, const Point& p)
{
    o<<"p("<<p.x<<","<<p.y<<","<<p.z<<")";
    return o;
}

inline std::ostream& operator<<(std::ostream& o, const Vector& v)
{
    o<<"v("<<v.x<<","<<v.y<<","<<v.z<<")";
    return o;
}

///@}
#endif
