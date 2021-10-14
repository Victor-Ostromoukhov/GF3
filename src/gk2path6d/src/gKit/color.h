
#ifndef _COLOR_H
#define _COLOR_H

#include <algorithm>

#include "real.h"


//! \addtogroup image
///@{

//! \file
//! manipulation de couleurs

//! representation d'une couleur (rgba) transparente ou opaque.
template< typename T >
struct ColorT
{
    //! constructeur par defaut.
    ColorT( ) : r(T(0)), g(T(0)), b(T(0)), a(T(1)) {}
    ColorT( const ColorT& color ) : r(color.r), g(color.g), b(color.b), a(color.a) {} 
    ColorT& operator= ( const ColorT& color ) { r= color.r; g= color.g; b= color.b; a= color.a; return *this; }
    explicit ColorT( const T _r, const T _g, const T _b, const T _a= T(1) ) : r(_r), g(_g), b(_b), a(_a) {}
    explicit ColorT( const T _value ) : r(_value), g(_value), b(_value), a(T(1)) {}
    
    //! cree une couleur avec les memes composantes que color, mais remplace sa composante alpha (color.r, color.g, color.b, alpha).
    ColorT( const ColorT& color, const T alpha ) : r(color.r), g(color.g), b(color.b), a(alpha) {}  // remplace alpha.
    
    operator ColorT<float>( ) const { return ColorT<float>(float(r), float(g), float(b), float(a)); }
    
    T power( ) const  { return std::max(r, std::max(g, std::max(b, T(0)))); }
    T grey( ) const  { return T(0.3) * r + T(0.6) * g + T(0.1) * b; }
    
    T r, g, b, a;
};

typedef ColorT<Float> Color;
typedef ColorT<float> Color32;


template< typename T>
ColorT<T> operator+ ( const ColorT<T>& a, const ColorT<T>& b )
{
    return ColorT<T>(a.r + b.r, a.g + b.g, a.b + b.b, a.a + b.a);
}

template< typename T>
ColorT<T> operator- ( const ColorT<T>& c )
{
    return ColorT<T>(-c.r, -c.g, -c.b, -c.a);
}

template< typename T>
ColorT<T> operator- ( const ColorT<T>& a, const ColorT<T>& b )
{
    return a + (-b);
}

template< typename T>
ColorT<T> operator* ( const ColorT<T>& a, const ColorT<T>& b )
{
    return ColorT<T>(a.r * b.r, a.g * b.g, a.b * b.b, a.a * b.a);
}

template< typename T>
ColorT<T> operator* ( const T k, const ColorT<T>& c )
{
    return ColorT<T>(c.r * k, c.g * k, c.b * k, c.a * k);
}

template< typename T>
ColorT<T> operator* ( const ColorT<T>& c, const T k )
{
    return k * c;
}

template< typename T>
ColorT<T> operator/ ( const ColorT<T>& a, const ColorT<T>& b )
{
    return ColorT<T>(a.r / b.r, a.g / b.g, a.b / b.b, a.a / b.a);
}

template< typename T>
ColorT<T> operator/ ( const T k, const ColorT<T>& c )
{
    return ColorT<T>(k / c.r, k / c.g, k / c.b, k / c.a);
}

template< typename T>
ColorT<T> operator/ ( const ColorT<T>& c, const T k )
{
    T kk= 1 / k;
    return kk * c;
}

//! utilitaire. renvoie une couleur noire.
Color Black( );// { return Color(0, 0, 0); }

//! utilitaire. renvoie une couleur blanche.
Color White( ); // { return Color(1, 1, 1); }

//! utilitaire. renvoie une couleur rouge.
Color Red( ); // { return Color(1, 0, 0); }

//! utilitaire. renvoie une couleur verte.
Color Green( ); // { return Color(0, 1, 0); }

//! utilitaire. renvoie une couleur bleue.
Color Blue( ); // { return Color(0, 0, 1); }

//! utilitaire. renvoie une couleur jaune.
Color Yellow( ); // { return Color(1, 1, 0); }

///@}
#endif
