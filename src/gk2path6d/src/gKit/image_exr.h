
#ifndef _IMAGE_EXR_H
#define _IMAGE_EXR_H

#include "image.h"

//! charge une image a partir d'un fichier .hdr. renvoie Image::error() en cas d'echec. a detruire avec image::release( ).
//! \param filemane nom de l'image .hdr a charger
Image read_image_exr( const char *filename );

//! enregistre une image dans un fichier .hdr.
int write_image_exr( const Image& image, const char *filename );

//! renvoie vrai si le nom de fichier se termine par .hdr.
bool is_exr_image( const char *filename );

#endif
