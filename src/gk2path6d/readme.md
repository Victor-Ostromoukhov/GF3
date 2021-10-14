# installation
dependances : sdl2 sdl2_image glew openexr

installation ubuntu : `apt-get install libsdl2-dev libsdl2-image-dev libglew-dev libopenexr-dev premake4`

installation mac : 
			      brew reinstall openexr
			      brew reinstall ilmbase

# documentation
base de code : http://perso.univ-lyon1.fr/jean-claude.iehl/Public/educ/M1IMAGE/html/index.html

# generer les projets / compiler

premake4 gmake

make help  --  pour lister les projets et les options de compilation.

make config=release64 -j16 path6d image_viewer scene_viewer

les executables seront dans `bin/`
