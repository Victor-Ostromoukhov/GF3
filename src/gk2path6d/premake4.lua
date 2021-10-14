solution "gKit2light"
	configurations { "debug", "release" }

	platforms { "x64" }

	includedirs { ".", "src/gKit" }

--~ 	defines { "FLOAT32" }

	defines { "REMAP_BOUNCE"}	-- utilise les memes dimensions pour source + brdf
--~ 	defines { "REMAP_BOUNCE1"}	-- idem uniquement 1 lobe brdf + 1 source == selecteur a 0 dans le cas BRDF_3D et SOURCE_3D (ou plutot !BRDF_2D !SOURCE_2D)
	defines { "BRDF_2D"}		-- 2d + normalisation au lieu de 3d
	defines { "SOURCE_2D" }		-- idem
-- 	defines { "SOURCE_QUADS" }	-- uniquement des sources quads
 	defines { "SOURCE_QUAD1" }	-- uniquement 1 source quad
	defines { "TRIANGLE_MAPPING" }
--~ 	defines { "LENS_SQUARE" }
--~ 	defines { "EXPORT_SAMPLES" }
	configuration "debug"
		defines { "DEBUG" }
		if _PREMAKE_VERSION >="5.0" then
			symbols "on"
		else
			flags { "Symbols" }
		end

	configuration "release"
--~ 		defines { "NDEBUG" }
--~ 		defines { "GK_RELEASE" }
		if _PREMAKE_VERSION >="5.0" then
			optimize "speed"
		else
			flags { "OptimizeSpeed" }
		end

	configuration "linux"
		buildoptions { "-mtune=native -march=native" }
		buildoptions { "-std=c++11" }
		buildoptions { "-W -Wall -Wextra -Wsign-compare -Wno-unused-parameter -Wno-unused-function -Wno-unused-variable", "-pipe" }
		buildoptions { "-flto"}
--~ 		buildoptions { "-Weffc++"}
		linkoptions { "-flto"}
		buildoptions { "-fopenmp" }
		linkoptions { "-fopenmp -L/usr/local/lib" }
		links { "GLEW", "SDL2", "SDL2_image", "SDL2_ttf", "GL", "IlmImf", "Imath", "Half" }
		defines { "EXR"}
		defines { "THREADS" }
		defines { "TTF" }

	configuration { "linux", "debug" }
		linkoptions { "-g"}	-- bugfix : premake4 ne genere pas le flag debug pour le linker...

	configuration { "windows" }
		defines { "WIN32", "_USE_MATH_DEFINES", "_CRT_SECURE_NO_WARNINGS" }
		defines { "NOMINMAX" } -- allow std::min() and std::max()

		includedirs { "extern/visual/include" }
		libdirs { "extern/visual/lib" }
		links { "opengl32", "glew32", "SDL2", "SDL2main", "SDL2_image" }

		buildoptions { "/openmp" }

		if _PREMAKE_VERSION >="5.0" then
			system "Windows"
			architecture "x64"
			flags { "MultiProcessorCompile", "NoMinimalRebuild" }
			disablewarnings { "4244", "4305" }
		end

	configuration { "windows", "release" }
		if _PREMAKE_VERSION >="5.0" then
			flags { "LinkTimeOptimization" }
		end

	configuration "macosx"
		frameworks= "-F /Library/Frameworks/"
		buildoptions { "-std=c++11" ,"-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include"}
		defines { "GK_MACOS","GL_SILENCE_DEPRECATION" }
		defines { "USE_32BITS_SEED" }
		buildoptions { frameworks }
		linkoptions { frameworks .. " -framework OpenGL"}
		defines { "EXR" }
		linkoptions { "-lSDL2 -lSDL2_image -L/usr/local/lib" }
		linkoptions { "-lIlmImf", "-lImath", "-lHalf","/usr/local/opt/libomp/lib/libomp.dylib" }


 -- description des fichiers communs
gkit_files = { "src/gKit/*.cpp", "src/gKit/*.h" }


project("path6d")
	language "C++"
	kind "ConsoleApp"
	targetdir "bin"
	includedirs { "path" }

	files ( gkit_files )
	files { "path/main.cpp", 
		"path/stbvh_builder.cpp", "path/bvh_builder.cpp", "path/stbvh.cpp", "path/bvh.cpp", 
		"path/framebuffer.cpp", "path/texturelib.cpp", "path/materials.cpp", "path/scene.cpp", 
		"path/sampler_s19.cpp", "path/sobol_spack.cpp", 
		"../Screenspace/ZSampler/zhash.cpp", 
		"../Screenspace/ZSampler/z.cpp", 
		"../Screenspace/ZSampler/zpp.cpp",
		"../Screenspace/ZSampler/art2x2-table.cpp"
		}
		
project("tile")
	language "C++"
	kind "ConsoleApp"
	targetdir "bin"
	includedirs { "path" }

	files ( gkit_files )
	files { "path/tile.cpp", 
		"path/sampler_s19.cpp", "path/sobol_spack.cpp", 
		"../Screenspace/ZSampler/zhash.cpp", 
		"../Screenspace/ZSampler/z.cpp", 
		"../Screenspace/ZSampler/zpp.cpp",
		"../Screenspace/ZSampler/art2x2-table.cpp"
		}
		
project("rank4d")
	language "C++"
	kind "ConsoleApp"
	targetdir "bin"
	includedirs { "path" }

	files ( gkit_files )
	files { "path/rank4d.cpp", 
		"path/stbvh_builder.cpp", "path/bvh_builder.cpp", "path/stbvh.cpp", "path/bvh.cpp", 
		"path/framebuffer.cpp", "path/texturelib.cpp", "path/materials.cpp", "path/scene.cpp", 
		"path/sampler_s19.cpp", "path/sobol_spack.cpp", 
		"../Screenspace/ZSampler/zhash.cpp", 
		"../Screenspace/ZSampler/z.cpp", 
		"../Screenspace/ZSampler/zpp.cpp",
		"../Screenspace/ZSampler/art2x2-table.cpp"
		}

project("blue6d")
	language "C++"
	kind "ConsoleApp"
	targetdir "bin"
	includedirs { "path" }

	files ( gkit_files )
	files { "path/blue.cpp", 
		"path/stbvh_builder.cpp", "path/bvh_builder.cpp", "path/stbvh.cpp", "path/bvh.cpp", 
		"path/framebuffer.cpp", "path/texturelib.cpp", "path/materials.cpp", "path/scene.cpp", 
		"path/sobol_spack.cpp" }

project("scene_viewer")
	language "C++"
	kind "ConsoleApp"
	targetdir "bin"
	includedirs { "path" }
	files ( gkit_files )
	files { "src/scene_viewer.cpp", 		
		"path/stbvh_builder.cpp", "path/bvh_builder.cpp", "path/stbvh.cpp", "path/bvh.cpp", 
		"path/framebuffer.cpp", "path/texturelib.cpp", "path/materials.cpp", "path/scene.cpp", 
		"path/sobol_spack.cpp" }

project("rank_viewer")
	language "C++"
	kind "ConsoleApp"
	targetdir "bin"
	includedirs { "path" }
	files ( gkit_files )
	files { "src/rank_viewer.cpp", 		
		"path/stbvh_builder.cpp", "path/bvh_builder.cpp", "path/stbvh.cpp", "path/bvh.cpp", 
		"path/framebuffer.cpp", "path/texturelib.cpp", "path/materials.cpp", "path/scene.cpp", 
		"path/sobol_spack.cpp" }

--~ project("mse")
--~ 	language "C++"
--~ 	kind "ConsoleApp"
--~ 	targetdir "bin"
--~
--~ 	files { "mse.cpp" }


 -- description des projets
projects = {
	"shader_kit",
	"image_viewer",
	"hdrtopng",
	"grid",
	"seeds",
	"stats",
	"ranks"
}

for i, name in ipairs(projects) do
	project(name)
		language "C++"
		kind "ConsoleApp"
		targetdir "bin"
		includedirs { "../evalRender" }
		files ( gkit_files )
		files { "src/" .. name..'.cpp' }
end
