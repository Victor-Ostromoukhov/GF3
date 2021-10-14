
#include <cstdio>
#include <vector>

#include "mesh.h"
#include "wavefront.h"
#include "mesh_cache.h"

#include "scene.h"


struct Object
{
    std::string name;
    std::string mesh;
    int mesh_start, mesh_stop;
    Transform transform;
    
    std::vector<Transform> frames;
    
    Object( ) : name(), mesh(), mesh_start(0), mesh_stop(0), transform(), frames() {}
    Object( const char *_name ) : name(_name), mesh_start(0), mesh_stop(0), transform(), frames() {}
};


Scene read_scene( const char *filename )
{
    Scene scene;
    
    FILE *in= fopen(filename, "rt");
    if(in == nullptr)
    {
        printf("[error] loading scene '%s'...\n", filename);
        return scene;
    }
    
    std::string path= normalize_path(pathname(filename));
    
    std::vector<Object> objects;
    Object *object= nullptr;
    
    char tmp[1024];
    char line_buffer[1024];
    bool error= true;
    for(;;)
    {
        // charge une ligne du fichier
        if(fgets(line_buffer, sizeof(line_buffer), in) == NULL)
        {
            error= false;       // fin du fichier, pas d'erreur detectee
            break;
        }
        
        // force la fin de la ligne, au cas ou
        line_buffer[sizeof(line_buffer) -1]= 0;
        
        // saute les espaces en debut de ligne
        char *line= line_buffer;
        while(*line && isspace(*line))
            line++;
        
        if(line[0] == '#')
            continue;
        
        if(line[0] == 'n')
        {
            if(sscanf(line, "newobj %[^\r\n]", tmp) == 1)
            {
                objects.emplace_back(tmp);
                object= &objects.back();
            }
        }
        
        if(!object)
            continue;
        
        if(line[0] == 'm')
        {
            if(sscanf(line, "mesh %[^\r\n]", tmp) == 1)
            {
                object->mesh= normalize_path(tmp);
                if(relative_path(object->mesh))
                    object->mesh= path + object->mesh;
            }
            
            //~ \todo matieres de l'instance...
            //~ if(sscanf(line, "materials %[^\r\n]", tmp) == 1)
        }
        else if(line[0] == 'f')
        {
            int start, stop;
            if(sscanf(line, "frames %d %d", &start, &stop) == 2)
            {
                object->mesh_start= start;
                object->mesh_stop= stop;
            }
            
            else if(sscanf(line, "frame %d", &start) == 1)
            {
                object->frames.push_back( object->transform );
                object->transform= Identity();
            }
        }
        else if(line[0] == 't')
        {
            double x, y, z;
            if(sscanf(line, "translation %lf %lf %lf", &x, &y, &z) == 3)
                object->transform= Translation(x, y, z) * object->transform;
        }
        else if(line[0] == 'r')
        {
            double a, x, y, z;
            if(sscanf(line, "rotation %lf %lf %lf %lf", &a, &x, &y, &z) == 4)       // degrees
                object->transform= Rotation(Vector(x, y, z), a) * object->transform;
        }
        else if(line[0] == 's')
        {
            double x, y, z;
            if(sscanf(line, "scale %lf %lf %lf", &x, &y, &z) == 3)
                object->transform= Scale(x, y, z) * object->transform;
        }
        
        else if(line[0] == 'u')
        {
            if(sscanf(line, "useobj %[^\r\n]", tmp) == 1)
            {
                for(int i= 0; i < int(objects.size()); i++)
                {
                    if(objects[i].name == std::string(tmp))
                    {
                        if(objects[i].mesh_start != objects[i].mesh_stop)
                        {
                            //~ printf("obj '%s', keyframes %d..%d\n", objects[i].name.c_str(), objects[i].mesh_start, objects[i].mesh_stop);
                            // key frames
                            scene.attach( SceneObject(objects[i].mesh.c_str(), objects[i].mesh_start, objects[i].mesh_stop, objects[i].transform) );
                        }
                        else
                        {
                            //~ printf("obj '%s', %d transforms\n", objects[i].name.c_str(), int(objects[i].frames.size()));
                            // animated transforms
                            if(objects[i].frames.size())
                            {
                                SceneObject move;
                                for(int f= 0; f < int(objects[i].frames.size()); f++)
                                    move.read_mesh(objects[i].mesh.c_str(), objects[i].frames[f]);
                                // \todo eventuellement composer avec la transformation courante... objects[i].transform
                                
                                scene.attach(move);
                            }
                            else
                                // static object
                                scene.attach( SceneObject(objects[i].mesh.c_str(), objects[i].transform) );
                        }
                        
                        break;
                    }
                }
            }
            
            else if(sscanf(line, "usecamera %[^\r\n]", tmp) == 1)
            {
                std::string filename= normalize_path(tmp);
                if(relative_path(filename))
                    filename= path + filename;
                
                scene.camera_filename= filename;
            }
            
            else if(sscanf(line, "uselens %[^\r\n]", tmp) == 1)
            {
                std::string filename= normalize_path(tmp);
                if(relative_path(filename))
                    filename= path + filename;
                
                scene.lens_filename= filename;
            }
        }
    }
    
    fclose(in);
    if(error)
        printf("[error] parsing line :\n%s\n", line_buffer);
    
    return scene;
}

// transforme les positions / normales des triangles du mesh.
Mesh transform( const Mesh& data, const Transform& m )
{
    
    std::vector<vec3> positions;
    positions.reserve(data.positions().size());    
    for(int i= 0; i < int(data.positions().size()); i++)
        positions.push_back( vec3(m(Point(data.positions()[i]))) );
        
    Transform mn= m.normal();
    std::vector<vec3> normals;
    normals.reserve(data.normals().size());    
    for(int i= 0; i < int(data.normals().size()); i++)
        normals.push_back( vec3(mn(Vector(data.normals()[i]))) );
    
    Mesh mesh(data.primitives());
    mesh.positions(positions);
    mesh.texcoords(data.texcoords());
    mesh.normals(normals);
    mesh.indices(data.indices());
    mesh.mesh_materials(data.mesh_materials());
    mesh.materials(data.materials());    
    return mesh;
}


bool Scene::intersect( Raydata& data ) const
{
    IntersectionAlpha filter(*this);
    bvh.intersect(data.ray, data.hit, filter);
    return data.intersect();
    
    // todo verifier que des textures semi transparentes existent... 
}


bool Scene::visible( const Float t, const Point& o, const Point& e ) const
{
    IntersectionAlpha filter(*this);
    return bvh.visible(t, o, e, filter);
}

bool Scene::visible( const Float t, const Point& o, const Vector& d ) const
{
    IntersectionAlpha filter(*this);
    return bvh.visible(t, o, d, filter);
}    
