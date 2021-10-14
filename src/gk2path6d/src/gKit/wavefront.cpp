
#include <cassert>
#include <cstdio>
#include <ctype.h>
#include <climits>

#include <string>
#include <algorithm>

#include "wavefront.h"


/*! renvoie le chemin d'acces a un fichier. le chemin est toujours termine par /
    pathname("path/to/file") == "path/to/"
    pathname("file") == "./"
 */
std::string pathname( const std::string& filename )
{
    std::string path= filename;
#ifndef WIN32
    std::replace(path.begin(), path.end(), '\\', '/');   // linux, macos : remplace les \ par /.
    size_t slash = path.find_last_of( '/' );
    if(slash != std::string::npos)
        return path.substr(0, slash +1); // inclus le slash
    else
        return "./";
#else
    std::replace(path.begin(), path.end(), '/', '\\');   // windows : remplace les / par \.
    size_t slash = path.find_last_of( '\\' );
    if(slash != std::string::npos)
        return path.substr(0, slash +1); // inclus le slash
    else
        return ".\\";
#endif
}


std::string normalize_path( const std::string& file )
{
    std::string tmp= file;
    
#ifndef WIN32
    std::replace(tmp.begin(), tmp.end(), '\\', '/');   // linux, macos : remplace les \ par /.
#else
    std::replace(tmp.begin(), tmp.end(), '/', '\\');   // windows : remplace les / par \.
#endif
    return tmp;
}

bool relative_path( const std::string& filename )
{
    if(filename.empty()) 
        return false;
    
    return (filename[0] != '/' && filename[0] != '\\' && filename[0]!='.');
}


Mesh read_mesh( const char *filename )
{
    FILE *in= fopen(filename, "rt");
    if(in == NULL)
    {
        printf("[error] loading mesh '%s'...\n", filename);
        return Mesh::error();
    }
    
    Mesh data(GL_TRIANGLES);
    
    printf("loading mesh '%s'...\n", filename);
    
    std::vector<vec3> positions;
    std::vector<vec2> texcoords;
    std::vector<vec3> normals;
    MaterialLib materials;
    int default_material_id= -1;
    int material_id= -1;
    
    std::vector<int> idp;
    std::vector<int> idt;
    std::vector<int> idn;
    
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
        
        if(line[0] == 'v')
        {
            double x, y, z;
            if(line[1] == ' ')          // position x y z
            {
                if(sscanf(line, "v %lf %lf %lf", &x, &y, &z) != 3)
                    break;
                positions.push_back( vec3(x, y, z) );
            }
            else if(line[1] == 'n')     // normal x y z
            {
                if(sscanf(line, "vn %lf %lf %lf", &x, &y, &z) != 3)
                    break;
                normals.push_back( vec3(x, y, z) );
            }
            else if(line[1] == 't')     // texcoord x y
            {
                if(sscanf(line, "vt %lf %lf", &x, &y) != 2)
                    break;
                texcoords.push_back( vec2(x, y) );
            }
        }
        
        else if(line[0] == 'f')         // triangle a b c, les sommets sont numerotes a partir de 1 ou de la fin du tableau (< 0)
        {
            idp.clear();
            idt.clear();
            idn.clear();
            
            int next= 0;
            for(line= line +1; ; line= line + next)
            {
                idp.push_back(0); 
                idt.push_back(0); 
                idn.push_back(0);         // 0: invalid index
                
                next= 0;
                if(sscanf(line, " %d/%d/%d %n", &idp.back(), &idt.back(), &idn.back(), &next) == 3) 
                    continue;
                else if(sscanf(line, " %d/%d %n", &idp.back(), &idt.back(), &next) == 2)
                    continue;
                else if(sscanf(line, " %d//%d %n", &idp.back(), &idn.back(), &next) == 2)
                    continue;
                else if(sscanf(line, " %d %n", &idp.back(), &next) == 1)
                    continue;
                else if(next == 0)      // fin de ligne
                    break;
            }
            
            // force une matiere par defaut, si necessaire
            if(material_id == -1)
            {
                if(default_material_id == -1)
                    // creer une matiere par defaut
                    default_material_id= data.mesh_material(Material());
                
                material_id= default_material_id;
                data.material(material_id);
                
                printf("usemtl default\n");
            }
            
            for(int v= 2; v +1 < (int) idp.size(); v++)
            {
                int idv[3]= { 0, v -1, v };
                for(int i= 0; i < 3; i++)
                {
                    int k= idv[i];
                    int p= (idp[k] < 0) ? int(positions.size()) + idp[k] : idp[k] -1;
                    int t= (idt[k] < 0) ? int(texcoords.size()) + idt[k] : idt[k] -1;
                    int n= (idn[k] < 0) ? int(normals.size())   + idn[k] : idn[k] -1;
                    
                    if(p < 0) break; // error
                    if(t >= 0) data.texcoord(texcoords[t]);
                    if(n >= 0) data.normal(normals[n]);
                    data.vertex(positions[p]);
                }
            }
        }
        
        else if(line[0] == 'm')
        {
           if(sscanf(line, "mtllib %[^\r\n]", tmp) == 1)
           {
               materials= read_materials( std::string(pathname(filename) + tmp).c_str() );
               // enregistre les matieres dans le mesh
               data.mesh_materials(materials.data);
           }
        }
        
        else if(line[0] == 'u')
        {
           if(sscanf(line, "usemtl %[^\r\n]", tmp) == 1)
           {
               material_id= -1;
               for(unsigned int i= 0; i < (unsigned int) materials.names.size(); i++)
                    if(materials.names[i] == tmp)
                        material_id= i;
                
                if(material_id == -1)
                {
                    // force une matiere par defaut, si necessaire
                    if(default_material_id == -1)
                        default_material_id= data.mesh_material(Material());
                    
                    material_id= default_material_id;
                }
                
                // selectionne une matiere pour le prochain triangle
                data.material(material_id);
           }
        }
    }
    
    fclose(in);
    
    if(error)
        printf("loading mesh '%s'...\n[error]\n%s\n\n", filename, line_buffer);
    
    return data;
}


MeshData read_mesh_data( const char *filename )
{
    FILE *in= fopen(filename, "rt");
    if(in == NULL)
    {
        printf("[error] loading mesh data '%s'...\n", filename);
        return MeshData();
    }
    
    printf("loading mesh data '%s'...\n", filename);
    
    MeshData data;
    MaterialLib materials;
    int default_material_id= -1;
    int material_id= -1;
    
    std::vector<int> idp;
    std::vector<int> idt;
    std::vector<int> idn;
    
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
        
        if(line[0] == 'v')
        {
            double x, y, z;
            if(line[1] == ' ')          // position x y z
            {
                if(sscanf(line, "v %lf %lf %lf", &x, &y, &z) != 3)
                    break;
                data.positions.push_back( vec3(x, y, z) );
            }
            else if(line[1] == 'n')     // normal x y z
            {
                if(sscanf(line, "vn %lf %lf %lf", &x, &y, &z) != 3)
                    break;
                data.normals.push_back( vec3(x, y, z) );
            }
            else if(line[1] == 't')     // texcoord x y
            {
                if(sscanf(line, "vt %lf %lf", &x, &y) != 2)
                    break;
                data.texcoords.push_back( vec2(x, y) );
            }
        }
        
        else if(line[0] == 'f')         // triangle a b c, les sommets sont numerotes a partir de 1 ou de la fin du tableau (< 0)
        {
            idp.clear();
            idt.clear();
            idn.clear();
            
            int next= 0;
            for(line= line +1; ; line= line + next)
            {
                idp.push_back(0); 
                idt.push_back(0); 
                idn.push_back(0);         // 0: invalid index
                
                next= 0;
                if(sscanf(line, " %d/%d/%d %n", &idp.back(), &idt.back(), &idn.back(), &next) == 3) 
                    continue;
                else if(sscanf(line, " %d/%d %n", &idp.back(), &idt.back(), &next) == 2)
                    continue;
                else if(sscanf(line, " %d//%d %n", &idp.back(), &idn.back(), &next) == 2)
                    continue;
                else if(sscanf(line, " %d %n", &idp.back(), &next) == 1)
                    continue;
                else if(next == 0)      // fin de ligne
                    break;
            }
            
            // force une matiere par defaut, si necessaire
            if(material_id == -1)
            {
                if(default_material_id == -1)
                {
                    // creer une matiere par defaut
                    default_material_id= int(data.materials.size());
                    data.materials.push_back(Material());
                }
                
                material_id= default_material_id;
                printf("usemtl default\n");
            }
            
            for(int v= 2; v +1 < int(idp.size()); v++)
            {
                int idv[3]= { 0, v -1, v };
                for(int i= 0; i < 3; i++)
                {
                    int k= idv[i];
                    int p= (idp[k] < 0) ? int(data.positions.size()) + idp[k] : idp[k] -1;
                    int t= (idt[k] < 0) ? int(data.texcoords.size()) + idt[k] : idt[k] -1;
                    int n= (idn[k] < 0) ? int(data.normals.size())   + idn[k] : idn[k] -1;
                    
                    if(p < 0) break; // error
                    
                    data.position_indices.push_back(p);
                    data.texcoord_indices.push_back(t);
                    data.normal_indices.push_back(n);
                }
                
                assert(material_id != -1);
                data.material_indices.push_back(material_id);
            }
        }
        
        else if(line[0] == 'm')
        {
           if(sscanf(line, "mtllib %[^\r\n]", tmp) == 1)
           {
               materials= read_materials( std::string(pathname(filename) + tmp).c_str() );
               data.materials= materials.data;
           }
        }
        
        else if(line[0] == 'u')
        {
           if(sscanf(line, "usemtl %[^\r\n]", tmp) == 1)
           {
               material_id= -1;
               for(unsigned int i= 0; i < (unsigned int) materials.names.size(); i++)
                    if(materials.names[i] == tmp)
                        material_id= i;
                
                if(material_id == -1)
                {
                    // force une matiere par defaut, si necessaire
                    if(default_material_id == -1)
                    {
                        // creer une matiere par defaut
                        default_material_id= int(data.materials.size());
                        data.materials.push_back(Material());
                    }
                    
                    material_id= default_material_id;
                }
           }
        }
    }
    
    fclose(in);
    
    if(error)
        printf("loading mesh data '%s'...\n[error]\n%s\n\n", filename, line_buffer);
    
    return data;
}


int write_mesh( const Mesh& mesh, const char *filename )
{
    if(mesh == Mesh::error())
        return -1;

    if(mesh.primitives() != GL_TRIANGLES)
        return -1;
    if(mesh.positions().size() == 0)
        return -1;
    if(filename == NULL)
        return -1;
    
    FILE *out= fopen(filename, "wt");
    if(out == NULL)
        return -1;
    
    printf("writing mesh '%s'...\n", filename);
    
    const std::vector<vec3>& positions= mesh.positions();
    for(unsigned int i= 0; i < (unsigned int) positions.size(); i++)
        fprintf(out, "v %lf %lf %lf\n", double(positions[i].x), double(positions[i].y), double(positions[i].z));
    fprintf(out, "\n");
    
    const std::vector<vec2>& texcoords= mesh.texcoords();
    bool has_texcoords= (texcoords.size() == positions.size());
    for(unsigned int i= 0; i < (unsigned int) texcoords.size(); i++)
        fprintf(out, "vt %lf %lf\n", double(texcoords[i].x), double(texcoords[i].y));
    fprintf(out, "\n");
    
    const std::vector<vec3>& normals= mesh.normals();
    bool has_normals= (normals.size() == positions.size());
    for(unsigned int i= 0; i < (unsigned int) normals.size(); i++)
        fprintf(out, "vn %lf %lf %lf\n", double(normals[i].x), double(normals[i].y), double(normals[i].z));
    fprintf(out, "\n");
    
    const std::vector<unsigned int>& indices= mesh.indices();
    bool has_indices= (indices.size() > 0);
    unsigned int n= has_indices ? (unsigned int) indices.size() : (unsigned int) positions.size();
    for(unsigned int i= 0; i +2 < n; i+= 3)
    {
        fprintf(out, "f");
        for(unsigned int k= 0; k < 3; k++)
        {
            unsigned int id= has_indices ? indices[i+k] +1 : i+k +1;
            fprintf(out, " %u", id);
            if(has_texcoords && has_normals)
                fprintf(out, "/%u/%u", id, id);
            else if(has_texcoords)
                fprintf(out, "/%u", id);
            else if(has_normals)
                fprintf(out, "//%u", id);
        }
        fprintf(out, "\n");
    }
    
    fclose(out);
    return 0;
}

int write_mesh_data( const MeshData& data, const char *filename )
{
    if(data.positions.size() == 0)
        return -1;
    if(filename == nullptr)
        return -1;
    
    FILE *out= fopen(filename, "wt");
    if(out == NULL)
        return -1;
    
    printf("writing mesh data '%s'...\n", filename);
    
    for(int i= 0; i < int(data.positions.size()); i++)
        fprintf(out, "v %lf %lf %lf\n", double(data.positions[i].x), double(data.positions[i].y), double(data.positions[i].z));
    fprintf(out, "\n");
    
    bool has_texcoords= (data.texcoords.size() == data.positions.size());
    for(int i= 0; i < int(data.texcoords.size()); i++)
        fprintf(out, "vt %lf %lf\n", double(data.texcoords[i].x), double(data.texcoords[i].y));
    fprintf(out, "\n");
    
    bool has_normals= (data.normals.size() == data.positions.size());
    for(int i= 0; i < int(data.normals.size()); i++)
        fprintf(out, "vn %lf %lf %lf\n", double(data.normals[i].x), double(data.normals[i].y), double(data.normals[i].z));
    fprintf(out, "\n");
    
    for(int i= 0; i +2 < int(data.position_indices.size()); i+= 3)
    {
        fprintf(out, "f");
        for(int k= 0; k < 3; k++)
        {
            int id= i+k;
            fprintf(out, " %u", data.position_indices[id] +1);
            if(has_texcoords && has_normals)
                fprintf(out, "/%u/%u", data.texcoord_indices[id] +1, data.normal_indices[id] +1);
            else if(has_texcoords)
                fprintf(out, "/%u", data.texcoord_indices[id] +1);
            else if(has_normals)
                fprintf(out, "//%u", data.normal_indices[id] +1);
        }
        fprintf(out, "\n");
    }
    
    fclose(out);
    return 0;
}


MaterialLib read_materials( const char *filename )
{
    MaterialLib materials;
    
    FILE *in= fopen(filename, "rt");
    if(in == NULL)
    {
        printf("[error] loading materials '%s'...\n", filename);
        return materials;
    }
    
    printf("loading materials '%s'...\n", filename);
    
    Material *material= NULL;
    std::string path= normalize_path(pathname(filename));
    
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
        
        if(line[0] == 'n')
        {
            if(sscanf(line, "newmtl %[^\r\n]", tmp) == 1)
            {
                materials.names.push_back( tmp );
                materials.data.push_back( Material(Black()) );
                material= &materials.data.back();
            }
        }
        
        if(material == NULL)
            continue;
        
        if(line[0] == 'K')
        {
            double r, g, b;
            if(sscanf(line, "Kd %lf %lf %lf", &r, &g, &b) == 3)
                material->diffuse= Color(r, g, b);
            else if(sscanf(line, "Ks %lf %lf %lf", &r, &g, &b) == 3)
                material->specular= Color(r, g, b);
            else if(sscanf(line, "Ke %lf %lf %lf", &r, &g, &b) == 3)
                material->emission= Color(r, g, b);
        }
        
        else if(line[0] == 'N')
        {
            double n;
            if(sscanf(line, "Ns %lf", &n) == 1)          // Ns, puissance / concentration du reflet, modele blinn phong
                material->ns= n;
        }
        
        else if(line[0] == 'm')
        {
            if(sscanf(line, "map_Kd %[^\r\n]", tmp) == 1)
            {
                material->diffuse_filename= normalize_path(tmp);
                if(relative_path(material->diffuse_filename))
                    material->diffuse_filename= normalize_path(path + tmp);
                    
                // pas propre : evite de concatener 2 fois le prefixe sur les noms de fichiers des textures relatifs...
                // a reprendre : lors de la relecture du cache, ne pas re-interpreter les noms de fichiers des textures...
                if(relative_path(material->diffuse_filename))
                    material->diffuse_filename= "./" + material->diffuse_filename;
            }
            if(sscanf(line, "map_Ks %[^\r\n]", tmp) == 1)
            {
                material->specular_filename= normalize_path(tmp);
                if(relative_path(material->specular_filename))
                    material->specular_filename= normalize_path(path + tmp);
                    
                // a reprendre : lors de la relecture du cache, ne pas re-interpreter les noms de fichiers des textures...
                if(relative_path(material->specular_filename))
                    material->specular_filename= std::string("./") + material->specular_filename;
            }
        }
    }
    
    fclose(in);
    if(error)
        printf("[error] parsing line :\n%s\n", line_buffer);
    
    return materials;
}
