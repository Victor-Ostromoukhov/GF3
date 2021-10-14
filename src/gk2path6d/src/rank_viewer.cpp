#include <string>

#include "app.h"
#include "widgets.h"

#include "vec.h"
#include "mat.h"
#include "orbiter.h"

#include "texture.h"
#include "mesh.h"
#include "draw.h"
#include "scene.h"

#include "brdf.h"
#include "points.h"
#include "sampler.h"
#include "sampler_sobol.h"

#include "options.h"

Options options;


class SceneViewer : public App
{
    Scene scene;
    UniformSource sources;

    Orbiter camera;
    int frame;

    Point scene_pmin;
    Point scene_pmax;

    int wireframe;

    Mesh source_edges;
    Mesh source_points;
    Mesh brdf_points;
    Mesh brdf_rays;
    Mesh source_rays;
    
public:
    SceneViewer( ) : App(options.width, options.height), frame(0), wireframe(0), m_sampler(nullptr)
    {
        if(options.scene_filename)
            // charge une description de scene, si possible
            scene= read_scene(options.scene_filename);
            
        else if(options.file_filename)
            // sinon, fabrique une scene avec un seul objet statique
            scene.attach( SceneObject(options.file_filename) );
            
        else
        {
            printf("[error] no scene...\n");
            return;
        }
        
        if(options.orbiter_filename)
            scene.camera_filename= options.orbiter_filename;
        
        if(scene.camera_filename.empty())
        {
            printf("[error] no camera...\n");
            return;
        }
        
        scene.build();
        
        //
        m_widgets= create_widgets();
        
        {
            // ensemble de sources de lumieres
            std::vector<Source> mesh_sources;
            scene.build_sources(mesh_sources);
            
            sources= UniformSource(mesh_sources);
            printf("%d sources.\n", sources.size());
            printf("%d quads.\n", int(sources.quads.size()));
            
            if(sources.size())
            {
                // construit les aretes des quads
                source_edges= Mesh(GL_LINES);
                source_points= Mesh(GL_POINTS);
                if(sources.quads.size())
                {
                    source_edges.default_color(Color(1, 0, 0));
                    source_points.default_color(Color(0.8, 0.4, 0));
                    for(auto const & quad : sources.quads)
                    {
                        source_edges.vertex(quad.a);
                        source_edges.vertex(quad.b);
                        
                        source_edges.vertex(quad.b);
                        source_edges.vertex(quad.c);
                        
                        source_edges.vertex(quad.c);
                        source_edges.vertex(quad.d);
                        
                        source_edges.vertex(quad.d);
                        source_edges.vertex(quad.a);
                    }
                }
                else
                {
                    // construit les aretes des triangles
                    source_edges.default_color(Color(0.8, 0.4, 0));
                    source_points.default_color(Color(1, 0, 0));
                    for(auto const & triangle : sources.sources)
                    {
                        source_edges.vertex(triangle.triangle.a);
                        source_edges.vertex(triangle.triangle.b);
                        
                        source_edges.vertex(triangle.triangle.b);
                        source_edges.vertex(triangle.triangle.c);
                        
                        source_edges.vertex(triangle.triangle.c);
                        source_edges.vertex(triangle.triangle.a);
                    }
                }
                
                // echantillonne les sources
                std::vector<Point> samples0;
                std::vector<Point> samples1;
                
                m_sampler= new SobolSampler(16, 1024, 0);
                
                for(unsigned int i= 0; i < m_sampler->size(); i++)
                {
                    m_sampler->index(i);
                    m_sampler->dimension(5);       // consomme les memes dimensions que path...
                    
                    Float pdf;
                    LightSample sample= sources.sample(*m_sampler, pdf);    // en fonction des options --source_2d, etc.
                    source_points.vertex(sample.s);
                    
            #if 0
                    //~ Float u1= sampler.sample1();    // selection source
                    //~ Float u2= sampler.sample2();    // position source
                    //~ Float u3= sampler.sample3();
                    //~ LightSample sample= points.sample(u1, u2, u3, pdf);
                    //~ LightSample sample= points.sample(u1, u2, pdf);
                    
                    //~ if(sample.source_id == 0)
                        //~ samples0.push_back(Point(u1, u2, u3));
                    //~ else if(sample.source_id == 1)
                        //~ samples1.push_back(Point(u1, u2, u3));
                    
                    //~ if(sample.source_id == 0)
                        //~ samples0.push_back(sample.s);
                    //~ else if(sample.source_id == 1)
                        //~ samples1.push_back(sample.s);
                }
                
                {
                    //~ FILE *out= fopen("samples0re.txt", "wt");
                    FILE *out= fopen("samples0.txt", "wt");
                    if(out)
                    {
                        for(int i= 0; i < int(samples0.size()); i++)
                            fprintf(out, "%f %f %f\n", samples0[i].x, samples0[i].y, samples0[i].z);
                        
                        fclose(out);
                    }
                }
                {
                    //~ FILE *out= fopen("samples1re.txt", "wt");
                    FILE *out= fopen("samples1.txt", "wt");
                    if(out)
                    {
                        for(int i= 0; i < int(samples1.size()); i++)
                            fprintf(out, "%f %f %f\n", samples1[i].x, samples1[i].y, samples1[i].z);
                        
                        fclose(out);
                    }
            #endif
                }
            }
        }
        
        // bbox
        scene.objects[0].frames[0].bounds(scene_pmin, scene_pmax);
        for(int i= 1; i < int(scene.objects.size()); i++)
        {
            Point bmin, bmax;
            scene.objects[i].frames[0].bounds(bmin, bmax);
            scene_pmin= min(scene_pmin, bmin);
            scene_pmax= max(scene_pmax, bmax);
        }

        // observe toute la scene
        camera.lookat(scene_pmin, scene_pmax);
    }

    int init( )
    {
        // etat openGL par defaut
        glCullFace(GL_BACK);
        glFrontFace(GL_CCW);
        //glEnable(GL_CULL_FACE); // n'affiche que les faces correctement orientees...
        glDisable(GL_CULL_FACE);    // les faces mal orientees sont affichees avec des hachures oranges...

        glDepthFunc(GL_LESS);
        glEnable(GL_DEPTH_TEST);
        glLineWidth(2);
        glPointSize(3);
        return 0;
    }

    int quit( )
    {
        delete m_sampler;
        return 0;
    }

    int render( )
    {
        if(wireframe)
        {
            glClearColor(1, 1, 1, 1);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        }
        else
        {
            glClearColor(0.2f, 0.2f, 0.2f, 1);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // recupere les mouvements de la souris
        int mx, my;
        unsigned int mb= SDL_GetRelativeMouseState(&mx, &my);
        int xmouse, ymouse;
        SDL_GetMouseState(&xmouse, &ymouse);

        // deplace la camera
        if(mb & SDL_BUTTON(1))
            camera.rotation(mx, my);      // tourne autour de l'objet
        else if(mb & SDL_BUTTON(3))
            camera.translation((float) mx / (float) window_width(), (float) my / (float) window_height()); // deplace le point de rotation
        //~ else if(mb & SDL_BUTTON(2))
            //~ camera.move(mx);           // approche / eloigne l'objet

        SDL_MouseWheelEvent wheel= wheel_event();
        if(wheel.y != 0)
        {
            clear_wheel_event();
            if(key_state(SDLK_LSHIFT))
                camera.move(1600.f * wheel.y);  // approche / eloigne l'objet
            else
                camera.move(16.f * wheel.y);  // approche / eloigne l'objet
        }

        if(key_state('w'))
        {
            clear_key_state('w');
            wireframe= !wireframe;
        }

        if(key_state('c'))
        {
            clear_key_state('c');
            camera.write_orbiter("orbiter.txt");
            if(!scene.camera_filename.empty())
                camera.write_orbiter(scene.camera_filename.c_str());
        }
        if(key_state('v'))
        {
            clear_key_state('v');
            if(!scene.camera_filename.empty())
            {
                if(camera.read_orbiter(scene.camera_filename.c_str()) < 0)
                    camera= Orbiter(scene_pmin, scene_pmax);
            }
            else
                if(camera.read_orbiter("orbiter.txt") < 0)
                    camera= Orbiter(scene_pmin, scene_pmax);

        }
        if(key_state('f'))
        {
            clear_key_state('f');
            camera= Orbiter(scene_pmin, scene_pmax);
        }
        
        // draw
        static int highlight= 1;
        if(key_state('h'))
        {
            clear_key_state('h');
            highlight= (highlight +1) % 2;
        }
        if(highlight)
        {
            for(int i= 0; i < int(scene.objects.size()); i++)
                draw(scene.objects[i].frames[0], camera);
        }
        
        // interface
        begin(m_widgets);
        int x= xmouse;
        int y= window_height() - ymouse -1;
        
        Transform view= camera.view();
        Transform projection= camera.projection(window_width(), window_height(), 45);
        Transform viewport= Viewport(window_width(), window_height());
        Transform T= Inverse(viewport * projection * view);
        Point o= camera.position();
        Point e= T(Point(x, y, 1));
        
        Raydata ray(o, e);
        scene.intersect(ray);
        Float t= 0;
        Float td= 0;
        if(ray.intersect())
        {
            t= ray.hit.t;
            td= distance(o, ray.hitp());
        }
        
        label(m_widgets, "pixel %d %d: t %f, %f", x, y, float(t), float(td));
        begin_line(m_widgets);
        
        static int sample= 0;
        button(m_widgets, "brdf", sample);
        
        static int show_rays= 1;
        button(m_widgets, "show rays", show_rays);
        static int show_shadows= 1;
        button(m_widgets, "show shadows", show_shadows);
        
        static int rays_n= 100;
        value(m_widgets, "samples", rays_n, 0, 1024, 16);
        
        if(key_state('b'))  // raccourci clavier...
        {
            clear_key_state('b');
            sample= (sample +1) % 2;
        }
        
        if(sample && ray.intersect())
        {
            std::vector<vec3> points;
            std::vector<vec3> rays;
            std::vector<vec3> shadows;
            
            const Material& material= scene.material(ray);
            Float d= material.diffuse.power();
            Float s= material.specular.power();
            if(options.force_diffuse_material) { d= 1; s= 0; }      // force matiere diffuse
            Float kd= d / (d+s);
            Float ks= s / (d+s);
            
            Blinn distribution(material.ns);
            Brdf brdf(kd, material.diffuse, ks, distribution, material.specular);
            
            Point p= ray.hitp();
            Vector pn= scene.normal(ray);
            if(dot(pn, ray.d()) > 0)
                pn= -pn;
            World world(pn);
            
            for(int i= 0; i < rays_n; i++)
            {
                m_sampler->index(i);
                m_sampler->dimension(5);
                
                Float brdf1= m_sampler->sample1();
                Float brdf2= m_sampler->sample2();
                Float source1= m_sampler->sample3();
                Float source2= m_sampler->sample4();
                
                Float brdf_pdf;
                Vector l= world(brdf.sample(brdf1, brdf2, world.inverse(-normalize(ray.d())), brdf_pdf));
                if(brdf_pdf > 0)
                {
                    Raydata bounce(p + pn * Float(.001), l, t);
                    if(scene.intersect(bounce))
                    {
                        Point q= bounce.hitp();
                        Vector qn= scene.normal(bounce);
                        if(dot(qn, bounce.d()) > 0)
                            qn= -qn;
                        
                        points.push_back(vec3(q));
                        
                        // teste aussi la visibilite d'une source
                        Float source_pdf;
                        LightSample sample= sources.sample(source1, source2, source_pdf);    // en fonction des options --source_2d, etc.
                        if(source_pdf > 0)
                        {
                            Raydata shadow(q + qn * Float(.001), sample.s + sample.sn * Float(.001));
                            if(!scene.intersect(shadow))
                            {
                                rays.push_back(vec3(p));
                                rays.push_back(vec3(q));
                                
                                shadows.push_back(vec3(q));
                                shadows.push_back(vec3(sample.s));
                            }
                        }
                    }
                }
            }
            
            label(m_widgets, "null samples %d/%d", rays_n - int(points.size()), rays_n);
            
            brdf_points.create(GL_POINTS);
            brdf_points.default_color(Red());
            brdf_points.positions(points);
            
            brdf_rays.create(GL_LINES);
            brdf_rays.default_color(Color(0.8, 0.4, 0));
            brdf_rays.positions(rays);
            
            source_rays.create(GL_LINES);
            source_rays.default_color(Color(0.4, 0.2, 0));
            source_rays.positions(shadows);
            
        }
        end(m_widgets);
        
        draw(m_widgets, window_width(), window_height());
        
        // source
        if(source_edges.vertex_count())
        {
            glDisable(GL_DEPTH_TEST);
            draw(source_edges, camera);
            draw(source_points, camera);
            glEnable(GL_DEPTH_TEST);
        }
        
        // brdf
        if(brdf_points.vertex_count())
        {
            glDisable(GL_DEPTH_TEST);
            draw(brdf_points, camera);
            glEnable(GL_DEPTH_TEST);
        }
        
        // brdf_rays
        if(show_rays && brdf_rays.vertex_count())
        {
            draw(brdf_rays, camera);
        }
        
        // shadow rays
        if(show_shadows && source_rays.vertex_count())
        {
            draw(source_rays, camera);
        }
        
        if(key_state('s'))
        {
            clear_key_state('s');
            static int calls= 1;
            printf("screenshot %d...\n", calls);
            screenshot("scene", calls++);
        }

        return 1;
    }
    
protected:
    Widgets m_widgets;
    SobolSampler *m_sampler;
};



int main( int argc, char *argv[] )
{
    if(argc == 1)
    {
        printf("usage: %s scene.txt\n", argv[0]);
        return 0;
    }

    options= Options((const char **) argv, (const char **) argv + argc);
    
    SceneViewer viewer;
    viewer.run();

    return 0;
}
