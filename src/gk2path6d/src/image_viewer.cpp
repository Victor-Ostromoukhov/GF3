
//! \file image_viewer.cpp permet de visualiser les images aux formats reconnus par gKit2 light bmp, jpg, tga, png, hdr, etc.

#include <cfloat>
#include <algorithm>

#include "app.h"
#include "widgets.h"

#include "image.h"
#include "image_io.h"
#include "image_hdr.h"
#include "image_exr.h"

#include "program.h"
#include "uniforms.h"
#include "texture.h"


struct ImageViewer : public App
{
    ImageViewer( std::vector<const char *>& _filenames ) : App(1024, 640), m_filenames(_filenames) {}
    
    void range( const Image& image )
    {
        constexpr int NBINS= 1000;
        
        int bins[NBINS] = {};
        float gmin= FLT_MAX;
        float gmax= 0;
        for(int y= 0; y < image.height(); y++)
        for(int x= 0; x < image.width(); x++)
        {
            Color32 color= image(x, y);
            float g= (color.r + color.g + color.b) / 3;
            //~ if(std::isinf(g) || std::isnan(g) || g == FLT_MAX || g == -FLT_MAX)
                //~ continue;
            
            if(g < gmin) gmin= g;
            if(g > gmax) gmax= g;
        }
        
        for(int y= 0; y < image.height(); y++)
        for(int x= 0; x < image.width(); x++)
        {
            Color32 color= image(x, y);
            float g= (color.r + color.g + color.b) / 3;
            //~ if(std::isinf(g) || std::isnan(g) || g == FLT_MAX || g == -FLT_MAX)
                //~ continue;
            
            int b= (g - gmin) * NBINS / (gmax - gmin);
            if(b >= NBINS) b= NBINS-1;
            if(b < 0) b= 0;
            bins[b]++;
        }
        
        float qbins= 0;
        for(int i= 0; i < NBINS; i++)
        {
            qbins= qbins + (float) bins[i] / (m_width * m_height);
            //~ if(qbins > .75f)
            if(qbins > .90f)
            {
                m_saturation= gmin + float(i + .5f) / NBINS * (gmax - gmin);
                //~ m_saturation_step= m_saturation / 40.f;
                m_saturation_step= gmin + float(i + .5f) / NBINS * (gmax - gmin) / NBINS;
                m_saturation_max= gmax;
                break;
            }
        }
        
        printf("range [%f..%f]\n", gmin, gmax);
        //~ for(int i= 0; i < NBINS; i++)
            //~ printf("%f ", ((float) bins[i] * NBINS / (m_width * m_height)));
        //~ printf("\n");
        printf("step %f\n", m_saturation_step);
        printf("saturation %f\n", m_saturation);
        float s= std::log2(1 / m_saturation);
        
        //~ m_saturation= std::log2(1 / m_saturation);
        //~ m_saturation_step= 0.05;
        m_compression= 2.2f;        
        //~ printf("saturation exponent %f (%f)\n", s, m_saturation * pow(float(2), s));
        //~ m_saturation= s;
    }
    
    void title( const int index )
    {
        char tmp[1024];
        sprintf(tmp, "buffer %02d: %s", index, m_filenames[index]);
        SDL_SetWindowTitle(m_window, tmp);        
    }
    
    int init( )
    {
        m_width= 0;
        m_height= 0;
        
    #if 0
        for(int i= 0; i < int(m_filenames.size()); i++)
        {
            printf("loading buffer %d...\n", i);
            
            Image image;
            if(is_exr_image(m_filenames[i]))
                image= read_image_exr(m_filenames[i]);
            else if(is_pfm_image(m_filenames[i]))
                image= read_image_pfm(m_filenames[i]);
            else if(is_hdr_image(m_filenames[i]))
                image= read_image_hdr(m_filenames[i]);
            else
                image= read_image(m_filenames[i]);
            
            if(image.size())
            {
                m_images.push_back(image);
                m_textures.push_back(make_texture(0, image));
                
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
                
                m_width= std::max(m_width, image.width());
                m_height= std::max(m_height, image.height());
            }
            else
            {
                m_filenames.erase(m_filenames.begin() + i);
                i= i -1;
            }
        }
    #else
        // parallelise le chargement des images
        std::vector<Image> images(m_filenames.size());
        std::vector<const char*> filenames(m_filenames.size());
        
    #pragma omp parallel for
        for(int i= 0; i < int(m_filenames.size()); i++)
        {
            if(is_exr_image(m_filenames[i]))
                images[i]= read_image_exr(m_filenames[i]);
            else if(is_pfm_image(m_filenames[i]))
                images[i]= read_image_pfm(m_filenames[i]);
            else if(is_hdr_image(m_filenames[i]))
                images[i]= read_image_hdr(m_filenames[i]);
            else
                images[i]= read_image(m_filenames[i]);
            
            if(images[i].size())
                filenames[i]= m_filenames[i];
        }
        
        m_filenames.clear();
        m_images.clear();
        m_textures.clear();
        for(int i= 0; i < int(filenames.size()); i++)
        {
            if(filenames[i] && images[i].size())
            {
                m_filenames.push_back(filenames[i]);
                
                m_images.push_back(images[i]);
                m_textures.push_back( make_texture(0, images[i]) );
                
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
                
                m_width= std::max(m_width, images[i].width());
                m_height= std::max(m_height, images[i].height());
            }
        }
    #endif
    
        if(m_images.empty())
        {
            printf("no image...\n");
            return -1;
        }
        
        assert(m_filenames.size() == m_images.size());
        assert(m_textures.size() == m_images.size());
        
        // change le titre de la fenetre
        title(0);
        
        // redminsionne la fenetre
        SDL_SetWindowSize(m_window, m_width, m_height);
        
        glGenVertexArrays(1, &m_vao);
        glBindVertexArray(m_vao);
        
        m_program= read_program( smart_path("data/shaders/tonemap.glsl") );
        program_print_errors(m_program);
        
        // 
        m_red= 1;
        m_green= 1;
        m_blue= 1;
        m_alpha= 1;
        m_gray= 0;
        m_smooth= 1;
        m_difference= 0;
        m_compression= 2.2f;
        m_saturation= 1;
        m_saturation_step= 1;
        m_saturation_max= 1000;
        m_index= 0;
        m_reference_index= -1;
        m_zoom= 4;
        m_graph= 0;
        
        // parametres d'exposition / compression
        range(m_images.back());
        
        //
        m_widgets= create_widgets();
        
        //
        glGenSamplers(1, &m_sampler_nearest);
        glSamplerParameteri(m_sampler_nearest, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glSamplerParameteri(m_sampler_nearest, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glSamplerParameteri(m_sampler_nearest, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
        glSamplerParameteri(m_sampler_nearest, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
        
        // etat openGL par defaut
        glUseProgram(0);
        glBindVertexArray(0);
        glBindTexture(GL_TEXTURE_2D, 0);
        
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);
        return 0;
    }
    
    int quit( )
    {
        glDeleteVertexArrays(1, &m_vao);
        glDeleteTextures(m_textures.size(), m_textures.data());
        
        release_program(m_program);
        release_widgets(m_widgets);
        return 0;
    }
    
    int render( )
    {
        // effacer l'image
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        if(key_state('r'))
        {
            clear_key_state('r');
            reload_program(m_program, smart_path("data/shaders/tonemap.glsl") );
            program_print_errors(m_program);
        }
        if(key_state('f'))
        {
            clear_key_state('f');
            SDL_SetWindowSize(m_window, m_width, m_height);
            glViewport(0, 0, m_width, m_height);
        }
        if(key_state('l'))
        {
            clear_key_state('l');
            m_compression= 1;
        }
        
        if(key_state(SDLK_LEFT))
        {
            clear_key_state(SDLK_LEFT);
            m_index= (m_index -1 + m_textures.size()) % m_textures.size();
            // change aussi le titre de la fenetre
            title(m_index);
        }
        
        if(key_state(SDLK_RIGHT))
        {
            clear_key_state(SDLK_RIGHT);
            m_index= (m_index +1 + m_textures.size()) % m_textures.size();
            // change aussi le titre de la fenetre
            title(m_index);
        }
        
        int xmouse, ymouse;
        unsigned int bmouse= SDL_GetMouseState(&xmouse, &ymouse);
        
        glBindVertexArray(m_vao);
        glUseProgram(m_program);
        
        // selection des buffers + filtrage
        GLuint sampler= 0;
        if(!m_smooth)
            sampler= m_sampler_nearest;
        
        program_use_texture(m_program, "image", 0, m_textures[m_index], sampler);
        if(m_reference_index == -1)
            program_use_texture(m_program, "image_next", 1, m_textures[(m_index +1) % m_textures.size()], sampler);
        else
            program_use_texture(m_program, "image_next", 1, m_textures[m_reference_index], sampler);
        
        // activer le split de l'ecran
        if(bmouse & SDL_BUTTON(1))
            program_uniform(m_program, "split", (int) xmouse);
        else
            program_uniform(m_program, "split", (int) window_width() +2);
        
        // parametres
        program_uniform(m_program, "channels", Color(m_red, m_green, m_blue, m_alpha));
        program_uniform(m_program, "gray", float(m_gray));
        program_uniform(m_program, "difference", float(m_difference));
        program_uniform(m_program, "compression", m_compression);
        program_uniform(m_program, "saturation", m_saturation);
        
        // zoom
        if(bmouse & SDL_BUTTON(3))
        {
            SDL_MouseWheelEvent wheel= wheel_event();
            if(wheel.y != 0)
            {
                m_zoom= m_zoom + float(wheel.y) / 4.f;
                if(m_zoom < .1f) m_zoom= .1f;
                if(m_zoom > 10.f) m_zoom= 10.f;
            }
        }
    
        program_uniform(m_program, "center", vec2( float(xmouse) / float(window_width()), float(window_height() - ymouse -1) / float(window_height())));
        if(bmouse & SDL_BUTTON(3))
            program_uniform(m_program, "zoom", m_zoom);
        else
            program_uniform(m_program, "zoom", 1.f);
        
        // graphes / courbes
        if(key_state('g'))
        {
            clear_key_state('g');
            m_graph= (m_graph +1) % 2;
        }

        program_uniform(m_program, "graph", int(m_graph));
        program_uniform(m_program, "line", vec2(float(window_height() - ymouse -1) / float(window_height()), float(window_height() - ymouse -1)));
        
        // dessine 1 triangle plein ecran
        glDrawArrays(GL_TRIANGLES, 0, 3);
        
        // actions
        if(key_state('c'))
        {
            clear_key_state('c');
            
            // change l'extension
            std::string file= m_filenames[m_index];
            size_t ext= file.rfind(".");
            if(ext != std::string::npos)
                file= file.substr(0, ext) + "-tone.png";
            
            printf("writing '%s'...\n", file.c_str());
            screenshot(file.c_str());
        }
        if(key_state('p'))
        {
            clear_key_state('p');
            
            // change l'extension
            std::string file= m_filenames[m_index];
            size_t ext= file.rfind(".");
            if(ext != std::string::npos)
                file= file.substr(0, ext) + "-export.pfm";
            
            printf("writing '%s'...\n", file.c_str());
            write_image_pfm(m_images[m_index], file.c_str());
        }
        if(key_state('h'))
        {
            clear_key_state('h');
            
            // change l'extension
            std::string file= m_filenames[m_index];
            size_t ext= file.rfind(".");
            if(ext != std::string::npos)
                file= file.substr(0, ext) + "-export.hdr";
            
            printf("writing '%s'...\n", file.c_str());
            write_image_hdr(m_images[m_index], file.c_str());
        }
        if(key_state('e'))
        {
            clear_key_state('e');
            
            // change l'extension
            std::string file= m_filenames[m_index];
            size_t ext= file.rfind(".");
            if(ext != std::string::npos)
                file= file.substr(0, ext) + "-export.exr";
            
            printf("writing '%s'...\n", file.c_str());
            write_image_exr(m_images[m_index], file.c_str());
        }
        
        begin(m_widgets);
            value(m_widgets, "saturation", m_saturation, 0.f, m_saturation_max*10, m_saturation_step);
            value(m_widgets, "compression", m_compression, .1f, 10.f, .1f);
        
            int reset= 0; 
            button(m_widgets, "reset", reset);
            if(reset) range(m_images[m_index]);

            int reload= 0; 
            button(m_widgets, "reload", reload);
            if(reload)
            {
                Image image;
                if(is_exr_image(m_filenames[m_index]))
                    image= read_image_exr(m_filenames[m_index]);
                else if(is_pfm_image(m_filenames[m_index]))
                    image= read_image_pfm(m_filenames[m_index]);
                else if(is_hdr_image(m_filenames[m_index]))
                    image= read_image_hdr(m_filenames[m_index]);
                else
                    image= read_image(m_filenames[m_index]);
                
                {
                    m_images[m_index]= image;
                    
                    // transfere la nouvelle version
                    glBindTexture(GL_TEXTURE_2D, m_textures[m_index]);
                    glTexImage2D(GL_TEXTURE_2D, 0,
                        GL_RGBA32F, image.width(), image.height(), 0,
                        GL_RGBA, GL_FLOAT, image.buffer());
                    
                    glGenerateMipmap(GL_TEXTURE_2D);                    
                }
            }
            
            int reference= (m_index == m_reference_index) ? 1 : 0;
            if(button(m_widgets, "reference", reference))
            {
                if(reference) m_reference_index= m_index;       // change de reference
                else m_reference_index= -1;     // deselectionne la reference
            }
        
        begin_line(m_widgets);
            button(m_widgets, "R", m_red);
            button(m_widgets, "G", m_green);
            button(m_widgets, "B", m_blue);
            button(m_widgets, "A", m_alpha);
            button(m_widgets, "gray", m_gray);
            button(m_widgets, "smooth", m_smooth);
            
            if(m_reference_index != -1)
                button(m_widgets, "diff to reference", m_difference);
            
        begin_line(m_widgets);
        {
            int x= xmouse;
            int y= window_height() - ymouse -1;
            int y1= ymouse;
            int px= float(x) / float(window_width()) * m_images[m_index].width();
            int py= float(y) / float(window_height()) * m_images[m_index].height();
            int py1= float(y1) / float(window_height()) * m_images[m_index].height();
            if(px >= 0 && px < m_images[m_index].width()
            && py >= 0 && py < m_images[m_index].height())
            {
                Color32 pixel= m_images[m_index](px, py);
                label(m_widgets, "pixel (%d %d / %d) : %f %f %f", px, py, py1, pixel.r, pixel.g, pixel.b);
            }
        }
        end(m_widgets);
        
        draw(m_widgets, window_width(), window_height());
        if(key_state('s'))
        {
            clear_key_state('s');
            
            static int calls= 0;
            screenshot("screenshot", ++calls);
            printf("screenshot %d...\n", calls);
        }
        
        if(key_state(SDLK_LCTRL) && key_state('w'))
        {
            // ferme l'image courante...
            clear_key_state('w');
            
            m_filenames.erase(m_filenames.begin() + m_index);
            m_images.erase(m_images.begin() + m_index);
            m_textures.erase(m_textures.begin() + m_index);
            if(m_reference_index == m_index)
                m_reference_index= -1;
            
            if(m_textures.empty())
                // derniere image fermee, quitter l'application
                return 0;
            
            m_index= m_index % int(m_textures.size());
            title(m_index);
        }
        return 1;
    }
    
protected:
    Widgets m_widgets;
    
    std::vector<const char *> m_filenames;
    std::vector<Image> m_images;
    std::vector<GLuint> m_textures;
    int m_width, m_height;

    GLuint m_program;
    GLuint m_vao;
    GLuint m_sampler_nearest;
    
    int m_red, m_green, m_blue, m_alpha, m_gray;
    int m_smooth;
    int m_difference;
    
    float m_compression;
    float m_saturation;
    float m_saturation_step;
    float m_saturation_max;

    float m_zoom;
    int m_index;
    int m_reference_index;
    int m_graph;
};


int main( int argc, char **argv )
{
    if(argc == 1)
    {
        printf("usage: %s image.[bmp|png|jpg|tga|hdr|pfm|exr]\n", argv[0]);
        return 0;
    }
    
    std::vector<const char *> options(argv +1, argv + argc);
    ImageViewer app(options);
    app.run();
    
    return 0;
}
