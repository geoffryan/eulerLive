#include <stdio.h>
#include "vis.h"


void vis_start()
{
   SDL_Init(SDL_INIT_TIMER|SDL_INIT_VIDEO|SDL_INIT_EVENTS); 
}

void canvas_init(canvas *c, int win_width, int win_height, int width,
                 int height)
{
    c->width = width;
    c->height = height;
    c->win_width = win_width;
    c->win_height = win_height;
    SDL_CreateWindowAndRenderer(win_width, win_height, SDL_WINDOW_BORDERLESS,
                                &(c->window), &(c->renderer));
    SDL_RenderSetLogicalSize(c->renderer, width, height);
    SDL_SetRenderDrawColor(c->renderer, 0, 0, 0, 255);
    SDL_RenderClear(c->renderer);
    SDL_RenderPresent(c->renderer);

    Uint32 pixFormatCode = SDL_PIXELFORMAT_ARGB8888;

    c->texture = SDL_CreateTexture(c->renderer, pixFormatCode,
                                    SDL_TEXTUREACCESS_STREAMING,
                                    width, height);
    c->fmt = SDL_AllocFormat(pixFormatCode);

    c->pixels = (Uint32 *)malloc(width*height*sizeof(Uint32));
    c->var = 0;
}

void canvas_free(canvas *c)
{
    free(c->pixels);
    SDL_FreeFormat(c->fmt);
    SDL_Quit();
}

void vis_quit()
{
    SDL_Quit();
}

int canvas_handleEvents(canvas *c)
{
    SDL_Event e;
    int code = 0;
    while(SDL_PollEvent(&e))
    {
        if(e.type == SDL_QUIT)
            code = 1;
        else if(e.type == SDL_KEYDOWN)
            code = canvas_handleKeyEvent(c, &e.key);
    }

    return code;
}

int canvas_handleKeyEvent(canvas *c, SDL_KeyboardEvent *key)
{
    int code = 0;

    if(key->keysym.sym == SDLK_q || key->keysym.sym == SDLK_ESCAPE)
        code = 1;
    else if(key->keysym.sym == SDLK_t)
        code = 0;
    else if(key->keysym.sym == SDLK_1)
        c->var = 0;
    else if(key->keysym.sym == SDLK_2)
        c->var = 1;
    else if(key->keysym.sym == SDLK_3)
        c->var = 2;
    else if(key->keysym.sym == SDLK_4)
        c->var = 3;
    else if(key->keysym.sym == SDLK_5)
        c->var = 4;

    return code;
}

void canvas_update(canvas *c, domain *dom)
{
    buff2pixels(dom, c);

    SDL_UpdateTexture(c->texture, NULL, c->pixels, c->width * sizeof(Uint32));
    SDL_RenderClear(c->renderer);
    SDL_RenderCopy(c->renderer, c->texture, NULL, NULL);
    SDL_RenderPresent(c->renderer);
}

void buff2pixels(domain *dom, canvas *c)
{
    int Ng = dom->Ng;
    int Nx = dom->Nx;
    int Ny = dom->Ny;
    int N = Nx*Ny;

    double *buff = dom->prim;
    if(c->var >= 0 && c->var < dom->Nq)
        buff += c->var * N;

    //Find max/min
    double min = 1.0e100;
    double max = -min;
    
    int i, j;
    for(j=Ng; j<Ny-Ng; j++)
        for(i=Ng; i<Nx-Ng; i++)
        {
            int k = Nx*j + i;
            min = buff[k] < min ? buff[k] : min;
            max = buff[k] > max ? buff[k] : max;
        }

    //printf("var %d, max: %.3e min: %.3e\n", c->var, max, min);

    int ilen = Nx-2*Ng;
    int jlen = Ny-2*Ng;
    if(ilen > c->width)
        ilen = c->width;
    if(jlen > c->height)
        jlen = c->height;

    double idv = 255 / (max-min);

    Uint8 rmin = 255;
    Uint8 rmax = 0;
    
    for(j=0; j<jlen; j++)
    {
        int offPix = (jlen-j-1)*ilen;
        int offBuff = Nx*(j+Ng) + Ng;
        for(i=0; i<ilen; i++)
        {
            Uint8 r = (Uint8) (idv*(buff[offBuff+i]-min));
            c->pixels[offPix + i] = SDL_MapRGB(c->fmt, r, 0, 0);
            rmin = r < rmin ? r : rmin;
            rmax = r > rmax ? r : rmax;
        }
    }
    //printf("  maxr: %hhu minr: %hhu\n", rmax, rmin);

}

void updateFrame(canvas *c, double t)
{
    int i, j;
    double om = 5.0;
    int h = c->height;
    int w = c->width;

    for(i=0; i<h; i++)
        for(j=0; j<w; j++)
        {
            double r1 = sqrt(i*i + j*j);
            double r2 = sqrt((i-h)*(i-h) + j*j);
            double r3 = sqrt(i*i + (j-w)*(j-w));
            double r4 = sqrt((i-h)*(i-h) + (j-w)*(j-w));
            Uint8 r = (Uint8) (256*(0.5*(sin(2*M_PI*r1/50 - om*t)+1)));
            Uint8 g = (Uint8) (256*(0.5*(sin(2*M_PI*r2/50 - om*t)+1)));
            Uint8 b = (Uint8) (256*(0.5*(sin(2*M_PI*r3/50 - om*t)+1)));
            double k = (0.5*(sin(2*M_PI*r4/50 - om*t)+1));
            c->pixels[w*i+j] = SDL_MapRGB(c->fmt, k*r, k*g, k*b);
        }

    SDL_UpdateTexture(c->texture, NULL, c->pixels, w * sizeof(Uint32));
    SDL_RenderClear(c->renderer);
    SDL_RenderCopy(c->renderer, c->texture, NULL, NULL);
    SDL_RenderPresent(c->renderer);
}
