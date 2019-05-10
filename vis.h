#ifndef EULERLIVE_VIS
#define EULERLIVE_VIS

#include <SDL2/SDL.h>
#include "euler.h"

struct canvas
{
    int width;
    int height;
    int win_width;
    int win_height;
    SDL_Window *window;
    SDL_Renderer *renderer;
    SDL_Texture *texture;
    SDL_PixelFormat *fmt;
    Uint32 *pixels;
    int var;
};
typedef struct canvas canvas;

void vis_start();
void vis_quit();

void canvas_init(canvas *c, int win_width, int win_height, int width,
                    int height);
void canvas_free(canvas *c);
int canvas_handleEvents();
int canvas_handleKeyEvent(canvas *c, SDL_KeyboardEvent *key);
void canvas_update(canvas *c, domain *dom);
void buff2pixels(domain *dom, canvas *c);
void updateFrame(canvas *c, double t);
#endif
