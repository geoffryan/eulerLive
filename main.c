#include <stdio.h>
#include <stdlib.h>
#include "euler.h"
#include "vis.h"

int main(int argc, char *argv[])
{
    int Nx = 150;
    int Ny = 300;
    int win_width = 800;
    int win_height = 800;

    canvas c;
    vis_start();
    canvas_init(&c, win_width, win_height, Nx, Ny);

    double PLM = 1.5;
    double CFL = 0.1;

    domain dom;
    domain_init(&dom, Nx, Ny, 2, -100, 100, 0.0, 400, 1, PLM, CFL);
    domain_initialize(&dom, &atmo);

    double tfinal = 10.0;


    int ret = 2;
    canvas_update(&c, &dom);
    while(dom.t < tfinal)
    {
        if(ret == 0)
            domain_step(&dom);
        ret = canvas_handleEvents(&c);
        if(ret == 1)
            break;
        canvas_update(&c, &dom);
    }

    canvas_free(&c);
    vis_quit();
    domain_free(&dom);

    return 0;
}
