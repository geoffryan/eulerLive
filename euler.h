#ifndef EULERLIVE_EULER
#define EULERLIVE_EULER

enum{RHO, PPP, UXX, UYY};
enum{DDD, TAU, SXX, SYY};

struct domain
{
    int Nx;
    int Ny;
    int Ng;
    int Nq;
    double xa;
    double xb;
    double ya;
    double yb;
    double dx;
    double dy;

    int aRHO;
    int aPPP;
    int aUXX;
    int aUYY;

    int aDDD;
    int aTAU;
    int aSXX;
    int aSYY;

    double t;
    double PLM;
    double CFL;

    double *prim;
    double *cons;
    double *RKcons;
};
typedef struct domain domain;


void domain_init(domain *dom, int Nx, int Ny, int Ng, double xa, double xb,
                 double ya, double yb, int Npassive, double PLM, double CFL);
void domain_free(domain *dom);

void domain_initialize(domain *dom,
                        void (*init_func)(double,double,double*,int));
void domain_prim2cons(domain *dom);
void domain_cons2prim(domain *dom);
double domain_dt(domain *dom);
void domain_flux_x(domain *dom, double dt);
void domain_flux_y(domain *dom, double dt);
void domain_source(domain *dom, double dt);
void domain_boundary(domain *dom);
void domain_substep(domain *dom, double dt);
void domain_step(domain *dom);


double minmod(double a, double b, double c, double PLM);
void wavespeeds_x(double primL[], double primR[], double *vm, double *vp);
void wavespeeds_y(double primL[], double primR[], double *vm, double *vp);
void prim2cons(double prim[], double U[], int Nq);
void flux_x(double prim[], double F[], int Nq);
void flux_y(double prim[], double F[], int Nq);

void atmo(double x, double y, double prim[], int Nq);

#endif
