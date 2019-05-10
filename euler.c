#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "euler.h"


void domain_init(domain *dom, int Nx, int Ny, int Ng, double xa, double xb,
                 double ya, double yb, int Npassive, double PLM, double CFL)
{
    int Nq = 4;
    if(Npassive > 0)
        Nq += Npassive;

    int Nx_glob = Nx+2*Ng;
    int Ny_glob = Ny+2*Ng;

    dom->Nq = Nq;
    dom->Ng = Ng;
    dom->Nx = Nx_glob;
    dom->Ny = Ny_glob;
    dom->xa = xa;
    dom->xb = xb;
    dom->ya = ya;
    dom->yb = yb;
    dom->dx = (xb-xa)/Nx;
    dom->dy = (yb-ya)/Ny;

    dom->aRHO = RHO * Nx_glob*Ny_glob;
    dom->aPPP = PPP * Nx_glob*Ny_glob;
    dom->aUXX = UXX * Nx_glob*Ny_glob;
    dom->aUYY = UYY * Nx_glob*Ny_glob;

    dom->aDDD = DDD * Nx_glob*Ny_glob;
    dom->aTAU = TAU * Nx_glob*Ny_glob;
    dom->aSXX = SXX * Nx_glob*Ny_glob;
    dom->aSYY = SYY * Nx_glob*Ny_glob;

    dom->PLM = PLM;
    dom->CFL = CFL;

    dom->t = 0.0;

    dom->prim = (double *)malloc(Nq*Nx_glob*Ny_glob * sizeof(double));
    dom->cons = (double *)malloc(Nq*Nx_glob*Ny_glob * sizeof(double));
    dom->RKcons = (double *)malloc(Nq*Nx_glob*Ny_glob * sizeof(double));
}

void domain_free(domain *dom)
{
    free(dom->prim);
    free(dom->cons);
    free(dom->RKcons);
}

void domain_initialize(domain *dom,
                        void (*init_func)(double,double,double*,int))
{
    int Ng = dom->Ng;
    int Nx = dom->Nx;
    int Ny = dom->Ny;

    int i, j;
    for(j=Ng; j<Ny-Ng; j++)
    {
        double y = dom->ya + (j-Ng)*dom->dy;
        for(i=Ng; i<Nx-Ng; i++)
        {
            double x = dom->xa + (i-Ng)*dom->dx;
            int k = j*Nx + i;

            double prim[dom->Nq];
            init_func(x, y, prim, dom->Nq);
            int q;
            for(q=0; q<dom->Nq; q++)
                dom->prim[q*Nx*Ny + k] = prim[q];
        }
    }
    domain_boundary(dom);
    domain_prim2cons(dom);
}

void domain_prim2cons(domain *dom)
{
    double *rho = dom->prim+dom->aRHO;
    double *P = dom->prim+dom->aPPP;
    double *ux = dom->prim+dom->aUXX;
    double *uy = dom->prim+dom->aUYY;
    double *den = dom->cons+dom->aDDD;
    double *E = dom->cons+dom->aTAU;
    double *sx = dom->cons+dom->aSXX;
    double *sy = dom->cons+dom->aSYY;

    int Ng = dom->Ng;
    int Nx = dom->Nx;
    int Ny = dom->Ny;

    int i, j;
    for(j=Ng; j<Ny-Ng; j++)
    {
        for(i=Ng; i<Nx-Ng; i++)
        {
            int k = j*Nx + i;

            double v2 = ux[k]*ux[k] + uy[k]*uy[k];

            den[k] = rho[k];
            sx[k] = rho[k]*ux[k];
            sy[k] = rho[k]*uy[k];
            E[k] = 0.5*rho[k]*v2 + 1.5*P[k];

            int q;
            for(q=4; q<dom->Nq; q++)
                dom->cons[q*Nx*Ny + k] = rho[k]*dom->prim[q*Nx*Ny+k];
        }
    }
}

void domain_cons2prim(domain *dom)
{
    double *rho = dom->prim+dom->aRHO;
    double *P = dom->prim+dom->aPPP;
    double *ux = dom->prim+dom->aUXX;
    double *uy = dom->prim+dom->aUYY;
    double *den = dom->cons+dom->aDDD;
    double *E = dom->cons+dom->aTAU;
    double *sx = dom->cons+dom->aSXX;
    double *sy = dom->cons+dom->aSYY;

    int Ng = dom->Ng;
    int Nx = dom->Nx;
    int Ny = dom->Ny;

    int i, j;
    for(j=Ng; j<Ny-Ng; j++)
    {
        for(i=Ng; i<Nx-Ng; i++)
        {
            //printf("c2p: %d %d\n", i, j);
            int k = j*Nx + i;

            //printf("P0: %.6lg %.6lg %.6lg %.6lg\n", rho[k],P[k],ux[k],uy[k]);
            //printf("C:  %.6lg %.6lg %.6lg %.6lg\n", den[k],E[k],sx[k],sy[k]);

            rho[k] = den[k];
            ux[k] = sx[k] / den[k];
            uy[k] = sy[k] / den[k];

            double v2 = ux[k]*ux[k] + uy[k]*uy[k];

            P[k] = 2.0*(E[k] - 0.5*den[k]*v2)/3.0;
            
            //printf("P1: %.6lg %.6lg %.6lg %.6lg\n", rho[k],P[k],ux[k],uy[k]);
            
            int q;
            for(q=4; q<dom->Nq; q++)
                dom->prim[q*Nx*Ny + k] = dom->cons[q*Nx*Ny+k]/den[k];
        }
    }
}

double domain_dt(domain *dom)
{
    double *rho = dom->prim+dom->aRHO;
    double *P = dom->prim+dom->aPPP;
    double *ux = dom->prim+dom->aUXX;
    double *uy = dom->prim+dom->aUYY;

    int Ng = dom->Ng;
    int Nx = dom->Nx;
    int Ny = dom->Ny;
    double dx = dom->dx;
    double dy = dom->dy;

    double dt = 1.0e100;

    int i, j;
    for(j=Ng; j<Ny-Ng; j++)
    {
        for(i=Ng; i<Nx-Ng; i++)
        {
            int k = j*Nx + i;

            double cs = sqrt(5*P[k]/(3*rho[k]));
            double dtx = dx / (fabs(ux[k]) + cs);
            double dty = dy / (fabs(uy[k]) + cs);

            if(dtx < dt)
                dt = dtx;
            if(dty < dt)
                dt = dty;
        }
    }

    return dom->CFL * dt;
}

double minmod(double a, double b, double c, double PLM)
{
    if(a*b <= 0.0 || b*c <= 0.0)
        return 0.0;

    double fa = fabs(PLM*a);
    double fb = fabs(b);
    double fc = fabs(PLM*c);

    if(fa < fb && fa < fc)
        return PLM*a;
    else if(fc < fa && fc < fb)
        return PLM*c;
    return b;
}

void wavespeeds_x(double primL[], double primR[], double *vm, double *vp)
{
    double csL = sqrt(5*primL[PPP]/(3*primL[RHO]));
    double csR = sqrt(5*primL[PPP]/(3*primL[RHO]));

    double cLm = primL[UXX] - csL;
    double cLp = primL[UXX] + csL;
    double cRm = primR[UXX] - csR;
    double cRp = primR[UXX] + csR;

    *vm = cLm < cRm ? cLm : cRm;
    *vp = cLp > cRp ? cLp : cRp;
}

void wavespeeds_y(double primL[], double primR[], double *vm, double *vp)
{
    double csL = sqrt(5*primL[PPP]/(3*primL[RHO]));
    double csR = sqrt(5*primL[PPP]/(3*primL[RHO]));

    double cLm = primL[UYY] - csL;
    double cLp = primL[UYY] + csL;
    double cRm = primR[UYY] - csR;
    double cRp = primR[UYY] + csR;

    *vm = cLm < cRm ? cLm : cRm;
    *vp = cLp > cRp ? cLp : cRp;
}

void prim2cons(double *prim, double *U, int Nq)
{
    double v2 = prim[UXX]*prim[UXX] + prim[UYY]*prim[UYY];
    U[DDD] = prim[RHO];
    U[SXX] = prim[RHO] * prim[UXX];
    U[SYY] = prim[RHO] * prim[UYY];
    U[TAU] = (0.5*prim[RHO]*v2 + 1.5*prim[PPP]);
    int q;
    for(q=4; q<Nq; q++)
        U[q] = prim[q] * U[DDD];
}

void flux_x(double *prim, double *F, int Nq)
{
    double v2 = prim[UXX]*prim[UXX] + prim[UYY]*prim[UYY];
    F[DDD] = prim[RHO] * prim[UXX];
    F[SXX] = prim[RHO] * prim[UXX]*prim[UXX] + prim[PPP];
    F[SYY] = prim[RHO] * prim[UXX]*prim[UYY];
    F[TAU] = (0.5*prim[RHO]*v2 + 2.5*prim[PPP])*prim[UXX];
    int q;
    for(q=4; q<Nq; q++)
        F[q] = prim[q] * F[DDD];
}

void flux_y(double *prim, double *F, int Nq)
{
    double v2 = prim[UXX]*prim[UXX] + prim[UYY]*prim[UYY];
    F[DDD] = prim[RHO] * prim[UYY];
    F[SXX] = prim[RHO] * prim[UXX]*prim[UYY];
    F[SYY] = prim[RHO] * prim[UYY]*prim[UYY] + prim[PPP];
    F[TAU] = (0.5*prim[RHO]*v2 + 2.5*prim[PPP])*prim[UYY];
    int q;
    for(q=4; q<Nq; q++)
        F[q] = prim[q] * F[DDD];
}

void domain_flux_x(domain *dom, double dt)
{
    int Ng = dom->Ng;
    int Nx = dom->Nx;
    int Ny = dom->Ny;
    int Nq = dom->Nq;
    int N = Nx*Ny;
    double dx = dom->dx;
    double dy = dom->dy;

    double *prim = dom->prim;
    double *cons = dom->cons;

    int i, j;
    const int LL = -1;
    const int L = 0;
    const int R = 1;
    const int RR = +2;
    for(j=Ng; j<Ny-Ng; j++)
    {
        for(i=Ng-1; i<Nx-Ng; i++)
        {
            int k = j*Nx + i;

            int q;
            double pL[Nq];
            double pR[Nq];
            double g;

            for(q=0; q<Nq; q++)
            {
                int kq = q*N+k;

                g = minmod(prim[kq+L]-prim[kq+LL],
                            0.5*(prim[kq+R]-prim[kq+LL]),
                            prim[kq+R]-prim[kq+L], dom->PLM);
                pL[q] = prim[kq+L] + 0.5*g;
                g = minmod(prim[kq+R]-prim[kq+L],
                            0.5*(prim[kq+RR]-prim[kq+L]),
                            prim[kq+RR]-prim[kq+R], dom->PLM);
                pR[q] = prim[kq+R] - 0.5*g;
            }

            double F[Nq];
            double vL, vR;
            wavespeeds_x(pL, pR, &vL, &vR);

            if(vL > 0)
                flux_x(pL, F, Nq);
            else if(vR < 0)
                flux_x(pR, F, Nq);
            else
            {
                double FL[Nq], FR[Nq], UL[Nq], UR[Nq];
                flux_x(pL, FL, Nq);
                flux_x(pR, FR, Nq);
                prim2cons(pL, UL, Nq);
                prim2cons(pR, UR, Nq);

                for(q=0; q<Nq; q++)
                    F[q] = (vR*FL[q] - vL*FR[q] + vL*vR*(UR[q]-UL[q]))
                            / (vR - vL);
            }

            for(q=0; q<Nq; q++)
            {
                int kq = q*N+k;
                dom->cons[kq+L] -= dt*F[q]/dx;
                dom->cons[kq+R] += dt*F[q]/dx;
            }
        }
    }
}

void domain_flux_y(domain *dom, double dt)
{
    int Ng = dom->Ng;
    int Nx = dom->Nx;
    int Ny = dom->Ny;
    int Nq = dom->Nq;
    int N = Nx*Ny;
    double dx = dom->dx;
    double dy = dom->dy;

    double *prim = dom->prim;
    double *cons = dom->cons;

    int i, j;
    const int LL = -Nx;
    const int L = 0;
    const int R = Nx;
    const int RR = +2*Nx;
    for(j=Ng-1; j<Ny-Ng; j++)
    {
        for(i=Ng; i<Nx-Ng; i++)
        {
            int k = j*Nx + i;

            int q;
            double pL[Nq];
            double pR[Nq];
            double g;

            for(q=0; q<Nq; q++)
            {
                int kq = q*N+k;

                g = minmod(prim[kq+L]-prim[kq+LL],
                            0.5*(prim[kq+R]-prim[kq+LL]),
                            prim[kq+R]-prim[kq+L], dom->PLM);
                pL[q] = prim[kq+L] + 0.5*g;
                g = minmod(prim[kq+R]-prim[kq+L],
                            0.5*(prim[kq+RR]-prim[kq+L]),
                            prim[kq+RR]-prim[kq+R], dom->PLM);
                pR[q] = prim[kq+R] - 0.5*g;
            }

            double F[Nq];
            double vL, vR;
            wavespeeds_y(pL, pR, &vL, &vR);

            if(vL > 0)
                flux_y(pL, F, Nq);
            else if(vR < 0)
                flux_y(pR, F, Nq);
            else
            {
                double FL[Nq], FR[Nq], UL[Nq], UR[Nq];
                flux_y(pL, FL, Nq);
                flux_y(pR, FR, Nq);
                prim2cons(pL, UL, Nq);
                prim2cons(pR, UR, Nq);

                for(q=0; q<Nq; q++)
                    F[q] = (vR*FL[q] - vL*FR[q] + vL*vR*(UR[q]-UL[q]))
                            / (vR - vL);
            }

            for(q=0; q<Nq; q++)
            {
                int kq = q*N+k;
                dom->cons[kq+L] -= dt*F[q]/dy;
                dom->cons[kq+R] += dt*F[q]/dy;
            }
        }
    }
}

void domain_source(domain *dom, double dt)
{
    int Ng = dom->Ng;
    int Nx = dom->Nx;
    int Ny = dom->Ny;
    int Nq = dom->Nq;
    int N = Nx*Ny;

    double *rho = dom->prim+dom->aRHO;
    double *P = dom->prim+dom->aPPP;
    double *ux = dom->prim+dom->aUXX;
    double *uy = dom->prim+dom->aUYY;
    double *den = dom->cons+dom->aDDD;
    double *E = dom->cons+dom->aTAU;
    double *sx = dom->cons+dom->aSXX;
    double *sy = dom->cons+dom->aSYY;

    const double g = -9.81;
    const double L = 1.0e3;
    const double R = 10;
    const double A = M_PI*L / (4*R*R*(0.5*M_PI-1));
    const double tend = 10.0;

    int i, j;
    for(j=Ng; j<Ny-Ng; j++)
    {
        double y = dom->ya + (j-Ng)*dom->dy;
        for(i=Ng; i<Nx-Ng; i++)
        {
            int k = j*Nx + i;
            double x = dom->xa + (i-Ng)*dom->dx;

            sy[k] += g*rho[k] * dt;
            E[k] += rho[k]*uy[k]*g * dt;

            if(x*x + y*y < R*R && dom->t < tend)
                E[k] += A*cos(M_PI*sqrt(x*x+y*y)/(2*R));
        }
    }
}

void domain_boundary(domain *dom)
{
    int Ng = dom->Ng;
    int Nx = dom->Nx;
    int Ny = dom->Ny;
    int Nq = dom->Nq;
    int N = Nx*Ny;

    double *ux = dom->prim+dom->aUXX;
    double *uy = dom->prim+dom->aUYY;

    int i,j,q;
    for(q=0; q<Nq; q++)
    {
        for(j=0; j<Ng; j++)
        {
            int jr = 2*Ng-j-1;
            for(i=0; i<Nx; i++)
            {
                int k = j*Nx + i;
                int kr = jr*Nx + i;
                dom->prim[q*N+k] = dom->prim[q*N+kr];
            }
        }
    }
    for(j=0; j<Ng; j++)
        for(i=0; i<Nx; i++)
            uy[j*Nx+i] *= -1;

    for(q=0; q<Nq; q++)
    {
        for(j=Ny-Ng; j<Ny; j++)
        {
            int jr = 2*(Ny-Ng)-j-1;
            for(i=0; i<Nx; i++)
            {
                int k = j*Nx + i;
                int kr = jr*Nx + i;
                dom->prim[q*N+k] = dom->prim[q*N+kr];
            }
        }
    }
    for(j=Ny-Ng; j<Ny; j++)
        for(i=0; i<Nx; i++)
            uy[j*Nx+i] *= -1;

    for(q=0; q<Nq; q++)
    {
        for(j=0; j<Ny; j++)
        {
            for(i=0; i<Ng; i++)
            {
                int ir = 2*Ng-i-1;
                int k = j*Nx + i;
                int kr = j*Nx + ir;
                dom->prim[q*N+k] = dom->prim[q*N+kr];
            }
        }
    }
    for(j=0; j<Ny; j++)
        for(i=0; i<Ng; i++)
            ux[j*Nx+i] *= -1;

    for(q=0; q<Nq; q++)
    {
        for(j=0; j<Ny; j++)
        {
            for(i=Nx-Ng; i<Nx; i++)
            {
                int ir = 2*(Nx-Ng)-i-1;
                int k = j*Nx + i;
                int kr = j*Nx + ir;
                dom->prim[q*N+k] = dom->prim[q*N+kr];
            }
        }
    }
    for(j=0; j<Ny; j++)
        for(i=Nx-Ng; i<Nx; i++)
            ux[j*Nx+i] *= -1;
}

void domain_substep(domain *dom, double dt)
{
    domain_flux_x(dom, dt);
    domain_flux_y(dom, dt);
    domain_source(dom, dt);
    domain_cons2prim(dom);
    domain_boundary(dom);
}

void domain_step(domain *dom)
{
    double dt = domain_dt(dom);
    printf("t: %.3le, dt: %.3le\n", dom->t, dt);
    domain_substep(dom, dt);
    dom->t += dt;
}

void atmo(double x, double y, double prim[], int Nq)
{
    const double cs0 = 343;
    const double rho0 = 1.225;
    const double g = 9.81;

    double cs20 = cs0*cs0;
    double cs2 = cs20 - 2*g*y/(3.0);

    prim[RHO] = rho0 * pow(cs2/cs20, 1.5);
    prim[PPP] = 0.6 * prim[RHO] * cs2;
    prim[UXX] = 0.0;
    prim[UYY] = 0.0;

    int q;
    for(q=4; q<Nq; q++)
        prim[q] = 0.0;
}


