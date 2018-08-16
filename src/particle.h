#ifndef PARTICLE_H
#define PARTICLE_H

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <gmp.h>

#ifndef ACCU_TEST
#define ACCU_TEST
static mpf_t G;
#endif

struct particle {
    double r;
    double m;
    double locx;
    double locy;
    bool is_large;
    double vx;
    double vy;
    double fx;
    double fy;

#ifdef ACCU_TEST
    mpf_t m_r;
    mpf_t m_m;
    mpf_t m_locx;
    mpf_t m_locy;
    mpf_t m_vx;
    mpf_t m_vy;
    mpf_t m_fx;
    mpf_t m_fy;
#endif

};


static inline void p_print(struct particle p)
{
#ifdef ACCU_TEST
    gmp_printf("Particle, %Ffkg, %Ffm @(%Ff, %Ff), v(%Ff,%Ff), f(%Ff, %Ff)\n",
            p.m_m, p.m_r, p.m_locx, p.m_locy, p.m_vx, p.m_vy, p.m_fx, p.m_fy);
#else
    printf("Particle: %fkg, %fm @(%f, %f), v=(%lf,%lf), f=(%lf, %lf)\n", p.m, p.r, p.locx, p.locy, p.vx, p.vy, p.fx, p.fy);
#endif
}

//Final a particle's locatin based on v and f.
void p_final_loc(struct particle*, double);
void make_sp(struct particle*, int, double, double);
//Update force on each other 
void update_force(struct particle*, struct particle*, int, int, int);

//Get collision time
double p_collide_time(struct particle*, struct particle*);

//Model Collision
void p_collide(struct particle*, struct particle*);

static inline void init_G(){
    mpf_init_set_str(G, "6.674E-11", 10);
}

static inline double distance(double x, double y, double a, double b) 
{
   return sqrt(pow(x-a, 2) + pow(y-b, 2)); 
}
static inline double angle(double x, double y, double a, double b) 
{
    return atan2(b-y, a-x);
}
#endif
