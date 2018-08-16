#include "particle.h"
#include "util.h"
#include <stdlib.h>
#include <math.h>

#define GRAVITY_G 6.674E-11
static inline double v_post_collide(double, double, double, double, double);
static inline double dot_product(double, double, double, double);
static inline double length(double,double, double, double);

void make_sp(struct particle* p, int size, double m, double r)
{
    p->m = m;
    p->r = r;
    p->locx = rand() % size;
    p->locy = rand() % size;
    p->vx = 0;
    p->vy = 0;
    p->fx = p->fy = 0;
    p->is_large = false;
}

void p_final_loc(struct particle* p, double d_t)
{
#ifdef ACCU_TEST
    mpf_t ax, ay,sx, sy, point_five, t, vtx, vty, t2, axt2, ayt2, dvx, dvy;
    mpf_inits(ax, ay, sx, sy, point_five, t, vtx, vty, t2, axt2, ayt2, dvx, dvy,NULL);

    mpf_set_d(point_five, 0.5);
    mpf_set_d(t, d_t);

    mpf_div(ax, p->m_fx, p->m_m);
    mpf_div(ay, p->m_fy, p->m_m);
    
    mpf_pow_ui(t2, t, 2);

    mpf_mul(vtx, p->m_vx, t);
    mpf_mul(axt2, ax, t2);
    mpf_mul(axt2, axt2, point_five);
    mpf_add(sx, axt2, vtx);
    
    mpf_mul(vty, p->m_vy, t);
    mpf_mul(ayt2, ay, t2);
    mpf_mul(ayt2, ayt2, point_five);
    mpf_add(sy, ayt2, vty);
    
    mpf_add(p->m_locx, p->m_locx, sx);
    mpf_add(p->m_locy, p->m_locy, sy);

    mpf_mul(dvx, ax, t);
    mpf_mul(dvy, ay, t);

    mpf_add(p->m_vx, p->m_vx, dvx);
    mpf_add(p->m_vy, p->m_vy, dvy);

    mpf_set_d(p->m_fx, 0);
    mpf_set_d(p->m_fy, 0);

    mpf_clears(ax, ay, sx, sy, point_five, t, vtx, vty, t2, axt2, ayt2, dvx, dvy,NULL);
#else
    double ax, ay;

    ax = p->fx / p->m;
    ay = p->fy / p->m;

    p->locx+=DISTANCE(p->vx, ax, d_t); 
    p->locy+=DISTANCE(p->vy, ay, d_t);

    //calculate new velocity
    p->vx += ax * d_t;
    p->vy += ay * d_t;

    //reset the force 
    p->fx = p->fy = 0;
#endif
}

void update_force(struct particle* p1, struct particle* p2, int myid, int x_off, int y_off)
{

#ifdef ACCU_TEST
    double dx_d, dy_d, theta_d;
    mpf_t d, theta, dx, dy, dx2, dy2, d2,r_sum, m2, f, fx, fy, m_cos, m_sin;
    mpf_inits(d, theta, dx, dy, dx2, dy2,d2, r_sum,m2,f,fx,fy, m_cos, m_sin, NULL);

    //Calculate Distance
    mpf_sub(dx, p2->m_locx, p1->m_locx);
    mpf_sub(dy, p2->m_locy, p1->m_locy);

    mpf_pow_ui(dx2, dx, 2);
    mpf_pow_ui(dy2, dy, 2);

    mpf_add(d2, dx2, dy2);
    mpf_sqrt(d, d2);

    //Calcualte Angle
    dx_d = mpf_get_d(dx);
    dy_d = mpf_get_d(dy);

    theta_d = atan2(dy_d, dx_d);
    mpf_set_d(theta, theta_d);

    mpf_add(r_sum, p1->m_r, p2->m_r);

    if(mpf_cmp(d, r_sum) <= 0) 
        goto clean_up;

    mpf_mul(m2, p1->m_m, p2->m_m);
    mpf_div(f, m2, d2);

    mpf_set_d(m_cos, cos(theta_d));
    mpf_set_d(m_sin, sin(theta_d));

    mpf_mul(fx, m_cos, f);
    mpf_mul(fx, fx, G);

    mpf_mul(fy, m_sin, f);
    mpf_mul(fy, fy, G);


    mpf_add(p1->m_fx, p1->m_fx, fx);
    mpf_add(p1->m_fy, p1->m_fy, fy);

    mpf_sub(p2->m_fx, p2->m_fx, fx);
    mpf_sub(p2->m_fy, p2->m_fy, fy);

clean_up:
    mpf_clears(d, theta, dx, dy, dx2, dy2,d2, r_sum,m2,f,fx,fy, m_cos, m_sin, NULL);

#else 
    double f, d, theta, fx, fy;

    d =  distance(p1->locx + x_off, p1->locy+y_off, p2->locx, p2->locy);
    theta = angle(p1->locx + x_off, p1->locy+y_off, p2->locx, p2->locy); 
    //Prevent division by zero
    if(d <= (p1->r + p2->r)) return;

    f =  p1->m * p2->m / (d*d);

    fx = cos(theta) * f * GRAVITY_G;
    fy = sin(theta) * f * GRAVITY_G;

    p1->fx += fx;
    p1->fy += fy;
    p2->fx -= fx;
    p2->fy -= fy;
#endif
}

double p_collide_time(struct particle* p1, struct particle *p2) 
{
    double a, b, c, dvx, dvy, dx, dy, r, t;

    dx = p1->locx - p2->locx;
    dy = p1->locy - p2->locy;
    dvx = p1->vx - p2->vx;
    dvy = p1->vy - p2->vy;
    r = p1->r + p2->r;

    a = dvx * dvx + dvy * dvy;
    b = -2 * (dx * dvx + dy * dvy);
    c = dx * dx + dy * dy - r * r;

    t = solve(a, b, c);
    return t;
}

void p_collide(struct particle* p1, struct particle* p2)
{
    double vx1, vy1, vx2, vy2, m1, m2, x1,y1, x2, y2, d1, d2;

    vx1 = p1->vx;
    vy1 = p1->vy;

    vx2 = p2->vx;
    vy2 = p2->vy;

    m1 = p1->m;
    m2 = p2->m;

    x1 = p1->locx;
    y1 = p1->locy;

    x2 = p2->locx;
    y2 = p2->locy;

    d1 = 2 * m2 /(m1 + m2) 
        *dot_product(vx1-vx2, vy1-vy2, x1-x2, y1-y2)
        /pow(length(x1, y1, x2, y2), 2);

    p1->vx = vx1 - d1 *(x1 - x2);
    p1->vy = vy1 - d1 *(y1 - y2);

    d2 = 2 * m1 /(m1 + m2) 
        *dot_product(vx2-vx1, vy2-vy1, x2-x1, y2-y1)
        /pow(length(x1, y1, x2, y2), 2);

    p2->vx = vx2 - d2 * (x2 - x1);
    p2->vy = vy2 - d2 * (y2 - y1);

    //Let the particles move for a while to avoid immediate collision again 
    p1->locx += p1->vx*30; 
    p1->locy += p1->vy*30; 
    p2->locx += p2->vx*30; 
    p2->locy += p2->vy*30; 
}

static inline double dot_product(double a, double b, double c, double d) 
{
    return a * c + b * d;
}

static inline double length(double x, double y, double a, double b) 
{
    return sqrt((x-a)*(x-a) + (y-b)*(y-b));
}


static inline double v_post_collide(double dm, 
        double m, double m_sum, double v, double v_other)
{
    return (v * dm + 2 * m * v_other) / m_sum;
}
