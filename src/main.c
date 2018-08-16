#define ACCU_TEST

#include "mpi.h"
#include "gmp.h"
#include "main.h"
#include "util.h"
#include "easyppm.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*  Static variables functions definition */
struct universe uni;
static struct region region;
static int myid, myrank, num_slaves;
static int nprocs;
static int proc_side;
static struct m_time *s_comm_timer;
static MPI_Comm slave_comm;

/*   Local helpers */
static void init_canvas(struct region*);
static int get_slave_id(int, int);
static void ps_final_loc(struct particle*,int, struct region*);
static void p_cross_region(struct region*, struct particle*);
static void add_swap(struct particle*, int, struct region*);
static int resize_ps(struct particle** list, int* cap, int size);
static bool init_sps(struct region*);
static void draw_circle_canvas(unsigned char**, struct particle*, int, int);
#ifdef ACCU_TEST
static void init_precise_ps(struct region*);
static void free_mpf(struct particle*, int);
#endif





#define BUF_FOR_SLAVE(res, id, size, row) ((res)+ (id)*(size)*(size)+(row)*(size))
#define DISTANCE(v, a, t) ((v)*(t) + 0.5*(a)*(t)*(t))
#define OFFSET_ROW_COL(x, delta) (((x) + (delta) + proc_side) % proc_side)
#define SLAVE_ID(r, c) (r * proc_side + c)
#define NEIGHBOR_IDX(x_off, y_off, l) ((x_off) * l + (y_off) + (l * l)/2)



int main(int argc,char ** argv) 
{
    struct m_time *total;

    srand(1);
    int ndims, reorder, dim_size[2], periods[2], coords[2];

    //Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    num_slaves = nprocs - 1;
    proc_side = sqrt(num_slaves);

    //Create new communicator for slaves
    slave_comm = MPI_COMM_WORLD;
    ndims = 2;
    dim_size[0] = dim_size[1] = proc_side;
    periods[0] = periods[1] = 1;
    reorder = 1;
    
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dim_size, periods, reorder, &slave_comm);
    coords[0] = myrank / proc_side;
    coords[1] = myrank % proc_side;

    if(myrank != MASTER_ID) {
      MPI_Cart_rank(slave_comm, coords, &myid);
    } else {
        myid = myrank;
    }

    if(myid == MASTER_ID) {
        total = init_timer("Total Time", myid);
        start_timer(total);
        master();
    } else {
        total = init_timer("Slave Time", myid);
        start_timer(total);
        slave();
    }

    MPI_Finalize();
    pause_timer(total);

    if(s_comm_timer != NULL) 
        print_timer(s_comm_timer);

    return 0;
}


void master()
{
    int i = 0;
    master_read_config(&uni);
    master_init(&uni);

    master_receive_region(&uni);
    master_draw_universe(&uni, i);

    master_free_universe(&uni);
}

/* 
 * Initialize resources for the master
 * */
void master_init(struct universe* uni) 
{
    //Allocate buffers
    int grid_size, img_size;
    grid_size = uni->gridsize;

    allocate_buf_1d(&uni->result_blue, grid_size*grid_size*(nprocs-1));
    allocate_buf_1d(&uni->result_red, grid_size*grid_size*(nprocs-1));

    //Utility buf used for receiving
    allocate_buf_1d(&uni->slave_buf_red, grid_size);
    allocate_buf_1d(&uni->slave_buf_blue, grid_size);

    //Allocate Image ppm
    img_size = (int)sqrt(nprocs-1) * grid_size;
    uni->ppm = easyppm_create(img_size, img_size, IMAGETYPE_PPM);
}

/*  
 *  Master integrates the canvas
 *  */
void master_draw_universe(struct universe* uni, int time)
{
    unsigned char* res_blue = uni->result_blue;
    unsigned char* res_red = uni->result_red;
    int grid_size, img_size, proc_size, proc_offset, x, y, i, j, r, c, idx, proc_num;
    PPM* ppm;
    char img_name[20];

    grid_size = uni->gridsize;
    proc_size = sqrt(nprocs-1);
    proc_offset = grid_size * grid_size;
    img_size =  grid_size * proc_size;

    ppm = &(uni->ppm);
    easyppm_clear(ppm, easyppm_rgb(0,0,0));

    sprintf(img_name, "img/image-%d.ppm", time);

    //Draw Small particles Red and Large blue
    for(x = 0; x<img_size; x++) {
        for(y =0; y < img_size; y++) {
            i = x / grid_size;
            j = y / grid_size;

            proc_num = j * proc_size +  i;

            r = x % grid_size;
            c = y % grid_size;

            idx = proc_num * proc_offset + r * grid_size + c;
            easyppm_set(ppm, x, y, easyppm_rgb(res_red[idx], 0, res_blue[idx]));
        }
    }

    easyppm_write(ppm, img_name);
}

/*  Master waits for canvas from all slaves */
void master_receive_region(struct universe* uni) 
{
    MPI_Status status;
    unsigned char* result_red, *slave_buf_red, *result_blue, *slave_buf_blue;
    int size, slave_id, row, offset;
    
    //Buffer pointers
    result_red = uni->result_red;
    slave_buf_red = uni->slave_buf_red;
    result_blue = uni->result_blue;
    slave_buf_blue = uni->slave_buf_blue;

    size = uni->gridsize;

    for(slave_id = 0; slave_id < num_slaves; slave_id++) {
        for(row = 0; row < size; row++) {
            offset = slave_id * size *size + row *size;
            //Receive red canvas from slave
            MPI_Recv(slave_buf_red, size, MPI_UNSIGNED_CHAR, slave_id,CANVAS_RED*size+row,MPI_COMM_WORLD, &status);
            //Copy to master's red canvas
            memcpy(result_red+offset, slave_buf_red, size);

            //Receive blue canvas from slave
            MPI_Recv(slave_buf_blue, size, MPI_UNSIGNED_CHAR, slave_id,CANVAS_BLUE*size+row,MPI_COMM_WORLD, &status);

            //Copy to master's blue canvas
            memcpy(result_blue+offset, slave_buf_blue, size);
        }
    }
}

/*  Free master's resources */
void master_free_universe(struct universe* uni)
{
    free(uni->result_red);
    free(uni->result_blue);
    free(uni->slave_buf_red);
    free(uni->slave_buf_blue);

    easyppm_destroy(&uni->ppm);
}



/* ***************************
*           Slave 
 * ****************************/
void slave()
{
    s_comm_timer = init_timer("S_COMM", myid);
    slave_read_config(&region);
    slave_init_region(&region);
    slave_simulate(&region);
    slave_free_region(&region);
}


/*
 * Slave convert particles to local canvas 
 * */
void slave_convert_canvas(struct region* rg)
{
    unsigned char **canvas_blue, **canvas_red;
    int i, x, y;
    struct particle  *p;

    canvas_blue = rg->canvas_blue;
    canvas_red = rg->canvas_red;

    //Conver particles to canvas
    for(i = 0; i<rg->num_large;i++) {
        p = rg->lgps+i;
        draw_circle_canvas(canvas_blue, p, rg->gridsize, CANVAS_BLUE);
    }
    for(i = 0; i<rg->num_small;i++) {
        p = rg->smps+i;
        x = round(p->locx);
        y = round(p->locy);
        if(canvas_red[x][y] < 255) {
            canvas_red[x][y]++;
        }
    }
}

/*  
 *  Slave send the local canvas to masters
 *  */
void slave_send_region(struct region* rg)
{
    int size, i;
    size = rg->gridsize;
    //Send two particles list over 
    start_timer(s_comm_timer);
    for(i = 0; i<size; i++) {
        MPI_Send(rg->canvas_blue[i], size, MPI_UNSIGNED_CHAR, MASTER_ID, CANVAS_BLUE*size+i, MPI_COMM_WORLD);
        MPI_Send(rg->canvas_red[i], size, MPI_UNSIGNED_CHAR, MASTER_ID, CANVAS_RED*size+i,MPI_COMM_WORLD);
    }
    pause_timer(s_comm_timer);
}

/*  Slav free resources */
void slave_free_region(struct region* rg)
{
    int i;

#ifdef ACCU_TEST
    free_mpf(rg->smps, rg->num_small);
    free_mpf(rg->lgps, rg->num_large);
#endif

    free(rg->smps);
    free(rg->lgps);

    for(i =0; i<8; i++) {
        free(rg->neigh_mov_ps[i]);
    }
    free(rg->neigh_mov_ps);

    free_buf_2d(rg->canvas_red, rg->gridsize);
    free_buf_2d(rg->canvas_blue, rg->gridsize);
}


/*  
 *  Slave running the simulation at each timeslot
 *  */
void slave_simulate(struct region* rg)
{
    int i;

    for(i = 0; i< rg->timeslots; i++) {
        slave_update_region(rg);
    }
    
    slave_convert_canvas(&region);
    slave_send_region(&region); 
    slave_reset_canvas(&region);
}

/*  Slave resetting canvas after sending 
 *  (Only for multiple drawings in a run) */
void slave_reset_canvas(struct region* rg) 
{
    unsigned char ** canvas_blue, ** canvas_red;
    int i, grid_size;

    canvas_blue = rg->canvas_blue;
    canvas_red = rg->canvas_red;
    grid_size = rg->gridsize;

    for(i = 0; i<grid_size; i++) {
        memset(canvas_red[i], 0, sizeof(unsigned char) * grid_size);
        memset(canvas_blue[i], 0, sizeof(unsigned char) * grid_size);
    }
}


void slave_update_region(struct region* rg)
{
    struct event *ev_collide, *ev_boarder;
    // Compute forces from local
    slave_compute_local(rg);

    // Send to neighbors and receive from neighbors to compute cross regions
    // forces
    slave_compute_neighbor(rg);
    
    // Bouncing off
    if(rg->collision_mode == MODE_BOUNCE) {
        //Compute the local minimal timestep  
        ev_collide = slave_compute_collison(rg);
        ev_boarder = slave_compute_boarder(rg);
        

        if(ev_collide != NULL || ev_boarder != NULL) {
            if(ev_collide == NULL) { /*  Boarder Event Happens */
                rg->event = ev_boarder;
            } else if(ev_boarder == NULL) { /*  Collide Event Happens */
                rg->event = ev_collide;
            } else if(ev_collide->timestamp < ev_boarder->timestamp){
                /*  Collide Happens First */
                rg->event = ev_collide;
            } else {
                /*  Boarder First  */
                rg->event = ev_boarder;
            }
        }
    }
    
    //Synchronize events with all slaves
    slave_sync_event(rg);
    
    // Finalize locations
    slave_final_loc(rg);

    // Handle event
    slave_handle_event(rg);

    //Boundary adjust
    slave_prepare_move(rg);

    // Move particles for next round 
    if(rg->event == NULL || rg->event->type != EVENT_BOARDER)
        slave_do_move(rg);

    slave_reset_round(rg);
}

/*  Preparing to move the particles by adding the particles to 
 *  the corresponding buffer */
void slave_prepare_move(struct region* rg)
{
    int i;
    struct particle* p;

    for(i = 0; i<rg->num_small+rg->num_large; i++) {
        p = (i < rg->num_small) ? rg->smps+i : rg->lgps+(i-rg->num_small);
        p_cross_region(rg, p);
    }
}

/*  Slave Synchronize on the timestep with all other */
void slave_sync_event(struct region* rg)
{
    int i;
    struct event * ev;
    double time_to_send, * times_recv, time_min;
    ev = rg->event;
    if(ev == NULL) { /*  Default timestep */
        time_to_send = rg->timestep;
    } else {
        time_to_send  = ev->timestamp;
    }

    times_recv = calloc(num_slaves, sizeof(double));
    if(times_recv == NULL) {
        fprintf(stderr, "Out of memroy: times_send/recv\n");
        exit(1);
    }

    //All-to-all broadcast
    MPI_Allgather(&time_to_send, 1, MPI_DOUBLE, times_recv,1 , MPI_DOUBLE, slave_comm);
    
    //Get the global minimal safe timestep
    time_min = rg->timestep;
    for(i = 0; i<num_slaves; i++) {
        time_min = MIN(time_min, times_recv[i]);
    }

    free(times_recv);
    if(time_min != time_to_send) { 
        /*  Someone other processes have events ealier */
        free(rg->event);
        rg->event = NULL;
    }
}

/*  Handle events:
 *  @Collision: calculate post-collision velocities
 *  @Crossing boarder: move to the correct buffer queue
 *  */
void slave_handle_event(struct region* rg)
{
    struct event * ev;
    struct particle **ps, *p;
    int grid_size, dir;

    ev  = rg->event;
    if(ev == NULL) return;
    grid_size = rg->gridsize;
    

    if(ev->type == EVENT_COLLIDE) { /*  Collision */
        rg->timestep = ev->timestamp;
        ps = ev->ps;
       
        //Particles collide and recalculate velocity
        p_collide(ps[0], ps[1]);
        
        free(rg->event);
        rg->event = NULL;
    } else if(ev->type  == EVENT_BOARDER) {
        rg->timestep = ev->timestamp;
        
        p = ev->ps[0];
        dir = ev->aux; 
       
        //Adjust the locations to out of bound so that 
        //move() will put the particle into the correct buffer for neighbor
        switch(dir){
            case DIR_UP:
                p->locy = -1;
            case DIR_DOWN:
                p->locy = grid_size;
            case DIR_LEFT:
                p->locx = -1;
            case DIR_RIGHT:
                p->locx = grid_size;
        }

        free(rg->event);
        rg->event = NULL;

    } else {
        fprintf(stderr, "Unsupported event %d \n", ev->type);
        exit(1);
    }
}


/*  Slave calculate the nearest collision timestep 
 *  by going through each particle pair, and solve 
 *  for the minimal time where two particles collide, 
 *  while updating the local minimal at the same time
 *  */
struct event* slave_compute_collison(struct region* rg) 
{
    struct event * ev;
    struct particle *p1, *p2;
    int num_small, num_large, i, j;
    double t_temp, t_collision;

    ev = malloc(sizeof(struct event));
    if(ev == NULL) {
        fprintf(stderr, "Out of memory: comp_collison\n");
    }
    ev->type = EVENT_COLLIDE;

    num_small = rg->num_small;
    num_large = rg->num_large;

    t_collision = rg->timestep;

    for(i = 0; i< num_small+num_large; i++) {
        p1 = (i < num_small) ? rg->smps+i : rg->lgps+(i-num_small);
        
        /*  For each particle pair */
        for(j = i+1; j< num_small+num_large; j++) {
            p2 = (j < num_small) ? rg->smps+j : rg->lgps+(j-num_small);

            // Calculate the collision time
            t_temp = p_collide_time(p1, p2);
            
            // Set the local minimal collision time
            if(t_temp > 0 && t_temp < t_collision) { 
                t_collision = t_temp;
                ev->ps[0] = p1;
                ev->ps[1] = p2;
                ev->timestamp = t_collision;
            }
        }
    }

    if(t_collision == rg->timestep) { /* No Collision detected */
        free(ev);
        return NULL;
    }

    return ev;
}

/*  Compute the nearest time for a crossing boarder event to occur
 *  Similar approach with slave_compute_collision()
 *  */
struct event* slave_compute_boarder(struct region* rg) 
{
    double x,y, t, t_cross, tx0, tx1, ty0, ty1, ax, ay;
    struct particle* p;
    struct event * ev;
    int num_small, num_large, i, grid_size;

    ev = malloc(sizeof(struct event));
    if(ev == NULL){
        fprintf(stderr, "Out of memory: comp_neighbor\n");
        exit(1);
    }
    ev->type = EVENT_BOARDER;
    ev->ps[0]= ev->ps[1] = NULL;

    num_small = rg->num_small;
    num_large = rg->num_large;
    grid_size = rg->gridsize;
    t = rg->timestep;
    t_cross = t;

    for(i = 0; i< num_small + num_large; i++) {
        p = ( i< num_small) ? rg->smps+i : rg->lgps+(i-num_small);

        x = p->locx;
        y = p->locy;

        ax = p->fx / p->m;
        ay = p->fy / p->m;

        x+=DISTANCE(p->vx, ax, t); 
        y+=DISTANCE(p->vy, ay, t);

        // Particles will be out of bound
        if(x < 0 || x > grid_size -1 || y < 0 || y > grid_size - 1) {
            //X direction - time it crosses the LEFT
            tx0 = solve(ax/2, p->vx, x);
            if(tx0 > 0 && tx0 < t_cross) {
                t_cross = tx0;
                ev->timestamp = t_cross;
                ev->ps[0] = p;
                ev->aux = DIR_LEFT;
            }

            //X direction - time it crosses the RIGHT
            tx1 = solve(ax/2, p->vx, x - grid_size +1);
            if(tx1 > 0 && tx1 < t_cross) {
                t_cross = tx1;
                ev->timestamp = t_cross;
                ev->ps[0] = p;
                ev->aux = DIR_RIGHT;
            }
            
            //Y direction - time it crosses TOP
            ty0 = solve(ay/2, p->vy, y);
            if(ty0 > 0 && ty0 < t_cross) {
                t_cross = ty0;
                ev->timestamp = t_cross;
                ev->ps[0] = p;
                ev->aux = DIR_UP;
            }
            
            //Y direction - time it crosses BOTTOM
            ty1 = solve(ay/2, p->vy, y - grid_size +1);
            if(ty1 > 0 && ty1 < t_cross) {
                t_cross = ty1;
                ev->timestamp = t_cross;
                ev->ps[0] = p;
                ev->aux = DIR_DOWN;
            }
        }
    }
    
    if(t_cross < rg->timestep) {
        return ev;
    } else {
        free(ev);
        return NULL;
    }
}


/*  Resetting the counting variables for each timeslot  */
void slave_reset_round(struct region* rg) 
{
    //Reset comm buf
    int i;

    for(i = 0; i<9; i++) {
        rg->neigh_mov_cnt[i] = 0;
    }

    //Reset Event and timesteps
    rg->timestep = rg->config_timestep;
}

/* 
 * Communication routine for slaves
 * */
int slave_converse(int dx, int dy, struct particle* data, 
        const int* send_cnt, struct particle** buf, int buf_size)
{ 
    MPI_Request req_cnt, req_ps;
    MPI_Status status_cnt, status_ps;
    int send_to, recv_from, recv_cnt;

    start_timer(s_comm_timer);
    
    //Calculate neighbors to converse with in this round
    send_to = get_slave_id(dx, dy);
    recv_from = get_slave_id(0-dx, 0-dy); 
    
    
    //Exchange the control information: how much data to be send
    MPI_Isend(send_cnt, 1, MPI_INT, send_to, 0, slave_comm, &req_cnt);
    MPI_Recv(&recv_cnt, 1, MPI_INT, recv_from, 0, slave_comm, &status_cnt);

    MPI_Wait(&req_cnt, &status_cnt);
    
    //Adjust the buffer if it is too small to have all data
    if(buf_size < recv_cnt) { /* Buf too small */
        *buf = realloc(*buf, recv_cnt * sizeof(struct particle));
        if(*buf == NULL) {
            fprintf(stderr, "OoM: resize_buf_comm\n");
            exit(1);
        }
    }
    
    //Exchange data:  particle arrays
    MPI_Isend(data, *send_cnt * sizeof(struct particle), MPI_BYTE, send_to, 0, slave_comm, &req_ps);
    MPI_Recv(*buf, recv_cnt * sizeof(struct particle), MPI_BYTE, recv_from, 0, slave_comm, &status_ps);
    MPI_Wait(&req_ps, &status_ps);

    pause_timer(s_comm_timer);
    return recv_cnt;
}


/*  
 *  Moving of the particles to the neighbors
 *  */
void slave_do_move(struct region* rg) 
{
    int dx, dy, recv_idx, cnt, i;

    struct particle* buf;

    //Buffer for storing particles from neighbors
    buf = malloc(sizeof(struct particle) * DEFAULT_BUF_SIZE);
    if(buf == NULL) {
        fprintf(stderr, "Out of memory - buf_move\n");
        exit(1);
    }

    for(dx = -1; dx <= 1; dx++) {
        for(dy = -1; dy <= 1; dy++) {
            if(dx == 0 && dy  == 0) continue; /*  Skip itself */

            recv_idx = NEIGHBOR_IDX(dx, dy, 3);

            // Call slave converse routine
            cnt = slave_converse(dx, dy, rg->neigh_mov_ps[recv_idx],
                    &rg->neigh_mov_cnt[recv_idx], &buf, DEFAULT_BUF_SIZE);

            //Combines particles from neighbors to itself
            for(i = 0; i< cnt; i++) {
                if(buf[i].is_large) {
                    rg->lgps[rg->num_large++] = buf[i];
                    resize_ps(&rg->lgps, &rg->cap_large, rg->num_large);
                } else {
                    rg->smps[rg->num_small++] = buf[i];
                    resize_ps(&rg->smps, &rg->cap_small, rg->num_small);
                }
            }

        }
    }
    free(buf);
}


/*  Finalize the locations of each particle
 *  */
void slave_final_loc(struct region* rg) 
{
    ps_final_loc(rg->smps, rg->num_small, rg);
    ps_final_loc(rg->lgps, rg->num_large, rg);
}

/*  Compute the forces acting from neighbors:
 *  It first receives and sends the particles.
 *  And calculate each pair's force.
 *  */
void slave_compute_neighbor(struct region* rg)
{

    int dx, dy, cnt, i,j, horizon, buf_size, buf_size_s, cnt_s;
    struct particle* buf, *p1, *p2, *buf_s;
    double x_off, y_off;
    
    /* For large and small paricles */
    buf_size = 0;
    buf_size_s = rg->num_small; 
    buf = NULL;
    buf_s = malloc(sizeof(struct particle) * buf_size_s);
    if(buf_s == NULL) {
        fprintf(stderr, "Out of memory: buf_s\n");
    }

    horizon = rg->horizon;

    for(dx = 0 - horizon; dx <= horizon; dx++) { 
        for(dy = 0 - horizon; dy <= horizon; dy++) {
            if(dx == 0 && dy == 0) continue;
            // Exchange large particles 
            cnt = slave_converse(dx, dy, rg->lgps, &rg->num_large, &buf, buf_size);
            buf_size = MAX(buf_size, cnt);

            // Exchange small particles
            cnt_s = slave_converse(dx, dy, rg->smps, &rg->num_small, &buf_s, buf_size_s);
            buf_size_s = MAX(buf_size_s, cnt_s);

            //Consolidate impacts from other processes
            /*  Relative position of processes  */
            x_off = (0-dx) * rg->gridsize; 
            y_off = (0-dy) * rg->gridsize;
    
            //Neighbors' Large parciles
            for(i = 0; i < cnt; i++) {
                p1 = &buf[i];
                for(j = 0; j<rg->num_small; j++) {
                    p2 = rg->smps+j;
                    update_force(p1, p2, myid, x_off, y_off);
                }
                for(j = 0; j<rg->num_large; j++) {
                    p2 = rg->lgps+j;
                    update_force(p1, p2, myid, x_off, y_off);
                }
            }

            //Neighnbor's Small particles
            for(i = 0; i<cnt_s; i++) {
                p1 = &buf_s[i];

                for(j = 0; j<rg->num_small; j++) {
                    p2 = rg->smps+j;
                    update_force(p1, p2, myid, x_off, y_off);
                }

                for(j = 0; j<rg->num_large; j++) {
                    p2 = rg->lgps+j;
                    update_force(p1, p2, myid, x_off, y_off);
                }
            }
        }
    }

    if(buf != NULL) 
        free(buf);
}

/*  
 *  Computing the forces of particles within a region
 *  */
void slave_compute_local(struct region* rg)
{
    struct particle *p1, *p2;
    int i,j;


    //update small with the rest
    for(i = 0; i<rg->num_small; i++) {
        p1 = rg->smps+i;
        for(j = i+1; j<rg->num_small; j++) {
            p2 = rg->smps+j;
            update_force(p1, p2, myid, 0, 0);
        }

        for(j =0; j<rg->num_large; j++) {
            p2 = rg->lgps+j; 
            update_force(p1, p2, myid, 0, 0);
        }
    }

    //update large among themselves
    for(i = 0; i<rg->num_large; i++) {

        p1 = rg->lgps+i;
        for(j = i+1; j<rg->num_large; j++) {

            p2 = rg->lgps+j;
            update_force(p1, p2, myid, 0, 0);
        }
    }

}


void slave_init_region(struct region* rg) 
{
    bool ok = false;
    int i;

    //Small particles (Large particles init in read_config
    ok = init_sps(rg);
    if(!ok) {
        printf("init_sp failed\n");
        exit(1);
    }

#ifdef ACCU_TEST
    init_precise_ps(rg);
#endif

    //Communication Channels
    rg->neigh_mov_ps = malloc(sizeof(struct particle*) * 9);
    if(rg->neigh_mov_ps == NULL) {
        fprintf(stderr, "OoM: neigh_mov_ps\n");
        exit(1);
    }
    for(i = 0; i<9; i++) {
        rg->neigh_mov_cnt[i] =0;
        rg->neigh_mov_cap[i] = BUF_SWAP_INIT_CAP;
        rg->neigh_mov_ps[i] = malloc(sizeof(struct particle) * rg->neigh_mov_cap[i]);
        if(rg->neigh_mov_ps[i] == NULL) {
            fprintf(stderr, "OoM: neigh_mav_ps[i]\n");
        }
    }

    //Canvas
    init_canvas(rg);
}

/*  Initizie the local canvas */
static void init_canvas(struct region* rg) 
{
    int grid_size, i;
    grid_size = rg->gridsize;

    rg->canvas_blue = (unsigned char**) malloc(sizeof(unsigned char*) * grid_size);
    rg->canvas_red = (unsigned char**) malloc(sizeof(unsigned char*) * grid_size);

    if(rg->canvas_blue == NULL || rg->canvas_red == NULL) {
        fprintf(stderr, "Out of memory: canvas \n");
        exit(1);
    }

    for(i=0; i<grid_size; i++) {
        rg->canvas_blue[i] =  (unsigned char*) calloc((size_t) grid_size, sizeof(unsigned char));
        rg->canvas_red[i] =  (unsigned char*) calloc((size_t) grid_size, sizeof(unsigned char));

        if(rg->canvas_blue[i] == NULL || rg->canvas_red[i] == NULL) {
            fprintf(stderr, "Out of memory: canvas \n");
            exit(1);
        }
    }
}

/*  Initizate the small particles */
static bool init_sps(struct region * rg)
{
    int i;
    rg->cap_small = 2 * rg->num_small;
    rg->smps = malloc(sizeof(struct particle) * rg->cap_small);
    if(rg->smps == NULL) {
        return false;
    }

    for(i = 0; i < rg->num_small; i++) {
        make_sp((rg->smps+i),rg->gridsize, rg->small_mass, rg->small_radius);
    }
    return true;
}

/*  Helper function to draw a large particle to canvas */
static void draw_circle_canvas(unsigned char** canvas, struct particle* p, int size, int color)
{
    int xs, xe, ys, ye, x, y;
    double r = p->r;
    xs = MAX(0, p->locx - r);
    xe = MIN(p->locx + r, size);
    ys = MAX(0, p->locy - r);
    ye = MIN(p->locy + r, size);


    for(x = xs; x < xe; x++) {
        for(y = ys; y<ye; y++) {
            if(distance(x,y, p->locx, p->locy) <= r) {
                canvas[x][y] = 255;
            }
        }
    }
}

/*  Finalize locations of each particle in the ps array */
static void ps_final_loc(struct particle* ps, int size, struct region* rg)
{
    int i;
    double t, ax, ay;
    struct particle *p;

    t= rg->timestep;

    for(i = 0; i<size; i++) {
        p=ps+i;
        //get acceleration
        ax = p->fx / p->m;
        ay = p->fy / p->m;

        //calculate locations
        p->locx+=DISTANCE(p->vx, ax, t); 
        p->locy+=DISTANCE(p->vy, ay, t);

        //calculate new velocity
        p->vx += ax * t;
        p->vy += ay * t;

        //reset the force 
        p->fx = p->fy = 0;
    }
}


/*  Move a particle to the buffer for neighbors,
 *  and delete it from the local particles array*/
static void p_cross_region(struct region* rg, struct particle* p)
{
    double x, y;
    int grid_size;

    x = p->locx;
    y = p->locy;
    grid_size = rg->gridsize;


#ifdef ACCU_TEST
    /* Simplification for accuracy testing (no crossing of boaders)*/
    if(p->locx < 0 || x > grid_size-1 || y < 0 || y > grid_size-1) {
        if(x < 0) 
            p->locx = grid_size-1;

        if(x > grid_size-1) 
            p->locx = 0;

        if(y < 0)
            p->locy = grid_size-1;

        if(y > grid_size-1)
            p->locy = 0;
    }
#else
    if(x < 0 || x > grid_size-1 || y < 0 || y > grid_size-1) {
        if(x < 0 && y < 0) { /* Top left */
            p->locx = grid_size-1;
            p->locy = grid_size-1;
            add_swap(p, NEIGHBOR_IDX(-1, -1, 3), rg);
        }
        else if(x < 0 && y > grid_size-1) { /*  Bottom Left */
            p->locx = grid_size-1;
            p->locy = 0;
            add_swap(p, NEIGHBOR_IDX(1, -1, 3), rg);
        }
        else if(x > grid_size-1 && y < 0) { /*  Top Right */
            p->locx = 0;
            p->locy = grid_size-1;
            add_swap(p, NEIGHBOR_IDX(-1, 1, 3), rg);
        } 
        else if(x > grid_size-1 && y > grid_size-1) { /* Bottom Right */
            p->locx = 0;
            p->locy = 0;
            add_swap(p, NEIGHBOR_IDX(1, 1, 3), rg);
        }
        else if(x < 0) { /*  Left */
            p->locx = grid_size-1;
            add_swap(p, NEIGHBOR_IDX(0,-1, 3), rg);
        } 
        else if(x > grid_size-1) { /*  Right */
            p->locx = 0;
            add_swap(p, NEIGHBOR_IDX(0, 1, 3), rg);
        }
        else if(y < 0) { /* Top */
            p->locy = grid_size-1;
            add_swap(p, NEIGHBOR_IDX(-1,0, 3), rg);
        } 
        else if(y> grid_size-1) { /*  Down  */
            p->locy = 0;
            add_swap(p, NEIGHBOR_IDX(1, 0, 3), rg);
        }
        else {
            fprintf(stderr, "Error....\n");
            exit(1);
        }
    }
#endif
}

/*  Add a particle to the neighbor buffer to be sent
 *  and swap the last particle in the particle array
 *  to replace the hole
 *  */
static void add_swap(struct particle* p, int dir, struct region* rg)
{
    int idx;
    int resized;
    idx = rg->neigh_mov_cnt[dir]++;

    //Copy to buffer (both small and large)
    rg->neigh_mov_ps[dir][idx] = *p;

    //Swap with last guy in the particle array
    if(p->is_large) {
        *p = rg->lgps[--rg->num_large];
        if(rg->num_large == 0) 
            memset(rg->lgps, 0, sizeof(struct particle));
    } else {
        *p = rg->smps[--rg->num_small];
        if(rg->num_small == 0) 
            memset(rg->smps, 0, sizeof(struct particle));
    }
    
    //Resize if needed
    resized = resize_ps(&rg->neigh_mov_ps[dir], &rg->neigh_mov_cap[dir], rg->neigh_mov_cnt[dir]);
    if(resized) {
        printf("Reszie: %d , %d\n", rg->neigh_mov_cap[dir], rg->neigh_mov_cnt[dir]);
    }

}

/*  Utility function to get the neighbor's slave id */
static int get_slave_id(int dx, int dy) 
{
    int r, c;

    r = myid / proc_side;
    c = myid % proc_side;

    return SLAVE_ID(OFFSET_ROW_COL(r, dx), OFFSET_ROW_COL(c, dy)); 
}

static int resize_ps(struct particle** list, int* cap, int size) 
{
    if(size >= *cap) {
        *cap = *cap * 2;
        printf("Reallocating on resize: %d -> %d\n", size, *cap);
        *list = realloc(*list, (*cap) * sizeof(struct particle));
        if(*list == NULL) {
            fprintf(stderr, "Out of memory: resizing\n");
            exit(1);
        }

        return true;
    } 

    return false;
}

#ifdef ACCU_TEST
static void init_precise_ps(struct region* rg)
{
    int i;
    struct particle *p;

    init_G();

    for(i = 0; i<rg->num_large; i++) {
        p = &rg->lgps[i];

        mpf_init_set_d(p->m_m, p->m);
        mpf_init_set_d(p->m_r, p->r);
        mpf_init_set_d(p->m_locx, p->locx);
        mpf_init_set_d(p->m_locy, p->locy);
        mpf_init_set_d(p->m_vx, p->vx);
        mpf_init_set_d(p->m_vy, p->vy);
        mpf_init_set_d(p->m_fx, p->fx);
        mpf_init_set_d(p->m_fy, p->fy);

    }

    for(i = 0; i<rg->num_small; i++) {
        p = &rg->smps[i];

        mpf_init_set_d(p->m_m, p->m);
        mpf_init_set_d(p->m_r, p->r);
        mpf_init_set_d(p->m_locx, p->locx);
        mpf_init_set_d(p->m_locy, p->locy);
        mpf_init_set_d(p->m_vx, p->vx);
        mpf_init_set_d(p->m_vy, p->vy);
        mpf_init_set_d(p->m_fx, p->fx);
        mpf_init_set_d(p->m_fy, p->fy);
    }
}

static void free_mpf(struct particle* ps, int size)
{
    int i;
    struct particle* p;

    for(i = 0; i<size; i++) {
        p = ps+i;

        mpf_clears(p->m_r,
                p->m_m,
                p->m_locx,
                p->m_locy,
                p->m_vx,
                p->m_vy,
                p->m_fx,
                p->m_fy,
                NULL);
    }
}
#endif

