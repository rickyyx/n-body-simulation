#ifndef MAIN_H
#define MAIN_H


#include "particle.h"
#include "easyppm.h"

#define MAX_NP 81

#define NUM_COLOR 2
#define MODE_BOUNCE 1
#define MODE_PASS 0

#define EVENT_COLLIDE 1
#define EVENT_BOARDER 0

#define DIR_UP 0
#define DIR_DOWN 1
#define DIR_LEFT 2
#define DIR_RIGHT 3

#define MASTER_ID num_slaves
#define CANVAS_RED 0
#define CANVAS_BLUE 1
#define DEFAULT_BUF_SIZE 100
#define BUF_SWAP_INIT_CAP 100

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

struct region {
    int timeslots;
    double timestep;
    int gridsize;

    struct particle* lgps;
    int num_large;
    int cap_large;

    struct particle* smps;
    int cap_small;
    int num_small;
    double small_mass;
    double small_radius;
    
    //Horizon
    int  horizon;
    
    //Communication
    struct particle** neigh_mov_ps;         /*  Buffer to store particles cross boarders */
    int neigh_mov_cnt[9];                   /*  Control for size of each neighbor buffer */
    int neigh_mov_cap[9];                   /*  Control for cap of the buffer */

    //Canvas to be sent to the master
    unsigned char** canvas_blue; 
    unsigned char** canvas_red;

    //Bouncing implementation
    struct event* event;
    int collision_mode;
    double config_timestep;
};

struct universe {
    //From Config File
    int gridsize;
    int timeslots;
    unsigned char* result_red;
    unsigned char* result_blue;
    unsigned char* slave_buf_red;
    unsigned char* slave_buf_blue;
    PPM  ppm;

    //Other members
    struct region* rgs;
};

struct event {
    int type;                                   /* EVENT_COLLISION or EVENT_BOARDER */
    struct particle* ps[2];                     /* Particles of the event */
    double timestamp;                           /* Future time event will happen */
    int aux;                                    /* Direction for boarder event */
};


/* **********************
 *         Master 
 * **********************/
void master();
void master_init(struct universe*);
void master_receive_region(struct universe*);
void master_draw_universe(struct universe*, int);
void master_free_universe(struct universe*);


/* **********************
 *         Slave 
 * **********************/
void slave();
void slave_init_region(struct region*);
void slave_free_region(struct region*);
void slave_simulate(struct region*);            /* Simulate for each timeslot */
void slave_convert_canvas(struct region*);      /* Convert particles to canvas */
void slave_send_region(struct region*);         /* Send canvas to the master*/
int slave_converse(int, int, struct particle*,
        const int*, struct particle**, int);    /* Converse routine between slaves */

//Update routine at each timestop
void slave_update_region(struct region*);       /* Update routine captures the below */
void slave_compute_local(struct region*);       /* Compute local forces  */
void slave_compute_neighbor(struct region*);    /* Compute forces cross regions */
struct event* slave_compute_boarder(struct region*); /* Compute collision time and event */
struct event* slave_compute_collison(struct region*); /* Compute crossing boarder event 
                                                         (in MODE_BOUNCE only) */
void slave_handle_event(struct region*);        /* Handle the event and do the updates */
void slave_sync_event(struct region*);          /* Sync events with all slaves */
void slave_prepare_move(struct region*);        /* Prepare to move the particles 
                                                   (in MODE_BOUNCE only) */
void slave_final_loc(struct region*);           /* Finalize the locations */
void slave_reset_round(struct region*);         /* Reset the round */
void slave_reset_canvas(struct region*);        /* Reset the canvas for each round */
void slave_do_move(struct region*);             /* Send the particles in the buffer to neighbors */


#endif
