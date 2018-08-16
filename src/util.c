#include "util.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define MAX_BUF_SIZE 256
static void parse(char*, char*, struct region *, FILE*);
static void parse_large_p(struct particle*, int, FILE*);

const char* config_file_name = "initialspec.txt";

/* ***********************
 * Utility Functions
 * ************************/

void allocate_buf_1d(unsigned char** ptr, size_t size) 
{
    *ptr = (unsigned char*)malloc(sizeof(unsigned char)*size);
    if(*ptr == NULL) {
        fprintf(stderr, "Out of memory: 1d\n");
        exit(1);
    }
}


void allocate_buf_2d(unsigned char*** ptr, size_t x_size, size_t y_size) 
{
    int i;
    *ptr = (unsigned char**) malloc(sizeof(unsigned char*) * x_size);
    if(*ptr == NULL) {
        fprintf(stderr, "Out of memory: 2d\n");
        exit(1);
    }

    for(i = 0; i<x_size; i++) {
        (*ptr)[i] = (unsigned char*) calloc(sizeof(unsigned char), y_size);
        if((*ptr)[i] == NULL) {
            fprintf(stderr, "Out of memory: 2d\n");
            exit(1);
        }
    }
}

void free_buf_2d(unsigned char** ptr, size_t size) 
{
    int i;

    for(i = 0; i < size; i++) {
        free(ptr[i]);
    }

    free(ptr);
}

void debug_region(struct region* rg) 
{
    int i = 0, j, k;
    printf("---------------REGION-----------------\n");
    printf(" P: %d\n", rg->num_small); 
    printf("Large P: %d\n", rg->num_large);

    for(i =0; i< MIN(rg->num_small, 5); i++) {
        p_print(rg->smps[i]);
    }

    for(i = 0; i<MIN(rg->num_large, 3); i++) {
        p_print(rg->lgps[i]);
    }

    //p_canvas(rg->canvas_blue, rg->gridsize);
    //p_canvas(rg->canvas_red, rg->gridsize);

    //Movement
    printf("------------Movement Counters--------\n");
    for(i =0, k = 0; i<3; i++) {
        for(j =0; j<3; j++) {
            if(i == 1 && j == 1){
                printf("(x,x),");
                continue;
            }
            printf("(%d %d),", rg->neigh_mov_cnt[k], rg->neigh_mov_cap[k]);
            k++;
        }
        printf("\n");
    }

    printf("--------------END REGION--------------\n");
}

void p_canvas(unsigned char** canvas, int size) 
{
    int i, j;
    unsigned char * data;
    printf("--------- Canvas ----------\n");
    for(i = 0; i<size; i++) {
        data = canvas[i];
        printf(" [");
        for(j = 0; j< size; j++) {
            printf("%u,", data[j]); 
        }
        printf(" ]\n");
    }
}

void debug_buf_1d(unsigned char* arr, int size) 
{
    int i;

    printf("[");
    for(i = 0; i<size; i++) {
        printf("%u, ", arr[i]);
    }
    printf("]\n");
}

void p_result(unsigned char* res, int size, int nproc) 
{
    int i, j, k;
    printf("---------------------- BUF -------------------\n");

    for(i = 0; i<nproc; i++) {
        printf("----- PROC: %d ---------\n\n", i);
        for(j = 0; j<size; j++) {
            for(k = 0; k<size; k++) {
                printf("%u, ", res[i*size*size+j*size+k]);
            }
            printf("\n");
        }
    }
}


void debug_event(struct event* ev) 
{
    if(ev != NULL) {
        printf("\n---------------EVENT--------------\n");
        printf("Event Type: %d\n", ev->type);
        printf("TimeStamp: %lf\n", ev->timestamp);
        
        printf("[Particles] \n");
        if(ev->ps[0] != NULL){
            p_print(*ev->ps[0]);
        }
        if(ev->ps[1] != NULL) {
            p_print(*ev->ps[1]);
        }
        printf("Event Aux:%d\n", ev->aux);
    }
}


void 
print_config(struct region * cfg) 
{
    int i;
    printf("Timeslots: %d\n", cfg->timeslots);
    printf("Timestep: %f\n", cfg->timestep);
    printf("Horizon: %d\n", cfg->horizon);
    printf("GridSize: %d\n", cfg->gridsize);
    printf("num_s_p: %d\n", cfg->num_small);
    printf("small_mass: %f\n", cfg->small_mass);
    printf("small_radius: %f\n", cfg->small_radius);
    printf("num_large_p: %d\n", cfg->num_large);

    for(i = 0; i<cfg->num_large; i++) {
        p_print(cfg->lgps[i]);
    }
}


void master_read_config(struct universe* cfg)
{
    char buf[MAX_BUF_SIZE], key[MAX_BUF_SIZE], value[MAX_BUF_SIZE];
    FILE *file;

    file = fopen(config_file_name, "r");

    while((fgets(buf, MAX_BUF_SIZE, file)) != NULL) {
        if(strlen(buf) <= 0) 
            continue;
        
        //Assume no comments, only key:value pair each line
        sscanf(buf, "%[^:]: %s", key, value);
        
        //Only need grid size
        if(strcmp(key, "GridSize") == 0) 
            cfg->gridsize = atoi(value);
        else if(strcmp(key, "TimeSlots") == 0) 
            cfg->timeslots = atoi(value);
    }
}

void
slave_read_config(struct region * cfg) 
{
    char buf[MAX_BUF_SIZE], key[MAX_BUF_SIZE], value[MAX_BUF_SIZE];
    FILE *file;

    file = fopen(config_file_name, "r");

    while((fgets(buf, MAX_BUF_SIZE, file)) != NULL) {
        if(strlen(buf) <= 0) 
            continue;
        
        //Assume no comments, only key:value pair each line
        sscanf(buf, "%[^:]: %s", key, value);
        parse(key, value, cfg, file);
    }
}

static void
parse(char* key, char* value, struct region* cfg, FILE* file)
{
    if(strcmp(key, "TimeSlots") == 0) {
        cfg->timeslots = atoi(value);
    } else if(strcmp(key, "TimeStep") == 0) {
        cfg->config_timestep = atof(value);
        cfg->timestep = cfg->config_timestep;
    } else if(strcmp(key, "Horizon") == 0 ) {
        cfg->horizon = atoi(value);
    } else if(strcmp(key, "GridSize") == 0) {
        cfg->gridsize = atoi(value);
    } else if(strcmp(key, "NumberOfSmallParticles") == 0) {
        cfg->num_small = atof(value);
    } else if(strcmp(key, "SmallParticleMass")==0) {
        cfg->small_mass = atof(value);
    } else if(strcmp(key, "SmallParticleRadius") == 0) {
        cfg->small_radius = atof(value);
    } else if(strcmp(key, "CollisionMode") == 0) {
        cfg->collision_mode = atoi(value);
    } else if(strcmp(key, "NumberOfLargeParticles") == 0) {
        cfg->num_large = atoi(value);
        cfg->cap_large = 2 * cfg->num_large;
        cfg->lgps = malloc(cfg->cap_large * sizeof(struct particle));
        if(cfg->lgps == NULL) {
            printf("Error: cannot allocate memory - cfg-lgps");
            exit(1);
        }
        parse_large_p(cfg->lgps, cfg->num_large, file);
    }
    else {
        printf("Unsupported key: %s, value %s\n", key, value);
    }
}

static void
parse_large_p(struct particle *lgps, int num_p, FILE* file) 
{
    char buf[MAX_BUF_SIZE];
    int i = 0;
   while(i<num_p){
       if((fgets(buf, MAX_BUF_SIZE, file))!= NULL) {
            sscanf(buf, "%lf %lf %lf %lf", &lgps[i].r, &lgps[i].m, &lgps[i].locx, &lgps[i].locy);
            lgps[i].vx = lgps[i].vy = lgps[i].fx = lgps[i].fy = 0.0;
            lgps[i].is_large = true;
       }
       i++;
   }
}


/**
 * Determines the current time
 * Credit: Lab 2
 **/
long long wall_clock_time()
{
#ifdef LINUX
	struct timespec tp;
	clock_gettime(CLOCK_REALTIME, &tp);
	return (long long)(tp.tv_nsec + (long long)tp.tv_sec * 1000000000ll);
#else
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (long long)(tv.tv_usec * 1000 + (long long)tv.tv_sec * 1000000000ll);
#endif
}


struct m_time*  init_timer(char name[32], int id) 
{
    struct m_time* t;

    t = malloc(sizeof(struct m_time));
    if(t == NULL) {
        fprintf(stderr, "Out of memory: timer\n");
    }
    t->id = id; 
    t->total_time  = t->start = 0;
    t->min_time = t->max_time = -1;
    t->sample = 0;
    strncpy(t->name, name, 32);

    return t;
}

void  start_timer(struct m_time* t)
{
    t->start = wall_clock_time(); 
}

void pause_timer(struct m_time* t)
{
    long long duration;
    duration = wall_clock_time() - t->start;
    t->sample++;
    t->total_time += duration;
}

void free_timer(struct m_time *t)
{
    free(t);
}

void print_timer(struct m_time* t) 
{
    printf("%6.3f\n",
            t->total_time / 1000000000.0);
}


double solve(double a, double b, double c)
{
    double delta, t1, t2;

    delta = b*b - 4*a*c;
    if(delta < 0) {
        return -1;
    }
    
    t1 = (-1*b + sqrt(delta)) / (2*a);
    t2 = (-1*b - sqrt(delta)) / (2*a);
    
    if(t1 < 0) {
        return t2;
    } else if (t2 < 0) {
        return t1;
    } else {
        return (t1 < t2) ? t1: t2;
    }
}



