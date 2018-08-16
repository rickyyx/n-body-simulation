#ifndef UTIL_H
#define UTIL_H

#include <time.h>
#include "main.h"


/*  Helper Utility */
void debug_event(struct event*);
void slave_read_config(struct region *);
void master_read_config(struct universe *);
void print_config(struct region *);
void debug_buf_1d(unsigned char*, int);
void debug_region(struct region*);
void p_canvas(unsigned char**, int);
void p_result(unsigned char*, int, int);
void free_buf_2d(unsigned char** ptr, size_t size);
void allocate_buf_1d(unsigned char** ptr, size_t size);
void allocate_buf_2d(unsigned char*** ptr, size_t x_size, size_t y_size);




/*  Time Utility  */
struct m_time {
    long long total_time;
    long long min_time;
    long long max_time;
    unsigned int sample;
    long long start;
    char name[32];
    int id;
};
long long wall_clock_time();
struct m_time* init_timer(char*, int);
void start_timer(struct m_time*);
void pause_timer(struct m_time*);
void free_timer(struct m_time*);
void print_timer(struct m_time*);

/*  Math Utility */
double solve(double, double, double);

#endif
