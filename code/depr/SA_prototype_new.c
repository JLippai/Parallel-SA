/*****************************************************************************/
// gcc -O1 -pthread -fopenmp -o test_SOR_omp test_SOR_omp.c -lrt -lm
// export NUM_OMP_THREADS=4

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>

#define GIG 1000000000
#define CPG 2.9           // Cycles per GHz -- Adjust to your computer

#define BASE  2
#define ITERS 1
#define DELTA 998

#define OPTIONS 10
#define BLOCK_SIZE 100     // TO BE DETERMINED

#define MINVAL   0.0
#define MAXVAL  10.0

#define TOL 0.00001
#define OMEGA 1.99      // TO BE DETERMINED

#define LAMBDA 20
#define MAX_DEVICE_HEIGHT 135000
#define MAX_DEVICE_WIDTH 80000
#define N_COMPONENTS 20
#define MOVES_PER_TEMP 2;
#define SPACING 0;
#define EDGE_DENSITY 0.2;

typedef int data_t;

/* Create abstract data type for vector -- here a 2D array */
typedef struct {
  long int len;
  data_t *data;
} vec_rec, *vec_ptr;

typedef struct {
  int neighbors[4];
  int x, y, width, height;
} component, *comp_ptr;

typedef struct {
  int v;
  component components[N_COMPONENTS];
} graph, *graph_ptr;

typedef struct {
  int grid[MAX_DEVICE_WIDTH/BLOCK_SIZE][MAX_DEVICE_HEIGHT/BLOCK_SIZE][N_COMPONENTS];
} layoutGrid, *grid_ptr;

int leftEdge(component c);

int NUM_THREADS = 4;

/* Number of bytes in a vector (SSE sense) */
#define VBYTES 16

/* Number of elements in a vector (SSE sense) */
#define VSIZE VBYTES/sizeof(data_t)

typedef data_t vec_t __attribute__ ((vector_size(VBYTES)));
typedef union {
  vec_t v;
  data_t d[VSIZE];
} pack_t;

/* used to pass parameters to worker threads */
struct thread_data{
  int thread_id;
  vec_ptr v;
  int* iterations;
  double* mean_change;
};

/* declare barrier and mutex lock*/
pthread_barrier_t barr;
pthread_mutex_t mean_lock;

/*****************************************************************************/
main(int argc, char *argv[])
{
  int OPTION;
  struct timespec diff(struct timespec start, struct timespec end);
  struct timespec time1, time2;
  struct timespec time_stamp[OPTIONS][ITERS+1];
  graph g;
  int convergence[OPTIONS][ITERS+1];
  vec_ptr new_vec(long int len);
  int set_vec_length(vec_ptr v, long int index);
  long int get_vec_length(vec_ptr v);
  int init_vector(vec_ptr v, long int len);
  int init_vector_rand(vec_ptr v, long int len);
  int init_device_rand(graph_ptr g, long int len);
  int print_vector(vec_ptr v);
  int *iterations;

  time_t t_seed;
  srand((unsigned) time(&t_seed));

  long int i, j, k;
  long int time_sec, time_ns;
  long int MAXSIZE = BASE+(ITERS+1)*DELTA;

  printf("\n Hello World -- SOR variations \n");

  // declare and initialize the vector structure
  vec_ptr v0 = new_vec(MAXSIZE);
  iterations = (int *) malloc(sizeof(int));
/*
  OPTION = 0;
  for (i = 0; i < ITERS; i++) {
    printf("\niter = %d length = %d OMEGA = %0.2f", i, BASE+(i+1)*DELTA, OMEGA);
    set_vec_length(v0, BASE+(i+1)*DELTA);
    init_vector_rand(v0, BASE+(i+1)*DELTA);
    //print_vector(v0);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    SOR(v0,iterations);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
    convergence[OPTION][i] = *iterations;
    //print_vector(v0);
  }

  OPTION++;
  for (i = 0; i < ITERS; i++) {
    printf("\niter = %d length = %d", i, BASE+(i+1)*DELTA);
    init_vector_rand(v0, BASE+(i+1)*DELTA);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    set_vec_length(v0,BASE+(i+1)*DELTA);
    SOR_ji(v0,iterations);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
    convergence[OPTION][i] = *iterations;
  }

  OPTION++;
  for (i = 0; i < ITERS; i++) {
    printf("\niter = %d length = %d", i, BASE+(i+1)*DELTA);
    init_vector_rand(v0, BASE+(i+1)*DELTA);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    set_vec_length(v0,BASE+(i+1)*DELTA);
    SOR_blocked(v0,iterations);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
    convergence[OPTION][i] = *iterations;
  }

  OPTION++;
  for (i = 0; i < ITERS; i++) {
    printf("\niter = %d length = %d", i, BASE+(i+1)*DELTA);
    init_vector_rand(v0, BASE+(i+1)*DELTA);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    set_vec_length(v0,BASE+(i+1)*DELTA);
    pt_SOR(v0,iterations);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
    convergence[OPTION][i] = *iterations;
  }

  OPTION++;
  for (i = 0; i < ITERS; i++) {
    printf("\niter = %d length = %d", i, BASE+(i+1)*DELTA);
    init_vector_rand(v0, BASE+(i+1)*DELTA);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    set_vec_length(v0,BASE+(i+1)*DELTA);
    pt_SOR_otherstripe(v0,iterations);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
    convergence[OPTION][i] = *iterations;
  }  

  OPTION++;
  for (i = 0; i < ITERS; i++) {
    printf("\niter = %d length = %d", i, BASE+(i+1)*DELTA);
    init_vector_rand(v0, BASE+(i+1)*DELTA);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    set_vec_length(v0,BASE+(i+1)*DELTA);
    pt_SOR_alternating(v0,iterations);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
    convergence[OPTION][i] = *iterations;
  }  

  NUM_THREADS = 2;
  OPTION++;
  for (i = 0; i < ITERS; i++) {
    printf("\niter = %d length = %d", i, BASE+(i+1)*DELTA);
    init_vector_rand(v0, BASE+(i+1)*DELTA);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    set_vec_length(v0,BASE+(i+1)*DELTA);
    pt_SOR_rb(v0,iterations);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
    convergence[OPTION][i] = *iterations;
  }  

  OPTION++;
  for (i = 0; i < ITERS; i++) {
    printf("\niter = %d length = %d", i, BASE+(i+1)*DELTA);
    init_vector_rand(v0, BASE+(i+1)*DELTA);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    set_vec_length(v0,BASE+(i+1)*DELTA);
    omp_SOR(v0,iterations);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
    convergence[OPTION][i] = *iterations;
  }    
*/
/*  OPTION++;
  for (i = 0; i < ITERS; i++) {
    printf("\niter = %d length = %d", i, BASE+(i+1)*DELTA);
    init_vector_rand(v0, BASE+(i+1)*DELTA);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    set_vec_length(v0,BASE+(i+1)*DELTA);
    omp_SOR_ji(v0,iterations);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
    convergence[OPTION][i] = *iterations;
  }     */ 
/*
  OPTION++;
  for (i = 0; i < ITERS; i++) {
    printf("\niter = %d length = %d", i, BASE+(i+1)*DELTA);
    init_vector_rand(v0, BASE+(i+1)*DELTA);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    set_vec_length(v0,BASE+(i+1)*DELTA);
    omp_SOR_blocked(v0,iterations);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
    convergence[OPTION][i] = *iterations;
  }    
*/
/*  OPTION++;
  for (i = 0; i < ITERS; i++) {
    printf("\niter = %d length = %d", i, BASE+(i+1)*DELTA);
    init_vector_rand(v0, BASE+(i+1)*DELTA);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    set_vec_length(v0,BASE+(i+1)*DELTA);
    omp_SOR_striped(v0,iterations);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
    convergence[OPTION][i] = *iterations;
  }
  init_device_rand(g);  */       
/*
  for (i = 0; i < ITERS; i++) {
    printf("\n%d, ", BASE+(i+1)*DELTA);
    for (j = 0; j < OPTIONS; j++) {
      if (j != 0) printf(", ");
      printf("%ld", (long int)((double)(CPG)*(double)
     (GIG * time_stamp[j][i].tv_sec + time_stamp[j][i].tv_nsec)));
      printf(", %d", convergence[j][i]);
    }
  }
*/
  printf("\n");
  
}/* end main */
/*********************************/
/*
public long calcCompOverlap(component randC) {
  long overlapSum = 0;
  int minX = leftEdge(randC) / blockSize;
  int maxX = rightEdge(randC) / blockSize;
  int minY = topEdge(randC) / blockSize;
  int maxY = bottomEdge(randC) / blockSize;
  for (int i = minX; i <= maxX; i++) {
    for (int j = minY; j <= maxY; j++) {
      long key = GridUtil.hash(i, j);
      List<Component> l = layoutGrid.get(key);
            //System.out.println("getting key " + key);
            //System.out.println(l);
      if (l == null) {
        continue;
      } else {
        for (Component c : l) {
          if (randC != c) {
            overlapSum += overlapArea(randC, c);
          }
        }
      }
    }

  }
  return overlapSum;
}

float calcOverlap() {

  float overlapSum = 0;
  for (List<Component> l : layoutGrid.values()) {
    for (int i = 0; i < l.size(); i++) {
      Component c1 = l.get(i);
      for (int j = i + 1; j < l.size(); j++) {
        Component c2 = l.get(j);
        overlapSum += overlapArea(c1, c2);
      }
    }
  }
  return overlapSum;
}
*/
int overlaps(component c1, component c2) {
  if (leftEdge(c1) > rightEdge(c2)) {
    return 0;
  } else if (rightEdge(c1) < leftEdge(c2)) {
    return 0;
  } else if (bottomEdge(c1) < topEdge(c2)) {
    return 0;
  } else if (topEdge(c1) > bottomEdge(c2)) {
    return 0;
  } else {
    return 1;
  }
}
/*
float calcCompOverlap(component randC) {
  float overlapSum = 0;
  int blockSize = BLOCK_SIZE;
  int minX = leftEdge(randC) / blockSize;
  int maxX = rightEdge(randC) / blockSize;
  int minY = topEdge(randC) / blockSize;
  int maxY = bottomEdge(randC) / blockSize;
  for (int i = minX; i <= maxX; i++) {
    for (int j = minY; j <= maxY; j++) {
      long key = GridUtil.hash(i, j);
      List<Component> l = layoutGrid.get(key);
        //System.out.println("getting key " + key);
        //System.out.println(l);
      if (l == null) {
        continue;
      } else {
        for (component c : l) {
          if (randC != c) {
            overlapSum += overlapArea(randC, c);
          }
        }
      }
    }
  }
  return overlapSum;
}
*/
int leftEdge(component c)
  return (c.x - SPACING);

int topEdge(component c)
  return c.y - SPACING;

int rightEdge(component c) 
  return c.x + c.width + SPACING;

int bottomEdge(component c)
  return c.y + c.height + SPACING;

int overlapArea(component c1, component c2) {
  if (!overlaps(c1, c2)) {
    return 0;
  }
  int left = max(leftEdge(c1), leftEdge(c2));
  int right = min(rightEdge(c1), rightEdge(c2));
  int top = max(topEdge(c1), topEdge(c2));
  int bottom = min(bottomEdge(c1), bottomEdge(c2));
        /*if ((bottom - top) <= routingSpacing || (right - left) <= routingSpacing) {
         return 0;
         }*/
  return (bottom - top) * (right - left);
}



void SAPlace(graph_ptr device) {
  float temp = 20 * average_cost(*device);
  float cost = cost(*device);
  float new_cost, rateAccept;
  int rangeX = MAX_DEVICE_WIDTH;
  int rangeY = MAX_DEVICE_HEIGHT;
  int i, j;
  int acceptCount = 0;
  while (temp > 0.005 * cost / N_COMPONENTS && temp > 2) {
    for (i = 0; i < N_COMPONENTS * MOVES_PER_TEMP; ++i) {
      j = i / MOVES_PER_TEMP;
      x_old = device->components[j].x;
      y_old = device->components[j].y;
      acceptCount += 1;
      if (rand() > 0.5*RAND_MAX)
        device->components[j].x = x_old + 2*(srand() - RAND_MAX/2)/RAND_MAX*rangeX;
      else
        device->components[j].y = y_old + 2*(srand() - RAND_MAX/2)/RAND_MAX*rangeY;
      new_cost = cost(*device);
      if (new_cost > cost)
        if (rand()/RAND_MAX >= exp((cost - new_cost)/temp)) {
          device->components[j].x = x_old;
          device->components[j].y = y_old;
          acceptCount -= 1;
        }
    }
    rateAccept = acceptCount/(MOVES_PER_TEMP * N_COMPONENTS);
    if (rateAccept > 0.96)
      temp *= 0.5;
    else if (rateAccept > 0.8)
      temp *= 0.9;
    else if (rateAccept > 0.15)
      temp *= 0.95;
    else
      temp *= 0.8;
    rangeX *= (1.0 - 0.44 + rateAccept);
    rangeY *= (1.0 - 0.44 + rateAccept);
  }
}

int overlapCost(graph_ptr device) {
  float totalCost = 0;
  int float cost_per_overlap = 10000;
  int i, j, k;

  for (i = 0; i < N_COMPONENTS; i++) {
    for (j = 0; j < 4; j++) {
      if (device->components[i].neighbors[j] > i) {
        totalCost += overlapArea(device->components[i], device->components[device->components[i].neighbors[j]])*10000;
      }
    }
  }
  return totalCost;
}

float channelPenalty(graph_ptr device)
{
   /* 
      As described in the Fluigi Thesis (Haiyao Huang 2016): 

      channelPenalty = channelLength*channelCost + numPenalty*overlapCost

      In other words, the total cost incurred in channel placment is the the sum of the 
      total cost of the channel lengths, plus the cost due to any overlapping channels. 

      The cost of a channel is given by the total length of channels, weighted by the cost per unit length. This is 2, as cited by Huang. 
      The cost of a overlap penalty is, reasonably, very large. Each overlap has a cost of 10000, as cited by Huang. 

      This function realizes the above. 
   */

    float cumulative_cost = 0; 
    float single_channel_cost = 0; 
    float single_overlap_cost = 0;
    float cost_per_unit_length = 2;
    int i,j,k;
    char meport;
    char youport;

    for (i = 0; i < N_COMPONENTS; i++)                              // Iterate over all components in device. 
    {
        for (j = 0; j < 4; j++)                                     // For each component, iterate over all ports. 
        {
            meport  = '';
            youport = '';
            single_channel_cost = 0;
            single_overlap_cost = 0;

            if ( (device->components[i]->neighbors[j]) > i )        // Check to make sure we have not already counted this components channel contribution.
            {
                int you_id = device->components[i]->neighbors[j];   // neighbors component id (to lookup/refference component list)
                int me_id  = i;                                     // current component id   (to lookup/reffernece component list)
                component me  = device->components[i];              // Current component 
                component you = device->components[you_id];         // Current neighbor component 

                switch (j)                                          // Find location of current components port. 
                {
                  case 0: // Top port 
                  {
                    int me_x = me.x + (me.width / 2);
                    int me_y = me.y;
                    meport = 't';
                  }
                  break;
                  case 1: // Left port
                  {
                    int me_x = me.x;
                    int me_y = me.y + (me.height / 2);
                    meport = 'l';
                  }
                  break;
                  case 2: // Right port
                  {
                    int me_x = me.x + me.width;
                    int me_y = me.y + (me.height / 2);
                    meport = 'r';
                  }
                  break;
                  case 3: // Bottom port 
                  {
                    int me_x = me.x + (me.width / 2);
                    int me_y = me.y + me.height;  
                    meport = 'b';
                  }
                  break;
                }

                for (k = 0; k < 4; k++)                       // Find location of neighbor component port. 
                {
                  if (you->neighbors[k] == me_id)
                  switch (k)
                  {
                    case 0: // Top port 
                    {
                      int you_x = you.x + (you.width / 2);
                      int you_y = you.y;
                      youport = 't'; 
                    }
                    break;
                    case 1: // Left port
                    {
                      int you_x = you.x;
                      int you_y = you.y + (you.height / 2);
                      youport = 'l';
                    }
                    break;
                    case 2: // Right port
                    {
                      int you_x = you.x + you.width;
                      int you_y = you.y + (you.height / 2);
                      youport = 'r';
                    }
                    break;
                    case 3: // Bottom port
                    {
                      int you_x = you.x + (you.width / 2);
                      int you_y = you.y + you.height;
                      youport = 'b';
                    }
                    break;
                  }
                }

                /* 
                    Recall: channelPenalty = channelLength*channelCost + numPenalty*overlapCost
                    Here, we calculate these seperate components.  
                */ 
                /* 
                    Calculate channel length cost: (channelLength*channelCost)
                    We measure the distance as a Manhattan Geometry, since channels cannot be routed diaganoly. 
                */
                single_channel_cost = (abs(you_x - me_x) + abs(you_y - me_y)) * cost_per_unit_length;

                /* 
                    Calculate overlap penalty cost: (numPenalty*overlapCost)
                    There are four types of overlaps, as described by Huang. These are checked here. 
                */
                if (meport == 't' && youport == 'b')
                  if (me_y > you_y)
                    single_overlap_cost = 10000;
                if (meport == 'b' && youport == 't')
                  if (me_y < you_y)
                    single_overlap_cost = 10000;
                if (meport == 'r' && youport == 'l')
                  if (me_x > you_x)
                    single_overlap_cost = 10000;
                if (meport == 'l' && youport == 'r')
                  if (me_x < you_x)
                    single_overlap_cost = 10000;

                cumulative_cost += (single_channel_cost + single_overlap_cost);
            }
        }
    }
    return cumulative_cost;
}

int cost(graph_ptr device) {
  return overlapCost(device) + channelPenalty(device);
}

void SAPlace_reorg(graph_ptr device_base) {
  float temp = 20 * average_cost(*device_base);
  float cost = cost(*device_base);
  float new_cost, rateAccept;
  int rangeX = MAX_DEVICE_WIDTH;
  int rangeY = MAX_DEVICE_HEIGHT;
  int i, j, k;
  int acceptCount = 0;
  graph device;
  while (temp > 0.005 * cost / N_COMPONENTS && temp > 2) {
    for (i = 0; i < N_COMPONENTS; ++i) {
      device = *device_base;
      acceptCount = 0;
      for (j = 0; j < MOVES_PER_TEMP; ++j) {
        x_old = device->components[i].x;
        y_old = device->components[i].y;
        acceptCount += 1;
        if (rand() > 0.5*RAND_MAX)
          device->components[i].x = x_old + 2*(rand() - RAND_MAX/2)/RAND_MAX*rangeX;
        else
          device->components[i].y = y_old + 2*(rand() - RAND_MAX/2)/RAND_MAX*rangeY;
        new_cost = cost(device);
        if (new_cost > cost)
          if (rand()/RAND_MAX >= exp((cost - new_cost)/temp)) {
            device->components[i].x = x_old;
            device->components[i].y = y_old;
            acceptCount -= 1;
          }
      }
      *device_base = device;
    }
    rateAccept = acceptCount/(MOVES_PER_TEMP * N_COMPONENTS);
    if (rateAccept > 0.96)
      temp *= 0.5;
    else if (rateAccept > 0.8)
      temp *= 0.9;
    else if (rateAccept > 0.15)
      temp *= 0.95;
    else
      temp *= 0.8;
    rangeX *= (1.0 - 0.44 + rateAccept);
    rangeY *= (1.0 - 0.44 + rateAccept);
  }
}

int init_device(graph_ptr g, long int len) {
  long int i;
  long int randc;
  long int randn;
  long int randp;
  if (len > 0) {
    g->v = len;
    for (i = 0; i < len; i++) {
      g->x = rand()%MAX_DEVICE_WIDTH;
      g->y = rand()%MAX_DEVICE_HEIGHT;
      g->width = MAX_DEVICE_WIDTH/N_COMPONENTS;
      g->height = MAX_DEVICE_HEIGHT/N_COMPONENTS;
      for (j = 0; j < 4; j++) {
        if (g->components[i].neighbors[j] == -1) {
          randc = rand()/RAND_MAX;
          if (randc < EDGE_DENSITY) {
            randn = rand() % N_COMPONENTS;
            g->components[i].neighbors[j] = randn;
            randp = rand() % 4;
            g->components[g->components[i].neighbors[j]].neighbors[randp] = i;
          }
        }
      }
    }
  }
  else return 0;
}