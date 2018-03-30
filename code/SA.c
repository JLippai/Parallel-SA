/*****************************************************************************/
// gcc -O1 -pthread -fopenmp -o test_SA SA_working.c -lrt -lm
// export NUM_OMP_THREADS=4

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>

#define GIG 1000000000
#define CPG 2.9           // Cycles per GHz -- Adjust to your computer

#define BASE  5
#define ITERS 20
#define DELTA 5

#define OPTIONS 3
#define BLOCK_SIZE 100     // TO BE DETERMINED

#define MINVAL   0.0
#define MAXVAL  10.0

#define TOL 0.00001
#define OMEGA 1.99      // TO BE DETERMINED

#define LAMBDA 20
#define MAX_DEVICE_HEIGHT 76200
#define MAX_DEVICE_WIDTH 76200
#define N_COMPONENTS 7
#define MOVES_PER_TEMP 2
#define SPACING 0
#define EDGE_DENSITY 0.2
#define OVERLAP_COST 10000

#define max(a,b) \
  ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b; })

#define min(a,b) \
  ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
    _a > _b ? _b : _a; })

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
  comp_ptr components;
} graph, *graph_ptr;

typedef struct {
  int width;//  MAX_DEVICE_WIDTH/BLOCK_SIZE;
  int height;// = MAX_DEVICE_HEIGHT/BLOCK_SIZE;
  int v;
  int grid[MAX_DEVICE_WIDTH/BLOCK_SIZE][MAX_DEVICE_HEIGHT/BLOCK_SIZE][N_COMPONENTS];
} layoutGrid, *grid_ptr;

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

  //g = (struct graph*) malloc(sizeof(struct graph));
  int convergence[OPTIONS][ITERS+1];
  vec_ptr new_vec(long int len);
  int set_vec_length(vec_ptr v, long int index);
  long int get_vec_length(vec_ptr v);
  int init_vector(vec_ptr v, long int len);
  int init_vector_rand(vec_ptr v, long int len);
  int init_device(graph_ptr g, long int len);
  long long int cost(graph_ptr g);
  int print_vector(vec_ptr v);
  void SAPlace_reorg(graph_ptr device_base, int* iterations);
  void SAPlace(graph_ptr device_base, int* iterations);
  void omp_SAPlace(graph_ptr device_base, int* iterations);
  graph_ptr new_device(long int len);
  void print_device(grid_ptr g);
  //int leftEdge(component c);

  printf("Version: 0\n");
  int *iterations;

  time_t t_seed;
  srand((unsigned) time(&t_seed));

  long int i, j, k;
  long int time_sec, time_ns;
  long int MAXSIZE = BASE+(ITERS+1)*DELTA;

  printf("\n Hello World -- SA variations \n");

  // declare and initialize the vector structure
  //vec_ptr v0 = new_vec(MAXSIZE);
  //iterations = (int *) malloc(sizeof(int));

  graph_ptr g0 = new_device(MAXSIZE);
  //init_device(g0,N_COMPONENTS);
  //print_device(g);
  /*
  printf("Step 1\n");  
  for (i = 0; i < N_COMPONENTS; ++i) {
    printf("%d: %d, %d\n",i, g->components[i].x, g->components[i].y);
    for (j = 0; j < 4; j++)
      printf("%d: %d\n", j, g->components[i].neighbors[j]);
  }
  */
/*
  SAPlace(g);
  long long int gc = cost(g);
  printf("\n%li\n",gc);*/
  /*
  for (i = 0; i < N_COMPONENTS; ++i) {
    printf("%d: %d, %d\n",i, g->components[i].x, g->components[i].y);
    for (j = 0; j < 4; j++)
      printf("%d: %d\n", j, g->components[i].neighbors[j]);
  }  
  */
  iterations = (int *) malloc(sizeof(int));
  OPTION = 0;
  for (i = 0; i < ITERS; i++) {
    printf("\niter = %d length = %d EDGE RHO = %0.2f", i, BASE+(i+1)*DELTA, EDGE_DENSITY);
    set_device_length(g0, BASE+(i+1)*DELTA);
    init_device(g0, BASE+(i+1)*DELTA);
    //print_vector(v0);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    //printf("\nCheck 0\n");
    SAPlace(g0,iterations);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
    convergence[OPTION][i] = *iterations;
    //print_vector(v0);
  }

  OPTION++;
  for (i = 0; i < ITERS; i++) {
    printf("\niter = %d length = %d EDGE RHO = %0.2f", i, BASE+(i+1)*DELTA, EDGE_DENSITY);
    set_device_length(g0, BASE+(i+1)*DELTA);
    init_device(g0, BASE+(i+1)*DELTA);
    //print_vector(v0);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    //printf("\nCheck 0\n");
    SAPlace_reorg(g0,iterations);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
    convergence[OPTION][i] = *iterations;
    //print_vector(v0);
  }

  OPTION++;
  for (i = 0; i < ITERS; i++) {
    printf("\niter = %d length = %d EDGE RHO = %0.2f", i, BASE+(i+1)*DELTA, EDGE_DENSITY);
    set_device_length(g0, BASE+(i+1)*DELTA);
    init_device(g0, BASE+(i+1)*DELTA);
    //print_vector(v0);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    //printf("\nCheck 0\n");
    omp_SAPlace(g0,iterations);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    time_stamp[OPTION][i] = diff(time1,time2);
    convergence[OPTION][i] = *iterations;
    //print_vector(v0);
  }    

  for (i = 0; i < ITERS; i++) {
    printf("\n%d, ", BASE+(i+1)*DELTA);
    for (j = 0; j < OPTIONS; j++) {
      if (j != 0) printf(", ");
      printf("%ld", (long int)((double)(CPG)*(double)
     (GIG * time_stamp[j][i].tv_sec + time_stamp[j][i].tv_nsec)));
      printf(", %d", convergence[j][i]);
    }
  }

  printf("\n");
  
}/* end main */
/*********************************/

struct timespec diff(struct timespec start, struct timespec end)
{
  struct timespec temp;
  if ((end.tv_nsec-start.tv_nsec)<0) {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return temp;
}

double fRand(double fMin, double fMax)
{
    double f = (double)random() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

/************************************/

vec_ptr new_vec(long int len)
{
  long int i;

  /* Allocate and declare header structure */
  vec_ptr result = (vec_ptr) malloc(sizeof(vec_rec));
  if (!result) return NULL;  /* Couldn't allocate storage */
  result->len = len;

  /* Allocate and declare array */
  if (len > 0) {
    data_t *data = (data_t *) calloc(len*len, sizeof(data_t));
    if (!data) {
    free((void *) result);
    printf("\n COULDN'T ALLOCATE STORAGE \n", result->len);
    return NULL;  /* Couldn't allocate storage */
  }
  result->data = data;
  }
  else result->data = NULL;

  return result;
}


int leftEdge(component c) {
  return (c.x - SPACING)/BLOCK_SIZE;
}

int topEdge(component c) {
  return (c.y - SPACING)/BLOCK_SIZE;
}

int rightEdge(component c) {
  return (c.x + c.width + SPACING)/BLOCK_SIZE;
}

int bottomEdge(component c) {
  return (c.y + c.height + SPACING)/BLOCK_SIZE;
}

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

int set_device_length(graph_ptr g, long int index)
{
  g->v = index;
  return 1;
}

int overlaps_boolarith(component c1, component c2) {
  return (!(leftEdge(c1) > rightEdge(c2)) && !(rightEdge(c1) < leftEdge(c2)) && !(bottomEdge(c1) < topEdge(c2)) && !(topEdge(c1) > bottomEdge(c2))) ? 1 : 0;
}

int overlapArea(component c1, component c2) {
  if (!overlaps_boolarith(c1, c2)) {
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


long long int overlapCost(graph_ptr device) {
  long int totalCost = 0;
  int i, j;
  for (i = 0; i < device->v; i++) {
    for (j = i+1; j < device->v; j++) {
      totalCost += overlapArea(device->components[i], device->components[j])*OVERLAP_COST;
        //printf("Contribution from overlap: %li\n", overlapArea(device->components[i], device->components[device->components[i].neighbors[j]])*OVERLAP_COST);
    }
  }
  return totalCost;
}

long long int overlapCost_unrolled(graph_ptr device) {
  long int totalCost, acc1;
  int i, j;
  for (i = 0; i < device->v; i++) {
    for (j = i+1; j < device->v -1; j+=2) {
      totalCost += overlapArea(device->components[i], device->components[j]);
      acc1 += overlapArea(device->components[i], device->components[j+1]);
        //printf("Contribution from overlap: %li\n", overlapArea(device->components[i], device->components[device->components[i].neighbors[j]])*OVERLAP_COST);
    }
    for (; j < device->v; j++)
      totalCost += overlapArea(device->components[i], device->components[j]);
  }
  totalCost = OVERLAP_COST*(totalCost + acc1);
  return totalCost + acc1;
}

long long int overlapCost_unrolled4(graph_ptr device) {
  long int totalCost = 0;
  long long int acc1 = 0;
  long long int acc2 = 0;
  long long int acc3 = 0;
  int i, j;
  for (i = 0; i < device->v; i++) {
    for (j = i+1; j < device->v -3; j+=4) {
      totalCost += overlapArea(device->components[i], device->components[j]);
      acc1 += overlapArea(device->components[i], device->components[j+1]);
      acc2 += overlapArea(device->components[i], device->components[j+2]);
      acc3 += overlapArea(device->components[i], device->components[j+3]);
        //printf("Contribution from overlap: %li\n", overlapArea(device->components[i], device->components[device->components[i].neighbors[j]])*OVERLAP_COST);
    }
    for (; j < device->v ; j++)
      totalCost += overlapArea(device->components[i], device->components[j]);
  }
  totalCost = OVERLAP_COST*(totalCost + acc1 + acc2 + acc3);
  return totalCost;
}

long int channelPenalty(graph_ptr device)
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

    long int cumulative_cost = 0; 
    int single_channel_cost = 0; 
    int single_overlap_cost = 0;
    int cost_per_unit_length = 2;
    int i,j,k;
    char meport;
    char youport;

    for (i = 0; i < device->v; i++)                              // Iterate over all components in device. 
    {
        for (j = 0; j < 4; j++)                                     // For each component, iterate over all ports. 
        {
            meport  = '0';
            youport = '0';
            single_channel_cost = 0;
            single_overlap_cost = 0;

            if ( (device->components[i].neighbors[j]) > i )        // Check to make sure we have not already counted this components channel contribution.
            {
                int you_id = device->components[i].neighbors[j];   // neighbors component id (to lookup/refference component list)
                int me_id  = i;                                     // current component id   (to lookup/reffernece component list)
                int me_x, you_x, me_y, you_y;
                component me  = device->components[i];              // Current component 
                component you = device->components[you_id];         // Current neighbor component 

                switch (j)                                          // Find location of current components port. 
                {
                  case 0: // Top port 
                  {
                    me_x = me.x + (me.width / 2);
                    me_y = me.y;
                    meport = 't';
                  }
                  break;
                  case 1: // Left port
                  {
                    me_x = me.x;
                    me_y = me.y + (me.height / 2);
                    meport = 'l';
                  }
                  break;
                  case 2: // Right port
                  {
                    me_x = me.x + me.width;
                    me_y = me.y + (me.height / 2);
                    meport = 'r';
                  }
                  break;
                  case 3: // Bottom port 
                  {
                    me_x = me.x + (me.width / 2);
                    me_y = me.y + me.height;  
                    meport = 'b';
                  }
                  break;
                }

                for (k = 0; k < 4; k++)                       // Find location of neighbor component port. 
                {
                  if (you.neighbors[k] == me_id)
                  switch (k)
                  {
                    case 0: // Top port 
                    {
                      you_x = you.x + (you.width / 2);
                      you_y = you.y;
                      youport = 't'; 
                    }
                    break;
                    case 1: // Left port
                    {
                      you_x = you.x;
                      you_y = you.y + (you.height / 2);
                      youport = 'l';
                    }
                    break;
                    case 2: // Right port
                    {
                      you_x = you.x + you.width;
                      you_y = you.y + (you.height / 2);
                      youport = 'r';
                    }
                    break;
                    case 3: // Bottom port
                    {
                      you_x = you.x + (you.width / 2);
                      you_y = you.y + you.height;
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

long long int cost(graph_ptr device) {
  //printf("\nWhoo\n");
  return overlapCost(device) + channelPenalty(device);
}

void SAPlace(graph_ptr device_base, int* iterations) {
  int i, j, x_old, y_old;
  float temp = 0;
  graph temp_init = *device_base;
  for (i = 0; i < 20; ++i) {
    rand_place(&temp_init);
    temp += cost(&temp_init)/20;
  }

  long long int old_cost = cost(device_base);
  long long int new_cost;
  float rateAccept;
  int rangeX = MAX_DEVICE_WIDTH;
  int rangeY = MAX_DEVICE_HEIGHT;
  int acceptCount;
  int tempCount = 0;
  while ( (temp > 0.005 * old_cost / device_base->v) && (old_cost > 2)) {
  	tempCount++;
    acceptCount = 0;
    for (i = 0; i < device_base->v * MOVES_PER_TEMP; ++i) {
      j = i / MOVES_PER_TEMP;
      x_old = device_base->components[j].x;
      y_old = device_base->components[j].y;
      acceptCount += 1;
      if (rand() > 0.5*RAND_MAX)
        device_base->components[j].x = x_old + 2*(float)(rand() - RAND_MAX/2)/RAND_MAX*rangeX;
      else
        device_base->components[j].y = y_old + 2*(float)(rand() - RAND_MAX/2)/RAND_MAX*rangeY;
      new_cost = cost(device_base);
      if (new_cost > old_cost) {
        if ((float)rand()/RAND_MAX >= exp((old_cost - new_cost)/temp)) {
          device_base->components[j].x = x_old;
          device_base->components[j].y = y_old;
          acceptCount -= 1;
        }
        else
          old_cost = new_cost;
      }
      else
        old_cost = new_cost;
    }
    rateAccept = (float)acceptCount/(MOVES_PER_TEMP * device_base->v);
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
  *iterations = tempCount;
}

void SAPlace_reorg(graph_ptr device_base, int* iterations) {
  int i, j, k, x_old, y_old;
  float temp = 0;
  graph temp_init = *device_base;
  for (i = 0; i < 20; ++i) {
    rand_place(&temp_init);
    temp += cost(&temp_init)/20;
  }
  long long int old_cost = cost(device_base);
  long long int new_cost;
  float rateAccept;
  int rangeX = MAX_DEVICE_WIDTH;
  int rangeY = MAX_DEVICE_HEIGHT;
  int acceptCount;
  int tempCount = 0;
  while ((temp > 0.005 * old_cost / device_base->v) && (old_cost > 2)) {
    tempCount++;
    acceptCount = 0;
    for (i = 0; i < device_base->v; ++i) {
      for (j = 0; j < MOVES_PER_TEMP; ++j) {
        x_old = device_base->components[i].x;
        y_old = device_base->components[i].y;
        acceptCount += 1;
        if (rand() > 0.5*RAND_MAX)
          device_base->components[i].x = x_old + 2*(float)(rand() - RAND_MAX/2)/RAND_MAX*rangeX;
        else
          device_base->components[i].y = y_old + 2*(float)(rand() - RAND_MAX/2)/RAND_MAX*rangeY;
        new_cost = cost(device_base);
        if (new_cost > old_cost) {
          if ((float)rand()/RAND_MAX >= exp((old_cost - new_cost)/temp)) {
            device_base->components[i].x = x_old;
            device_base->components[i].y = y_old;
            acceptCount -= 1;
          }
          else {
            old_cost = new_cost;
          }
        }
        else {
          old_cost = new_cost;
        }
      }
    }
    rateAccept = (float)acceptCount/(MOVES_PER_TEMP * device_base->v);
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
  *iterations = tempCount;
}

graph_ptr new_device(long int len) {
  long int i;

  /* Allocate and declare header structure */
  graph_ptr result = (graph_ptr) malloc(sizeof(graph));
  if (!result) return NULL;
  result->v = len;

  /* Allocate and declare array */
  if (len > 0) {
    comp_ptr components = (comp_ptr) calloc(len, sizeof(component));
    if (!components) {
      free((void *) result);
      printf("\n COULDN'T ALLOCATE STORAGE \n", result-> v);
      return NULL;
    }
    result->components = components;
  }
  else result->components = NULL;
  return result;
}

void omp_SAPlace(graph_ptr device_base, int* iterations) {
  int i, j, x_old, y_old;
  float temp = 0;
  graph temp_init = *device_base;
  for (i = 0; i < 20; ++i) {
    rand_place(&temp_init);
    temp += cost(&temp_init)/20;
  }

  long long int old_cost = cost(device_base);
  long long int new_cost;
  float rateAccept;
  int rangeX = MAX_DEVICE_WIDTH;
  int rangeY = MAX_DEVICE_HEIGHT;
  int acceptCount;
  int tempCount = 0;
  int acceptChange;
  int perturbX, perturbY;
  graph device_copy;

  while ( (temp > 0.005 * old_cost / device_base->v) && (old_cost > 2)) {
    tempCount++;
    acceptCount = 0;
    acceptChange = 0;
    device_copy = *device_base;
    old_cost = cost(device_base);
    omp_set_num_threads(4);
    #pragma omp parallel for firstprivate(acceptChange, old_cost, device_copy) private(j, new_cost, perturbX, perturbY)
    for (i = 0; i < device_base->v * MOVES_PER_TEMP; ++i) {
      j = i / MOVES_PER_TEMP;
      acceptChange += 1;
      perturbX = (rand() > 0.5) ? 2*(float)(rand() - RAND_MAX/2)/RAND_MAX*rangeX : 0;
      perturbY = (rand() > 0.5) ? 0 : 2*(float)(rand() - RAND_MAX/2)/RAND_MAX*rangeY;
      device_copy.components[j].x += perturbX;
      device_copy.components[j].y += perturbY;
      new_cost = cost(&device_copy);
      if (new_cost > old_cost) {
        if ((float)rand()/RAND_MAX >= exp((old_cost - new_cost)/temp)) {
          perturbX = 0;
          perturbY = 0;
          acceptChange -= 1;
        }
      }
      #pragma omp critical
        acceptCount += acceptChange;
      #pragma omp critical
        device_base->components[j].x += perturbX;
      #pragma omp critical
        device_base->components[j].y += perturbY;
      
    }
    rateAccept = (float)acceptCount/(MOVES_PER_TEMP * device_base->v);
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
  *iterations = tempCount;
}

int init_device(graph_ptr g, long int len) {
  long int i, j;
  float randc;
  long int randn;
  long int randp;
  if (len > 0) {
    g->v = len;
    for (i = 0; i < len; i++) {
      g->components[i].x = rand()%MAX_DEVICE_WIDTH;
      g->components[i].y = rand()%MAX_DEVICE_HEIGHT;
      g->components[i].width = MAX_DEVICE_WIDTH/len;
      g->components[i].height = MAX_DEVICE_HEIGHT/len;
      for (j = 0; j < 4; j++)
        g->components[i].neighbors[j] = -1;
    }
    for (i = 0; i < len; ++i) {
      for (j = 0; j < 4; ++j) {
        if (g->components[i].neighbors[j] == -1) {
          randc = ((float)rand())/RAND_MAX;
          if (randc < EDGE_DENSITY) {
            //printf("Booyah!\n");
            randn = rand() % len;
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

int rand_place(graph_ptr g) {
  long int i;
  for (i = 0; i < g->v; i++) {
    g->components[i].x = rand()%MAX_DEVICE_WIDTH;
    g->components[i].y = rand()%MAX_DEVICE_HEIGHT;
  }
  
}


void print_device(grid_ptr g) 
{
  int i,j,k;
  int tally;
  for (i = 0; i < MAX_DEVICE_WIDTH/BLOCK_SIZE; i++)
  {
    printf("\n");
    for(j = 0; j < MAX_DEVICE_HEIGHT/BLOCK_SIZE; j++)
    {
      tally = 0;
      for (k = 0; k < N_COMPONENTS; k++)
      {
        if (g->grid[i][j][k] == 1)
          tally++;
      }
        if (tally == 0)
          printf("%s", ".");
        else 
          printf("%d", tally);
      printf("%s"," ");
    }
  }
} 

void fill_grid(graph_ptr g, grid_ptr grid)
{
  int i,j,k;
  int x, y, w, h;

  for (i = 0; i < N_COMPONENTS; i++)
  {
    x = g->components[i].x;
    y = g->components[i].y;
    w = g->components[i].width;
    h = g->components[i].height;

    grid->grid[x/BLOCK_SIZE][y/BLOCK_SIZE][i] = 1;
    for (k = 0; k < w/BLOCK_SIZE; k++)
      grid->grid[x/BLOCK_SIZE + k][y/BLOCK_SIZE][i] = 1;
    for (j = 0; j < h/BLOCK_SIZE; j++)
      grid->grid[x/BLOCK_SIZE][y/BLOCK_SIZE + h][i] = 1;
  }
}