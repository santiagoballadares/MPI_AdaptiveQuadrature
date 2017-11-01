#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>
#include "stack.h"

#define EPSILON 1e-3
#define F(arg)  cosh(arg)*cosh(arg)*cosh(arg)*cosh(arg)
#define A 0.0
#define B 5.0

#define FARMER 0

#define TAG_NO_MORE_TASKS  99
#define TAG_NEW_TASK        0
#define TAG_RESULT          1
#define TAG_DO_WORK         2

#define SLEEPTIME 1


int *tasks_per_process;

double  farmer(int);

void    worker(int);


int main(int argc, char **argv) {
  int i, my_id, num_procs;
  double area, a, b;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD,&my_id);

  if(num_procs < 2) {
    fprintf(stderr, "ERROR: Must have at least 2 processes to run\n");
    MPI_Finalize();
    exit(1);
  }

  if (my_id == FARMER) { // Farmer
    // init counters
    tasks_per_process = (int *) malloc(sizeof(int)*(num_procs));
    for (i=0; i<num_procs; ++i) {
      tasks_per_process[i] = 0;
    }
  }

  if (my_id == FARMER) { // Farmer
    area = farmer(num_procs);
  } else { //Workers
    worker(my_id);
  }

  if(my_id == FARMER) { // Farmer
    fprintf(stdout, "Area=%lf\n", area);
    fprintf(stdout, "\nTasks Per Process\n");
    for (i=0; i<num_procs; i++) {
      fprintf(stdout, "%d\t", i);
    }
    fprintf(stdout, "\n");
    for (i=0; i<num_procs; i++) {
      fprintf(stdout, "%d\t", tasks_per_process[i]);
    }
    fprintf(stdout, "\n");
    free(tasks_per_process);
  }
  MPI_Finalize();
  return 0;
}


double farmer(int num_procs)
{
  MPI_Status status;

  int i;
  double total = 0.0;

  // Workers array: idle: 0 | Busy: 1
  // num_idle keeps track of the number of idle workers
  int num_workers = num_procs - 1;
  int num_idle = num_workers;
  int *workers = (int *) malloc(sizeof(int)*num_workers);
  for (i=0; i<num_workers; ++i)
  {
    workers[i] = 0;
  }

  // Bag of tasks (Stack structure) & initial task
  stack *mybag = new_stack();
  double task[2] = {A, B};
  push(task, mybag);

  // Start main loop
  while ( !is_empty(mybag) || num_idle != num_workers )
  {
	// Distribute tasks
    i = 0;
    while ( !is_empty(mybag) && num_idle > 0 )
    {
	  // Look for idle workers
      if ( workers[i] == 0 )
      {
        // Send task
        MPI_Ssend(pop(mybag), 2, MPI_DOUBLE, i+1, TAG_DO_WORK, MPI_COMM_WORLD);

        // Update workers array, idle counter, and tasks_per_process array
        workers[i] = 1;
        num_idle--;
        tasks_per_process[i+1]++;
      }

	  // Update index for workers array
      i = (i < num_workers) ? i+1 : 0;
    }

    // Receive result or new task 1
    MPI_Recv(task, 2, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    if ( status.MPI_TAG == TAG_RESULT )
    {
      // Accumulate result
      total += task[0];
    }
    else if ( status.MPI_TAG == TAG_NEW_TASK )
    {
      // Push task 1
      push(task, mybag);

      // Use status.MPI_SOURCE and TAG_NEW_TASK to get task 2
      MPI_Recv(task, 2, MPI_DOUBLE, status.MPI_SOURCE, TAG_NEW_TASK, MPI_COMM_WORLD, &status);

      // Push task 2
      push(task, mybag);
    }

    // Update workers array and idle counter
    workers[status.MPI_SOURCE - 1] = 0;
    num_idle++;
  }

  // Tell workers to stop execution
  for (i=0; i<num_workers; ++i)
  {
    MPI_Ssend(task, 2, MPI_DOUBLE, i+1, TAG_NO_MORE_TASKS, MPI_COMM_WORLD);
  }
  
  // Free memory
  free(workers);

  return total;
}


void worker(int mypid)
{
  MPI_Status status;

  int tag;
  double task[2], left, right, fleft, fright, lrarea, mid, fmid, larea, rarea;
  task[0] = task[1] = left = right = fleft = fright = lrarea = mid = fmid = larea = rarea = 0.0;

  // Receive first task
  MPI_Recv(task, 2, MPI_DOUBLE, FARMER, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  tag = status.MPI_TAG;

  while ( tag != TAG_NO_MORE_TASKS )
  {
    // Calculate the Adaptive Quadrature
    left = task[0];
    right = task[1];
    fleft = F(left);
    fright = F(right);
    lrarea = (fleft + fright) * (right-left) / 2.0;

    mid = (left + right) / 2.0;
    fmid = F(mid);
    larea = (fleft + fmid) * (mid - left) / 2.0;
    rarea = (fmid + fright) * (right - mid) / 2.0;

    // Call sleep to simulate processing in multiple cores
    usleep(SLEEPTIME);

    // Evaluate approximation
    if ( fabs((larea + rarea) - lrarea) > EPSILON )
    {
      // Send new task 1
      task[0] = left;
      task[1] = mid;
      MPI_Ssend(task, 2, MPI_DOUBLE, FARMER, TAG_NEW_TASK, MPI_COMM_WORLD);

      // Send new task 2
      task[0] = mid;
      task[1] = right;
      MPI_Ssend(task, 2, MPI_DOUBLE, FARMER, TAG_NEW_TASK, MPI_COMM_WORLD);
    }
    else
    {
      // Send result
      task[0] = larea + rarea;
      task[1] = 0.0;
      MPI_Ssend(task, 2, MPI_DOUBLE, FARMER, TAG_RESULT, MPI_COMM_WORLD);
    }

    // Receive another task or terminate message
    MPI_Recv(task, 2, MPI_DOUBLE, FARMER, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    tag = status.MPI_TAG;
  }
}
