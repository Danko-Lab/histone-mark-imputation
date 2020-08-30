/*
 *
 *  Monte Carlo simulation of Pol II
 *
 *  Parameters: 
 *  - Initiation.
 *  - Release from pause.
 *  - Gene length.
 *  - Number of cells.
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>

const int tu_size = 25000;
const int duration = 1000; // seconds.
const int offset = 1000000;

const int uas_size   = 250;  // TU size.
const int uas_offset = 110;  // How far downstream is the UAS TSS?
const double uas_ratio  = 0.75; // Sense / UAS.

double p_i = 0.01;  // Probability of initiation in each time-step (1/s).
double p_r = 0.1;   // Probability of pause release (1/s).
const double p_t = 0.05;  // Probability to terminate after passing PAS.

const int k_min   = 8;   // Elongation rate. 8-50 bp / s (0.5-3 kb / min).
const int k_max   = 50;

struct polii {
 int position;
 int cell;
 int paused; // 1 if paused. Position is 30 +/- 5 bp.
 char strand;
 struct polii *next;
};

// Initiate a new Pol II.
void new_polii(struct polii **first, int cell) {
  struct polii *new = (struct polii*)malloc(sizeof(struct polii));

  // Decide the strand.
  double draw = (double)rand()/ (double)RAND_MAX;

  if (draw < uas_ratio) {
   new->position = 30; // Paused Pol II is at 30bp.
   new->cell = cell;
   new->paused = 1;
   new->strand = '+';
   new->next = *first;
  }
  else {
   new->position = 0-30-uas_offset; // Paused Pol II is at 30bp.
   new->cell = cell;
   new->paused = 1;
   new->strand = '-';
   new->next = *first;
  }

  *first = new; // Set *first to the new element.
}

// Terminate an existing Pol II.
void terminate_polii(struct polii *remove, struct polii *previous, struct polii **first) {
  if(remove == NULL || previous == NULL)
    return;

  // What happens if remove is the first Pol II?
  if (previous == NULL) { // Then we are popping the first element. This should be vanishingly rare.
   if(*first == NULL ) 
     return;

   *first = (*first)->next;
   free(remove);
  }
  else {
   previous->next = remove->next; // Set *previous.next to point to *remove.next.
   free(remove);
  }
}

void move_polii(struct polii *current) {
 double draw = (double)rand()/ (double)RAND_MAX;
 int k = k_min + (int)(draw*((double)k_max - (double)k_min));

 if(current->strand == '+')
   current->position += k; // Move Pol II.
 else 
   current->position -= k;
}

// Release a Pol II into productive elongation. 
void release_polii(struct polii *current) {
  // Release Pol II.
  current->paused=0;
  move_polii(current);
}

// Designed to print out the position of all RNA polymerase.
void in_silico_proseq(struct polii *first) {
  struct polii *current = first;
  while(current != NULL) {
    printf("chr6_ssto_hap7\t%d\t%d\tn\t1\t%c\n", current->position+offset, current->position+1+offset, current->strand);
    current = current->next;
  }
} 

int mcmc_polii() {
 srand(time(0));

 int ncells= 10000, i, t;
 struct polii *first = NULL, *previous = NULL, *current = NULL;

 double draw;

 int *paused = (int *)calloc(ncells, sizeof(int)); // Use calloc to initalize to 0.

 for(t=0;t<duration;t++) {
  // For each cell, determine whether a new Pol II is added.
  for(i=0;i<ncells;i++) {
    draw = (double)rand()/ (double)RAND_MAX;
    if(draw < p_i && paused[i] == 0) { // Don't initiatlize unless the pause position is free.
      new_polii(&first, i); // Create a new pol ii, place it in the pause.
      paused[i] = 1;
    }
  }
 
  current = first;
  previous = NULL;
  while(current != NULL) {
   if(current->paused) { // If paused, check whether to proceed.
    draw = (double)rand()/ (double)RAND_MAX;
    if(draw < p_r)
      release_polii(current); // Draw a number and ask whether this Pol II releases.
      paused[current->cell] = 0;
   }
   else {
    move_polii(current); // If not paused, move the Pol II.
 
    // If position > tu_size, terminate.
    if((current->position > tu_size && current->strand == '+') || (current->position < (uas_size - uas_offset) && current->strand == '-')) {
      draw = (double)rand()/ (double)RAND_MAX;
      if(draw < p_t) {
        terminate_polii(current, previous, &first);
        current = previous->next; // Have to adjust the pointer. Otherwise previous will point to a terminated polii.
        continue;
      }
    }
   }
   
   previous = current;
   current = current->next;
  }
 }

 in_silico_proseq(first);

 return(0);
}

int main(int argc, char *argv[]) {

 if(argc != 3) {
  printf("Usage: ./mcmc_polii $INITIATION[1/s] $RELEASE[1/s]\n");
  return(1);
 }

 p_i= (double)atof(argv[1]);
 p_r= (double)atof(argv[2]);
 //printf("Running MCMC: %f\t%f\n", p_i, p_r);

 return(mcmc_polii());
}


