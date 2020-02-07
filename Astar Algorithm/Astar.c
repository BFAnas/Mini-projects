/*
 * Astar.c
 * 
 * Copyright 2020 Anas <anas@anas-MACH-WX9>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265359
#define R 6.3781e6 //earth radius in meters

unsigned long nnodes = 23895681UL;   // 3472620, 23895681

typedef char Queue;
enum whichQueue {NONE, OPEN, CLOSED};

/* Data structure for the status of the node */
typedef struct {
double g, h;
unsigned long parent;
Queue whq;
} AStarStatus;

/* Data structure for the nodes */
typedef struct {
unsigned long id; // Node identification
char *name;
double lat, lon; // Node position
unsigned short nsucc; // Number of node successors; i. e. length of successors
unsigned long *successors; // List of successors
} node;



/* linked list to store the open list nodes in a sorted way */
struct linked_list{
	unsigned long index; // internal id of the node
	double f;
	struct linked_list * next;
	};

typedef struct linked_list open_list;

node * nodes;
open_list *Open;
AStarStatus *Status;


	
/* function to insert a new_node in a list in a sorted way */
void sortedInsert(open_list ** head_ref, open_list * new_node) 
{ 
	open_list * current; 
	/* Special case for the head end */
	if (*head_ref == NULL || (*head_ref)->f >= new_node->f) { 
		new_node->next = *head_ref; 
		*head_ref = new_node; 
	} 
	else{ 
		/* Locate the node before the point of insertion */
		current = *head_ref; 
		while (current->next!=NULL && 
					 current->next->f < new_node->f) 
		{ 
				current = current->next; 
		} 
		new_node->next = current->next; 
		current->next = new_node; 
	} 
}

/* Given a reference (pointer to pointer) to the head of a list 
   and a key, deletes the first occurrence of key in linked list */
void deleteNode(open_list **head_ref, unsigned long key) 
{ 
	// Store head node 
	open_list* temp = *head_ref, *prev; 

	// If head node itself holds the key to be deleted 
	if (temp != NULL && temp->index == key) { 
			*head_ref = temp->next;   // Changed head 
			free(temp);               // free old head 
			return; 
	} 

	// Search for the key to be deleted, keep track of the 
	// previous node as we need to change 'prev->next' 
	while (temp != NULL && temp->index != key) { 
			prev = temp; 
			temp = temp->next; 
	} 

	// If key was not present in linked list 
	if (temp == NULL) return; 

	// Unlink the node from linked list 
	prev->next = temp->next; 

	free(temp);  // Free memory 
} 

double Haversine(node *init, node *prev)
{	
	// Haversine formula 
	double phi1= init->lat/180*PI;
	double phi2= prev->lat/180*PI;
	double lambda1=init->lon/180*PI;
	double lambda2=prev->lon/180*PI;
	double dif_lat= fabs(phi1-phi2);
	double dif_lon= fabs(lambda1-lambda2);

	double a = sin(dif_lat/2)*sin(dif_lat/2)+cos(phi1)*cos(phi2)*sin(dif_lon/2)*sin(dif_lon/2);
	double c= 2*atan2(sqrt(a),sqrt(1-a));
	
	return c*R;
}

/* Function for binary search */
long bi_search(node *Vector, unsigned long Vector_length, unsigned long target) {
	unsigned long imin, imax, imid;
	imin = 0; imax = Vector_length-1;
	while (imax >= imin) {
		imid = (imax + imin)*0.5;
		if (Vector[imid].id == target) {
			return imid;
		}
		else if (Vector[imid].id < target) {
			imin = imid + 1;
		}
		else {
			imax = imid - 1;
		}
	}
	return -1;
}


void reading()
{
	int i;
	static unsigned long * allsuccessors;
	unsigned long total_nsucc=0UL;

	FILE * fin;

	fin=fopen("binary.bin","rb");
	
//GLOBAL DATA HEADER
	if ((fread(&nnodes,sizeof(unsigned long),1,fin)+
	fread(&total_nsucc,sizeof(unsigned long),1,fin))!=2) 
	printf("when reading the header of the binary data file\n");
	// GETTING MEMORY FOR ALL THE DATA
	if ((nodes= (node *) malloc(sizeof(node)*nnodes))==NULL) 
	printf("when allocating memory for the nodes vector\n");
	if ((allsuccessors= (unsigned long*) 
	malloc(sizeof(unsigned long)*total_nsucc))==NULL) 
	printf("when allocating memory for the edges vector\n");
//READING ALL DATA FROM FILE
	if (fread(nodes,sizeof(node),nnodes,fin)!=nnodes) 
	printf("when reading nodes from the binary data file\n");
	if (fread(allsuccessors,sizeof(unsigned long),total_nsucc,fin)!=total_nsucc) 
	printf("when reading successors from the binary data file\n");

	fclose(fin);
//SETTING POINTERS TO SUCCESSORS 
	for (i=0;i<nnodes; i++) if (nodes[i].nsucc){
		nodes[i].successors=allsuccessors;
		allsuccessors+=nodes[i].nsucc;
	}
}


int AStar (unsigned long start, unsigned long goal)
{
	node *current , *successor;
	
	int i, nelements=0;
	open_list *succ; 
	double successor_current_cost;
	
	
	unsigned long start_index = bi_search(nodes, nnodes, start);
	unsigned long goal_index = bi_search(nodes, nnodes, goal);
	
	/* Allocating memory for the status of the nodes */
	Status = (AStarStatus *) calloc(nnodes, sizeof(AStarStatus));
	
	Open = (open_list *) calloc(1, sizeof(open_list));
	succ = (open_list *) malloc(1* sizeof(open_list));
	
	/* Put node_start in the OPEN list */
	Open->index = start_index;
	
	/* The number of elements in the open list now is nelements = 1 */
	nelements += 1;

	/* Change the Queue status for node start to open */
	Status->whq = 1;
	
	/* Set g and h for the start node */
	Status->g = 0;
	Status->h =  Haversine(&nodes[start_index], &nodes[goal_index]);
	
	/* f(start) = h(start) */
	Open->f = Status->h;
	//printf("%f \n", Open->f);
	
	while (nelements>0){
		current = &nodes[Open->index]; //takes first node in the open list: the one with lowest f distance 
		unsigned long current_index = bi_search(nodes, nnodes, current->id); //find its position in the nodes array
		
		if ( current->id == goal) return 0; //if node_current is node_goal we have found the solution; break
		if (current->nsucc > 0) {
			/*Generate each state node_successor that come after node_current*/
			for ( i = 0; i < (current->nsucc); i++){ 
				unsigned long successor_index = bi_search(nodes, nnodes, *(current->successors+i));
				successor = &nodes[successor_index];
				successor_current_cost = Status[current_index].g + Haversine(successor, current); 
				if ( Status[successor_index].whq == 1){
						if ( Status[successor_index].g <= successor_current_cost) continue; 
					}

				else if ( Status[successor_index].whq == 2){
						if ( Status[successor_index].g <= successor_current_cost) continue;
						Status[successor_index].whq = 1;
						succ = (open_list *) malloc(1* sizeof(open_list));
						succ->index = successor_index;
						succ->f = successor_current_cost + Haversine(successor, &nodes[goal_index]);
						nelements += 1; 
						sortedInsert(&Open, succ);
				}

				else{
					Status[successor_index].whq = 1;
					succ = (open_list *) malloc(1* sizeof(open_list));
					succ->index = successor_index;
					
					succ->f = successor_current_cost + Haversine(successor, &nodes[goal_index]); 
					
					sortedInsert(&Open, succ);
					nelements += 1;
					Status[successor_index].h = Haversine(successor, &nodes[goal_index]); 
				}
				
				Status[successor_index].g = successor_current_cost;	
				Status[successor_index].parent = current_index;
			}
		}

	Status[current_index].whq = 2;
	deleteNode (&Open, current_index); 
	nelements -= 1;
	}

	return 1;
}


int main ()
{	
	unsigned long start = 240949599; //771979683  240949599
	unsigned long goal = 195977239; //429854583  195977239
	time_t ttime, ttime1; 
	ttime = clock();
	int counter=0;
	
	reading();
	printf("Total time spent reading : %fs \n",((clock()-(double)ttime))/CLOCKS_PER_SEC);
	
	//unsigned long start_index = bi_search(nodes, nnodes, start);
	unsigned long goal_index = bi_search(nodes, nnodes, goal);	
	
	/* Allocating memory for the status of the nodes */
	Status = (AStarStatus *) malloc(nnodes * sizeof(AStarStatus));
	Open = (open_list *) malloc(1* sizeof(open_list));
	
	ttime1 = clock();
	AStar(start, goal);
	printf("Total time spent in A* : %fs \n",((clock()-(double)ttime1))/CLOCKS_PER_SEC);
	
	AStarStatus* curr;
	long n;
		
	FILE * output;
	output=fopen("output_a_star.txt","w");	

	printf("Iteration |  Node  Id  | Distance  |    Lat    |   Long \n");
	
	n=0;
	curr = &Status[goal_index];
	unsigned long id = goal_index;
	while (curr->parent != 0){
		printf("     %4.ld | %10.lu | %6.2f | %2.6f | %2.6f \n", n, nodes[id].id, curr->g, nodes[id].lat, nodes[id].lon);
		fprintf(output,"     %4.ld | %10.lu | %6.2f | %2.6f | %2.6f \n", n, nodes[id].id, curr->g, nodes[id].lat, nodes[id].lon);
		id = curr->parent;
		curr = &Status[curr->parent];
		n++;
	}
	
	printf("     %4.ld | %10.lu | %6.2f | %2.6f | %2.6f \n", n, nodes[id].id, curr->g, nodes[id].lat, nodes[id].lon);
	fprintf(output,"     %4.ld | %10.lu | %6.2f | %2.6f | %2.6f \n", n, nodes[id].id, curr->g, nodes[id].lat, nodes[id].lon);
	
	for (int i=0; i<nnodes; i++){
		if (Status[i].whq != 0) counter++;
	} 
	
	printf("Total number of visited nodes : %d \n", counter);
	printf("Total time spent : %fs ",((clock()-(double)ttime))/CLOCKS_PER_SEC);
	
	return 0;
}
