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


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//23895681  3472620
#define FILENAME "spain.csv"
unsigned long nnodes = 23895681UL;
int nfields = 2009;

/* Data structure for the nodes */
typedef struct{
unsigned long id; // Node identification
char *name;
double lat, lon; // Node position
unsigned short nsucc; // Number of node successors; i. e. length of successors
unsigned long *successors; // List of successors
} node;

node *nodes;

/* Function for binary search */
unsigned long imin, imax, imid;
long bi_search(node *Vector, unsigned long Vector_length, unsigned long target){
	imin = 0; imax = Vector_length-1;
	while (imax >= imin){
		imid = (imax + imin)*0.5;
		if (Vector[imid].id == target){
			return imid;
		}
		else if (Vector[imid].id < target){
			imin = imid + 1;
		}
		else {
			imax = imid - 1;
		}
	}
	return -1;
}

/* Function for writing the graph in a binary file*/
void writing(node *nodes){
	FILE *fin;
	int i;
	unsigned long total_nsucc=0UL;
	for(i=0;i<nnodes; i++) total_nsucc+=nodes[i].nsucc;
	if((fin=fopen("binary.bin","wb"))==NULL) 
	printf("the output binary data file cannot be opened\n");
	/* Global data---header */
	 if((fwrite(&nnodes,sizeof(unsigned long),1,fin)+
	 	fwrite(&total_nsucc,sizeof(unsigned long),1,fin))!=2) 
	 	printf("when initializing the output binary data file\n");
	 /* Writting all nodes */
	 if( fwrite(nodes,sizeof(node),nnodes,fin)!=nnodes) 
	 printf ("when writing nodes to the output binary data file\n");
	//WRITING SUCCESORS IN BLOCKS
	for (i=0;i<nnodes; i++) 
		if(nodes[i].nsucc){
			if(fwrite(nodes[i].successors,sizeof(unsigned long),
			nodes[i].nsucc,fin)!=nodes[i].nsucc)
			printf("when writing edges to the output binary data file\n");	
		}
fclose(fin);
}


int main(int argc, char **argv){
	time_t ttime; 
	ttime = clock();
	/* Allocating memory for the vector of nodes */
	nodes = (node*) malloc(nnodes*sizeof(node));
	/* Allocating  memory for nodes successors and initializing
	 * node ids and number of successors to zero */
	for(int i=0;i<nnodes;i++){
		nodes[i].id = 0;
		nodes[i].nsucc = 0;
		nodes[i].successors = malloc(sizeof(unsigned long)*1);
	}	
	/* Open the file for reading */
  char *line_buf = NULL;
  size_t line_buf_size = 0;
  int line_count = 0;
  ssize_t line_size;
  FILE *fp = fopen(FILENAME, "r");
  if (!fp)
  {
    fprintf(stderr, "Error opening file '%s'\n", FILENAME);
    return EXIT_FAILURE;
  }
  /* Get the first line of the file. */
  line_size = getline(&line_buf, &line_buf_size, fp);
  /* nodenum variable for given the order number of the node */
  unsigned long nodenum = 0;
  while (line_size >= 0) {
    /* Go line by line and get nodes lines*/
    if (line_buf[0] == 'n') {
			/* Go through the line field by field and get the node info*/
			int field = 0;
			char *token;
			while ((token = strsep(&line_buf, "|")) != NULL) {
				if (field == 1) {
					/*We get the node id from the first field*/
					nodes[nodenum].id = strtoul(token, NULL, 10);
				}
				else if (field == 2) {
					/*From the second field we get the node name*/
					nodes[nodenum].name = token;
				}
				else if (field == 9) {
					/*From the 9th field we get the node latitude*/
					nodes[nodenum].lat = atof(token);
				}
				else if (field == 10) {
					/*From the 10th field we get the node longitude*/
					nodes[nodenum].lon = atof(token);
				}
				field += 1;
			}
			nodenum += 1;
		}		
		
		
		/* Go line by line and get number of neighbours*/
    else if (line_buf[0] == 'w')
    {
			/* Go through the line field by field and count the neighbours*/
			int field = 0;
			char *token;
			char *oneway;
			unsigned long token0;
			long index; 
			long index0;
			
			while ((token = strsep(&line_buf, "|")) != NULL) {
				if (field == 7) {
					oneway = token;
				}
				else if (field > 8) {
					index = bi_search(nodes, nnodes, strtoul(token, NULL, 10));
					if (field == 9) index0 = index;

					else if (index != -1 && index0 != -1) {
						nodes[index0].nsucc += 1;
						nodes[index0].successors = realloc(nodes[index0].successors,
						sizeof(unsigned long)*nodes[index0].nsucc);
						nodes[index0].successors[nodes[index0].nsucc-1] = strtoul(token, NULL, 10);
					
						if ((strcmp(oneway,"oneway") != 0) && (field > 9)) {
							nodes[index].nsucc += 1;
							nodes[index].successors = realloc(nodes[index].successors,
							sizeof(unsigned long)*nodes[index].nsucc);
							nodes[index].successors[nodes[index].nsucc-1] = token0;
						}
					}
					token0 = strtoul(token, NULL, 10);
					index0 = index;					
				}
				field += 1;
			}
		}
		
	/* Increment the line count */
    line_count++;

	/* Get the next line */
    line_size = getline(&line_buf, &line_buf_size, fp);
	}
	
	printf("Total time spent before writing : %fs ",((clock()-(double)ttime))/CLOCKS_PER_SEC);
	
	writing(nodes);
	free(nodes); 
	fclose(fp);
	
	printf("Total time spent : %fs ",((clock()-(double)ttime))/CLOCKS_PER_SEC);
	return 0;
}
