#include <stdio.h>
#include <stdlib.h>
#include "list.h"
#include "proteins.h"


/* Global variables */
    vertex skeleton[41][41];
    amino_acid_chain *chain;

/* Function declarations */
void print_structure(vertex thing[41][41]);
void DFS(int i, int j);

/* This is where the fun begins */
int main(int argc, char **argv){
    chain = create_chain();

    int i,j;
    for(i=0; i < 41; i++){
        for (j = 0; j < 41; j++){
            skeleton[i][j].x = i;
            skeleton[i][j].y = j;
        }
    }
    DFS(20, 20);
    print_structure(skeleton);
    
    return 0;
}


/*
 * Performs a depth-first search.  Starts with a specific
 * "vertex" and recursively does a depth-first search on 
 * each vertex around this one.  
 */
void DFS(int i, int j){
    //TODO: put amino acids in vertices
    //
    if(list_empty(&chain->amino_acid_list)){
            //do stuff to end it like:
            //take a snapshot of the system
            //calculate its energy
            //if need be, save it in another list
            return;
    }
    /* Get the next amino acid to be placed into the skeleton */
    struct list_elem *e = list_pop_front(&chain->amino_acid_list);
    amino_acid *ami_aci = list_entry(e, struct amino_acid, elem);
    /*Set the attributes for this part of the skeleton */
    skeleton[i][j].amino = ami_aci;
    skeleton[i][j].next = NULL;
    /* Start decision making for the next part of the chain */
    if((i-1 >= 0) && (skeleton[i-1][j].amino == NULL)){
        skeleton[i][j].next = &(skeleton[i-1][j]);
        skeleton[i-1][j].prev = &(skeleton[i][j]);
        DFS(i-1, j);
        list_push_front(&chain->amino_acid_list,&(skeleton[i][j].amino)->elem);   //put it back on the stack/queue to be reused
        skeleton[i][j].amino = NULL;    //remove the amino acid that used to be here
        //so that other DFS searches can use this vertex
    }
    if((j+1 <= 40) && (skeleton[i][j+1].amino == NULL)){
        skeleton[i][j].next = &(skeleton[i][j+1]);
        skeleton[i][j+1].prev = &(skeleton[i][j]);
        DFS(i, j+1);
        list_push_front(&chain->amino_acid_list,&(skeleton[i][j].amino)->elem);   //put it back on the stack/queue to be reused
        skeleton[i][j].amino = NULL;    //remove the amino acid that used to be here
        //so that other DFS searches can use this vertex
    }
    if((i+1 <= 40) && (skeleton[i+1][j].amino == NULL)){
        skeleton[i][j].next = &(skeleton[i+1][j]);
        skeleton[i+1][j].prev = &(skeleton[i][j]);
        DFS(i+1, j);
        list_push_front(&chain->amino_acid_list,&(skeleton[i][j].amino)->elem);   //put it back on the stack/queue to be reused
        skeleton[i][j].amino = NULL;    //remove the amino acid that used to be here
        //so that other DFS searches can use this vertex
    }
    if((j-1 >= 0) && (skeleton[i][j-1].amino == NULL)){
        skeleton[i][j].next = &(skeleton[i][j-1]);
        skeleton[i][j-1].prev = &(skeleton[i][j]);
        DFS(i, j-1);
        list_push_front(&chain->amino_acid_list,&(skeleton[i][j].amino)->elem);   //put it back on the stack/queue to be reused
        skeleton[i][j].amino = NULL;    //remove the amino acid that used to be here
        //so that other DFS searches can use this vertex
    }
    

}

void print_structure(vertex thing[41][41]){
    int i, j;
    for(i = 0; i < 8; i++){
        for(j = 0; j < 8; j++){
            printf("(%d, %d)\t", thing[i][j].x, thing[i][j].y);
        }
        printf("\n");
    }
    printf("\n\n Done.\n");
}











