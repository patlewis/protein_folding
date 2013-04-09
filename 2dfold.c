#include <stdio.h>
#include <stdlib.h>
#include "list.h"
#include "proteins.h"


/* Global variables */
    vertex skeleton[41][41];
    amino_acid_chain *chain;
    two_d_protein *proteins;//used as an array
    int num_proteins;
    double min_energy;
    double max_energy;

/* Function declarations */
void print_structure(vertex thing[41][41]);
void DFS(int i, int j);
double calculate_energy(void);
void sort(void);

/* This is where the fun begins */
int main(int argc, char **argv){
    proteins = (two_d_protein *)malloc(100 * sizeof(two_d_protein));
    chain = create_chain();
    min_energy = 0;
    max_energy = 0;
    num_proteins = 0;

    int i,j;
    for(i=0; i < 41; i++){
        for (j = 0; j < 41; j++){
            skeleton[i][j].x = i;
            skeleton[i][j].y = j;
        }
    }
    DFS(20, 20);
    print_structure(skeleton);
    free_chain(chain);    
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
    if(list_empty(&chain->amino_acid_list)){//if all protein is used
        double en = calculate_energy();
        int i;
        two_d_protein pro;
        if(num_proteins < 100){
            pro = two_d_protein_create(skeleton, en);
            proteins[num_proteins] = pro;
            num_proteins++;
            sort();
            return;
        }
        for(i = 0; i < num_proteins; i++){
            if (en < proteins[i].energy){   //if it belongs in arraya
                two_d_protein_free(&(proteins[99]));
                pro = two_d_protein_create(skeleton, en);
                proteins[99] = pro;
                sort();
                return;
            }
        }
        //if we got here, then it wasn't added to the array
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

/*
 * Prints a 2d representation of the structure.
 * What does it sound like it does?
 */
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


/*
 * Calculates the energy of the completed structure.  Counts everything
 * twice because I didn't have enough time to think of a good solution.
 * It divides by 2 at the end to get to the right answer.  If you want 
 * to change the values of the interactions, go ahead and do so.  I just
 * used the energies that the other people used in their paper.
 */
double calculate_energy(void){
    double energy = 0;
    const double EHH = -2.3;
    const double EHP = -1;
    const double EPP = 0;
    int i, x, y;
    x = 20;
    y = 20;
    vertex current = skeleton[x][y];
    for(i = 0; i < 20; i++){
        x = current.x;
        y = current.y;
        /* Decision block for seeing if there's something we can caclulate the energy of */
        if(skeleton[x+1][y].amino != NULL && skeleton[x+1][y].amino != (current.next)->amino){
            if((current.amino)->hydro && (skeleton[x+1][y].amino)->hydro){ energy += EHH; continue; }
            else if ((current.amino)->hydro || (skeleton[x+1][y].amino)->hydro){ energy += EHP; continue; }
            else { energy += EPP; continue; }
        }
        if(skeleton[x-1][y].amino != NULL && skeleton[x-1][y].amino != (current.next)->amino){
            if((current.amino)->hydro && (skeleton[x-1][y].amino)->hydro){ energy += EHH; continue; }
            else if ((current.amino)->hydro || (skeleton[x-1][y].amino)->hydro){ energy += EHP; continue; }
            else { energy += EPP; continue; }
        }
        if(skeleton[x][y+1].amino != NULL && skeleton[x][y+1].amino != (current.next)->amino){
            if((current.amino)->hydro && (skeleton[x][y+1].amino)->hydro){ energy += EHH; continue; }
            else if ((current.amino)->hydro || (skeleton[x][y+1].amino)->hydro){ energy += EHP; continue; }
            else { energy += EPP; continue; }
        }
        if(skeleton[x][y-1].amino != NULL && skeleton[x][y-1].amino != (current.next)->amino){
            if((current.amino)->hydro && (skeleton[x][y-1].amino)->hydro){ energy += EHH; continue; }
            else if ((current.amino)->hydro || (skeleton[x][y-1].amino)->hydro){ energy += EHP; continue; }
            else { energy += EPP; continue; }
        }
    current = *(current.next);
    }
    //end for loop
    return energy/2;
}


/* Sorting algorithm based in insertion sort.
 * Insertion sort runs in O(n) in the best case,
 * and all of these cases are close to the best
 * case since we're really only finding a spot
 * for one element, so it's pretty cool.
 */
void sort(void){
    int i,j;
    for(i = 0; i < num_proteins; i++){
        double v = proteins[i].energy;
        for(j = i-1; j >= 0; j--){
            if(proteins[j].energy <= v) break;
            proteins[j+1] = proteins[j];           
        }
        proteins[j+1].energy = v;
    }
}




