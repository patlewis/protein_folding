#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>


#include "list.h"
#include "proteins.h"

#define MAX_PROTEINS 10

/* Function declarations */
void print_structure(vertex thing[41][41][41]);
void DFS(int i, int j, int k, int previous);
double calculate_energy(void);
int sort(const void *a, const void *b);
void print_current_structure(void);
static bool chain_is_finished(void);




/* Global variables */
    vertex skeleton[41][41][41];
    amino_acid_chain *chain;
    three_d_protein proteins[MAX_PROTEINS];
    int num_proteins;
    double min_energy;
    double max_energy;
    unsigned long total_structures;
    size_t protein_size;


/* This is where the fun begins */
int main(int argc, char **argv){
    chain = create_chain();
    protein_size = list_size(&chain->amino_acid_list);
    min_energy = 0;
    max_energy = 0;
    num_proteins = 0;
    total_structures = 0;
    int i,j,k;

    for(i=0; i < 41; i++){
        for (j = 0; j < 41; j++){
            for (k = 0; k < 41; k++){
                skeleton[i][j][k].x = i;
                skeleton[i][j][k].y = j;
                skeleton[i][j][k].z = k;
            }
        }
    }
    DFS(20, 20, 20, 0);
    for(i = 0; i < MAX_PROTEINS; i++){
//        printf("Printing...\n");
        printf("====================================================================================================================\n");
        print_protein_3(proteins[i], (int)protein_size);
        printf("Structure number %d\n", i+1);
    }
    free_chain(chain);    
    printf("\n\n\tTotal number of structures: %lu\n\n", total_structures);
  //  printf("Max number of structures we can keep track of: %lu\n", ULONG_MAX);
    return 0;
}


/*
 * Performs a depth-first search.  Starts with a specific
 * "vertex" and recursively does a depth-first search on 
 * each vertex around this one.  
 *
 * 'prev' doesn't do anything so I should probably take
 * it out, but right now I just want to finish this
 */
void DFS(int i, int j, int k, int prev){
    
    /* Get the next amino acid to be placed into the skeleton */
    struct list_elem *e = list_pop_front(&chain->amino_acid_list);
    amino_acid *ami_aci = list_entry(e, struct amino_acid, elem);
    /*Set the attributes for this part of the skeleton */
    skeleton[i][j][k].amino = ami_aci;
    skeleton[i][j][k].next = NULL;
    //see if stop condition reached
    if(chain_is_finished()){
        list_push_front(&chain->amino_acid_list,&(skeleton[i][j][k].amino)->elem);   //put it back on the stack/queue to be reused
        skeleton[i][j][k].amino = NULL;    //remove the amino acid that used to be here
        skeleton[i][j][k].next = NULL;  //it was already null so I don't remember why I put this in
        return;
    }
    /* Start decision making for the next part of the chain */
    if((i-1 >= 0) && (skeleton[i-1][j][k].amino == NULL)){//if "North" spot is available, go there
        skeleton[i][j][k].next = &(skeleton[i-1][j][k]);
        skeleton[i-1][j][k].prev = &(skeleton[i][j][k]);
        DFS(i-1, j, k, 1);
        skeleton[i][j][k].next = NULL;
    }
    if((j+1 <= 40) && (skeleton[i][j+1][k].amino == NULL)){//if "East" is available...
        skeleton[i][j][k].next = &(skeleton[i][j+1][k]);
        skeleton[i][j+1][k].prev = &(skeleton[i][j][k]);
        DFS(i, j+1, k, 2);
        skeleton[i][j][k].next = NULL;
    }
    if((i+1 <= 40) && (skeleton[i+1][j][k].amino == NULL)){//if "South"...
        skeleton[i][j][k].next = &(skeleton[i+1][j][k]);
        skeleton[i+1][j][k].prev = &(skeleton[i][j][k]);
        DFS(i+1, j, k, 3);
        skeleton[i][j][k].next = NULL;
    }
    if((j-1 >= 0) && (skeleton[i][j-1][k].amino == NULL)){//if "West"...
        skeleton[i][j][k].next = &(skeleton[i][j-1][k]);
        skeleton[i][j-1][k].prev = &(skeleton[i][j][k]);
        DFS(i, j-1, k, 4);
        skeleton[i][j][k].next = NULL;
    }
    if((k+1 <= 40) && (skeleton[i][j][k+1].amino == NULL)){//if "UP"...
        skeleton[i][j][k].next = &(skeleton[i][j][k+1]);
        skeleton[i][j][k+1].prev = &(skeleton[i][j][k]);
        DFS(i, j, k+1, 5);
        skeleton[i][j][k].next = NULL;
    }
    if((k-1 >= 0) && (skeleton[i][j][k-1].amino == NULL)){//if "DOWN"...
        skeleton[i][j][k].next = &(skeleton[i][j][k-1]);
        skeleton[i][j][k-1].prev = &(skeleton[i][j][k]);
        DFS(i, j, k-1, 6);
        skeleton[i][j][k].next = NULL;
    }
    //once this node has no more branches to make:
    list_push_front(&chain->amino_acid_list,&(skeleton[i][j][k].amino)->elem);   //put it back on the stack/queue to be reused
    skeleton[i][j][k].amino = NULL;    //remove the amino acid that used to be here
                                    //so that other DFS searches can use this vertex
    

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
    static const double EHH = -2.3;
    static const double EHP = -1;
    static const double EPP = 0;
    static const double SPECIAL = -3.5; //for ARG and ASP
    int i, x, y, z;
    x = 20;
    y = 20;
    z = 20;
    vertex current = skeleton[x][y][z];
    for(i = 0; i < protein_size; i++){
        x = current.x;
        y = current.y;
        z = current.z;
        vertex *nxt = current.next;
        vertex *prv = current.prev;
        char *cname = (current.amino)->name;
        /* Decision blocks for seeing if there's something we can caclulate the energy of */
        if(nxt == NULL){ //last block in the chain
            //If amino acid to the "south" and it is not directly "next to" in the original chain
            if(skeleton[x+1][y][z].amino != NULL && skeleton[x+1][y][z].amino != prv->amino){
                if(strcmp(cname, "ARG") == 0 && strcmp((skeleton[x+1][y][z].amino)->name, "ASP") == 0)
                {energy += SPECIAL;} //if special pair
                else if(strcmp(cname, "ASP") == 0 && strcmp((skeleton[x+1][y][z].amino)->name, "ARG") == 0)
                {energy += SPECIAL;} //if special pair
                else if((current.amino)->hydro && (skeleton[x+1][y][z].amino)->hydro){ energy += EHH;      }  //if both hydrophilic
                else if ((current.amino)->hydro || (skeleton[x+1][y][z].amino)->hydro){ energy += EHP;}  // if different types
                else{ energy += EPP;                                                               }  //if both hydrophobic
            }
            //If amino acid to the "north" and it is not directly "next to" in the original chain
            if(skeleton[x-1][y][z].amino != NULL && skeleton[x-1][y][z].amino != prv->amino){
                if(strcmp(cname, "ARG") == 0 && strcmp((skeleton[x-1][y][z].amino)->name, "ASP") == 0)
                {energy += SPECIAL;} //if special pair
                else if(strcmp(cname, "ASP") == 0 && strcmp((skeleton[x-1][y][z].amino)->name, "ARG") == 0)
                {energy += SPECIAL;} //if special pair
                else if((current.amino)->hydro && (skeleton[x-1][y][z].amino)->hydro){ energy += EHH;      }  //if both hydrophilic
                else if ((current.amino)->hydro || (skeleton[x-1][y][z].amino)->hydro){ energy += EHP;}  // if different types
                else{ energy += EPP;                                                               } //if both hydrophobic    
            }
            //If amino acid to the "east" and it is not directly "next to" in the original chain
            if(skeleton[x][y+1][z].amino != NULL && skeleton[x][y+1][z].amino != prv->amino){
                if(strcmp(cname, "ARG") == 0 && strcmp((skeleton[x][y+1][z].amino)->name, "ASP") == 0)
                {energy += SPECIAL;} //if special pair
                else if(strcmp(cname, "ASP") == 0 && strcmp((skeleton[x][y+1][z].amino)->name, "ARG") == 0)
                {energy += SPECIAL;} //if special pair
                else if((current.amino)->hydro && (skeleton[x][y+1][z].amino)->hydro){ energy += EHH;      }  //if both hydrophilic
                else if ((current.amino)->hydro || (skeleton[x][y+1][z].amino)->hydro){ energy += EHP;}  // if different types
                else{ energy += EPP;                                                               }  //if both hydrophobic          
            }
            //If amino acid to the "west" and it is not directly "next to" in the original chain
            if(skeleton[x][y-1][z].amino != NULL && skeleton[x][y-1][z].amino != prv->amino){
                if(strcmp(cname, "ARG") == 0 && strcmp((skeleton[x][y-1][z].amino)->name, "ASP") == 0)
                {energy += SPECIAL;} //if special pair
                else if(strcmp(cname, "ASP") == 0 && strcmp((skeleton[x][y-1][z].amino)->name, "ARG") == 0)
                {energy += SPECIAL;} //if special pair
                else if((current.amino)->hydro && (skeleton[x][y-1][z].amino)->hydro){ energy += EHH;      }  //if both hydrophilic
                else if ((current.amino)->hydro || (skeleton[x][y-1][z].amino)->hydro){ energy += EHP;}  // if different types
                else{ energy += EPP;                                                               }  //if both hydrophobic
            }
            //If amino acid to the "up" and it is not directly "next to" in the original chain
            if(skeleton[x][y][z+1].amino != NULL && skeleton[x][y][z+1].amino != prv->amino){
                if(strcmp(cname, "ARG") == 0 && strcmp((skeleton[x][y][z+1].amino)->name, "ASP") == 0)
                {energy += SPECIAL;} //if special pair
                else if(strcmp(cname, "ASP") == 0 && strcmp((skeleton[x][y][z+1].amino)->name, "ARG") == 0)
                {energy += SPECIAL;} //if special pair
                else if((current.amino)->hydro && (skeleton[x][y][z+1].amino)->hydro){ energy += EHH;      }  //if both hydrophilic
                else if ((current.amino)->hydro || (skeleton[x][y][z+1].amino)->hydro){ energy += EHP;}  // if different types
                else{ energy += EPP;                                                               }  //if both hydrophobic
            }
            //If amino acid to the "down" and it is not directly "next to" in the original chain
            if(skeleton[x][y][z-1].amino != NULL && skeleton[x][y][z-1].amino != prv->amino){
                if(strcmp(cname, "ARG") == 0 && strcmp((skeleton[x][y][z-1].amino)->name, "ASP") == 0)
                {energy += SPECIAL;} //if special pair
                else if(strcmp(cname, "ASP") == 0 && strcmp((skeleton[x][y][z-1].amino)->name, "ARG") == 0)
                {energy += SPECIAL;} //if special pair
                else if((current.amino)->hydro && (skeleton[x][y][z-1].amino)->hydro){ energy += EHH;      }  //if both hydrophilic
                else if ((current.amino)->hydro || (skeleton[x][y][z-1].amino)->hydro){ energy += EHP;}  // if different types
                else{ energy += EPP;                                                               } //if both hydrophobic    
            }
            continue; //if this is the last one in the chain
        }
        //If there is no "previous" link
        if(prv == NULL)
        {
            //If amino acid to the "south" and it is not directly "next to" in the original chain
            if(skeleton[x+1][y][z].amino != NULL && skeleton[x+1][y][z].amino != nxt->amino){
                if(strcmp(cname, "ARG") == 0 && strcmp((skeleton[x+1][y][z].amino)->name, "ASP") == 0)
                {energy += SPECIAL;} //if special pair
                else if(strcmp(cname, "ASP") == 0 && strcmp((skeleton[x+1][y][z].amino)->name, "ARG") == 0)
                {energy += SPECIAL;} //if special pair
                else if((current.amino)->hydro && (skeleton[x+1][y][z].amino)->hydro){ energy += EHH;      }  //if both hydrophilic
                else if ((current.amino)->hydro || (skeleton[x+1][y][z].amino)->hydro){ energy += EHP;}  // if different types
                else{ energy += EPP;                                                               }  //if both hydrophobic
            }
            //If amino acid to the "north" and it is not directly "next to" in the original chain
            if(skeleton[x-1][y][z].amino != NULL && skeleton[x-1][y][z].amino != nxt->amino){
                if(strcmp(cname, "ARG") == 0 && strcmp((skeleton[x-1][y][z].amino)->name, "ASP") == 0)
                {energy += SPECIAL;} //if special pair
                else if(strcmp(cname, "ASP") == 0 && strcmp((skeleton[x-1][y][z].amino)->name, "ARG") == 0)
                {energy += SPECIAL;} //if special pair
                else if((current.amino)->hydro && (skeleton[x-1][y][z].amino)->hydro){ energy += EHH;      }  //if both hydrophilic
                else if ((current.amino)->hydro || (skeleton[x-1][y][z].amino)->hydro){ energy += EHP;}  // if different types
                else{ energy += EPP;                                                               } //if both hydrophobic    
            }
            //If amino acid to the "east" and it is not directly "next to" in the original chain
            if(skeleton[x][y+1][z].amino != NULL && skeleton[x][y+1][z].amino != nxt->amino){
                if(strcmp(cname, "ARG") == 0 && strcmp((skeleton[x][y+1][z].amino)->name, "ASP") == 0)
                {energy += SPECIAL;} //if special pair
                else if(strcmp(cname, "ASP") == 0 && strcmp((skeleton[x][y+1][z].amino)->name, "ARG") == 0)
                {energy += SPECIAL;} //if special pair
                else if((current.amino)->hydro && (skeleton[x][y+1][z].amino)->hydro){ energy += EHH;      }  //if both hydrophilic
                else if ((current.amino)->hydro || (skeleton[x][y+1][z].amino)->hydro){ energy += EHP;}  // if different types
                else{ energy += EPP;                                                               }  //if both hydrophobic          
            }
            //If amino acid to the "west" and it is not directly "next to" in the original chain
            if(skeleton[x][y-1][z].amino != NULL && skeleton[x][y-1][z].amino != nxt->amino){
                if(strcmp(cname, "ARG") == 0 && strcmp((skeleton[x][y-1][z].amino)->name, "ASP") == 0)
                {energy += SPECIAL;} //if special pair
                else if(strcmp(cname, "ASP") == 0 && strcmp((skeleton[x][y-1][z].amino)->name, "ARG") == 0)
                {energy += SPECIAL;} //if special pair
                else if((current.amino)->hydro && (skeleton[x][y-1][z].amino)->hydro){ energy += EHH;      }  //if both hydrophilic
                else if ((current.amino)->hydro || (skeleton[x][y-1][z].amino)->hydro){ energy += EHP;}  // if different types
                else{ energy += EPP;                                                               }  //if both hydrophobic
            }
            //If amino acid to the "up" and it is not directly "next to" in the original chain
            if(skeleton[x][y][z+1].amino != NULL && skeleton[x][y][z+1].amino != nxt->amino){
                if(strcmp(cname, "ARG") == 0 && strcmp((skeleton[x][y][z+1].amino)->name, "ASP") == 0)
                {energy += SPECIAL;} //if special pair
                else if(strcmp(cname, "ASP") == 0 && strcmp((skeleton[x][y][z+1].amino)->name, "ARG") == 0)
                {energy += SPECIAL;} //if special pair
                else if((current.amino)->hydro && (skeleton[x][y][z+1].amino)->hydro){ energy += EHH;      }  //if both hydrophilic
                else if ((current.amino)->hydro || (skeleton[x][y][z+1].amino)->hydro){ energy += EHP;}  // if different types
                else{ energy += EPP;                                                               }  //if both hydrophobic
            }
            //If amino acid to the "down" and it is not directly "next to" in the original chain
            if(skeleton[x][y][z-1].amino != NULL && skeleton[x][y][z-1].amino != nxt->amino){
                if(strcmp(cname, "ARG") == 0 && strcmp((skeleton[x][y][z-1].amino)->name, "ASP") == 0)
                {energy += SPECIAL;} //if special pair
                else if(strcmp(cname, "ASP") == 0 && strcmp((skeleton[x][y][z-1].amino)->name, "ARG") == 0)
                {energy += SPECIAL;} //if special pair
                else if((current.amino)->hydro && (skeleton[x][y][z-1].amino)->hydro){ energy += EHH;      }  //if both hydrophilic
                else if ((current.amino)->hydro || (skeleton[x][y][z-1].amino)->hydro){ energy += EHP;}  // if different types
                else{ energy += EPP;                                                               } //if both hydrophobic    
            }
            current = *(current.next);
            continue;
        }
        //If there is a "previous" link we can check then we do this:
        //If amino acid to the "south" and it is not directly "next to" in the original chain
        if(skeleton[x+1][y][z].amino != NULL && \
            skeleton[x+1][y][z].amino != nxt->amino && \
            skeleton[x+1][y][z].amino != prv->amino){
                if(strcmp(cname, "ARG") == 0 && strcmp((skeleton[x+1][y][z].amino)->name, "ASP") == 0)
                {energy += SPECIAL;} //if special pair
                else if(strcmp(cname, "ASP") == 0 && strcmp((skeleton[x+1][y][z].amino)->name, "ARG") == 0)
                {energy += SPECIAL;} //if special pair
                else if((current.amino)->hydro && (skeleton[x+1][y][z].amino)->hydro){ energy += EHH;      }  //if both hydrophilic
                else if ((current.amino)->hydro || (skeleton[x+1][y][z].amino)->hydro){ energy += EHP;}  // if different types
                else{ energy += EPP;                                                               }  //if both hydrophobic 
        }
        //If amino acid to the "north" and it is not directly "next to" in the original chain
        if(skeleton[x-1][y][z].amino != NULL && \
            skeleton[x-1][y][z].amino != nxt->amino && \
            skeleton[x-1][y][z].amino != prv->amino){
                if(strcmp(cname, "ARG") == 0 && strcmp((skeleton[x-1][y][z].amino)->name, "ASP") == 0)
                {energy += SPECIAL;} //if special pair
                else if(strcmp(cname, "ASP") == 0 && strcmp((skeleton[x-1][y][z].amino)->name, "ARG") == 0)
                {energy += SPECIAL;} //if special pair
                else if((current.amino)->hydro && (skeleton[x-1][y][z].amino)->hydro){ energy += EHH;      }  //if both hydrophilic
                else if ((current.amino)->hydro || (skeleton[x-1][y][z].amino)->hydro){ energy += EHP;}  // if different types
                else{ energy += EPP;                                                               } //if both hydrophobic
        }
        //If amino acid to the "east" and it is not directly "next to" in the original chain
        if(skeleton[x][y+1][z].amino != NULL && \
            skeleton[x][y+1][z].amino != nxt->amino && \
            skeleton[x][y+1][z].amino != prv->amino){
                if(strcmp(cname, "ARG") == 0 && strcmp((skeleton[x][y+1][z].amino)->name, "ASP") == 0)
                {energy += SPECIAL;} //if special pair
                else if(strcmp(cname, "ASP") == 0 && strcmp((skeleton[x][y+1][z].amino)->name, "ARG") == 0)
                {energy += SPECIAL;} //if special pair
                else if((current.amino)->hydro && (skeleton[x][y+1][z].amino)->hydro){ energy += EHH;      }  //if both hydrophilic
                else if ((current.amino)->hydro || (skeleton[x][y+1][z].amino)->hydro){ energy += EHP;}  // if different types
                else{ energy += EPP;                                                               }  //if both hydrophobic
        }
        //If amino acid to the "west" and it is not directly "next to" in the original chain
        if(skeleton[x][y-1][z].amino != NULL && \
            skeleton[x][y-1][z].amino != nxt->amino && \
            skeleton[x][y-1][z].amino != prv->amino){
                if(strcmp(cname, "ARG") == 0 && strcmp((skeleton[x][y-1][z].amino)->name, "ASP") == 0)
                {energy += SPECIAL;} //if special pair
                else if(strcmp(cname, "ASP") == 0 && strcmp((skeleton[x][y-1][z].amino)->name, "ARG") == 0)
                {energy += SPECIAL;} //if special pair
                else if((current.amino)->hydro && (skeleton[x][y-1][z].amino)->hydro){ energy += EHH;      }  //if both hydrophilic
                else if ((current.amino)->hydro || (skeleton[x][y-1][z].amino)->hydro){ energy += EHP;}  // if different types
                else{ energy += EPP;                                                               }  //if both hydrophobic
        }
        //If amino acid to the "up" and it is not directly "next to" in the original chain
        if(skeleton[x][y][z+1].amino != NULL && \
            skeleton[x][y][z+1].amino != nxt->amino && \
            skeleton[x][y][z+1].amino != prv->amino){
                if(strcmp(cname, "ARG") == 0 && strcmp((skeleton[x][y][z+1].amino)->name, "ASP") == 0)
                {energy += SPECIAL;} //if special pair
                else if(strcmp(cname, "ASP") == 0 && strcmp((skeleton[x][y][z+1].amino)->name, "ARG") == 0)
                {energy += SPECIAL;} //if special pair
                else if((current.amino)->hydro && (skeleton[x][y][z+1].amino)->hydro){ energy += EHH;      }  //if both hydrophilic
                else if ((current.amino)->hydro || (skeleton[x][y][z+1].amino)->hydro){ energy += EHP;}  // if different types
                else{ energy += EPP;                                                               }  //if both hydrophobic
        }
        //If amino acid to the "down" and it is not directly "next to" in the original chain
        if(skeleton[x][y][z-1].amino != NULL && \
            skeleton[x][y][z-1].amino != nxt->amino && \
            skeleton[x][y][z-1].amino != prv->amino){
                if(strcmp(cname, "ARG") == 0 && strcmp((skeleton[x][y][z-1].amino)->name, "ASP") == 0)
                {energy += SPECIAL;} //if special pair
                else if(strcmp(cname, "ASP") == 0 && strcmp((skeleton[x][y][z-1].amino)->name, "ARG") == 0)
                {energy += SPECIAL;} //if special pair
                else if((current.amino)->hydro && (skeleton[x][y][z-1].amino)->hydro){ energy += EHH;      }  //if both hydrophilic
                else if ((current.amino)->hydro || (skeleton[x][y][z-1].amino)->hydro){ energy += EHP;}  // if different types
                else{ energy += EPP;                                                               }  //if both hydrophobic
        }
        current = *(current.next);
    } //end for loop
    return energy/2;
}


/* Sorting algorithm based in insertion sort.
 * Insertion sort runs in O(n) in the best case,
 * and all of these cases are close to the best
 * case since we're really only finding a spot
 * for one element, so it's pretty cool.
 */
int sort(const void *p1, const void *p2){
    three_d_protein *pro1 = (three_d_protein *)p1;
    three_d_protein *pro2 = (three_d_protein *)p2;
    if(pro1->energy < pro2->energy) return -1;
    else if (pro2->energy < pro1->energy) return 1;
    else return 0;
}


/*
 * For debugging purposes only.  Don't use it.
 *
 * Seriously.  Don't.
 */
void print_current_structure(){
    int i, j, k;
    printf("========================================================================================================================\n");
    for(i = 0; i < 27; i++){
        for(j = 0; j < 27; j++){
            for(k = 0; k < 27; k++){
                if(skeleton[i][j][k].amino == NULL) printf("    ");
                else printf(" %s", (skeleton[i][j][k].amino)->name);
            }
        }
        printf("\n");
    }
    printf("Structure number %lu\n", ++total_structures);
    printf("========================================================================================================================\n");
}


/*
 * Checks to see of there are no more amino acids left to 
 * put in the protein.  If there are, it immediately returns
 * false.  Otherwise, it calculates the energy of this protein
 * structure and stores informations as necessary.
 */
static bool chain_is_finished(){
    if(list_empty(&chain->amino_acid_list)){//if all protein is used
        total_structures++;
        double energy = calculate_energy();
        int count;
        //if still populating the list
        if(num_proteins < MAX_PROTEINS){
            proteins[num_proteins++] = three_d_protein_create(skeleton, energy);
//            print_protein(proteins[num_proteins-1]);
            return true;
        }
        //if list is full, see if this should go in
        for(count = 0; count < MAX_PROTEINS; count++){
            if (energy < proteins[count].energy){   //if it belongs in array
                proteins[MAX_PROTEINS-1] = three_d_protein_create(skeleton, energy);
                qsort(proteins, MAX_PROTEINS, sizeof(two_d_protein), sort);
                return true;
            }
        }
        //if we got here, then it wasn't added to the array, but the chain 
        //was finished
        return true;
    }
    else return false;
}
