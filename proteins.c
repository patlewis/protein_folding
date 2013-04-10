#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "proteins.h"
#include "list.h"

#define SIZE 41

static const char *names[] = {"ALA", "ARG", "ASN", "ASP", "ASX", "CYS", "GLU", \
    "GLN", "GLX", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};

static const char *philic[] = {"ASN", "ASP", "ASX", "GLU", "GLN", "GLX"};

amino_acid_chain * create_chain(void){
    amino_acid_chain *chain = (amino_acid_chain *)malloc(sizeof(amino_acid_chain));
    char line[10000];
    printf("Reading input...");
    if(fgets(line, 9999, stdin) == NULL){
        printf("Some kind of error with the input.\n");
        exit(0);
    }
    /* Create list for temporary storage of amino acids */
    list_init(&chain->amino_acid_list);
    /* Find size of chain and do basic error checking to make sure it's valid */
    int i;
    int size = 1;
    char *token = strtok(line, "- \n");
    bool matchfound;
    while(token != NULL){
        //error checking
        matchfound = false;
        for(i = 0; i < 22; i++){
            if(strcmp(token, names[i]) == 0){   //if match found, good
                matchfound = true;
                list_push_back(&chain->amino_acid_list, &((new_amino_acid(token))->elem));
                token = strtok(NULL, "- \n");
                break;
            }
        }
        if(!matchfound){
            printf("\nError: could not rectify input \"%s\" at position %d.\n", token, size);
            exit(0);
        }
        size++;
    }
    printf("Valid amino acid chain.\n");
    /* Now that chain is created empty, fill it */
    return chain;
}

amino_acid *new_amino_acid(char *label){                                                                                                  
    amino_acid *acid = (amino_acid *)malloc(sizeof(amino_acid));
    strncpy(acid->name, label, (size_t)3);
    int i;
    bool phil = false;
    for(i = 0; i < 6; i++){
        if(strcmp(label, *(philic+i)) == 0){
            acid->hydro = 1;
            phil = true;
        }
    }
    if(!phil) acid->hydro = 0;
    return acid;
}

void free_chain(amino_acid_chain *chain){
    struct list_elem *e;
    struct amino_acid *acid;
    while(!list_empty(&chain->amino_acid_list)){
        e = list_pop_front(&chain->amino_acid_list);
        acid = list_entry(e, struct amino_acid, elem); 
        free(acid);
    }
    free(chain);
}

two_d_protein two_d_protein_create(struct vertex array[SIZE][SIZE], double energy){
    struct two_d_protein *pro = (two_d_protein *)malloc(sizeof(two_d_protein));
    amino_acid *ac;
    int i, j;
    for(i = 0; i < 41; i++){
        for(j = 0; j < 41; j++){
            ac = array[i][j].amino;
            if(ac != NULL){
                pro->structure[i][j] = *(new_amino_acid(ac->name));
            }
        }
    }
    pro->energy = energy;
    return *pro;
}


void two_d_protein_free(two_d_protein *pro){
    if(pro != NULL) free(pro);
}


void print_protein(two_d_protein pro){
    int i,j;
    for(i = 0; i < 10; i++){
        for(j = 0; j < 10; j++){
            if(strcmp((pro.structure[i][j]).name, "") != 0){
                printf("%s ", (pro.structure[i][j]).name);
            }
            else printf("    ");
        }
        printf("\n");
    }
    printf("Energy: %f\n\n\n", pro.energy);
}



