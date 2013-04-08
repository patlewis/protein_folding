#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "proteins.h"
#include "list.h"

static const char *names[] = {"ALA", "ARG", "ASN", "ASP", "ASX", "CYS", "GLU", \
    "GLN", "GLX", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};


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
