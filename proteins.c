#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "proteins.h"

static const char *names[] = {"ALA", "ARG", "ASN", "ASP", "ASX", "CYS", "GLU", "GLN", "GLX", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};

amino_acid * create_amino_acid(char* label){
    amino_acid *acid = malloc(sizeof(amino_acid));
    strcpy(acid->namelabel);
    return acid;
}

amino_acid_chain * create_chain(char *list){
    /* Find size of chain and do basic error checking to make sure it's valid */
    char str1[4000]; 
    char str2[4000]; 
    strcpy(str1, list);
    strcpy(str2, list);
    str1[3999] = '\0';
    str2[3999] = '\0';
    int size = 0;
    int i = 0;
    amino_acid_chain *chain;
    char *token = strtok(str1, "-");
    bool matchfound;
    while(token != NULL){
        //error checking
        matchfound = false;
        for(i = 0; i < 22; i++){
            if(strcmp(token, names[i]) == 0){   //if match found, good
                matchfound = true;
                size++;
                token = strtok(str1);
                break;
            }
        }
        if(!matchfound){
            printf("Error: could not rectify input %s at position %d.\n");
            exit(0);
        }    
    }
    /* Now that chain is created empty, fill it */
    chain = malloc(sizeof(amino_acid_chain) + size * sizeof(amino_acid));
    for(i = 0; i < size; i++){
        token = strtok(str2, "-");
        add_to_chain(create_amino_acid(token), chain);
    }
    free(str1);
    free(str2);
    return chain;
}

static void add_to_chain(amino_acid *current, amino_acid_chain *current_chain)
{
    //TODO
}


