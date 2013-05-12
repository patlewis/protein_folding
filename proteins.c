#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "proteins.h"
#include "list.h"

#define SIZE 41

static const char *names[] =  {"ALA","ARG","ASN","ASP","ASX","CYS","GLU","GLN", \
    "GLX","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","TYR","TRP"};

static const char *philic[] = {"ASN", "ARG", "ASP", "LYS", "GLU", "GLN", "SER" };

amino_acid_chain * create_chain(void){
    amino_acid_chain *chain = (amino_acid_chain *)malloc(sizeof(amino_acid_chain));
    char line[10000];
    printf("Reading input...\n");
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
        printf("%s-", token);
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
    two_d_protein *pro = (two_d_protein *)malloc(sizeof(two_d_protein));
    list_init(&pro->am_ac_list);
    vertex current = array[20][20];
    amino_acid *acid = current.amino;
    while(current.next != NULL){ //does NOT add last amino acid to chain
        //printf("current.amino = %s\n", acid->name);
        //create an amino acid
        amino_acid *ac = new_amino_acid(acid->name);
        ac->x = current.x;
        ac->y = current.y;
        //add it to the list
        list_push_back(&pro->am_ac_list, &ac->elem);
        //set up next iteration
        current = *(current.next);
        acid = current.amino;
    }
    //add last amino acid
    amino_acid *ac = new_amino_acid(acid->name);
    ac->x = current.x;
    ac->y = current.y;
    //add it to the list
    list_push_back(&pro->am_ac_list, &ac->elem);
    pro->energy = energy;
    return *pro;
}


void two_d_protein_free(two_d_protein *pro){
    if(pro != NULL) free(pro);
}

void print_protein(two_d_protein pro, int max){
    amino_acid *acid;
    struct list_elem *e;
    //char **pnames = (char**) malloc(4*20);//list of names
    int i = 0;
    printf("Coordinates: {");
    for(e = list_begin(&(pro.am_ac_list)); e != list_end(&(pro.am_ac_list)) && i < max; e = list_next(e)){
        acid = list_entry(e, struct amino_acid, elem);
        printf("{%d,%d}", acid->x, acid->y); 
      //  *(names + i) = acid->name;
        i++;
        if(i < max) printf(",");
    }
    printf("}\n");    
    //print all except last
    //printf("Names: {");
    //for(i = 0; i < max; i++){
    //    printf("\"%s\"", *(pnames+i));
    //   if(i < (max-1)) printf(",");
    //}
    //printf("}\n");
    printf("Energy: %f\n\n\n", pro.energy);
}

/* DEPRECATED-this structure no longer supported
void print_protein(two_d_protein pro){
    amino_acid current = pro.structure[20][20];
    while(current.next != NULL){
        printf("%d, %d, %s\n", current.x, current.y, current.name);
        current = *(current.next);
    }
    // Print this out because the last one didn't print 
    printf("%d, %d, %s\n", current.x, current.y, current.name);
    printf("Energy: %f\n\n\n", pro.energy);

}
*/


bool less_than(two_d_protein pro1, two_d_protein pro2){
    if (pro1.energy < pro2.energy) return true;
    else return false;
}


