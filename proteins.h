#include "list.h"


typedef struct amino_acid {
    char name[4];           /*the amino acid's name     */
    struct amino_acid *next;       /*next one in the chain     */
    struct amino_acid *previous;   /*previous one in the chain */
    struct list_elem elem;
} amino_acid;

typedef struct amino_acid_chain{
    struct list amino_acid_list;
} amino_acid_chain;

typedef struct two_d_protein {
    amino_acid structure[41][41];
    double energy;
    struct list_elem *elem;
} two_d_protein;

typedef struct three_d_protein {
    amino_acid structure[41][41][41];
    double energy;
    struct list_elem *elem;
} three_d_protein;


/*
 * Creates a new amino acid with the given label.  
 */
amino_acid *new_amino_acid(char *label);

/*
 * Creates a textual representation of an amino acid chain
 * based on text from the user input (standard in).  Does
 * some error checking to make sure that the input is formatted
 * corectly.
 */
amino_acid_chain * create_chain(void);

/*
 * Frees all the memory assiociated with this amino acid chain
 */
void free_chain(amino_acid_chain *chain);
