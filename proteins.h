#include "list.h"

#define SIZE 41


typedef struct amino_acid {
    char name[4];           /*the amino acid's name     */
    int hydro;              /*0 for hydrophobic, 1 for hydrophilic */
    int x;  //coordinates in the lattice
    int y;
    struct list_elem elem;
    struct amino_acid *next;//for printing purposes
} amino_acid;

typedef struct amino_acid_chain{
    struct list amino_acid_list;
} amino_acid_chain;

typedef struct vertex{
    amino_acid *amino;
    struct vertex *prev;
    struct vertex *next;
    //coordinates
    int x;
    int y;
} vertex;

typedef struct two_d_protein {
    struct list am_ac_list; 
    double energy;
    struct list_elem *elem;
} two_d_protein;

typedef struct three_d_protein {
    vertex structure[SIZE][SIZE][SIZE];
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

/*
 * Creates a two-dimensional protein object for ease
 * of representation.
 */
two_d_protein two_d_protein_create(struct vertex array[SIZE][SIZE], double energy);


/*
 * Prints out a text-based representation of the protein structure.  At this 
 * point in time, it's a list of coordinates for each amino acid.
 * Has to take a maximum size as a parameter because for some reason the 
 * lists are a little funky
 */
void print_protein(two_d_protein pro, int max);

/*
 * Frees all the memory associated with a protein structure
 */
void two_d_protein_free(two_d_protein *);
/*
 * Used to sort protein structures for qsort
 */
bool less_than(two_d_protein pro1, two_d_protein pro2);
