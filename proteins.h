
typedef struct amino_acid {
    char name[4];           /*the amino acid's name     */
    struct amino_acid *next;       /*next one in the chain     */
    struct amino_acid *previous;   /*previous one in the chain */
} amino_acid;

typedef struct amino_acid_chain{
    amino_acid *first;      /* First amino acid in the chain        */
    amino_acid *last;       /* Last amino acid in the chain         */
} amino_acid_chain;

typedef struct protein {
    amino_acid structure[41][41];
    double energy;
} protein;


