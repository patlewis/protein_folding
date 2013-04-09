#include <stdio.h>
#include <stdlib.h>
#include "proteins.h"

int main(int argc, char **argv){
    amino_acid_chain *chain;
    chain = create_chain();
    two_d_protein *pro = (two_d_protein *)malloc(sizeof(struct two_d_protein));
    int i, j;
    j = 9;
    for(i = 0; i < 10; i++){
        pro->structure[i][i] = *(new_amino_acid("LEU"));
        if(i != j) pro->structure[j][i] = *(new_amino_acid("ALA"));
        j--;
    }
    pro->energy = -7.171718191017;
    print_protein(*pro);
    free_chain(chain);
    return 0;
}
