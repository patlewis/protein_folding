#include <stdio.h>
#include "proteins.h"

int main(int argc, char **argv){
    amino_acid_chain *chain;
    chain = create_chain();
    free_chain(chain);
    return 0;
}
