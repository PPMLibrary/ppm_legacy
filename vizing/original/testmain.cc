#include <iostream>
#include "vizing_wrap.cc"

main() {
    int *p, *q;
    int edges[10];
    int i;
    int coloring[15];

    p = new int (4);
    q = new int (5);
    edges[0] = 1; edges[1] = 2;
    edges[2] = 1; edges[3] = 4;
    edges[4] = 1; edges[5] = 3;
    edges[6] = 2; edges[7] = 4;
    edges[8] = 2; edges[9] = 3;
    vizing_coloring_(p,q,edges,coloring);

    for(i=0; i<3*(*q); i+=3) cout << coloring[i] << " <--> " << coloring[i+1] << ": " << coloring[i+2] << endl;
    
    return 0;
}
