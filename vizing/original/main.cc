#include <iostream.h>
#include "color.h"
#include "graph.h"
main( )
{
    int v, e;
    cin >> v;
    cin >> e;
    Array<Edge> A(e);
   

    for (int i = 0; i < e; i++) {
	int x, y;
	cin >> x >> y;
	Edge tmp(x, y);
	A.insert(i, tmp);
    }


    ColorGraph G(v, A);
    cout << v << endl << e << endl;

    cout << G << endl;

}


