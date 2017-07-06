//This program provides a way to verify that the graph-edge-coloring
//program is correct!


#include <iostream.h>

enum Boolean {FALSE=0, TRUE=1};


Boolean Verify()
{

    int v, e;
    cin >> v;
    cin >> e;

    int A[e], B[e], color[e];

    for (int i=0; i<e; i++) 
	cin >> A[i] >> B[i] >> color[i];





    for (int j=0; j<e; j++) 
	for (int k=0; k<e; k++) {
	    if (A[j]==A[k] && B[j]!=B[k] && color[j] == color[k])
		return FALSE;
	    if (A[j] != A[k] && B[j]==B[k] && color[j] == color[k])
		return FALSE;
	    if (A[j] == B[k] && A[k] != B[j] && color[j] == color[k])
		return FALSE;
	}
    
    return TRUE;
}


main()
{

    if (Verify())
	cout << "Correct!"<<endl;
    else
	cout <<"Wrong!"<<endl;
}
	



