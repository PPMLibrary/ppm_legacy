/*****************************************************************************/
/*                                                                           */
/* This program finds a valid edge coloring for graphs of class one or       */
/* determines that a graph is of class two.                                  */
/*                                                                           */
/* INPUT: A set of lines ending with a EOF character (^D) specifying a       */
/* graph in the following manner:                                            */
/*                                                                           */
/* c followed by a comment no longer than 255 characters                     */
/* p followed by the number of vertices,                                     */
/* q followed by the number of edges, and                                    */
/* e followed by two numbers from [0, p) specifying the vertices of an edge. */
/*                                                                           */
/* It is assumed that unique p and q lines precede any e lines, and that the */
/* number of e lines is exactly q.                                           */
/*                                                                           */
/* OUTPUT: A line specifying the class of the graph, the numbers of vertices */
/* and edges, and the maximum degree,                                        */
/* a set of lines specifying the reordered edges, and for class one graphs,  */
/* their colors, and                                                         */
/* for class one graphs, a line specifying the backtracking coordinate of    */
/* the solution.                                                             */
/*                                                                           */
/*****************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/times.h>


/*****************************************************************************/
/* Bit Array Definitions                                                     */
/*****************************************************************************/

#define bitsInWord 32
#define logBitsInWord 5

#define allocate(a, bits) a=(int*)calloc(bits/bitsInWord+1, sizeof(int))
#define reset(a, bit) a[bit>>logBitsInWord]&=~(1<<(bit&(bitsInWord-1)))
#define set(a, bit) a[bit>>logBitsInWord]|=1<<(bit&(bitsInWord-1))
#define get(a, bit) (a[bit>>logBitsInWord]&(1<<(bit&(bitsInWord-1))))


/*****************************************************************************/
/* Variable Declarations                                                     */
/*****************************************************************************/

struct vertex {
  unsigned int degree;
};

struct edge {
  unsigned int vertex1;
  unsigned int vertex2;
  unsigned int color;
  unsigned int coordinate;
  unsigned int dblMaxChoic;
};

main()
{
  vertex *v, *vi;
  edge *e, *ei, *ej, *ek, *next;
  unsigned int p, q, Delta, firstVertex, color, i, j, k, t, v1, v2, bit1, bit2;
  int *Used;
  char c, comment[256];
  tms start, finish;
  clock_t time;


  /***************************************************************************/
  /* Reading the Input                                                       */
  /***************************************************************************/

  while(scanf("%c", &c)!=EOF) {
    if(c=='e') {
      scanf("%d%d", &v1, &v2); 
      (*ek).vertex1=v1; (*ek).vertex2=v2; ek++;
      v[v1].degree++; v[v2].degree++;
    }
    else if(c=='c') gets(comment);
    else if(c=='q') { scanf("%d", &q); ek=e=(edge*)calloc(q, sizeof(edge)); }
    else if(c=='p') { scanf("%d", &p); v=(vertex*)calloc(p, sizeof(vertex)); }
  }


  /***************************************************************************/
  /* Initialization                                                          */
  /***************************************************************************/

  times(&start);
  Delta=0; vi=v;
  for(i=0; i<p; i++) {
    if((*vi).degree>Delta) { Delta=(*vi).degree; firstVertex=i; }
    vi++;
  }
  allocate(Used, (p*Delta));
  t=2*Delta; ek=e; for(k=0; k<q; k++) { (*ek).dblMaxChoic=t; ek++; }
  i=0; ek=ei=e;
  while(i<Delta) {
    if((*ek).vertex1==firstVertex || (*ek).vertex2==firstVertex) {
      t=(*ek).vertex1; (*ek).vertex1=(*ei).vertex1; (*ei).vertex1=t;
      t=(*ek).vertex2; (*ek).vertex2=(*ei).vertex2; (*ei).vertex2=t;
      (*ek).dblMaxChoic=(*ei).dblMaxChoic; (*ei).dblMaxChoic=0;
      (*ei).color=i;
      bit1=(*ei).vertex1*Delta+i; set(Used, bit1);
      bit2=(*ei).vertex2*Delta+i; set(Used, bit2);
      ej=&e[i+1];
      for(j=i+1; j<q; j++) {
	if((*ei).vertex1==(*ej).vertex1 || (*ei).vertex1==(*ej).vertex2 ||
	   (*ei).vertex2==(*ej).vertex1 || (*ei).vertex2==(*ej).vertex2)
	  (*ej).dblMaxChoic--;
	ej++;
      }
      i++; ei++;
    }
    ek++;
  }


  /***************************************************************************/
  /* Rordering Edges                                                         */
  /***************************************************************************/

  ek=&e[Delta];
  for(k=Delta; k<q; k++) {
    next=ek; ej=&e[k+1];
    for(j=k+1; j<q; j++) {
      if((*ej).dblMaxChoic<(*next).dblMaxChoic) next=ej;
      ej++;
    }
    t=(*next).vertex1; (*next).vertex1=(*ek).vertex1; (*ek).vertex1=t;
    t=(*next).vertex2; (*next).vertex2=(*ek).vertex2; (*ek).vertex2=t;
    t=(*next).dblMaxChoic; (*next).dblMaxChoic=(*ek).dblMaxChoic;
    (*ek).dblMaxChoic=t;
    ej=&e[k+1];
    for(j=k+1; j<q; j++) {
      if((*ek).vertex1==(*ej).vertex1 || (*ek).vertex1==(*ej).vertex2 ||
	 (*ek).vertex2==(*ej).vertex1 || (*ek).vertex2==(*ej).vertex2)
	(*ej).dblMaxChoic--;
      ej++;
    }
    ek++;
  }


  /***************************************************************************/
  /* Backtracking                                                            */
  /***************************************************************************/

  k=Delta; ek=&e[k];
  while(Delta<=k && k<q) {
    color=(*ek).color; 
    bit1=(*ek).vertex1*Delta+color; bit2=(*ek).vertex2*Delta+color;
    while((get(Used, bit1) || get(Used, bit2)) && color<Delta) {
      color++; bit1++; bit2++;
    }
    if(color<Delta) {
      (*ek).color=color; set(Used, bit1); set(Used, bit2);
      k++; ek++;
    }
    else {
      (*ek).color=0; (*ek).coordinate=0;
      k--; ek--;
      bit1=(*ek).vertex1*Delta+(*ek).color; reset(Used, bit1);
      bit2=(*ek).vertex2*Delta+(*ek).color; reset(Used, bit2);
      (*ek).color++; (*ek).coordinate++;
    }
  }
  times(&finish);

  /***************************************************************************/
  /* Output                                                                  */
  /***************************************************************************/

  if(k==q) {
    printf("Class one graph. p=%d, q=%d, Delta=%d.\n", p, q, Delta);
    ei=e;
    for(i=0; i<q; i++) {
      printf("e %d %d\t%d\n", (*ei).vertex1, (*ei).vertex2, (*ei).color);
      ei++;
    }
    printf("coordinate");
    ei=e; for(i=0; i<q; i++) { printf(" %d", (*ei).coordinate); ei++; }
    printf("\n");
  }
  else {
    printf("Class two graph. p=%d, q=%d, Delta=%d.\n", p, q, Delta);
    ei=e;
    for(i=0; i<q; i++) { 
      printf("e %d %d\n", (*ei).vertex1, (*ei).vertex2); 
      ei++;
    }
  }
  time=finish.tms_utime-start.tms_utime+finish.tms_stime-start.tms_stime;
  printf("elapsed time %d.%d sec\n", time/60, (time%60)*100/60);
}

