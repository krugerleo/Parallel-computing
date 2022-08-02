#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#ifndef N
#define N 5
#endif
#ifndef FS
#define FS 38
#endif

struct node {
   int data;
   int fibdata;
   struct node* next;
};

int fib(int n) {
   int x, y;
   if (n < 2) {
      return (n);
   } else {
      x = fib(n - 1);
      y = fib(n - 2);
	  return (x + y);
   }
}

void processwork(struct node* p) 
{
    if(p == NULL)
        return;
    int n;
    n = p->data;
    p->fibdata = fib(n);
    }


struct node* init_list(struct node* p) {
    int i;
    struct node* head = NULL;
    struct node* temp = NULL;
    
    head = malloc(sizeof(struct node));
    p = head;
    p->data = FS;
    p->fibdata = 0;
    for (i=0; i< N; i++) {
       temp  =  malloc(sizeof(struct node));
       p->next = temp;
       p = temp;
       p->data = FS + i + 1;
       p->fibdata = i+1;
    }
    p->next = NULL;
    return head;
}
struct node* nextNode(int n, struct node* p)
{
    printf("%d\n",n);
    int i = 0;    
    while((i < n) && (p != NULL)){
        p = p->next;
        i++;
    } 
    return p;
};

#define NUM_THREADS 2

int main(int argc, char *argv[]) {
    double start, end;
    struct node *p=NULL;
    struct node *temp=NULL;
    struct node *head=NULL;
    
    printf("Process linked list\n");
    printf("  Each linked list node will be processed by function 'processwork()'\n");
    printf("  Each ll node will compute %d fibonacci numbers beginning with %d\n",N,FS);      

    p = init_list(p);
    head = p;
    int x =0;

    start = omp_get_wtime();
    #pragma omp parallel num_threads(NUM_THREADS) firstprivate(x,p)
    {
        int id,nthrds;
        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();
        while (p != NULL) {
            printf("%d : %d : %d\n",p->data, id, x);
            p = nextNode(id+x,p);
            processwork(p);                        
            x+=2;
        }
    }

    end = omp_get_wtime();
    p = head;
    while (p != NULL) {
    printf("%d : %d\n",p->data, p->fibdata);
    temp = p->next;
    free (p);
    p = temp;
    }  
    free (p);

    printf("Compute Time: %f seconds\n", end - start);

    return 0;
}

