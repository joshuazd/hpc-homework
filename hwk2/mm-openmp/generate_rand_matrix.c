#include <stdio.h>
#include <stdlib.h> // srand and rand
#include <time.h>

int main()
{
   time_t t;
   /* Intializes random number generator */
   srand((unsigned) time(&t));

   int i,r,j,n=4000;
   for ( i=0 ; i<n ; i++ ) {
      for ( j=0; j<n ; j++ ) {
         r = rand() % 200;
         printf("%d ",r);
      }
      printf("\n");
   }
 
   return 0;
}
