#include <stdlib.h>
#include "proxSortedL1.h"

/* ----------------------------------------------------------------------- */
int evaluateProx(double *y, double *lambda, double *x, size_t n, int *order)
/* ----------------------------------------------------------------------- */
{  double   d;
   double  *s     = NULL;
   double  *w     = NULL;
   size_t  *idx_i = NULL;
   size_t *idx_j = NULL;
   size_t  i,j,k;
   int      result = 0;
   
   /* Allocate memory */
   s     = (double *)malloc(sizeof(double) * n);
   w     = (double *)malloc(sizeof(double) * n);
   idx_i = (size_t *)malloc(sizeof(size_t) * n);
   idx_j = (size_t *)malloc(sizeof(size_t) * n);
   
   if ((s != NULL) && (w != NULL) && (idx_i != NULL) && (idx_j != NULL))
   {
      k = 0;
      for (i = 0; i < n; i++)
      {
         idx_i[k] = i;
         idx_j[k] = i;
         s[k]     = y[i] - lambda[i];
         w[k]     = s[k];
         
         while ((k > 0) && (w[k-1] <= w[k]))
         {  k --;
            idx_j[k] = i;
            s[k]    += s[k+1];
            w[k]     = s[k] / (i - idx_i[k] + 1);
         }
         
         k++;
      }
      
      if (order == NULL)
      {  for (j = 0; j < k; j++)
         {  d = w[j]; if (d < 0) d = 0;
            for (i = idx_i[j]; i <= idx_j[j]; i++)
            {  x[i] = d;
            }
         }
      }
      else
      {  for (j = 0; j < k; j++)
         {  d = w[j]; if (d < 0) d = 0;
            for (i = idx_i[j]; i <= idx_j[j]; i++)
            {  x[order[i]] = d;
            }
         }
      }
   }
   else
   {  result = -1;
   }
   
   /* Deallocate memory */
   if (s     != NULL) free(s);
   if (w     != NULL) free(w);
   if (idx_i != NULL) free(idx_i);
   if (idx_j != NULL) free(idx_j);
   
   return result;
}
