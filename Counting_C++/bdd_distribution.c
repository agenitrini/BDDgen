#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <sys/time.h>
#include <omp.h>
/*-------------------------------------------*/
int bit_length(slong k) {
  int b = 0;
  while (k != 0) {
    b++;
    k = k >> 1;
  }
  return b;
}
/*-------------------------------------------*/
int max_size(int k) {
  //  2**(k-(k-k.bit_length()+1).bit_length()+1)+2**2**((k-k.bit_length()+1).bit_length()-1)-1 
  unsigned int theta = (unsigned int) ((int) k - bit_length(k)+1);
  unsigned int L = (unsigned int) ((int) k-bit_length(theta)+1);
  int M = (unsigned int) ((int) (1<< L)+(1<<(1<< ((int) bit_length(theta)-1))) -3);
  return M;
}
/*--------------------------------------------*/
void set_max_profile(int *p, int k) {
  int N = max_size(k);
  int top = 0;
  int bottom = N+2;
  int i = 0;
  while (top < N) {
    int m = FLINT_MIN(top+1, bottom-ceil(sqrt(bottom)));
    p[i] = m;
    top += m;
    bottom -=m;
    i++;
  }
}
/*---------------------------------------------*/
void print_array(int *A, int n) {
  for (int i =0; i < n; i++){
    printf("%d ", A[i]);
  }
}
/* ---------------------------------------------------------------------*/
void next_bis(fmpz_poly_t *src, int m, int nb_nodes) {
  // tradeoff speed for memory: next iteration is computed in place
  fmpz_poly_t x, res; // temporary variable
  fmpz_poly_init(x);
  fmpz_poly_init(res);
  int deg_src = FLINT_MIN(2*m-2, m-1+ nb_nodes);
  int deg_dest = FLINT_MIN(2*m, m+nb_nodes);
  /* int Delta = deg_dest - deg_src; */
  /* Delta is 1 or 2 */
  for (int i=deg_dest; i > 0; i--) { 
      fmpz_poly_zero(res);
    if (i < deg_src+1) {
      fmpz_poly_derivative(x, src[i]); // x = d/du src[1]
      fmpz_poly_shift_left(res, x, 1); // dest[1] = u d/du src[1]
    }
    if (i>=1 && i-1 < deg_src+1) {
      fmpz_poly_add(res, res, src[i-1]);
    }
    fmpz_poly_shift_left(x, res, 1); // x = u * dest[i]
    fmpz_poly_sub(res, res, x); // dest[i] = dest[i] - u dest[i]
    if (i >= 2 && i-2 < deg_src+1) {
      fmpz_poly_shift_left(x, src[i-2], 1); //
      fmpz_poly_add(res, res, x);
    }
    fmpz_poly_set_trunc(src[i], res, nb_nodes+1);
  }
  fmpz_poly_zero(src[0]);
  fmpz_poly_clear(x);
  fmpz_poly_clear(res);
}
/*--------------------------------------------*/
void count(int k, fmpz_poly_t R) {
  int N = max_size(k);
  int MAX_PROFILE[k];
  set_max_profile(MAX_PROFILE, k);
  printf("max size=%d, max profile: [", N);
  print_array(MAX_PROFILE, k);
  printf("]\n");
  fmpz_poly_t *previous, *current;

  previous = malloc((N+2) * sizeof(fmpz_poly_t));
  int len_prev = N+2;
  for (int r=0; r < N+2; r++) {
    fmpz_poly_init(previous[r]);
  }
  fmpz_t f;
  fmpz_init(f);
  fmpz_set_si(f, 1);
  for (int m=0; m < N+2; m++) {
    fmpz_poly_set_fmpz(previous[m], f);
    fmpz_mul_ui(f, f, 2);
  }
  fmpz_clear(f);

  fmpz_poly_t *src;
  
  src = malloc((N+2) * sizeof(fmpz_poly_t));
   for (int i=0; i < N+2; i++)
    fmpz_poly_init(src[i]);
  
  int T=N, B=0;
  for (int exponent=1; exponent < k+1; exponent++) {
    struct timeval stop, start;
    gettimeofday(&start, NULL);
    printf("level %d", exponent);
    T -= MAX_PROFILE[k-exponent];
    B += MAX_PROFILE[k-exponent];
    int nb_nodes = MAX_PROFILE[k-exponent];
    printf(" (T=%d, B=%d, nb_nodes=%d) ", T, B, nb_nodes);
    fmpz_poly_one(src[0]); // initiate
    current = malloc((T+2) * sizeof(fmpz_poly_t));
    for (int i =0; i < T+2; i++) {
      fmpz_poly_init(current[i]);
      fmpz_poly_zero(current[i]);
    }
    fmpz_poly_zero(current[0]);
    for (int m=1; m < T+2; m++) {
      if (T < 20 || (m % ((T+1) /10) == 1)) {
	printf(".");
      }
      int dX = m+FLINT_MIN(m, nb_nodes);
      next_bis(src, m, nb_nodes);
      fmpz_poly_t partial_Sum;
      fmpz_poly_t x; /* Pu */
#pragma omp parallel private(partial_Sum, x) shared(current)
      {
	fmpz_poly_init(partial_Sum);
	fmpz_poly_zero(partial_Sum);
	fmpz_poly_init(x);
#pragma omp for nowait 
	  for (int j=0; j < dX+1; j++) {
	    fmpz_poly_mullow(x, src[j], previous[j], B+1); // Multiplication keeping only the lowest B+1 coefficients of the product 
	    fmpz_poly_add(partial_Sum, partial_Sum, x);
	  }
#pragma omp critical
	{
	  fmpz_poly_add(current[m], current[m], partial_Sum);
	}
	fmpz_poly_clear(x);
	fmpz_poly_clear(partial_Sum);
      }
    }
    for (int i=0; i < len_prev; i++) {
      fmpz_poly_clear(previous[i]);
    }
    free(previous);
    previous = current;
    len_prev = T+2;
    gettimeofday(&stop, NULL);
    long unsigned int delta_us = (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec;
    printf(" in %g seconds (approx)\n", (float) delta_us / 1000000.); 
  }
  fmpz_poly_set(R, current[1]);

  for (int i=0; i < len_prev; i++) {
    fmpz_poly_clear(previous[i]);
  }
  free(previous);
  free(src);
}
/*---------------------------------------------------*/
int main(int argc, char **argv) {
  int k = atol(argv[1]);
#if defined(_OPENMP)
  if (argc > 2) {
    int nb_threads = atol(argv[2]);
    //omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(nb_threads);
    printf("Setting to %d threads\n", nb_threads);
  }
#endif
  setvbuf(stdout, NULL, _IONBF, 0);

  fmpz_poly_t M;
  fmpz_poly_init(M);

  count(k, M);
  fmpz_t coeff;
  fmpz_init(coeff);
  for (int i = 0; i < fmpz_poly_length(M); i++) {
    fmpz_poly_get_coeff_fmpz(coeff, M, i);
    printf("%i\t", i);
    fmpz_print(coeff);
    printf("\n");
  }
  fmpz_clear(coeff);
  fmpz_poly_clear(M);
  return 0;
}
