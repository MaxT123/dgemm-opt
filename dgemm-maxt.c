const char* dgemm_desc = "Grid search for best block size.";

/*
#if !defined(BLOCK_SIZE)
  #define BLOCK_SIZE 41
#endif
*/

#define min(a,b) (((a) < (b))? (a) : (b))

/*
 * This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N.
 */
static void do_block (int lda, int M, int N, int K, double* A, double* B, double* C)
{
  // For each row i of A
  for (int i = 0; i < M; ++i) {
    //For each column j of B
    for (int j = 0; j < N; ++j) {
      // Compute C(i,j)
      
      /*
      //(?)
      //double cij0 = C[i+j*lda];
      double cij0 = 0, cij1 = 0, cij2 = 0, cij3 = 0;
      // Loop unrolling
      int k=0;
      for (k = 0; k < (K-3); k+=4) {
        cij0 += A[i+k*lda] * B[k+j*lda];
        cij1 += A[i+(k+1)*lda] * B[(k+1)+j*lda];
        cij2 += A[i+(k+2)*lda] * B[(k+2)+j*lda];
        cij3 += A[i+(k+3)*lda] * B[(k+3)+j*lda];
      }
      for (; k<K; ++k) {
          cij0 += A[i+k*lda] * B[k+j*lda];
      C[i+j*lda] += cij0 + cij1 + cij2 + cij3;
      }
      */
      
      
      double cij0=0, cij1=0, cij2=0, cij3=0;
      int k=1;
      for (k=0;k<(K-3);k+=4)
      {
          cij0+=A[i+k*lda] * B[k+j*lda];
          cij1+=A[i+(k+1)*lda] * B[(k+1)+j*lda];
          cij2+=A[i+(k+2)*lda] * B[(k+2)+j*lda];
          cij3+=A[i+(k+3)*lda] * B[(k+3)+j*lda];
      }
      for (;k<K;k++)
          cij0+=A[i+k*lda] * B[k+j*lda];
      C[i+j*lda]+=cij0+cij1+cij2+cij3;
      
    }
  }
}

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format. 
 * On exit, A and B maintain their input values. */  
void square_dgemm (int BLOCK_SIZE, int lda, double* A, double* B, double* C)
{
  // For each block-row of A
  for (int i = 0; i < lda; i += BLOCK_SIZE) {
    // For each block-column of B
    for (int j = 0; j < lda; j += BLOCK_SIZE) {
      // Accumulate block dgemms into block of C
      for (int k = 0; k < lda; k += BLOCK_SIZE) {
        // Correct block dimensions if block "goes off edge of" the matrix
        int M = min (BLOCK_SIZE, lda-i);
        int N = min (BLOCK_SIZE, lda-j);
        int K = min (BLOCK_SIZE, lda-k);
        // Perform individual block dgemm
        do_block(lda, M, N, K, A + i + k*lda, B + k + j*lda, C + i + j*lda);
      }
    }
  }
}
