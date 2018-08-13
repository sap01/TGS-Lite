/*
 * This source code is inspired by
 * that of R package 'bnstruct'
 * version 1.0.2
 */
#include <R.h>
#include <Rinternals.h>

SEXP fbp( SEXP aflml )
{
	R_len_t ni, si, bi, pos, pos2; 
	R_len_t n_nodes = nrows(aflml); 
	R_len_t pow_nodes = ncols(aflml);
	R_len_t mask = 1; 
	int *bps; // force to 32bit integer
	double *ml = REAL(aflml), *bss;
	SEXP res;
	
	// allocate output
	PROTECT( res = allocMatrix(INTSXP, n_nodes, pow_nodes) );
	bps = INTEGER(res);
	memset( bps, 0, sizeof(int) * n_nodes * pow_nodes );
	
	// allocate internal structures 
	bss = (double *) R_alloc( n_nodes * pow_nodes, sizeof(double) );
	memset( bss, 0, sizeof(double) * n_nodes * pow_nodes );
	
	for( ni = 0; ni < n_nodes; ni++ )
		for( si = 0; si < pow_nodes; si++ )
		{
			/* a node can not be parent of itself */
			if( si & (mask<<ni) )
				continue;
			
			pos = si*n_nodes + ni;
			bps[pos] = si+1;
			bss[pos] = ml[pos];
			
			/* cycle through all parents and remove one in turn */
			for( bi = 0; bi < n_nodes; bi++ )
				if( si & (mask<<bi) )
				{
					pos2 = ( si ^ (mask<<bi) )*n_nodes + ni; /* remove bi-th parent */
					if( bss[pos2] > bss[pos] )
					{
						bps[pos] = bps[pos2];
						bss[pos] = bss[pos2];
					}
				}
			}
	
	UNPROTECT(1);
	return( res );
}

SEXP fbs( SEXP bps, SEXP aflml )
{
	R_len_t si, sink, nop, pos;
	R_len_t n_nodes = nrows(aflml); 
	R_len_t pow_nodes = ncols(aflml);
	R_len_t mask = 1;
	double skore;
	int *sinks; // force to 32bit integer
	int *ps = INTEGER(bps);
	double *ml = REAL(aflml), *scores;
	SEXP res;
	
	// allocate output
	PROTECT( res = allocVector(INTSXP, pow_nodes) );
	sinks = INTEGER(res);
	for( si = 0; si < pow_nodes; si++ ) sinks[si] = -1;
		
	// allocate internal structures 
	scores = (double *) R_alloc( pow_nodes, sizeof(double) );
	memset( scores, 0, sizeof(double) * pow_nodes );
	
	for( si = 0; si < pow_nodes; si++ )
		for( sink = 0; sink < n_nodes; sink++ )
			if( si & (mask<<sink) )
			{
				nop = si ^ (mask<<sink);
				skore = scores[nop] + ml[(ps[nop*n_nodes+sink]-1)*n_nodes+sink];
				
				if( sinks[si] == -1 || skore > scores[si] )
				{
					scores[si] = skore;
					sinks[si] = sink+1;
				}
			}
			
	UNPROTECT(1);
	return( res );
}

SEXP na_rows_int( SEXP mat )
{
	int i,j,stride, *res;
	int n_rows = nrows( mat );
	int n_cols = ncols( mat );
	int * m = INTEGER( mat );
	SEXP(result);
	
	// allocate oputput
	PROTECT(result = allocVector(INTSXP, n_rows));
	res = INTEGER(result);
	memset( res, 0, sizeof(int) * n_rows );
	
	for( j = 0; j < n_cols; j++ )
	{
		stride = n_rows * j;
		for( i = 0; i < n_rows; i++ )
			res[i] |= (int)(NA_INTEGER == (m[i+stride]));
	}
	
	UNPROTECT(1);
	return( result );
}

SEXP fumt_mask( SEXP n_elements, SEXP pattern )
{
	int i, times, ind, stride, start;
	int bitmask = 1;
	int n_el = asInteger(n_elements);
	PROTECT( pattern = coerceVector(pattern,INTSXP) );
	int * pat = INTEGER(pattern);
	int l_pat = LENGTH(pattern);
	int pow_el = bitmask<<n_el; // faster than power for int
	SEXP result;
	
	// allocate output
	PROTECT(result = allocVector(INTSXP, pow_el));
	int * mask = INTEGER(result);
	memset( mask, 0, sizeof(int) * pow_el );
	
	// assign pattern
	for( i = 0; i < l_pat; i++ )
		mask[ bitmask << ((int)pat[i]-1) ] = 1;
	
	// compute fumt
	for( i = 0; i < n_el; i++ )
	{
		stride = bitmask<<i;
		for( times = 0; times < (bitmask<<(n_el-i-1)); times++ )
		{
			start = times * (bitmask<<(i+1));
			for( ind = 0; ind < stride; ind++ )
				mask[ start + ind + stride ] += mask[ start + ind ];
		}
	}
	
	UNPROTECT(2);
	return(result);
}

SEXP all_fam_log_marg_lik( SEXP data, SEXP node_sizes, SEXP imp_fam_mask, SEXP iss, SEXP func ) {
	
	unsigned int i,j,n_pa,pos;
	
	// Begin: Get inputs
	unsigned int * d = INTEGER(data);
	unsigned int n_nodes = ncols(data);
	
	// No. of obs.
	unsigned int n_ex = nrows(data);	 
	
	unsigned int * ns = INTEGER(node_sizes);
	
	// Impossible family mask
	unsigned int * ifm = INTEGER(imp_fam_mask);
	unsigned int pow_nodes = ncols(imp_fam_mask);
	
	double alpha = *(REAL(iss));
  
  /*
   * 0 for BDeu
   * 1 for AIC
   * 2 for BIC
   */
  unsigned int scoring_func = *INTEGER(func);
	
	unsigned int * pa = (unsigned int *) R_alloc( n_nodes, sizeof(unsigned int) );
	// End: Get inputs
	
	// allocate and initialize output
	SEXP result;
	PROTECT( result = allocMatrix(REALSXP, n_nodes, pow_nodes) );
	double * aflml = REAL(result);
	for( i = 0; i < n_nodes*pow_nodes; i++ )
		aflml[i] = R_NegInf;
	
	// compute log likelihood
	for( i = 0; i < n_nodes; i++ )
		for( j = 0; j < pow_nodes; j++ ) {
			pos = j*n_nodes + i;
			if( ifm[pos] ) {
				/* 
				 * 'n_pa' = No. of parents of node 'i'.
				 * 'pa' = Arr of node indices of parents in the j^th parent set. 
				 * k^th elt of the arr contains node index of the k^th parent. 
				 * The first 'n_pa' number of 
				 * elts. contain node indices of parents.
				 */
				// Rprintf("get bits\n");
				n_pa = get_bits( j, pa, n_nodes );

			  /*
			   * Calc score when the i-th node is parented by the j-th parent set.
			   * The j^th parent set is represented by 'pa'.
			   */			  
				// Rprintf("log lik, node %d, n parents %d\n",i,n_pa);
				aflml[pos] = score_node_1(d, n_nodes, n_ex, ns, i, pa, n_pa, scoring_func, alpha);
        // bdeu_score( d, n_nodes, n_ex, ns, i, pa, n_pa, alpha );
				// Rprintf("end\n");
			}
		}
	
	UNPROTECT(1);
	return result;	
}

/* 
 * Input:
 * 'word': Parent set index in decimal.
 * 'bits':  Unsigned int dyn arr of length total no. of nodes in the data. 
 * 'size': Number of noeds in the data.
 * 
 * Output:
 * Returns 'count': No. of parents.
 * Modifies 'bits': Arr of node indices of parents. i^th elt of the arr
 * contains node index of the i^th parent. The first 'count' number of 
 * elts. contain node indices of parents.
 */
unsigned int get_bits( unsigned int word, unsigned int * bits, unsigned int size ) {
	unsigned int i, count = 0, bitmask = 1;
	for( i = 0; i < size; i++ ) {
	  
	  /* 
	   * Bit-wise binary AND, and
	   * bit-wise binary 'i' no. of left-shifts.
	   * Does not modify 'bitmask'.
	   */
	  if( word & bitmask<<i ) {
	    bits[count++] = i;
	  }
	}
	
	return count;
}

double score_node_1( int* data, int ncols_data, int nrows_data, int* node_sizes, unsigned int ni, int* pars, int length_pars, int func, double ess )
{
  double score;
  
  switch (func)
  {
  case 0 : score = bdeu_score( data, ncols_data, nrows_data, node_sizes,
                               ni, pars, length_pars, ess );
    break;
    
  case 1 : score = log_likelihood( data, ncols_data, nrows_data, node_sizes,
                                   ni, pars, length_pars, 0.5*log(nrows_data) );
    break;
    
  case 2 : score = log_likelihood( data, ncols_data, nrows_data, node_sizes,
                                   ni, pars, length_pars, 1.0 );
    break;
  }
  
  return score;
}


