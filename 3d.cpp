//#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <string.h>
#include <float.h>
//#include <alloc.h>
using namespace std;
fstream fp_out;
FILE *fp_rec,*fp_tst;

#define FMAX(a,b) ((a)>(b)? (a):(b))
#define FMIN(a,b) ((a)<(b)? (a):(b))
#define TINY 1e-30
//#define TINY 1.0e-10

double PI = acos(-1.);

void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double b[]);
void nrerror(char error_text[]);
double rel_err(double x,double x0);

////////////////////////////////////////////////////////////////////////////////////////////////
int ** imatrix(long nrow, long ncol);
double ** dmatrix(long nrow, long ncol);
double ***dMallocArray3D(int dim1, int dim2, int dim3);
int ***iMallocArray3D(int dim1, int dim2, int dim3);
void free2D(double **x);
void free2D(int **x);
void free3D(double ***x);
void free3D(int ***x);

void P_avg(int nx,int ny,int nz,double **P,double *Pavg);

int mji2k(int m, int j, int i, int nj, int ni);
void val_Ab_sp(int i0,int j0,int m0,int nx,int ny,int nz,double rmax,double *P0,double D,double lca,double dt,double **A,double *b);

int main()
{
	double *P0,*P;
	double ***Pbig;
	double *swap;
	double dt;
	double lca;
	int nx,ny,nz;
	int m,j,i,m0,j0,i0;
	int k,nk,nbig;
	//int filnum;
	//char fname[20],fname_tmp[10];
	//double Pavg;
	double D;	
	double rmax;
	int ntmax;

	int * indx;	
	double d;
	double ** A;


	clock_t time_used;
	time_used=clock();

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//nx = 1+3+5*2;	计算时间稍长
	nx = 1+3+5;	
	ny = nx;
	nz = nx;
	nk = nx*ny*nz;

	nbig = (nx-1+2)/5;
	ntmax = nx-2;

	i0 = 0;
	j0 = 0;
	m0 = 0;
	

	rmax = (double) ntmax;

	P0 = (double *) malloc((size_t)(nk*sizeof(double)));
	P = (double *) malloc((size_t)(nk*sizeof(double)));
	Pbig = dMallocArray3D(nbig, nbig, nbig);

	indx = (int *) malloc((size_t) (nk*sizeof(int)));
	A = dmatrix(nk, nk);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	dt = 1.5e-6;
	lca = 1e-7;
	D = 1e-9;
	//n_loop_max = n_fine_cell/2 + n_fine_cell;	//31
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout<<"i0 = "<<i0<<"  j0 = "<<j0<<"  m0 = "<<m0<<endl;
	cout<<"nx = "<<nx<<"  ny = "<<ny<<"  nz = "<<nz<<endl;
	//cout<<"n_loop_max = "<<n_loop_max<<endl;


	if((fp_rec = fopen("test_rslt.txt","w+" )) == NULL)
		printf( "The file 'test_rslt.txt' was not opened\n" );
	else
		printf( "The file 'test_rslt.txt' was opened\n" );


	for(k=0;k<nk;k++){
		P[k]  = 0.;
		P0[k] = 0.;
	}
	//P0[j0][i0] = 0.0085*pow(n_fine_cell,2.);
	k = mji2k(m0,j0,i0,ny,nx);
	P0[k] = 1.;
	//P0[k] = 10.;


	dt /= 11.;
	lca /= 11.;

	for(int nt=0;nt<ntmax;nt++){

		for(j=0;j<nk;j++){
			for(i=0;i<nk;i++){
				A[j][i] = 0.;
			}
			A[j][j] = 1.;
		}

		val_Ab_sp(i0,j0,m0,nx,ny,nz,rmax,P0,D,lca,dt,A,P);
		ludcmp(A,nk,indx,&d);	
		lubksb(A,nk,indx,P);	

		swap = P0; P0 = P; P = swap;
		cout<<"nt = "<<nt<<endl;
	}


	for(m=0;m<nbig;m++){
		for(j=0;j<nbig;j++){
			for(i=0;i<nbig;i++){
				Pbig[m][j][i] = 0.;
			}
		}
	}

	int ib,jb,mb;

	for(m=0;m<nz-1;m++){
		mb = (m+2)/5;
		for(j=0;j<ny-1;j++){
			jb = (j+2)/5;
			for(i=0;i<nx-1;i++) {
				ib = (i+2)/5;
				k = mji2k(m,j,i,ny,nx);
				if(mb==0&&jb==0&&ib==0){
					if(m==0&&j==0&&i==0) Pbig[mb][jb][ib] += P0[k];
					else if((m==0&&j==0)||(m==0&&i==0)||(j==0&&i==0)) Pbig[mb][jb][ib] += 2.*P0[k];
					else if(m==0||j==0||i==0) Pbig[mb][jb][ib] += 4.*P0[k];
					else Pbig[mb][jb][ib] += 8.*P0[k];
				}
				else if(mb==0&&jb==0){
					if(m==0&&j==0) Pbig[mb][jb][ib] += P0[k];
					else if(m==0||j==0) Pbig[mb][jb][ib] += 2.*P0[k];
					else Pbig[mb][jb][ib] += 4.*P0[k];
				}
				else if(mb==0&&ib==0){
					if(m==0&&i==0) Pbig[mb][jb][ib] += P0[k];
					else if(m==0||i==0) Pbig[mb][jb][ib] += 2.*P0[k];
					else Pbig[mb][jb][ib] += 4.*P0[k];
				}
				else if(jb==0&&ib==0){
					if(j==0&&i==0) Pbig[mb][jb][ib] += P0[k];
					else if(j==0||i==0) Pbig[mb][jb][ib] += 2.*P0[k];
					else Pbig[mb][jb][ib] += 4.*P0[k];
				}
				else if(mb==0){
					if(m==0) Pbig[mb][jb][ib] += P0[k];
					else Pbig[mb][jb][ib] += 2.*P0[k];
				}
				else if(jb==0){
					if(j==0) Pbig[mb][jb][ib] += P0[k];
					else Pbig[mb][jb][ib] += 2.*P0[k];
				}
				else if(ib==0){
					if(i==0) Pbig[mb][jb][ib] += P0[k];
					else Pbig[mb][jb][ib] += 2.*P0[k];
				}
				else Pbig[mb][jb][ib] += P0[k];	
			}
		}
	}

for(m=0;m<nz;m++){
	for(j=0;j<ny;j++){
		for(i=0;i<nx;i++) {
			k = mji2k(m,j,i,ny,nx);
			fprintf(fp_rec,"%e\t",P0[k]);
		}
		fprintf(fp_rec,"\n");
	}
}

	for(m=0;m<nbig;m++){
		for(j=0;j<nbig;j++){
			for(i=0;i<nbig;i++) {
				fprintf(fp_rec,"%e\t",Pbig[m][j][i]/125.);
			}
			fprintf(fp_rec,"\n");
		}
		fprintf(fp_rec,"\n");
	}

	free(P);
	free(P0);
	free2D(A);
	free3D(Pbig);
	free(indx);

	fclose(fp_rec);

	time_used=clock()-time_used;
	double time_u = (double)time_used/CLOCKS_PER_SEC;
	cout<<endl<<"Total time taken: "<<time_u<<" s"<<endl;
	getchar();
	return 0;
}


void val_Ab_sp(int i0,int j0,int m0,int nx,int ny,int nz,double rmax,double *P0,double D,double lca,double dt,double **A,double *b)
{
	int m,j,i,k;
	double r;

	double dtx2 = lca*lca/D/dt;
	double dtx2p6 = dtx2+6.;

	for(m=1;m<nz-1;m++){
		for(j=1;j<ny-1;j++){
			for(i=1;i<nx-1;i++) {
				k = mji2k(m,j,i,ny,nx);
				r = sqrt(1.*(i-i0)*(i-i0)+1.*(j-j0)*(j-j0)+1.*(m-m0)*(m-m0));
				if(r<=rmax){
					A[k][k] = dtx2p6;
					A[k][k-1] = -1.;
					A[k][k+1] = -1.;
					A[k][mji2k(m,j+1,i,ny,nx)] = -1.;
					A[k][mji2k(m,j-1,i,ny,nx)] = -1.;
					A[k][mji2k(m+1,j,i,ny,nx)] = -1.;
					A[k][mji2k(m-1,j,i,ny,nx)] = -1.;
					b[k] = dtx2*P0[k];
				}
				else {
					b[k] = P0[k];
				}
			}
		}
	}

    m = 0;
	for(j=1;j<ny-1;j++) {
		for(i=1;i<nx-1;i++) {
			k = mji2k(m,j,i,ny,nx);
			r = sqrt(1.*(i-i0)*(i-i0)+1.*(j-j0)*(j-j0));
			if(r<=rmax){
				A[k][k] = dtx2p6;
				A[k][k-1] = -1.;
				A[k][k+1] = -1.;
				A[k][mji2k(m,j+1,i,ny,nx)] = -1.;
				A[k][mji2k(m,j-1,i,ny,nx)] = -1.;
				A[k][mji2k(m+1,j,i,ny,nx)] = -2.;
				b[k] = dtx2*P0[k];
			}
			else {
				b[k] = P0[k];
			}
		}
	}
	m = nz-1;
    for(j=0;j<ny;j++) {
        for(i=0;i<nx;i++) {
            k = mji2k(m,j,i,ny,nx);
            b[k] = P0[k];
        }
    }
    j = 0;
	for(m=1;m<nz-1;m++) {
		for(i=1;i<nx-1;i++) {
			k = mji2k(m,j,i,ny,nx);
			r = sqrt(1.*(i-i0)*(i-i0)+1.*(m-m0)*(m-m0));
			if(r<=rmax){
				A[k][k] = dtx2p6;
				A[k][k-1] = -1.;
				A[k][k+1] = -1.;
				A[k][mji2k(m,j+1,i,ny,nx)] = -2.;
				A[k][mji2k(m-1,j,i,ny,nx)] = -1.;
				A[k][mji2k(m+1,j,i,ny,nx)] = -1.;
				b[k] = dtx2*P0[k];
			}
			else {
				b[k] = P0[k];
			}
		}
	}
	j = ny-1;
    for(m=0;m<nz;m++) {
        for(i=0;i<nx;i++) {
            k = mji2k(m,j,i,ny,nx);
            b[k] = P0[k];
        }
    }
    i = 0;
	for(m=1;m<nz-1;m++) {
		for(j=1;j<ny-1;j++) {
			k = mji2k(m,j,i,ny,nx);
			r = sqrt(1.*(j-j0)*(j-j0)+1.*(m-m0)*(m-m0));
			if(r<=rmax){
				A[k][k] = dtx2p6;
				A[k][k+1] = -2.;
				A[k][mji2k(m,j+1,i,ny,nx)] = -1.;
				A[k][mji2k(m,j-1,i,ny,nx)] = -1.;
				A[k][mji2k(m-1,j,i,ny,nx)] = -1.;
				A[k][mji2k(m+1,j,i,ny,nx)] = -1.;
				b[k] = dtx2*P0[k];
			}
			else {
				b[k] = P0[k];
			}
		}
	}
	i = nx-1;
    for(m=0;m<nz;m++) {
        for(j=0;j<ny;j++) {
            k = mji2k(m,j,i,ny,nx);
            b[k] = P0[k];
        }
    }
    m = 0;
    j = 0;
	for(i=1;i<nx-1;i++) {
		k = mji2k(m,j,i,ny,nx);
		A[k][k] = dtx2p6;
		A[k][k-1] = -1.;
		A[k][k+1] = -1.;
		A[k][mji2k(m,j+1,i,ny,nx)] = -2.;
		A[k][mji2k(m+1,j,i,ny,nx)] = -2.;
		b[k] = dtx2*P0[k];
	}
    m = 0;
    i = 0;
	for(j=1;j<ny-1;j++) {
		k = mji2k(m,j,i,ny,nx);
		A[k][k] = dtx2p6;
		A[k][k+1] = -2.;
		A[k][mji2k(m,j-1,i,ny,nx)] = -1.;
        A[k][mji2k(m,j+1,i,ny,nx)] = -1.;
		A[k][mji2k(m+1,j,i,ny,nx)] = -2.;
		b[k] = dtx2*P0[k];
	}
    j = 0;
    i = 0;
	for(m=1;m<nz-1;m++) {
		k = mji2k(m,j,i,ny,nx);
		A[k][k] = dtx2p6;
		A[k][k+1] = -2.;
		A[k][mji2k(m,j+1,i,ny,nx)] = -2.;
        A[k][mji2k(m-1,j,i,ny,nx)] = -1.;
		A[k][mji2k(m+1,j,i,ny,nx)] = -1.;
		b[k] = dtx2*P0[k];
	}
    m = 0;
    j = 0;
	i = 0;
	k = mji2k(m,j,i,ny,nx);
	A[k][k] = dtx2p6;
	A[k][k+1] = -2.;
	A[k][mji2k(m,j+1,i,ny,nx)] = -2.;
	A[k][mji2k(m+1,j,i,ny,nx)] = -2.;
	b[k] = dtx2*P0[k];
}


int mji2k(int m,int j,int i,int nj,int ni)
{
	int k;
	k = m*nj*ni+j*ni+i;
	return k;
}

void P_avg(int nx,int ny,double **P,double *Pavg)
{
	int i,j;
	*Pavg = 0.;
		for(j=0;j<ny;j++){
			for(i=0;i<nx;i++){
				*Pavg += P[j][i];
			}
		}
	*Pavg /= (double) nx*ny;
}


//=================================================================
// ����ͨ�ó���
double ** dmatrix(long nrow, long ncol)
// allocate a double matrix with subscript range m[nrow][ncol]
{
	long i;
	double **m;
	// allocate pointers to rows
	m=(double **) malloc((size_t)(nrow*sizeof(double*)));
	// allocate rows and set pointers to them
	m[0]=(double *) malloc((size_t)(nrow*ncol*sizeof(double)));
	for(i=1;i<nrow;i++) m[i]=m[i-1]+ncol;
	// return pointer to array of pointers to rows
	return m;
}

int ** imatrix(long nrow, long ncol)
// allocate a double matrix with subscript range m[nrow][ncol]
{
	long i;
	int **m;
	// allocate pointers to rows
	m=(int **) malloc((size_t)(nrow*sizeof(int*)));
	// allocate rows and set pointers to them
	m[0]=(int *) malloc((size_t)(nrow*ncol*sizeof(int)));
	for(i=1;i<nrow;i++) m[i]=m[i-1]+ncol;
	// return pointer to array of pointers to rows
	return m;
}

void free2D(double **x)
//�˶γ���ɱ�֤�ڴ�����
{
	if(x != NULL){
		free(x[0]);
		free(x);
		x = NULL;
	}
}

void free2D(int **x)
//�˶γ���ɱ�֤�ڴ�����
{
	if(x != NULL){
		free(x[0]);
		free(x);
		x = NULL;
	}
}

double ***dMallocArray3D(int dim1, int dim2, int dim3)
//�˶γ���ɱ�֤�ڴ��ַ����
{
	double *** x;
	x = new double **[dim1];
	x[0] = new double *[dim1*dim2];
	x[0][0] = new double[dim1*dim2*dim3];
	for (int i = 0; i < dim1; i++) {
		if (i > 0) {
			x[i] = x[i-1] + dim2;
			x[i][0] = x[i-1][0] + (dim2*dim3);
		}
		for (int j = 1; j < dim2; j++) { 
			x[i][j] = x[i][j-1] + dim3;
		} 
	}
	return x;
}

int ***iMallocArray3D(int dim1, int dim2, int dim3)
//�˶γ���ɱ�֤�ڴ��ַ����
{
	int *** x;
	x = new int **[dim1];
	x[0] = new int *[dim1*dim2];
	x[0][0] = new int[dim1*dim2*dim3];
	for (int i = 0; i < dim1; i++) {
		if (i > 0) {
			x[i] = x[i-1] + dim2;
			x[i][0] = x[i-1][0] + (dim2*dim3);
		}
		for (int j = 1; j < dim2; j++) { 
			x[i][j] = x[i][j-1] + dim3;
		} 
	}
	return x;
}

void free3D(int ***x)
//�˶γ���ɱ�֤�ڴ�����
{
	if(x != NULL){
		free(x[0][0]);
		free(x[0]);
		free(x);
		x = NULL;
	}
}

void free3D(double ***x)
//�˶γ���ɱ�֤�ڴ�����
{
	if(x != NULL){
		free(x[0][0]);
		free(x[0]);
		free(x);
		x = NULL;
	}
}



//////////////////////////////////////////////////////////////
void lubksb(double **a, int n, int *indx, double b[])
/*********************************************************************************************************
  Solves the set of n linear equations A . X = B. Here a[1..n][1..n] is input, not as the matrix A but
  rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input as the
  permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector B, and returns
  with the solution vector X. a, n, and indx are not modified by this routine and can be left in place
  for successive calls with different right-hand sides b. This routine takes into account the possibility
  that b will begin with many zero elements, so it is eficient for use in matrix inversion.
*********************************************************************************************************/
{
int i,ii=-1,ip,j;
   double sum;
   for (i=0;i<n;i++) {
   	/* When ii is set to a positive value, it will become the index of
         the first nonvanishing element of b. We now do the forward
         substitution, equation (2.3.6). The only new wrinkle is to
         unscramble the permutation as we go. */
      ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii>-1) for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i; 	/* A nonzero element was encountered, so from now on we
					will have to do the sums in the loop above. */
		b[i]=sum;
   }
   for (i=n-1;i>=0;i--) { 	/* Now we do the backsubstitution, equation (2.3.7). */
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i]; 	/* Store a component of the solution vector X. */
   } 								/* All done! */
}


void ludcmp(double **a, int n, int *indx, double *d)
/************************************************************************************************************
Given a matrix a[0..n-1][0..n-1], this routine replaces it by the LU decomposition of a rowwise permutation
of itself. a and n are input. a is output, arranged as in equation (2.3.14) above; indx[0..n-1] is an output
vector that records the row permutation effected by the partial pivoting; d is output as +/- 1 depending
on whether the number of row interchanges was even or odd, respectively. This routine is used in
combination with lubksb to solve linear equations or invert a matrix.
************************************************************************************************************/
{
   int i,imax,j,k;
   double big,dum,sum,temp;
   double *vv; 			/* vv stores the implicit scaling of each row. */
   //vv=dvector(n);
   vv=(double *)malloc((size_t) (n*sizeof(double)));
   *d=1.0; 			/* No row interchanges yet. */
   for (i=0;i<n;i++) { 	/* Loop over rows to get the implicit scaling information.  */
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
	   if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
	   /* No nonzero largest element. */
      vv[i]=1.0/big; 		/* Save the scaling. */
	}

	for (j=0;j<n;j++) { 		/* This is the loop over columns of Crout's method. */
	   for (i=0;i<j;i++) { 	/* This is equation (2.3.12) except for i = j. */
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
   	}
		big=0.0; 			/* Initialize for the search for largest pivot element. */
		for (i=j;i<n;i++) { 	/* This is i = j of equation (2.3.12) and i = j +1 : : : N
						of equation (2.3.13). */
			sum=a[i][j];
			for (k=0;k<j;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
		   	/* Is the figure of merit for the pivot better than the best so far? */
				big=dum;
		   	imax=i;
      	}
   	}
   	if (j != imax) { 		/* Do we need to interchange rows? */
			for (k=0;k<n;k++) { 	/* Yes, do so... */
	   		dum=a[imax][k];
	      	a[imax][k]=a[j][k];
	      	a[j][k]=dum;
      	}
      	*d = -(*d); 		/* ...and change the parity of d. */
	   	vv[imax]=vv[j]; 	/* Also interchange the scale factor. */
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		/*
	   If the pivot element is zero the matrix is singular (at least to the precision of the
	   algorithm). For some applications on singular matrices, it is desirable to substitute
	   TINY for zero.
		*/
		if (j != (n-1)) { 			/* Now, finally, divide by the pivot element. */
			dum=1.0/(a[j][j]);
	   	for (i=j+1;i<n;i++) a[i][j] *= dum;
   	}
	} 					/* Go back for the next column in the reduction. */
   free(vv);
}

double rel_err(double x,double x0)
{
	return fabs(x-x0)/FMAX(fabs(x0),TINY);
}


void nrerror(char error_text[])
// Numerical Recipes standard error handler
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}


