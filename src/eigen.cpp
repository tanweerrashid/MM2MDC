#include "eigen.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
//#include "Octree.h"

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

// for reducing two upper triangular systems of equations into 1
void qr ( float *mat1, float *mat2, float *rvalue )
{
	int i, j;
	float temp1 [ 8 ] [ 4 ];

	for ( i = 0; i < 4; i++ )
	{
		for ( j = 0; j < i; j++ )
		{
			temp1 [ i ] [ j ] = 0;
			temp1 [ i + 4 ] [ j ] = 0;
		}
		for ( j = i; j < 4; j++ )
		{
			temp1 [ i ] [ j ] = mat1 [ ( 7 * i - i * i ) / 2 + j ];
			temp1 [ i + 4 ] [ j ] = mat2 [ ( 7 * i - i * i ) / 2 + j ];
		}
	}

	qr ( temp1, 8, rvalue );
}

void qr ( float *mat, float *eqn )
{
	int i, j, k;
	float a, b, mag, temp;
	float eqs [ 4 ] [ 4 ];

	for ( i = 0; i < 4; i++ )
	{
		for ( j = 0; j < i; j++ )
		{
			eqs [ i ] [ j ] = 0;
		}
		for ( j = i; j < 4; j++ )
		{
			eqs [ i ] [ j ] = mat [ ( 7 * i - i * i ) / 2 + j ];
		}
	}

	for ( i = 0; i < 4; i++ )
	{
		a = eqs [ i ] [ i ];
		b = eqn [ i ];

		if ( fabsf ( a ) > 0.000001f || fabsf ( b ) > 0.000001f )
		{
			mag = sqrtf ( a * a + b * b );
			a /= mag;
			b /= mag;

			for ( k = 0; k < 4; k++ )
			{
				temp = a * eqs [ i ] [ k ] + b * eqn [ k ];
				eqn [ k ] = b * eqs [ i ] [ k ] - a * eqn [ k ];
				eqs [ i ] [ k ] = temp;
			}
		}

		for ( j = i - 1; j >= 0; j-- )
		{
			if ( eqs [ j ] [ j ] < 0.000001f && eqs [ j ] [ j ] > -0.000001f )
			{
				a = eqs [ i ] [ i ];
				b = eqs [ j ] [ i ];

				if ( fabsf ( a ) > 0.000001f || fabsf ( b ) > 0.000001f )
				{
					mag = sqrtf ( a * a + b * b );
					a /= mag;
					b /= mag;

					for ( k = 0; k < 4; k++ )
					{
						temp = a * eqs [ i ] [ k ] + b * eqs [ j ] [ k ];
						eqs [ j ] [ k ] = b * eqs [ i ] [ k ] - a * eqs [ j ] [ k ];
						eqs [ i ] [ k ] = temp;
					}
				}
			}
		}
	}

	mat [ 0 ] = eqs [ 0 ] [ 0 ];
	mat [ 1 ] = eqs [ 0 ] [ 1 ];
	mat [ 2 ] = eqs [ 0 ] [ 2 ];
	mat [ 3 ] = eqs [ 0 ] [ 3 ];
	mat [ 4 ] = eqs [ 1 ] [ 1 ];
	mat [ 5 ] = eqs [ 1 ] [ 2 ];
	mat [ 6 ] = eqs [ 1 ] [ 3 ];
	mat [ 7 ] = eqs [ 2 ] [ 2 ];
	mat [ 8 ] = eqs [ 2 ] [ 3 ];
	mat [ 9 ] = eqs [ 3 ] [ 3 ];
}

// WARNING: destroys eqs in the process
void qr ( float eqs[][4], int num, float *rvalue )
{
	int i, j, k;

	qr ( eqs, num );
	for ( i = 0; i < 10; i++ )
	{
		rvalue [ i ] = 0;
	}

	k = 0;
	for ( i = 0; i < num && i < 4; i++ )
	{
		for ( j = i; j < 4; j++ )
		{
			rvalue [ k++ ] = eqs [ i ] [ j ];
		}
	}
}

void qr ( float eqs[][4], int num )
{
	int i, j, k;
	float a, b, mag, temp;

	for ( i = 0; i < 4 && i < num; i++ )
	{
		for ( j = i + 1; j < num; j++ )
		{
			a = eqs [ i ] [ i ];
			b = eqs [ j ] [ i ];

			if ( fabsf ( a ) > 0.000001f || fabsf ( b ) > 0.000001f )
			{
				mag = sqrtf ( a * a + b * b );
				a /= mag;
				b /= mag;

				for ( k = 0; k < 4; k++ )
				{
					temp = a * eqs [ i ] [ k ] + b * eqs [ j ] [ k ];
					eqs [ j ] [ k ] = b * eqs [ i ] [ k ] - a * eqs [ j ] [ k ];
					eqs [ i ] [ k ] = temp;
				}
			}
		}
		for ( j = i - 1; j >= 0; j-- )
		{
			if ( eqs [ j ] [ j ] < 0.000001f && eqs [ j ] [ j ] > -0.000001f )
			{
				a = eqs [ i ] [ i ];
				b = eqs [ j ] [ i ];

				if ( fabsf ( a ) > 0.000001f || fabsf ( b ) > 0.000001f )
				{
					mag = sqrtf ( a * a + b * b );
					a /= mag;
					b /= mag;

					for ( k = 0; k < 4; k++ )
					{
						temp = a * eqs [ i ] [ k ] + b * eqs [ j ] [ k ];
						eqs [ j ] [ k ] = b * eqs [ i ] [ k ] - a * eqs [ j ] [ k ];
						eqs [ i ] [ k ] = temp;
					}
				}
			}
		}
	}

}

void jacobi ( float u[][3], float d[], float v[][3] )
{
	int j, iq, ip, i;
	float tresh, theta, tau, t, sm, s, h, g, c, b [ 3 ], z [ 3 ];
	float a [ 3 ] [ 3 ];

	a [ 0 ] [ 0 ] = u [ 0 ] [ 0 ];
	a [ 0 ] [ 1 ] = u [ 0 ] [ 1 ];
	a [ 0 ] [ 2 ] = u [ 0 ] [ 2 ];
	a [ 1 ] [ 0 ] = u [ 1 ] [ 0 ];
	a [ 1 ] [ 1 ] = u [ 1 ] [ 1 ];
	a [ 1 ] [ 2 ] = u [ 1 ] [ 2 ];
	a [ 2 ] [ 0 ] = u [ 2 ] [ 0 ];
	a [ 2 ] [ 1 ] = u [ 2 ] [ 1 ];
	a [ 2 ] [ 2 ] = u [ 2 ] [ 2 ];

	for ( ip = 0; ip < 3; ip++ ) 
	{
		for ( iq = 0; iq < 3; iq++ )
		{
			v [ ip ] [ iq ] = 0.0f;
		}
		v [ ip ] [ ip ] = 1.0f;
	}

	for ( ip = 0; ip < 3; ip++ )
	{
		b [ ip ] = a [ ip ] [ ip ];
		d [ ip ] = b [ ip ];
		z [ ip ] = 0.0f;
	}

	for ( i = 1; i <= 50; i++ )
	{
		sm = 0.0f;
		for ( ip = 0; ip < 2; ip++ )
		{
			for ( iq = ip + 1; iq < 3; iq++ )
			{
				sm += (float)fabs ( a [ ip ] [ iq ] );
			}
		}

		if ( sm == 0.0f )
		{
			// sort the stupid things and transpose
			a [ 0 ] [ 0 ] = v [ 0 ] [ 0 ];
			a [ 0 ] [ 1 ] = v [ 1 ] [ 0 ];
			a [ 0 ] [ 2 ] = v [ 2 ] [ 0 ];
			a [ 1 ] [ 0 ] = v [ 0 ] [ 1 ];
			a [ 1 ] [ 1 ] = v [ 1 ] [ 1 ];
			a [ 1 ] [ 2 ] = v [ 2 ] [ 1 ];
			a [ 2 ] [ 0 ] = v [ 0 ] [ 2 ];
			a [ 2 ] [ 1 ] = v [ 1 ] [ 2 ];
			a [ 2 ] [ 2 ] = v [ 2 ] [ 2 ];

			if ( fabs ( d [ 0 ] ) < fabs ( d [ 1 ] ) )
			{
				sm = d [ 0 ];
				d [ 0 ] = d [ 1 ];
				d [ 1 ] = sm;

				sm = a [ 0 ] [ 0 ];
				a [ 0 ] [ 0 ] = a [ 1 ] [ 0 ];
				a [ 1 ] [ 0 ] = sm;
				sm = a [ 0 ] [ 1 ];
				a [ 0 ] [ 1 ] = a [ 1 ] [ 1 ];
				a [ 1 ] [ 1 ] = sm;
				sm = a [ 0 ] [ 2 ];
				a [ 0 ] [ 2 ] = a [ 1 ] [ 2 ];
				a [ 1 ] [ 2 ] = sm;
			}
			if ( fabs ( d [ 1 ] ) < fabs ( d [ 2 ] ) )
			{
				sm = d [ 1 ];
				d [ 1 ] = d [ 2 ];
				d [ 2 ] = sm;

				sm = a [ 1 ] [ 0 ];
				a [ 1] [ 0 ] = a [ 2 ] [ 0 ];
				a [ 2 ] [ 0 ] = sm;
				sm = a [ 1 ] [ 1 ];
				a [ 1 ] [ 1 ] = a [ 2 ] [ 1 ];
				a [ 2 ] [ 1 ] = sm;
				sm = a [ 1 ] [ 2 ];
				a [ 1 ] [ 2 ] = a [ 2 ] [ 2 ];
				a [ 2 ] [ 2 ] = sm;
			}
			if ( fabs ( d [ 0 ] ) < fabs ( d [ 1 ] ) )
			{
				sm = d [ 0 ];
				d [ 0 ] = d [ 1 ];
				d [ 1 ] = sm;

				sm = a [ 0 ] [ 0 ];
				a [ 0 ] [ 0 ] = a [ 1 ] [ 0 ];
				a [ 1 ] [ 0 ] = sm;
				sm = a [ 0 ] [ 1 ];
				a [ 0 ] [ 1 ] = a [ 1 ] [ 1 ];
				a [ 1 ] [ 1 ] = sm;
				sm = a [ 0 ] [ 2 ];
				a [ 0 ] [ 2 ] = a [ 1 ] [ 2 ];
				a [ 1 ] [ 2 ] = sm;
			}

			v [ 0 ] [ 0 ] = a [ 0 ] [ 0 ];
			v [ 0 ] [ 1 ] = a [ 0 ] [ 1 ];
			v [ 0 ] [ 2 ] = a [ 0 ] [ 2 ];
			v [ 1 ] [ 0 ] = a [ 1 ] [ 0 ];
			v [ 1 ] [ 1 ] = a [ 1 ] [ 1 ];
			v [ 1 ] [ 2 ] = a [ 1 ] [ 2 ];
			v [ 2 ] [ 0 ] = a [ 2 ] [ 0 ];
			v [ 2 ] [ 1 ] = a [ 2 ] [ 1 ];
			v [ 2 ] [ 2 ] = a [ 2 ] [ 2 ];

			return;
		}

		if ( i < 4 )
		{
			tresh = 0.2f * sm / 9;
		}
		else
		{
			tresh = 0.0f;
		}

		for ( ip = 0; ip < 2; ip++ )
		{
			for ( iq = ip + 1; iq < 3; iq++ ) 
			{
				g = 100.0f * (float)fabs ( a [ ip ] [ iq ] );
				if ( i > 4 && (float)( fabs ( d [ ip ] ) + g ) == (float)fabs ( d [ ip ] )
					&& (float)( fabs ( d [ iq ] ) + g ) == (float)fabs ( d [ iq ] ) )
				{
					a [ ip ] [ iq ] = 0.0f;
				}
				else
				{
					if ( fabs ( a [ ip ] [ iq ] ) > tresh )
					{
						h = d [ iq ] - d [ ip ];
						if ( (float)( fabs ( h ) + g ) == (float)fabs ( h ) )
						{
							t = ( a [ ip ] [ iq ] ) / h;
						}
						else
						{
							theta = 0.5f * h / ( a [ ip ] [ iq ] );
							t = 1.0f / ( (float)fabs ( theta ) + (float)sqrt ( 1.0f + theta * theta ) );
							if ( theta < 0.0f ) 
							{
								t = -1.0f * t;
							}
						}

						c = 1.0f / (float)sqrt ( 1 + t * t );
						s = t * c;
						tau = s / ( 1.0f + c );
						h = t * a [ ip ] [ iq ];
						z [ ip ] -= h;
						z [ iq ] += h;
						d [ ip ] -= h;
						d [ iq ] += h;
						a [ ip ] [ iq ] = 0.0f;
						for ( j = 0; j <= ip - 1; j++ )
						{
							ROTATE ( a, j, ip, j, iq )
						}
						for ( j = ip + 1; j <= iq - 1; j++ )
						{
							ROTATE ( a, ip, j, j, iq )
						}
						for ( j = iq + 1; j < 3; j++ )
						{
							ROTATE ( a, ip, j, iq, j )
						}
						for ( j = 0; j < 3; j++ )
						{
							ROTATE ( v, j, ip, j, iq )
						}
					}
				}
			}
		}

		for ( ip = 0; ip < 3; ip++ )
		{
			b [ ip ] += z [ ip ];
			d [ ip ] = b [ ip ];
			z [ ip ] = 0.0f;
		}
	}

	v [ 0 ] [ 0 ] = a [ 0 ] [ 0 ];
	v [ 0 ] [ 1 ] = a [ 0 ] [ 1 ];
	v [ 0 ] [ 2 ] = a [ 0 ] [ 2 ];
	v [ 1 ] [ 0 ] = a [ 1 ] [ 0 ];
	v [ 1 ] [ 1 ] = a [ 1 ] [ 1 ];
	v [ 1 ] [ 2 ] = a [ 1 ] [ 2 ];
	v [ 2 ] [ 0 ] = a [ 2 ] [ 0 ];
	v [ 2 ] [ 1 ] = a [ 2 ] [ 1 ];
	v [ 2 ] [ 2 ] = a [ 2 ] [ 2 ];

	printf ( "too many iterations in jacobi\n" );
	//exit ( 1 );
}

int estimateRank ( float *a )
{
	float w [ 3 ];
	float u [ 3 ] [ 3 ];
	float mat [ 3 ] [ 3 ];
	int i;

	mat [ 0 ] [ 0 ] = a [ 0 ] * a [ 0 ];
	mat [ 0 ] [ 1 ] = a [ 1 ] * a [ 0 ];
	mat [ 0 ] [ 2 ] = a [ 2 ] * a [ 0 ];
	mat [ 1 ] [ 0 ] = a [ 1 ] * a [ 0 ];
	mat [ 1 ] [ 1 ] = a [ 1 ] * a [ 1 ] + a [ 4 ] * a [ 4 ];
	mat [ 1 ] [ 2 ] = a [ 1 ] * a [ 2 ] + a [ 4 ] * a [ 5 ];
	mat [ 2 ] [ 0 ] = a [ 2 ] * a [ 0 ];
	mat [ 2 ] [ 1 ] = a [ 1 ] * a [ 2 ] + a [ 5 ] * a [ 4 ];
	mat [ 2 ] [ 2 ] = a [ 2 ] * a [ 2 ] + a [ 5 ] * a [ 5 ] + a [ 7 ] * a [ 7 ];

	jacobi ( mat, w, u );

	if ( w [ 0 ] == 0.0f )
	{
		return 0;
	}
	else
	{
		for ( i = 1; i < 3; i++ )
		{
			if ( w [ i ] < 0.1f )
			{
				return i;
			}
		}

		return 3;
	}
}

void matInverse ( float mat[][3], float midpoint[], float rvalue[][3], float w[], float u[][3] )
{
	// there is an implicit assumption that mat is symmetric and real
	// U and V in the SVD will then be the same matrix whose rows are the eigenvectors of mat
	// W will just be the eigenvalues of mat
//	float w [ 3 ];
//	float u [ 3 ] [ 3 ];
	int i;

	jacobi ( mat, w, u );

	if ( w [ 0 ] == 0.0f )
	{
//		printf ( "error: largest eigenvalue is 0!\n" );
	}
	else
	{
		for ( i = 1; i < 3; i++ )
		{
			if ( w [ i ] < 0.1f ) // / w [ 0 ] < TOLERANCE )
			{
					w [ i ] = 0;
			}
			else
			{
				w [ i ] = 1.0f / w [ i ];
			}
		}
		w [ 0 ] = 1.0f / w [ 0 ];
	}

	rvalue [ 0 ] [ 0 ] = w [ 0 ] * u [ 0 ] [ 0 ] * u [ 0 ] [ 0 ] +
					w [ 1 ] * u [ 1 ] [ 0 ] * u [ 1 ] [ 0 ] +
					w [ 2 ] * u [ 2 ] [ 0 ] * u [ 2 ] [ 0 ];
	rvalue [ 0 ] [ 1 ] = w [ 0 ] * u [ 0 ] [ 0 ] * u [ 0 ] [ 1 ] +
					w [ 1 ] * u [ 1 ] [ 0 ] * u [ 1 ] [ 1 ] +
					w [ 2 ] * u [ 2 ] [ 0 ] * u [ 2 ] [ 1 ];
	rvalue [ 0 ] [ 2 ] = w [ 0 ] * u [ 0 ] [ 0 ] * u [ 0 ] [ 2 ] +
					w [ 1 ] * u [ 1 ] [ 0 ] * u [ 1 ] [ 2 ] +
					w [ 2 ] * u [ 2 ] [ 0 ] * u [ 2 ] [ 2 ];
	rvalue [ 1 ] [ 0 ] = w [ 0 ] * u [ 0 ] [ 1 ] * u [ 0 ] [ 0 ] +
					w [ 1 ] * u [ 1 ] [ 1 ] * u [ 1 ] [ 0 ] +
					w [ 2 ] * u [ 2 ] [ 1 ] * u [ 2 ] [ 0 ];
	rvalue [ 1 ] [ 1 ] = w [ 0 ] * u [ 0 ] [ 1 ] * u [ 0 ] [ 1 ] +
					w [ 1 ] * u [ 1 ] [ 1 ] * u [ 1 ] [ 1 ] +
					w [ 2 ] * u [ 2 ] [ 1 ] * u [ 2 ] [ 1 ];
	rvalue [ 1 ] [ 2 ] = w [ 0 ] * u [ 0 ] [ 1 ] * u [ 0 ] [ 2 ] +
					w [ 1 ] * u [ 1 ] [ 1 ] * u [ 1 ] [ 2 ] +
					w [ 2 ] * u [ 2 ] [ 1 ] * u [ 2 ] [ 2 ];
	rvalue [ 2 ] [ 0 ] = w [ 0 ] * u [ 0 ] [ 2 ] * u [ 0 ] [ 0 ] +
					w [ 1 ] * u [ 1 ] [ 2 ] * u [ 1 ] [ 0 ] +
					w [ 2 ] * u [ 2 ] [ 2 ] * u [ 2 ] [ 0 ];
	rvalue [ 2 ] [ 1 ] = w [ 0 ] * u [ 0 ] [ 2 ] * u [ 0 ] [ 1 ] +
					w [ 1 ] * u [ 1 ] [ 2 ] * u [ 1 ] [ 1 ] +
					w [ 2 ] * u [ 2 ] [ 2 ] * u [ 2 ] [ 1 ];
	rvalue [ 2 ] [ 2 ] = w [ 0 ] * u [ 0 ] [ 2 ] * u [ 0 ] [ 2 ] +
					w [ 1 ] * u [ 1 ] [ 2 ] * u [ 1 ] [ 2 ] +
					w [ 2 ] * u [ 2 ] [ 2 ] * u [ 2 ] [ 2 ];
}

float calcError ( float mat[], float point[] )
{
	float rvalue = mat [ 9 ] * mat [ 9 ];
	float tmp;

	tmp = mat [ 0 ] * point [ 0 ] + mat [ 1 ] * point [ 1 ] + mat [ 2 ] * point [ 2 ] - mat [ 3 ];
	rvalue += tmp * tmp;

	tmp = mat [ 4 ] * point [ 1 ] + mat [ 5 ] * point [ 2 ] - mat [ 6 ];
	rvalue += tmp * tmp;

	tmp = mat [ 7 ] * point [ 2 ] - mat [ 8 ];
	rvalue += tmp * tmp;

	return rvalue;
}

float calcPoint ( float midpoint[], float rvalue[], float *mat )
{
	float newB [ 3 ];
	float a [ 3 ] [ 3 ];
	float inv [ 3 ] [ 3 ];
	float w [ 3 ];
	float u [ 3 ] [ 3 ];
	float b [ 3 ];

	a [ 0 ] [ 0 ] = mat [ 0 ] * mat [ 0 ];
	a [ 0 ] [ 1 ] = mat [ 1 ] * mat [ 0 ];
	a [ 0 ] [ 2 ] = mat [ 2 ] * mat [ 0 ];
	a [ 1 ] [ 0 ] = mat [ 1 ] * mat [ 0 ];
	a [ 1 ] [ 1 ] = mat [ 1 ] * mat [ 1 ] + mat [ 4 ] * mat [ 4 ];
	a [ 1 ] [ 2 ] = mat [ 1 ] * mat [ 2 ] + mat [ 4 ] * mat [ 5 ];
	a [ 2 ] [ 0 ] = mat [ 2 ] * mat [ 0 ];
	a [ 2 ] [ 1 ] = mat [ 1 ] * mat [ 2 ] + mat [ 5 ] * mat [ 4 ];
	a [ 2 ] [ 2 ] = mat [ 2 ] * mat [ 2 ] + mat [ 5 ] * mat [ 5 ] + mat [ 7 ] * mat [ 7 ];

	b [ 0 ] = mat [ 0 ] * mat [ 3 ];
	b [ 1 ] = mat [ 1 ] * mat [ 3 ] + mat [ 4 ] * mat [ 6 ];
	b [ 2 ] = mat [ 2 ] * mat [ 3 ] + mat [ 5 ] * mat [ 6 ] + mat [ 7 ] * mat [ 8 ];


	matInverse ( a, midpoint, inv, w, u );

	newB [ 0 ] = b [ 0 ] - a [ 0 ] [ 0 ] * midpoint [ 0 ] - a [ 0 ] [ 1 ] * midpoint [ 1 ] - a [ 0 ] [ 2 ] * midpoint [ 2 ];
	newB [ 1 ] = b [ 1 ] - a [ 1 ] [ 0 ] * midpoint [ 0 ] - a [ 1 ] [ 1 ] * midpoint [ 1 ] - a [ 1 ] [ 2 ] * midpoint [ 2 ];
	newB [ 2 ] = b [ 2 ] - a [ 2 ] [ 0 ] * midpoint [ 0 ] - a [ 2 ] [ 1 ] * midpoint [ 1 ] - a [ 2 ] [ 2 ] * midpoint [ 2 ];

	rvalue [ 0 ] = inv [ 0 ] [ 0 ] * newB [ 0 ] + inv [ 1 ] [ 0 ] * newB [ 1 ] + inv [ 2 ] [ 0 ] * newB [ 2 ] + midpoint [ 0 ];
	rvalue [ 1 ] = inv [ 0 ] [ 1 ] * newB [ 0 ] + inv [ 1 ] [ 1 ] * newB [ 1 ] + inv [ 2 ] [ 1 ] * newB [ 2 ] + midpoint [ 1 ];
	rvalue [ 2 ] = inv [ 0 ] [ 2 ] * newB [ 0 ] + inv [ 1 ] [ 2 ] * newB [ 1 ] + inv [ 2 ] [ 2 ] * newB [ 2 ] + midpoint [ 2 ];

/*	switch ( Octree::method )
	{
	case 0:
		matInverse ( a, midpoint, inv, w, u );

		newB [ 0 ] = b [ 0 ] - a [ 0 ] [ 0 ] * midpoint [ 0 ] - a [ 0 ] [ 1 ] * midpoint [ 1 ] - a [ 0 ] [ 2 ] * midpoint [ 2 ];
		newB [ 1 ] = b [ 1 ] - a [ 1 ] [ 0 ] * midpoint [ 0 ] - a [ 1 ] [ 1 ] * midpoint [ 1 ] - a [ 1 ] [ 2 ] * midpoint [ 2 ];
		newB [ 2 ] = b [ 2 ] - a [ 2 ] [ 0 ] * midpoint [ 0 ] - a [ 2 ] [ 1 ] * midpoint [ 1 ] - a [ 2 ] [ 2 ] * midpoint [ 2 ];

		rvalue [ 0 ] = inv [ 0 ] [ 0 ] * newB [ 0 ] + inv [ 1 ] [ 0 ] * newB [ 1 ] + inv [ 2 ] [ 0 ] * newB [ 2 ] + midpoint [ 0 ];
		rvalue [ 1 ] = inv [ 0 ] [ 1 ] * newB [ 0 ] + inv [ 1 ] [ 1 ] * newB [ 1 ] + inv [ 2 ] [ 1 ] * newB [ 2 ] + midpoint [ 1 ];
		rvalue [ 2 ] = inv [ 0 ] [ 2 ] * newB [ 0 ] + inv [ 1 ] [ 2 ] * newB [ 1 ] + inv [ 2 ] [ 2 ] * newB [ 2 ] + midpoint [ 2 ];
		break;
	case 1:
		rvalue [ 0 ] = midpoint [ 0 ];
		rvalue [ 1 ] = midpoint [ 1 ];
		rvalue [ 2 ] = midpoint [ 2 ];
	}
*/
	return calcError ( mat, rvalue );
}
