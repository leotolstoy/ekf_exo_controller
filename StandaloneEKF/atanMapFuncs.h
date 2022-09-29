#include <math.h>

double MAXIMUM_STRIDELENGTH = 4;

double arctanMap(double x)
{
     return (MAXIMUM_STRIDELENGTH/M_PI)*(atan( (M_PI/(MAXIMUM_STRIDELENGTH)) * (x ) )) + MAXIMUM_STRIDELENGTH/2;
}

double dAarctanMap(double x)
{
     return (1/(1 + (((M_PI/MAXIMUM_STRIDELENGTH)*x)*((M_PI/MAXIMUM_STRIDELENGTH)*x)) ));

}

double invArctanMap(double y)
{
     return ((MAXIMUM_STRIDELENGTH)/M_PI)*(tan( (M_PI/MAXIMUM_STRIDELENGTH)* (y-MAXIMUM_STRIDELENGTH/2)));
}

	

