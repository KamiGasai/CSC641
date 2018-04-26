/////////////////////////////////////////////////////////////////////
//  statlib.cpp - Elementary Statistics Library                    //
//                                                                 //
//  1. Random number generators                                    // 
//  2. Basic statistics                                            //
//  3. Correlation                                                 //
//  4. Regression                                                  //
//  5. Visualization                                               //
//                                                                 //
//  Jozo Dujmovic 9/29/2006                                        //
/////////////////////////////////////////////////////////////////////

#include<iostream.h>
#include<math.h>
#include<iomanip.h>
#include<string.h>
#include<stdlib.h>




/////////////////////////////////////////////////////////////////////
//  1. Random number generators                                    // 
/////////////////////////////////////////////////////////////////////

//---------------------------------------------------------------
//  Machine-dependent standard uniform random number generator 
//  from the standard library
//---------------------------------------------------------------
double urn(void)   
{
	return double(rand())/double(RAND_MAX);
}


//------------------------------------------------------
//  Fast Fibonacci-style machine-independent standard
//  uniform random number generator. The quality of 
//  randomness is not tested.
//  Jozo Dujmovic, Dec 2002
//------------------------------------------------------
double FiboURN(void)
{
	static double p=sqrt(2.)-1., q=sqrt(31.)-5.; 
	static int flag = 1;

	if(flag)
	{
		p += q;
		if(p>1.) p -= 1.;
		flag = 0;
		return p;
	}
	else
	{
		q += p;
		if(q>1.) q -= 1.;		
		flag = 1;
		return q;
	}
}


//------------------------------------------------------------------
// Machine-independent additive generator of standard uniform random 
// numbers. Good quality of randomness.
//                 r(i) := (r(i-1) + r(i-17)) mod xmod
// Jozo Dujmovic, 1998
//------------------------------------------------------------------
double uniform(void)
{
      static unsigned long r[18]={131071,43691,262657,649657,274177,
                                  524287,121369,61681,179951,513239,
                                  333667,909091,1777,8617,87211,
                                  174763,65537,0},
                           xmod = 1048573;
      static double rnmax=xmod;
      static int i=17,j=16,k=0,n=18;
      double rn;

      r[i] = (r[j] + r[k]) % xmod;
      rn   = r[i]/rnmax;
      i    = (i+1)%n;
      j    = (j+1)%n;
      k    = (k+1)%n;
      return rn;
 }



//-----------------------------------------------------------------------
//     MacLaren - Marsaglia Shuffler for a random number generator
//     definied as urn( ).
//     Jozo Dujmovic, 1998
//-----------------------------------------------------------------------
double rng(void)
{
      static int n=200, firstcall=1;
      static double table[200];
      double rnumber;
      int i, itable;

      if (firstcall)
      {
         for(i=0; i<n; i++) table[i]=urn( ); // Any desired distribution
         firstcall = 0;
      }

      itable        = int(n * urn( ));       // Uniform selection from table
      rnumber       = table[itable];
      table[itable] = urn( );                // Any desired distribution
      return rnumber;
}


//----------------------------------------------------------------
//  Standard uniform random number selector function
//----------------------------------------------------------------
double URN(int selector)
{
	if(selector == 1)  return urn();       // Standard Library
	if(selector == 2)  return FiboURN();   // Fibonacci (2 values)
	if(selector == 3)  return uniform();   // Fibonacci (array)
	if(selector == 4)  return rng();       // Shuffler
	return urn();
}



//------------------------------------------------------------------
// Generator of a pair of random numbers (x and y) that have
// a desired degree of (positive) linear correlation r
// Conditions: O <= r <= 1,  O <= x <= 1,  O <= y <= 1
//------------------------------------------------------------------
void CorrelatedPair(double& x, double& y, double r)
// x = standard random number
// y = correlated random number computed from x 
// r = parameter similar to desired correlation coefficient
//
// Method: Generate random x. Generate random y as x + noise.
//         Noise is a small random value denoted as delta:
//         delta = (1-r^0.72)(URN - 1/2),    -0.5 < delta < 0.5
//         y = x + delta         if 0 < x+delta < 1
//           = x + delta -1      if x+delta > 1
//           = x + delta +1      if x+delta < 0
//         
// Jozo Dujmovic, June 2002 (Foz de Iguazu conference)
//------------------------------------------------------------------
{
	x = URN(1);         // R = resulting correlation between x and y
	if(r < 0.)        y = 1.-x ;                             // R=-1
	else if (r == 0.) y = URN(1);                            // R= 0
	else if (r < 1.)  y = x + (1.-pow(r,0.72))*(URN(1)-0.5); // R= r
	else              y = x ;                                // R= 1

	if(y < 0.)        y += 1.;
	else if(y > 1.)   y -= 1.;
}








/////////////////////////////////////////////////////////////////////
//  2. Basic statistics                                            //
/////////////////////////////////////////////////////////////////////

// Mean of a sequence of n values
double mean(double x[], int n)     
{
  double sx=0.;
  for(int i=0; i<n; i++) sx += x[i];
  return sx/n;
}


// Standard deviation of a sequence of n values
double sigma(double x[], int n) 
{
  double sx=0., sxx=0.;
  for(int i=0; i<n; i++)
  {
    sx  += x[i];
    sxx += x[i]*x[i];
  }
  sx  /= n;
  sxx /= n;
  return sqrt(sxx - sx*sx);
}


// Mean and standard deviation of a sequence of n values
void meansigma(double x[], int n, double& mean, double& sigma)
{
  double sx=0., sxx=0.;
  for(int i=0; i<n; i++)
  {
    sx  += x[i];
    sxx += x[i]*x[i];
  }
  mean = sx / n;
  sxx /= n;
  sigma = sqrt(sxx - mean*mean);
}




/////////////////////////////////////////////////////////////////////
//  3. Correlation                                                 //
/////////////////////////////////////////////////////////////////////

// Coefficient of correlation between sequences x[] and y[]
double r(double x[], double y[], int n)  
{
  double sxy=0.;
  for(int i=0; i<n; i++) sxy += x[i]*y[i];
  sxy /= n;
  return (sxy - mean(x,n)*mean(y,n))/(sigma(x,n)*sigma(y,n));
}

void show(double x[], int n)
{
  for (int i = 0; i<n; i++) cout << x[i] << ' ';
  cout << '\n';
}

void sort(double a[], int n)
{ double t;
  int i,j;
  for(i=0; i<n-1; i++)
	  for(j=i+1; j<n; j++)
		  if(a[i] > a[j]) {t=a[i]; a[i]=a[j]; a[j]=t;}
}

double sum(double a[], int from, int how_many)
{
	double s=0;
	for(int i=0; i<how_many; i++) s += a[from++];
	return s;
}

double mean(double a[], int from, int how_many)
{
	double s=0;
	for(int i=0; i<how_many; i++) s += a[from++];
	return s/how_many;
}

//---------------------------------------------------------------
// Achieved correlation as a function of desired correlation
//---------------------------------------------------------------
void CorreTestMean(void) 
{
	const int N=100, K=1000;
	double x[N], y[N], R, Rave, Erel, meanErel, Eabs, meanEabs;
	int i,k, cnt;

    cout.setf(ios::fixed, ios::floatfield);
    cout.setf(ios::showpoint);
    cout << '\n' << setw(10) << "Desired r" << setw(15) << "Achieved Rave" 
		         << setw(16) << "Abs error" 
				 << setw(17) << "Rel Error\n";
	meanErel = 0.; cnt = 0; meanEabs=0.;
	for(R=0.; R<=1.00001; R+=0.05)
	{
		Rave = 0.; 
		for(k=0; k<K; k++)
		{
			for(i=0; i<N; i++) CorrelatedPair(x[i],y[i],R);
			Rave += r(x,y,N);
		}
		Rave /= K;

		Eabs = fabs(Rave-R);
        meanEabs += Eabs;

		if(R>0.)
		{
			Erel = 100.*fabs(Rave-R)/R;
			meanErel += Erel; cnt++;

			cout << setw(10) << R << setw(15) << Rave << setw(15) 
				 << 100.*Eabs <<'%' << setw(15) << Erel << "%";
			if(R<0.9999) cout << '\n';
			else cout << "\nAverage Eabs = " << 100.*meanEabs/(cnt+1)
			     	  <<"%,   Average Erel = " << meanErel/cnt << "%\n";
		}
		else
			cout << setw(10) << R << setw(15) << Rave << setw(15) << 100.*Eabs <<"%\n";
	}
        
}




/////////////////////////////////////////////////////////////////////
//  4. Regression                                                  //
/////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////
// Linear Regression and Correlation
// Finds a and b for the regression line y = ax+b and the correlation
// coefficient r
// Jozo Dujmovic  3/3/2001
///////////////////////////////////////////////////////////////////////
void lire(double x[],double y[],int nxy,double& a,double& b,double& r)
{ double sx, sy, sxx, syy, sxy;
  sx=sy=sxx=syy=sxy=0.;
  for(int i=0; i<nxy; i++)
  {  sx  += x[i];
     sy  += y[i];
     sxx += x[i]*x[i];
     syy += y[i]*y[i];
     sxy += x[i]*y[i];
  }
  a = (nxy*sxy - sx*sy)/(nxy*sxx - sx*sx);
  b = (sy - a*sx)/nxy;
  r = (nxy*sxy-sx*sy)/(sqrt(nxy*sxx-sx*sx)*sqrt(nxy*syy-sy*sy));
}

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  FMIN1  - Minimum of function F(X)                              //
//  Jozo Dujmovic, 11/18/2001                                      //
//                                                                 //
//    F     = the function to be minimized                         //
//    STEP  = initial value of step                                //
//    ERROR = maximum absolute error of XMIN                       //
//    XMIN  = an initial estimate of local minimum which is        //
//            replaced by the resulting coordinate of minimum      //
//    FMIN  = F(XMIN) = the resulting minimum of function F        //
/////////////////////////////////////////////////////////////////////
                                                                   //
void FMIN1(double (*F)(double), double STEP,double ERROR,          //
           double& XMIN, double& FMIN)                             //
{ double XNEXT, FNEXT;                                             //
      FMIN = F(XMIN);                                              //
      while (fabs(STEP) >= ERROR*0.1353)                           //
      {  STEP  = -STEP;                                            //
         XNEXT = XMIN + STEP;                                      //
         FNEXT = F(XNEXT);                                         //
         while (FNEXT < FMIN)                                      //
         {  XMIN  = XNEXT;                                         //
            FMIN  = FNEXT;                                         //
            XNEXT = XNEXT + STEP;                                  //
            FNEXT = F(XNEXT);                                      //
         }                                                         //
         STEP = STEP/2.7183;                                       //
      }                                                            //
}                                                                  //
/////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////
//  5. Visualization                                               //
/////////////////////////////////////////////////////////////////////

//------------------
// Set parecision
//------------------
void setprec(int n)  
{ 
   cout << setiosflags(ios::fixed | ios::showpoint) << setprecision(n);
}

//----------------------------
// Show an array of numbers
//----------------------------
void ShowArray(double x[], int n)
{
  setprec(4);
  for (int i = 0; i<n; i++) 
	  cout << x[i] << ((i+1)%10 ? ' ' : '\n');
  cout << '\n';
}
void ShowArray(int x[], int n)
{
  cout << '\n';
  for (int i = 0; i<n; i++) 
	  cout << setw(6) << x[i] << ((i+1)%10 ? ' ' : '\n');
  cout << '\n';
}


/////////////////////////////////////////////////////////////////////
// Print-plot (x[i],y[i])-pairs for i=0,...,n-1.  x>0 and y>0
// xname, yname =  names of x and y axis
// nx, ny = picture size = number characters per x and y axis 
// Suggested values: nx = multiple of 10 (e.g. nx=50), 20<=ny=100
// Jozo Dujmovic, 12/20/2002
/////////////////////////////////////////////////////////////////////
void xyplot(double x[], double y[], int n, char xname[], char yname[], 
			int nx, int ny)
{
	double xmax, ymax, X1,X2,Y1,Y2;
	int i,j,k;
	char c;
	cout << "\n\n         " << yname << setprecision(2);
    xmax=x[0]; ymax=y[0];
	for(i=1; i<n; i++)
	{
		if(x[i]>xmax) xmax=x[i];
		if(y[i]>ymax) ymax=y[i];
	}
    for(j=ny; j>0; j--)                  // y-axis loop
	{
		Y1 = (j-1)*(ymax/ny); Y2 = j*(ymax/ny); 
		cout << '\n' << setw(8) << Y2 << " |";
		for(i=1; i<=nx; i++)             // x-axis loop
		{
			X1 = (i-1)*(xmax/nx); X2 = i*(xmax/nx); 
			for(c=' ', k=0; k<n; k++) 
				if(X1<x[k] && x[k]<=X2 && Y1<y[k] && y[k]<=Y2) c='*';
			cout << c;
		}
	}
	cout << '\n' << setw(8) << 0 << " +";
    for(i=1; i<=nx; i++) cout << (i%10 ? '-' : '+'); 
	cout << '\n' << setw(10) << 0 << "  ";
	for(i=1; i<=nx/10; i++) cout << setw(10) << 10.*i*xmax/nx;
	cout << "\n";
	for(i=0; i<int(10+nx-strlen(xname)); i++) 
		cout << ' '; cout<< xname << "\n\n";
}


//---------------------------------------------------------------
//  Frequency distribution
//  x[0 .. nx-1]        = array of analyzed values
//  number_of_intervals = intervals for distribution plot <= 100
//  xmin, xmax          = range of x values
//  if xmin<xmax then the program uses the given values; otherwise
//  the program finds minimum and maximum in the array x[]
//---------------------------------------------------------------
void distribution(double x[], int nx, double xmin, double xmax, int number_of_intervals)
{
	int i,j, F[100], Fmax=0, Fsum=0,  // F[]=frequency distribution
		FMAX=50;                      // Maximum number of asterisks
	                                  // displayed in the histogram
	if(xmin >= xmax)
	{
		xmin=xmax=x[0];
		for(i=1; i<nx; i++)               // Find the range of values
		{                                 // (xmin <= x[i] <= xmax)
			if(x[i]>xmax)      xmax=x[i];
			else if(x[i]<xmin) xmin=x[i];
		}
	}

	cout << "\nFREQUENCY DISTRIBUTION\n";
	cout << "\n     xmin = " << xmin << "    xmax = " << xmax ;
                                          // Compute the frequencies
	if(number_of_intervals > 100) number_of_intervals = 100 ;
	for(i=0; i < number_of_intervals; i++) F[i]=0;
    for(i=0; i<nx; i++)
	   ++F[int(0.9999*number_of_intervals*(x[i]-xmin)/(xmax-xmin))];

	ShowArray(F,number_of_intervals);

    for(i=0; i<number_of_intervals; i++)  // Find the maximum frequency
	{									  // and limit the maximum to
		Fsum += F[i];                     // the desired maximum FMAX
		if(F[i]>Fmax) Fmax=F[i];    
	}
	cout << "Total = " << Fsum << "   Max frequency = " << Fmax << endl;

	if(Fmax>FMAX)                         
	{	for(i=0; i<number_of_intervals; i++)
			F[i] = int(0.499+(FMAX*F[i])/Fmax);
	    Fmax = FMAX;
	}
	
	setprec(3);
	cout << "\n" << setw(10) << xmin << "+";
	for(i=0; i<(1+(Fmax>50 ? 10 : Fmax/5)); i++) cout << "----+";
	
	for(i=0; i<number_of_intervals; i++)
	{
		cout << '\n' << setw(10) 
			 << (xmin + (i+1.)*(xmax-xmin)/number_of_intervals) << '|';
		for(j=1; j<=F[i]; j++) cout << '*'; 
		cout << "  " << F[i];
	}

	cout << "\n          +";
	for(i=0; i<(1+(Fmax>50 ? 10 : Fmax/5)); i++) cout << "----+";	
    cout << "\n\n";
}


void main(void)
{

  double f[10]={0,0,0,0,0,0,0,0,0,0};
  long int m, i, n=0, k=10000;

  setprec(3);
  cout << "\n\nTest of uniform distribution of the rng generator\n";
  while(n<90000)
  {
         for(i=0; i<k; i++)
            f[int(10*rng())] += 1.;
         n += k;
         cout<< '\n'<< n << ' ';
         for(i=0; i<10; i++) cout << f[i]/n << ' ';
  }

  cout << "\n\nTest of correlation between arrays of numbers\n";
  double x[10000] = {0,1,2,3,4,5,6,7,8,9},
         y[10000] = {9,8,7,6,5,4,3,2,1,0},
         z[10000] = {1,2,3,4,5,6,7,8,9,10},
         a[10000] = {0,2,1,3,5,4,7,6,10,9};
  cout << "\nx = "; show(x,10);
  cout << "y = "; show(y,10);
  cout << "z = "; show(z,10);
  cout << "a = "; show(a,10);
  cout << "\nmean  x = " << mean(x,10)
       << "\nsigma x = " << sigma(x,10)
       << "\nr(x,x)  = " << r(x,x,10)
       << "\nr(x,y)  = " << r(x,y,10)
       << "\nr(x,z)  = " << r(x,z,10)
       << "\nr(x,a)  = " << r(x,a,10)
       << "\nr(y,a)  = " << r(y,a,10)
       << "\nr(z,a)  = " << r(z,a,10)
       << "\n\n";


  cout << "\nTest of the generator of 1000 correlated numbers\n";
  for(i=0; i<1000; i++) CorrelatedPair(x[i], y[i], -1.);
  cout << "\nDesired r = -1      Achieved R = " << r(x,y,1000);
  for(i=0; i<1000; i++) CorrelatedPair(x[i], y[i], 0.);
  cout << "\nDesired r =  0      Achieved R = " << r(x,y,1000);
  for(i=0; i<1000; i++) CorrelatedPair(x[i], y[i], 1.);
  cout << "\nDesired r =  1      Achieved R = " << r(x,y,1000);
  for(i=0; i<1000; i++) CorrelatedPair(x[i], y[i], 0.5);
  cout << "\nDesired r =  0.5    Achieved R = " << r(x,y,1000);

  cout << "\n\nSample graph of correlated x,y values\n";
  xyplot(x,y,200,"x", "y", 50, 40);


  cout << "\n\nAverage desired and achieved correlations\n";
  CorreTestMean();

  
  n=1000;
  for(m=1; m<=4; m++)
  {
     cout << "\n\nAn array of " << n << " uniform random numbers\n"
	   << "Type of distribution = " << m << endl;
     for(i=0; i<n; i++) x[i]= URN(m);
  // ShowArray(x,n);
     distribution(x,n, 0., 1., 10);

  }
  cout << endl;
}