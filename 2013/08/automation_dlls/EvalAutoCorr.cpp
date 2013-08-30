/* Read and evaluate FAST P7888 *.lst file */
/* Originally by Wolfgang
 * Modified and commented by Hiroki
 * v1 25/01/2011
 * v2 25/04/2012
 * 
 * Re-made to be an autocorrelation evaluation function by Nic
 * v1 27/08/2013
 */

#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <exception>

//#define DEBUGGING

#include <windows.h>   // for messagebox
#include <sstream>

using namespace std;

/* Function declarations */
int EvalAutoCorr(unsigned int CorrelationHistogram[], unsigned long long gtime[], 
		  unsigned int total_counts, unsigned int HistogramLength, int taumin, int taumax, double *t_res);

int EvalAutoCorr(unsigned long CorrelationHistogram[], unsigned long long gtime[], 
		  unsigned long total_counts, unsigned int HistogramLength, int taumin, int taumax, double *t_res)
{
	
  //Make an auto-correlation histogram
  long tau;	//time difference
  double taures;	//time resolution of the histogram
  taures = (double) (taumax-taumin) / HistogramLength;
  *t_res = taures;
  
  #ifdef DEBUGGING
  stringstream msg;
  msg.str("");
  msg << "total counts = " << total_counts <<endl
	  << "taumax = " << taumax
      << "\ntaumin = " << taumin 
	  << "\ntaures = " << taures 
	  << "\nHistogram length = " << HistogramLength << endl;
  MessageBox(NULL, msg.str().c_str(), "EvalAutoCorr", MB_OK);
  #endif
	
  //initialize histogram
  for (unsigned int i = 0; i < HistogramLength; i++) CorrelationHistogram[i] = 0;
	
  unsigned long j1, j2, j2end, lim;
  long i[5];
  lim = total_counts;
  j2 = 0;
  try
  {
	  for (j1=0; j1 < total_counts; j1++) {
		i[0] = (long) j1;
		if (j1 > lim){
			throw i;
		}
		//Given j1-th gtime, find all events 
		//within gtime+taumin < gtime < gtime+taumax (taumin < 0)
			
		//First find the left most events
		while( ( (long) (gtime[j2] - gtime[j1]) < taumin) && (j2 < total_counts-1) ) 
		{
				j2++;
				i[1] = (long) j2;
				if (j2 > lim){
					throw i;
				}
		}
		j2end = j2;
				
		//Calculate time difference until it becomes larger than taumax
		while( ( (long) (gtime[j2end]-gtime[j1]) < taumax) && (j2end < total_counts-1) ) 
		{
			tau = (long) (gtime[j2end]-gtime[j1]);
			i[4] = tau;
			j2end++;
			i[2] = (long) j2end;
			if (j2end > lim){
				throw i;
			}
			i[3] = (long) ((tau - taumin)/taures);
			if ((unsigned int) i[3] > HistogramLength){
				throw i;
			}
			if(i[3] < HistogramLength && j2end != j1){
				CorrelationHistogram[(unsigned int) ((tau - taumin)/taures) ]++;
			}
		}
	  }
  }
  catch(long i[])
  {
	stringstream msg;
	msg.str("");
	msg << "Exception: j1, j2, jend, m, tau_last: " << i[0] << ", " << i[1] 
		<< ", " << i[2] << ", " << i[3] << ", " << i[4]
	    << "\ntotal_counts = " << total_counts << ", HistogramLength = "<< HistogramLength << endl;
	MessageBox(NULL, msg.str().c_str(), "EvalAutoCorr", MB_OK); 
  }
  
  #ifdef DEBUGGING
  msg.str("");
  msg << "Succesfully calculated autocorrelation!" << endl;
  MessageBox(NULL, msg.str().c_str(), "EvalCrossCorr", MB_OK);
  #endif	
	
  return(0);
}
