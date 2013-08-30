#ifndef _DLL_H_
#define _DLL_H_

#if BUILDING_DLL
# define DLLIMPORT __declspec (dllexport)
#else /* Not BUILDING_DLL */
# define DLLIMPORT __declspec (dllimport)
#endif /* Not BUILDING_DLL */

DLLIMPORT long ReadP7888(char filepath[], unsigned long CorrelationHistogram[], 
	unsigned long HistogramLength, long taumin, long taumax, double Frequency[]);

#endif /* _DLL_H_ */
