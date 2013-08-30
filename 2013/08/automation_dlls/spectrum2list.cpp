/* Read and conver FAST P7888 spectrum data to list of timestamps */
/* N. Seymour-Smith 29/08/2013 */

#include <cmath>

/* Function declarations */
int sum_spectrum(unsigned int Spectrum[], unsigned int SpectrumLength, unsigned int * TotalCounts);

int spectrum2list(unsigned int Spectrum[], unsigned int SpectrumLength, unsigned long long List[]);

int sum_spectrum(unsigned int Spectrum[], unsigned int SpectrumLength, unsigned int * TotalCounts)
{
	int i, sum = 0;
	for (i = 0; i < SpectrumLength; i++){
		sum += Spectrum[i];
	}
	*TotalCounts = sum;
}

int spectrum2list(unsigned int Spectrum[], unsigned int SpectrumLength, unsigned long long List[])
{
	unsigned long long i = 0;
	int j, k = 0;
	for (i = 0; i < SpectrumLength; i++){
		if (Spectrum[i]){
			for(j = 0; j < Spectrum[i]; j++){
				List[k] = i+1;
				k++;
			}
		}
	}
}
