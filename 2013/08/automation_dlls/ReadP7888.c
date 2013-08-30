/* Read and evaluate FAST P7888 *.lst file */

#include <fstream>
#include <string>
#include <math.h>

// #define DEBUGGING

#ifdef DEBUGGING
       #include <windows.h>   // for messagebox
       #include <sstream>
#endif

#include "ReadP7888.h"


using namespace std;

void four1(double data[], unsigned long nn);

long ReadP7888(char filepath[], unsigned long CorrelationHistogram[], 
	unsigned long HistogramLength, long taumin, long taumax, double Frequency[])
	{
    
    int k;		
    errno = 0;
      ifstream  binfile(filepath,  ios::binary); // ios::nocreate |
    if (errno>0)
    {  // perror(filepath);
       return(errno);
    }
	  
    unsigned long filelength;
    
    // get length of file:
    binfile.seekg (0, ios::end);
    filelength = binfile.tellg();
    binfile.seekg (0, ios::beg);

    string newline;
    string Body = "";
    string test;
    char cline[100] ;
    bool neof;

    bool ascii;
    long range;
    long headerpos;
    unsigned long datasize;
		
    ascii = false;
    do
    {
      neof=binfile.getline(cline,100);
      newline = cline;
      test = strtok(cline,"=");
      if (test.compare("fmt")==0)
      {    test = strtok(NULL,"");
           if (test.compare("asc")==0) ascii = true;
      }
      else
      if (test.compare("range")==0)
      {    test = strtok(NULL,"");
           range = atol(test.c_str());
      }
      Body += newline;
    }
    while (neof & (newline.compare(0,6,"[DATA]")!=0));

#ifdef DEBUGGING
    MessageBox(NULL, Body.c_str(), "File header", MB_OK);
#endif

    headerpos = binfile.tellg();
    datasize  = (filelength-headerpos)/4;

  // unsigned long maxdatasize;
  // if ((datasize>maxdatasize) && (maxdatasize>0)) datasize = maxdatasize;

    unsigned char *channel;
    unsigned long *time;
    unsigned __int64 *gtime;

    channel = new unsigned char[datasize];
    time    = new unsigned long[datasize];
    gtime   = new unsigned __int64[datasize];

    /* read actual data */

    unsigned long i = 0;

    if (ascii)
    {
        do
        {
              neof = binfile.getline(cline,100);
              time[i] = atol(cline);
              i++;
        }
        while (neof);
        datasize = i;
   }
   else
   {
        binfile.read ((char*)time, sizeof (unsigned long) * datasize);
   }

   errno = 0;
   binfile.close( );
   if (errno>0)
   { // perror(filepath);
     return(errno);
   }

  // extract time stamp (in ns) from data
  // the highest 32th bit represents the channel

    int chwidth = 1;    // channel width in ns = 2 for 4 channel readout
    unsigned long mask = 0x7fffffff;    // 0x3fffffff; for 4 channel readout
    unsigned long  timemin, timemax;
    timemin = timemax  = ( time[0] & mask ) * chwidth  ;

    unsigned int nch[2] = {0};
    int lap = 0;
    unsigned long lasttime = 0;

    unsigned __int64 timerange =  2147483648ull  * chwidth;  // ull = ui64, 
                                                             // 1073741824ull for 4 channel readout
     
    for (unsigned long i=0; i<datasize; i++)
    {
        channel[i]  = time[i] >> 31;    // 30 for 4 channel readout
        time[i]     = ( time[i] & mask ) * chwidth ;
        if (time[i] > timemax) { timemax = time[i]; };
        if (time[i] < timemin) { timemin = time[i]; };
        if (time[i] < lasttime) lap++;   
        if (time[i] > timerange)  
        {                     
#ifdef DEBUGGING
           std::stringstream errmessage;
           errmessage << "time[" << i << "] = " << time[i] << 
                         "\nlap = " << lap << 
                         "\ntimemax = " << timemax << endl;
           MessageBox(NULL, errmessage.str().c_str(), "Error", MB_OK);
#endif        
           return(1001);
        }
        gtime[i] = timerange * lap + time[i];  
        lasttime = time[i];
        nch[channel[i]]++;
    }

    // store indices of channels separately
    unsigned int *chnr[2];
    unsigned int ch_index[2] = {0};
    for( k=0; k<2; k++)  chnr[k] = new unsigned int[nch[k]]; 

    for (unsigned long i=0; i<datasize; i++)
    {
        k = channel[i] ; 
        chnr[k][ch_index[k]++] = i;
    }
 
    if (timemax-64*chwidth > timerange)
    {
 #ifdef DEBUGGING
       std::stringstream errmessage;
       errmessage << "Timing error:\n timemax = " << timemax << 
                     "\n timerange = " << timerange << endl;
       MessageBox(NULL, errmessage.str().c_str(), "Error", MB_OK);
#endif
       return(1002);
    }

#ifdef DEBUGGING
       float  totaltime = gtime[datasize-1] * 1E-9;   // final time in seconds 
    double rate[2];
    for (int k=0; k<2; k++)
       rate[k] = (double) nch[k] / totaltime ;

    {  std::stringstream msg;
       msg << "Total time = " << totaltime << "s" << 
              "\nRate channel 0 =" << rate[0] << "/s" <<
              "\nRate channel 1 =" << rate[1] << "/s" <<           
              endl;
       MessageBox(NULL, msg.str().c_str(), "Info", MB_OK);
    }
#endif

    long d;
    double taures;

    taures = (double) (taumax-taumin) / HistogramLength;

#ifdef DEBUGGING
    {
      std::stringstream msg;
      msg << "taures = " << taures << "ns" << 
          "\ntaumin =" << taumin << "ns" <<
          "\ntaumax =" << taumax << "ns" <<           
          endl;
      MessageBox(NULL, msg.str().c_str(), "Info", MB_OK);
    }
#endif
  
    
//    int *hist;      
//    hist = new int[HistogramLength];   
//    for (unsigned int i=0; i<HistogramLength; i++) hist[i]=0;

    for (unsigned int i=0; i<HistogramLength; i++) CorrelationHistogram[i]=0;

#ifdef DEBUGGING
    int check=0;
#endif    

    unsigned long i1, i2, i2end, k2;
    for(int k1=0; k1<2; k1++)
    {
            k2 = (k1+1) % 2;            
            i2 = 0;
            
            for (i1=0; i1<nch[k1]; i1++)
            {             
                while( ( (long) (gtime[chnr[k2][i2]]-gtime[chnr[k1][i1]])<taumin) && (i2<nch[k2]-1) ) i2++;
                           
                i2end = i2;
            
                while( ( (long) (gtime[chnr[k2][i2end]]-gtime[chnr[k1][i1]])<taumax) && (i2end<nch[k2]-1) )
                {
                        d = gtime[chnr[k2][i2end]]-gtime[chnr[k1][i1]];
                        i2end++;                        
#ifdef DEBUGGING                        
                        if ((d>=taumax) || (d<taumin))
                        {
                           std::stringstream msg;
                           msg << "d = " << d << "ns" << 
                               "\ntaumin =" << taumin << "ns" <<
                               "\ntaumax =" << taumax << "ns" <<
                               "\ni1 =" << i1 << "\ni2 =" << i2 << 
                               "\ni2end =" << i2end <<           
                               endl;
                           MessageBox(NULL, msg.str().c_str(), "Info", MB_OK);
                           check++;
                           if (check>20) return(1003);
                        }
#endif                        
                        CorrelationHistogram[(unsigned int) ((d-taumin)/taures) ]++;   
                }             
            }
    }

    if (channel!=NULL) delete [] channel;
    if (time!=NULL)    delete [] time;
    if (gtime!=NULL)    delete [] gtime;
    for( k=0; k<2; k++)  
         if (chnr[k]!=NULL) delete [] chnr[k]; 

    double *doublehist;      
    doublehist = new double[2*HistogramLength];   

    for (unsigned long j=0,jj=0;j<HistogramLength;j++,jj+=2) 
    { 
        doublehist[jj]    = CorrelationHistogram[j];
        doublehist[jj+1]  = 0.0;                           // real input data
    }

    // FFT
    four1(doublehist-1,HistogramLength);

    doublehist[0]=0;

    for (unsigned long j=0;j<2*HistogramLength;j+=2) 
        doublehist[j] = sqrt( doublehist[j]*doublehist[j]+doublehist[j+1]*doublehist[j+1] ); 

    for (unsigned long j=0,jj=0;j<HistogramLength/2;j++,jj+=2)  Frequency[j] = doublehist[jj];
    
    if (doublehist!=NULL)   delete [] doublehist;      

    return(0);
		
}	
	
 
