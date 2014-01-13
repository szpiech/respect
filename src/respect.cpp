/*
respect - a program for calculating resampled uni- and bi-variate minor allele
          frequency spectrums

    Copyright (C) 2011  Zachary A Szpiech (szpiechz@umich.edu)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <pthread.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_combination.h>
#include <ctime>
#include "binom.h"
#include "eps_primitives.h"

using namespace std;

/**********EPS creation methods**********/
/*
 * x,y - coordinates in 1/72 of an inch to place the lower left corner of the plot
 * size - length in 1/72 of an inch of the side of the plot 
 * data - 2d data matrix (assumed square)
 * dataSize - dim of matrix
 */
void plotMAFTriangle(ostream &o, int x, int y, int size,
		     double** data, int dataSize, int numBins, long max,
		     bool logScale, string pop1, string pop2);

void plotMAFSpectrum(ostream &o, int x, int y, int size,
		     double* data, int dataSize, double ymin,
		     double ymax, string pop);
void plotAxesTicks(ostream &o, int x, int y, int size, int num);
void plotDiscreteColorScale(ostream &o, int x, int y,int height, int fontSize,
			    double* binColors, string *binLabels, int numBins);
double data2color(long data, long max, int numBins);

/****************************************/

const string ARG_NDR = "--ndr"; //int
const string ARG_NDC = "--ndc"; //int
const string ARG_NLOCI = "-l"; //int
const string ARG_NROW = "-r"; //int
const string ARG_NRESAMP = "-n"; //int
const string ARG_SORT = "-s"; //int
const string ARG_THREADS = "--threads";
const string ARG_DATAFILE = "-f"; //string
const string ARG_OUTPUT = "-o"; //string
const string ARG_PAIRWISE = "--bivariate-all"; //bool
const string ARG_BIPOP1 = "--bivariate-pop1"; //string
const string ARG_BIPOP2 = "--pop2"; //string
const string ARG_HELP = "--help"; //bool
const string ARG_SEED = "--seed"; //int
const string ARG_BINMAX = "--overflow"; //int
const string ARG_NUMBINS = "--num-bins"; //int

const int MISSING = -9;
const int DERIVED = 1;
const int ANCESTRAL = 0;

map<string,bool> argb;
map<string,long int> argi;
map<string,string> args;

typedef struct
{
  short first;
  short second;
} Pair;

typedef struct
{
  int first_index;
  int last_index;
  int numPops;
  int N;
  int id;
  Pair popPair;
} work_order_t;

void parameterDefaults()
{
  argi[ARG_SEED] = time(NULL);
  argb[ARG_PAIRWISE] = false;
  argb[ARG_HELP] = false;
  argi[ARG_NDR] = 0;
  argi[ARG_NDC] = 0;
  argi[ARG_NLOCI] = 0;
  argi[ARG_NROW] = 0;
  argi[ARG_NRESAMP] = 0;
  argi[ARG_SORT] = 0;
  argi[ARG_THREADS] = 1;
  argi[ARG_BINMAX] = -1;
  argi[ARG_NUMBINS] = 6;
  args[ARG_OUTPUT] = "none";
  args[ARG_DATAFILE] = "none";
  args[ARG_BIPOP1] = "none";
  args[ARG_BIPOP2] = "none";
  return;
}

void parse_cmd_line(int argc, char* argv[]);
void printHelp();
bool parameterCheckFail();
void calc_spectrum(void* work_order);
void calc_bispectrum(void* order);
void readData(ifstream &fin, int ndr, int ndc, int nloci, int nrows, int sortBy, 
	      map<string,int> &pop2index, map<int,string> &index2pop);
int drawAllele(double freq, int id);
double** FREQ_DATA;
double** SPECTRUM;
double** BISPECTRUM;

gsl_rng ** r;
const gsl_rng_type * T;



pthread_mutex_t mutex_spectrum = PTHREAD_MUTEX_INITIALIZER;


int main(int argc, char* argv[])
{
  parameterDefaults();
  parse_cmd_line(argc,argv);

  if(argb[ARG_HELP])
    {
      printHelp();
      return 0;
    }

  if(parameterCheckFail())
    {
      return 0;
    }

  bool ALL_BIVARIATE = argb[ARG_PAIRWISE];
  int ndr = argi[ARG_NDR];
  int ndc = argi[ARG_NDC];
  int nloci = argi[ARG_NLOCI];
  int nrows = argi[ARG_NROW];
  int N = argi[ARG_NRESAMP];
  int sortBy = argi[ARG_SORT]-1; //shifts to start at 0
  int num_threads = argi[ARG_THREADS];
  string outputBase = args[ARG_OUTPUT];
  string datafile = args[ARG_DATAFILE];
  string biPop1 = args[ARG_BIPOP1];
  string biPop2 = args[ARG_BIPOP2];
  int numBins = argi[ARG_NUMBINS]-2; //always plot at least 2, transform to be the number of boxes between the extremes
  long max = argi[ARG_BINMAX];
  bool SELF_MAX = (max <= 0);

  ifstream fin;
  fin.open(datafile.c_str());
  if(fin.fail())
    {
      cerr << "Failed to open " << datafile << " for reading.\n";
      return 0;
    }

  map<string,int> pop2index; 
  map<int,string> index2pop;
  
  readData(fin,ndr,ndc,nloci,nrows,sortBy,pop2index,index2pop);

  int numPops = pop2index.size();
   
  work_order_t *order;
  unsigned long int *NUM_PER_THREAD = new unsigned long int[num_threads];
  unsigned long int div = nloci/num_threads;

  gsl_rng_env_setup();
  T = gsl_rng_default;  
  r = new gsl_rng*[num_threads];
  unsigned long int s = argi[ARG_SEED];
  
  for(unsigned long int i = 0; i < num_threads; i++)
    {
      r[i] = gsl_rng_alloc(T);
      gsl_rng_set(r[i],s+i);
    }

  for(int i = 0; i < num_threads; i++)
    {
      NUM_PER_THREAD[i] = 0;
      NUM_PER_THREAD[i] += div;
    }

  for(int i = 0; i < nloci%num_threads;i++)
    {
      NUM_PER_THREAD[i]++;
    }
 
  pthread_t *peer = new pthread_t[num_threads];
  unsigned long int prev_index = 0;
  Pair *popPair;
  short numPairs;
  bool BIVARIATE = false;

  //If we're calculating bivariate MAF spectrums
  //get a list of the pairs of populations
  if(ALL_BIVARIATE)
    {
      numPairs = nCk(numPops,2);
      popPair = new Pair[numPairs];
      
      gsl_combination * c;
      
      c = gsl_combination_calloc (numPops, 2);
      int count = 0;
      do
	{
	  popPair[count].first = gsl_combination_get(c,0);
	  popPair[count].second = gsl_combination_get(c,1);
	  count++;
	}
      while (gsl_combination_next (c) == GSL_SUCCESS);
      gsl_combination_free(c);
      BIVARIATE = true;
    }
  else if(biPop1.compare("none") != 0 && biPop2.compare("none") != 0) //only one pair specified on CMDLN
    {
      bool fail = 0;
      if(pop2index.count(biPop1) == 0)
	{
	  cerr << "ERROR: " << biPop1 << " not in datafile.\n";
	  fail |= 1;
	}
      if(pop2index.count(biPop2) == 0)
	{
	  cerr << "ERROR: " << biPop2 << " not in datafile.\n";
	  fail |= 1;
	}
      if(fail) return 0;
      
      numPairs = 1;
      popPair = new Pair;
      popPair->first = pop2index[biPop1];
      popPair->second = pop2index[biPop2];
      BIVARIATE = true;
    }

  //Do the bivariate claculations
  //This could be more efficient, currently recalculates for each pair
  //Could do all simultaneously, though would take more memory
  if(BIVARIATE) 
    {
      for(int p = 0; p < numPairs; p++) //each pair of pops
	{
	  BISPECTRUM = new double*[N+1];
	  for(int j = 0; j < N+1; j++)
	    {
	      BISPECTRUM[j] = new double[N+1-j];
	      for(int n = 0; n < N+1-j; n++)
		{
		  BISPECTRUM[j][n] = 0;
		}
	    }


	  if(!peer) peer = new pthread_t[num_threads];
	  prev_index = 0;
	  for(int i = 0; i < num_threads; i++)
	    {
	      order = new work_order_t;
	      order->first_index = prev_index;
	      order->last_index = prev_index+NUM_PER_THREAD[i];
	      prev_index += NUM_PER_THREAD[i];
	      order->numPops = numPops;
	      order->id = i;
	      order->N = N;
	      order->popPair = popPair[p];
	      pthread_create(&(peer[i]),
			     NULL,
			     (void *(*)(void*))calc_bispectrum,
			     (void *)order);
	    }
	  
	  for(int i = 0; i < num_threads; i++)
	    {
	      pthread_join(peer[i],NULL);
	    }

	  delete [] peer;
	  peer = NULL;
	  
	  //output files setup
	  string outputFilename = outputBase;
	  string epsFilename;
	  outputFilename += ".";
	  outputFilename += index2pop[popPair[p].first];
	  outputFilename += ".";
	  outputFilename += index2pop[popPair[p].second];
	  outputFilename += ".bi.spec";
	  epsFilename = outputFilename;
	  epsFilename += ".eps";

	  ofstream fout, epsout;

	  long local_max = max;
	  
	  fout.open(outputFilename.c_str());
	  
	  if(fout.fail())
	    {
	      cerr << "ERROR: Could not open " << outputFilename << " for writing.\n";
	      return 0;
	    }
	  
	  for(int n = 0; n < N+1; n++)
	    {
	      fout << index2pop[popPair[p].first] << " ";
	    }
	  fout << endl;

	  for(int n = 0; n < N+1; n++)
	    {
	      fout << double(n)/double(N) << " ";
	    }
	  fout << endl;
	  
	  for(int i = 0; i <  N+1; i++)
	    {
	      fout << index2pop[popPair[p].second] << " " 
		   << double(i)/double(N) << " ";
	      for(int n = 0; n < N+1-i; n++)
		{
		  fout << BISPECTRUM[i][n] << " ";
		  if(SELF_MAX && local_max < BISPECTRUM[i][n]) local_max = BISPECTRUM[i][n];
		}
	      fout << endl;
	    }
      	  
	  
	  fout.close();

	  //A hack to guess the best overflow level
	  //currently takes twice the mean number of counts in all cells
	  //except the one with the most (this will probably be huge)
	  local_max = 2*int((nloci-local_max)/double((N)*(N+1)/2 + N)+0.5);

	  //cout << local_max << "\n";

	  epsout.open(epsFilename.c_str());
	  
	  if(epsout.fail())
	    {
	      cerr << "ERROR: Could not open " << epsFilename << " for writing.\n";
	      return 0;
	    }

	  time_t rawtime;
	  struct tm * timeinfo;
	  string timeStr, creatorStr = argv[0];
	  
	  time ( &rawtime );
	  timeinfo = localtime ( &rawtime );
	  timeStr = asctime(timeinfo);
	  creatorStr += " written by Zachary A Szpiech";

	  setUpEPS(epsout,0,0,612,792,creatorStr,epsFilename,timeStr);
	  plotMAFTriangle(epsout,1.5*72,3.5*72,5*72,BISPECTRUM,N+1,numBins,local_max,0,
			  index2pop[popPair[p].first],index2pop[popPair[p].second]);
	  finalizeEPS(epsout);


	  epsout.close();
	  
	  for(int j = 0; j < N+1; j++)
	    {
	      delete [] BISPECTRUM[j];
	    }
	  delete [] BISPECTRUM;
	}

      delete [] popPair;
    }
  else
    {
      
      SPECTRUM = new double*[numPops];
      for(int i = 0; i < numPops; i++)
	{
	  SPECTRUM[i] = new double[N/2+1];
	  for(int n = 0; n < N/2+1; n++)
	    {
	      SPECTRUM[i][n] = 0;
	    }
	}
      
      prev_index = 0;
      for(int i = 0; i < num_threads; i++)
	{
	  order = new work_order_t;
	  order->first_index = prev_index;
	  order->last_index = prev_index+NUM_PER_THREAD[i];
	  prev_index += NUM_PER_THREAD[i];
	  order->numPops = numPops;
	  order->id = i;
	  order->N = N;
	  pthread_create(&(peer[i]),
			 NULL,
			 (void *(*)(void*))calc_spectrum,
			 (void *)order);
	}
      
      for(int i = 0; i < num_threads; i++)
	{
	  pthread_join(peer[i],NULL);
	}
      
      string outputFilename = outputBase;
      outputFilename += ".uni.spec";
      string epsFilename;

      ofstream fout, epsout;

      fout.open(outputFilename.c_str());

      if(fout.fail())
	{
	  cerr << "ERROR: Could not open " << outputFilename << " for writing.\n";
	  return 0;
	}

      for(int n = 0; n < N/2+1; n++)
	{
	  fout << double(n)/double(N) << " ";
	}
      fout << endl;

      for(int i = 0; i < numPops; i++)
	{
	  fout << index2pop[i] << " ";
	  int local_max = max;
	  for(int n = 0; n < N/2+1; n++)
	    {
	      fout << SPECTRUM[i][n] << " ";
	      if(SELF_MAX && local_max < SPECTRUM[i][n]) local_max = SPECTRUM[i][n];
	    }
	  fout << endl;

	  local_max = 3*int((nloci-local_max)/double(N/2)+0.5);

	  epsFilename = outputBase;
	  epsFilename += ".";
	  epsFilename += index2pop[i];
	  epsFilename += ".uni.spec.eps";
	  
	  epsout.open(epsFilename.c_str());

	  if(epsout.fail())
	    {
	      cerr << "ERROR:  Failed to open " << epsFilename << " for writing.\n";
	      return 0;
	    }
	  
	  time_t rawtime;
	  struct tm * timeinfo;
	  string timeStr, creatorStr = argv[0];
	  
	  time ( &rawtime );
	  timeinfo = localtime ( &rawtime );
	  timeStr = asctime(timeinfo);
	  creatorStr += " written by Zachary A Szpiech";
	  
	  setUpEPS(epsout,0,0,612,792,creatorStr,epsFilename,timeStr);
	  plotMAFSpectrum(epsout,1.5*72,3.5*72,5*72,&(SPECTRUM[i][0]),N/2+1,0,local_max,index2pop[i]);
	  finalizeEPS(epsout);
	  
	  epsout.close();

	}
      fout.close();
      

     
      
      for(int i = 0; i < numPops; i++)
	{
	  delete [] SPECTRUM[i];
	}

      delete [] SPECTRUM;





    }

  delete [] NUM_PER_THREAD;


  return 0;
}

void calc_spectrum(void* order)
{
  work_order_t *p = (work_order_t*)order;
  int first, last, numPops, N, id;
  first = p->first_index;
  last = p->last_index;
  numPops = p->numPops;
  id = p->id;
  N = p->N;

  double** LOCAL_SPECTRUM = new double*[numPops];
  for(int i = 0; i < numPops; i++)
    {
      LOCAL_SPECTRUM[i] = new double[N/2+1];
      for(int n = 0; n < N/2+1; n++)
	{
	  LOCAL_SPECTRUM[i][n] = 0;
	}
    }

  for(int locus = first; locus < last; locus++)
    {
      for(int pop = 0; pop < numPops; pop++)
	{
	  unsigned long int count = 0;
	  for(int n = 0; n < N; n++)
	    {
	      count += drawAllele(FREQ_DATA[pop][locus],id); 
	    }
	  if(count > N/2) count = N - count;
	  LOCAL_SPECTRUM[pop][count]++;
	}
    }
  
  pthread_mutex_lock(&mutex_spectrum);
  for(int i = 0; i < numPops; i++)
    {
      for(int n = 0; n < N/2+1; n++)
	{
	  SPECTRUM[i][n] += LOCAL_SPECTRUM[i][n];
	}
    }
  pthread_mutex_unlock(&mutex_spectrum);
  delete p;
  for(int i = 0; i < numPops; i++)
    {
      delete [] LOCAL_SPECTRUM[i];
    }
  delete [] LOCAL_SPECTRUM;
  return;
}


void calc_bispectrum(void* order)
{
  work_order_t *p = (work_order_t*)order;
  int first, last, numPops, N, id;
  Pair popPair = p->popPair;
  first = p->first_index;
  last = p->last_index;
  numPops = p->numPops;
  id = p->id;
  N = p->N;

  double** LOCAL_BISPECTRUM = new double*[N+1];
  for(int i = 0; i < N+1; i++)
    {
      LOCAL_BISPECTRUM[i] = new double[N+1-i];
      for(int n = 0; n < N+1-i; n++)
	{
	  LOCAL_BISPECTRUM[i][n] = 0;
	}
    }

  //cerr << popPair.first << " " << popPair.second << "\n";

  // cerr << "CHECK 1\n";

  for(int locus = first; locus < last; locus++)
    {
      unsigned long int count1 = 0;
      unsigned long int count2 = 0;
      for(int n = 0; n < N; n++)
	{
	  count1 += drawAllele(FREQ_DATA[popPair.first][locus],id); 
	  count2 += drawAllele(FREQ_DATA[popPair.second][locus],id); 
	}

      //cerr << "1:" << count1 << " " << count2 << endl;

      if(count1+count2 > N) //then 0 is minor allele
	{
	  count1 = N - count1;
	  count2 = N - count2;
	}
      else if(count1+count2 == N)
	{
	  if(gsl_rng_uniform(r[id]) < 0.5)
	    {
	      count1 = N - count1;
	      count2 = N - count2;
	    }
	}
      //cerr << "" << count1 << " " << count2 << endl;
      //cerr << FREQ_DATA[popPair.first][locus] << " " << FREQ_DATA[popPair.second][locus] << endl;
      LOCAL_BISPECTRUM[count1][count2]++;
      
    }
  //cerr << "CHECK 2\n";
  pthread_mutex_lock(&mutex_spectrum);
  for(int i = 0; i < N+1; i++)
    {
      for(int n = 0; n < N+1-i; n++)
	{
	  BISPECTRUM[i][n] += LOCAL_BISPECTRUM[i][n];
	}
    }
  pthread_mutex_unlock(&mutex_spectrum);
  delete p;

  for(int i = 0; i < N+1; i++)
    {
     delete [] LOCAL_BISPECTRUM[i];
    }
  delete [] LOCAL_BISPECTRUM;
  return;
}




int drawAllele(double freq, int id)
{
  double rand = gsl_rng_uniform(r[id]);
  return (rand < freq);
}

void readData(ifstream &fin, int ndr, int ndc, int nloci, int nrows, int sortBy, 
	      map<string,int> &pop2index, map<int,string> &index2pop)
{
  string str;
  string pop;

  for(int i = 0; i < ndr; i++) getline(fin,str);

  streampos begin = fin.tellg();
  int index = 0;
  for(int row = 0; row < nrows; row++)
    {
      for(int i = 0; i <= sortBy; i++)
	{
	  fin >> pop;
	}
      if(pop2index.count(pop) == 0)
	{
	  pop2index[pop] = index;
	  index2pop[index] = pop;
	  index++;
	}
      getline(fin,str);
    }

  fin.seekg(begin);

  int numPops = pop2index.size();

  FREQ_DATA = new double*[numPops];
  double** TOTAL = new double*[numPops];

  for(int i = 0; i < numPops; i++)
    {
      FREQ_DATA[i] = new double[nloci];
      TOTAL[i] = new double[nloci];
      for(int j = 0; j < nloci; j++)
	{
	  FREQ_DATA[i][j] = 0;
	  TOTAL[i][j] = 0;
	}
    }

  int allele;

  for(int row = 0; row < nrows; row++)
    {
      for(int i = 0; i < ndc; i++)
	{
	  fin >> str;
	  if(i == sortBy) pop = str;
	}

      for(int locus = 0; locus < nloci; locus++)
	{
	  fin >> allele;
	  if(allele == DERIVED)
	    {
	      FREQ_DATA[pop2index[pop]][locus]++;
	      TOTAL[pop2index[pop]][locus]++;
	    }
	  else if (allele == ANCESTRAL)
	    {
	      TOTAL[pop2index[pop]][locus]++;
	    }
	}
    }

  for(int i = 0; i < numPops; i++)
    {
      for(int j = 0; j < nloci; j++)
	{
	  FREQ_DATA[i][j] /= TOTAL[i][j];
	}
      
    }
  
  for(int i = 0; i < numPops; i++)
    {
      delete [] TOTAL[i];
    }
  delete [] TOTAL;
  
 return;
}

bool parameterCheckFail()
{
  bool fail = 0;  
  if(argi[ARG_NDR] <= 0)
    {
      cerr << "Non-data rows must be > 0\n";
      fail |= 1;
    }
  if(argi[ARG_NDC] <= 0)
    {
      cerr << "Non-data cols must be > 0\n";
      fail |= 1;
    }
  if(argi[ARG_NLOCI] <= 0)
    {
      cerr << "Number of loci must be > 0\n";
      fail |= 1;
    }
  if(argi[ARG_NROW] <= 0)
    {
      cerr << "Number of rows of data must be > 0\n";
      fail |= 1;
    }
  if(argi[ARG_NRESAMP] <= 0)
    {
      cerr << "Number of alleles to resample must be > 0\n";
      fail |= 1;
    }
  if(argi[ARG_THREADS] < 1)
    {
      cerr << "Number of threads must be >= 1\n";
      fail |= 1;
    }
  if(argi[ARG_SORT] <= 0 || argi[ARG_SORT] > argi[ARG_NDC])
    {
      cerr << "Non-data column index to distinguish populations must be > 0 and <" 
	   << argi[ARG_NDC] << "\n";
      fail |= 1;
    }
  if(args[ARG_DATAFILE].compare("none") == 0)
    {
      cerr << "Datafile must be named something other than \"none\" \n";
      fail |= 1;
    }
  if(argi[ARG_NRESAMP]%2 == 1)
    {
      argi[ARG_NRESAMP]++;
      cerr << "Resample number must be even.  Changed to " 
	   << argi[ARG_NRESAMP] << "\n";
    }
  if(args[ARG_BIPOP1].compare("none") != 0 && args[ARG_BIPOP2].compare("none") == 0)
    {
      cerr << "Must define " << ARG_BIPOP2 << " when " << ARG_BIPOP1 << " is defined.\n";
      fail |= 1;
    }
  if(args[ARG_BIPOP2].compare("none") != 0 && args[ARG_BIPOP1].compare("none") == 0)
    {
      cerr << "Must define " << ARG_BIPOP1 << " when " << ARG_BIPOP2 << " is defined.\n";
      fail |= 1;
    }
  if(argb[ARG_PAIRWISE] && 
     (args[ARG_BIPOP2].compare("none") != 0 && args[ARG_BIPOP1].compare("none") != 0))
    {
      cerr << "WARNING: " << ARG_PAIRWISE << " flag is overridden by setting " 
	   << ARG_BIPOP1 << " and " << ARG_BIPOP2 << endl;
      argb[ARG_PAIRWISE] = false;
    }
  if(argi[ARG_NUMBINS] < 3)
    {
      cerr << "ERROR: fewer than 3 color bins not allowed.\n";
      fail |= 1;
    }
  return fail;
}

void printHelp()
{
  cout << ARG_NDR << " <int> number of non-data rows that precede the data\n"
       << ARG_NDC << " <int> number of non-data columns that precede the data\n"
       << ARG_NLOCI << " <int> number of loci\n"
       << ARG_NROW << " <int> the number of chromosomes sampled\n"
       << "\t(for N individuals, this is 2N for diploid)\n"
       << ARG_SORT << " <int> the non-data column index by which populations are determined\n"
       << ARG_NRESAMP << " <int> number of alleles to resample\n"
       << ARG_THREADS << " <int> number of threads to spawn\n"
       << ARG_SEED << " <int> seed for the RNG, an RNG is allocated for each thread\n"
       << "\tand is seeded with seed + thread id (starting at 0)\n"
       << ARG_DATAFILE << " <string> datafile\n"
       << ARG_BINMAX << " <int> if any bins spill over this number, they are colored in the largest bin color\n"
       << ARG_NUMBINS << " <int> number of bin colors to use in the bivariate plots\n"
       << ARG_OUTPUT << " <string> uses this as the base name for output files\n"
       << ARG_PAIRWISE << " when present the bivariate spectrum is \n"
       << "\tcalculated pairwise between all populations\n"
       << ARG_BIPOP1 << " <string> calculates the bivariate spectrum between the population defined here\n"
       << "\tand the one defined by " << ARG_BIPOP2 << endl
       << ARG_BIPOP2 << " <string> companion command to " << ARG_BIPOP1 << "\n"
       << ARG_HELP << " when present prints this help\n";
}

void parse_cmd_line(int argc, char* argv[])
{
 for(int i = 1; i < argc; i++)
    {
      if(argb.count(argv[i]) > 0)
	{
	  argb[argv[i]] = true;
	}
      else if(argi.count(argv[i]) > 0)
	{
	  argi[argv[i]] = atoi(argv[i+1]);
	  i++;
	}
      else if(args.count(argv[i]) > 0)
	{
	  args[argv[i]] = argv[i+1];
	  i++;
	}
      else
	{
	  cerr << "WARNING: Did not recognize " << argv[i] << " as a valid argument.\n"
	       << "\tYou may encounter unexpected behavior.\n";
	}
    }
 
 if(args[ARG_OUTPUT].compare("none") == 0)
   {
     args[ARG_OUTPUT] = args[ARG_DATAFILE];
   }
 
 return;
}

/**EPS METHODS*/

void plotMAFSpectrum(ostream &o, int x, int y, int size,
		     double* data, int dataSize, double ymin, double ymax, string pop)
{
  int boxSize = size/dataSize;
  size = boxSize*dataSize;
  int fontSize = size/20;

  plotBox(o,x,y,size,0);
  
  double min, max, range;

  /*
  min = data[0];
  max = data[0];

  for(int i = 1; i < dataSize; i++)
    {
      if(min > data[i]) min = data[i];
      if(max < data[i]) max = data[i];
    }
  */
  range = ymax - ymin;

  double height;
  bool mark = false;
  string markText;
  int markTextCount = 0;
  char buffer[50];
  for(int i = 0; i < dataSize; i++)
    {
      height = double(size)*((data[i]-ymin)/ymax);
      if(height > size)
	{
	  height = size;
	  mark = true;
	}
      plotRect(o,x+i*boxSize,y,boxSize-1,height,0.5);
      if(mark)
	{
	  plotTextCenter(o,x+(i+0.5)*boxSize,y+height,fontSize,0,"*");
	  markText = "* ";
	  //sprintf(buffer,"%.2f",double(i)/double(dataSize));
	  //markText += buffer;
	  //markText += ",";
	  sprintf(buffer,"%ld",long(data[i]));
	  markText += buffer;
	  //markText += ")";
	  plotText(o,x,y-(4+markTextCount)*fontSize,fontSize,0,markText);
	  mark = false;
	  markTextCount++;
	}
    }

  
  string xaxis = "MAF (";
  xaxis += pop;
  xaxis += ")";
  string yaxis = "Count";
  string ytick1, ytick2, ytick3;

  sprintf(buffer,"%d",int(ymin));
  ytick1 = buffer;
  sprintf(buffer,"%d",int((range/2.0)+ymin));
  ytick2 = buffer;
  sprintf(buffer,"%d",int(ymax));
  ytick3 = buffer;

  //x-axis
  plotTextCenter(o,x+size/2,y-3*fontSize,fontSize,0,xaxis);
  plotTextCenter(o,x+boxSize/2,y-fontSize,fontSize,0,"0");
  plotTextCenter(o,x+((dataSize-1)/2)*boxSize+boxSize/2,y-fontSize,fontSize,0,"0.25");
  plotTextCenter(o,x+(dataSize-1)*boxSize+boxSize/2,y-fontSize,fontSize,0,"0.5");
  //y-axis 
  plotTextCenter(o,x-3*fontSize,y+(size/2),fontSize,90,yaxis);  
  plotTextRight(o,x-fontSize/2,y-fontSize/4,fontSize,0,ytick1);
  plotTextRight(o,x-fontSize/2,y+(size/2)-fontSize/4,fontSize,0,ytick2);
  plotTextRight(o,x-fontSize/2,y+(size)-fontSize/4,fontSize,0,ytick3);


  //plotAxesTicks(o,x,y,size,4);
}


void plotAxesTicks(ostream &o, int x, int y, int size, int num)
{
  o << "newpath\n"
    << "\t0 setgray\n";

  for(int i = 0; i < num+1; i++)
    {
      o << "\t" << x+(double(size)/double(num))*i << " " << y << " moveto\n"
	<< "\t" << 0 << " " << size/80 << " rmoveto\n"
	<< "\t" << 0 << " " << -size/40 << " rlineto\n"
	<< "\t1 setlinewidth\n"
	<< "stroke\n";
    }
  
  o << "\t" << x << " " << y << " moveto\n";

  for(int i = 0; i < num+1; i++)
    {
      o << "\t" << x << " " << y+(double(size)/double(num))*i << " moveto\n"
	<< "\t" << size/80 << " " << 0 << " rmoveto\n"
	<< "\t" << -size/40 << " " << 0 << " rlineto\n"
	<< "\t1 setlinewidth\n"
	<< "stroke\n";
    }
    
}


double data2color(long data, long max, int numBins)
{
  double color;
  if(data == 0) return data;
  
  for(int i = 1; i <= numBins; i++)
    {
      if(data <= i*max/numBins) return i/double(numBins+1);
    }
  return 1;
}

/*
 * x,y - coordinates in 1/72 of an inch to place the lower left corner of the plot
 * size - length in 1/72 of an inch of the side of the plot 
 * data - pointer to 2d data matrix (assumed square)
 * dataSize - dim of matrix
 */
void plotMAFTriangle(ostream &o, int x, int y, int size,
		     double** data, int dataSize, int numBins, long max, 
		     bool logScale, string pop1, string pop2)
{
  int boxSize = size/dataSize;
  size = boxSize*dataSize; 

  //Bounding triangle
  plotLLTriangle(o,x,y,size,1.0);
 
  //plot color boxes for each data point
  for(int i = 0; i < dataSize; i++)
    {
      for(int j = 0; j < dataSize-i; j++)
	{
	  if(i == dataSize-1-j)
	    {
	      plotLLTriangle(o,i*boxSize+x,j*boxSize+y,boxSize,data2color(data[i][j],max,numBins));
	    }
	  else plotBox(o,i*boxSize+x,j*boxSize+y,boxSize,data2color(data[i][j],max,numBins));
	}
    }
  
  int fontSize = size/20;

  string xaxis = "MAF (";
  xaxis += pop1;
  xaxis += ")";
  string yaxis = "MAF (";
  yaxis += pop2;
  yaxis += ")";

  //x-axis
  plotTextCenter(o,x+(size/2),y-3*fontSize,fontSize,0,xaxis);
  plotTextCenter(o,x+boxSize/2,y-fontSize,fontSize,0,"0");
  plotTextCenter(o,x+((dataSize-1)/2)*boxSize+boxSize/2,y-fontSize,fontSize,0,"0.5");
  plotTextCenter(o,x+(dataSize-1)*boxSize+boxSize/2,y-fontSize,fontSize,0,"1.0");
  //y-axis 
  plotTextCenter(o,x-3*fontSize,y+(size/2),fontSize,90,yaxis);  
  plotTextRight(o,x-fontSize/2,y+fontSize/4,fontSize,0,"0");
  plotTextRight(o,x-fontSize/2,y+((dataSize-1)/2)*boxSize+fontSize/4,fontSize,0,"0.5");
  plotTextRight(o,x-fontSize/2,y+(dataSize-1)*boxSize+fontSize/4,fontSize,0,"1.0");

  //colorscale
  // int numBins = 5;
  double *binColors = new double[numBins+2];
  string *binLabels = new string[numBins+2];
  char buffer[50];

  binColors[0] = 0;
  binLabels[0] = "0";
  binLabels[1] = "1-";

  for(int i = 1; i <= numBins; i++)
    {      
      sprintf(buffer,"%ld",i*max/numBins);
      binLabels[i]+=buffer;
      sprintf(buffer,"%ld",i*max/numBins+1);
      binLabels[i+1]+=buffer;
      binLabels[i+1]+="-";
      binColors[i] = i/double(numBins+1);
    }
  sprintf(buffer,"%ld",max);
  binLabels[numBins+1] = ">";
  binLabels[numBins+1] += buffer;
  binColors[numBins+1] = 1;
  

  /*  for(int i = 0; i < max; i++)
    {
      binColors[i] = double(i)/double(max);
      sprintf(buffer,"%.1f",binColors[i]);
      binLabels[i] = buffer;
    }
  */
  plotDiscreteColorScale(o,x+1.1*size,y,size,fontSize,binColors,binLabels,numBins+2);
  
  delete [] binLabels;
  delete [] binColors;
  return;
}

void plotDiscreteColorScale(ostream &o, int x, int y,int height, int fontSize,
			    double* binColors, string *binLabels, int numBins)
{
  int width = double(height)*0.05;
  int boxHeight = height/numBins;
  height = boxHeight*numBins; 
 
  plotRect(o,x,y,width,height,1);

  for(int i = 0; i < numBins; i++)
    {
      plotRect(o,x,y+boxHeight*i,width,boxHeight,binColors[i]);
      plotText(o,x+width+fontSize/2,y+boxHeight*(i+0.5)-fontSize/4,fontSize,0,binLabels[i]);
    }

  return;
}
