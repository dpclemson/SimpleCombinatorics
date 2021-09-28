#include <iostream>
#include <sys/time.h>

#include "TComplex.h"
#include "TProfile.h"
#include "TFile.h"
#include "TRandomGen.h"
#include "Recursion.C"

using namespace std;

// --- correlation function (generates particles, calls recursion)
void get_corr(int,int,double,int);

// --- recursion function (uses vector of angles to do calculations)
void do_recursion(vector<double>&);

// --- gets system time, executes get_corr inside of over sequences/events
void execute(int,int,int,double,int,unsigned int);
void execute(int,int,int,double);
void execute(int,int,int);

void RanGenX()
{
  int howmany = 10;
  double space = 0.1; // default space between correlated ntuples...
  execute(howmany,700,2,space);
  execute(howmany,700,4,space);
  execute(howmany,700,6,space);
  execute(howmany,700,8,space);
}

void execute(int sequences, int nparticles, int ntuple)
{
  execute(sequences,nparticles,ntuple,0.1);
}

void execute(int sequences, int nparticles, int ntuple, double space)
{
  execute(sequences,nparticles,ntuple,space,-1,0);
}

void execute(int sequences, int nparticles, int ntuple, double space, int sequence_id, unsigned int seed)
{

  int stop = nparticles/ntuple;

  struct timeval Time;

  gettimeofday(&Time,0);
  int begintime = Time.tv_sec;
  //cout<<"begintime is "<<begintime<<endl;

  Init(nparticles); // initialize histograms

  // --- generate the correlations
  for ( int j = 0; j < sequences; ++j )
    {
      if ( j % 10 == 0 ) cout << "Executing sequence j = " << j << endl;
      for ( int i = 1; i < stop; ++i )
	{
	  get_corr(i,ntuple,space,seed);
	}
    }

  // --- Use a time struct to get the current date and time to append to the output file names
  time_t t = time(NULL);
  struct tm tm = *localtime(&t);

  char timestamp[100];
  char daypart[100];
  char timepart[100];

  printf("now: %02d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
  sprintf(timepart,"%02d%02d", tm.tm_hour, tm.tm_min);
  sprintf(daypart,"%02d%02d%02d",tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday);
  sprintf(timestamp,"%s-%s",daypart,timepart);
  cout<<"time stamp is "<<timestamp<<endl;

  char outfilename[100];
  sprintf(outfilename,"OutputFiles/OutFile_%s_k%d.root",timestamp,ntuple);
  if ( sequence_id >= 0 ) sprintf(outfilename,"CondorOutput/OutFile_%04d_k%d.root",sequence_id,ntuple);

  //--- make an output file to write the histograms
  TFile* HistFile = new TFile(outfilename,"recreate");
  HistFile->cd();
  // --- write recursion histo
  for ( int cs = 0; cs < 2; ++cs )
    {
      for(int c = 0; c < maxCorrelator; ++c )
        {
	  hmult_recursion[cs][c]->Write();
        }
    }
  HistFile->Write();
  delete HistFile;

  Delete(); // delete histograms

  gettimeofday(&Time,0);
  int endtime = Time.tv_sec;
  //cout<<"endtime is "<<endtime<<endl;

  int tdiff = endtime-begintime;

  cout<<"End of program."<<endl;
  cout<<"Execution time: "<<tdiff<<" seconds"<<endl;

} //end of void RanGen2()



void get_corr(int nparticles, int ntuple, double space, int seed)
{

  if ( seed < 0 ) seed = 0;
  TRandomMT64 angle(seed);

  vector <double> ang; // inserting pairs into single object

  // int = integer, i can always be at 0,
  // i < 100 means go up to 100 events, ++i checks if it's less
  // than 100 and adds 1 -- if more than 100, than the program stops

  //int stop = nparticles/ntuple;

  for ( int i = 0; i < nparticles; ++i )
    {
      //double means an exact number, phi1 is the name of the double for particle 1.
      double phi1 = angle.Uniform(-3.1415926535,3.1415926535);
      double phi2 = 0; //same as above, but for particle 2
      if ( phi1 > 0 ) phi2 = phi1 - 3.1415926535;
      if ( phi1 < 0 ) phi2 = phi1 + 3.1415926535;

      ang.push_back(phi1); //push_back to stack numbers on return
      ang.push_back(phi2);

      if ( ntuple <= 2 ) continue;

      double phi3 = 0;
      double phi4 = 0;
      if ( phi1 > 0 ) phi3 = phi1 - space;
      if ( phi1 < 0 ) phi3 = phi1 + space;
      if ( phi3 > 0 ) phi4 = phi3 - 3.1415926535;
      if ( phi3 < 0 ) phi4 = phi3 + 3.1415926535;

      ang.push_back(phi3);
      ang.push_back(phi4);

      if ( ntuple <= 4 ) continue;

      double phi5 = 0;
      double phi6 = 0;
      if ( phi3 > 0 ) phi5 = phi3 - space;
      if ( phi3 < 0 ) phi5 = phi3 + space;
      if ( phi5 > 0 ) phi6 = phi5 - 3.1415926535;
      if ( phi5 < 0 ) phi6 = phi5 + 3.1415926535;

      ang.push_back(phi5);
      ang.push_back(phi6);

      if ( ntuple <= 6 ) continue;

      double phi7 = 0;
      double phi8 = 0;
      if ( phi5 > 0 ) phi7 = phi5 - space;
      if ( phi5 < 0 ) phi7 = phi5 + space;
      if ( phi7 > 0 ) phi8 = phi7 - 3.1415926535;
      if ( phi7 < 0 ) phi8 = phi7 + 3.1415926535;

      ang.push_back(phi7);
      ang.push_back(phi8);

    } // end of nparticles for loop

  do_recursion(ang);

  return;

} //end of get_corr



void do_recursion(vector<double>& ang)
{

  for(int h=0;h<maxHarmonic;h++)
    {
      for(int w=0;w<maxPower;w++)
        {
          Qvector[h][w] = TComplex(0.,0.); // initialize Q-vector to zero
        } // end of loop over Qvector
    } // end of loop over maxHarmonic

  int mult = ang.size();

  for ( int i = 0; i < mult; ++i )
    {
      //cout << ang[i] << " " ;
      for(int h=0;h<maxHarmonic;h++)
        {
          double phi = ang[i];
          // do the summation for the Q-vectors
          for(int w=0;w<maxPower;w++)
            {
              Qvector[h][w] += TComplex(cos(h*phi),sin(h*phi));
            } // end of loop over powers
        } // end of loop over harmonics
    } // end of loop over ang vector

  // --- from generic formulas ----------------------------------------------------------------------------
  //  2-p correlations
  int harmonics_Two_Num[2] = {2,-2}; // 2, -2
  int harmonics_Two_Den[2] = {0,0}; // recursion gives the right combinatorics
  TComplex twoRecursion = Recursion(2,harmonics_Two_Num)/Recursion(2,harmonics_Two_Den).Re();
  //double spwTwoRecursion = Recursion(2,harmonics_Two_Den).Re();
  double wTwoRecursion = 1.0;
  hmult_recursion[0][0]->Fill(mult,twoRecursion.Re(),wTwoRecursion);
  hmult_recursion[1][0]->Fill(mult,twoRecursion.Im(),wTwoRecursion);
  //  4-p correlations
  int harmonics_Four_Num[4] = {2,2,-2,-2};
  int harmonics_Four_Den[4] = {0,0,0,0};
  TComplex fourRecursion = Recursion(4,harmonics_Four_Num)/Recursion(4,harmonics_Four_Den).Re();
  //double spwFourRecursion = Recursion(4,harmonics_Four_Den).Re();
  double wFourRecursion = 1.0;
  hmult_recursion[0][2]->Fill(mult,fourRecursion.Re(),wFourRecursion);
  hmult_recursion[1][2]->Fill(mult,fourRecursion.Im(),wFourRecursion);
  //  6-p correlations:
  int harmonics_Six_Num[6] = {2,2,2,-2,-2,-2};
  int harmonics_Six_Den[6] = {0,0,0,0,0,0};
  TComplex sixRecursion = Recursion(6,harmonics_Six_Num)/Recursion(6,harmonics_Six_Den).Re();
  //double spwSixRecursion = Recursion(6,harmonics_Six_Den).Re();
  double wSixRecursion = 1.0;
  hmult_recursion[0][4]->Fill(mult,sixRecursion.Re(),wSixRecursion);
  hmult_recursion[1][4]->Fill(mult,sixRecursion.Im(),wSixRecursion);
  //  8-p correlations
  int harmonics_Eight_Num[8] = {2,2,2,2,-2,-2,-2,-2};
  int harmonics_Eight_Den[8] = {0,0,0,0,0,0,0,0};
  TComplex eightRecursion = Recursion(8,harmonics_Eight_Num)/Recursion(8,harmonics_Eight_Den).Re();
  //double spwEightRecursion = Recursion(8,harmonics_Eight_Den).Re();
  double wEightRecursion = 1.0;
  hmult_recursion[0][6]->Fill(mult,eightRecursion.Re(),wEightRecursion);
  hmult_recursion[1][6]->Fill(mult,eightRecursion.Im(),wEightRecursion);

  // --- print statements for diagnostic purposes
  //cout << twoRecursion.Re() << endl;
  //cout << fourRecursion.Re() << endl;
  //cout << sixRecursion.Re() << endl;
  //cout << eightRecursion.Re() << endl;

} // end do_recursion functions

