#include <iostream>

#include "RanGenX.C"

using namespace std;

int main(int argc, char *argv[])
{

  // --- preamble

  cout << "Processing with arguments ";
  for ( int i = 0; i < argc; ++i )
    {
      cout << argv[i] << " ";
    }
  cout << endl;

  // ---

  int sequence = -1;
  if ( argc > 1 ) sequence = atoi(argv[1]);

  int nevents = 10;
  if ( argc > 2 ) nevents = atoi(argv[2]);

  int get_seed = 0;
  unsigned int seed = 0;
  if ( argc > 3 ) get_seed = atoi(argv[3]);
  if ( get_seed > 0 ) seed = get_seed * 1e6; // worth considering how far apart random seeds need to be

  double space = 0.1;
  if ( argc > 4 ) space = atof(argv[4]);

  cout << "Sequence " << sequence << endl;
  cout << "Number of events " << nevents << endl;
  cout << "Random seed " << seed << endl;
  cout << "Tuplet space " << space << endl;

  // --- nevents, nparticles, correlator, harmonic, sequence, seed
  execute(nevents,250,2,2,space,sequence,seed); // nevents, max 250 tracks, 2-p correlation, 2nd harmonic, sequence number, seed
  execute(nevents,250,3,2,space,sequence,seed); // nevents, max 250 tracks, 3-p correlation, 2nd harmonic, sequence number, seed
  execute(nevents,250,4,2,space,sequence,seed); // nevents, max 250 tracks, 4-p correlation, 2nd harmonic, sequence number, seed
  execute(nevents,250,5,2,space,sequence,seed); // nevents, max 250 tracks, 5-p correlation, 2nd harmonic, sequence number, seed
  execute(nevents,250,6,2,space,sequence,seed); // nevents, max 250 tracks, 6-p correlation, 2nd harmonic, sequence number, seed
  execute(nevents,250,7,2,space,sequence,seed); // nevents, max 250 tracks, 7-p correlation, 2nd harmonic, sequence number, seed
  execute(nevents,250,8,2,space,sequence,seed); // nevents, max 250 tracks, 8-p correlation, 2nd harmonic, sequence number, seed
  execute(nevents,250,9,2,space,sequence,seed); // nevents, max 250 tracks, 9-p correlation, 2nd harmonic, sequence number, seed
  execute(nevents,250,10,2,space,sequence,seed); // nevents, max 250 tracks, 10-p correlation, 2nd harmonic, sequence number, seed
  execute(nevents,250,11,2,space,sequence,seed); // nevents, max 250 tracks, 11-p correlation, 2nd harmonic, sequence number, seed
  execute(nevents,250,12,2,space,sequence,seed); // nevents, max 250 tracks, 12-p correlation, 2nd harmonic, sequence number, seed
  execute(nevents,250,13,2,space,sequence,seed); // nevents, max 250 tracks, 13-p correlation, 2nd harmonic, sequence number, seed
  execute(nevents,250,14,2,space,sequence,seed); // nevents, max 250 tracks, 14-p correlation, 2nd harmonic, sequence number, seed

  execute(nevents,250,2,3,space,sequence,seed); // nevents, max 250 tracks, 2-p correlation, 3rd harmonic, sequence number, seed
  execute(nevents,250,3,3,space,sequence,seed); // nevents, max 250 tracks, 3-p correlation, 3rd harmonic, sequence number, seed
  execute(nevents,250,4,3,space,sequence,seed); // nevents, max 250 tracks, 4-p correlation, 3rd harmonic, sequence number, seed
  execute(nevents,250,5,3,space,sequence,seed); // nevents, max 250 tracks, 5-p correlation, 3rd harmonic, sequence number, seed
  execute(nevents,250,6,3,space,sequence,seed); // nevents, max 250 tracks, 6-p correlation, 3rd harmonic, sequence number, seed
  execute(nevents,250,7,3,space,sequence,seed); // nevents, max 250 tracks, 7-p correlation, 3rd harmonic, sequence number, seed
  execute(nevents,250,8,3,space,sequence,seed); // nevents, max 250 tracks, 8-p correlation, 3rd harmonic, sequence number, seed
  execute(nevents,250,9,3,space,sequence,seed); // nevents, max 250 tracks, 9-p correlation, 3rd harmonic, sequence number, seed
  execute(nevents,250,10,3,space,sequence,seed); // nevents, max 250 tracks, 10-p correlation, 3rd harmonic, sequence number, seed
  execute(nevents,250,11,3,space,sequence,seed); // nevents, max 250 tracks, 11-p correlation, 3rd harmonic, sequence number, seed
  execute(nevents,250,12,3,space,sequence,seed); // nevents, max 250 tracks, 12-p correlation, 3rd harmonic, sequence number, seed
  execute(nevents,250,13,3,space,sequence,seed); // nevents, max 250 tracks, 13-p correlation, 3rd harmonic, sequence number, seed
  execute(nevents,250,14,3,space,sequence,seed); // nevents, max 250 tracks, 14-p correlation, 3rd harmonic, sequence number, seed

  return 0;

}
