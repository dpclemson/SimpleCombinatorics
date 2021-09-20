const int maxCorrelator = 12;
const int maxHarmonic = 10;
const int maxPower = 9;
TComplex Qvector[maxHarmonic][maxPower];

TProfile* hmult_recursion[2][maxCorrelator];

TComplex Q(int, int);
TComplex Recursion(int, int*);
TComplex Recursion(int, int*, int, int);

void Init();

TComplex Recursion(int n, int* harmonic)
{
  return Recursion(n,harmonic,1,0); // 1 and 0 are defaults from above
}

TComplex Recursion(int n, int* harmonic, int mult, int skip)
{
 // Calculate multi-particle correlators by using recursion (an improved faster version) originally developed by
 // Kristjan Gulbrandsen (gulbrand@nbi.dk).

  int nm1 = n-1;
  TComplex c(Q(harmonic[nm1], mult));
  if (nm1 == 0) return c;
  c *= Recursion(nm1, harmonic);
  if (nm1 == skip) return c;

  int multp1 = mult+1;
  int nm2 = n-2;
  int counter1 = 0;
  int hhold = harmonic[counter1];
  harmonic[counter1] = harmonic[nm2];
  harmonic[nm2] = hhold + harmonic[nm1];
  TComplex c2(Recursion(nm1, harmonic, multp1, nm2));
  int counter2 = n-3;
  while (counter2 >= skip) {
    harmonic[nm2] = harmonic[counter1];
    harmonic[counter1] = hhold;
    ++counter1;
    hhold = harmonic[counter1];
    harmonic[counter1] = harmonic[nm2];
    harmonic[nm2] = hhold + harmonic[nm1];
    c2 += Recursion(nm1, harmonic, multp1, counter2);
    --counter2;
  }
  harmonic[nm2] = harmonic[counter1];
  harmonic[counter1] = hhold;

  if (mult == 1) return c-c2;
  return c-double(mult)*c2;

}

TComplex Q(int n, int p)
{
  // Using the fact that Q{-n,p} = Q{n,p}^*.
  if(n>=0){return Qvector[n][p];}
  return TComplex::Conjugate(Qvector[-n][p]);
} // TComplex Q(int n, int p)
// -------------------------------------------------------------------------------


void Init()
{
  for ( int cs = 0; cs < 2; ++cs )
    {
      for(int c = 0; c < maxCorrelator; ++c )
        {
          hmult_recursion[cs][c] = new TProfile(Form("hmult_recursion_%d_%d",cs,c),"",700,-0.5,699.5,-1.1,1.1);
        }
    }
}

void Delete()
{
  for ( int cs = 0; cs < 2; ++cs )
    {
      for(int c = 0; c < maxCorrelator; ++c )
        {
          delete hmult_recursion[cs][c];
        }
    }
}

