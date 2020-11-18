#include "main.h"
#include "NuclearDensityCalculator.h"
#include "savers/3Dsavers.h"

//#define DEBUG
using namespace std;

int main()
{
  double bound = 7;
  uint n = 3;
  arma::mat m = SolverSchrodinger::solve1D(-bound, bound, n);

#ifdef DEBUG
  m.print("The Matrix");
  auto *checker = new OrthogonalityChecker();
  for (uint i = 0; i < HERM_QUADRA_N_MAX; i++)
    {
      for (uint j = 0; j < HERM_QUADRA_N_MAX; j++)
        {
          cout << "params: " << i << " " << j << "res: " << checker->checkFor(i, j) << '\n';
        }
    }
  delete checker;
#endif
  arma::rowvec z = arma::regspace(-bound, STEP, bound).as_row();
#ifdef DEBUG
    bool test = SolverSchrodinger::test1DSolution(z, m);
    cout << to_string(test);
    if (!test){
      return 1;
    }
#endif
    Saver::saveToCSV(z, m);
  //  cout << "file saved!" << endl;

    NuclearDensityCalculator nuclearDensityCalculator;
  //  nuclearDensityCalculator.printRhoDefs();
    arma::mat res = nuclearDensityCalculator.naive_method(arma::regspace(-10,1, 10), arma::regspace(-20,1, 20));
    arma::cube cube = arma::zeros(41,41,41);
    int r = -20;
    for (int k = 0 ; k < 41 ; k++) {
      for (int t = 0 ; t<40  ; t++) {
        double theta = PI*t/20;
        int x = (r*cos(theta) + 20);
        int y = (r*sin(theta) + 20);
        cube.tube(x,y) = res.row(k);
      }
      r ++;
    }
    //cube.slice(0) = res;
    //cube.slice(1)= res;
   // cube.slice(2) = res;
  //  cube.slice(3)=res;
    std::cout << cubeToDf3(cube);

   // res.print();
    return 0;
}
