#include "main.h"
#include "NuclearDensityCalculator.h"
#include "savers/3Dsavers.h"

using namespace std;

int main() {
    NuclearDensityCalculator nuclearDensityCalculator;
    int r_bound = 10;
    int r_step = 1;
    int r_len = 2 * r_bound / r_step + 1;
    int z_bound = 20;
    int z_step = 2;
    int z_len = 2 * z_bound / z_step + 1;

    /* Chrono start*/
    auto start = std::chrono::high_resolution_clock::now();
    arma::mat res = nuclearDensityCalculator.naive_method(arma::regspace(-r_bound, r_step, r_bound), arma::regspace(-z_bound, z_step, z_bound));
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cerr << "time " << duration.count() / 1000000.0 << "s" << std::endl;
    /* Chrono end */


    arma::cube cube = arma::zeros(r_len, r_len, z_len);
    int r = -r_bound;
    for (int k = 0; k < r_len; k++) {
        for (int t = 0; t < 40; t++) {
            double theta = PI * t / 20;
            int x = (r * cos(theta) + r_bound);
            int y = (r * sin(theta) + r_bound);
            cube.tube(x, y) = res.row(k);
        }
        r++;
    }

    std::cout << cubeToDf3(cube);
    return 0;
}
