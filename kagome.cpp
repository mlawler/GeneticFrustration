#include <iostream>
#include <string>
#include <complex>
#include <random>
#include <unordered_set>
#include <set>
#include <tuple>
#include <boost/range/irange.hpp> // to make python-like irange(0,10) objects
#include <armadillo> // to use matrices and linear algebra!

using namespace std;
using namespace boost;
using namespace arma;

typedef complex<double> dcomplex;
	
// System size:
const unsigned int N1 = 10;
const unsigned int N2 = 10;


//Mat<dcomplex> 
cx_mat* RandomGenerate() {
	cx_mat ones{9,9,fill::ones};
	cx_mat *J = new cx_mat[N1,N2];
    
//	for (auto i : irange(0,N1)) {
//		for (auto j : irange(0,N2)) {
//			J[N1][N2] = cx_mat{9,9,fill::randu} - 0.5*ones;
//		}
//	}
	
	return J;
}

int main() {
	cx_mat* J;
	J = RandomGenerate();
	
	cout << J[1][1] << endl;
	delete [] J;
	return 0;
}