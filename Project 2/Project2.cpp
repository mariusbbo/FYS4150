#include <iostream>
#include <armadillo>
#include <fstream>
#include <math.h>
#include <assert.h>

#define pi 3.14159265359



// function creating tridiagonal matrix
// input: - size of matrix A, N
//        - value of main diagonal, d
//        - value of sub-diagonal, a
//        - value of super-diagonal, e
// returns tridiagonal matrix A
arma::mat create_tridiagonal(int N, double d, double a, double e){
  arma::mat A = arma::mat(N, N, arma::fill::zeros); // define matrix A filled with zeros
  // insert values into the first and last row of A
  A(0,0) = d;
  A(0,1) = a;
  A(N-1,N-1) = d;
  A(N-1,N-2) = e;
  // insert values into the rest of A
  for (int i = 1; i < N-1; i++){
    A(i, i) = d; // main diagonal
    A(i, i-1) = a; // sub-diagonal
    A(i, i+1) = e; // super-diagonal
  }
  return A;
}

// function computing eigenvalues and eigenvectors analytically
void analytical_eigvals_and_eigvecs(int N, double d, double a){
  arma::vec analytical_eigvals = arma::vec(N, arma::fill::zeros); // array for storing analytical eigenvalues
  arma::mat analytical_eigvecs = arma::mat(N, N, arma::fill::zeros); // array for storing analytical eigenvectors
  // use analytical expressions from project description
  for (int i = 0; i < N; i++){
    analytical_eigvals(i) = d + 2*a*cos((i+1)*pi/(N+1)); // calculate eigenvalues analytically
    for (int j = 0; j < N; j++){
      analytical_eigvecs(j,i) = sin((j+1)*(i+1)*pi/(N+1)); // calculate eigenvectors analytically
    }
  }
  // normalise eigenvectors
  arma::mat eigvecs_a = arma::normalise(analytical_eigvecs);

  // compare with numerical results if N=6 for Problem 3
  if (N == 6){
    // output analytical eigenvalues and eigenvectors
    std::cout << "Analytical eigenvalues:" << "\n";
    for (int i = 0; i < N; i++){
    std::cout << analytical_eigvals(i) << "\n";
    }
    std::cout << "\n";
    std::cout << "Analytical eigenvectors" << "\n";
    std::cout << eigvecs_a << "\n";
  }

  // write analytical eigenvectors corresponding to three lowest eigenvalues to file
  // write only if N=9 (n=10) or N=99 (n=100) for Problem 7
  if (N == 9 || N == 99){
    std::ofstream myfile;
    myfile.open ("analytical_eigenvectors_n10.txt");
    for (int i = 0; i < N; i++){
      myfile << eigvecs_a(i,0) << " " << eigvecs_a(i,1) << " " << eigvecs_a(i,2) << "\n";
    }
    myfile.close();
  }
}

// function for calculating eigenvalues and eigenvectors with
// eig_sym() from armadillo for Problem 3
void compare_eigvals_and_eigvecs(arma::mat A, int N, double a, double d){
  arma::vec eigvals_num; // define array for storing eigenvalues
  arma::mat eigvecs_num; // define array for storing eigenvectors
  arma::eig_sym(eigvals_num, eigvecs_num, A); // find eigenvalues and eigenvectors


  // output numerical eigenvalues and eigenvectors
  std::cout << "Compare eigenvalues and eigenvectors from armadillo with analytical ones." << "\n";
  std::cout << "Some eigenvectors have to be scaled." << "\n\n";
  std::cout << "Armadillo eigenvalues:" << "\n";
  for (int i = 0; i < N; i++){
    std::cout << eigvals_num(i) << "\n";
  }
  std::cout << "\n";
  std::cout << "Armadillo eigenvectors:" << "\n";
  std::cout << eigvecs_num << "\n";

  analytical_eigvals_and_eigvecs(N, d, a);
}

// function finding the maximum off-diagonal element of matrix A
// - update indices k and l of max off-diagonal element
// - returns max off-diagonal value
double max_offdiag_symmetric(const arma::mat& A, int& k, int& l){
  int N = A.n_rows; // get size of the matrix A
  assert(N > 1 && A.is_square()); // check if A is larger than 1x1 and if A is a square matrix

  double maxval = 0; // initialize maximum value of the off-diagonal elements in A
  // loop through only the elements above the main diagonal because of symmetric matrix A
  for (int i = 0; i < N-1; i++){
    for (int j = i+1; j < N; j++){
      // check if the value of the element of A is larger than the present maximum value
      if (abs(A(i,j)) > maxval){
        maxval = abs(A(i,j)); // update maximum value
        k = i; // update row number k for new maximum value
        l = j; // update column number l for new maximum value
      }
    }
  }
  return maxval; // returns maximum value of off-diagonal elements
}

// function for testing the function max_offdiag_symmetric() for Problem 4b
void test_max_off_diag(){
  // declare indices k and l
  int k, l;

  // create matrix to test the function max_offdiag_symmetric
  arma::mat test_matrix = arma::mat(4, 4, arma::fill::eye);
  test_matrix(0,3) = 0.5;
  test_matrix(1,2) = -0.7;
  test_matrix(2,1) = -0.7;
  test_matrix(3,0) = 0.5;

  // output test matrix and maximum off-diagonal value in absolute value
  std::cout << "Test the function finding the maximum off-diagonal element of A." << "\n\n";
  std::cout << "Test matrix:" << "\n\n";
  std::cout << test_matrix << "\n";
  double test_max = max_offdiag_symmetric(test_matrix, k, l);
  std::cout << "Maximum value (in absolute value):" << " " << test_max << "\n";

}

// function performing Jacobis rotation method
// - modifies the input matrices A and R
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l){
  int N = A.n_rows; // size of matrix A

  // declare tan(theta), t, cos(theta), c, and sin(theta), s
  double t;
  double c;
  double s;

  // set values of tan, cos and sin if a_k,l = 0
  if (A(k,l) == 0){
    t = 0;
    c = 1;
    s = 0;
  }
  else{ // compute value, tau, needed to find angle
    double tau = (A(l,l) - A(k,k)) / (2 * A(k,l));

  // find tan(theta), t
  // choose smallest value depending on sign of tau
  if (tau >= 0){
    t = 1/(tau + sqrt(1 + pow(tau, 2)));
  }
  else if (tau < 0){
    t = -1/(-tau + sqrt(1 + pow(tau, 2)));
  }

  c = 1/sqrt(1 + pow(t, 2)); // find cos(theta), c
  s = c*t; // find sin(theta), s
  }
  // set temporary values of elements to not confuse A^(m+1) with A^(m)
  double a_kk, a_ll, a_ik, a_il;

  // update elements of A with both indices k and l
  a_kk = A(k,k)*pow(c,2) - 2*A(k,l)*c*s + A(l,l)*pow(s,2);
  a_ll = A(l,l)*pow(c,2) + 2*A(k,l)*c*s + A(k,k)*pow(s,2);
  A(k,k) = a_kk;
  A(l,l) = a_ll;
  A(k,l) = 0;
  A(l,k) = 0;

  // update elements of A, a_ik, a_il, a_ki, a_li,
  // where i is not equal to k and l
  for (int i = 0; i < N; i++){
    if (i != k && i != l){
      a_ik = A(i,k)*c - A(i,l)*s;
      a_il = A(i,l)*c + A(i,k)*s;
      A(i,k) = a_ik;
      A(k,i) = a_ik;
      A(i,l) = a_il;
      A(l,i) = a_il;
    }
  }
  // set temporary values of elements to not confuse R^(m+1) with R^(m)
  double r_ik, r_il;

  // update elements of R
  for (int i = 0; i < N; i++){
    r_ik = R(i,k)*c - R(i,l)*s;
    r_il = R(i,l)*c + R(i,k)*s;
    R(i,k) = r_ik;
    R(i,l) = r_il;
  }
}

// Jacobi method eigensolver:
// - Runs jacobi_rotate until max off-diagonal element < eps
// - Writes the eigenvalues as entries in the vector "eigenvalues"
// - Writes the eigenvectors as columns in the matrix "eigenvectors"
// - Stops if the number of iterations reaches "maxiter"
// - Writes the number of iterations to the integer "iterations"
// - Sets the bool reference "converged" to true if convergence was reached before hitting maxiter
void jacobi_eigensolver(arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors,
  const int maxiter, int& iterations, bool& converged){
    int N = A.n_rows; // size of matrix A

    // initialize matrix R as the identity matrix
    arma::mat R = arma::mat(N, N, arma::fill::eye);

    // declare indices k and l
    int k, l;

    // run function to find maximum off-diagonal element
    double maxval = max_offdiag_symmetric(A, k, l);

    // runs Jacobis rotation method as long as the maximum off-diagonal
    // element of A is smaller than "epsilon"
    while (maxval > eps){
      // stops the code if the number of iterations become larger
      // than the maximum number of iterations, "maxiter"
      if (iterations > maxiter){
        std::cout << "The number of iterations is larger than maxiter!" << "\n";
        std::exit(0);
      }
      else{
        maxval = max_offdiag_symmetric(A, k, l); // set maximum value of off-diagonal elements of A
        jacobi_rotate(A, R, k, l); // runs Jacobi's rotation method
        iterations += 1; // count number of iterations
      }
    }

    // set converged to true if convergence was reached before the
    // number of iterations reached "maxiter"
    if (iterations < maxiter){
      converged = true;
    }
    else{
      converged = false;
    }

    // define vector and matrix for storing eigenvalues and eigenvectors temporary
    arma::vec eigvals = arma::vec(A.n_rows, arma::fill::zeros);
    arma::mat eigvecs = arma::mat(A.n_rows, A.n_rows, arma::fill::zeros);

    // add eigenvalues from diagonal of A to eigenvalue vector
    for (int i = 0; i < N; i++){
      eigvals(i) = A(i,i);
    }

    arma::uvec eigvals_index = arma::sort_index(eigvals); // get index of eigenvalues in ascending order
    for (int i = 0; i < N; i++){
      eigenvalues(i) = eigvals(eigvals_index(i)); // sort eigenvalues in ascending order
      for (int j = 0; j < N; j++){
        eigvecs(j,i) = R(j,eigvals_index(i)); // sort eigenvectors corresponding to sorted eigenvalues
      }
    }
    // normalise eigenvectors
    eigenvectors = arma::normalise(eigvecs);
    // output eigenvalues and eigenvectors if N=6 for Problem 5
    if (N == 6){
      std::cout << "Compare eigenvalues and eigenvectors computed" << "\n";
      std::cout << "by the rotation method with analytical ones." << "\n";
      std::cout << "Some eigenvectors might have to be scaled." << "\n\n";
      std::cout << "Rotation method eigenvalues:" << "\n";
      std::cout << std::fixed;
      for (int i = 0; i < N; i++){
        std::cout <<  eigenvalues(i) << "\n";
      }
      std::cout << "\n";
      std::cout << "Rotation method eigenvectors:" << "\n";
      std::cout << eigenvectors << "\n";
      double h = 1./N;
      double d = 2/pow(h,2);
      double a = -1/pow(h,2);
      analytical_eigvals_and_eigvecs(N, d, a);
    }
  }

// function finding number of iterations for each N
// - writes N and number of iterations to file
void estimate_iterations(){
  // define lists for N and number of iterations
  arma::vec N_list = arma::linspace(5, 150, 30);
  arma::vec iterations_list = arma::vec(N_list.n_rows);

  // runs Jacobis rotation method for each N in N_list
  for (int i = 0; i < N_list.n_rows; i++){
    int N = N_list(i);
    double h = 1./N; // stepsize
    double d = 2/std::pow(h, 2); // value of the main diagonal
    double a = -1/std::pow(h, 2); // value of the sub-diagonal and super-diagonal
    arma::mat A = create_tridiagonal(N, d, a, a);
    //arma::mat A = arma::mat(N, N, arma::fill::randu); // generate random NxN matrix
    //A = arma::symmatu(A); // make matrix symmetric
    double epsilon = 1e-8; // convergence criterion max off-diagonal value < epsilon
    int maxiter = 100000; // maximum number of iterations
    int iterations = 0; // initialize iteration counter
    bool converged; // bool telling if the algorithm converged
    arma::vec eigenvalues = arma::vec(N, arma::fill::zeros); // vector for storing eigenvalues
    arma::mat eigenvectors = arma::mat(N, N, arma::fill::zeros); // matrix for storing eigenvectors

    // runs Jacobis rotation method
    jacobi_eigensolver(A, epsilon, eigenvalues, eigenvectors,
      maxiter, iterations, converged);

    // append number of iterations for each N to list
    iterations_list(i) = iterations;
  }

  // write N and number of iterations to file
  std::ofstream myfile;
  myfile.open ("Nvals_iterations_random.txt");
  for (int i = 0; i < N_list.n_rows; i++){
    myfile << N_list(i) << " " << iterations_list(i) << "\n";
  }
  myfile.close();
}

// function computing eigenvectors for N=9 (n=10) or N=99 (n=100) and write to file
void eigvecs_diffeq(){
  int N = 9; // size of matrix
  double h = 1./N; // stepsize
  double d = 2/std::pow(h, 2); // value of the main diagonal
  double a = -1/std::pow(h, 2); // value of the sub-diagonal and super-diagonal
  arma::mat A = create_tridiagonal(N, d, a, a); // function to create tridiagonal A
  double epsilon = 1e-8; // convergence criterion max off-diagonal value < epsilon
  int maxiter = 100000; // maximum number of iterations
  int iterations = 0; // initialize iteration counter
  bool converged; // bool telling if the algorithm converged
  arma::vec eigenvalues = arma::vec(N, arma::fill::zeros); // vector for storing eigenvalues
  arma::mat eigenvectors = arma::mat(N, N, arma::fill::zeros); // matrix for storing eigenvectors

  // compute eigenvalues and eigenvectors analytically and with
  // Jacobis rotation method for comparison
  analytical_eigvals_and_eigvecs(N, d, a);
  jacobi_eigensolver(A, epsilon, eigenvalues, eigenvectors, maxiter, iterations, converged);

  // write eigenvectors corresponding to three lowest eigenvalues to file
  std::ofstream myfile;
  myfile.open ("numerical_eigenvectors_n10.txt");
  for (int i = 0; i < N; i++){
    myfile << eigenvectors(i,0) << " " << eigenvectors(i,1) << " " << eigenvectors(i,2) << "\n";
  }
  myfile.close();
}


int main(){
  int N = 6; // size of matrix
  double h = 1./N; // stepsize
  double d = 2/std::pow(h, 2); // value of the main diagonal
  double a = -1/std::pow(h, 2); // value of the sub-diagonal and super-diagonal
  arma::mat A = create_tridiagonal(N, d, a, a); // function to create tridiagonal A
  double epsilon = 1e-8; // convergence criterion max off-diagonal value < epsilon
  int maxiter = 1000; // maximum number of iterations
  int iterations = 0; // initialize iteration counter
  bool converged; // bool telling if the algorithm converged
  arma::vec eigenvalues = arma::vec(N, arma::fill::zeros); // vector for storing eigenvalues
  arma::mat eigenvectors = arma::mat(N, N, arma::fill::zeros); // matrix for storing eigenvectors

  // Problem 3: uncomment below to compare eigvals and eigvecs from armadillo
  //            with analytical eigvals and eigvecs
  //compare_eigvals_and_eigvecs(A, N, a, d);

  // Problem 4b: uncomment below to test the function max_offdiag_symmetric
  //             which finds the maximum off-diagonal element
  //test_max_off_diag();

  // Problem 5b: uncomment below to compare eigvals and eigvecs computed with rotation method
  //             with analytical eigvals and eigvecs
  //jacobi_eigensolver(A, epsilon, eigenvalues, eigenvectors, maxiter, iterations, converged);

  // Problem 6: the call below finds the number of iterations for different N
  //            and writes N and number of iterations to file
  //estimate_iterations();

  // Problem 7: the call below finds eigvals and eigvecs with rotation method and
  //            analytically for n = 10 (can be changed to n = 100) and writes
  //            eigenvectors to file
  eigvecs_diffeq();

  return 0;
}
