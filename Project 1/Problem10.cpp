#include <iostream>
#include <armadillo>
#include <fstream>

double u(double x); // function computing u(x)
double f(double x); // function computing f(x)

int main(){
    double n = 1000.; // number of x values
    double h = 1/n; // stepsize
    arma::vec x = arma::linspace(0, 1, n); //create array with x values
    arma::vec func = arma::vec(n); //define array for function values of u(x)
    for (int i = 0; i < n; i++){ // compute u(x) for each x value
      func(i) = u(x(i));
    }

    double n_time = 100.; // number of times to run timing test
    double count_general = 0.; // counters for timing tests
    // loop for running timing test
    for (int i = 0; i < n_time; i++){

      // Start measuring time
      auto t1_g = std::chrono::high_resolution_clock::now();

      // general algorithm

      // define arrays
      arma::vec g = arma::vec(n).fill(0); // g
      arma::vec a = arma::vec(n).fill(-1); // subdiagonal
      arma::vec b = arma::vec(n).fill(2); // main diagonal
      arma::vec c = arma::vec(n).fill(-1); // superdiagonal
      arma::vec b_tilde = arma::vec(n).fill(0); // new notation for main diagonal values
      arma::vec g_tilde = arma::vec(n).fill(0); // new notation for g
      arma::vec v = arma::vec(n); // v

      for (int i = 1; i < n-1; i++){ // compute values of g from g_i = h^2*f_i
        g(i) = std::pow(h, 2) * f(x(i));
      }
      // boundary conditions
      v(0) = 0;
      v(n-1) = 0;
      // set known values
      b_tilde(1) = b(1);
      g_tilde(1) = g(1);
      g(1) = g(1) + v(0);
      g(n-2) = g(n-2) + v(n-1);

      // compute new variables b_tilde and g_tilde
      for (int i = 2; i < n-1; i++){
        b_tilde(i) = b(i) - a(i)/b_tilde(i-1)*c(i-1);
        g_tilde(i) = g(i) - a(i)/b_tilde(i-1)*g_tilde(i-1);
      }
      // compute last value of v before boundary value
      v(n-2) = g_tilde(n-2)/b_tilde(n-2);
      // compute v_i
      for (int i = 3; i < n; i++){
        v(n-i) = (g_tilde(n-i) - c(n-i) * v(n-i+1))/b_tilde(n-i);
          }

      // Stop measuring time
      auto t2_g = std::chrono::high_resolution_clock::now();
      // Calculate the elapsed time
      double duration_seconds_general = std::chrono::duration<double>(t2_g - t1_g).count();
      count_general += duration_seconds_general;
    }

      double avg_time_general = count_general/n_time;
      printf ("Average time general algorithm with n=%.0e: %.3e seconds \n", n, avg_time_general);
      double count_special = 0.;
      for (int i = 0; i < n_time; i++){

      // special algorithm

      // Start measuring time
      auto t1_s = std::chrono::high_resolution_clock::now();

      // define arrays
      arma::vec g = arma::vec(n).fill(0); // g
      arma::vec b_tilde = arma::vec(n).fill(0); // new notation for main diagonal values
      arma::vec g_tilde = arma::vec(n).fill(0); // new notation for g
      arma::vec v = arma::vec(n); // v

      // compute values of g from g_i = h^2*f_i
      for (int i = 1; i < n-1; i++){
        g(i) = std::pow(h, 2) * f(x(i));
      }
      // boundary conditions
      v(0) = 0;
      v(n-1) = 0;
      // set known values
      b_tilde(1) = 2;
      g_tilde(1) = g(1);
      g(1) = g(1) + v(0);
      g(n-2) = g(n-2) + v(n-1);

      // compute new variables for main diagonal and g
      for (int i = 2; i < n-1; i++){
        b_tilde(i) = 2 - 1/b_tilde(i-1);
        g_tilde(i) = g(i) + 1/b_tilde(i-1)*g_tilde(i-1);
      }
      // compute last value of v before boundary value
      v(n-2) = g_tilde(n-2)/b_tilde(n-2);
      // compute v
      for (int i = 3; i < n; i++){
        v(n-i) = (g_tilde(n-i) + v(n-i+1))/b_tilde(n-i);
          }
      // Stop measuring time
      auto t2_s = std::chrono::high_resolution_clock::now();
      // Calculate the elapsed time
      double duration_seconds_special = std::chrono::duration<double>(t2_s - t1_s).count();
      count_special += duration_seconds_special;
    }

    double avg_time_special = count_special/n_time;
    printf ("Average time special algorithm with n=%.0e: %.3e seconds \n", n, avg_time_special);

  return 0;
}

// function computing u(x)
double u(double x){
  return 1 - (1-std::exp(-10))*x - std::exp(-10*x);
}

// function computing f(x)
double f(double x){
  return 100*std::exp(-10*x);
}
