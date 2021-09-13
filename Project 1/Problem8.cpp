#include <iostream>
#include <armadillo>
#include <fstream>

double u(double x); // function computing u(x)
double f(double x); // function computing f(x)

int main(){
    double n = 1000.; // number of x values
    arma::vec x = arma::linspace(0, 1, n); //create array with x values
    arma::vec func = arma::vec(n); //define array for function values of u(x)
    for (int i = 0; i < n; i++){ // compute u(x) for each x value
      func(i) = u(x(i));
    }
    // write u(x) and x to file
    std::ofstream myfile;
    myfile.open ("x_u_vals_n1000.txt");
    scientific(myfile).precision(5);
    for (int i = 0; i < n; i++){
      myfile << x(i) << " " << func(i) << "\n";
    }
    myfile.close();

    double h = 1/n; // stepsize
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
    // write v and x to file
    myfile.open ("x_v_vals_n1000.txt");
    scientific(myfile).precision(5);
    for (int i = 0; i < n; i++){
      myfile << x(i) << " " << v(i) << "\n";
    }
    myfile.close();
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
