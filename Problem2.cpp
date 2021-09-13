#include <iostream>
#include <armadillo>
#include <fstream>

double u(double x);

int main(){
  double n = 1000.; // number of x values
  arma::vec x = arma::linspace(0, 1, n); //create array with x values
  arma::vec func = arma::vec(n); //define array for function values of u(x)
  for (int i = 0; i < n; i++){ // compute u(x) for each x value
    func(i) = u(x(i));
  }
  // write values of u and x to file
  std::ofstream myfile;
  myfile.open ("x_u_vals.txt");
  scientific(myfile).precision(5);
  for (int i = 0; i < n; i++){ // write u(x) and x to file
    myfile << x(i) << " " << func(i) << "\n";
  }
  myfile.close();
  return 0;
}

// function computing u(x)
double u(double x){
  return 1 - (1-std::exp(-10))*x - std::exp(-10*x);
}
