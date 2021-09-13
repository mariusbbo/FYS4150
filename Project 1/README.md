# Project 1

Info on how to go through the code for Project 1.

Problem 2:
---------
Problem2.cpp
----
C++ code that computes the exact solution u(x). 
Outputs values of u(x) and x to a file u_x_vals.txt. 

Build command: g++ Problem2.cpp -o Problem2.exe
Run command: ./Problem2.exe

Problem2.py
----
Python script that reads the file x_u_vals.txt and plot u(x) as function of x.
Plot is saved as pdf file.

Run command: python3 Problem2.py

Problem 7:
--------
Problem7.cpp
----
C++ code that uses the general algorithm to solve the matrix equation Av=g for different values of the number of grid points n.
The value of n has to be changed manually. The code outputs v and x to a file x_v_vals_n100.txt. The number after n in the 
filename represents the value of n used in the code. The filename also has to be updated manually.

Build command: g++ Problem7.cpp -o Problem7.exe
Run command: ./Problem7.exe

Problem7.py
----
Python script that reads the file x_v_vals_n100.txt for n from 10 to 10^4 and make a plot that compares v with u(x) for different n. 

Run command: python3 Problem7.py
