# Project2

In this folder there are three programs to find eigenvalues and eigenvectors of a matrix A. We suppose to compile and run the program from terminal
- Jacobi:
main.cpp: it's the general code about finding eigenvalues. It's developed to work on a loop. We did this to write to a file running times, number of iterations and in general to get data for plots without running it every time the code. 
		 Input values by terminal are:
		* n:  number of mesh points er of points of the problem, the matrix for solving the eigenvalue problem will be (n-1)X(n-1).
		* loop: number of loop's repetitions, to perform one loop insert 2. If you insert loop=3, the first run'll be with n=200, the second run with n=200*2=400. So if you want multiple runs, you should put n=5: in this case you can obtain a wide range of data. 
			Indeed the loop looks like 
			...
			for(int ii = 1; ii < loop; ii++){
        		   int t=n*ii;
			... //calculations
			}
		* rho_min;
		* rho_max;
		* interacting: insert 0 if you want non-interacting case, insert 1 if you want the interacting one;
		* eps: tolerance;
		* method: 0 if you want brute-force Jacobi, 1 if you want cyclic jacobi;
		* omega: pulsation for harmonic potentials.

Example: do you want the eigenvalues of matrix A in the interacting case using brute-force Jacobi, rho_min=0, rho_max=25, eps=1e-8, omega=0.25 for just n=200, no more n? To compile: c++ -o p.x main.cpp jacobi.cpp define_matrices.cpp	To run: 
./p.x 200 2 0 25 1 1e-8 0 0.25 
	If you want to write to a file the first eigenvalue on n or first eigenvalue on rho or the elapsed time and interations it is possible to uncomment that specific part in the program.

Note: you should insert an omega even if you are in the non-interacting case, in that situation you can insert 0.

With: 200 2 0 25 0 1e-8 0 0.25
You should get:
Eigenvalues:
2.99511
6.9755
10.9401

With: 200 2 0 25 1 1e-8 0 0.25
You should get:
Eigenvalues
1.2497
2.18873
3.1471

WIth: 200 2 0 25 0 1e-8 1 0.25
You should get:
Eigenvalues
3.289
7.4095
11.4783

With: 200 2 0 25 1 1e-8 1 0.25
You should get:
Eigenvalues
1.25995
2.21053
3.17974

Inside you can find the function jacobi.cpp, which is actually the core of the program, and the function define_matrices.cpp, which is responsible for creating the right matrix A depending on which case are you considering (non-interacting or interacting).
Moreover you can find unit_test.cpp and unit_test_main.cpp. These perform unit tests on jacobi.cpp and define_matrices.cpp. These tests include a comparison between the eigenvalues calculated by the program and the eigenvalues from the theory and verification of the conservation of Frobenius norm in both non-interacting and interacting case with both brute-force Jacobi and cyclic Jacobi. You have just to run these files, 
watch out to not run together main.cpp and unit_test_main.cpp because it would be a conflict between main. To compile: c++ -o p.x unit_test_main.cpp unit_test.cpp jacobi.cpp define_matrices.cpp (only debug mode in Qt, it takes some time :\ )

- eigenpairs_with_lib: 
main.cpp: it's the code, which finds eigenvalues of a matrix A using function tqli in lib.h ("Numerical recipe"). We use this to prove our results. 
		       It's set with n=200, rho_min=0, rho_max=25, omega=0.25 and the non interacting case. To find eigenvalues for the interacting case, please comment loop line 49 and uncomment loop line 53 please. It prints out the run time and the first three eigenvalues.
Non-interacting case n=200:
Eigenvalues:
3.27423
7.36418
11.3968

Non-interacting case n=1000:
Eigenvalues:
3.05655
7.08375
11.1032

Interacting case n=200:
Eigenvalues:
1.25971
2.20967
3.17794

Interacting case n=1000:
1.25343
2.19745
3.16142

To compile: c++ -o p.x main.cpp lib.cpp	
 To run: ./p.x
-omega analysis: 
main.cpp: it's the code, which finds eigenvalues and eigenvectors of a matrix A, specialized to our analysis of the physics of the solutions of Schr\"odinger equation in the interacting case: omega[] and rho_n[] (which represents rho_max) have a set of values (and eps=1e-10) for that reasons.   
	Input values by terminal are:
	* n: number of mesh points of the problem, the matrix for solving the eigenvalue problem will be (n-1)X(n-1).
	If you want you can change the values inside omega, or rho_n, but you need to declare still 4 for omega or 5 for rho_n. If you declare less, you obtain strange values, because the element missed in the array will be filled with a random number.
	However the values already set are suitable for our discussion.
	
The programm creates 5 files. It fills the files following the scheme:  first column -> variable rho // second column -> eigenvalue // third column -> first eigenvector // fourth column -> second eigenvector // and so on.
The first 4 files are about the interacting case with with pulsations declared in vector omega and rho_max declared in the first 4 elements of rho_n. The last refers to the non interacting case with rho_max=rho_n[4].

To compile: c++ -o p.x main.cpp
To run: 
./p.x 200  