In this folder there are three programs to solve linear systems. We suppose to compile and run the program from terminal
- p1general.cpp: it's the general code about Gaussian Elimination. It's developed to work on a loop. We did this to write to a file running times, and errors without running it every time the code. 
		 Input values by terminal are:
		* n: number of initial points;
		* s: step to increase the number of points at the end of the loop by: n=s*n (if the product of n*s is too big it could lead to an overflow)
		* q: numbers of loop's repetitions;
		* bb: value of the diagonal elements of A
		* aa: value of the lower diagonal elements of A
		* cc: value of the upper diagonal elements of A

	Example: do you want to solve Poisson equation with n=1000 ? 
	After compiling the program with armadillo like:
	c++ -o p.x p1general.cpp -larmadillo -llapack -lblas -I/usr/local/include -L/usr/local/lib
	Run the program: 
./p.x 1000 1 1 2 -1 -1 
	If you want to write to a file the analytical/numerical solutions or the elapsed time it is possible to uncomment that specific part in the program.

- p1specialized.cpp: it's the specialized code about Gaussian elimination, i.e. we assume the matrix is filled with aa=cc=-1, bb=2. It works on a loop as p1.general.cpp 
		 Input values by terminal are:
		* n: number of initial points;
		* s: step to increase the number of points in the loop, following the expression: n=s*n;
		* q: number of time you loop's repetitions;

		Example:  do you want to solve Poisson equation with n=1000? 
		After compiling the program with armadillo like:
	      c++ -o p.x p1specialized.cpp -larmadillo -llapack -lblas -I/usr/local/include -L/usr/local/lib
		Run the program: 
./p.x 1000 1 1 2 -1 -1
		Then you should decide if you want to have printed the analytical/numerical solutions/time/relative errors and comment the other things.

-LU.cpp: It gives the solution of a system of linear equation (Av=r) with LU decomposition (only need to declare the matrix A and the vector r). It works on a loop as the other cases.
	Imput values by terminal are:
	* n: number of initial points;
	* s: step to increase the number of points in the loop, following the expression: n=s*n; (if the product of n*s is too big it could lead to overflow)
	* q: number of time you loop's repetitions;
	In this case you can print the numerical solution and the running time.

All the programs are developed to end if s*n isn't an integer (you aren't discretizing the interval (0,1)) OR s<0 OR q<0 OR n<0.
	
