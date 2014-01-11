Repo: matmul
---------------
Features
--------------
- *C*
- ``pthread`` *POSIX threading library*
Description
-------------
As part of a class project to implement matrix multiplication using threads,
I wrote this code to practice using the pthreads library and to demonstrate
their use in matrix multiplication.  The code creates two tables A and B of 
size determined by command line arguments and populates them with random
values.  It then distributes the rows of the result table to the number of
threads given as a command line argument, creates that many threads, and calls
on each thread to calculate the values for all rows for which it is
responsible.
