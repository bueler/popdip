# popdip

POPDIP: a POsitive-variables Primal-Dual Interior Point method

This algorithm, implemented in Matlab/Octave, solves the positivity-constrained,
n-dimensional problem

$$\begin{matrix} \min & f(x) \\\\ \text{subject to} & x \ge 0 \end{matrix}$$

This project is just for fun, and not research.  That is, nothing here is an
original contribution beyond well-known algorithms.  Specifically, I regard
this algorithm as a special case of Algorithm 16.1 in section 16.7 of the
textbook Griva, Nash, Sofer, "Linear and Nonlinear Optimization", 2nd ed., SIAM
Press 2009.

## running the small example

Go to the directory `matlab/`, and start Matlab or Octave.  The two examples
are in `small.m` and `obstacle.m`.  The main algorithm is in `popdip.m`.

To test, run

        >> small
        >> obstacle

## documentation

Read `doc/doc.pdf`.

## history

The first version of this project was in my Math 661 Optimization course in
Fall 2018.  See [bueler.github.io/M661F18/](https://bueler.github.io/M661F18/index.html)]
and look among the Matlab/Octave codes.