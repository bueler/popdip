# POPDIP: a POsitive-variables Primal-Dual Interior Point method

This algorithm, implemented in Matlab/Octave, solves the positivity-constrained,
n-dimensional problem

$$\begin{matrix} \min & f(x) \\\\ \text{subject to} & x \ge 0 \end{matrix}$$

This project is just for fun, and not research.  That is, nothing here is an
original contribution beyond well-known algorithms.  Specifically, I regard
this algorithm as a special case of Algorithm 16.1 in section 16.7 of the
textbook Griva, Nash, Sofer, _Linear and Nonlinear Optimization_, 2nd ed., SIAM
Press 2009.

## running examples

Go to the directory `matlab/`, and start Matlab or Octave.  The main algorithm
is in `popdip.m`.

        >> help popdip

The two functioning examples are in `small.m` and `obstacle.m`:

        >> small
        >> obstacle

## documentation

Read `doc/doc.pdf` after doing `make` in `doc/`.

## history

The first version of this project was in my Math 661 Optimization course in
Fall 2018.  See [bueler.github.io/M661F18/](https://bueler.github.io/M661F18/index.html)
and look among the Matlab/Octave codes.
