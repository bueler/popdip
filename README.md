# POPDIP: a POsitive-variables Primal-Dual Interior Point method

This algorithm, implemented in Matlab/Octave, solves the following positivity-constrained n-dimensional problem with linear equality constraints:

$$\begin{matrix} \min & f(x) \\\\ \text{subject to} & A x = b  \\\\ & x \ge 0 \end{matrix}$$

This project is just for fun, and not research.  That is, nothing here is an original contribution beyond well-known algorithms.  Specifically, POPDIP is a special case of Algorithm 16.1 in section 16.7 of the book Griva, Nash, & Sofer, _Linear and Nonlinear Optimization_, 2nd ed., SIAM Press 2009.

## examples

Go to the directory `matlab/`, and start Matlab or Octave.  The main algorithm is in `popdip.m`.

        >> help popdip

To use `popdip.m` effectively, you will need a driver program.  Three examples are `small.m`, `linear.m`, and `obstacle.m`.  Give them a try:

        >> small
        >> linear
        >> obstacle

## documentation

Please read [the PDF documentation](doc.pdf) and/or check out the [slides](slides.pdf).

These LaTeX documents can be rebuilt, using `pdflatex` and `bibtex`, by doing `make` in the `doc/` directory.

## history

The first version of this project was demonstrated in my Math 661 Optimization course in Fall 2018.  (See [bueler.github.io/M661F18/](https://bueler.github.io/M661F18/index.html) and look among the Matlab/Octave codes.)  The current version 1.2 is an example in the Fall 2024 offering of the same course; see [bueler.github.io/opt/](https://bueler.github.io/opt/).
