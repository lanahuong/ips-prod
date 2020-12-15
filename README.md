# Introduction to scientific programming project 2

## Subject

The goal was to compute and visualise the local density of a nuclear system given the density matrix.
This uses some code ips-dev to compute de basis of functions. Visualisation is done using Povray and
gnuplot.

Unit tests and a good management of the git repository were also required.

## Dependencies

* gnuplot
* armadillo
* povray

## Usage

To build the whole project, go to the root of the project and run :

```
make
```

To run the unit tests you can run `bin/tests --gtest_color=yes`.

To build only the tests use :

```
make tests
```

To generate the documentation, run from the root :

```
make doc
```

The documentation is then available in HTML format at `doc/html/index.html`.

To get the solutions, run :

```
bin/nuclearDensity
```

To generate the visual of the density as a function of _z_ and _r_ you can run :

```
gnuplot scripts/density-r-z.p
```

The resulting image is `tmp/density-r-z.png`.

To genarate the visual of the density in 3 dimensions you can run from the root of the project :

```
povray visu.pov
```

The resulting image is `visu.png`

To clear the project build, run :

```
make clean
```