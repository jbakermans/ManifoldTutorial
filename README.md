# Manifold Tutorial

<!-- TABLE OF CONTENTS -->
## Table of Contents

* [About](#about)
* [Exercises](#exercises)
	* [Part 1](#part-1)
	* [Part 2](#part-2)
	* [Part 3](#part-3)	
* [Contact](#contact)


## About

This is a Matlab tutorial to accompany [Sara Solla's World Wide Neuro talk](https://www.crowdcast.io/e/sara-sollas-world-wide) from April 29, 2020. It contains Matlab scripts that implement and visualise the more technical sections of the talk in three different parts. The Matlab scripts contain a number of blanks; solving the exercise by filling out the blanks provides hands-on experience with the methods explained in the talk. Probably the best approach is to first watch the part of the talk that comes with an exercise, then solve the corresponding exercise, then check the solution against the provided solutions, and then move on to the next part of the talk.


## Exercises

There are three Matlab exercises, corresponding to three technical parts of the talk: ```part1.m```, ```part2.m```, and ```part3.m```. Each script contains a number of blanks like this: ```'___'``` (triple underscore). Try filling them out and then running the result section by section. Press cmd + Enter to have Matlab run the current section, so it won't produce errors from sections that haven't been solved yet.

### Part 1

Script ```part1.m``` goes with minutes 15 to 25 of the [talk](https://www.crowdcast.io/e/sara-sollas-world-wide). Load synthetic data and find low-dimensional manifold. Project data on neural modes; see how latent trajectory generates data. Solutions can be found in ```part1_solutions.m```.

### Part 2

Script ```part2.m``` goes with minutes 45 to 55 of the [talk](https://www.crowdcast.io/e/sara-sollas-world-wide). Load synthetic data from two different groups of cells, to simulate different recording days, then find how latent dynamics across days are oriented relative to each other. Solutions can be found in ```part2_solutions.m```.

### Part 3

Script ```part3.m``` goes with minutes 55 to 60 of the [talk](https://www.crowdcast.io/e/sara-sollas-world-wide). Use basis vectors for different days in common ambient space to align latent dynamics via Canonical Correlation Analysis, using QR-decomposition. Solutions can be found in ```part3_solutions.m```.

<!-- CONTACT -->
## Contact

[Jacob Bakermans](http://users.ox.ac.uk/~phys1358/) - jacob.bakermans [at] gmail.com

Project Link: [https://github.com/jbakermans/ManifoldTutorial](https://github.com/jbakermans/ManifoldTutorial)
