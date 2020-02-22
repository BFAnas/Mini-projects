In this exercise we are going to use a genetic algorithm for a medical
optimal control problem, i.e. finding the best curative treatment for a tumor if possible, or a
palliative treatment otherwise.
We are going to give a statement of the problem, an overview of the genetic algorithm and it's
adaptation to this particular problem, finally we are going to present our implementation of the
algorithm in C, discuss the results and the improvements performed.

Command for compiling : gcc -o Genetic_main Genetic_main.c RKF78.c ran1.c CrossoverAndMutation.c Fitness.c -lm

