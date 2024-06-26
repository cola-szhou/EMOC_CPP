// This file provides some random number utility functions
#pragma once
#ifndef RAND_MAX
#define RAND_MAX 0x7fff
#endif



void advance_random();
void warmup_random(double seed);
void randomize();
double randomperc();
int rnd(int low, int high);
double rndreal(double low, double high);
void random_permutation(int *perm, int size);
double g_rand(double mi, double sigma);
