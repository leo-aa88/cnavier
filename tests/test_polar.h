// Tests of the polar grid (test_polar.c), run by test_solver.c

#ifndef TEST_POLAR_H_INCLUDED
#define TEST_POLAR_H_INCLUDED

// The suite's check(): passes when value <= limit
typedef void (*check_fn)(const char *name, double value, double limit);

void test_polar(check_fn check);

#endif // TEST_POLAR_H_INCLUDED
