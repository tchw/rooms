#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <upc_relaxed.h>

#define N 20

shared float dislikes[N * N];
shared float initial_cost[THREADS];
shared float cost[THREADS];
shared int assignment[THREADS * N];

float rand_interval_float(float max) {
	return ((float) rand() / (float) RAND_MAX) * max;
}

unsigned int rand_interval(unsigned int min, unsigned int max) {
	int r;
	const unsigned int range = 1 + max - min, buckets = RAND_MAX / range,
			limit = buckets * range;

	do {
		r = rand();
	} while (r >= limit);

	return min + (r / buckets);
}

void swap(shared int* a, shared int* b) {
	int t = *b;
	*b = *a;
	*a = t;
}

float calculate_cost(int t) {
	float result = 0.0f;
	for (int i = 0; i < N / 2; ++i) {
		result = result
				+ dislikes[N * assignment[N * t + 2 * i]
						+ assignment[N * t + 2 * i + 1]];
	}
	return result;
}

int main() {
	srand(time(0) * (MYTHREAD + 1));

	shared float* dislikes = (float*)upc_all_alloc()

	// generate dislikes
	if (MYTHREAD == 0) {
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) {
				dislikes[N * i + j] = rand_interval_float(10.0f);
				dislikes[N * j + i] = dislikes[N * i + j];
			}
		}

//		for (int i = 0; i < N; ++i) {
//			for (int j = 0; j <= i; ++j) {
//				printf("%.2f  ", dislikes[N * i + j]);
//			}
//			printf("\n");
//		}
	}

	// prevent other threads from proceeding if dislikes are not yet generated
	upc_barrier;

	sleep(MYTHREAD);
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j <= i; ++j) {
			printf("%.2f  ", dislikes[N * i + j]);
		}
		printf("\n");
	}

	// generate initial random solution
	for (int i = 0; i < N; ++i) {
		assignment[N * MYTHREAD + i] = i;
	}
	for (int i = 0; i < N - 1; ++i) {
		int j = rand_interval(i, N - 1);
		swap(&assignment[N * MYTHREAD + i], &assignment[N * MYTHREAD + j]);
	}

//	printf("%%%%%%%%%%%%%%%%%%%%%%%\n");
//	for (int i = 0; i < N / 2; ++i) {
//		printf("[%d], %-3d %-3d\n", MYTHREAD, assignment[MYTHREAD][i][0],
//				assignment[MYTHREAD][i][1]);
//	}
//	printf("%%%%%%%%%%%%%%%%%%%%%%%\n");

// calculate initial cost
	cost[MYTHREAD] = calculate_cost(MYTHREAD);
	initial_cost[MYTHREAD] = cost[MYTHREAD];

	// algorithm
	float T = 1.0f;
	int i = 1000;
	while (i > 0) {
		// generate new solution
		shared int* a = &assignment[N * MYTHREAD + rand_interval(0, 19)];
		shared int* b = &assignment[N * MYTHREAD + rand_interval(0, 19)];
		swap(a, b);

		float cost_ = calculate_cost(MYTHREAD);
		if (cost_ < cost[MYTHREAD]
				|| rand_interval_float(1.0f)
						<= exp((cost[MYTHREAD] - cost_) / T)) {
			cost[MYTHREAD] = cost_;
		} else {
			// restore previous solution
			swap(a, b);
		}

		T *= 0.999f;
		i -= 1;
	}

	// wait for all solutions
	upc_barrier;

	if (MYTHREAD == 0) {
		// find winner
		int winner = 0;
		for (int i = 1; i < THREADS; ++i) {
			if (cost[i] < cost[winner]) {
				winner = i;
			}
		}

		// print gains from all threads
		for (int i = 0; i < THREADS; ++i) {
			printf("[%d] %f -> %f\n", i, initial_cost[i], cost[i]);
		}

		// printf winner assignment
		printf("winner: %d\n", winner);
		for (int i = 0; i < N / 2; ++i) {
			printf("%-3d %-3d\n", assignment[N * winner + 2 * i],
					assignment[N * winner + 2 * i + 1]);
		}
		printf("\n");

	}

	return 0;
}
