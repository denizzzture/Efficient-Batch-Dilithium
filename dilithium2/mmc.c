#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#define D_MMC 3
#define MAX_ENTRIES 100000

typedef struct {
    int curBest;
    int curPrediction;
    int postfixes[256]; 
} PostfixDictionary;

void initializeDictionary(PostfixDictionary *dict) {
    dict->curBest = 0;
    dict->curPrediction = 0;
    memset(dict->postfixes, 0, sizeof(dict->postfixes));
}

int predict(PostfixDictionary *dict, int *count) {
	assert(count > 0);
    *count = dict->curBest;
    return dict->curPrediction;
}

bool incrementPostfix(PostfixDictionary *dict, int in, bool makeNew) {
    int *curp = &(dict->postfixes[in]);
    int curCount;
    bool newEntry = false;

    if (*curp > 0) {
        curCount = ++(*curp);
    } else if (makeNew) {
        newEntry = true;
        curCount = dict->postfixes[in] = 1;
    } else {
        return false;
    }
    if ((curCount > dict->curBest) || ((curCount == dict->curBest) && (in > dict->curPrediction))) {
        dict->curPrediction = in;
        dict->curBest = curCount;
    }

    return newEntry;
}

double multi_mmc_test(int *data, long len, int alph_size, const int verbose, const char *label) {
    int winner, cur_winner;
    long *entries = malloc(D_MMC * sizeof(long));
    int x[D_MMC] = {0};
	long scoreboard[D_MMC] = {0};
	long i, d, N, C, run_len, max_run_len;

    PostfixDictionary (*M)[MAX_ENTRIES] = malloc(D_MMC * sizeof(PostfixDictionary[MAX_ENTRIES]));

    if (entries == NULL || M == NULL) {
        printf("Memory allocation failed.\n");
        return -1.0;
    }

    if (len < 3) {
        printf("\t*** Warning: not enough samples to run multiMMC test (need more than %d) ***\n", 3);
        free(entries);
        free(M);
        return -1.0;
    }

    N = len - 2;
    winner = 0;
    C = 0;
    run_len = 0;
    max_run_len = 0;

    for (d = 0; d < D_MMC; d++) {
        if (d < N) {
            memcpy(x, data, (d + 1) * sizeof(int));
            initializeDictionary(&M[d][0]);
            M[d][0].curPrediction = data[d+1];
            M[d][0].curBest = 1;
            entries[d] = 1;
        }
    }

    for (i = 2; i < len; i++) {
        bool found_x = false;
        cur_winner = winner;
        memset(x, 0, D_MMC * sizeof(int));

        for (d = 0; (d < D_MMC) && (i - 2 >= d); d++) {
            PostfixDictionary *curp;

            if ((d == 0) || found_x) {
                memcpy(x, data + i - d - 1, (d + 1) * sizeof(int));
                curp = &M[d][0];
                if(curp == NULL) found_x = false;
				else found_x = true;
            }

            printf("Debug: curp->curBest: %d, curp->curPrediction: %d, data[i]: %d, winner: %d, C: %ld\n",
                   curp->curBest, curp->curPrediction, data[i], winner, C);

            if (found_x) {
                int count;
                int prediction = predict(curp, &count);

                if (count > 0 && prediction == data[i]) {
                    if (scoreboard[d] >= scoreboard[winner-1]) winner = d+1;
                    if (d == cur_winner) {
                        C++;
                        if (++run_len > max_run_len) max_run_len = run_len;
                    }
                } else if (d == cur_winner) {
                    run_len = 0;
                }

                if (count == 0 || prediction != data[i]) {
                    if (entries[d] < MAX_ENTRIES) {
                        memcpy(x, data + i - d - 1, (d + 1) * sizeof(int));
                        initializeDictionary(&M[d][entries[d]]);
                        M[d][entries[d]].curPrediction = data[i];
                        M[d][entries[d]].curBest = 1;
                        entries[d]++;
                    }
                } else {
                    incrementPostfix(curp, data[i], entries[d] < MAX_ENTRIES);
                }
            } else if (entries[d] < MAX_ENTRIES) {
                memcpy(x, data + i - d - 1, (d + 1) * sizeof(int));
                initializeDictionary(&M[d][entries[d]]);
                M[d][entries[d]].curPrediction = data[d];
                M[d][entries[d]].curBest = 1;
                entries[d]++;
            }
        }
    }
	free(entries);
    free(M);
    printf("%lf\n", (double) C / (double) N);
    return 0.0;
}




int main() {
    int S[] = {2,1,3,2,1,3,1,3,1};
	multi_mmc_test(S, 9, 3, 1, "Example");
    return 0;
}