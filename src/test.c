#include "hmm.h"

int main(int argc, char *argv[]) {

    char *modelList = argv[1];
	char *test_seq = argv[2];
	char *result = argv[3];
	char *compare = "data/test_lbl.txt";
	HMM hmm_initial[5];
	int count[5];
	int max_num = 5;

	start_testing(modelList, &hmm_initial, max_num, test_seq, result); // load output的hmm
	double accuracy = calculate_accuracy(compare, result);
	printf("accracy = %f%%\n", accuracy);
    return 0;
}