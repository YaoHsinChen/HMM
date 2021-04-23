#include "hmm.h"
#include <time.h>

int main(int argc, char *argv[]) {

	char *p;
	long conv = strtol(argv[1], &p, 10);
	int iteration = conv;
	char *model_init = argv[2];
	char *train_seq = argv[3];
	char *model = argv[4]; // output
	long trainingTime = CLOCKS_PER_SEC * 60;
	HMM hmm_initial;

	loadHMM( &hmm_initial, model_init );
	dumpHMM( stderr, &hmm_initial );
	printf("==============================\n");
	printf("System Running...\n");
	printf("training HMM...\n");
	printf("==============================\n");
	int count = 0;
	while((clock()/CLOCKS_PER_SEC < 60) && count < iteration) 
	{
		hmm_initial = *(trainingModel(train_seq, &hmm_initial));
		count++;
	}
	FILE *fp = open_or_die( model, "w");
	dumpHMM( fp, &hmm_initial);
	fclose(fp);

    return 0;
}