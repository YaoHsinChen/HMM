#ifndef HMM_HEADER_
#define HMM_HEADER_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#ifndef MAX_STATE
#	define MAX_STATE	10
#endif

#ifndef MAX_OBSERV
#	define MAX_OBSERV	26
#endif

#ifndef MAX_SEQ
#	define	MAX_SEQ		200
#endif

#ifndef MAX_LINE
#	define MAX_LINE 	256
#endif

typedef struct{
   char *model_name;
   int state_num;					//number of state
   int observ_num;					//number of observation
   double initial[MAX_STATE];			//initial prob.
   double transition[MAX_STATE][MAX_STATE];	//transition prob.
   double observation[MAX_OBSERV][MAX_STATE];	//observation prob.
} HMM;

static FILE *open_or_die( const char *filename, const char *ht )
{
   FILE *fp = fopen( filename, ht );
   if( fp == NULL ){
      perror( filename);
      exit(1);
   }

   return fp;
}

static void loadHMM( HMM *hmm, const char *filename )
{
   int i, j;
   FILE *fp = open_or_die( filename, "r");

   hmm->model_name = (char *)malloc( sizeof(char) * (strlen( filename)+1));
   strcpy( hmm->model_name, filename );

   char token[MAX_LINE] = "";
   while( fscanf( fp, "%s", token ) > 0 )
   {
      if( token[0] == '\0' || token[0] == '\n' ) continue;

      if( strcmp( token, "initial:" ) == 0 ){
         fscanf(fp, "%d", &hmm->state_num );

         for( i = 0 ; i < hmm->state_num ; i++ )
            fscanf(fp, "%lf", &( hmm->initial[i] ) );
      }
      else if( strcmp( token, "transition:" ) == 0 ){
         fscanf(fp, "%d", &hmm->state_num );

         for( i = 0 ; i < hmm->state_num ; i++ )
            for( j = 0 ; j < hmm->state_num ; j++ )
               fscanf(fp, "%lf", &( hmm->transition[i][j] ));
      }
      else if( strcmp( token, "observation:" ) == 0 ){
         fscanf(fp, "%d", &hmm->observ_num );

         for( i = 0 ; i < hmm->observ_num ; i++ )
            for( j = 0 ; j < hmm->state_num ; j++ )
               fscanf(fp, "%lf", &( hmm->observation[i][j]) );
      }
   }
}

static void dumpHMM( FILE *fp, HMM *hmm )
{
   int i, j;

   //fprintf( fp, "model name: %s\n", hmm->model_name );
   fprintf( fp, "initial: %d\n", hmm->state_num );
   for( i = 0 ; i < hmm->state_num - 1; i++ )
      fprintf( fp, "%.5lf ", hmm->initial[i]);
   fprintf(fp, "%.5lf\n", hmm->initial[ hmm->state_num - 1 ] );

   fprintf( fp, "\ntransition: %d\n", hmm->state_num );
   for( i = 0 ; i < hmm->state_num ; i++ ){
      for( j = 0 ; j < hmm->state_num - 1 ; j++ )
         fprintf( fp, "%.5lf ", hmm->transition[i][j] );
      fprintf(fp,"%.5lf\n", hmm->transition[i][hmm->state_num - 1]);
   }

   fprintf( fp, "\nobservation: %d\n", hmm->observ_num );
   for( i = 0 ; i < hmm->observ_num ; i++ ){
      for( j = 0 ; j < hmm->state_num - 1 ; j++ )
         fprintf( fp, "%.5lf ", hmm->observation[i][j] );
      fprintf(fp,"%.5lf\n", hmm->observation[i][hmm->state_num - 1]);
   }
}

double viterbi(HMM *hmm, char *sequence, int T) {
   double delta[T][hmm->state_num];
   // int track[T][hmm->state_num];

   for(int t = 0; t < T; t++) {
      for(int j = 0; j < hmm->state_num; j++) {
         delta[t][j] = 0;
         // track[t][j] = 0;
      }
   }

   for(int t = 0; t < T; t++) {
      for(int j = 0; j < hmm->state_num; j++) {
         if(t == 0) {
            delta[t][j] = hmm->initial[j] * hmm->observation[sequence[t]-65][j];
            // printf("%l5f\n", delta[t][j]);
         }
         else {
            double p = -1e9;
            for(int i = 0; i < hmm->state_num; ++i) {
               double max = delta[t-1][i] * hmm->transition[i][j];//printf("%l5f\n", max);
               if(max > p) {
                  p = max;
                  // track[t][j] = i;
               }
            }
            delta[t][j] = p * hmm->observation[sequence[t]-65][j];
            // printf("%l5f\n", delta[t][j]);
         }
      }
   }
   double p_star = -1e9;
   for(int j = 0; j < hmm->state_num; j++) {
      if(delta[T-1][j] > p_star) {
         p_star = delta[T-1][j];
      }
   }
   // printf("%l5f\n", p_star);
   return p_star;
}

static void start_testing( const char *listname, HMM *hmm, const int max_num, const char *test_seq, const char *result)
{
   /*load models*/
   FILE *fp = open_or_die( listname, "r" );
   int count = 0;
   double p_star = 0;
   char filename[MAX_LINE] = "";
   while( fscanf(fp, "%s", filename) == 1 ){
      // printf("%s\n",filename);
      loadHMM( &hmm[count], filename );
      // printf("%l5f\n", hmm[3].initial[3]);
      count ++;

      if( count >= max_num ){
         break;
      } 
   }
   fclose(fp);
   /*-------------------------------------*/
   /*load sequence, run viterbi and output*/
   FILE *seq = open_or_die( test_seq, "r");
   FILE *re = open_or_die( result, "w" );

   char seq_token[MAX_LINE] = "";
   int T = 0;
   int modelNum = 0;
   double prob_of_model[5];
   double max_prob = -1e9;

   for (int i = 0; i < 5; ++i)
      prob_of_model[i] = 0;

   while(1) {
      if(fscanf( seq, "%s", seq_token ) > 0 ) {
         T = strlen(seq_token);
         for(int i = 0 ; i < 5; i++) {
            prob_of_model[i] = viterbi(&hmm[i], seq_token, T);
            // printf("%l5f\n", prob_of_model[i]);
            if(max_prob < prob_of_model[i]) {
               max_prob = prob_of_model[i];
               modelNum = i;
               // printf("%d\n", modelNum);
            }
            // printf("%s\n", hmm[modelNum].model_name);
         }
         fprintf(re, "%s", hmm[modelNum].model_name);
         fprintf(re, " %e\n", max_prob);
         modelNum = 0;
         max_prob = -1e9;
      }
      else
         break;
   }
   fclose(seq);
   fclose(re);
}

double **calculateAlpha(HMM *hmm, char *sequence, int T) 
{
   //printf("calculate Alpha...\n");
   double **alpha = (double**)malloc(sizeof(double*)*T);
   for(int i = 0; i < T; i++)
      alpha[i] = (double*)calloc(hmm->state_num, sizeof(double));

   for(int t = 0; t < T; t++) {
      for(int j = 0; j < hmm->state_num; j++) {
         if(t == 0)
            alpha[t][j] = hmm->initial[j] * hmm->observation[sequence[t]-65][j];
         else {
            double p = 0;
            for (int i = 0; i < hmm->state_num; i++)
               p += alpha[t-1][i] * hmm->transition[i][j];
            alpha[t][j] = p * hmm->observation[sequence[t]-65][j];
         }
      }
   }
   // printf("alpha========================\n");
   // for(int i = 0; i < T; i++){
   //    for(int j = 0; j < hmm->state_num - 1; j++)
   //       printf("%.5lf ", alpha[i][j]);
   //    printf("%.5lf\n", alpha[i][hmm->state_num - 1]);
   // }
   return alpha;
}

double **calculateBeta(HMM *hmm, char *sequence, int T) 
{
   //printf("calculate Beta...\n");
   double **beta = (double**)malloc(sizeof(double*)*T);//[hmm->state_num][hmm->T];
    for(int i = 0; i < T; i++)
      beta[i] = (double*)malloc(sizeof(double)*hmm->state_num);

   for(int t = T-1; t >= 0; --t) {
      for(int i = 0; i < hmm->state_num; i++) {
         if(t == T-1)
            beta[t][i] = 1.0;
         else {
            double p = 0;
            for(int j = 0; j < hmm->state_num; j++)
               p += hmm->transition[i][j] * hmm->observation[sequence[t+1]-65][j] * beta[t+1][j];
            beta[t][i] = p;
         }
      }
   }
   // printf("beta========================\n");
   // for(int i = 0; i < T; i++){
   //    for(int j = 0; j < hmm->state_num - 1; j++)
   //       printf("%.5lf ", beta[i][j]);
   //    printf("%.5lf\n", beta[i][hmm->state_num - 1]);
   // }
    return beta;
}

double **calculateGama(double **alpha, double **beta, HMM *hmm, int T)
{
   //printf("calculate Gama...\n");
   double **gama = (double**)malloc(sizeof(double*)*T);
   for(int i = 0; i < T; i++)
      gama[i] = (double*)malloc(sizeof(double)*hmm->state_num);

   for(int t = 0; t < T; t++) {
      for(int i = 0; i < hmm->state_num; i++) {
         double sum = 0.0;
         for(int j = 0; j < hmm->state_num; j++){
            sum = sum + alpha[t][j] * beta[t][j];
         }
         gama[t][i] = (alpha[t][i]*beta[t][i]) / sum;
      }
   }
   // printf("gama========================\n");
   // for(int i = 0; i < T; i++){
   //    for(int j = 0; j < hmm->state_num - 1; j++)
   //       printf("%.5lf ", gama[i][j]);
   //    printf("%.5lf\n", gama[i][hmm->state_num - 1]);
   // }
    return gama;
}

double ***calculateEpsilon(double **alpha, double **beta, HMM *hmm, char *sequence, int T)
{
   //printf("calculate Epsilon...\n");
   double  ***epsilon = (double***)malloc(sizeof(double**)*(T-1));
   for(int i = 0; i < T-1; i++) {
      epsilon[i] = (double**)malloc(sizeof(double*)*hmm->state_num);
      for(int j = 0; j < hmm->state_num; j++) {
         epsilon[i][j] = (double*)malloc(sizeof(double)*hmm->state_num);
      }
   }
   for(int t = 0; t < T-1; ++t) {
      double p = 0;
      for(int i = 0; i < hmm->state_num; i++) {
         for( int j = 0; j < hmm->state_num; j++) {
            p += alpha[t][i] * hmm->transition[i][j] * hmm->observation[sequence[t+1]-65][j] * beta[t+1][j];
         }
      }
      assert(p != 0);
      for(int i = 0; i < hmm->state_num; i++) {
         for(int j = 0; j < hmm->state_num; j++) {
            epsilon[t][i][j] = alpha[t][i] * hmm->transition[i][j] * hmm->observation[sequence[t+1]-65][j] * beta[t+1][j] / p;
         }
      }
   }
   return epsilon;
}

HMM* trainingModel(const char *filename, HMM *hmm_initial) 
{
   FILE *fp = open_or_die (filename, "r");
   char token[MAX_LINE] = "";
   double lines = 0;
   int T;
   double pi_sum[hmm_initial->state_num];
   double gama_sum_A[hmm_initial->state_num];
   double gama_sum_B[hmm_initial->state_num];
   double epsilon_sum[hmm_initial->state_num][hmm_initial->state_num];
   double observation_sum[hmm_initial->observ_num][hmm_initial->state_num];

   for(int i = 0; i < hmm_initial->state_num; i++) {
      pi_sum[i] = 0;
      gama_sum_A[i] = 0;
      gama_sum_B[i] = 0;
   }
   for(int i = 0; i < hmm_initial->state_num; i++) {
      for(int j = 0; j < hmm_initial->state_num; j++) {
         epsilon_sum[i][j] = 0;
      }
   }
   for(int i = 0; i < hmm_initial->observ_num; i++) {
      for(int j = 0; j < hmm_initial->state_num; j++) {
         observation_sum[i][j] = 0;
      }
   }
   while(1) {
      if(fscanf( fp, "%s", token ) > 0 ) {
         double **alpha;
         double **beta;
         double **gama;
         double ***epsilon;
         T = strlen(token);
         alpha = calculateAlpha(hmm_initial, token, T); 
         beta = calculateBeta(hmm_initial, token, T); 
         gama = calculateGama(alpha, beta, hmm_initial, T); 
         epsilon = calculateEpsilon(alpha, beta, hmm_initial, token, T); 

         for(int i  = 0; i < hmm_initial->state_num; i++) {
            pi_sum[i] += gama[0][i];
         }//printf("pi_sum\n");
         for(int t = 0; t < T-1; t++) {
            for(int i = 0; i < hmm_initial->state_num; i++) {
               gama_sum_A[i] += gama[t][i];   
            }
         }//printf("gama_sum_A\n");
         for(int t = 0; t < T; t++) {
            for(int i = 0; i < hmm_initial->state_num; i++) {
               gama_sum_B[i] += gama[t][i];   
            }
         }//printf("gama_sum_B\n");
         for(int t = 0; t < T-1; t++) {
            for(int i = 0; i < hmm_initial->state_num; i++) {
               for(int j = 0; j < hmm_initial->state_num; j++) {
                  epsilon_sum[i][j] += epsilon[t][i][j];
               }
            }
         }//printf("epsilon_sum\n");
         for(int t = 0 ; t < T; t++) {
            int i = token[t]-'A';
            for(int j = 0; j < hmm_initial->state_num; j++) {
               observation_sum[i][j] += gama[t][j];
            }
         }//printf("observation_sum\n");
         lines++;
         // printf("finished\n");
         //printf("=============================\n");
         /* free  all space */
         for(int i = 0; i < T; i++) {
            free(alpha[i]);
            free(beta[i]);
            free(gama[i]);
         }
         free(alpha);
         free(beta);
         free(gama);
         for(int i = 0; i < T-1; i++) {
            for(int j = 0; j < hmm_initial->state_num; j++) {
               free(epsilon[i][j]);
            }
            free(epsilon[i]);
         }
         free(epsilon);
      } // end of if
      else 
         break;
   } // end of while


   // printf("re-estimate all parameters for lamda...\n");
   for(int i = 0; i < hmm_initial->state_num; i++) {
      hmm_initial->initial[i] = pi_sum[i] / lines;
   }//printf("end_of_pi\n");

   for(int i = 0; i < hmm_initial->state_num; i++) {
      for(int j = 0; j < hmm_initial->state_num; j++) {
         hmm_initial->transition[i][j] = epsilon_sum[i][j] / gama_sum_A[i];
      }
   }//printf("end_of_A\n");

   for(int i = 0; i < hmm_initial->observ_num; i++) {
      for(int j = 0; j < hmm_initial->state_num; j++) {
         hmm_initial->observation[i][j] = observation_sum[i][j] / gama_sum_B[j];
      }
   }//printf("end_of_B\n");
   fclose(fp);
   return hmm_initial;
}

double calculate_accuracy(const char *compare, const char *result) {
   double accuracy = 0;
   FILE *com = open_or_die (compare, "r");
   FILE *res = open_or_die (result, "r");
   char com_token[MAX_LINE] = "";
   char res_token[MAX_LINE] = "";
   char prob[MAX_LINE] = "";
   double identical = 0;

   while(1) {
      if((fscanf( com, "%s", com_token ) > 0 ) && fscanf( res, "%s %s", res_token, prob ) > 0 ) {
         if(!strcmp(com_token, res_token)) {
            identical++;
         }
      }
      else
         break;
   }
   accuracy = identical / 2500 * 100;
   return accuracy;
}

#endif
