/**
 * File: main.cpp
 * Date: Dec 2014
 * Author: Qing Wang
 * Description: Main entry of LDA program
 */
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include "lda.cpp"

int main(int argc, char* argv[])
{
	if(argc<20) {
		printf(" Usage: %s\n",argv[0]);
		printf("    input_file_name\n");
		printf("    out_topics_file_name\n");
		printf("    out_position_file_name\n");
		printf("    out_assign_file_name\n");
		printf("    out_theta_file_name\n");
		printf("    use_pos\n");
		printf("    n_word\n");
		printf("    K\n");
		printf("    alpha\n");
		printf("    beta\n");
		printf("    mu0_x\n");
		printf("    mu0_y\n");
		printf("    kappa0\n");
		printf("    Lambda0_xx\n");
		printf("    Lambda0_xy\n");
		printf("    Lambda0_yy\n");
		printf("    nu0\n");
		printf("    n_iter\n");
		printf("    seed\n");
		printf("  (See main.cpp for details)\n");
		return 0;
	}
	// argv[1] is input, argv[2] is out_topics,
	// argv[3] is out_pos, argv[4] is out_assign,
	// argv[5] is out_theta,
	uint32_t _ = 6; 
	uint32_t use_pos =  strtol(argv[_++],NULL,10);
	uint32_t n_word =   strtol(argv[_++],NULL,10);
	uint32_t K =        strtol(argv[_++],NULL,10);
	double alpha =      strtod(argv[_++],NULL);
	double beta =       strtod(argv[_++],NULL);
	double mu0_x =      strtod(argv[_++],NULL);
	double mu0_y =      strtod(argv[_++],NULL);
	double kappa0 =     strtod(argv[_++],NULL);
	double Lambda0_xx = strtod(argv[_++],NULL);
	double Lambda0_xy = strtod(argv[_++],NULL);
	double Lambda0_yy = strtod(argv[_++],NULL);
	double nu0 =        strtod(argv[_++],NULL);
	uint32_t n_iter =   strtol(argv[_++],NULL,10);
	uint64_t seed =     strtol(argv[_++],NULL,10);

	lcg64(seed);
	LDA<uint32_t,uint32_t,uint8_t> lda;
	// Order matters: 1. Read data. 2. set K. 3. config prior
	//lda.read_lda_c("../hdp/ap/ap.dat");
	lda.read_data(argv[1]);
	lda.n_word = n_word;
	printf("Read file done: %u doc, %u word\n",lda.dat.size(),lda.n_word);
	lda.set_K(K);
	printf("Set K = %u\n",lda.topic.size());
	// alpha, beta, mu0_x, mu0_y, kappa0, Lambda0_xx, Lambda0_xy, Lambda0_yy, nu0
	lda.config({alpha,beta,mu0_x,mu0_y,kappa0,Lambda0_xx,Lambda0_xy,Lambda0_yy,nu0},65536*8);
	lda.use_position = use_pos>0?true:false;

	//lda.print();
	lda.init();
	for(uint32_t i=0;i<n_iter;i++) {
		lda.gibbs();
		if((i+1)%10==0) {
			printf("iter: %u\n",i+1);
			fflush(stdout);
		}
	}
	FILE* f_topic = fopen(argv[2],"w");
	lda.output_topics(f_topic);
	fclose(f_topic);
	FILE* f_pos = fopen(argv[3],"w");
	lda.output_positions(f_pos);
	fclose(f_pos);
	FILE* f_assign = fopen(argv[4],"w");
	lda.output_assignments(f_assign);
	fclose(f_assign);
	FILE* f_theta = fopen(argv[5],"w");
	lda.output_theta(f_theta);
	fclose(f_theta);
	return 0;
}

