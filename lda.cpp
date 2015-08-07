/**
 * File: lda.cpp
 * Date: Dec 2014
 * Author: Qing Wang
 * Description: Data Structure for LDA with position
 */
#pragma once

#include <cstdint>
#include <cstdio>
#include <vector>
#include <cstring>
#include "rng.cpp"
#include "pct.cpp"

using namespace std;

/**
 * Data Structure representing a topic
 */
class Topic {
public:
	/** c_sum: how many words (including duplication) appeared in this topic*/
	uint32_t c_sum;
	/** c[w]: how many times word w appeared in this topic*/
	uint32_t * c;
	/** n[p]: how many times player p appeared in this topic */
	uint32_t * n;
	/** x[p]: sum of x for this topic */
	float* x;
	/** y[p]: sum of y for this topic */
	float* y;
	/** xx[p]: sum of xx for this topic */
	float* xx;
	/** xy[p]: sum of xy for this topic */
	float* xy;
	/** yy[p]: sum of yy for this topic */
	float* yy;

	Topic():c_sum(0),c(NULL),n(NULL),x(NULL),y(NULL),xx(NULL),xy(NULL),yy(NULL) {}

	~Topic() { destroy(); }

	void make(uint32_t n_word, uint32_t n_player)
	{
		destroy();
		c_sum = 0;
		c = new uint32_t[n_word]();
		n = new uint32_t[n_player]();
		x = new float[n_player]();
		y = new float[n_player]();
		xx= new float[n_player]();
		xy= new float[n_player]();
		yy= new float[n_player]();
	}

	void destroy()
	{
		delete[] c;
		delete[] n;
		delete[] x;
		delete[] y;
		delete[] xx;
		delete[] xy;
		delete[] yy;
		memset(this,0,sizeof(*this));
	}
};

/**
 * Structure for Prior parameter
 */
struct Prior {
	/** Prior for document param \theta_d */
	double alpha;
	/** Prior for topic param \phi_k */
	double beta;
	/** NIW prior \mu_0 */
	double mu0_x;
	double mu0_y;
	/** NIW prior \kappa_0 */
	double kappa0;
	/** NIW prior \Lambda_0 */
	double Lambda0_xx;
	double Lambda0_xy;
	double Lambda0_yy;
	/** NIW prior \nu_0 */
	double nu0;
};

/**
 * Structure for LDA with position data
 */
template <typename word_t, typename player_t, typename topic_t>
class LDA
{
public:
	bool use_position;
	/** Size of vocabulary */
	uint32_t n_word;
	/** Size of squads */
	uint32_t n_player;

	Prior prior;

	/** Data Structure for each token in a doc */
	struct Token {
		word_t w;
		player_t p;
		topic_t z;
		float x;
		float y;
	};

	/** Data and assignments. dat[doc][index] = {z,w,p,x,y} */
	vector<vector<Token> > dat;

	/** Topics */
	vector<Topic> topic;

	/** Statistics for Theta, theta[k][doc] = #{z_d,i = k} */
	vector<vector<uint32_t> > theta;

	/** PreComputedTable for speed. */
	PCT pct_log_a, pct_log_b, pct_log_Nb;

	LDA()
	{
		n_word = n_player = 0;
		prior = {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0};
	};

	/**
	 * Read data in row-format.
	 */
	void read_data(const char * filename)
	{
		FILE * f = fopen(filename,"r");
		if(!f)
			printf("Cannot open %s for reading.\n",filename);
		uint32_t d;
		word_t w;
		player_t p;
		float x,y;
		while(true) {
			if(5>fscanf(f,"%u %u %u %f %f\n",&d,&w,&p,&x,&y))
				return;
			n_word = n_word>w?n_word:(w+1);
			n_player = n_player>p?n_player:(p+1);
			while(dat.size()<=d)
				dat.push_back(vector<Token>());
			dat[d].push_back({w,p,0,x,y});
		}
	}

	/**
	 * Read data in lda-c format (without position info)
	 */
	void read_lda_c(const char * filename)
	{
		FILE * f = fopen(filename,"r");
		if(!f)
			printf("Cannot open %s for reading.\n",filename);
		while(true) {
			uint32_t n,w,m;
			if(1>fscanf(f,"%u ",&n))
				return;
			dat.push_back(vector<Token>());
			for(uint32_t i=0;i<n;i++) {
				if(2>fscanf(f,"%u:%u",&w,&m))
					printf("Incomplete data format.\n");
				for(uint32_t j=0;j<m;j++) {
					n_word = n_word>w?n_word:(w+1);
					dat[dat.size()-1].push_back({w,0,0,0,0});
				}
			}
		}
	}

	/**
	 * Set num of topic of LDA
	 * Should be called after reading data.
	 * @param K num of topic
	 */
	void set_K(topic_t K)
	{
		topic.clear();
		topic.resize(K);
		theta.clear();
		for(topic_t k=0;k<K;k++) {
			theta.push_back(vector<uint32_t>());
			theta[k].resize(dat.size());
			topic[k].make(n_word,n_player);
		}
	}

	/**
	 * Config prior and make PCT
	 * Should be called after set_K()
	 */
	void config(const Prior& p, uint32_t buffer_size)
	{
		prior = p;
		pct_log_a.make(buffer_size,log,prior.alpha);
		pct_log_b.make(buffer_size,log,prior.beta);
		pct_log_Nb.make(buffer_size,log,n_word*prior.beta);
	}

	/**
	 * Initialize self
	 * Should be called after config()
	 */
	void init()
	{
		gibbs(true);
	}

	/**
	 * Do gibbs sampling
	 */
	void gibbs(bool first_run=false)
	{
		static vector<uint32_t> x;
		if(x.size()!=dat.size()) {
			x.clear();
			for(uint32_t d=0;d<dat.size();d++)
				x.push_back(d);
		}
		shuffle(&(x[0]),dat.size());
		for(uint32_t d=0;d<dat.size();d++)
			for(uint32_t i=0;i<dat[x[d]].size();i++)
				reassign(x[d],i,first_run);
	}

	/**
	 * Reassign a token. This is called by gibbs().
	 */
	void reassign(uint32_t d, uint32_t i, bool first_run=false)
	{
		static vector<double> prob;
		static vector<double> qprob;
		// Remove statistics
		if(not first_run) {
			topic_t k = dat[d][i].z;
			word_t w = dat[d][i].w;
			player_t p = dat[d][i].p;
			float x = dat[d][i].x;
			float y = dat[d][i].y;

			topic[k].c[w]--;
			topic[k].c_sum--;
			topic[k].n[p]--;
			topic[k].x[p]-=x;
			topic[k].y[p]-=y;
			topic[k].xx[p]-=x*x;
			topic[k].xy[p]-=x*y;
			topic[k].yy[p]-=y*y;
			theta[k][d]--;
		}
		prob.resize(topic.size());
		qprob.resize(topic.size());
		// P(z_{d,i}=k | z_{-(d,i)})
		for(topic_t k=0;k<topic.size();k++)
			prob[k] = pct_log_a(theta[k][d]);
		// P(w_{d,i} | z_{d,i}=k)
		for(topic_t k=0;k<topic.size();k++)
			prob[k] += pct_log_b(topic[k].c[dat[d][i].w])-pct_log_Nb(topic[k].c_sum);
		// P(x_{d,i} | z_{d,i}=k)
		for(topic_t k=0;k<topic.size();k++)
			if(use_position) {
				player_t p=dat[d][i].p;
				uint32_t n = topic[k].n[p];
				double x = topic[k].x[p];
				double y = topic[k].y[p];
				double xx = topic[k].xx[p];
				double xy = topic[k].xy[p];
				double yy = topic[k].yy[p];
				double mu0_x = prior.mu0_x;
				double mu0_y = prior.mu0_y;
				double kappa0 = prior.kappa0;
				double nu0 = prior.nu0;
				double Lambda0_xx = prior.Lambda0_xx;
				double Lambda0_xy = prior.Lambda0_xy;
				double Lambda0_yy = prior.Lambda0_yy;

				double kappa_n = kappa0 + n;
				double mu_x = (kappa0*mu0_x + x)/kappa_n;
				double mu_y = (kappa0*mu0_y + y)/kappa_n;
				double nu_n = nu0 + n;
				double Lambda_xx = Lambda0_xx;
				double Lambda_xy = Lambda0_xy;
				double Lambda_yy = Lambda0_yy;
				if(n>0) {
					Lambda_xx += xx - x*x/n + kappa0*n/kappa_n*(x/n-mu0_x)*(x/n-mu0_x);
					Lambda_xy += xy - x*y/n + kappa0*n/kappa_n*(x/n-mu0_x)*(y/n-mu0_y);
					Lambda_yy += yy - y*y/n + kappa0*n/kappa_n*(y/n-mu0_y)*(y/n-mu0_y);
				}
				double Sigma_xx = Lambda_xx*(kappa_n+1)/kappa_n/(nu_n-2+1);
				double Sigma_xy = Lambda_xy*(kappa_n+1)/kappa_n/(nu_n-2+1);
				double Sigma_yy = Lambda_yy*(kappa_n+1)/kappa_n/(nu_n-2+1);
				double det_Sigma = Sigma_xx*Sigma_yy - Sigma_xy*Sigma_xy;
				double SigmaInv_xx = Sigma_yy/det_Sigma;
				double SigmaInv_xy = -Sigma_xy/det_Sigma;
				double SigmaInv_yy = Sigma_xx/det_Sigma;
				double nu = nu_n - 2 + 1;
				double pp = 0; // prob
				pp += lgamma((nu+2)/2);
				pp -= lgamma(nu/2);
				pp -= log(nu);
				//pp -= log(3.1416);
				pp -= 0.5*log(det_Sigma);
				double dx = (dat[d][i].x-mu_x);
				double dy = (dat[d][i].y-mu_y);
				double s_norm = SigmaInv_xx*dx*dx + SigmaInv_yy*dy*dy + 2*SigmaInv_xy*dx*dy;
				pp -= (nu+2)/2*log(1+s_norm/nu);
				prob[k] += pp;
			}
		// Calculate probability and draw rmultinorm
		prop_exp(&prob[0],topic.size());
		dat[d][i].z = rmultinorm(&prob[0],&qprob[0],topic.size());
		// Update statistics
		if(true) {
			topic_t k = dat[d][i].z;
			word_t w = dat[d][i].w;
			player_t p = dat[d][i].p;
			float x = dat[d][i].x;
			float y = dat[d][i].y;

			topic[k].c[w]++;
			topic[k].c_sum++;
			topic[k].n[p]++;
			topic[k].x[p]+=x;
			topic[k].y[p]+=y;
			topic[k].xx[p]+=x*x;
			topic[k].xy[p]+=x*y;
			topic[k].yy[p]+=y*y;
			theta[k][d]++;
		}
	}

	void output_topics(FILE* fo)
	{
		for(topic_t k=0;k<topic.size();k++) {
			for(word_t w=0;w<n_word-1;w++)
				fprintf(fo,"%u\t",topic[k].c[w]);
			fprintf(fo,"%u\n",topic[k].c[n_word-1]);
		}
	}

	void output_assignments(FILE* fo)
	{
		for(uint32_t d=0;d<dat.size();d++) {
			for(uint32_t i=0;i<dat[d].size()-1;i++)
				fprintf(fo,"%u ",dat[d][i].z);
			fprintf(fo,"%u\n",dat[d][dat[d].size()-1].z);
		}
	}

	void output_theta(FILE* fo)
	{
		for(uint32_t d=0;d<dat.size();d++) {
			for(topic_t k=0;k<topic.size()-1;k++)
				fprintf(fo,"%u\t",theta[k][d]);
			fprintf(fo,"%u\n",theta[topic.size()-1][d]);
		}
	}

	void output_positions(FILE* fo)
	{
		double kappa0 = prior.kappa0;
		double mu0_x = prior.mu0_x;
		double mu0_y = prior.mu0_y;
		for(topic_t k=0;k<topic.size();k++) {
			for(player_t p=0;p<n_player;p++) {
				uint32_t n = topic[k].n[p];
				fprintf(fo,"%le\t%le",
						(kappa0*mu0_x+topic[k].x[p])/(kappa0+n),
						(kappa0*mu0_y+topic[k].y[p])/(kappa0+n));
				if(p!=n_player-1)
					fprintf(fo,"\t");
				else
					fprintf(fo,"\n");
			}
		}
	}

};

