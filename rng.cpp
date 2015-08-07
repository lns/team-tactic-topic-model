/**
 * File: rng.cpp
 * Date: Oct 2013, Dec 2014
 * Author: Qing Wang
 * Description: LCG64 Random Number Generator and some facilities
 */
#pragma once

#include <cstdint>
#define A_Default (18145460002477866997ull)
#define A_LEcuyer1 (2862933555777941757ull)
#define A_LEcuyer2 (3202034522624059733ull)
#define A_LEcuyer3 (3935559000370003845ull)


static uint64_t __lcg64_r = 0;

inline uint64_t lcg64(void)
{
	uint64_t * const r = &__lcg64_r;
	const uint64_t a = A_Default;
	const uint64_t b = 1;
	(*r) = ((*r)*a + b);
	return (*r);
}

inline void lcg64(uint64_t seed)
{
	__lcg64_r = seed;
}

inline double drand(void)
{
	uint64_t M = ~0ull;
	return((double)lcg64()/M);
}

template <typename T>
void shuffle(T* x, uint32_t len)
{
	for(uint32_t i=0;i<len-1;i++)
	{
		uint32_t j = i+lcg64()%(len-1-i);
		// swap x[i] & x[j]
		T t = x[i];
		x[i] = x[j];
		x[j] = t;
	}
}

/**
 * Proportional (for numerical stability) Exponential function
 */
#include <cmath>
template <typename T>
inline void prop_exp(T* x, uint32_t len)
{
	T mean = 0;
	for(uint32_t i=0;i<len;i++)
		if(std::isfinite(x[i]))
			mean += x[i];
	mean /= len;
	for(uint32_t i=0;i<len;i++)
		x[i] = exp(x[i]-mean);
}

/**
 * Input:
 * 		Probability (without normalization) in p
 * 		length of p and cum: len
 * Output:
 * 		return in [0,len-1]
 * 		Cumulative Prob in (without normalization) in cum
 */
template <typename T>
inline uint32_t rmultinorm(T* p, T* cum, uint32_t len, bool cal_cum=true)
// Sequential search is faster for small size prob array
#define SEQ_SEARCH
{
	if(cal_cum)
	{
		cum[0] = p[0];
		for(uint32_t i=1;i<len;i++)
			cum[i] = cum[i-1] + p[i];
	}
	double r = drand()*cum[len-1];
//	for(uint32_t i=0;i<len;i++)
//		printf("p[%u]: %le  q[%u]: %le\n",i,p[i],i,cum[i]);
#ifndef SEQ_SEARCH
	uint32_t med = len/2;
	uint32_t low=0;
	uint32_t high=len-1;
	do {
		if(cum[med]>r)
			high = med;
		else
			low = med;
		med = (high+low)/2;
	} while(high>1+low);
	if(r<cum[low])
		return(low);
	else
		return(high);
#else
	for(uint32_t i=0;i<len;i++)
		if(cum[i]>=r)
			return(i);
#endif
	return(~0u);
}

