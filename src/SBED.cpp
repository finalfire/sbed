#include <iostream>
#include <string>
#include <bitset>
#include <algorithm>
#include <map>
#include <vector>
#include "MatchingSchema.h"
#include "Alignment.h"
#include "Utility.h"
#include "EditDistance.h"

/* Definitions */
#define endl '\n'
#define _ASCII_START 32
#define _DEBUG false
#define _MAX_ARGS_NUM 4

/* Consts */
const unsigned _ASCII_LEN = 255 - 0;

/* Functions */
void extract_sigma(std::string&, std::string&);
void define_mapping(std::string&, std::map<char,int>&);
void print_alignment(const Alignment<int>&, const std::string&, const std::string&, const matching_schema<bool>&, int&, bool);
// for edit distance

/* Solver functions */
int bruteforce(const unsigned*, const unsigned*, const size_t&, const size_t&, const unsigned*, const unsigned*, const size_t&, const size_t&, matching_schema<bool>&);
int hill_climbing(const unsigned*, const unsigned*, const size_t&, const size_t&, const unsigned*, const unsigned*, const size_t&, const size_t&, const size_t&, matching_schema<bool>&);


/* main */
int main(int argc, char *argv[]) {
	std::ios_base::sync_with_stdio(false);

	// arguments: p1 p2 [y|(* \ {y})] [true|false]
	size_t p1 = argc > 2 ? fast_atoi(argv[1]) : 1;		// default = 1
	size_t p2 = argc > 2 ? fast_atoi(argv[2]) : 1;		// default = 1

	bool has_constraints = false;						// if there are constraints
	if (argc > 3) has_constraints = strcmp(argv[3],"1") || strcmp(argv[3], "y");
	/* if default_constraints_mode is true:
	 * 		constraints in input mean ONLY THESE MATCH
	 * otherwise:
	 * 		constraints in input mean ALL EXCEPT THESE MATCH
	 */
	bool default_constraints_mode = false;				// the semantic of the constraints
	//if (has_constraints) default_constraints_mode = strcmp(argv[4], "t");


	/* (1) read strings, extract sigmas and constraints too */
	std::string s1, s2;
	std::string sigma1(""), sigma2("");

	read_stdin(s1, s2);
	extract_sigma(s1, sigma1);
	extract_sigma(s2, sigma2);

	size_t sigma1l = sigma1.size();
	size_t sigma2l = sigma2.size();
	size_t s1l = s1.size();
	size_t s2l = s2.size();

	std::vector<p_constr> constraints;
	if (has_constraints)
		read_constraints(constraints);


	/* define the mapping from char -> int */
	std::map<char, int> map1;
	std::map<char, int> map2;

	define_mapping(sigma1, map1);
	define_mapping(sigma2, map2);

	/* integer representations of strings and sigmas */
	unsigned* s1i = new unsigned[s1l];
	unsigned* s2i = new unsigned[s2l];
	unsigned* sigma1i = new unsigned[sigma1l];
	unsigned* sigma2i = new unsigned[sigma2l];

	for (size_t i = 0; i < sigma1l; ++i) sigma1i[i] = i;		// sigma1
	for (size_t i = 0; i < sigma2l; ++i) sigma2i[i] = i;		// sigma2
	for (size_t i = 0; i < s1l; ++i) s1i[i] = map1[s1[i]];		// s1
	for (size_t i = 0; i < s2l; ++i) s2i[i] = map2[s2[i]];		// s2


	if (_DEBUG) {
		std::cout << "strings: " << s1 << ", " << s2 << endl;
		std::cout << "sigmas: " << sigma1 << ", " << sigma2 << endl;
		std::cout << "int rep (s1): ";
		for (size_t i = 0; i < s1l; ++i) std::cout << s1i[i]; std::cout << endl;
		std::cout << "int rep (s2): ";
		for (size_t i = 0; i < s2l; ++i) std::cout << s2i[i]; std::cout << endl;
		std::cout << "int rep (sigma1): ";
		for (size_t i = 0; i < sigma1l; ++i) std::cout << sigma1i[i]; std::cout << endl;
		std::cout << "int rep (sigma2): ";
		for (size_t i = 0; i < sigma2l; ++i) std::cout << sigma2i[i]; std::cout << endl;
	}


	/* identity (classical) matching schema */
	matching_schema<bool> ms(sigma1l, sigma2l, p1, p2, true, default_constraints_mode);
	ms.set_general(sigma1, sigma2, false);

	if (has_constraints)
		ms.set_constraints(map1, map2, constraints, !default_constraints_mode);


	if (_DEBUG) {
		if (has_constraints) {
			std::cout << "Constraints: ";
			for(std::vector<p_constr>::size_type i = 0; i < constraints.size(); ++i)
				std::cout << constraints[i].first << ", " << constraints[i].second << endl;
		}
		ms.print_matching_schema(sigma1, sigma2);
	}


	//int d = bruteforce(s1i, s2i, s1l, s2l, sigma1i, sigma2i, sigma1l, sigma2l, ms);
	int d = hill_climbing(s1i, s2i, s1l, s2l, sigma1i, sigma2i, sigma1l, sigma2l, p1, ms);
	std::cout << d << endl;

	//Alignment<int> a = compute_alignment(s1i, s2i, s1l, s2l, ms);
	//std::cout << distance_from_alignment(a, sigma1, sigma2, ms, false) << endl;


	return 0;
}


/* === */


/** SOLVER FUNCTIONS **/
int bruteforce(const unsigned* s1, const unsigned* s2,  const size_t& s1l, const size_t& s2l, const unsigned* sig1, const unsigned* sig2, const size_t& sig1l, const size_t& sig2l, matching_schema<bool>& m) {

	int distance = edit_distance_matching_schema(s1, s2, s1l, s2l, m);
	int current = distance;

	unsigned* perm1 = new unsigned[sig1l]; for (unsigned i = 0; i < sig1l; ++i) perm1[i] = i;
	unsigned* perm2 = new unsigned[sig2l]; for (unsigned i = 0; i < sig2l; ++i) perm2[i] = i;
	unsigned* best_perm1 = new unsigned[sig1l]; for (unsigned i = 0; i < sig1l; ++i) perm1[i] = i;
	unsigned* best_perm2 = new unsigned[sig2l]; for (unsigned i = 0; i < sig2l; ++i) perm2[i] = i;

	//do {
		do {
			current = edit_distance_matching_schema_enhanced(s1, s2, s1l, s2l, perm1, perm2, sig1l, sig2l, m);

			if (current < distance) {
				distance = current;
				std::copy(perm1, perm1 + sig1l, best_perm1);
				std::copy(perm2, perm2 + sig2l, best_perm2);
			}
		} while(std::next_permutation(perm2, perm2+sig2l));
	//} while(std::next_permutation(perm1, perm1+sig1l));

	if (_DEBUG) {
		std::cout << "best matching schema (brute force):" << endl;
		for (size_t i = 0; i < sig1l; ++i) std::cout << best_perm1[i]; std::cout << endl;
		for (size_t i = 0; i < sig2l; ++i) std::cout << best_perm2[i]; std::cout << endl;
	}

	return distance;
}

int hill_climbing(const unsigned* s1, const unsigned* s2,  const size_t& s1l, const size_t& s2l, const unsigned* sig1, const unsigned* sig2, const size_t& sig1l, const size_t& sig2l, const size_t& p1, matching_schema<bool>& m) {

	int d = edit_distance_matching_schema(s1, s2, s1l, s2l, m);
	int minDist = d;
	int minMinDist = minDist;

	// for swap
	unsigned* sigma1_o = new unsigned[sig1l]; for (unsigned i = 0; i < sig1l; ++i) sigma1_o[i] = i;
	unsigned* sigma2_o = new unsigned[sig2l]; for (unsigned i = 0; i < sig2l; ++i) sigma2_o[i] = i;
	unsigned* sigma1_t = new unsigned[sig1l]; for (unsigned i = 0; i < sig1l; ++i) sigma1_t[i] = i;
	unsigned* sigma2_t = new unsigned[sig2l]; for (unsigned i = 0; i < sig2l; ++i) sigma2_t[i] = i;
	// for fixpoint on min
	unsigned* sigma1_min = new unsigned[sig1l]; for (unsigned i = 0; i < sig1l; ++i) sigma1_min[i] = i;
	unsigned* sigma2_min = new unsigned[sig2l]; for (unsigned i = 0; i < sig2l; ++i) sigma2_min[i] = i;
	unsigned* sigma1_min_min = new unsigned[sig1l]; for (unsigned i = 0; i < sig2l; ++i) sigma2_min[i] = i;
	unsigned* sigma2_min_min = new unsigned[sig2l]; for (unsigned i = 0; i < sig2l; ++i) sigma2_min_min[i] = i;

	size_t attempts = 1, shuffle_tries = 10;
	unsigned temp = 0, tries = 0, k_shuffle = 0;


	bool improved = true;
	while (improved) {
		improved = false;

		//if (!isValid(sigma1_o, sgl1, this->mped->getP1()))
		//	cout << "! not valid" << endl;


		for (size_t ip = 0; ip < sig1l; ip++) {
			for (size_t jp = ip; jp < sig1l; jp++) {

				// here comes the swap for sigma1
				std::copy(sigma1_t, sigma1_t + sig1l, sigma1_o);
				temp = sigma1_o[ip]; sigma1_o[ip] = sigma1_o[jp]; sigma1_o[jp] = temp;

				if (isValid(sigma1_o, sig1l, p1)) {

					for (size_t ipp = 0; ipp < sig2l; ipp++) {
						for (size_t jpp = ipp; jpp < sig2l; jpp++) {

							// here comes the swap for sigma2
							std::copy(sigma2_t, sigma2_t + sig2l, sigma2_o);
							temp = sigma2_o[ipp]; sigma2_o[ipp] = sigma2_o[jpp]; sigma2_o[jpp] = temp;

							d = edit_distance_matching_schema_enhanced(s1, s2, s1l, s2l, sigma1_o, sigma2_o, sig1l, sig2l, m);

							if (d < minDist) {
								minDist = d;

								improved = true;
								std::copy(sigma1_o, sigma1_o + sig1l, sigma1_min);
								std::copy(sigma2_o, sigma2_o + sig2l, sigma2_min);
							}

						}

						std::copy(sigma2_t, sigma2_t + sig2l, sigma2_o);
					}

				}
			}
		}

		if (improved) {
			// copy sigmaMin to sigmaOrig
			std::copy(sigma1_min, sigma1_min + sig1l, sigma1_o);
			std::copy(sigma2_min, sigma2_min + sig2l, sigma2_o);
			// copy sigmaOrig to sigmaTmp
			std::copy(sigma1_o, sigma1_o + sig1l, sigma1_t);
			std::copy(sigma2_o, sigma2_o + sig2l, sigma2_t);
		} else {
			if (minDist < minMinDist) {
				minMinDist = minDist;

				// copy sigmaMin to sigmaMinMin
				std::copy(sigma1_min, sigma1_min + sig1l, sigma1_min_min);
				std::copy(sigma2_min, sigma2_min + sig2l, sigma2_min_min);

				improved = true;
				tries = 0;
			}

			if (tries < attempts) {
				improved = true;
				tries++;

				// random swap
				// for sigma1, we try _SHUFFLE_TRIES times, then if is still not valid, we retry with the original one
				for (k_shuffle = 0; k_shuffle < shuffle_tries && !isValid(sigma1_t, sig1l, p1); ++k_shuffle)
					shuffle(sigma1_t, sig1l);
				if (k_shuffle == shuffle_tries) std::copy(sigma1_o, sigma1_o + sig1l, sigma1_t);
				// no constraints on the shuffle for sigma2
				shuffle(sigma2_t, sig2l);

				std::copy(sigma1_t, sigma1_t + sig1l, sigma1_o);
				std::copy(sigma2_t, sigma2_t + sig2l, sigma2_o);

				minDist = edit_distance_matching_schema_enhanced(s1, s2, s1l, s2l, sigma1_o, sigma2_o, sig1l, sig2l, m);
			} //else {
				// here's the final step: improved is false so we don't continue

				/* here we return the optimal matching schema if we need to */
				/*
				this->computed_sigma1 = new unsigned short[sgl1];
				copy(sigma1_min_min, sigma1_min_min + sgl1, this->computed_sigma1);
				this->computed_sigma2 = new unsigned short[sgl2];
				copy(sigma2_min_min, sigma2_min_min + sgl2, this->computed_sigma2);
				cout << ">> MPED SCHEMA: ";
				for (size_t i = 0; i < sgl1; i++)
					cout << this->mped->get_Sigma1()[computed_sigma1[i]];
				cout << " - ";
				for (size_t i = 0; i < sgl2; i++)
					cout << this->mped->get_Sigma2()[computed_sigma2[i]];
				cout << endl;
				this->result = min_dist;*/
			//}
		}
	}

	return minMinDist;
}


/**
 * extract_sigmas works with a bitset of length _ASCII_LEN = 255 - _ASCII_START
 * in order to index the symbols we are going to read.
 */
void extract_sigma(std::string& s, std::string& e) {
	std::bitset<_ASCII_LEN> symbols;

	for (size_t i = 0; i < s.size(); ++i)
		symbols[(int) s[i]] = (symbols[(int) s[i]] || 1);

	for (size_t i = 0; i < _ASCII_LEN; ++i)
		if (symbols[i])
			e += (char) i;
}

void define_mapping(std::string& s, std::map<char,int>& m) {
	for (size_t i = 0; i < s.size(); ++i)
		m.insert(std::pair<char,int>(s[i], i));
}

void print_alignment(const Alignment<int>& alignment, const std::string& sigma1, const std::string& sigma2, const matching_schema<bool>& m, int& distance, bool is_identity) {
	// print first string
	for (size_t i = 0; alignment.a[i] != _END_ALIGNMENT; i++)
		std::cout << (alignment.a[i] != _GAP_FLAG ? sigma1[alignment.a[i]] : '-');
	std::cout << endl;

	// print second string
	for (size_t i = 0; alignment.b[i] != _END_ALIGNMENT; i++)
		std::cout << (alignment.b[i] != _GAP_FLAG ? sigma2[alignment.b[i]] : '-');
	std::cout << endl;

	// print the alignment (a string containing * or space)
	size_t index = 0;
	size_t temp = 0;
	for (; alignment.a[index] != _END_ALIGNMENT; index++) {
		// if both chars are not gaps
		if (alignment.a[index] != _GAP_FLAG && alignment.b[index] != _GAP_FLAG) {
			if (	!m.ms[alignment.a[index]][alignment.b[index]] || 									// if they are in match or
					(is_identity && (sigma1[alignment.a[index]] == sigma2[alignment.b[index]]))		) {	// the self_identity is active and they are the same char
				std::cout << '*';
				temp++;
			} else
				std::cout << ' ';
		} else
			std::cout << ' ';
	}
	std::cout << endl;

	// here we return the new edit distance including the self_identity case
	distance = index - temp;
}

