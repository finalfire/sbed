/*
 * Utility.h
 *
 *  Created on: 29/mag/2016
 *      Author: finalfire
 */

#ifndef UTILITY_H_
#define UTILITY_H_

// a pair representing a constraint pair
typedef std::pair<char, char> p_constr;

/* Utility functions */
inline unsigned fast_atoi(const char* str) {
	unsigned val = 0;
    while( *str ) val = val*10 + (*str++ - '0');
    return val;
}

inline const int min(const int a, const int b, const int c) {
	if ((a<=b) && (a<=c)) return a;
	else if ((b<=a) && (b<=c)) return b;
	else return c;
}

inline void read_stdin(std::string& a, std::string& b) {
	getline(std::cin, a);
	getline(std::cin, b);
}

template <typename T>
inline void print(T x) { std::cout << x << std::endl; }

void read_constraints(std::vector<p_constr>& l) {
	std::string s;
	while (getline(std::cin, s)) {
		p_constr t(s.substr(0,1)[0], s.substr(2,3)[0]);
		l.push_back(t);
	}
}



#endif /* UTILITY_H_ */
