/*
 * Utility.h
 *
 *  Created on: 29/mag/2016
 *      Author: finalfire
 */

#ifndef UTILITY_H_
#define UTILITY_H_

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
inline void print(T x) { std::cout << x << endl; }



#endif /* UTILITY_H_ */
