
#ifndef UTILITY_H_
#define UTILITY_H_

// a pair representing a constraint pair
typedef std::pair<char, char> p_constr;

/* Utility functions */
inline unsigned fast_atoi(const char* str) {
	unsigned val = 0;
    while (*str) val = val*10 + (*str++ - '0');
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

inline void printX(int x) {
	std::cout << x << std::endl;
}

void read_constraints(std::vector<p_constr>& l) {
	std::string s;
	while (getline(std::cin, s)) {
		// TODO: this works cause we have only symbols
		// we need a refactor in order to use tokens
		p_constr t(s.substr(0,1)[0], s.substr(2,3)[0]);
		l.push_back(t);
	}
}


/* USED FOR HILL CLIMBING HEURISTIC */

inline int index_of(const unsigned& x, const unsigned* s, const size_t& n) {
	for (size_t i = 0; i < n; ++i)
		if (s[i] == x)
			return i;
	return -1;
}

// Fisher-Yates shuffle
void shuffle(unsigned* a, size_t s) {
	size_t j = 0;
	unsigned short t = 0;
	for (size_t i = s-1; i > 0; i--) {
		j = rand() % (i+1);
		t = a[j]; a[j] = a[i]; a[i] = t;
	}
}

bool isValid(unsigned* s, int s_size, int p) {
	int currentMin = INT_MAX-1;
	int lastMin = INT_MAX;

	for (int i = 0; i < s_size; i++) {
		if (currentMin > s[i])
			currentMin = s[i];

		if (i >= p)
			if (lastMin > currentMin)
				return false;

		if ((i % p) == (p-1)) {
			lastMin = currentMin;
			currentMin = INT_MAX;
		}
	}
	return true;
}

#endif /* UTILITY_H_ */
