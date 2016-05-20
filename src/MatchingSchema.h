#ifndef MATCHINGSCHEMA_H_
#define MATCHINGSCHEMA_H_

template <typename T>
struct matching_schema {
	size_t n, m;
	size_t p1, p2;
	T** ms;

	matching_schema(size_t a, size_t b, size_t c, size_t d, const T x) : n(a), m(b), p1(c), p2(d) {
		ms = new T*[n];
		for (size_t i = 0; i < n; ++i) ms[i] = new T[m];
		init(x);
	}

	~matching_schema() {
		for (size_t i = 0; i < n; ++i)
			delete[] ms[i];
		delete[] ms;
	}

	void init(T b) {
		for (size_t i = 0; i < n; ++i)
			for (size_t j = 0; j < m; ++j)
				ms[i][j] = b;
	}

	void set_identity(std::string& sigma1, std::string& sigma2, T value) {
		for (size_t i = 0; i < n; ++i)
			for (size_t j = 0; j < m; ++j)
				if (sigma1[i] == sigma2[j])
					ms[i][j] = value;
	}

	void set_general(const std::string& sigma1, const std::string& sigma2, T value) {
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				if (p1 !=0 && p2 !=0 && (i/p2 == j/p1))
					ms[i][j] = value;
				else
					ms[i][j] = !value;
	}

	void print_matching_schema(const std::string& sigma1, const std::string& sigma2) {
		std::cout << "\t";
		for (int j = 0; j < m; ++j)
			std::cout << sigma2[j] << "\t";
		std::cout << std::endl;
		for (int i = 0; i < n; ++i) {
			std::cout << sigma1[i] << "\t";
			for (int j = 0; j < m; ++j)
				std::cout << ms[i][j] << "\t";
			std::cout << std::endl;
		}
	}
};



#endif /* MATCHINGSCHEMA_H_ */
