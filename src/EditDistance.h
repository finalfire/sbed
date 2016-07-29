#ifndef EDITDISTANCE_H_
#define EDITDISTANCE_H_

#include "Utility.h"
#include "Alignment.h"
#include "MatchingSchema.h"

const short _END_ALIGNMENT = -2;
const short _GAP_FLAG = -1;

unsigned edit_distance_matching_schema_enhanced(const unsigned* a, const unsigned* b, const size_t& al, const size_t& bl,
		const unsigned* sg1, const unsigned* sg2, const size_t& sg1l, const size_t& sg2l, const matching_schema<bool>& m) {

	unsigned** d = new unsigned*[al+1];
	for (size_t i = 0; i < al + 1; ++i) d[i] = new unsigned[bl+1];

	// (1) first row and first column
	for (size_t i = 0; i < al+1; ++i) d[i][0] = i;
	for (size_t j = 0; j < bl+1; ++j) d[0][j] = j;

	// (2) fill the matrix
	for (size_t i = 1; i < al+1; ++i)
		for (size_t j = 1; j < bl+1; ++j) {
			d[i][j] = min(
					d[i-1][j] + 1,																		// deletion
					d[i][j-1] + 1,																		// insertion
					d[i-1][j-1] + 1 * m.ms[index_of(a[i-1], sg1, sg1l)][index_of(b[j-1], sg2, sg2l)]	// if in the matching schema there's a false, they match
			);
		}

	// (3) the computed distance
	unsigned my_dist = d[al][bl];

	for (size_t i = 0; i < al; ++i) delete[] d[i];
	delete[] d;

	return my_dist;
}

unsigned edit_distance_matching_schema(const unsigned* a, const unsigned* b, const size_t& al, const size_t& bl, const matching_schema<bool>& m) {

	unsigned** d = new unsigned*[al+1];
	for (size_t i = 0; i < al + 1; ++i) d[i] = new unsigned[bl+1];

	// (1) first row and first column
	for (size_t i = 0; i < al+1; ++i) d[i][0] = i;
	for (size_t j = 0; j < bl+1; ++j) d[0][j] = j;

	// (2) fill the matrix
	for (size_t i = 1; i < al+1; ++i)
		for (size_t j = 1; j < bl+1; ++j) {
			d[i][j] = min(
					d[i-1][j] + 1,							// deletion
					d[i][j-1] + 1,							// insertion
					d[i-1][j-1] + 1 * m.ms[a[i-1]][b[j-1]]	// if in the matching schema there's a false, they match
			);
		}

	// (3) the computed distance
	unsigned my_dist = d[al][bl];

	for (size_t i = 0; i < al; ++i) delete[] d[i];
	delete[] d;

	return my_dist;
}

Alignment<int> compute_alignment(const unsigned* a, const unsigned* b, const size_t& al, const size_t& bl, const matching_schema<bool>& m) {
	unsigned** d = new unsigned*[al+1];
	for (size_t i = 0; i < al + 1; ++i) d[i] = new unsigned[bl+1];
	char** path = new char*[al+1];
	for (size_t i = 0; i < al + 1; ++i) path[i] = new char[bl+1];

	// (1) first row and first column
	for (size_t i = 0; i <= al; ++i) { d[i][0] = i; path[i][0] = 'n'; }
	for (size_t j = 0; j <= bl; ++j) { d[0][j] = j; path[0][j] = 'o'; }

	// (2) fill the matrix
	for (size_t i = 1; i <= al; ++i) {
		for (size_t j = 1; j <= bl; ++j) {
			d[i][j] = min(
					d[i-1][j] + 1,								// deletion
					d[i][j-1] + 1,								// insertion
					d[i-1][j-1] + 1 * m.ms[a[i-1]][b[j-1]]		// if in the matching schema there's a false, they match
			);


			// (3) annotating the path for the backtrace
			if (d[i][j] == (d[i-1][j-1] + 1 * m.ms[a[i-1]][b[j-1]]))
				path[i][j] = 'd';
			else if (d[i][j] == d[i-1][j] + 1)
				path[i][j] = 'n';
			else
				path[i][j] = 'o';
		}
	}

	// (3) the computed distance
	unsigned my_dist = d[al][bl];

	// (4) here I start the backtrace
	unsigned maxl = (al > bl) ? al : bl;
	int x = 2 * maxl - 1;

	std::vector<int> all1(2*maxl, -1);
	std::vector<int> all2(2*maxl, -1);

	int i = al;
	int j = bl;

	while (i >= 0 || j >= 0) {
		if (i == 0 && j == 0) {
			break;
		} else {
			if (path[i][j] == 'n') {
				all1[x] = a[i-1];
				all2[x] = _GAP_FLAG;
				i = i-1;
				x--;
			} else if (path[i][j] == 'd') {
				all1[x] = a[i-1];
				all2[x] = b[j-1];
				i = i-1;
				j = j-1;
				x--;
			} else {
				all1[x] = _GAP_FLAG;
				all2[x] = b[j-1];
				j = j-1;
				x--;
			}
		}
	}

	// (5) here I define the alignment
	Alignment<int> alignment(my_dist, 2*maxl-x);

	size_t k = 0;
	for (k = x+1; k <= 2*maxl-1; k++)
		alignment.a[k-x-1] = all1[k];
	alignment.a[k-x-1] = -2;
	for (k = x+1; k <= 2*maxl-1; k++)
		alignment.b[k-x-1] = all2[k];
	alignment.b[k-x-1] = -2;


	// (6) deallocate everything
	for (size_t i = 0; i < al; ++i) delete[] path[i];
	delete[] path;
	for (size_t i = 0; i < al; ++i) delete[] d[i];
	delete[] d;

	return alignment;
}

int distance_from_alignment(const Alignment<int>& alignment, const std::string& sigma1, const std::string& sigma2, const matching_schema<bool>& m, bool is_identity) {
	// print the alignment (a string containing * or space)
	size_t index = 0;
	size_t temp = 0;
	for (; alignment.a[index] != _END_ALIGNMENT; index++) {
		// if both chars are not gaps
		if (alignment.a[index] != _GAP_FLAG && alignment.b[index] != _GAP_FLAG) {
			if (	!m.ms[alignment.a[index]][alignment.b[index]] || 									// if they are in match or
					(is_identity && (sigma1[alignment.a[index]] == sigma2[alignment.b[index]]))		) {	// the self_identity is active and they are the same char
				temp++;
			}
		}
	}

	// here we return the new edit distance including the self_identity case
	return index - temp;
}

#endif /* EDITDISTANCE_H_ */
