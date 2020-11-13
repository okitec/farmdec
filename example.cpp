// C++ farmdec example
//
// g++ -Wall example.cpp -L -l:libfarmdec.a

#include <cstdint>
#include <iostream>

#include "farmdec.h"

int main() {
	uint32_t in = 0xd10083ff; // sub sp, sp, #32
	farmdec::Inst out;
	int n = decode(&in, 1, &out);
	std::cout << n;
	return 0;
}
