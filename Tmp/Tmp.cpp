// Tmp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "../algtop/mymath.h"

int tmp()
{
	array gen_degs = { 1, 1, 1, 1 }; /* x, y, z, w */
	Poly rel0 = { {{0, 4}} }; /* x^4 */
	Poly rel1 = { {{0, 2}, {1, 4}} }; /* x^2y^4 */
	Poly rel2 = { {{0, 1}, {2, 2}}, {{0, 2}, {1, 1}} }; /* x^2y + xz^2 */
	Poly rel3 = { {{3, 4}}, {{2, 4}} }; /* z^4 + w^4 */
	Poly x = { {{0, 1}} };
	Poly z = { {{2, 1}} };
	Poly w = { {{3, 1}} };
	Poly1d gb;
	
	add_rels(gb, { rel0, rel1, rel2, rel3 }, gen_degs, -1);

	Mon mon = { {0, 1}, {1, 2} };
	std::cout << FnGetDeg{ gen_degs }(mon) << '\n';
	return 0;
}

int main()
{
	array a = { 0, 1, 2, 3 };
	for (auto p1 = a.begin(); p1 != a.end(); ++p1)
		for (auto p2 = p1 + 1; p2 != a.end(); ++p2)
			std::cout << *p1 << ' ' << *p2 << '\n';
	return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
