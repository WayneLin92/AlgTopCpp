// tmp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "mymath.h"
#include "myparser.h"
#include <fstream>
#include <sstream>

int main(int argc, char** argv)
{
	/*test();
	return 0;*/
	std::ifstream file;
	if (argc > 1) {
		char* filename = argv[1];
		file.open(filename);
	}
	else
		file.open("relations-7.txt");

	if (file.fail()) {
		std::cout << "File can not open." << std::endl;
		return 0;
	}
	std::string line;
	std::vector<std::string> gen_names;
	std::vector<int> gen_degs;
	Relations rel_gens;
	file >> "gen_names:";
	load_vector(gen_names, file, "[", ",", "]", load_py_str);
	file >> std::ws >> "gen_degs:";

	load_vector(gen_degs, file, "[", ",", "]");
	file >> std::ws >> "relations:" >> rel_gens;
	file.close();
	std::cout << "gen_names: ";
	dump_vector(gen_names, std::cout, "[", ", ", "]\n");
	std::cout << "gen_degs: " << gen_degs << '\n';
	std::cout << "relations: " << rel_gens << '\n';

	Relations relations;
	for (Poly& rel : rel_gens) {
		sort(rel);
		std::cout << "adding " << rel << '\n';
		add_rel(rel, relations, &gen_degs);
	}
	/*for (Poly& rel : relations)
		std::cout << rel << '\n';*/
	std::cout << "rels size:" << relations.size() << '\n';

	return 0;
}