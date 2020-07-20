#ifndef DATABSE_H
#define DATABSE_H

#include "sqlite3/sqlite3.h"
#include <iostream>
#include <vector>

void execute_cmd(sqlite3* conn, const char* cmd);
std::vector<int> str_to_mon(const char* str_mon);
inline const char* sqlite3_column_str(sqlite3_stmt* stmt, int iCol) { return reinterpret_cast<const char*>(sqlite3_column_text(stmt, iCol)); }
/* num_mons[t] is the number of monomials up to degree t */
void load_mon_indices(sqlite3* conn, const char* table_name, std::vector<int>& num_mons);
std::string mon_to_str(const std::vector<int>& mon);
int main_generate_basis(int argc, char** argv);

int main_generate_diff(int argc, char** argv);

#endif /* DATABSE_H */