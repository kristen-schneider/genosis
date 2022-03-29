#ifndef preprocessing_h    // To make sure you don't declare the function more than once by including the header multi>
#define preprocessing_h
//declare creating db and refs function
#include <vector>
#include <map>
#include "node.h"
#include <libpmemkv.hpp>
#include <libpmemobj++/container/concurrent_hash_map.hpp>
#include <libpmemobj++/p.hpp>
#include <libpmemobj++/persistent_ptr.hpp>
#include <libpmemobj++/pool.hpp>

std::map <int, std::vector<std::vector<uint8_t>> >  createRefs(int n , string filepath );
std::vector<std::vector<uint8_t>> createRefs_spike(int n , string filepath );
std::vector<std::vector<uint8_t>> createdb(int n , std::vector<string> );
std::vector<std::vector<uint8_t>> createdb_OHE(int n , const char* filepath );
int dfs2pmem(pmem::obj::pool<root_tree_pmem> *pop1, Node* startNode, int _max_lvl, int N);
std::pair< int, Node*> PmemCreatTree(pmem::obj::pool<root_tree_pmem> *pop1);

void printChildren( Node* parent);

#endif