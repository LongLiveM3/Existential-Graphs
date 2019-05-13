// Copyright 2019 Luca Istrate, Danut Matei
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <set>
#include <map>
#include <utility>
#include <cassert>
#include "./aegraph.h"

std::string strip(std::string s) {
    // deletes whitespace from the beginning and end of the string
    s.erase(0, s.find_first_not_of(" \n\r\t"));
    s.erase(s.find_last_not_of(" \n\r\t")+1);
    return s;
}

std::pair<std::string, std::string> split_first(std::string s,
    char delimiter = ',') {
    // returns a pair that contains: <first_cut, rest_of_graph>

    int numOpen = 0;

    int len = s.size();
    for (int i = 0; i < len; i++) {
        char c = s[i];
        if (c == delimiter && numOpen == 0)
            return std::make_pair(
                    strip(s.substr(0, i)), strip(s.substr(i+1, len)));
        if (c == '[')
            numOpen += 1;
        if (c == ']')
            numOpen -= 1;
    }

    return {strip(s), std::string()};
}


std::vector<std::string> split_level(std::string s, char delimiter = ',') {
    // splits 's' into separate entities (atoms, subgraphs)

    std::vector<std::string> result;
    auto snd = s;
    while (true) {
        auto p = split_first(snd, delimiter);
        auto fst = p.first;
        snd = p.second;

        result.push_back(fst);

        if (snd.empty())
            return result;
    }
}


int AEGraph::num_subgraphs() const {
    return subgraphs.size();
}


int AEGraph::num_atoms() const {
    return atoms.size();
}


int AEGraph::size() const {
    return num_atoms() + num_subgraphs();
}


bool AEGraph::operator<(const AEGraph& other) const {
    return this->repr() < other.repr();
}

bool AEGraph::operator==(const AEGraph& other) const {
    return this->repr() == other.repr();
}

bool AEGraph::operator!=(const AEGraph& other) const {
    return this->repr() != other.repr();
}

AEGraph AEGraph::operator[](const int index) const {
    // offers an easier way of accessing the nested graphs
    if (index < num_subgraphs()) {
        return subgraphs[index];
    }

    if (index < num_subgraphs() + num_atoms()) {
        return AEGraph('(' + atoms[index - num_subgraphs()] + ')');
    }

    return AEGraph("()");
}

std::ostream& operator<<(std::ostream &out, const AEGraph &g) {
    out << g.repr();
    return out;
}

AEGraph::AEGraph(std::string representation) {
    // constructor that creates an AEGraph structure from a
    // serialized representation
    char left_sep = representation[0];
    char right_sep = representation[representation.size() - 1];

    assert((left_sep == '(' && right_sep == ')')
        || (left_sep == '[' && right_sep == ']'));

    // if the left separator is '(' then the AEGraph is the entire
    // sheet of assertion
    if (left_sep == '(') {
        is_SA = true;
    } else {
        is_SA = false;
    }

    // eliminate the first pair of [] or ()
    representation = representation.substr(1, representation.size() - 2);

    // split the graph into separate elements
    auto v = split_level(representation);
    // add them to the corresponding vector
    for (auto s : v) {
        if (s[0] != '[') {
            atoms.push_back(s);
        } else {
            subgraphs.push_back(AEGraph(s));
        }
    }

    // also internally sort the new graph
    this->sort();
}

std::string AEGraph::repr() const {
    // returns the serialized representation of the AEGraph

    std::string left, right;
    if (is_SA) {
        left = '(';
        right = ')';
    } else {
        left = '[';
        right = ']';
    }

    std::string result;
    for (auto subgraph : subgraphs) {
        result += subgraph.repr() + ", ";
    }

    int len = atoms.size();
    if (len != 0) {
        for (int i = 0; i < len - 1; i++) {
            result += atoms[i] + ", ";
        }
        result += atoms[len - 1];
    } else {
        if (subgraphs.size() != 0)
            return left + result.substr(0, result.size() - 2) + right;
    }

    return left + result + right;
}


void AEGraph::sort() {
    std::sort(atoms.begin(), atoms.end());

    for (auto& sg : subgraphs) {
        sg.sort();
    }

    std::sort(subgraphs.begin(), subgraphs.end());
}

bool AEGraph::contains(const std::string other) const {
    // checks if an atom is in a graph
    if (find(atoms.begin(), atoms.end(), other) != atoms.end())
        return true;

    for (const auto& sg : subgraphs)
        if (sg.contains(other))
            return true;

    return false;
}

bool AEGraph::contains(const AEGraph& other) const {
    // checks if a subgraph is in a graph
    if (find(subgraphs.begin(), subgraphs.end(), other) != subgraphs.end())
        return true;

    for (const auto& sg : subgraphs)
        if (sg.contains(other))
            return true;

    return false;
}

std::vector<std::vector<int>> AEGraph::get_paths_to(const std::string other)
    const {
    // returns all paths in the tree that lead to an atom like <other>
    std::vector<std::vector<int>> paths;

    int len_atoms = num_atoms();
    int len_subgraphs = num_subgraphs();

    for (int i = 0; i < len_atoms; i++) {
        if (atoms[i] == other && size() > 1) {
            paths.push_back({i + len_subgraphs});
        }
    }

    for (int i = 0; i < len_subgraphs; i++) {
        if (subgraphs[i].contains(other)) {
            auto r = subgraphs[i].get_paths_to(other);
            for (auto& v : r)
                v.insert(v.begin(), i);
            copy(r.begin(), r.end(), back_inserter(paths));
        }
    }

    return paths;
}

std::vector<std::vector<int>> AEGraph::get_paths_to(const AEGraph& other)
    const {
    // returns all paths in the tree that lead to a subgraph like <other>
    std::vector<std::vector<int>> paths;
    int len_subgraphs = num_subgraphs();

    for (int i = 0; i < len_subgraphs; i++) {
        if (subgraphs[i] == other && size() > 1) {
            paths.push_back({i});
        } else {
            auto r = subgraphs[i].get_paths_to(other);
            for (auto& v : r)
                v.insert(v.begin(), i);
            copy(r.begin(), r.end(), back_inserter(paths));
        }
    }
    return paths;
}

void possibleDCR(AEGraph g, int anterior,
	std::vector<std::vector<int>>& res, std::vector<int>& path){
	// daca nodul trecut are un subgraf + daca nodul curent are un subgraf
	// si nodul curent nu are atomi
	if(anterior >= 1 && g.num_subgraphs() == 1 && g.num_atoms() == 0){
		res.push_back(path);
	}
	// daca nodul curent are mai multe subgrafuri
	if(g.num_subgraphs() >= 1){
		for(int i = 0; i < g.num_subgraphs(); ++i){
			// se adauga la path
			path.push_back(i);
			possibleDCR(g.subgraphs[i], g.num_subgraphs(), res, path);
			path.pop_back();
		}
	}
}

std::vector<std::vector<int>> AEGraph::possible_double_cuts() const {
    // 10p
    std::vector<std::vector<int>> res = {};
    std::vector<int> path = {};
    AEGraph g(repr());
    possibleDCR(g, -1, res, path);
    return res;
}

void doubleCutR(AEGraph& g, std::vector<int>& where, int& ok){
		for(int i = 0; i < g.num_subgraphs(); ++i){
			// daca drumul e acelasi cu cel dat in where
			if(where.front() == i){
		    	if(where.size() > 1){
		    		where.erase(where.begin());
		    		doubleCutR(g.subgraphs[i], where, ok);
		    	}
		  		if(ok && where.size() == 1){
		  			// eliminarea primului set de paranteze patrate(prima taietura)
			  		g.subgraphs[where[0]] =
			  			g.subgraphs[where[0]].subgraphs[0];
			  		// trecerea prin toti atomii nivelului si copierea lor
					for(auto a : g.subgraphs[where[0]].atoms){
						g.atoms.push_back(a);
					}
					g.subgraphs[where[0]].atoms.clear();
					// eliminarea celui de-al doilea set de paranteze(a doua taietura)
					g.subgraphs.erase(g.subgraphs.begin() + where[0]);
					ok = 0;
					doubleCutR(g, where, ok);
				}
			}
	    }
}

AEGraph AEGraph::double_cut(std::vector<int> where) const {
    // 10p
    AEGraph g(repr());
    int ok = 1;
    doubleCutR(g, where, ok);
    return g;
}

void possibleER(AEGraph g, std::vector<std::vector<int>>& res,
	std::vector<int>& path, int ok){
	// verificarea conditiilor de nivel (nivel initial sau nivel par
	// cu mai multi frati in taietura)
    if (path.size() == 1){
      	res.push_back(path);
  	} else if (path.size() % 2 == 1 && ok){
      	res.push_back(path);
    }
    ok = 0;
    if (g.size() >= 1){
    	// parcurgerea subgrafurilor si atomilor
  		for(int i = 0; i < g.size(); ++i){
  			path.push_back(i);
  			// verificarea conditiei de mai multi frati in taietura
	        if (g.size() > 1)
	        	ok = 1;
	        // in cazul pentru subgrafuri se apeleaza recursiv functia
	        if (i < g.num_subgraphs()){
	  			possibleER(g[i], res, path, ok);
	  		// in cazul pentru atomi sunt verificate conditiile de nivel
	        } else if ((path.size() % 2 == 1 && ok) || path.size() == 1){
	          	res.push_back(path);
	        }
  			path.pop_back();
  	  	}
    }
}

std::vector<std::vector<int>> AEGraph::possible_erasures(int level) const {
    // 10p
    std::vector<std::vector<int>> res = {};
    std::vector<int> path = {};
    int ok = level;
    AEGraph g(repr());
    if (size() != 0) {
    	possibleER(g, res, path, ok);
    }
	return res;
}

void eraseR(AEGraph& g, std::vector<int>& where, int& ok){
	for(int i = 0; i < g.size(); ++i){
		// daca drumul e acelasi cu cel dat in where
		if(where.front() == i){
			if(where.size() > 1){
				where.erase(where.begin());
				eraseR(g.subgraphs[i], where, ok);
			}
	      	if(ok && where.size() == 1){
	      		// eliminarea subgrafului daca path-ul din where e subgraf
	        	if (i < g.num_subgraphs()){
	          		g.subgraphs.erase(g.subgraphs.begin() + where[0]);
	        	} else {  // eliminarea atomului daca path-ul din where e atom
	          		g.atoms.erase(g.atoms.begin()
	          			+ where[0] - g.num_subgraphs());
	        	}
	        	ok = 0;
	      	}
		}
	}
}

AEGraph AEGraph::erase(std::vector<int> where) const {
    // 10p
    AEGraph g(repr());
    int ok = 1;
    eraseR(g, where, ok);
    return g;
}

std::vector<std::vector<int>> AEGraph::possible_deiterations() const {
    // 20p
	std::vector<std::vector<int>> res = {};
    std::vector<std::vector<int>> path;
	std::string at;
    AEGraph g(repr());
    // trecerea prin toti atomii si toate subgrafurile nodului curent
	int i = 0;
	while(i < g.size()){
		// daca e subgraf
		if(i < g.num_subgraphs()){
			path = g.get_paths_to(g[i]);
		} else {  // daca e atom
			at = g[i].repr();
			at = at[1];
			path = g.get_paths_to(at);
		}
		// copierea inafara deiterarii cautate
		for(auto j : path){
			if(j[0] != i){
				res.push_back(j);
			}
		}
		++i;
	}
    return res;
}

void deiterateR(AEGraph g1, AEGraph& g2, std::vector<int>& where,
  std::vector<std::vector<int>>& path, int& ok){
  	// se parcurge recursiv graful pana la aflarea grafului dat de where
	for(int i = 0; i < g1.size(); ++i){
		if(where.front() == i){
			if(where.size() > 1){
				where.erase(where.begin());
				deiterateR(g1.subgraphs[i], g2, where, path, ok);
			}
			// salvarea drumurilor pentru deiterare pana la graful
	    	// aflat anterior
	    	if (ok && where.size() == 1){
	    		// cazul pentru subgrafuri
	        	if (i < g1.num_subgraphs()){
	          		path = g2.get_paths_to(g1.subgraphs[where[0]]);
	        	// cazul pentru atomi
	        	} else {
	    			path = g2.get_paths_to(g1.atoms[where[0] -
	            	g1.num_subgraphs()]);
	        	}
	        	ok = 0;
	      	}
    	}
	}
}

AEGraph AEGraph::deiterate(std::vector<int> where) const {
    // 10p
    // graful initial prin care se ajunge la graful dat de where
    AEGraph g1(repr());
    AEGraph g2(repr());  // graful final dupa stergeri
    int ok = 1;
    std::vector<std::vector<int>> path;  // partial path
    deiterateR(g1, g2, where, path, ok);
    // se sterg subgrafurile si grafurile date de path-urile
    // salvate in vectorul de drumuri
    for (unsigned int i = 0; i < path.size(); i++){
      g2 = g1.erase(path[i]);
    }
	return g2;
}
