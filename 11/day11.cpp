#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <chrono>

using namespace std;

// Simple hash for 3-letter strings
uint16_t s3hash(string s) {
    return s[0] - 'a' + (s[1] - 'a') * 26 + (s[2] - 'a') * 26 * 26;
}

// Parse input file into a map of <input -> vector<outputs>>
const unordered_map<uint16_t, vector<uint16_t>> parse_connections(string filename) {
    unordered_map<uint16_t, vector<uint16_t>> connections;

    ifstream infile(filename);
    string line;
    while(getline(infile, line)) {
        istringstream iss(line);

        string source;
        getline(iss, source, ':');
        uint16_t source_idx = s3hash(source);

        connections[source_idx] = vector<uint16_t>();

        string dest;
        while (iss >> dest) {
            connections[source_idx].push_back(s3hash(dest));
        }

    }
    return connections;
}

uint64_t dfs_memo(
    uint16_t node, 
    const unordered_map<uint16_t, vector<uint16_t>>& connections,
    unordered_map<uint16_t, uint64_t>* paths_to_out 
) {
    if (paths_to_out->contains(node)) {
        return paths_to_out->at(node);
    }
    if (!connections.contains(node)) {
        paths_to_out->insert({node, 0});
        return 0;
    }
    uint64_t count = 0;
    for (auto dest : connections.at(node)) {
        count += dfs_memo(dest, connections, paths_to_out);
    }
    paths_to_out->insert({node, count});
    return count;
}

uint64_t count_paths(
    string start,
    string end,
    const unordered_map<uint16_t, vector<uint16_t>>& connections
) {
    // We can do better than a brute force search here.
    // We'll run a DFS until we find `end`. When we find `end`, we know that the nodes on the path to `end` can all reach `end`. 
    // When we fully explore a node (on return from its DFS), we'll record the number of times that it reached `out`.
    // When we reach a node for which we've recorded, then we can avoid re-exploring that tail and simply inherit its count.
    // Basically, we're memoizing fully evaluated tails. 
    unordered_map<uint16_t, uint64_t> paths_to_out;    // The number of times a given node can reach `out`
    paths_to_out[s3hash(end)] = 1;    // Our base condition - there is one way to reach `out`, because we're already there :)

    return dfs_memo(s3hash(start), connections, &paths_to_out);
}

uint64_t part1() {
    const auto connections = parse_connections("input1.txt");
    return count_paths("you", "out", connections);
}

uint64_t part2() {
    const auto connections = parse_connections("input2.txt");

    // We can reuse our solution for part 1.
    // The number of paths from `svr` to `out` that pass through both `fft` and `dac` is:
    // number of paths from `svr` -> `fft` -> `dac` -> `out`
    //         + paths from `svr` -> `dac` -> `fft` -> `out`

    // ==        paths from (`svr` -> `fft`) * (`fft` -> `dac`) * (`dac -> out`)
    //         + paths from (`svr` -> `dac`) * (`dac` -> `fft`) * (`fft` -> out`)
    return count_paths("svr", "fft", connections) * count_paths("fft", "dac", connections) * count_paths("dac", "out", connections)
        + count_paths("svr", "dac", connections) * count_paths("dac", "fft", connections) * count_paths("fft", "out", connections);
}

int main() {

    auto start = chrono::high_resolution_clock::now();
    auto count = part1();
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double, std::milli> ms_double = end - start;
    cout << "Part 1: " << count;
    cout << " in " << ms_double << " ms." << endl;

    start = chrono::high_resolution_clock::now();
    count = part2();
    end = chrono::high_resolution_clock::now();

    ms_double = end - start;
    cout << "Part 2: " << count;
    cout << " in " << ms_double << " ms." << endl;
}
