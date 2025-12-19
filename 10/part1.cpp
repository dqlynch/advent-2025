#include <vector>
#include <unordered_map>
#include <iostream>
#include <tuple>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#include <queue>
#include <stack>

using namespace std;
using mask16_t = uint16_t;

vector<tuple<vector<char>, vector<vector<size_t>>, vector<uint16_t>>> parse_to_raw(const string& filename) {
    vector<tuple<vector<char>, vector<vector<size_t>>, vector<uint16_t>>> inputs;

    ifstream file(filename);
    string line;
    while (getline(file, line)) {
        vector<char> light_pattern;
        vector<vector<size_t>> buttons;
        vector<uint16_t> joltages;

        size_t pos = 0;

        // Parse pattern in square brackets [.##.]
        size_t start = line.find('[', pos);
        size_t end = line.find(']', start);
        string pattern_str = line.substr(start + 1, end - start - 1);
        for (char c : pattern_str) {
            light_pattern.push_back(c);
        }
        pos = end + 1;

        // Parse groups in parentheses
        while ((start = line.find('(', pos)) != string::npos) {
            end = line.find(')', start);
            string group_str = line.substr(start + 1, end - start - 1);

            vector<size_t> group;
            stringstream ss(group_str);
            string num;
            while (getline(ss, num, ',')) {
                group.push_back(stoi(num));
            }
            buttons.push_back(group);
            pos = end + 1;

            // Check if we've reached the curly braces
            if (line.find('{', pos) < line.find('(', pos)) {
                break;
            }
        }

        // Parse final numbers in curly braces {3,5,4,7}
        start = line.find('{', pos);
        end = line.find('}', start);
        string final_str = line.substr(start + 1, end - start - 1);
        stringstream ss(final_str);
        string num;
        while (getline(ss, num, ',')) {
            joltages.push_back(stoi(num));
        }

        inputs.push_back(make_tuple(light_pattern, buttons, joltages));
    }

    return inputs;
}

void print_bitmasks(const vector<pair<mask16_t, vector<mask16_t>>>& bitmasks) {
    for (auto line : bitmasks) {
        cout << line.first << ", [ ";
        for (auto button : line.second) {
            cout << button << " ";
        }
        cout << "]" << endl;
    }
}

vector<pair<mask16_t, vector<mask16_t>>> to_bitmasks(
    const vector<tuple<vector<char>, vector<vector<size_t>>, vector<uint16_t>>>& raw_input
) {
    vector<pair<mask16_t, vector<mask16_t>>> bitmasks;

    // Convert raw input to bitmasks:
    // We have a max pattern length of 10 (for input 1)
    // [.##...#.] -> 0b0000000001100010 = 98
    // (0, 2, 3, 4) = [#.###] -> 0b0000000000010111
    mask16_t one = 1;

    for (const auto& line : raw_input) {
        const auto& pattern = get<0>(line);
        mask16_t lights = 0;
        for (char light : pattern) {
            lights <<= 1;       // bitshift by 1
            if (light == '#') {
                lights |= one;  // and mark bit if light is on
            }
        }

        vector<mask16_t> buttons;
        size_t offset = pattern.size() - 1;
        for (const auto& group : get<1>(line)) {
            mask16_t button = 0;
            for (size_t pos : group) {
                // Light 0 is (num_lights-1) bits from the right
                button |= (one << (offset - pos));
            }
            buttons.push_back(button);
        }

        // Ignore joltage for now
        bitmasks.push_back(make_pair(lights, buttons));
    }
    return bitmasks;
}

int bfs(mask16_t lights, const vector<mask16_t>& buttons) {
    queue<pair<int, mask16_t>> queue;    // <depth, current lights>
    queue.push(make_pair(0, 0));

    // Our lights are max length 10, so we can fit all permutations in a 1024-length array,
    // using the bitmask as an index directly
    char explored[1024];
    memset(explored, 0, 1024);

    while (!queue.empty()) {
        int depth = queue.front().first;
        mask16_t current_lights = queue.front().second;
        queue.pop();

        if (current_lights == lights) {
            // We've reached our goal permutation
            return depth;
        }

        if (explored[current_lights]) {
            // Skip light permutations we've already explored - they cannot lead to an optimal solution
            continue;
        }

        for (auto button : buttons) {
            mask16_t next = button ^ current_lights;
            queue.push(make_pair(depth + 1, next));
        }

        explored[current_lights] = true;
    }

    // We've exhausted all permutations and have not found our goal
    exit(1);
}

void part1() {
    const auto raw_input = parse_to_raw("input1.txt");
    auto start = chrono::high_resolution_clock::now();

    const auto bitmasks = to_bitmasks(raw_input);

    // Simple BFS until we find our end position, with some early
    // return conditions to prevent cycling
    int count = 0;
    for (const auto& line : bitmasks) {
        count += bfs(line.first, line.second);
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, std::milli> ms_double = end - start;
    cout << "Part 1: " << count << " in " << ms_double << endl;
}

int main() {
    part1();

    return 0;
}
