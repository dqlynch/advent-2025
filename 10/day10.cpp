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

const double EPS = 1e-9;
const double MAX = 1e32; 

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

using Tableau = vector<vector<double>>;

void print_tableau(const Tableau& tableau, const vector<int>& basis) {
    cout << "\n";
    for (size_t y = 0; y < tableau.size(); ++y) {
        const auto& row = tableau[y];
        cout << "[";
        for (auto x : row) {
            cout << ((x > 9 || x < 0) ? " " : "  ") << x;
        }
        cout << "  ]    ";
        if (y < basis.size()) {
            cout << ((basis[y] > 9 || basis[y] < 0) ? " " : "  ") << basis[y];
        }
        cout << endl;
    }
}

Tableau prepare_initial_tableau(
    tuple<vector<char>, vector<vector<size_t>>, vector<uint16_t>> raw_line
) {
    const auto& buttons = get<1>(raw_line);
    const auto& joltages = get<2>(raw_line);

    size_t h = joltages.size() + 1;                  // One per equation, plus one for objective func
    size_t w = buttons.size() + joltages.size() + 1; // One per variable, plus one per artificial variable, plus one for RHS

    auto tableau = Tableau(h, vector<double>(w));

    for (size_t i = 0; i < buttons.size(); ++i) {
        const auto& button = buttons.at(i);
        int column_total = 0;
        for (size_t idx : button) {
            tableau[idx][i] = 1;                // Equation constraint LHS; i.e. button press effect
            ++column_total;
        }
        tableau[joltages.size()][i] = column_total;
    }
    int column_total = 0;
    for (size_t i = 0; i < joltages.size(); ++i) {
        tableau[i][w - 1] = joltages.at(i);     // Equation constraint RHS; i.e. joltage amount
        column_total += joltages.at(i);
        tableau[i][i + buttons.size()] = 1;     // Slack variable identity matrix
    }
    tableau[h - 1][w - 1] = column_total;       // W function RHS = sum(RHS)


    return tableau;
}

double get_obj_rhs(const Tableau& tableau) {
    return tableau[tableau.size() - 1][tableau[0].size() - 1];
}

double get_max_obj_coef(const Tableau& tableau) {
        auto& obj_row = tableau.at(tableau.size() - 1);
        auto max_it = max_element(obj_row.begin(), obj_row.end() - 1);

        return *max_it;
}

void perform_pivot(Tableau* tableau, vector<int>* basis, size_t pivot_row_idx, size_t pivot_col_idx) {
    const size_t h = tableau->size();
    const size_t w = tableau->at(0).size();

    // Normalize our pivot row, such that the pivot value is 1
    auto& pivot_row = tableau->at(pivot_row_idx);
    double pivot_val = pivot_row[pivot_col_idx];
    for (size_t x = 0; x < w; ++x) {
        pivot_row[x] /= pivot_val;
    }

    // Subtract pivot row from each other row, normalizing such that the pivot column index end up as 0,
    // leaving an identity column for our chose pivot column - i.e. Gaussian elimination
    for (size_t y = 0; y < h; ++y) {
        // Skip the pivot row itself
        if (y == pivot_row_idx) {
            continue;
        }

        auto& row = tableau->at(y);
        double norm = row[pivot_col_idx];
        for (size_t x = 0; x < w; ++x) {
            row[x] -= pivot_row[x] * norm;
        }
    }

    // Update basis now that we've pivoted.
    basis->at(pivot_row_idx) = pivot_col_idx;
}

void pivot(Tableau* tableau, vector<int>* basis, bool contains_slack_vars) {
    const size_t h = tableau->size();
    const size_t w = tableau->at(0).size();
    const size_t num_constraints = h - 1;
    const size_t num_vars = contains_slack_vars ? w - 1 - num_constraints : w - 1;

    // 1. Identify pivot column - maximum objective value of real vars only
    const auto& obj_row = tableau->at(h - 1);
    size_t pivot_col_idx = max_element(obj_row.begin(), obj_row.begin() + num_vars) - obj_row.begin();

    // 2. Identify pivot row
    vector<double> ratios(num_constraints);
    for (size_t y = 0; y < num_constraints; ++y) {
        if (tableau->at(y)[pivot_col_idx] > EPS) {
            ratios[y] = tableau->at(y)[w - 1] / tableau->at(y)[pivot_col_idx];
        }
        else {
            ratios[y] = MAX;
        }
    }
    size_t pivot_row_idx = min_element(ratios.begin(), ratios.end()) - ratios.begin();

    perform_pivot(tableau, basis, pivot_row_idx, pivot_col_idx);
}

void flush_artificial_bases(Tableau* tableau, vector<int>* basis) {
    const size_t h = tableau->size();
    const size_t w = tableau->at(0).size();
    const size_t num_constraints = h - 1;
    const size_t num_vars = w - 1 - num_constraints;

    for (size_t y = 0; y < num_constraints; ++y) {
        if (basis->at(y) == -1) {
            // This row is still owned by an artifical variable. 
            // If there is a non-zero real variable in this row, we can pivot on it so it becomes the basis instead.
            for (size_t x = 0; x < num_vars; ++x) {
                if (abs(tableau->at(y)[x]) > EPS) {
                    perform_pivot(tableau, basis, y, x);
                    break;
                }
            }

            // If it's all zeros, then it's truly redundant - 0 + ... + 0 = 0
        }
    }
}

void remove_artificial_vars(Tableau* tableau) {
    const size_t h = tableau->size();
    const size_t w = tableau->at(0).size();
    const size_t num_constraints = h - 1;
    const size_t num_vars = w - 1 - num_constraints;

    for (auto& row : *tableau) {
        row.erase(row.begin() + num_vars, row.end() - 1);
    }
}

void calculate_objective(Tableau* tableau, vector<int>* basis) {
    const size_t h = tableau->size();
    const size_t w = tableau->at(0).size();

    // Start with a naive objective function - 1 coefficients everywhere...
    auto& obj_row = tableau->at(h - 1);
    for (size_t x = 0; x < w - 1; ++x) {
        obj_row[x] = -1;
    }

    // ...and then substitute out our basis variables, so that they remain basis variables.
    for (size_t y = 0; y < h - 1; ++y) {
        int basis_col_idx = basis->at(y);
        if (basis_col_idx == -1) {
            continue;
        }

        auto& basis_row = tableau->at(y);
        const double norm = obj_row[basis_col_idx];
        if (abs(norm) < EPS) {
            continue;
        }
        
        for (size_t x = 0; x < w; ++x) {
            obj_row[x] -= basis_row[x] * norm;
        }
    }
    
}


uint64_t solve(Tableau* tableau) {
    // The basis of rows 0 -> n-2; excludes W/Z
    // -1 basis is implied to be an artificial variable
    vector<int> basis(tableau->size() - 1, -1);     

    // Phase 1: Optimize our fake objective function to remove
    // slack variables
    cout << "Phase 1" << endl;
    while (get_obj_rhs(*tableau) > EPS) {
        pivot(tableau, &basis, true);
        print_tableau(*tableau, basis);
    }


    // First, flush out any remaninig artificial basis - if 
    // we skip these in the following step, then we may be ignoring
    // valid constraints.
    cout << "Phase 1.5" << endl;
    flush_artificial_bases(tableau, &basis);
    print_tableau(*tableau, basis);

    remove_artificial_vars(tableau);
    print_tableau(*tableau, basis);
    calculate_objective(tableau, &basis);
    print_tableau(*tableau, basis);

    cout << "Phase 2" << endl;
    // Phase 2: Optimize actual objective Z
    // Now we have _a_ solution - optimize it until we have no more positive
    // objective coefficients.
    while (get_max_obj_coef(*tableau) > EPS) {
        pivot(tableau, &basis, false);
        print_tableau(*tableau, basis);
    }

    //print_tableau(*tableau, basis);

    return llround(get_obj_rhs(*tableau));
}

void part2() {
    const auto raw_input = parse_to_raw("line12.txt");
    auto start = chrono::high_resolution_clock::now();

    /*
    1. Data format
        Consider that we have buttons [a, b, c, d]:
            [(0,1,2,3,4) (0,3,4) (0,1,2,4,5) (1,2)]
        and want to reach the target joltages:
            {10,11,11,5,10,5}

        We can phrase this as a system of linear equations:
            a + b + c      = 10
            a     + c + d  = 11
            a     + c + d  = 11
            a + b          = 5
            a + b + c      = 10
                    c      = 5

        with an objective function:
            Z - a - b - c - d - e = 0
        (We are trying to minimize Z = sum(a..e), i.e. minimize the 
        number of button presses to satisfy our system.)

        We'll can write this as a matrix,
            [[ 1,  1,  1,  0] | 10]
             [ 1,  0,  1,  1] | 11]
             [ 1,  0,  1,  1] | 11]
             [ 1,  1,  0,  0] | 5 ]
             [ 1,  1,  1,  0] | 10]
             [ 0,  0,  1,  0] | 5 ]
             [-1, -1, -1, -1] | 0 ]]
        where each row represents an equation from above, and each column represents a variable, or button.

        We could solve this with Gaussian Elimination and then optimize
        the resulting free variables, but this is tedious. We also have additional
        constraints: our variables must be non-negative, and we're optimizing an
        objective function. Instead, we can use the Simplex Algorithm
        - https://en.wikipedia.org/wiki/Simplex_algorithm

        Specifically, we'll use two phase Simplex:
        - http://www.lokminglui.com/lpch3.pdf (Two phase Simplex)

        Phase 1:
        We can't optimize with the starting point (a, ... , e) = (0, ..., 0), as this violates our constraints.
        We need a "basic feasible solution" (BFS) to start optimizing.
        We'll need to introduce artificial "slack" variables: one per equation, to get a starting (but fake) solution.
        We introduce (s_1, ..., s_6), with W = (s1 + ... + s6)
        We'll first optimize these to 0, replacing our objective function, which leaves us with a BFS using only real variables. 

        In Phase 2, when we've optimized to W = 0, we'll remove the artificial variables and substitute back in
        our real objective function. 

        So, adding our slack variables, our system of equations looks like:
            a + b + c     + s1 = 10
            a     + c + d + s2 = 11
            a     + c + d + s3 = 11
            a + b         + s4 = 5
            a + b + c     + s5 = 10
                    c     + s6 = 5
        with objective: minimize `W = s1 + s2 + s3 + s4 + s5 + s6`

        With some basic algebra (e.g. substitute `s1 = 10 - a - b - c`), we can write W in terms
        our RHS values and our actual variables:
                W = 52 - (5a + 3b + 5c + 2d)
            =>  W + 5a + 3b + 5c + 2d = 52
        Note that the coefficients are simply the sum of each column, so in practice our 
        initial W row values are simple to calculate.
            
               a   b   c   d  s1  s2  s3  s4  s5  s6   RHS    basis
            [[ 1,  1,  1,  0,  1,  0,  0,  0,  0,  0 | 10]     s1
             [ 1,  0,  1,  1,  0,  1,  0,  0,  0,  0 | 11]     s2
             [ 1,  0,  1,  1,  0,  0,  1,  0,  0,  0 | 11]     s3
             [ 1,  1,  0,  0,  0,  0,  0,  1,  0,  0 | 5 ]     s4
             [ 1,  1,  1,  0,  0,  0,  0,  0,  1,  0 | 10]     s5
             [ 0,  0,  1,  0,  0,  0,  0,  0,  0,  1 | 5 ]     s6
             [ 5,  3,  5,  2,  0,  0,  0,  0,  0,  0 | 52]]    W

            

    2. Pivot - Phase 1
        Now, we "pivot" our matrix - minimizing our artificial objective function W.
        * Pick the highest value of the bottom row. This is our pivot column. (ex: column `c` - abritrary between `a` and `c`)
            * This is the highest coefficient of our W, so the "most value" to reduce.
        * Pick the row with the smallest value of RHS / (pivot column), skipping 0 denominators. This is our pivot row. (ex: s6)
        * Pivot: normalize the pivot row and subtract it from all other rows such that pivot column is all 0
            * (This is standard Gaussian elimination)

            normalized pivot row:
             [ 0,  0,  1,  0,  0,  0,  0,  0,  0,  1 | 5 ]     s6
         
               a   b   c   d  s1  s2  s3  s4  s5  s6   RHS
            [[ 1,  1,  0,  0,  1,  0,  0,  0,  0, -1 | 5 ]     s1
             [ 1,  0,  0,  1,  0,  1,  0,  0,  0, -1 | 6 ]     s2
             [ 1,  0,  0,  1,  0,  0,  1,  0,  0, -1 | 6 ]     s3
             [ 1,  1,  0,  0,  0,  0,  0,  1,  0,  0 | 5 ]     s4
             [ 1,  1,  0,  0,  0,  0,  0,  0,  1, -1 | 5 ]     s5
             [ 0,  0,  1,  0,  0,  0,  0,  0,  0,  1 | 5 ]     c
             [ 5,  3,  0,  2,  0,  0,  0,  0,  0, -5 | 27]]    W

        Continue until W is 0...

            [[ 1,  1,  0,  0,  1,  0,  0,  0,  0, -1 | 5 ]     s1

               a   b   c   d  s1  s2  s3  s4  s5  s6   RHS
            [[ 1,  1,  0,  0,  1,  0,  0,  0,  0, -1 | 5 ]     a
             [ 0, -1,  0,  1, -1,  1,  0,  0,  0,  0 | 1 ]     s2
             [ 0, -1,  0,  1, -1,  0,  1,  0,  0,  0 | 1 ]     s3
             [ 0,  0,  0,  0, -1,  0,  0,  1,  0,  1 | 0 ]     s4
             [ 0,  0,  0,  0, -1,  0,  0,  0,  1,  0 | 0 ]     s5
             [ 0,  0,  1,  0,  0,  0,  0,  0,  0,  1 | 5 ]     c
             [ 0, -2,  0,  2, -5,  0,  0,  0,  0,  0 | 2]]     W


             [ 0, -1,  0,  1, -1,  1,  0,  0,  0,  0 | 1 ]     s2

               a   b   c   d  s1  s2  s3  s4  s5  s6   RHS
            [[ 1,  1,  0,  0,  1,  0,  0,  0,  0, -1 | 5 ]     a
             [ 0, -1,  0,  1, -1,  1,  0,  0,  0,  0 | 1 ]     d
             [ 0,  0,  0,  0,  0, -1,  1,  0,  0,  0 | 0 ]     s3
             [ 0,  0,  0,  0, -1,  0,  0,  1,  0,  1 | 0 ]     s4
             [ 0,  0,  0,  0, -1,  0,  0,  0,  1,  0 | 0 ]     s5
             [ 0,  0,  1,  0,  0,  0,  0,  0,  0,  1 | 5 ]     c
             [ 0,  0,  0,  0, -3, -2,  0,  0,  0,  0 | 0]]     W

        W is now 0. All of our slack variables are either out of the basis, or have a value of 0.

    3. Pivot - Phase 3
        Now, we can remove our artificial slack variables. We'll restore our original objective function:
          minimize Z = a + b + c + d
          (We can substitute basis variables (a, c, d) using Gaussian elimination)

               a   b   c   d    RHS
            [[ 1,  1,  0,  0, | 5 ]     a
             [ 0, -1,  0,  1, | 1 ]     d
             [ 0,  0,  0,  0, | 0 ]     null        <-- these rows represent redundant constraints
             [ 0,  0,  0,  0, | 0 ]     null        <--
             [ 0,  0,  0,  0, | 0 ]     null        <--
             [ 0,  0,  1,  0, | 5 ]     c
             [ 0, -1,  0,  0, | 11]]    Z

        We have no positive value in Z, so our tableau is already optimal.

    */

    uint64_t count = 0;
    for (const auto& line : raw_input) {
        auto tableau = prepare_initial_tableau(line);
        auto result = solve(&tableau);
        cout << result << endl;
        count += result;
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, std::milli> ms_double = end - start;
    //cout << "Part 2: " << count << " in " << ms_double << endl;
}

int main() {
    //part1();
    part2();

    return 0;
}