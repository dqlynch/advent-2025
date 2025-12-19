#include <vector>
#include <iostream>
#include <tuple>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#include <queue>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <cstring>

using namespace std;
using mask16_t = uint16_t;

const double EPS = 1e-9;
const double INT_TOL = 1e-5; 
const double MAX = 1e32; 

// Global debug flag
bool DEBUG_LOG = false;

// --- Parser ---

vector<tuple<vector<char>, vector<vector<size_t>>, vector<uint16_t>>> parse_input(istream& input) {
    vector<tuple<vector<char>, vector<vector<size_t>>, vector<uint16_t>>> inputs;
    string line;
    while (getline(input, line)) {
        if (line.empty()) continue;
        vector<char> light_pattern;
        vector<vector<size_t>> buttons;
        vector<uint16_t> joltages;
        size_t pos = 0;
        size_t start = line.find('[', pos);
        size_t end = line.find(']', start);
        if (start == string::npos || end == string::npos) continue;
        string pattern_str = line.substr(start + 1, end - start - 1);
        for (char c : pattern_str) light_pattern.push_back(c);
        pos = end + 1;
        while ((start = line.find('(', pos)) != string::npos) {
            end = line.find(')', start);
            string group_str = line.substr(start + 1, end - start - 1);
            vector<size_t> group;
            stringstream ss(group_str);
            string num;
            while (getline(ss, num, ',')) group.push_back(stoi(num));
            buttons.push_back(group);
            pos = end + 1;
            if (line.find('{', pos) < line.find('(', pos)) break;
        }
        start = line.find('{', pos);
        end = line.find('}', start);
        string final_str = line.substr(start + 1, end - start - 1);
        stringstream ss(final_str);
        string num;
        while (getline(ss, num, ',')) joltages.push_back(stoi(num));
        inputs.push_back(make_tuple(light_pattern, buttons, joltages));
    }
    return inputs;
}

// --- Simplex / ILP Implementation ---

using Tableau = vector<vector<double>>;
enum ConstraintType { LE, GE }; 
struct Constraint {
    int var_idx;
    ConstraintType type;
    int val;
};

struct SimplexResult {
    bool feasible;
    double objective_value;
    vector<double> variables;
};

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
    auto& pivot_row = tableau->at(pivot_row_idx);
    double pivot_val = pivot_row[pivot_col_idx];
    for (size_t x = 0; x < w; ++x) pivot_row[x] /= pivot_val;
    for (size_t y = 0; y < h; ++y) {
        if (y == pivot_row_idx) continue;
        auto& row = tableau->at(y);
        double norm = row[pivot_col_idx];
        if (abs(norm) > EPS) {
            for (size_t x = 0; x < w; ++x) row[x] -= pivot_row[x] * norm;
        }
    }
    basis->at(pivot_row_idx) = pivot_col_idx;
}

void pivot(Tableau* tableau, vector<int>* basis, size_t max_col_idx) {
    const size_t h = tableau->size();
    const size_t w = tableau->at(0).size();
    const size_t num_constraints = h - 1;
    const auto& obj_row = tableau->at(h - 1);
    
    // Maximization of Coefficient => Minimization of RHS (due to Gaussian elimination)
    int pivot_col_idx = -1;
    double max_val = EPS;
    for (size_t x = 0; x < max_col_idx; ++x) {
        if (obj_row[x] > max_val) {
            max_val = obj_row[x];
            pivot_col_idx = x;
        }
    }
    if (pivot_col_idx == -1) return; 

    int pivot_row_idx = -1;
    double min_ratio = MAX;
    for (size_t y = 0; y < num_constraints; ++y) {
        double val = tableau->at(y)[pivot_col_idx];
        if (val > EPS) {
            double rhs = tableau->at(y)[w - 1];
            if (abs(rhs) < EPS) rhs = 0.0; // Avoid negative zero issues
            double ratio = rhs / val;
            if (ratio >= -EPS && ratio < min_ratio) {
                min_ratio = ratio;
                pivot_row_idx = y;
            }
        }
    }
    if (pivot_row_idx != -1) perform_pivot(tableau, basis, pivot_row_idx, pivot_col_idx);
}

SimplexResult solve_simplex_with_constraints(
    const tuple<vector<char>, vector<vector<size_t>>, vector<uint16_t>>& raw_line,
    const vector<Constraint>& constraints,
    int depth
) {
    const auto& buttons = get<1>(raw_line);
    const auto& joltages = get<2>(raw_line);

    size_t num_vars = buttons.size();
    size_t num_base_eq = joltages.size();
    size_t num_extra = constraints.size();
    size_t total_constraints = num_base_eq + num_extra;
    
    // Column Mapping
    size_t col_idx = 0;
    size_t start_vars = col_idx; col_idx += num_vars;
    
    size_t start_le_slacks = col_idx;
    size_t count_le = 0; for(const auto& c : constraints) if(c.type == LE) count_le++;
    col_idx += count_le;

    size_t start_ge_surplus = col_idx;
    size_t count_ge = 0; for(const auto& c : constraints) if(c.type == GE) count_ge++;
    col_idx += count_ge;

    size_t start_artificials = col_idx;
    size_t num_artificials = num_base_eq + count_ge;
    // Track where artificials start to exclude them in Phase 2 pivot search
    size_t artificial_start_idx = col_idx;
    
    col_idx += num_artificials;

    size_t rhs_col = col_idx;
    size_t width = col_idx + 1;
    size_t height = total_constraints + 1; 

    Tableau tableau(height, vector<double>(width, 0.0));
    vector<int> basis(total_constraints, -1);

    // --- Fill Tableau ---
    size_t current_row = 0;
    size_t current_art = 0;
    size_t current_le_slack = 0;
    size_t current_ge_surplus = 0;

    // Base Equations
    for (size_t i = 0; i < num_base_eq; ++i) {
        for(size_t j=0; j<buttons.size(); ++j) {
            bool affects = false;
            for(size_t affected : buttons[j]) if(affected == i) affects = true;
            if(affects) tableau[current_row][start_vars + j] = 1.0;
        }
        tableau[current_row][start_artificials + current_art] = 1.0;
        basis[current_row] = start_artificials + current_art;
        
        tableau[current_row][rhs_col] = joltages[i];
        
        // Phase 1 Obj: Add row coeffs to Obj.
        for(size_t col=0; col<width-1; ++col) tableau[height-1][col] += tableau[current_row][col];
        tableau[height-1][start_artificials + current_art] = 0.0; // Basic var must be 0 in obj
        // Add RHS to Obj RHS
        tableau[height-1][rhs_col] += tableau[current_row][rhs_col];

        current_art++; current_row++;
    }

    // Constraints
    for (const auto& c : constraints) {
        tableau[current_row][start_vars + c.var_idx] = 1.0;
        tableau[current_row][rhs_col] = c.val;

        if (c.type == LE) {
            tableau[current_row][start_le_slacks + current_le_slack] = 1.0;
            basis[current_row] = start_le_slacks + current_le_slack;
            current_le_slack++;
        } else {
            tableau[current_row][start_ge_surplus + current_ge_surplus] = -1.0;
            tableau[current_row][start_artificials + current_art] = 1.0;
            basis[current_row] = start_artificials + current_art;
            
            // Phase 1 Obj
            for(size_t col=0; col<width-1; ++col) tableau[height-1][col] += tableau[current_row][col];
            tableau[height-1][start_artificials + current_art] = 0.0;
            tableau[height-1][rhs_col] += tableau[current_row][rhs_col];
            
            current_ge_surplus++; current_art++;
        }
        current_row++;
    }

    // Phase 1 Optimization
    int p1_iters = 0;
    while (get_max_obj_coef(tableau) > EPS) {
        pivot(&tableau, &basis, artificial_start_idx); 
        p1_iters++;
        if (p1_iters > 50000) break;
    }

    // --- CRITICAL FIX: DRIVE ARTIFICIALS OUT OF BASIS ---
    // If Phase 1 ended with 0 objective, but artificial variables remain in the basis (at 0),
    // we must swap them out or they will corrupt Phase 2 logic.
    for(size_t r = 0; r < total_constraints; ++r) {
        int col_idx = basis[r];
        // Check if basis var is Artificial
        if (col_idx >= (int)start_artificials && col_idx < (int)(start_artificials + num_artificials)) {
            // Check value. If non-zero, it's infeasible.
            double val = tableau[r][rhs_col];
            if (abs(val) > 1e-4) {
                 return { false, 0.0, {} };
            }
            
            // It is zero (degenerate). Try to pivot it out.
            // Look for a non-artificial column with non-zero coeff in this row.
            bool pivoted_out = false;
            for (size_t c = 0; c < artificial_start_idx; ++c) {
                if (abs(tableau[r][c]) > EPS) {
                    perform_pivot(&tableau, &basis, r, c);
                    pivoted_out = true;
                    break;
                }
            }
            
            // If we couldn't pivot it out, the row is redundant (0 = 0).
            // We can leave it, but it's safer to ensure it doesn't mess up objective calcs.
            // In Simplex, a redundant row with an ArtVar basis just stays 0=0. 
            // Phase 2 won't pick this row for pivot because min_ratio test handles 0s.
        }
    }
    // ----------------------------------------------------

    double p1_obj = get_obj_rhs(tableau);
    if (p1_obj > 1e-4) return { false, 0.0, {} };

    // Prepare Phase 2
    auto& obj_row = tableau[height-1];
    fill(obj_row.begin(), obj_row.end(), 0.0);
    for(size_t i=0; i<num_vars; ++i) obj_row[i] = -1.0; 

    // Substitute basis
    for(size_t r=0; r<total_constraints; ++r) {
        int basis_col = basis[r];
        if (basis_col == -1) continue; 
        double coeff = obj_row[basis_col]; // -1.0
        if (abs(coeff) > EPS) {
            for(size_t c=0; c<width; ++c) obj_row[c] -= tableau[r][c] * coeff;
        }
    }

    // Phase 2 Optimization
    int p2_iters = 0;
    while (get_max_obj_coef(tableau) > EPS) {
        pivot(&tableau, &basis, artificial_start_idx);
        p2_iters++;
        if (p2_iters > 50000) break;
    }

    double p2_obj = get_obj_rhs(tableau); 
    
    if (p2_obj < -1e-4) {
         return { false, 0.0, {} };
    }

    vector<double> vars(num_vars, 0.0);
    for(size_t r=0; r<total_constraints; ++r) {
        if(basis[r] < (int)num_vars && basis[r] >= 0) {
            vars[basis[r]] = tableau[r][rhs_col];
        }
    }

    // --- INTERNAL HARD VERIFICATION ---
    // Final check for floating point drift.
    for(size_t i = 0; i < num_base_eq; ++i) {
        double row_sum = 0;
        for(size_t j = 0; j < num_vars; ++j) {
            bool affects = false;
            for(size_t affected : buttons[j]) { 
                if(affected == i) { 
                    affects = true; 
                    break; 
                } 
            }
            if (affects) row_sum += vars[j];
        }
        if (abs(row_sum - joltages[i]) > 0.1) {
            return { false, 0.0, {} };
        }
    }
    // ----------------------------------
    
    return { true, p2_obj, vars };
}

// --- Branch and Bound ---

long long global_min_obj = -1;

void solve_bnb(
    const tuple<vector<char>, vector<vector<size_t>>, vector<uint16_t>>& raw_line,
    vector<Constraint> constraints,
    int depth
) {
    SimplexResult res = solve_simplex_with_constraints(raw_line, constraints, depth);

    if (!res.feasible) return;

    if (global_min_obj != -1 && res.objective_value >= global_min_obj - INT_TOL) return;

    int best_var = -1;
    double max_dist = -1.0; 
    double best_val = 0;

    for(size_t i=0; i<res.variables.size(); ++i) {
        double val = res.variables[i];
        double dist = abs(val - round(val));
        if (dist > INT_TOL) {
            if (dist > max_dist) {
                max_dist = dist;
                best_var = i;
                best_val = val;
            }
        }
    }

    if (best_var == -1) {
        long long int_obj = llround(res.objective_value);
        if (global_min_obj == -1 || int_obj < global_min_obj) {
            global_min_obj = int_obj;
        }
        return;
    }

    bool try_floor_first = (best_val - floor(best_val)) < 0.5;

    vector<Constraint> first_branch = constraints;
    vector<Constraint> second_branch = constraints;

    if (try_floor_first) {
        first_branch.push_back({best_var, LE, (int)floor(best_val)});
        second_branch.push_back({best_var, GE, (int)ceil(best_val)});
        solve_bnb(raw_line, first_branch, depth + 1);
        solve_bnb(raw_line, second_branch, depth + 1);
    } else {
        first_branch.push_back({best_var, GE, (int)ceil(best_val)});
        second_branch.push_back({best_var, LE, (int)floor(best_val)});
        solve_bnb(raw_line, first_branch, depth + 1);
        solve_bnb(raw_line, second_branch, depth + 1);
    }
}

void part2() {
    ifstream file("input1.txt");
    if (!file.is_open()) {
        cout << "Failed to open input1.txt" << endl;
        return;
    }
    auto inputs = parse_input(file);
    if (inputs.empty()) return;

    // cout << "Solving for " << inputs.size() << " lines in input1.txt..." << endl;
    auto start = chrono::high_resolution_clock::now();

    long long total_score = 0;
    
    for (size_t i = 0; i < inputs.size(); ++i) {
        global_min_obj = -1; // Reset for each line
        solve_bnb(inputs[i], {}, 0);

        if (global_min_obj != -1) {
            total_score += global_min_obj;
            cout << global_min_obj << endl;
        } else {
            cout << "Line " << i + 1 << ": No feasible integer solution found." << endl;
        }
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, std::milli> ms_double = end - start;
    
    cout << "Total Score: " << total_score << endl;
    // cout << "Total Time: " << ms_double.count() << "ms" << endl;
}

int main() {
    part2();
    return 0;
}