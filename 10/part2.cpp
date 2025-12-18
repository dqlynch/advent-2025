#include <vector>
#include <iostream>
#include <tuple>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <algorithm>

using namespace std;

const double EPS = 1e-9;
const double INT_TOL = 1e-5;
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
        // for (char c : pattern_str) {
        //     //light_pattern.push_back(c);  // We don't need this for part 2
        // }
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
    // leaving an identity column for our chosen pivot column - i.e. Gaussian elimination
    for (size_t y = 0; y < h; ++y) {
        // Skip the pivot row itself
        if (y == pivot_row_idx) {
            continue;
        }

        auto& row = tableau->at(y);
        double norm = row[pivot_col_idx];
        if (abs(norm) > EPS) {
            for (size_t x = 0; x < w; ++x) {
                row[x] -= pivot_row[x] * norm;
            }
        }
    }

    // Update basis now that we've pivoted.
    basis->at(pivot_row_idx) = pivot_col_idx;
}

void pivot(Tableau* tableau, vector<int>* basis, size_t max_col_idx) {
    const size_t h = tableau->size();
    const size_t w = tableau->at(0).size();
    const size_t num_constraints = h - 1;

    // 1. Identify pivot column - maximum objective value of real vars only
    const auto& obj_row = tableau->at(h - 1);
    int pivot_col_idx = -1;
    double max_val = EPS;
    for (size_t x = 0; x < max_col_idx; ++x) {
        if (obj_row[x] > max_val) {
            max_val = obj_row[x];
            pivot_col_idx = x;
        }
    }

    if (pivot_col_idx == -1) {
        return;
    }

    // 2. Identify pivot row - minimum ratio test
    int pivot_row_idx = -1;
    double min_ratio = MAX;
    for (size_t y = 0; y < num_constraints; ++y) {
        double val = tableau->at(y)[pivot_col_idx];
        if (val > EPS) {
            double rhs = tableau->at(y)[w - 1];
            if (abs(rhs) < EPS) {
                rhs = 0.0; // Avoid negative zero issues
            }
            double ratio = rhs / val;
            if (ratio >= -EPS && ratio < min_ratio) {
                min_ratio = ratio;
                pivot_row_idx = y;
            }
        }
    }

    if (pivot_row_idx != -1) {
        perform_pivot(tableau, basis, pivot_row_idx, pivot_col_idx);
    }
}

bool check_and_flush_artificial_bases(
    Tableau* tableau,
    vector<int>* basis,
    size_t total_constraints,
    size_t num_vars,
    size_t start_artificials,
    size_t num_artificials,
    size_t rhs_col
) {
    // If Phase 1 ended with 0 objective, but artificial variables remain in the basis (at 0),
    // we must swap them out or they will corrupt Phase 2 logic.
    for (size_t r = 0; r < total_constraints; ++r) {
        int col_idx = basis->at(r);

        // Check if basis var is artificial
        if (col_idx >= (int)start_artificials && col_idx < (int)(start_artificials + num_artificials)) {
            // Check RHS value - if non-zero, it's infeasible
            double val = tableau->at(r)[rhs_col];
            if (abs(val) > 1e-4) {
                return false; // Infeasible
            }

            // It is zero (degenerate). Try to pivot it out.
            // Look for a non-artificial column with non-zero coeff in this row.
            for (size_t c = 0; c < num_vars; ++c) {
                if (abs(tableau->at(r)[c]) > EPS) {
                    perform_pivot(tableau, basis, r, c);
                    break;
                }
            }

            // If we couldn't pivot it out, the row is redundant (0 = 0): we don't have to 
            // do anything, and this row will be ignored. 
        }
    }

    return true; // Feasible
}

struct TableauDimensions {
    size_t num_vars;
    size_t num_base_eq;
    size_t num_extra;
    size_t total_constraints;
    size_t count_le;
    size_t count_ge;
    size_t num_artificials;

    size_t start_vars;
    size_t start_le_slacks;
    size_t start_ge_surplus;
    size_t start_artificials;
    size_t rhs_col;

    size_t width;
    size_t height;
};

TableauDimensions calculate_tableau_dimensions(
    size_t num_vars,
    size_t num_base_eq,
    const vector<Constraint>& constraints
) {
    TableauDimensions dim;

    dim.num_vars = num_vars;
    dim.num_base_eq = num_base_eq;
    dim.num_extra = constraints.size();
    dim.total_constraints = num_base_eq + dim.num_extra;

    // Count constraint types
    dim.count_le = 0;
    dim.count_ge = 0;
    for (const auto& c : constraints) {
        if (c.type == LE) {
            dim.count_le++;
        } else {
            dim.count_ge++;
        }
    }

    // Column mapping: [variables | LE slack vars | GE surplus vars | artificial vars | RHS]
    size_t col_idx = 0;
    dim.start_vars = col_idx;
    col_idx += num_vars;

    dim.start_le_slacks = col_idx;
    col_idx += dim.count_le;

    dim.start_ge_surplus = col_idx;
    col_idx += dim.count_ge;

    dim.start_artificials = col_idx;
    dim.num_artificials = num_base_eq + dim.count_ge;
    col_idx += dim.num_artificials;

    dim.rhs_col = col_idx;
    dim.width = col_idx + 1;
    dim.height = dim.total_constraints + 1;

    return dim;
}

void initialize_tableau(
    Tableau* tableau,
    vector<int>* basis,
    const TableauDimensions& dim,
    const vector<vector<size_t>>& buttons,
    const vector<uint16_t>& joltages,
    const vector<Constraint>& constraints
) {
    size_t current_row = 0;
    size_t current_le_slack = 0;
    size_t current_ge_surplus = 0;

    // Step 1: Add base equality constraints (structure only, no artificial vars yet)
    for (size_t i = 0; i < dim.num_base_eq; ++i) {
        for (size_t j = 0; j < buttons.size(); ++j) {
            bool affects = false;
            for (size_t affected : buttons[j]) {
                if (affected == i) {
                    affects = true;
                    break;
                }
            }
            if (affects) {
                (*tableau)[current_row][dim.start_vars + j] = 1.0;
            }
        }
        (*tableau)[current_row][dim.rhs_col] = joltages[i];
        current_row++;
    }

    // Step 2: Add inequality constraints from Branch and Bound
    for (const auto& c : constraints) {
        (*tableau)[current_row][dim.start_vars + c.var_idx] = 1.0;
        (*tableau)[current_row][dim.rhs_col] = c.val;

        if (c.type == LE) {
            // x <= k  =>  x + s = k  (slack variable)
            (*tableau)[current_row][dim.start_le_slacks + current_le_slack] = 1.0;
            (*basis)[current_row] = dim.start_le_slacks + current_le_slack;
            current_le_slack++;
        } else {
            // x >= k  =>  x - s = k  (surplus variable, artificial added later)
            (*tableau)[current_row][dim.start_ge_surplus + current_ge_surplus] = -1.0;
            current_ge_surplus++;
        }
        current_row++;
    }

    // Step 3: Add artificial variables and build Phase 1 objective
    // Artificial variables are needed for:
    // - All base equality constraints (no slack/surplus to use as initial basis)
    // - All GE constraints (surplus variable is negative, can't be in initial basis)
    size_t current_art = 0;

    // Add artificial variables for base equality constraints
    for (size_t i = 0; i < dim.num_base_eq; ++i) {
        (*tableau)[i][dim.start_artificials + current_art] = 1.0;
        (*basis)[i] = dim.start_artificials + current_art;

        // Add this row's contribution to Phase 1 objective
        for (size_t col = 0; col < dim.width - 1; ++col) {
            (*tableau)[dim.height - 1][col] += (*tableau)[i][col];
        }
        (*tableau)[dim.height - 1][dim.start_artificials + current_art] = 0.0; // Basic var must be 0 in obj
        (*tableau)[dim.height - 1][dim.rhs_col] += (*tableau)[i][dim.rhs_col];

        current_art++;
    }

    // Add artificial variables for GE constraints
    for (size_t i = 0; i < constraints.size(); ++i) {
        if (constraints[i].type == GE) {
            size_t row = dim.num_base_eq + i;
            (*tableau)[row][dim.start_artificials + current_art] = 1.0;
            (*basis)[row] = dim.start_artificials + current_art;

            // Add this row's contribution to Phase 1 objective
            for (size_t col = 0; col < dim.width - 1; ++col) {
                (*tableau)[dim.height - 1][col] += (*tableau)[row][col];
            }
            (*tableau)[dim.height - 1][dim.start_artificials + current_art] = 0.0;
            (*tableau)[dim.height - 1][dim.rhs_col] += (*tableau)[row][dim.rhs_col];

            current_art++;
        }
    }
}

void calculate_objective(Tableau* tableau, vector<int>* basis, size_t num_vars) {
    const size_t h = tableau->size();
    const size_t w = tableau->at(0).size();

    // Start with a naive objective function - 1 coefficients everywhere...
    auto& obj_row = tableau->at(h - 1);
    fill(obj_row.begin(), obj_row.end(), 0.0);
    for (size_t x = 0; x < num_vars; ++x) {
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

SimplexResult solve_simplex_with_constraints(
    const tuple<vector<char>, vector<vector<size_t>>, vector<uint16_t>>& raw_line,
    const vector<Constraint>& constraints
) {
    /*

        Phase 1:
        We can't optimize with the starting point (a, ... , e) = (0, ..., 0), as this violates our constraints.
        We need a "basic feasible solution" (BFS) to start optimizing.
        We'll need to introduce artificial variables: one per equation, to get a starting (but fake) solution.
        We introduce (a_1, ..., a_6), with W = (a1 + ... + a6)
        We'll first optimize these to 0, replacing our objective function, which leaves us with a BFS using only real variables.

        In Phase 2, when we've optimized to W = 0, we'll remove the artificial variables and substitute back in
        our real objective function.

        So, adding our artificial variables, our system of equations looks like:
            a + b + c     + a1 = 10
            a     + c + d + a2 = 11
            a     + c + d + a3 = 11
            a + b         + a4 = 5
            a + b + c     + a5 = 10
                    c     + a6 = 5
        with objective: minimize `W = a1 + a2 + a3 + a4 + a5 + a6`

        With some basic algebra (e.g. substitute `a1 = 10 - a - b - c`), we can write W in terms
        our RHS values and our actual variables:
                W = 52 - (5a + 3b + 5c + 2d)
            =>  W + 5a + 3b + 5c + 2d = 52
        Note that the coefficients are simply the sum of each column, so in practice our
        initial W row values are simple to calculate.

               a   b   c   d  s1  s2  s3  s4  s5  s6   RHS    basis
            [[ 1,  1,  1,  0,  1,  0,  0,  0,  0,  0 | 10]     a1
             [ 1,  0,  1,  1,  0,  1,  0,  0,  0,  0 | 11]     a2
             [ 1,  0,  1,  1,  0,  0,  1,  0,  0,  0 | 11]     a3
             [ 1,  1,  0,  0,  0,  0,  0,  1,  0,  0 | 5 ]     a4
             [ 1,  1,  1,  0,  0,  0,  0,  0,  1,  0 | 10]     a5
             [ 0,  0,  1,  0,  0,  0,  0,  0,  0,  1 | 5 ]     a6
             [ 5,  3,  5,  2,  0,  0,  0,  0,  0,  0 | 52]]    W



    2. Pivot - Phase 1
        Now, we "pivot" our matrix - minimizing our artificial objective function W.
        * Pick the highest value of the bottom row. This is our pivot column. (ex: column `c` - arbitrary between `a` and `c`)
            * This is the highest coefficient of our W, so the "most value" to reduce.
        * Pick the row with the smallest value of RHS / (pivot column), skipping 0 denominators. This is our pivot row. (ex: s6)
        * Pivot: normalize the pivot row and subtract it from all other rows such that pivot column is all 0
            * (This is standard Gaussian elimination)

            normalized pivot row:
             [ 0,  0,  1,  0,  0,  0,  0,  0,  0,  1 | 5 ]     a6

               a   b   c   d  a1  a2  a3  a4  a5  a6   RHS
            [[ 1,  1,  0,  0,  1,  0,  0,  0,  0, -1 | 5 ]     a1
             [ 1,  0,  0,  1,  0,  1,  0,  0,  0, -1 | 6 ]     a2
             [ 1,  0,  0,  1,  0,  0,  1,  0,  0, -1 | 6 ]     a3
             [ 1,  1,  0,  0,  0,  0,  0,  1,  0,  0 | 5 ]     a4
             [ 1,  1,  0,  0,  0,  0,  0,  0,  1, -1 | 5 ]     a5
             [ 0,  0,  1,  0,  0,  0,  0,  0,  0,  1 | 5 ]     c
             [ 5,  3,  0,  2,  0,  0,  0,  0,  0, -5 | 27]]    W

        Continue until W is 0...

            [[ 1,  1,  0,  0,  1,  0,  0,  0,  0, -1 | 5 ]     a1

               a   b   c   d  a1  a2  a3  a4  a5  a6   RHS
            [[ 1,  1,  0,  0,  1,  0,  0,  0,  0, -1 | 5 ]     a
             [ 0, -1,  0,  1, -1,  1,  0,  0,  0,  0 | 1 ]     a2
             [ 0, -1,  0,  1, -1,  0,  1,  0,  0,  0 | 1 ]     a3
             [ 0,  0,  0,  0, -1,  0,  0,  1,  0,  1 | 0 ]     a4
             [ 0,  0,  0,  0, -1,  0,  0,  0,  1,  0 | 0 ]     a5
             [ 0,  0,  1,  0,  0,  0,  0,  0,  0,  1 | 5 ]     c
             [ 0, -2,  0,  2, -5,  0,  0,  0,  0,  0 | 2]]     W


             [ 0, -1,  0,  1, -1,  1,  0,  0,  0,  0 | 1 ]     a2

               a   b   c   d  a1  a2  a3  a4  a5  a6   RHS
            [[ 1,  1,  0,  0,  1,  0,  0,  0,  0, -1 | 5 ]     a
             [ 0, -1,  0,  1, -1,  1,  0,  0,  0,  0 | 1 ]     d
             [ 0,  0,  0,  0,  0, -1,  1,  0,  0,  0 | 0 ]     a3
             [ 0,  0,  0,  0, -1,  0,  0,  1,  0,  1 | 0 ]     a4
             [ 0,  0,  0,  0, -1,  0,  0,  0,  1,  0 | 0 ]     a5
             [ 0,  0,  1,  0,  0,  0,  0,  0,  0,  1 | 5 ]     c
             [ 0,  0,  0,  0, -3, -2,  0,  0,  0,  0 | 0]]     W

        W is now 0. All of our artificial variables are either out of the basis, or have a value of 0.

    3. Pivot - Phase 2
        Now, we can remove our artificial variables. We'll restore our original objective function:
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

        We have no positive value in Z, so our tableau is already optimal. Typically, we would pivot
        with the same logic as above to optimize our objective function. We're attempting to minimize Z, so we
        want to pivot until there are no remaining positive coefficients in our objective row. 

    N.B. - Additional Constraints
        We'll need to also handle integer constraints. Simplex optimizes in continous space, so we can get results
        with fractional button presses. Fortunately, we already optimize for positive values.

        We handle integer requirements by allowing the calling function to add arbitrary additional <= or >= constraints.
        For example, if we return a solution with an (x = 3.7) value for a variable, the calling function will force it to the
        nearest integer with one of two additional constraints: (x <= 3), or (x >= 4). 

        We model inequalities using slack or surplus variables. 
        - For x <= k: Add slack variable s, where x + s = k
        - For x >= k: Add surplus variable s and artificial variable a, where x - s + a = k
          (The artificial variable is needed to get an initial basic feasible solution)
    */

    const auto& buttons = get<1>(raw_line);
    const auto& joltages = get<2>(raw_line);

    // Calculate tableau dimensions and setup
    TableauDimensions dim = calculate_tableau_dimensions(buttons.size(), joltages.size(), constraints);

    Tableau tableau(dim.height, vector<double>(dim.width, 0.0));
    vector<int> basis(dim.total_constraints, -1);

    // Set up an initial tableau with:
    // - Our equality constraints from the initial problem
    // - Any additional LE/GE constraints for integer forcing
    // - Our artificial variables, necessary to reach a BFS
    initialize_tableau(&tableau, &basis, dim, buttons, joltages, constraints);

    // Phase 1: Optimize to remove artificial variables
    int phase1_iters = 0;
    while (get_max_obj_coef(tableau) > EPS) {
        pivot(&tableau, &basis, dim.start_artificials);
        phase1_iters++;
        if (phase1_iters > 50000) {
            break;
        }
    }

    // Check and flush artificial bases
    if (!check_and_flush_artificial_bases(&tableau, &basis, dim.total_constraints, dim.num_vars,
                                          dim.start_artificials, dim.num_artificials, dim.rhs_col)) {
        return { false, 0.0, {} };
    }

    // Check feasibility: Phase 1 objective should be 0
    double phase1_obj = get_obj_rhs(tableau);
    if (phase1_obj > 1e-4) {
        return { false, 0.0, {} };
    }

    // Phase 2: Optimize real objective function
    calculate_objective(&tableau, &basis, dim.num_vars);

    int phase2_iters = 0;
    while (get_max_obj_coef(tableau) > EPS) {
        pivot(&tableau, &basis, dim.start_artificials);
        phase2_iters++;
        if (phase2_iters > 50000) {
            break;
        }
    }

    double phase2_obj = get_obj_rhs(tableau);

    if (phase2_obj < -1e-4) {
        return { false, 0.0, {} };
    }

    // Extract variable values from basis
    vector<double> vars(dim.num_vars, 0.0);
    for (size_t r = 0; r < dim.total_constraints; ++r) {
        if (basis[r] < (int)dim.num_vars && basis[r] >= 0) {
            vars[basis[r]] = tableau[r][dim.rhs_col];
        }
    }

    // Verification: check that our solution satisfies base constraints
    for (size_t i = 0; i < dim.num_base_eq; ++i) {
        double row_sum = 0;
        for (size_t j = 0; j < dim.num_vars; ++j) {
            bool affects = false;
            for (size_t affected : buttons[j]) {
                if (affected == i) {
                    affects = true;
                    break;
                }
            }
            if (affects) {
                row_sum += vars[j];
            }
        }
        if (abs(row_sum - joltages[i]) > 0.1) {
            return { false, 0.0, {} };
        }
    }

    return { true, phase2_obj, vars };
}

void solve_bnb(
    const tuple<vector<char>, vector<vector<size_t>>, vector<uint16_t>>& raw_line,
    vector<Constraint> constraints,
    int64_t* best_obj,
    int depth
) {
    /*
    1. Modeling Our Problem
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

    2. Solving the System of Equations

         We can use the Simplex Algorithm to solve this system of equations and optimize for our objective function.
        - https://en.wikipedia.org/wiki/Simplex_algorithm

        Phrasing our system of equations as a matrix...
            [[ 1,  1,  1,  0] | 10]
             [ 1,  0,  1,  1] | 11]
             [ 1,  0,  1,  1] | 11]
             [ 1,  1,  0,  0] | 5 ]
             [ 1,  1,  1,  0] | 10]
             [ 0,  0,  1,  0] | 5 ]
             [-1, -1, -1, -1] | 0 ]]
        ...where each row represents an equation from above, and each column represents a variable, or button.

        Specifically, we'll use two phase Simplex:
        - http://www.lokminglui.com/lpch3.pdf (Two phase Simplex)

        We'll first find a basic feasible solution (BFS), and then optimize it according to our objective function Z.
        For details, see `solve_simplex_with_constraints`. (This is where the vast majority of effort is!)

    4. Dealing with Integer Constraints
        The simplex algorithm gives us continuous solutions, but we need integer ones - we can't press a button 3.5 times.

        Branch and Bound works by:
        1. Solving the relaxed Linear Programming problem (allowing fractional values)
        2. If all variables are integers, we're done
        3. If not, pick a fractional variable (e.g., x = 3.7) and create two branches to force
           the fractional variable to its neighboring integer values:
           - Branch 1: Add constraint x <= 3
           - Branch 2: Add constraint x >= 4
        4. Recursively solve each branch
        5. Prune branches where:
           - The LP is infeasible (no solution exists) - This is possible because we're 
             introducing arbitrary constraints!
           - The LP objective is worse than the best integer solution found so far

        For additional constraints from branching:
        - For x <= k: Add slack variable s, where x + s = k
        - For x >= k: Add surplus variable s and artificial variable a, where x - s + a = k
          (The artificial variable is needed to get an initial basic feasible solution)


    */
    SimplexResult fractional_result = solve_simplex_with_constraints(raw_line, constraints);

    // Prune if infeasible...
    if (!fractional_result.feasible) {
        return;
    }

    // ...or worse than current best.
    if (*best_obj != -1 && fractional_result.objective_value >= *best_obj - INT_TOL) {
        return;
    }

    // Find the most fractional variable
    int best_var = -1;
    double max_dist = -1.0;
    double best_val = 0;

    for (size_t i = 0; i < fractional_result.variables.size(); ++i) {
        double val = fractional_result.variables[i];
        double dist = abs(val - round(val));
        if (dist > INT_TOL) {
            if (dist > max_dist) {
                max_dist = dist;
                best_var = i;
                best_val = val;
            }
        }
    }

    // All variables are integers - we have a solution!
    if (best_var == -1) {
        int16_t int_obj = llround(fractional_result.objective_value);
        if (*best_obj == -1 || int_obj < *best_obj) {
            *best_obj = int_obj;
        }
        return;
    }

    // Branch: create two subproblems
    // Try the closer branch first (heuristic for faster convergence)
    bool try_floor_first = (best_val - floor(best_val)) < 0.5;

    vector<Constraint> first_branch = constraints;
    vector<Constraint> second_branch = constraints;

    if (try_floor_first) {
        first_branch.push_back({best_var, LE, (int)floor(best_val)});
        second_branch.push_back({best_var, GE, (int)ceil(best_val)});
        solve_bnb(raw_line, first_branch, best_obj, depth + 1);
        solve_bnb(raw_line, second_branch, best_obj, depth + 1);
    } else {
        first_branch.push_back({best_var, GE, (int)ceil(best_val)});
        second_branch.push_back({best_var, LE, (int)floor(best_val)});
        solve_bnb(raw_line, first_branch, best_obj, depth + 1);
        solve_bnb(raw_line, second_branch, best_obj, depth + 1);
    }
}

void part2() {
    const auto raw_input = parse_to_raw("input1.txt");
    auto start = chrono::high_resolution_clock::now();

    uint64_t count = 0;
    for (const auto& line : raw_input) {
        int64_t best_obj = -1; // Start with no solution found
        solve_bnb(line, {}, &best_obj, 0);

        if (best_obj != -1) {
            count += best_obj;
            //cout << best_obj << endl;
        } else {
            cout << "No feasible integer solution found." << endl;
        }
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, std::milli> ms_double = end - start;
    cout << "Part 2: " << count << " in " << ms_double << endl;
}

int main() {
    part2();

    return 0;
}
