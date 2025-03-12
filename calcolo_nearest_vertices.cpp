#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <random>
#include <list>

using namespace std;

constexpr int NUM_VERTS_OVERS = 38920; // numero vertici triangolazione aumentata

constexpr int NUM_VERTS = 26850; // numero vertici triangolazione originale
constexpr int NUM_TRIANGLES = 53696; // numero triangoli triangolazione originale

constexpr double cutoff_range_pert = 1; // range massimo per calcolo vertici vicini PER PERTURBAZIONE
constexpr double cutoff_range_smooth = 6; // range massimo per calcolo vertici vicini PER SMOOTHING GEODETICO

constexpr int NUM_DIMENSIONS = 3;

void cross(const double a[NUM_DIMENSIONS], const double b[NUM_DIMENSIONS], double * res) {
    res[0] = a[1] * b[2] - a[2] * b[1];
    res[1] = a[2] * b[0] - a[0] * b[2];
    res[2] = a[0] * b[1] - a[1] * b[0];
}

int read_dim_1(const string& file_name) {
    ifstream file(file_name);

    string line;
    getline(file, line);
    stringstream ss(line);
    string cell;
    getline(ss, cell, ',');

    int dim1 = stoi(cell);
    return dim1;
}

void read_matrix_from_csv(const string& file_name, int dim1, int dim2, double** mat, bool read_dim_1 = false, int skip = 0) {
    ifstream file(file_name);

    string line;

    if (read_dim_1) {
        getline(file, line);
        stringstream ss(line);

        string cell;
        getline(ss, cell, ',');

        dim1 = stoi(cell);
    }

    int row = 0;
    while (getline(file, line) && row < dim1) {
        stringstream ss(line);
        string cell;
        int col = 0;
        while (getline(ss, cell, ',') && col < dim2) {
            if(skip <= row) {
                mat[row - skip][col] = stod(cell);
            }
            col++;
        }
        row++;
    }
    file.close();
}

void read_matrix_from_csv(const string& file_name, int dim1, int dim2, int** mat, bool read_dim_1 = false) {
    ifstream file(file_name);

    string line;

    if (read_dim_1) {
        getline(file, line);
        stringstream ss(line);

        string cell;
        getline(ss, cell, ',');

        dim1 = stoi(cell);
    }

    int row = 0;
    while (getline(file, line) && row < dim1) {
        stringstream ss(line);
        string cell;
        int col = 0;
        while (getline(ss, cell, ',') && col < dim2) {
            mat[row][col] = stoi(cell);
            col++;
        }
        row++;
    }
    file.close();
}

double closest_point_on_triangle_dist(const double a[NUM_DIMENSIONS], const double b[NUM_DIMENSIONS],
    const double c[NUM_DIMENSIONS], const double p[NUM_DIMENSIONS], double* proj_point){
    double ab[NUM_DIMENSIONS];
    ab[0] = b[0] - a[0];
    ab[1] = b[1] - a[1];
    ab[2] = b[2] - a[2];
    double ac[NUM_DIMENSIONS];
    ac[0] = c[0] - a[0];
    ac[1] = c[1] - a[1];
    ac[2] = c[2] - a[2];
    double ap[NUM_DIMENSIONS];
    ap[0] = p[0] - a[0];
    ap[1] = p[1] - a[1];
    ap[2] = p[2] - a[2];

    const double d1 = ab[0] * ap[0] + ab[1] * ap[1] + ab[2] * ap[2]; // dot(ab, ap)
    const double d2 = ac[0] * ap[0] + ac[1] * ap[1] + ac[2] * ap[2]; // dot(ac, ap)
    if (d1 <= 0.0 && d2 <= 0.0){
        proj_point[0] = a[0];
        proj_point[1] = a[1];
        proj_point[2] = a[2];
    } else {
        double bp[NUM_DIMENSIONS];
        bp[0] = p[0] - b[0];
        bp[1] = p[1] - b[1];
        bp[2] = p[2] - b[2];

        const double d3 = ab[0] * bp[0] + ab[1] * bp[1] + ab[2] * bp[2]; // dot(ab, bp);
        const double d4 = ac[0] * bp[0] + ac[1] * bp[1] + ac[2] * bp[2]; // dot(ac, bp);
        if (d3 >= 0.0 && d4 <= d3){
            proj_point[0] = b[0];
            proj_point[1] = b[1];
            proj_point[2] = b[2];
        } else {
            double cp[NUM_DIMENSIONS];
            cp[0] = p[0] - c[0];
            cp[1] = p[1] - c[1];
            cp[2] = p[2] - c[2];
            const double d5 = ab[0] * cp[0] + ab[1] * cp[1] + ab[2] * cp[2]; // dot(ab, cp)
            const double d6 = ac[0] * cp[0] + ac[1] * cp[1] + ac[2] * cp[2]; // dot(ac, cp)
            if (d6 >= 0.0 && d5 <= d6){
                proj_point[0] = c[0];
                proj_point[1] = c[1];
                proj_point[2] = c[2];
            } else {
                const double vc = d1 * d4 - d3 * d2;
                if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0){
                    const double v = d1 / (d1 - d3);

                    proj_point[0] = a[0] + v * ab[0];
                    proj_point[1] = a[1] + v * ab[1];
                    proj_point[2] = a[2] + v * ab[2];
                } else {
                    const double vb = d5 * d2 - d1 * d6;
                    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0){
                        const double v = d2 / (d2 - d6);

                        proj_point[0] = a[0] + v * ac[0];
                        proj_point[1] = a[1] + v * ac[1];
                        proj_point[2] = a[2] + v * ac[2];
                    } else {
                        const double va = d3 * d6 - d5 * d4;
                        if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0){
                            const double v = (d4 - d3) / ((d4 - d3) + (d5 - d6));

                            double bc[NUM_DIMENSIONS];
                            bc[0] = c[0] - b[0];
                            bc[1] = c[1] - b[1];
                            bc[2] = c[2] - b[2];

                            proj_point[0] = b[0] + v * bc[0];
                            proj_point[1] = b[1] + v * bc[1];
                            proj_point[2] = b[2] + v * bc[2];
                        } else {
                            const double den = 1.0 / (va + vb + vc);
                            const double v = vb * den;
                            const double w = vc * den;
                            proj_point[0] = a[0] + v * ab[0] + w * ac[0];
                            proj_point[1] = a[1] + v * ab[1] + w * ac[1];
                            proj_point[2] = a[2] + v * ab[2] + w * ac[2];
                        }
                    }
                }
            }
        }
    }

    double d_tmp = (p[0] - proj_point[0])*(p[0] - proj_point[0]) +
            (p[1] - proj_point[1])*(p[1] - proj_point[1]) +
                (p[2] - proj_point[2])*(p[2] - proj_point[2]);
    return d_tmp;
}

class kd_tree {
    public:
        kd_tree(double** verts, int** triangles, const int depth, const int* current_indices, const int size, const int current_depth = 0) {
            if(depth - current_depth > 0) {
                this->size = size;

                if(this->size > 1) {
                    auto* col_tmp = new double[this->size];
                    auto* sel_tmp_right = new bool[this->size];
                    auto* sel_tmp_left = new bool[this->size];

                    int num_right = 0;
                    int num_left = 0;

                    int c_tmp = 0;
                    this->split_dimension = (current_depth % NUM_DIMENSIONS) - 1;

                    do{
                        c_tmp++;
                        this->split_dimension = (this->split_dimension + 1) % NUM_DIMENSIONS;

                        for (int i = 0; i < this->size; i++) {
                            double s_tmp = 0;
                            for(int j = 0; j < 3; j++) {
                                s_tmp += verts[triangles[current_indices[i]][j]][this->split_dimension];
                            }
                            col_tmp[i] = s_tmp/3; // ATTENZIONE: sarebbe forse meglio usare il massimo per ogni triangolo (forse non "divido" solo un triangolo in meno?)...

                            /*
                            for (int j = 0; j < 3; j++) {
                                col_tmp[i * 3 + j] = verts[triangles[current_indices[i]][j]][this->split_dimension];
                            }
                            */
                        }
                        sort(col_tmp, col_tmp + this->size);

                        if (this->size % 2 == 0) {
                            this->split_val = (col_tmp[this->size / 2 - 1] + col_tmp[this->size / 2]) / 2;
                        }
                        else {
                            this->split_val = col_tmp[(this->size-1) / 2];
                        }

                        this->current_depth = current_depth;
                        this->remaining_depth = depth - current_depth;

                        num_right = 0;
                        num_left = 0;
                        for (int i = 0; i < this->size; i++) {
                            // sel_tmp[i] = verts[current_indices[i]][this->split_dimension] > this->split_val;

                            double tmp_max = verts[triangles[current_indices[i]][0]][this->split_dimension];
                            double tmp_min = verts[triangles[current_indices[i]][0]][this->split_dimension];
                            for(int j = 1; j < 3; j++) {
                                double tmp = verts[triangles[current_indices[i]][j]][this->split_dimension];
                                if(tmp_max < tmp) {
                                    tmp_max = tmp;
                                }
                                if(tmp_min > tmp) {
                                    tmp_min = tmp;
                                }
                            }

                            sel_tmp_right[i] = tmp_max > this->split_val;
                            sel_tmp_left[i] = tmp_min <= this->split_val;

                            num_right += sel_tmp_right[i];
                            num_left += sel_tmp_left[i];
                        }
                    }while( (num_right == this->size || num_left == this->size) && c_tmp < 3);
                    delete[] col_tmp;

                    if(num_right == this->size || num_left == this->size) {
                        delete[] sel_tmp_left;
                        delete[] sel_tmp_right;

                        this->is_leaf = true;
                        this->remaining_depth = 0;
                        this->current_depth = current_depth;
                        this->split_val = 0;
                        this->split_dimension = 0;

                        auto* inds_tmp = new int[this->size];
                        for (int i = 0; i < this->size; i++) {
                            inds_tmp[i] = current_indices[i];
                        }
                        this->indices = inds_tmp;
                    }
                    else {
                        auto* right_indices = new int[num_right];
                        auto* left_indices = new int[num_left];
                        int rc = 0;
                        int lc = 0;
                        for (int i = 0; i < this->size; i++) {
                            if(sel_tmp_right[i]) {
                                right_indices[rc] = current_indices[i];
                                rc++;
                            }
                            if(sel_tmp_left[i]) {
                                left_indices[lc] = current_indices[i];
                                lc++;
                            }
                        }

                        delete[] sel_tmp_left;
                        delete[] sel_tmp_right;

                        this->right = new kd_tree(verts, triangles, depth, right_indices, num_right, current_depth + 1);
                        delete[] right_indices;
                        this->left = new kd_tree(verts, triangles, depth, left_indices, num_left, current_depth + 1);
                        delete[] left_indices;
                    }
                }
                else {
                    this->is_leaf = true;
                    this->remaining_depth = 0;
                    this->current_depth = current_depth;

                    auto* inds_tmp = new int[this->size];
                    for (int i = 0; i < this->size; i++) {
                        inds_tmp[i] = current_indices[i];
                    }
                    this->indices = inds_tmp;
                }
            }
            else {
                this->is_leaf = true;
                this->size = size;
                this->remaining_depth = depth - current_depth;
                this->current_depth = current_depth;

                auto* inds_tmp = new int[this->size];
                for (int i = 0; i < this->size; i++) {
                    inds_tmp[i] = current_indices[i];
                }
                this->indices = inds_tmp;
            }

        }

        bool is_leaf = false;
        int split_dimension = 0;
        double split_val = 0;
        int size = 0;
        int current_depth = 0;
        int remaining_depth = 0;
        kd_tree* left = nullptr;
        kd_tree* right = nullptr;
        int* indices = nullptr;

        static void delete_kd_tree(const kd_tree* node) {
            if(!node->is_leaf) {
                delete_kd_tree(node->left);
                delete_kd_tree(node->right);
            }

            delete[] node->indices;
            delete node;
        }

        static double get_nearest_dist(const double p[NUM_DIMENSIONS], const kd_tree* node, double** verts, int** triangles, unsigned int* nearest_tr_ind,
            uint8_t* nearest_tr_count, bool* already_checked, double* proj_pt, double current_dist = 10000, const bool stop_if_lower = false, double x_dist = 0, double y_dist = 0, double z_dist = 0) {
            if(node->is_leaf) {
                for(int i = 0; i < node->size; i++) {
                    if(!already_checked[node->indices[i]]) {
                        already_checked[node->indices[i]] = true;

                        double proj_pt_tmp[NUM_DIMENSIONS];
                        double d_tmp = closest_point_on_triangle_dist(verts[triangles[node->indices[i]][0]],
                            verts[triangles[node->indices[i]][1]],
                            verts[triangles[node->indices[i]][2]],
                            p, proj_pt_tmp);

                        if(const double diff = d_tmp - current_dist; diff < (1e-10)) {
                            if(diff > -1e-10) {
                                nearest_tr_ind[nearest_tr_count[0]] = node->indices[i];
                                nearest_tr_count[0] += 1;
                            }
                            else {
                                if(stop_if_lower) {
                                    current_dist = -10;
                                    return current_dist;
                                }

                                nearest_tr_ind[0] = node->indices[i];
                                for(int k = 1; k < nearest_tr_count[0]; k++) {
                                    nearest_tr_ind[k] = 999999;
                                }
                                nearest_tr_count[0] = 1;

                                current_dist = d_tmp;
                                proj_pt[0] = proj_pt_tmp[0];
                                proj_pt[1] = proj_pt_tmp[1];
                                proj_pt[2] = proj_pt_tmp[2];
                            }
                        }
                    }
                }

                return current_dist;
            }
            else {
                double ch = p[node->split_dimension] - node->split_val;

                if(ch > 0) {
                    current_dist = get_nearest_dist(p, node->right, verts, triangles, nearest_tr_ind, nearest_tr_count, already_checked, proj_pt, current_dist, stop_if_lower, x_dist, y_dist, z_dist);
                    if(stop_if_lower && current_dist < -5) return current_dist;

                    if(node->split_dimension == 0) {
                        x_dist = ch;
                    } else if(node->split_dimension == 1) {
                        y_dist = ch;
                    } else {
                        z_dist = ch;
                    }

                    if(current_dist + (1e-10) > (x_dist*x_dist + y_dist*y_dist + z_dist*z_dist)) {
                        current_dist = get_nearest_dist(p, node->left, verts, triangles, nearest_tr_ind, nearest_tr_count, already_checked, proj_pt, current_dist, stop_if_lower, x_dist, y_dist, z_dist);
                        if(stop_if_lower && current_dist < -5) return current_dist;
                    }
                }
                else {
                    current_dist = get_nearest_dist(p, node->left, verts, triangles, nearest_tr_ind, nearest_tr_count, already_checked, proj_pt, current_dist, stop_if_lower, x_dist, y_dist, z_dist);
                    if(stop_if_lower && current_dist < -5) return current_dist;

                    if(node->split_dimension == 0) {
                        x_dist = -ch;
                    } else if(node->split_dimension == 1) {
                        y_dist = -ch;
                    } else {
                        z_dist = -ch;
                    }

                    if(current_dist + (1e-10) > (x_dist*x_dist + y_dist*y_dist + z_dist*z_dist)) {
                        current_dist = get_nearest_dist(p, node->right, verts, triangles, nearest_tr_ind, nearest_tr_count, already_checked, proj_pt, current_dist, stop_if_lower, x_dist, y_dist, z_dist);
                        if(stop_if_lower && current_dist < -5) return current_dist;
                    }
                }
            }

            return current_dist;
        }

        static void get_inside_range(const double p[NUM_DIMENSIONS], const kd_tree* node, double** verts, int** triangles, const double range,
            double* dists, bool* is_inside_range, double x_dist = 0, double y_dist = 0, double z_dist = 0) {

            if(node->is_leaf) {
                for(int i = 0; i < node->size; i++) {
                    double proj_pt_tmp[NUM_DIMENSIONS];
                    double d_tmp = closest_point_on_triangle_dist(verts[triangles[node->indices[i]][0]],
                        verts[triangles[node->indices[i]][1]],
                        verts[triangles[node->indices[i]][2]],
                        p, proj_pt_tmp);

                    dists[node->indices[i]] = d_tmp;

                    if(d_tmp < (range*range + 1e-10)) {
                        is_inside_range[node->indices[i]] = true;
                    }
                }
            }
            else {
                double ch = p[node->split_dimension] - node->split_val;
                if(ch > 0) {
                    get_inside_range(p, node->right, verts, triangles, range, dists, is_inside_range, x_dist, y_dist, z_dist);

                    if(node->split_dimension == 0) {
                        x_dist = ch;
                    } else if(node->split_dimension == 1) {
                        y_dist = ch;
                    } else {
                        z_dist = ch;
                    }

                    if(range*range + (1e-10) > (x_dist*x_dist + y_dist*y_dist + z_dist*z_dist)) {
                        get_inside_range(p, node->left, verts, triangles, range, dists, is_inside_range, x_dist, y_dist, z_dist);
                    }
                }
                else {
                    get_inside_range(p, node->left, verts, triangles, range, dists, is_inside_range, x_dist, y_dist, z_dist);

                    if(node->split_dimension == 0) {
                        x_dist = -ch;
                    } else if(node->split_dimension == 1) {
                        y_dist = -ch;
                    } else {
                        z_dist = -ch;
                    }

                    if(range*range + (1e-10) > (x_dist*x_dist + y_dist*y_dist + z_dist*z_dist)) {
                        get_inside_range(p, node->right, verts, triangles, range, dists, is_inside_range, x_dist, y_dist, z_dist);
                    }
                }
            }
        }
};


int main() {
    double** verts = nullptr;
    verts = new double*[NUM_VERTS];

    for (int h = 0; h < NUM_VERTS; h++) {
        verts[h] = new double[NUM_DIMENSIONS];
        for (int w = 0; w < NUM_DIMENSIONS; w++) {
            verts[h][w] = 0;
        }
    }

    int** triangles = nullptr;
    triangles = new int*[NUM_TRIANGLES];

    for (int h = 0; h < NUM_TRIANGLES; h++) {
        triangles[h] = new int[3];
        for (int w = 0; w < 3; w++) {
            triangles[h][w] = 0;
        }
    }

    double** verts_overs = nullptr;
    verts_overs = new double*[NUM_VERTS_OVERS];

    for (int h = 0; h < NUM_VERTS_OVERS; h++) {
        verts_overs[h] = new double[NUM_DIMENSIONS];
        for (int w = 0; w < NUM_DIMENSIONS; w++) {
            verts_overs[h][w] = 0;
        }
    }

    read_matrix_from_csv("CAD DATASETS/CAD_verts_old.csv", NUM_VERTS, NUM_DIMENSIONS, verts);
    read_matrix_from_csv("CAD DATASETS/CAD_triangles_old.csv", NUM_TRIANGLES, 3, triangles);
    read_matrix_from_csv("CAD DATASETS/CAD_verts_overs.csv", NUM_VERTS_OVERS, 3, verts_overs);


    int indices[NUM_TRIANGLES];
    for (int i = 0; i < NUM_TRIANGLES; i++) {
        indices[i] = i;
    }
    auto* my_tree = new kd_tree(verts, triangles, 18, indices, NUM_TRIANGLES);

    list<int> neigh_verts_list[NUM_TRIANGLES];

    for (int i = 0; i < NUM_VERTS_OVERS; i++) {
        bool is_inside_range[NUM_TRIANGLES];
        double dists[NUM_TRIANGLES];
        for (int h = 0; h < NUM_TRIANGLES; h++) {
            is_inside_range[h] = false;
            dists[h] = -1.0;
        }

        kd_tree::get_inside_range(verts_overs[i], my_tree, verts, triangles, cutoff_range_pert, dists, is_inside_range);

        for (int j = 0; j < NUM_TRIANGLES; j++) {
            if (is_inside_range[j]) {
                neigh_verts_list[j].push_back(i);
            }
        }

        if(i % 1000 == 0) {
            cout << i << endl;
        }
    }

    ofstream output_file ("neigh_verts_list_pert.txt");
    if (output_file.is_open()) {
        for(int i1 = 0; i1 < NUM_TRIANGLES; i1++) {

            for (auto const& element : neigh_verts_list[i1]) {
                output_file << element << " ";
            }

            output_file << "\n";
        }
    }
    else cout << "Unable to open file";
    output_file.close();
    
    
    list<int> neigh_verts_list_2[NUM_TRIANGLES];

    for (int i = 0; i < NUM_VERTS_OVERS; i++) {
        bool is_inside_range[NUM_TRIANGLES];
        double dists[NUM_TRIANGLES];
        for (int h = 0; h < NUM_TRIANGLES; h++) {
            is_inside_range[h] = false;
            dists[h] = -1.0;
        }

        kd_tree::get_inside_range(verts_overs[i], my_tree, verts, triangles, cutoff_range_smooth, dists, is_inside_range);

        for (int j = 0; j < NUM_TRIANGLES; j++) {
            if (is_inside_range[j]) {
                neigh_verts_list_2[j].push_back(i);
            }
        }

        if(i % 1000 == 0) {
            cout << i << endl;
        }
    }

    ofstream output_file2 ("neigh_verts_list_smooth.txt");
    if (output_file2.is_open()) {
        for(int i1 = 0; i1 < NUM_TRIANGLES; i1++) {

            for (auto const& element : neigh_verts_list_2[i1]) {
                output_file2 << element << " ";
            }

            output_file2 << "\n";
        }
    }
    else cout << "Unable to open file";
    output_file2.close();
    

    kd_tree::delete_kd_tree(my_tree);

    for (int h = 0; h < NUM_TRIANGLES; h++) {
        delete[] triangles[h];
    }
    delete[] triangles;
    for (int h = 0; h < NUM_VERTS; h++) {
        delete[] verts[h];
    }
    delete[] verts;
    for (int h = 0; h < NUM_VERTS_OVERS; h++) {
        delete[] verts_overs[h];
    }
    delete[] verts_overs;

    return 0;
}


