#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>
// #include <omp.h>

using namespace std;

constexpr int NUM_VERTS = 38920;
constexpr int NUM_TRIANGLES = 67362;

constexpr int NUM_VERTS_FATHER = 26850;
constexpr int NUM_TRIANGLES_FATHER = 53696;

int NUM_TRIANGLES_BATCH = 100;
constexpr int NUM_DIMENSIONS = 3;

constexpr double sigma2_boot = 0.04;// 0.0;

constexpr int W = 25;

int NUM_VERTS_SCAN = 0;
int NUM_TRIANGLES_SCAN = 0;
constexpr int NUM_BOOT_PER_CENTR = 75;
constexpr double AREA_BASELINE = 62002.07; // AREA CAD
int NUM_POINTS_OVERSAMP = 0;

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

void read_lists_from_csv(const string& file_name_count, const string& file_name_lists, int dim, short* counts, vector<int>* lists) {
    ifstream file(file_name_count);
    ifstream file2(file_name_lists);

    string line;
    string line2;
    int row = 0;
    while (row < dim && getline(file, line) && getline(file2, line2)) {
        stringstream ss2(line2);
        string cell2;

        counts[row] = static_cast<short>(stoi(line));

        int col = 0;
        while (getline(ss2, cell2, ',') && col < counts[row]) {
            lists[row].push_back(stoi(cell2));

            col++;
        }
        row++;
    }

    file.close();
    file2.close();
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

    const double d1 = ab[0] * ap[0] + ab[1] * ap[1] + ab[2] * ap[2];
    const double d2 = ac[0] * ap[0] + ac[1] * ap[1] + ac[2] * ap[2];
    if (d1 <= 0.0 && d2 <= 0.0){
        proj_point[0] = a[0];
        proj_point[1] = a[1];
        proj_point[2] = a[2];
    } else {
        double bp[NUM_DIMENSIONS];
        bp[0] = p[0] - b[0];
        bp[1] = p[1] - b[1];
        bp[2] = p[2] - b[2];

        const double d3 = ab[0] * bp[0] + ab[1] * bp[1] + ab[2] * bp[2];
        const double d4 = ac[0] * bp[0] + ac[1] * bp[1] + ac[2] * bp[2];
        if (d3 >= 0.0 && d4 <= d3){
            proj_point[0] = b[0];
            proj_point[1] = b[1];
            proj_point[2] = b[2];
        } else {
            double cp[NUM_DIMENSIONS];
            cp[0] = p[0] - c[0];
            cp[1] = p[1] - c[1];
            cp[2] = p[2] - c[2];
            const double d5 = ab[0] * cp[0] + ab[1] * cp[1] + ab[2] * cp[2];
            const double d6 = ac[0] * cp[0] + ac[1] * cp[1] + ac[2] * cp[2];
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

        static double get_nearest_dist(const double p[NUM_DIMENSIONS], const kd_tree* node, double** verts, int** triangles, vector<unsigned int>* nearest_tr_ind,
            vector<double>* nearest_tr_dist, bool* already_checked, double* proj_pt, double current_dist = 10000, double x_dist = 0, double y_dist = 0, double z_dist = 0) {
            if(node->is_leaf) {
                for(int i = 0; i < node->size; i++) {
                    if(!already_checked[node->indices[i]]) {
                        already_checked[node->indices[i]] = true;

                        double proj_pt_tmp[NUM_DIMENSIONS];
                        double d_tmp = closest_point_on_triangle_dist(verts[triangles[node->indices[i]][0]],
                            verts[triangles[node->indices[i]][1]],
                            verts[triangles[node->indices[i]][2]],
                            p, proj_pt_tmp);

                        if(const double diff = d_tmp - current_dist; diff < 1e-10) {
                            if(diff >= 0) {
                                nearest_tr_ind[0].push_back(node->indices[i]);
                                nearest_tr_dist[0].push_back(d_tmp);
                            }
                            else {
                                current_dist = d_tmp;

                                for (int j = 0; j < nearest_tr_dist[0].size(); j++) {
                                    if (nearest_tr_dist[0][j] > (current_dist + 1e-10)) {
                                        nearest_tr_dist[0].erase(nearest_tr_dist[0].begin() + j);
                                        nearest_tr_ind[0].erase(nearest_tr_ind[0].begin() + j);
                                        j--;
                                    }
                                }

                                nearest_tr_ind[0].push_back(node->indices[i]);
                                nearest_tr_dist[0].push_back(d_tmp);

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
                    current_dist = get_nearest_dist(p, node->right, verts, triangles, nearest_tr_ind, nearest_tr_dist, already_checked, proj_pt, current_dist, x_dist, y_dist, z_dist);

                    if(node->split_dimension == 0) {
                        x_dist = ch;
                    } else if(node->split_dimension == 1) {
                        y_dist = ch;
                    } else {
                        z_dist = ch;
                    }

                    if(current_dist + (1e-10) > (x_dist*x_dist + y_dist*y_dist + z_dist*z_dist)) {
                        current_dist = get_nearest_dist(p, node->left, verts, triangles, nearest_tr_ind, nearest_tr_dist, already_checked, proj_pt, current_dist, x_dist, y_dist, z_dist);
                    }
                }
                else {
                    current_dist = get_nearest_dist(p, node->left, verts, triangles, nearest_tr_ind, nearest_tr_dist, already_checked, proj_pt, current_dist, x_dist, y_dist, z_dist);

                    if(node->split_dimension == 0) {
                        x_dist = -ch;
                    } else if(node->split_dimension == 1) {
                        y_dist = -ch;
                    } else {
                        z_dist = -ch;
                    }

                    if(current_dist + (1e-10) > (x_dist*x_dist + y_dist*y_dist + z_dist*z_dist)) {
                        current_dist = get_nearest_dist(p, node->right, verts, triangles, nearest_tr_ind, nearest_tr_dist, already_checked, proj_pt, current_dist, x_dist, y_dist, z_dist);
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

// ATTENZIONE: POTREI ANCHE CALCOLARE PRIMA TUTTE LE INFO NECESSARIE PER LA FUNZIONE closest_point_on_triangle_dist (cioè ab, ac, bc...)

int main() {
    double** verts = nullptr;
    verts = new double*[NUM_VERTS];
    for (int h = 0; h < NUM_VERTS; h++) {
        verts[h] = new double[NUM_DIMENSIONS];
        for (int w = 0; w < NUM_DIMENSIONS; w++) {
            verts[h][w] = 0;
        }
    }

    double** verts_father = nullptr;
    verts_father = new double*[NUM_VERTS_FATHER];
    for (int h = 0; h < NUM_VERTS_FATHER; h++) {
        verts_father[h] = new double[NUM_DIMENSIONS];
        for (int w = 0; w < NUM_DIMENSIONS; w++) {
            verts_father[h][w] = 0;
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

    int** triangles_father = nullptr;
    triangles_father = new int*[NUM_TRIANGLES_FATHER];
    for (int h = 0; h < NUM_TRIANGLES_FATHER; h++) {
        triangles_father[h] = new int[3];
        for (int w = 0; w < 3; w++) {
            triangles_father[h][w] = 0;
        }
    }

    read_matrix_from_csv("CAD DATASETS/CAD_verts_overs.csv", NUM_VERTS, NUM_DIMENSIONS, verts);
    read_matrix_from_csv("CAD DATASETS/CAD_triangles_overs.csv", NUM_TRIANGLES, 3, triangles);

    read_matrix_from_csv("CAD DATASETS/CAD_verts_old.csv", NUM_VERTS_FATHER, NUM_DIMENSIONS, verts_father);
    read_matrix_from_csv("CAD DATASETS/CAD_triangles_old.csv", NUM_TRIANGLES_FATHER, 3, triangles_father);

    auto* indices = new int[NUM_TRIANGLES_FATHER];
    for (int i = 0; i < NUM_TRIANGLES_FATHER; i++) {
        indices[i] = i;
    }
    auto* my_tree = new kd_tree(verts_father, triangles_father, 15, indices, NUM_TRIANGLES_FATHER);

    auto* geod_neighbors = new vector<int>[NUM_TRIANGLES];
    auto* geod_neighbors_count = new short[NUM_TRIANGLES];
    read_lists_from_csv("CAD DATASETS/geod_neig_tr_count.txt",
        "CAD DATASETS/geod_neig_tr_list.txt", NUM_TRIANGLES, geod_neighbors_count, geod_neighbors);

    auto* children = new vector<int>[NUM_TRIANGLES_FATHER];
    auto* children_count = new short[NUM_TRIANGLES_FATHER];
    read_lists_from_csv("CAD DATASETS/child_tr_count.txt",
        "CAD DATASETS/child_tr_list.txt", NUM_TRIANGLES_FATHER, children_count, children);

    // bool* already_checked = new bool[NUM_TRIANGLES_FATHER];

    // #pragma omp parallel for num_threads(4)
    int egg_ind = 0;
    for (int egg_ind2 = 0; egg_ind2 < 5000; egg_ind2++) {
        if(egg_ind2 % 50 == 0){
            egg_ind += 100;
        }
        egg_ind++;
        
        cout << "Indice corrente: " << egg_ind << endl;

        NUM_VERTS_SCAN = read_dim_1("/Volumes/EXTERNAL_USB/parts_simulated_SIM2.51/SCAN_verts_" + to_string(egg_ind) + ".csv");
        NUM_TRIANGLES_SCAN = read_dim_1("/Volumes/EXTERNAL_USB/parts_simulated_SIM2.51/SCAN_triangles_" + to_string(egg_ind) + ".csv");

        double** verts_scan = nullptr;
        verts_scan = new double*[NUM_VERTS_SCAN];

        for (int h = 0; h < NUM_VERTS_SCAN; h++) {
            verts_scan[h] = new double[NUM_DIMENSIONS];
            for (int w = 0; w < NUM_DIMENSIONS; w++) {
                verts_scan[h][w] = 0;
            }
        }

        int** triangles_scan = nullptr;
        triangles_scan = new int*[NUM_TRIANGLES_SCAN];

        for (int h = 0; h < NUM_TRIANGLES_SCAN; h++) {
            triangles_scan[h] = new int[3];
            for (int w = 0; w < 3; w++) {
                triangles_scan[h][w] = 0;
            }
        }

        NUM_POINTS_OVERSAMP = NUM_TRIANGLES_SCAN*NUM_BOOT_PER_CENTR;

        read_matrix_from_csv("/Volumes/EXTERNAL_USB/parts_simulated_SIM2.51/SCAN_verts_" + to_string(egg_ind) + ".csv", NUM_VERTS_SCAN, NUM_DIMENSIONS, verts_scan, true);
        read_matrix_from_csv("/Volumes/EXTERNAL_USB/parts_simulated_SIM2.51/SCAN_triangles_" + to_string(egg_ind) + ".csv", NUM_TRIANGLES_SCAN, 3, triangles_scan, true);

        double tr_areas_scan[NUM_TRIANGLES_SCAN];
        double tr_areas_sum = 0.0;
        for (int h = 0; h < NUM_TRIANGLES_SCAN; h++) {
            double tmp1[NUM_DIMENSIONS];
            double tmp2[NUM_DIMENSIONS];
            double tmp3[NUM_DIMENSIONS];

            for (int w = 0; w < NUM_DIMENSIONS; w++) {
                tmp1[w] = verts_scan[triangles_scan[h][1]][w] - verts_scan[triangles_scan[h][0]][w];
                tmp2[w] = verts_scan[triangles_scan[h][2]][w] - verts_scan[triangles_scan[h][0]][w];
            }
            
            cross(tmp1, tmp2, tmp3);
            

            tr_areas_scan[h] = 0;
            for (double w : tmp3) {
                tr_areas_scan[h] += w*w;
            }
            tr_areas_scan[h] = sqrt(tr_areas_scan[h])/2;

            tr_areas_sum += tr_areas_scan[h];
        }
        double mult = tr_areas_sum / AREA_BASELINE;
        

        auto* weights = new double[NUM_TRIANGLES];
        auto* weighted_dists = new double[NUM_TRIANGLES];
        for (int i = 0; i < NUM_TRIANGLES; i++) {
            weights[i] = 0.0;
            weighted_dists[i] = 0.0;
        }

        random_device rd{};
        mt19937 gen{rd()};
        normal_distribution<> norm_dist_sc{0, sqrt(sigma2_boot)};
        uniform_real_distribution<> unif_dist{0, 1};

        int centr_generator_range[NUM_TRIANGLES_SCAN];
        centr_generator_range[0] = static_cast<int>(ceil(exp(log(tr_areas_scan[0]) - log(tr_areas_sum) + log(NUM_POINTS_OVERSAMP) + log(mult))));

        for (int i = 1; i < NUM_TRIANGLES_SCAN; i++) {
            centr_generator_range[i] = centr_generator_range[i - 1] +
                static_cast<int>(ceil(exp(log(tr_areas_scan[i]) - log(tr_areas_sum) + log(NUM_POINTS_OVERSAMP) + log(mult))));
        }

        NUM_POINTS_OVERSAMP = centr_generator_range[NUM_TRIANGLES_SCAN - 1];

        int i_tr = 0;
        bool first_iteration = true;

        for(int b = 0; b < NUM_POINTS_OVERSAMP; b++) {
            int i = b % (NUM_TRIANGLES_BATCH * 1000);

            if (first_iteration) {
                for (int j = 0; j < NUM_TRIANGLES_SCAN; j++) {
                    if (centr_generator_range[j] > b) {
                        i_tr = j;
                        break;
                    }
                }
                first_iteration = false;
            }
            else {
                if(b == centr_generator_range[i_tr]) {
                    i_tr++;
                }
            }

            double tmp1 = norm_dist_sc(gen);
            double tmp2 = norm_dist_sc(gen);
            double tmp3 = norm_dist_sc(gen);

            double tmp4 = unif_dist(gen);
            double tmp5 = unif_dist(gen);
            if((tmp4 + tmp5) > 1) {
                tmp4 = 1 - tmp4;
                tmp5 = 1 - tmp5;
            }

            double sp[NUM_DIMENSIONS];
            sp[0] = verts_scan[triangles_scan[i_tr][0]][0] +
                (verts_scan[triangles_scan[i_tr][1]][0] - verts_scan[triangles_scan[i_tr][0]][0]) * tmp4 +
                    (verts_scan[triangles_scan[i_tr][2]][0] - verts_scan[triangles_scan[i_tr][0]][0]) * tmp5;
            sp[1] = verts_scan[triangles_scan[i_tr][0]][1] +
                (verts_scan[triangles_scan[i_tr][1]][1] - verts_scan[triangles_scan[i_tr][0]][1]) * tmp4 +
                    (verts_scan[triangles_scan[i_tr][2]][1] - verts_scan[triangles_scan[i_tr][0]][1]) * tmp5;
            sp[2] = verts_scan[triangles_scan[i_tr][0]][2] +
                (verts_scan[triangles_scan[i_tr][1]][2] - verts_scan[triangles_scan[i_tr][0]][2]) * tmp4 +
                    (verts_scan[triangles_scan[i_tr][2]][2] - verts_scan[triangles_scan[i_tr][0]][2]) * tmp5;

            double oversampled_grid[NUM_DIMENSIONS];
            oversampled_grid[0] = sp[0] + tmp1;
            oversampled_grid[1] = sp[1] + tmp2;
            oversampled_grid[2] = sp[2] + tmp3;

            //unsigned int nearest_triangle_ind[12];
            //for(unsigned int & k : nearest_triangle_ind) {
            //    k = 999999;
            //}
            //uint8_t nearest_triangle_count[1];
            //nearest_triangle_count[0] = 0;
            vector<unsigned int> nearest_triangle_ind[1];
            vector<double> nearest_triangle_dist[1];

            bool already_checked[NUM_TRIANGLES_FATHER];
            for (int ac = 0; ac < NUM_TRIANGLES_FATHER; ac++) {
                already_checked[ac] = false;
            }
            double proj_pt_cad[NUM_DIMENSIONS];

            kd_tree::get_nearest_dist(oversampled_grid, my_tree, verts_father, triangles_father, nearest_triangle_ind, nearest_triangle_dist, already_checked, proj_pt_cad);

            vector<unsigned int> nearest_triangle_ind_children;
            for (int nti = 0; nti < nearest_triangle_ind[0].size(); nti++) {
                for (int child_ind = 0; child_ind < children_count[nearest_triangle_ind[0][nti]]; child_ind++) {
                    double proj_tmp[NUM_DIMENSIONS];
                    double d_tmp = closest_point_on_triangle_dist(verts[triangles[children[nearest_triangle_ind[0][nti]][child_ind]][0]],
                            verts[triangles[children[nearest_triangle_ind[0][nti]][child_ind]][1]],
                            verts[triangles[children[nearest_triangle_ind[0][nti]][child_ind]][2]],
                            proj_pt_cad, proj_tmp);

                    if (d_tmp < 1e-10) {
                        nearest_triangle_ind_children.push_back(children[nearest_triangle_ind[0][nti]][child_ind]);
                    }
                }
            }

            double cur_dist = 100000.0;
            for (int nti = 0; nti < nearest_triangle_ind_children.size(); nti++) {
                for (int neigh_ind = 0; neigh_ind < geod_neighbors_count[nearest_triangle_ind_children[nti]]; neigh_ind++) {
                    double proj_pt_scan[NUM_DIMENSIONS];
                    double d_tmp = closest_point_on_triangle_dist(verts[triangles[geod_neighbors[nearest_triangle_ind_children[nti]][neigh_ind]][0]],
                            verts[triangles[geod_neighbors[nearest_triangle_ind_children[nti]][neigh_ind]][1]],
                            verts[triangles[geod_neighbors[nearest_triangle_ind_children[nti]][neigh_ind]][2]],
                            sp, proj_pt_scan);
                    if (cur_dist > d_tmp) {
                        cur_dist = d_tmp;
                    }
                }
            }

            for(int i_id = 0; i_id < nearest_triangle_ind_children.size(); i_id++) {
                weighted_dists[nearest_triangle_ind_children[i_id]] += sqrt(cur_dist) / static_cast<double>(nearest_triangle_ind_children.size());
                weights[nearest_triangle_ind_children[i_id]] += 1.0 / static_cast<double>(nearest_triangle_ind_children.size());
            }

            //if(i == (NUM_TRIANGLES_BATCH * 1000 - 1) | b == (NUM_POINTS_OVERSAMP - 1)) {
            //    cout << "Num punti " << b + 1 << " su " << NUM_POINTS_OVERSAMP << endl;
            //}
        }

        ofstream my_file2_2_fin ("/Volumes/EXTERNAL_USB/RESULTS SIMULAZIONE_SIM2.51/weights_" + to_string(egg_ind) + ".txt");
        ofstream my_file2_3_fin ("/Volumes/EXTERNAL_USB/RESULTS SIMULAZIONE_SIM2.51/weighted_dists_" + to_string(egg_ind) + ".txt");
        if (my_file2_2_fin.is_open()) {
            for(int i1 = 0; i1 < NUM_TRIANGLES; i1++) {
                my_file2_2_fin << weights[i1] << "\n";
                my_file2_3_fin << weighted_dists[i1] << "\n";
            }
        }
        else cout << "Unable to open file";
        my_file2_2_fin.close();
        my_file2_3_fin.close();


        for (int h = 0; h < NUM_TRIANGLES_SCAN; h++) {
            delete[] triangles_scan[h];
        }
        delete[] triangles_scan;
        for (int h = 0; h < NUM_VERTS_SCAN; h++) {
            delete[] verts_scan[h];
        }
        delete[] verts_scan;

        delete[] weights;
        delete[] weighted_dists;
    }

    kd_tree::delete_kd_tree(my_tree);
    delete[] indices;

    for (int h = 0; h < NUM_TRIANGLES; h++) {
        delete[] triangles[h];
    }
    delete[] triangles;
    for (int h = 0; h < NUM_VERTS; h++) {
        delete[] verts[h];
    }
    delete[] verts;

    for (int h = 0; h < NUM_TRIANGLES_FATHER; h++) {
        delete[] triangles_father[h];
    }
    delete[] triangles_father;
    for (int h = 0; h < NUM_VERTS_FATHER; h++) {
        delete[] verts_father[h];
    }
    delete[] verts_father;

    // delete[] already_checked;

    delete[] geod_neighbors;
    delete[] geod_neighbors_count;

    return 0;
}

