#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <random>
#include <filesystem>

using namespace std;

constexpr int NUM_COLS = 16840;
constexpr int NUM_ROWS = 150;
constexpr int NUM_SAMPLES = 10;

constexpr int minwin = 1;
constexpr int maxwin = 10;
constexpr int kmax = 50;
constexpr int m0 = 100;

constexpr int num_perm = 25000;

constexpr int num_stats = 4 * 2; // numero statistiche di controllo

constexpr double lambda = 0.01;
constexpr double alpha_err = 0.005; // ARL0 = 200

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
            ++col;
        }
        ++row;
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
            ++col;
        }
        ++row;
    }
    file.close();
}

int main() {
    double** data_mat1 = nullptr;
    data_mat1 = new double*[NUM_ROWS];
    for (int h = 0; h < NUM_ROWS; ++h) {
        data_mat1[h] = new double[NUM_COLS];
        for (int w = 0; w < NUM_COLS; ++w) {
            data_mat1[h][w] = 0;
        }
    }
    double** data_mat2 = nullptr;
    data_mat2 = new double*[NUM_ROWS];
    for (int h = 0; h < NUM_ROWS; ++h) {
        data_mat2[h] = new double[NUM_COLS];
        for (int w = 0; w < NUM_COLS; ++w) {
            data_mat2[h][w] = 0;
        }
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1e-9); // Small noise distribution

    int win[kmax];
    for (int i = 0; i < kmax; ++i) {
        win[i] = i + 1;
        if (win[i] < minwin) {
            win[i] = minwin;
        }
        if (win[i] > maxwin) {
            win[i] = maxwin;
        }
    }
    double working_weight[maxwin];
    for (int i = 0; i < maxwin; ++i) {
        working_weight[i] = maxwin-1.0-i;
        working_weight[i] = pow(1.0-lambda, working_weight[i]);
    }

    double** tmp2_1 = nullptr;
    tmp2_1 = new double*[kmax];
    for (int h = 0; h < kmax; ++h) {
        tmp2_1[h] = new double[NUM_COLS];
        for (int w = 0; w < NUM_COLS; ++w) {
            tmp2_1[h][w] = 0;
        }
    }
    double** tmp2_2 = nullptr;
    tmp2_2 = new double*[kmax];
    for (int h = 0; h < kmax; ++h) {
        tmp2_2[h] = new double[NUM_COLS];
        for (int w = 0; w < NUM_COLS; ++w) {
            tmp2_2[h][w] = 0;
        }
    }

    for(int s_ind = 0; s_ind < NUM_SAMPLES; ++s_ind) {
        cout << s_ind + 1 << endl;
        double** data_mat = nullptr;
        data_mat = new double*[NUM_ROWS];
        for (int h = 0; h < NUM_ROWS; ++h) {
            data_mat[h] = new double[NUM_COLS];
            for (int w = 0; w < NUM_COLS; ++w) {
                data_mat[h][w] = 0;
            }
        }

        read_matrix_from_csv("area_warp_smooth_mats/area_warp_mat_smooth_" + to_string(s_ind+1) + ".csv", NUM_ROWS, NUM_COLS, data_mat);
        for (int r = 0; r < NUM_ROWS; ++r) {
            for (int c = 0; c < NUM_COLS; ++c) {
                data_mat1[r][c] = data_mat[r][c] + dis(gen);
            }
        }

        read_matrix_from_csv("backward_dist_mats/backward_dist_smooth_" + to_string(s_ind+1) + ".csv", NUM_ROWS, NUM_COLS, data_mat);
        for (int r = 0; r < NUM_ROWS; ++r) {
            for (int c = 0; c < NUM_COLS; ++c) {
                data_mat2[r][c] = data_mat[r][c] + dis(gen);
            }
        }

        read_matrix_from_csv("forward_dist_mats/forward_dist_smooth_" + to_string(s_ind+1) + ".csv", NUM_ROWS, NUM_COLS, data_mat);
        for (int r = 0; r < NUM_ROWS; ++r) {
            for (int c = 0; c < NUM_COLS; ++c) {
                double tmp = dis(gen);
                if (data_mat2[r][c] < data_mat[r][c] + tmp) {
                    data_mat2[r][c] = data_mat[r][c] + tmp;
                }
            }
        }

        for (int h = 0; h < NUM_ROWS; ++h) {
            delete[] data_mat[h];
        }
        delete[] data_mat;

        for (int row = 0; row < kmax; ++row) {
            for (int col = 0; col < NUM_COLS; ++col) {
                tmp2_1[row][col] = 0.0;
                tmp2_2[row][col] = 0.0;
            }
        }

        double t_stat[num_stats][kmax];
        double lo_lim[num_stats][kmax];
        double up_lim[num_stats][kmax];
        for (int ns = 0; ns < num_stats; ++ns) {
            for (int k = 0; k < kmax; ++k) {
                t_stat[ns][k] = 0;
                lo_lim[ns][k] = 0;
                up_lim[ns][k] = 0;
            }
        }

        for (int i = 0; i < kmax; ++i) {
            int uwl_tmp = m0 + i;
            int lwl_tmp = m0 + 1 + i - win[i];

            double samp_mean1[NUM_COLS];
            double samp_mean_sq1[NUM_COLS];
            double samp_mean2[NUM_COLS];
            double samp_mean_sq2[NUM_COLS];
            for (int col = 0; col < NUM_COLS; ++col) {
                samp_mean1[col] = 0;
                samp_mean_sq1[col] = 0;
                samp_mean2[col] = 0;
                samp_mean_sq2[col] = 0;
            }
            for (int owid = 0; owid <= uwl_tmp; ++owid) {
                for (int col = 0; col < NUM_COLS; ++col) {
                    samp_mean1[col] += data_mat1[owid][col];
                    samp_mean2[col]  += data_mat2[owid][col];
                    samp_mean_sq1[col] += data_mat1[owid][col] * data_mat1[owid][col];
                    samp_mean_sq2[col] += data_mat2[owid][col] * data_mat2[owid][col];
                }
            }

            double samp_mean1_tmp[NUM_COLS];
            double samp_mean_sq1_tmp[NUM_COLS];
            double samp_mean2_tmp[NUM_COLS];
            double samp_mean_sq2_tmp[NUM_COLS];
            for (int col = 0; col < NUM_COLS; ++col) {
                samp_mean1_tmp[col] = samp_mean1[col];
                samp_mean2_tmp[col] = samp_mean2[col];
                samp_mean_sq1_tmp[col] = samp_mean_sq1[col];
                samp_mean_sq2_tmp[col] = samp_mean_sq2[col];
            }

            for (int wid = uwl_tmp; wid >= lwl_tmp; --wid) {
                for (int col = 0; col < NUM_COLS; ++col) {
                    samp_mean1_tmp[col] -= data_mat1[wid][col];
                    samp_mean2_tmp[col] -= data_mat2[wid][col];
                    samp_mean_sq1_tmp[col] -= data_mat1[wid][col] * data_mat1[wid][col];
                    samp_mean_sq2_tmp[col] -= data_mat2[wid][col] * data_mat2[wid][col];
                }
            }
            
            double samp_sd1[NUM_COLS];
            double samp_sd2[NUM_COLS];
            for (int col = 0; col < NUM_COLS; ++col) {
                samp_mean1_tmp[col] /= lwl_tmp;
                samp_mean_sq1_tmp[col] /= lwl_tmp;
                samp_mean2_tmp[col] /= lwl_tmp;
                samp_mean_sq2_tmp[col] /= lwl_tmp;

                samp_sd1[col] = sqrt(samp_mean_sq1_tmp[col] - samp_mean1_tmp[col]*samp_mean1_tmp[col]);
                samp_sd2[col] = sqrt(samp_mean_sq2_tmp[col] - samp_mean2_tmp[col]*samp_mean2_tmp[col]);
            }
            for (int wid = uwl_tmp; wid >= lwl_tmp; --wid) {
                for (int col = 0; col < NUM_COLS; ++col) {

                    tmp2_1[i][col] += (data_mat1[wid][col] - samp_mean1_tmp[col]) / samp_sd1[col] * working_weight[maxwin - (m0 + i - wid + 1)];
                    tmp2_2[i][col] += (data_mat2[wid][col] - samp_mean2_tmp[col]) / samp_sd2[col] * working_weight[maxwin - (m0 + i - wid + 1)];
                }
            }

            double m11 = 0;
            double m21 = 0;
            double m31 = 0;
            double m41 = 0;
            double m12 = 0;
            double m22 = 0;
            double m32 = 0;
            double m42 = 0;
            for (int col = 0; col < NUM_COLS; ++col) {
                m11 += tmp2_1[i][col];
                m21 += tmp2_1[i][col]*tmp2_1[i][col];
                m31 += tmp2_1[i][col]*tmp2_1[i][col]*tmp2_1[i][col];
                m41 += tmp2_1[i][col]*tmp2_1[i][col]*tmp2_1[i][col]*tmp2_1[i][col];
                m12 += tmp2_2[i][col];
                m22 += tmp2_2[i][col]*tmp2_2[i][col];
                m32 += tmp2_2[i][col]*tmp2_2[i][col]*tmp2_2[i][col];
                m42 += tmp2_2[i][col]*tmp2_2[i][col]*tmp2_2[i][col]*tmp2_2[i][col];
            }
            m11 /= NUM_COLS;
            m21 /= NUM_COLS;
            m31 /= NUM_COLS;
            m41 /= NUM_COLS;
            m12 /= NUM_COLS;
            m22 /= NUM_COLS;
            m32 /= NUM_COLS;
            m42 /= NUM_COLS;

            t_stat[0][i] = m11;
            t_stat[1][i] = m21 - m11*m11;
            t_stat[2][i] = (m31 - 3*m11*t_stat[1][i] - m11*m11*m11) / pow(t_stat[1][i], 1.5);
            t_stat[3][i] = (m41 - 4*m11*m31 + 6*m11*m11*t_stat[1][i] + 3*m11*m11*m11*m11) / (t_stat[1][i] * t_stat[1][i]);

            t_stat[4][i] = m12;
            t_stat[5][i] = m22 - m12*m12;
            t_stat[6][i] = (m32 - 3*m12*t_stat[5][i] - m12*m12*m12) / pow(t_stat[5][i], 1.5);
            t_stat[7][i] = (m42 - 4*m12*m32 + 6*m12*m12*t_stat[5][i] + 3*m12*m12*m12*m12) / (t_stat[5][i] * t_stat[5][i]);

            int sample[m0 + i + 1];
            for (int c = 0; c < m0 + i + 1; ++c) {
                sample[c] = c;
            }

            int count = 0;
            double t_perm[num_stats][num_perm];

            int end_id = max(m0+i-win[i]+1, m0);
            int st_id = m0 + i - 1;
            int min_id = end_id - win[end_id - m0] + 1;
            // m0 + i -> min_id ... m0 + i - min_id + 1
            int working_size = (m0 + i)-min_id+1;
            int working_st_id = st_id - min_id;
            int working_end_id = end_id - min_id;
            double working_tmp1[working_size][NUM_COLS];
            double working_tmp2[working_size][NUM_COLS];

            while (count < num_perm) {
                shuffle(sample, sample + st_id + 2, gen);
                int *working_samp = sample + min_id;
                bool perm_ok = true;

                if (i > 0) {
                    for (int j = 0; j < working_size; ++j) {
                        for (int c = 0; c < NUM_COLS; ++c) {
                            working_tmp1[j][c] = data_mat1[working_samp[j]][c];
                            working_tmp2[j][c] = data_mat2[working_samp[j]][c];
                        }
                    }
                    
                    double samp_mean1_cpy[NUM_COLS];
                    double samp_mean2_cpy[NUM_COLS];
                    double samp_mean_sq1_cpy[NUM_COLS];
                    double samp_mean_sq2_cpy[NUM_COLS];
                    for (int c = 0; c < NUM_COLS; ++c) {
                        samp_mean1_cpy[c] = samp_mean1[c];
                        samp_mean2_cpy[c] = samp_mean2[c];
                        samp_mean_sq1_cpy[c] = samp_mean_sq1[c];
                        samp_mean_sq2_cpy[c] = samp_mean_sq2[c];
                    }

                    for (int j = working_st_id; j >= working_end_id; --j) {
                        // j -> (j - win[i - (st_id - j + 1)] + 1)

                        int tmp_i = i - (working_st_id - j + 1);
                        int tmp_win_lim = j - win[tmp_i] + 1;

                        double t_stat_tmp[num_stats];

                        m11 = 0;
                        m21 = 0;
                        m31 = 0;
                        m41 = 0;
                        m12 = 0;
                        m22 = 0;
                        m32 = 0;
                        m42 = 0;
                        double tmp2_tmp1[NUM_COLS];
                        double tmp2_tmp2[NUM_COLS];
                        for (int c = 0; c < NUM_COLS; ++c) {
                            tmp2_tmp1[c] = 0.0;
                            tmp2_tmp2[c] = 0.0;

                            samp_mean1_cpy[c] -= working_tmp1[j+1][c];
                            samp_mean2_cpy[c] -= working_tmp2[j+1][c];
                            samp_mean_sq1_cpy[c] -= working_tmp1[j+1][c] * working_tmp1[j+1][c];
                            samp_mean_sq2_cpy[c] -= working_tmp2[j+1][c] * working_tmp2[j+1][c];
                            
                            samp_mean1_tmp[c] = samp_mean1_cpy[c];
                            samp_mean2_tmp[c] = samp_mean2_cpy[c];
                            samp_mean_sq1_tmp[c] = samp_mean_sq1_cpy[c];
                            samp_mean_sq2_tmp[c] = samp_mean_sq2_cpy[c];
                        }
                        int ww_ind = maxwin;
                        double ww_tmp;
                        for (int wid = j; wid >= (tmp_win_lim); --wid) {
                            ww_tmp = working_weight[--ww_ind];
                            double* p1 = working_tmp1[wid];
                            double* p2 = working_tmp2[wid];
                            for (int col = 0; col < NUM_COLS; ++col) {
                                samp_mean1_tmp[col] -= p1[col];
                                samp_mean2_tmp[col] -= p2[col];
                                samp_mean_sq1_tmp[col] -= p1[col] * p1[col];
                                samp_mean_sq2_tmp[col] -= p2[col] * p2[col];
                            }
                        }
                        for (int col = 0; col < NUM_COLS; ++col) {
                            samp_mean1_tmp[col] /= (tmp_win_lim + min_id);
                            samp_mean_sq1_tmp[col] /= (tmp_win_lim + min_id);
                            samp_mean2_tmp[col] /= (tmp_win_lim + min_id);
                            samp_mean_sq2_tmp[col] /= (tmp_win_lim + min_id);

                            samp_sd1[col] = sqrt(samp_mean_sq1_tmp[col] - samp_mean1_tmp[col]*samp_mean1_tmp[col]);
                            samp_sd2[col] = sqrt(samp_mean_sq2_tmp[col] - samp_mean2_tmp[col]*samp_mean2_tmp[col]);
                        }

                        ww_ind = maxwin;
                        for (int wid = j; wid >= (tmp_win_lim); --wid) {
                            ww_tmp = working_weight[--ww_ind];
                            double* p1 = working_tmp1[wid];
                            double* p2 = working_tmp2[wid];
                            for (int col = 0; col < NUM_COLS; ++col) {
                                tmp2_tmp1[col] += (p1[col] - samp_mean1_tmp[col]) / samp_sd1[col]  * ww_tmp;
                                tmp2_tmp2[col] += (p2[col] - samp_mean2_tmp[col]) / samp_sd2[col]  * ww_tmp;
                            }
                        }

                        for (int c = 0; c < NUM_COLS; ++c) {
                            m11 += tmp2_tmp1[c];
                            m21 += tmp2_tmp1[c]*tmp2_tmp1[c];
                            m31 += tmp2_tmp1[c]*tmp2_tmp1[c]*tmp2_tmp1[c];
                            m41 += tmp2_tmp1[c]*tmp2_tmp1[c]*tmp2_tmp1[c]*tmp2_tmp1[c];

                            m12 += tmp2_tmp2[c];
                            m22 += tmp2_tmp2[c]*tmp2_tmp2[c];
                            m32 += tmp2_tmp2[c]*tmp2_tmp2[c]*tmp2_tmp2[c];
                            m42 += tmp2_tmp2[c]*tmp2_tmp2[c]*tmp2_tmp2[c]*tmp2_tmp2[c];
                        }
                        m11 /= NUM_COLS;
                        m21 /= NUM_COLS;
                        m31 /= NUM_COLS;
                        m41 /= NUM_COLS;
                        m12 /= NUM_COLS;
                        m22 /= NUM_COLS;
                        m32 /= NUM_COLS;
                        m42 /= NUM_COLS;

                        t_stat_tmp[0] = m11;
                        t_stat_tmp[1] = m21 - m11*m11;
                        t_stat_tmp[2] = (m31 - 3*m11*t_stat_tmp[1] - m11*m11*m11) / pow(t_stat_tmp[1], 1.5);
                        t_stat_tmp[3] = (m41 - 4*m11*m31 + 6*m11*m11*t_stat_tmp[1] + 3*m11*m11*m11*m11) / (t_stat_tmp[1] * t_stat_tmp[1]);

                        t_stat_tmp[4] = m12;
                        t_stat_tmp[5] = m22 - m12*m12;
                        t_stat_tmp[6] = (m32 - 3*m12*t_stat_tmp[5] - m12*m12*m12) / pow(t_stat_tmp[5], 1.5);
                        t_stat_tmp[7] = (m42 - 4*m12*m32 + 6*m12*m12*t_stat_tmp[5] + 3*m12*m12*m12*m12) / (t_stat_tmp[5] * t_stat_tmp[5]);

                        for (int sid = 0; sid < num_stats; ++sid) {
                            if (t_stat_tmp[sid] < lo_lim[sid][tmp_i] || t_stat_tmp[sid] > up_lim[sid][tmp_i]) {
                                perm_ok = false;
                                break;
                            }
                        }

                        if (!perm_ok) {
                            break;
                        }
                    }
                }


                if (perm_ok) {
                    m11 = 0;
                    m21 = 0;
                    m31 = 0;
                    m41 = 0;
                    m12 = 0;
                    m22 = 0;
                    m32 = 0;
                    m42 = 0;
                    double tmp2_tmp1[NUM_COLS];
                    double tmp2_tmp2[NUM_COLS];
                    for (int c = 0; c < NUM_COLS; ++c) {
                        tmp2_tmp1[c] = 0.0;
                        tmp2_tmp2[c] = 0.0;

                        samp_mean1_tmp[c] = samp_mean1[c];
                        samp_mean2_tmp[c] = samp_mean2[c];
                        samp_mean_sq1_tmp[c] = samp_mean_sq1[c];
                        samp_mean_sq2_tmp[c] = samp_mean_sq2[c];
                    }
                    int ww_ind = maxwin;
                    double ww_tmp;
                    for (int wid = m0 + i; wid >= (lwl_tmp); --wid) {
                        ww_tmp = working_weight[--ww_ind];
                        double* p1 = data_mat1[sample[wid]];
                        double* p2 = data_mat2[sample[wid]];
                        for (int col = 0; col < NUM_COLS; ++col) {
                            samp_mean1_tmp[col] -= p1[col];
                            samp_mean2_tmp[col] -= p2[col];
                            samp_mean_sq1_tmp[col] -= p1[col] * p1[col];
                            samp_mean_sq2_tmp[col] -= p2[col] * p2[col];
                        }
                    }
                    for (int col = 0; col < NUM_COLS; ++col) {
                        samp_mean1_tmp[col] /= lwl_tmp;
                        samp_mean_sq1_tmp[col] /= lwl_tmp;
                        samp_mean2_tmp[col] /= lwl_tmp;
                        samp_mean_sq2_tmp[col] /= lwl_tmp;

                        samp_sd1[col] = sqrt(samp_mean_sq1_tmp[col] - samp_mean1_tmp[col]*samp_mean1_tmp[col]);
                        samp_sd2[col] = sqrt(samp_mean_sq2_tmp[col] - samp_mean2_tmp[col]*samp_mean2_tmp[col]);
                    }

                    ww_ind = maxwin;
                    for (int wid = m0 + i; wid >= (lwl_tmp); --wid) {
                        ww_tmp = working_weight[--ww_ind];
                        double* p1 = data_mat1[sample[wid]];
                        double* p2 = data_mat2[sample[wid]];
                        for (int col = 0; col < NUM_COLS; ++col) {
                            tmp2_tmp1[col] += (p1[col] - samp_mean1_tmp[col]) / samp_sd1[col] * ww_tmp;
                            tmp2_tmp2[col] += (p2[col] - samp_mean2_tmp[col]) / samp_sd2[col] * ww_tmp;
                        }
                    }

                    for (int c = 0; c < NUM_COLS; ++c) {
                        m11 += tmp2_tmp1[c];
                        m21 += tmp2_tmp1[c]*tmp2_tmp1[c];
                        m31 += tmp2_tmp1[c]*tmp2_tmp1[c]*tmp2_tmp1[c];
                        m41 += tmp2_tmp1[c]*tmp2_tmp1[c]*tmp2_tmp1[c]*tmp2_tmp1[c];

                        m12 += tmp2_tmp2[c];
                        m22 += tmp2_tmp2[c]*tmp2_tmp2[c];
                        m32 += tmp2_tmp2[c]*tmp2_tmp2[c]*tmp2_tmp2[c];
                        m42 += tmp2_tmp2[c]*tmp2_tmp2[c]*tmp2_tmp2[c]*tmp2_tmp2[c];
                    }
                    m11 /= NUM_COLS;
                    m21 /= NUM_COLS;
                    m31 /= NUM_COLS;
                    m41 /= NUM_COLS;
                    m12 /= NUM_COLS;
                    m22 /= NUM_COLS;
                    m32 /= NUM_COLS;
                    m42 /= NUM_COLS;

                    t_perm[0][count] = m11;
                    t_perm[1][count] = m21 - m11*m11;
                    t_perm[2][count] = (m31 - 3*m11*t_perm[1][count] - m11*m11*m11) / pow(t_perm[1][count], 1.5);
                    t_perm[3][count] = (m41 - 4*m11*m31 + 6*m11*m11*t_perm[1][count] + 3*m11*m11*m11*m11) / (t_perm[1][count] * t_perm[1][count]);

                    t_perm[4][count] = m12;
                    t_perm[5][count] = m22 - m12*m12;
                    t_perm[6][count] = (m32 - 3*m12*t_perm[5][count] - m12*m12*m12) / pow(t_perm[5][count], 1.5);
                    t_perm[7][count] = (m42 - 4*m12*m32 + 6*m12*m12*t_perm[5][count] + 3*m12*m12*m12*m12) / (t_perm[5][count] * t_perm[5][count]);

                    if (count % 1000 == 0) cout << count << endl;
                    ++count;
                }
            }

            int num_ok = num_perm;
            int alpha = num_perm - 1;
            double t_perm_sort[num_stats][num_perm];
            for (int sid = 0; sid < num_stats; ++sid) {
                for (int p = 0; p < num_perm; ++p) {
                    t_perm_sort[sid][p] = t_perm[sid][p];
                }
                sort(t_perm_sort[sid], t_perm_sort[sid] + (num_perm));
            }

            double lims_up_prec[num_stats];
            double lims_down_prec[num_stats];
            int num_ok_prec;
            int alpha_prec;

            double lims_up[num_stats];
            double lims_down[num_stats];
            for (int sid = 0; sid < num_stats; ++sid) {
                lims_up[sid] = t_perm_sort[sid][alpha] + 1e-12;
                lims_down[sid] = t_perm_sort[sid][num_perm - alpha - 1] - 1e-12;
            }

            while(num_ok > (num_perm)*(1-alpha_err)){
                alpha_prec = alpha;
                --alpha;
                for (int sid = 0; sid < num_stats; ++sid) {
                    lims_up_prec[sid] = lims_up[sid];
                    lims_down_prec[sid] = lims_down[sid];

                    lims_up[sid] = t_perm_sort[sid][alpha] + 1e-12;
                    lims_down[sid] = t_perm_sort[sid][num_perm - alpha - 1] - 1e-12;
                }

                num_ok_prec = num_ok;
                num_ok = 0;
                for (int c = 0; c < num_perm; ++c) {
                    bool is_ok = true;
                    for (int sid = 0; sid < num_stats; ++sid) {
                        if (t_perm[sid][c] < lims_down[sid] || t_perm[sid][c] > lims_up[sid]) {
                            is_ok = false;
                            break;
                        }
                    }
                    if (is_ok) {
                        ++num_ok;
                    }
                }
            }

            for (int sid = 0; sid < num_stats; ++sid) {
                lo_lim[sid][i] = lims_down_prec[sid];
                up_lim[sid][i] = lims_up_prec[sid];
            }

            bool is_ok = true;
            for (int sid = 0; sid < num_stats; ++sid) {
                if (t_stat[sid][i] < lo_lim[sid][i] || t_stat[sid][i] > up_lim[sid][i]) {
                    is_ok = false;
                    break;
                }
            }

            cout << i << " " << num_ok_prec*1.0/(num_perm) << " " << (alpha_prec+1)*1.0/(num_perm) << endl;
            cout << lo_lim[0][i] << " " << lo_lim[1][i] << " " << lo_lim[2][i] << " " << lo_lim[3][i] << endl;
            cout << t_stat[0][i] << " " << t_stat[1][i] << " " << t_stat[2][i] << " " << t_stat[3][i] << endl;
            cout << up_lim[0][i] << " " << up_lim[1][i] << " " << up_lim[2][i] << " " << up_lim[3][i] << endl;
            cout  << endl;
            cout << lo_lim[4][i] << " " << lo_lim[5][i] << " " << lo_lim[6][i] << " " << lo_lim[7][i] << endl;
            cout << t_stat[4][i] << " " << t_stat[5][i] << " " << t_stat[6][i] << " " << t_stat[7][i] << endl;
            cout << up_lim[4][i] << " " << up_lim[5][i] << " " << up_lim[6][i] << " " << up_lim[7][i] << endl;
            cout << endl;

            if (!is_ok) {
                break;
            }
        }

        std::filesystem::create_directory("res_surv_new");
        ofstream output_file ("res_surv_new/up_lim_" + to_string(s_ind+1) + ".txt");
        ofstream output_file2 ("res_surv_new/lo_lim_" + to_string(s_ind+1) + ".txt");
        ofstream output_file3 ("res_surv_new/t_stat_" + to_string(s_ind+1) + ".txt");
        if (output_file.is_open()) {
            for(int i1 = 0; i1 < num_stats; ++i1) {
                for (int i2 = 0; i2 < kmax; ++i2) {
                    output_file << up_lim[i1][i2] << " ";
                    output_file2 << lo_lim[i1][i2] << " ";
                    output_file3 << t_stat[i1][i2] << " ";
                }

                output_file << "\n";
                output_file2 << "\n";
                output_file3 << "\n";
            }
        }
        else cout << "Unable to open file";
        output_file.close();
        output_file2.close();
        output_file3.close();

        ofstream output_file4 ("res_surv_new/tmp2_1_" + to_string(s_ind+1) + ".txt");
        if (output_file4.is_open()) {
            for(int i1 = 0; i1 < kmax; ++i1) {
                for (int i2 = 0; i2 < NUM_COLS; ++i2) {
                    output_file4 << tmp2_1[i1][i2] << " ";
                }
                output_file4 << "\n";
            }
        }
        else cout << "Unable to open file";
        output_file4.close();
        ofstream output_file5 ("res_surv_new/tmp2_2_" + to_string(s_ind+1) + ".txt");
        if (output_file5.is_open()) {
            for(int i1 = 0; i1 < kmax; ++i1) {
                for (int i2 = 0; i2 < NUM_COLS; ++i2) {
                    output_file5 << tmp2_2[i1][i2] << " ";
                }
                output_file5 << "\n";
            }
        }
        else cout << "Unable to open file";
        output_file5.close();
    }

    for (int h = 0; h < kmax; ++h) {
        delete[] tmp2_1[h];
    }
    delete[] tmp2_1;
    for (int h = 0; h < kmax; ++h) {
        delete[] tmp2_2[h];
    }
    delete[] tmp2_2;

    for (int h = 0; h < NUM_ROWS; ++h) {
        delete[] data_mat1[h];
    }
    delete[] data_mat1;
    for (int h = 0; h < NUM_ROWS; ++h) {
        delete[] data_mat2[h];
    }
    delete[] data_mat2;

    return 0;
}














