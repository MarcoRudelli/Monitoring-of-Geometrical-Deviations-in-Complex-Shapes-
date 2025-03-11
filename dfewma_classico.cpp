#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <random>

using namespace std;

constexpr string sim_name = "2.51";
constexpr string in_file_name = "LBS"; //CV_quantiles //CV_moments //LBS
constexpr string out_file_name = "lbs"; //qt //mm //lbs
constexpr int NUM_COLS = 15; //10; //8; //15;

constexpr int NUM_ROWS = 150;
constexpr int NUM_SAMPLES = 39;

constexpr int minwin = 1;
constexpr int maxwin = 10;
constexpr int kmax = 50;
constexpr int m0 = 100;

constexpr int num_perm = 25000;

constexpr bool useMean = false;

constexpr double lambda = 0.01;
constexpr double alpha_err = 0.005;

void rank_vector(const double *vec, int *ranks, const int length) {
    int inds[length];
    for (int i = 0; i < length; ++i) {
        inds[i] = i;
    }

    std::sort(inds, inds + length, [&vec](size_t i1, size_t i2) {
        return vec[i1] < vec[i2];
    });

    for (int i = 0; i < length; ++i) {
        ranks[inds[i]] = i + 1;
    }
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
    //int rl[NUM_SAMPLES];

    double** t_data_mat1 = nullptr;
    t_data_mat1 = new double*[NUM_COLS];
    for (int h = 0; h < NUM_COLS; ++h) {
        t_data_mat1[h] = new double[NUM_ROWS];
        for (int w = 0; w < NUM_ROWS; ++w) {
            t_data_mat1[h][w] = 0;
        }
    }

    int** tmp1 = nullptr;
    tmp1 = new int*[NUM_ROWS];
    for (int h = 0; h < NUM_ROWS; ++h) {
        tmp1[h] = new int[NUM_COLS];
        for (int w = 0; w < NUM_COLS; ++w) {
            tmp1[h][w] = 0;
        }
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 1e-9); // Small noise distribution

    int win[kmax];
    double mu[kmax];
    double sd[kmax];
    for (int i = 0; i < kmax; ++i) {
        win[i] = i + 1;
        if (win[i] < minwin) {
            win[i] = minwin;
        }
        if (win[i] > maxwin) {
            win[i] = maxwin;
        }

        int sampleNum = i + 1 + m0;
        mu[i] = (1.0 - pow(1.0 - lambda, win[i])) / lambda * (sampleNum + 1.0) / 2.0;

        constexpr double lam = 1.0-lambda;
        double temp = pow(lam, win[i]);

        sd[i] = (1.0 - temp*temp) * (sampleNum + 1.0) * (sampleNum - 1.0) / ((2.0 * lambda - lambda * lambda)*12);
        sd[i] = sd[i] - ( (lam - temp)/(lambda*lambda) - (lam*lam - temp*temp)/(2.0*lambda*lambda-lambda*lambda*lambda) ) * (sampleNum+1)/6.0;
        sd[i] = sqrt(sd[i]);
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

    for(int s_ind = 0; s_ind < NUM_SAMPLES; ++s_ind) {
        double** data_mat = nullptr;
        data_mat = new double*[NUM_ROWS];
        for (int h = 0; h < NUM_ROWS; ++h) {
            data_mat[h] = new double[NUM_COLS];
            for (int w = 0; w < NUM_COLS; ++w) {
                data_mat[h][w] = 0;
            }
        }
        read_matrix_from_csv("/Volumes/EXTERNAL_USB/sim_" + sim_name + "/results_SIM" + sim_name + "_" + in_file_name + "/res_SIM" + sim_name + "_" + to_string(s_ind+1) + ".csv", NUM_ROWS, NUM_COLS, data_mat);
        for (int r = 0; r < NUM_ROWS; ++r) {
            for (int c = 0; c < NUM_COLS; ++c) {
                t_data_mat1[c][r] = data_mat[r][c] + dis(gen);
            }
        }

        for (int h = 0; h < NUM_ROWS; ++h) {
            delete[] data_mat[h];
        }
        delete[] data_mat;

        for (int row = 0; row < kmax; ++row) {
            for (int col = 0; col < NUM_COLS; ++col) {
                tmp2_1[row][col] = 0.0;
            }
        }

        double t_stat[kmax];
        double lim[kmax];
        for (int k = 0; k < kmax; ++k) {
            t_stat[k] = 0;
            lim[k] = 0;
        }

        for (int i = 0; i < kmax; ++i) {
            //int win_ids[win[i]];
            //for (int j = 0; j < win[i]; j++) {
            //    win_ids[win[i] - j - 1] = m0 + i - j;
            //} // (m0 + i - 0) -> (m0 + i - win[i] + 1)
            int uwl_tmp = m0 + i;
            int lwl_tmp = m0 + 1 + i - win[i];

            int rank_v_tmp[m0 + 1 + i];
            for (int col = 0; col < NUM_COLS; ++col) {
                rank_vector(t_data_mat1[col], rank_v_tmp, uwl_tmp + 1);
                for (int row = 0; row < uwl_tmp + 1; ++row) {
                    tmp1[row][col] = rank_v_tmp[row];
                }
            }

            for (int wid = uwl_tmp; wid >= lwl_tmp; --wid) {
                for (int col = 0; col < NUM_COLS; ++col) {
                    tmp2_1[i][col] += tmp1[wid][col] * working_weight[maxwin - (m0 + i - wid + 1)];
                }
            }

            double m11 = 0;

            for (int col = 0; col < NUM_COLS; ++col) {
                tmp2_1[i][col] = (tmp2_1[i][col] - mu[i]) / sd[i];
                tmp2_1[i][col] *= tmp2_1[i][col];

                if (useMean) m11 += tmp2_1[i][col];
                else m11 = m11 < tmp2_1[i][col] ? tmp2_1[i][col] : m11;
            }
            if (useMean) m11 /= NUM_COLS;

            t_stat[i] = m11;

            int sample[m0 + i + 1];
            for (int c = 0; c < m0 + i + 1; ++c) {
                sample[c] = c;
            }

            int count = 0;
            //double t_perm[num_perm + 1];
            double t_perm[num_perm];
            //t_perm[num_perm] = t_stat[i];

            // to check ids: (m0 + i - 1) -> (m0 + i - win[i] + 1)
            int end_id = max(m0+i-win[i]+1, m0);
            int st_id = m0 + i - 1;
            int min_id = end_id - win[end_id - m0] + 1;
            // m0 + i -> min_id ... m0 + i - min_id + 1
            int working_size = (m0 + i)-min_id+1;
            int working_st_id = st_id - min_id;
            int working_end_id = end_id - min_id;
            int working_tmp1[working_size][NUM_COLS];

            while (count < num_perm) {
                shuffle(sample, sample + st_id + 2, gen);
                int *working_samp = sample + min_id;
                bool perm_ok = true;

                if (i > 0) {
                    for (int j = 0; j < working_size; ++j) {
                        for (int c = 0; c < NUM_COLS; ++c) {
                            working_tmp1[j][c] = tmp1[working_samp[j]][c];
                        }
                    }

                    for (int j = working_st_id; j >= working_end_id; --j) {
                        // j -> (j - win[i - (st_id - j + 1)] + 1)

                        for (int k = j; k >= 0; --k) {
                            for (int c = 0; c < NUM_COLS; ++c) {
                                if (working_tmp1[j+1][c] < working_tmp1[k][c]) {
                                    --working_tmp1[k][c];
                                }
                                //else {
                                //    if (working_tmp1[j+1][c] == working_tmp1[k][c]) {
                                //        working_tmp1[k][c] -= 0.5;
                                //    }
                                //}
                            }
                        }

                        int tmp_i = i - (working_st_id - j + 1);
                        int tmp_win_lim = j - win[tmp_i] + 1;

                        double t_stat_tmp;

                        m11 = 0;
                        double tmp2_tmp1[NUM_COLS];
                        for (double & c : tmp2_tmp1) {
                            c = 0.0;
                        }
                        int ww_ind = maxwin;
                        double ww_tmp;
                        for (int wid = j; wid >= (tmp_win_lim); --wid) {
                            ww_tmp = working_weight[--ww_ind];
                            int* p1 = working_tmp1[wid];
                            for (int col = 0; col < NUM_COLS; ++col) {
                                tmp2_tmp1[col] += p1[col] * ww_tmp;
                            }
                        }
                        for (double & c : tmp2_tmp1) {
                            c = (c - mu[tmp_i]) / sd[tmp_i];
                            c *= c;

                            if (useMean) m11 += c;
                            else m11 = m11 < c ? c : m11;
                        }
                        if (useMean) m11 /= NUM_COLS;

                        t_stat_tmp = m11;


                        if (t_stat_tmp > lim[tmp_i]) {
                            perm_ok = false;
                            break;
                        }
                    }
                }


                if (perm_ok) {
                    m11 = 0;
                    double tmp2_tmp1[NUM_COLS];
                    for (double & c : tmp2_tmp1) {
                        c = 0.0;
                    }
                    int ww_ind = maxwin;
                    double ww_tmp;
                    for (int wid = m0 + i; wid >= (lwl_tmp); --wid) {
                        ww_tmp = working_weight[--ww_ind];
                        int* p1 = tmp1[sample[wid]];
                        for (int col = 0; col < NUM_COLS; ++col) {
                            tmp2_tmp1[col] += p1[col] * ww_tmp;
                        }
                    }
                    for (double & c : tmp2_tmp1) {
                        c = (c - mu[i]) / sd[i];
                        c *= c;

                        if (useMean) m11 += c;
                        else m11 = m11 < c ? c : m11;
                    }
                    if (useMean) m11 /= NUM_COLS;

                    t_perm[count] = m11;

                    if (count % 1000 == 0) cout << count << endl;
                    ++count;
                }
            }

            sort(t_perm, t_perm + (num_perm)); //+ 1));

            //lim[i] = t_perm[static_cast<int>(ceil((num_perm + 1) * (1 - alpha_err)))] + 1e-12;
            lim[i] = t_perm[static_cast<int>(ceil((num_perm) * (1 - alpha_err)))] + 1e-12;

            cout << i << endl;
            cout << lim[i] << " " << t_stat[i] << endl;
            if (t_stat[i] > lim[i]) {
                break;
            }
        }


        ofstream output_file ("/Volumes/EXTERNAL_USB/sim_" + sim_name + "/res_surv_old_" + out_file_name + "/lim_" + to_string(s_ind+1) + ".txt");
        ofstream output_file2 ("/Volumes/EXTERNAL_USB/sim_" + sim_name + "/res_surv_old_" + out_file_name + "/t_stat_" + to_string(s_ind+1) + ".txt");
        if (output_file.is_open()) {
            for (int i2 = 0; i2 < kmax; ++i2) {
                output_file << lim[i2] << "\n";
                output_file2 << t_stat[i2] << "\n";
            }
        }
        else cout << "Unable to open file";
        output_file.close();
        output_file2.close();

        ofstream output_file4 ("/Volumes/EXTERNAL_USB/sim_" + sim_name + "/res_surv_old_" + out_file_name + "/tmp2_" + to_string(s_ind+1) + ".txt");
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
    }

    for (int h = 0; h < kmax; ++h) {
        delete[] tmp2_1[h];
    }
    delete[] tmp2_1;

    for (int h = 0; h < NUM_COLS; ++h) {
        delete[] t_data_mat1[h];
    }
    delete[] t_data_mat1;

    for (int h = 0; h < NUM_ROWS; ++h) {
        delete[] tmp1[h];
    }
    delete[] tmp1;

    return 0;
}













