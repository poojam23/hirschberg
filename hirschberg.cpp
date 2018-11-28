#include<bits/stdc++.h>
using namespace std;

void printResult(map<int, vector<int> > &result) {
    for (map<int, vector<int> >::iterator iter = result.begin(); iter != result.end(); ++iter) {
        cout << iter->first << ": ";
        for (int i = 0; i < iter->second.size(); ++i)
            cout << iter->second[i] << ",";
        cout << endl;
    }
}

void printVector(vector<int> &vec) {
    for (int i = 0; i < vec.size(); ++i)
        cout << vec[i] << ",";
    cout << endl;
}

pair<vector<int>, vector<int> > prefix(int i1, int j1, int i2, int j2, string &seq1, string &seq2, map<string, int> &scoring) {
    cout << "prefix: " << i1 << ", " << j1 << ", " << i2 << ", " << j2 << endl; 
    // 2-column solution
    int col_length = i2 - i1 + 1;
    vector<int> col1, col2(col_length);
    for (int i = 0; i < col_length; ++i)
        col1.push_back(-i);
    for (int j = j1 + 1; j <= j2; j++) {
        for (int i = 0; i < col_length; ++i) {
           int max = col1[i] + scoring["indel"];
           if (i > 0) { // j > j1 always
                max = max > (col2[i-1] + scoring["indel"]) ? max : (col2[i-1] + scoring["indel"]);
                int character_compare_score = 0;
                if (seq1[i1+i] == seq2[j])
                    character_compare_score = scoring["match"];
                else
                    character_compare_score = scoring["mismatch"];
                max = max > (col1[i-1] + character_compare_score) ? max : (col1[i-1] + character_compare_score);
           }
           col2[i] = max;
        }
        if (j < j2)
            col1 = col2;
    }
    printVector(col2);
    //vector<vector<int> > return_vec(2);
    //return_vec[0] = col1;
    //return_vec[1] = col2;
    //return return_vec;
    return make_pair(col1, col2);
}

vector<int> suffix(int i1, int j1, int i2, int j2, string &seq1, string &seq2, map<string, int> &scoring) {
    cout << "suffix: " << i1 << ", " << j1 << ", " << i2 << ", " << j2 << endl; 
    // 2-column solution
    int col_length = i2 - i1 + 1;
    vector<int> col1(col_length), col2(col_length);
    for (int i = col_length - 1; i >= 0; --i)
        col1[col_length - 1 - i] = -i;
    for (int j = j2 - 1; j >= j1; --j) {
        for (int i = col_length - 1; i >= 0; --i) {
            int max = col1[i] + scoring["indel"];
            if (i < col_length - 1) {
                max = max > (col2[i+1] + scoring["indel"]) ? max : (col2[i+1] + scoring["indel"]);
                int character_compare_score = 0;
                if (seq1[i1+i+1] == seq2[j+1])
                    character_compare_score = scoring["match"];
                else
                    character_compare_score = scoring["mismatch"];
                max = max > (col1[i+1] + character_compare_score) ? max : (col1[i+1] + character_compare_score);
            }
            col2[i] = max;
        }
        col1 = col2;
    }
    printVector(col2);
    return col2;
}

void hirschberg(int i1, int j1, int i2, int j2, string &seq1, string &seq2, map<string, int> &scoring, map<int, vector<int> > &result) {
    cout << "hirsch: " << i1 << ", " << j1 << ", " << i2 << ", " << j2 << endl; 
    if (j1 >= j2-1)
        return;
    vector<int> pre = prefix(i1, j1, i2, j1+((j2-j1)/2), seq1, seq2, scoring).second;
    vector<int> suf = suffix(i1, j1+((j2-j1)/2), i2, j2, seq1, seq2, scoring);
    int max = pre[0] + suf[0];
    vector<int> i_star_list;
    i_star_list.push_back(i1);
    for (int i = 1; i < pre.size(); ++i) {
        if (max < pre[i]+suf[i]) {
            max = pre[i] + suf[i];
            i_star_list.clear();
            i_star_list.push_back(i + i1);
        } else if (max == pre[i] + suf[i]) {
            i_star_list.push_back(i + i1);
        }
    }
    int i_star = i_star_list[0];
    result[(j1+((j2-j1)/2))] = i_star_list;
    printResult(result);
    hirschberg(i1, j1, i_star, (j1+((j2-j1)/2)), seq1, seq2, scoring, result);
    hirschberg(i_star, (j1+((j2-j1)/2)), i2, j2, seq1, seq2, scoring, result);
}

void printAlignment(map<int, vector<int> > &result, string &seq1, string &seq2, map<string, int> &scoring) {
    cout << "printing alignment..." << endl;
    string seq1_align, seq2_align;
    int prev_col_last_cell = -1;
    for (map<int, vector<int> >::iterator iter = result.begin(); iter != result.end(); ++iter) {
        int col_no = iter->first;
        for (int i = 0; i < iter->second.size(); ++i) {
            if (i == 0) {
                if (prev_col_last_cell == -1) {
                    if (iter->second[i] == 0)
                        seq1_align.push_back('-');
                    else
                        seq1_align.push_back(seq1[iter->second[i]]);
                    seq2_align.push_back(seq2[col_no]);
                } else if (prev_col_last_cell == iter->second[i]) {
                    seq1_align.push_back('-');
                    seq2_align.push_back(seq2[col_no]);
                } else if (prev_col_last_cell == iter->second[i]-1) {
                    seq1_align.push_back(seq1[iter->second[i]]);
                    seq2_align.push_back(seq2[col_no]);
                }
            } else {
                seq1_align.push_back(seq1[iter->second[i]]);
                seq2_align.push_back('-');
            }
        }
        prev_col_last_cell = iter->second.back();
        cout << prev_col_last_cell << endl;
    }
    if (prev_col_last_cell == seq1.size()-1) {
        seq1_align.push_back('-');
        seq2_align.push_back(seq2.back());
    } else {
        cout << "last column stuff: " << prev_col_last_cell << endl;
        pair<vector<int>, vector<int> > prefix_pair = prefix(0, 0, seq1.size()-1, seq2.size()-1, seq1, seq2, scoring);
        vector<int> col1 = prefix_pair.first, col2 = prefix_pair.second;
        cout << "col1: ";
        printVector(col1);
        cout << "col2: ";
        printVector(col2);
        int second_last_value = col1[prev_col_last_cell];
        int diagonal_down = col2[prev_col_last_cell+1];
        cout << second_last_value << " " << diagonal_down << endl;
        int seq1_start = 0;
        if ((seq2.back() == seq1[prev_col_last_cell+1]) && (diagonal_down == second_last_value+1)) {
            cout << "picking 1" << endl;
            seq1_align.push_back(seq1[prev_col_last_cell+1]);
            seq2_align.push_back(seq2.back());
            seq1_start = prev_col_last_cell + 2;
        } else if (diagonal_down == second_last_value-1) {
            cout << "picking 2" << endl;
            seq1_align.push_back(seq1[prev_col_last_cell+1]);
            seq2_align.push_back(seq2.back());
            seq1_start = prev_col_last_cell + 2;
        } else if ((diagonal_down == col2[prev_col_last_cell]-1) && (col2[prev_col_last_cell] == second_last_value-1)) {
            cout << "picking 3" << endl;
            seq1_align.push_back('-');
            seq2_align.push_back(seq2.back());
            seq1_start = prev_col_last_cell + 1;
        }
        for (int i = seq1_start; i < seq1.size(); ++i) {
            seq1_align.push_back(seq1[i]);
            seq2_align.push_back('-');
        }
    }
    cout << seq1_align << endl;
    cout << seq2_align << endl;
}

int main() {
    string sequence1, sequence2;
    cin >> sequence1 >> sequence2;
    map<string, int> scoring;
    // cin >> scoring["indel"] >> scoring["mismatch"] >> scoring["match"];
    scoring["indel"] = -1;
    scoring["mismatch"] = -1;
    scoring["match"] = 1;
    sequence1 = "*" + sequence1;
    sequence2 = "*" + sequence2;
    if (sequence2.size() < sequence1.size()) // keeping sequence1 to be the shorter sequence
        swap(sequence1, sequence2);
    cout << "seq1: " << sequence1 << endl;
    cout << "seq2: " << sequence2 << endl;
    map<int, vector<int> > result;
    hirschberg(0, 0, sequence1.size()-1, sequence2.size()-1, sequence1, sequence2, scoring, result);
    for (map<int, vector<int> >::iterator iter = result.begin(); iter != result.end(); ++iter) {
        cout << "for column " << iter->first << endl;
        for (int i = 0; i < iter->second.size(); ++i) {
            cout << iter->second[i] << ", ";
        }
        cout << endl;
    }
    printAlignment(result, sequence1, sequence2, scoring);
}
