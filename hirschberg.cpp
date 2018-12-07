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

int getAlignmentScore(string &seq1_align, string &seq2_align, map<string, int> &scoring) {
    int score = 0;
    for (int i = 0; i < seq1_align.size(); ++i) {
        if (seq1_align[i] == '-' || seq2_align[i] == '-')
            score += scoring["indel"];
        else if (seq1_align[i] != seq2_align[i])
            score += scoring["mismatch"];
        else
            score += scoring["match"];
    }
    return score;
}

vector<int> prefix(int i1, int j1, int i2, int j2, string &seq1, string &seq2, map<string, int> &scoring) {
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
    return col2;
}

vector<int> suffix(int i1, int j1, int i2, int j2, string &seq1, string &seq2, map<string, int> &scoring) {
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
    return col2;
}

void hirschberg(int i1, int j1, int i2, int j2, string &seq1, string &seq2, map<string, int> &scoring, map<int, vector<int> > &result) {
    if (j1 >= j2-1)
        return;
    vector<int> pre = prefix(i1, j1, i2, j1+((j2-j1)/2), seq1, seq2, scoring);
    vector<int> suf = suffix(i1, j1+((j2-j1)/2), i2, j2, seq1, seq2, scoring);
    int max = pre[0] + suf[0], list_begin_ix = 0;
    for (int i = 1; i < pre.size(); ++i) {
        if (max < pre[i] + suf[i]) {
            max = pre[i] + suf[i];
            list_begin_ix = i;
        }
    }
    vector<int> i_star_list;
    i_star_list.push_back(i1 + list_begin_ix);
    for (int i = list_begin_ix+1; i < pre.size(); ++i) {
        if ((pre[i] == pre[i-1] + scoring["indel"]) && max == pre[i]+suf[i])
            i_star_list.push_back(i1 + i);
    }
    int i_star = i_star_list[0];
    result[(j1+((j2-j1)/2))] = i_star_list;
    hirschberg(i1, j1, i_star, (j1+((j2-j1)/2)), seq1, seq2, scoring, result);
    hirschberg(i_star, (j1+((j2-j1)/2)), i2, j2, seq1, seq2, scoring, result);
}

void addFirstLastColumnToResult(map<int, vector<int> > &result, string &seq1, string &seq2, map<string, int> &scoring) {
    vector<int> i_star_first_col, i_star_last_col;
    vector<int> first_col_suffix = suffix(0, 0, seq1.size()-1, seq2.size()-1, seq1, seq2, scoring);
    vector<int> last_col_prefix = prefix(0, 0, seq1.size()-1, seq2.size()-1, seq1, seq2, scoring);
    int first_col_best_score = first_col_suffix[0], last_col_best_score = last_col_prefix.back();
    cout << "correct alignment score: " << last_col_best_score << endl;
    for (int i = 0; i < seq1.size(); ++i) {
        int first_col_score = first_col_suffix[i] - i;
        int last_col_score = last_col_prefix[seq1.size()-i-1] - i;
        if (i == 0 || first_col_score == first_col_best_score) {
            i_star_first_col.push_back(i);
        }
        if (i == 0 || last_col_score == last_col_best_score) {
            i_star_last_col.push_back(seq1.size() - i - 1);
        }
    }
    reverse(i_star_last_col.begin(), i_star_last_col.end());
    result[0] = i_star_first_col;
    result[seq2.size()-1] = i_star_last_col;
}

void printAlignment(map<int, vector<int> > &result, string &seq1, string &seq2, map<string, int> &scoring) {
    string seq1_align, seq2_align;
    addFirstLastColumnToResult(result, seq1, seq2, scoring);
    int next_col_first_cell = seq1.size();
    map<int, vector<int> >::reverse_iterator riter = result.rbegin();
    for (; riter != result.rend(); ++riter) {
        int col_no = riter->first;
        bool col_started = false;
        for (int i = riter->second.size()-1; i >= 0; --i) {
            if (!col_started) {
                if (riter->second[i] > next_col_first_cell)
                    continue;
                if (riter->second[i] == next_col_first_cell) {
                    seq1_align.push_back('-');
                    seq2_align.push_back(seq2[col_no+1]);
                } else if (riter->second[i] == next_col_first_cell-1) {
                    if (riter != result.rbegin()) {
                        seq1_align.push_back(seq1[next_col_first_cell]);
                        seq2_align.push_back(seq2[col_no+1]);
                    }
                }
                col_started = true;
            } else {
                seq1_align.push_back(seq1[riter->second[i+1]]);
                seq2_align.push_back('-');
            }
        }
        next_col_first_cell = riter->second[0];
    }
    reverse(seq1_align.begin(), seq1_align.end());
    reverse(seq2_align.begin(), seq2_align.end());
    cout << seq1_align << endl;
    cout << seq2_align << endl;
    int alignment_score = getAlignmentScore(seq1_align, seq2_align, scoring);
    cout << "our alignment score: " << alignment_score << endl;
}

int main() {
    string sequence1, sequence2;
    cin >> sequence1 >> sequence2;
    map<string, int> scoring;
    scoring["indel"] = -1;
    scoring["mismatch"] = -1;
    scoring["match"] = 1;
    sequence1 = "*" + sequence1;
    sequence2 = "*" + sequence2;
    if (sequence2.size() < sequence1.size()) // keeping sequence1 as the shorter sequence
        swap(sequence1, sequence2);
    cout << "seq1: " << sequence1 << endl;
    cout << "seq2: " << sequence2 << endl;
    map<int, vector<int> > result;
    hirschberg(0, 0, sequence1.size()-1, sequence2.size()-1, sequence1, sequence2, scoring, result);
    //printResult(result);
    printAlignment(result, sequence1, sequence2, scoring);
}
