#pragma GCC optimize(2)
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <iterator>
#include <algorithm>
#include <ctime>
#include <omp.h>
using namespace std;
clock_t start_time, end_time;

const int N = 200;
vector<vector<vector<int>>> cost_table(N, vector<vector<int>>(N, vector<int>(N, 1000)));
vector<vector<vector<int>>> route_table(N, vector<vector<int>>(N, vector<int>(N, -1)));
vector<int> path;

int pair_cost(char a, char b)
{
    if (a == b)
        return 0;
    else
    {
        if (a == '-' || b == '-')
            return 2;
        else
            return 3;
    }
}

int cost(char a, char b, char c)
{
    return (pair_cost(a, b) + pair_cost(a, c) + pair_cost(b, c));
}

void build_table(string seq1, string seq2, string seq3, int min_val)
{
    seq1.insert(0, "-");
    seq2.insert(0, "-");
    seq3.insert(0, "-");
    vector<int> move(7, 1000);
    int len1 = seq1.length();
    int len2 = seq2.length();
    int len3 = seq3.length();
    // for (int i = 0; i < len1; i++)
    // {
    //     cost_table[i][0][0] = 4 * i;
    //     route_table[i][0][0] = 1;
    // }
    // for (int i = 0; i < len2; i++)
    // {
    //     cost_table[0][i][0] = 4 * i;
    //     route_table[0][i][0] = 2;
    // }

    // for (int i = 0; i < len3; i++)
    // {
    //     cost_table[0][0][i] = 4 * i;
    //     route_table[0][0][i] = 3;
    // }

    for (int i = 0; i < len1; i++)
    {
        for (int j = 0; j < len2; j++)
        {
            for (int k = 0; k < len3; k++)
            {
                if (i == 0 && j == 0)
                {
                    cost_table[i][j][k] = 4 * k;
                    route_table[i][j][k] = 3;
                    // continue;
                }
                else if (j == 0 && k == 0)
                {
                    cost_table[i][j][k] = 4 * i;
                    route_table[i][j][k] = 1;
                    // continue;
                }
                else if (i == 0 && k == 0)
                {
                    cost_table[i][j][k] = 4 * j;
                    route_table[i][j][k] = 2;
                    // continue;
                }
                else if (i == 0)
                {
                    move[0] = 1000;
                    move[1] = 1000;
                    move[2] = cost_table[i][j - 1][k] + cost('-', seq2[j], '-');
                    move[3] = cost_table[i][j][k - 1] + cost('-', '-', seq3[k]);
                    move[4] = 1000;
                    move[5] = 1000;
                    move[6] = cost_table[i][j - 1][k - 1] + cost('-', seq2[j], seq3[k]);

                    int min_cost = *min_element(move.begin(), move.end());
                    cost_table[i][j][k] = min_cost;
                    route_table[i][j][k] = min_element(move.begin(), move.end()) - move.begin();
                }
                else if (j == 0)
                {
                    move[0] = 1000;
                    move[1] = cost_table[i - 1][j][k] + cost(seq1[i], '-', '-');
                    move[2] = 1000;
                    move[3] = cost_table[i][j][k - 1] + cost('-', '-', seq3[k]);
                    move[4] = 1000;
                    move[5] = cost_table[i - 1][j][k - 1] + cost(seq1[i], '-', seq3[k]);
                    move[6] = 1000;

                    int min_cost = *min_element(move.begin(), move.end());
                    cost_table[i][j][k] = min_cost;
                    route_table[i][j][k] = min_element(move.begin(), move.end()) - move.begin();
                }
                else if (k == 0)
                {
                    move[0] = 1000;
                    move[1] = cost_table[i - 1][j][k] + cost(seq1[i], '-', '-');
                    move[2] = cost_table[i][j - 1][k] + cost('-', seq2[j], '-');
                    move[3] = 1000;
                    move[4] = cost_table[i - 1][j - 1][k] + cost(seq1[i], seq2[j], '-');
                    move[5] = 1000;
                    move[6] = 1000;

                    int min_cost = *min_element(move.begin(), move.end());
                    cost_table[i][j][k] = min_cost;
                    route_table[i][j][k] = min_element(move.begin(), move.end()) - move.begin();
                }
                else if (i > 0 && j > 0 && k > 0)
                {
                    move[0] = cost_table[i - 1][j - 1][k - 1] + cost(seq1[i], seq2[j], seq3[k]);
                    move[1] = cost_table[i - 1][j][k] + cost(seq1[i], '-', '-');
                    move[2] = cost_table[i][j - 1][k] + cost('-', seq2[j], '-');
                    move[3] = cost_table[i][j][k - 1] + cost('-', '-', seq3[k]);
                    move[4] = cost_table[i - 1][j - 1][k] + cost(seq1[i], seq2[j], '-');
                    move[5] = cost_table[i - 1][j][k - 1] + cost(seq1[i], '-', seq3[k]);
                    move[6] = cost_table[i][j - 1][k - 1] + cost('-', seq2[j], seq3[k]);

                    int min_cost = *min_element(move.begin(), move.end());
                    cost_table[i][j][k] = min_cost;
                    route_table[i][j][k] = min_element(move.begin(), move.end()) - move.begin();
                    // if (cost_table[i][j][k] > min_val)
                    // {
                    //     cost_table[len1 - 1][len2 - 1][len3 - 1] = 1000;
                    //     return;
                    // }
                }
            }
        }
    }
}

void build_path(string seq1, string seq2, string seq3)
{

    path.clear();

    int i = seq1.length();
    int j = seq2.length();
    int k = seq3.length();

    // int pos = 0;

    while ((i > 0) || (j > 0) || (k > 0))
    {
        switch (route_table[i][j][k])
        {
        case 0:
        {
            path.push_back(0);
            i--;
            j--;
            k--;
            break;
        }

        case 1:
        {
            path.push_back(1);
            i--;
            break;
        }

        case 2:
        {
            path.push_back(2);
            j--;
            break;
        }
        case 3:
        {
            path.push_back(3);
            k--;
            break;
        }

        case 4:
        {
            path.push_back(4);
            i--;
            j--;
            break;
        }

        case 5:
        {
            path.push_back(5);
            i--;
            k--;
            break;
        }
        case 6:
        {
            path.push_back(6);
            j--;
            k--;
            break;
        }
        default:
            break;
        }
    }
    reverse(path.begin(), path.end());
}

void reconstruct_seq(string seq1, string seq2, string seq3)
{
    int i = 0, j = 0, k = 0;

    string new_seq1;
    string new_seq2;
    string new_seq3;
    int pos = 0;

    while (pos < path.size())
    {
        if (path[pos] == 0)
        {
            new_seq1.append(1, seq1[i++]);
            new_seq2.append(1, seq2[j++]);
            new_seq3.append(1, seq3[k++]);
        }
        else if (path[pos] == 1)
        {
            new_seq1.append(1, seq1[i++]);
            new_seq2.append(1, '-');
            new_seq3.append(1, '-');
        }
        else if (path[pos] == 2)
        {
            new_seq1.append(1, '-');
            new_seq2.append(1, seq2[j++]);
            new_seq3.append(1, '-');
        }
        else if (path[pos] == 3)
        {
            new_seq1.append(1, '-');
            new_seq2.append(1, '-');
            new_seq3.append(1, seq3[k++]);
        }
        else if (path[pos] == 4)
        {
            new_seq1.append(1, seq1[i++]);
            new_seq2.append(1, seq2[j++]);
            new_seq3.append(1, '-');
        }
        else if (path[pos] == 5)
        {
            new_seq1.append(1, seq1[i++]);
            new_seq2.append(1, '-');
            new_seq3.append(1, seq3[k++]);
        }
        else if (path[pos] == 6)
        {
            new_seq1.append(1, '-');
            new_seq2.append(1, seq2[j++]);
            new_seq3.append(1, seq3[k++]);
        }
        pos++;
    }
    cout << new_seq1 << "\n"
         << new_seq2 << "\n"
         << new_seq3 << "\n";
}

int main()
{
    // omp_set_num_threads(4);
    ifstream query("MSA_query.txt");
    ifstream data("MSA_database.txt");
    string tmp;
    getline(query, tmp);
    vector<string> base;
    vector<string> Query;
    for (int i = 0; i < 8; i++)
    {
        getline(query, tmp);
        Query.push_back(tmp);
    }
    for (int i = 0; i < 100; i++)
    {
        getline(data, tmp);
        base.push_back(tmp);
    }
    string seq1;
    string seq2;
    string seq3;
    string best_align_1;
    string best_align_2;
    int min_cost[2] = {10000, 10000};
    start_time = clock();
    for (int i = 0; i < 2; i++)
    {
        seq1 = Query[i + 6];
        int len1 = seq1.length();
        for (int j = 0; j < 100; j++)
        {
            seq2 = base[j];
            int len2 = seq2.length();
// #pragma omp parallel for
            for (int k = j + 1; k < 100; k++)
            {
                seq3 = base[k];
                int len3 = seq3.length();
                int min_val = min_cost[i];
                build_table(seq1, seq2, seq3, min_val);
                if (cost_table[len1][len2][len3] < min_cost[i])
                {
                    min_cost[i] = cost_table[len1][len2][len3];
                    best_align_1 = seq2;
                    best_align_2 = seq3;
                    build_path(seq1, seq2, seq3);
                }
                // cout << k << "th done" << endl;
                // Sleep(10);
            }
        }
        printf("%d \n", min_cost[i]);
        reconstruct_seq(seq1, best_align_1, best_align_2);
    }
    end_time = clock();
    cout << double(end_time - start_time) / CLOCKS_PER_SEC;

    // seq1 = "IPZJJLMLTKJULOSTKTJOGLKJOBLTXGKTPLUWWKOMOYJBGALJUKLGLOSVHWBPGWSLUKOBSOPLOOKUKSARPPJ";
    // seq2 = "IPZJJPLLTHUULOSTXTJOGLKJGBLLMMPJPLUWGKOMOYJBZAYUKOFLOSZHGBPHXPLXKJBXKJLAUUOJHWTWWPQ";
    // seq3 = "IPMJJLLLTHOULOSTMAJIGLKJPVLLXGKTPLTWWKOMOYJBZPYUKLILOSZHGBPGWXLZKJBSWJLPJUUMHKTRAP";
    // build_table(seq1, seq2, seq3, 1000);
    // build_path(seq1, seq2, seq3);
    // reconstruct_seq(seq1, seq2, seq3);
    // printf("%d", cost_table[seq1.length()][seq2.length()][seq3.length()]);
    return 0;
}
