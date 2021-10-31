#pragma GCC optimize(2)
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <iterator>
#include <algorithm>
#include <ctime>
using namespace std;
clock_t start_time,end_time;

const int N = 200;
vector<vector<int>> cost_table(N, vector<int>(N, 0));
vector<vector<int>> route_table(N, vector<int>(N, -1));
vector<int> path;

int cost(char a, char b)
{
    if ((a == '-') || (b == '-'))
    {
        return 2;
    }
    else
    {
        if (a == b)
        {
            return 0;
        }
        else
        {
            return 3;
        }
    }
}

void build_table(string seq1, string seq2)
{
    seq1.insert(0, "-");
    seq2.insert(0, "-");
    int len1 = seq1.length();
    int len2 = seq2.length();
    int i = 0, j = 0;
    for (int i = 0; i < len1; i++)
    {
        for (int j = 0; j < len2; j++)
        {
            if (i == 0)
            {
                cost_table[i][j] = 2 * j;
                route_table[i][j] = 2;
                continue;
            }
            if (j == 0)
            {
                cost_table[i][j] = 2 * i;
                route_table[i][j] = 1;
                continue;
            }
            int move_xy = cost_table[i - 1][j - 1] + cost(seq1[i], seq2[j]);
            int move_x = cost_table[i - 1][j] + cost(seq1[i], '-');
            int move_y = cost_table[i][j - 1] + cost('-', seq2[j]);
            int min_cost = min(min(move_xy, move_x), min(move_x, move_y));
            cost_table[i][j] = min_cost;
            if (min_cost == move_xy)
                route_table[i][j] = 0;
            else if (min_cost == move_x)
                route_table[i][j] = 1;
            else
                route_table[i][j] = 2;
        }
    }
}

void build_path(string seq1, string seq2)
{
    path.clear();

    int i = seq1.length();
    int j = seq2.length();

    int pos = 0;
    while ((i > 0) || (j > 0))
    {
        if (route_table[i][j] == 0)
        {
            path.push_back(0);
            i--;
            j--;
        }
        else if (route_table[i][j] == 1)
        {
            path.push_back(1);
            i--;
        }
        else if (route_table[i][j] == 2)
        {
            path.push_back(2);
            j--;
        }
    }
    reverse(begin(path), end(path));
}

void reconstruct_seq(string seq1, string seq2)
{
    int pos = 0;
    int i = 0,j = 0;
    string new_seq1,new_seq2;
    while (pos<path.size())
    {
        if (path[pos] == 0)
        {
            new_seq1.append(1,seq1[i++]);
            new_seq2.append(1,seq2[j++]);
            pos++;
        }
        else if (path[pos] == 1)
        {
            new_seq1.append(1,seq1[i++]);
            new_seq2.append(1,'-');
            pos++;
        }
        else if (path[pos] == 2)
        {
            new_seq1.append(1,'-');
            new_seq2.append(1,seq2[j++]);
            pos++;
        }
    }
    cout << new_seq1<<"\n"<<new_seq2<<endl;
}

int main()
{
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
    // vector<string> best_align(5);
    string best_align;
    int min_cost[5] = {10000, 10000, 10000, 10000, 10000};
    getline(query, seq1);
    int start_time = clock();
    start_time = clock();
    for (int i = 0; i < 5; i++)
    {
        seq1 = Query[i];
        int len1 = seq1.length();
        for (int j = 0; j < 100; j++)
        {
            seq2 = base[j];
            int len2 = seq2.length();
            build_table(seq1, seq2);
            if (cost_table[len1][len2] < min_cost[i])
            {
                min_cost[i] = cost_table[len1][len2];
                best_align = seq2;
                build_path(seq1,seq2);
            }
            // for (auto n : cost_table)
            // {
            //     n.clear();
            // }
        }
        // build_table(seq1, best_align);
        // build_path(seq1, best_align);
        cout<<"Minimum cost = "<<min_cost[i]<<endl;
        reconstruct_seq(seq1, best_align);
        // for (auto n : cost_table)
        // {
        //     n.clear();
        // }
        // for (auto n : route_table)
        // {
        //     n.clear();
        // }
        // fill(begin(path), end(path), -1);
        // data.seekg(0);
    }
    end_time = clock();
    cout<<double(end_time - start_time)/CLOCKS_PER_SEC;
    return 0;
}
