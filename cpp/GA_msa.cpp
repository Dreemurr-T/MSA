#pragma GCC optimize(2)
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <ctime>
using namespace std;
clock_t start_time,end_time;
// const int MAX_LENGTH = 140;
const int N = 100;
const int num = N * (N - 1) / 2;
const int ROW = 3;
#define Cross_rate 0.8
#define Mutation_rate 0.1
#define Iter 500

class GA
{
private:
    double F1[num]; //fitness for initial population
    double F2[num]; //fitness for mutate population
    double P1[num]; //probability for selection
    double P2[num];
    string initial_population[num][ROW] = {};
    string new_population[num][ROW] = {};
    string rand_population[num][ROW] = {}; //population for rand
    string cross_population[num][ROW] = {};
    string mutation_population[num][ROW] = {};
    string best_alignment[Iter][ROW];
    int min_cost[Iter];

public:
    int times = 0; //counter for iteration
    // vector<vector<string>> initial_population(100,vector<string>(2));
    void Genetic_A(vector<string> seq1, string seq2);
    void Init_Population(vector<string> seq1, string seq2);
    int cost(char a, char b);
    void Cal_fitness_1(); //calculate fitness for initial population
    void Cal_fitness_2(); //calculate fitness for mutation population
    void Selection();
    void Crossover();
    void Mutation();
    void Best_Solution();
    void Mixing_population();
};

void GA::Genetic_A(vector<string> seq1, string seq2)
{
    Init_Population(seq1, seq2);
    Cal_fitness_1();
    for (times = 0; times < Iter; times++)
    {
        Selection();
        Crossover();
        Mutation();
        Cal_fitness_2();
        Best_Solution();
        Mixing_population();
        // cout << times + 1 << " iteration done" << endl;
    }
    cout << "Minimun cost = " << min_cost[Iter - 1] << "\n"
         << best_alignment[Iter - 1][0] << "\n"
         << best_alignment[Iter - 1][1] << "\n"
         << best_alignment[Iter - 1][2] << endl;
}

void GA::Init_Population(vector<string> seq1, string seq2)
{
    string tmp1, tmp2, tmp3 = seq2;
    int count = 0;
    int len1 = seq2.length();
    for (int i = 0; i < N; i++)
    {
        tmp1 = seq1[i];
        int len2 = seq1[i].length();
        for (int j = i + 1; j < N; j++)
        {
            tmp2 = seq1[j];
            int len3 = seq1[j].length();
            int len = max(max(len1, len2), max(len2, len3));
            while (tmp1.length() < 1.5 * len)
            {
                int pos = rand() % tmp1.length();
                tmp1.insert(pos, "-");
            }
            initial_population[count][0] = tmp1;
            while (tmp2.length() < 1.5 * len)
            {
                int pos = rand() % tmp2.length();
                tmp2.insert(pos, "-");
            }
            initial_population[count][1] = tmp2;
            while (tmp3.length() < 1.5 * len)
            {
                int pos = rand() % tmp3.length();
                tmp3.insert(pos, "-");
            }
            initial_population[count][2] = tmp3;
            for (int pos = 1.5 * len - 1; pos >= 0; pos--)
            {
                if (initial_population[count][0][pos] == '-' && initial_population[count][1][pos] == '-' && initial_population[count][2][pos] == '-')
                {
                    initial_population[count][0].erase(pos, 1);
                    initial_population[count][1].erase(pos, 1);
                    initial_population[count][2].erase(pos, 1);
                }
            }
            tmp1 = seq1[i];
            tmp3 = seq2;
            count++;
        }
    }

    // for (int i = 0; i < N; i++)
    // {
    //     cout << initial_population[i][0] << "\n"
    //          << initial_population[i][1] << "\n"
    //          << initial_population[i][2] << endl;
    // }
}

int GA::cost(char a, char b)
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

void GA::Cal_fitness_1()
{
    memset(F1, 0, sizeof(F1));
    int pos = 0;
    for (int i = 0; i < num; i++)
    {
        for (pos = 0; pos < initial_population[i][0].length(); pos++)
        {
            F1[i] += cost(initial_population[i][0][pos], initial_population[i][1][pos]);
            F1[i] += cost(initial_population[i][1][pos], initial_population[i][2][pos]);
            F1[i] += cost(initial_population[i][0][pos], initial_population[i][2][pos]);
        }
        // F1[i] = (double)1 / F1[i];
    }
}

void GA::Selection()
{
    double sum1 = 0;
    double sum2 = 0;
    for (int i = 0; i < num; i++)
        sum1 += (double)1 / F1[i];
    for (int i = 0; i < num; i++)
    {
        sum2 += (double)1 / F1[i];
        P1[i] = (double)sum2 / sum1;
    }
    for (int i = 0; i < num; i++)
    {
        P2[i] = (double)i / num;
    }
    for (int i = 0; i < 0.8 * num; i++)
    {
        double r = (double)rand() / RAND_MAX;
        for (int j = 0; j < num; j++)
        {
            if (P1[j] >= r)
            {
                for (int k = 0; k < ROW; k++)
                    rand_population[i][k] = initial_population[j][k];
                break;
            }
        }
    }
    for (int i = 0.8 * num; i < num; i++)
    {
        double r = (double)rand() / RAND_MAX;
        for (int j = 0; j < num; j++)
        {
            if (P2[j] >= r)
            {
                for (int k = 0; k < ROW; k++)
                    rand_population[i][k] = initial_population[j][k];
                break;
            }
        }
    }
}

void GA::Crossover()
{
    int chromo1, chromo2;
    // string parent1[2];
    // string parent2[2];
    int pos = 0;
    for (int i = 0; i < num / 2; i++)
    {
        chromo1 = rand() % num;
        chromo2 = rand() % num;
        string tmp1 = rand_population[chromo1][1];
        string tmp2 = rand_population[chromo2][1];

        while (chromo1 == chromo2)
        {
            chromo2 = rand() % num;
        }

        // for (int j = 0; j < ROW; j++)
        // {
        //     parent1[j] = rand_population[chromo1][j];
        //     parent2[j] = rand_population[chromo2][j];
        // }

        double p = (double)rand() / RAND_MAX;
        if (Cross_rate >= p)
        {
            cross_population[pos][0] = rand_population[chromo1][0];
            cross_population[pos][1] = rand_population[chromo1][1];
            cross_population[pos][2] = rand_population[chromo2][2];

            cross_population[pos + 1][0] = rand_population[chromo2][0];
            cross_population[pos + 1][1] = rand_population[chromo2][1];
            cross_population[pos + 1][2] = rand_population[chromo1][2];

            while (cross_population[pos][0].length() < cross_population[pos][2].length())
            {
                int pos1 = rand() % cross_population[pos][0].length();
                int pos2 = rand() % cross_population[pos][1].length();
                int pos3 = rand() % cross_population[pos + 1][2].length();
                cross_population[pos][0].insert(pos1, "-");
                cross_population[pos][1].insert(pos2, "-");
                cross_population[pos + 1][2].insert(pos3, "-");
            }

            while (cross_population[pos][0].length() > cross_population[pos][2].length())
            {
                int pos1 = rand() % cross_population[pos][2].length();
                int pos2 = rand() % cross_population[pos + 1][0].length();
                int pos3 = rand() % cross_population[pos + 1][1].length();
                cross_population[pos][2].insert(pos1, "-");
                cross_population[pos + 1][0].insert(pos2, "-");
                cross_population[pos + 1][1].insert(pos3, "-");
            }

            pos += 2;
        }
        else
        {
            // cross_population[pos][0] = rand_population[chromo1][0];
            // cross_population[pos][1] = rand_population[chromo1][1];
            // cross_population[pos][2] = rand_population[chromo1][2];
            // cross_population[pos + 1][0] = rand_population[chromo2][0];
            // cross_population[pos + 1][1] = rand_population[chromo2][1];
            for (int j = 0; j < ROW; j++)
            {
                cross_population[pos][j] = rand_population[chromo1][j];
                cross_population[pos + 1][j] = rand_population[chromo2][j];
            }
            pos += 2;
        }
    }
}

void GA::Mutation()
{
    for (int i = 0; i < num; i++)
    {
        for (int j = 0; j < ROW; j++)
            mutation_population[i][j] = cross_population[i][j];
    }
    int point1, point2, point3;
    for (int i = 0; i < num; i++)
    {
        double r = (double)rand() / RAND_MAX;
        if (Mutation_rate >= r)
        {
            int x = rand() % 2;
            if (x % 2 == 0)
            {
                point1 = mutation_population[i][0].find("-");
                point2 = mutation_population[i][1].find("-");
                point3 = mutation_population[i][2].find("-");
                if (point1 != mutation_population[i][0].npos && point2 != mutation_population[i][1].npos && point3 != mutation_population[i][2].npos)
                {
                    mutation_population[i][0].erase(point1, 1);
                    mutation_population[i][1].erase(point2, 1);
                    mutation_population[i][2].erase(point3,1);
                }
            }
            else
            {
                point1 = rand() % mutation_population[i][0].length();
                point2 = rand() % mutation_population[i][1].length();
                point3 = rand() % mutation_population[i][2].length();
                mutation_population[i][0].insert(point1, "-");
                mutation_population[i][1].insert(point2, "-");
                mutation_population[i][2].insert(point3, "-");
            }
        }
        for (int j = mutation_population[i][0].length() - 1; j >= 0; j--)
        {
            if (mutation_population[i][0][j] == '-' && mutation_population[i][1][j] == '-' && mutation_population[i][2][j] == '-')
            {
                mutation_population[i][0].erase(j, 1);
                mutation_population[i][1].erase(j, 1);
                mutation_population[i][2].erase(j, 1);
            }
        }
    }
}

void GA::Cal_fitness_2()
{
    memset(F2, 0, sizeof(F2));
    int pos = 0;
    for (int i = 0; i < num; i++)
    {
        for (pos = 0; pos < mutation_population[i][0].length(); pos++)
        {
            F2[i] += cost(mutation_population[i][0][pos], mutation_population[i][1][pos]);
            F2[i] += cost(mutation_population[i][1][pos], mutation_population[i][2][pos]);
            F2[i] += cost(mutation_population[i][0][pos], mutation_population[i][2][pos]);
        }
    }
}

void GA::Best_Solution()
{
    string all_population[2 * num][ROW];
    double all_fitness[2 * num];
    for (int i = 0; i < 2 * num; i++)
    {
        if (i < num)
        {
            for (int j = 0; j < ROW; j++)
            {
                all_population[i][j] = initial_population[i][j];
            }
            all_fitness[i] = F1[i];
        }
        else
        {
            for (int j = 0; j < ROW; j++)
            {
                all_population[i][j] = mutation_population[i - num][j];
            }
            all_fitness[i] = F2[i - num];
        }
    }
    int max_index = 0;
    min_cost[times] = 10000;
    for (int i = 0; i < 2 * num; i++)
    {
        if (all_fitness[i] <= min_cost[times])
        {
            min_cost[times] = all_fitness[i];
            max_index = i;
        }
    }
    best_alignment[times][0] = all_population[max_index][0];
    best_alignment[times][1] = all_population[max_index][1];
    best_alignment[times][2] = all_population[max_index][2];
    // cout << "The " << times << "th best solution is\n"
    //      << all_population[max_index][0] << "\n"
    //      << all_population[max_index][1] << endl;
    // cout << "Cost is " << max_fitness << endl;
}

void GA::Mixing_population()
{
    int mut_num = 0.8 * num;
    int init_num = 0.2 * num;
    // int new_num = 0.1*N;
    int mix_index = 0;
    double F1_copy[num];
    copy(F1, F1 + num, F1_copy);
    double F2_copy[num];
    copy(F2, F2 + num, F2_copy);
    double F3[num];

    double F1_sort[num];
    int F1_preindex[num];
    copy(F1, F1 + num, F1_sort);
    sort(F1_sort, F1_sort + num);
    // reverse(F1_sort,F1_sort + N);
    for (int i = 0; i < num; i++)
    {
        F1_preindex[i] = (find(F1_copy, F1_copy + num, F1_sort[i]) - F1_copy);
        F1_copy[F1_preindex[i]] = 0;
    }
    for (int i = 0; i < init_num; i++)
    {
        for (int j = 0; j < ROW; j++)
        {
            new_population[mix_index][j] = initial_population[F1_preindex[i]][j];
        }
        F3[mix_index] = F1[F1_preindex[i]];
        mix_index++;
    }

    double F2_sort[num];
    int F2_preindex[num];
    copy(F2, F2 + num, F2_sort);
    sort(F2_sort, F2_sort + num);
    // reverse(F2_sort,F2_sort+N);
    for (int i = 0; i < num; i++)
    {
        F2_preindex[i] = (find(F2_copy, F2_copy + num, F2_sort[i]) - F2_copy);
        F2_copy[F2_preindex[i]] = 0;
    }
    for (int i = 0; i < mut_num; i++)
    {
        for (int j = 0; j < ROW; j++)
        {
            new_population[mix_index][j] = mutation_population[F2_preindex[i]][j];
        }
        F3[mix_index] = F2[F2_preindex[i]];
        mix_index++;
    }
    for (int i = 0; i < num; i++)
    {
        for (int j = 0; j < ROW; j++)
            initial_population[i][j] = new_population[i][j];
        F1[i] = F3[i];
    }
}

int main()
{
    srand((unsigned)time(NULL));
    GA test;
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
    // string seq1 = "ABCD";
    // string seq2 = "KJXXJAJKPXKJJXJKPXKJXXJAJKPXKJJXJKPXKJXXJAJKPXKJXXJAJKHXKJXXJAJKPXKJXXJAJKHXKJXX";
    start_time = clock();
    for (int i = 6; i < 8; i++)
    {
        test.Genetic_A(base, Query[i]);
    }
    end_time = clock();
    cout<<double(end_time-start_time)/CLOCKS_PER_SEC<<endl;
    // cout<<test.initial_population[0];
    return 0;
}