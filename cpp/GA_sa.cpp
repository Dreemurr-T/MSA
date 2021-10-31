#pragma GCC optimize(2)
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <omp.h>
#include <ctime>
using namespace std;

clock_t start_time,end_time;
// const int MAX_LENGTH = 140;
const int N = 100;
const int ROW = 2;
#define Cross_rate 0.8
#define Mutation_rate 0.1
#define Iter 500

class GA
{
private:
    double F1[N]; //fitness for initial population
    double F2[N]; //fitness for mutate population
    double P1[N]; //probability for selection
    double P2[N];
    string initial_population[N][ROW] = {};
    string new_population[N][ROW] = {};
    string rand_population[N][ROW] = {}; //population for rand
    string cross_population[N][ROW] = {};
    string mutation_population[N][ROW] = {};
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
    }
    cout << "Minimun cost = " << min_cost[Iter - 1] << "\n"
         << best_alignment[Iter - 1][0] << "\n"
         << best_alignment[Iter - 1][1] << endl;
}

void GA::Init_Population(vector<string> seq1, string seq2)
{
    string new_seq2[N];
    for (int i = 0; i < N; i++)
    {
        new_seq2[i] = seq2;
    }
    int len2 = seq2.length();
    for (int i = 0; i < N; i++)
    {
        int len1 = seq1[i].length();
        int len = max(len1, len2);
        while (seq1[i].length() < 1.1 * len)
        {
            int pos = rand() % seq1[i].length();
            seq1[i].insert(pos, "-");
        }
        initial_population[i][0] = seq1[i];
        while (new_seq2[i].length() < 1.1 * len)
        {
            int pos = rand() % new_seq2[i].length();
            new_seq2[i].insert(pos, "-");
        }
        initial_population[i][1] = new_seq2[i];
        for (int pos = 1.1 * len - 1; pos >= 0; pos--)
        {
            if (initial_population[i][0][pos] == '-' && initial_population[i][1][pos] == '-')
            {
                initial_population[i][0].erase(pos, 1);
                initial_population[i][1].erase(pos, 1);
            }
        }
    }

    // for (int i = 0; i < N; i++)
    // {
    //     cout << initial_population[i][0] << "\n"
    //          << initial_population[i][1] << endl;
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
    for (int i = 0; i < N; i++)
    {
        for (pos = 0; pos < initial_population[i][0].length(); pos++)
        {
            F1[i] += cost(initial_population[i][0][pos], initial_population[i][1][pos]);
        }
        // F1[i] = (double)1 / F1[i];
    }
}

void GA::Selection()
{
    double sum1 = 0;
    double sum2 = 0;
    for (int i = 0; i < N; i++)
        sum1 += (double)1 / F1[i];
    for (int i = 0; i < N; i++)
    {
        sum2 += (double)1 / F1[i];
        P1[i] = (double)sum2 / sum1;
    }
    for (int i = 0; i < N; i++)
    {
        P2[i] = (double)i / N;
    }
    for (int i = 0; i < 0.8 * N; i++)
    {
        double r = (double)rand() / RAND_MAX;
        for (int j = 0; j < N; j++)
        {
            if (P1[j] >= r)
            {
                for (int k = 0; k < ROW; k++)
                    rand_population[i][k] = initial_population[j][k];
                break;
            }
        }
    }
    for (int i = 0.8 * N; i < N; i++)
    {
        double r = (double)rand() / RAND_MAX;
        for (int j = 0; j < N; j++)
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
    for (int i = 0; i < N / 2; i++)
    {
        chromo1 = rand() % N;
        chromo2 = rand() % N;
        while (chromo1 == chromo2)
            chromo2 = rand() % N;
        // for (int j = 0; j < ROW; j++)
        // {
        //     parent1[j] = rand_population[chromo1][j];
        //     parent2[j] = rand_population[chromo2][j];
        // }

        double p = (double)rand() / RAND_MAX;
        if (Cross_rate >= p)
        {
            cross_population[pos][0] = rand_population[chromo1][0];
            cross_population[pos][1] = rand_population[chromo2][1];

            cross_population[pos + 1][0] = rand_population[chromo2][0];
            cross_population[pos + 1][1] = rand_population[chromo1][1];

            while (cross_population[pos][0].length() < cross_population[pos][1].length())
            {
                int pos1 = rand() % cross_population[pos][0].length();
                int pos2 = rand() % cross_population[pos + 1][1].length();
                cross_population[pos][0].insert(pos1, "-");
                cross_population[pos + 1][1].insert(pos2, "-");
            }

            while (cross_population[pos][0].length() > cross_population[pos][1].length())
            {
                int pos1 = rand() % cross_population[pos][1].length();
                int pos2 = rand() % cross_population[pos + 1][0].length();
                cross_population[pos][1].insert(pos1, "-");
                cross_population[pos + 1][0].insert(pos2, "-");
            }

            pos += 2;
        }
        else
        {
            cross_population[pos][0] = rand_population[chromo1][0];
            cross_population[pos][1] = rand_population[chromo1][1];
            cross_population[pos + 1][0] = rand_population[chromo2][0];
            cross_population[pos + 1][1] = rand_population[chromo2][1];
            pos += 2;
        }
    }
}

void GA::Mutation()
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < ROW; j++)
            mutation_population[i][j] = cross_population[i][j];
    }
    int point1, point2;
    for (int i = 0; i < N; i++)
    {
        double r = (double)rand() / RAND_MAX;
        if (Mutation_rate >= r)
        {
            // int x = rand() % 2;
            // if (x % 2 == 0)
            // {
            //     point1 = mutation_population[i][0].find("-");
            //     point2 = mutation_population[i][1].find("-");
            //     if (point1 != mutation_population[i][0].npos)
            //         mutation_population[i][0].erase(point1, 1);
            //     if (point2 != mutation_population[i][1].npos)
            //         mutation_population[i][1].erase(point2, 1);
            // }
            // else
            // {
            point1 = rand() % mutation_population[i][0].length();
            point2 = rand() % mutation_population[i][1].length();
            mutation_population[i][0].insert(point1, "-");
            mutation_population[i][1].insert(point2, "-");
            // }
        }
        for (int j = mutation_population[i][0].length() - 1; j >= 0; j--)
        {
            if (mutation_population[i][0][j] == '-' && mutation_population[i][1][j] == '-')
            {
                mutation_population[i][0].erase(j, 1);
                mutation_population[i][1].erase(j, 1);
            }
        }
    }
}

void GA::Cal_fitness_2()
{
    memset(F2, 0, sizeof(F2));
    int pos = 0;
    for (int i = 0; i < N; i++)
    {
        for (pos = 0; pos < mutation_population[i][0].length(); pos++)
        {
            F2[i] += cost(mutation_population[i][0][pos], mutation_population[i][1][pos]);
        }
    }
}

void GA::Best_Solution()
{
    string all_population[2 * N][ROW];
    double all_fitness[2 * N];
    for (int i = 0; i < 2 * N; i++)
    {
        if (i < N)
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
                all_population[i][j] = mutation_population[i - N][j];
            }
            all_fitness[i] = F2[i - N];
        }
    }
    int max_index = 0;
    min_cost[times] = 1000;
    for (int i = 0; i < 2 * N; i++)
    {
        if (all_fitness[i] <= min_cost[times])
        {
            min_cost[times] = all_fitness[i];
            max_index = i;
        }
    }
    best_alignment[times][0] = all_population[max_index][0];
    best_alignment[times][1] = all_population[max_index][1];
    // cout << "The " << times << "th best solution is\n"
    //      << all_population[max_index][0] << "\n"
    //      << all_population[max_index][1] << endl;
    // cout << "Cost is " << max_fitness << endl;
}

void GA::Mixing_population()
{
    int mut_num = 0.8 * N;
    int init_num = 0.2 * N;
    // int new_num = 0.1*N;
    int mix_index = 0;
    double F1_copy[N];
    copy(F1, F1 + N, F1_copy);
    double F2_copy[N];
    copy(F2, F2 + N, F2_copy);
    double F3[N];

    double F1_sort[N];
    int F1_preindex[N];
    copy(F1, F1 + N, F1_sort);
    sort(F1_sort, F1_sort + N);
    // reverse(F1_sort,F1_sort + N);
    for (int i = 0; i < N; i++)
    {
        F1_preindex[i] = (find(F1_copy, F1_copy + N, F1_sort[i]) - F1_copy);
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

    double F2_sort[N];
    int F2_preindex[N];
    copy(F2, F2 + N, F2_sort);
    sort(F2_sort, F2_sort + N);
    // reverse(F2_sort,F2_sort+N);
    for (int i = 0; i < N; i++)
    {
        F2_preindex[i] = (find(F2_copy, F2_copy + N, F2_sort[i]) - F2_copy);
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
    for (int i = 0; i < N; i++)
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
    for (int i = 0; i < 5; i++)
    {
        test.Genetic_A(base, Query[i]);
    }
    end_time = clock();
    cout<<double(end_time-start_time)/CLOCKS_PER_SEC<<endl;
    // cout<<test.initial_population[0];
    return 0;
}