#pragma GCC optimize(2)
#include "omp.h"
#include <iostream>
#include <fstream>
#include <queue>
#include <stdio.h>
#include <tuple>
#include <algorithm>
#include <time.h>
using namespace std;

clock_t start_time;
clock_t end_time;

struct Node
{
    int x;
    int y;
    int z;
    int f;
    bool operator<(const Node other) const
    {
        return f > other.f;
    }
    bool operator!=(const Node other) const
    {
        return (x != other.x || y != other.y || z != other.z);
    }
    Node() : x(), y(), z(), f() {}
    Node(int x, int y, int z) : x(x), y(y), z(z), f(0) {}
};

const int N = 200;
// vector<vector<vector<int>>> cost_table(N, vector<vector<int>>(N, vector<int>(N, 0)));
vector<vector<vector<int>>> cost_table(3, vector<vector<int>>(N, vector<int>(N, 0)));
vector<vector<vector<Node>>> came_from(N, vector<vector<Node>>(N, vector<Node>(N)));
vector<int> path;

class Graph
{
private:
    string seq1;
    string seq2;
    string seq3;

public:
    Graph(string seq1, string seq2, string seq3) : seq1(seq1), seq2(seq2), seq3(seq3) {}
    vector<Node> get_neighbours(Node current);
    int pair_cost(Node current, Node neighbour);
    int cost(Node current, Node neighbour);
};

vector<Node> Graph::get_neighbours(Node current)
{
    vector<Node> neighbours;
    // neighbours.reserve(256);
    if (current.x != seq1.length() - 1 && current.y != seq2.length() - 1 && current.z != seq3.length() - 1)
    {
        Node neighbour(current.x + 1, current.y + 1, current.z + 1);
        neighbours.push_back(neighbour);
    }
    if (current.x != seq1.length() - 1 && current.y != seq2.length() - 1)
    {
        Node neighbour(current.x + 1, current.y + 1, current.z);
        neighbours.push_back(neighbour);
    }
    if (current.y != seq2.length() - 1 && current.z != seq3.length() - 1)
    {
        Node neighbour(current.x, current.y + 1, current.z + 1);
        neighbours.push_back(neighbour);
    }
    if (current.x != seq1.length() - 1 && current.z != seq3.length() - 1)
    {
        Node neighbour(current.x + 1, current.y, current.z + 1);
        neighbours.push_back(neighbour);
    }
    if (current.x != seq1.length() - 1)
    {
        Node neighbour(current.x + 1, current.y, current.z);
        neighbours.push_back(neighbour);
    }
    if (current.y != seq2.length() - 1)
    {
        Node neighbour(current.x, current.y + 1, current.z);
        neighbours.push_back(neighbour);
    }
    if (current.z != seq3.length() - 1)
    {
        Node neighbour(current.x, current.y, current.z + 1);
        neighbours.push_back(neighbour);
    }
    return neighbours;
}

int Graph::pair_cost(Node current, Node neighbour)
{
    int cost[3];
    if (neighbour.x != current.x && neighbour.y != current.y)
    {
        if (seq1[neighbour.x] != seq2[neighbour.y])
            cost[0] = 3;
        else
            cost[0] = 0;
    }
    else
    {
        if (neighbour.x == current.x && neighbour.y == current.y)
            cost[0] = 0;
        else
            cost[0] = 2;
    }

    if (neighbour.x != current.x && neighbour.z != current.z)
    {
        if (seq1[neighbour.x] != seq3[neighbour.z])
            cost[1] = 3;
        else
            cost[1] = 0;
    }
    else
    {
        if (neighbour.x == current.x && neighbour.z == current.z)
            cost[1] = 0;
        else
            cost[1] = 2;
    }

    if (neighbour.y != current.y && neighbour.z != current.z)
    {
        if (seq3[neighbour.z] != seq2[neighbour.y])
            cost[2] = 3;
        else
            cost[2] = 0;
    }
    else
    {
        if (neighbour.y == current.y && neighbour.z == current.z)
            cost[2] = 0;
        else
            cost[2] = 2;
    }

    return (cost[0] + cost[1] + cost[2]);
}

int Graph::cost(Node current, Node neighbour)
{
    // int pair1_cost = pair_cost(current.x, current.y, neighbour.x, neighbour.y);
    // int pair2_cost = pair_cost(current.x, current.z, neighbour.x, neighbour.z);
    // int pair3_cost = pair_cost(current.y, current.z, neighbour.y, neighbour.z);
    // return (pair1_cost + pair2_cost + pair3_cost);
    return pair_cost(current, neighbour);
}

void build_cost_table(string seq1, string seq2, string seq3, int key)
{
    // seq1.insert(0, "-");
    // seq2.insert(0, "-");
    reverse(seq1.begin() + 1, seq1.end());
    reverse(seq2.begin() + 1, seq2.end());
    reverse(seq3.begin() + 1, seq3.end());
    vector<int> move(3);
    if (key == 1)
    {
        for (int i = 0; i < seq1.length(); i++)
        {
            for (int j = 0; j < seq2.length(); j++)
            {
                if (i == 0)
                {
                    cost_table[0][i][j] = 2 * j;
                    continue;
                }
                if (j == 0)
                {
                    cost_table[0][i][j] = 2 * i;
                    continue;
                }
                if (seq1[i] == seq2[j])
                    move[0] = cost_table[0][i - 1][j - 1];
                else
                    move[0] = cost_table[0][i - 1][j - 1] + 3;

                move[1] = cost_table[0][i - 1][j] + 2;
                move[2] = cost_table[0][i][j - 1] + 2;
                int min_cost = *min_element(move.begin(), move.end());
                cost_table[0][i][j] = min_cost;
            }
        }
    }
    if (key == 2)
    {
        for (int i = 0; i < seq1.length(); i++)
        {
            for (int j = 0; j < seq3.length(); j++)
            {
                if (i == 0)
                {
                    cost_table[1][i][j] = 2 * j;
                    continue;
                }
                if (j == 0)
                {
                    cost_table[1][i][j] = 2 * i;
                    continue;
                }
                if (seq1[i] == seq3[j])
                    move[0] = cost_table[1][i - 1][j - 1];
                else
                    move[0] = cost_table[1][i - 1][j - 1] + 3;

                move[1] = cost_table[1][i - 1][j] + 2;
                move[2] = cost_table[1][i][j - 1] + 2;
                int min_cost = *min_element(move.begin(), move.end());
                cost_table[1][i][j] = min_cost;
            }
        }
    }
    if (key == 3)
    {
        for (int i = 0; i < seq2.length(); i++)
        {
            for (int j = 0; j < seq3.length(); j++)
            {
                if (i == 0)
                {
                    cost_table[2][i][j] = 2 * j;
                    continue;
                }
                if (j == 0)
                {
                    cost_table[2][i][j] = 2 * i;
                    continue;
                }
                if (seq2[i] == seq3[j])
                    move[0] = cost_table[2][i - 1][j - 1];
                else
                    move[0] = cost_table[2][i - 1][j - 1] + 3;

                move[1] = cost_table[2][i - 1][j] + 2;
                move[2] = cost_table[2][i][j - 1] + 2;
                int min_cost = *min_element(move.begin(), move.end());
                cost_table[2][i][j] = min_cost;
            }
        }
    }
}

int Heuristic(Node current, Node dest)
{
    int sum = 0;
    sum += cost_table[0][dest.x - current.x][dest.y - current.y];
    sum += cost_table[1][dest.x - current.x][dest.z - current.z];
    sum += cost_table[2][dest.y - current.y][dest.z - current.z];
    return sum;
}

int Astar_search(string seq1, string seq2, string seq3, Node start, Node dest, int min_val)
{
    priority_queue<Node, vector<Node>> pq;
    vector<vector<vector<int>>> real_cost(seq1.length(), vector<vector<int>>(seq2.length(), vector<int>(seq3.length(), 1000)));
    Graph graph(seq1, seq2, seq3);
    pq.push(start);
    real_cost[start.x][start.y][start.z] = 0;
    came_from[start.x][start.y][start.z] = start;
    while (!pq.empty())
    {
        Node current = pq.top();
        pq.pop();
        // printf("%d \n", real_cost[current.x][current.y][current.z]);
        if (current.x == dest.x && current.y == dest.y && current.z == dest.z)
            break;

        for (auto neighbour : graph.get_neighbours(current))
        {
            int new_cost = real_cost[current.x][current.y][current.z] + graph.cost(current, neighbour);
            if (new_cost > min_val)
                break;
            if (real_cost[neighbour.x][neighbour.y][neighbour.z] == 1000 || real_cost[neighbour.x][neighbour.y][neighbour.z] > new_cost)
            {
                real_cost[neighbour.x][neighbour.y][neighbour.z] = new_cost;
                neighbour.f = new_cost + Heuristic(neighbour, dest);
                // if (neighbour.f > cost_table[dest.x][dest.y])
                //     continue;
                // else
                // {
                pq.push(neighbour);
                came_from[neighbour.x][neighbour.y][neighbour.z] = current;
                // }
                // pq.push(neighbour);
                // came_from[neighbour.x][neighbour.y] = current;
            }
        }
    }
    // while (!pq.empty())
    //     pq.pop();
    // int tmp = real_cost[dest.x][dest.y][dest.z];
    // vector<vector<int>>().swap(real_cost);
    return real_cost[dest.x][dest.y][dest.z];
}

void build_path(Node start, Node dest)
{
    path.clear();
    Node current = dest;
    while (came_from[current.x][current.y][current.z] != start)
    {
        Node father = came_from[current.x][current.y][current.z];
        if (current.x - father.x == 1 && current.y - father.y == 1 && current.z - father.z == 1)
            path.push_back(0);
        else if (current.x - father.x == 1 && current.y == father.y && current.z == father.z)
            path.push_back(1);
        else if (current.x == father.x && current.y - father.y == 1 && current.z == father.z)
            path.push_back(2);
        else if (current.x == father.x && current.y == father.y && current.z - father.z == 1)
            path.push_back(3);
        else if (current.x - father.x == 1 && current.y - father.y == 1 && current.z == father.z)
            path.push_back(4);
        else if (current.x - father.x == 1 && current.y == father.y && current.z - father.z == 1)
            path.push_back(5);
        else if (current.x == father.x && current.y - father.y == 1 && current.z - father.z == 1)
            path.push_back(6);
        current = father;
    }
}

void rebuild_align(string seq1, string seq2, string seq3)
{
    reverse(path.begin(), path.end());
    string new_seq1;
    string new_seq2;
    string new_seq3;
    seq1.erase(0, 1);
    seq2.erase(0, 1);
    seq3.erase(0, 1);
    int pos_1 = 0;
    int pos_2 = 0;
    int pos_3 = 0;
    for (int i = 0; i < path.size(); i++)
    {
        switch (path[i])
        {
        case 0:
            new_seq1.append(1, seq1[pos_1++]);
            new_seq2.append(1, seq2[pos_2++]);
            new_seq3.append(1, seq3[pos_3++]);
            break;

        case 1:
            new_seq1.append(1, seq1[pos_1++]);
            new_seq2.append(1, '-');
            new_seq3.append(1, '-');
            break;
        case 2:
            new_seq1.append(1, '-');
            new_seq2.append(1, seq2[pos_2++]);
            new_seq3.append(1, '-');
            break;
        case 3:
            new_seq1.append(1, '-');
            new_seq2.append(1, '-');
            new_seq3.append(1, seq3[pos_3++]);
            break;
        case 4:
            new_seq1.append(1, seq1[pos_1++]);
            new_seq2.append(1, seq2[pos_2++]);
            new_seq3.append(1, '-');
            break;
        case 5:
            new_seq1.append(1, seq1[pos_1++]);
            new_seq2.append(1, '-');
            new_seq3.append(1, seq3[pos_3++]);
            break;
        case 6:
            new_seq1.append(1, '-');
            new_seq2.append(1, seq2[pos_2++]);
            new_seq3.append(1, seq3[pos_3++]);
            break;
        default:
            break;
        }
    }
    cout << new_seq1 << endl;
    cout << new_seq2 << endl;
    cout << new_seq3 << endl;
}

int main()
{
    ///////////////////////////////////////////////////////////////
    // string seq1;
    // string seq2;
    // string seq3;
    // seq1 = "IPZJJLMLTKJULOSTKTJOGLKJOBLTXGKTPLUWWKOMOYJBGALJUKLGLOSVHWBPGWSLUKOBSOPLOOKUKSARPPJ";
    // seq2 = "IPZJJPLLTHUULOSTXTJOGLKJGBLLMMPJPLUWGKOMOYJBZAYUKOFLOSZHGBPHXPLXKJBXKJLAUUOJHWTWWPQ";
    // seq3 = "IPMJJLLLTHOULOSTMAJIGLKJPVLLXGKTPLTWWKOMOYJBZPYUKLILOSZHGBPGWXLZKJBSWJLPJUUMHKTRAP";
    // seq1.insert(0, "-");
    // seq2.insert(0, "-");
    // seq3.insert(0, "-");
    // Node start(0, 0, 0);
    // Node dest(seq1.length() - 1, seq2.length() - 1, seq3.length() - 1);
    // for (int i = 1; i <= 3; i++)
    //     build_cost_table(seq1, seq2, seq3, i);
    // int real_cost = Astar_search(seq1, seq2, seq3, start, dest, 1000);
    // cout << real_cost << endl;
    // build_path(start,dest);
    // rebuild_align(seq1,seq2,seq3);
    ////////////////////////////////////////////////////////////////
    start_time = clock();
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
    Node start(0, 0, 0);

    int min_cost[2] = {10000, 10000};
    for (int i = 0; i < 2; i++)
    {
        seq1 = Query[i + 6];
        seq1.insert(0, "-");
        int len1 = seq1.length();
        for (int j = 0; j < 100; j++)
        {
            seq2 = base[j];
            seq2.insert(0, "-");
            int len2 = seq2.length();
            build_cost_table(seq1, seq2, seq3, 0);
            if (cost_table[0][seq1.length() - 1][seq2.length() - 1] > (min_cost[i]) / 2)
                continue;
            for (int k = j + 1; k < 100; k++)
            {
                seq3 = base[k];
                seq3.insert(0, "-");
                int len3 = seq3.length();
                Node dest(len1 - 1, len2 - 1, len3 - 1);
                for (int key = 2; key <= 3; key++)
                    build_cost_table(seq1, seq2, seq3, key);
                if (cost_table[0][dest.x][dest.y] + cost_table[1][dest.x][dest.z] + cost_table[2][dest.y][dest.z] > min_cost[i])
                    continue;
                int real_cost = Astar_search(seq1, seq2, seq3, start, dest, min_cost[i]);
                // cout << k << "th done" << endl;
                // printf("%d \n",real_cost);
                if (real_cost < min_cost[i])
                {
                    min_cost[i] = real_cost;
                    build_path(start, dest);
                    best_align_1 = seq2;
                    best_align_2 = seq3;
                }
            }
        }
        rebuild_align(seq1, best_align_1, best_align_2);
        printf("Minimum cost = %d \n", min_cost[i]);
        printf("\n");
    }
    end_time = clock();
    cout << double(end_time - start_time) / CLOCKS_PER_SEC << endl;
    return 0;
}