#include <iostream>
#include <fstream>
#include <queue>
#include <stdio.h>
#include <tuple>
#include "omp.h"
#include <algorithm>
#include <time.h>
using namespace std;

clock_t start_time;
clock_t end_time;

struct Node
{
    int x;
    int y;
    int f;
    bool operator<(const Node other) const
    {
        return f > other.f;
    }
    bool operator!=(const Node other) const
    {
        return (x != other.x || y != other.y);
    }
    Node() : x(), y(), f() {}
    Node(int x, int y) : x(x), y(y), f() {}
};

const int N = 190;
vector<vector<int>> cost_table(N, vector<int>(N, 0));
vector<vector<Node>> came_from(N, vector<Node>(N));
vector<int> path;

class Graph
{
private:
    string seq1;
    string seq2;

public:
    Graph(string seq1, string seq2) : seq1(seq1), seq2(seq2) {}
    vector<Node> get_neighbours(Node current);
    int cost(Node current, Node neighbour);
};

vector<Node> Graph::get_neighbours(Node current)
{
    vector<Node> neighbours;
    neighbours.reserve(256);
    if (current.x != seq1.length() - 1 && current.y != seq2.length() - 1)
    {
        Node neighbour(current.x + 1, current.y + 1);
        neighbours.push_back(neighbour);
    }
    if (current.x != seq1.length() - 1)
    {
        Node neighbour(current.x + 1, current.y);
        neighbours.push_back(neighbour);
    }
    if (current.y != seq2.length() - 1)
    {
        Node neighbour(current.x, current.y + 1);
        neighbours.push_back(neighbour);
    }
    return neighbours;
}

int Graph::cost(Node current, Node neighbour)
{
    if (neighbour.x != current.x && neighbour.y != current.y)
    {
        if (seq1[neighbour.x] != seq2[neighbour.y])
            return 3;
        else
            return 0;
    }
    else
        return 2;
}

void build_cost_table(string seq1, string seq2)
{
    // seq1.insert(0, "-");
    // seq2.insert(0, "-");
    reverse(seq1.begin() + 1, seq1.end());
    reverse(seq2.begin() + 1, seq2.end());
    vector<int> move(3);
    for (int i = 0; i < seq1.length(); i++)
    {
        for (int j = 0; j < seq2.length(); j++)
        {
            if (i == 0)
            {
                cost_table[i][j] = 2 * j;
                continue;
            }
            if (j == 0)
            {
                cost_table[i][j] = 2 * i;
                continue;
            }
            if (seq1[i] == seq2[j])
                move[0] = cost_table[i - 1][j - 1];
            else
                move[0] = cost_table[i - 1][j - 1] + 3;

            move[1] = cost_table[i - 1][j] + 2;
            move[2] = cost_table[i][j - 1] + 2;
            int min_cost = *min_element(move.begin(), move.end());
            cost_table[i][j] = min_cost;
        }
    }
}

int Heuristic(Node current, Node dest)
{
    return cost_table[dest.x - current.x][dest.y - current.y];
}

int Astar_search(string seq1, string seq2, Node start, Node dest, int min_val)
{
    priority_queue<Node, vector<Node>> pq;
    vector<vector<int>> real_cost(seq1.length(), vector<int>(seq2.length(), 1000));
    Graph graph(seq1, seq2);
    pq.push(start);
    real_cost[start.x][start.y] = 0;
    came_from[start.x][start.y] = start;
    while (!pq.empty())
    {
        Node current = pq.top();
        pq.pop();
        if (current.x == dest.x && current.y == dest.y)
            break;

        for (auto neighbour : graph.get_neighbours(current))
        {
            int new_cost = real_cost[current.x][current.y] + graph.cost(current, neighbour);
            if (new_cost > min_val)
                break;
            if (real_cost[neighbour.x][neighbour.y] == 1000 || real_cost[neighbour.x][neighbour.y] > new_cost)
            {
                real_cost[neighbour.x][neighbour.y] = new_cost;
                neighbour.f = new_cost + Heuristic(neighbour, dest);
                if (neighbour.f > cost_table[dest.x][dest.y])
                    continue;
                else
                {
                    pq.push(neighbour);
                    came_from[neighbour.x][neighbour.y] = current;
                }
                // pq.push(neighbour);
                // came_from[neighbour.x][neighbour.y] = current;
            }
        }
    }
    // while (!pq.empty())
    //     pq.pop();
    int tmp = real_cost[dest.x][dest.y];
    // vector<vector<int>>().swap(real_cost);
    return tmp;
}

void build_path(Node start, Node dest)
{
    path.clear();
    Node current = dest;
    while (came_from[current.x][current.y] != start)
    {
        Node father = came_from[current.x][current.y];
        if (current.x - father.x == 1 && current.y - father.y == 1)
            path.push_back(0);
        else if (current.x - father.x == 1)
        {
            path.push_back(1);
        }
        else if (current.y - father.y == 1)
        {
            path.push_back(2);
        }
        current = father;
    }
}

void rebuild_align(string seq1, string seq2, vector<int> path)
{
    reverse(path.begin(), path.end());
    string new_seq1;
    string new_seq2;
    seq1.erase(0, 1);
    seq2.erase(0, 1);
    int pos_1 = 0;
    int pos_2 = 0;
    for (int i = 0; i < path.size(); i++)
    {
        switch (path[i])
        {
        case 0:
            new_seq1.append(1, seq1[pos_1++]);
            new_seq2.append(1, seq2[pos_2++]);
            // pos_1++;
            // pos_2++;
            break;

        case 1:
            new_seq1.append(1, seq1[pos_1++]);
            new_seq2.append(1, '-');
            // pos_1++;
            break;
        case 2:
            new_seq1.append(1, '-');
            new_seq2.append(1, seq2[pos_2++]);
            // pos_2++;
            break;

        default:
            break;
        }
    }
    cout << new_seq1 << endl;
    cout << new_seq2 << endl;
}

int main()
{
    // string seq1;
    // string seq2;
    // seq1 = "KJXXJAJKPXKJJXJKPXKJXXJAJKPXKJJXJKPXKJXXJAJKPXKJXXJAJKHXKJXXJAJKPXKJXXJAJKHXKJXX";
    // seq2 = "XHAPXJAJXXXJAJXDJAJXXXJAPXJAHXXXJAJXDJAJXXXJAJXXXJPPXJAJXXXHAPXJAJXXXJAJXXXJAJXXXJA";
    // seq1.insert(0, "-");
    // seq2.insert(0, "-");
    // Node start(0, 0);
    // Node dest(seq1.length() - 1, seq2.length()-1);
    // build_cost_table(seq1, seq2);
    // int real_cost = Astar_search(seq1, seq2, start, dest);
    // cout << real_cost << endl;
    // return 0;
    // omp_set_num_threads(8);
    start_time = clock();
    ifstream query("MSA_query.txt");
    ifstream data("MSA_database.txt");
    string tmp;
    getline(query, tmp);
    vector<string> base;
    vector<string> Query;
    for (int i = 0; i < 5; i++)
    {
        getline(query, tmp);
        Query.push_back(tmp);
    }
    for (int i = 0; i < 100; i++)
    {
        getline(data, tmp);
        base.push_back(tmp);
    }
    // for (int i = 0; i < 5; i++)
    // {
    //     cout<<Query[i]<<endl;
    // }
    string seq1;
    string seq2;
    string best_align;

    int min_cost[5] = {10000, 10000, 10000, 10000, 10000};

    for (int i = 0; i < 5; i++)
    {
        seq1 = Query[i];
        seq1.insert(0, "-");
        int len1 = seq1.length();
        // #pragma omp parallel for
        for (int j = 0; j < 100; j++)
        {
            seq2 = base[j];
            seq2.insert(0, "-");
            int len2 = seq2.length();
            Node start(0, 0);
            Node dest(len1 - 1, len2 - 1);
            build_cost_table(seq1, seq2);
            if (min_cost[i] < cost_table[dest.x][dest.y])
                continue;
            int real_cost = Astar_search(seq1, seq2, start, dest, min_cost[i]);
            // printf("%d \n",real_cost);
            min_cost[i] = real_cost;
            best_align = seq2;
            build_path(start, dest);
            // best_align = seq2;
        }
        printf("Minimum cost = %d \n", min_cost[i]);
        rebuild_align(seq1, best_align, path);
        printf("\n");
        // for (auto n : cost_table)
        // {
        //     n.clear();
        // }
    }
    end_time = clock();
    cout << double(end_time - start_time) / CLOCKS_PER_SEC << endl;
    return 0;
}