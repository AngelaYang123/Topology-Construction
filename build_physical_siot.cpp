#include<bits/stdc++.h>
using namespace std;

int node_num;
vector<int> *Virtual;
double **CT;
int **random;

class Node{
    public:
        double x = 0, y = 0;
}; Node *node;

void calculate_distance(){
    int i, j;
    double tmp;

    for(i=0; i<node_num; i++)
        for(j=i; j<node_num; j++){
            if(i == j){
                CT[i][j] = 0;
            }
            else{
                tmp = pow(node[i].x - node[j].x, 2) + pow(node[i].y-node[j].y, 2);
                CT[i][j] = sqrt(tmp);
                CT[j][i] = CT[i][j];
            }
        }

    return;       
}

int *parent;

void update_parent(int a, int b){
    int target, temp;

    if(a > b){
        temp = parent[a];
        target = a;
        while(temp != -1){
            target = temp;
            temp = parent[temp];
        }
        parent[target] = b;
    }
    else{
        temp = parent[b];
        target = b;
        while(temp != -1){
            target = temp;
            temp = parent[temp];
        }
        parent[target] = a;
    }

    return;
}

bool iscycle(int a, int b){
    int temp, temp1, temp2;

    temp = parent[a];
    while(temp != -1){
        if(temp == b)
            return true;
        temp2 = temp;
        temp = parent[temp];
    }

    temp1 = parent[b];
    while(temp1 != -1){
        if(temp1 == a)
            return true;
        if(temp2 == temp1)
            return true;
        temp1 = parent[temp1];
    }
    
    return false;
}

bool isneighbor(int a, int b){
    int i;

    for(i=0; i<Virtual[a].size(); i++)
        if(Virtual[a][i] == b)
            return true;

    for(i=0; i<Virtual[b].size(); i++)
        if(Virtual[b][i] == a)
            return true;

    return false;
}

bool isconnected(){
    int i, neighbor;
    queue<int> Q;
    bool visited[node_num] = {false};

    Q.push(0);
    visited[0] = true;

    while(!Q.empty()){
        for(i=0; i<Virtual[Q.front()].size(); i++){
            neighbor = Virtual[Q.front()][i];
            if(!visited[neighbor]){
                visited[neighbor] = true;
                Q.push(neighbor);
            }
        }
        Q.pop();
    }

    for(i=0; i<node_num; i++)
        if(!visited[i])
            return false;
    
    return true;
}

void add_edge(int a, int b){
    Virtual[a].push_back(b);
    Virtual[b].push_back(a);

    return;
}

void matrix_to_vector(){
    int i, j;

    for(i=0; i<node_num; i++)
        for(j=i+1; j<node_num; j++)
            if(random[i][j] == 1){
                Virtual[i].push_back(j);
                Virtual[j].push_back(i);
            }

    return;
}

void clean_vector(){
    int i;

    for(i=0; i<node_num; i++)
        while(!Virtual[i].empty())
            Virtual[i].pop_back();
    
    return;
}

void get_graph(){
    int i, j, neighbor;

    fstream output;
    output.open("SIoT_16_physical.txt", ios::out);

    for(i=0; i<node_num; i++){
        sort(Virtual[i].begin(), Virtual[i].end());
        for(j=0; j<Virtual[i].size(); j++){
            neighbor = Virtual[i][j];
            if(i < neighbor)
                output << i << "\t" << neighbor << "\t" << CT[i][neighbor] << endl;
        }
    }

    output.close();
    return;
}

void random_build(){
    int i, j, k;

    for(i=0; i<node_num; i++)
        for(j=i+1; j<node_num; j++){
            k = rand()%(1-0+1)+0;
            random[i][j] = k;
            random[j][i] = k;
        }
    
    matrix_to_vector();

    if(isconnected())
        get_graph();
    else
        clean_vector();
    
    return;
}

void build_graph(){
    int i, j;
    int left, right;
    double min;

    while(1){
        min = 10e20;
        for(i=0; i<node_num; i++)
            for(j=i; j<node_num; j++){
                if(i == j || isneighbor(i, j) || iscycle(i, j)) 
                    continue;

                if(min > CT[i][j]){
                    min = CT[i][j];
                    left = i;
                    right = j;
                }
            }

        update_parent(left, right);
        if(min == 10e20)
            break;
        add_edge(left, right);
    }

    return;
}

void initial(){
    int i, j;

    node = new Node[node_num];
    Virtual = new vector<int> [node_num];
    parent = new int [node_num];

    CT = new double*[node_num];
    random = new int*[node_num];
    for(i=0; i<node_num; i++){
        CT[i] = new double[node_num];
        random[i] = new int[node_num];
        parent[i] = -1;
        for(j=0; j<node_num; j++)
            random[i][j] = 0;
    }
    
    return;
}

void end(){
    int i;

    delete[] node;
    delete[] Virtual;
    delete[] parent;

    for(i=0; i<node_num; i++){
        delete[] CT[i];
        delete[] random[i];
    }
    
    delete[] CT;
    delete[] random;

    return;
}

int main(){
    int i;
    double a, b;

    fstream input;
    input.open("SIoT_16.txt", ios::in);

    input >> node_num;
    initial();

    for(i=0; i<node_num; i++){
        input >> a >> b;
        node[i].x = a * 5;
        node[i].y = b * 5;
    }

    input.close();

    calculate_distance();

    // for(i=0; i<node_num; i++){
    //     for(int j=0; j<node_num; j++)
    //         cout << CT[i][j] << " ";
    //     cout << endl;
    // }

    random_build();
    // build_graph();
    // get_graph();

    end();
    return 0;
}