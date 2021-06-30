#include<bits/stdc++.h>
using namespace std;

int node_num;
vector<int> *Virtual, *Physical, *Social;
double **metropolis, **before, **after, **HT;
int **CT;

class Node{
    public:
        bool DP = false;
}; Node *node;

/////////////// H ////////////////////
void calculate_metropolis_matrix(){
    int i, j, neighbor;
    double sum = 0;
    
    for(i=0; i<node_num; i++){
        for(j=0; j<Virtual[i].size(); j++){
            neighbor = Virtual[i][j];

            metropolis[i][neighbor] = (double)1 / (2 * max(Virtual[i].size(), Virtual[neighbor].size()));
            sum += metropolis[i][neighbor]; 

            before[i][neighbor] = metropolis[i][neighbor] * -1;
            after[i][neighbor] = before[i][neighbor];
        }

        metropolis[i][i] = (double)1 - sum;
        sum = 0;
        before[i][i] = (double)1 - metropolis[i][i];
                        
        before[i][node_num] = 1;
        after[i][i] = before[i][i];
        after[i][node_num] = 1;  
    }
}

void gauss_elimination(){
    int i, j, x = 0, y = 0;
    double zero = 0.00000001;  
    double temp;  
  
    while(x < node_num+1 && y < node_num){  
        while(x < node_num+1 && after[y][x] < zero && after[y][x] > (-1*zero)){  
            j = y+1;  
            while(j < node_num && after[y][x] < zero && after[y][x] > (-1*zero)) 
                j++;  
   
            if(j >= node_num){ 
                x++; 
                continue; 
            }  
  
            for(i=x; i<node_num+1; i++){  
                temp = after[j][i];  
                after[j][i] = after[y][i];  
                after[y][i] = temp;  
            }  
             
            break;  
        }  

        if(x >= node_num+1)
            break;  
   
        for(i=node_num+1-1; i>x; i--){  
            after[y][i] /= after[y][x];  
        }  
        after[y][x] = 1;  
  
        for(j=0; j<node_num; j++){  
            if(j == y)
                continue; 
   
            for(i=node_num+1-1; i>=x; i--)
                after[j][i] -= after[y][i]*after[j][x];   
        }  

        x++;y++;  
    }  
}

bool isreach(int a, int b){
    int i, neighbor;
    queue<int> Q;
    bool visited[node_num] = {false};

    Q.push(a);
    visited[a] = true;

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

    if(visited[b])
        return true;
    else
        return false;
}

void calculate_HT(){
    int i, j, k;
            
    calculate_metropolis_matrix();

    for(i=0; i<node_num; i++){

        for(j=0; j<node_num+1; j++){
            if(i == j)
                after[i][j] = 1;
            else
                after[i][j] = 0;
        }        

        gauss_elimination();

        for(j=0; j<node_num; j++)
            HT[j][i] = after[j][node_num];
        
        for(j=0; j<node_num; j++)
            for(k=0; k<node_num+1; k++)
                after[j][k] = before[j][k];

    } 
    
    for(i=0; i<node_num; i++)
        for(j=0; j<node_num; j++)
            if(i != j && !isreach(i, j))
                HT[i][j] = 10e10;
}
//////////////////////////////////////

void add_edge(int a, int b){
    Virtual[a].push_back(b);
    Virtual[b].push_back(a);
}

int add_DP(int DP_num, int a, int b){

    if(node[a].DP == false){
        node[a].DP = true;
        DP_num++;
    }

    if(node[b].DP == false){
        node[b].DP = true;
        DP_num++;
    }

    return DP_num;
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

bool isfriend(int a, int b){
    int i;

    for(i=0; i<Social[a].size(); i++)
        if(Social[a][i] == b)
            return true;

    for(i=0; i<Social[b].size(); i++)
        if(Social[b][i] == a)
            return true;

    return false;
}

void get_graph(int edge_num, int constraint, double min_HTCT, double min_CT, double min_comp){
    int i, j;

    fstream output;
    string filename;

    filename = ".//32_node_for_regression//OURS_" + to_string(constraint) + "_" + to_string(min_HTCT) + "_" + to_string(min_CT-min_comp) + "_" + to_string(min_comp) + ".txt";
    output.open(filename, ios::out);

    output << node_num << endl;

    for(i=0; i<node_num; i++){
        sort(Virtual[i].begin(), Virtual[i].end());
        for(j=0; j<Virtual[i].size(); j++){
            if(i < Virtual[i][j] && isfriend(i, Virtual[i][j]))
                output << i << "\t" << Virtual[i][j] << "\t" << 0 << endl;
            else if(i < Virtual[i][j] && !isfriend(i, Virtual[i][j]))
                output << i << "\t" << Virtual[i][j] << "\t" << 1 << endl;
        }
    }

    output.close();
}      

void process(int constraint){
    int i, j;
    int left, right;

    int edge_num = 0;
    int DP_num = 0, DP_edge_num = 0, check = 0;

    int soc_degree, vir_degree;

    double max1 = -1, max2 = -1, max3 = -1, min = 10e20;
    double min_HTCT = 10e20, min_HT, min_CT;
    int min_edge;
    double score;

    fstream output;

    while(1){

        calculate_HT();

        max1 = -1; max2 = -1;

        for(i=0; i<node_num; i++){
            for(j=0; j<node_num; j++)
                if(i != j && max1 < HT[i][j])
                    max1 = HT[i][j];

            for(j=0; j<Virtual[i].size(); j++)
                if(max2 < CT[i][Virtual[i][j]]){
                    max2 = CT[i][Virtual[i][j]];
                }
        }

        double objective = (0.3165 * max1 + 121.6365 * (double)DP_num/node_num + 18.6251) * max2; //max1 * max2; 
        cout << objective << endl;

        if(isconnected()){
            if(min_HTCT > objective){
                min_HTCT = objective;
                min_HT = max1;
                min_CT = max2;
                min_edge = edge_num;
                // if(min_edge == 64)
                //     get_graph(min_edge, constraint, min_HTCT, min_CT);
            }
        }
        else{
            min = 10e20;
            for(i=0; i<node_num; i++)
                for(j=0; j<node_num; j++){
                    check = 0;
                    score = 0;

                    if(i == j || isneighbor(i, j)) 
                        continue;

                    if(!isfriend(i, j))
                        score = 10000000;

                    double formula = CT[i][j] + ((Virtual[i].size() + Virtual[j].size()) * 10000000000) + score;
                    // cout << "*" << i << " " << j << " " << formula << endl;
                    if(min > formula && DP_num <= constraint){
                        
                        if(!isfriend(i, j) && node[i].DP == false)
                            check++;
                        if(!isfriend(i, j) && node[j].DP == false)   
                            check++;
                        if(DP_num + check > constraint)
                            continue; 

                        min = formula;
                        left = i;
                        right = j;  
                    }
                }

            goto line;              
        }

        max1 = -1; 

        if(DP_num != constraint){
            for(i=0; i<node_num; i++)
                for(j=0; j<node_num; j++){
                    check = 0;
                    soc_degree = Social[i].size() * Social[j].size();

                    if(isneighbor(i, j) || i == j)
                        continue;

                    if(max1 < HT[i][j]/CT[i][j]/soc_degree && DP_num < constraint){
                    
                        if(!isfriend(i, j) && node[i].DP == false)
                            check++;
                        if(!isfriend(i, j) && node[j].DP == false)   
                            check++;
                        if(DP_num + check > constraint)
                            continue;                  

                        max1 = HT[i][j]/CT[i][j]/soc_degree;
                        left = i;
                        right = j;
                    }
                }
        }

        else{
            for(i=0; i<node_num; i++)
                for(j=0; j<node_num; j++){
                    vir_degree = Virtual[i].size() * Virtual[j].size();

                    if(isneighbor(i, j) || i == j)
                        continue;

                    if(max1 < HT[i][j]/CT[i][j]/vir_degree && (isfriend(i, j) || !isfriend(i, j) && node[i].DP == true && node[j].DP == true)){
                        max1 = HT[i][j]/CT[i][j]/vir_degree;
                        left = i;
                        right = j;
                    }
                }
        }

        if(max1 == -1)
            break;

        line:
        if(min == 10e20)
            break;

        add_edge(left, right);
        edge_num++;
        
        if(!isfriend(left, right)){
            DP_num = add_DP(DP_num, left, right);
        }

        DP_edge_num = 0;
        for(i=0; i<node_num; i++)
            if(node[i].DP == true)
                DP_edge_num += Virtual[i].size();

        cout << edge_num << "  " << left << " " << right << "  " << DP_num << "  ";
    }

    cout << endl;
    cout << min_edge << " " << min_HT << " * " << min_CT << " = " << min_HTCT << endl; 
}

void initial(){
    int i, j;
    
    node = new Node[node_num];

    Virtual = new vector<int> [node_num];
    Physical = new vector<int> [node_num];
    Social = new vector<int> [node_num];

    metropolis = new double*[node_num];
    before = new double*[node_num];
    after = new double*[node_num];
    HT = new double*[node_num];

    CT = new int*[node_num];

    for(i=0; i<node_num; i++){
        metropolis[i] = new double[node_num];
        before[i] = new double[node_num+1];
        after[i] = new double[node_num+1];
        HT[i] = new double[node_num];

        CT[i] = new int[node_num];
    }
   
    for(i=0; i<node_num; i++){
        for(j=0; j<node_num; j++){
            if(i == j){
                metropolis[i][j] = 1;
                before[i][j] = 0;
                after[i][j] = 0;
                HT[i][j] = 0;

                CT[i][j] = 0;
            }
            else{
                metropolis[i][j] = 0;
                before[i][j] = 0;
                after[i][j] = 0;
                HT[i][j] = -1;

                CT[i][j] = -1;
            }
        }
    }
}

void end(){
    int i;

    delete[] node;

    delete[] Virtual;
    delete[] Physical;
    delete[] Social;
    
    for(i=0; i<node_num; i++){
        delete[] metropolis[i];
        delete[] before[i];
        delete[] after[i];
        delete[] HT[i];

        delete[] CT[i];
    }
    delete[] metropolis;
    delete[] before;
    delete[] after;
    delete[] HT;
    
    delete[] CT;
}

int main(int argc, char *argv[]){
    int i, j;
    int socLink, phyLink;
    int constraint;
    int left, right;
    int weight;

    fstream soc, phy;

    soc.open(argv[1], ios::in);
    phy.open(argv[2], ios::in);

    soc >> node_num >> socLink;
    phy >> node_num >> phyLink;
    
    constraint = atoi(argv[3]); 
    
    initial();

    for(i=0; i<socLink; i++){
        soc >> left >> right;
        Social[left].push_back(right);
        Social[right].push_back(left);
    }

    for(i=0; i<phyLink; i++){
        phy >> left >> right >> weight;
        CT[left][right] = weight;
        CT[right][left] = weight;
    }        

    soc.close();
    phy.close();

    process(constraint);
    end();
    
    return 0;
}