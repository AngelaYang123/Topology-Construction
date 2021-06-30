#include<bits/stdc++.h>
using namespace std;

int node_num, flag = 0, *total_edge;
double *obj;
vector<int> *Virtual, *Physical, *Social;
double **metropolis, **before, **after, **HT;
double **CT;
string mode;

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

    return;
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

    return; 
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
    
    return;
}
//////////////////////////////////////

void CT_adjustment(int seed, int computation){
    default_random_engine gen(seed);
    uniform_real_distribution<double> dis(0, 1);
    for(int i = 0; i < node_num; i++)
    {
        double comp;
        double rand = dis(gen);
        double interval[5];
        
        switch (computation)
        {
            case 1:
                interval[0] = 0;
                interval[1] = 0;
                interval[2] = 0;
                interval[3] = 0;
                interval[4] = 1;
                break;
            case 2:
                interval[0] = 0;
                interval[1] = 0;
                interval[2] = 0;
                interval[3] = 1;
                interval[4] = 0;
                break;
            case 3:
                interval[0] = 0;
                interval[1] = 0;
                interval[2] = 1;
                interval[3] = 0;
                interval[4] = 0;
                break;
            case 4:
                interval[0] = 0;
                interval[1] = 1;
                interval[2] = 0;
                interval[3] = 0;
                interval[4] = 0;
                break;
            case 5:
                interval[0] = 1;
                interval[1] = 0;
                interval[2] = 0;
                interval[3] = 0;
                interval[4] = 0;
                break;
        }
        if(rand < interval[0])
            comp = 2.532212885;
        else if(rand < interval[1])
            comp = 6.746268657;
        else if(rand < interval[2])
            comp = 11.16049383;
        else if(rand < interval[3])
            comp = 14.125;
        else if(rand < interval[4])
            comp = 21.02325581;
        
        comp *= 1000;
        for(int j = 0; j < node_num; j++)
            if(i != j){
                CT[i][j] += comp;
            }
    }

    return;
}

void add_edge(int a, int b){
    Virtual[a].push_back(b);
    Virtual[b].push_back(a);

    return;
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

void get_graph(int constraint){
    int i, j;

    fstream output;
    string filename;

    filename = to_string(constraint) + "_"  + mode + ".txt";
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
    return;
}      

void process1(int constraint, int seed, int density, int computation){
    int i, j;
    int left, right;

    int edge_num = 0;
    int DP_num = 0, DP_edge_num = 0, check = 0;

    double max1 = -1, max2 = -1, min = 10e20;
    double min_HTCT = 10e20, min_HT, min_CT;
    int min_edge;
    double score;

    fstream output;
    string filename;

    while(1){
        
        calculate_HT();

        max1 = -1; max2 = -1;

        for(i=0; i<node_num; i++){
            for(j=0; j<node_num; j++)
                if(i != j && max1 < HT[i][j])
                    max1 = HT[i][j];

            for(j=0; j<Virtual[i].size(); j++)
                if(max2 < CT[i][Virtual[i][j]])
                    max2 = CT[i][Virtual[i][j]];
        }
        
        double objective = (0.3165 * max1 + 121.6365 * (double)DP_num/node_num + 18.6251) * max2;

        if(isconnected()){
            if(min_HTCT > objective){
                min_HTCT = objective;
                min_HT = max1;
                min_CT = max2;
                min_edge = edge_num;
                if(flag == 1 && min_edge == total_edge[constraint]){
                    get_graph(constraint);
                    return;
                }
            }
        }
        else{
            min = 10e20;
            for(i=0; i<node_num; i++)
                for(j=0; j<node_num; j++){
                    check = 0;

                    if(i == j || isneighbor(i, j) || iscycle(i, j)) 
                        continue;

                    if(min > CT[i][j] && DP_num <= constraint){
                        
                        if(!isfriend(i, j) && node[i].DP == false)
                            check++;
                        if(!isfriend(i, j) && node[j].DP == false)   
                            check++;
                        if(DP_num + check > constraint)
                            continue; 

                        min = CT[i][j];
                        left = i;
                        right = j;
                    }
                }

            update_parent(left, right);
            goto line;              
        }

        max1 = -1;

        for(i=0; i<node_num; i++)
            for(j=0; j<node_num; j++){
                check = 0;

                if(isneighbor(i, j) || i == j)
                    continue;

                if(max1 < HT[i][j] && DP_num <= constraint){
                    
                    if(!isfriend(i, j) && node[i].DP == false)
                        check++;
                    if(!isfriend(i, j) && node[j].DP == false)   
                        check++;
                    if(DP_num + check > constraint)
                        continue;                  

                    max1 = HT[i][j];
                    left = i;
                    right = j;
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
    }

    cout << constraint << " " << min_edge << " " << min_HT << " " << min_CT << " " << min_HTCT << endl; 
    obj[constraint] = min_HTCT;
    total_edge[constraint] = min_edge;
    filename = to_string(node_num) + "_node_regression_for_HOA.txt";
    output.open(filename, ios::out | ios::app);
    output << constraint << "\t" << min_edge << "\t" << min_HT << "\t" << min_CT << "\t" << min_HTCT << endl; 

    output.close();
    return;
}

void process2(int constraint, int seed, int density, int computation){
    int i, j;
    int left, right;

    int edge_num = 0;
    int DP_num = 0, DP_edge_num = 0, check = 0;

    double max1 = -1, max2 = -1, min = 10e20;
    double min_HTCT = 10e20, min_HT, min_CT;
    int min_edge;
    double score;

    fstream output;
    string filename;

    while(1){
        
        calculate_HT();

        max1 = -1; max2 = -1;

        for(i=0; i<node_num; i++){
            for(j=0; j<node_num; j++)
                if(i != j && max1 < HT[i][j])
                    max1 = HT[i][j];

            for(j=0; j<Virtual[i].size(); j++)
                if(max2 < CT[i][Virtual[i][j]])
                    max2 = CT[i][Virtual[i][j]];
        }
        
        double objective = (0.3165 * max1 + 121.6365 * (double)DP_num/node_num + 18.6251) * max2;

        if(isconnected()){
            if(min_HTCT > objective){
                min_HTCT = objective;
                min_HT = max1;
                min_CT = max2;
                min_edge = edge_num;
                if(flag == 1 && min_edge == total_edge[constraint]){
                    get_graph(constraint);
                    return;
                }
            }
        }
        else{
            min = 10e20;
            for(i=0; i<node_num; i++)
                for(j=0; j<node_num; j++){
                    check = 0;

                    if(i == j || isneighbor(i, j) || iscycle(i, j)) 
                        continue;

                    if(min > CT[i][j] && DP_num <= constraint){
                        
                        if(!isfriend(i, j) && node[i].DP == false)
                            check++;
                        if(!isfriend(i, j) && node[j].DP == false)   
                            check++;
                        if(DP_num + check > constraint)
                            continue; 

                        min = CT[i][j];
                        left = i;
                        right = j;
                    }
                }

            update_parent(left, right);
            goto line;              
        }

        max1 = -1;

        for(i=0; i<node_num; i++)
            for(j=0; j<node_num; j++){
                check = 0;

                if(isneighbor(i, j) || i == j)
                    continue;

                if(max1 < HT[i][j]/CT[i][j] && DP_num <= constraint){
                    
                    if(!isfriend(i, j) && node[i].DP == false)
                        check++;
                    if(!isfriend(i, j) && node[j].DP == false)   
                        check++;
                    if(DP_num + check > constraint)
                        continue;                  

                    max1 = HT[i][j]/CT[i][j];
                    left = i;
                    right = j;
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
    }

    cout << constraint << " " << min_edge << " " << min_HT << " " << min_CT << " " << min_HTCT << endl; 
    obj[constraint] = min_HTCT;
    total_edge[constraint] = min_edge;
    filename = to_string(node_num) + "_node_regression_for_GLS.txt";
    output.open(filename, ios::out | ios::app);
    output << constraint << "\t" << min_edge << "\t" << min_HT << "\t" << min_CT << "\t" << min_HTCT << endl; 

    output.close();
    return;
}

void process3(int constraint, int seed, int density, int computation){
    int i, j;
    int left, right;

    int edge_num = 0;
    int DP_num = 0, DP_edge_num = 0, check = 0;

    int soc_degree, vir_degree;

    double max1 = -1, max2 = -1, min = 10e20;
    double min_HTCT = 10e20, min_HT, min_CT;
    int min_edge;
    double score;

    fstream output;
    string filename;

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

        double objective = (0.3165 * max1 + 121.6365 * (double)DP_num/node_num + 18.6251) * max2;

        if(isconnected()){
            if(min_HTCT > objective){
                min_HTCT = objective;
                min_HT = max1;
                min_CT = max2;
                min_edge = edge_num;
                if(flag == 1 && min_edge == total_edge[constraint]){
                    get_graph(constraint);
                    return;
                }
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
    }

    cout << constraint << " " << min_edge << " " << min_HT << " " << min_CT << " " << min_HTCT << endl; 
    obj[constraint] = min_HTCT;
    total_edge[constraint] = min_edge;
    filename = to_string(node_num) + "_node_regression_for_AutoTag.txt";
    output.open(filename, ios::out | ios::app);
    output << constraint << "\t" << min_edge << "\t" << min_HT << "\t" << min_CT << "\t" << min_HTCT << endl; 

    output.close();
    return;
}

int find_suitable_DP_constraint(){
    int i, min_DP;
    double min = 10e20;

    fstream output;
    string filename;

    for(i=2; i<=node_num; i++)
        if(min > obj[i]){
            min = obj[i];
            min_DP = i;
        }
    
    filename = to_string(node_num) + "_node_regression_for_" + mode + ".txt";
    output.open(filename, ios::out | ios::app);
    output << endl;
    output << "suitable DP constraint = " << min_DP << endl; 

    output.close();
    return min_DP;
}

void initial(){
    int i, j;

    node = new Node[node_num];
    parent = new int [node_num];
    Virtual = new vector<int> [node_num];
    metropolis = new double*[node_num];
    before = new double*[node_num];
    after = new double*[node_num];
    HT = new double*[node_num];

    for(i=0; i<node_num; i++){
        metropolis[i] = new double[node_num];
        before[i] = new double[node_num+1];
        after[i] = new double[node_num+1];
        HT[i] = new double[node_num];
    }
   
    for(i=0; i<node_num; i++){
        parent[i] = -1;
        for(j=0; j<node_num; j++){
            if(i == j){
                metropolis[i][j] = 1;
                before[i][j] = 0;
                after[i][j] = 0;
                HT[i][j] = 0;
            }
            else{
                metropolis[i][j] = 0;
                before[i][j] = 0;
                after[i][j] = 0;
                HT[i][j] = -1;
            }
        }
    }

    return;
}

void end(){
    int i;

    delete[] node;
    delete[] parent;
    delete[] Virtual;
     
    for(i=0; i<node_num; i++){
        delete[] metropolis[i];
        delete[] before[i];
        delete[] after[i];
        delete[] HT[i];
    }
    delete[] metropolis;
    delete[] before;
    delete[] after;
    delete[] HT;

    return;
}

int main(int argc, char *argv[]){
    int i, j;
    int id, left, right;
    int socLink;
    int dp, constraint, density, seed, computation;
    double num;

    fstream soc, phy;
    node_num = atoi(argv[1]);
    string filename1 = argv[2];
    string filename2 = argv[3];
    mode = argv[4];
    
    Physical = new vector<int> [node_num];
    Social = new vector<int> [node_num];
    obj = new double[node_num+1];
    total_edge = new int[node_num+1];

    CT = new double*[node_num];  
    for(i=0; i<node_num; i++)
        CT[i] = new double[node_num];

    soc.open(filename1, ios::in); 
    phy.open(filename2, ios::in); 

    soc >> node_num >> socLink;
    phy >> node_num;
    
    density = socLink * 2 / node_num;
    seed = atoi(argv[5]);
    computation = atoi(argv[6]);

    for(i=0; i<socLink; i++){
        soc >> left >> right;
        Social[left].push_back(right);
        Social[right].push_back(left);
    }

    for(i=0; i<node_num; i++)
        for(j=0; j<node_num; j++){
            phy >> num;
            CT[i][j] = num;
        }

    CT_adjustment(seed, computation);

    soc.close();
    phy.close();

    for(dp=2; dp<=node_num; dp++){
        initial();
        if(mode == "HOA")
            process1(dp, seed, density, computation);
        else if(mode == "GLS")
            process2(dp, seed, density, computation);
        else
            process3(dp, seed, density, computation);
        end();   
    }

    constraint = find_suitable_DP_constraint();
    flag = 1;
    initial();
    if(mode == "HOA")
        process1(constraint, seed, density, computation);
    else if(mode == "GLS")
        process2(constraint, seed, density, computation);
    else
        process3(constraint, seed, density, computation);
    end(); 

    delete[] Physical;
    delete[] Social;
    delete[] obj;
    delete[] total_edge;

    for(i=0; i<node_num; i++)
        delete[] CT[i];
    delete[] CT;

    return 0;
}