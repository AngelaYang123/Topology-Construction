#include<bits/stdc++.h>
using namespace std;

#include "Eigen\Eigen"
using namespace Eigen;

int node_num = 16;
int regular = 4;
int link_num = node_num * regular / 2;
int graph_num = 10000, graph_count = 0;
vector<int> *Virtual;
double **metropolis, **before, **after, **HT;
int **M;

fstream result1;
fstream result2;

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

void matrix_to_vector(){
    int i, j;

    for(i=0; i<node_num; i++)
        for(j=0; j<node_num; j++)
            if(M[i][j] == 1 && i < j){
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

double calculate_var(){
    int i, j;
    double var;
    
    for(i=0; i<node_num; i++){
        j = Virtual[i].size() - regular;
        var += abs(j);
    }
    var /= node_num;
    return var;
}

void get_HT_and_SG(){
    int i, j, k;
    double max1 = -1, max2 = -1, var;

    fstream output;
    string filename;

    filename = ".//graph//graph_" + to_string(graph_count) + ".txt";
    output.open(filename, ios::out);

    calculate_HT();
    Matrix<double, 16, 16> A;

    for(i=0; i<node_num; i++){
        for(j=0; j<node_num; j++){
            cout << HT[i][j] << " ";
            A(i, j) = metropolis[i][j];
            if(i != j && max1 < HT[i][j])
                max1 = HT[i][j];
        }
        cout << endl;
    }

    EigenSolver<Matrix<double, 16, 16>> es(A);
    Matrix<double, 16, 16> D = es.pseudoEigenvalueMatrix();

    for(i=0; i<node_num; i++)
        if(A(i, i) != 1 && max2 < A(i, i))
            max2 = A(i, i);

    var = calculate_var();

    output << "HT = " << max1 << endl;
    output << "SG = " << max2 << endl;
    output << "var = " << var << endl; 

    output << node_num << "\t" << link_num << endl;

    for(i=0; i<node_num; i++)
        for(j=0; j<Virtual[i].size(); j++){
            k = Virtual[i][j];
            if(i < k)
                output << i << "\t" << k << endl;
        }

    result1 << var << "\t" << max1 << endl;
    result2 << var << "\t" << max2 << endl;
    output.close();
    return;
}

void initial(){
    int i, j;
   
    for(i=0; i<node_num; i++){
        for(j=0; j<node_num; j++){
            if(i == j){
                metropolis[i][j] = 1;
                before[i][j] = 0;
                after[i][j] = 0;
                HT[i][j] = 0;
                M[i][j] = 0;
            }
            else{
                metropolis[i][j] = 0;
                before[i][j] = 0;
                after[i][j] = 0;
                HT[i][j] = -1;
                M[i][j] = 0;
            }
        }
    }
}

void end(){
    int i;

    delete[] Virtual;
    
    for(i=0; i<node_num; i++){
        delete[] metropolis[i];
        delete[] before[i];
        delete[] after[i];
        delete[] HT[i];
        delete[] M[i];
    }
    delete[] metropolis;
    delete[] before;
    delete[] after;
    delete[] HT;
    delete[] M;
}

int main(int argc, char *argv[]){
    int i, j;
    int r, link_count = 0;

    Virtual = new vector<int> [node_num];

    metropolis = new double*[node_num];
    before = new double*[node_num];
    after = new double*[node_num];
    HT = new double*[node_num];
    M = new int*[node_num];

    for(i=0; i<node_num; i++){
        metropolis[i] = new double[node_num];
        before[i] = new double[node_num+1];
        after[i] = new double[node_num+1];
        HT[i] = new double[node_num];
        M[i] = new int[node_num];
    }
    
    result1.open(".//graph//result1.txt", ios::out);
    result2.open(".//graph//result2.txt", ios::out);

    result1 << "var" << "\t" << "HT" << endl;
    result2 << "var" << "\t" << "SG" << endl;

    while(graph_count < graph_num){

        link_count = 0;
        initial();
        clean_vector(); 

        for(i=0; i<node_num-1; i++){
            for(j=i+1; j<node_num; j++){
                r = rand()%(1-0+1)+0;
                if(r == 1){
                    link_count++;
                    cout << link_count << endl;
                    M[i][j] = r;
                    M[j][i] = r;
                }
                if(link_count >= link_num)
                    break;
            }
            if(link_count >= link_num)
                break;
        }

        if(link_count >= link_num){
            matrix_to_vector();
            if(isconnected() == 1){
                get_HT_and_SG();
                graph_count++;               
            }
            else
                continue;
        }
        else
            continue;

    }

    end();
    result1.close();
    result2.close();
    return 0;
}