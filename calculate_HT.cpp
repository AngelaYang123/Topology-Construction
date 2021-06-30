#include<bits/stdc++.h>
using namespace std;

#include "Eigen\Eigen"
using namespace Eigen;

int node_num;
vector<int> *Virtual;
double **metropolis, **before, **after, **HT;

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

void initial(){
    int i, j;

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
}

void end(){
    int i;

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
}

int main(int argc, char *argv[]){
    int i, j;
    int link_num;
    int left, right;

    fstream input;

    input.open(argv[1], ios::in);
    
    input >> node_num >> link_num;
        
    initial();

    for(i=0; i<link_num; i++){
        input >> left >> right;
        Virtual[left].push_back(right);
        Virtual[right].push_back(left);
    }

    input.close();
    
    double max1 = -1, max2 = -1;

    calculate_HT();
    Matrix<double, 10, 10> A;

    for(i=0; i<node_num; i++){
        for(j=0; j<node_num; j++){
            cout << HT[i][j] << " ";
            A(i, j) = metropolis[i][j];
            if(i != j && max1 < HT[i][j])
                max1 = HT[i][j];
        }
        cout << endl;
    }

    EigenSolver<Matrix<double, 10, 10>> es(A);
    Matrix<double, 10, 10> D = es.pseudoEigenvalueMatrix();

    for(i=0; i<node_num; i++)
        if(A(i, i) != 1 && max2 < A(i, i))
            max2 = A(i, i);

    cout << "HT = " << max1 << endl;
    cout << "SG = " << max2 << endl;

    end();
    return 0;
}