#include<bits/stdc++.h>
using namespace std;

int main(){
    int i, count = 0, data = 10000;
    queue<double> var;
    double a, b, first, sum = 0;

    fstream result1, output1;
    fstream result2, output2;
    result1.open("result1.txt", ios::in);
    result2.open("result2.txt", ios::in);
    output1.open(".//graph//average1.txt", ios::out);
    output2.open(".//graph//average2.txt", ios::out);

    result1 >> a >> b;
    var.push(a);
    sum = b;
    count++;

    for(i=1; i<data; i++){
        result1 >> a >> b;
        first = var.front();
        if(a == first){
            sum += b;
            count++;
            if(i == data-1){
                output1 << first << "\t" << sum / count << endl;
                var.pop();
            }
        }
        else{
            output1 << first << "\t" << sum / count << endl;
            var.pop();
            var.push(a);
            sum = b;
            count = 1;
            if(i == data-1){
                output1 << a << "\t" << b << endl;
                var.pop();
            }
        }
    }
    
    count = 0;
    result2 >> a >> b;
    var.push(a);
    sum = b;
    count++;

    for(i=1; i<data; i++){
        result2 >> a >> b;
        first = var.front();
        if(a == first){
            sum += b;
            count++;
            if(i == data-1){
                output2 << first << "\t" << sum / count << endl;
                var.pop();
            }
        }
        else{
            output2 << first << "\t" << sum / count << endl;
            var.pop();
            var.push(a);
            sum = b;
            count = 1;
            if(i == data-1){
                output2 << a << "\t" << b << endl;
                var.pop();
            }
        }
    }

    result1.close();
    result2.close();
    output1.close();
    output2.close();
    return 0;
}