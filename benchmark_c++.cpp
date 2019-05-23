#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
using namespace std;


const int n = 5000;
vector< vector< float> > base, truth, sample, partials, B;
float margin_factor = 2 * sqrt(log (n));

void start(){
    for(int i=0; i<n; i++){
        vector<float> base_part, truth_part, sample_part, B_part, partials_part;
        for(int j=0; j<n; j++){
            base_part.push_back(0);
            truth_part.push_back(0);
            sample_part.push_back(0);
            B_part.push_back(0);
            partials_part.push_back(0);
        }
        partials_part.push_back(0);
        B_part.push_back(0);
        base.push_back(base_part);
        truth.push_back(truth_part);
        sample.push_back(sample_part);
        partials.push_back(partials_part);
        B.push_back(B_part);
    }
}

void print(vector< vector<float > > A){
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            cout<<A[i][j]<<" ";
        }
        cout<<endl;
    }
}

//Initial matrix
void build_matrix(){
    for(int i=0; i<n ; i++){
        for(int j=0; j<n; j++){
            float level = floor((4*(i+j))/(2*n));
            base[i][j] = (2*level + 1)/8;
        }
    }

}

//Randomly shuffle and add noise
void initialize(){
    vector<int> perm1, perm2;
    for(int i=0; i<n; i++){
        perm1.push_back(i);
        perm2.push_back(i);
    }
    random_shuffle(perm1.begin(),perm1.end());
    random_shuffle(perm2.begin(), perm2.end());


    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            truth[perm1[i]][perm2[j]] = base[i][j];
            float test = ((float) rand()/RAND_MAX);
            if(test  < base[i][j])
                sample[perm1[i]][perm2[j]] = 1;
            else
                sample[perm1[i]][perm2[j]] = 0;
        }
    }    
}


//Function for drawing fresh sample
void resample(){
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            float test = ((float) rand()/RAND_MAX);
            if(test < truth[i][j])
                sample[i][j] = 1;
            else
                sample[i][j] = 0;
        }
    }
}

//Compute Frobenius Squared Error
float sumsq(vector< vector<float > > A, vector< vector<float > > B){
    float sum = 0;
    for(int i=0; i<n; i++){
        for(int j=0; j< n; j++){
            sum += (A[i][j] - B[i][j]) * (A[i][j] - B[i][j]);
        }
    }
    return sum;
}

//Permute rows of the matrices
void perm_rows(vector<int> perm){
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            B[i][j] = truth[i][j];
        }
    }
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            truth[i][j] = B[perm[i]][j];
        }
    }

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            B[i][j] = sample[i][j];
        }
    }
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            sample[i][j] = B[perm[i]][j];
        }
    }


    for(int i=0; i<n; i++){
        for(int j=0; j<n + 1; j++){
            B[i][j] = partials[i][j];
        }
    }
    for(int i=0; i<n; i++){
        for(int j=0; j<n + 1; j++){
            partials[i][j] = B[perm[i]][j];
        }
    }
}


//Sort rows by sum
vector<int> sum_sort(){
    vector<pair<float, int> > sums;
    for(int i=0; i<n ;i++){
        float sum = 0;
        for(int j=0; j<n; j++){
            sum += sample[i][j];
        }
        sums.push_back(make_pair(sum, i));
    }

    sort(sums.begin(), sums.end());
    vector<int> perm;
    for(int i=0; i<n; i++){
       perm.push_back(sums[i].second);
    }
    return perm;
}

//Compute partials sums to make later comparisons more efficient
void compute_partials(){
    for(int i=0; i<n; i++){
        partials[i][0] = 0;
    }
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            partials[i][j + 1] = partials[i][j] + sample[i][j];
        }
    }
}




//Compare rows based on partial sums on blocks of size ~0.5 sqrt(n) log n
int sum_comp(int a, int b){
    int s =  2 * floor(sqrt(n / log(n)));
    int left = 0;
    int right = n - 1;
    if(partials[a][right + 1] - partials[a][left] > margin_factor * sqrt(right + 1 - left) + partials[b][right + 1] - partials[b][left])
        return 1;
    else if(partials[b][right + 1] - partials[b][left] > margin_factor * sqrt(right + 1 - left) + partials[a][right + 1] - partials[a][left])
        return 2;
    for(int i=0; i<s; i++){
        left = floor((float) n * i / s);
        right = floor((float) n * (i + 1) / s) - 1;
        if(partials[a][right + 1] - partials[a][left] > margin_factor * sqrt(right + 1 - left) + partials[b][right + 1] - partials[b][left])
            return 1;
        else if(partials[b][right + 1] - partials[b][left] > margin_factor * sqrt(right + 1 - left) + partials[a][right + 1] - partials[a][left])
            return 2;
    }
    return 0;
}



//Transpose matrix
void transpose(){
    for(int i=0; i<n; i++){
        for(int j=i + 1; j<n; j++){
            float a = truth[i][j];
            truth[i][j] = truth[j][i];
            truth[j][i] = a;

            a = sample[i][j];
            sample[i][j] = sample[j][i];
            sample[j][i] = a;
        }
    }
}


//Graph structure for topological sorting
class Graph { 
    int V;
    list<int>* adj; 
    void topologicalSortUtil(int v, bool visited[], stack<int>& Stack); 
  
public: 
    Graph(int V); 
    void addEdge(int v, int w); 
    stack<int> topologicalSort(); 
}; 
  
Graph::Graph(int V) 
{ 
    this->V = V; 
    adj = new list<int>[V]; 
} 
  
void Graph::addEdge(int v, int w) 
{ 
    adj[v].push_back(w); 
} 
  
void Graph::topologicalSortUtil(int v, bool visited[], 
                                stack<int>& Stack) 
{ 
    visited[v] = true; 
    list<int>::iterator i; 
    for (i = adj[v].begin(); i != adj[v].end(); ++i) 
        if (!visited[*i]) 
            topologicalSortUtil(*i, visited, Stack); 
    Stack.push(v); 
} 
  
stack<int> Graph::topologicalSort() 
{ 
    stack<int> Stack; 
  
    bool* visited = new bool[V]; 
    for (int i = 0; i < V; i++) 
        visited[i] = false; 
    for (int i = 0; i < V; i++) 
        if (visited[i] == false) 
            topologicalSortUtil(i, visited, Stack); 
    return Stack;
} 



//Construct graph and sort rows
vector<int> row_sort(){
    compute_partials();
    Graph g(n);
    for(int i=0; i<n; i++){
        for(int j = i+1; j < n; j++){
            if(sum_comp(i,j) == 1)
                g.addEdge(j,i);
            else if(sum_comp(i,j) == 2)
                g.addEdge(i,j);
        }
    }
    stack<int> verts = g.topologicalSort();
    vector<int> perm;
    for(int i=0; i<n; i++){
        perm.push_back(verts.top());
        verts.pop();
    }

    return perm;
}

//Computes inverse of a permutation
vector<int> inverse(vector<int> perm){
    int inv[n];
    for(int i=0; i<n; i++){
        inv[perm[i]] = i;
    }
    vector<int> inverse;
    for(int i=0; i<n; i++){
        inverse.push_back(inv[i]);
    }
    return inverse;
}



int main() 
{
    start();
    cout<<"constructed"<<endl;
    srand(time(NULL));
    build_matrix();
    initialize();

    //Note to permute the columns we transpose the matrix, permute the rows, and then transpose back
    vector<int> pre, inv_pre, col_perm, row_perm;
    pre= sum_sort();
    inv_pre = inverse(pre);
    perm_rows(pre);
    transpose();
    resample();
    col_perm = row_sort();
    transpose();
    perm_rows(inv_pre);


    transpose();
    resample();
    pre = sum_sort();
    inv_pre = inverse(pre);
    perm_rows(pre);
    transpose();
    resample();
    row_perm = row_sort();
    transpose();
    perm_rows(inv_pre);
    transpose();


    perm_rows(row_perm);
    transpose();
    perm_rows(col_perm);
    transpose();


    
    //Write matrices to output file
    ofstream myfile;
    myfile.open (to_string(n) + "benchmark.csv");
    string s;
    for(int i=0; i<n; i++){
        s = "";
        for(int j=0; j<n; j++){
            s += (to_string(round(sample[i][j])) + ",");
        }
        s += "\n";
        myfile << s;
    }
    myfile.close();

    ofstream myfile2;
    myfile2.open (to_string(n) + "benchmark_truth.csv");
    for(int i=0; i<n; i++){
        s = "";
        for(int j=0; j<n; j++){
            s += (to_string(round(8 * truth[i][j])) + ",");
        }
        s += "\n";
        myfile2 << s;
    }
    myfile2.close();

    cout<<sumsq(truth, base)<<endl;
    


    return 0;
}
