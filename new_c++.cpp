#include <bits/stdc++.h>
#include <iostream>
using namespace std;


const int n = 12000;
vector< vector< float> > base, truth, sample, partials, B;
int min_horizontal_interval = 50;
int min_vertical_interval = 50;
int top = 0;
int bottom = n-1;
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

void print(vector< vector<float> > A){
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

//Compute Frobenius Squared Error
float sumsq(vector< vector< float> > A, vector< vector< float> > B){
    float sum = 0;
    for(int i=0; i<n; i++){
        for(int j=0; j< n; j++){
            sum += (A[i][j] - B[i][j]) * (A[i][j] - B[i][j]);
        }
    }
    return sum;
}


//Permute rows of the matrices
void perm_rows(int perm[][2], int l){
    for(int i=0; i<l; i++){
        for(int j=0; j<n; j++){
            B[perm[i][0]][j] = truth[perm[i][0]][j];
        }
    }
    for(int i=0; i<l; i++){
        for(int j=0; j<n; j++){
            truth[perm[i][0]][j] = B[perm[i][1]][j];
        }
    }

    for(int i=0; i<l; i++){
        for(int j=0; j<n; j++){
            B[perm[i][0]][j] = sample[perm[i][0]][j];
        }
    }
    for(int i=0; i<l; i++){
        for(int j=0; j<n; j++){
            sample[perm[i][0]][j] = B[perm[i][1]][j];
        }
    }


    for(int i=0; i<l; i++){
        for(int j=0; j<n + 1; j++){
            B[perm[i][0]][j] = partials[perm[i][0]][j];
        }
    }
    for(int i=0; i<l; i++){
        for(int j=0; j<n + 1; j++){
            partials[perm[i][0]][j] = B[perm[i][1]][j];
        }
    }
}


//Pivoting step for fixed test interval
void row_step(int left, int right){
    if(bottom - top < min_vertical_interval || right - left < min_horizontal_interval)
        return;
    vector<pair<float, int> > to_sort;
    for(int i=top; i<= bottom; i++){
        to_sort.push_back(make_pair(partials[i][right + 1] - partials[i][left], i));
    }
    sort(to_sort.begin(), to_sort.end());
    int progress = 0;
    while(partials[to_sort[bottom - top - progress].second][right + 1] - partials[to_sort[bottom - top - progress].second][left] >
        partials[to_sort[progress].second][right + 1] - partials[to_sort[progress].second][left] + margin_factor * sqrt(right + 1 - left) ){
            progress ++;
       }
    if(progress > 0){
        int test_top = progress - 1;
        int test_bottom = bottom - top - progress + 1;


        int perm_front[test_top + 1], perm_middle[test_bottom - test_top - 1], perm_back[bottom - top + 1 - test_bottom];
        int perm[bottom - top + 1][2];
        for(int i=0; i<test_top + 1; i++){
            perm_front[i] = to_sort[i].second;
        }

        for(int i=test_top + 1; i < test_bottom; i++){
            perm_middle[i - test_top - 1] = to_sort[i].second;
        }
        for(int i= test_bottom; i < bottom - top + 1; i++){
            perm_back[i - test_bottom] = to_sort[i].second;
        }

        sort(perm_front, perm_front + sizeof(perm_front)/sizeof(perm_front[0]));
        sort(perm_middle, perm_middle + sizeof(perm_middle)/sizeof(perm_middle[0]));
        sort(perm_back, perm_back + sizeof(perm_back)/sizeof(perm_back[0]));
        for(int i=0; i<test_top + 1; i++){
            perm[i][0] = top + i;
            perm[i][1] =  perm_front[i];
        }
        for(int i=test_top + 1; i < test_bottom; i++){
            perm[i][0] = top + i;
            perm[i][1] =  perm_middle[i - test_top - 1];
        }
        for(int i= test_bottom; i < bottom - top + 1; i++){
            perm[i][0] = top + i;
            perm[i][1] = perm_back[i - test_bottom];
        }


        perm_rows(perm, bottom - top + 1);

        top = test_top + top + 1;
        bottom = test_bottom + top - 1;
    }

}

//Pivoting step over all test intervals
void compound_row_step(int left, int right){
    if(bottom - top < min_vertical_interval || right - left < min_horizontal_interval)
        return;
    row_step(left, right);
    compound_row_step(left, (left + right)/2);
    compound_row_step((left + right)/2 + 1 ,right);
}

//Full row pass
void secondary_row_step(){
    if(bottom - top < min_vertical_interval)
        return;
    int t = top;
    int b = bottom;
    compound_row_step(0,n-1);
    top = t;
    bottom = (t+b)/2;
    secondary_row_step();
    top = (t+b)/2 + 1;
    bottom = b;
    secondary_row_step();
}

//Compute partials sums to make comparisons more efficient
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


//Fully sort the rows
void full_row_step(){
    compute_partials();
    secondary_row_step();
    top = 0;
    bottom = n-1;
}

//Transpose the matrix.  Note to sort the columns we transpose the matrix, sort the rows using the above and then transpose back.
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

int main()
{
    start();
    cout<<"constructed"<<endl;
    srand(time(NULL));
    build_matrix();
    initialize();
    int num_iterations = 5;
    for(int i=0; i< num_iterations; i++){
        cout<<i<<endl;
        full_row_step();
        transpose();
        full_row_step();
        transpose();
        cout<<sumsq(truth, base)<<endl;
    }

    //Write matrices to output file
    ofstream myfile;
    myfile.open(to_string(n) + "new.csv");
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
    myfile2.open (to_string(n) + "new_truth.csv");
    for(int i=0; i<n; i++){
        s = "";
        for(int j=0; j<n; j++){
            s += (to_string(round(8 * truth[i][j])) + ",");
        }
        s += "\n";
        myfile2 << s;
    }
    myfile2.close();




    return 0;
}
