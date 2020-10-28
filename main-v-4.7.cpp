// list1 and list2 in cmbr mergd into list
// global cumualative instance sum variable added
// tepmporary data structure distributed over grid cells
// grid approach changed. Added more dimentions
// isDeleted added for cmbr_arr and cmbr_map
// grid update at the all-to-all matching level
// created a new data structure to save both cbr_map amd cmbr_arr info together 

#include <iostream>
#include <vector>
#include <array>
#include <bitset>
#include<string>
#include <chrono>
#include <set>
#include <bits/stdc++.h> 

// distance threshold
// float const DIST = 5.0;
float const DIST = 200.0;

// number features
int const FMAX = 13;
// int const FMAX = 3;

#include "readFille-v-4.1.cpp"

using namespace std::chrono;

// prevalence threshold
// double const PI = 0.3;
double const PI = 0.25;
// double const PI = 0;
int prev_size = 0;
// new id for cmbr
int CMBR_ID = FMAX;
// max combination id in a cell
int MAX_COMB = FMAX+1;
int prevalent_count = 0; 

// saves all counts for features
vector<int> fcount(FMAX);
// feature id list
// int feature_ids[FMAX] = {1, 5, 8, 9, 10, 14, 20, 24, 28, 39, 40, 42, 43};

// data structure to save MBR info
struct mbr {
    float x1;
    float y1;
    float x2;
    float y2;
    bool empty = true;
};

// data structure hold cmbr info
struct cmbr {
    bitset<FMAX> combination; 
    int featureCount = -1;
    bool isDeleted = false;
    bool isAccessed = false;
    vector<mbr> cmbr_array;
    vector<int> instnace_list;
    vector<set<int>> inst_array;
};


vector<set<int>> tmp_inst_array(FMAX);

// saves all mbr information
vector<vector<vector<vector<mbr>>>>  mbr_array;
vector<vector<vector<int>>> instance_sum; 
vector<vector<vector<vector<cmbr>>>>  cmbr_map;


bool tp_flag = false;
bool mp_flag = false;
bool le_flag = false;

// print time info
void print_time(string str){
    if (tp_flag)
    {
        cout << "Time -> " << str << endl;        
    }
}

// print messge
void print_message(string str) {
    if (mp_flag)
    {
        cout << str << endl;        
    }
}

// print less messge
void print_message_less(string str) {
    if (le_flag)
    {
        cout << str << endl;        
    }
}

// calculate MBR for a given datapoint
mbr getMBR(float px, float py) {
    mbr box;
    box.x1 = px - DIST;
    box.y1 = py - DIST;
    box.x2 = px + DIST;
    box.y2 = py + DIST;

    return box;
}

// returns MBR list for a selected feature
void getMBRList(struct table_row *data) {

    int j = 0, k = 0;
    int row, col;
    mbr temp;

    print_message("Grid set Start...");
    for (int i = 0; i < ROWS; ++i) { 
        // check if feature id changes
        if (j != (data[i].id - 1)) {
            j = (data[i].id - 1);
            k++;
        }

        temp = getMBR(data[i].x, data[i].y);
        col = floor((temp.x1 - GRID_MIN_X) / (DIST * 2));
        row = GRID_ROWS - 1 - floor((temp.y2 - GRID_MIN_Y) / (DIST * 2));

        // calculate MBR using the getMBR() and assign it to the relavant feature instance
        mbr_array[row][col][k].push_back(temp);
        fcount[k] += 1; 
    }
    print_message("Grid set...");
    // calculate cumulative sums
    for (int i = 0; i < GRID_ROWS; ++i)
    {
        for (int j = 0; j < GRID_COLS; ++j)
        {
            for (int k = 0; k < FMAX; ++k)
            {                          
                if (j + 1 < GRID_COLS)
                { 
                    instance_sum[i][j + 1][k] = mbr_array[i][j][k].size() + instance_sum[i][j][k];          
                } else if(i + 1 < GRID_ROWS) { 
                    instance_sum[i + 1][0][k] = mbr_array[i][j][k].size() + instance_sum[i][j][k];          
                }
            }
        }
    }

}

float getMin(float a, float b) {
    if (a < b)
    {
        return a;
    }
    return b;
}

float getMax(float a, float b) {
    if (a > b)
    {
        return a;
    }
    return b;
}

// returns CMBR for a given two MBRs
mbr calculateCMBR(float ax1, float ay1, float ax2, float ay2, float bx1, float by1, float bx2, float by2) {
    mbr c;
    //if (ax2 > bx1 && ax1 < bx2 && ay2 > by1 && ay1 < by2) {}      
    c.x1 = getMax(ax1, bx1);
    c.y1 = getMax(ay1, by1);
    c.x2 = getMin(ax2, bx2);
    c.y2 = getMin(ay2, by2);
    if(!(c.x1 > c.x2 || c.y1 > c.y2))   
    {
        c.empty = false;
    }
    return c;
}

// returns a cmbr structure for selected 2 MBRs with the count of CMBRs
void getCMBRLayerWCount2(int fid1, int fid2, int crow, bool cmbrFlag, bitset<FMAX> comb, int fCount, int s) {

    // temperary array sizes
    int arr_size = 20000;
    // 1D array to hold layer CMBRs
    vector<mbr> arr(arr_size);
    // arr.resize(arr_size);
    int arr_loc = 0;

    // temperary array sizes
    int li_size = 20000;
    // 1D array to hold layer CMBRs
    vector<int> li(arr_size);
    // arr.resize(arr_size);
    int li_loc = 0;

    // temporary cmbr reference
    cmbr ret;

    // temporary varibale to track of calculated CMBR
    mbr cmbr_v;

    bool ZeroFlag = false;

    int c, r, n1row, n1col,  n2row, n2col, destR, destC;

    int n1[7] = {0, 0, 0, 0, 1, 2, 3};
    int n2[7] = {0, 1, 2, 3, 0, 0, 0};

    vector<mbr>* box1;
                        
    auto start = high_resolution_clock::now();

    // traverse grid
    for (int row = 0; row < GRID_ROWS; ++row)
    {
        for (int col = 0; col < GRID_COLS; ++col)
        {
            // check if we have anything to match in the given cell
            if (cmbrFlag)
            {
                if (cmbr_map[row][col][crow].size() >= fid1 +1 )
                {
                    c =  cmbr_map[row][col][crow][fid1].cmbr_array.size();                    
                } else {
                    c = -1;
                }
            } else {
                c = mbr_array[row][col][fid1].size();
            }
            // if both features do not exists in the cell, we may skip the cell
            if ( c <= 0 && mbr_array[row][col][fid2].size() <= 0)
            {
                continue;
            }


            n1row = row;
            n1col = col;
            for (int ii = 0, iii = 0; ii < 7; ++ii, ++iii)
            {
                if (ii == 4)
                {
                    iii = 1;
                }
                if (ii < 4)
                {
                    n2row = row + iii/2;
                    n2col = col + iii%2;
                } else {
                    n2row = row;
                    n2col = col;
                    n1row = row + iii/2;
                    n1col = col + iii%2;
                }

                // destination row and column
                destR = row + iii/2;
                destC = col + iii%2;
                
                if(n1row >= GRID_ROWS || n1col >= GRID_COLS || n2row >= GRID_ROWS || n2col >= GRID_COLS) 
                {
                    continue;
                }

                // change MBR or CMBR for first box when checking for CMBR
                if (cmbrFlag)
                {
                    // if there are no CMBRs in the selected cell, skip calculation and move to next
                    if (cmbr_map[n1row][n1col][crow].size() < fid1 + 1)
                    {
                        continue;
                    }
                    box1 = &cmbr_map[n1row][n1col][crow][fid1].cmbr_array;
                } else {
                    box1 = &mbr_array[n1row][n1col][fid1];
                }

                                
                for (int i = 0; i < (*box1).size(); ++i)
                {
                                        
                    for (int j = 0; j < mbr_array[n2row][n2col][fid2].size(); ++j)
                    {
                        

                        cmbr_v = calculateCMBR((*box1)[i].x1, (*box1)[i].y1, (*box1)[i].x2, (*box1)[i].y2, mbr_array[n2row][n2col][fid2][j].x1, mbr_array[n2row][n2col][fid2][j].y1, mbr_array[n2row][n2col][fid2][j].x2, mbr_array[n2row][n2col][fid2][j].y2);
                        if (!cmbr_v.empty)
                        {  
                            if (arr_loc >= arr_size)
                            {
                                // if arr is full, we may increaase the space
                                arr_size *= 3;
                                arr.resize(arr_size);                                 
                            }

                            if (li_loc+fCount-1 >= li_size)
                            {
                                li_size *=3;
                                li.resize(li_size);
                            }
                            // saving instance ids
                            arr[arr_loc++] = cmbr_v; 
                            li[li_loc++] = j + instance_sum[n2row][n2col][fid2];    

                            if (cmbrFlag)
                            {
                                for (int ft1 = 0; ft1 < fCount-1; ++ft1)
                                {
                                    li[li_loc] = cmbr_map[n1row][n1col][crow][fid1].instnace_list[i*(fCount-1) + ft1];  
                                    tmp_inst_array[ft1+1].insert(li[li_loc++]);                                  
                                }   
                            } else {
                                li[li_loc++] = i + instance_sum[n1row][n1col][fid1];
                                tmp_inst_array[1].insert(i + instance_sum[n1row][n1col][fid1]);
                            }                   
                         
                            tmp_inst_array[0].insert(j + instance_sum[n2row][n2col][fid2]);                           
                            if (destR == 0 && destC == 0)
                            {
                                ZeroFlag = true;
                            }
                        } 
                    }
                }
                ret.combination = comb;  
                ret.featureCount = fCount; 

                if (arr_loc > 0)
                { 
                    // if cmbr_map local sizes are not synch with latest, we may resize fid2 refers to step 1,2,...
                    if (cmbr_map[destR][destC][fid2-1].size() < pow(2, (fid2)))
                    {
                        cmbr_map[destR][destC][fid2-1].resize(pow(2, (fid2)));
                    } 

                    if (cmbr_map[destR][destC][fid2-1][s].isAccessed)
                    {

                        cmbr_map[destR][destC][fid2-1][s].cmbr_array.insert(cmbr_map[destR][destC][fid2-1][s].cmbr_array.end(), arr.begin(), arr.begin() + arr_loc);
                        cmbr_map[destR][destC][fid2-1][s].instnace_list.insert(cmbr_map[destR][destC][fid2-1][s].instnace_list.end(), li.begin(), li.begin() + li_loc);
                        
                        for (int fiid = 0; fiid < fCount; ++fiid)
                        {
                            merge(cmbr_map[destR][destC][fid2-1][s].inst_array[fiid].begin(), cmbr_map[destR][destC][fid2-1][s].inst_array[fiid].end(), tmp_inst_array[fiid].begin(), tmp_inst_array[fiid].end(), inserter(cmbr_map[destR][destC][fid2-1][s].inst_array[fiid], cmbr_map[destR][destC][fid2-1][s].inst_array[fiid].begin()));
                        }
                    } else {
                        ret.cmbr_array.insert(ret.cmbr_array.begin(), arr.begin(), arr.begin() + arr_loc);
                        ret.instnace_list.insert(ret.instnace_list.begin(), li.begin(), li.begin() + li_loc);

                        ret.inst_array = tmp_inst_array;
                        ret.isAccessed = true;

                        cmbr_map[destR][destC][fid2-1][s] = ret;
                    }               
                    tmp_inst_array.clear();
                    tmp_inst_array.resize(12);

                    arr_loc = 0;
                    li_loc = 0;
                }             
                ret = {}; //reset ret
            }          
        }
    }
    // if cell 0 did not get any CMBRs, we may add a null to referenece it
    if (!ZeroFlag)
    {
        ret.combination = comb;  
        ret.featureCount = fCount;
        cmbr_map[0][0][fid2-1].push_back(ret);
        ret = {};
    } 
}

//erasing filtered out cmbr_maps layerwise
void erase_combination(int k, int comb)
{
    auto start = high_resolution_clock::now();

    cmbr_map[0][0][k][comb].isDeleted = true;
    cmbr_map[0][0][k][comb].cmbr_array.clear();
    cmbr_map[0][0][k][comb].instnace_list.clear();
    cmbr_map[0][0][k][comb].inst_array.clear();
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start); 
    print_time("Function: erase_cmbr_map " + to_string(duration.count()));
}

void cmbr_filter_combinationwise(int k, int comb)
{
    vector<set<int>> grid_pr(FMAX);
    int cmbr_count = 0;

    int fc = 0;

    auto start = high_resolution_clock::now();
    //for merging all the grid cells to cell 0
    for (int row = 0; row < GRID_ROWS; row++)
    {
        for (int col = 0; col < GRID_COLS; col++)
        {   
            if(cmbr_map[row][col][k].size() > 0 && cmbr_map[row][col][k][comb].isAccessed)
            {
                cmbr_count += cmbr_map[row][col][k][comb].cmbr_array.size();                 
                for(int ft = 0; ft <  cmbr_map[row][col][k][comb].featureCount; ft++)                   
                {
                    merge(grid_pr[ft].begin(), grid_pr[ft].end(), cmbr_map[row][col][k][comb].inst_array[ft].begin(), cmbr_map[row][col][k][comb].inst_array[ft].end(), inserter(grid_pr[ft], grid_pr[ft].begin()));    
                }//ft
            }// is delete if check cmbr_map
            
        }//col
    }//row
    // check if CMBR count is > 0
    if (cmbr_count > 0)
    {
        double pr = 1.0;
        bitset<FMAX> bit_comb = cmbr_map[0][0][k][comb].combination;
        string result = "";
        
        int pos = FMAX;
        int ft_loop = 0;

        for(int ft = grid_pr.size()-1; ft >= 0 ; ft--)    
        {
            if(grid_pr[ft].size() > 0)
            {

                if (ft == 0)
                {
                    fc++;
                }
                ft_loop += 1;

                //for printing the bit comb features only 1st time
                if(ft_loop == 1){
                    result += "(";
                    int flag =0;
                    for(int n = FMAX-1; n>=0; n--){
                        if(bit_comb[n]==1 && flag==1){
                            result += ", "; 
                        }
                        if(bit_comb[n]==1){
                            result += to_string(FMAX-1-n);
                            flag=1;
                        }       
                    }
                    result += ") -> (";
                    // result += ", ";
                }
                double total_instances = 0;
                for(int i = pos-1; i >= 0; i--){
                    if(bit_comb[i] == 1)
                    {
                        pos = i;
                        break;
                    }
                }
                total_instances = fcount[FMAX-1-pos];
                result += to_string(grid_pr[ft].size()) + ":" + to_string(fcount[FMAX-1-pos]);
                if(ft != 0){
                    result += ", ";
                }
                else{
                    result += ") -> " + to_string(cmbr_count);
                    // result += ", " + to_string(cmbr_count[comb]) + "\n";
                }
                double temp_pr = grid_pr[ft].size()/total_instances;
                if(temp_pr < pr){
                    pr = temp_pr;               
                }
            }//only those features with size > 0            
        }//ft
        if( pr < PI ){
            erase_combination(k, comb);
            result = "";
        }
        else{
            prevalent_count++;
            print_message(result);// << endl;
        }
        // }
    }//isdeleted check

    grid_pr.clear();

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start); 
    print_time("Function: cmbr_filter_layerwise " + to_string(duration.count()));

}

// check whether to proceed with the new combination or not
bool isCombinationPrevalent(int k, bitset<FMAX> combination, int fCount) {
    cmbr tmp;
    tmp.combination = combination;  
    tmp.featureCount = fCount;

    int count;


    for (int ncomb = 0; ncomb < cmbr_map[0][0][k].size(); ++ncomb)
    {
        if (cmbr_map[0][0][k][ncomb].isDeleted)
        {
            count = cmbr_map[0][0][k][ncomb].featureCount;

            for (int i = 0; i < FMAX; ++i)
            {
                if (combination[i] == 1 &&  cmbr_map[0][0][k][ncomb].combination[i] == 1 && i != FMAX-2-k)
                {
                    if (count != 2)
                    {
                        count--;                        
                    } else {                        
                        tmp.isDeleted = true;
                        cmbr_map[0][0][k].push_back(tmp);
                        return false;
                    }
                }
            }
        }
    }
    return true; 
}


// traverse through all the steps.  Maintains count for each CMBR
void buildCMBRList() {
    // numbers layers need to consider. (#features-1)
    int layers = FMAX-1; 
    bitset<FMAX> temp_comb; 

    int fc = 0, comb_loc;

    auto start = high_resolution_clock::now();
    for (int k = 0; k < layers; ++k)
    {  
        comb_loc = 0;

        // layer 1. Instance wise CMBR only.No previos layer
        if (k == 0) {

            // get combination bit pattern
            temp_comb.reset();               
            temp_comb[FMAX-1-k] = 1;
            temp_comb[FMAX-2-k] = 1;

            getCMBRLayerWCount2(k, k+1, 0, false, temp_comb, 2, comb_loc);
            cmbr_filter_combinationwise(k, comb_loc++); 

        } else {
            // find all the CMBRs from 1st feature to K+1 feature 
            auto start_i = high_resolution_clock::now();
            for (int i = 0; i <= k; ++i)
            {  
                // get combination bit pattern
                temp_comb.reset();
                temp_comb[FMAX-1-i] = 1;
                temp_comb[FMAX-2-k] = 1;

                getCMBRLayerWCount2(i, k+1, 0, false, temp_comb, 2, comb_loc);                  
                cmbr_filter_combinationwise(k, comb_loc++); 
                
            }

            // find CMBRs with K+1 feature and previous layer CMBRs
            for (int i = 0; i < k ; ++i)
            {   
                auto start_jj = high_resolution_clock::now();  
                for(int jj=0; jj< cmbr_map[0][0][i].size(); ++jj)
                {
                    if (!cmbr_map[0][0][i][jj].isDeleted)
                    { 
                        // get combination bit pattern
                        temp_comb.reset();
                        temp_comb = cmbr_map[0][0][i][jj].combination; // take combintion id from previous step
                        temp_comb[FMAX-2-k] = 1; 
    
                        if (isCombinationPrevalent(k, temp_comb, cmbr_map[0][0][i][jj].featureCount+1))
                        {
                            getCMBRLayerWCount2(jj, k+1, i, true, temp_comb, cmbr_map[0][0][i][jj].featureCount+1, comb_loc);                 
                            cmbr_filter_combinationwise(k, comb_loc); 
                        }
                        comb_loc++;
                    } 

                }
                auto stop_jj = high_resolution_clock::now();
                auto duration_jj = duration_cast<microseconds>(stop_jj - start_jj); 
                print_time("Function: buildCMBR jj to i loop " + to_string(duration_jj.count()));                               
            }

            auto stop_i = high_resolution_clock::now();
            auto duration_i = duration_cast<microseconds>(stop_i - start_i); 
            print_time("Function: buildCMBR k to i loop " + to_string(duration_i.count()));
        }       

        print_message("Layer " + to_string(k) + " Built. Prevalent count: " + to_string(prevalent_count)); 
        prevalent_count = 0;      
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start); 
    print_message("Function: buildCMBRList " + to_string(duration.count()));
    return ;
} 

int main(int argc, char* argv[])
{
    
    // freopen ("out_v_4.2_1676.txt","w",stdout);  
    auto start = high_resolution_clock::now();

    if (argc > 1 && strcmp(argv[1], "true") == 0)
    {
        le_flag = true;
    }
    if (argc > 2 && strcmp(argv[2], "true") == 0)
    {
        mp_flag = true;
    }
    if (argc > 3 && strcmp(argv[2], "true") == 0)
    {
        tp_flag = true;
    }

    // read data into a table_row structure type 1D array
    struct table_row *dat;
    // dat = createArray("data/Seattle2012_1.csv");
    // dat = createArray("data/newData/Seattle2012_1676.csv");
    dat = createArray("data/newData/Point_Of_Interest_modified.csv");

    // print_message(to_string(mbr_array.max_size()));
    // print_message(to_string(instance_sum.max_size()));
    // print_message(to_string(cmbr_map.max_size()));

    print_message("Total rows: " + to_string(ROWS));


    print_message("Initialize start..");

    // initialize sizes for main 3 data structures 
    mbr_array.resize(GRID_ROWS, vector<vector<vector<mbr>>>(GRID_COLS, vector<vector<mbr>>(FMAX)));
    print_message("mbr_array initialized..");
    instance_sum.resize(GRID_ROWS, vector<vector<int>>(GRID_COLS, vector<int>(FMAX, 0))); 
    print_message("instance_sum initialized..");
    cmbr_map.resize(GRID_ROWS, vector<vector<vector<cmbr>>>(GRID_COLS, vector<vector<cmbr>>(FMAX-1)));
    print_message("cmbr_map initialized..");

    // test ----- START ----

    // number of rows in the data file
    // ROWS = 6;

    // // grid origin
    // GRID_MIN_X = 0.0, GRID_MIN_Y = 0.0;

    // // grid number of rows
    // GRID_ROWS = ceil((80 - GRID_MIN_X)/(DIST * 2)) + 2;

    // // grid number of columns
    // GRID_COLS = ceil((80 - GRID_MIN_Y)/(DIST * 2)) + 2;


    // struct table_row dat[14] = {{1, 50,50}, {1, 15,15}, {1, 82, 25}, {1, 13, 20},
    //                                  {2, 55, 50}, {2, 56, 50}, {2, 80, 80}, {2, 11, 11},
    //                                  {3, 11, 15}, {3, 65, 65},
    //                                  {4, 35, 25}, {4, 15, 15}, {4, 50, 50},
    //                                  {5, 40, 40}};

    // struct table_row dat[6] = {{1, 65,65}, {1, 13, 80}, {1, 25,25},
    //                             {2, 80, 80}, {2, 30, 25},
    //                             {3, 28, 28}};


    // test ----- END -----

    // calculate MBR for all the datapoints. 
    // returns a 2D array. 1st-D : Features, 2nd-D: instances per each feature 
    getMBRList(dat); 

    print_message("mbr array constructed");
    // cout << GRID_ROWS << " " << GRID_COLS << endl;

    // // print mbr array
    // for (int r = 0; r < GRID_ROWS; ++r)
    // {
    //     for (int c = 0; c < GRID_COLS; ++c)
    //     {
    //         for (int i = 0; i < mbr_array[r][c].size(); ++i)
    //         {
    //             for (int ii = 0; ii < mbr_array[r][c][i].size(); ++ii)
    //             {
    //                 // cout << "Row: " << r << " Column: " << " Feature id: " << i << " Count: " << mbr_array[r][c][i].size() << endl;
    //                 cout << r << "," << c << "," << i << "," << mbr_array[r][c][i][ii].x1 << "," << mbr_array[r][c][i][ii].y1 << "," << mbr_array[r][c][i][ii].x2 << "," << mbr_array[r][c][i][ii].y2 << endl;
    //                 cout << instance_sum[r][c][i] << endl;
    //             }                
    //         }
    //     }
    // }

    // build CMBR tree 
    buildCMBRList();
    print_message("cmbr layers constructed");

    // print cmbr array
    // for (int r = 0; r < GRID_ROWS; ++r)
    // {
    //    for (int c = 0; c < GRID_COLS; ++c)
    //      {
    //          for (int i = 0; i < cmbr_map[r][c].size(); ++i)
    //          {
    //              for (int j = 0; j < cmbr_map[r][c][i].size(); ++j)
    //              {
    //                 for (int ii = 0; ii < cmbr_map[r][c][i][j].cmbr_array.size(); ++ii)
    //                 {                        
    //                     cout << r << ", " << c << ", " << i << ", " << j << " -* " << cmbr_map[r][c][i][j].cmbr_array[ii].x1 << ", "  << cmbr_map[r][c][i][j].cmbr_array[ii].y1 << ", " << cmbr_map[r][c][i][j].cmbr_array[ii].x2 << ", " << cmbr_map[r][c][i][j].cmbr_array[ii].y2 << endl;                           
    //                 }
    //                 // cout << "List1" << endl;
    //                 for (int ii = 0; ii < cmbr_map[r][c][i][j].list1.size(); ++ii)
    //                 {
    //                     cout << "l1[";
    //                     for (int jj = 0; jj < cmbr_map[r][c][i][j].list1[ii].size(); ++jj)
    //                     {
    //                         cout << cmbr_map[r][c][i][j].list1[ii][jj] << ", ";
    //                     }
    //                     cout << "]" << endl;
    //                 }
    //                 // cout << "List2" << endl;
    //                 for (int ii = 0; ii < cmbr_map[r][c][i][j].list2.size(); ++ii)
    //                 {
    //                     cout << "l2[";
    //                     for (int jj = 0; jj < cmbr_map[r][c][i][j].list2[ii].size(); ++jj)
    //                     {
    //                         cout << cmbr_map[r][c][i][j].list2[ii][jj] << ", ";
    //                     }
    //                     cout << "]" << endl;
    //                 }
    //                  // if (!cmbr_map[0][0][i][j].isDeleted && cmbr_map[r][c][i][j].cmbr_array.size() > 0)
    //                  // {
    //                  //     cout << cmbr_map[r][c][i][j].combination << "[" << cmbr_map[r][c][i][j].cmbr_array.size() << "] ";                
    //                  // }
    //              }
    //              //cout << endl;
    //          }
    //      }
    //  }

    // clear global vectors
    mbr_array.clear();
    instance_sum.clear();
    cmbr_map.clear();
    tmp_inst_array.clear();
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start); 
    print_message_less("Function: main " + to_string(duration.count()));


    // fclose(stdout);

    return 0;
}