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
float const DIST = 5.0;

// number features
int const FMAX = 13;
// int const FMAX = 3;

#include "readFille-v-4.1.cpp"

using namespace std::chrono;

// prevalence threshold
double const PI = 0.3;
// double const PI = 0;
int prev_size = 0;
// new id for cmbr
int CMBR_ID = FMAX;
// max combination id in a cell
int MAX_COMB = FMAX+1;
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
    int featureCount = 0;
    bool isDeleted = false;
    vector<mbr> cmbr_array;
    vector<vector<int>> list1;
    vector<vector<int>> list2;
};

int GROWS = 802;
int GCOLS = 602;

// saves all mbr information
vector<vector<vector<vector<mbr>>>>  mbr_array(GROWS, vector<vector<vector<mbr>>>(GCOLS, vector<vector<mbr>>(FMAX)));

// data structure saves instance cumulative sum 
vector<vector<vector<int>>> instance_sum(GROWS, vector<vector<int>>(GCOLS, vector<int>(FMAX, 0))); 

// // data structure saves MBR infor per cell
// struct mbr_cell {
//     vector<vector<mbr>>  mbrs(FMAX);
//     vector<int> instance_sum(FMAX, 0);
// };

// // saves all mbr information
// vector<vector<mbr_cell>>  mbr_array(801, vector<mbr_cell>(601));

// 2D vector to keep track of all the combinations and instances realted
vector<vector<vector<vector<cmbr>>>>  cmbr_map(GROWS, vector<vector<vector<cmbr>>>(GCOLS, vector<vector<cmbr>>(FMAX-1)));

// intermediate data structure to hold unique instances ids for cmbrs in a perticular step in each cell
vector<vector<vector<vector<set<int>>>>> instance_array;

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

    for (int i = 0; i < ROWS; ++i) { 
        // check if feature id changes
        if (j != (data[i].id - 1)) {
            j = (data[i].id - 1);
            k++;
        }

        temp = getMBR(data[i].x, data[i].y);
        col = floor((temp.x1 - GRID_MIN_X) / (DIST * 2));
        row = GRID_ROWS - 1 - floor((temp.y2 - GRID_MIN_Y) / (DIST * 2));
        // cout << "row " << row << " col " << col << " fid " << k << " point " << floor(temp.y2 - GRID_MIN_Y) / (DIST * 2) << ", " << temp.y2 << endl;

        // calculate MBR using the getMBR() and assign it to the relavant feature instance
        mbr_array[row][col][k].push_back(temp);
        // cout << "row " << row << " col " << col << " fid " << k << " point " << data[i].x << ", " << data[i].y << endl;        
        fcount[k] += 1; 
    }

    // calculate cumulative sums
    for (int i = 0; i < mbr_array.size(); ++i)
    {
        for (int j = 0; j < mbr_array[i].size(); ++j)
        {
            for (int k = 0; k < FMAX; ++k)
            {                          
                if (j + 1 < GRID_COLS)
                { 
                    instance_sum[i][j + 1][k] = mbr_array[i][j][k].size() + instance_sum[i][j][k];      
                    // if (mbr_array[i][j][k].size() > 0)
                    // {
                    //     cout << "ins: " << i << ", " << j << ", " << k << "- " << instance_sum[i][j + 1][k] << endl;
                    // }    
                } else if(i + 1 < GRID_ROWS) { 
                    instance_sum[i + 1][0][k] = mbr_array[i][j][k].size() + instance_sum[i][j][k];    
                    // if (mbr_array[i][j][k].size() > 0)
                    // {
                    //     cout << "ins: " << i << ", " << j << ", " << k << "- " << instance_sum[i + 1][j][k] << endl;
                    // }         
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

void print_time(string str){
    // cout << "Time -> " << str << endl;
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

// read previous stel cmbr instance combination and save it to list1
vector<vector<int>> instanceCombinationBuild(vector<vector<int>> list1, vector<vector<int>> map_l1, vector<vector<int>> map_l2, int k, int r, int c) {
    vector<int> ttlist1;
    vector<vector<int>> ttlist2;
    int tt=0;
    auto start = high_resolution_clock::now();
    // cout << "===" << list1.size() << endl;
    // cout << "===" << map_l1.size() << endl;

    // cout << "list 1 " << endl; 
    // for (int i = 0; i < list1.size(); ++i)
    // {
    //     for (int j = 0; j < list1[i].size(); ++j)
    //     {
    //         cout << list1[i][j] << ", ";
    //     }
    //     cout << endl;
    // }

    // cout << "map 1 " << endl; 
    // for (int i = 0; i < map_l1.size(); ++i)
    // {
    //     for (int j = 0; j < map_l1[i].size(); ++j)
    //     {
    //         cout << map_l1[i][j] << ", ";
    //     }
    //     cout << endl;
    // }

    // cout << "map 1 " << endl; 
    // for (int i = 0; i < map_l2.size(); ++i)
    // {
    //     for (int j = 0; j < map_l2[i].size(); ++j)
    //     {
    //         cout << map_l2[i][j] << ", ";
    //     }
    //     cout << endl;
    // }

    for (int aa = 0; aa < list1.size(); ++aa)
    {
        // goes through list2 of previos k-1 step to find the insatnce which made new CMBRs
        for (int bb = 0; bb < map_l2.size(); ++bb)
        {
            if(bb == 0) {
                // track the size of list2 first row
                tt = map_l2[bb].size();
            }
            else {         
                // if current row does not have search index, move to next row
                tt += map_l2[bb].size();
            }
            // cout << "===***" << list1[aa][0] << endl;
            // cout << "===***--" << tt << endl;

            // check if current row contains the index in the returned list1
            if ( tt > list1[aa][0])
            {
                // assign list1 row into new row
                ttlist1.insert(ttlist1.end(), map_l1[bb].begin(), map_l1[bb].end());
                tt = map_l2 [bb].size() - tt + list1[aa][0];
                
                // add list2 cell value to the same row. We have the old CMBR instance now
                ttlist1.push_back(map_l2[bb][tt]);
                
                // push the new CMBR instance combination into a 2D array
                ttlist2.push_back(ttlist1);
                // cout << r << " " <<  c << " " << k << " tt size: " << ttlist1.size() << endl;
                // update global temporary cmbr unique instance set
                // cout << "bu-Push: i: ";
                for (int ii = 0, j = ttlist1.size()-1; ii < ttlist1.size(); ++ii, --j)
                {
                    // cout << "ii " << ii << endl;
                    instance_array[r][c][k][ii+1].insert(ttlist1[j]);
                    // cout << ttlist1[j] << "-";
                }               
                // instance_array[k].insert(ttlist1.begin(), ttlist1.end()); //test this not sure

                //clear temporary 1D array. Ready for next combination creation
                ttlist1.clear();
                // stop looking for more. We have only 1 row-cell combination at a time
                break;
            }                    
        }
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start); 
    print_time("Function: instanceCombinationBuild " + to_string(duration.count()));
    return ttlist2;                 
}

// returns a cmbr structure for selected 2 MBRs with the count of CMBRs
void getCMBRLayerWCount2(int fid1, int fid2, int crow, bool cmbrFlag, bitset<FMAX> comb, int fCount) {

    // 1D array to hold layer CMBRs
    vector<mbr> arr;
    arr.reserve(20000);
    // temporary cmbr reference
    cmbr ret;

    // temporary varibale to track of calculated CMBR
    mbr cmbr_v;

    // l1 and l2 maintains the lists from mbr1 and mbr2. l1[0] l2[0] will give 0th CMBR instances
    vector<vector<int>> l1;
    vector<vector<int>> l2;
    l1.reserve(2000);
    l2.reserve(2000);

    // temporary arrays to hel build l1 and l2 rows
    vector<int> t1;
    vector<int> t2;
    t1.reserve(2);
    t2.reserve(1000);

    int c, r, n1row, n1col,  n2row, n2col, s, t;

    int n1[7] = {0, 0, 0, 0, 1, 2, 3};
    int n2[7] = {0, 1, 2, 3, 0, 0, 0};

    vector<mbr> box1;

    int init_layer_count = cmbr_map[0][0][fid2-1].size() +1 ;  

    // cout <<  "first insert cmbr_map: " << init_layer_count << endl;                      


    auto start = high_resolution_clock::now();

    // traverse grid
    for (int row = 0; row < GRID_ROWS; ++row)
    {
        for (int col = 0; col < GRID_COLS; ++col)
        {
            if (mbr_array[row][col].size() <= 0)
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
                
                if(n1row >= GRID_ROWS || n1col >= GRID_COLS || n2row >= GRID_ROWS || n2col >= GRID_COLS) 
                {
                    continue;
                }

                // cout << "r " << row << " c " << col << endl;
                // change MBR or CMBR for first box when checking for CMBR
                if (cmbrFlag)
                {
                    // cout << "fid1 " << fid1 << endl; 
                    // cout << "id1 " << id1 << endl; 
                    // if (cmbr_map[n1row][n1col][crow][fid1].cmbr_array.size() > 0)
                    // {
                    //     cout << "fid2 " << fid2 << " fid1 " << fid1 << " crow " << crow << " new fid1=" << cmbr_map[n1row][n1col][crow][fid1].cmbr_array.size() << endl;                
                        
                    // }
                    box1 = cmbr_map[n1row][n1col][crow][fid1].cmbr_array;
                } else {
                    box1 = mbr_array[n1row][n1col][fid1];
                }
                                
                for (int i = 0; i < box1.size(); ++i)
                {
                    // cout << "i " << i << " j size: " << mbr_array[n2row][n2col][fid2].size() << endl;
                    
                    for (int j = 0; j < mbr_array[n2row][n2col][fid2].size(); ++j)
                    {
                        // cout << "j " << j << endl;

                        cmbr_v = calculateCMBR(box1[i].x1, box1[i].y1, box1[i].x2, box1[i].y2, mbr_array[n2row][n2col][fid2][j].x1, mbr_array[n2row][n2col][fid2][j].y1, mbr_array[n2row][n2col][fid2][j].x2, mbr_array[n2row][n2col][fid2][j].y2);
                        if (!cmbr_v.empty)
                        {  
                            c = floor((cmbr_v.x1 - GRID_MIN_X) / (DIST * 2));
                            r = GRID_ROWS - 1 - floor((cmbr_v.y2 - GRID_MIN_Y) / (DIST * 2));
                            // cout << row + iii/2 << " " << col + iii%2 << " " << instance_array[row + iii/2][col + iii%2][init_layer_count].size() << " " << init_layer_count << endl;

                            arr.push_back(cmbr_v); 
                            t2.push_back(j + instance_sum[n2row][n2col][fid2]); 
                            instance_array[row + iii/2][col + iii%2][init_layer_count - 1][0].insert(j + instance_sum[n2row][n2col][fid2]);                           
                            // cout << "Push: j: " << j << " " << instance_sum[n2row][n2col][fid2] << endl;
                            // cout << comb << " -- " << row + iii/2 << " " << col + iii%2 << " " << cmbr_v.x1 << " " << cmbr_v.y1 << " " << cmbr_v.x2 << " " << cmbr_v.y2 << endl;
                        } 
                    }
                    // if there are any CMBRs for i, add ID i to l1
                    if (t2.size() > 0)
                    {
                        if (!cmbrFlag)
                        {
                            instance_array[row + iii/2][col + iii%2][init_layer_count - 1][1].insert(i + instance_sum[n1row][n1col][fid1]);
                            // cout << "Push: i: " << i << " " << instance_sum[n1row][n1col][fid1] << endl;                                            
                            t1.push_back(i + instance_sum[n1row][n1col][fid1]);
                        } else {
                            t1.push_back(i);
                        }
                        l1.push_back(t1); // push 1D arrays to 2D array
                        l2.push_back(t2); // push 1D arrays to 2D array
                        t1.clear(); //clear 1D array
                        t2.clear(); //clear 1D array   
                        t1.reserve(2);
                        t2.reserve(1000);                          
                    }
                }
                ret.combination = comb;  
                ret.featureCount = fCount;
                s = cmbr_map[row + iii/2][col + iii%2][fid2-1].size();  
                // cout << "init " << init_layer_count << " s " << s << endl;

                if (arr.size() > 0)
                {                    
                    // check if the CMBR pattern already exists in the selected cell
                    if (init_layer_count == s)
                    {

                        cmbr_map[row + iii/2][col + iii%2][fid2-1][s-1].cmbr_array.insert(cmbr_map[row + iii/2][col + iii%2][fid2-1][s-1].cmbr_array.end(), arr.begin(), arr.end());

                        cmbr_map[row + iii/2][col + iii%2][fid2-1][s-1].list2.insert(cmbr_map[row + iii/2][col + iii%2][fid2-1][s-1].list2.end(), l2.begin(), l2.end());

                        if (cmbrFlag)
                        {
                            // cout << "222 bqqqtttend inner " << endl;

                            // ret.list1 = instanceCombinationBuild(l1, cmbr_map[n1row][n1col][crow][fid1].list1, cmbr_map[n1row][n1col][crow][fid1].list2, init_layer_count - 1, row + iii/2, col + iii%2);                            
                            ret.list1 = instanceCombinationBuild(l1, cmbr_map[n1row][n1col][crow][fid1].list1, cmbr_map[n1row][n1col][crow][fid1].list2, init_layer_count - 1, row + iii/2, col + iii%2);                            
                            // cout << "4444 qqqtttend inner " << endl;
                            
                            cmbr_map[row + iii/2][col + iii%2][fid2-1][s-1].list1.insert(cmbr_map[row + iii/2][col + iii%2][fid2-1][s-1].list1.end(), ret.list1.begin(), ret.list1.end());                            
                        } else {
                            cmbr_map[row + iii/2][col + iii%2][fid2-1][s-1].list1.insert(cmbr_map[row + iii/2][col + iii%2][fid2-1][s-1].list1.end(), l1.begin(), l1.end());
                        }
                        // cout << "tttend inner " << endl;
                    } else {
                        ret.cmbr_array = arr;
                        ret.list1 = l1;
                        // check if matching is CMBR vs MBR                        
                        if (cmbrFlag)
                        {
                            ret.list1 = instanceCombinationBuild(l1, cmbr_map[n1row][n1col][crow][fid1].list1, cmbr_map[n1row][n1col][crow][fid1].list2, init_layer_count - 1, row + iii/2, col + iii%2);                            
                        } else {
                            ret.list1 = l1;
                        }
                        // cout << "rrrrrrrtttend inner " << endl;

                        // insert new feature combination if it does not exist in the cell 
                        cmbr_map[row + iii/2][col + iii%2][fid2-1].push_back(ret); 
                        // cout <<  "in else n insert cmbr_map: " << cmbr_map[row + iii/2][col + iii%2][fid2-1].size() << endl;                      
                    }
                    l1.clear();
                    l2.clear();
                    arr.clear();
                    l1.reserve(2000);
                    l2.reserve(2000);
                    arr.reserve(20000); 

                // insert new feature combination if it does not exist in the cell even if there are no cmbrs found in this cell
                } else if (init_layer_count != s) {
                    cmbr_map[row + iii/2][col + iii%2][fid2-1].push_back(ret);
                    // cout <<  "out else n insert cmbr_map: " << cmbr_map[row + iii/2][col + iii%2][fid2-1].size() << endl;                      
                }            
                ret = {}; //reset ret
            }          
        }
    }
}

//erasing filtered out cmbr_maps layerwise
void erase_cmbr_map(int k, vector<int> erase_list)
{
    //before erase
    //cout << "CMBR MAP before erase " << endl;
    //for(int i = 0; i<cmbr_map[k].size(); i++){
        //cout << cmbr_map[k][i].combination << " ";  
    //}
    //cout << endl;   
    // erase code last to first
    auto start = high_resolution_clock::now();
    for(int ii=erase_list.size()-1; ii >= 0 ; ii--)
    {
        // cout << "erasing " << cmbr_map[k][erase_list[ii]].combination << " index is " << erase_list[ii] << " k is " << k << endl;
        // cmbr_map[k].erase(cmbr_map[k].begin() + erase_list[ii]);
        // cmbr_arr[k].erase(cmbr_arr[k].begin() + erase_list[ii]);

        cmbr_map[0][0][k][erase_list[ii]].isDeleted = true;
        cmbr_map[0][0][k][erase_list[ii]].cmbr_array.clear();
        cmbr_map[0][0][k][erase_list[ii]].list1.clear();
        cmbr_map[0][0][k][erase_list[ii]].list2.clear();
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start); 
    print_time("Function: erase_cmbr_map " + to_string(duration.count()));
    //after erase
    //cout << "CMBR MAP after erase " << endl;
    //for(int i = 0; i<cmbr_map[k].size(); i++){
        //cout << cmbr_map[k][i].combination << " ";  
    //}
    //cout << endl;
    return;
}

//cmbr filter layerwise where k is the layer number
//vector<vector<vector<vector<set<int>>>>> instance_array; row col comb feature instance
void cmbr_filter_layerwise(int k)
{
    //vector<vector<int>> grid_pr;//(cmbr_map[0][0][k].size());
    vector<int> erase_list;
    auto start = high_resolution_clock::now();
    cout << "CMBR Filter called for Layer " << k << endl;
    //for merging all the grid cells to cell 0
    for (int row = 0; row < GRID_ROWS; row++)
    {
        for (int col = 0; col < GRID_COLS; col++)
        {   
            //prev_size += cmbr_map[row][col][k].size();
            for(int comb = 0; comb < instance_array[row][col].size(); comb++)
            {
                //each combination like AB, abc
                // if (!cmbr_map[row][col][k][comb].isDeleted)
                {
                    for(int ft = 0; ft < instance_array[row][col][comb].size(); ft++)                   
                    {
                        // cout << ft << ", ";
                        if(instance_array[row][col][comb][ft].size() > 0){                          
                            // cout << "Before size is " << instance_array[0][0][comb][ft].size() << " Combination is " << comb << " Feature is " << ft << endl;
                            // for(auto const &e: instance_array[0][0][comb][ft]){
                            //     cout << e << " ";
                            // } cout << endl;   

                            // cout << "cmbr Before size is " << cmbr_map[row][col][k][comb].cmbr_array.size() << " Combination is " << comb << " Feature is " << ft << endl;
                            // for (int i = 0; i < cmbr_map[row][col][k][comb].cmbr_array.size(); ++i)
                            // {
                            //     cout << i << " ";  
                            // }
                            // cout << endl;  

                            merge(instance_array[0][0][comb][ft].begin(), instance_array[0][0][comb][ft].end(), instance_array[row][col][comb][ft].begin(), instance_array[row][col][comb][ft].end(), inserter(instance_array[0][0][comb][ft], instance_array[0][0][comb][ft].begin()));
                            // cout << "After size is " << instance_array[0][0][comb][ft].size() << " Combination is " << comb << " Feature is " << ft << endl;
                            // for(auto const &e: instance_array[0][0][comb][ft]){
                            //     cout << e << " ";
                            // } cout << endl;
                        }// for merging and printing only
                        
                    }//ft
                }// is delete if check cmbr_map
            }//comb
        }//col
    }//row
    cout << "out" << endl;
    cout << "****" << endl;
    //now instance_array[0][0] will have the union info of all the cells
    for(int comb = 0; comb<instance_array[0][0].size(); comb++) 
    {
        if (!cmbr_map[0][0][k][comb].isDeleted)
        {
            // if(instance_array[0][0].size() > 0){
                double pr = 1.0;
                bitset<FMAX> bit_comb = cmbr_map[0][0][k][comb].combination;
                int pos = FMAX;
                for(int ft = instance_array[0][0][comb].size()-1; ft >= 0 ; ft--)   
                {
                    if(instance_array[0][0][comb][ft].size() > 0){
                        cout << "Comb is " << bit_comb << endl;
                        cout << "Feature is " << ft << " size is " << instance_array[0][0][comb][ft].size() << endl;
                        double total_instances = 0;
                        for(int i = pos-1; i >= 0; i--){
                            if(bit_comb[i] == 1)
                            {
                                pos = i;
                                break;
                            }
                        }
                        total_instances = fcount[FMAX-1-pos];
                        double temp_pr = instance_array[0][0][comb][ft].size()/total_instances;
                        if(temp_pr < pr){
                            pr = temp_pr;               
                        }
                    }//only those features with size > 0            
                }//ft
                if( pr < PI ){
                    cout << "Erasing Combination " << bit_comb << endl;
                    erase_list.push_back(comb);
                }
            // }
        }
    }//comb

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start); 
    print_time("Function: cmbr_filter_layerwise " + to_string(duration.count()));
    
    erase_cmbr_map(k, erase_list);
    erase_list.clear();

    return;
}


// check for selected CMBR
bool isCombinationValid(thirteenBits const c) {
    bool flag;

    for (int i = 0; i < seletedFeatures.size(); ++i)
    {
        flag = true;
        for (int j = 0; j < FMAX; ++j)
        {
          if (c[j] != seletedFeatures[i][j])
            {
                flag = false;
                break;
            }  
        }
        if (flag)
        {
            return true;
        }
    }
  
    return false;
}


// returns a 2D vector of all the CMBRs of the features.  Maintains count for each CMBR
void buildCMBRList() {
    // numbers layers need to consider. (#features-1)
    int layers = FMAX-1; 
    // int layers = 2; 

    bitset<FMAX> temp_comb; 

    int fc = 0;

    auto start = high_resolution_clock::now();
    for (int k = 0; k < layers; ++k)
    {  
        //cout <<"Layer " << k << " Building ..." << endl;
        //b = features[k+1]-1;

        // layer 1. Instance wise CMBR only.No previos layer
        if (k == 0) {
            //a = features[k]-1;      

            //cout << "Feature: " << k << " with feature: " << k+1 << endl;
            //cout << "ptr[a]= " << fcount[k] << " ptr[b]= " << fcount[k+1] << endl;

            // get combination bit pattern
            temp_comb.reset();               
            temp_comb[FMAX-1-k] = 1;
            temp_comb[FMAX-2-k] = 1;

            cout << temp_comb << endl;

            // resize instance array for the new feature combination
            // instance_array.resize(instance_array.size()+1, vector<set<int>>(2));
            // cout << instance_array.size() << endl;
            instance_array.resize(GRID_ROWS, vector<vector<vector<set<int>>>>(GRID_COLS, vector<vector<set<int>>>(1, vector<set<int>>(2))));
            // features[k]-1 returns the feature id -1 value of the kth feature
            getCMBRLayerWCount2(k, k+1, 0, false, temp_comb, 2); 
            fc++;

        } else {
            // find all the CMBRs from 1st feature to K+1 feature 
            auto start_i = high_resolution_clock::now();
            for (int i = 0; i <= k; ++i)
            {  
                //cout << "Feature: degree2 ->" << i << " " << k+1 << endl;

                // get combination bit pattern
                temp_comb.reset();
                temp_comb[FMAX-1-i] = 1;
                temp_comb[FMAX-2-k] = 1;
                cout << temp_comb << endl;

                if (i == 0)
                {
                    fc += (k+1); 
                    // resize instance array for the new feature combination
                    // instance_array.resize(instance_array.size()+1, vector<set<int>>(2));
                    instance_array.resize(GRID_ROWS, vector<vector<vector<set<int>>>>(GRID_COLS, vector<vector<set<int>>>(fc, vector<set<int>>(12))));      //*******change manual 12**********                                 
                }                

                getCMBRLayerWCount2(i, k+1, 0, false, temp_comb, 2); 
            //}

            // find CMBRs with K+1 feature and previous layer CMBRs
            //for (int i = 0; i < k ; ++i)
            //{   
                auto start_jj = high_resolution_clock::now();  
                for(int jj=0; jj< cmbr_map[0][0][i].size() && i<k; ++jj)
                {
                    //cout << "Layer: " << i << " loc: " << jj << endl;

                    // temp = getCMBRLayerWCount(cmbr_arr[i][jj].cmbr_array, mbrs[k+1]); 
                    if (!cmbr_map[0][0][i][jj].isDeleted)
                    { 
                        // get combination bit pattern
                        temp_comb.reset();
                        temp_comb = cmbr_map[0][0][i][jj].combination; // take combintion id from previous step
                        temp_comb[FMAX-2-k] = 1; 
                        cout << temp_comb << endl;

                        // if (jj == 0)
                        // {
                        //     cout << "feature count " << fc << endl;
                        //     // resize instance array for the new feature combination
                        //     instance_array.resize(GRID_ROWS, vector<vector<vector<set<int>>>>(GRID_COLS, vector<vector<set<int>>>(fc, vector<set<int>>(cmbr_map[0][0][i][jj].featureCount+1))));
                            //cout << "c: " << instance_array[0][0].size()<< " " << instance_array[0][0][0].size() << " " << instance_array[0][0][0][0].size() << endl;
                        // }      
                        
                        getCMBRLayerWCount2(jj, k+1, i, true, temp_comb, cmbr_map[0][0][i][jj].featureCount+1); 
                        fc++;
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
        
        cmbr_filter_layerwise(k); 

        // clear temporary set structure
        instance_array.clear();

        //cout <<"Layer " << k << " Built Successfully!!!" << endl;       
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start); 
    print_time("Function: buildCMBRList " + to_string(duration.count()));
    return ;
} 

int main()
{
    
    // freopen ("out_v_4.2_1676.txt","w",stdout);   

    // read data into a table_row structure type 1D array
    struct table_row *dat;
    dat = createArray("data/newData/Seattle2012_1676.csv");
    // dat = createArray("data/Point_Of_Interest_modified.csv");

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

    cout << "mbr array constructed" << endl;
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
    cout << "cmbr layers constructed" << endl;

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
    

    // fclose(stdout);

    return 0;
}