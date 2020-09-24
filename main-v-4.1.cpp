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
#include <bits/stdc++.h> 

// distance threshold
float const DIST = 5.0;

// number features
int const FMAX = 13;
// int const FMAX = 5;

#include "readFille.cpp"

using namespace std::chrono;


// prevalence threshold
double const PI = 0.3;
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
    bool isDeleted = false;
    vector<mbr> cmbr_array;
    vector<vector<int>> list1;
    vector<vector<int>> list2;
};

// saves all mbr information
vector<vector<vector<vector<mbr>>>>  mbr_array(801, vector<vector<vector<mbr>>>(601, vector<vector<mbr>>(FMAX)));

// 2D vector to keep track of all the combinations and instances realted
vector<vector<vector<vector<cmbr>>>>  cmbr_map(801, vector<vector<vector<cmbr>>>(601, vector<vector<cmbr>>(FMAX-1)));

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
    cout << "Time -> " << str << endl;
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
vector<vector<int>> instanceCombinationBuild(vector<vector<int>> list1, vector<vector<int>> map_l1, vector<vector<int>> map_l2) {
    vector<int> ttlist1;
    vector<vector<int>> ttlist2;
    int tt=0;
    auto start = high_resolution_clock::now();
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
void getCMBRLayerWCount2(int fid1, int fid2, int crow, bool cmbrFlag, bitset<FMAX> comb) {

    // 1D array to hold layer CMBRs
    vector<mbr> arr;
    // temporary cmbr reference
    cmbr ret;

    // temporary varibale to track of calculated CMBR
    mbr cmbr_v;

    // l1 and l2 maintains the lists from mbr1 and mbr2. l1[0] l2[0] will give 0th CMBR instances
    vector<vector<int>> l1;
    vector<vector<int>> l2;

    // temporary arrays to hel build l1 and l2 rows
    vector<int> t1;
    vector<int> t2;

    int c, r, nrow, ncol, s;

    vector<mbr> box1;

    int init_layer_count = cmbr_map[0][0][fid2-1].size();    

    auto start = high_resolution_clock::now();

    // traverse grid
    for (int row = 0; row < GRID_ROWS; ++row)
    {
        for (int col = 0; col < GRID_COLS; ++col)
        {
            // cout << "r " << row << " c " << col << endl;
            // change MBR or CMBR for first box when checking for CMBR
            if (cmbrFlag)
            {
                // cout << "fid1 " << fid1 << endl; 
                // cout << "id1 " << id1 << endl; 
                // cout << "id1 " << id1 << " fid1 " << fid1 << " new fid1=" << cmbr_map[row][col][crow][id1].cmbr_array.size() << endl;                
                box1 = cmbr_map[row][col][crow][fid1].cmbr_array;
            } else {
                box1 = mbr_array[row][col][fid1];
            }
            for (int ii = 0; ii < 4; ++ii)
            {
                nrow = row + ii/2;
                ncol = col + ii%2;
                if(nrow >= GRID_ROWS || ncol >= GRID_COLS) 
                {
                    continue;
                }

                for (int i = 0; i < box1.size(); ++i)
                {
                    // cout << "i " << i << endl;
                    
                    for (int j = 0; j < mbr_array[nrow][ncol][fid2].size(); ++j)
                    {
                        // cout << "j " << j << endl;

                        cmbr_v = calculateCMBR(box1[i].x1, box1[i].y1, box1[i].x2, box1[i].y2, mbr_array[nrow][ncol][fid2][j].x1, mbr_array[nrow][ncol][fid2][j].y1, mbr_array[nrow][ncol][fid2][j].x2, mbr_array[nrow][ncol][fid2][j].y2);
                        if (!cmbr_v.empty)
                        {                              
                            c = floor((cmbr_v.x1 - GRID_MIN_X) / (DIST * 2));
                            r = GRID_ROWS - 1 - floor((cmbr_v.y2 - GRID_MIN_Y) / (DIST * 2));

                            arr.push_back(cmbr_v); 
                            t2.push_back(j);                            
                        } 
                    }
                    // if there are any CMBRs for i, add ID i to l1
                    if (t2.size() > 0)
                    {
                        t1.push_back(i);
                        l1.push_back(t1); // push 1D arrays to 2D array
                        l2.push_back(t2); // push 1D arrays to 2D array
                        t1.clear(); //clear 1D array
                        t2.clear(); //clear 1D array                
                    }
                }
                ret.combination = comb;  
                s = cmbr_map[nrow][ncol][fid2-1].size();                
                if (arr.size() > 0)
                {
                    // check if the CMBR pattern already exists in the selected cell
                    if (init_layer_count < s)
                    {
                        cmbr_map[nrow][ncol][fid2-1][s-1].cmbr_array.insert(cmbr_map[nrow][ncol][fid2-1][s-1].cmbr_array.end(), arr.begin(), arr.end());
                        cmbr_map[nrow][ncol][fid2-1][s-1].list1.insert(cmbr_map[nrow][ncol][fid2-1][s-1].list1.end(), l1.begin(), l1.end());
                        cmbr_map[nrow][ncol][fid2-1][s-1].list1.insert(cmbr_map[nrow][ncol][fid2-1][s-1].list1.end(), l2.begin(), l2.end());
                    } else {
                        ret.cmbr_array = arr;
                        ret.list1 = l1;
                        // check if matching is CMBR vs MBR                        
                        if (cmbrFlag)
                        {
                            ret.list2 = instanceCombinationBuild(l2, cmbr_map[row][col][crow][fid1].list1, cmbr_map[row][col][crow][fid1].list2);                            
                        } else {
                            ret.list2 = l2;
                        }
                    }
                    l1.clear();
                    l2.clear();
                    arr.clear();                        
                }
                // insert new feature combination if it does not exist in the cell 
                if (init_layer_count >= s)
                {
                    cmbr_map[nrow][ncol][fid2-1].push_back(ret);
                }            
                ret = {}; //reset ret
            }           
        }
    }
}

//erasing filtered out cmbr_maps layerwise
// void erase_cmbr_map(int k, vector<int> erase_list)
// {
//     //before erase
//     //cout << "CMBR MAP before erase " << endl;
//     //for(int i = 0; i<cmbr_map[k].size(); i++){
//         //cout << cmbr_map[k][i].combination << " ";  
//     //}
//     //cout << endl;   
//     // erase code last to first
//     auto start = high_resolution_clock::now();
//     for(int ii=erase_list.size()-1; ii >= 0 ; ii--)
//     {
//         // cout << "erasing " << cmbr_map[k][erase_list[ii]].combination << " index is " << erase_list[ii] << " k is " << k << endl;
//         // cmbr_map[k].erase(cmbr_map[k].begin() + erase_list[ii]);
//         // cmbr_arr[k].erase(cmbr_arr[k].begin() + erase_list[ii]);

//         cmbr_map[k][erase_list[ii]].isDeleted = true;

//         cmbr_map[k][erase_list[ii]].cmbr_array.clear();
//         cmbr_map[k][erase_list[ii]].list1.clear();
//         cmbr_map[k][erase_list[ii]].list2.clear();
//     }
//     auto stop = high_resolution_clock::now();
//     auto duration = duration_cast<microseconds>(stop - start); 
//     print_time("Function: erase_cmbr_map " + to_string(duration.count()));
//     //after erase
//     //cout << "CMBR MAP after erase " << endl;
//     //for(int i = 0; i<cmbr_map[k].size(); i++){
//         //cout << cmbr_map[k][i].combination << " ";  
//     //}
//     //cout << endl;
//     return;
// }

// //cmbr filter layerwise where k is the layer number
// void cmbr_filter_layerwise(int k)
// {
//     cout << "CMBR Filter Layer is " << k << " and size is " << cmbr_map[k].size() <<endl;
//     prev_size += cmbr_map[k].size();
//     vector<int> erase_list;

//     auto start = high_resolution_clock::now();
//     for(int ii = 0; ii<cmbr_map[k].size(); ii++)
//     {
//         if (!cmbr_map[k][ii].isDeleted)
//         {
            
//             string result = "";
//             double pr = 1.0;
//             double total_instances= 0;
//             bitset<FMAX> bit_comb = cmbr_map[k][ii].combination;
//             result += "(";
//             int flag =0;
//             for(int n = FMAX-1; n>=0; n--){
//                 if(bit_comb[n]==1 && flag==1){
//                     result += ", "; 
//                 }
//                 if(bit_comb[n]==1){
//                     result += to_string(FMAX-1-n);
//                     flag=1;
//                 }       
//             }
//             result += ") -> (";
//             //cout << result << endl;
//             // list_1_each_cell_item_count
//             int L1_item_count;
//             //cout << "ii is " << ii << " combination is " << bit_comb << " and count is " << cmbr_map[k][ii].count << endl;

//             vector<int> list_1;
//             for(int x=0; x<cmbr_map[k][ii].list1.size(); x++)
//             {
//                 L1_item_count = 0;
//                 int tmp = 0;
//                 for (int y=0; y<cmbr_map[k][ii].list1[x].size(); y++)
//                 {
//                     ////cout << " --- List 1 --- > " << cmbr_map[k][ii].list1[x][y] ;
//                     tmp = cmbr_map[k][ii].list1[x][y] ;
//                     L1_item_count++;

//                 }
//                 ////cout << endl;
//                 ////cout << "L1 item count is " << L1_item_count << endl;
//                 if (L1_item_count == 1)
//                 {
//                     list_1.push_back(tmp);
//                 }
//             }

//             if (L1_item_count == 1)
//             {
//                 sort(list_1.begin(), list_1.end());
//                 list_1.erase(unique(list_1.begin(), list_1.end()), list_1.end());
//                 //cout << "list 1 unique val " << list_1.size() << endl;
//                 int first_pos;
//                 for(first_pos = FMAX-1; first_pos >= 0 ; first_pos--)
//                 {
//                     if(bit_comb[first_pos]==1)
//                     break;  
//                 }
//                 total_instances = fcount[FMAX-1-first_pos];
//                 //cout << "test pr before if list 1 " << list_1.size()/total_instances << endl;
//                 result += to_string(list_1.size()) + ":" + to_string(fcount[FMAX-1-first_pos]) + ", ";
//                 if((list_1.size()/total_instances) < pr)
//                 {   
//                     pr = list_1.size()/total_instances;
//                 }
//             }

//             vector<vector<int>> list_1_2d(L1_item_count);

//             if(L1_item_count > 1)
//             {
//                 ////cout << "inside l1 count > 1" << endl;
//                 for(int j = 0; j<L1_item_count; j++)
//                 {
//                 ////cout << "j is -> " << j << endl;
//                     for(int z = 0; z<cmbr_map[k][ii].list1.size(); z++)
//                     {
//                         ////cout << endl << "z is -> " << z << endl;
//                         list_1_2d[j].push_back(cmbr_map[k][ii].list1[z][j]);
//                         ////cout << "  ^^^   " << endl;
//                         ////cout << "  " << cmbr_map[k][ii].list1[z][j] << "  " ;
//                     }
//                     ////cout << endl;
//                 }
//                 int pos = FMAX;
//                 for(int j = 0; j<L1_item_count; j++)
//                 {
//                     sort(list_1_2d[j].begin(), list_1_2d[j].end());
//                     list_1_2d[j].erase(unique(list_1_2d[j].begin(), list_1_2d[j].end()), list_1_2d[j].end());
//                     //cout << "list 1 2D unique val " << list_1_2d[j].size() << endl;

//                     for(int p = pos-1 ; p >= 0 ; p--)
//                     {
//                         if(bit_comb[p]==1)
//                         {
//                             pos = p;
//                             break;
//                         }   
//                     }
//                     //{4, 4, 2, 3, 1}
//                     total_instances = fcount[FMAX-1-pos];
//                     //cout << "total_instances " << total_instances << " pos " << pos << endl;
//                     //cout << "test pr before if list 1 2d " << list_1_2d[j].size()/total_instances << endl;
//                     result += to_string(list_1_2d[j].size()) + ":" + to_string(fcount[FMAX-1-pos]) + ", ";
//                     if(((list_1_2d[j].size())/total_instances) < pr)
//                     {
//                         pr = (list_1_2d[j].size())/total_instances;
//                     }
//                 }
//             }
            
//             //cout << "List 1 pr -> " << pr << endl;      
            
//             vector<int> list_2;
//             for(int x=0; x<cmbr_map[k][ii].list2.size(); x++)
//             {
//                 for (int y=0; y<cmbr_map[k][ii].list2[x].size(); y++)
//                 {
//                     list_2.push_back(cmbr_map[k][ii].list2[x][y]);
//                     ////cout << " --- List 2 --- > " << cmbr_map[k][ii].list2[x][y];
//                 }    
//                 ////cout << endl;
//             }
//             ////cout << endl;
//             sort(list_2.begin(), list_2.end());
//             list_2.erase(unique(list_2.begin(), list_2.end()), list_2.end());
//             //cout << "list 2 unique val " << list_2.size() << endl;

//             int last_pos;
//             for(last_pos = 0; last_pos < FMAX ; last_pos++)
//             {
//                 if(bit_comb[last_pos]==1)
//                     break;  
//             }
//             total_instances = fcount[FMAX-1-last_pos];
//             //cout << "test pr before if list 2 " << list_2.size()/total_instances << endl;
//             result += to_string(list_2.size()) + ":" + to_string(fcount[FMAX-1-last_pos]) + ") -> " + to_string(cmbr_map[k][ii].cmbr_array.size());
//             if((list_2.size()/total_instances) < pr)
//             {
//                 pr = list_2.size()/total_instances;
//             }       
//             //cout << "List 2 -> pr " << pr << endl;

//             if (pr < PI)
//             {
//                 //cout << "to be removed comb is " << cmbr_map[k][ii].combination << " index is " << ii<< " and K is " << k << endl;
//                 //cout << "pr is " << pr << " PI is " << PI << endl;
//                 erase_list.push_back(ii);
//                 pr = 1.0;
//                 result = "";
//             }
//             else
//             {
//                 //cout << "survivor cmbr is --->  " << cmbr_map[k][ii].combination << endl;
//                 //Debug format for comparison
//                 cout << result << endl;
//             }
//         }
//     }//end ii
//     auto stop = high_resolution_clock::now();
//     auto duration = duration_cast<microseconds>(stop - start); 
//     print_time("Function: cmbr_filter_layerwise " + to_string(duration.count()));

//     erase_cmbr_map(k, erase_list);
//     erase_list.clear();
//     return;
// }

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

    bitset<FMAX> temp_comb; 


    auto start = high_resolution_clock::now();
    for (int k = 0; k < layers; ++k)
    {           
        //cout <<"Layer " << k << " Building ..." << endl;
        //b = features[k+1]-1;

        // layer 1. Instance wise CMBR only.No previos layer
        if (k == 0) {
            //a = features[k]-1;      

            cout << "Feature: " << k << " with feature: " << k+1 << endl;
            //cout << "ptr[a]= " << fcount[k] << " ptr[b]= " << fcount[k+1] << endl;

            // get combination bit pattern
            temp_comb.reset();               
            temp_comb[FMAX-1-k] = 1;
            temp_comb[FMAX-2-k] = 1;
            // features[k]-1 returns the feature id -1 value of the kth feature
            getCMBRLayerWCount2(k, k+1, 0, false, temp_comb);            

        } else {
            // find all the CMBRs from 1st feature to K+1 feature 
            auto start_i = high_resolution_clock::now();
            for (int i = 0; i <= k; ++i)
            {  
                cout << "Feature: degree2 ->" << i << " " << k+1 << endl;

                // get combination bit pattern
                temp_comb.reset();
                temp_comb[FMAX-1-i] = 1;
                temp_comb[FMAX-2-k] = 1;

                getCMBRLayerWCount2(i, k+1, 0, false, temp_comb); 
            //}

            // find CMBRs with K+1 feature and previous layer CMBRs
            //for (int i = 0; i < k ; ++i)
            //{   
                auto start_jj = high_resolution_clock::now();  
                for(int jj=0; jj< cmbr_map[0][0][i].size() && i<k; ++jj)
                {
                    cout << "Layer: " << i << " loc: " << jj << endl;

                    // temp = getCMBRLayerWCount(cmbr_arr[i][jj].cmbr_array, mbrs[k+1]); 
                    if (!cmbr_map[0][0][i][jj].isDeleted)
                    { 
                        // get combination bit pattern
                        temp_comb.reset();
                        temp_comb = cmbr_map[0][0][i][jj].combination; // take combintion id from previous step
                        temp_comb[FMAX-2-k] = 1; 

                        getCMBRLayerWCount2(jj, k+1, i, true, temp_comb);        
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
        
        // cmbr_filter_layerwise(k); 

        //cout <<"Layer " << k << " Built Successfully!!!" << endl;       
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start); 
    print_time("Function: buildCMBRList " + to_string(duration.count()));
    return ;
} 

int main()
{
    
    // freopen ("out.txt","w",stdout);   

    // read data into a table_row structure type 1D array
    struct table_row *dat;
    dat = createArray("data/newData/Seattle2012_tt_1.csv");
    // dat = createArray("data/Point_Of_Interest_modified.csv");

    // test ----- START ----

    // number of rows in the data file
    // ROWS = 14;

    // // grid origin
    // GRID_MIN_X = 0.0, GRID_MIN_Y = 0.0;

    // // grid number of rows
    // GRID_ROWS = ceil((1010 - GRID_MIN_X)/DIST * 2);

    // // grid number of columns
    // GRID_COLS = ceil((1010 - GRID_MIN_Y)/DIST * 2);


    // struct table_row dat[14] = {{1, 500,500}, {1, 1005,1005}, {1, 825, 325}, {1, 130, 200},
    //                                  {2, 505, 500}, {2, 506, 500}, {2, 830, 250}, {2, 101, 101},
    //                                  {3, 1010, 1010}, {3, 515, 515},
    //                                  {4, 1005, 1005}, {4, 135, 205}, {4, 509, 506},
    //                                  {5, 400, 400}};



    // test ----- END -----

    // calculate MBR for all the datapoints. 
    // returns a 2D array. 1st-D : Features, 2nd-D: instances per each feature 
    getMBRList(dat); 

    cout << "mbr array constructed" << endl;

    // print mbr array
    // for (int r = 0; r < GRID_ROWS; ++r)
    // {
    //     for (int c = 0; c < GRID_COLS; ++c)
    //     {
    //         for (int i = 0; i < mbr_array[r][c].size(); ++i)
    //         {
    //             if (mbr_array[r][c][i].size() > 0)
    //             {
    //                 // cout << "Row: " << r << " Column: " << " Feature id: " << i << " Count: " << mbr_array[r][c][i].size() << endl;
    //                 cout << r << "," << c << "," << i << "," << mbr_array[r][c][i].size() << endl;                    
    //             }
    //         }
    //     }
    // }

    // build CMBR tree 
    buildCMBRList();
    cout << "cmbr layers constructed" << endl;

    // // print cmbr array
    // for (int r = 0; r < GRID_ROWS; ++r)
    // {
    //     for (int c = 0; c < GRID_COLS; ++c)
    //     {
    //         for (int i = 0; i < cmbr_map[r][c].size(); ++i)
    //         {
    //             for (int j = 0; j < cmbr_map[r][c][i].size(); ++j)
    //             {
    //                 if (!cmbr_map[r][c][i][j].isDeleted && cmbr_map[r][c][i][j].cmbr_array.size() > 0)
    //                 {
    //                     cout << cmbr_map[r][c][i][j].combination << "[" << cmbr_map[r][c][i][j].cmbr_array.size() << "] ";                
    //                 }
    //             }
    //             // cout << endl;
    //         }
    //     }
    // }
    

    // fclose(stdout);

    return 0;
}