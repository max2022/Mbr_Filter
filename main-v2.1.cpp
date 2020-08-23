#include <iostream>
#include <vector>
#include <array>
#include <bitset>

#include "readFille.cpp"

// number features
int const FMAX = 13;
// distance threshold
float const DIST = 5.0;
// prevalence threshold
double const PI = 0.3;
int prev_size = 0;
// saves all counts for features
vector<int> fcount(FMAX);

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
    int count = 0;
    vector<mbr> cmbr_array;
    vector<vector<int>> list1;
    vector<vector<int>> list2;
};

// data structure to hold combination and count
struct cmbr_comb {
    bitset<FMAX> combination;
    int count = 0;
    vector<vector<int>> list1;
    vector<vector<int>> list2;
};

// 2D vector to keep track of all the combinations and counts
vector<vector<cmbr_comb>> cmbr_map(FMAX-1);
// 2D vector for holding all the cmbr info
vector<vector<cmbr> > cmbr_arr(FMAX-1);

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
vector<vector<mbr>> getMBRList(struct table_row *data) {
    // output 2D mbr. 
    vector<vector<mbr>> arr(FMAX);
    int j = 0, k = 0;

    for (int i = 0; i < ROWS; ++i) { 
        // check if feature id changes
        if (j != (data[i].id - 1)) {
            j = (data[i].id - 1);
            k++;
        }

        // calculate MBR using the getMBR() and assign it to the relavant feature instance
        arr[k].push_back(getMBR(data[i].x, data[i].y));
        fcount[k] += 1; 
    }
    return arr;
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

// read previous stel cmbr instance combination and save it to list1
vector<vector<int>> instanceCombinationBuild(vector<vector<int>> list1, vector<vector<int>> map_l1, vector<vector<int>> map_l2) {
    vector<int> ttlist1;
    vector<vector<int>> ttlist2;
    int tt=0;

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
    return ttlist2;                 
}

// returns a cmbr structure for selected 2 MBRs with the count of CMBRs
cmbr getCMBRLayerWCount(vector<mbr> mbrs1,  vector<mbr> mbrs2) {

    // 1D array to hold layer CMBRs
    vector<mbr> arr;

    // controls next available sapce in the layer array
    int insid = 0;

    // temporary varibale to track of calculated CMBR
    mbr cmbr_v;

    // l1 and l2 maintains the lists from mbr1 and mbr2. l1[0] l2[0] will give 0th CMBR instances
    vector<vector<int>> l1;
    vector<vector<int>> l2;

    // temporary arrays to hel build l1 and l2 rows
    vector<int> t1;
    vector<int> t2;

    //nested loop to check all the CMBRs for all combinations of instances 
    for (int i = 0; i < mbrs1.size(); ++i)
    { 
        for (int j = 0; j < mbrs2.size(); ++j)
        {
            cmbr_v= calculateCMBR(mbrs1[i].x1, mbrs1[i].y1, mbrs1[i].x2, mbrs1[i].y2, mbrs2[j].x1, mbrs2[j].y1, mbrs2[j].x2, mbrs2[j].y2);
            if (!cmbr_v.empty)
            {
                insid++;
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

    // create return structure
    struct cmbr ret;
    if (insid > 0)
    {
        // ret((insid - 1), arr);
        ret.count = insid;
        ret.cmbr_array = arr;
        ret.list1 = l1;
        ret.list2 = l2;
        arr.clear();
        l1.clear();
        l2.clear();
    }
    cout << "count: " << ret.count << endl;
    cout << "Size: " << ret.cmbr_array.size() << endl;
    //return 2D array with CMBRs
    return ret;
}

//erasing filtered out cmbr_maps layerwise
void erase_cmbr_map(int k, vector<int> erase_list)
{
    //before erase
    cout << "CMBR MAP before erase " << endl;
    for(int i = 0; i<cmbr_map[k].size(); i++){
        cout << cmbr_map[k][i].combination << " ";  
    }
    cout << endl;   
    // erase code last to first
    for(int ii=erase_list.size()-1; ii >= 0 ; ii--)
    {
        //cout << "erasing " << cmbr_map[k][erase_list[ii]].combination << " index is " << erase_list[ii] << " k is " << k << endl;
        cmbr_map[k].erase(cmbr_map[k].begin() + erase_list[ii]);
        cmbr_arr[k].erase(cmbr_arr[k].begin() + erase_list[ii]);
    }

    //after erase
    cout << "CMBR MAP after erase " << endl;
    for(int i = 0; i<cmbr_map[k].size(); i++){
        cout << cmbr_map[k][i].combination << " ";  
    }
    cout << endl;
    return;
}

//cmbr filter layerwise where k is the layer number
void cmbr_filter_layerwise(int k)
{
    cout << "CMBR_MAP K is " << k << " size is " << cmbr_map[k].size() <<endl;
    prev_size += cmbr_map[k].size();
    vector<int> erase_list;
    for(int ii = 0; ii<cmbr_map[k].size(); ii++)
    {
        double pr = 1.0;
        double total_instances= 0;
        bitset<FMAX> bit_comb = cmbr_map[k][ii].combination;
        // list_1_each_cell_item_count
        int L1_item_count;
        cout << "ii is " << ii << " combination is " << bit_comb << " and count is " << cmbr_map[k][ii].count << endl;

        vector<int> list_1;
        for(int x=0; x<cmbr_map[k][ii].list1.size(); x++)
        {
            L1_item_count = 0;
            int tmp = 0;
            for (int y=0; y<cmbr_map[k][ii].list1[x].size(); y++)
            {
                //cout << " --- List 1 --- > " << cmbr_map[k][ii].list1[x][y] ;
                tmp = cmbr_map[k][ii].list1[x][y] ;
                L1_item_count++;

            }
            //cout << endl;
            //cout << "L1 item count is " << L1_item_count << endl;
            if (L1_item_count == 1)
            {
                list_1.push_back(tmp);
            }
        }

        if (L1_item_count == 1)
        {
            sort(list_1.begin(), list_1.end());
            list_1.erase(unique(list_1.begin(), list_1.end()), list_1.end());
            cout << "list 1 unique val " << list_1.size() << endl;
            int first_pos;
            for(first_pos = FMAX-1; first_pos >= 0 ; first_pos--)
            {
                if(bit_comb[first_pos]==1)
                break;  
            }
        
            total_instances = fcount[FMAX-1-first_pos];
            cout << "test pr before if list 1 " << list_1.size()/total_instances << endl;
            if((list_1.size()/total_instances) < pr)
            {   
                pr = list_1.size()/total_instances;
            }
        }

        vector<vector<int>> list_1_2d(L1_item_count);

        if(L1_item_count > 1)
        {
            //cout << "inside l1 count > 1" << endl;
            for(int j = 0; j<L1_item_count; j++)
            {
            //cout << "j is -> " << j << endl;
                for(int z = 0; z<cmbr_map[k][ii].list1.size(); z++)
                {
                    //cout << endl << "z is -> " << z << endl;
                    list_1_2d[j].push_back(cmbr_map[k][ii].list1[z][j]);
                    //cout << "  ^^^   " << endl;
                    //cout << "  " << cmbr_map[k][ii].list1[z][j] << "  " ;
                }
                //cout << endl;
            }
            int pos = FMAX;
            for(int j = 0; j<L1_item_count; j++)
            {
                sort(list_1_2d[j].begin(), list_1_2d[j].end());
                list_1_2d[j].erase(unique(list_1_2d[j].begin(), list_1_2d[j].end()), list_1_2d[j].end());
                cout << "list 1 2D unique val " << list_1_2d[j].size() << endl;

                for(int p = pos-1 ; p >= 0 ; p--)
                {
                    if(bit_comb[p]==1)
                    {
                        pos = p;
                        break;
                    }   
                }
                //{4, 4, 2, 3, 1}
                total_instances = fcount[FMAX-1-pos];
                cout << "total_instances " << total_instances << " pos " << pos << endl;
                cout << "test pr before if list 1 2d " << list_1_2d[j].size()/total_instances << endl;
                if(((list_1_2d[j].size())/total_instances) < pr)
                {
                    pr = (list_1_2d[j].size())/total_instances;
                }
            }
        }
        
        cout << "List 1 pr -> " << pr << endl;      
        
        vector<int> list_2;
        for(int x=0; x<cmbr_map[k][ii].list2.size(); x++)
        {
            for (int y=0; y<cmbr_map[k][ii].list2[x].size(); y++)
            {
                list_2.push_back(cmbr_map[k][ii].list2[x][y]);
                //cout << " --- List 2 --- > " << cmbr_map[k][ii].list2[x][y];
            }    
            //cout << endl;
        }
        //cout << endl;
        sort(list_2.begin(), list_2.end());
        list_2.erase(unique(list_2.begin(), list_2.end()), list_2.end());
        cout << "list 2 unique val " << list_2.size() << endl;

        int last_pos;
        for(last_pos = 0; last_pos < FMAX ; last_pos++)
        {
            if(bit_comb[last_pos]==1)
                break;  
        }
        total_instances = fcount[FMAX-1-last_pos];
        cout << "test pr before if list 2 " << list_2.size()/total_instances << endl;
        if((list_2.size()/total_instances) < pr)
        {
            pr = list_2.size()/total_instances;
        }       
        cout << "List 2 -> pr " << pr << endl;

        if (pr < PI)
        {
            cout << " to be removed comb is " << cmbr_map[k][ii].combination << " index is " << ii<< " and K is " << k << endl;
            cout << "pr is " << pr << " PI is " << PI << endl;
            erase_list.push_back(ii);
            pr = 1.0;
        }
        else
        {
            cout << "survivor cmbr is --->  " << cmbr_map[k][ii].combination << endl;
        }
    }//end ii
    erase_cmbr_map(k, erase_list);
    erase_list.clear();
    return;
}

// returns a 2D vector of all the CMBRs of the features.  Maintains count for each CMBR
void buildCMBRList(vector<vector<mbr>> mbrs, int *features) {
    // numbers layers need to consider. (#features-1)
    int layers = FMAX-1; 

    // temporary 1D array to hold the CMBRs of a layer
    cmbr temp;

    // varible to temporariy hold feature IDs
    int a, b;
    // varibale to temporarily hold combination data and the count of CMBRs for that combination
    struct cmbr_comb comb;

    for (int k = 0; k < layers; ++k)
    {  
        cout <<"Layer " << k << " Building ..." << endl;
        b = features[k+1]-1;

        // layer 1. Instance wise CMBR only.No previos layer
        if (k == 0) {
            a = features[k]-1;

            cout << "Feature: " << a+1 << " with feature: " << b+1 << endl;
            cout << "ptr[a]= " << fcount[k] << " ptr[b]= " << fcount[k+1] << endl;
            // features[k]-1 returns the feature id -1 value of the kth feature
            temp = getCMBRLayerWCount(mbrs[k], mbrs[k+1]);            
            // append all the calculated CMBRs to the layer 1 if CMBRs exists
            if (temp.count > 0) 
            {
                // update combinations global array
                comb.combination.reset();               
                comb.combination[FMAX-1-k] = 1;
                comb.combination[FMAX-2-k] = 1;
                comb.count = temp.count;
                comb.list1 = temp.list1;
                comb.list2 = temp.list2;
                cmbr_map[k].push_back(comb); 
                cmbr_arr[k].push_back(temp);
            }
        } else {
            // find all the CMBRs from 1st feature to K+1 feature 
            for (int i = 0; i <= k; ++i)
            {    
                a = features[i]-1;
                cout << "Feature: " << a+1 << " with feature: (K != 0) " << b+1 << endl;

                temp = getCMBRLayerWCount(mbrs[k], mbrs[k+1]); 
                // if CMBRs exists, add to the layer
                if (temp.count > 0)
                {
                    // update combinations global array
                    comb.combination.reset();
                    comb.combination[FMAX-1-i] = 1;
                    comb.combination[FMAX-2-k] = 1;
                    comb.count = temp.count;
                    comb.list1 = temp.list1;
                    comb.list2 = temp.list2;
                    cmbr_map[k].push_back(comb); 
                    cmbr_arr[k].push_back(temp);                    
                }
                // arr[k].insert( arr[k].end(), temp.begin(), temp.end());
            //}

            // find CMBRs with K+1 feature and previous layer CMBRs
            //for (int i = 0; i < k ; ++i)
            //{    
                for(int jj=0; jj< cmbr_arr[i].size() && i<k ; ++jj)
                    {
                    cout << "Layer: " << i << " loc: " << jj << " with feature: " << b+1 << endl;
                    temp = getCMBRLayerWCount(cmbr_arr[i][jj].cmbr_array, mbrs[k+1]);        
                    
                    // if CMBRs exists, add to the layer
                    if (temp.count > 0)
                    {
                        // update combinations global array
                        // change the required bit related to featured id into 1
                        comb.combination.reset();
                        comb.combination = cmbr_map[i][jj].combination; // take combintion id from previous step
                        comb.combination[FMAX-2-k] = 1;
                        
                        comb.count = temp.count;    

                        // add created CMBR instance list to object
                        comb.list1 =  instanceCombinationBuild(temp.list1, cmbr_map[i][jj].list1, cmbr_map[i][jj].list2);
                        // add list2 returned from the method
                        comb.list2 = temp.list2;
                        // push created feature combination to output array
                        cmbr_map[k].push_back(comb); 
                        // push created CMBR list and other info to CMBR output array 
                        cmbr_arr[k].push_back(temp);  
                    }               
                }                               
            }
        }
        
        cmbr_filter_layerwise(k); 
        cout <<"Layer " << k << " Built Successfully!!!" << endl;       
    }
    return ;
} 

int main()
{
    freopen ("mbr_filter_v2.txt","w",stdout); 
    //for batch
    //freopen ("mbr_filter_v2.1.txt","w",stdout);   

    // feature id list
    static int feature_ids[FMAX] = {1, 5, 8, 9, 10, 14, 20, 24, 28, 39, 40, 42, 43};

    // read data into a table_row structure type 1D array
    struct table_row *dat;
    dat = createArray("Seattle2012_tt_1.csv");

    // calculate MBR for all the datapoints. 
    // returns a 2D array. 1st-D : Features, 2nd-D: instances per each feature 
    vector<vector<mbr>>  mbr_array = getMBRList(dat); 

    cout << "mbr array constructed" << endl;
   
    // build CMBR tree 
    buildCMBRList(mbr_array, feature_ids);
    cout << "cmbr layers constructed" << endl;

    // print bitmap array
    
    // testing getMBR() START
    /*struct table_row test_dat[14] = {{1, 500,500}, {1, 1005,1005}, {1, 825, 325}, {1, 130, 200},
                                     {2, 510, 500}, {2, 511, 500}, {2, 830, 250}, {2, 101, 101},
                                     {3, 1010, 1010}, {3, 515, 515},
                                     {4, 1005, 1005}, {4, 135, 205}, {4, 509, 506},
                                     {5, 400, 400}};
    
    static int test_feature_sizes[5] = {0};
    static int test_feature_ids[FMAX] = {1, 2, 3, 4, 5};

    mbr**  test_mbr_array = getMBRList(test_dat, 14, test_feature_sizes);  


    //vector< vector<cmbr> > test_cmbr_layers =
buildCMBRList(test_mbr_array, test_feature_sizes, test_feature_ids);
*/
    cout << "map: " << endl;
    int curr_size = 0;
    for (int i = 0; i < cmbr_map.size(); ++i)
    {
        cout << "CMBR layer is " << i << endl;
        for (int j = 0; j < cmbr_map[i].size(); ++j)
        {       
            cout << cmbr_map[i][j].combination << "[" << cmbr_map[i][j].count << "]";
            curr_size++;
        }
        cout << "\n";
    }
    cout << "Size before = " << prev_size << " Size after = " << curr_size << endl;
    fclose(stdout);

    return 0;
}
