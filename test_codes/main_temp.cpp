#include <iostream>
#include <deque>
#include <string>
#include <vector>
#include <array>
#include <bitset>

#include <boost/assign/std/vector.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

#include "readFille.cpp"

//#define INST 45000//10325 //max instance per feature
//small data
//#define INST 50
//#define FEATURES 43
//#define FMAX 13


// for testing
#define INST 10
#define FEATURES 5
#define FMAX 5

//#define PI 0.5
#define DIST 10

using namespace std;

typedef boost::geometry::model::d2::point_xy<float> point_xy;
typedef boost::geometry::model::polygon<point_xy > polygon;

// data structure hold cmbr info
struct cmbr {
	int count = 0;
	std::vector<polygon> cmbr_array;
	std::vector<std::vector<int>> list1;
	std::vector<std::vector<int>> list2;
};

// data structure to hold combination and count
struct cmbr_comb {
	std::bitset<FMAX> combination;
	int count = 0;
	std::vector<std::vector<int>> list1;
	std::vector<std::vector<int>> list2;
};

// 2D vector to keep track of all the combinations and counts
std::vector<std::vector<cmbr_comb>> cmbr_map(FMAX-1);

std::vector< std::vector<cmbr> > cmbr_arr(FMAX-1);

double const PI = 0.3;
int prev_size = 0;
//int fcount[FMAX] = {4060, 1404, 1899, 1367, 6617, 2579, 1695, 671, 1528, 940, 10325, 10315, 606};
// feature serial 1, 5, 8, 9, 10, 14, 20, 24, 28, 39, 40, 42, 43
//small data
//int fcount[FMAX] = {50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50};
int fcount[FMAX] = {5, 5, 2, 3, 1};

// returns MBR for a given data point
polygon getMBR(float px, float py) {
    using namespace boost::assign;    

    // return MBR
    polygon ret;
    // threshhold constant detrmines the MBR size for an instance
    float const d = DIST;

    // Create points to represent a rectagle.
    std::vector<point_xy> points;
    points +=
      point_xy(px-d,py-d),
      point_xy(px-d,py+d),
      point_xy(px+d,py+d),
      point_xy(px+d,py-d),
      point_xy(px-d,py-d);

    // Create a polygon object and assign the points to it.
    boost::geometry::assign_points(ret, points);

    return ret;
}

// returns MBR for a selected feature

polygon** getMBRList(struct table_row *data, int li, int *ptr) {
    // output 2D polygon. 
    polygon** arr = 0;

    // Creating 1st-D
    arr = new polygon *[FEATURES];
	//cout << " inside get mbr list " << endl;

    for (int i = 0; i < li; ++i)
    { 
         // only creates the 2nd dimension for the selected indices,
        // iff the  (feature id -1) equals to 1st-D index
        // helps to not assign 2nd-D to the unwanted indices
        if (ptr[data[i].id-1] == 0)
         {
            arr[data[i].id-1] = new polygon [INST];
			//cout << "id is " << data[i].id-1 << endl;
         }
		//cout << " inside for loop " << i << " ptr[data[i].id-1] = " << ptr[data[i].id-1] << endl;
        // calculate MBR using the getMBR() and assign it to the relavant feature instance
        arr[data[i].id-1][ptr[data[i].id-1]++] = getMBR(data[i].x, data[i].y);
    }

    return arr;
}

// returns CMBR for a given two MBRs
polygon getCMBR(polygon a, polygon b) {
    // Dequeue to save all the intersections from the input polygons
    // We may only have rectangles as input and output always wil be a one rectangle
    std::deque<polygon> output;
    // output rectangle
    polygon ret;

    // check intersection between a and b
    boost::geometry::intersection(a, b, output);
    // if an intersection exists, it will be assigned to the output. 
    // At max, we may have only one intersection for two given rectangles.
    if (!output.empty())
    {
        ret = output.front();        
    }

    return ret;
}

// read previous stel cmbr instance combination and save it to list1
std::vector<std::vector<int>> instanceCombinationBuild(std::vector<std::vector<int>> list1, std::vector<std::vector<int>> map_l1, std::vector<std::vector<int>> map_l2) {
    std::vector<int> ttlist1;
    std::vector<std::vector<int>> ttlist2;
    int tt=0;

    for (int aa = 0; aa < list1.size(); ++aa)
    {
        
        // goes through list2 of previos k-1 step to find the insatnce which made new CMBRs
        for (int bb = 0; bb < map_l2.size(); ++bb)
        {
            if(bb == 0) {
                // track the size of list2 first row
                tt = map_l2[bb].size();
            }else {         
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

// check if there can be a CMBR
bool isIntersection(float x1, float x2, float y1, float y2) {
        // std::cout << x1 << " " << x2 << std::endl;

    if ((abs(x1 - x2) < DIST*2) && (abs(y1 - y2) < DIST*2))
    {
        return true;
    }
    return false;
}

// returns a cmbr structure for selected 2 MBRs with the count of CMBRs
cmbr getCMBRLayerWCount(polygon *mbrs1,  polygon *mbrs2, int a, int b) {

    // 1D array to hold layer CMBRs
    std::vector<polygon> arr;

    // controls next available sapce in the layer array
    int insid = 0;

    // temporary varibale to track of calculated CMBR
    polygon cmbr_v;

    // l1 and l2 maintains the lists from mbr1 and mbr2. l1[0] l2[0] will give 0th CMBR instances
    std::vector<std::vector<int>> l1;
    std::vector<std::vector<int>> l2;

    // temporary arrays to hel build l1 and l2 rows
	std::vector<int> t1;
	std::vector<int> t2;

    //nested loop to check all the CMBRs for all combinations of instances 
    for (int i = 0; i < a; ++i)
    { 
        for (int j = 0; j < b; ++j)
        {
			auto poly1 = boost::begin(boost::geometry::exterior_ring(mbrs1[i]));
            auto poly2 = boost::begin(boost::geometry::exterior_ring(mbrs2[j]));
            if (isIntersection(boost::geometry::get<0>(*poly1), boost::geometry::get<0>(*poly2), boost::geometry::get<1>(*poly1), boost::geometry::get<1>(*poly2)))
            {
		        cmbr_v= getCMBR(mbrs1[i], mbrs2[j]);
		         
		        if (boost::geometry::num_points(cmbr_v) > 0)
		        {
		        	insid++;
		            arr.push_back(cmbr_v); 

		            t2.push_back(j);
		            // std::cout << i << ": " << j << " ";  
		        }
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

    std::cout << "count: " << ret.count << std::endl;
    std::cout << "Size: " << ret.cmbr_array.size() << std::endl;
    

    //return 2D array with CMBRs
    return ret;
}

//erasing filtered out cmbr_maps layerwise
void erase_cmbr_map(int k, vector<int> erase_list)
{
	// erase code last to first
	for(int ii=erase_list.size()-1; ii >= 0 ; ii--)
	{
		//cout << "erasing " << cmbr_map[k][erase_list[ii]].combination << " index is " << erase_list[ii] << " k is " << k << endl;
		cmbr_map[k].erase(cmbr_map[k].begin() + erase_list[ii]);
		cmbr_arr[k].erase(cmbr_arr[k].begin() + erase_list[ii]);
	}
	return;
}

//cmbr filter layerwise where k is the layer number
void cmbr_filter_layerwise(int k)
{
	//Iqra
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
			//cmbr_map[k].erase(cmbr_map[k].begin() + ii);
			//cmbr_arr[k].erase(cmbr_arr[k].begin() + ii);
			erase_list.push_back(ii);
			pr = 1.0;

		}
		else
		{
			cout << "survivor cmbr is --->  " << cmbr_map[k][ii].combination << endl;
		}
		//cout << "--------------------- " << ii << endl; 

	}//end ii
	erase_cmbr_map(k, erase_list);
	erase_list.clear();
	//iqra 
	return;
}


// returns a 2D vector of all the CMBRs of the features.  Maintains count for each CMBR
//std::vector<std::vector<cmbr>> 
void buildCMBRList(polygon **mbrs, int *ptr, int *features) {
    // numbers layers need to consider. (#features-1)
    int layers = FMAX-1;
    // 2D array to hold all CMBRS. outer D represents the layers and inner D 
    // represents CMBRs of that layer
    //std::vector< std::vector<cmbr> > cmbr_arr(layers); 

    // size of each layer. This can be dynamicaally calculated using ptr[a]*ptr[b]. 
    // But it can exceed memory alloation
    int x = 1;

    // temporary 1D array to hold the CMBRs of a layer
    cmbr temp;

    // varible to temporariy hold feature IDs
    int a, b;
    // varibale to temporarily hold combination data and the count of CMBRs for that combination
	struct cmbr_comb comb;

    for (int k = 0; k < layers; ++k)
    {  
        std::cout <<"Layer " << k << " Building ..." << std::endl;
        // define 2nd-D size 
        // arr[k] = std::vector<cmbr>(x); 

        b = features[k+1]-1;

        // layer 1. Instance wise CMBR only.No previos layer
        if (k == 0)	{
        	a = features[k]-1;

	        std::cout << "Feature: " << a+1 << " with feature: " << b+1 << std::endl;
			cout << "ptr[a]= " << ptr[a] << " ptr[b]= " << ptr[b] << endl;
            // features[k]-1 returns the feature id -1 value of the kth feature
            temp = getCMBRLayerWCount(mbrs[a], mbrs[b], ptr[a], ptr[b]);            
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
            // find all the CMBRs from kth feature to K+1 feature 
            for (int i = 0; i <= k; ++i)
            {    
        		a = features[i]-1;
	        	std::cout << "Feature: " << a+1 << " with feature: (K != 0) " << b+1 << std::endl;

                temp = getCMBRLayerWCount(mbrs[a], mbrs[b], ptr[a], ptr[b]); 
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

					//std::cout << FMAX-1-i << std::endl;

                	cmbr_arr[k].push_back(temp);                	
                }
                // arr[k].insert( arr[k].end(), temp.begin(), temp.end());
            }

            // find CMBRs with K+1 feature and previous layer CMBRs
            for (int i = 0; i < k-1; ++i)
            {    
				for(int jj=0; jj< cmbr_arr[i].size(); ++jj)
					{
			    	std::cout << "Layer: " << i << " loc: " << jj << " with feature: " << b+1 << std::endl;

		    		temp = getCMBRLayerWCount(&cmbr_arr[i][jj].cmbr_array[0], mbrs[b], cmbr_arr[i][jj].count, ptr[b]);        
		        	
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
        
		//Iqra
        cmbr_filter_layerwise(k);
        //iqra
        // arr[k].insert( arr[k].end(), temp.begin(), temp.end()); 
        std::cout <<"Layer " << k << " Built Successfully!!!" << std::endl;       
    }
	/*for (int k=layers-1 ; k>=0; k--)
	{
		cmbr_filter_layerwise(k);
	}*/
    return ; //cmbr_arr;
} 


int main()
{
   	freopen ("mbr_filter_temp.txt","w",stdout);   
    //array to hold the number of instances for each feature
    /*static int feature_sizes[FMAX] = {0};

    // feature id list
    static int feature_ids[FMAX] = {1, 5, 8, 9, 10, 14, 20, 24, 28, 39, 40, 42, 43};

    // read data into a table_row structure type 1D array
    struct table_row *dat;
    dat = createArray("Seattle2012_1.csv");
	//small data
	//dat = createArray("data_iqr.csv");

    // calculate MBR for all the datapoints. 
    // returns a 2D array. 1st-D : Features, 2nd-D: instances per each feature 
    polygon**  mbr_array = getMBRList(dat, ROWS, feature_sizes);   
	//cout << "mbr array constructed" << endl;

    // build CMBR tree 
    //std::vector< std::vector<cmbr> > cmbr_layers = 
	buildCMBRList(mbr_array, feature_sizes, feature_ids);
	//cout << "cmbr layers constructed" << endl;
    // print bitmap array
	*/
	// testing getMBR() START
    struct table_row test_dat[16] = {{1, 500,500}, {1, 700,700}, {1, 702,700}, {1, 825, 325}, {1, 130, 200},
                                    {2, 510, 500}, {2, 1000, 1000}, {2, 705, 700}, {2, 706, 700}, {2, 101, 101},
                                    {3, 100, 100}, {3, 800, 800},
                                    {4, 1005, 1005}, {4, 135, 205}, {4, 20, 20},
                                    {5, 703, 701}};
    
    static int test_feature_sizes[5] = {0};
    static int test_feature_ids[FMAX] = {1, 2, 3, 4, 5};

    polygon**  test_mbr_array = getMBRList(test_dat, 16, test_feature_sizes);  


    //std::vector< std::vector<cmbr> > test_cmbr_layers =
buildCMBRList(test_mbr_array, test_feature_sizes, test_feature_ids);


    int curr_size = 0;
    for (int i = 0; i < cmbr_map.size(); ++i)
    {
		cout << "CMBR layer is " << i << endl;
    	for (int j = 0; j < cmbr_map[i].size(); ++j)
    	{		
    		std::cout << cmbr_map[i][j].combination << "[" << cmbr_map[i][j].count << "]";
        	curr_size++;
    	}
    	std::cout << "\n";
    }
    cout << "Size before = " << prev_size << " Size after = " << curr_size << endl;
	fclose(stdout);
    return 0;
}
