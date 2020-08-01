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
#include <boost/geometry/geometries/adapted/boost_polygon.hpp>
#include <boost/geometry/geometries/linestring.hpp>

#include <boost/foreach.hpp>

#include "readFille.cpp"

#define INST 40000
#define FEATURES 50
#define FMAX 13

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
	int count;
	std::vector<std::vector<int>> list1;
	std::vector<std::vector<int>> list2;
};

// 2D vector to keep track of all the combinations and counts
std::vector<std::vector<cmbr_comb>> cmbr_map(12);

// returns MBR for a given data point
polygon getMBR(float px, float py) {
    using namespace boost::assign;    

    // return MBR
    polygon ret;
    // threshhold constant detrmines the MBR size for an instance
    float const d = 2.0;

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

    for (int i = 0; i < li; ++i)
    { 
         // only creates the 2nd dimension for the selected indices,
        // iff the  (feature id -1) equals to 1st-D index
        // helps to not assign 2nd-D to the unwanted indices
        if (ptr[data[i].id-1] == 0)
         {
            arr[data[i].id-1] = new polygon [INST];
         }

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

// returns a 1D vector of CMBRs for a selected layer 
std::vector<polygon> getCMBRLayer(polygon *mbrs1,  polygon *mbrs2, int a, int b) {
    // defining temporary size for the vector
    int ss = 5;

    // 1D array to hold layer CMBRs
    std::vector<polygon> arr(ss);

    // controls next available sapce in the layer array
    int insid = 0;

    // temporary varibale to track of calculated CMBR
    polygon cmbr;

    //nested loop to check all the CMBRs for all combinations of instances 
    for (int i = 0; i < a; ++i)
    { 
        for (int j = 0; j < b; ++j)
        {
            cmbr= getCMBR(mbrs1[i], mbrs2[j]);
            // std::cout << !boost::geometry::is_empty(cmbr) <<std::endl;              
            // checking if the CMBR exists
            if (!boost::geometry::is_empty(cmbr))    
            {
                arr[insid++] = cmbr; 
                std::cout << i << ": " << j << " ";  
            } 
            // memory safe condition
            // if we reach max memory for the layer combination allocation, 
            // it may stop assigning further CMBRS
            if (ss <= insid)
            {           
                break;
            }       
        }
        // memory safe condition
        // if we reach max memory for the layer combination allocation, 
        // it may stop assigning further CMBRS
        if (ss <= insid)
        {
            std::cout << "Error: CMBR array not enough!!!" <<std::endl;              

            break;
        }
    }

    //return 2D array with CMBRs
    return arr;
}

// returns a 2D vector of all the CMBRs of the features
std::vector< std::vector<polygon> > getCMBRList(polygon **mbrs, int *ptr, int *features) {
    // numbers layers need to consider. (#features-1)
    int layers = 12;
    // 2D array to hold all CMBRS. outer D represents the layers and inner D 
    // represents CMBRs of that layer
    std::vector< std::vector<polygon> > arr(layers); 


    // size of each layer. This can be dynamicaally calculated using ptr[a]*ptr[b]. 
    // But it can exceed memory alloation
    int x = 100;

    // temporary 1D array to hold the CMBRs of a layer
    std::vector<polygon> temp(x);

    for (int k = 0; k < layers; ++k)
    {  
        std::cout <<"Layer " << k << " Building ..." << std::endl;
        // define 2nd-D size 
        arr[k] = std::vector<polygon>(x);   

        // layer 1. Instance wise CMBR only.No previos layer
        if (k == 0)
        {
            // features[k]-1 returns the feature id -1 value of the kth feature
            temp = getCMBRLayer(mbrs[features[k]-1], mbrs[features[k+1]-1], ptr[features[k]-1], ptr[features[k+1]-1]);            
        } else {
            // find all the CMBRs from kth feature to K+1 feature 
            for (int i = 0; i <= k; ++i)
            {    
                temp = getCMBRLayer(mbrs[features[k]-1], mbrs[features[k+1]-1], ptr[features[k]-1], ptr[features[k+1]-1]); 
                arr[k].insert( arr[k].end(), temp.begin(), temp.end());
            }

            // find CMBRs with K+1 feature and previous layer CMBRs
            temp = getCMBRLayer(&arr[k-1][0], mbrs[features[k+1]-1], arr[k-1].size(), ptr[features[k+1]-1]);        
        }
        // append all the calculated CMBRs to the realted layer 
        arr[k].insert( arr[k].end(), temp.begin(), temp.end()); 
        std::cout <<"Layer " << k << " Built Successfully!!!" << std::endl;       
    }

    return arr;
} 

// returns a cmbr structure for selected 2 MBRs with the count of CMBRs
cmbr getCMBRLayerWCount(polygon *mbrs1,  polygon *mbrs2, int a, int b) {
    // defining temporary size for the vector
    // int ss = 150000;
    int ss = 10;

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
    // for (int i = 0; i < a; ++i)
    for (int i = 0; i < 2000; ++i)
    { 
        // for (int j = 0; j < b; ++j)
        for (int j = 0; j < 2000; ++j)
        {
            cmbr_v= getCMBR(mbrs1[i], mbrs2[j]);
            // std::cout << !boost::geometry::is_empty(cmbr_v) <<std::endl;              
            // checking if the CMBR exists
            if (!boost::geometry::is_empty(cmbr_v))    
            {
            	insid++;
                arr.push_back(cmbr_v); 

                // update l1 and l2
                t2.push_back(j);
                // std::cout << i << ": " << j << " ";  
            } 
            // memory safe condition
            // if we reach max memory for the layer combination allocation, 
            // it may stop assigning further CMBRS
            if (ss <= insid)
            {           
                break;
            }       
        }

        // if there are any CMBRs for i, add ID i to l1
        if (t2.size() > 0)
        {
        	t1.push_back(i);
        		
        }

        // push t1 and t2 into l1 and l2
        if (t1.size() > 0)
        {
        	l1.push_back(t1);
        	t1 = {};
        }
        if (t2.size() > 0)
        {
        	l2.push_back(t2);
        	t2 = {};
        }
        // memory safe condition
        // if we reach max memory for the layer combination allocation, 
        // it may stop assigning further CMBRS
        if (ss <= insid)
        {
            std::cout << "Error: CMBR array not enough!!!" <<std::endl;              

            break;
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
    }

    std::cout << "count: " << ret.count << std::endl;
    std::cout << "Size: " << ret.cmbr_array.size() << std::endl;
    

    //return 2D array with CMBRs
    return ret;
}

// returns a 2D vector of all the CMBRs of the features.  Maintains count for each CMBR
std::vector<std::vector<cmbr>> buildCMBRList(polygon **mbrs, int *ptr, int *features) {
    // numbers layers need to consider. (#features-1)
    int layers = 3;
    // 2D array to hold all CMBRS. outer D represents the layers and inner D 
    // represents CMBRs of that layer
    std::vector< std::vector<cmbr> > arr(layers); 



    // size of each layer. This can be dynamicaally calculated using ptr[a]*ptr[b]. 
    // But it can exceed memory alloation
    int x = 1;

    // temporary 1D array to hold the CMBRs of a layer
    cmbr temp;

    // varible to temporariy hold feature IDs
    int a, b;
    // varibale to temporarily hold combination data and the count of CMBRs for that combination
	struct cmbr_comb comb;

	int tt=0;

	std::vector<int> ttlist1;
	std::vector<std::vector<int>> ttlist2;


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

	        	arr[k].push_back(temp);
	        }
        } else {
            // find all the CMBRs from kth feature to K+1 feature 
            for (int i = 0; i <= k; ++i)
            {    
        		a = features[i]-1;
	        	std::cout << "Feature: " << a+1 << " with feature: " << b+1 << std::endl;

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

					std::cout << FMAX-1-i << std::endl;

                	arr[k].push_back(temp);                	
                }
                // arr[k].insert( arr[k].end(), temp.begin(), temp.end());
            }

            // find CMBRs with K+1 feature and previous layer CMBRs
            for (int i = 0; i < arr[k-1].size(); ++i)
            {    
	        	std::cout << "Layer: " << k-1 << " loc: " << i << " with feature: " << b+1 << std::endl;

        		temp = getCMBRLayerWCount(&arr[k-1][i].cmbr_array[0], mbrs[b], arr[k-1][i].count, ptr[b]);        
            	// if CMBRs exists, add to the layer
                if (temp.count > 0)
                {
                	// update combinations global array
	        		comb.combination.reset();
                	comb.combination = cmbr_map[k-1][i].combination; // take combintion id from previous step
					comb.combination[FMAX-2-k] = 1;
					comb.count = temp.count;	


					for (int a = 0; a < temp.list1.size(); ++a)
					{
						tt = cmbr_map[k-1][i].list2[0].size();

						for (int b = 1; b < cmbr_map[k-1][i].list2[b].size(); ++b)
						{
							if ( tt < temp.list1[a][0])
							{
								// assign
								ttlist1.insert(ttlist1.end(), cmbr_map[k-1][i].list1[b].begin(), cmbr_map[k-1][i].list1[b].end());
								tt = tt + cmbr_map[k-1][i].list2[b].size()-1 - temp.list1[a][0];
								// tt=0;
								ttlist1.push_back(cmbr_map[k-1][i].list2[b][tt]);

								ttlist2.push_back(ttlist1);
								ttlist1.clear();

								break;
							}
							tt += cmbr_map[k-1][i].list2[b].size();					
						}
					}				

					comb.list1 = ttlist2;
					ttlist1.clear();
					comb.list2 = temp.list2;
					cmbr_map[k].push_back(comb); 
						
                	arr[k].push_back(temp);                	
                }           		           	
            }
        }
        
        // arr[k].insert( arr[k].end(), temp.begin(), temp.end()); 
        std::cout <<"Layer " << k << " Built Successfully!!!" << std::endl;       
    }

    return arr;
} 

int main()
{
    //array to hold the number of instances for each feature
    static int feature_sizes[FEATURES] = {0};

    // feature id list
    static int feature_ids[13] = {1, 5, 8, 9, 10, 14, 20, 24, 28, 39, 40, 42, 43};

    // read data into a table_row structure type 1D array
    struct table_row *dat;
    dat = createArray("Seattle2012.csv");

    // calculate MBR for all the datapoints. 
    // returns a 2D array. 1st-D : Features, 2nd-D: instances per each feature 
    polygon**  mbr_array = getMBRList(dat, ROWS, feature_sizes);   

    // calculate CMBR list 
    // std::vector< std::vector<polygon> > cmbr_array = getCMBRList(mbr_array, feature_sizes, feature_ids);

    // build CMBR tree 
    std::vector< std::vector<cmbr> > cmbr_layers = buildCMBRList(mbr_array, feature_sizes, feature_ids);

    // print bitmap array
    for (int i = 0; i < cmbr_map.size(); ++i)
    {
    	for (int j = 0; j < cmbr_map[i].size(); ++j)
    	{		
    		std::cout << cmbr_map[i][j].combination << "[" << cmbr_map[i][j].count << "]";
    	}
    	std::cout << "\n";
    }

    // prints the feature sizes 
    // for (int i = 0; i < FEATURES; ++i)
    // {
    //     std::cout << feature_sizes[i] << std::endl;
    // }

    // testing getMBR()
    // polygon x = getMBR(5,5);
    // polygon y = getMBR(2,2);
    // polygon cmbr;
    // cmbr = getCMBR(x, y);
    // std::cout << boost::geometry::wkt(cmbr) << std::endl;

    

    // std::cout <<  << std::endl;

    // print all CMBRs
    // int i = 0;
    // std::cout << "green && blue:" << std::endl;
    // BOOST_FOREACH(polygon const& p, cmbrs)
    // {
    //     std::cout << i++ << ": " << boost::geometry::area(p) << std::endl;
    //     std::cout << i++ << ": " << boost::geometry::wkt(p) << std::endl;
    // }


    return 0;
}