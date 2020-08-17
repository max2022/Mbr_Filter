#include <iostream>
#include <fstream>
#include <deque>
#include <string>
#include <vector>
#include <array>
#include <bitset>


#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

#include <boost/assign/std/vector.hpp>

// #include <boost/geometry/geometries/adapted/boost_polygon.hpp>
// #include <boost/geometry/geometries/linestring.hpp>
// #include <boost/type_traits/is_empty.hpp>

// #include <boost/foreach.hpp>

#include "readFille.cpp"

#define INST 45000
#define FEATURES 50
#define FMAX 13
#define DIST 10

// for testing
// #define INST 10
// #define FEATURES 5
// #define FMAX 5

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
std::vector<std::vector<cmbr_comb>> cmbr_map(FMAX - 1);

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
    polygon tmp;

    // writing data to file. Visualize purpose
    // std::ofstream myfile;
    // myfile.open ("MBRAll.csv");

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
        tmp = getMBR(data[i].x, data[i].y);
        arr[data[i].id-1][ptr[data[i].id-1]++] = tmp;
        // myfile << data[i].id-1 << "," << boost::geometry::dsv(tmp) <<"\n";
    }
    // myfile.close();


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
cmbr getCMBRLayerWCount(polygon *mbrs1,  polygon *mbrs2, int a, int b, int kk) {
    // defining temporary size for the vector
    // int ss = 150000;
    // int ss = 11000;

    // selects set of limited instances
    // int max_rows = 25000000;
    // int max_rows = 1000;

    // 1D array to hold layer CMBRs
    std::vector<polygon> arr;

    // controls next available sapce in the layer array
    int insid = 0;

    // temporary varibale to track of calculated CMBR
    polygon cmbr_v;

    // writing data to file. Visualize purpose
    // std::ofstream myfile;
    // myfile.open ("CMBRAll.csv", std::ios_base::app);

    // l1 and l2 maintains the lists from mbr1 and mbr2. l1[0] l2[0] will give 0th CMBR instances
    std::vector<std::vector<int>> l1;
    std::vector<std::vector<int>> l2;

    // temporary arrays to hel build l1 and l2 rows
    std::vector<int> t1;
    std::vector<int> t2;

    // for testing protection 
    // if (max_rows <= a)
    // {
    //     a = max_rows;
    // }
    // if (max_rows <= b)
    // {
    //    b = max_rows; 
    // }

    //nested loop to check all the CMBRs for all combinations of instances 
    for (int i = 0; i < a; ++i)
    // for (int i = 0; i < 50; ++i)
    { 
        for (int j = 0; j < b; ++j)
        // for (int j = 0; j < 50; ++j)
        {
            auto poly1 = boost::begin(boost::geometry::exterior_ring(mbrs1[i]));
            auto poly2 = boost::begin(boost::geometry::exterior_ring(mbrs2[j]));
            if (isIntersection(boost::geometry::get<0>(*poly1), boost::geometry::get<0>(*poly2), boost::geometry::get<1>(*poly1), boost::geometry::get<1>(*poly2)))
            {  
                // std::cout << i << ": " << j << std::endl;  
            
                cmbr_v= getCMBR(mbrs1[i], mbrs2[j]);
                // std::cout << boost::geometry::num_points(cmbr_v) <<std::endl; 
                // std::cout << boost::geometry::wkt(cmbr_v) <<std::endl;              

                // std::cout << !boost::geometry::is_empty(cmbr_v) <<std::endl;              
                // checking if the CMBR exists
                // if (!boost::geometry::is_empty(cmbr_v)) 
                if (boost::geometry::num_points(cmbr_v) > 0)
                {
                    insid++;
                    arr.push_back(cmbr_v); 

                    // update  l2
                    t2.push_back(j);
                    // std::cout << i << ": " << j << std::endl;  
                    // std::cout << boost::geometry::dsv(mbrs1[i]) <<std::endl;              
                    // std::cout << boost::geometry::dsv(mbrs2[j]) <<std::endl;              
                    // std::cout << boost::geometry::dsv(cmbr_v) <<std::endl;     
                    // myfile << kk << "," << i << "," << j << "," << boost::geometry::dsv(cmbr_v) <<"\n"; 

                } 
                // memory safe condition
                // if we reach max memory for the layer combination allocation, 
                // it may stop assigning further CMBRS
                // if (ss <= insid)
                // {           
                //     break;
                // }   
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
        // memory safe condition
        // if we reach max memory for the layer combination allocation, 
        // it may stop assigning further CMBRS
        // if (ss <= insid)
        // {
        //     std::cout << "Error: CMBR array not enough!!!" <<std::endl;              

        //     break;
        // }
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
        l1.clear();
        l2.clear();
    }

    std::cout << "count: " << ret.count << std::endl;
    std::cout << "Size: " << ret.cmbr_array.size() << std::endl;
    // myfile.close();
    

    //return 2D array with CMBRs
    return ret;
}

// temporary method to vector print 
void printToFile(std::vector<std::vector<int>> arr) {
    
    for(int i=0; i<arr.size(); ++i) {
        for (int j = 0; j < arr[i].size(); ++j)
        {
            std::cout << arr[i][j] << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << "-------------------" << std::endl;
}

void printVectorToFile(std::vector<std::vector<int>> arr, std::vector<std::vector<int>> arr2, std::bitset<FMAX> combination, std::ofstream& f) {
    for (int k = 0; k < arr.size(); ++k)
    {
        f << combination << "," << k << ",";
        for (int l = 0; l < arr[k].size(); ++l)
        {
            f << arr[k][l] <<",";   
        }
        f << "\n"; 
        f << combination << "," << k << ",";          
        for (int l = 0; l < arr2[k].size(); ++l)
        {
            f << arr2[k][l] <<", ";   
        }
        f << "\n";           
    }
}

// returns a 2D vector of all the CMBRs of the features.  Maintains count for each CMBR
std::vector<std::vector<cmbr>> buildCMBRList(polygon **mbrs, int *ptr, int *features) {
    // numbers layers need to consider. (#features-1)
    int layers = FMAX - 1;
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

    for (int k = 0; k < layers; ++k)
    {  
        std::cout <<"Layer " << k << " Building ..." << std::endl;
        // define 2nd-D size 
        // arr[k] = std::vector<cmbr>(x); 

        b = features[k+1]-1;

        // layer 1. Instance wise CMBR only.No previos layer
        if (k == 0) {
            a = features[k]-1;

            std::cout << "Feature: " << a+1 << " with feature: " << b+1 << std::endl;
            // features[k]-1 returns the feature id -1 value of the kth feature
            temp = getCMBRLayerWCount(mbrs[a], mbrs[b], ptr[a], ptr[b], k);            
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

                // std::cout << comb.combination << std::endl;
                // std::cout << "List 1" << std::endl;
                // printToFile(comb.list1);                    
                // std::cout << "List 2" << std::endl;
                // printToFile(comb.list2);                    

                arr[k].push_back(temp);
            }
        } else {
            // find all the CMBRs from kth feature to K+1 feature 
            for (int i = 0; i <= k; ++i)
            {    
                a = features[i]-1;
                std::cout << "Feature: " << a+1 << " with feature: " << b+1 << std::endl;

                temp = getCMBRLayerWCount(mbrs[a], mbrs[b], ptr[a], ptr[b], k); 
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

                    // std::cout << comb.combination << std::endl;
                    // std::cout << "List 1" << std::endl;
                    // printToFile(comb.list1);                    
                    // std::cout << "List 2" << std::endl;
                    // printToFile(comb.list2);

                    arr[k].push_back(temp);                 
                }
                // arr[k].insert( arr[k].end(), temp.begin(), temp.end());
            }

            // find CMBRs with K+1 feature and previous layer CMBRs
            for (int i = 0; i < arr[k-1].size(); ++i)
            {    
                std::cout << "Layer: " << k-1 << " loc: " << i << " with feature: " << b+1 << std::endl;

                temp = getCMBRLayerWCount(&arr[k-1][i].cmbr_array[0], mbrs[b], arr[k-1][i].count, ptr[b], k);        
                
                // if CMBRs exists, add to the layer
                if (temp.count > 0)
                {
                    // update combinations global array
                    // change the required bit related to featured id into 1
                    comb.combination.reset();

                    comb.combination = cmbr_map[k-1][i].combination; // take combintion id from previous step
                    
                    comb.combination[FMAX-2-k] = 1;
                    
                    comb.count = temp.count;    

                    // add created CMBR instance list to object
                    comb.list1 = instanceCombinationBuild(temp.list1, cmbr_map[k-1][i].list1, cmbr_map[k-1][i].list2);
                    // // clear temporary 2D array. ready for next feature combination from k-1 step
                    // ttlist1.clear();
     //                ttlist2.clear();
                    // add list2 returned from the method
                    comb.list2 = temp.list2;
                    // push created feature combination to output array
                    cmbr_map[k].push_back(comb); 
                    // push created CMBR list and other info to CMBR output array 
                    arr[k].push_back(temp);  

                    // std::cout << comb.combination << std::endl;
                    // std::cout << "List 1" << std::endl;
                    // printToFile(comb.list1);                    
                    // std::cout << "List 2" << std::endl;
                    // printToFile(comb.list2);                
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
    static int feature_sizes[43] = {0};

    // feature id list
    static int feature_ids[FMAX] = {1, 5, 8, 9, 10, 14, 20, 24, 28, 39, 40, 42, 43};

    // std::ofstream out("out.txt");
    // std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
    // std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!

    // read data into a table_row structure type 1D array
    struct table_row *dat;
    dat = createArray("Seattle2012_1.csv");

    // calculate MBR for all the datapoints. 
    // returns a 2D array. 1st-D : Features, 2nd-D: instances per each feature 
    polygon**  mbr_array = getMBRList(dat, ROWS, feature_sizes);   

    // for (int ii = 0; ii < 43; ++ii)
    // {
    //     std::cout << feature_sizes[ii] << std::endl;
    // }

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

    std::cout << "-----" << std::endl;

    // // writing data to file. Visualize purpose
    // std::ofstream myfile;
    // myfile.open ("CMBR_map.csv", std::ios_base::app);

    // testing getMBR() START
    // struct table_row test_dat[14] = {{1, 500,500}, {1, 700,700}, {1, 825, 325}, {1, 130, 200},
    //                                 {2, 510, 500}, {2, 1000, 1000}, {2, 830, 250}, {2, 101, 101},
    //                                 {3, 100, 100}, {3, 515, 515},
    //                                 {4, 1005, 1005}, {4, 135, 205}, {4, 509, 506},
    //                                 {5, 400, 400}};
    
    // static int test_feature_sizes[5] = {0};
    // static int test_feature_ids[FMAX] = {1, 2, 3, 4, 5};

    // polygon**  test_mbr_array = getMBRList(test_dat, 14, test_feature_sizes);  

    // std::cout << "MBR List: " <<std::endl;
    // for (int i = 0; i < 6; ++i)
    // {
    //     for (int j = 0; j < test_feature_sizes[i]; ++j)
    //     {
    //         std::cout << boost::geometry::dsv(test_mbr_array[i][j]) << " " <<std::endl;
    //     }
    //     std::cout << std::endl;
    // }


    // std::vector< std::vector<cmbr> > test_cmbr_layers = buildCMBRList(test_mbr_array, test_feature_sizes, test_feature_ids);

    // std::cout << "map: " << std::endl;

    // for (int i = 0; i < cmbr_map.size(); ++i)
    // {
    //     for (int j = 0; j < cmbr_map[i].size(); ++j)
    //     {       
    //         std::cout << cmbr_map[i][j].combination << "[" << cmbr_map[i][j].count << "]";
    //     }
    //     std::cout << "\n";
    // }
    // testing END

    // for (int i = 0; i < cmbr_map.size(); ++i)
    // {
    //     for (int j = 0; j < cmbr_map[i].size(); ++j)
    //     { 
    //         printVectorToFile(cmbr_map[i][j].list1, cmbr_map[i][j].list2, cmbr_map[i][j].combination, myfile);          
    //     }
    // }


    // myfile.close();

    // polygon x = getMBR(5,5);
    // polygon y = getMBR(2,2);
    // polygon cmbr;
    // cmbr = getCMBR(x, y);
    // // std::cout << boost::geometry::dsv(x) << std::endl;
    // auto it = boost::begin(boost::geometry::exterior_ring(x));
    // std::cout << boost::geometry::get<0>(*it) << std::endl;
    // std::cout << boost::geometry::get<0>(x.getpoint()) << std::endl;
    // std::cout << boost::geometry::dsv(y) << std::endl;
    // std::cout << boost::geometry::dsv(cmbr) << std::endl;


    //getting the vertices back
    // for(auto it = boost::begin(boost::geometry::exterior_ring(x)); it != boost::end(boost::geometry::exterior_ring(x)); ++it)
    // {
    //     double x = boost::geometry::get<0>(*it);
    //     double y = boost::geometry::get<1>(*it);
    //     std::cout << x << " " << y << std::endl;
    //     //use the coordinates...
    // }

    // std::cout.rdbuf(coutbuf); //reset to standard output again
    // std::cout << "File writing Completed!" << std::endl;

    return 0;
}