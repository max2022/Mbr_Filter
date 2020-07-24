#include <iostream>
#include <deque>
#include <string>
#include <vector>
#include <array>

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

typedef boost::geometry::model::d2::point_xy<float> point_xy;
typedef boost::geometry::model::polygon<point_xy > polygon;

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

    // calculate CMBR for 2 selected features. Works as the layer 1 output
    std::vector< std::vector<polygon> > cmbr_array = getCMBRList(mbr_array, feature_sizes, feature_ids);

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