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
         // if (data[i].id-1 == 8)
         //  {
         //    std::cout << data[i].id << ": " << data[i].x << ": " <<data[i].y << ": " << ptr[data[i].id-1] <<std::endl;              
         //  } 

        // calculate MBR using the getMBR() and assign it to the relavant feature instance
        arr[data[i].id-1][ptr[data[i].id-1]++] = getMBR(data[i].x, data[i].y);
    }

    return arr;
    // return my2DArray;
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

polygon** getCMBRList(polygon **mbrs, int a, int b, int *ptr) 
{
    // 2D array to hold all CMBRS. Should represents levels
    polygon** arr = 0;

    // level controls the layer
    int level = 0;
    // controls next available sapce in the layer array
    int insid = 0;
    // memory safe variable
    // controls max number of CMBRs for each layer
    int ss = 100;  

    // Create 1st-D for the output. Assumed we have only 12 features
    arr = new polygon *[12]; 

    // allocate 2nd-D for the layer 1
    arr[level] = new polygon [ss]; 

    // temporary varibale to track of calculated CMBR
    polygon cmbr;

    //nested loop to check all the CMBRs for all combinations of instances 
    for (int i = 0; i < ptr[a]; ++i)
    { 
        for (int j = 0; j < ptr[b]; ++j)
        {
            cmbr= getCMBR(mbrs[a][i], mbrs[b][j]);
            // std::cout << !boost::geometry::is_empty(cmbr) <<std::endl;              
            // checking if the CMBR exists
            if (!boost::geometry::is_empty(cmbr))    
            {
                arr[level][insid++] = getCMBR(mbrs[a][i], mbrs[b][j]);   
            } 
            // memory safe condition
            // if we reach max memory for the layer combination allocation, 
            // it may stop assigning further CMBRS
            if (ss <= insid)
            {           
                break;
            }       
        }
        // prints the outer feature instance id-1
        std::cout << i  << std::endl;

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

int main()
{
    //array to hold the number of instances for each feature
    static int feature_sizes[FEATURES] = {0};

    // read data into a table_row structure type 1D array
    struct table_row *dat;
    dat = createArray("Seattle2012.csv");

    // calculate MBR for all the datapoints. 
    // returns a 2D array. 1st-D : Features, 2nd-D: instances per each feature 
    polygon**  mbr_array = getMBRList(dat, ROWS, feature_sizes);

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


    // calculate CMBR for 2 selected features. Works as the layer 1 output
    polygon**  cmbr_array = getCMBRList(mbr_array, 0, 7, feature_sizes);

     

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