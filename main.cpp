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

    polygon ret;
    // threshhold constant
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
polygon** getMBRList(struct table_row *data, int li) {
    int ptr[FEATURES] = {0};
    polygon** arr = 0;

    arr = new polygon *[FEATURES];

    for (int i = 0; i < li; ++i)
    {  
        if (ptr[data[i].id-1] == 0)
         {
            arr[data[i].id-1] = new polygon [INST];
         }
         // if (data[i].id-1 == 8)
         //  {
            // std::cout << data[i].id << ": " << data[i].x << ": " <<data[i].y << ": " << ptr[data[i].id-1] <<std::endl;              
         //  } 

        arr[data[i].id-1][ptr[data[i].id-1]++] = getMBR(data[i].x, data[i].y);
    }

    return arr;
    // return my2DArray;
}

// returns CMBR for a given two MBRs
polygon getCMBR(polygon a, polygon b) {
    std::deque<polygon> output;
    polygon ret;
    boost::geometry::intersection(a, b, output);
    ret = output.front();

    return ret;
}

polygon** getCMBRList(polygon **mbrs, int a, int b) 
{
    // 2D array to hold all CMBRS. Should represents levels
    polygon** arr = 0;

    int level = 0, insid = 0, len = 10000;
    int ss = 100000;  
    arr = new polygon *[12]; 
    arr[level] = new polygon [ss]; 

    //nested loop to check all the combinations of instances and then check CMBRs
    for (int i = 0; i < len; ++i)
    { 
        for (int j = 0; j < len; ++j)
        {
            // if (mbrs[a][i] != 0 && mbrs[b][j] != 0)
            // {
                arr[level][insid++] = getCMBR(mbrs[a][i], mbrs[b][j]);   
                std::cout << i << ": " << j << std::endl;             
            // }   
            if (ss <= insid)
            {
                break;
            }       
        }
        if (ss <= insid)
        {
            break;
        }
    }

    //return 2D array with CMBRs
    return arr;
}

int main()
{
    // read data into a 2D array
    struct table_row *dat;
    dat = createArray("Seattle2012.csv");
    // std::cout << dat[0].x << std::endl;

    // calculate MBR for all the datapoints
    polygon**  mbr_array = getMBRList(dat, ROWS);
    // getMBRList(dat, ROWS);
    std::cout <<  boost::geometry::wkt(mbr_array[41][0]) << std::endl;


    // calculate CMBR recursively for each of the features
    polygon x = getMBR(5,5);
    polygon y = getMBR(2,2);
    polygon cmbr;
    cmbr = getCMBR(x, y); 
    polygon**  cmbr_array = getCMBRList(mbr_array, 0, 7);


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