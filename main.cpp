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

#define INST 3000
#define FEATURES 50

typedef boost::geometry::model::d2::point_xy<double> point_xy;
typedef boost::geometry::model::polygon<point_xy > polygon;


// returns MBR for a given data point
polygon getMBR(double px, double py) {
    using namespace boost::assign;    

    polygon ret;
    // threshhold constant
    double const d = 2.0;

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
void getMBRList(struct table_row *data, int li) {
    int ptr[li+1] = {0};
    // polygon** arr;
    std::array<std::array<polygon, INST>, FEATURES> arr;

    // arr = new polygon *[INST];

    for (int i = 0; i < li; ++i)
    {  
        // if (ptr[i] == 0)
        //  {
        //     arr[i] = new polygon [li];
        //  } 
        std::cout << data[i].id <<std::endl;
        arr[data[i].id-1][ptr[data[i].id-1]++] = getMBR(data[i].x, data[i].y);   
        std::cout << boost::geometry::wkt(getMBR(data[i].x, data[i].y)) <<std::endl;
        std::cout << boost::geometry::wkt(arr[data[i].id-1][data[i].id-1]) <<std::endl;
        if (true)
        {
            std::cout << data[i].id <<std::endl;
            std::cout << i <<std::endl;

            break;
        }

    }

    // return arr;
}

// returns CMBR for a given two MBRs
polygon getCMBR(polygon a, polygon b) {
    std::deque<polygon> output;
    polygon ret;
    boost::geometry::intersection(a, b, output);
    ret = output.front();

    return ret;
}


int main()
{

    // polygon green, blue;
    // polygon* mbr_array;
    // std::array<std::array<polygon, INST>, FEATURES> mbr_array;

    // boost::geometry::read_wkt("POLYGON((0 0,0 10, 10 10, 10 0, 0 0))", green);
  

    // boost::geometry::read_wkt("POLYGON((5 0,5 5,15 5,15 0,5 0))", blue);
    
    // std::cout << boost::geometry::wkt(green) << std::endl;
    // std::cout << boost::geometry::wkt(blue) << std::endl;

  

    // std::deque<polygon> routput;
    // boost::geometry::intersection(g, b, routput);
    // std::cout << boost::geometry::wkt(routput) << std::endl;
    

    // read data into a 2D array
    struct table_row *dat;
    dat = createArray("Seattle2012.csv");
    // std::cout << dat[0].x << std::endl;

    // calculate MBR for all the datapoints
    // mbr_array = getMBRList(dat, FEATURES);
    getMBRList(dat, FEATURES);

    // calculate CMBR recursively for each of the features
    polygon x = getMBR(5,5);
    polygon y = getMBR(2,2);
    polygon cmbr;
    cmbr = getCMBR(x, y); 

    std::cout << boost::geometry::wkt(cmbr) << std::endl;

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