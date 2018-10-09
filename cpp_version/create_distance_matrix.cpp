#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cassert>

typedef std::vector<double> Point;

double dist(const Point& a, const Point& b)
{

    double result = 0.0;

    assert(a.size() == b.size());
    for (size_t i = 0; i < a.size(); i++) {
        result += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return sqrt(result);
}

Point get_point(std::string s)
{
    std::stringstream sstr(s);

    Point p;
    double next;
    sstr >> next;
    while (sstr.good()) {
        p.push_back(next);
        sstr >> next;
    }
    return p;
}

int main(int argc, char** argv)
{

    if (argc < 2) {
        std::cerr << "Arguments" << std::endl;
        return 0;
    }

    std::vector<Point> points;

    std::ifstream ifstr;

    ifstr.open(argv[1]);

    if (not iftsr.good()) {
        std::cerr << "Cannot open file " << argv[1] << std::endl;
        return 1;
    }

    std::string next_line;

    std::getline(ifstr, next_line);

    Point p = get_point(next_line);

    int n = 0;

    while (ifstr.good()) {
        points.push_back(p);
        n++;
        std::getline(ifstr, next_line);
        p = get_point(next_line);
    }

    std::cout << n << std::endl;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << dist(points[i], points[j]) << " ";
        }
        std::cout << std::endl;
    }

    ifstr.close();

    int no_queries = 0;
    std::vector<Point> queries;

    if (argc >= 3) {

        ifstr.open(argv[2]);
        std::getline(ifstr, next_line);
        p = get_point(next_line);

        while (ifstr.good()) {
            queries.push_back(p);
            no_queries++;
            std::getline(ifstr, next_line);
            p = get_point(next_line);
        }
        ifstr.close();
    }

    std::cout << no_queries << std::endl;
    for (int i = 0; i < no_queries; i++) {

        Point& query = queries[i];

        for (int j = 0; j < n; j++) {
            std::cout << dist(query, points[j]) << " ";
        }
        std::cout << std::endl;
    }

}
