#include <iostream>
#include <cmath>
#include <vector>

struct Point {
  double x;
  double y;
  Point(double x,double y) : x(x), y(y) {}
};

double dist(const Point& a, const Point& b) {
  return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));
}

int main() {

  std::vector<Point> points;

  double x,y;
  std::cin >> x >> y;
  int n=0;

  while(std::cin.good()) {

    points.push_back(Point(x,y));
    n++;
    std::cin >> x >> y;
  }

  std::cout << n << std::endl;

  for(int i=0;i<n;i++) {
    for(int j=0;j<n;j++) {
      std::cout << dist(points[i],points[j]) << " ";
    }
    std::cout << std::endl;
  }
}
