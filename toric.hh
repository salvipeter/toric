#pragma once

#include <geometry.hh>
#include <map>

class Toric {
public:
  static Toric load(std::string filename);
  Geometry::Point3D eval(const Geometry::Point2D &u) const;
  Geometry::TriMesh tessellate(size_t resolution) const;

private:
  using Index = std::array<int, 2>;

  template <typename T> double L(size_t k, const T &u) const;
  double coefficient(const Index &v) const;
  void updateLines();
  Geometry::Point2DVector parameters(size_t resolution) const;

  size_t n;
  std::vector<Index> domain;
  std::map<Index, Geometry::Point3D> cnet;
  std::map<Index, double> weights;
  std::vector<Geometry::Vector3D> lines;
};

template <typename T>
double Toric::L(size_t k, const T &u) const {
  return lines[k][0] * u[0] + lines[k][1] * u[1] + lines[k][2];
}
