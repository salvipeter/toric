#include "toric.hh"

#include <algorithm>
#include <fstream>
#include <cmath>
#include <numeric>

using namespace Geometry;

static size_t binomial(size_t n, size_t k) {
  if (k > n)
    return 0;
  size_t result = 1;
  for (size_t d = 1; d <= k; ++d, --n)
    result = result * n / d;
  return result;
}

// 1 for interior control points
// C(dk, i) for the i-th control point on the k-th side (of degree dk, i.e., dk+1 CP)
double Toric::coefficient(const Index &v) const {
  double result = 1;
  for (size_t i = 0; i < n; ++i)
    if (L(i, v) == 0) {
      size_t im = (i + n - 1) % n;
      auto p = Point2D(v[0], v[1]);
      auto p0 = Point2D(domain[im][0], domain[im][1]);
      auto p1 = Point2D(domain[i][0], domain[i][1]);
      size_t d = std::round((p1 - p0).norm() / Vector2D(lines[i][0], lines[i][1]).norm());
      size_t k = std::round((p - p0).norm() / Vector2D(lines[i][0], lines[i][1]).norm());
      result = binomial(d, k);
      break;
    }
  return result;
}

Point3D Toric::eval(const Point2D &u) const {
  auto blend = [&](const Index &v) {
    double result = coefficient(v);
    for (size_t i = 0; i < n; ++i) {
      double lu = L(i, u);
      size_t lv = std::round(L(i, v));
      result *= std::pow(lu, lv);
    }
    return result;
  };

  Point3D result(0, 0, 0);
  double bsum = 0;
  for (const auto &[v, p] : cnet) {
    auto b = blend(v) * weights.at(v);
    result += p * b;
    bsum += b;
  }
  return result / bsum;
}

// - side i between vertices i-1 and i
// - assumes that the domain points are given in CCW order
void Toric::updateLines() {
  lines.clear();
  for (size_t i = 0; i < n; ++i) {
    size_t im = (i + n - 1) % n;
    int u = domain[im][1] - domain[i][1];
    int v = domain[i][0] - domain[im][0];
    int w = std::gcd(u, v);
    u /= w; v /= w;
    lines.emplace_back(u, v, -(domain[i][0] * u + domain[i][1] * v));
  }
}


// TRP file loader

Toric Toric::load(std::string filename) {
  Toric result;
  std::ifstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  size_t n_cpts;
  f >> result.n;
  for (size_t i = 0; i < result.n; ++i) {
    Index v;
    f >> v[0] >> v[1];
    result.domain.push_back(v);
  }
  f >> n_cpts;
  for (size_t i = 0; i < n_cpts; ++i) {
    Index v;
    Point3D p;
    double w;
    f >> v[0] >> v[1] >> p[0] >> p[1] >> p[2] >> w;
    result.cnet[v] = p;
    result.weights[v] = w;
  }
  result.updateLines();
  return result;
}


// Mesh generation

static size_t meshSize(size_t n, size_t resolution) {
  if (n == 3)
    return (resolution + 1) * (resolution + 2) / 2;
  if (n == 4)
    return (resolution + 1) * (resolution + 1);
  return 1 + n * resolution * (resolution + 1) / 2;
}

Point2DVector Toric::parameters(size_t resolution) const {
  Point2DVector parameters;
  parameters.reserve(meshSize(domain.size(), resolution));

  Point2DVector vertices;
  std::transform(domain.begin(), domain.end(), std::back_inserter(vertices),
                 [](const Index &v) { return Point2D(v[0], v[1]); });

  if (n == 3) {
    for (size_t j = 0; j <= resolution; ++j) {
      double u = (double)j / resolution;
      auto p = vertices[0] * u + vertices[2] * (1 - u);
      auto q = vertices[1] * u + vertices[2] * (1 - u);
      for (size_t k = 0; k <= j; ++k) {
        double v = j == 0 ? 1.0 : (double)k / j;
        parameters.push_back(p * (1 - v) + q * v);
      }
    }
  } else if (n == 4) {
    for (size_t j = 0; j <= resolution; ++j) {
      double u = (double)j / resolution;
      auto p = vertices[0] * (1 - u) + vertices[1] * u;
      auto q = vertices[3] * (1 - u) + vertices[2] * u;
      for (size_t k = 0; k <= resolution; ++k) {
        double v = (double)k / resolution;
        parameters.push_back(p * (1 - v) + q * v);
      }
    }
  } else { // n > 4
    Point2D center(0, 0);
    for (size_t i = 0; i < n; ++i)
      center += vertices[i];      // probably should weight by edge length?
    center /= n;

    parameters.push_back(center);
    for (size_t j = 1; j <= resolution; ++j) {
      double u = (double)j / (double)resolution;
      for (size_t k = 0; k < n; ++k)
        for (size_t i = 0; i < j; ++i) {
          double v = (double)i / (double)j;
          Point2D ep = vertices[(k+n-1)%n] * (1.0 - v) + vertices[k] * v;
          Point2D p = center * (1.0 - u) + ep * u;
          parameters.push_back(p);
        }
    }
  }
  return parameters;
}

static TriMesh meshTopology(size_t n, size_t resolution) {
  TriMesh mesh;
  mesh.resizePoints(meshSize(n, resolution));

  if (n == 3) {
    size_t prev = 0, current = 1;
    for (size_t i = 0; i < resolution; ++i) {
      for (size_t j = 0; j < i; ++j) {
        mesh.addTriangle(current + j, current + j + 1, prev + j);
        mesh.addTriangle(current + j + 1, prev + j + 1, prev + j);
      }
      mesh.addTriangle(current + i, current + i + 1, prev + i);
      prev = current;
      current += i + 2;
    }
  } else if (n == 4) {
    for (size_t i = 0; i < resolution; ++i)
      for (size_t j = 0; j < resolution; ++j) {
        size_t index = i * (resolution + 1) + j;
        mesh.addTriangle(index, index + resolution + 1, index + 1);
        mesh.addTriangle(index + 1, index + resolution + 1, index + resolution + 2);
      }
  } else { // n > 4
    size_t inner_start = 0, outer_vert = 1;
    for (size_t layer = 1; layer <= resolution; ++layer) {
      size_t inner_vert = inner_start, outer_start = outer_vert;
      for (size_t side = 0; side < n; ++side) {
        size_t vert = 0;
        while(true) {
          size_t next_vert = (side == n - 1 && vert == layer - 1) ? outer_start : (outer_vert + 1);
          mesh.addTriangle(inner_vert, outer_vert, next_vert);
          ++outer_vert;
          if (++vert == layer)
            break;
          size_t inner_next = (side == n - 1 && vert == layer - 1) ? inner_start : (inner_vert + 1);
          mesh.addTriangle(inner_vert, next_vert, inner_next);
          inner_vert = inner_next;
        }
      }
      inner_start = outer_start;
    }
  }
  return mesh;
}

TriMesh Toric::tessellate(size_t resolution) const {
  auto mesh = meshTopology(domain.size(), resolution);
  auto params = parameters(resolution);
  for (size_t i = 0; i < params.size(); ++i)
    mesh[i] = eval(params[i]);
  return mesh;
}
