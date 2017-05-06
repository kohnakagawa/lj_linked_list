#include <iostream>
#include <chrono>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <random>
#include <cmath>
#include <cassert>
#include "force_linked_list.hpp"

struct double4 { double x, y, z, w; };
enum {X = 0, Y, Z, W, D};

const double density = 1.0;
const int N = 400000;
const double L[] {50.0, 50.0, 50.0};
const double dt = 0.001;
const double CUTOFF_LENGTH = 3.0;
const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;

void add_particle(double x, double y, double z,
                  double4* q,
                  int& particle_number) {
  static std::mt19937 mt(2);
  std::uniform_real_distribution<double> ud(0.0, 0.1);
  q[particle_number].x = x + ud(mt);
  q[particle_number].y = y + ud(mt);
  q[particle_number].z = z + ud(mt);
  particle_number++;
}

void init(double4* q,
          double4* p,
          int& particle_number) {
  const double s = 1.0 / std::pow(density * 0.25, 1.0 / 3.0);
  const double hs = s * 0.5;
  int sx = static_cast<int>(L[0] / s);
  int sy = static_cast<int>(L[1] / s);
  int sz = static_cast<int>(L[2] / s);
  for (int iz = 0; iz < sz; iz++) {
    for (int iy = 0; iy < sy; iy++) {
      for (int ix = 0; ix < sx; ix++) {
        double x = ix * s;
        double y = iy * s;
        double z = iz * s;
        add_particle(x       ,y      ,z,      q, particle_number);
        add_particle(x       ,y + hs ,z + hs, q, particle_number);
        add_particle(x + hs  ,y      ,z + hs, q, particle_number);
        add_particle(x + hs  ,y + hs ,z,      q, particle_number);
      }
    }
  }
  for (int i = 0; i < particle_number; i++) {
    p[i].x = 0.0;
    p[i].y = 0.0;
    p[i].z = 0.0;
  }
}

void show_momentum(double4* p,
                   const int pn) {
#if 1
  for (int i = 0; i < 5; i++) {
    std::cerr << p[i].x << " " << p[i].y << " " << p[i].z << std::endl;
  }
  for (int i = pn - 5; i < pn; i++) {
    std::cerr << p[i].x << " " << p[i].y << " " << p[i].z << std::endl;
  }
#else
  for (int i = 0; i < pn; i++) {
    std::cerr << p[i].x << " " << p[i].y << " " << p[i].z << std::endl;
  }
#endif
}

#define MEASURE(func)                                                   \
  do {                                                                  \
    using namespace std::chrono;                                        \
    func;                                                               \
    auto beg = system_clock::now();                                     \
    for (int i = 0; i < 100; i++) func;                                 \
    auto end = system_clock::now();                                     \
    std::cout << duration_cast<milliseconds>(end - beg).count() << " [ms]\n"; \
  } while (0)


int main(int argc, char* argv[]) {
  double4* q = new double4[N];
  double4* p = new double4[N];

  // create particles
  int particle_number = 0;
  init(q, p, particle_number);

  std::cout << "# N = " << particle_number << " ";

  ForceLinkedList<double4> f_linked_list(CUTOFF_LENGTH, L[X], L[Y], L[Z], dt, particle_number);
  f_linked_list.Initialize();
  f_linked_list.RegisterCellIdx(q);

  MEASURE(f_linked_list.CalculateForce(q, p));
  show_momentum(p, particle_number);

  delete [] q;
  delete [] p;
}
