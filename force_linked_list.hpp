#pragma once

#include <cassert>
#include <vector>
#include <numeric>
#include "simd_util.hpp"

template <typename Vec>
struct ForceLinkedList {
  int32_t particle_number_ = -1;
  int32_t cell_numb_[3], all_cell_ = -1;
  Vec cell_leng_, inv_cell_leng_;
  double cutoff_length_  = 0.0, cutoff_length2_ = 0.0;
  double dt_ = 0.0;

  int32_t *cell_id_of_ptcl_ = nullptr, *neigh_cell_id_ = nullptr;
  int32_t *number_in_cell_ = nullptr;
  int32_t *cell_pointer_ = nullptr, *cell_pointer_buf_ = nullptr;
  int32_t *ptcl_id_in_cell_ = nullptr;
  int32_t *next_dst_ = nullptr;

  std::vector<int32_t>* ptcl_id_of_neigh_cell_ = nullptr;

  enum : int32_t {
    MAX_PARTNERS = 100,
    SORT_FREQ = 50,
    NUM_NEIGH_CELL = 13,
    NUM_PTCL_IN_NEIGH_CELL = 500,
  };

  int32_t GenHash(const int32_t* idx) const {
    const auto ret = idx[0] + (idx[1] + idx[2] * cell_numb_[1]) * cell_numb_[0];
#ifdef DEBUG
    assert(ret >= 0);
    assert(ret < all_cell_);
#endif
    return ret;
  }

  int32_t GenHash(const Vec& q) const {
    int32_t idx[] = {
      static_cast<int32_t>(q.x * inv_cell_leng_.x),
      static_cast<int32_t>(q.y * inv_cell_leng_.y),
      static_cast<int32_t>(q.z * inv_cell_leng_.z)
    };
    ApplyPBC(idx);
    return GenHash(idx);
  }

  void ApplyPBC(int32_t* idx) const {
    for (int i = 0; i < 3; i++) {
      if (idx[i] < 0) idx[i] += cell_numb_[i];
      if (idx[i] >= cell_numb_[i]) idx[i] -= cell_numb_[i];
    }
  }

  void Allocate() {
    cell_id_of_ptcl_       = new int32_t [particle_number_];
    neigh_cell_id_         = new int32_t [NUM_NEIGH_CELL * all_cell_];
    number_in_cell_        = new int32_t [all_cell_];
    cell_pointer_          = new int32_t [all_cell_ + 1];
    cell_pointer_buf_      = new int32_t [all_cell_ + 1];
    ptcl_id_in_cell_       = new int32_t [particle_number_];
    next_dst_              = new int32_t [particle_number_];
    ptcl_id_of_neigh_cell_ = new std::vector<int32_t> [all_cell_];
    for (int32_t i = 0; i < all_cell_; i++) {
      ptcl_id_of_neigh_cell_[i].resize(NUM_PTCL_IN_NEIGH_CELL);
    }
  }

  void Deallocate() {
    delete [] cell_id_of_ptcl_;
    delete [] neigh_cell_id_;
    delete [] number_in_cell_;
    delete [] cell_pointer_;
    delete [] cell_pointer_buf_;
    delete [] ptcl_id_in_cell_;
    delete [] next_dst_;
    delete [] ptcl_id_of_neigh_cell_;
  }

  void MakeNeighCellId() {
    int32_t icell_id = 0;
    for (int32_t iz = 0; iz < cell_numb_[2]; iz++)
      for (int32_t iy = 0; iy < cell_numb_[1]; iy++)
        for (int32_t ix = 0; ix < cell_numb_[0]; ix++) {
          int32_t jcell_id = 0;
          for (int32_t jz = -1; jz < 2; jz++)
            for (int32_t jy = -1; jy < 2; jy++)
              for (int32_t jx = -1; jx < 2; jx++) {
                int32_t idx[] = { ix + jx, iy + jy, iz + jz };
                ApplyPBC(idx);
                neigh_cell_id_[NUM_NEIGH_CELL * icell_id + jcell_id] = GenHash(idx);
                jcell_id++;
                if (jcell_id == NUM_NEIGH_CELL) goto OUT;
              }
        OUT:
          icell_id++;
        }
#ifdef DEBUG
    assert(icell_id == all_cell_);
    for (int i = 0; i < all_cell_ * NUM_NEIGH_CELL; i++) {
      assert(neigh_cell_id_[i] >= 0);
      assert(neigh_cell_id_[i] < all_cell_);
    }
#endif
  }

  void MakeCellidOfPtcl(const Vec* q) {
    std::fill(number_in_cell_,
              number_in_cell_ + all_cell_,
              0);
    for (int32_t i = 0; i < particle_number_; i++) {
      const auto hash = GenHash(q[i]);
      number_in_cell_[hash]++;
      cell_id_of_ptcl_[i] = hash;
    }
  }

  void MakeNextDest() {
    cell_pointer_[0] = cell_pointer_buf_[0] = 0;
    for (int32_t i = 0; i < all_cell_; i++) {
      const auto g_ptr = cell_pointer_[i] + number_in_cell_[i];
      cell_pointer_[i + 1] = g_ptr;
      cell_pointer_buf_[i + 1] = g_ptr;
    }

    for (int32_t i = 0; i < particle_number_; i++) {
      const auto hash = cell_id_of_ptcl_[i];
      const auto dst = cell_pointer_buf_[hash];
      next_dst_[i] = dst;
      ptcl_id_in_cell_[dst] = i;
      cell_pointer_buf_[hash]++;
    }

#ifdef DEBUG
    assert(cell_pointer_[all_cell_] == particle_number_);
#endif
  }

  void CheckSorted(const Vec* q) const {
    for (int32_t cell = 0; cell < all_cell_; cell++) {
      const auto beg = cell_pointer_[cell    ];
      const auto end = cell_pointer_[cell + 1];
      for (int32_t i = beg; i < end; i++) {
        const auto hash = GenHash(q[ptcl_id_in_cell_[i]]);
        if (hash != cell) {
          std::cerr << "particle data is not correctly sorted.\n";
          std::exit(1);
        }
      }
    }
  }

  void MakeNeighCellPtclId() {
    for (int32_t icell = 0; icell < all_cell_; icell++) {
      ptcl_id_of_neigh_cell_[icell].clear();
      const auto icell_beg = cell_pointer_[icell    ];
      const auto icell_end = cell_pointer_[icell + 1];
      ptcl_id_of_neigh_cell_[icell].insert(ptcl_id_of_neigh_cell_[icell].end(),
                                           &ptcl_id_in_cell_[icell_beg],
                                           &ptcl_id_in_cell_[icell_end]);
      for (int32_t k = 0; k < NUM_NEIGH_CELL; k++) {
        const auto jcell = neigh_cell_id_[NUM_NEIGH_CELL * icell + k];
        const auto jcell_beg = cell_pointer_[jcell    ];
        const auto jcell_end = cell_pointer_[jcell + 1];
        ptcl_id_of_neigh_cell_[icell].insert(ptcl_id_of_neigh_cell_[icell].end(),
                                             &ptcl_id_in_cell_[jcell_beg],
                                             &ptcl_id_in_cell_[jcell_end]);
      }
    }
  }

  void CalculateForceBruteForce(const Vec* q,
                                Vec* p) {
    for (int i = 0; i < particle_number_ - 1; i++) {
      for (int j = i + 1; j < particle_number_; j++) {
        const auto dx = q[j].x - q[i].x;
        const auto dy = q[j].y - q[i].y;
        const auto dz = q[j].z - q[i].z;
        const auto r2 = dx * dx + dy * dy + dz * dz;
        if (r2 > cutoff_length2_) continue;
        const auto r6 = r2 * r2 * r2;
        const auto df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt_;
        p[i].x += df * dx;
        p[i].y += df * dy;
        p[i].z += df * dz;
        p[j].x -= df * dx;
        p[j].y -= df * dy;
        p[j].z -= df * dz;
      }
    }
  }

  void CalculateForceNaive(const Vec* q,
                           Vec* p) {
    for (int32_t icell = 0; icell < all_cell_; icell++) {
      const auto icell_beg  = cell_pointer_[icell    ];
      const auto icell_size = number_in_cell_[icell];
      const int32_t* pid_of_neigh_cell_loc = &ptcl_id_of_neigh_cell_[icell][0];
      const int32_t num_of_neigh_cell      = ptcl_id_of_neigh_cell_[icell].size();
      for (int32_t l = 0; l < icell_size; l++) {
        const auto i  = ptcl_id_in_cell_[icell_beg + l];
        const auto qi = q[i];
        double pfx = 0.0, pfy = 0.0, pfz = 0.0;
        for (int32_t k = l + 1; k < num_of_neigh_cell; k++) {
          const auto j = pid_of_neigh_cell_loc[k];
          const auto dx = q[j].x - qi.x;
          const auto dy = q[j].y - qi.y;
          const auto dz = q[j].z - qi.z;
          const auto r2 = (dx * dx + dy * dy + dz * dz);
          if (r2 > cutoff_length2_) continue;
          const auto r6 = r2 * r2 * r2;
          const auto df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt_;
          pfx += df * dx;
          pfy += df * dy;
          pfz += df * dz;
          p[j].x -= df * dx;
          p[j].y -= df * dy;
          p[j].z -= df * dz;
        }
        p[i].x += pfx;
        p[i].y += pfy;
        p[i].z += pfz;
      }
    }
  }

  void CalculateForce1x4(const Vec* q,
                         Vec* p) {
    const auto vzero = _mm256_setzero_pd();
    const auto vcl2  = _mm256_set1_pd(cutoff_length2_);
    const auto vc24  = _mm256_set1_pd(24.0 * dt_);
    const auto vc48  = _mm256_set1_pd(48.0 * dt_);

    for (int32_t icell = 0; icell < all_cell_; icell++) {
      const auto icell_beg  = cell_pointer_[icell];
      const auto icell_size = number_in_cell_[icell];
      const int32_t* pid_of_neigh_cell_loc = &ptcl_id_of_neigh_cell_[icell][0];
      const int32_t num_of_neigh_cell      = ptcl_id_of_neigh_cell_[icell].size();
      for (int32_t l = 0; l < icell_size; l++) {
        const auto i  = ptcl_id_in_cell_[icell_beg + l];
        const auto vqi = _mm256_loadu_pd(&q[i].x);

        auto vpf = _mm256_setzero_pd();
        const auto num_loop = num_of_neigh_cell - (l + 1);
        for (int32_t k = 0; k < (num_loop / 4) * 4; k += 4) {
          const auto ja = pid_of_neigh_cell_loc[k + l + 1];
          const auto jb = pid_of_neigh_cell_loc[k + l + 2];
          const auto jc = pid_of_neigh_cell_loc[k + l + 3];
          const auto jd = pid_of_neigh_cell_loc[k + l + 4];

          auto vqja = _mm256_loadu_pd(&q[ja].x);
          auto vqjb = _mm256_loadu_pd(&q[jb].x);
          auto vqjc = _mm256_loadu_pd(&q[jc].x);
          auto vqjd = _mm256_loadu_pd(&q[jd].x);

          auto vdq_a = _mm256_sub_pd(vqja, vqi);
          auto vdq_b = _mm256_sub_pd(vqjb, vqi);
          auto vdq_c = _mm256_sub_pd(vqjc, vqi);
          auto vdq_d = _mm256_sub_pd(vqjd, vqi);

          __m256d dvx, dvy, dvz;
          transpose_4x4(vdq_a, vdq_b, vdq_c, vdq_d,
                        dvx, dvy, dvz);

          auto vr2 = _mm256_fmadd_pd(dvz, dvz,
                                     _mm256_fmadd_pd(dvy, dvy,
                                                     _mm256_mul_pd(dvx, dvx)));
          auto mask = _mm256_cmp_pd(vr2, vcl2, _CMP_LE_OS);
          const int32_t hash = _mm256_movemask_pd(mask);
          if (hash == 0) continue;

          auto vr6 = _mm256_mul_pd(_mm256_mul_pd(vr2, vr2), vr2);
          auto vdf = _mm256_div_pd(_mm256_fmsub_pd(vc24, vr6, vc48),
                                   _mm256_mul_pd(_mm256_mul_pd(vr6, vr6), vr2));
          vdf = _mm256_blendv_pd(vzero, vdf, mask);

          auto vdf_a = _mm256_permute4x64_pd(vdf, 0x00);
          auto vdf_b = _mm256_permute4x64_pd(vdf, 0x55);
          auto vdf_c = _mm256_permute4x64_pd(vdf, 0xaa);
          auto vdf_d = _mm256_permute4x64_pd(vdf, 0xff);

          vpf = _mm256_fmadd_pd(vdf_a, vdq_a, vpf);
          vpf = _mm256_fmadd_pd(vdf_b, vdq_b, vpf);
          vpf = _mm256_fmadd_pd(vdf_c, vdq_c, vpf);
          vpf = _mm256_fmadd_pd(vdf_d, vdq_d, vpf);

          auto vpja = _mm256_loadu_pd(&p[ja].x);
          vpja = _mm256_fnmadd_pd(vdf_a, vdq_a, vpja);
          _mm256_storeu_pd(&p[ja].x, vpja);

          auto vpjb = _mm256_loadu_pd(&p[jb].x);
          vpjb = _mm256_fnmadd_pd(vdf_b, vdq_b, vpjb);
          _mm256_storeu_pd(&p[jb].x, vpjb);

          auto vpjc = _mm256_loadu_pd(&p[jc].x);
          vpjc = _mm256_fnmadd_pd(vdf_c, vdq_c, vpjc);
          _mm256_storeu_pd(&p[jc].x, vpjc);

          auto vpjd = _mm256_loadu_pd(&p[jd].x);
          vpjd = _mm256_fnmadd_pd(vdf_d, vdq_d, vpjd);
          _mm256_storeu_pd(&p[jd].x, vpjd);
        }
        auto vpi = _mm256_loadu_pd(&p[i].x);
        vpi = _mm256_add_pd(vpi, vpf);
        _mm256_storeu_pd(&p[i].x, vpi);

        const auto qi = q[i];
        double pfx = 0.0, pfy = 0.0, pfz = 0.0;
        for (int32_t k = (num_loop / 4) * 4; k < num_loop; k++) {
          const auto j = pid_of_neigh_cell_loc[k + l + 1];
          const auto dx = q[j].x - qi.x;
          const auto dy = q[j].y - qi.y;
          const auto dz = q[j].z - qi.z;
          const auto r2 = (dx * dx + dy * dy + dz * dz);
          if (r2 > cutoff_length2_) continue;
          const auto r6 = r2 * r2 * r2;
          const auto df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt_;
          pfx += df * dx;
          pfy += df * dy;
          pfz += df * dz;
          p[j].x -= df * dx;
          p[j].y -= df * dy;
          p[j].z -= df * dz;
        }
        p[i].x += pfx;
        p[i].y += pfy;
        p[i].z += pfz;
      }
    }
  }

  void CalculateForcePair(const Vec* q,
                          Vec* p,
                          const int i,
                          const int j) {
    const auto dx = q[j].x - q[i].x;
    const auto dy = q[j].y - q[i].y;
    const auto dz = q[j].z - q[i].z;
    const auto r2 = (dx * dx + dy * dy + dz * dz);
    if (r2 > cutoff_length2_) return;
    const auto r6 = r2 * r2 * r2;
    const auto df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt_;
    p[i].x += df * dx;
    p[i].y += df * dy;
    p[i].z += df * dz;
    p[j].x -= df * dx;
    p[j].y -= df * dy;
    p[j].z -= df * dz;
  }

  void CalculateForce4x1(const Vec* q,
                         Vec* p) {
    const auto vzero = _mm256_setzero_pd();
    const auto vcl2  = _mm256_set1_pd(cutoff_length2_);
    const auto vc24  = _mm256_set1_pd(24.0 * dt_);
    const auto vc48  = _mm256_set1_pd(48.0 * dt_);

    for (int32_t icell = 0; icell < all_cell_; icell++) {
      const auto icell_beg  = cell_pointer_[icell];
      const auto icell_size = number_in_cell_[icell];
      const int32_t* pid_of_neigh_cell_loc = &ptcl_id_of_neigh_cell_[icell][0];
      const int32_t num_of_neigh_cell      = ptcl_id_of_neigh_cell_[icell].size();
      for (int32_t l = 0; l < (icell_size / 4) * 4; l += 4) {
        const auto ia = ptcl_id_in_cell_[icell_beg + l    ];
        const auto ib = ptcl_id_in_cell_[icell_beg + l + 1];
        const auto ic = ptcl_id_in_cell_[icell_beg + l + 2];
        const auto id = ptcl_id_in_cell_[icell_beg + l + 3];

        auto vqia = _mm256_loadu_pd(&q[ia].x);
        auto vqib = _mm256_loadu_pd(&q[ib].x);
        auto vqic = _mm256_loadu_pd(&q[ic].x);
        auto vqid = _mm256_loadu_pd(&q[id].x);

        auto vpfa = _mm256_setzero_pd();
        auto vpfb = _mm256_setzero_pd();
        auto vpfc = _mm256_setzero_pd();
        auto vpfd = _mm256_setzero_pd();
        for (int32_t k = l + 4; k < num_of_neigh_cell; k++) {
          const auto j = pid_of_neigh_cell_loc[k];

          auto vqj = _mm256_loadu_pd(&q[j].x);

          auto vdq_a = _mm256_sub_pd(vqj, vqia);
          auto vdq_b = _mm256_sub_pd(vqj, vqib);
          auto vdq_c = _mm256_sub_pd(vqj, vqic);
          auto vdq_d = _mm256_sub_pd(vqj, vqid);

          __m256d dvx, dvy, dvz;
          transpose_4x4(vdq_a, vdq_b, vdq_c, vdq_d,
                        dvx, dvy, dvz);

          auto vr2 = _mm256_fmadd_pd(dvz, dvz,
                                     _mm256_fmadd_pd(dvy, dvy,
                                                     _mm256_mul_pd(dvx, dvx)));
          auto mask = _mm256_cmp_pd(vr2, vcl2, _CMP_LE_OS);
          const int32_t hash = _mm256_movemask_pd(mask);
          if (hash == 0) continue;

          auto vr6 = _mm256_mul_pd(_mm256_mul_pd(vr2, vr2), vr2);
          auto vdf = _mm256_div_pd(_mm256_fmsub_pd(vc24, vr6, vc48),
                                   _mm256_mul_pd(_mm256_mul_pd(vr6, vr6), vr2));
          vdf = _mm256_blendv_pd(vzero, vdf, mask);

          auto vdf_a = _mm256_permute4x64_pd(vdf, 0x00);
          auto vdf_b = _mm256_permute4x64_pd(vdf, 0x55);
          auto vdf_c = _mm256_permute4x64_pd(vdf, 0xaa);
          auto vdf_d = _mm256_permute4x64_pd(vdf, 0xff);

          vpfa = _mm256_fmadd_pd(vdf_a, vdq_a, vpfa);
          vpfb = _mm256_fmadd_pd(vdf_b, vdq_b, vpfb);
          vpfc = _mm256_fmadd_pd(vdf_c, vdq_c, vpfc);
          vpfd = _mm256_fmadd_pd(vdf_d, vdq_d, vpfd);

          auto vpj = _mm256_loadu_pd(&p[j].x);
          vpj = _mm256_fnmadd_pd(vdf_a, vdq_a, vpj);
          vpj = _mm256_fnmadd_pd(vdf_b, vdq_b, vpj);
          vpj = _mm256_fnmadd_pd(vdf_c, vdq_c, vpj);
          vpj = _mm256_fnmadd_pd(vdf_d, vdq_d, vpj);
          _mm256_storeu_pd(&p[j].x, vpj);
        }
        auto vpia = _mm256_loadu_pd(&p[ia].x);
        auto vpib = _mm256_loadu_pd(&p[ib].x);
        auto vpic = _mm256_loadu_pd(&p[ic].x);
        auto vpid = _mm256_loadu_pd(&p[id].x);

        vpia = _mm256_add_pd(vpia, vpfa);
        vpib = _mm256_add_pd(vpib, vpfb);
        vpic = _mm256_add_pd(vpic, vpfc);
        vpid = _mm256_add_pd(vpid, vpfd);

        _mm256_storeu_pd(&p[ia].x, vpia);
        _mm256_storeu_pd(&p[ib].x, vpib);
        _mm256_storeu_pd(&p[ic].x, vpic);
        _mm256_storeu_pd(&p[id].x, vpid);

        CalculateForcePair(q, p, ia, ib);
        CalculateForcePair(q, p, ia, ic);
        CalculateForcePair(q, p, ia, id);
        CalculateForcePair(q, p, ib, ic);
        CalculateForcePair(q, p, ib, id);
        CalculateForcePair(q, p, ic, id);
      }
      for (int32_t l = (icell_size / 4) * 4; l < icell_size; l++) {
        const auto i  = ptcl_id_in_cell_[icell_beg + l];
        const auto qi = q[i];
        double pfx = 0.0, pfy = 0.0, pfz = 0.0;
        for (int32_t k = l + 1; k < num_of_neigh_cell; k++) {
          const auto j = pid_of_neigh_cell_loc[k];
          const auto dx = q[j].x - qi.x;
          const auto dy = q[j].y - qi.y;
          const auto dz = q[j].z - qi.z;
          const auto r2 = (dx * dx + dy * dy + dz * dz);
          if (r2 > cutoff_length2_) continue;
          const auto r6 = r2 * r2 * r2;
          const auto df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2)) * dt_;
          pfx += df * dx;
          pfy += df * dy;
          pfz += df * dz;
          p[j].x -= df * dx;
          p[j].y -= df * dy;
          p[j].z -= df * dz;
        }
        p[i].x += pfx;
        p[i].y += pfy;
        p[i].z += pfz;
      }
    }
  }

  ForceLinkedList(const double cutoff_length,
                  const double Lx,
                  const double Ly,
                  const double Lz,
                  const double dt,
                  const int particle_number) {
    cell_numb_[0]   = static_cast<int32_t>(Lx / cutoff_length);
    cell_numb_[1]   = static_cast<int32_t>(Ly / cutoff_length);
    cell_numb_[2]   = static_cast<int32_t>(Lz / cutoff_length);
    all_cell_       = cell_numb_[0] * cell_numb_[1] * cell_numb_[2];

    cell_leng_.x    = Lx / cell_numb_[0];
    cell_leng_.y    = Ly / cell_numb_[1];
    cell_leng_.z    = Lz / cell_numb_[2];

    cutoff_length_  = cutoff_length;
    cutoff_length2_ = cutoff_length * cutoff_length;
    dt_             = dt;

    particle_number_ = particle_number;
  }
  ~ForceLinkedList() {
    Deallocate();
  }

  // disable copy
  const ForceLinkedList<Vec>& operator = (const ForceLinkedList<Vec>& obj) = delete;
  ForceLinkedList<Vec>(const ForceLinkedList<Vec>& obj) = delete;

  // disable move
  ForceLinkedList<Vec>& operator = (ForceLinkedList<Vec>&& obj) = delete;
  ForceLinkedList<Vec>(ForceLinkedList<Vec>&& obj) = delete;

  void Initialize() {
    inv_cell_leng_.x = 1.0 / cell_leng_.x;
    inv_cell_leng_.y = 1.0 / cell_leng_.y;
    inv_cell_leng_.z = 1.0 / cell_leng_.z;

    Allocate();
    MakeNeighCellId();
  }

  void RegisterCellIdx(Vec* q) {
    MakeCellidOfPtcl(q);
    MakeNextDest();
#ifdef DEBUG
    CheckSorted(q);
#endif
    MakeNeighCellPtclId();
  }

  void CalculateForce(const Vec* q,
                      Vec* p) {
#ifdef SIMD1x4
    CalculateForce1x4(q, p);
#elif SIMD4x1
    CalculateForce4x1(q, p);
#elif BRUTE_FORCE
    CalculateForceBruteForce(q, p);
#else
    CalculateForceNaive(q, p);
#endif
  }
};
