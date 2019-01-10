#include <cmath>
#include <ctime>
#include <vector>

class Wolff2D {
  // Number of sites
  int N_, L_;
  // Coefficient
  float K_;
  // Current configuration
  std::vector<int> conf_, confbin_;
  // Current cluster (indices)
  std::vector<int> cluster_;
  // Remember if the spin is in cluster
  std::vector<bool> not_in_cluster_;
  // Nearest neighbor map (right, left, down, up)
  std::vector<std::vector<int>> nnmap_;
  // Addition probability
  float P_add_;

  std::vector<int> FindNeighbors(const int site) {
    std::vector<int> nn;
    int x = site % L_, y = site / L_;
    // Right
    if (x == L_ - 1) {
      nn.push_back(site + 1 - L_);
    } else {
      nn.push_back(site + 1);
    }

    // Left
    if (x == 0) {
      nn.push_back(site - 1 + L_);
    } else {
      nn.push_back(site - 1);
    }

    // Down
    if (y == L_ - 1) {
      nn.push_back(x);
    } else {
      nn.push_back(site + L_);
    }

    // Up
    if (y == 0) {
      nn.push_back(L_ * (L_ - 1) + x);
    } else {
      nn.push_back(site - L_);
    }

    return nn;
  }

  int ToBinary(const int x) {
    if (x == 1) {
      return 1;
    } else {
      return 0;
    }
  }

  void CreateCluster() {
    // Start by adding random site
    int site = rand() % N_, checked_ind = 0;
    int cluster_direction = conf_[site];
    cluster_.push_back(site);
    not_in_cluster_[site] = false;

    while (checked_ind < cluster_.size()) {
      for (std::size_t i = 0; i < nnmap_[cluster_[checked_ind]].size(); i++) {
        site = nnmap_[cluster_[checked_ind]][i];
        if (not_in_cluster_[site] && (cluster_direction * conf_[site] > 0)) {
          // Random float in [0, 1]
          float r = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
          if (P_add_ > r) {
            // Add site to cluster
            not_in_cluster_[site] = false;
            cluster_.push_back(site);
          }
        }
      }
      checked_ind++;
    }
  }

  void FlipCluster() {
    for (std::size_t i = 0; i < cluster_.size(); i++) {
      conf_[cluster_[i]] *= -1;
      confbin_[cluster_[i]] = ToBinary(conf_[cluster_[i]]);
      // Reset cluster
      not_in_cluster_[cluster_[i]] = true;
    }
    cluster_.clear();
  }

 public:
  int Nsites() { return N_; }

  void Init(const int L, const float K) {
    // Set random seed
    srand(time(NULL));

    // Define parameters
    L_ = L;
    N_ = L * L;
    K_ = K;
    P_add_ = 1 - exp(-2.0 * K_);

    // Initialize configuration and nn map
    for (int i = 0; i < N_; i++) {
      conf_.push_back(2 * (rand() % 2) - 1);
      confbin_.push_back(ToBinary(conf_[i]));
      nnmap_.push_back(FindNeighbors(i));
      not_in_cluster_.push_back(true);
    }
  }

  std::vector<std::vector<int>> Run(const int Nsamples, const int Ncorr,
                                    const int Nburn) {
    std::vector<std::vector<int>> data;

    // Burn in flips
    for (int i = 0; i < Nburn; i++) {
      CreateCluster();
      FlipCluster();
    }

    std::cout << "Burn-in completed." << std::endl;

    // Sample flips
    for (int i = 0; i < Nsamples; i++) {
      for (int j = 0; j < Ncorr; j++) {
        CreateCluster();
        FlipCluster();
      }
      data.push_back(confbin_);

      if (i % (Nsamples / 10) == 0) {
        std::cout << "i=" << i << std::endl;
      }
    }

    return data;
  }
};
