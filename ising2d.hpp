#include <cmath>
#include <vector>

class Ising2D {
  // Number of sites
  int N_, L_;
  // Coefficient
  float K_;
  // Current configuration
  std::vector<int> conf_, confbin_;
  // Nearest neighbor map (right, left, down, up)
  std::vector<std::vector<int>> nnmap_;
  // Current energy
  float energy_;

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

 public:
  int Nsites() { return N_; }

  void Init(const int L, const float K) {
    // Set random seed
    srand(time(NULL));

    // Initialize configuration and nn map
    L_ = L;
    N_ = L * L;
    K_ = K;
    conf_.push_back(2 * (rand() % 2) - 1);
    confbin_.push_back(ToBinary(conf_[0]));
    nnmap_.push_back(FindNeighbors(0));
    for (int i = 1; i < N_; i++) {
      conf_.push_back(2 * (rand() % 2) - 1);
      confbin_.push_back(ToBinary(conf_[i]));
      nnmap_.push_back(FindNeighbors(i));
    }

    // Calculate initial energy
    // (add only right (0) and down (2) interactions
    // to avoid repetition)
    energy_ = 0;
    for (int i = 0; i < N_; i++) {
      energy_ += -K_ * conf_[i] * (conf_[nnmap_[i][0]] + conf_[nnmap_[i][2]]);
    }
  }

  void SingleFlip() {
    // Generate random site to flip
    int site = rand() % N_;
    // Calculate difference in energy
    float denergy = 0.0;
    for (std::size_t j = 0; j < nnmap_[site].size(); j++) {
      denergy += conf_[nnmap_[site][j]];
    }
    denergy *= 2.0 * K_ * conf_[site];

    // Random float in [0, 1]
    float r = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    // Check acceptance
    if (exp(-denergy) > r) {
      conf_[site] *= -1;
      confbin_[site] = ToBinary(conf_[site]);
      energy_ += denergy;
    }
  }

  std::vector<std::vector<int>> Run(const int Nsamples, const int Ncorr,
                                    const int Nburn) {
    std::vector<std::vector<int>> data;

    // Burn in sweeps
    for (int i = 0; i < Nburn; i++) {
      SingleFlip();
    }

    // Sample sweeps
    for (int i = 0; i < Nsamples; i++) {
      for (int j = 0; j < Ncorr; j++) {
        SingleFlip();
      }
      data.push_back(confbin_);
    }

    return data;
  }
};
