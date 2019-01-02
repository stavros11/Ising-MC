#include <cmath>
#include <vector>

class Ising1D {
  // Number of sites
  int N_;
  // Coefficient
  float K_;
  // Current configuration
  std::vector<int> conf_, confbin_;
  // Current energy
  float energy_;

 public:
  void Init(const int N, const float K) {
    // Set random seed
    srand(time(NULL));

    // Initialize configuration and calculate energy
    N_ = N;
    K_ = K;
    energy_ = 0;
    conf_.push_back(2 * (rand() % 2) - 1);
    confbin_.push_back(ToBinary(conf_[0]));
    for (int i = 1; i < N_; i++) {
      conf_.push_back(2 * (rand() % 2) - 1);
      confbin_.push_back(ToBinary(conf_[i]));
      energy_ += -K_ * conf_[i - 1] * conf_[i];
    }
    // Add PBC term to energy
    energy_ += -K_ * conf_[0] * conf_[N_ - 1];
  }

  int ToBinary(const int x) {
    if (x == 1) {
      return 1;
    } else {
      return 0;
    }
  }

  void SingleFlip() {
    // Generate random site to flip
    int site = rand() % N_;
    // Calculate difference in energy
    float denergy = 2.0 * K_ *
                    (conf_[(site - 1) % N_] + conf_[(site + 1) % N_]) *
                    conf_[site];
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
