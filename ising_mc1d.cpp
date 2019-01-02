#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

class Ising {
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

int main() {
  int N, Nsamples, Ncorr, Nburn, n_temperatures;
  float T;
  std::vector<float> Tlist;
  Ising sampler;
  std::ofstream file;

  // User gives parameters
  std::cout << "N = ";
  std::cin >> N;
  std::cout << "n_temperatures = ";
  std::cin >> n_temperatures;
  for (int i = 0; i < n_temperatures; i++) {
    std::cout << "T = ";
    std::cin >> T;
    Tlist.push_back(T);
  }
  std::cout << "Nsamples = ";
  std::cin >> Nsamples;
  std::cout << "Ncorrelation = ";
  std::cin >> Ncorr;
  std::cout << "Nburn = ";
  std::cin >> Nburn;

  // Sample
  std::vector<std::vector<std::vector<int>>> data;
  for (int iT = 0; iT < n_temperatures; iT++) {
    sampler.Init(N, 1.0 / Tlist[iT]);
    data.push_back(sampler.Run(Nsamples, Ncorr, Nburn));
  }

  // Write samples to file
  file.open("confs.txt");
  for (std::size_t i = 0; i < data.size(); i++) {
    for (std::size_t j = 0; j < data[i].size(); j++) {
      for (int k = 0; k < N; k++) {
        file << data[i][j][k];
      }
      file << "\n";
    }
  }
  file.close();

  // Write params to file
  file.open("params.txt");
  file << "N=" << N << "\n";
  for (int iT = 0; iT < n_temperatures; iT++) {
    file << "T=" << Tlist[iT] << ", ";
  }
  file << "\n";
  file << "Nsamples=" << Nsamples << "\n";
  file << "Ncorrelation=" << Ncorr << "\n";
  file << "Nburn=" << Nburn << "\n";

  return 0;
}
