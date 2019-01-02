#include <fstream>
#include <iostream>
#include <vector>
#include "ising1d.hpp"

int main() {
  int N, Nsamples, Ncorr, Nburn, n_temperatures;
  float T;
  std::vector<float> Tlist;
  Ising1D sampler;
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
