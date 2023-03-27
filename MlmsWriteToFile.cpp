/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include "MlmsWriteToFile.h"
#include "Mlms.h"

void writeToFile(const matrix &arr, const std::string name) {
  std::ofstream file;
  file.open(name, std::ofstream::out | std::ofstream::trunc);
  for (std::size_t i = 0; i < arr.shape[0]; i++) {
    for (std::size_t j = 0; j < arr.shape[1]; j++) {
      file << arr(i, j);
      file << "\t";
    }
    file << std::endl;
  }
  file.close();
}
