/* Copyright 2023 <Ilia Fedotov @ Uni Freiburg> */
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "MlmsWriteToFile.h"

void writeVecToFile(std::vector<double> arr) {
  std::ofstream file;
  file.open("vectorResult.txt");
  int size = sqrt(arr.size());
  int row = 0;
  for (int i = 0; i < arr.size(); i++) {
    file << (i%size)*size+(i/size)*size;
  }
  file << "Writing this to a file." << std::endl;
  file.close();
  return;
}

void writeTimeToFile(std::vector<double> arr) {
  std::ofstream file;
  file.open("timeToCompute.txt");
  file << "Writing this to a file." << std::endl;
  file.close();
  return;
}
