#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <fstream>

// based on https://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c

inline bool exists_test3 (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}