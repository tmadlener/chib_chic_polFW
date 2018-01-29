#include "utils.h"

#include <fstream>

#include <sys/types.h>
#include <sys/stat.h>


//#define WINDOWS

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif
#include<iostream>


// Function implementations

bool file_exists(const std::string &fileName)
{
  std::ifstream infile(fileName);
  return infile.good();
}

bool is_directory(const std::string &dirName)
{
  struct stat st;
  return (stat(dirName.c_str(), &st) == 0) && (st.st_mode & S_IFDIR);
}

std::string get_current_workingdir() {
  char buff[FILENAME_MAX];
  GetCurrentDir(buff, FILENAME_MAX);
  return std::string(buff);
}


