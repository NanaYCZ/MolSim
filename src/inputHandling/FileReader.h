#pragma once

#include <list>
#include <optional>

#include "particleModel/updating/CellCalculator.h"


class FileReader {
 public:
  FileReader();
  virtual ~FileReader();

  struct SphereData;

  struct CuboidData;

  struct MembraneData;

  struct ProgramArgs; 


  /**
   * @brief read all program arguments
   * 
   * Reads all the Arguments that are specified within the ProgramArgs struct.
   * That means delta_t and end_t for the simulation, how often output files
   * should be written, the basename of these files and lastly an arbitrary
   * amount of cuboids and then an arbitrary amount of spheres. This data
   * should be specified in a certain xsd format (that is specified within 
   * 'project_root_dir'/input/parameters.xsd) and is then read into a struct that
   * is returned
   * 
   * 
   * @param filename name from which the arguments are read (has to be in xml format)
   * 
   * @return struct contains all necessary info for the program
   * 
  */
  ProgramArgs readProgramArguments(std::string filename);


  static void initializeCorrectInitialTemp(FileReader::ProgramArgs& args);


};
