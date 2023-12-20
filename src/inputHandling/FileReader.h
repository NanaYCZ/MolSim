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

  /**
   * @brief Reads Cuboids of Particles from a file and returns a list of
   * CuboidData structs
   *
   * reads cuboids from a file specified by the given filename(file has to have
   * specific format). A cuboid in the file is just a list of the parameters
   * that are in the CuboidData struct. Turns the parameters given in the file
   * into a struct containing the parameters and creates a list of CuboidData
   * structs based on the read data.
   *
   * @param particleContainer reference to the ParticleContainer to add to
   * @param filename Filename of the file containing CuboidData
   *
   *
   * @return Returns a list of structs that contain the data of the Cuboids that
   * were read. Mainly used for testing and debugging purposes.
   *
   */
  std::list<CuboidData> readCuboidFile(std::string filename);
};
