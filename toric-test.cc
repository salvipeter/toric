#include "toric.hh"

int main(int argc, char **argv) {
  if (argc != 3 && argc != 4) {
    std::cerr << "Usage: " << argv[0] << " <input.trp> <output.obj> [resolution]" << std::endl;
    return 1;
  }
  size_t resolution = 50;
  if (argc == 4)
    resolution = std::atoi(argv[3]);

  Toric::load(argv[1]).tessellate(resolution).writeOBJ(argv[2]);
}
