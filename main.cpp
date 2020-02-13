#include "Axis.h"
#include "Grid.h"
#include "MetaDynamics.h"

#include <string>

int main(int argc, char* argv[]) {
    MetaDynamics meta;
    meta.readFromJson(argv[1]);
    meta.compute();
    meta.writeToFile("result.pmf");
    meta.writeGradients("result.grad");
    meta.writeFakeCount("result.count");
    return 0;
}
