#pragma once
// Non-convex models of the samples (teapot, bunny, ring), decomposed once into convex parts and shared
// by every body that uses them.
#include "rigid/Shapes.h"

#include <memory>

namespace rf {

std::shared_ptr<const CompoundShape> teapotShape();
std::shared_ptr<const CompoundShape> bunnyShape();
// A steel ring of a chain: 10 cm across, its tube 24 mm thick.
std::shared_ptr<const CompoundShape> ringShape();

} // namespace rf
