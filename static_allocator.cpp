#include "static_allocator.h"

namespace {
// Force initialize memory at the begining of program. FastPage<n> is allocated before any usage.
const auto& STRING_PAGES = dma::StringMemoryTraits::Instance();
const auto& VECTOR_PAGES = dma::VectorMemoryTraits<>::Instance();
}  // namespace
