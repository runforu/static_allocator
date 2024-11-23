#ifndef STATIC_ALLOCATOR_H
#define STATIC_ALLOCATOR_H
#include <array>
#include <atomic>
#include <cstdarg>
#include <cstddef>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

// direct memory allocation
namespace dma {

// there are 8192 elements by default for the same size.
static constexpr uint16_t kMaxCellCount = 8192;

template<uint16_t MaxCellCount = kMaxCellCount>
class FastPage {
    static constexpr uint64_t kAllBitOne = 0xFFFFFFFFFFFFFFFF;
    static constexpr uint64_t kOne = 1;
    static constexpr uint64_t kMaskSize = MaxCellCount / 64;

public:
    FastPage(uint16_t byte_per_cell, uint16_t cell_count) : byte_per_cell_(byte_per_cell), cell_count_(cell_count) {
        memory_ = std::unique_ptr<uint8_t[]>(new uint8_t[byte_per_cell * cell_count] {0});
        flags_ = std::unique_ptr<std::atomic_flag[]>(new std::atomic_flag[cell_count]);
        // seems no need to clear.
        for (uint16_t i = 0; i < cell_count_; ++i) {
            flags_[i].clear(std::memory_order_relaxed);
        }

        // setup bitmap
        valid_mask_size_ = Divide64(cell_count_ - static_cast<uint16_t>(1)) + static_cast<uint16_t>(1);
        auto mask = kAllBitOne;
        std::fill_n(masks_.begin(), valid_mask_size_, mask);
        if (auto shift = Modulo64(cell_count_)) {
            masks_[valid_mask_size_ - 1] = ~(mask << shift);
        }
    }

    FastPage(const FastPage& page) = delete;
    FastPage(FastPage&& page) = delete;
    FastPage operator=(const FastPage& page) = delete;
    FastPage operator=(FastPage&& page) = delete;
    ~FastPage() noexcept = default;

    static uint16_t Divide64(uint16_t value) noexcept {
        static constexpr unsigned int kDivide64Bits = 6;
        return static_cast<uint16_t>(value >> kDivide64Bits);
    }

    static uint16_t Multiply64(uint16_t value) noexcept {
        static constexpr unsigned int kMultiply64Bits = 6;
        return static_cast<uint16_t>(value << kMultiply64Bits);
    }

    static uint16_t Modulo64(uint16_t value) noexcept {
        static constexpr unsigned int kModulo64Mask = 0x3f;
        return static_cast<uint16_t>(value & kModulo64Mask);
    }

    static uint16_t GetFirstOne(uint64_t offset) noexcept {
        // clang-format off
        // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers, readability-magic-numbers)
        static constexpr std::array<uint8_t, 256> bit_map {
            /* 00 */ 0, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            /* 10 */ 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            /* 20 */ 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            /* 30 */ 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            /* 40 */ 6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            /* 50 */ 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            /* 60 */ 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            /* 70 */ 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            /* 80 */ 7, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            /* 90 */ 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            /* A0 */ 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            /* B0 */ 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            /* C0 */ 6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            /* D0 */ 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            /* E0 */ 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
            /* F0 */ 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
        };
        // clang-format on

        uint16_t value = 0;
        if ((offset & 0xFFFFFFFFU) == 0) {
            value += 32U;
            offset >>= 32U;
        }
        if ((offset & 0xFFFFU) == 0) {
            value += 16U;
            offset >>= 16U;
        }
        if ((offset & 0xFFU) == 0) {
            value += 8U;
            offset >>= 8U;
        }
        return bit_map[offset & 0xFFU] + value;
        // NOLINTEND(cppcoreguidelines-avoid-magic-numbers, readability-magic-numbers)
    }

    void* Allocate(std::size_t /*size*/) noexcept {
        static uint16_t next_ {};
        auto mask_position = next_++;
        if (next_ >= valid_mask_size_) {
            next_ = 0;
        }
        if (mask_position >= valid_mask_size_) {
            mask_position = 0;
        }

        for (uint16_t i = 0; i < valid_mask_size_; ++i) {
            auto mask = masks_[mask_position].load(std::memory_order_acquire);
            while (mask != 0) {
                auto lowest_one_index = GetFirstOne(mask);
                auto index = lowest_one_index + Multiply64(mask_position);
                if (!flags_[index].test_and_set(std::memory_order_acquire)) {
                    masks_[mask_position].fetch_and(~(kOne << lowest_one_index), std::memory_order_release);
                    return &memory_[index * byte_per_cell_];
                }
                mask &= (~(kOne << lowest_one_index));  // next item
            }
            ++mask_position;
            if (mask_position >= valid_mask_size_) {
                mask_position = 0;
            }
        }
        return nullptr;
    }

    bool Release(void* pointer) noexcept {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
        auto index = static_cast<uint16_t>(reinterpret_cast<uint8_t*>(pointer) - memory_.get()) / byte_per_cell_;
        auto mask_position = Divide64(index);
        auto offset_in_mask = Modulo64(index);
        masks_[mask_position].fetch_or(kOne << offset_in_mask, std::memory_order_release);
        flags_[index].clear(std::memory_order_release);
        return true;
    }

    uint16_t BytePerCell() {
        return byte_per_cell_;
    }

    uint16_t CellCount() noexcept {
        return cell_count_;
    }

    bool Contains(void* pointer) noexcept {
        const void* begin = memory_.get();
        const void* end = memory_.get() + static_cast<unsigned int>(byte_per_cell_ * cell_count_);
        return pointer >= begin && pointer < end;
    }

    // For diagnosing usage
    uint8_t* StartPoint() noexcept {
        return memory_.get();
    }

private:
    // 1 of the bit of 64 bits means the cell is not allocated, 0 means the cell is allocated.
    static constexpr unsigned int k64BitsAlignment = 8;
    alignas(k64BitsAlignment) std::array<std::atomic_uint_fast64_t, kMaskSize> masks_;
    alignas(k64BitsAlignment) std::unique_ptr<uint8_t[]> memory_;
    alignas(k64BitsAlignment) std::unique_ptr<std::atomic_flag[]> flags_;

    uint16_t byte_per_cell_ {};
    uint16_t cell_count_ {};
    uint16_t valid_mask_size_ {};
};

class StringMemoryTraits {
private:
    // Not exposing StringMemoryTraits::Iterator
    using Iterator = std::vector<std::unique_ptr<FastPage<>>>::iterator;

public:
    static StringMemoryTraits& Instance() noexcept {
        static StringMemoryTraits instance;
        return instance;
    }

    StringMemoryTraits(const StringMemoryTraits&) = delete;
    StringMemoryTraits(StringMemoryTraits&&) = delete;
    StringMemoryTraits operator=(const StringMemoryTraits&) = delete;
    StringMemoryTraits operator=(StringMemoryTraits&&) = delete;

    ~StringMemoryTraits() = default;

    std::size_t MaxSize() noexcept {
        return pages_.back()->BytePerCell();
    }

    Iterator begin() {
        return pages_.begin();
    }

    Iterator end() {
        return pages_.end();
    }

private:
    // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers, readability-magic-numbers)
    StringMemoryTraits() {
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {32, 2048}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {48, 2048}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {64, 2048}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {96, 1024}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {128, 1024}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {160, 1024}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {192, 1024}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {224, 1024}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {256, 1024}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {320, 768}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {384, 768}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {448, 512}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {512, 256}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {448, 256}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {512, 256}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {640, 160}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {768, 64}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {896, 64}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {1024, 32}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {1280, 24}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {1536, 16}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {1792, 16}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {2048, 16}));
        pages_.emplace_back(std::unique_ptr<FastPage<>>(new FastPage<> {4096, 8}));
    }
    // NOLINTEND(cppcoreguidelines-avoid-magic-numbers, readability-magic-numbers)

    std::vector<std::unique_ptr<FastPage<>>> pages_;
};

// For Vector, there are 128 elements for the same size.
static constexpr uint16_t kMaxCellCountForVector = 128;
template<uint16_t MaxCellCountForVector = kMaxCellCountForVector>
class VectorMemoryTraits {
private:
    // Not exposing VectorMemoryTraits::Iterator VectorMemoryTraits::VectorFastPage
    using VectorFastPage = FastPage<MaxCellCountForVector>;
    using Iterator = typename std::vector<std::unique_ptr<VectorFastPage>>::iterator;

public:
    static VectorMemoryTraits& Instance() noexcept {
        static VectorMemoryTraits instance;
        return instance;
    }

    VectorMemoryTraits(const VectorMemoryTraits&) = delete;
    VectorMemoryTraits(VectorMemoryTraits&&) = delete;
    VectorMemoryTraits operator=(const VectorMemoryTraits&) = delete;
    VectorMemoryTraits operator=(VectorMemoryTraits&&) = delete;

    ~VectorMemoryTraits() = default;

    std::size_t MaxSize() noexcept {
        return pages_.back()->BytePerCell();
    }

    Iterator begin() {
        return pages_.begin();
    }

    Iterator end() {
        return pages_.end();
    }

private:
    VectorMemoryTraits() noexcept {
        static constexpr uint16_t kMaxCellSize = 2048;
        static constexpr uint16_t kCellCountDecrement = 2;
        static constexpr uint16_t kStartingSize = 64;
        static constexpr uint16_t kSizeIncrement = 64;
        for (uint16_t i = kStartingSize, j = MaxCellCountForVector; i <= kMaxCellSize;
             i += kSizeIncrement, j -= kCellCountDecrement) {
            pages_.emplace_back(new VectorFastPage {i, j});
        }
    }

    std::vector<std::unique_ptr<VectorFastPage>> pages_;
};

template<typename MemoryTraits>
class MemoryHelper {
public:
    static MemoryHelper& Instance() {
        static MemoryHelper instance;
        return instance;
    }

    MemoryHelper(const MemoryHelper&) = delete;
    MemoryHelper(MemoryHelper&&) = delete;
    MemoryHelper operator=(const MemoryHelper&) = delete;
    MemoryHelper operator=(MemoryHelper&&) = delete;
    ~MemoryHelper() noexcept = default;

    void* allocate(std::size_t size) noexcept {
        for (auto& page : pages_) {
            if (page->BytePerCell() >= size) {
                return page->Allocate(size);
            }
        }
        return nullptr;
    }

    void deallocate(void* const pointer, std::size_t /* size */) noexcept {
        for (auto& page : pages_) {
            if (page->Contains(pointer) && page->Release(pointer)) {
                break;
            }
        }
    }

    std::size_t MaxSize() noexcept {
        return pages_.MaxSize();
    }

private:
    MemoryHelper() = default;
    MemoryTraits& pages_ = MemoryTraits::Instance();
};

template<typename MemoryTraits>
class MemoryHelperWithLog {
public:
    static MemoryHelperWithLog& Instance() {
        static MemoryHelperWithLog instance;
        return instance;
    }

    MemoryHelperWithLog(const MemoryHelperWithLog&) = delete;
    MemoryHelperWithLog(MemoryHelperWithLog&&) = delete;
    MemoryHelperWithLog operator=(const MemoryHelperWithLog&) = delete;
    MemoryHelperWithLog operator=(MemoryHelperWithLog&&) = delete;
    ~MemoryHelperWithLog() noexcept = default;

    void* allocate(std::size_t size) noexcept {
        for (auto& page : pages_) {
            if (page->BytePerCell() >= size) {
                if (auto* pointer = page->Allocate(size)) {
                    statistics_.OnAllocate(page->BytePerCell());
                    return pointer;
                }
            }
        }
        statistics_.OnNoMemory(size);
        return nullptr;
    }

    void deallocate(void* const pointer, std::size_t /* size */) noexcept {
        for (auto& page : pages_) {
            if (page->Contains(pointer) && page->Release(pointer)) {
                statistics_.OnRelease(page->BytePerCell());
                break;
            }
        }
    }

    std::size_t MaxSize() noexcept {
        return pages_.MaxSize();
    }

private:
    struct Statistics {
        void OnAllocate(std::size_t size) {
            auto& flag = flags_[size];
            while (flag.test_and_set(std::memory_order_acquire)) {
                ;
            }
            auto& info = informations_[size];
            ++info.alloctions_;
            auto value = info.alloctions_ - info.releases_;
            if (value > info.max_load_) {
                info.max_load_ = value;
            }
            flag.clear(std::memory_order_release);
        }

        void OnRelease(std::size_t size) {
            auto& flag = flags_[size];
            while (flag.test_and_set(std::memory_order_acquire)) {
                ;
            }
            auto& info = informations_[size];
            ++info.releases_;
            auto value = info.alloctions_ - info.releases_;
            if (value > info.max_load_) {
                info.max_load_ = value;
            }
            flag.clear(std::memory_order_release);
        }

        void OnNoMemory(std::size_t size) {
            auto& flag = flags_[size];
            while (flag.test_and_set(std::memory_order_acquire)) {
                ;
            }
            auto& info = informations_[size];
            ++info.no_memory_;
            flag.clear(std::memory_order_release);
        }

        ~Statistics() {
            int index = 0;
            for (auto& info : informations_) {
                if (info.alloctions_ != 0 || info.releases_ != 0 || info.max_load_ != 0) {
                    // Show the statistics of memory usage when program ends.
                    std::cout << "Cell size: " << index << " \t\tallocate times: " << info.alloctions_
                              << " \t\t release times: " << info.releases_ << " \t\t max load: " << info.max_load_
                              << " \t\t no memory: " << info.no_memory_ << std::endl;
                }
                ++index;
            }
        }

        Statistics() noexcept {
            for (auto& flag : flags_) {
                flag.clear(std::memory_order_relaxed);
            }
        }

        Statistics(const Statistics& statistics) = delete;
        Statistics(Statistics&& statistics) noexcept = delete;
        Statistics& operator=(const Statistics& statistics) = delete;
        Statistics& operator=(Statistics&& statistics) noexcept = delete;

    private:
        struct Information {
            uint64_t alloctions_;
            uint64_t releases_;
            uint64_t max_load_;
            uint64_t no_memory_;
        };
        static constexpr unsigned int kMaximumAllocations = 4096;
        std::array<Information, kMaximumAllocations> informations_ {};
        std::array<std::atomic_flag, kMaximumAllocations> flags_ {};
    };

    MemoryHelperWithLog() = default;
    MemoryTraits& pages_ = MemoryTraits::Instance();
    Statistics statistics_;
};

}  // namespace dma

template<typename Type, typename Memory>
class Allocator {
public:
    using value_type = Type;
    using size_type = size_t;
    using difference_type = ptrdiff_t;
    using pointer = Type*;
    using const_pointer = const Type*;
    using reference = Type&;
    using const_reference = const Type&;

    template<class U>
    struct rebind {
        using other = Allocator<U, Memory>;
    };

    constexpr Allocator() noexcept = default;
    Allocator(const Allocator&) noexcept = default;
    Allocator(Allocator&&) noexcept = default;
    Allocator& operator=(const Allocator&) noexcept = default;
    Allocator& operator=(Allocator&&) noexcept = default;
    ~Allocator() noexcept = default;
    template<typename Other>
    explicit Allocator(Allocator<Other, Memory> const& other) noexcept {}

    Type* allocate(std::size_t size) noexcept {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
        return reinterpret_cast<Type*>(Memory::Instance().allocate(size * sizeof(Type)));
    }

    void deallocate(Type* const pointer, std::size_t size) noexcept {
        Memory::Instance().deallocate(pointer, size);
    }

    Type* allocate(std::size_t size, const void* /* unused */) noexcept {
        return allocate(size);
    }

    void construct(Type* pointer, const Type& value) noexcept {
        new (pointer) Type(value);
    }

    void destroy(Type* pointer) noexcept {
        pointer->~Type();
    }

    std::size_t max_size() const noexcept {
        return Memory::Instance().MaxSize();
    }

    template<class Other>
    bool operator==(const Allocator<Other, Memory>& /* unused */) const noexcept {
        return true;
    }

    template<class Other>
    bool operator!=(const Allocator<Other, Memory>& /* unused */) const noexcept {
        return false;
    }
};

#endif  // STATIC_ALLOCATOR_H
