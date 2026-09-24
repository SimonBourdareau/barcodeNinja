// FastqReadPool.h

#pragma once

#include <vector>
#include <memory>

// Template class for a memory pool
template <typename T>
class ObjectPool {
private:
    std::vector<std::unique_ptr<T>> pool;  // Pool of pre-allocated objects
    std::size_t nextAvailable = 0;          // Index of the next available object in the pool

public:
    /**
     * @brief Constructs a memory pool with a specified number of objects.
     *
     * @param size The number of objects to pre-allocate.
     */
    ObjectPool(std::size_t size) {
        pool.reserve(size);
        for (std::size_t i = 0; i < size; ++i) {
            pool.push_back(std::make_unique<T>());
        }
    }

    /**
     * @brief Allocates an object from the pool.
     *
     * @return A pointer to the allocated object.
     * @throws std::bad_alloc if no objects are available in the pool.
     */
    T* allocate() {
        if (nextAvailable >= pool.size()) {
            throw std::bad_alloc();  // No more objects available in the pool
        }
        return pool[nextAvailable++].get();
    }

    /**
     * @brief Deallocates an object and returns it to the pool.
     *
     * @param ptr A pointer to the object to deallocate.
     */
    void deallocate(T* ptr) {
        if (ptr == nullptr) {
            return;  // Ignore null pointers
        }

        // Reset the object for reuse
        ptr->~T();
        new (ptr) T();  // Reconstruct the object in-place

        // Mark the object as available
        --nextAvailable;
    }

    /**
     * @brief Resets the entire pool, marking all objects as available.
     */
    void resetPool() {
        nextAvailable = 0;  // Reset the index to the start of the pool
    }

    /**
     * @brief Returns the size of the pool.
     *
     * @return The size of the pool.
     */
    std::size_t  size() {
        return pool.size();
    }

    /**
     * @brief Returns the number of available objects in the pool.
     *
     * @return The number of available objects.
     */
    std::size_t available() const {
        return pool.size() - nextAvailable;
    }

    /**
     * @brief Clears the pool, freeing memory.
     */
    void clearPool() {
      // Iterate through the pool and deallocate each object
      for (auto& ptr : pool) {
        deallocate(ptr.get());
      }
      pool.clear();
      nextAvailable = 0;
    }
};
