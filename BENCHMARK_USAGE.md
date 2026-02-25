# Benchmarking atn_add_n Implementations

This document explains how to use the benchmarking functionality to compare the performance of `atn_add_n` (standard) and `atn_add_n_intrinsic` (compiler intrinsics) implementations.

## Usage

### From C++ Code

Include the header and call the benchmark function:

```cpp
#include "Athena.hpp"

int main()
{
    // Benchmark with default parameters (100 limbs, 10000 iterations)
    Athena::benchmark_add_n_implementations();
    
    // Or customize the parameters
    // benchmark_add_n_implementations(limbCount, iterations)
    Athena::benchmark_add_n_implementations(50, 100000);   // 50 limbs, 100k iterations
    Athena::benchmark_add_n_implementations(200, 5000);     // 200 limbs, 5k iterations
    
    return 0;
}
```

### From Sandbox/main.cpp

You can add the benchmark call to your existing main.cpp:

```cpp
#include "Athena.hpp"
#include <iostream>

int main()
{
    std::cout << "Running benchmark..." << std::endl;
    
    // Test with different sizes
    Athena::benchmark_add_n_implementations(10, 100000);    // Small: 10 limbs
    Athena::benchmark_add_n_implementations(100, 10000);    // Medium: 100 limbs
    Athena::benchmark_add_n_implementations(1000, 1000);    // Large: 1000 limbs
    
    return 0;
}
```

## Output Example

```
=== Benchmarking atn_add_n implementations ===
Limb count: 100
Iterations: 10000
Limb size: 64 bits

atn_add_n (standard):
  Average: 245.32 ns
  Min:     198.00 ns
  Max:     3421.00 ns

atn_add_n_intrinsic:
  Average: 187.45 ns
  Min:     154.00 ns
  Max:     2893.00 ns

Performance comparison:
  Intrinsic version is 1.31x faster
  Speedup: 30.9%

Correctness: Both implementations produce the same carry value.
```

## Parameters

- **limbCount**: Number of limbs (multi-precision integers) to add. 
  - Each limb is typically 32 or 64 bits depending on your system.
  - Larger values test performance with bigger numbers.
  - Default: 100

- **iterations**: Number of times to repeat the operation for averaging.
  - Higher values give more accurate results but take longer.
  - Default: 10000

## Tips

1. **Warm-up**: Run the benchmark a few times to warm up the CPU cache.
2. **Consistent Environment**: Close other applications for more consistent results.
3. **Multiple Sizes**: Test with different limb counts to see how performance scales.
4. **Release Build**: For accurate performance testing, use Release build configuration with optimizations enabled.
5. **Compiler Intrinsics**: The intrinsic version uses:
   - MSVC: `_addcarry_u64` / `_addcarry_u32`
   - GCC/Clang: `__builtin_add_overflow`

## Enabling Optimizations

### Quick Method (Visual Studio):
1. Change configuration dropdown to **`x64-release`**
2. Build solution

### Command Line:
```sh
# Configure for release
cmake --preset x64-release

# Build
cmake --build out/build/x64-release --config Release
```

### Optimization Flags Applied:
**MSVC (Release):**
- `/O2` - Maximum optimization (speed)
- `/Ob2` - Inline function expansion
- `/Oi` - Enable intrinsic functions
- `/Ot` - Favor fast code
- `/GL` - Whole program optimization
- `/LTCG` - Link-time code generation
- `/arch:AVX2` - AVX2 instructions (if supported)

**GCC/Clang (Release):**
- `-O3` - Aggressive optimizations
- `-march=native` - Optimize for your CPU
- `-DNDEBUG` - Disable assertions

### Verify Optimizations:
Run the benchmark - it will show:
```
Build: RELEASE (optimized)
Compiler: MSVC 1930
```

If it shows "DEBUG (not optimized)", you're not building in Release mode.
