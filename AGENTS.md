# TBKOSTER Development Guidelines for Agentic Coding

This document provides essential information for agentic coding agents working on the TBKOSTER codebase, a Fortran-based tight-binding magnetic materials simulation package.

## Build System and Commands

### Basic Build Commands
```bash
# Configure build (from project root)
cmake -B build

# Build the project
cmake --build build

# Clean build artifacts
cmake --build build --target clean
cmake --build build --target distclean  # removes generated .f90 files

# Install (if configured)
cmake --build build --target install
```

### Advanced Build Options
```bash
# Release build with optimizations
cmake -DCMAKE_BUILD_TYPE=Release -B build

# Debug build with checks
cmake -DCMAKE_BUILD_TYPE=Debug -B build

# Enable optional features
cmake -DENABLE_OPENMP=ON -DENABLE_TESTS=ON -B build

# Intel compiler with MKL
BLA_VENDOR=Intel10_64lp_seq FC=ifort cmake -B build

# OpenBLAS with gfortran
BLA_VENDOR=OpenBLAS FC=gfortran cmake -B build
```

### Test Commands
```bash
# Run all tests
ctest --test-dir build

# Run specific test
ctest --test-dir build -R atom_test
ctest --test-dir build -R math_test

# Verbose test output
ctest --test-dir build --verbose

# Run test manually
./build/unit_tests/math/math_test
```

## Code Architecture

### Module Structure
- **Core library**: Static library in `src/` containing all computational modules
- **Main executable**: `TBKOSTER.x` - main simulation program
- **Tools**: Post-processing executables in `tools/` (bands.x, pdos.x, mft.x, etc.)
- **Unit tests**: Individual test programs in `unit_tests/` subdirectories

### Key Modules
- `precision_mod`: Centralized numerical precision control
- `calculation_mod`: Main calculation orchestrator
- `atom_mod`, `element_mod`: Atomic and elemental properties
- `hamiltonian_tb_mod`: Tight-binding Hamiltonian implementation
- `mesh_mod`: k-point mesh generation and management
- `energy_mod`, `forces_mod`: Energy and force calculations

## Code Style Guidelines

### Module Organization
```fortran
! Standard module structure with consistent naming
module example_mod
   use, intrinsic :: iso_fortran_env
   use other_required_mod, only: specific_routines
   implicit none
   private
   
   ! Public interface
   public :: public_subroutine, public_function
   
   ! Private implementations
   private :: internal_helper
   
contains
   
subroutine public_subroutine(arg1, arg2)
   ! Implementation
end subroutine public_subroutine

end module example_mod
```

### Naming Conventions
- **Modules**: End with `_mod` suffix (e.g., `precision_mod`)
- **Subroutines/Functions**: Lowercase with underscores (`calculate_energy`)
- **Variables**: Descriptive names, lowercase (`n_atoms`, `k_vector`)
- **Constants**: All uppercase with parameters (`PI`, `BOHR_RADIUS`)
- **Types**: Lowercase ending with `_t` (`atom_type`, `lattice_type`)

### Precision Management
Always use the centralized precision module:
```fortran
use precision_mod, only: rp, ip  ! real and integer precision
real(rp) :: variable_name       ! Working precision real
integer(ip) :: count_var        ! Working precision integer
```

### Import Statement Guidelines
```fortran
! Preferred order and grouping
use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
use precision_mod, only: rp, ip, rl
use math_mod, only: matrix_multiply, cross_product
use string_mod, only: int2str, lower_case
```

### Documentation Style
Use FORD-compatible comments:
```fortran
!> Main calculation orchestrator for TBKOSTER simulations
!> 
!> This module coordinates the overall workflow including input reading,
!> calculation execution, and output generation. It manages the
!> different calculation types (SCF, MD, SD) and their
!> associated post-processing steps.
!> 
!> @author Cyrille Barreteau, Mathieu Cesar, Pascal Thibaudeau
!> @date 2020
```

### Error Handling
```fortran
! Check file existence before operations
inquire(file=filename, exist=file_exists)
if(.not. file_exists) then
   write(error_unit, '(a)') 'Error: File not found: ' // trim(filename)
   stop
endif

! Validate array bounds
if(size(array) /= expected_size) then
   write(error_unit, '(a)') 'Error: Array size mismatch'
   stop
endif
```

### Array and Memory Management
```fortran
! Automatic arrays preferred over allocatable when size known
real(rp) :: temp_array(n_atoms)

! Explicit allocation with error checking
allocate(large_array(n_points), stat=alloc_stat)
if(alloc_stat /= 0) then
   write(error_unit, '(a)') 'Memory allocation failed'
   stop
endif

! Always deallocate in reverse order of allocation
if(allocated(large_array)) deallocate(large_array)
```

### Compiler Directives and Preprocessor
```fortran
! Use preprocessor for optional features
#if defined(OpenMP_Fortran_FOUND)
   use omp_lib
   !$omp parallel do
#endif

#if defined(Enable_flush)
   call TBKOSTER_flush(output_unit)
#endif
```

## File I/O Conventions

### Input File Format
- Main input: `in_master.txt` in namelist format
- Configuration files: `in_*.txt` in subdirectories
- Use consistent namelist structure with clear variable naming

### Output File Format
- Main output: `out_log.txt` for execution log
- Results: `out_*.txt` in namelist format for reusability
- Tools: Generate `.dat` files compatible with gnuplot/xmgrace

### File Path Handling
```fortran
use string_mod, only: sl, trim

character(len=sl) :: filename
filename = trim(input_dir) // '/' // trim(base_name)
```

## Testing Guidelines

### Unit Test Structure
```fortran
program test_module
   use module_under_test, only: routine_to_test
   implicit none
   
   ! Test data
   real(rp), parameter :: tolerance = 1.0e-12_rp
   
   call test_routine1()
   call test_routine2()
   
contains

subroutine test_routine1()
   ! Arrange
   real(rp) :: input_value, expected_result, actual_result
   
   ! Act
   input_value = 2.0_rp
   expected_result = 4.0_rp
   actual_result = routine_to_test(input_value)
   
   ! Assert
   if(abs(actual_result - expected_result) < tolerance) then
      print *, 'PASS: test_routine1'
   else
      print *, 'FAIL: test_routine1 - Expected:', expected_result, 'Got:', actual_result
   endif
end subroutine test_routine1
```

### Test Compilation
Each test has its own CMakeLists.txt and compiles to a standalone executable. Use CTest for batch execution.

## Performance Considerations

### Optimization Flags
- Release: `-Ofast` with bounds checking enabled for safety
- Debug: `-O0` with full checking (`-fcheck=all`)
- OpenMP parallelization over k-points when enabled

### Memory Efficiency
- Prefer automatic arrays over allocatable where possible
- Use `intent(in)`, `intent(out)`, `intent(inout)` explicitly
- Leverage BLAS/LAPACK for linear algebra operations

## Git and Development Workflow

### Commit Message Format
```
component: brief description

Detailed explanation of changes, including:
- What was modified
- Why it was necessary
- Any breaking changes or side effects

Resolves: issue-number (if applicable)
```

### Branch Strategy
- `main`: Stable releases
- `dev`: Development with experimental features
- Feature branches: `feature/description` or `fix/description`

## Platform-Specific Notes

### Intel MKL Integration
- Requires `MKLROOT` environment variable
- Supports both `lp64` and `ilp64` integer models
- Automatic library detection and linking

### OpenMP Usage
- Conditional compilation with `OpenMP_Fortran_FOUND` macro
- Parallelization focused on k-point loops
- Thread-safe implementations required

### Debugging Support
- Full runtime checking in Debug builds
- Stack traces on gfortran (`-fbacktrace`)
- Bounds checking always enabled (`-fcheck=all`)

This document should be updated as new features are added or coding standards evolve. All contributors should follow these guidelines to maintain code quality and consistency across the TBKOSTER project.