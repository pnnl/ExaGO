/**
 * @file testBase.hpp
 *
 * @author Asher Mancinelli <asher.mancinelli@pnnl.gov>, PNNL
 * @author Slaven Peles <slaven.peles@pnnl.gov>, PNNL
 *
 */
#pragma once

/*
  When checking for nulls, we often want to throw the same exception each time.
  We are using a macro to avoid changing the error message in many places should the need
  arise.
*/
#define THROW_NULL_DEREF throw std::runtime_error("Attempted null dereference")

#include <iostream>
#include <limits>
#include <cmath>
#include <type_traits>

namespace exago
{

namespace tests
{

using LocalOrdinalType = int;
using GlobalOrdinalType = long long int;
using RealType = double;

/**
  @brief Base class for all testing classes. Each child class will call the
  static methods to report failures/successes.
 */
class TestBase
{

  static constexpr RealType zero = 0.0;
  static constexpr RealType quarter = 0.25;
  static constexpr RealType half = 0.5;
  static constexpr RealType one = 1.0;
  static constexpr RealType two = 2.0;
  static constexpr RealType three = 3.0;
  static constexpr RealType eps =
      10*std::numeric_limits<RealType>::epsilon();
  static constexpr LocalOrdinalType SKIP_TEST = -1;

  /*
   * must be const pointer and const dest for
   * const string declarations to pass
   * -Wwrite-strings
   */
  static constexpr const char * const  RED       = "\033[1;31m";
  static constexpr const char * const  GREEN     = "\033[1;32m";
  static constexpr const char * const  YELLOW    = "\033[1;33m";
  static constexpr const char * const  CLEAR     = "\033[0m";

protected:

  /**
   * @brief Returns true if the value agrees with the reference to within
   * relative tolerance. Uses absolute tolerancew when the reference is zero.
   * 
   */
  [[nodiscard]]
  static
  LocalOrdinalType isEqual(const RealType value, const RealType reference)
  {
    return (std::abs(value - reference)/(one + std::abs(reference)) < eps);
  }

  /// Prints error output for each rank
  static void printMessage(const LocalOrdinalType fail, const char* funcname, const LocalOrdinalType rank)
  {
    if(fail > 0)
    {
      std::cout << RED << "--- FAIL: Test " << funcname << " on rank " << rank << CLEAR << "\n";
    }
    else if (fail == SKIP_TEST)
    {
      if(rank == 0)
      {
        std::cout << YELLOW << "--- SKIP: Test " << funcname << CLEAR << "\n";
      }
    }
    else
    {
      if(rank == 0)
      {
        std::cout << GREEN << "--- PASS: Test " << funcname << CLEAR << "\n";
      }
    }
  }

};

}} // namespace exago::tests
