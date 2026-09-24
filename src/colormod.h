// colormod.h

#pragma once

#include <ostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>

namespace Color {

/**
 * @brief This file defines a set of utilities for handling color codes in terminal output.
 *
 * This header file provides an enum class `Code` representing color codes, a `Modifier` class
 * for applying color codes to output, and various utility functions and macros for creating
 * colored exceptions and combining color codes.
 *
 * AUTHOR : Simon BOURDAREAU
 * INSTITUTION : Stowers Institute for Medical Research
 * DATE : August 24
 * PORTABILITY : C++17
 */

/**
 * @brief Enum class representing color codes.
 */
enum class Code {
  // Standard foreground colors
  FG_BLACK             = 30,
    FG_RED               = 31,
    FG_GREEN             = 32,
    FG_YELLOW            = 33,
    FG_BLUE              = 34,
    FG_MAGENTA           = 35,
    FG_CYAN              = 36,
    FG_WHITE             = 37,
    FG_DEFAULT           = 39,
    FG_DEFAULT_NO_BOLD   = 0,

    // Bright foreground colors
    FG_BRIGHT_BLACK      = 90,
    FG_BRIGHT_RED        = 91,
    FG_BRIGHT_GREEN      = 92,
    FG_BRIGHT_YELLOW     = 93,
    FG_BRIGHT_BLUE       = 94,
    FG_BRIGHT_MAGENTA    = 95,
    FG_BRIGHT_CYAN       = 96,
    FG_BRIGHT_WHITE      = 97,

    // Standard background colors
    BG_BLACK             = 40,
    BG_RED               = 41,
    BG_GREEN             = 42,
    BG_YELLOW            = 43,
    BG_BLUE              = 44,
    BG_MAGENTA           = 45,
    BG_CYAN              = 46,
    BG_WHITE             = 47,
    BG_DEFAULT           = 49,

    // Bright background colors
    BG_BRIGHT_BLACK      = 100,
    BG_BRIGHT_RED        = 101,
    BG_BRIGHT_GREEN      = 102,
    BG_BRIGHT_YELLOW     = 103,
    BG_BRIGHT_BLUE       = 104,
    BG_BRIGHT_MAGENTA    = 105,
    BG_BRIGHT_CYAN       = 106,
    BG_BRIGHT_WHITE      = 107,

    // Text styles
    BOLD                 = 1,
    FAINT                = 2,
    ITALIC               = 3,
    UNDERLINE            = 4,
    BLINK                = 5,
    REVERSE              = 7,
    CONCEALED            = 8,
    STRIKETHROUGH        = 9,
};

/**
 * @brief Class representing a color modifier.
 */
class Modifier {
  std::vector<Code> codes;
public:
  /**
   * @brief Constructs a Modifier with a variadic list of codes.
   *
   * @tparam Codes The types of the codes.
   * @param pCodes The codes to initialize the Modifier with.
   */
  template<typename... Codes>
  Modifier(Codes... pCodes) : codes{pCodes...} {}

  /**
   * @brief Constructs a Modifier with a range of codes.
   *
   * @tparam Iterator The type of the iterator.
   * @param begin The beginning of the range.
   * @param end The end of the range.
   */
  template<typename Iterator>
  Modifier(Iterator begin, Iterator end) : codes(begin, end) {}

  /**
   * @brief Outputs the color codes to an output stream.
   *
   * @param os The output stream.
   * @param mod The Modifier to output.
   * @return std::ostream& The output stream.
   */
  friend std::ostream& operator<<(std::ostream& os, const Modifier& mod) {
    os << "\033[";
    for (size_t i = 0; i < mod.codes.size(); ++i) {
      os << static_cast<int>(mod.codes[i]);
      if (i < mod.codes.size() - 1) {
        os << ";";
      }
    }
    return os << "m";
  }

  /**
   * @brief Predefined color modifiers.
   */
  static const Modifier red;
  static const Modifier green;
  static const Modifier blue;
  static const Modifier defaultColor;
  static const Modifier reset;
  static const Modifier bgRed;
  static const Modifier bgGreen;
  static const Modifier bgBlue;
  static const Modifier bgDefault;
  static const Modifier bold;
  static const Modifier italic;

  static const Modifier black;
  static const Modifier yellow;
  static const Modifier magenta;
  static const Modifier cyan;
  static const Modifier white;
  static const Modifier brightRed;
  static const Modifier brightGreen;
  static const Modifier brightYellow;
  static const Modifier brightBlue;
  static const Modifier brightMagenta;
  static const Modifier brightCyan;
  static const Modifier brightWhite;

  static const Modifier bgBlack;
  static const Modifier bgYellow;
  static const Modifier bgMagenta;
  static const Modifier bgCyan;
  static const Modifier bgWhite;
  static const Modifier bgBrightBlack;
  static const Modifier bgBrightRed;
  static const Modifier bgBrightGreen;
  static const Modifier bgBrightYellow;
  static const Modifier bgBrightBlue;
  static const Modifier bgBrightMagenta;
  static const Modifier bgBrightCyan;
  static const Modifier bgBrightWhite;

  // New style modifiers
  static const Modifier faint;
  static const Modifier underline;
  static const Modifier blink;
  static const Modifier reverse;
  static const Modifier concealed;
  static const Modifier strikethrough;

  /**
   * @brief Gets the codes of the Modifier.
   *
   * @return const std::vector<Code>& The codes of the Modifier.
   */
  const std::vector<Code>& getCodes() const { return codes; }
};

// Define the static members inline
inline const Modifier Modifier::red{Code::FG_RED};
inline const Modifier Modifier::green{Code::FG_GREEN};
inline const Modifier Modifier::blue{Code::FG_BLUE};
inline const Modifier Modifier::defaultColor{Code::FG_DEFAULT};
inline const Modifier Modifier::reset{Code::FG_DEFAULT_NO_BOLD};
inline const Modifier Modifier::bgRed{Code::BG_RED};
inline const Modifier Modifier::bgGreen{Code::BG_GREEN};
inline const Modifier Modifier::bgBlue{Code::BG_BLUE};
inline const Modifier Modifier::bgDefault{Code::BG_DEFAULT};
inline const Modifier Modifier::bold{Code::BOLD};
inline const Modifier Modifier::italic{Code::ITALIC};

inline const Modifier Modifier::black{Code::FG_BLACK};
inline const Modifier Modifier::yellow{Code::FG_YELLOW};
inline const Modifier Modifier::magenta{Code::FG_MAGENTA};
inline const Modifier Modifier::cyan{Code::FG_CYAN};
inline const Modifier Modifier::white{Code::FG_WHITE};
inline const Modifier Modifier::brightRed{Code::FG_BRIGHT_RED};
inline const Modifier Modifier::brightGreen{Code::FG_BRIGHT_GREEN};
inline const Modifier Modifier::brightYellow{Code::FG_BRIGHT_YELLOW};
inline const Modifier Modifier::brightBlue{Code::FG_BRIGHT_BLUE};
inline const Modifier Modifier::brightMagenta{Code::FG_BRIGHT_MAGENTA};
inline const Modifier Modifier::brightCyan{Code::FG_BRIGHT_CYAN};
inline const Modifier Modifier::brightWhite{Code::FG_BRIGHT_WHITE};

inline const Modifier Modifier::bgBlack{Code::BG_BLACK};
inline const Modifier Modifier::bgYellow{Code::BG_YELLOW};
inline const Modifier Modifier::bgMagenta{Code::BG_MAGENTA};
inline const Modifier Modifier::bgCyan{Code::BG_CYAN};
inline const Modifier Modifier::bgWhite{Code::BG_WHITE};
inline const Modifier Modifier::bgBrightBlack{Code::BG_BRIGHT_BLACK};
inline const Modifier Modifier::bgBrightRed{Code::BG_BRIGHT_RED};
inline const Modifier Modifier::bgBrightGreen{Code::BG_BRIGHT_GREEN};
inline const Modifier Modifier::bgBrightYellow{Code::BG_BRIGHT_YELLOW};
inline const Modifier Modifier::bgBrightBlue{Code::BG_BRIGHT_BLUE};
inline const Modifier Modifier::bgBrightMagenta{Code::BG_BRIGHT_MAGENTA};
inline const Modifier Modifier::bgBrightCyan{Code::BG_BRIGHT_CYAN};
inline const Modifier Modifier::bgBrightWhite{Code::BG_BRIGHT_WHITE};

inline const Modifier Modifier::faint{Code::FAINT};
inline const Modifier Modifier::underline{Code::UNDERLINE};
inline const Modifier Modifier::blink{Code::BLINK};
inline const Modifier Modifier::reverse{Code::REVERSE};
inline const Modifier Modifier::concealed{Code::CONCEALED};
inline const Modifier Modifier::strikethrough{Code::STRIKETHROUGH};

/**
 * @brief Combines multiple color codes into a single Modifier.
 *
 * @tparam Args The types of the arguments.
 * @param args The arguments to combine.
 * @return Modifier The combined Modifier.
 */
template<typename... Args>
Modifier colorCombiner(Args&&... args) {
  std::vector<Code> codes;
  (processArg(codes, std::forward<Args>(args)), ...);
  return Modifier(codes.begin(), codes.end());
}

/**
 * @brief Processes a single code argument.
 *
 * @param codes The vector of codes to add to.
 * @param code The code to add.
 */
inline void processArg(std::vector<Code>& codes, Code code) {
  codes.push_back(code);
}

/**
 * @brief Processes a Modifier argument.
 *
 * @param codes The vector of codes to add to.
 * @param mod The Modifier to add.
 */
static void processArg(std::vector<Code>& codes, const Modifier& mod) {
  codes.insert(codes.end(), mod.getCodes().begin(), mod.getCodes().end());
}

/**
 * @brief Template class for colored exceptions.
 *
 * @tparam ExceptionType The type of the base exception.
 */
template<typename ExceptionType>
class ColoredException : public ExceptionType {
public:
  /**
   * @brief Constructs a ColoredException with a message and a color.
   *
   * @param message The exception message.
   * @param color The color of the exception message.
   */
  ColoredException(const std::string& message, const Modifier& color)
    : ExceptionType(message), color(color) {
    updateWhatStr();
  }

  /**
   * @brief Gets the exception message.
   *
   * @return const char* The exception message.
   */
  const char* what() const noexcept override {
    return what_str.c_str();
  }

private:
  Modifier color;
  std::string what_str;

  /**
   * @brief Updates the what string with the colored message.
   */
  void updateWhatStr() {
    std::ostringstream oss;
    oss << color << ExceptionType::what() << Modifier(Code::FG_DEFAULT_NO_BOLD);
    what_str = oss.str();
  }
};

/**
 * @brief Creates a colored runtime_error.
 *
 * @param message The exception message.
 * @param color The color of the exception message.
 * @return ColoredException<std::runtime_error> The colored runtime_error.
 */
inline ColoredException<std::runtime_error> makeRuntimeError(const std::string& message, const Modifier& color = Modifier{Code::FG_RED}) {
  return ColoredException<std::runtime_error>(message, color);
}

/**
 * @brief Creates a colored invalid_argument.
 *
 * @param message The exception message.
 * @param color The color of the exception message.
 * @return ColoredException<std::invalid_argument> The colored invalid_argument.
 */
inline ColoredException<std::invalid_argument> makeInvalidArgument(const std::string& message, const Modifier& color = Modifier{Code::FG_RED}) {
  return ColoredException<std::invalid_argument>(message, color);
}

/**
 * @brief Macro to create a runtime_error with optional color.
 */
#define GET_MACRO(_1, _2, NAME, ...) NAME
#define color_runtime_error(...) GET_MACRO(__VA_ARGS__, color_runtime_error_2, color_runtime_error_1)(__VA_ARGS__)
#define color_invalid_argument(...) GET_MACRO(__VA_ARGS__, color_invalid_argument_2, color_invalid_argument_1)(__VA_ARGS__)

/**
 * @brief Macro to create a runtime_error with a message.
 */
#define color_runtime_error_1(msg) makeRuntimeError(msg)

/**
 * @brief Macro to create a runtime_error with a message and a color.
 */
#define color_runtime_error_2(msg, color) makeRuntimeError(msg, color)

/**
 * @brief Macro to create an invalid_argument with a message.
 */
#define color_invalid_argument_1(msg) makeInvalidArgument(msg)

/**
 * @brief Macro to create an invalid_argument with a message and a color.
 */
#define color_invalid_argument_2(msg, color) makeInvalidArgument(msg, color)

} // namespace Color

/**
 * @example
 *
 * // Example usage of the Modifier class
 * #include <iostream>
 * #include "color.h"
 *
 * int main() {
 *   std::cout << Color::Modifier::red << "This text is red." << Color::Modifier::reset << std::endl;
 *   std::cout << Color::Modifier::bgBlue << "This text has a blue background." << Color::Modifier::reset << std::endl;
 *   std::cout << Color::Modifier::bold << "This text is bold." << Color::Modifier::reset << std::endl;
 *
 *   // Combining multiple modifiers
 *   std::cout << Color::colorCombiner(Color::Code::FG_GREEN, Color::Code::BOLD) << "This text is green and bold." << Color::Modifier::reset << std::endl;
 *
 *   // Using colored exceptions
 *   try {
 *     throw color_runtime_error("This is a colored runtime error.");
 *   } catch (const std::runtime_error& e) {
 *     std::cerr << e.what() << std::endl;
 *   }
 *
 *   return 0;
 * }
 */
