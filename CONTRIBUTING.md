# Contributing Guidelines

This document defines the mandatory code formatting standards and submission requirements for this project. All contributions must conform to the specified formatting rules to be accepted.

## Code Formatting

This project uses **clang-format** to enforce consistent code style. All C++ code must be formatted according to our style configuration before submission.

### Formatting Style

Our code style is based on **Google C++ Style Guide** with the following custom modifications:

- **Brace Style**: Custom brace wrapping
  - Braces are placed on new lines for classes, functions, structs, unions, namespaces, enums
  - Braces are placed on new lines for control statements (if, for, while, etc.)
  - Braces are placed on new lines before catch and else statements
- **Indentation**: 2 spaces (no tabs)
- **Column Limit**: 80 characters per line

### Before Submitting a Pull Request

**1. Install clang-format**

We use version 13 of clang-format.

```bash
# Ubuntu/Debian
sudo apt install clang-format-13

# macOS
brew install clang-format@13

# Python (cross-platform)
#   - clang only 
pip install clang-format==13
#   - all dev requirements at once
pip install -r requirements-dev.txt

# Or check if already installed
clang-format --version
```

**2. Format your code**

Before committing any changes, run clang-format on all modified C++ files:

```bash
# Format a single file
clang-format -i path/to/your/file.cpp

# Format all C++ files recursively
find . \( -name "*.h" -o -name "*.hpp" -o -name "*.cc" -o -name "*.cpp" -o -name "*.cxx" \) -exec clang-format -i {} \;
```

**3. Verify formatting**

Check that your files are properly formatted:

```bash
# Check if files need formatting (should show no output if properly formatted)
find . \( -name "*.h" -o -name "*.hpp" -o -name "*.cc" -o -name "*.cpp" -o -name "*.cxx" \) -exec clang-format -i {} \;
git diff --exit-code
```

### Configuration Details

The project's code formatting is defined in [`.clang-format`](.clang-format):

Thank you for helping maintain code quality!
