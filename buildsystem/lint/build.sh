#!/usr/bin/env bash
fail=0
for dir in src tests include; do
  find ${CMAKE_SOURCE_DIR}/$dir -name '*.c' -o -name '*.cpp' -o -name '*.h' -o -name '*.hpp' \
    | while read f
      do
        ${CLANGFORMAT_EXECUTABLE} -style=file --dry-run -Werror $f
        if [ $? -ne 0 ]; then
          echo Found incorrectly formatted file: $f
          fail=$((fail + 1))
        fi
      done
done
echo Got $fail failures.
exit $fail
