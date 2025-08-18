# Revelio — starter Makefile
# References:
# - GNU Make manual (automatic variables like $@, $<): https://www.gnu.org/s/make/manual/ (see docs)
# - GCC options: optimize (-O3) and tune for host (-march=native)
# - MSVC users: prefer the CLI example in README or use CMake/NMake separately.

CXX ?= g++
PREFIX ?= /usr/local
SRC := revelio.cpp
BIN := revelio

# Default flags for release builds
CXXFLAGS ?= -std=c++17 -O3 -march=native -Wall -Wextra -Wpedantic -DNDEBUG
LDFLAGS ?=

.PHONY: all release debug clean install uninstall format test help

all: release

release: $(BIN)

$(BIN): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) $(LDFLAGS) -o $@

# Address/UB sanitizers for quick runtime checks (GCC/Clang)
debug: CXXFLAGS := -std=c++17 -O0 -g -Wall -Wextra -Wpedantic -fsanitize=address,undefined -fno-omit-frame-pointer
debug: LDFLAGS := -fsanitize=address,undefined
debug: $(BIN)

test: $(BIN)
	./$(BIN) --help || true

install: $(BIN)
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 0755 $(BIN) $(DESTDIR)$(PREFIX)/bin/$(BIN)

uninstall:
	rm -f $(DESTDIR)$(PREFIX)/bin/$(BIN)

# Optional source formatting (requires clang-format)
format:
	@command -v clang-format >/dev/null 2>&1 && clang-format -i $(SRC) || echo "clang-format not found; skipping"

clean:
	rm -f $(BIN)

help:
	@echo "Targets:"
	@echo "  release (default) - build with -O3"
	@echo "  debug              - build with sanitizers"
	@echo "  test               - run a lightweight self-check"
	@echo "  install            - install to \`$$(PREFIX)/bin\` (default PREFIX=$(PREFIX))"
	@echo "  uninstall          - remove the installed binary"
	@echo "  format             - clang-format the source (if available)"
	@echo "  clean              - remove build artifacts"
