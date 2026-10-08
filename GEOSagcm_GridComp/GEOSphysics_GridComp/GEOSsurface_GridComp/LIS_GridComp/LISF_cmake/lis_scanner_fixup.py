#!/usr/bin/env python3
"""Recover Fortran module dependencies that CMake's scanner drops.

cmFortranParser permanently loses sync on a line-continuation ``&`` that is
immediately followed by a preprocessor directive.  Every MODULE defined after
that point in the file is missing from the target's ``provides`` list, so no
ordering edge reaches the consumers and parallel make compiles them first.

This script over-approximates: any source containing that pattern is treated
as suspect, and *every* module it defines gets an explicit edge to *every*
source that USEs it.  The extra edges are harmless where CMake already got it
right, and the scan is re-run at configure time so the fixups stay correct
across LISF bumps.

Usage: lis_scanner_fixup.py <sources-list-file> <output-file>

Output is one record per consumer:  ``<consumer>|<provider>,<provider>,...``
"""

import re
import sys

CONTINUATION = re.compile(r"&\s*(!.*)?$")
DIRECTIVE = re.compile(r"^\s*#")
MODULE_DEF = re.compile(r"^\s*module\s+([A-Za-z_]\w*)\s*(!.*)?$", re.IGNORECASE)
USE_STMT = re.compile(
    r"^\s*use\s*(?:,\s*intrinsic\s*)?(?:::)?\s*([A-Za-z_]\w*)", re.IGNORECASE
)
FORTRAN_EXT = (".F90", ".f90", ".F", ".f")


def read_lines(path):
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        return handle.read().splitlines()


def scanner_desyncs(lines):
    """True if a continuation line is followed by a preprocessor directive."""
    for index, line in enumerate(lines):
        if not CONTINUATION.search(line):
            continue
        for following in lines[index + 1:]:
            if not following.strip():
                continue
            if DIRECTIVE.match(following):
                return True
            break
    return False


def main():
    sources_file, output_file = sys.argv[1], sys.argv[2]
    sources = [s for s in read_lines(sources_file) if s.strip()]

    defines = {}  # source -> {MODULE, ...}, suspect sources only
    uses = {}  # source -> {MODULE, ...}

    for source in sources:
        if not source.endswith(FORTRAN_EXT):
            continue
        lines = read_lines(source)

        used = {m.group(1).upper() for m in map(USE_STMT.match, lines) if m}
        if used:
            uses[source] = used

        if not scanner_desyncs(lines):
            continue
        defined = {m.group(1).upper() for m in map(MODULE_DEF.match, lines) if m}
        defined.discard("PROCEDURE")
        if defined:
            defines[source] = defined

    records = []
    for consumer, used in sorted(uses.items()):
        providers = sorted(
            provider
            for provider, modules in defines.items()
            if provider != consumer and used & modules
        )
        if providers:
            records.append("{}|{}".format(consumer, ",".join(providers)))

    with open(output_file, "w", encoding="utf-8") as handle:
        handle.write("\n".join(records))
        if records:
            handle.write("\n")


if __name__ == "__main__":
    main()
