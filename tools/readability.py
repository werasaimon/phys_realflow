#!/usr/bin/env python3
"""Readability check of the SDK sources: no function longer than 60 lines, every file opens with a
comment paragraph, and the comment density per file. Prints a table and exits with 1 if a function
is too long or a file has no head paragraph.

A "function" is a definition line at column 0 - `Type Class::method(...) {` or `Type name(...) {` -
up to the next `}` at column 0. Usage: python tools/readability.py [root] [--limit 60]
"""
import os
import re
import sys

ROOT = next((a for a in sys.argv[1:] if not a.startswith('--')), '.')
LIMIT = 60
if '--limit' in sys.argv:
    LIMIT = int(sys.argv[sys.argv.index('--limit') + 1])
DIRS = ('src', 'samples', 'verification')
DEFINITION = re.compile(r'^[A-Za-z_][A-Za-z0-9_:<>,\s\*&]*\s\*?&?[A-Za-z_][A-Za-z0-9_:]*\s*\([^;]*\)\s*(const)?\s*(noexcept)?\s*\{\s*$')


def functions(lines):
    """Yields (name, length) for every column-0 function definition."""
    i = 0
    while i < len(lines):
        line = lines[i]
        if DEFINITION.match(line) and not line.startswith(('if', 'for', 'while', 'switch', 'else')):
            j = i + 1
            while j < len(lines) and lines[j] != '}':
                j += 1
            name = line.split('(')[0].split()[-1]
            yield name, j - i + 1
            i = j
        i += 1


def head_paragraph(lines):
    """True if the file starts (after #pragma once / blank lines) with a comment."""
    for line in lines[:4]:
        s = line.strip()
        if s.startswith('//') or s.startswith('/*'):
            return True
        if s and s != '#pragma once':
            return False
    return False


rows, long_functions, headless = [], [], []
for d in DIRS:
    for dirpath, _, names in os.walk(os.path.join(ROOT, d)):
        for name in sorted(names):
            if not name.endswith(('.cpp', '.h')):
                continue
            path = os.path.join(dirpath, name)
            with open(path, encoding='utf-8', errors='replace') as f:
                lines = f.read().split('\n')
            comments = sum(1 for l in lines if l.strip().startswith('//'))
            longest = 0
            for fname, length in functions(lines):
                longest = max(longest, length)
                if length > LIMIT:
                    long_functions.append((length, path, fname))
            if not head_paragraph(lines):
                headless.append(path)
            rows.append((path, len(lines), comments, longest))

print('%-52s %6s %8s %8s' % ('file', 'lines', 'comment', 'longest'))
for path, n, c, longest in rows:
    print('%-52s %6d %7d%% %8d' % (path.replace('\\', '/'), n, 100 * c // max(n, 1), longest))
print('\nfiles: %d, lines: %d, comment lines: %d' % (len(rows), sum(r[1] for r in rows), sum(r[2] for r in rows)))
if long_functions:
    print('\nFUNCTIONS LONGER THAN %d LINES:' % LIMIT)
    for length, path, fname in sorted(long_functions, reverse=True):
        print('%5d  %-45s %s' % (length, path.replace('\\', '/'), fname))
if headless:
    print('\nFILES WITHOUT A HEAD PARAGRAPH:')
    for path in headless:
        print('   ' + path.replace('\\', '/'))
if not long_functions and not headless:
    print('\nOK: every function is at most %d lines and every file explains itself.' % LIMIT)
sys.exit(1 if long_functions or headless else 0)
