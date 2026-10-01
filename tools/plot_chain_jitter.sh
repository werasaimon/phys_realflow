#!/bin/sh
# Rebuild the recorded chain figures and normalize generated SVG whitespace for Git.
# Run from the repository root; requires Gnuplot and Graphviz, no simulation or Python.
set -eu
gnuplot tools/plot_chain_jitter.gnuplot
dot -Tsvg docs/images/chain-jitter-20261001/split-energy.dot \
    -o docs/images/chain-jitter-20261001/split-energy.svg
for svg_file in docs/images/chain-jitter-20261001/*.svg; do
    awk '
        { sub(/[ \t]+$/, "") }
        /^$/ { pending++; next }
        { for (i = 0; i < pending; i++) print ""; pending = 0; print }
    ' "$svg_file" > "$svg_file.tmp"
    mv -- "$svg_file.tmp" "$svg_file"
done
