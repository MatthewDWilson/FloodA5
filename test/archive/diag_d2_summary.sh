#!/usr/bin/env bash
# test/diag_d2_summary.sh — damping diagnostics (D2, 2026-10): summarise a set of
# runs named <cfg>_c07, <cfg>_twin (rainfall perturbed by 1e-10), <cfg>_c05 and
# <cfg>_c035 (Courant 0.5 / 0.35) in DIR: twin divergence at the last frame and
# its growth, stored volume per Courant number, and the rms depth difference
# from the Courant-0.7 run. Run from the repository root:
#   bash test/diag_d2_summary.sh DIR CFG...
D=$1; shift; T=test
for c in "$@"; do
  [ -f "$D/${c}_c035.log" ] && grep -q "Simulation finished" "$D/${c}_c035.log" || { echo "$c: incomplete"; continue; }
  tw=$(julia --project=. "$T/diag_divergence.jl" "$D/${c}_c07.h5" "$D/${c}_twin.h5" twin 2>/dev/null)
  last=$(echo "$tw" | grep -E "^ +[0-9]+ " | tail -1); gr=$(echo "$tw" | grep growth)
  vol() { grep "domain_vol=" "$D/$1.log" | tail -1 | grep -oE 'domain_vol=[0-9.e+]+' | cut -d= -f2; }
  c5=$(julia --project=. "$T/diag_divergence.jl" "$D/${c}_c07.h5" "$D/${c}_c05.h5" c05 2>/dev/null | grep -E "^ +[0-9]+ " | tail -1)
  c35=$(julia --project=. "$T/diag_divergence.jl" "$D/${c}_c07.h5" "$D/${c}_c035.h5" c035 2>/dev/null | grep -E "^ +[0-9]+ " | tail -1)
  printf "%-8s twin@6h: %s | %s\n" $c "$(echo $last | awk '{printf "rms %s max %s >1cm %s", $3, $4, $5}')" "$gr"
  printf "%-8s storage C0.7/0.5/0.35: %s / %s / %s Mm3 | vs C0.7 rms dh: C0.5 %s  C0.35 %s\n" "" \
    $(echo "$(vol ${c}_c07) $(vol ${c}_c05) $(vol ${c}_c035)" | awk '{printf "%.3f %.3f %.3f", $1/1e6, $2/1e6, $3/1e6}') \
    "$(echo $c5 | awk '{print $3}')" "$(echo $c35 | awk '{print $3}')"
done
