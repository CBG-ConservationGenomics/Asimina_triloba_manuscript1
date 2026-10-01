#!/bin/bash
set -e #stops the script if any command fails

BUSCO_DIR="busco_results"
OUTFILE="busco_summary_table.tsv"

echo -e "Species\tC_pct\tS_pct\tD_pct\tF_pct\tM_pct\tn\tC_count\tS_count\tD_count\tF_count\tM_count\tTotal_count" > "$OUTFILE" #write headers

for dir in "$BUSCO_DIR"/*; do
  [ -d "$dir" ] || continue
  species=$(basename "$dir")

  summary=$(ls -t "$dir"/short_summary*.txt 2>/dev/null | head -n 1)
  if [ -z "$summary" ]; then
    echo "=== No BUSCO summary for $species (skipping) ==="
    continue
  fi

  awk -v sp="$species" '
    /C:[0-9.]+%/ && !got_pct {
      got_pct=1
      match($0, /C:([0-9.]+)%/, a); C=a[1]
      match($0, /S:([0-9.]+)%/, a); S=a[1]
      match($0, /D:([0-9.]+)%/, a); D=a[1]
      match($0, /F:([0-9.]+)%/, a); F=a[1]
      match($0, /M:([0-9.]+)%/, a); M=a[1]
      match($0, /n:([0-9]+)/,  a); n=a[1]
    }
    /Complete BUSCOs \(C\)/{Ccount=$1}
    /Complete and single-copy BUSCOs \(S\)/{Scount=$1}
    /Complete and duplicated BUSCOs \(D\)/{Dcount=$1}
    /Fragmented BUSCOs \(F\)/{Fcount=$1}
    /Missing BUSCOs \(M\)/{Mcount=$1}
    /Total BUSCO groups searched/{Tcount=$1}
    END {
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
        sp,C,S,D,F,M,n,Ccount,Scount,Dcount,Fcount,Mcount,Tcount
    }
  ' "$summary" >> "$OUTFILE"
done

echo "=== Wrote: $OUTFILE ==="
