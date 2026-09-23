#!/usr/bin/env bash
# Download dbSNP compressed VCF + index 
#
# Use: bash ./download_dbSNP.sh [GRCh38|GRCh37] [output_folder]
#
# Warning: heavy file. If interrupted, download will resume by rerruning the command
set -euo pipefail

BUILD="${1:-GRCh38}"
DEST="${2:-reference}"

case "${BUILD}" in 
  GRCh38) FILE="GCF_000001405.40.gz" ;;   
  GRCh37) FILE="GCF_000001405.25.gz" ;;
  *) echo "Build must be GRCh38 or GRCh37" >&2; exit 2 ;;
esac

URL="https://ftp.ncbi.nih.gov/snp/latest_release/VCF"
mkdir -p "${DEST}"

fetch() {
  if command -v wget >/dev/null; then wget -c -O "$2" "$1"
  else curl -fL -C - -o "$2" "$1"; fi
}

for ext in ".md5" ".tbi" ""; do
  echo ">> ${FILE}${ext}"
  fetch "${URL}/${FILE}${ext}" "${DEST}/${FILE}${ext}
done

echo ">> conferindo md5 (demora alguns minutos)"
expected=$(awk '{print $1}' "$DEST/$FILE.md5")
if command -v md5sum >/dev/null; then actual=$(md5sum "$DEST/$FILE" | awk '{print $1}')
else actual=$(md5 -q "$DEST/$FILE"); fi
if [ "$expected" != "$actual" ]; then
  echo "md5 não confere! Apague $DEST/$FILE e baixe de novo." >&2; exit 1
fi

touch "${DEST}/${FILE}.tbi"

echo "Download complete"
