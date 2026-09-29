#!/bin/sh
set -eu

if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    echo "usage: $0 input.xpm [output.jpg]" >&2
    exit 1
fi

in="$1"
out="${2:-${in%.*}.jpg}"

tmp="$(mktemp --suffix=.jpg)"
trap 'rm -f "$tmp"' EXIT

magick "$in" -crop 3600x1800+0+0 +repage "$tmp"

exiftool \
  -overwrite_original \
  -XMP-GPano:UsePanoramaViewer=True \
  -XMP-GPano:ProjectionType=equirectangular \
  -XMP-GPano:FullPanoWidthPixels=3600 \
  -XMP-GPano:FullPanoHeightPixels=1800 \
  -XMP-GPano:CroppedAreaImageWidthPixels=3600 \
  -XMP-GPano:CroppedAreaImageHeightPixels=1800 \
  -XMP-GPano:CroppedAreaLeftPixels=0 \
  -XMP-GPano:CroppedAreaTopPixels=0 \
  "$tmp" >/dev/null

mv "$tmp" "$out"
trap - EXIT

echo "wrote $out"
