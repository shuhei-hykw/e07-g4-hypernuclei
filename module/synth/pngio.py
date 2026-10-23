"""Minimal 8-bit grayscale PNG writer (no external dependency)."""

import struct
import zlib

import numpy as np

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
BIT_DEPTH = 8
COLOR_GRAY = 0
FILTER_NONE = 0
ZLIB_LEVEL = 6


def _chunk(tag, data):
  body = tag + data
  return (struct.pack(">I", len(data)) + body
          + struct.pack(">I", zlib.crc32(body) & 0xffffffff))


def write_gray_png(path, image):
  """Write a 2-D uint8 array as a grayscale PNG."""
  image = np.ascontiguousarray(image, dtype=np.uint8)
  height, width = image.shape
  header = struct.pack(">IIBBBBB", width, height, BIT_DEPTH, COLOR_GRAY,
                       0, 0, 0)
  rows = np.empty((height, width + 1), dtype=np.uint8)
  rows[:, 0] = FILTER_NONE
  rows[:, 1:] = image
  data = zlib.compress(rows.tobytes(), ZLIB_LEVEL)
  with open(path, "wb") as f:
    f.write(PNG_SIGNATURE + _chunk(b"IHDR", header) + _chunk(b"IDAT", data)
            + _chunk(b"IEND", b""))
