"""Crop a PNG sequence with a common bounding box.

The crop is computed from the union of foreground pixels across all frames, so
the resulting images keep a stable camera frame when assembled as a video.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from PIL import Image


def foreground_mask(image: Image.Image) -> np.ndarray:
    """Return a mask for hand-drawn foreground strokes.

    The source frames have a white background and a light gray page border. A
    simple "not white" mask would include the border, so this mask keeps dark
    ink and saturated colored strokes while ignoring pale background elements.
    """

    rgba = np.asarray(image.convert("RGBA"))
    rgb = rgba[..., :3].astype(np.int16)
    alpha = rgba[..., 3]

    darkness = np.min(rgb, axis=-1) < 135
    saturation = np.max(rgb, axis=-1) - np.min(rgb, axis=-1) > 80
    return (alpha > 0) & (darkness | saturation)


def bbox_from_mask(mask: np.ndarray) -> tuple[int, int, int, int]:
    ys, xs = np.nonzero(mask)
    if len(xs) == 0 or len(ys) == 0:
        raise ValueError("no foreground pixels found")
    return int(xs.min()), int(ys.min()), int(xs.max()) + 1, int(ys.max()) + 1


def pad_bbox(
    bbox: tuple[int, int, int, int],
    size: tuple[int, int],
    padding: int,
) -> tuple[int, int, int, int]:
    left, top, right, bottom = bbox
    width, height = size
    return (
        max(0, left - padding),
        max(0, top - padding),
        min(width, right + padding),
        min(height, bottom + padding),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("images"))
    parser.add_argument("--output-dir", type=Path, default=Path("images/cropped"))
    parser.add_argument("--pattern", default="bg-*.png")
    parser.add_argument("--padding", type=int, default=80)
    args = parser.parse_args()

    paths = sorted(args.input_dir.glob(args.pattern))
    if not paths:
        raise SystemExit(f"No images found for {args.input_dir / args.pattern}")

    images = [Image.open(path).convert("RGBA") for path in paths]
    sizes = {image.size for image in images}
    if len(sizes) != 1:
        raise SystemExit(f"All images must have the same size, got: {sorted(sizes)}")

    boxes = [bbox_from_mask(foreground_mask(image)) for image in images]
    left = min(box[0] for box in boxes)
    top = min(box[1] for box in boxes)
    right = max(box[2] for box in boxes)
    bottom = max(box[3] for box in boxes)
    crop_box = pad_bbox((left, top, right, bottom), images[0].size, args.padding)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for path, image in zip(paths, images, strict=True):
        image.crop(crop_box).save(args.output_dir / path.name)

    print(f"Frames: {len(paths)}")
    print(f"Input size: {images[0].size[0]} x {images[0].size[1]}")
    print(
        "Crop box: "
        f"left={crop_box[0]}, top={crop_box[1]}, "
        f"right={crop_box[2]}, bottom={crop_box[3]}"
    )
    print(f"Output size: {crop_box[2] - crop_box[0]} x {crop_box[3] - crop_box[1]}")
    print(f"Output dir: {args.output_dir}")


if __name__ == "__main__":
    main()
