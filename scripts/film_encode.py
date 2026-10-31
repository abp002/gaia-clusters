"""The film's frames (scripts/film.py), with titles, into an mp4.

uv run scripts/film_encode.py <region> [--preview]    # figures/film-<region>.mp4 (-preview.mp4)

Titles follow film.py's shots, in seconds, and every number in them comes from the region's meta.json.
Fonts are macOS's Charter and SF Mono (New York's default instance has no hyphen for ffmpeg).
"""
import argparse
import json
import subprocess
from pathlib import Path

LY = 3.2616  # light years per parsec
SERIF, MONO = "/System/Library/Fonts/Supplemental/Charter.ttc", "/System/Library/Fonts/SFNSMono.ttf"
FADE = 0.6

args = argparse.ArgumentParser()
args.add_argument("region")
args.add_argument("--preview", action="store_true")
a = args.parse_args()
root = Path(__file__).resolve().parent.parent
meta = json.loads((root / "viewer" / "data" / a.region / "meta.json").read_text())
frames = root / "data" / "film" / (a.region + ("-preview" if a.preview else ""))
out = root / "figures" / f"film-{a.region}{'-preview' if a.preview else ''}.mp4"

recovered = [c for c in meta["clusters"] if c["recovered"]]
big = max((c for c in recovered if c["dist"] < 1000), key=lambda c: c["n"])
ly = lambda pc: f"{round(pc * LY, -3 if pc * LY > 10_000 else -2):,.0f}"
TITLES = [  # (from, to, title, line under it)
    (1.0, 5.2, f"{meta['stars'] / 1e6:.1f} million stars, measured by Gaia",
     f"a slice of the Milky Way towards {meta['title']}"),
    (8.0, 13.0, f"Out to {ly(meta['reach'])} light years from the Sun", "each star at its own distance"),
    (21.0, 27.0, "Stars born together move together", "three million years, at the motions Gaia measured"),
    (37.5, 44.5, f"{meta['recovered']} of {meta['recoverable']} known clusters recovered",
     f"in red, the {meta['recoverable'] - meta['recovered']} that HDBSCAN missed · scored against Hunt & Reffert (2023)"),
    (50.5, 54.8, big["name"].replace("_", " "), f"{ly(big['dist'])} light years away"),
    (56.0, 60.0, "gaia-clusters", "alejandroborrego.dev/gaia-clusters"),
]


def text(t, size, y, font, start, end, colour="white"):
    t = t.replace("\\", "\\\\").replace(":", "\\:").replace("'", "’").replace(",", "\\,")
    alpha = (f"if(lt(t\\,{start})\\,0\\,if(lt(t\\,{start + FADE})\\,(t-{start})/{FADE}\\,"
             f"if(lt(t\\,{end - FADE})\\,1\\,if(lt(t\\,{end})\\,({end}-t)/{FADE}\\,0))))")
    return (f"drawtext=fontfile={font}:text='{t}':fontsize=h*{size}:fontcolor={colour}:alpha='{alpha}'"
            f":x=w*0.07:y=h*{y}")


filters = []
for start, end, title, under in TITLES:
    filters += [text(title, 0.052, 0.78, SERIF, start, end), text(under, 0.021, 0.86, MONO, start + 0.3, end, "0xb8c4dd")]
fade = "fade=t=out:st=54.5:d=1.2"  # to black under the last title
rate = 10 if a.preview else 30
subprocess.run(["ffmpeg", "-v", "error", "-y", "-framerate", str(rate), "-pattern_type", "glob", "-i", str(frames / "*.png"),
                "-vf", ",".join([fade, *filters, "format=yuv420p"]), "-c:v", "libx264", "-crf", "20", "-preset", "slow",
                "-movflags", "+faststart", str(out)], check=True)
print(f"wrote {out}")
