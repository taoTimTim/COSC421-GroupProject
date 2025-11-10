import os, shutil

# dataset root (inside the data folder)
root = "./data/ds003969"

# destination for .bdf files
dest = "./bdf_data"
os.makedirs(dest, exist_ok=True)

copied = 0

for dirpath, _, filenames in os.walk(root):
    for f in filenames:
        if f.endswith(".bdf"):
            src = os.path.join(dirpath, f)
            dst = os.path.join(dest, f)
            print(f"Copying: {src} -> {dst}")
            shutil.copy2(src, dst)
            copied += 1

print(f"\n✅ Done! Copied {copied} .bdf files into bdf_data/")
