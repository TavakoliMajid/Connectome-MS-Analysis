# save as: inspect_all_connectomes.py
import os, re, pathlib, re

BASE = r"C:\\Majid\\image_processing\\INsIDER_Subj\\INsIDER_Subj\\Processed_Connectomes"  # <-- set your path

FILE_RE = re.compile(r"^(INsIDER_[^_]+)_(ACT|TREKKER)\.csv$", re.IGNORECASE)

def tokenize(line: str):
    s = line.strip().strip(",;")
    if not s:
        return []
    s = s.replace("\t", ",").replace(";", ",")
    s = s.replace("[", "").replace("]", "").replace('"', '').replace("'", "")
    parts = [t for t in s.split(",") if t != ""]
    if len(parts) == 1:
        parts = [t for t in re.split(r"\s+", s) if t != ""]
    return parts

def analyze_file(path: pathlib.Path):
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    lengths = [len(tokenize(ln)) for ln in lines if ln.strip()]
    if not lengths:
        return {"empty": True}
    unique = sorted(set(lengths))
    return {
        "n_lines": len(lines),
        "unique_lengths": unique,
        "min_len": min(lengths),
        "max_len": max(lengths),
        "is_square": (len(unique) == 1 and unique[0] == len(lines))
    }

def main():
    folder = pathlib.Path(BASE)
    files = [f for f in folder.iterdir() if f.is_file() and FILE_RE.match(f.name)]
    print(f"Scanning {len(files)} CSV files in {folder}")
    problems = []
    for f in files:
        info = analyze_file(f)
        if "empty" in info:
            print(f"[EMPTY] {f.name}")
            problems.append((f.name, "empty"))
            continue
        if not info["is_square"]:
            print(f"[PROBLEM] {f.name}: unique row lengths={info['unique_lengths']} "
                  f"(lines={info['n_lines']})")
            problems.append((f.name, f"row_lengths={info['unique_lengths']}"))
        else:
            print(f"[OK] {f.name}: {info['n_lines']}x{info['unique_lengths'][0]}")
    print("="*60)
    print(f"Total files with problems: {len(problems)}")
    for fn, why in problems[:10]:
        print(f"  - {fn}: {why}")
    if len(problems) > 10:
        print("  ... (more)")
        
if __name__ == "__main__":
    main()
