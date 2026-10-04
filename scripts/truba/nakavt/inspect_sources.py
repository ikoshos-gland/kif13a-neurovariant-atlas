"""C ekseni: indirilen ek dosyalari duz metne/TSV'ye cevirir ve hedef genleri arar.

Compute node'da calisir (SLURM). Excel sayfalarini interim/sheets/<dosya>/<sayfa>.tsv olarak yazar,
PDF ve DOCX'i interim/text/ altina duz metin olarak cikarir, Chundru zip'ini acar, sonra her sayfanin
basligini ve KIF13A / KIF1A / KIF13B satirlarini log'a basar.
"""
from __future__ import annotations

import csv
import gzip
import re
import sys
import zipfile
from pathlib import Path

BASE = Path("/arf/scratch/avural/kif13a_nakavt")
sys.path.insert(0, str(BASE / "pylib"))

import openpyxl  # noqa: E402

RAW = BASE / "raw"
SHEETS = BASE / "interim" / "sheets"
TEXT = BASE / "interim" / "text"
GENES = re.compile(r"(?<![A-Z0-9])(KIF13A|KIF1A|KIF13B)(?![A-Z0-9])")


def xlsx_to_tsv(path: Path) -> list[Path]:
    out_dir = SHEETS / path.stem
    out_dir.mkdir(parents=True, exist_ok=True)
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    outs = []
    for ws in wb.worksheets:
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", ws.title)[:80]
        out = out_dir / f"{safe}.tsv"
        n = 0
        with out.open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            for row in ws.iter_rows(values_only=True):
                w.writerow(["" if v is None else str(v).replace("\t", " ").replace("\n", " ") for v in row])
                n += 1
        print(f"  sheet {path.name} :: {ws.title!r} -> {out.name} ({n} rows)", flush=True)
        outs.append(out)
    return outs


def docx_text(path: Path) -> Path:
    with zipfile.ZipFile(path) as z:
        xml = z.read("word/document.xml").decode("utf-8", "replace")
    xml = re.sub(r"</w:p>", "\n", xml)
    xml = re.sub(r"</w:tc>", "\t", xml)
    txt = re.sub(r"<[^>]+>", "", xml)
    out = TEXT / f"{path.stem}.txt"
    out.write_text(txt, encoding="utf-8")
    print(f"  docx -> {out.name} ({len(txt):,} chars)", flush=True)
    return out


def pdf_text(path: Path) -> Path | None:
    try:
        from pypdf import PdfReader
    except ImportError:
        print("  pypdf yok, PDF atlandi", flush=True)
        return None
    reader = PdfReader(str(path))
    parts = []
    for i, page in enumerate(reader.pages, start=1):
        parts.append(f"\n=== page {i} ===\n" + (page.extract_text() or ""))
    out = TEXT / f"{path.stem}.txt"
    out.write_text("".join(parts), encoding="utf-8")
    print(f"  pdf -> {out.name} ({len(reader.pages)} pages)", flush=True)
    return out


def show_sheet(tsv: Path, max_hits: int = 40) -> None:
    with tsv.open(encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    print(f"\n### {tsv.parent.name}/{tsv.name}: {len(rows)} rows", flush=True)
    for r in rows[:4]:
        print("   HEAD", " | ".join(c[:60] for c in r[:30]))
    hits = [r for r in rows if any(GENES.search(c) for c in r)]
    print(f"   target-gene rows: {len(hits)}")
    for r in hits[:max_hits]:
        print("   HIT ", " | ".join(c[:80] for c in r[:40]))


def main() -> None:
    SHEETS.mkdir(parents=True, exist_ok=True)
    TEXT.mkdir(parents=True, exist_ok=True)
    tsvs = []
    for p in sorted(RAW.glob("*.xlsx")):
        tsvs += xlsx_to_tsv(p)
    for p in sorted(RAW.glob("*.docx")):
        docx_text(p)
    for p in sorted(RAW.glob("*.pdf")):
        pdf_text(p)
    for p in sorted(RAW.glob("chundru_supp_data*.zip")):
        dest = BASE / "interim" / "chundru_data"
        dest.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(p) as z:
            for info in z.infolist():
                print(f"  zip member {info.filename} {info.file_size:,} B", flush=True)
            z.extractall(dest)
        for f in sorted(dest.rglob("*")):
            if "__MACOSX" in f.parts:  # zip icindeki macOS kaynak catali dosyalari Excel degildir
                continue
            if f.is_file() and f.suffix.lower() in {".xlsx"}:
                tsvs += xlsx_to_tsv(f)
            elif f.is_file() and f.suffix.lower() in {".tsv", ".txt", ".csv", ".gz"}:
                opener = gzip.open if f.suffix == ".gz" else open
                with opener(f, "rt", encoding="utf-8", errors="replace") as fh:
                    head = [next(fh, "") for _ in range(3)]
                    fh.seek(0)
                    hits = [line for line in fh if GENES.search(line)]
                print(f"\n### {f.relative_to(dest)}: head", flush=True)
                for h in head:
                    print("   HEAD", h.rstrip()[:400])
                print(f"   target-gene lines: {len(hits)}")
                for h in hits[:60]:
                    print("   HIT ", h.rstrip()[:600])
    for t in tsvs:
        show_sheet(t)
    for p in sorted(RAW.glob("ddg2p*.csv.gz")):
        with gzip.open(p, "rt", encoding="utf-8") as fh:
            rd = csv.reader(fh)
            head = next(rd)
            print(f"\n### {p.name} header: {head}")
            for r in rd:
                if r and GENES.search(",".join(r[:3])):
                    print("   HIT ", " | ".join(c[:120] for c in r))


if __name__ == "__main__":
    main()
