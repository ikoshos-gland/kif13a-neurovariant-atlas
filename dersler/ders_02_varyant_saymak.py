r"""Ders 2: Bir genin butun varyantlarini cekmek, sinifina gore saymak, ve bir tuzaga dusmek.

Calistirmak icin:
    python dersler\ders_02_varyant_saymak.py

Ogrenecegimiz uc sey:
  1. gnomAD'ye "bu transkriptteki butun varyantlari ver" demek
  2. consequence alanina gore saymak (2.671 missense buradan geliyor)
  3. AYNI geni iki farkli sekilde sorunca neden iki farkli cevap geldigi
"""

from collections import Counter

import requests

API = "https://gnomad.broadinstitute.org/api"
MANE = "ENST00000259711"   # KIF13A'nin MANE Select transkripti

# ---------------------------------------------------------------- 1. transkript sorgusu
TRANSKRIPT_SORGUSU = """
{
  transcript(transcript_id: "%s", reference_genome: GRCh38) {
    transcript_id
    variants(dataset: gnomad_r4) {
      variant_id
      consequence
      hgvsp
      joint { ac an filters }
    }
  }
}
"""

# ---------------------------------------------------------------- 2. gen sorgusu (TUZAK)
GEN_SORGUSU = """
{
  gene(gene_symbol: "%s", reference_genome: GRCh38) {
    variants(dataset: gnomad_r4) {
      variant_id
      consequence
      transcript_id
    }
  }
}
"""


def sor(sorgu):
    c = requests.post(API, json={"query": sorgu}, timeout=180)
    c.raise_for_status()
    v = c.json()
    if "errors" in v:
        raise RuntimeError(v["errors"])
    return v["data"]


def guvenilir(varyant):
    """Bu varyantin gercekten tasiyicisi var mi?

    AC = allele count, yani bu varyanti tasiyan kromozom sayisi.
    AC0 filtresi = 'tasiyici gorundu ama hepsi kalite kontrolunden dustu'.
    Ikisi birlikte 'gercekten gorulmus' tanimini veriyor.
    """
    j = varyant.get("joint") or {}
    ac = j.get("ac") or 0
    filtreler = j.get("filters") or []
    return ac >= 1 and "AC0" not in filtreler


if __name__ == "__main__":
    print("1) MANE transkriptindeki butun varyantlar cekiliyor...")
    tx = sor(TRANSKRIPT_SORGUSU % MANE)["transcript"]
    varyantlar = tx["variants"]
    print(f"   {len(varyantlar):,} varyant geldi\n")

    print("2) consequence alanina gore sayim (en sik 10 sinif):")
    sayim = Counter(v["consequence"] for v in varyantlar)
    for sinif, n in sayim.most_common(10):
        print(f"   {sinif:32s} {n:6,}")

    missense = [v for v in varyantlar if v["consequence"] == "missense_variant"]
    print(f"\n3) missense sayisi: {len(missense):,}")

    saglam = [v for v in missense if guvenilir(v)]
    ac0 = [v for v in missense if not guvenilir(v)]
    print(f"   bunlardan gercekten tasiyicisi olan : {len(saglam):,}")
    print(f"   AC0 / sifir tasiyici                : {len(ac0):,}")

    print("\n4) TUZAK: ayni geni bir de 'gen sorgusu' ile soralim")
    gen = sor(GEN_SORGUSU % "KIF13A")["gene"]
    print(f"   gen sorgusu     : {len(gen['variants']):,} varyant")
    print(f"   transkript sorgusu: {len(varyantlar):,} varyant")
    transkriptler = Counter(v.get("transcript_id") for v in gen["variants"]
                            if v["consequence"] == "missense_variant")
    print(f"\n   gen sorgusundaki missense'ler {len(transkriptler)} farkli transkripte gore numaralanmis:")
    for t, n in transkriptler.most_common():
        isaret = "  <- MANE" if t and t.startswith(MANE) else ""
        print(f"     {t} {n:6,}{isaret}")
    print("\n   Iste bu yuzden proje kurali: HER ZAMAN transkript sorgusu.")
