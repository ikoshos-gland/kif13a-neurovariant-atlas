r"""Ders 1b: LOEUF ve missense Z neden var, ham oran neden yetmiyor.

Calistirmak icin:
    python dersler\ders_01b_loeuf_ve_z.py

Fikir: ayni sayilari uc farkli sekilde okuyup farki gormek.
  1. ham oran         = gozlenen / beklenen
  2. LOEUF            = o oranin guven araliginin UST siniri
  3. naif Z           = (beklenen - gozlenen) / sqrt(beklenen)   <- kendimiz hesapliyoruz
  4. gnomAD'nin Z'si  = ayni fikir ama ampirik olarak kalibre edilmis
"""

import math
import requests

API = "https://gnomad.broadinstitute.org/api"

SORU = """
{
  gene(gene_symbol: "%s", reference_genome: GRCh38) {
    symbol
    gnomad_constraint {
      obs_lof exp_lof oe_lof oe_lof_upper
      obs_mis exp_mis oe_mis mis_z
    }
  }
}
"""

# Kucuk genler bilerek listede: LOEUF'un neden var oldugunu onlar gosteriyor.
GENLER = ["INS", "UBL5", "DEFB1", "KIF13A", "KIF1A", "TCF4", "ANKRD11", "SCN1A"]


def cek(sembol):
    c = requests.post(API, json={"query": SORU % sembol}, timeout=60)
    c.raise_for_status()
    v = c.json()
    if "errors" in v:
        raise RuntimeError(v["errors"])
    g = v["data"]["gene"]
    return g["symbol"], g["gnomad_constraint"]


def naif_z(gozlenen, beklenen):
    """Sayma verisinde gurultunun kabaca sqrt(beklenen) oldugu varsayimiyla."""
    if not beklenen:
        return float("nan")
    return (beklenen - gozlenen) / math.sqrt(beklenen)


if __name__ == "__main__":
    print(f"{'gen':9s} {'goz LoF':>8s} {'bek LoF':>8s} {'o/e':>6s} {'LOEUF':>6s} "
          f"| {'goz mis':>8s} {'bek mis':>8s} {'o/e':>6s} {'naif Z':>7s} {'gnomAD Z':>9s}")
    print("-" * 95)
    for sembol in GENLER:
        ad, k = cek(sembol)
        if k is None:
            print(f"{ad:9s} kisitlilik verisi yok")
            continue
        print(f"{ad:9s} {k['obs_lof']:8.0f} {k['exp_lof']:8.1f} {k['oe_lof']:6.3f} {k['oe_lof_upper']:6.3f} "
              f"| {k['obs_mis']:8.0f} {k['exp_mis']:8.1f} {k['oe_mis']:6.3f} "
              f"{naif_z(k['obs_mis'], k['exp_mis']):7.2f} {k['mis_z']:9.2f}")

    print()
    print("Bak: kucuk genlerde ham o/e ile LOEUF arasindaki makas aciliyor.")
    print("Bak: naif Z her zaman gnomAD'nin Z'sinden buyuk. Neden oldugunu derste konusuyoruz.")
