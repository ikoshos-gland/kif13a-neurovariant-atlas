r"""Ders 1c: Bir gende varyantlarin uc sinifi ve neden ayri ayri olculuyor.

Calistirmak icin:
    python dersler\ders_01c_uc_sinif.py

DNA uclu harf gruplari (kodon) halinde okunur, her kodon bir aminoasit soyler.
Tek harflik bir degisim uc sonuctan birini verir:

  sessiz (synonymous)  AAA -> AAG   ikisi de Lys. Protein AYNI. Secilim gormez.
  missense             AAA -> GAA   Lys yerine Glu. Protein var ama farkli.
  gen susturucu (LoF)  AAA -> TAA   "dur" emri. Protein yarim kalir, hucre mRNA'yi yok eder.

gnomAD gozlenen/beklenen hesabini bu uc sinif icin AYRI AYRI yapar, cunku
her sinif farkli bir soruya cevap verir. Sessiz sinif ayrica bir denetimdir:
secilim onu gormedigi icin orani 1'e yakin CIKMALIDIR. Cikmiyorsa model bozuktur.
"""

import requests

API = "https://gnomad.broadinstitute.org/api"

SORU = """
{
  gene(gene_symbol: "%s", reference_genome: GRCh38) {
    symbol
    gnomad_constraint {
      obs_syn exp_syn oe_syn
      obs_mis exp_mis oe_mis mis_z
      obs_lof exp_lof oe_lof oe_lof_upper
    }
  }
}
"""

GENLER = ["KIF13A", "KIF1A", "ANKRD11", "TCF4"]


def cek(sembol):
    c = requests.post(API, json={"query": SORU % sembol}, timeout=60)
    c.raise_for_status()
    v = c.json()
    if "errors" in v:
        raise RuntimeError(v["errors"])
    g = v["data"]["gene"]
    return g["symbol"], g["gnomad_constraint"]


if __name__ == "__main__":
    for sembol in GENLER:
        ad, k = cek(sembol)
        print(f"\n{ad}")
        print(f"  {'sinif':16s} {'gozlenen':>9s} {'beklenen':>9s} {'o/e':>7s}   okuma")
        print(f"  {'sessiz':16s} {k['obs_syn']:9.0f} {k['exp_syn']:9.1f} {k['oe_syn']:7.3f}   "
              f"denetim, 1'e yakin olmali")
        print(f"  {'missense':16s} {k['obs_mis']:9.0f} {k['exp_mis']:9.1f} {k['oe_mis']:7.3f}   "
              f"Z = {k['mis_z']:.2f}")
        print(f"  {'gen susturucu':16s} {k['obs_lof']:9.0f} {k['exp_lof']:9.1f} {k['oe_lof']:7.3f}   "
              f"LOEUF = {k['oe_lof_upper']:.3f}")

    print("\nDikkat: her gen IKI ayri hukum aliyor, biri missense icin biri gen susturucu icin.")
    print("Ayni gen bir tarafta cok kisitli, oteki tarafta hic kisitsiz olabilir.")
