# KIF13A NeuroVariant Atlas

İnsan KIF13A genindeki missense varyantların nörolojik etkisini, paralog KIF1A'nın klinik bilgisini
(KAND) ve 2026 yapı/genom modellerini birleştirerek önceliklendiren, yeniden üretilebilir analiz hattı.

**Durum (15 Eyl 2026):** veri erişim katmanı kuruldu ve doğrulandı; WP1 çekirdek veri dondurması alındı.

## Hızlı başlangıç

```powershell
# sağlık kontrolü: 14 MCP sunucusu, her biri gerçek bir KIF13A sorgusuyla
mcp\atlas_server\.venv\Scripts\python.exe mcp\smoke_test.py

# WP1 ham verisini yeniden çek (gnomAD, ClinVar, AlphaMissense, AlphaFold, UniProt, PDB, ifade)
mcp\atlas_server\.venv\Scripts\python.exe scripts\pull_wp1_core.py
```

Sunucular kullanıcı kapsamında kayıtlı; herhangi bir klasörde açılan yeni Claude Code oturumunda yüklenir. Kurulum, kaynak
eşlemesi ve bilinen tuzaklar: [docs/MCP_KURULUM.md](docs/MCP_KURULUM.md).

## Referans kimlikler

Tüm protein pozisyonları MANE Select numaralandırmasıdır; kimlikler [config/targets.yaml](config/targets.yaml).

| Gen | MANE transkript | UniProt (MANE izoformu) | Uzunluk |
|---|---|---|---|
| KIF13A | ENST00000259711.11 / NM_022113.6 | Q9H1H9 | 1805 aa |
| KIF1A | ENST00000498729.9 / NM_001244008.2 | Q12756-3 (kanonik Q12756 = 1690 aa) | 1791 aa |
| KIF13B | ENST00000524189.6 / NM_015254.4 | Q9NQT8 | 1826 aa |
| KIF1B | ENST00000676179.1 / NM_001365951.3 | O60333 | 1816 aa |
| KIF1C | ENST00000320785.10 / NM_006612.6 | O43896 | 1103 aa |

## İlk proje planına düzeltmeler

- KIF13A referans proteini 1770 değil **1805 aa**; doygunluk matrisi 1805 × 19 = **34.295** değişim.
- KIF1A klinik numaralandırması (1791 aa) UniProt kanoniğinden (1690 aa) 394. kalıntıdan sonra kayar.
- **8H3T bir KIF1A yapısı değil.** KIF1A motor yapıları: 8UTN-8UTS, P305L 8UTT-8UTW, P364L 8UTY,
  R350G 9YA5/9YA7, R350W 9YAB/9YAI.
- KIF1A.org'un herkese açık varyant tablosu yok; DECIPHER API anahtarı gerektiriyor.

## Klasörler

`config/` doğrulanmış kimlikler · `mcp/` sunucular, kurulum, test · `scripts/` veri çekme ·
`data/raw` ham veri + `PROVENANCE.jsonl` · `data/interim` türetilmiş haritalar · `workflow/` boru hattı ·
`results/`, `figures/` çıktılar · `docs/` belgeler

Lisans: MIT (repo açıldığında eklenecek).
