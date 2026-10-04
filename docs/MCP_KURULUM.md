# MCP kurulumu ve veri erişimi

Son doğrulama: 15 Eyl 2026, Windows 11, Python 3.12.14 (uv), Node 24.17. Her sunucu gerçek bir KIF13A
sorgusuyla test edildi (`mcp/smoke_test.py`).

## 1. Mimari: üç katman

1. **Proje veri katmanı, `kif13a-atlas`** (bu repoda, `mcp/atlas_server`). Boru hattının dayandığı her
   veriyi çeker, `data/raw` altına yazar, her dosyaya kaynak URL, parametre, UTC zaman ve sha256 içeren
   bir `.prov.json` ekler. Figürler yalnızca buradan gelen dosyalardan üretilir.
2. **Keşif sunucuları** (üçüncü taraf, sürümleri sabit). Soru sormak, hızlı bakmak, literatür taramak için.
3. **MCP olmayan hesaplama** (TRUBA). Model çalıştırma işleri veri çekme değildir; aşağıda 5. bölüm.

## 2. Planın veri tablosu ve karşılığı

| Plan katmanı | Kaynak | Sunucu / araç | Durum |
|---|---|---|---|
| Dizi & varyant | UniProt Q9H1H9 | `uniprot`, `kif13a-atlas.uniprot_features` | çalışıyor |
| | Ensembl MANE | `ensembl`, `config/targets.yaml` | çalışıyor |
| | gnomAD v4.1 | `kif13a-atlas.gnomad_variants` / `gnomad_constraint`, `gnomad` | çalışıyor |
| | ClinVar | `kif13a-atlas.clinvar_variants`, `gnomad.gnomad_search_clinvar` | çalışıyor |
| | DECIPHER | yok | istemci anahtarı gerekiyor (403) |
| | KIF1A.org varyant tablosu | yok | herkese açık tablo/API yok; ClinVar + literatür ekleri kullanılacak |
| | OMIM #601255 | `tooluniverse` (OMIM_API_KEY ile) | anahtar bekliyor |
| Yapı | AlphaFold DB | `kif13a-atlas.alphafold_model`, `protein` | çalışıyor |
| | AlphaFold3 Server | yok | web arayüzü, API yok; elle iş gönderilir |
| | OpenFold3 / ESMFold2 | yok | TRUBA GPU işi |
| | PDB KIF1A cryo-EM | `kif13a-atlas.pdb_structures`, `pdbe` | çalışıyor (8H3T KIF1A değil, bkz. 6) |
| | Foldseek | `protein.protein_find_similar`, `tooluniverse` | çalışıyor |
| Varyant skorları | AlphaMissense | `kif13a-atlas.alphamissense`, `vep_annotate` | çalışıyor |
| | CADD v1.7 | `kif13a-atlas.cadd_region`, `vep_annotate` | çalışıyor |
| | SpliceAI, REVEL, ESM1b, PrimateAI, MetaRNN | `kif13a-atlas.vep_annotate` | çalışıyor |
| | PrimateAI-3D | yok | lisanslı; dbNSFP'deki PrimateAI (3D değil) var |
| | ESM1b/ESM-C tam doygunluk | yok | TRUBA işi; ESM Forge için ESM_API_KEY |
| | Evo 2 40B | `kif13a-atlas.evo2_score` (NVIDIA NIM forward) | çalışıyor; hedefli setler (çağrı başına en fazla 50 SNV), gen çapı doygunluk TRUBA |
| Ekspresyon & ağ | CELLxGENE | `kif13a-atlas.cellxgene_expression` (toplu) | çalışıyor; hücre düzeyi Census TRUBA'da |
| | Allen Brain Cell Atlas | yok | mevcut analiz: `ikoshos-gland/kif13a-allenbrain` |
| | HPA | `kif13a-atlas.hpa_summary`, `tooluniverse` | çalışıyor |
| | GTEx v10 | `kif13a-atlas.gtex_expression`, `tooluniverse` | çalışıyor |
| | STRING | `string` (resmi) | çalışıyor |
| | BioGRID | `tooluniverse` | BIOGRID_ACCESS_KEY bekliyor |
| Genetik epidemiyoloji | Open Targets | `opentargets` (resmi) | çalışıyor |
| | GWAS Catalog | `biomcp`, `tooluniverse` | çalışıyor |
| | PGC özet istatistikleri | yok | GB boyutunda dosyalar; WP3'te TRUBA'da indirilecek |
| Literatür | PubMed, bioRxiv | `pubmed`, `biorxiv` | çalışıyor |
| | Europe PMC, Semantic Scholar, OpenAlex | `tooluniverse` | OpenAlex anahtarsız: günlük 1.000 kredi, arama 10 kredi |
| | Google Scholar | `google-scholar` (kullanıcı kapsamı) | önceden kurulu |
| | Zotero | `zotero` | Zotero masaüstü açıkken |
| Ajan | Biomni | yok | Windows'ta gerçekçi değil (conda, 13-30 GB); gerekirse TRUBA/WSL |

## 3. Sıfırdan kurulum (başka makinede tekrar etmek için)

```powershell
python -m pip install --user uv
uv python install 3.12
cd mcp\atlas_server;  uv sync;  cd ..\..
cd mcp\node;  npm install;  cd ..\..
uv tool install --python 3.12 biomcp-cli==0.8.25
uv tool install --python 3.12 pdbe-mcp-server==1.1.6
uv tool install --python 3.12 zotero-mcp-server==0.12.4
uv tool install --python 3.12 tooluniverse==1.4.1
copy .env.example .env      # anahtarları doldur (hepsi isteğe bağlı)
mcp\atlas_server\.venv\Scripts\python.exe mcp\setup_mcp.py   # 14 sunucuyu kullanıcı kapsamında kaydeder
mcp\atlas_server\.venv\Scripts\python.exe mcp\smoke_test.py
```

`setup_mcp.py` sunucuları **kullanıcı kapsamında** kaydeder (`claude mcp add-json --scope user`, kayıt
`~/.claude.json` içinde). Hangi klasörde olursa olsun **yeni açılan her Claude Code oturumunda** yüklenirler;
çalışan bir oturum yeni sunucuları görmez. Yalnızca bu klasöre bağlamak istersen `--scope project`.
Kontrol: `claude mcp list`.

## 4. Anahtarlar

Tanımlı (15 Eyl 2026, ikisi de canlı doğrulandı): `NCBI_API_KEY` (E-utilities sınırı 3'ten 10 istek/s'ye
çıktı), `NVIDIA_API_KEY` (Evo 2 40B NIM). Yenisini eklemek için `.env` dosyasına yaz,
`mcp\atlas_server\.venv\Scripts\python.exe mcp\setup_mcp.py` çalıştır, yeni oturum aç. Microsoft Store
`python` komutunu kullanma: `%APPDATA%\uv` klasörünü göremediği için her sunucuyu "kurulu değil" sayar.
`.env` gitignore'da; anahtarlar repoya girmez, `~/.claude.json` içine gömülür.

## 5. Neden bazı şeyler MCP değil

MCP bir soru-cevap kanalıdır; araç çağrısı saniyeler, en fazla birkaç dakika sürmelidir. 34.295 missense
değişimi için ESM skorlama, AlphaFold3 kompleksleri, 100-200 ns MD veya milyonlarca hücrelik Census analizi
saatler ile günler sürer. Bunlar TRUBA'da SLURM işi olarak koşar, çıktıları `data/processed` altına gelir,
MCP ile yalnızca sonuç sorgulanır.

## 6. Kurulumda bulunan tuzaklar

- **8H3T KIF1A değil** (AlpH X-ışını yapısı). Doğru KIF1A cryo-EM yapıları: 8UTN-8UTS (yabanıl tip),
  8UTT-8UTW (P305L), 8UTY (P364L), 9YA5/9YA7 (R350G), 9YAB/9YAI (R350W). 8UTX aralığın içinde ama KIF1A değil.
- **gnomAD gen sorgusu MANE'ye göre numaralamıyor.** 2.823 missense satırı 5 transkripte dağılıyor; MANE
  transkript sorgusu 2.671 missense veriyor.
- **ClinVar alan adı değişti:** `germline_classification`. Eski sunucular (ör. openpharma clinvar) bu yüzden boş dönüyor.
- **koido/gnomad-mcp** v4 kısıtlılık alanlarını döndürmüyor; yerine cyanheads sunucusu kuruldu.
- **BioMCP Rust ile yeniden yazıldı**; PyPI paketi `biomcp-cli` (eski `biomcp-python` değil).
- **Windows Akıllı Uygulama Denetimi** `zotero-mcp.exe` başlatıcısını engelliyor (WinError 4551); sunucu
  `python -m zotero_mcp.cli serve` ile başlatılıyor. Güvenlik ayarına dokunulmadı.
- **Windows'ta stdio kilitlenmesi:** MCP sunucusu stdin'i okurken bir araç çağrısında numpy ilk kez içe
  aktarılırsa süreç kilitleniyor. `structure.py` Biopython'u modül yüklenirken içe aktarıyor.
- **mcp SDK 2.x** API'si farklı; özel sunucu `mcp<2` ile sabit.
- **uv**, Microsoft Store Python'unu yorumlayıcı olarak kullanamıyor; uv'nin kendi Python 3.12'si kurulu.
  "minor version link" uyarısı zararsız.
- **CELLxGENE Census** (`tiledbsoma`) Windows'ta yok; toplu ifade için Where's My Gene API kullanılıyor.
  `pc` alanı belgede 0-100 yazsa da gerçekte 0-1 oran; araç `n/total` kontrol sütunu ekliyor.
- **GTEx v10** GENCODE v39 kimliği ister (`ENSG00000137177.20`); v26 kimliğiyle boş döner.
- **OpenAlex anahtar istemiyor.** Yanıt başlıklarına göre anahtarsız günlük sınır 1.000 kredi, bir
  arama 10 kredi (yaklaşık 100 arama/gün). Araştırma ajanının "100 kredi" bilgisi yanlıştı.
- **Evo 2 NIM** logitleri `unembed` katmanından base64 npz olarak döndürüyor; sözlük bayt kodları
  (A=65, C=67, G=71, T=84). 8.192 nt pencere yaklaşık 5 s, yanıt 1,3 MB. Aynı girdi aynı skoru veriyor.
- **Evo 2 NIM tam nedensel davranmıyor.** Tek baz değişince varyanttan ÖNCEKİ bazların log-olasılığı da
  hafifçe kayıyor (4.096 bazda p99 0,03, en fazla 0,36, toplam yaklaşık 0,4 nat). `evo2_score` bu yüzden
  skoru varyant bazından itibaren topluyor. Referans panel (`scripts/evo2_control_panel.py`): klinik ve
  kanonik splice pozitifleri -17 ile -68, yaygın intronik varyantlar -2,5 ile +3,8; |skor| < 3 ayırt edilemez.
