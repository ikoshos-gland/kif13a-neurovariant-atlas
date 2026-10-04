# KIF13A hücre tipi eQTL kolokalizasyon ekseni (A ekseni)

**Durum: ÖN SONUÇ.** Analizi yapan ajan, boru hattı TRUBA'da bittikten sonra ama bu belgeyi yazmadan
durduruldu (4 Ekim 2026). Bu belge ana oturumda, ajanın bıraktığı tablolardan yazıldı; sayılar
`data/interim/eqtl_hucre_tipi/` altındaki dosyalardan okundu, yeniden hesaplanmadı. Script'ler
`scripts/truba/eqtl/` (SLURM dosyaları dahil), tam çalışma alanı TRUBA'da
`/arf/scratch/avural/kif13a_eqtl` (33 GB, ana iş günlüğü `logs/eqtl_pipe_6440189.out`).

## 0. Hüküm: oligodendrosit ve OPC için karşı kanıt

Sorulan soru (H_A): ak madde GWAS haplotipi (rs144963738 / rs115824472 / rs75428660, %3,2 frekans)
KIF13A ifadesini özellikle oligodendrositte ya da OPC'de değiştiriyor ve bu iki sinyal aynı nedensel
varyanttan geliyor mu?

**Cevap hayır.** KIF13A'nın oligodendrosit ve OPC'de güçlü eQTL'leri var, ama bunlar ak madde
sinyalinden **farklı varyantlar**. İki büyük tek çekirdek atlasında (Zeng 2026 ve SingleBrain 2025), dört
ak madde GWAS'ı ve kortikal kıvrım GWAS'ı için, oligodendrosit ve OPC'de PP.H3 (iki ayrı nedensel
varyant) 0,78 ile 1,00 arasında, PP.H4 (ortak varyant) en fazla 0,015. eQTL öncüleri GWAS öncüsüyle LD
içinde değil (r² en fazla 0,054), ve GWAS öncüsündeki eQTL p değeri oligodendrositte 0,78.

Tek ipucu mikroglia'da ve kırılgan: Zeng bağışıklık hücrelerinde üçüncü koşullu sinyal (rs73373199)
EUR LD paneliyle PP.H4 0,55-0,62 veriyor, tüm 1000G paneliyle 0,002-0,03'e düşüyor. Örnek içi LD ile
doğrulanmadan bulgu sayılmamalı.

## 1. Veri kaynakları

| Kaynak | Ne | Etkin N (hücre tipi başına, ortanca) | Erişim |
|---|---|---|---|
| Zeng ... Roussos, Nat Genet 2026, PMID 42778721 | 1.384 donör, 8 sınıf + 27 alt sınıf tek çekirdek cis-eQTL | yaklaşık 3.000-3.500 | Synapse `syn61929919`; ücretsiz hesap ve **download** yetkili token |
| SingleBrain 2025 meta-analizi | tek çekirdek eQTL meta-analizi | 1.000-12.000 | Zenodo 14908182; yalnız KIF13A penceresi HTTP Range ile çekildi |
| Bryois ve ark. 2022 | 8 beyin hücre tipi cis-eQTL | yaklaşık 200-1.000 | Zenodo 7276971 |
| BIG40 (UK Biobank) | L3, FA, ICVF sol ön korona radiata (UKB 25319, 25079, 25367) | GWAS, N yaklaşık 31.000 | açık, GRCh37 |
| GWAS Catalog | UDIP-FA ak madde (GCST90860937), sulkus derinliği (GCST90095129) | | açık |
| 1000 Genomes faz 3 | LD paneli (EUR ve tümü) | | açık |

Pencere: KIF13A ±1 Mb. Komşu genler (NUP153, NHLRC1, TPMT, APH1B, NDUFS2) aynı boru hattından geçti.

## 2. Doğrulama ve pozitif kontroller (geçti)

- Python ile yazılan coloc, R `coloc` 5.2.3 ile aynı sonucu veriyor (en büyük fark 2,9e-05;
  `validation_python_vs_R.tsv`).
- Yayımlanmış kolokalizasyonlar birebir yeniden üretildi (`posctrl.tsv`): Bryois Tablo S5, Alzheimer
  GWAS'ı ile oligodendrositte CR1 PP.H4 0,993, NDUFS2 0,926, ADAMTS4 0,768 (yayın 0,993 / 0,926 / 0,768).
  Zeng `coloc_class.tsv`, CR1 0,958 (yayın 0,958), ADAMTS4 0,917-0,991 (yayın 0,944).
- BIG40'ın GRCh37 olduğu test edildi (`build_test.tsv`): konum ve allel eşleşmesi GRCh37'de %70,2,
  bir baz kaydırınca %0,3.
- LD uyumu: r² ≥ 0,9 olan varyant çiftlerinde eQTL etki işaretleri neredeyse hiç çelişmiyor
  (`ld_uyum_tanisi.tsv`, oran en fazla 0,0008).

Boru hattı gerçek bir kolokalizasyonu görebiliyor. Dolayısıyla KIF13A'daki negatif sonuç, yöntemin
kör olmasından kaynaklanmıyor.

## 3. Ana sonuç (`coloc_main.tsv`)

| Kaynak | Hücre tipi | KIF13A eQTL öncüsü (p) | GWAS öncüsünde eQTL p | PP.H3 | PP.H4 |
|---|---|---|---|---|---|
| Zeng 2026 | Oligo | rs113815188 (2,7e-14) | 0,78 (BIG40 L3) | 0,98 | 0,0007 |
| Zeng 2026 | OPC | rs760849 (3,6e-37) | 0,0045 marjinal, koşullu yaklaşık 0,9 | 0,98 | 0,00008 |
| Zeng 2026 | Astrosit | rs6459581 (2,7e-09) | 0,65 | 0,98 | 0,0001 |
| Zeng 2026 | Uyarıcı nöron | rs7769287 (2,0e-17) | 0,07 | 0,97 | 0,011 |
| SingleBrain 2025 | Oligo | rs13211637 (4,7e-07) | 0,46 | 0,95 | 0,002 |
| SingleBrain 2025 | OPC | rs760849 (4,4e-17) | 0,11 | 0,99 | 0,00009 |

Değerler BIG40 L3 sol ön korona radiata içindir; UDIP-FA, ICVF, FA ve sulkus derinliği aynı tabloyu
verir (oligodendrosit ve OPC'de PP.H4 en fazla 0,015). OPC öncüsü rs760849 iki atlasta da aynı varyant,
yani OPC eQTL'i sağlam ama ak madde sinyaliyle ilgisiz.

## 4. Koşullu analiz ve mikroglia ipucu (`coloc_conditional.tsv`)

Bir hücre tipinde birden fazla eQTL sinyali varsa, en güçlüsü baskın gelip zayıf bir ortak sinyali
gizleyebilir. Bu yüzden sinyaller tek tek koşullanarak (dış LD paneliyle) yeniden test edildi.

- Oligodendrosit ve OPC: koşullu sinyallerin hiçbiri PP.H4 0,001'i geçmiyor.
- **Zeng bağışıklık (mikroglia), 3. koşullu sinyal rs73373199** (GWAS öncüsüyle r² 0,30): EUR paneliyle
  PP.H4 0,55 (BIG40 L3) ve 0,62 (UDIP-FA). Tüm 1000G paneliyle koşullama sırası değişiyor ve PP.H4 0,002
  ile 0,034'e düşüyor.
- SingleBrain bağışıklık: koşullu sinyal doğrudan rs75428660 (öncüyle r² 0,97), PP.H4 0,17, ve boru
  hattının LD uyumsuzluğu bayrağı açık.
- Bryois OPC ile sulkus derinliği PP.H4 0,535 veriyor, ama aynı test Zeng ve SingleBrain'de 1e-5
  düzeyinde. Küçük örneklemin ürünü olarak okunmalı.

Okuma: mikroglia'da ak madde haplotipine yakın bir KIF13A eQTL sinyali **olabilir**, ama sonuç
kullanılan LD paneline bağlı. Örnek içi LD ya da yazarların ince haritalama setleri olmadan iddia
edilemez.

## 5. Güç

- Bryois 2022 bu soru için yetersiz: N 1.100'de ve simüle edilen en büyük etkide (allel başına 0,388 SD)
  PP.H4 ≥ 0,8'e ulaşma gücü 0,32-0,34 (`power_sim.tsv`). Bryois'deki düşük PP.H4 değerleri bilgi taşımıyor.
- Zeng ve SingleBrain'de etkin N 3.000'in üzerinde ve oligodendrosit/OPC eQTL'leri 1e-14 ile 1e-37
  arasında; buradaki yüksek PP.H3 değerleri güç eksikliği değil, iki sinyalin gerçekten ayrı olduğunu
  gösteriyor.
- rs76091407 (UDIP-FA öncüsü) Zeng'in varyant kümesinde yok; UDIP-FA testleri pencerenin kendi
  öncüsüyle yapıldı.

## 6. Aşırı iddia uyarıları

1. **"KIF13A ak maddeyi oligodendrosit üzerinden etkiliyor."** Bu veriyle desteklenmiyor; tersine iki
   büyük atlas iki sinyalin ayrı olduğunu söylüyor.
2. **"Mikroglia'da kolokalizasyon bulundu."** Hayır; LD paneli değişince kayboluyor.
3. **"Ak madde haplotipinin KIF13A ile ilgisi yok."** Fazla. MetaBrain korteks eQTL bağlantısı (bulk doku,
   r² 0,97) duruyor; bu analiz yalnız ölçülen tek çekirdek hücre tiplerinde ortak varyant olmadığını
   gösteriyor. Etki ölçülmeyen bir hücre durumunda, gelişim döneminde ya da başka bir gen üzerinden
   olabilir.

## 7. Sonraki adım

Mikroglia ipucunu örnek içi LD ile test etmek ya da kapatmak: yazarların ince haritalama (credible set)
çıktıları varsa onlarla, yoksa LD gerektirmeyen bir yaklaşımla. Sonuç ne olursa olsun bu belgeye
eklenmeli.

## 8. Bu eksende öğrenilen tuzaklar

1. Synapse token'ı okuyabilip indiremeyebilir: izin kaydı `canDownload: true` derken her dosya 403
   döndü; etkin token'da yalnız `view` yetkisi vardı. Yetkileri `GET auth/v1/personalAccessToken` ile
   kontrol et (yalnız ad ve yetki listesini yazdır).
2. Dış LD paneliyle koşullu kolokalizasyon kırılgandır; PP.H4 panel değişince değişiyorsa hipotezdir.
3. Kendi yazdığın coloc'u hem referans uygulamayla hem yayımlanmış bir sonuçla doğrula; doğrulanmamış bir
   boru hattının negatif sonucu yorumlanamaz.
4. SingleBrain dosyaları büyük; tamamını işleyen iş bellek yetersizliğinden düştü (iş 6439632).
   Yalnız gen penceresini HTTP Range ile çekmek (`s06_singlebrain_stream.py`) çözdü.
