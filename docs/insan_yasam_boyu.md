# KIF13A insan yaşam boyu ekseni (B ekseni)

**Tarih:** 3-4 Ekim 2026. Hesaplar TRUBA'da SLURM işleri olarak çalıştı (`/arf/scratch/avural/kif13a_insan_yas`),
login düğümünde hesap yapılmadı. Script'ler `scripts/truba/insan_yas/`, küçük sonuç tabloları
`data/interim/insan_yas/` altında. Yalnız canlı doğrulanmış PMID, DOI ve erişim numarası kullanıldı
(PubMed ve Europe PMC üzerinden, 3 Ekim 2026). Hiçbir veri kullanım sözleşmesi kabul edilmedi, hesap açılmadı.

## 0. Hüküm: fare yaşlanma bulgusu insanda tekrarlanmadı

**H_B(1), yaşlanma: tekrarlanmadı.** Farede olgun oligodendrositte görülen +0,250 log2 (1,19x) yaş artışı insan
olgun oligodendrositinde yok. Asıl test, Nature 2026 yaşam boyu atlasının donörleri (PsychAD "Aging_Cohort",
CELLxGENE üzerinden açık): 232 yetişkin donörde KIF13A psödobulk ifadesi on yılda **-0,018 log2** (SE 0,007,
p = 0,017); genç yetişkin (20-39 yaş, n = 55) ile yaşlı (60 yaş ve üstü, n = 82) arasında fark **-0,059 log2
(0,96x)**, Mann-Whitney p = 0,035, Cliff's delta -0,21. Yön farenin tersi. Bağımsız üç kohortun meta-analizi
(PsychAD, SEA-AD MTG, Seeker ak madde; 335 donör) aynı küçük düşüşü verir: on yılda -0,018 (p = 0,002, I² = 0),
yaşlı eksi genç -0,057 log2, %95 GA [-0,096, -0,017]. **Fare büyüklüğündeki artış (+0,25) bu güven aralığının
çok dışında ve %80 güçle saptanabilir en küçük etki 0,085 log2:** bu bir güç eksikliği değil, gerçek bir
negatiftir. Gözlenen küçük düşüş ise KIF13A'ya özgü değil: aynı modelde ifade edilen 18.985 genin arasında
empirik p = 0,078, ifade düzeyi eşleştirilmiş genler arasında 0,26; 90 yaş üstü 9 donör çıkarılınca p = 0,10'a
düşüyor. Bunu "insanda KIF13A yaşla azalır" diye okumak aşırı iddia olur; doğru ifade "artmıyor" dur.
OPC'de fare gibi değişim yok (on yılda -0,009, p = 0,48; meta p = 0,41).

**H_B(2), oligodendrosit soyu birinciliği: tekrarlandı, bir farkla.** Üç insan veri setinde (PsychAD 26 alt
sınıf, ak madde 10 tip, postnatal multiome 10 tip) oligodendrosit soyu psödobulk, hücre ortalaması ve derinliğe
göre standartlaştırılmış tespit oranında ilk iki sırayı alıyor. Fark: insanda **OPC birinci, olgun oligodendrosit
ikinci** (farede tersi). Astrosite karşı olgun oligodendrositte tespit oranı derinlik eşitlenince 1,25x, tespit
edilen hücrede düzey 2,0x (fare: 1,8x ve 2,0x). Düzey bileşeni birebir tekrarlanıyor, tespit bileşeni insanda
daha küçük.

**H_B(3), gelişim: OPC'de erken çocuklukta artış var ama özgül değil.** KIF13A, OPC'de bebeklikten genç
yetişkinliğe yükseliyor ve bu üç bağımsız kaynakta görülüyor (PsychAD Spearman rho = 0,43, p = 0,001, n = 56;
Velmeshev rho = 0,42, p = 0,023, n = 29; Herring rho = 0,82, p < 0,001, n = 20). Değişimin yarısı yaklaşık
3-5 yaşa kadar tamamlanıyor, yani erken postnatal OPC olgunlaşmasıyla (PDGFRA düşüşü, MOG/OPALIN/BCAS1 artışı
aynı pencerede) örtüşüyor, ergenliğe uzanan geç prefrontal dönemle değil. Ama aynı artış KIF13B'de ve
astrositte KIF13A'da da var; olgun oligodendrositte hiçbir kohortta yok. Bu, miyelinleşmeye özgü bir program
değil, genel bir glial olgunlaşma izi.

**Erişim:** Atlasın kendi yaş istatistikleri (dreamlet özet tabloları, psödobulk `.rds`, PMI içeren donör
tablosu) Synapse'te sözleşmesiz (`OPEN`) ama indirmek için Synapse hesabı gerekiyor; ham FASTQ sözleşmeli. Ham
sayımlar CELLxGENE'de açık olduğundan analiz DUA olmadan yapılabildi; eksik kalan PMI düzeltmesi ve yazarların
kendi KIF13A satırıyla karşılaştırma (bkz. 8).

---

## 1. Bu eksen neden var

HANDOFF.md 3.6 ve `docs/ak_madde_goruntuleme.md` 1: KIF13A'nın hücre tipi izi oligodendrosit soyunu gösteriyor
ve farede olgun oligodendrositte yaşla artıyor (`C:\Users\pc\Claude\kif13a-allenbrain`, 53 genç ve 44 yaşlı
fare, fare başına tek değer: +0,250 log2, p = 0,0002, Cliff's delta +0,435; OPC +0,131, p = 0,22; ifade eden
hücre oranı değişmiyor). İnsan delesyon vakalarının ortak fenotipi konuşma ve dil gecikmesi olduğundan,
gelişimsel seyir de soruldu. Üç hipotez:

- **H_B(1)** İnsanda da KIF13A olgun oligodendrositte yaşla artar, OPC'de artmaz.
- **H_B(2)** Oligodendrosit soyu birinciliği tespit/düzey ayrımı yapıldığında da durur.
- **H_B(3)** (keşif) KIF13A'nın OPC ve oligodendrositteki bebeklik-çocukluk-ergenlik seyri miyelinleşme
  dönemiyle örtüşür mü.

## 2. Veri setleri

Hepsi CELLxGENE Discover'dan, şema 7.1.0, son revizyon 11 Haziran 2026, indirme 3 Ekim 2026
(`data/interim/insan_yas/DOWNLOADS.tsv`, boyutlar uzak sunucuyla bayt bayt eşleşti). Her dosyada ham sayımlar
tamsayı olarak doğrulandı (`pseudobulk_meta/<set>/obs_summary.txt`): PsychAD'de ham sayım `X`'te, diğerlerinde
`raw/X`'te; diğer setlerin `X`'i normalize (0,16-10,6 arası ondalık), kullanılmadı.

| Set | Yayın (canlı doğrulandı) | CELLxGENE dataset_version_id | Bölge, yöntem | Donör (analiz) | Yaş |
|---|---|---|---|---|---|
| PsychAD yaşam boyu ("Aging_Cohort") | Yang 2026 Nature, PMID 42778698, doi 10.1038/s41586-026-10271-7; koleksiyon: Lee 2024 medRxiv, PMID 41602098, doi 10.1101/2024.10.31.24316513 | 13e8f1dd-962f-47b7-9cf7-c71d2b21b8a5 | DLPFC, 10x 3' v3, çekirdek | 291 (yetişkin 232; <20 yaş OPC 56, Oligo 43) | 1 ay-90+ |
| SEA-AD MTG Oligo / OPC | Gabitto 2024 Nat Neurosci, PMID 39402379, doi 10.1038/s41593-024-01774-5 | 021a2e09-19c6-42eb-9e52-3ed2cf3d6032 / 01a8996c-8817-43c5-ac45-d66342ddb97e | MTG, 10x v3 + multiome | 83 / 84 (46'sı 90+) | 65-90+ |
| SEA-AD DLPFC Oligo / OPC | aynı | d847bbb8-9be5-473c-9f3a-4af707a2cd05 / 92762ef7-6cb4-4946-ac81-f54c640a72a4 | DLPFC, aynı donörler | 80 / 80 | 68-90+ |
| Ak madde ("Brain matters") | Seeker 2023 Acta Neuropathol Commun, PMID 37217978, doi 10.1186/s40478-023-01568-z | a08ce2b5-813f-4596-832b-088f02efb55d | BA4 ak madde, serebellum ak madde, servikal omurilik; 10x v3 | 20 (10 genç 34-45, 10 yaşlı 61-74) | 34-74 |
| Postnatal multiome | Clarence 2025 Nat Genet, PMID 39962241, doi 10.1038/s41588-025-02083-8 | e6ef2a07-1b8e-49a8-a771-15b81971eac7 | 4 bölge, 10x multiome | 10 (yalnız tanımlayıcı) | bebek-62 |
| Velmeshev entegrasyonu | Velmeshev 2023 Science, PMID 37824647, doi 10.1126/science.adf0834 | 2001fbc7-de84-4e84-9a76-4290bcd65415 | korteks, 4 kaynak çalışma | "Velmeshev" alt kümesi 54, "Herring" alt kümesi 27 (gelişim) | doğum öncesi-54 |

**PsychAD verisinin açık olduğunun kanıtı.** Atlasın GitHub deposu (DiseaseNeuroGenomics/DLPFC_Lifespan_Atlas,
README) etkileşimli görüntüleme için `cellxgene.cziscience.com/e/4442d412-91cb-4261-acca-8adf5fa04c11.cxg/`
bağlantısını veriyor; bu, CELLxGENE'nin herkese açık dizinindeki "Aging_Cohort" veri setidir (dataset_id
4442d412-..., indirilen sürüm 13e8f1dd-...). CELLxGENE sürümünde 291 donör ve 1.332.155 çekirdek var; makale
kalite kontrolünden sonra 284 donör ve 1.307.674 çekirdek bildiriyor. Fark (7 donör) açıklanmadı; analiz
CELLxGENE sürümüne aittir.

**Kullanılamayan / tekrar seti olamayanlar.**
- PsychAD HBCC kohortu (300 donör) indirildi ve psödobulk'u çıkarıldı, ama tanısız ("normal") 146 donörün
  **146'sı da** yaşam boyu kohortunda. Kalan 122 donör psikiyatrik tanılı. Bağımsız tekrar seti olamaz.
- PsychAD MSSM ve RADC kohortları (36 GB ve 6 GB) indirilmedi: yaşam boyu kohortu zaten bunların alt kümesi
  (`is_primary_data = False`) ve çoğu donör tanılı.
- Makalenin "306 nörotipik donörlük tekrar seti" üç yayımlanmış veri setinin yazarlarca entegrasyonu;
  indirilebilir bir nesne olarak bulunmadı.
- SEA-AD'nin 5 "Neurotypical reference" donörü (29-60 yaş) PMI'leri "Reference" olarak verildiği ve ayrı bir
  kohorttan geldiği için yetişkin modelinden düştü; SEA-AD yalnız 65 yaş üstü eğimi test ediyor.

## 3. Yöntem

**Psödobulk** (`pseudobulk.py`): h5ad diskten parça parça okunur, ham sayımlar donör x hücre tipi (çok bölgeli
setlerde x bölge, SEA-AD'de x yöntem, Velmeshev'de x kaynak çalışma) gruplarına toplanır. Grup başına: tüm
genler için toplam sayım ve sayımı > 0 olan hücre sayısı; hücre başına UMI ve gen sayısı; odak genler için
hücre düzeyinde sayım. Bir donör bir hücre tipinden en az **20 çekirdek** varsa teste girer.

**Dört ölçü** (her odak gen için, donör başına tek değer):
- `pb`: psödobulk log2(CPM+1) = log2(toplam sayım / toplam UMI x 1e6 + 1). **Birincil ölçü.**
- `cellmean`: hücrelerin log2(CPM+1) ortalaması, sıfırlar dahil. Fare analizinin metriği.
- `det`: tespit oranı (en az bir sayım alan çekirdek oranı).
- `lvl`: yalnız tespit edilen çekirdeklerde ortalama log2(CPM+1) (tespit edilen hücrede düzey).
- Ek: Poisson beklentisi `det_expected` = ortalama(1 - exp(-p x UMI_hücre)), p grubun psödobulk oranı;
  `det_obs_over_exp` gözlenen/beklenen tespit.

**Hücre tipi tanımları (her setin kendi anotasyonu).** PsychAD: yazar `subclass` "Oligo" (alt tipler ör.
Oligo_OPALIN) ve "OPC"; karşılaştırma için "Astro", "Micro". SEA-AD: Allen taksonomisi `Subclass`
"Oligodendrocyte" ve "OPC" (süpertipler OPC_1, OPC_2, OPC_2_1-SEAAD, OPC_2_2-SEAAD). Ak madde: yazar
`broad_cell_type` "Oligo", "OPC", "Astrocyte", "Microglia_Macrophages". Postnatal ve Velmeshev: gönderenin CL'ye
eşlediği `cell_type` "oligodendrocyte", "oligodendrocyte precursor cell". Farede karşılıkları "327 Oligo NN" ve
"326 OPC NN".

**Yaş.** `development_stage` etiketinden yıl; "80 year-old and over stage" (90 yaş ve üstü) 90 olarak kodlandı
ve her sette bu donörler çıkarılarak duyarlılık analizi yapıldı; "unknown" dışlandı; "infant stage" 0,5.
Genç/yaşlı tanımı atlasın kendi gruplarından: genç yetişkin 20-39, yaşlı 60 ve üstü (ak madde setinde yazarların
tasarımı 30-45 ve 60-75). Fare karşılığı varsayımdır: 2 aylık fare genç yetişkin, 18 aylık fare yaşlı insan.

**Modeller.** Yetişkin (20 yaş ve üstü) donörlerde EKK: ölçü ~ yaş (on yıl) + cinsiyet + ortalama log(gen
sayısı) [+ PMI, varsa] [+ set özgü ortak değişken]. PsychAD'de PMI h5ad'de yok, beyin bankası (Source: HBCC/MSSM)
ortak değişken. SEA-AD'de PMI aralık orta noktası ve ADNC (0-3). Ak madde setinde PMI `tissue_handling_interval`.
Çok bölgeli setlerde değerler bölge içinde merkezlenip donör başına ortalanır (bölge sabit etkisi); test birimi
her zaman donör. İkili karşılaştırma: Mann-Whitney, Welch t, Cliff's delta. Meta-analiz: ters varyans sabit etki
ve DerSimonian-Laird rastgele etki; SEA-AD MTG ve DLPFC aynı donörler olduğu için ana meta'da MTG, duyarlılıkta
DLPFC.

## 4. Sonuçlar

### 4.1 Yaşlanma, fareyle yan yana (olgun oligodendrosit ve OPC, KIF13A)

| | Fare (Allen, 53 vs 44) | İnsan PsychAD (55 vs 82) | İnsan ak madde (10 vs 10) | İnsan SEA-AD MTG (65-90+, eğim) | Meta (bağımsız 3 set) |
|---|---|---|---|---|---|
| **Oligo, yaşlı eksi genç** | **+0,250 log2 (1,19x)**, p = 0,0002, Cliff +0,435 | **-0,059 (0,96x)**, p = 0,035, Cliff -0,21 | -0,055 (0,96x), p = 0,076, Cliff -0,48 | yok (genç grup yok) | **-0,057**, GA [-0,096, -0,017], p = 0,005 |
| Oligo, on yıllık eğim | | -0,018, p = 0,017 | -0,020, p = 0,072 | -0,006, p = 0,87 | **-0,018**, p = 0,002, I² = 0 |
| **OPC, yaşlı eksi genç** | +0,131 (1,09x), p = 0,22, Cliff +0,146 | -0,045, p = 0,28 | +0,075, p = 0,73 | yok | -0,026, p = 0,58 |
| OPC, on yıllık eğim | | -0,009, p = 0,48 | +0,012, p = 0,73 | -0,095, p = 0,13 | -0,009, p = 0,41, I² = 0,12 |
| Oligo, 90+ çıkarılınca | | eğim -0,013, p = 0,10; fark -0,043, p = 0,10 | (90+ yok) | -0,029, p = 0,54 (n = 37) | |

İnsan değerleri `pb` ölçüsü (psödobulk log2 CPM), kovaryat düzeltilmiş eğim; ikili fark ham ortalama farkı.
Kaynak: `summary/combined/kif13a_age_table.tsv`, `meta_age.tsv`.

Fare metriğiyle (`cellmean`) PsychAD olgun oligodendrositte fark -0,026 (p = 0,51), OPC'de +0,032 (p = 0,96):
yön yine artış değil, ama bu ölçü çekirdek verisinde çok gürültülü (bkz. 6).

### 4.2 Tespit oranı ve tespit edilen hücrede düzey (yaşlanma)

PsychAD olgun oligodendrosit, genç (n = 55) ve yaşlı (n = 82):
- Tespit oranı 0,852 ve 0,856: **değişmiyor** (p = 0,89), farede olduğu gibi.
- Tespit edilen hücrede düzey 9,03 ve 8,96 log2: **-0,073**, p = 0,0009. Farede bu bileşen artıyordu; insanda
  hafif azalıyor.
- Gözlenen/Poisson beklenen tespit oranı on yılda -0,003 (p < 0,001): yaşla tespitteki değişim neredeyse
  tamamen düzey ve derinlikle açıklanıyor, ek bir "açılma/kapanma" yok.

OPC'de her iki bileşen de değişmiyor (tespit +0,010, p = 0,20; düzey -0,066, p = 0,18).

### 4.3 Meta-analiz

`meta_age.tsv`. Ana set PsychAD + SEA-AD MTG + ak madde. Olgun oligodendrosit KIF13A eğimi -0,0180 (SE 0,0058),
p = 0,0021, heterojenlik yok (Q p = 0,92). SEA-AD DLPFC ile değiştirilince sonuç aynı (-0,0180, p = 0,0021). OPC
-0,0093 (SE 0,0113), p = 0,41. Paralog KIF13B olgun oligodendrositte daha güçlü azalıyor (-0,034, p = 1,3e-05,
I² = %55; yaşlı eksi genç -0,153), OPC'de değişmiyor (+0,021, p = 0,15).

### 4.4 Hücre tipi sıralaması, tespit/düzey ayrımı (H_B(2))

Yalnız yetişkin donörler, donör başına ortalanmış. "Derinlik standardize" = tespit oranı ya da düzey her UMI
ondalık diliminde ayrı hesaplanıp tüm hücrelerin UMI dağılımına göre ağırlıklandırıldı (`ranking_all.tsv`).

| Set | Tip | Sıra (pb / cellmean / derinlik std tespit) | Astrosite karşı pb | Naif cellmean farkı | Tespit (ham) | Tespit (derinlik std) | Düzey (derinlik std) |
|---|---|---|---|---|---|---|---|
| PsychAD (26 alt sınıf) | OPC | 1 / 1 / 1 | 4,8x | 16,3x | 1,43x | 1,33x | 3,4x |
| PsychAD | Oligo | 2 / 2 / 2 | 2,6x | 5,1x | 1,27x | 1,25x | 2,0x |
| Ak madde (10 tip) | OPC | 1 / 1 / 1 | 4,6x | 42,9x | 2,33x | 1,64x | 2,2x |
| Ak madde | Oligo | 2 / 2 / 2 | 2,1x | 7,6x | 1,80x | 1,40x | 1,4x |
| Postnatal (10 tip, 4 yetişkin) | OPC | 1 / 1 / 1 | 4,1x | 13,4x | 1,50x | 1,51x | 2,6x |
| Postnatal | Oligo | 2 / 2 / 2 | 2,8x | 5,9x | 1,36x | 1,42x | 1,9x |
| Fare (Allen, 265 alt sınıf) | Oligo | 1 | | ~16x | 1,8x | | 2,0x |

Fare deposunun uyarısı insanda da geçerli: "naif" fark (log ortalamalarının farkından üs alma) 5-43x gösterirken
gerçek ayrım tespitte 1,3-1,6x, düzeyde 1,4-3,4x. Derinlik etkisi gerçek: PsychAD'de nöronların medyan UMI'si
yaklaşık 25.000, olgun oligodendrositin 5.097; ham tespit oranında iki nöron alt sınıfı olgun oligodendrositin
önüne geçiyor (sıra 4), derinlik eşitlenince oligodendrosit 2. sıraya dönüyor.

### 4.5 Gelişim (H_B(3))

PsychAD, KIF13A psödobulk log2 CPM medyanı (donör sayısı):

| Yaş | 0-1 | 1-5 | 6-11 | 12-19 | 20-39 | 40-59 | 60+ |
|---|---|---|---|---|---|---|---|
| OPC | 9,78 (11) | 9,97 (7) | 9,98 (4) | 10,03 (34) | 10,04 (54) | 10,01 (95) | 9,96 (83) |
| Oligo | yok* | 9,16 (6) | 9,08 (4) | 9,11 (33) | 9,09 (55) | 9,07 (95) | 9,01 (82) |

\* 1 yaş altında hiçbir donörde 20 olgun oligodendrosit çekirdeği yok (bu yaşta "Oligo" alt sınıfının çekirdek
oranı medyanı 0).

20 yaş altı Spearman korelasyonu ve LOWESS eğrisinde bebeklikten genç yetişkin düzeyine olan değişimin yarısının
tamamlandığı yaş (`timing_all.tsv`):

| Gen, tip | PsychAD rho (p), yarı yaş | Velmeshev rho (p), yarı yaş | Herring rho (p), yarı yaş |
|---|---|---|---|
| KIF13A, OPC | 0,43 (0,001), 2,8 | 0,42 (0,023), 0,6 | 0,82 (<0,001), 4,8 |
| KIF13A, Oligo | -0,11 (0,50) | 0,37 (0,068) | 0,14 (0,64) |
| KIF13A, Astro | 0,42 (0,002), 2,9 | 0,29 (0,13) | 0,78 (<0,001), 1,3 |
| KIF13B, OPC | 0,30 (0,026) | 0,53 (0,003) | 0,83 (<0,001) |
| PDGFRA, OPC | -0,48 (<0,001), 3,9 | -0,17 (0,38) | -0,09 (0,72) |
| MOG, OPC | 0,45 (0,001), 2,4 | | |
| OPALIN, OPC | 0,28 (0,040), 4,4 | 0,67 (<0,001), 4,4 | 0,72 (<0,001), 4,5 |
| BCAS1, OPC | 0,39 (0,003), 1,7 | -0,30 (0,11) | -0,14 (0,57) |

Okuma: OPC'de KIF13A artışı üç kaynakta da var, PsychAD'de bebeklikten genç yetişkine +0,25 log2 (1,19x).
Zamanlama erken çocukluk (yarısı 3-5 yaşa kadar), atlasın kendi bileşim analizinde bulduğu 24 yaş geçiş
noktasından çok önce. Olgun oligodendrositte gelişimsel değişim yok. Astrositte ve KIF13B'de benzer artış
olduğu için bu bir "miyelinleşme geni" imzası değil.

## 5. Kontroller

### 5.1 Pozitif kontroller (boru hattı bilinen yaş etkisini buluyor mu)

Seeker 2023 (PMID 37217978), ak maddede: yaşlı OPC'lerde PLP1 ve PDGFRA düşük; yaşlı mikrogliada CD74, NEAT1,
HLA-C yüksek, CX3CR1 ve P2RY12 düşük. Not: o çalışma hücre düzeyinde test etti; burada donör düzeyinde
tekrarlandı. Ak madde setinde bu kontrol döngüseldir (aynı veri), bağımsız test PsychAD'dir.

| Kontrol | PsychAD eğim / on yıl (p) | PsychAD yaşlı eksi genç (p) | Ak madde yaşlı eksi genç (p) | Beklenen yön tuttu mu |
|---|---|---|---|---|
| OPC PDGFRA | -0,101 (2,8e-08) | -0,299 (1,1e-04) | -0,291 (0,045) | evet, meta -0,298, p = 2,3e-05 |
| OPC PLP1 | -0,026 (0,36) | -0,253 (0,028) | -0,207 (0,38) | ikili testte evet |
| Mikroglia CD74 | +0,129 (6,1e-08) | +0,635 (1,9e-08) | +0,318 (0,19) | evet |
| Mikroglia NEAT1 | +0,108 (0,0012) | +0,662 (1,7e-06) | +0,344 (0,19) | evet |
| Mikroglia HLA-C | +0,159 (5,4e-05) | +0,808 (1,9e-06) | +0,628 (0,19) | evet |
| Mikroglia CX3CR1 | -0,182 (1,9e-04) | -0,884 (7,7e-06) | -0,493 (0,11) | evet |
| Mikroglia P2RY12 | -0,139 (4,3e-04) | -0,649 (5,5e-05) | -0,542 (0,13) | evet |

Yedi kontrolün yedisi PsychAD'de beklenen yönde; ikili testte yedisi, yaş eğiminde altısı p < 0,05. Boru hattı, KIF13A'dan çok daha küçük
örneklemle bile gerçek yaş etkilerini yakalıyor. Fareden türetilen yaşlanma oligodendrosit işaretleri (Kaya 2022,
PMID 36280798, doi 10.1038/s41593-022-01183-6: C4b, Serpina3n; çoğunlukla ak maddede) insan korteks
oligodendrositinde değişmedi (C4A, C4B, STAT1 p > 0,5; B2M +0,043, p = 0,11; SERPINA3 PsychAD matrisinde yok).
Bu tür-ötesi bir öngörüydü, pozitif kontrol olarak çalışmadı ve öyle sunulmamalı.

Gelişim ekseni için pozitif kontrol: atlas olgun oligodendrosit oranının yaşla logaritmik arttığını, OPC oranının
azaldığını bildiriyor (Extended Data Fig. 1). Burada aynı veride Spearman (log2(yaş+1) ile) oligo oranı
rho = +0,33 (p = 9e-09), OPC oranı rho = -0,39 (p = 4e-12), n = 290 donör: tekrarlandı.

### 5.2 Negatif kontrol: tüm genlerin yaş dağılımı

Aynı model, aynı donörler, ifade edilen tüm genler (ortalama CPM >= 1, donörlerin >= %80'inde sayım > 0);
ayrıca KIF13A'ya ±0,5 log2 ifade düzeyinde eşleştirilmiş en fazla 1000 rastgele gen (`controls_negative_null.tsv`).

| Set, tip | Gen sayısı | KIF13A t | Empirik p (tüm genler) | Empirik p (eşleşmiş) | p < 0,05 olan gen oranı |
|---|---|---|---|---|---|
| PsychAD Oligo | 18.985 | -2,41 | 0,078 | 0,26 (221 gen) | %13,3 |
| PsychAD OPC | 19.061 | -0,71 | 0,59 | 0,77 | %17,5 |
| SEA-AD MTG Oligo | 17.219 | -0,17 | 0,92 | 0,94 | %17,3 |
| Ak madde Oligo | 12.967 | -1,94 | 0,044 | 0,034 | %2,8 |

PsychAD'de genlerin %13'ü yaşla nominal anlamlı değişiyor; KIF13A'nın küçük düşüşü bu genel kaymanın içinde.

### 5.3 Paralog KIF13B

Olgun oligodendrositte KIF13B de azalıyor ve daha güçlü (meta -0,034 / on yıl; PsychAD yaşlı eksi genç -0,163,
p = 0,0001). OPC'de değişmiyor. Gelişimde OPC'de KIF13A ile aynı yönde artıyor. Yani yaşlanmada da gelişimde de
KIF13A'ya özgü bir davranış görülmüyor.

### 5.4 Derinlik (hücre başına yakalanan gen sayısı)

Her modelde ortalama log(gen sayısı) ortak değişken; ayrıca modeller onsuz da kuruldu (`beta_no_depth`).
PsychAD olgun oligodendrositte derinlik yaşla ilişkili değil (rho = 0,095, p = 0,15) ve KIF13A eğimi derinlik
düzeltmesiyle ve düzeltmesiz neredeyse aynı (-0,0177 ve -0,0158). OPC, astrosit ve mikrogliada derinlik yaşla
artıyor (rho 0,19-0,29), bu yüzden `det` ve `lvl` ölçüleri bu tiplerde derinlik düzeltmesine duyarlı (bkz.
`depth_checks_all.tsv`); birincil `pb` ölçüsü değil.

## 6. Güç

Fare etkisi: genç yetişkin ile yaşlı arasında +0,25 log2. Fare olgun oligodendrositinde tespit oranı 0,93 ve yaşla
değişmediği için bu fark pratikte ifade düzeyinde 1,19x'lik bir artıştır; insanda karşılığı psödobulk log2 CPM'de
+0,25'tir. Doğrusal birikim varsayılırsa PsychAD'in genç ve yaşlı grupları arasındaki 44,7 yıllık farkla on yılda
+0,056 log2 eder. Gözlenen artık varyansla ve gerçek yaş dağılımıyla (`power_vs_mouse.tsv`):

| Set, tip | n | Fare etkisini görme gücü (ikili, +0,25) | %80 güçle en küçük fark | Eğim gücü (+0,056/on yıl) | %80 güçle en küçük eğim |
|---|---|---|---|---|---|
| PsychAD Oligo | 232 | **>0,99** | 0,085 log2 | >0,99 | 0,021 |
| PsychAD OPC | 232 | >0,99 | 0,143 | >0,99 | 0,034 |
| Ak madde Oligo | 20 | >0,99 | 0,075 | >0,99 | 0,028 |
| Ak madde OPC | 20 | 0,51 | 0,330 | 0,46 | 0,099 |
| SEA-AD MTG Oligo | 83 | (genç yok) | | 0,75 | 0,095 |
| SEA-AD MTG OPC | 84 | | | 0,31 | 0,175 |
| Meta Oligo | 335 | >0,99 | | >0,99 | |

Olgun oligodendrosit için negatif sonuç güç eksikliğiyle açıklanamaz: fare büyüklüğünde bir artış olsaydı
PsychAD tek başına neredeyse kesin görürdü ve meta güven aralığının üst sınırı (-0,017) +0,25'in çok altında.
Fare metriğiyle (`cellmean`) güç yalnız 0,54 (en küçük fark 0,34): çekirdek verisinde hücre başına 1-2 sayımla
bu metrik derinlik gürültüsüne boğuluyor, o yüzden birincil ölçü olarak kullanılmadı.

## 7. Aşırı iddia uyarıları

1. **"İnsanda KIF13A yaşla azalır" denemez.** Düşüş küçük (0,96x), genom geneli dağılımda sıradan (empirik
   p = 0,078 ve eşleşmişte 0,26) ve 90 yaş üstü 9 donör çıkarılınca p = 0,10. Söylenebilen: fare tipi bir artış
   insanda yok.
2. **Doku ve hücre farkı.** Fare: tüm beyin bölgelerinden tek hücre (10x v3, sitoplazma dahil). İnsan: çoğunlukla
   prefrontal ve temporal gri madde, çekirdek RNA'sı. KIF13A mRNA'sının uzantılardaki ya da sitoplazmadaki payı
   çekirdek verisinde görünmez.
3. **Yaş eşlemesi varsayımdır.** 18 aylık fare ile 60 yaş üstü insan eşleştirmesi bir tercih; SEA-AD yalnız
   65-90+ aralığını, ak madde seti yalnız 34-74'ü kapsıyor.
4. **Eksik ortak değişkenler.** PsychAD'de PMI ve RIN h5ad'de yok (atlasın kendi modeli PMI içeriyor); 90 yaş
   üstü yaşlar sansürlü (90 kodlandı).
5. **CPM görelidir.** Hücre başına toplam mRNA yaşla değişirse oran tek başına kayar.
6. **Gelişim bulgusu özgül değil ve az donöre dayanıyor.** 1-11 yaş arasında PsychAD'de 4-7 donör; bebeklikteki
   "OPC" sınıfı olgunlaşmamış durumları da içerir, sınıf içi bileşim kayması ifade kayması gibi görünebilir.
   Astrosit ve KIF13B aynı eğrileri çiziyor.
7. **Sıralama tanımlayıcıdır.** Donör başına ortalanmış ama hücre tipleri arasında resmi test yok; "düzey"
   ölçüsü düşük sayım rejiminde derinlikle iç içe, derinlik dilimlemesi bunu yalnız yaklaşık düzeltir.
8. **SEA-AD bir Alzheimer kohortudur.** ADNC düzeltildi ama yaş ile nöropatoloji iç içe; donörlerin %55'i 90+.

## 8. Sözleşme ya da hesap olmadan yapılamayanlar

Synapse REST API ile 3 Ekim 2026'da sorgulandı (anonim):

| Synapse nesnesi | İçerik | Kısıt düzeyi | Anonim indirme |
|---|---|---|---|
| syn52396927 | Proje (DLPFC Lifespan Atlas) | | görüntülenebilir |
| syn77643981 | `aging_umap40pca_age_synapse_clean.h5ad` | OPEN, erişim şartı yok | hayır (Synapse girişi gerekir) |
| syn62064718 | `lifespan_pseudobulk.rds` | OPEN | hayır |
| syn72107266 | Supplementary Data 1, donör metadatası (PMI dahil) | OPEN | hayır |
| syn62062069, syn62063803 | Supplementary Data 4 ve 7, dreamlet yaş DEG özet istatistikleri | OPEN | hayır |
| syn57373409 | Ham FASTQ (1.494 PsychAD donörü) | CONTROLLED_BY_ACT, şartlar 5592528 (ADKP Clickwrap) ve 9603055 (AD Knowledge Portal ACT onaylı) | hayır |
| syn68498949 | Donör adı dönüştürücü | CONTROLLED_BY_ACT, aynı iki şart | hayır |

GitHub deposunda (DiseaseNeuroGenomics/DLPFC_Lifespan_Atlas) yalnız kod var, gen düzeyinde sonuç tablosu yok.

Sonuç: bu soru için **DUA gerekmiyor**, ham sayımlar CELLxGENE'de açık. Eksik kalan iki şey yalnız **ücretsiz
bir Synapse hesabı** ister (sözleşme değil, Synapse kullanım koşulları): (a) Supplementary Data 1'den PMI ile
modeli yazarların kovaryat setine yaklaştırmak, (b) Supplementary Data 4/7'deki yazarların kendi dreamlet
sonucunda KIF13A'nın Oligo ve OPC satırlarını bu analizle karşılaştırmak. Ham FASTQ gerekmiyor.

## 9. Yeniden üretim

TRUBA'da, `BASE=/arf/scratch/avural/kif13a_insan_yas`, script'ler `$BASE/scripts` (yerel kopya
`scripts/truba/insan_yas/`, birebir aynı). Her adım bir SLURM işi; login düğümünde yalnız `sbatch`.

```bash
cd /arf/scratch/avural/kif13a_insan_yas
# 1. indirme (yalnız debug bölümünde internet var), anahtarlar iki nokta ile ayrılır
sbatch --export=ALL,KEYS=all scripts/download.slurm
# 2. psödobulk, set başına (GROUP'ta da ayırıcı iki nokta)
sbatch --export=ALL,KEY=psychad_aging,GROUP=donor_id:subclass,WORKERS=6 scripts/pseudobulk.slurm
sbatch --export=ALL,KEY=wm_all,GROUP=donor_id:broad_cell_type:tissue,WORKERS=4 --cpus-per-task=8 --mem=16G scripts/pseudobulk.slurm
for k in seaad_mtg_oligo seaad_mtg_opc seaad_dlpfc_oligo seaad_dlpfc_opc; do
  sbatch --export=ALL,KEY=$k,GROUP=donor_id:Subclass:assay,WORKERS=6 --cpus-per-task=16 --mem=32G scripts/pseudobulk.slurm; done
sbatch --export=ALL,KEY=postnatal,GROUP=donor_id:cell_type:tissue,WORKERS=6 --cpus-per-task=16 --mem=32G scripts/pseudobulk.slurm
sbatch --export=ALL,KEY=velmeshev,GROUP=donor_id:cell_type:dataset,WORKERS=4 --cpus-per-task=8 --mem=16G scripts/pseudobulk.slurm
sbatch --export=ALL,KEY=psychad_hbcc,GROUP=donor_id:subclass,WORKERS=6 scripts/pseudobulk.slurm   # örtüşme kontrolü için
# 3. donör düzeyi istatistik
sbatch --export=ALL,KEYS=psychad_aging:wm_all:seaad_mtg_oligo:seaad_mtg_opc:seaad_dlpfc_oligo:seaad_dlpfc_opc:postnatal:velmeshev_velmeshev:velmeshev_herring --mem=60G --cpus-per-task=28 scripts/analyze.slurm
# 4. birleştirme, güç, meta-analiz
sbatch --export=ALL,SCRIPT=combine.py --mem=8G --cpus-per-task=4 scripts/analyze.slurm
```

Ortam: `conda activate allenbrain` (anndata 0.12.19, h5py 3.16.0, numpy 2.4.6, pandas 2.3.3, scipy 1.17.1); yeni
paket kurulmadı. Süreler: PsychAD psödobulk 7 dakika (28 çekirdek), diğerleri 1 dakikanın altında; disk 38 GB
(`raw/` 37 GB).

### Dosyalar

- `scripts/truba/insan_yas/`: `datasets.tsv`, `download.slurm`, `pseudobulk.py`, `pseudobulk.slurm`,
  `analyze.py`, `stats_lib.py`, `analyze.slurm`, `combine.py`.
- `data/interim/insan_yas/summary/combined/`: `kif13a_age_table.tsv`, `meta_age.tsv`, `power_vs_mouse.tsv`,
  `controls_positive.tsv`, `controls_negative_null.tsv`, `depth_checks_all.tsv`, `ranking_all.tsv`,
  `development_all.tsv`, `timing_all.tsv`, `composition_all.tsv`, `all_age_tests.tsv`.
- `data/interim/insan_yas/summary/<set>/`: `donor_values.tsv` (donör başına tüm ölçüler), `age_tests.tsv`,
  `null_age.tsv.gz` (tüm genlerin yaş eğimi).
- `data/interim/insan_yas/pseudobulk_meta/<set>/`: grup ve donör tabloları, `obs_summary.txt`, `meta.json`.
- Büyük psödobulk matrisleri (`pb_counts.npz`, `pb_detect.npz`, `focal_cells.npz`) yalnız TRUBA'da; scratch
  dönemsel silindiğinden gerekirse 2. adımdan yeniden üretilir.

### Bu eksende öğrenilen tuzaklar

1. **`sbatch --export` virgülle böler.** `KEYS=a,b,c` yalnız `a`'yı geçirir, gerisi sessizce kaybolur. Script'ler
   iki nokta ayırıcı kabul ediyor.
2. **CELLxGENE'de ham sayımın yeri sete göre değişir.** PsychAD'de `X` ham, `raw` yok; diğerlerinde `X`
   normalize, ham sayım `raw/X`'te. `pseudobulk.py` tamsayı olduğunu doğrulamadan devam etmez.
3. **"80 year-old and over stage" 90 yaş ve üstü demektir** ve SEA-AD'de donörlerin yarısından fazlası bu
   etikette. Sürekli yaş modeline koyarken sansürlü olduğunu unutma; duyarlılık analizi şart.
4. **PsychAD kohortları iç içe.** Yaşam boyu kohortu `is_primary_data = False`, yani HBCC/MSSM/RADC'deki
   çekirdeklerin kopyası; HBCC'yi "tekrar seti" sanmak aynı donörleri iki kez saymak olur.
5. **Fare metriği çekirdek verisinde gürültülü.** Hücre başına log2(CPM+1) ortalaması, tek sayımlık çekirdeklerde
   neredeyse log2(1e6/UMI) olur; tespit edilen hücrede "düzey" de aynı nedenle derinliğe bağlı. Yaş testi için
   psödobulk, tipler arası karşılaştırma için derinlik dilimlemesi kullan.
6. **Login düğümü.** Hesap içeren hiçbir adım login'de çalışmadı; donör örtüşmesi gibi küçük kontroller `cut`,
   `sort`, `comm` ile, gerisi SLURM işinde yapıldı. (Bir kez boş girdili bir `python3` çağrısı yanlışlıkla login'de
   başladı ve hiçbir şey yapmadan kapandı; tekrarlanmadı.)
