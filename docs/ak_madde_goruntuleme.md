# KIF13A ak madde görüntüleme ekseni

**Tarih:** 17 Eylül 2026. Bu belgedeki her sayı o gün canlı çekildi. Üretim: `scripts/imaging_gwas_axis.py`,
veriler `data/interim/imaging_gwas/`, her indirmenin `.prov.json` yan dosyası ve
`data/raw/PROVENANCE.jsonl` kaydı var. Yalnız canlı doğrulanmış PMID, DOI ve erişim numarası
kullanıldı.

## 0. Hüküm: zayıf destek

**Zayıf destek.** Hücre tipi ekseninin öngörüsü şuydu: KIF13A bir oligodendrosit geniyse, insandaki
izi ak maddede olmalı. Tasarlanan ve en güçlü test olan UK Biobank BIG40 taramasında (33.224 kişiye
kadar, 3.935 görüntüleme fenotipi, GRCh37) **KIF13A ±100 kb penceresinde ak madde mikroyapısında
5e-8'i geçen hiçbir ilişki yok; sınıfın en iyi p değeri 3,5e-07** (rs144963738, chr6:18.012.583, sol
ön korona radiatada L3, UKB 25319, β = -0,11, SE = 0,021, N = 31.356). Bu değer tek başına hiçbir şey
söylemez, çünkü beyin fenotipi beklenmeyen negatif kontrol lokusunda (LCT) aynı sınıfın en iyi p
değeri 7,6e-07, yani yalnızca 2,2 kat daha kötü. Ayrıca BIG40'ın kendi yerel tepe tablosunda
(-log10 p > 7,5, 10.889 tepe, tüm IDP'ler) ±250 kb penceresinde **sıfır** kayıt var. Buna karşılık
üç bağımsız gözlem öngörüyü tam olarak yok saymayı engelliyor: (i) BIG40'ta bu lokusun 3,6 % frekanslı
haplotipine dokunan 8 varyantın **hepsi** aynı ak madde fenotipini gösteriyor ve p ≤ 0,01 olan ~200
IDP'nin yarısı ak madde mikroyapısı; (ii) aynı haplotipteki rs76091407, bağımsız ve ak maddeye özel
bir GWAS'ta (PMID 42236496, N = 11.736) **p = 4e-09 ile genom çapında anlamlı bir ak madde
mikroyapı ilişkisi** taşıyor ve en yakın geni KIF13A (32,6 kb); (iii) bu haplotip, MetaBrain'in
insan korteksindeki **KIF13A cis-eQTL'inin üçüncü bağımsız sinyaliyle r² = 0,969** (rs115824472,
p = 5,1e-06, N = 2.646) ve NUP153, NHLRC1 ya da TPMT'nin hiçbir eQTL öncü varyantıyla LD içinde
değil (r² < 0,05). **Güç ifadesi:** N = 31.355'te, %80 güçle saptanabilir alt sınır tek bir IDP'nin
varyansının **%0,126'sı** (5e-8'de), %0,130'u (BIG40'ın 7,5 eşiğinde) ve %0,188'i (BIG40'ın kendi
11,1 Bonferroni eşiğinde). Haplotipte **gözlenen** etki %0,086, yani saptanabilir sınırın **altında**:
BIG40'ın negatif sonucu bu büyüklükte gerçek bir etkiyi dışlamıyor, onu görmeye gücü yok. Aynı etkiyi
5e-8'de %80 güçle görmek için yaklaşık **46.000 kişi** gerekir. En önemli tek olumsuz sonuç ise
spesifik öngörünün çürütülmesidir: AZ PheWAS'ın n = 44'lük singulum OD izi, N = 31.355'te aynı
fenotipte ±100 kb'de en iyi p = 1,4e-04 veriyor, yani gürültü.

---

## 1. Bu eksen neden var

HANDOFF.md 3.6: iki bağımsız hücre tipi analizi KIF13A'nın oligodendrosit soyu geni olduğunda
birleşiyor (insan CELLxGENE beyin: oligodendrosit soyu ortalaması 1,49'a karşı nöron tiplerinde 0,16;
fare Allen Brain Cell Atlas, 2.341.350 hücre: olgun oligodendrosit birinci, yaşlanmayla artış yalnız
orada). Oligodendrositler miyelin yapar. Dolayısıyla KIF13A varyasyonunun insanda bir beyin
karşılığı varsa, beklenen fenotip **ak madde**, KIF1A'nın yaptığı nöronal aksonal taşıma hastalığı
değil. Bu öngörü bugüne kadar insan verisine karşı hiç test edilmedi. Hastalık genetiğinden tamamen
bağımsız bir yön olduğu için, buradan gelecek olumlu ya da olumsuz sonuç projenin başka hiçbir
ekseninde bulunmayan bir bilgi taşır.

## 2. Veri kaynağı: BIG40, ne buldum ve nasıl eriştim

| | |
|---|---|
| Kaynak | Oxford Brain Imaging Genetics Server, BIG40 |
| Ana sayfa | `https://open.win.ox.ac.uk/ukbiobank/big40/` (kaydedildi, sha256 87ad7795...) |
| Sürüm | release2; sayfadaki sürüm günlüğünün son kaydı **BIG40 26/03/21** (33k birleşik özet istatistikleri 15/10/20'de eklenmiş) |
| Yayın | Smith SM ve ark., *Nat Neurosci* 2021;24:737-745, **PMID 33875891**, PMC7610742, **doi 10.1038/s41593-021-00826-4** (tam metin canlı çekildi) |
| Yapı | **GRCh37/hg19**, iki yerden doğrulandı: ana sayfada "The reference human assembly is GRCH37/hg19", PheWeb `/about` sayfasında "All positions are on GRCh37" |
| Fenotip sayısı | **3.935** (3.929 UKB IDP'si + 6 türetilmiş bağlantı özelliği), 17 kategori |
| Örneklem | keşif 22.138, tekrar 11.086, **birleşik en fazla 33.224**; IDP başına N değişiyor (dMRI için 31.355-31.356, T1 hacimleri 33.224, WMH hacmi 32.114) |
| Varyant | 17.103.079 SNP; özet istatistiklerde MAF ≥ 0,001, INFO ≥ 0,3 |
| Etki ölçeği | Fenotipler birim varyansa ölçekli, β bir alleldeki bir standart sapma değişimi |

**Arayüz envanteri (hepsi canlı denendi):**

1. **PheWeb 33k sorgu arayüzü** (`/pheweb33k/`, PheWeb 1.1.21). Birincil yol. İki uç nokta işe yaradı:
   - `api/region/<IDP>/lz-results/?filter=...`: bir IDP için verilen penceredeki tüm varyantlar.
   - `api/variant/<chr-pos-ref-alt>`: bir varyant için **tüm 3.935 IDP**.
   Yaklaşık 10 çağrı/s, 3.935 IDP'lik tam tarama ~7 dakika ve ~1 MB. Kütlesel indirmeye gerek yok.
2. **Gen sayfası** (`/pheweb33k/gene/<SYMBOL>`): gen ±100 kb için IDP başına en iyi p, p < 1e-6
   olanlar (en az ilk 4). Bir çağrı, bağımsız çapraz kontrol olarak kullanıldı.
3. **`release2/Table1.txt`**: BIG40'ın kendi yerel tepe tablosu, -log10 p > 7,5, **10.889 satır**,
   22k keşif + 11k tekrar. Tüm IDP'leri tek indirmede kapsar (0,9 MB).
4. **`release2/stats33k/NNNN.txt.gz`**: IDP başına tam özet istatistik, **her biri 328-344 MB**, tek
   parçalı gzip (BGZF senkron noktası yok, test edildi -> rastgele erişim imkânsız) ve ~0,1-0,45 MB/s.
   Yalnız **üç** IDP için kullanıldı (aşağıda gerekçesi), akış öneki okunup pencere geçilince
   durduruldu (~100-120 MB indirme). Diske yalnız çıkarılan pencere yazıldı.
5. `BIG40-IDPs_v4/IDPs.html`: IDP adları, UKB alan kimlikleri, 17 kategori, IDP başına N ve
   kalıtsallık. Sınıflandırmanın kaynağı.

**PheWeb kopyasının gerçek içeriği (varsayılmadı, ölçüldü).** `--stage probe`: dört farklı IDP'nin
yanıtında en büyük p değeri sırasıyla 0,0099 / 0,0084 / 0,0099 / 0,010 ve en küçük MAF 0,010-0,014.
Yani BIG40'ın PheWeb kopyası yalnız **p ≤ 0,01 ve MAF ≥ 0,01** olan ilişkileri tutuyor. Bunun iki
sonucu var: (a) bir penceredeki **minimum p tam olarak geri kazanılabilir**, hiçbir şey dönmeyen IDP
için minimum p > 0,01 demektir; (b) MAF 0,001-0,01 aralığı görünmez, bu yüzden kritik IDP'lerde tam
özet istatistiğe inildi.

**Sorgu yolunun tam özet istatistikle doğrulanması.** İki IDP'de pencere iki yoldan bağımsız
çıkarıldı ve birebir tuttu:

| IDP | PheWeb bölge sorgusu | Tam özet istatistik (MAF ≥ 0,001) | Öncü |
|---|---|---|---|
| wg rh intensity-contrast medialorbitofrontal (1415) | 1,1e-08 | **1,11e-08** | rs143422697 |
| OD cingulum cingulate gyrus L (2012) | 2,1e-05 | **2,09e-05** | rs13194542 |

Tam istatistikte ±250 kb penceresinde 4.945, ±100 kb'de 2.879 varyant var; PheWeb'de görünmeyen
nadir varyantların hiçbiri bu iki IDP'de öncüyü geçmiyor. Sorgu yolu güvenilir.

## 3. Fenotip sınıflandırması (p değerine bakmadan önce yazıldı)

Sınıflandırma `scripts/imaging_gwas_axis.py` içindeki `classify()` fonksiyonunda kodludur, yalnız IDP
adı ve BIG40 kategorisini kullanır, `--stage classify` ile `data/interim/imaging_gwas/idp_classification.tsv`
dosyasına yazıldı ve çıkarım aşaması bu dosyayı **okur**, yeniden türetmez.

| Sınıf | IDP | Kural |
|---|---|---|
| **ak madde mikroyapı** | **675** | BIG40 kategorileri: WM tract FA (75), WM tract diffusivity (300 = MD, L1, L2, L3), WM tract MO (75), WM tract OD (75), WM tract ICVF (75), WM tract ISOVF (75). Traktlar ve FA iskeleti üzerinde difüzyon MRG |
| **ak madde makroyapı** | **8** | Ad listesiyle: `IDP_T2_FLAIR_BIANCA_WMH_volume`, `aseg_global_volume_WM-hypointensities`, SIENAX ak madde hacmi (normalize / normalize değil), `aseg_lh/rh_volume_CerebralWhiteMatter`, `aseg_lh/rh_volume_Cerebellum-White-Matter` |
| **gri madde ve diğer (negatif kontrol sınıfı)** | **3.236** | rfMRI bağlantısallığı (1.701), bölgesel ve doku hacmi (647 eksi 7 ak madde hacmi = 640), kortikal alan (372), kortikal kalınlık (306), rfMRI düğüm genliği (76), **kortikal gri-ak kontrastı (70)**, bölgesel ve doku yoğunluğu (41), tfMRI aktivasyonu (16), bölgesel T2* (14) |
| QC | 16 | Görüntüleme kalite ölçütleri, fenotip değil; ayrı tutuldu |

**Tartışmalı tek atama, ve önemi büyük: kortikal gri-ak kontrastı (70 IDP).** Literatürde bu ölçüt
kısmen korteks içi miyelinasyonun göstergesi sayılır, yani bir eleştirmen onu miyelin fenotipleriyle
birlikte sınıflamak isteyebilir. Analiz planının kendi tanımları kesin: ak madde mikroyapısı difüzyon
MRG'dir, ak madde makroyapısı hacim ve lezyon yüküdür. Gri-ak kontrastı ikisi de değildir, bu yüzden
kontrol sınıfında kaldı. **Bu atama, lokusun 5e-8'i geçen tek hitini kontrol sınıfına koyuyor**;
taşınırsa ak madde tarafı "5e-8'i geçen bir hit var" diyebilir. İki okuma da aşağıda ayrı ayrı
verildi ki okuyucu kendi kararını verebilsin. Kararın çıkarımdan önce yazılmış olması dosya
zaman damgasıyla kayıtlı, ama dürüst olmak gerekirse sınıflandırmayı sonlandırdığım sırada gen
sayfasındaki en iyi hitin bu kategoride olduğunu **görmüştüm**; bu yüzden sonucu tek bir okumaya
dayandırmıyorum.

## 4. Eşikler (varsayılmadı, yayından canlı okundu)

`--stage thresholds` Smith 2021'in PMC tam metnini çekip eşik cümlelerini ayrıştırıyor:

| Eşik | -log10 p | p | Kaynak |
|---|---|---|---|
| Geleneksel genom çapı | 7,30 | **5e-08** | alan standardı |
| BIG40'ın tek fenotip eşiği | **7,5** | **3,16e-08** | Smith 2021: "standard GWAS threshold of -Log10(P) >= 7.5"; Manhattan grafikleri ve Table1 bununla filtrelenmiş |
| **BIG40'ın IDP'ler arası çoklu test eşiği** | **11,1** | **7,94e-12** | Smith 2021: "Bonferroni corrected level of -Log10(P) >= 11.1", tanım olarak "7.5 divided by the number of IDPs" |
| Bağımsız yeniden hesap | 11,095 | 8,04e-12 | 7,5 + log10(3.935); yayınla birebir tutuyor |

## 5. Çıkarım: KIF13A lokusu

Koordinatlar Ensembl **GRCh37** REST'ten canlı: KIF13A chr6:17.759.414-17.987.854 (eksi iplik),
NUP153 chr6:17.615.269-17.706.656, NHLRC1 chr6:18.120.718-18.122.851, TPMT chr6:18.128.542-18.155.305.
Pencereler gen gövdesi ±100 kb (chr6:17.659.414-18.087.854) ve ±250 kb (chr6:17.509.414-18.237.854).
**±250 kb penceresi TPMT'yi ve NHLRC1'i içine alır, ±100 kb almaz;** bölüm 5.3 bunun neden önemli
olduğunu gösteriyor.

### 5.1 Sınıf başına minimumlar (3.935 IDP'nin tamamı tarandı)

| Pencere | Sınıf | IDP | p ≤ 0,01 olan IDP | en iyi p | en iyi IDP | öncü varyant | 5e-8 geçen | 3,16e-8 geçen | 7,94e-12 geçen |
|---|---|---|---|---|---|---|---|---|---|
| ±100 kb | **ak madde mikroyapı** | 675 | 631 | **3,5e-07** | L3, sol ön korona radiata (UKB 25319) | rs144963738 chr6:18.012.583 | **0** | **0** | **0** |
| ±100 kb | **ak madde makroyapı** | 8 | 8 | 2,2e-04 | SIENAX ak madde hacmi (UKB 25007) | rs114157579 chr6:17.872.329 | 0 | 0 | 0 |
| ±100 kb | gri madde ve diğer | 3.236 | 3.017 | **1,1e-08** | gri-ak kontrastı, sağ medial orbitofrontal (UKB 27037) | rs143422697 chr6:17.828.055 | **1** | **1** | **0** |
| ±100 kb | QC | 16 | 14 | 4,1e-04 | T1 doğrusal olmayan hizalama bükülmesi | rs71556153 | 0 | 0 | 0 |
| ±250 kb | ak madde mikroyapı | 675 | 669 | 3,5e-07 | aynı | aynı | 0 | 0 | 0 |
| ±250 kb | ak madde makroyapı | 8 | 8 | 2,2e-04 | aynı | aynı | 0 | 0 | 0 |
| ±250 kb | gri madde ve diğer | 3.236 | 3.186 | 1,1e-08 | aynı | aynı | 1 | 1 | 0 |
| ±250 kb | QC | 16 | 16 | 5,6e-05 | rfMRI ters tSNR | rs6911008 | 0 | 0 | 0 |

Pencere genişletmek hiçbir sınıfın minimumunu değiştirmiyor: her iki sınıfın öncüsü de ±100 kb
içinde. Tam tablo `data/interim/imaging_gwas/summary_by_class.tsv`, sınıf başına ilk altı hit
`top_hits.tsv`.

### 5.2 Üst hitler

**(A) Tek 5e-8 geçen hit, kontrol sınıfında: rs143422697.** chr6:17.828.055 (GRCh37) /
17.827.824 (GRCh38), T>G, KIF13A intronu. BIG40 af 0,018; gnomAD v4 NFE 0,0140. Fenotip UKB 27037,
sağ medial orbitofrontal gri-ak yoğunluk kontrastı, N = 31.968, h² = 0,268.
**β = +0,18113, SE = 0,031696, p = 1,11e-08** (tam özet istatistikten, PheWeb ile birebir).
-log10 p = 7,96: 5e-8'i ve BIG40'ın 7,5 eşiğini geçiyor, **BIG40'ın kendi 11,1 eşiğinin 3,1 büyüklük
mertebesi altında kalıyor.** Yalnız değil: aynı varyantta sağ frontal kutup 1,9e-05, sol frontal
kutup 2,6e-05, sol rostral orta frontal 2,6e-04, sağ pars orbitalis 4,3e-04, hepsi aynı yönde. Yani
bir frontal gri-ak kontrast kümesi, ama bunlar birbiriyle korelasyonlu ölçütler, tek bir sinyal
sayılmalı.

**(B) Ak madde mikroyapı sınıfının öncüsü: rs144963738.** chr6:18.012.583, G>C, af 0,037, KIF13A'nın
5' ucundan 24,7 kb ötede. En iyi 8 fenotipi:

| IDP | UKB | sınıf | β | SE | p |
|---|---|---|---|---|---|
| L3 sol ön korona radiata | 25319 | ak madde mikroyapı | -0,11 | 0,021 | **3,5e-07** |
| L3 sağ ön korona radiata | 25320 | ak madde mikroyapı | -0,10 | 0,021 | 1,5e-06 |
| ICVF sol ön korona radiata | 25367 | ak madde mikroyapı | +0,10 | 0,021 | 2,0e-06 |
| ICVF sağ ön korona radiata | 25366 | ak madde mikroyapı | +0,098 | 0,021 | 3,3e-06 |
| FA sol ön korona radiata | - | ak madde mikroyapı | +0,094 | 0,021 | 7,4e-06 |
| MD sol ön korona radiata | - | ak madde mikroyapı | -0,092 | 0,021 | 1,3e-05 |
| MD korpus kallozum genusu | - | ak madde mikroyapı | -0,089 | 0,021 | 2,5e-05 |
| L3 korpus kallozum genusu | - | ak madde mikroyapı | -0,088 | 0,021 | 3,0e-05 |

Desen iç tutarlı ve miyelin yönünde: radyal difüzivite (L3) ve ortalama difüzivite aşağı, FA ve ICVF
yukarı, iki yarımkürede simetrik, tek bir anatomik yapıda (ön korona radiata ve korpus kallozum
genusu) toplanmış. p ≤ 0,01 olan 205 IDP'nin **102'si** ak madde mikroyapısı (beklenen 35,2); ilk
20'nin 14'ü. Aynı haplotipteki en güçlü 8 varyant (chr6:18.008.324-18.020.433) aynı fenotipi
gösteriyor, p 3,5e-07 ile 6,9e-07 arasında.

**(C) Ak madde makroyapı: hiçbir şey.** 8 IDP'nin en iyisi 2,2e-04 (SIENAX ak madde hacmi). WMH
hacminde (UKB 25781, N = 32.114) pencerede kayda değer hiçbir şey yok. Pozitif kontroller aynı
fenotipte 1,1e-29 ve 4,9e-47 verdiği için bu, gücün değil sinyalin yokluğudur.

### 5.3 BIG40'ın kendi tepe tablosu ve pencere seçiminin tuzağı

`release2/Table1.txt`, 10.889 yerel tepe (-log10 p > 7,5, 22k keşif, 11k tekrar ile birlikte):
**chr6:17.509.414-18.237.854 penceresinde sıfır kayıt.** En yakın kayıtlar 16,42 Mb (2), 17,21 Mb,
18,98 Mb ve 19,38 Mb'de. Yani BIG40'ın kendi anlamlılık ölçütüne göre bu lokusta hiçbir IDP için bir
tepe yok. (33k birleşik analizde 1,1e-08'in çıkması çelişki değil: Table1 22k keşif örneklemine
dayanır.)

**Pencere tuzağı.** ±250 kb penceresinde, p ≤ 0,01 ile en az 60 IDP'ye dokunan varyantlar arasında ak
madde mikroyapı **payı** en yüksek olanlar KIF13A haplotipi değil, chr6:18.142.205-18.165.244
aralığındaki **yaygın** varyantlar (af 0,20-0,39, pay 0,83-0,87, en iyi ak madde p = 7,5e-07,
rs9297056, sağ arka korona radiatada L1). Bu aralık **NHLRC1 ve TPMT** bölgesidir, KIF13A'nın 5'
ucundan 155-177 kb ötede. LD kontrolü bunları kesin olarak ayırıyor: rs9297056'nın TPMT korteks eQTL
öncüsü rs4716224 ile r² = 0,672, NHLRC1 korteks eQTL öncüsü rs3930695 ile r² = 0,265, ama KIF13A
haplotipinin hiçbir üyesiyle **r² < 0,05**. Dolayısıyla ±250 kb penceresinde "ak madde zenginleşmesi"
diyen biri KIF13A'yı değil TPMT/NHLRC1'i görüyor olabilir. **KIF13A'ya atıf için doğru pencere
±100 kb'dir.**

## 6. Kontroller

### 6.1 Pozitif kontroller: boru hattı sinyal bulabiliyor mu

Literatürden canlı seçildi. Zhao B ve ark. 2021, *Science*, **PMID 34140357**,
**doi 10.1126/science.abf3736** ("Common genetic variation influencing human white matter
microstructure", 43.802 kişi, 109 lokus) üç adayın hepsini ak madde bağlamında anıyor; **VCAN** ayrıca
Elliott LT ve ark. 2018 *Nature*, **PMID 30305740** ile birlikte UKB difüzyon MRG bağlamında geçiyor;
**TRIM47/TRIM65** ve **EFEMP1** ak madde hiperintensitesi ve serebral küçük damar hastalığı
literatürünün kanonik lokusları (Europe PMC tam metin sorgularıyla doğrulandı). Aynı çıkarım, aynı
683 ak madde IDP'si, aynı ±250 kb:

| Lokus | ak madde makroyapı en iyi p | ak madde mikroyapı en iyi p | 5e-8 geçen | 7,94e-12 geçen | Öncü fenotip |
|---|---|---|---|---|---|
| **VCAN** (5q14.3) | 2,1e-18 (WMH hacmi) | **6,8e-106** | **434** / 675 | **369** | ICVF, alt longitudinal fasikulus |
| **TRIM47** (17q25.1) | **4,9e-47** (WMH hacmi) | 1,6e-14 | 48 + 4 | 3 + 4 | WMH hacmi, ICVF tapetum |
| **EFEMP1** (2p16.1) | **1,1e-29** (WMH hacmi) | 1,6e-14 | 20 + 2 | 4 + 2 | WMH hacmi, ISOVF arka talamik radyasyon |
| KIF13A | 2,2e-04 | 3,5e-07 | **0** | **0** | - |

**Üç pozitif kontrolün üçü de geçti**, hem de BIG40'ın en katı eşiğini fazlasıyla. Boru hattı
çalışıyor: ak madde sinyali varsa 1e-14 ile 1e-106 arasında görünüyor. KIF13A'nın en iyi ak madde
p değeri bu aralıktan 7 ile 99 büyüklük mertebesi uzakta.

### 6.2 Negatif kontroller ve nöronal karşılaştırma

| Lokus | rol | ak madde mikroyapı en iyi p | ak madde makroyapı | diğer sınıf en iyi p | herhangi bir eşiği geçen |
|---|---|---|---|---|---|
| **KIF1A** (chr2) | nöronal kinezin | 6,3e-06 (OD forniks/stria terminalis) | 2,6e-04 | 3,9e-08 (rfMRI bağlantısallığı) | 5e-8'de 1, BIG40 7,5'te **0** |
| **LCT** (2q21) | negatif kontrol | **7,6e-07** (OD forniks/stria terminalis) | 3,3e-03 | - | 0 |
| **MYH7** (14q11.2) | negatif kontrol | bkz. bölüm 6.3 | | | |

**KIF1A karşılaştırması bilerek kurulmuş bir sınavdı ve öngörüyü nicel olarak destekliyor ama zayıf
biçimde:** nöronal kinezinin ak madde mikroyapı minimumu 6,3e-06, KIF13A'nın 3,5e-07, yani 18 kat
daha iyi. Ayrıca KIF1A'nın en iyi ak madde varyantlarının fenotipleri birbirinden kopuk (forniks,
üst fronto-oksipital, kortikospinal, üst korona radiata) ve p ≤ 0,01 olan IDP'lerdeki ak madde payı
0,23-0,31, KIF13A haplotipinde 0,49-0,50. Ama iki lokus türsel olarak aynı yerde: her ikisinde de
"diğer" sınıfında 5e-8'i geçen **bir** hit var (1,1e-08 ve 3,9e-08), ikisi de BIG40'ın 7,5 eşiğini
geçemiyor ya da sınırda kalıyor, ve **hiçbirinde hiçbir sınıfta 11,1 eşiğini geçen hiçbir şey yok.**

**Asıl soğuk duş LCT'dir.** Beyin fenotipi beklenmeyen bir lokusta, aynı pencere boyutunda, ak madde
mikroyapı minimumu 7,6e-07 çıkıyor. KIF13A'nın 3,5e-07'si bunun yalnız 2,2 katı. Yani **minimum p
istatistiği tek başına KIF13A ile rastgele bir lokusu ayırmıyor.** (LCT'nin "diğer" sınıfı taranmadı,
ama gen sayfası orada 2,3e-13'lük bir rfMRI bağlantısallığı hiti gösteriyor; 3.935 korelasyonlu
fenotipte, beyin fenotipi beklenmeyen bir lokusta bile, 1e-13 düzeyinde bir hit çıkabiliyor. Bu,
KIF13A'nın 1,1e-08'ini ayrıca bağlamına oturtuyor.)

### 6.3 Rastgele eşleşmiş pencereler: ampirik boş dağılım

RASTGELE_NULL_PLACEHOLDER

## 7. Lokus mu gen mi: bu boşluğu kapatma denemesi

### 7.1 İki ayrı haplotip

LD, Ensembl GRCh37 REST, 1000 Genomes faz 3 EUR:

| Öncü | Konum (GRCh37) | BIG40 af | gnomAD v4 NFE | r² ≥ 0,2 ortak | LD aralığı | KIF13A gövdesindeki ortak sayısı |
|---|---|---|---|---|---|---|
| **rs143422697** | 17.828.055 | 0,018 | 0,0140 | 6 | 17.827.542-18.074.064 (247 kb) | 3 / 6 |
| **rs75428660** | 17.886.045 | 0,043 | 0,0363 | 185 | 17.717.338-18.099.929 (383 kb) | 122 / 184 |
| **rs76091407** | 18.020.433 | 0,036 | 0,0313 | 177 | 17.776.213-18.226.776 (451 kb) | 107 / 176 |

rs143422697 ile rs75428660 arasında **r² < 0,05**: lokusta **iki bağımsız haplotip** var.
rs75428660 ile rs76091407 r² = 0,667, rs144963738 ile rs75428660 r² = 0,969. Yani (B) ve (C)
bölümlerindeki ak madde sinyali ile yayımlanmış kortikal morfoloji sinyali **aynı** haplotiptir.

**Komşular:** üç haplotipin hiçbiri NUP153 gövdesine bir varyantla girmiyor (NUP153,
rs75428660 bloğunun başlangıcından 10,7 kb önce bitiyor) ve hiçbirinin LD kümesinde NUP153'ün
korteks eQTL öncüsü rs77458411 yok (r² < 0,05). **Projenin "dürüst varsayılan açıklaması" olan
NUP153 bu lokusta elenir.** Buna karşılık rs76091407 bloğu konum olarak NHLRC1'i (2,1 kb'lik küçük
bir gen) kapsıyor ve TPMT'ye kadar uzanıyor, dolayısıyla kredibl küme KIF13A ile sınırlı **değil**;
ama NHLRC1 ve TPMT'nin eQTL öncü varyantlarının hiçbiri bu haplotipte değil (aşağıda).

### 7.2 Beyin eQTL'i: sinyal KIF13A ifadesine dokunuyor mu

**MetaBrain cis-eQTL tarayıcısı** (`metabrain.nl/cis-eqtls.html`, arka uç
`eqtlgen.molgeniscloud.org/api/v2/metabrain_cis_eqtl_v2`, yayın sürümü 2023-07-18; GRCh38 konumları).
Erişilebilen tek beyin eQTL sorgu arayüzü buydu: eQTL Catalogue REST API'si bugün **HTTP 410**
("This API is no longer available") veriyor, PsychENCODE yalnız kütlesel indirme sunuyor, Brain xQTL
Serve bağlantı zaman aşımına düştü.

| Gen | Beyin bölgesi | Sıra | Öncü | p | β | N |
|---|---|---|---|---|---|---|
| **KIF13A** | Korteks EUR | 1 | rs594115 (17.860.176) | **2,92e-46** | -0,486 | 2.682 |
| **KIF13A** | Korteks EUR | 2 | rs1924466 (17.772.533) | 9,16e-12 | -0,186 | 2.681 |
| **KIF13A** | Korteks EUR | **3** | **rs115824472** (18.009.434) | **5,14e-06** (q = 8,8e-03) | -0,388 | 2.646 |
| KIF13A | Serebellum EUR | 1 | rs6910353 (18.102.327) | 4,29e-10 | -0,384 | 491 |
| NUP153 | Korteks EUR | 1 | rs77458411 (17.657.491) | 1,47e-05 | -0,183 | 2.683 |
| NHLRC1 | Korteks EUR | 1 | rs3930695 (18.149.437) | 1,65e-12 | -0,219 | 2.436 |
| TPMT | Korteks EUR | 1 | rs113340706 (18.163.848) | 2,71e-59 | -1,058 | 2.276 |

**Kritik bağlantı:** rs115824472, KIF13A'nın **üçüncü bağımsız korteks eQTL sinyalinin öncüsüdür** ve
görüntüleme haplotipiyle **r² = 0,969** (rs75428660 ve rs144963738 ile). Buna karşılık aynı haplotip
NUP153, NHLRC1 ve TPMT'nin **hiçbir** eQTL öncüsüyle LD içinde değil (hepsi r² < 0,05) ve KIF13A'nın
kendi iki güçlü eQTL sinyaliyle de (rs594115, rs1924466) LD içinde değil.

**GTEx v8** (API canlı): rs75428660 KIF13A için **yedi** dokuda anlamlı eQTL, biri beyin:
Brain_Putamen_basal_ganglia p = 2,6e-05 NES = -0,221; ayrıca tiroit 1,4e-07 (NES -0,365), karaciğer
8,1e-07 (NES -0,703), deri 1,7e-07 (NES **+0,338**), fibroblast 1,1e-05, tibial arter 2,4e-05. Ek
olarak deride KIF13A için bir **sQTL** (intron kümesi chr6:17.951.429-17.987.054, p = 2,2e-06) ve
**NHLRC1** için iki deri eQTL'i.
Dürüst okuma: (i) 13 GTEx beyin dokusunun **yalnız biri** anlamlı ve p = 2,6e-05, eQTL ölçeğinde
zayıf; (ii) etki **yönü dokular arasında değişiyor** (deride +, tiroit, karaciğer ve beyinde -), bu
kararsız bir eQTL'in işaretidir; (iii) aynı varyant NHLRC1'e de dokunuyor, yani GTEx tek başına
KIF13A'yı seçmiyor. rs143422697 ve rs4716216 GTEx'te `maf01 = false`, yani **hiç test edilmemiş**:
1,4-1,8 % frekanslı bir varyant için GTEx'in gücü yok, bu yüzden oradaki eQTL yokluğu bilgi
taşımıyor.

### 7.3 Open Targets kredibl kümeleri ve L2G

Open Targets Platform GraphQL (canlı, kısıtlılık sayıları kullanılmadı, proje kuralı):
rs75428660'ın **8 kredibl kümesi** var, ikisi GWAS: **GCST90095129 "Vertex-wise sulcal depth"
p = 1e-09** ve **GCST90095131 "Vertex-wise cortical thickness" p = 5e-09**, N = 33.748, ikisi de
**L2G = 0,622 ile KIF13A'yı** gösteriyor. Diğer altısı eQTL/pQTL: GTEx kas ve fibroblast, TwinsUK
yağ dokusu (×2), HipSci iPSC (ENSG00000137177 = KIF13A) ve plazma SCARA5 pQTL'i
(UKB_PPP_EUR, p = 1,3e-13; SCARA5 chr8'de olduğu için bu bir *trans* etkidir). **Bu kredibl
kümelerin hiçbiri beyin dokusundan değil.**

> **Proje kaydına düzeltme.** `docs/vaka_arama.md` 7e bölümü "KIF13A'ya atanmış tek bir GWAS
> credible set'i yok" diyor. 17 Eyl 2026'da canlı sorguda KIF13A'ya **L2G 0,622 ile iki GWAS kredibl
> kümesi** atanmış durumda (sulkus derinliği ve kortikal kalınlık). Eski ifade gen sayfasından gelen
> `studyType` dağılımına dayanıyordu; varyant tarafından sorgulandığında tablo değişiyor. 7e
> güncellenmelidir.

rs143422697 Open Targets'ın varyant dizininde **hiç yok** (çok nadir), dolayısıyla BIG40'ın 5e-8
geçen tek hiti için L2G kanıtı yok.

## 8. Peşine düşülen iki iz ve ortaya çıkan üçüncüsü

### 8.1 Singulum orientasyon dağılımı: öngörü çürütüldü

`docs/vaka_arama.md` 7c, AstraZeneca PheWAS'ın en iyi ikinci KIF13A hitini kaydediyor:
p = 2,654e-05, *"Mean OD in cingulum cingulate gyrus on FA skeleton (left)"*, **n = 44**, beş
collapsing modelinde tekrarlıyor ve orada gürültü sayılmış. Bu fenotip tam olarak
oligodendrosit hipotezinin öngördüğü sınıftadır, bu yüzden en spesifik tek sınavdı.

BIG40'ta aynı fenotip **`IDP_dMRI_TBSS_OD_Cingulum_cingulate_gyrus_L`, UKB alanı 25427, pheno 2012,
N = 31.355, h² = 0,265**; tanımı birebir aynı ("Mean OD (orientation dispersion index) in cingulum
cingulate gyrus (left) on FA skeleton"). Tam özet istatistiği indirildi (MAF ≥ 0,001, tam kapsam):

| Pencere | varyant | en iyi p | öncü |
|---|---|---|---|
| ±100 kb | 2.879 | **1,43e-04** | rs186536699 chr6:17.873.162, β = +0,480, SE = 0,126 (çok nadir) |
| ±250 kb | 4.945 | 2,09e-05 | rs13194542 chr6:17.551.996, β = +0,0353, SE = 0,0083 (gen gövdesinin dışında, 3' ucundan 207 kb ötede, NUP153'ün de ötesinde) |

2.879 varyantta beklenen en küçük p bağımsızlık varsayımıyla 3,5e-04, LD nedeniyle etkin test sayısı
daha az; gözlenen 1,43e-04 tam olarak gürültü beklentisidir. **713 kat daha büyük bir örneklemde
aynı fenotip, aynı lokus: hiçbir şey.** AZ PheWAS izi kapandı, doğru karar zaten gürültü demekti.

### 8.2 rs75428660: doğrulandı, KIF13A'yı gösteriyor, ama fenotipi ak madde değil

`vaka_arama.md` 7f'nin "tek genom çapında anlamlı nörolojik hikaye" dediği varyant. Canlı
doğrulandı (GWAS Catalog REST + Ensembl iki yapıda):

| | |
|---|---|
| Konum | **GRCh38 chr6:17.885.814, GRCh37 chr6:17.886.045**, A>G, KIF13A intronu, uzaklık 0 |
| Frekans | BIG40 af 0,043; gnomAD v4 NFE 0,0363; 1000G küresel 0,0092 |
| **Sulkus derinliği** | p = **1e-09**, β = 6,09 z, **GCST90095129**, PMID **34910505** ("The genetic architecture of human cortical folding"), N = 33.748 Avrupa kökenli |
| **Kortikal kalınlık** | p = **5e-09**, β = 5,85 z, **GCST90095131**, aynı yayın |
| **Kortikal kalınlık (bağımsız)** | p = **2e-08**, **GCST90091061**, PMID **34560273** ("Vertex-wise multivariate genome-wide association study identifies 780 unique genetic loci"), N = 35.657 |
| Geni işaret ediyor mu | **Evet, KIF13A'yı:** L2G 0,622 (iki GWAS kredibl kümesi), MetaBrain korteks KIF13A eQTL sıra 3 ile r² = 0,969, GTEx KIF13A eQTL 7 doku (biri putamen) + sQTL |
| Komşuyu işaret ediyor mu | NUP153 hayır (LD dışı, eQTL öncüsü LD dışı). NHLRC1 kısmen: GTEx'te aynı varyant iki deri dokusunda NHLRC1 eQTL'i, ama NHLRC1'in beyin eQTL öncüleri LD dışı |
| BIG40'taki profili | En iyi fenotipleri **ak madde mikroyapısı**: L3 sol/sağ ön korona radiata 1,2e-05 ve 1,4e-05, FA sol ön korona radiata 5,1e-05, ICVF sağ/sol ön korona radiata 7,7e-05 ve 7,8e-05, L3 korpus kallozum genusu 1,3e-04. p ≤ 0,01 olan 139 IDP'nin **51'i** ak madde mikroyapısı (beklenen 23,8) |

**Hüküm:** iz gerçek ve KIF13A'yı gösteriyor, ama yayımlanmış ve iki bağımsız çalışmada genom çapında
anlamlı olan fenotipi **kortikal morfoloji** (sulkus derinliği, kortikal kalınlık), ak madde değil.
BIG40'ta bu varyantın kendi ak madde ilişkileri 1e-05 düzeyinde, yani eşiğin çok altında. Bu iz
oligodendrosit hipotezinin lehine değil, "lokusun bir beyin fenotipi var ama beklediğimiz fenotip
değil" tarafına yazılır.

### 8.3 Aramadığım ama bulduğum iz: rs76091407, ak maddede genom çapında anlamlı

GWAS Catalog'un ±250 kb bölge sorgusu (GRCh38 chr6:17.508.198-18.238.274; 135 SNP, 196 ilişki)
penceredeki **tek** ak madde fenotipini ortaya çıkardı:

| | |
|---|---|
| Varyant | **rs76091407**, GRCh38 chr6:18.020.202 / GRCh37 chr6:18.020.433, G>A |
| Konum | GWAS Catalog'un hesabına göre en yakın gen **KIF13A, 32.567 bp** (Catalog gen sonunu 17.987.635 alıyor; Ensembl GRCh38 sınırıyla 31.928 bp, GRCh37'de 32.579 bp), NHLRC1 100.238 bp. KIF13A eksi iplik olduğu için bu mesafe genin **5' ucundan**, yani promotor tarafındandır |
| İlişki | **p = 4e-09**, EFO "white matter microstructure measurement" |
| Çalışma | **GCST90860937**, PMID **42236496**, Zhao X ve ark., *Nat Commun* 2026, doi **10.1038/s41467-026-73996-z**, "Genetic architecture of white matter microstructure captured by unsupervised deep representation learning of fractional anisotropy maps", **N = 11.736** İngiliz kökenli, fenotip **UDIP-FA** (gözetimsiz derin görüntüleme fenotipi) |
| BIG40'taki profili | En iyi fenotip L3 sol ön korona radiata **6,9e-07**; p ≤ 0,01 olan 201 IDP'nin **100'ü** ak madde mikroyapısı (beklenen 34,5), ilk 20'nin 15'i |
| Haplotip | rs75428660 ile r² = 0,667, rs144963738 ile r² = 0,689, rs115824472 (KIF13A korteks eQTL sıra 3) ile r² = 0,689; rs143422697 ile r² < 0,05 |

**Bu, bu eksenin ürettiği en olumlu tek bulgudur ve hemen üç çekince ister.** (i) Fenotip gözetimsiz
derin öğrenmeyle türetilmiş gizli bir bileşendir, doğrudan yorumlanabilir bir trakt ölçütü değil.
(ii) Çalışma **N = 11.736 ile 586 lokus ve 939 öncü SNP** bildiriyor; karşılaştırma için Zhao B 2021
*Science* (PMID 34140357) **N = 43.802 ile 109 lokus** bildiriyor. Dört kat küçük bir örneklemde beş
kat fazla lokus, lokus başına önsel güvenilirliği düşürür ve bu hit o 939'dan biridir. (iii) Öncü
varyant KIF13A gövdesinin dışındadır ve LD bloğu konum olarak NHLRC1'i kapsayıp TPMT'ye uzanır.
Yine de, bu haplotipin BIG40'taki alt eşik deseniyle birebir aynı fenotipi göstermesi ve KIF13A
korteks eQTL'i ile r² = 0,69 olması rastlantı olarak açıklanması kolay bir şey değil.

## 9. Güç: bu negatif sonuç neyi dışlıyor

Fenotipler birim varyansa ölçekli ve kodlama toplamsal olduğu için ki-kare ≈ N × R²; dolayısıyla %80
güçle saptanabilir **açıklanan varyans** frekanstan bağımsızdır, frekans yalnız bunun hangi β'ya
karşılık geldiğini belirler. `--stage report` bunu gözlenen standart hatalarla da çapraz kontrol
ediyor (dMRI IDP'lerinde SE ortancası 0,0279, en bilgili %5'te 0,0082).

| N | Eşik | Saptanabilir R² | β (MAF 0,05) | β (MAF 0,25) |
|---|---|---|---|---|
| 31.355 (dMRI) | 5e-08 | **%0,126** | 0,115 SD | 0,058 SD |
| 31.355 (dMRI) | BIG40 7,5 (3,16e-08) | %0,130 | 0,117 SD | 0,059 SD |
| 31.355 (dMRI) | BIG40 Bonferroni 11,1 (7,94e-12) | %0,188 | 0,141 SD | 0,071 SD |
| 33.224 (T1 hacim) | 5e-08 | %0,119 | 0,112 SD | 0,056 SD |
| 33.224 (T1 hacim) | BIG40 Bonferroni 11,1 | %0,178 | 0,137 SD | 0,069 SD |

**Doğru ifade budur:** BIG40'ın ak madde negatifi, tek bir ak madde IDP'sinin varyansının
**%0,126'sından büyük** bir etkiyi 5e-8'de dışlar (ve %0,188'den büyüğünü BIG40'ın kendi eşiğinde).
Daha küçük etkiler dışlanmaz. Haplotipte **gözlenen** etki, rs144963738'in β = -0,11 ve af = 0,037
değerleriyle R² = 2 × 0,037 × 0,963 × 0,11² = **%0,086**, yani saptanabilir sınırın altında. Sonuç:
gözlenen alt eşik p değeri, bu büyüklükte **gerçek** bir etkinin bu örneklemde üreteceği şeyin tam
olarak kendisidir. Aynı etkiyi 5e-8'de %80 güçle görmek için gereken örneklem
N = 39,6 / 0,00086 ≈ **46.000** kişidir; UK Biobank görüntüleme kohortu bu büyüklüğe ulaştığı için
bu soru gelecekteki bir BIG40 sürümüyle kesin olarak kapatılabilir.

Ak madde **makroyapı** tarafında güç ifadesi daha keskin: WMH hacminde (N = 32.114) pozitif
kontroller 1,1e-29 ve 4,9e-47 verirken KIF13A penceresinde en iyi değer 2,2e-04. TRIM47 ve EFEMP1
ölçeğinde bir WMH etkisi olsaydı 20 ile 40 büyüklük mertebesi arayla görülürdü. **Ak madde lezyon
yükü tarafı gerçekten boştur ve bu boşluk güç eksikliğiyle açıklanamaz.**

## 10. Bu eksen ne kadar ileri götürülebilir

**Lisanslı en güçlü cümle:**

> UK Biobank beyin görüntüleme genetiğinde KIF13A lokusunda ak madde mikroyapısı yönünde tutarlı ama
> eşik altı bir iz vardır: 33 bin kişilik BIG40 taramasında lokusun 3,6 % frekanslı haplotipi, iki
> yarımkürede ön korona radiata ve korpus kallozum genusunda radyal difüzivite ve ortalama
> difüzivitede azalma, FA ve ICVF'te artış gösteriyor (en iyi p = 3,5e-07, hiçbir eşiği geçmiyor);
> aynı haplotipteki bir varyant ak maddeye özel bağımsız bir GWAS'ta genom çapında anlamlıdır
> (p = 4e-09, N = 11.736) ve haplotip, insan korteksindeki KIF13A cis-eQTL'inin üçüncü bağımsız
> sinyaliyle r² = 0,97 içindedir. BIG40'ın bu etkiyi doğrulamaya gücü yoktur: %80 güçle
> saptanabilir alt sınır varyansın %0,126'sı, gözlenen etki %0,086.

**Kaçınılması gereken aşırı iddia** (ve neden):

1. **"UK Biobank görüntülemesi KIF13A'nın bir ak madde geni olduğunu doğruluyor."** Yanlış. Ak madde
   sınıflarında 5e-8'i geçen hiçbir şey yok, BIG40'ın kendi tepe tablosunda pencerede sıfır kayıt
   var, ve negatif kontrol lokusunun (LCT) ak madde minimumu yalnız 2,2 kat daha kötü.
2. **"Lokus genom çapında anlamlı bir hit taşıyor, demek ki KIF13A beyin fenotipiyle ilişkili."**
   Yarı doğru ve yanlış yönde: 5e-8'i geçen tek hit gri-ak yoğunluk kontrastıdır (1,1e-08), ayrı bir
   haplotiptedir, Open Targets'ta hiç yoktur ve BIG40'ın kendi 11,1 eşiğinin 3,1 büyüklük mertebesi
   altındadır. Ayrıca aynı veri tabanında, beyin fenotipi beklenmeyen bir lokusta 2,3e-13 çıkıyor.
3. **"±250 kb penceresi ak madde fenotipleri için zenginleşmiş."** Bu doğru ama nedeni KIF13A değil:
   en zengin varyantlar (pay 0,87) TPMT/NHLRC1 bölgesinin yaygın haplotipidir ve KIF13A haplotipiyle
   r² < 0,05'tir. ±100 kb kullanılmalıdır.
4. **"Singulum orientasyon dağılımı hipotezi doğrulandı."** Çürütüldü: 713 kat daha büyük örneklemde
   aynı fenotipte ±100 kb'de en iyi p = 1,4e-04.
5. **"KIF13A ifadesinin beyin eQTL'i sinyali açıklıyor."** Fazla ileri: beyin eQTL bağlantısı tek bir
   MetaBrain korteks sinyalinin üçüncü sırası (p = 5,1e-06, q = 8,8e-03) ve tek bir GTEx beyin dokusu
   (putamen, p = 2,6e-05, dokular arası yön tutarsız). Kolokalizasyon yapılmadı; yapılması gereken
   iş budur.
6. **"Fenotip aktarılamaz kuralı çürütüldü/doğrulandı."** Ne biri ne öteki: KIF1A ve KIF13A BIG40'ta
   türsel olarak aynı yerde duruyor (her ikisinde de "diğer" sınıfında tek bir 5e-8 hiti, hiçbir
   sınıfta 11,1 eşiğini geçen hiçbir şey yok). Ayrımı yapan tek şey KIF13A'nın ak madde
   minimumunun 18 kat daha iyi ve fenotipçe tutarlı olmasıdır, ki bu tek başına zayıf bir ayrımdır.

**Sıradaki en verimli üç iş:** (i) UDIP-FA özet istatistiklerini (GCST90860937) indirip KIF13A
lokusundaki sinyalin BIG40'ın ön korona radiata fenotipiyle **kolokalize** olup olmadığını test etmek;
(ii) rs115824472/rs75428660 haplotipi ile MetaBrain KIF13A korteks eQTL'i arasında resmi
kolokalizasyon (coloc/SMR); (iii) BIG40'ın bir sonraki sürümünü ya da 45-50 bin kişilik başka bir UKB
görüntüleme GWAS'ını beklemek, çünkü gereken örneklem hesaplanmış durumda (~46.000).

## 11. Yeniden üretim

Sistem `python` 3.13 ile, depo kökünden (proje venv'i Windows tarafından bloklu, CLAUDE.md):

```powershell
python -u scripts\imaging_gwas_axis.py --stage sources
python -u scripts\imaging_gwas_axis.py --stage probe
python -u scripts\imaging_gwas_axis.py --stage classify
python -u scripts\imaging_gwas_axis.py --stage thresholds
python -u scripts\imaging_gwas_axis.py --stage scan
python -u scripts\imaging_gwas_axis.py --stage sumstats --phenos 2012,1415,1479 --chrom 6 --lo 17509414 --hi 18237854
python -u scripts\imaging_gwas_axis.py --stage followup
python -u scripts\imaging_gwas_axis.py --stage report
python -u scripts\imaging_gwas_axis.py --stage figure
```

Çıktılar: `data/interim/imaging_gwas/` (sınıflandırma, lokus başına tarama `.tsv.gz`, sınıf
özetleri, üst hitler, varyant sınıf profilleri, güç tablosu, eşikler), `data/raw/imaging_gwas/`
(kaynak sayfalar, IDP tablosu, Table1, Smith 2021 tam metni, Ensembl koordinatları, varyant
PheWAS'ları), figür `figures/KIF13A_ak_madde_goruntuleme.svg`.

### Bu eksende öğrenilen tuzaklar

1. **BIG40 GRCh37'dir, GWAS Catalog REST GRCh38 döndürür** ve bu lokusta iki yapı arasındaki fark
   yalnız **231 bp**'dir. rs75428660'ın konumu her iki yapıda da KIF13A gövdesinin içine düşer, yani
   yanlış yapı sessizce doğru görünür. Her varyant iki yapıda ayrı ayrı çözüldü.
2. **BIG40'ın PheWeb kopyası yalnız p ≤ 0,01 ve MAF ≥ 0,01 tutar.** Bu bir kusur değil avantajdır
   (minimum p tam geri kazanılır) ama bilinmezse "bu varyant yok" diye okunur. MAF 0,001-0,01
   aralığı yalnız tam özet istatistikte vardır.
3. **BIG40 PheWeb'in bölge uç noktası bir IDP için kalıcı HTTP 500 veriyor:**
   `IDP_dMRI_TBSS_FA_Posterior_corona_radiata_L` (pheno 1479, UKB 25083), ki bu bir **ak madde**
   fenotipidir. Sessizce atlanırsa ak madde taraması delik kalır. Script bu adı
   `scan_<lokus>_failed.txt` dosyasına yazıyor ve rapor aşaması onu tam özet istatistikten
   dolduruyor.
4. **Ensembl'in `MAF` alanı 1000 Genomes *küresel* frekansıdır.** rs75428660 için 0,0092 döner ama
   gnomAD v4 NFE 0,0363 ve BIG40 (UKB) af 0,043'tür. 4,7 katlık bu fark "PheWeb yanlış varyantı
   döndürdü" gibi okunabilir; okunmamalı, popülasyona özel frekansa bakılmalı.
5. **BIG40 özet istatistik dosyaları tek parçalı gzip'tir** (BGZF senkron noktası yok, test edildi),
   yani HTTP Range ile rastgele erişim mümkün değil; chr6:18 Mb'ye ulaşmak için dosyanın ~%30'unu
   (100-120 MB) akıtmak gerekir. 3.935 IDP için bu yol kapalıdır, sorgu arayüzü kullanılmalıdır.
6. **BIG40'ın IDP tablosu ile PheWeb aynı IDP'yi farklı yazar** (boşluk, `+`, parantez sanitize
   edilir). `[^A-Za-z0-9]+` -> `_` normalizasyonu 3.933/3.935 eşleştirir; `ThalamNuclei_lh/rh_volume_MV(Re)`
   PheWeb'de `...MVRe` olur ve elle eşlenmelidir. Bu yapılmazsa 1.876 IDP "PheWeb'de yok" görünür.
7. **PheWeb'in `/api/variant/` yanıtı `variant` anahtarı taşımaz**, alanlar (`chrom`, `pos`, `ref`,
   `alt`, `rsids`, `nearest_genes`, `phenos`) kökte durur.
8. **Ensembl LD `pairwise` uç noktası bir varyant panelde yoksa da boş liste döndürür**, yani "boş"
   ile "r² düşük" ayrılamaz. Doğru yol pencere sorgusudur (`/ld/human/<rsid>/<pop>?window_size=500&r2=0.05`)
   ve hedef varyantın dönen kümede olup olmadığına bakmaktır.
9. **eQTL Catalogue REST API'si kapandı** (HTTP 410, "This API is no longer available"). Erişilebilen
   beyin eQTL sorgu arayüzü MetaBrain'in DataTables arka ucudur
   (`eqtlgen.molgeniscloud.org/api/v2/metabrain_cis_eqtl_v2`, `?q=GeneSymbol==KIF13A`).
10. **GTEx'te `maf01: false` "eQTL değil" anlamına gelmez, "test edilmedi" anlamına gelir.**
    rs143422697 için GTEx'in boş dönmesi bilgi taşımıyor.
11. **Pencere genişliği sonucu değiştirir.** ±250 kb TPMT'yi (chr6:18.128.542-18.155.305) kapsar ve
    orada KIF13A ile LD içinde olmayan, ak madde açısından daha zengin ayrı bir haplotip vardır.
    Gen atfı için ±100 kb kullanılmalıdır.
12. **BIG40'ın `Table1.txt`'si 22k keşif örneklemindendir**, 33k birleşik analizden değil; 33k'da
    7,5'i geçen bir hit Table1'de olmayabilir. Çelişki sanılmamalı.
