# Bölgesel kısıtlılık: KIF13A motor domaini yerel olarak ölçüldüğünde missense'i kısıtlıyor mu

Ölçüm tarihi: 16 Eylül 2026. Script: `scripts/regional_constraint.py`.
Veri: `data/interim/regional_constraint/` ve `data/raw/rmc/`.
Figür: `figures/KIF13A_bolgesel_kisitlilik.svg`.

Soru: gen düzeyindeki missense Z = 2,79 ve gözlenen/beklenen = 0,861 sayısı 1805 aminoasitlik bir
proteinin ortalamasıdır; motor domain bunun yalnız 352'si. Uzun ve toleranslı bir kuyruk, kısıtlı bir
motor domaini ortalamada yutabilir. Bu belge kısıtlılığı **yerel olarak** ölçüyor.

---

## 0. Hüküm

**Zayıf kısıtlı.** gnomAD v4.1 ekzomlarında (730 bin ekzom, MANE transkripti) KIF13A motor domaini
(5-352) için missense gözlenen/beklenen **0,814** (%90 Poisson aralığı 0,745-0,892; 333 gözlenen,
409,0 beklenen; tükenme p = 9,5x10⁻⁵). gnomAD'nin kendi ölçeğine çevrildiğinde **0,78**; buna karşılık
kuyruk (520-1805) 0,89 ve gen bütünü 0,86. Yani **evet, gen düzeyindeki ortalama bir sinyali
gizliyordu, ama gizlediği sinyal küçük**: motor domain kuyruktan yalnızca yaklaşık %12 daha tükenmiş
ve bu fark sınırda (motor/kuyruk oranı 0,879, p = 0,036). Ölçümün hipotezi kurtarmaya en çok
yaklaşan parçası bir kat daha aşağıda: motor domainin **yüzeyi hiç kısıtlı değil** (o/e 1,008),
gömülü çekirdeği 0,678, mikrotübül arayüzü 0,746, **ATP cebi ise 0,523** (0,362-0,776; 18 gözlenen,
34,4 beklenen; p = 0,0016) ve katman A adaylarının 23'ü tam o ATP cebinde oturuyor. Tavanı belirleyen
iki karşılaştırma ise sert: (i) aynı boru hattı KIF1A motorunda 0,450, KIF1A'nın ATP cebinde **0,140**
veriyor, yani KIF13A her yapısal sınıfta KIF1A'dan 2 ila 4 kat daha sığ; (ii) aynı ölçüm KIF13A'nın
CC1+FHA bölgesinde (381-519) **0,785** veriyor, yani motor domain KIF13A'nın en kısıtlı bölgesi bile
değil. Bunu bağımsız olarak dört yayımlanmış kaynak da söylüyor: gnomAD'nin bölgesel missense
kısıtlılık algoritması KIF13A'yı **arayıp hiçbir kırılma noktası bulamamış**, MTR ve CCR ise
KIF13A'nın tek tükenmiş bölgesi olarak motor domaini değil **353-519'u** gösteriyor.

Tek cümleyle: **motor domainde gerçek ama sığ bir seçilim izi var; ortalama onu gizliyordu ama
ortaya çıkan şey KIF1A'nınkinin dörtte biri kadar ve KIF13A'ya özgü bile değil.**

| Bölge | Etiket |
|---|---|
| motor domain bütünü (5-352) | **zayıf kısıtlı** (o/e 0,814, aralık üst sınırı 0,892) |
| motor ATP cebi (25 kalıntı) | **orta kısıtlı** (o/e 0,523, üst sınır 0,776) |
| motor yüzeyi (153 kalıntı) | **hayır, tolere ediyor** (o/e 1,008) |
| kuyruk (520-1805) | **hayır, tolere ediyor** (o/e 0,927) |

---

## 1. Yayımlanmış bölgesel kısıtlılık ne diyor

Kendi hesabımızdan önce hazır kaynaklar canlı olarak tarandı. Tam URL, HTTP durumu ve dosya listesi
8. bölümde; burada yalnız sonuçlar.

### 1.1 gnomAD bölgesel missense kısıtlılığı (RMC): KIF13A için açık bir negatif

Brifte verilen `fordist_constraint_official_regional_missense_constraint_transcript_level.txt.gz`
adresi **404**. Gerçek dosyalar gnomAD'nin genel GCS kovasında, v2.1.1 sürümünde duruyor ve
indirildi. gnomAD bu analizi **18.629 transkript** üzerinde çalıştırmış; 5.127'sinde bölgesel kırılma
bulmuş (17.842 bölge satırı), 13.502'sinde bulamamış.

**KIF13A ENST00000259711, kırılma bulunamayanlar listesinde.** `..._transcripts_with_rmc.tsv`
dosyasında sıfır satırı var; `..._transcripts_without_rmc.tsv` dosyasında tek satır:
obs 774, exp 951,68, **o/e 0,8133**, z 2,0467. Yani bu tam olarak bizim aradığımız şeyi arayan,
amaca özel yazılmış bir algoritma KIF13A'ya bakmış ve motor domain dahil hiçbir yerde sınır
bulamamış. Bu bir veri boşluğu değil, **ölçülmüş bir negatif sonuçtur.**

**KIF1A ENST00000498729, 11 bölge** (aynı dosyadan birebir):

| aa | gözlenen | beklenen | o/e | ki-kare | p |
|---|---|---|---|---|---|
| 1-84 | 20 | 55,35 | 0,3614 | 22,57 | 2,0x10⁻⁶ |
| 84-128 | 0 | 23,68 | 0,0000 | 23,68 | 1,1x10⁻⁶ |
| 128-154 | 16 | 21,50 | 0,7441 | 1,41 | 0,235 |
| 154-177 | 0 | 16,36 | 0,0000 | 16,36 | 5,2x10⁻⁵ |
| 177-298 | 26 | 77,86 | 0,3339 | 34,55 | 4,2x10⁻⁹ |
| **298-369** | **3** | **43,69** | **0,0687** | **37,89** | **7,5x10⁻¹⁰** |
| 370-433 | 40 | 45,81 | 0,8731 | 0,74 | 0,390 |
| 434-691 | 50 | 155,33 | 0,3219 | 71,42 | 2,9x10⁻¹⁷ |
| 691-1445 | 330 | 476,20 | 0,6930 | 44,88 | 2,1x10⁻¹¹ |
| 1446-1499 | 4 | 24,36 | 0,1642 | 17,01 | 3,7x10⁻⁵ |
| 1499-1791 | 151 | 187,22 | 0,8065 | 7,01 | 8,1x10⁻³ |

Motoru kapsayan bölgelerin toplamı (yalnız aritmetik, yeniden türetme yok): aa 1-369 için
65 / 238,43 = **0,2726**; motor sonrası aa 370-1791 için 575 / 888,92 = 0,6469. Toplam
640 / 1127,35 = 0,5677, yani transkript düzeyindeki v2.1.1 değerini birebir veriyor.

**KIF13B ENST00000524189, 3 bölge:** 1-1511 o/e 0,8162; 1511-1629 o/e 0,1994; 1629-1826 o/e 0,9973.
Motor domain 1-1511 bloğunun içinde eriyor, yani RMC KIF13B için motor çözünürlüğü **vermiyor**.

**ExAC RMC (eski dosya):** yalnız KIF1A var (1-689 o/e 0,3249; 690-1792 o/e 0,7532).
KIF13A ve KIF13B bu dosyada hiç yok.

**v4 için RMC yayımlanmamış.** `release/4.1/regional_missense_constraint` boş gövdeyle 200 dönüyor;
`release/4.1/constraint/` yalnız transkript düzeyi metrikleri içeriyor.

**Sürüm tuzağı doğrulandı.** GraphQL'de `Gene` tipinde `gnomad_v2_regional_missense_constraint` ve
`exac_regional_missense_constraint_regions` alanları var ama **yalnız `reference_genome: GRCh37`
ile dolu geliyor**; GRCh38 ile her üç gen için `null`. Aynı şekilde gen düzeyi kısıtlılık GRCh38'de
v4.1, GRCh37'de v2.1.1 değerlerini veriyor:

| | gnomAD v4.1 (GRCh38) | gnomAD v2.1.1 (GRCh37) |
|---|---|---|
| KIF13A o/e missense | 0,861 (1901 / 2209,0) | 0,813 (774 / 951,7) |
| KIF13A missense Z | 2,79 | 2,05 |
| KIF1A o/e missense | 0,677 (1683 / 2485,9) | 0,568 |
| KIF1A missense Z | 6,87 | 5,16 |
| KIF13B o/e missense | 0,913 (2173 / 2379,0) | 0,818 |

**Transkript kimliği sorunu yok.** ExAC/gnomAD-v2 kanonik transkript listesi indirildi:
KIF13A -> ENST00000259711, KIF1A -> ENST00000498729, KIF13B -> ENST00000524189. Üçü de bugünün MANE
Select'i ile **aynı**; yalnız sürüm numarası farklı (v2/GRCh37'de .6 / .2 / .1, MANE'de .11 / .9 / .6)
ve protein uzunlukları aynı (1805 / 1791 / 1826). Yani v2 ve v4 sayıları aynı transkript üzerinde
karşılaştırılabilir.

### 1.2 MTR (Missense Tolerance Ratio) v2.0

Eski adres `biosig.unimelb.edu.au` yeni adrese 301 yönlendiriyor; kalıntı başına tablo
(34,7 milyon satır, GRCh37 / Ensembl v95, ~220 bin ekzom) indirildi. Tabloda `aapos` doğrudan
verildiği için koordinat eşlemesi gerekmiyor. Düşük MTR = daha az tolere edilen.

| Gen | Motor | Bağlayıcı 353/355-519 | Kuyruk 520-son |
|---|---|---|---|
| **KIF13A** | 0,9062 (FDR<0,05 olan kalıntı %1,4) | **0,8003 (%9,6)** | 0,8946 (%1,4) |
| **KIF1A** | **0,5258 (%61,1)** | 0,7199 (%19,4) | 0,7981 (%7,5) |
| KIF13B | 0,8544 (%4,9) | 0,8719 (%0) | 0,8722 (%3,8) |

KIF1A'nın motoru dramatik biçimde intoleran. **KIF13A'nın motoru kuyruğundan ayırt edilemiyor**;
tek tükenmiş bölgesi 353-519.

### 1.3 MPC

Eski düz dosya artık yok (404); MPC yalnız Hail tabloları olarak duruyor. dbNSFP üzerinden
MyVariant.info ile kurtarıldı (KIF13A 11.903, KIF1A 11.728 olası missense; bunlar gözlenen değil
**olası** varyantlar, yani seçim yanlılığı yok).

| Bölge | KIF13A | KIF1A |
|---|---|---|
| motor | ortalama **0,813**, MPC >= 2 olan **%0** (n=2274) | ortalama **2,263**, MPC >= 2 olan **%72,9** (n=2331) |
| bağlayıcı | 0,771, %0 (n=1100) | 2,135, %64,3 (n=1067) |
| kuyruk | 0,624, %0 (n=8510) | 1,130, %6,7 (n=8312) |
| tüm protein | 0,674, **en yüksek 1,184** | 1,448, %25,2, en yüksek 3,237 |

**KIF13A'nın proteinin herhangi bir yerindeki en yüksek MPC değeri 1,18'dir; tek bir varyant bile
klasik MPC >= 2 eşiğine ulaşmıyor.** Doğrudan sonuç: katman A'nın 52 adayının hiçbiri bir MPC
filtresinden geçemez, çünkü KIF13A'nın hiçbir varyantı geçemez.

Uyarı: MPC, RMC'nin üzerine kurulmuştur. Yani 1.3, 1.1'den **bağımsız bir kanıt değil**, aynı v2
hesabının aşağı akımıdır.

### 1.4 CCR (Constrained Coding Regions v2)

BED dosyası indirildi (8,19 milyon satır, GRCh37; yapı ampirik olarak doğrulandı: KIF1A CCR aralığı
GRCh37 CDS ile örtüşüyor, GRCh38 ile örtüşmüyor).

| Gen | Bölge | Ortalama CCR yüzdeliği | >= 90 | >= 95 | >= 99 |
|---|---|---|---|---|---|
| **KIF13A** | motor 5-352 | 65,76 | %25,6 | %4,3 | **%0** |
| | bağlayıcı 353-519 | 61,95 | %28,1 | %28,1 | **%17,4** |
| | kuyruk 520-1805 | 51,42 | %9,2 | %5,4 | **%0** |
| **KIF1A** | motor 5-354 | **83,39** | %58,3 | %43,4 | **%14,3** |
| | bağlayıcı 355-519 | 59,96 | %20,6 | %0 | %0 |
| | kuyruk 520-1791 | 60,97 | %27,8 | %22,7 | %7,5 |

Yine aynı tablo: KIF13A'nın motor domaininde 99. yüzdeliğe ulaşan **tek bir kalıntı yok**,
bağlayıcıda kalıntıların %17,4'ü var.

### 1.5 Bu bölümün sonucu

Dört kaynak da aynı şeyi söylüyor: **KIF13A motor domaini yayımlanmış bölgesel ölçütlerde dikkat
çekmiyor, dikkat çeken bölge 353-519.** Bu, 2. bölümdeki kendi hesabımızın neden gerekli olduğunu
da açıklıyor: dördü de gnomAD v2 / ExAC çağının veri hacmiyle kurulmuş (KIF13A için 774 gözlenen
missense), bizimki v4.1 ile (1901 gözlenen missense) yaklaşık 2,5 kat daha güçlü ve ayrıca
domain **içinde** 25 kalıntılık sınıflara inebiliyor; yayımlanmış kaynakların hiçbiri bunu yapmıyor.

Makine tarafından erişilemeyenler: MetaDome (iş tabanlı servis, tarayıcı gerekiyor) ve MTR3D
(genom çapında kalıntı tablosu yayımlanmıyor).

---

## 2. Yöntem

Samocha 2014 çerçevesi, genin kendi sinonim varyantlarıyla kalibre edilmiş hali.

### 2.1 Mutasyon hızı tablosu

**Kaynak:** Samocha KE ve ark., *Nature Genetics* 46:944-950 (2014), PMID **25086666**,
doi **10.1038/ng.3050**. Künye Europe PMC üzerinden canlı doğrulandı. Trinükleotid tablosu
(64 bağlam x 3 alternatif = 192 satır) denovonear paketinin (MIT) yeniden dağıtımından indirildi:

```
https://raw.githubusercontent.com/jeremymcrae/denovonear/master/src/denovonear/data/rates.txt
```

Yerel kopya: `data/raw/mutation_rates/samocha2014_trinucleotide_rates_20260916.txt` (+ `.prov.json`).

Tablo tam simetrik değil: bir bağlam ile ters tümleyeni arasındaki göreli sapmanın medyanı 0,015,
en büyüğü 0,129. Germ hattı mutasyonu şerit simetriktir, bu sapma sonlu intergenik veriden gelen
kestirim gürültüsüdür. Bu yüzden **hızlar ters tümleyen çiftleri üzerinden ortalandı**; bu aynı
zamanda geni hangi şeritten okuduğumuzun sonucu değiştirmesini engelliyor (KIF13A, KIF1A, KIF13B
eksi şerit; KIF1B, KIF1C, KIF5A, KIF5C artı şerit). `--raw-rates` ile kapatılabilir, sonuç
değişmiyor (6. bölüm).

### 2.2 CDS ve trinükleotid bağlamı

Her gen için Ensembl REST'ten MANE transkript yapısı (`lookup/id/{ENST}?expand=1;utr=1`) ve genin
GRCh38 artı şerit genomik dizisi (`sequence/region/human/{kromozom}:{baş}..{son}:1`) alındı. CDS
transkript yönünde kuruldu ve **her kodlayan bazın trinükleotid bağlamı genomdan** okundu, spliceli
diziden değil. Yani ekzon-ekzon birleşme noktasında bölünen bir kodon gerçek intronik komşusunu
alıyor, bir sonraki ekzonun ilk bazını değil.

Doğrulama: her genin CDS'i ATG ile başlamalı, tek bir sonlandırıcıyla bitmeli, içinde erken durdurma
olmamalı ve uzunluğu `config/targets.yaml`'daki protein uzunluğunu vermeli. KIF13A: 5418 nt -> 1805 aa,
39 kodlayan ekzon, eksi şerit. KIF1A: 5376 nt -> 1791 aa, 48 kodlayan ekzon. Hepsi geçti.

### 2.3 Olası her değişimin sayımı

CDS'in her bazı için 3 alternatif baz, standart genetik kodla sinonim / missense / nonsense olarak
sınıflandırıldı (kodon 1'de ayrıca `start_lost`). Sonlandırıcı kodon dahil edilmedi. KIF13A için
16.245 olası SNV, bunların 2.274'ü motor domaininde missense.

### 2.4 Gözlenen varyantlar

gnomAD v4.1 GraphQL **transkript** sorgusu (gen sorgusu değil, tuzak 3), `hgvsc` üzerinden CDS
konumuna yerleştirildi. Sayılan varyant: **ekzom**, PASS, AC >= 1, AF <= %0,1.

Genom verisi bilerek dışarıda: gnomAD'nin v4 kısıtlılık hattı ekzomlar üzerinde çalışıyor ve bu seçim
gözlenen sayıyı **birebir** yeniden üretiyor (3. bölüm). KIF13A'nın 3.914 kodlayan SNV'sinden
1.061'inin ekzom AC'si 0, 53'ü filtreli, 21'i frekans eşiğinin üstünde.

**Bağımsız doğrulama:** gnomAD'nin bildirdiği referans bazı ile Ensembl CDS'inden okuduğumuz baz
3.914 varyantın **hepsinde** uyuştu (script uyuşmazlıkta durur). Bu, CDS kurgusunun ve konum
eşlemesinin doğrudan kanıtıdır.

### 2.5 Kalibrasyon: beklenti neden sinonimlerden gelmeli

Naif yaklaşım "beklenen = sabit x mutasyon hızı"dır ve 730 bin ekzomluk bir veri setinde yanlıştır:
yüksek hızlı CpG geçişleri **doygunluğa** ulaşır, hızın iki katına çıkması gözlenme olasılığını iki
katına çıkarmaz. Beklenti bu yüzden hızın doğrusal bir katı olarak değil, bir **fonksiyonu** olarak
kuruldu:

> Panelin tüm genlerindeki **olası sinonim** değişimler mutasyon hızına göre 24 nicelik kutusuna
> bölünür; her kutuda gerçekten gözlenmiş olanların oranı, o hızdaki bir değişimin gnomAD'de
> görülme olasılığıdır. Bir bölgenin beklenen missense sayısı = o bölgedeki olası missense
> değişimlerin bu olasılıklarının toplamı.

Panel: 11 otozomal kinezin, **31.410 olası sinonim değişim**, bunların **7.378'i** gözlenmiş (%23,5).
Kutuların gözlenme olasılığı %4,8'den %68,7'ye çıkıyor, doygunluk açıkça görülüyor.
Tablo: `data/interim/regional_constraint/calibration_bins_20260916.tsv`.

Sinonim üzerinden kalibre etmenin anlamı: dizileme derinliği, örneklem büyüklüğü ve demografik geçmiş
missense ile sinonimi aynı şekilde etkiler, oranda sadeleşir.

**X kromozomu neden panel dışı:** ilk çalıştırmada panele KIF4A (chrX) da konmuştu. Hemizigot
erkeklerden ötürü X'in allel sayısı farklı olduğundan KIF4A'nın sinonim o/e'si 0,56 çıktı ve ortak
kalibrasyonu aşağı çekerek bütün otozomal genlerin beklentisini birkaç yüzde eksiltti. Panel artık
otozomlarla sınırlı ve script otozom dışı bir gen görürse duruyor.

### 2.6 İstatistik

- **Güven aralığı:** o/e için %90 Poisson olabilirlik aralığı, gnomAD'nin LOEUF için kullandığı ızgara
  yöntemiyle (o/e 0-3 arası 0,001 adım, olabilirlik normalize edilip %5 ve %95 yüzdelikleri alınır).
  Aralığın **üst sınırı** missense için LOEUF'un karşılığıdır; hüküm ölçeği buna bakar, böylece az
  gözlemli bir bölge geniş aralıkla "kısıtlı" ilan edilemez.
- **Tükenme p değeri:** tek yanlı Poisson P(X <= gözlenen).
- **İki bölgenin karşılaştırılması:** koşullu binom testi. Toplam gözlemin A bölgesine düşen payı,
  sıfır hipotezinde Binom(n, beklenen_A/(beklenen_A+beklenen_B)) dağılır. Bu test **ölçekten
  bağımsızdır**, yani mutlak kalibrasyondaki hata iki bölgeyi aynı ölçüde etkilediği sürece sonucu
  bozmaz.
- **Sinonim normalize sütun:** bir bölgenin beklentisi kendi sinonim o/e'siyle çarpılarak yeniden
  hesaplanır; yerel bir kalibrasyon sapması (kapsama, CpG metilasyonu) varsa bu sütun onu siler.

### 2.7 Hüküm ölçeği

| o/e %90 aralığının üst sınırı | Etiket |
|---|---|
| >= 0,95 | hayır, tolere ediyor |
| 0,80 - 0,95 | zayıf kısıtlı |
| 0,60 - 0,80 | orta kısıtlı |
| < 0,60 | güçlü kısıtlı |

---

## 3. Sağlama kontrolleri

| # | Kontrol | Sonuç | Durum |
|---|---|---|---|
| 1 | **KIF1A motor domaini pozitif kontrol** | o/e **0,450** (0,402-0,505), 207 gözlenen / 459,9 beklenen -> güçlü kısıtlı | GEÇTİ |
| 2 | **gnomAD gen düzeyini yeniden üretme** | KIF13A gözlenen missense **1901 vs 1901** (%0,0); beklenen 2124,8 vs 2209,0 (-%3,8); o/e 0,895 vs 0,861 (+%4,0) | GEÇTİ |
| 3 | **Sinonim o/e her pencerede 1,0 civarında mı (KIF13A)** | medyan 0,958; örtüşmeyen 18 pencerede ki-kare/sd **1,21**; en sapan pencere 1-100, o/e 1,42 | GEÇTİ |
| 4 | **Aynısı KIF1A'da** | medyan 1,084; ki-kare/sd **1,12**; en sapan pencere 1001-1100, o/e 1,23 | GEÇTİ |
| 5 | **Motor domainde kaç gözlenen missense var** | **333** gözlenen, 409,0 beklenen; %90 aralık genişliği 0,147 | GEÇTİ |
| 6 | **Motor domainde sinonim o/e (yerel kalibrasyon)** | **0,999** (0,879-1,139), 160 gözlenen / 160,2 beklenen | GEÇTİ |
| 7 | **Motor kendi kuyruğundan anlamlı olarak tükenmiş mi** | motor/kuyruk 0,891, **p = 0,060** (353-1805'e karşı); 0,879, p = 0,036 (520-1805'e karşı). KIF1A'da aynı test 0,575, p = 5,7x10⁻¹⁴ | KALDI (sınırda) |
| 8 | **Ölçüm klinik gerçekliği izliyor mu** | kısıtlılık sırası ile motor domaindeki ClinVar P/LP missense sayısı: KIF1A 0,45 -> 132; KIF5A 0,46 -> 36; KIF5C 0,48 -> 7; KIF13A 0,81 -> 0; KIF13B 0,90 -> 0 | GEÇTİ |
| 9 | **Nonsense iç tutarlılık** | KIF13A motor nonsense o/e 0,410 (10 / 24,4); KIF1A 0,206 (4 / 19,4). Fonksiyon kaybı ekseninin "motor domainde tolere edilen kesik varyantlar var" bulgusuyla uyumlu | GEÇTİ |

Bir ek dış doğrulama: boru hattımız KIF1A'nın **en tükenmiş 50 kodonluk penceresini 301-350'de**
buluyor (12 gözlenen / 63,4 beklenen, o/e 0,189). gnomAD'nin RMC algoritması, bambaşka bir veri
sürümünde (v2.1.1) ve bambaşka bir yöntemle (değişim noktası araması), aynı bölgeyi **298-369**
olarak seçmiş (o/e 0,069). İki hesap birbirini görmeden aynı yere düşüyor.

### Kontrol 2 hakkında: kalan %4 nereden geliyor

Gözlenen sayı 11 genin 10'unda gnomAD ile **birebir** aynı (KIF13A 1901, KIF13B 2173, KIF1B 1622,
KIF1C 1315, KIF5A 861, KIF5C 780, KIF3A 581, KIF21A 1677, KIF16B 1421, KIF14 1751; KIF1A'da 1685'e
karşı 1683, iki varyant fark). Sistematik fark **beklenen** tarafında ve yönü her gende aynı: bizim
beklentimiz %3-6 daha düşük, dolayısıyla o/e'miz %3-6 daha yüksek.

Nedeni yapısaldır ve kusur değildir: bizim kalibrasyonumuz **panelin sinonim o/e'sini tanım gereği
1,0'a** oturtur, oysa gnomAD'nin genom çapındaki modeli bu genleri ortalama 0,95 sinonim o/e'de
bırakır. Aradaki fark tam olarak gözlenen %3-6'dır. Bu yüzden her bölgesel sayı için **gnomAD
ölçeğine çevrilmiş** bir sütun da veriliyor: gen düzeyindeki oranımızı gnomAD'nin yayımlanan
değerine oturtan tek bir katsayı (KIF13A için x0,962, KIF1A için x0,944). Gen içi ve genler arası
karşılaştırmalar bu katsayıdan etkilenmez.

### Kontrol 3 hakkında: neden sabit bant değil ki-kare

100 kodonluk bir pencerede yalnız 45 civarı sinonim varyant beklenir, yani sırf Poisson gürültüsü
pencerelerin üçte birini 0,80-1,25 bandının dışına atar. Doğru soru dağılımın Poisson'dan **fazla**
olup olmadığıdır. Örtüşmeyen pencerelerde ki-kare/serbestlik derecesi 1,21 (KIF13A) ve 1,12 (KIF1A),
yani gözlenen yayılım gürültüyle uyumlu.

**Ama tek bir gerçek yerel sapma var ve saklanmıyor:** KIF13A'nın ilk 80 kodonunda sinonim o/e
**1,486**. Orada beklentimiz olduğundan düşük, büyük olasılıkla promotöre bitişik CpG adası yüzünden
(hipometile CpG'ler tablonun ortalama hızından yavaş mutasyona uğrar). Sonuç: o bölgenin missense
o/e'si (1,020) yukarı doğru şişkindir; sinonim normalize değeri 0,686. Bu yüzden aşağıda
"motor domainin N yarısı toleranslı" okuması **güvenli değil** diye işaretlendi. Motor domainin
**bütününde** sinonim o/e 0,999 olduğu için domain düzeyindeki rakam bu sapmadan etkilenmiyor.

### Kapsama

gnomAD v4.1 ekzom kapsaması, kodlayan bazlarda 20x üzeri örnek oranı: KIF13A medyan **0,990**
(5418 bazın 270'i 0,80'in altında), KIF1A medyan **0,980** (5376 bazın 20'si altında).
Ham: `data/raw/gnomad/KIF13A_exome_coverage_20260916.tsv`.

---

## 4. Sonuçlar

### 4.1 KIF13A, domain başına

| Bölge | aa | Gözlenen missense | Beklenen | o/e | %90 aralık | gnomAD ölçeğinde | Sinonim o/e | p (tükenme) | Kuyruğa karşı p | Hüküm |
|---|---|---|---|---|---|---|---|---|---|---|
| **motor** | 5-352 | **333** | 409,0 | **0,814** | 0,745-0,892 | **0,783** | 0,999 | 9,5x10⁻⁵ | 0,036 | zayıf kısıtlı |
| boyun bağlayıcı | 353-380 | 32 | 33,8 | 0,946 | 0,714-1,270 | 0,910 | 0,956 | 0,42 | 0,98 | tolere ediyor |
| **CC1** | 381-468 | 73 | 98,9 | **0,738** | 0,611-0,897 | 0,710 | 0,821 | 0,0039 | 0,065 | zayıf kısıtlı |
| FHA | 469-519 | 53 | 61,5 | 0,862 | 0,691-1,083 | 0,829 | 1,374 | 0,15 | 0,65 | tolere ediyor |
| **kuyruk** | 520-1805 | **1406** | 1517,5 | **0,927** | 0,887-0,968 | **0,891** | 0,969 | 0,0022 | - | tolere ediyor |
| gen bütünü | 1-1805 | 1901 | 2124,8 | 0,895 | 0,860-0,931 | 0,861 | 0,980 | - | - | zayıf kısıtlı |

Ek olarak, 1. bölümdeki yayımlanmış kaynakların işaret ettiği bölgeler aynı ölçekte:

| Bölge | Gözlenen / beklenen | o/e | %90 aralık | Sinonim o/e | Kuyruğa karşı p |
|---|---|---|---|---|---|
| boyun+CC1+FHA 353-519 | 158 / 194,2 | 0,813 | 0,715-0,928 | 1,013 | 0,13 |
| CC1+FHA 381-519 | 126 / 160,4 | **0,785** | 0,680-0,911 | 1,028 | 0,083 |

Okunuşu: gen düzeyindeki 0,86, motor domaindeki 0,78 ile kuyruktaki 0,89'un ağırlıklı ortalamasıdır;
kuyruk proteinin %71'i olduğu için ortalamayı kuyruk belirliyor. **Yani evet, gen düzeyindeki sayı
bir sinyali gizliyordu. Ama gizlediği sinyal küçük ve motor domaine özgü değil:** CC1+FHA bölgesi
(0,785) motor domainden (0,814) biraz daha tükenmiş durumda. 50 kodonluk ince pencerelerde KIF13A'nın
iki en tükenmiş noktası **296-345** (o/e 0,574, aralık 0,435-0,768, sinonim 0,916, p = 3,3x10⁻⁴) ve
**376-425** (o/e 0,633, aralık 0,480-0,847, sinonim 0,941, p = 0,0031); birincisi motor domainin
mikrotübül arayüzü, ikincisi CC1. MTR ve CCR'nin bağımsız olarak işaret ettiği yer bu ikinci nokta.

### 4.2 Motor domainin içi: asıl bulgu

"Ortalama sinyali gizler" argümanı bir kat daha aşağıda çok daha sert geçerli.

| Yapısal sınıf (1-380) | Kalıntı | Gözlenen missense | Beklenen | o/e | %90 aralık | Sinonim o/e | p | Katman A adayı |
|---|---|---|---|---|---|---|---|---|
| **ATP cebi** | 25 | 18 | 34,4 | **0,523** | 0,362-0,776 | 1,498 | 0,0016 | **23** |
| gömülü çekirdek | 124 | 99 | 146,0 | 0,678 | 0,576-0,801 | 0,876 | 2,3x10⁻⁵ | 9 |
| mikrotübül arayüzü | 50 | 41 | 55,0 | 0,746 | 0,581-0,968 | 1,129 | 0,030 | **18** |
| boyun bağlayıcı | 28 | 32 | 33,8 | 0,946 | 0,714-1,270 | 0,956 | 0,42 | 1 |
| **yüzey** | 153 | 179 | 177,7 | **1,008** | 0,892-1,140 | 0,972 | 0,56 | 1 |

Yapısal sınıflar `data/interim/structure/KIF13A_motor_residue_structure_20260916.tsv` dosyasından,
yani yapı ekseninin KIF1A 8UTN üzerine yerleştirmesinden geliyor (`docs/yapi_ekseni.md`).

Alt bölge kesitleri:

| Bölge | Gözlenen / beklenen | o/e | %90 aralık | Sinonim o/e | Sinonim normalize | Not |
|---|---|---|---|---|---|---|
| motor N yarısı 5-180 | 187 / 205,8 | 0,909 | 0,807-1,026 | 1,134 | 0,801 | kalibrasyon burada şüpheli |
| motor C yarısı 181-352 | 146 / 203,2 | **0,719** | 0,628-0,824 | 0,872 | 0,824 | p = 1,5x10⁻⁵ |
| MT arayüzü / anahtar bölgesi 270-360 | 66 / 105,2 | **0,627** | 0,514-0,770 | 0,876 | 0,716 | p = 2,8x10⁻⁵ |
| ilk 80 kodon (5' CpG adası) | 99 / 97,1 | 1,020 | 0,867-1,205 | **1,486** | 0,686 | **kullanılmaz**, 3. bölüm |

### 4.3 Yorum güvenliği: "motor domainler zaten kısıtlıdır" itirazı

Bu proje daha önce yalnız konumla açıklanabilen bir ölçümle (0,995 AUROC, "motor domainde mi")
yanıldı. Buradaki ayna görüntüsü risk şu: motor domainler yoğun ve gömülü olduğu için herhangi bir
globüler domain bir miktar kısıtlılık gösterir. İki bağımsız ölçümle çerçevelendi.

**a) Kinezin paneli.** Aynı hesap 11 otozomal kinezinin motor domaininde.

| Gen | Motor aa | o/e | %90 aralık | gnomAD ölç. | Gözlenen | Kuyruk o/e | Motor/kuyruk | p | Motorda ClinVar P/LP missense |
|---|---|---|---|---|---|---|---|---|---|
| KIF1A | 5-354 | **0,450** | 0,402-0,505 | 0,425 | 207 | 0,782 | 0,575 | 5,7x10⁻¹⁴ | **132** |
| KIF5A | 44-327 | 0,462 | 0,407-0,527 | 0,444 | 163 | 0,775 | 0,596 | 2,9x10⁻⁹ | **36** |
| KIF5C | 16-327 | 0,482 | 0,427-0,544 | 0,467 | 184 | 0,788 | 0,611 | 5,0x10⁻⁹ | 7 |
| KIF1B | 11-354 | 0,638 | 0,577-0,707 | 0,612 | 262 | 0,771 | 0,828 | 0,0057 | 1 |
| KIF3A | 20-345 | 0,668 | 0,604-0,740 | 0,642 | 260 | 0,714 | 0,936 | 0,46 | 0 |
| KIF21A | 15-371 | 0,771 | 0,704-0,846 | 0,747 | 320 | 0,890 | 0,867 | 0,023 | 1 |
| **KIF13A** | **5-352** | **0,814** | **0,745-0,892** | **0,783** | **333** | **0,914** | **0,891** | **0,060** | **0** |
| KIF14 | 364-701 | 0,864 | 0,790-0,945 | 0,835 | 339 | 0,935 | 0,923 | 0,22 | 1 |
| KIF1C | 11-348 | 0,895 | 0,823-0,974 | 0,847 | 383 | 0,948 | 0,944 | 0,36 | 7 |
| KIF16B | 9-358 | 0,900 | 0,828-0,979 | 0,865 | 382 | 0,920 | 0,978 | 0,73 | 0 |
| KIF13B | 11-353 | 0,905 | 0,832-0,986 | 0,860 | 374 | 0,970 | 0,933 | 0,23 | 0 |

Motor domain sınırları: KIF13A ve KIF1A için projenin kendi tablosu; diğerleri için Ensembl'in
Pfam PF00225 ek açıklaması, **çeviri numaralandırmasında** (UniProt izoform kayması riski yok).
KIF14'ün motor domaini N ucunda değil (364-701); onun dışında hepsi N terminal.

Sonuç: **motor domainler genel olarak kısıtlı değildir.** Aralık 0,45'ten 0,91'e yayılıyor.
KIF13A 11 gen içinde 7. sırada, KIF13B / KIF16B / KIF1C ile aynı yarıda, KIF1A / KIF5A / KIF5C
grubundan açıkça ayrı.

Klinik eşlik: kısıtlılık sıralamasının tepesi ClinVar'daki patojenik missense yükünü izliyor
(KIF1A 132, KIF5A 36, KIF5C 7). KIF1C bir istisna (o/e 0,895 ama 7 P/LP) ve nedeni tutarlı:
KIF1C hastalığı çekiniktir, heterozigot taşıyıcıya karşı seçilim beklenmez. KIF21A'nın 11 P/LP
kaydının yalnız 1'i motor domainde. Sayılar
`{GEN}[gene] AND "missense variant"[molecular consequence] AND ("clinsig pathogenic"[Properties] OR "clinsig likely pathogenic"[Properties])`
sorgusuyla alındı; tuzak 14'teki `pathogenic[clinsig]` biçimi kullanılmadı. Kalıntı konumları kayıt
başlığındaki transkriptten okundu ve bu transkriptler ilgili genin MANE RefSeq'i çıktı (KIF1A için
NM_001244008.2, KIF5A için NM_004984.4, KIF1C için NM_006612.6).

**b) Aynı yapısal sınıflar KIF1A'da.** Yapı ekseninin süperpozisyon eşleşmesi kullanılarak aynı
sınıflar KIF1A kalıntılarına taşındı.

| Yapısal sınıf | KIF13A o/e | KIF1A o/e | KIF13A / KIF1A |
|---|---|---|---|
| yüzey | **1,008** (0,892-1,140) | 0,613 (0,528-0,715) | 1,6 kat |
| gömülü çekirdek | 0,678 (0,576-0,801) | 0,335 (0,269-0,421) | 2,0 kat |
| mikrotübül arayüzü | 0,746 (0,581-0,968) | 0,338 (0,242-0,483) | 2,2 kat |
| **ATP cebi** | **0,523** (0,362-0,776) | **0,140** (0,073-0,294) | 3,7 kat |

Gradyan **aynı yönde**: yüzey > gömülü çekirdek ~ mikrotübül arayüzü > ATP cebi. Yani KIF13A'nın
motor domaini bir motor domaini gibi davranıyor, işlevsel merkezde daha kısıtlı. Ama iki fark
belirleyici: (1) her sınıfta 2 ila 4 kat daha sığ; (2) KIF13A'nın motor **yüzeyi hiç kısıtlı değil**
(1,008), KIF1A'nınki belirgin biçimde kısıtlı (0,613). "Gömülülük zaten kısıtlılık üretir" itirazı
ölçüldü ve kısmen doğru çıktı; ama KIF13A'nın ATP cebi salt gömülülüğün açıkladığından daha fazla
tükenmiş (0,523'e karşı gömülü çekirdeğin 0,678'i), dolayısıyla orada işlevsel bir sinyal var.

### 4.4 Kayan pencere profili

`figures/KIF13A_bolgesel_kisitlilik.svg`

- **Panel A:** KIF13A, 100 kodonluk pencere, 10 kodon adım, %90 aralık bandı, kontrol olarak sinonim
  çizgisi ve gen düzeyindeki 0,861 referans çizgisi. Motor domain gölgeli.
- **Panel B:** aynısı KIF1A için (pozitif kontrol). Motor domaindeki derin çukur açıkça görünüyor.
- **Panel C:** 11 kinezinin motor domain o/e değeri, aralıklarıyla ve motor domaindeki ClinVar P/LP
  sayısıyla.
- **Panel D:** motor domainin içinde yapısal sınıfa göre kısıtlılık, KIF13A ve KIF1A yan yana.

Panel A'da göze çarpan iki şey: (i) sinonim çizgisinin proteinin ilk ~80 kodonunda 1,4'ün üstüne
fırlaması (3. bölümdeki CpG adası sapması, bilerek gizlenmedi), (ii) missense çizgisinin motor
domainin **ikinci yarısında** (yaklaşık 250-350) ve hemen ardından CC1'de tabana inmesi.

Ham tablolar: `sliding_window_100aa_20260916.tsv` (171 KIF13A penceresi) ve
`sliding_window_50aa_20260916.tsv`.

Dürüstlük notu: KIF13A'nın **tüm proteindeki en tükenmiş 100 kodonluk penceresi motor domaininde
değil**, kuyrukta (981-1080, o/e 0,607). O pencerede sinonim o/e 0,757, yani yerel kalibrasyon
sapmalı ve sinonim normalize değer 0,80'e geliyor; yine de kuyrukta da bakılmaya değer bir bölge
olduğu not edilmeli.

---

## 5. Bu, 52 katman A adayına ne yapıyor

Katman A adaylarının hepsi 11-370 aralığında; 51'i motor domaininde, 1'i boyun bağlayıcıda (L370P).
Yapısal dağılımı: **23 ATP cebi, 18 mikrotübül arayüzü, 9 gömülü çekirdek, 1 yüzey, 1 boyun
bağlayıcı.**

Bölgesel ölçüm bu listenin önsel olasılığını **yukarı** çekiyor, ama ölçülebilir kadar ve tavanı
belli:

1. **Eski itiraz zayıfladı.** "KIF13A missense'i tolere ediyor (o/e 0,86), dolayısıyla motor
   domaindeki bir missense muhtemelen zararsızdır" cümlesi artık kurulamaz. Adayların oturduğu
   bölgelerin gerçek değerleri 0,52 (ATP cebi), 0,68 (gömülü çekirdek) ve 0,75 (mikrotübül arayüzü);
   0,86 değil. Katman A'nın 50/52'si tükenmiş sınıflarda, yalnız 1'i tükenmemiş sınıfta (yüzey).
   Liste bölgesel kısıtlılıkla **kendiliğinden hizalı** kurulmuş; bu, yöntemin bağımsız bir
   doğrulamasıdır.
2. **Ama tavan alçak ve iki yönden geliyor.** (i) Aynı ölçüm KIF1A'nın ATP cebinde 0,140 veriyor;
   KIF13A'daki seçilim baskısı KIF1A'nınkinin yaklaşık dörtte biri. Tam penetranslı, ağır baskın bir
   hastalık mekanizması bu kalıntılarda çalışıyor olsaydı gnomAD'de daha derin bir çukur bırakırdı.
   (ii) KIF13A'nın **hiçbir varyantı** MPC >= 2 eşiğine ulaşmıyor (proteinin en yüksek MPC değeri
   1,18), çünkü MPC bölgesel kısıtlılıktan türer ve gnomAD KIF13A'da bölgesel kısıtlılık bulamamıştır.
   Yani katman A, MPC kullanan herhangi bir dış filtrede **tamamen** elenir. Desteklenen ifade
   "muhtemelen patojenik" değil, **"popülasyon verisinde zayıf ama ölçülebilir bir seçilim izi
   taşıyan kalıntılar"**dır.
3. **Liste içinde ayrışma.** Katman A'nın 16'sı 270-360 aralığında (o/e 0,627), 27'si 5-180
   aralığında. 5-180'in ölçümü kalibrasyon sapması yüzünden güvenilmez (sinonim o/e 1,134, ilk 80
   kodonda 1,486), o yüzden o 27 aday için bölgesel kısıtlılık **kanıt sunmuyor**, ne lehte ne
   aleyhte. Bölgesel kanıtın lehte olduğu adaylar ATP cebindeki 23 ile 270-360'taki 16'dır; bu iki
   küme büyük ölçüde örtüşür.
4. **Yeni bir soru açıldı: CC1.** Hem bizim ölçümümüz (381-468, o/e 0,738, p = 0,0039; 376-425
   penceresi 0,633) hem MTR hem CCR, KIF13A'nın en tükenmiş bölgesi olarak CC1'i gösteriyor.
   Katman A'da oradan **tek bir aday bile yok**, çünkü liste motor domain + paralog kanıtı üzerine
   kurulmuştu ve KIF1A'nın CC1'i hizalanabilir değil. Bu, aday listesinin bilinen bir kör noktası
   olarak kayda geçmeli.
5. **De novo ekseniyle aynı yöne bakıyor.** `docs/denovo_ekseni.md` 31.058 nörogelişimsel bozukluk
   triosunda KIF13A için 2 de novo missense (beklenen 3,16) buluyor, KIF1A için 53 (beklenen 4,81).
   Bölgesel kısıtlılığın "sığ seçilim izi" tablosu ile bu bağımsız olarak uyumlu: motor domaininde
   ağır baskın missense hastalığı olan bir gen hasta triolarında de novo missense zenginleşmesi
   gösterir, KIF13A göstermiyor. İki eksen aynı sonuca farklı veriden varıyor.
6. **Fonksiyon kaybı ekseniyle çelişki yok.** Motor domainde nonsense o/e 0,410 çıkıyor, yani
   `docs/fonksiyon_kaybi_ekseni.md`'nin "motor domaini erken kesmek heterozigot halde tolere
   ediliyor" bulgusu bağımsız bir hesapla yeniden üretilmiş oldu. Bu, katman A'nın mekanizma
   gerekçesini değiştirmiyor: doz kaybı değil, baskın negatif hipotezi hâlâ tek tutarlı seçenek ve
   hâlâ ölçülmedi.

---

## 6. Duyarlılık ve sınırlar

KIF13A motor domaini (5-352) için o/e, dört alternatif parametre setinde:

| Çalıştırma | Motor o/e | %90 aralık | gnomAD ölç. | KIF1A motor o/e |
|---|---|---|---|---|
| ana çalıştırma (simetrik hızlar, AF <= %0,1, 24 kutu) | 0,814 | 0,745-0,892 | 0,783 | 0,450 |
| `--raw-rates` (hız tablosu simetrikleştirilmedi) | 0,817 | 0,747-0,894 | 0,784 | 0,450 |
| `--max-af 0.0001` (AF <= %0,01) | 0,817 | 0,746-0,895 | 0,780 | 0,461 |
| `--bins 12` | 0,811 | 0,741-0,888 | 0,781 | 0,459 |
| `--bins 40` | 0,819 | 0,749-0,897 | 0,785 | 0,449 |

Aralık 0,811-0,819; hüküm hiçbir varyantta değişmiyor. Çıktılar `*_rawrates`, `*_af1e4`,
`*_bins12`, `*_bins40` sonekli dosyalarda.

### Bilinen sınırlar

1. **Metilasyon modellenmedi.** gnomAD'nin kendi modeli CpG'leri metilasyon düzeyine göre ayırır,
   bizimki ayırmaz. Ölçülen tek somut sonucu KIF13A'nın ilk ~80 kodonundaki sinonim şişmesidir
   (o/e 1,486). Motor domainin bütünü etkilenmiyor (sinonim o/e 0,999), ama N yarısının okuması
   güvenli değil.
2. **ATP cebi 25 kalıntı.** 18 gözlenen missense ile aralık geniş (0,362-0,776). Sinonim o/e 1,498
   de yalnız 21 gözlem üzerinde, o yüzden sinonim normalize değer (0,349) gürültülüdür; konuşurken
   ham 0,523 kullanılmalı.
3. **Motor-kuyruk farkı sınırda.** p = 0,036 ile 0,060 arası, karşılaştırma sınırına göre değişiyor.
   Tek başına "motor domain kuyruktan kısıtlı" manşeti bu sayıyla kurulamaz.
4. **Yapısal sınıflar KIF1A yapısına yerleştirmeden geliyor.** `docs/yapi_ekseni.md` katı cisim
   yerleştirmenin çakışma ürettiğini ve mesafe değerlerinin güvenilmez olduğunu söylüyor; sınıf
   ataması ("temas var mı") güvenilir sayılıyor. Sınıf sınırları değişirse 4.2'deki sayılar değişir.
5. **gnomAD ekzon ± 75 baz.** Derin intronik ve düzenleyici bölge yok; bu analiz yalnız CDS için.
6. **1. bölümdeki dış kaynaklar farklı yapıdadır** (MTR ve CCR GRCh37, RMC ve MPC gnomAD v2.1.1).
   Transkript kimliği doğrulandı ve üçü de MANE ile aynı; ama sürüm farkı yüzünden o sayılar bizim
   v4.1 sayılarımızla **doğrudan** aynı ölçekte değildir, yön karşılaştırması için kullanılabilirler.

---

## 7. Bu eksen ne kadar ileri götürülebilir

### En güçlü lisanslı cümle

> gnomAD v4.1'de (730 bin ekzom, MANE transkripti) KIF13A motor domaini, gen düzeyindeki 0,86
> ortalamasının gizlediği gerçek fakat sığ bir missense tükenmesi taşır: gözlenen/beklenen 0,78
> (%90 aralık 0,75-0,89; 333 gözlenen, 409 beklenen; p = 9,5x10⁻⁵). Tükenme domain boyunca düzgün
> dağılmamıştır: motor yüzeyi hiç tükenmemiş (o/e 1,01), gömülü çekirdek 0,68, mikrotübül arayüzü
> 0,75, ATP cebi ise 0,52'dir ve katman A adaylarının 23'ü tam bu ATP cebinde oturur. Aynı hesap
> KIF1A'nın ATP cebinde 0,14 vermekte, yani KIF13A'daki seçilim izi karşılaştırma genindekinin
> yaklaşık dörtte biri kadardır; ayrıca aynı ölçüm KIF13A'nın CC1 bölgesini (0,74) motor domainden
> biraz daha tükenmiş göstermektedir.

Bu cümlenin her parçası bu depoda bir dosyaya bağlıdır ve pozitif kontrolle, gen düzeyi yeniden
üretimiyle ve pencere başına sinonim kontrolüyle desteklenmiştir. Bağımsız üç eksen aynı yöne
bakıyor: popülasyon kısıtlılığı sığ, hasta triolarında de novo missense zenginleşmesi yok
(`docs/denovo_ekseni.md`), ClinVar'da sınıflandırılmış tek bir patojenik missense yok.

### Kaçınılacak aşırı iddia

- **"KIF13A motor domaini kısıtlıdır, dolayısıyla KIF1A gibi bir hastalık genidir."** Yanlış. Aynı
  boru hattı KIF1A motorunda 0,45, KIF13A motorunda 0,81 veriyor; kinezin panelinde KIF13A 11 gen
  içinde 7. sırada ve kısıtlılık göstermeyen KIF13B / KIF16B / KIF1C ile aynı kümede.
- **"Motor domain proteinin kısıtlı bölgesidir."** Değil. CC1+FHA (0,785) motor domainden (0,814)
  biraz daha tükenmiş ve üç bağımsız yayımlanmış kaynak (RMC'nin negatifi, MTR, CCR) aynı yeri
  gösteriyor.
- **"Motor domain kuyruktan anlamlı olarak kısıtlıdır."** p = 0,036-0,060, karşılaştırmanın nasıl
  kurulduğuna bağlı. Manşet yapılamaz.
- **"ATP cebi güçlü kısıtlıdır."** Aralığın üst sınırı 0,776, yani ölçek gereği en fazla
  "orta kısıtlı". 18 gözlem üzerinde nokta tahmini tek başına anılmamalı.
- **"Motor domainin N yarısı missense'i tolere ediyor."** O bölgede kalibrasyon sapmalı
  (sinonim o/e 1,134, ilk 80 kodonda 1,486). Veri o ifadeyi desteklemiyor, sadece ölçemiyor.
- **"gnomAD'nin bölgesel kısıtlılık bulamaması bir veri boşluğudur."** Değil, ölçülmüş bir negatif:
  KIF13A 18.629 transkriptlik analize dahil edilmiş ve kırılma bulunamayanlar listesine yazılmış.
- **Herhangi bir tekil kalıntı için bu sayıların kullanılması.** En küçük güvenilir birim yaklaşık
  25-50 kalıntılık bir sınıftır; tek kalıntı için popülasyon verisinde hiç güç yok.

### Sonraki ölçülebilir adım

1. **Metilasyon farkındalıklı hızlar.** gnomAD'nin bağlam + metilasyon tablosu kullanılırsa 5' uçtaki
   sapma kapanır ve motor domainin N yarısı ilk kez okunabilir hale gelir. Katman A'nın 27 adayı
   orada oturduğu için bu doğrudan aday listesini ilgilendiriyor.
2. **Değişim noktası araması.** gnomAD'nin RMC algoritmasını v4.1 sayımları üzerinde çalıştırmak:
   v2'de 774 gözlenen missense ile kırılma bulunamamıştı, v4.1'de 1901 var. Sınırları elle
   koymak yerine veriye bulduran bir hesap, CC1 ile motor C yarısı arasındaki seçimi yapabilir.
3. **Hizalanmış konumlarda kinezin ailesi havuzu.** 11 kinezinin motor domainleri hizalanıp aynı
   hizalama sütununda kısıtlılık toplanırsa, tek gende 25 kalıntıyla elde edilemeyen güç elde edilir
   ve "bu kalıntı ailede kısıtlı mı" sorusu kalıntı çözünürlüğünde sorulabilir.
4. **CC1 kör noktası.** Aday listesi CC1'i hiç kapsamıyor. Eğer bölgesel kanıt izlenecekse oraya da
   bir aday üretimi gerekir; ama paralog aktarımı orada çalışmaz, o yüzden başka bir kanıt ekseni
   gerekir.
5. **Baskın negatif testi hâlâ ölçülmedi** ve bu eksen onu ölçemez. Bölgesel kısıtlılık seçilim
   olduğunu gösterir, mekanizmasını göstermez.

---

## 8. Yayımlanmış kaynakların taranması: erişim kaydı

Hepsi 16 Eylül 2026'da canlı çekildi; dosyalar `data/raw/rmc/` altında, her biri `.prov.json`
yan dosyasıyla (kaynak, URL, UTC zaman damgası, sha256, kayıt sayısı).

| Kaynak | Durum | Adres |
|---|---|---|
| gnomAD RMC v2.1.1, kırılma bulunanlar | 200, 2.184.736 B, 17.842 bölge / 5.127 transkript | `https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/regional_missense_constraint/gnomAD_v2.1.1_transcripts_with_rmc.tsv` |
| gnomAD RMC v2.1.1, kırılma bulunamayanlar | 200, 1.022.909 B, 13.502 transkript | `.../gnomAD_v2.1.1_transcripts_without_rmc.tsv` |
| Kova listelemesi | 200 | `https://storage.googleapis.com/storage/v1/b/gcp-public-data--gnomad/o?prefix=release%2F2.1.1%2Fregional_missense_constraint%2F&delimiter=%2F` |
| ExAC RMC (eski) | 200, 647.497 B, 6.502 satır | `https://storage.googleapis.com/gcp-public-data--gnomad/legacy/exac_browser/regional_missense_constraint.tsv` |
| Briften gelen dosya adı | **404** | `.../fordist_constraint_official_regional_missense_constraint_transcript_level.txt.gz` |
| gnomAD v4 RMC | **yok** (boş gövdeyle 200; `release/4.1/constraint/` yalnız transkript metrikleri) | `.../release/4.1/regional_missense_constraint` |
| gnomAD GraphQL RMC alanları | 200; yalnız GRCh37'de dolu | `https://gnomad.broadinstitute.org/api` |
| ExAC/gnomAD-v2 kanonik transkript listesi | 200, 465.936 B | `https://storage.googleapis.com/gcp-public-data--gnomad/legacy/exac_browser/canonical_transcripts.txt.gz` |
| MTR-Viewer v2.0 kalıntı tablosu | 200, 538.549.548 B, 34.701.624 satır | `https://biosig.lab.uq.edu.au/mtr-viewer/static/mtrtable_2.0_ens95.txt.gz` |
| MTR eski barındırma | 301 yönlendirme | `http://biosig.unimelb.edu.au/mtr-viewer/` |
| MTR FTP | port 21 zaman aşımı | `ftp://mtr-viewer.mdhs.unimelb.edu.au` |
| MPC düz dosyası | **404** (yalnız Hail tablosu olarak var) | `.../legacy/exac_browser/fordist_constraint_official_mpc_values_v2.txt.gz` |
| MPC, dbNSFP üzerinden | 200; KIF13A 11.903, KIF1A 11.728 varyant | `https://myvariant.info/v1/query` (`dbnsfp.mpc.score`) |
| CCR v2 BED | 200, 153.510.897 B, 8.188.409 satır | `https://s3.us-east-2.amazonaws.com/ccrs/ccrs/ccrs.autosomes.v2.20180420.bed.gz` |
| MetaDome | transkript uç noktası 200, kalıntı profili iş tabanlı (tarayıcı gerekiyor) | `https://stuart.radboudumc.nl/metadome/api/get_transcripts/KIF13A` |
| MTR3D / MTR3D-AF2 | genom çapında kalıntı tablosu yayımlanmıyor | `https://biosig.lab.uq.edu.au/mtr3d/datasets` |

Not: MTR, MPC ve CCR'nin özgün yayın künyeleri bu çalışmada **ayrıca doğrulanmadı**; yukarıda
yalnız canlı çekilen adresler ve dosyalar anılmaktadır. Canlı doğrulanan tek künye Samocha 2014'tür
(PMID 25086666, doi 10.1038/ng.3050).

`data/raw/rmc/` içeriği: `regional_missense_constraint_gnomAD_v2.1.1_transcripts_with_rmc.tsv`,
`..._without_rmc.tsv`, `exac_legacy_regional_missense_constraint.tsv`,
`rmc_regions_KIF13A_KIF1A_KIF13B.tsv` (17 hedef satır, birebir),
`gnomad_graphql_rmc_KIF13A_KIF1A_KIF13B.json` (iki yapı x üç gen),
`mtrtable_2.0_ens95.txt.gz`, `mtr_{KIF13A,KIF1A,KIF13B}_*.csv`,
`ccrs.autosomes.v2.20180420.bed.gz`, `ccr_KIF13A_KIF1A_KIF13B.bed`,
`exac_canonical_transcripts.txt.gz`, `ensembl_*_expand.json` (GRCh38),
`ensembl37_*_expand.json` (GRCh37), `region_summary_motor_vs_tail.tsv` (354 satır, uzun biçim).

---

## 9. Yeniden üretim

**Yorumlayıcı uyarısı.** Projenin venv'i (`mcp\atlas_server\.venv\Scripts\python.exe`) bu oturumda
Windows Uygulama Denetimi tarafından reddedildi: uv'nin indirdiği CPython
(`%APPDATA%\uv\python\cpython-3.12.14-windows-x86_64-none\python.exe`) engelli
("Uygulama Denetimi ilkesi bu dosyayı engelledi"), venv'in kendi `python.exe`'si de onu
başlatamıyor. Bu yüzden aşağıdaki komutlar Anaconda yorumlayıcısıyla çalıştırıldı. Script yalnız
`httpx`, `pyyaml` ve `numpy` istiyor; venv kullanılabilir hale geldiğinde ikisi de çalışır.

```powershell
cd C:\Users\pc\Claude\kif13a-neurovariant-atlas

# ana çalıştırma (tüm tablolar + figür)
C:\Users\pc\anaconda3\python.exe -u scripts\regional_constraint.py

# venv düzeldiğinde eşdeğeri
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\regional_constraint.py

# her indirmeyi tazele
C:\Users\pc\anaconda3\python.exe -u scripts\regional_constraint.py --refresh

# duyarlılık çalıştırmaları
C:\Users\pc\anaconda3\python.exe -u scripts\regional_constraint.py --raw-rates     --tag rawrates --no-figure
C:\Users\pc\anaconda3\python.exe -u scripts\regional_constraint.py --max-af 0.0001 --tag af1e4    --no-figure
C:\Users\pc\anaconda3\python.exe -u scripts\regional_constraint.py --bins 12       --tag bins12   --no-figure
C:\Users\pc\anaconda3\python.exe -u scripts\regional_constraint.py --bins 40       --tag bins40   --no-figure
```

Süre: önbellek boşken yaklaşık 12 dakika (gnomAD transkript sorguları ve ekzon başına kapsama
sorguları baskın), önbellek doluyken 1 dakikanın altında.

### Üretilen dosyalar

`data/interim/regional_constraint/`, her biri `.prov.json` yan dosyasıyla ve
`data/raw/PROVENANCE.jsonl`'a bir satırla:

| Dosya | İçerik |
|---|---|
| `regional_constraint_summary_20260916.json` | hüküm, yöntem, tüm tablolar, sağlamalar |
| `per_domain_20260916.tsv` | KIF13A domainleri + panelin motor/kuyruk kırılımı |
| `motor_domain_panel_20260916.tsv` | 11 kinezinin motor domaini, ClinVar P/LP sayılarıyla |
| `KIF13A_motor_subregions_20260916.tsv` | alt bölgeler ve yapısal sınıflar (KIF13A + KIF1A) |
| `sliding_window_100aa_20260916.tsv` | 100 kodonluk pencereler, tüm panel |
| `sliding_window_50aa_20260916.tsv` | 50 kodonluk pencereler |
| `gene_level_reproduction_20260916.tsv` | bizim sayılarımız vs gnomAD v4.1 ve v2.1.1 |
| `calibration_bins_20260916.tsv` | hız kutuları ve gözlenme olasılıkları |
| `sanity_checks_20260916.tsv` | 9 sağlama kontrolü ve sonuçları |
| `KIF13A_possible_snv_constraint_20260916.tsv` | 16.245 olası SNV: bağlam, hız, beklenti, gözlendi mi |

`data/raw/` altındaki indirmeler: `mutation_rates/samocha2014_trinucleotide_rates_20260916.txt`,
`ensembl/{GEN}_{ENST}_genomic_20260916.txt`, `ensembl/{ENST}_structure_20260916.json`,
`ensembl/{GEN}_{ENSP}_pfam_20260916.json`, `gnomad/{GEN}_{ENST}_gnomad_r4_codingsnv_20260916.tsv`,
`gnomad/{GEN}_exome_coverage_20260916.tsv`, `gnomad/regional_constraint_reference_20260916.json`,
`clinvar/{GEN}_clinvar_plp_missense_20260916.tsv`, ve 8. bölümdeki `rmc/` dosyaları.

### Bu çalışmadan çıkan yeni tuzaklar

24. **Kısıtlılık kalibrasyonuna cinsiyet kromozomu geni karıştırılmaz.** KIF4A (chrX) panele
    konduğunda hemizigot allel sayısı yüzünden sinonim o/e 0,56 çıkıyor ve ortak kalibrasyonu
    bozarak bütün otozomal genlerin beklentisini eksiltiyor. Script artık otozom dışı bir gende
    duruyor.
25. **Sinonim o/e kontrolü sabit bantla yapılmaz.** 100 kodonluk pencerede ~45 sinonim beklenir;
    Poisson gürültüsü tek başına pencerelerin üçte birini 0,80-1,25 bandının dışına atar. Doğru
    test örtüşmeyen pencerelerde ki-kare/serbestlik derecesidir.
26. **Genin 5' ucunda kalibrasyon sapar.** KIF13A'nın ilk 80 kodonunda sinonim o/e 1,486; promotör
    CpG adasının hipometilasyonu trinükleotid tablosunun ortalama hızından sapıyor. O bölgenin
    missense okuması kullanılmamalı.
27. **gnomAD kısıtlılığını yeniden üretmek için gözlenen taraf yalnız ekzomdan alınır.** Ekzom
    AC >= 1, PASS, AF <= %0,1 kuralı 11 genin 10'unda gnomAD'nin yayımladığı gözlenen missense
    sayısını birebir veriyor; joint (ekzom + genom) sayılar şişiriyor.
28. **Kendi kalibrasyonunuz panelin sinonim o/e'sini 1,0'a oturtur, gnomAD'ninki oturtmaz.** Bu tek
    başına o/e değerlerinizi gnomAD'ninkinden %3-6 yüksek gösterir. Mutlak sayı vermeden önce gen
    düzeyindeki oranı gnomAD'nin yayımlanan değerine bağlayan katsayıyla çevirin.
29. **Ensembl'in `overlap/translation/{ENSP}` uç noktası domain sınırlarını çevirinin kendi
    numaralandırmasında verir**, yani UniProt izoform kayması (KIF1A tuzağı) devreye girmez. Panel
    genlerinin motor domain sınırları için doğru kaynak budur.
30. **gnomAD bölgesel missense kısıtlılığı yalnız v2.1.1 için yayımlanmıştır ve GraphQL'de yalnız
    `reference_genome: GRCh37` ile döner**; GRCh38 ile `null` gelir. Brifte dolaşan
    `fordist_constraint_official_regional_missense_constraint_transcript_level.txt.gz` adı 404;
    gerçek dosyalar `gnomAD_v2.1.1_transcripts_with_rmc.tsv` ve `..._without_rmc.tsv`.
31. **Bir genin RMC dosyasında olmaması "veri yok" demek değildir.** 18.629 transkriptin 13.502'si
    "kırılma bulunamadı" listesindedir ve KIF13A oradadır. Bu ölçülmüş bir negatiftir.
32. **MPC bağımsız bir kanıt değildir**, RMC'nin üzerine kuruludur. KIF13A'nın hiçbir varyantının
    MPC >= 2 olmaması ile gnomAD'nin KIF13A'da bölgesel kısıtlılık bulamaması aynı olgunun iki
    yüzüdür, iki ayrı kanıt gibi sayılmamalıdır.
