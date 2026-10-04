# Devir belgesi: KIF13A NeuroVariant Atlas

Son güncelleme: 4 Eki 2026 (3.11-3.15 eklendi). Bu dosya projeyi devralan herkes (insan ya da başka
bir model) için tek giriş noktasıdır.

> **Bugünkü durum, tek paragraf:** Her hastalık hipotezi kendi testinden geçti ve hiçbiri kurtarılamadı.
> Haployetmezlik de novo ekseninde düştü; missense yalnız ATP cebinde ve KIF1A'nın dörtte biri derinlikte
> kısıtlı, popEVE KIF1A ile birebir aynı patojenik değişimlerin hiçbirini KIF13A'da ağır saymıyor; çekinik
> soru açık verilerle karara bağlanamıyor; faredeki yaşla artış insanda yok; ak madde haplotipi
> oligodendrosit ve OPC'deki KIF13A eQTL'iyle kolokalize olmuyor. Ayakta kalan: iki türde oligodendrosit
> soyu ifadesi, sığ ama gerçek bir ATP cebi kısıtlılığı, hiç ölçülmemiş baskın-negatif mekanizma ve
> örneklem dışında (NSD1/NSD3) çalışan bir yöntem. Buradaki her sayı bu depodaki bir dosyadan gelir; kaynak ve tarih `data/raw/PROVENANCE.jsonl`
içinde kayıtlıdır.

## 1. Proje ne yapıyor

KIF13A geninde hangi missense varyantların zararlı olabileceğini sıralıyor. Sorun şu: KIF13A'da klinik
veri yok, dolayısıyla tahmin araçlarının bu gende çalışıp çalışmadığı ölçülemez. Çözüm: aynı aileden
gelen ve yüzlerce kesinleşmiş hastalık varyantı olan KIF1A üzerinde araçları ölçmek, sonra ölçülmüş
eşikleri ve kalıntı düzeyindeki kanıtı KIF13A'ya taşımak.

Taşınabilen şey protein mekanizmasıdır, çünkü iki protein aynı motoru kullanır. Taşınamayan şey
fenotiptir, çünkü iki gen farklı hücrelerde çalışır. Bu ayrım yöntemin belkemiğidir.

## 2. Referans kimlikler (15 Eyl 2026'da canlı doğrulandı)

| Gen | MANE transkript | RefSeq | UniProt (MANE izoformu) | Uzunluk |
|---|---|---|---|---|
| KIF13A | ENST00000259711.11 | NM_022113.6 | Q9H1H9 | 1805 aa |
| KIF1A | ENST00000498729.9 | NM_001244008.2 | Q12756-3 | 1791 aa |

KIF1A'nın UniProt kanoniği (Q12756, 1690 aa) klinik numaralandırma **değildir**. Kanonikten MANE'ye:
1-393 aynı, 395-848 arası +9, 849-1690 arası +101. Ayrıntı: `config/targets.yaml`.

KIF13A domainleri (MANE): motor 5-352, boyun bağlayıcı 353-380, CC1 381-468, FHA 469-519, kuyruk 520-1805.

## 3. Ölçülmüş sonuçlar

### 3.1 Kısıtlılık (gnomAD v4.1)

| | LOEUF | pLI | missense Z | gözlenen/beklenen missense |
|---|---|---|---|---|
| KIF13A | 0,449 | 0,9998 | 2,79 | 0,86 (1901/2209) |
| KIF1A | 0,259 | 1,0000 | 6,87 | 0,68 (1683/2486) |

İkisi de fonksiyon kaybına dayanıksız, ama KIF13A missense değişimleri büyük ölçüde tolere ediyor.
Bu, projenin ilk varsayımına karşı en güçlü argüman.

### 3.2 Kalibrasyon (KIF1A motor domaini, 110 pozitif vs 61 popülasyon kontrolü)

Kalıntı-gruplu bootstrap, %95 güven aralığı ile:

| Model | AUROC | Aralık | %90 özgüllükte duyarlılık | Eşik |
|---|---|---|---|---|
| MetaRNN (etiketle eğitilmiş, sadece referans) | 0,949 | 0,913-0,977 | 0,85 | 0,928 |
| **AlphaMissense** | **0,925** | 0,878-0,964 | 0,79 | **0,968** |
| REVEL (etiketle eğitilmiş, sadece referans) | 0,920 | 0,868-0,959 | 0,71 | 0,855 |
| ESM1b | 0,851 | 0,774-0,917 | 0,45 | -14,4 |
| PrimateAI | 0,801 | 0,720-0,875 | 0,37 | 0,883 |
| CADD v1.7 | 0,781 | 0,696-0,858 | 0,55 | 31 |
| Evo 2 40B | 0,742 | 0,643-0,839 | 0,32 | -22,7 |
| Referans: KIF13A'da özdeş kalıntı | 0,744 | | | |
| Referans: "motor domainde mi" (tüm protein) | **0,995** | | | |

Son satır kritik: tüm protein üzerinde ClinVar benign setiyle ölçüm yapılırsa konum tek başına her modeli
geçiyor. Bu yüzden ana metrik motor domain içi ölçümdür. Eşit ağırlıklı konsensüs 0,911 ile AlphaMissense'i
geçemiyor; Evo 2 eklemek 0,887'ye düşürüyor.

Pratik sonuç: AlphaMissense'in genel eşiği 0,564 bu gende 982 varyantı zararlı sayardı, kalibre edilmiş
0,968 eşiği 295 varyant veriyor.

Ayrıntı: `docs/kalibrasyon_olcumleri.md`, setlerin kuralları `docs/KIF1A_ClinVar_kurasyonu.md`.

### 3.3 Aday listesi

`data/processed/KIF13A_candidates_20260915.tsv` (gnomAD'de görülen 2.671 missense).
Katmanlar (parantezde güvenilir gözlem): A 52 (26), B 73 (43), C 170 (110), D 33 (23).

- A: eşiği geçiyor + kalıntı KIF13A'da özdeş + KIF1A'da o kalıntıda hastalık varyantı var
- B: eşiği geçiyor + özdeş kalıntı + motor domainde
- C: eşiği geçiyor, başka yerde (çoğu kuyrukta, paralog desteği yok)
- D: eşiğin altında ama paralog desteği var (modelin yanılmış olabileceği grup)

### 3.4 Doygunluk tablosu (16 Eyl 2026, kör noktayı kapatır)

`data/processed/KIF13A_saturation_lookup_20260916.tsv`, 34.295 satır: **olası her aminoasit değişimi**.
Neden gerekliydi: KIF1A'nın patojenik varyantlarının %89'u gnomAD'de yok, yani hasta varyantları
popülasyon verisine dayalı listeye hiç girmiyor.

| | Sayı |
|---|---|
| Tek baz değişimiyle ulaşılabilir | 10.627 |
| Kalibre eşiği geçen (tümü / ulaşılabilir) | 8.038 / 1.499 |
| Eşiği geçen ulaşılabilir olanlardan gnomAD'de güvenilir taşıyıcısı olan | yalnız 176 |
| gnomAD'de hiç görülmemiş | 31.680 |
| Katman A (tümü / ulaşılabilir) | 929 / 234 |
| Katman A ulaşılabilir olanlardan **eski listede olmayan** | 182 |
| KIF1A'da **birebir aynı değişim** patojenik olan | 68 (46'sı gnomAD'de hiç görülmemiş) |
| Bunlardan katman A | 52 (35'i gnomAD'de hiç görülmemiş) |

Yani popülasyona dayalı eski liste, birebir eşleşen kanıta sahip adayların çoğunu kaçırıyordu. Örnekler:
K108T (ATP cebi, AM 0,997), C97R ve C97Y, G120V, A90D, Y94D, G83R. Hiçbiri bugüne kadar bir insanda
görülmedi; görülürse ne diyeceğimiz artık tabloda hazır.

Dokuz sağlama kontrolünün hepsi geçti; bunlar arasında mevcut 2.671 adayın katmanlarının birebir yeniden
üretilmesi (52/52 katman A dahil), fenotipli tek hasta varyantının (K1122E) artık tabloda bulunması ve
onun `c.3364A>G` gösteriminin ClinVar kaydıyla birebir tutması var (CDS'in bağımsız doğrulaması).

Yapısal sütunlar yalnız 1-380 arası pozisyonlar için doludur ve `structure_assessed` bayrağı taşır; kuyruk
için yapı değerlendirmesi yapılmış gibi okunmamalıdır. Belge: `docs/missense_sonuc.md`, figür:
`figures/KIF13A_missense_saturation_blindspot.svg`.

### 3.5 Yapı ekseni

KIF13A'nın deneysel yapısı yok. AlphaFold modeli KIF1A'nın mikrotübüle bağlı cryo-EM yapısına (8UTN)
oturtuldu: 311 kalıntı çiftinde 1,03 Å RMSD. Katman A'nın dağılımı: 23 ATP cebi, 18 mikrotübül arayüzü,
9 gömülü çekirdek, 1 yüzey, 1 boyun bağlayıcı. Yani **%79'u işlevsel bir bölgede**. Katman B tam tersi
(çoğu gömülü çekirdek), yani ayrı bir mekanizma sınıfı.

Sınır: katı cisim yerleştirme çakışma üretiyor, "temas var mı" güvenilir ama mesafe değeri değil.
Ayrıntı: `docs/yapi_ekseni.md`. Doğru KIF1A yapıları: 8UTN-8UTS yabanıl, 8UTT-8UTW P305L, 8UTY P364L,
9YA5/9YA7 R350G, 9YAB/9YAI R350W. **8UTX ve 8H3T KIF1A değildir.**

### 3.6 Hücre tipi ekseni

İnsan (CELLxGENE, beyin): KIF13A oligodendrosit soyunda en yüksek (OPC birinci), oligodendrosit soyu
ortalaması 1,49 iken nöron tiplerinde 0,16. KIF1A tam tersi: beyin sapı motor nöronu birinci.

Fare (Allen Brain Cell Atlas, `C:\Users\pc\Claude\kif13a-allenbrain`, 2,34 milyon hücre): olgun
oligodendrosit birinci, yaşlanmayla artış yalnız olgun oligodendrositte (+0,250 log2, donör düzeyinde
p=0,0002), motor nöron hipotezi negatif. O deponun kendi uyarısı geçerli: sıralama tespit oranıyla ifade
düzeyini karıştırıyor.

Sonuç: mekanizma KIF1A'dan aktarılabilir, fenotip aktarılamaz. KIF13A varyantları hastalık yapıyorsa
beklenen tablo ak madde ve miyelin ile ilgilidir, KIF1A'nın nöronal aksonal taşıma hastalığı değil.

### 3.7 İnsanda KIF13A hakkında bilinen her şey

- ClinVar'da fenotip taşıyan **tek** kayıt: VCV002570672, p.Lys1122Glu, VUS, HudsonAlpha araştırma
  gönderimi; apraksi, nöbet, konuşma apraksisi, orta düzey zihinsel yetersizlik, Chiari malformasyonu.
  AlphaMissense 0,804, kuyrukta, KIF1A'da karşılığı yok.
- VCV000442594: 253 kb delesyon, içindeki tek protein kodlayan gen KIF13A, etkilenmiş çocuk, konuşma
  ve dil gecikmesi, 2017'den beri gerekçesiz VUS.
- Lam 2021 (PMID 33526817): mozaik çerçeve kayması c.2357dupA, kuyruk bölgesi.
- Kurumsal statü: ClinGen, Gene2Phenotype, PanelApp ve OMIM fenotipinde yok; yalnız SysNDD'de "Limited".
- Düzeltme: A254V ve R305L kayıtları hasta değil, TCGA tümör splicing anotasyonlarıdır.

Ayrıntı ve vaka arama yolları: `docs/vaka_arama.md`.

### 3.8 Fonksiyon kaybı ekseni (16 Eyl 2026, `docs/fonksiyon_kaybi_ekseni.md`)

gnomAD v4.1, MANE transkripti, güvenilir gözlemli varyantlar:

| | KIF13A | KIF1A | KIF13B |
|---|---|---|---|
| Güvenilir LoF varyantı | 204 | 107 | 301 |
| NMD kaçışı (son ekzon, 50 nt kuralı) | 45 | 1 | 18 |
| **Gerçekten yıkıma uğrayabilecek (NMD-competent, LOFTEE HC)** | **158** | 102 | 283 |
| LoF allel yükünün NMD kaçış bölgesindeki payı | **%59,6** (765'te 456) | %0,1 | %5,4 |
| gözlenen/beklenen LoF | 0,372 | 0,202 | 0,572 |

KIF13A'nın yüksek görünen LoF yükü büyük ölçüde yapaydır: son ekzon proteinin %15,5'ini kodluyor ve en sık
dört "LoF" alleli NMD'den kaçıyor. En sık allel p.Trp1800Ter (AC 180, tek homozigot ondadır): 1805
aminoasidin son 5'ini kesiyor, yani işlev kaybı sayılmaz. İlk gerçek aday p.Arg1395Ter (AC 29).

**En önemli tek bulgu:** motor domainde, sağlıklı insanlarda görülen ve NMD ile yıkılabilecek 29-31 kesik
protein varyantı var (R268*, R173*, R207*, A356fs gibi, AC 2-5). Yani motor domaini erken kesmek
heterozigot halde tolere ediliyor. Bu, katman A'nın gerekçesini zayıflatır: oradaki bir missense ağır
baskın bir hastalık yapıyorsa mekanizma doz kaybı değil, baskın-negatif etki olmalıdır. Bu da ölçülmedi.

Kısıtlılık bağlamı: ClinGen'de haployetmezlik skoru 3 olan 12 nörogelişimsel gene karşı **KIF13A hepsinden
daha az kısıtlı** (LOEUF 0,449'a karşı 0,094-0,369). Ayrıca düzeltme: **KIF1A bir haployetmezlik geni
değildir** (ClinGen HI = 0, "kanıt yok"), dolayısıyla karşılaştırma ölçütü olarak kullanılamaz. Mütevazı
missense kısıtlılığı da LoF aleyhine kanıt değildir (ANKRD11 missense Z = -1,21 ama HI = 3).

İnsan ve fare kanıtı: ClinVar'da **sınıflandırılmış tek bir kesik KIF13A varyantı yok**; LoF etiketli 6
kayıt TCGA splicing gönderimleridir. 20 ClinVar CNV'sinden yalnız VCV000442594 (253 kb) tek gen olarak
KIF13A içerir; patojenik olanların hepsi 23-237 gen kapsıyor. IMPC'de 21 anlamlı fenotip çağrısının
**hepsi homozigot**; heterozigotlar yaşayabilirlik dışında hiç fenotiplenmemiş (1 test). Yani
haployetmezliğin belirleyici fare testi mevcut değil.

Hüküm: LoF, missense'ten daha iyi uyuyor ama yalnızca göreli olarak. Kanıtın desteklediği ifade "olası doz
duyarlılığı, eksik penetransla", kanıtlanmış haployetmezlik değil. Yedi çürütme testi belgenin 6.3
bölümünde.

### 3.9 De novo ekseni (16 Eyl 2026, `docs/denovo_ekseni.md`)

Doz duyarlılığını kesen test: norogelişimsel bozukluk triolarında de novo varyant yükü.

| | KIF13A | KIF1A | KIF13B | ANKRD11 | CHD8 |
|---|---|---|---|---|---|
| Kaplanis 2020 (31.058 DD triosu), gözlenen de novo LoF | **0** | 0 | 2 | 111 | 19 |
| beklenen LoF (2 × N × mu_lof) | 0,367 | 0,298 | 0,337 | 0,308 | 0,497 |
| gözlenen de novo missense | 2 | **53** | 6 | 17 | 9 |
| beklenen missense | 3,163 | 4,806 | 3,763 | 7,065 | 5,138 |
| denovoWEST p | 0,709 | 5,2e-21 | 0,286 | 1e-14 | 4,3e-16 |
| Tüm kataloglar toplamı (~97.700 etkin trio), LoF | **0** | 2 (anotasyon hatası) | 2 | 144 | 74 |
| SFARI Gene 2026 Q2 | **listede değil** | 2S | 2 | 1S | 1S |

**Hüküm: KARŞI KANIT.** Etkilenmiş hiçbir bireyde tek bir de novo LoF yok; bulunan tek de novo LoF
(p.N1649Kfs\*3) sağlıklı bir kardeşte ve son ekzonda (NMD kaçışı), yani iki ayrı nedenle doz kanıtı
değil. Boru hattı çalışıyor: ANKRD11 148 kat, CHD8 47 kat zenginleşme veriyor, KIF1A da beklenen
deseni (53 missense, 0 LoF) gösteriyor. Bu, ANKRD11/CHD8 tipi haployetmezliği **dışlıyor**, ama
%95 güvenle yalnız **2,6 katın üzerini**; hafif doz duyarlılığı hâlâ dışlanamıyor.

Yan bulgular: Kaplanis koordinatları **GRCh37**'dir (test edildi, 4/4), tablo bunu yazmaz.
Gene4Denovo2'nin `AAChange_refGene` sütunu KIF13A için **13 kalıntı düşük** üçüncü bir izoform
kullanır; doğrusu `AAChange_mrefGene` (NM_022113). **VCV000442594'ün ClinVar `Origin` alanı
"unknown"**, yani 253 kb'lik delesyonun de novo olduğu kayıtlı değil ve de novo argümanına giremez.
Motor domainde üç de novo missense var (D65E, R173Q, Q237H) ama biri ALS triosundan ve sayılar
şanstan ayrılmıyor (uzunluk oranına göre en az bir isabet olasılığı zaten 0,855).

### 3.10 Bölgesel kısıtlılık (16 Eyl 2026, `docs/bolgesel_kisitlilik.md`)

3.1'deki missense Z = 2,79 / o/e = 0,86 bir **ortalamaydı** ve 1805 aminoasidin yalnız 352'si motor
domain. Kısıtlılık yerel olarak ölçüldü (`scripts/regional_constraint.py`): Samocha 2014 trinükleotid
hızları (PMID 25086666), gnomAD v4.1 ekzomları, panelin sinonim varyantlarıyla kalibre edilmiş
beklenti, %90 Poisson aralığı.

| KIF13A bölgesi | aa | o/e missense | %90 aralık | gnomAD ölçeğinde |
|---|---|---|---|---|
| **motor** | 5-352 | **0,814** | 0,745-0,892 | 0,783 |
| CC1 | 381-468 | 0,738 | 0,611-0,897 | 0,710 |
| kuyruk | 520-1805 | 0,927 | 0,887-0,968 | 0,891 |
| motor: ATP cebi (25 kalıntı) | | **0,523** | 0,362-0,776 | 0,503 |
| motor: mikrotübül arayüzü (50) | | 0,746 | 0,581-0,968 | 0,718 |
| motor: gömülü çekirdek (124) | | 0,678 | 0,576-0,801 | 0,652 |
| motor: **yüzey** (153) | | **1,008** | 0,892-1,140 | 0,969 |

Hüküm: **zayıf kısıtlı.** Ortalama gerçekten bir sinyal gizliyordu ama küçük bir sinyal. Asıl bulgu
bir kat aşağıda: motor yüzeyi hiç kısıtlı değil, ATP cebi iki kat tükenmiş ve **katman A'nın 23
adayı tam o cepte** (18'i mikrotübül arayüzünde, 50/52'si tükenmiş sınıflarda). Tavan iki yönden
belli: (i) aynı hesap KIF1A motorunda 0,450, KIF1A'nın ATP cebinde **0,140**, yani KIF13A her
yapısal sınıfta 2-4 kat daha sığ; (ii) 11 otozomal kinezin panelinde KIF13A 7. sırada, KIF13B /
KIF16B / KIF1C ile aynı kümede. Motor domain KIF13A'nın en kısıtlı bölgesi bile değil: CC1+FHA
(381-519) 0,785.

Yayımlanmış kaynaklar (hepsi canlı çekildi, `data/raw/rmc/`): **gnomAD'nin bölgesel missense
kısıtlılık algoritması KIF13A'yı aramış ve hiçbir kırılma noktası bulamamış** (18.629 transkriptin
13.502'si gibi); MTR ve CCR, KIF13A'nın tek tükenmiş bölgesi olarak motoru değil **353-519**'u
gösteriyor; **KIF13A'nın hiçbir varyantı MPC >= 2 eşiğine ulaşmıyor** (proteinin en yüksek MPC
değeri 1,18), yani katman A herhangi bir MPC filtresinde tamamen elenir.

Sağlamalar: KIF1A pozitif kontrolü 0,450 ile geçti; gözlenen missense sayısı 11 genin 10'unda
gnomAD'nin yayımladığı sayıyla **birebir** (KIF13A 1901/1901); motor domainde sinonim o/e 0,999.
Sınırda kalan tek kontrol motor-kuyruk farkı (p = 0,036-0,060).

### 3.11 Ak madde görüntüleme ekseni (17 Eyl 2026, `docs/ak_madde_goruntuleme.md`)

Hüküm **zayıf destek**. UK Biobank BIG40'ta (33 bine kadar kişi, 3.935 fenotip) KIF13A ±100 kb içinde ak
madde sınıflarında 5e-8'i geçen bir şey yok; en iyi p 3,5e-07 (rs144963738), negatif kontrol lokusundan
yalnız 2,2 kat iyi. Gözlenen haplotip etkisi (varyansın %0,086'sı) saptanabilir sınırın (%0,126) altında;
görmek için yaklaşık 46.000 kişi gerekir. Lehte üç gözlem: haplotipteki 8 varyantın hepsi aynı ak madde
profilini veriyor; rs76091407 bağımsız bir ak madde GWAS'ında genom çapında anlamlı (UDIP-FA,
GCST90860937, p 4e-09, N 11.736); haplotip MetaBrain korteks KIF13A eQTL'inin 3. sinyaliyle r² 0,969.
rs75428660'ın genom çapında anlamlı fenotipi kortikal kıvrım ve kalınlık, ak madde değil.

### 3.12 İnsan yaşam boyu ekseni (3 Eki 2026, `docs/insan_yasam_boyu.md`)

TRUBA'da, CELLxGENE'deki açık ham sayımlarla (Nature 2026 yaşam boyu atlası PsychAD "Aging_Cohort",
SEA-AD, ak madde ve gelişim setleri), donör düzeyinde.

- **Faredeki yaşla artış insanda tekrarlanmadı.** 335 donör meta-analizi: olgun oligodendrositte yaşlı
  eksi genç -0,057 log2 (%95 GA -0,096 ile -0,017); fare +0,250. Fare büyüklüğündeki etkiyi görme gücü
  0,99'un üzerinde, yani gerçek bir negatif. Küçük düşüş KIF13A'ya özgü değil: "artmıyor" denmeli,
  "azalıyor" değil. OPC değişmiyor.
- **Oligodendrosit soyu birinciliği tekrarlandı**, ama insanda OPC birinci, olgun oligodendrosit ikinci
  (farede tersi). Astrosite karşı tespit 1,25x, düzey 2,0x.
- **Gelişim:** OPC'de KIF13A bebeklikten artıyor (üç kaynak), yarısı 3-5 yaşta tamam; aynı artış KIF13B'de
  ve astrositte de var, yani genel glial olgunlaşma.

### 3.13 Hücre tipi eQTL kolokalizasyonu (3-4 Eki 2026, `docs/eqtl_hucre_tipi.md`, ÖN SONUÇ)

Zeng 2026 (1.384 donör, Synapse), SingleBrain 2025 ve Bryois 2022 tek çekirdek eQTL'leri; BIG40 ve
UDIP-FA ak madde GWAS'ları. Boru hattı R coloc ile ve yayımlanmış kolokalizasyonlarla doğrulandı.
**Oligodendrosit ve OPC'de karşı kanıt:** KIF13A'nın bu hücrelerde güçlü eQTL'leri var (OPC öncüsü
rs760849, p 3,6e-37) ama ak madde sinyalinden farklı varyantlar; PP.H3 0,78-1,00, PP.H4 en fazla 0,015.
Tek ipucu mikroglia'da (koşullu sinyal rs73373199, PP.H4 0,55-0,62) ve LD paneli değişince kayboluyor.
Ajan belgeyi yazmadan durduruldu; belge tablolardan ana oturumda yazıldı.

### 3.14 Nakavt ekseni ve popEVE (3 Eki 2026, `docs/nakavt_ekseni.md`)

Hüküm **bilgi vermiyor**: aynı tasarım, iki kopya kaybıyla kesin hastalık yapan KIF1A'yı bile ayırt
edemiyor. Genes & Health (44.028 otozigot yetişkin) 0 nakavt, beklenen 0,17; RGC-ME (983.578) 2 homozigot,
büyük olasılıkla NMD'den kaçan son intron varyantı c.4581+1G>C (çıkarım, doğrulanmadı); Chundru 2024
(29.745 DD triosu) 0, nötr beklenti 0,0063. Nakavt sayıları Poisson değil, ampirik dağılımla
karşılaştırıldı. **popEVE:** KIF13A ne 442 genlik listede ne 123 yeni adayda; KIF1A'da patojenik olan
birebir aynı 68 değişimin **hiçbiri** KIF13A'da ağır değil (KIF1A'nın kendi pozitiflerinde %77,5).
İnsan işi: RGC-ME tarayıcısında lisans onayı verip 6-17771113-C-G ve 6-17764129-C-T satırlarındaki
homozigot sayısına bakmak.

### 3.15 AlphaGenome Atlas ve AVI skoru (4 Eki 2026, değerlendirildi, çalıştırılmadı)

DeepMind, 8 Eyl 2026; önbaskı PMID 42779988. AVI, gnomAD frekansıyla eğitilmiş (ClinVar yalnız
değerlendirmede), yani kalibrasyon setine girebilir; ama missense için protein sinyali AlphaMissense ve
gen kısıtlılığını bilmiyor, bu yüzden motor domainde büyük kazanç beklenmiyor. Asıl kullanımı splice
ekseni (26 KIF1A splice pozitifi, K1122E) ve bir kolokalizasyon başarılı olursa nedensel varyant seçimi.
Doku listesinde oligodendrosit yok. Erişim için kullanıcının API anahtarı alması gerekiyor.

## 4. Açık soru: mekanizma hangisi

> 4 Eki 2026 notu: aşağıdaki üç hipotezin durumu 3.11-3.14 ile güncellendi; özet en üstteki paragrafta.
> Baskın-negatif hâlâ tek tutarlı mekanizma ve hâlâ ölçülmedi. Dürüst yön, KIF1A/KIF13A'yı çoğunlukla
> negatif bir vaka çalışması olarak kullanan çok çiftli yöntem makalesi.

Proje motor domain missense üzerine kuruldu. 16 Eyl 2026 ölçümlerinden sonra durum şu:

- **Missense hipotezi** ayakta ama zayıf, ve artık **yerel olarak ölçüldü** (3.10): gen düzeyindeki
  0,86 gerçekten bir ortalamaydı, ama altından çıkan sinyal küçük. Motor domain o/e 0,814
  (%90 aralık 0,745-0,892, gnomAD ölçeğinde 0,78) ve kuyruk 0,927; fark sınırda (p = 0,036-0,060).
  İşe yarayan kırılım yapısaldır: motor **yüzeyi hiç kısıtlı değil** (1,008), **ATP cebi 0,523** ve
  katman A'nın 23 adayı tam o ceptedir. Tavan iki ölçümle belli: aynı hesap KIF1A'nın ATP cebinde
  **0,140** veriyor ve KIF13A 11 kinezinlik panelde 7. sırada. Ayrıca gnomAD'nin bölgesel kısıtlılık
  algoritması KIF13A'yı arayıp kırılma noktası **bulamamış**, MTR ve CCR de tükenmiş bölge olarak
  motoru değil 353-519'u gösteriyor. ClinVar'da hâlâ tek bir patojenik missense yok. Desteklenen
  ifade: "işlevsel çekirdekte gerçek ama KIF1A'nınkinin dörtte biri kadar derin bir seçilim izi".
- **Haployetmezlik hipotezi** artık **büyük ölçüde elendi**: LoF yükünün çoğu NMD kaçışı artefaktı,
  gen ClinGen'in haployetmezlik genlerinden daha az kısıtlı, sınıflandırılmış kesik varyant yok,
  farede heterozigot hiç test edilmemiş, ve (3.9) **yaklaşık 97.700 triodan hiçbirinde etkilenmiş bir
  bireyde de novo LoF yok**. De novo ekseni bunun üstüne nicel bir sınır koyuyor: doz duyarlılığı
  kalacaksa %95 güvenle **2,6 kattan küçük** olmak zorunda, yani ClinGen haployetmezlik genlerinde
  ölçtüğümüz 47-148 kat aralığından iki büyüklük mertebesi uzakta. Desteklenen en güçlü ifade
  "hafif veya eksik penetranslı doz duyarlılığı olabilir", kanıtlanmış haployetmezlik değil.
- **Baskın-negatif hipotezi** şu an katman A'nın tek tutarlı gerekçesi ve de novo ekseninden sonra
  **göreli olarak güçlendi**: motor domainde erken kesilme tolere ediliyor, de novo LoF hiç yok, ama
  tam boy hatalı protein zehirleyici olabilir. Hiç ölçülmedi. KIF1A'da bu mekanizma gösterilmiştir
  (Anazawa 2022, PMID 35917346), dolayısıyla test edilebilir bir hipotezdir. Ayrıca de novo ekseni
  bu hipoteze sınanabilir üç gerçek insan gözlemi verdi: motor domainde **D65E, R173Q, Q237H**
  (biri ALS triosundan), doygunluk tablosuyla kesiştirilmeyi bekliyor.

Üç iş paralel yürütüldü:

| İş | Çıktı | Durum |
|---|---|---|
| Fonksiyon kaybı ekseni | `scripts/lof_axis.py`, `data/interim/lof/`, `docs/fonksiyon_kaybi_ekseni.md` | bitti |
| De novo ekseni | `scripts/denovo_axis.py`, `data/interim/denovo/`, `docs/denovo_ekseni.md` | bitti |
| Bölgesel kısıtlılık | `scripts/regional_constraint.py`, `data/interim/regional_constraint/`, `data/raw/rmc/`, `docs/bolgesel_kisitlilik.md`, `figures/KIF13A_bolgesel_kisitlilik.svg` | bitti |
| Missense sonucu ve doygunluk tablosu | `scripts/build_saturation_table.py`, tablo, `docs/missense_sonuc.md` | bitti |
| Çok çiftli ölçekleme fizibilitesi | `scripts/paralog_pairs_scan.py`, `docs/olcekleme_fizibilite.md` | bitti |

### 4.1 Ölçekleme taraması sonucu (16 Eyl 2026)

ClinVar'da en az 30 patojenik missense taşıyan 445 gen bulundu; kardeşi klinik olarak boş olan (en fazla 2
patojenik missense) 3.198 çift puanlandı.

| | Sonuç |
|---|---|
| Dört zayıflık filtresini geçen çift | 575 (%18), 183 zengin gen |
| Aktarım kazancının (lift) medyanı | 0,048 |
| Kazancı sıfır ya da negatif olan çift | %35,5 |
| Kardeşi hiç kısıtlı olmayan çift | %54 |
| **KIF1A/KIF13A sırası** | **152 / 3.198** (tekilleştirilmiş 71 / 337), sıfır zayıflık bayrağı |

En iyi beş çift: NSD1/NSD3, COMP/THBS1, PRKCG/PRKCE, MTM1/MTMR4, PCDH19/PCDH10.

Dürüst hüküm: **paralog aktarımı genel bir özellik değil**, yaklaşık altı çiftten birinde işe yarıyor.
Dolayısıyla yayınlanabilir iddia "paralog aktarımı çalışır" değil, "ucuz ölçütlerle hangi çiftte çalışacağı
önceden kestirilebilir" olmalı. KIF1A/KIF13A üst çeyrekte ve hiç zayıflık bayrağı yok, yani vaka
çalışmamız şanslı bir seçim değil. Kazancının (0,128) düşük görünmesinin nedeni KIF13A'nın motor
domaininin zaten %62-70 özdeş olması: yüksek arka plan, kazancı bastırıyor ama aktarımı kolaylaştırıyor.

Sonraki çift olarak NSD1/NSD3 ve PCDH19/PCDH10 öneriliyor (ikisi de nörogelişimsel, sıcak noktaları dar,
kardeşleri kısıtlı); her biri yaklaşık bir hafta ve süre çoğunlukla elle ClinVar kürasyonuna gidiyor.
RYR1/RYR3 ölçütleri en iyi olan çift ama 5.038 aminoasitlik proteinler yüzünden uygulanabilir değil.

### 4.2 NSD1/NSD3 doğrulaması (17 Eyl 2026, `docs/nsd1_nsd3_dogrulama.md`)

Taramanın 1. sıradaki çifti (fizibilite 0,931, sıfır zayıflık bayrağı) tahmin satırı kapalı tutularak
baştan sona çalıştırıldı: `scripts/nsd_pair_calibration.py` (10 aşama), veri `data/interim/nsd_pair/`,
aday listesi `data/processed/NSD3_candidates_20260917.tsv`. Kimlikler `config/targets.yaml` içinde
(NSD1 ENST00000439151.7 / NM_022455.5 / Q96L73, 2696 aa; NSD3 ENST00000317025.13 / NM_023034.2 /
Q9BZ95, 1437 aa; **ikisinde de izoform kayması yok**).

| | NSD1/NSD3 | KIF1A/KIF13A |
|---|---|---|
| taramadaki sıra | 1 / 3.198 | 152 / 3.198 |
| öngörülen kazanç (lift) | 0,250 | 0,128 |
| paylaşılan domain | PHD1-CysHis modülü (1539-2210), özdeşlik 0,646 | motor domain, 0,62-0,70 |
| çekirdek pozitif (domain içi) | 120 (82 kalıntı) | 110 (70 kalıntı) |
| popülasyon negatifi AC >= 3 | 72 | 61 |
| **AlphaMissense AUROC** | **0,973** (0,942-0,995) | 0,925 (0,878-0,964) |
| %90 özgüllükte duyarlılık / eşik | 0,95 / **0,9732** | 0,79 / 0,968 |
| CADD / ESM1b / PrimateAI | 0,820 / 0,902 / 0,857 | 0,781 / 0,851 / 0,801 |
| paralog-özdeşlik referansı | 0,740 | 0,744 |
| "domain içinde mi" (tüm protein) | 0,964 | 0,995 |
| öğretmen / öğrenci domain o/e missense | 0,433 / 0,760 | 0,450 / 0,814 |
| genel eşikle / kalibre eşikle çağrılan | 399 / 101 (1.336 missense) | 982 / 295 (2.671) |
| katman A/B/C/D | 23 / 45 / 33 / 30 | 52 / 73 / 170 / 33 |

**Hüküm: öngörü tuttu**, iki uyarıyla. (i) Taramanın sıralamada kullandığı kazanç ölçütü, ölçülen
paralog katkısını kestirmedi: kazanç iki kat yüksek ama özdeşlik referans çizgisi ayırt edilemez
(0,740 ve 0,744). (ii) Kardeşin kısıtlılığı taramada gen düzeyinde ölçülüyor; bölgesel hesap NSD3'ün
paylaşılan modülünün proteinin geri kalanından daha kısıtlı **olmadığını** gösteriyor (o/e missense
0,760'a karşı modül dışı 0,774; NSD1'de modül 0,433, modül dışı 1,032), yani aktarılan eşiğin gerçek
bir hastalık mekanizmasına oturduğu gösterilmedi. Kestirim iddiasının asıl sınavı da yapılmadı: iki
çift de taramanın üst %5'inde, "başarısız olacak" denen bir çift hiç çalıştırılmadı. Sonraki deney o
yanlış-pozitif kontrolü olmalıdır.

## 5. Dosya haritası (okuma sırası)

1. `HANDOFF.md` (bu dosya) ve `CLAUDE.md` (çalışma kuralları)
2. `config/targets.yaml` kimlikler ve izoform tuzağı
3. `docs/KIF1A_ClinVar_kurasyonu.md` kalibrasyon setlerinin kuralları
4. `docs/kalibrasyon_olcumleri.md` model ölçümleri
5. `docs/yapi_ekseni.md`, `docs/fonksiyon_kaybi_ekseni.md`, `docs/denovo_ekseni.md`,
   `docs/bolgesel_kisitlilik.md`, `docs/vaka_arama.md`
5b. `docs/olcekleme_fizibilite.md` (çok çiftli tarama) ve `docs/nsd1_nsd3_dogrulama.md`
   (yöntemin örneklem dışı doğrulaması, ikinci vaka çalışması)
5c. `docs/ak_madde_goruntuleme.md`, `docs/insan_yasam_boyu.md`, `docs/eqtl_hucre_tipi.md`,
   `docs/nakavt_ekseni.md`. Son üçü TRUBA'da çalıştı; script'leri `scripts/truba/`, küçük sonuçları
   `data/interim/{insan_yas,eqtl_hucre_tipi,nakavt}/`, tam çalışma alanları
   `/arf/scratch/avural/kif13a_{insan_yas,eqtl,nakavt}` (scratch dönemsel silinir)
5d. `skill/kif13a-atlas/`: Claude Code skill'inin kopyası (yöntem, durum, 77 tuzak, sonraki adımlar);
   kendi makinende kullanmak için `~/.claude/skills/` altına kopyala
6. `data/processed/` aday listesi ve doygunluk tablosu
7. `progress.md` oturum günlüğü, `data/raw/PROVENANCE.jsonl` her indirmenin kaydı

## 6. Yeniden üretim

Python ortamı: `mcp\atlas_server\.venv\Scripts\python.exe` (Microsoft Store `python`'u kullanma,
`%APPDATA%\uv` klasörünü göremez).

```powershell
mcp\atlas_server\.venv\Scripts\python.exe scripts\pull_wp1_core.py
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\kif1a_clinvar_curation.py
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\kif1a_curation_extras.py
mcp\atlas_server\.venv\Scripts\python.exe scripts\kif1a_build_calibration_sets.py
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\evaluate_calibration_sets.py
mcp\atlas_server\.venv\Scripts\python.exe scripts\build_kif13a_candidates.py
mcp\atlas_server\.venv\Scripts\python.exe scripts\build_saturation_table.py
mcp\atlas_server\.venv\Scripts\python.exe scripts\structure_axis.py
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\lof_axis.py
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\denovo_axis.py
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\regional_constraint.py
```

> **16 Eyl 2026 uyarısı:** venv yorumlayıcısı bu makinede artık çalışmıyor; Windows Smart App Control
> uv'nin yönettiği `%APPDATA%\uv\python\cpython-3.12.14-windows-x86_64-none\python.exe` dosyasını
> blokluyor ("Uygulama Denetimi ilkesi bu dosyayı engelledi"). De novo ekseni bu yüzden sistem
> Python 3.13 ile üretildi (`python -u scripts\denovo_axis.py`); `atlas_mcp.core` orada sorunsuz
> içe aktarılıyor ve provenans yazımı normal çalışıyor. Bölgesel kısıtlılık ekseni aynı nedenle
> Anaconda ile üretildi (`C:\Users\pc\anaconda3\python.exe -u scripts\regional_constraint.py`).
> Gereken tek şey `httpx`, `pyyaml` ve `numpy` taşıyan bir Python 3.12+.

Veri erişimi MCP sunucuları üzerinden kuruldu (`docs/MCP_KURULUM.md`), ama yukarıdaki script'ler doğrudan
HTTP kullanır; başka bir model ya da ortam MCP olmadan da çalıştırabilir. Anahtarlar `.env` içinde
(git dışı): `NCBI_API_KEY`, `NVIDIA_API_KEY`.

## 7. Devralanın bilmesi gereken tuzaklar

1. Protein uzunluğu 1805, 1770 değil (o başka bir izoform).
2. KIF1A klinik numaralandırması UniProt kanoniğinden kayar; kuyruk varyantları yanlış kalıntıya düşer.
3. gnomAD'nin gen sorgusu MANE'ye göre numaralamaz; transkript sorgusu kullanılmalı.
4. gnomAD listesi yalnız ekzon ± 75 bazı kapsar, derin intron yoktur.
5. AC0 filtreli 596 varyantın güvenilir taşıyıcısı yoktur.
6. ClinVar'ın sınıflandırma alanı `germline_classification`; eski alan boş döner.
7. Bir karşılaştırmanın iki sınıfı asla farklı skor kaynağından beslenmemeli (CADD 0,914'e karşı 0,781).
8. Negatifler pozitiflerle aynı bölgeden seçilmeli, yoksa konum tek başına kazanır (0,995).
9. Sınıflandırması in silico kanıta dayanan varyantlar ayrı tutulmalı; REVEL ve MetaRNN etiketle
   eğitildiği için konsensüse giremez.
10. Bootstrap ve çapraz doğrulama kalıntı bazlı olmalı, varyant bazlı değil.
11. VEP sıralanmamış VCF girdisini reddeder; dbNSFP alanları transkript listesi olarak gelir.
12. Dosya globları `.prov.json` yan dosyalarını da yakalar.
13. Evo 2 barındırılan servisi tam nedensel değildir; skor varyant bazından itibaren toplanmalı ve
    mutlak değeri 3'ün altındaki skorlar gürültü sayılmalı.
14. **ClinVar Entrez sorgusunda `pathogenic[clinsig]` diye bir alan yoktur.** Entrez bunu sessizce
    `pathogenic[All Fields]` olarak çevirir ve şişmiş sayı döndürür (KIF1A için 692). Doğrusu
    `("clinsig pathogenic"[Properties] OR "clinsig likely pathogenic"[Properties])`; tek gen ve missense
    filtreleriyle 135 verir ve bu bizim kürasyonumuzla birebir tutar.
15. Fonksiyon kaybı sayılırken NMD kaçışı ayrılmalı: KIF13A'nın LoF allel yükünün %60'ı son ekzondaki
    kaçan varyantlardan gelir ve bunlar işlev kaybı sayılmaz.
16. **Kaplanis 2020'nin ek tablosu GRCh37'dir ve bunu hiçbir yerde yazmaz.** GRCh38 sanılırsa dört
    KIF13A satırının da sonuç sınıfı değişir. Varsayma, test et (her satırı iki yapıda anotla).
17. **Gene4Denovo2'nin `AAChange_refGene` sütunu KIF13A için NM_001105567/8 kullanır ve MANE'den 13
    kalıntı düşük numara verir.** Doğrusu `AAChange_mrefGene` (NM_022113). Bu, 1770 tuzağından
    farklı, **üçüncü** bir izoform kayması.
18. **denovo-db her varyant için her transkriptten bir satır yazar**; satır saymak sayıyı dört katına
    kadar şişirir (KIF13A'nın "101 satırı" 25 varyanttır).
19. **Gene4Denovo2'nin LoF etiketi tek başına güvenilmez**: KIF1A'da 3 bazlık bir çerçeve içi
    insersiyon `stopgain` etiketli; Kaplanis'in kendi tablosunda `inframe_insertion`.
20. **`1 - CDF` formundaki Poisson p değeri pozitif kontrollerde taşar** ve sahte `p = 0` basar;
    üst kuyruk log uzayında toplanmalı.
21. **VCV000442594'ün `Origin` alanı "unknown"**; 253 kb'lik delesyon de novo kanıtı olarak
    kullanılamaz.
22. Sıfır gözlem sıfır etki değildir: KIF13A'da beklenen de novo LoF zaten 0,37 (Kaplanis) ve 1,16
    (toplam), dolayısıyla yalnız 8,2x ve 2,6x üzeri dışlanır.
23. Proje venv'i (`mcp\atlas_server\.venv\Scripts\python.exe`) bu makinede **çalışmıyor**: Windows
    Smart App Control uv'nin yönettiği yorumlayıcıyı blokluyor. Sistem Python 3.13 `atlas_mcp.core`
    dahil her şeyi çalıştırıyor.
24. **NSD3'ün eski sembolü WHSC1L1'dir ve Ensembl `lookup/symbol` bu sembole 400 döner.** Ayrıca
    NSD3'ün literatürde çok çalışılan kısa izoformu Q9BZ95-3 (NM_060248, 645 aa) SET domaini hiç
    taşımaz; o izoform üzerinden verilmiş bir kalıntı numarası MANE ile karşılaştırılamaz.
25. `pathogenic[clinsig]` şişmesi NSD1'de bağımsız olarak yeniden üretildi: doğru Properties biçimi
    158, yanlış biçim 371. Yani tuzak 14 gene özgü değil, sorgu biçimine özgüdür.
26. Ensembl REST 17 Eyl 2026'da aralıklı olarak 500 ve 503 döndürdü; `core.request`'in beş denemesi
    yetmedi, Ensembl çağrıları dokuz denemeye çıkarıldı (`nsd_pair_calibration.py::ens_get`).

Uzun liste ve nedenleri: `~/.claude/skills/kif13a-atlas/references/tuzaklar.md` (Claude skill'i) ve
ilgili docs dosyaları.

## 8. Yapılmamış işler

Splice ekseni (26 patojenik splice varyantı hiç skorlanmadı), ağırlıklı konsensüs, vaka arama başvuruları
(insan hesabı gerekiyor), mikroglia kolokalizasyon ipucunun örnek içi LD ile testi, baskın-negatif
testinin tasarımı, taramanın "başarısız olur" dediği bir çiftin denenmesi, Nextflow/DVC paketlemesi ve
Zenodo DOI. Depo 4 Eki 2026'da ilk kez commit edilip GitHub'a gönderildi; `data/raw` (831 MB) provenans
kaydından yeniden üretilebildiği için dışarıda.
