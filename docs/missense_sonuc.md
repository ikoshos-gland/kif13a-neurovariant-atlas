# Missense ekseni: sonuçlar

Tarih: 16 Eyl 2026.
Üreticiler: `scripts/build_kif13a_candidates.py` (gnomAD'de görülen 2.671 missense), `scripts/structure_axis.py`
(yapı ekseni), `scripts/build_saturation_table.py` (bütün değişim uzayı).
Çıktılar: `data/processed/KIF13A_candidates_20260915.tsv`,
`data/processed/KIF13A_saturation_lookup_20260916.tsv` (34.295 satır),
`data/processed/KIF13A_saturation_summary_20260916.json`,
`figures/KIF13A_missense_saturation_blindspot.svg`.
Kalibrasyon: `docs/kalibrasyon_olcumleri.md`. Kürasyon: `docs/KIF1A_ClinVar_kurasyonu.md`.
Yapı: `docs/yapi_ekseni.md`. İnsan kayıtları: `docs/vaka_arama.md`.

## 0. Bu belgenin iddiası ve iddia etmediği şey

**İddia etmediğimiz şey:** "aşağıdaki varyantlar hastalığa yol açar". Bugün bunu söyleyecek tek bir insan
gözlemi bile yok.

**İddia ettiğimiz şey:** *eğer* bir KIF13A missense bozukluğu varsa, en iyi adaylar bunlardır; ve bu
"eğer"in kendisinin lehine ve aleyhine olan kanıt aşağıdadır. Öncülün aleyhindeki kanıt (bölüm 3) en az
adayların kendisi kadar yer tutuyor, çünkü şu an daha ağır basan taraf o.

Belgenin üçüncü işi, listeyi üretirken yaptığımız bir yöntem hatasını düzeltmek ve hatanın ne kadar
büyük olduğunu göstermek (bölüm 2).

## 1. Üç kanıt ekseni ve nerede buluşuyorlar

Üç eksen de birbirinden bağımsız girdilere dayanıyor. Hiçbiri KIF13A klinik etiketiyle eğitilmedi, zaten
öyle bir etiket yok.

| Eksen | Girdi | KIF1A'da ölçülen gücü | KIF13A'ya aktarım |
|---|---|---|---|
| Model skoru | AlphaMissense, Q9H1H9 kanonik (= MANE) | motor içi AUROC **0,925** (%95 GA 0,878-0,964), sens@90%spec 0,79 | gen-kalibre eşik **0,968** (genel eşik 0,564) |
| Paralog aktarımı | KIF1A-KIF13A çift hizalaması + kürate KIF1A P/LP seti | tek başına AUROC **0,744** | kürate 111 pozitifin 98'i KIF13A'da **özdeş** kalıntıya düşüyor (70 kalıntının 60'ı) |
| Yapı | KIF13A AlphaFold modeli, 8UTN KIF1A-mikrotübül kompleksine bindirildi (budanmış RMSD 1,03 Å) | modelden bağımsız, sağlama kontrolleriyle doğrulandı (T104, R258, P303, R348) | ATP cebi / MT arayüzü / gömülü çekirdek / yüzey sınıfı, yalnız 1-380 |

**Buluşma noktası ölçülebilir bir şey.** Model ve paralog ekseni bağımsız olarak seçilen katman A
adaylarının (eski listede 52 varyant) **41'i (%79) ya ATP cebinde ya da mikrotübül arayüzünde**
(`docs/yapi_ekseni.md`, bölüm 4). Yüzeyde duran tek aday M214R ve onun pLDDT'si 42, yani zaten
güvenilmez bölge. Bu, üçüncü eksenin ilk ikisini doğruladığı anlamına geliyor: eğer katman A rastgele
yüksek skorlu varyantlardan oluşsaydı, yapısal dağılımının motor alanın genel dağılımını (ATP cebi 25,
arayüz 50, gömülü 124, yüzey 149) taklit etmesi gerekirdi, etmiyor.

Kalibrasyonun bize öğrettiği iki şey de burada duruyor:

- **Genel eşik kullanılamaz.** AlphaMissense'in kendi "likely pathogenic" sınırı 0,564. Saturasyon
  tablosunda bu sınırı **18.459** değişim geçiyor, gen-kalibre 0,968 sınırını ise **8.038**. Aradaki
  10.421 değişim, KIF1A'da %90 özgüllüğü tutturmak için feda edilmesi gereken yanlış pozitiflerdir.
- **Konsensüs yardım etmiyor.** Etiketle eğitilmemiş dört modelin sıra ortalaması AUROC 0,911, tek
  başına AlphaMissense 0,925. Evo 2 40B bu missense setinde 0,742 ve karıştırınca performansı düşürüyor.
  Bu yüzden missense ekseninin skoru tek modeldir, konsensüs değil.

## 2. Kapatılan kör nokta: gnomAD bir hastalık kataloğu değil

### 2.1 Hata

Eski aday listesi gnomAD v4.1'in KIF13A'da gördüğü **2.671 missense varyantla** başlıyordu. Yani liste
şu soruyu soruyordu: "popülasyonda görülmüş değişimlerden hangileri zararlı görünüyor?" Doğru soru ise
şu: "olası bütün değişimlerden hangileri zararlı görünüyor?"

Fark akademik değil, ölçtük:

- Kürate **111 KIF1A patojenik missense varyantının 99'u (%89) gnomAD'de hiç yok.** Baskın, de novo,
  ağır bir pediatrik bozukluğun alleli için beklenen tam da budur.
- KIF13A'nın fenotipi kayıtlı **tek** hasta varyantı olan p.Lys1122Glu (VCV002570672, HudsonAlpha;
  apraksi, konuşma apraksisi, nöbet, orta düzey zihinsel yetersizlik, Chiari malformasyonu) da gnomAD'de
  yok. Bu yüzden 2.671'lik listeye **hiç girmedi**.

Yani listenin kuruluşu, aradığımız şeyin türünü yapısal olarak dışarıda bırakıyordu.

### 2.2 Düzeltme

`scripts/build_saturation_table.py` bütün uzayı kuruyor: **1805 aa x 19 = 34.295 satır**, her satır bir
olası aminoasit değişimi. gnomAD artık bir filtre değil, sadece bir sütun.

Satır başına: AlphaMissense skoru ve sınıfı, iki eşiği geçip geçmediği; MANE CDS'ten (ENST00000259711,
bir kez indirilip önbelleğe alındı, 5.418 nt, 1805 aa'ya çevrildiği doğrulandı) tek nükleotid değişimiyle
erişilebilirlik, kodon değişimi ve `c.` gösterimi; gnomAD durumu (güvenilir taşıyıcı / AC0 elenmiş / hiç
görülmemiş); varsa ClinVar kaydı; paralog kanıtı (denk KIF1A kalıntısı, özdeşlik, o kalıntıdaki KIF1A
patojenik ve tolere edilen varyantları, **aynı aminoasidin KIF1A'da patojenik olup olmadığı**); 1-380
için yapı ekseni sütunları; ve mevcut betikten değiştirilmeden alınan katman kuralı.

### 2.3 Uzayın şekli

| | n | gnomAD'de güvenilir taşıyıcı | gnomAD'de yok / AC0 |
|---|---|---|---|
| Bütün olası değişimler | 34.295 | 2.045 (%6) | 32.250 (%94) |
| Tek nükleotidle erişilebilir | 10.627 | 2.045 (%19) | 8.582 (%81) |
| AlphaMissense ≥ 0,968 ve erişilebilir | **1.499** | 176 (%12) | 1.323 (%88) |
| Katman A ve erişilebilir | **234** | 26 (%11) | 208 (%89) |
| Katman A + birebir paralog eşi | **52** | 12 (%23) | 40 (%77) |
| Referans: KIF1A kürate patojenik missense | 111 | 12 (%11) | 99 (%89) |

Şekil: `figures/KIF13A_missense_saturation_blindspot.svg`.

Katmanların tamamı (erişilebilirlik koşulu olmadan): A 929, B 1.944, C 5.165, D 211, geri kalan 26.046.
Erişilebilir olanlar: A 234, B 403, C 862, D 110.

**Son satırdaki oranın, kendi adaylarımızdaki oranla neredeyse aynı olması** (%89 ve %89) bu belgedeki en
öğretici sayıdır: KIF1A'nın gerçek hastalık varyantları gnomAD'de ne kadar görünmezse, bizim en iyi
adaylarımız da o kadar görünmez. Eski liste tam olarak bu görünmez kısmı atıyordu.

### 2.4 Kör noktanın büyüklüğü

- Katman A'daki 234 erişilebilir değişimin **52'si** eski listede vardı, **182'si yoktu**.
- KIF1A'da **birebir aynı** aminoasit değişiminin patojenik olduğu, tek nükleotidle erişilebilir 68
  KIF13A değişimi var (52'si katman A). Bunların yalnız **15'i** gnomAD'de güvenilir biçimde görülmüş.
  **35'i gnomAD'de hiç yok** ve eski listeye hiç girmemişti.
- ClinVar'daki 232 KIF13A missense kaydının **41'i** gnomAD'de yok; K1122E bunlardan biri.

### 2.5 Sağlama kontrolleri (hepsi geçti)

| Kontrol | Sonuç |
|---|---|
| Satır sayısı 1805 x 19 | 34.295 |
| AlphaMissense referans kalıntısı MANE dizisiyle uyuşuyor mu | 34.295/34.295, ofset yok |
| CDS 1805 aa'ya çevriliyor ve stop ile bitiyor mu | evet, 5.418 nt, TGA |
| p.Lys1122Glu tabloda, AlphaMissense 0,8044 | evet |
| p.Lys1122Glu paralog eşlemesi yok | evet (kuyruk, hizalama dışı) |
| p.Lys1122Glu gnomAD'de yok, ClinVar VCV002570672 taşıyor | evet |
| p.Lys1122Glu tek nükleotidle erişilebilir, c.3364A>G | evet (ClinVar kaydındaki `c.` ile birebir aynı, CDS'in bağımsız doğrulaması) |
| p.Lys1122Glu yapısal olarak değerlendirilmemiş | evet, `structure_assessed` = False |
| Eski listedeki 2.671 satırın katmanı yeniden üretiliyor mu | 2.671/2.671, 0 uyuşmazlık |
| Eski katman A'nın 52 üyesi katman A'da kalıyor mu | 52/52 |
| gnomAD'deki 2.615 missense protein değişiminin hepsi tabloda mı | 2.615/2.615 |

Katman kuralı yeniden yazılmadı: eşikler, alan sınırları ve kural `build_kif13a_candidates.py`'den
içe aktarıldı, ve yukarıdaki iki satır kuralın aynı sonucu verdiğini kanıtlıyor.

## 3. Öncülün aleyhindeki argüman

Bu bölüm listeyi zayıflatmak için değil, listenin ne kadarına güvenilebileceğini söylemek için var.

### 3.1 KIF13A missense'i KIF1A'dan çok daha iyi tolere ediyor

gnomAD v4.1 kısıtlama (`data/raw/gnomad/*_constraint_gnomad_r4_20260915.json`):

| Ölçüm | KIF13A | KIF1A | Okunuşu |
|---|---|---|---|
| missense Z | **2,79** | 6,87 | KIF13A'da missense kısıtlaması mütevazı |
| gözlenen/beklenen missense | **0,86** (GA 0,83-0,89) | 0,68 | beklenen missense'in %86'sı gerçekten görülüyor |
| gözlenen/beklenen missense, ham | 1.901 / 2.209 | | eksik olan yalnız ~308 varyant |
| LOEUF | 0,449 | | LoF'a karşı **güçlü** kısıtlama |
| pLI | 0,9998 | | aynı |
| gözlenen/beklenen LoF, ham | 79 / 212,2 | | beklenen LoF'un %63'ü eksik |

Bu profil tek bir şey söylüyor: **KIF13A'nın seçilim baskısı LoF üzerinde, missense üzerinde değil.**
KIF1A'da durum tersidir ve KIF1A'nın hastalığı da gerçekten missense ağırlıklıdır. Yani paralogdan
aktardığımız şey, paraloğun kendi hastalık mekanizmasıdır ve o mekanizmanın KIF13A'da geçerli olduğunu
gösteren bağımsız bir kanıtımız yok.

### 3.2 Bugüne kadarki bütün insan gözlemleri LoF biçimli ya da motor alanın dışında

| Gözlem | Ne | Nerede | Mekanizma |
|---|---|---|---|
| Lam ve ark. 2021, PMID 33526817 | `NM_022113.6:c.2357dupA`, de novo postzigotik mozaik çerçeve kayması, NMD öngörülüyor | ~786. kodon, **motor alanın dışında** | **LoF** |
| VCV000442594 | yaklaşık 253 kb heterozigot delesyon, tek protein kodlayan geni KIF13A, etkilenmiş çocuk, konuşma ve dil gecikmesi | tüm gen | **LoF / haployetmezlik** |
| VCV002570672 | p.Lys1122Glu, VUS, araştırma gönderimi, fenotipli tek missense | **1122, kuyruk**, motor alanın dışında | bilinmiyor |
| 6p22.3 delesyon sendromu, PMID 34303382 | 2,5 Mb de novo, 13 gen | bitişik gen | **LoF**, KIF13A'ya atfedilemez |

Üçünün de ortak noktası: hiçbiri motor alan missense'i değil. SysNDD KIF13A'yı yalnız "Limited" olarak
ve yalnız Lam 2021'e dayanarak kaydediyor.

### 3.3 ClinVar'da KIF13A missense için sıfır patojenik etiket

232 missense kaydın 215'i VUS (tek gönderen, 1 yıldız), 6'sı B/LB, 11'inde germline sınıflandırması yok.
**P/LP sayısı sıfır.** Yani listeyi doğrulayabilecek bir dış etiket kümesi yok; bu eksenin değerlendirmesi
KIF1A üzerinden dolaylı yapılmak zorunda ve bu, ölçtüğümüz 0,925 AUROC'un KIF13A'ya taşındığı varsayımına
dayanıyor.

### 3.4 Dürüst özet

Kanıt bugün **haployetmezlik** hipotezini destekliyor, missense hipotezini değil. Missense hipotezi hâlâ
ayakta çünkü hiç test edilmedi: KIF13A motor alanında de novo bir missense taşıyan fenotipli bir hasta
literatürde yok, ClinVar'da yok. Yokluk, yokluğun kanıtı değil, ama "henüz kimse bakmadı" ile "baktılar,
bulunmadı" arasındaki farkı da bir kanıt olarak kullanamayız.

Buna karşılık bir şey var: eğer KIF13A'da missense bozukluğu olsaydı, ilk göreceğimiz varyantların
**ATP cebinde ve mikrotübül arayüzünde**, gnomAD'de **görülmemiş** ve **KIF1A'da birebir eşi patojenik**
olması beklenirdi. Bölüm 5.2'deki 40 satır tam olarak bu tanıma uyuyor ve bugüne kadar hiçbir listede
yoktu.

## 4. Hipotezi ne doğrular, ne öldürür

Bu bölümün amacı, hipotezi çürütülebilir kılmak.

### 4.1 Doğrulayacak gözlem

**En güçlü tek gözlem:** nörogelişimsel fenotipi olan bir probandda, KIF13A motor alanında (5-352)
**de novo** bir missense, ebeveyn örnekleriyle doğrulanmış, ve tercihen bölüm 5'teki kalıntılardan
birinde. Bu tek vaka ACMG'de PS2 (de novo) + PM1 (fonksiyonel alan) + PM2 (popülasyonda yok) + PP3
(gen-kalibre in silico) verir ve eksenin bütün öngörüsünü doğrular.

**İkinci en güçlü:** aynı kalıntıda, ilişkisiz iki probandda bağımsız de novo missense (PS4 yönünde).
R348 (KIF1A R350) bunun için en iyi adres: iki farklı değişimi (R348G, R348Q) zaten gnomAD'de düşük
frekansta duruyor, arayüzde, pLDDT 93 ve KIF1A'da R350G/R350W'nin 2026 yapıları var.

**Fonksiyonel destek (tek başına yeterli değil):** R348Q, E257K, G283V ya da L253P için tek moleküllü
hareketlilik ölçümünde hız veya işlemsellik kaybı. Dikkat: KIF1A'da mekanizma heterojen (Boyle 2021
hareketlilik kaybı, Anazawa 2022 baskın negatif, Chiba 2019 aşırı aktivasyon), bu yüzden "hareketlilik
değişmedi" sonucu varyantı temize çıkarmaz.

### 4.2 Öldürecek gözlem

- **Bir sağlıklı erişkin homozigot ya da bileşik heterozigot**, katman A kalıntılarından birinde. Şu an
  katman A'nın hiçbir üyesinde gnomAD'de homozigot yok ve en yüksek AF 9,9e-06 (P303R), yani hiçbiri
  frekansla elenmiyor. Bir homozigot çıkarsa o kalıntı düşer.
- **Büyük bir NDD kohortunda motor alan missense zenginleşmesinin olmaması.** DDD, GeneDx, SPARK ya da
  benzeri bir kohortta KIF13A motor alan de novo missense sayısı mutasyon oranı beklentisini aşmıyorsa,
  hipotez pratik olarak ölür. Bu, tek bir sorguyla yapılabilecek en ucuz ve en belirleyici testtir.
- **K1122E'nin ailede ayrışmaması** (etkilenmemiş ebeveynde bulunması) tek missense ipucunu da götürür;
  bu hipotezi öldürmez ama "KIF13A missense'i nörogelişimle ilişkili" iddiasının elindeki tek fenotipli
  dayanağı alır.
- **Mozaik/LoF vakalarının artması, missense vakalarının artmaması.** Önümüzdeki yıllarda ClinVar'a
  KIF13A LoF vakaları eklenip motor missense eklenmiyorsa, gen NDD geni olarak doğrulanır ama **missense
  ekseni yanlış eksen çıkar**. Bu, bu proje için en olası kötü sonuçtur ve bölüm 3.1'deki kısıtlama
  profili zaten bunu işaret ediyor.

## 5. Adaylar

Ortak sütunlar: `AM` AlphaMissense; `ATP / MT` nükleotide ve tübüline en kısa ağır-atom uzaklığı (Å),
yalnız bindirilmiş motor alan için; `pLDDT` AlphaFold güveni; `KIF1A` o kalıntıda kürate patojenik
varyantlar (`*` = fonksiyonel literatürde çalışılmış). "birebir" = KIF1A'da **aynı** aminoasit değişimi
patojenik. "yapı zayıf" = pLDDT < 70 veya bindirme çifti budandı, yani yapı etiketi tek başına kanıt
sayılmaz. Tam tablo: `data/processed/KIF13A_saturation_lookup_20260916.tsv`.

### 5.1 gnomAD'de güvenilir biçimde gözlenmiş katman A (26)

Bunlar gerçekten var olan alleller. Hiçbiri baskın bir bozukluk için fazla yaygın değil (hepsinin
AF < 1e-05, homozigot yok), yani hiçbiri frekansla elenmiyor; ama görülmüş olmaları patojenik olmadıkları
yönünde zayıf bir kanıttır.

| Değişim | c. | AM | Yapısal sınıf | ATP / MT | pLDDT | KIF1A | Not |
|---|---|---|---|---|---|---|---|
| L370P | c.1109T>C | 0,9999 | boyun bağlayıcı | 35,1 / 4,8 | 79 | L372V | AC 1, yapı zayıf |
| D306H | c.916G>C | 0,9996 | MT arayüzü | 20,3 / 2,1 | 90 | D308G | AC 1 |
| G283V | c.848G>T | 0,9992 | MT arayüzü | 15,7 / 3,6 | 90 | G279R | AC 3 |
| H222R | c.665A>G | 0,9991 | gömülü çekirdek | 7,8 / 11,2 | 87 | H218R | AC 2, **birebir** |
| E257K | c.769G>A | 0,9988 | ATP cebi | 4,9 / 2,6 | 75 | E253D, E253K* | AC 4, **birebir** |
| R258S | c.774A>C | 0,9988 | MT arayüzü | 10,8 / 0,6 | 67 | R254G, R254P, R254Q*, R254W* | AC 1, yapı zayıf |
| A254V | c.761C>T | 0,9980 | ATP cebi | 4,3 / 9,7 | 86 | A250P, A250S, A250D | AC 6 |
| R348G | c.1042C>G | 0,9974 | MT arayüzü | 19,4 / 1,4 | 93 | R350W*, R350G* | AC 1, **birebir** |
| R173G | c.517C>G | 0,9972 | MT arayüzü | 24,8 / 2,4 | 87 | R169W, R169T*, R169K | AC 1 |
| T310I | c.929C>T | 0,9971 | gömülü çekirdek | 17,0 / 8,8 | 94 | T312I, T312P | AC 1, **birebir** |
| S278L | c.833C>T | 0,9968 | gömülü çekirdek | 12,3 / 6,5 | 85 | S274L* | AC 1, **birebir** |
| M214R | c.641T>G | 0,9967 | yüzey | 5,4 / 8,6 | 42 | M210R | AC 1, **birebir**, yapı zayıf |
| R13L | c.38G>T | 0,9963 | ATP cebi | 3,1 / 15,1 | 94 | R13L, R13C, R13S, R13H* | AC 1, **birebir** |
| A73P | c.217G>C | 0,9960 | gömülü çekirdek | 5,5 / 24,2 | 90 | A68E | AC 1 |
| P303R | c.908C>G | 0,9959 | MT arayüzü | 21,5 / 4,9 | 90 | P305L* | AC 16 |
| R220C | c.658C>T | 0,9956 | ATP cebi | 4,3 / 7,3 | 66 | R216H, R216C | AC 3, **birebir**, yapı zayıf |
| H175R | c.524A>G | 0,9948 | MT arayüzü | 29,2 / 3,0 | 76 | H171P | AC 5 |
| R220H | c.659G>A | 0,9935 | ATP cebi | 4,3 / 7,3 | 66 | R216H, R216C | AC 3, **birebir**, yapı zayıf |
| S219N | c.656G>A | 0,9932 | ATP cebi | 2,1 / 11,5 | 62 | S215N, S215R | AC 14, **birebir**, yapı zayıf |
| R11W | c.31C>T | 0,9925 | ATP cebi | 3,1 / 15,4 | 97 | R11L, R11Q*, R11W | AC 1, **birebir** |
| E19A | c.56A>C | 0,9916 | gömülü çekirdek | 8,0 / 15,7 | 86 | E19K | AC 1 |
| R305Q | c.914G>A | 0,9915 | MT arayüzü | 22,6 / 2,4 | 92 | R307P, R307G, R307Q | AC 2, **birebir** |
| D306N | c.916G>A | 0,9914 | MT arayüzü | 20,3 / 2,1 | 90 | D308G | AC 1 |
| Q75H | c.225A>C | 0,9896 | gömülü çekirdek | 5,8 / 23,4 | 93 | Q70K | AC 2 |
| G272S | c.814G>A | 0,9874 | MT arayüzü | 11,8 / 4,3 | 64 | G268W, G268R, G268E | AC 1, yapı zayıf |
| R348Q | c.1043G>A | 0,9750 | MT arayüzü | 19,4 / 1,4 | 93 | R350W*, R350G* | AC 14 |

### 5.2 Kuramsal: gnomAD'de güvenilir biçimde hiç görülmemiş (yeni)

Katman A + **birebir paralog eşi** + yapısal kanıt güvenilir (pLDDT ≥ 70, bindirme çifti korundu).
Böyle 32 değişim var; en yüksek 16'sı aşağıda. 32'nin **27'si eski aday listesinde hiç yoktu**; kalan
5'i (E19K, T104K, R11L, P14T, R11Q) listede yalnızca AC0 etiketiyle, yani güvenilir taşıyıcısı olmadan
duruyordu. **Bir hastada bunlardan biri görülürse, bu belgenin en güçlü öngörüsü doğrulanmış olur.**

| Değişim | c. | AM | Yapısal sınıf | ATP / MT | pLDDT | KIF1A | gnomAD |
|---|---|---|---|---|---|---|---|
| L253P | c.758T>C | 0,9998 | ATP cebi | 4,7 / 12,3 | 92 | L249M, L249P, L249Q | hiç görülmemiş |
| G255R | c.763G>C | 0,9998 | ATP cebi | 3,7 / 6,7 | 81 | G251R* | hiç görülmemiş |
| L279P | c.836T>C | 0,9998 | gömülü çekirdek | 8,2 / 6,0 | 86 | L275Q, L275P | hiç görülmemiş |
| G283R | c.847G>A | 0,9998 | MT arayüzü | 15,7 / 3,6 | 90 | G279R | hiç görülmemiş |
| L253Q | c.758T>A | 0,9996 | ATP cebi | 4,7 / 12,3 | 92 | L249M, L249P, L249Q | hiç görülmemiş |
| R305P | c.914G>C | 0,9994 | MT arayüzü | 22,6 / 2,4 | 92 | R307P, R307G, R307Q | hiç görülmemiş |
| I275S | c.824T>G | 0,9993 | gömülü çekirdek | 8,6 / 6,3 | 72 | I271S | hiç görülmemiş |
| E19K | c.55G>A | 0,9992 | gömülü çekirdek | 8,0 / 15,7 | 86 | E19K | AC0 elenmiş |
| L160P | c.479T>C | 0,9992 | gömülü çekirdek | 11,4 / 16,8 | 86 | L157F, L157P | hiç görülmemiş |
| G203R | c.607G>A | 0,9992 | gömülü çekirdek | 11,9 / 20,9 | 77 | G199R | hiç görülmemiş |
| E151K | c.451G>A | 0,9987 | gömülü çekirdek | 9,2 / 14,0 | 88 | E148K, E148D | hiç görülmemiş |
| R305G | c.913C>G | 0,9987 | MT arayüzü | 22,6 / 2,4 | 92 | R307P, R307G, R307Q | hiç görülmemiş |
| D306G | c.917A>G | 0,9986 | MT arayüzü | 20,3 / 2,1 | 90 | D308G | hiç görülmemiş |
| C97R | c.289T>C | 0,9983 | gömülü çekirdek | 16,2 / 18,2 | 96 | C92R, C92Y | hiç görülmemiş |
| C97Y | c.290G>A | 0,9974 | gömülü çekirdek | 16,2 / 18,2 | 96 | C92R, C92Y | hiç görülmemiş |
| K108T | c.323A>C | 0,9973 | ATP cebi | 2,4 / 10,7 | 92 | K103T | hiç görülmemiş |

**Önemli uyarı: katman A bir kalıntı iddiasıdır, bir değişim iddiası değil.** 929 katman A satırı yalnız
**59 kalıntıya** dağılıyor (kalıntı başına ortalama 15,7 satır). Örneğin R305'in (KIF1A R307) 19
değişiminin tamamı katman A'dadır, çünkü kalıntı hem özdeş hem KIF1A'da patojeniktir ve AlphaMissense o
kalıntıda her değişimi eşiğin üstünde skorlar. Tek nükleotid süzgeci bu şişmeyi kırıyor: R305'te 19
değişimden yalnız 4'ü (R305G, R305L, R305P, R305Q), D306'da 7'si tek bazla erişilebilir. Katman A'nın
erişilebilir 234 satırı 56 kalıntıya karşılık geliyor. Tablodaki satırlar bu yüzden bağımsız hipotezler
değil, kalıntı başına kümelenmiş hipotezlerdir.

### 5.3 Tek fenotipli hasta varyantı

| Değişim | c. | AM | Alan | Paralog | gnomAD | ClinVar |
|---|---|---|---|---|---|---|
| K1122E | c.3364A>G | 0,8044 | kuyruk (520-1805) | **eşleşme yok** | hiç görülmemiş | VCV002570672, VUS, 1 yıldız |

AlphaMissense 0,8044 genel eşiği (0,564) geçiyor ama gen-kalibre eşiği (0,968) geçmiyor, bu yüzden
katmanı `-`. Paralog hizalaması bu pozisyona ulaşmıyor ve kuyruk yapısal olarak **hiç değerlendirilmedi**.
Yani üç eksenden ikisi bu varyant hakkında konuşamıyor, üçüncüsü de "sınırın altında" diyor. Bunu bir
başarısızlık olarak değil, eksenin kapsamının dürüst bir ifadesi olarak yazıyoruz: **bu eksen K1122E'yi
değerlendiremez.** Değerlendirebilecek olan şey aile ayrışması ve fonksiyonel testtir
(`docs/vaka_arama.md`, Lead 1).

## 6. Kısıtlar

- **Eksen bir patojenite kanıtı değil, bir önceliklendirmedir.** 0,968 eşiği KIF1A'da %90 özgüllük
  veriyor, yani her 10 popülasyon kontrolünden 1'i eşiği geçiyor. 1.499 erişilebilir değişimin eşiği
  geçmesi, 1.499 hastalık varyantı olduğu anlamına gelmez.
- **Paralog aktarımı KIF1A'nın mekanizmasını varsayar.** KIF1A'da mekanizma tek yönlü değil (hareketlilik
  kaybı, baskın negatif, aşırı aktivasyon), KIF13A'da hangisinin geçerli olduğu bilinmiyor.
- **Yapı ekseni yalnız 1-380.** Kuyruk (520-1805, protein uzunluğunun %71'i) yapısal olarak
  değerlendirilmedi ve tabloda o sütunlar boştur. Boşluk "yapısal olarak önemsiz" demek değildir.
- **Yapı bir KIF13A yapısı değil**, KIF1A çerçevesine bindirilmiş bir AlphaFold modelidir; yan zincirler
  yeniden optimize edilmedi ve arayüz sayısını gerçek KIF1A ayak izine göre yaklaşık %30 şişiriyor.
  Ayrıntı ve sağlama kontrolleri: `docs/yapi_ekseni.md`.
- **Tek nükleotid erişilebilirliği bir öncelik ölçütü değil, bir gerçekçilik ölçütüdür.** 34.295 değişimin
  yalnız 10.627'si (%31) tek baz değişimiyle olabilir; geri kalanı de novo olarak pratikte görülmez. Ama
  erişilebilir olmak bir varyantı daha zararlı yapmaz.
- **gnomAD'de görülmüş olmak zayıf bir negatif kanıttır, kesin değil.** KIF1A'da p.Ile119Thr 1.165 allel
  ile gnomAD'de duruyor ve AlphaMissense 0,992. Bölüm 5.1'deki "görülmüş" etiketi bu yüzden eleme
  gerekçesi değil, bir bağlam sütunudur.
- **AC0 ayrı bir kategoridir.** 570 protein değişimi gnomAD'de var ama bütün taşıyıcı genotipleri kalite
  filtresinden düştü; bunlar ne "görülmüş" ne "görülmemiş" sayılabilir, tabloda ayrı tutuldular.
- **Hiçbir KIF13A klinik etiketi yok**, bu yüzden bu eksenin KIF13A üzerindeki gerçek performansı
  ölçülemedi ve ölçülemez. Bildirilen bütün AUROC değerleri KIF1A'ya aittir.

## 7. Yeniden üretme

```powershell
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\build_kif13a_candidates.py
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\structure_axis.py
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\build_saturation_table.py
```

Üçüncü betik Ensembl'den MANE CDS'ini bir kez indirir
(`data/raw/sequences/KIF13A_ENST00000259711_cds.fasta`), sonraki çalışmalar önbelleği kullanır.
Bütün çıktılar `atlas_mcp.core.save_tsv` / `save_json` ile yazıldığı için yanlarında `.prov.json`
dosyaları ve `data/raw/PROVENANCE.jsonl` kaydı vardır. Betik sağlama kontrollerinden biri bile geçmezse
sıfırdan farklı bir çıkış kodu döndürür.

## 8. Sırada

- Büyük NDD kohortlarında (DDD, GeneDx, SPARK, Genomics England) KIF13A motor alan de novo missense
  sayımı: bölüm 4.2'deki en belirleyici ve en ucuz test.
- HudsonAlpha ile K1122E üzerine temas ve ISCA site 1 ile 253 kb delesyon üzerine temas
  (`docs/vaka_arama.md`), ikisi de missense ekseninin dışında ama genin durumunu belirler.
- R348Q / E257K / G283V için tek moleküllü hareketlilik testi tasarımı; mekanizma heterojenliği nedeniyle
  hem hız hem işlemsellik hem de aşırı aktivasyon okunmalı.
- Saturasyon tablosunun kuyruk kısmı için bir yapısal eksen: FHA ve CC bölgeleri AlphaFold'da makul
  güvenle modelleniyor, K1122 çevresi ise modellenmiyor.
