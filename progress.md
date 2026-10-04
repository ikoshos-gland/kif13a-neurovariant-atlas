# Öğrenme ilerlemesi (protein-bioinfo öğretme modu)

## 2026-09-15 · Modül 3 + 11 (veritabanları, genetik varyantlar)

Soru: "KIF13A'nın kaç varyantı var, bunu nasıl anladık?"

Plan (her adım bir tur, sonunda bir soru):
1. [anlatıldı] gnomAD nedir, KIF13A için nasıl soruldu, 11.728 nereden geldi (gen sorgusu, satır sayımı)
2. [anlatıldı] Hepsi proteini değiştirmiyor: sonuç türleri (intron 6.089, missense 2.671, sinonim 1.062; MANE listesi 10.936).
   Öğrencinin cevabı: sinonim (doğru, ama "korunma mekanizması" değil kodon dejenerasyonu) + fonksiyonel tolerans
   (farklı bir soru: etki). Eksik: yarıdan fazlası intronda. Düzeltme: gnomAD listesi yalnız ekzon ±75 baz.
   Sorulan: HGVS c. gösteriminden intronik varyant nasıl anlaşılır (c.3545-45A>G vs c.5364G>T)
   Öğrenci "bilmiyorum" dedi -> HGVS c. numaralandırması anlatıldı (c.1 = ATG'nin A'sı, yalnız kodlayan bazlar
   numara alır; intron = numara + işaret + uzaklık). Sorulan: c.3544+10G>A mı c.-14C>T mi intronik (tuzak: UTR eksi işareti)
   Öğrenci "devam et anlatmaya soru sorma" dedi: bundan sonra soru yok, kesintisiz anlatım.
2b. [anlatıldı] Gen listesi 11.728 = 5 transkriptin ekzon ±75 bazı; MANE listesi 10.936; fark 792 (153 missense)
3. [anlatıldı] ClinVar: 313 kayıt, 293 tek-gen, 217 VUS, 61 sınıflandırmasız, 15 B/LB, 0 P/LP
4. [anlatıldı] Teorik: 1805 x 19 = 34.295 (AlphaMissense satırı); CDS'te 16.245 olası SNV, bunlarla ulaşılabilen
   10.621 benzersiz missense; gnomAD'de görülen 2.615 benzersiz missense (%24,6)

Açık sorular: -

## 2026-09-15 · Ölçüm turu

Kalibrasyon setleri ölçüldü (`docs/kalibrasyon_olcumleri.md`): motor içi AUROC AlphaMissense 0,925,
ESM1b 0,851, CADD 0,781; eşit ağırlıklı konsensüs 0,911. AlphaMissense için gen-özgü eşik 0,968.
Sırada: Evo 2 ölçümü (342 NIM çağrısı, onay bekliyor), splice negatif seti, ağırlıklı konsensüs.

## 2026-09-15 · KIF13A aday listesi

`scripts/build_kif13a_candidates.py` + `scripts/figure_candidates.py`:
2.671 missense (2.075 güvenilir gözlem, 596 AC0 filtreli). Kalibre eşik 0,968 -> 295 varyant
(genel eşik 0,564 olsa 982 olurdu). Katmanlar (güvenilir): A 26, B 43, C 110, D 23.
Çıktı: `data/processed/KIF13A_candidates_20260915.tsv`, `figures/KIF13A_candidates_landscape.svg`.
Sırada: Evo 2 ölçümü, yapı ekseni (ATP/mikrotübül arayüzüne uzaklık).

## 2026-09-16 · Beş ajanın sonuçları işlendi

Beş paralel ajan tamamlandı, sonuçları doğrulanıp `HANDOFF.md` ve skill durum dosyasına işlendi.

1. **Yapı ekseni** (`docs/yapi_ekseni.md`): KIF13A AlphaFold modeli KIF1A 8UTN yapısına oturtuldu,
   311 çiftte 1,03 Å RMSD. Katman A: 23 ATP cebi, 18 MT arayüzü, 9 gömülü, 1 yüzey, 1 boyun.
   Sıralamaya girmemiş üçüncü bir kanıt katmanı olarak aynı varyantları işaret ediyor.
2. **Vaka arama** (`docs/vaka_arama.md`): ClinVar'da fenotipli tek kayıt p.Lys1122Glu (VCV002570672),
   ve KIF13A'yı tek protein kodlayan gen olarak içeren 253 kb delesyon (VCV000442594). KIF13A ClinGen,
   G2P, PanelApp ve OMIM'de yok; SysNDD'de "Limited". A254V ve R305L hasta değil, TCGA splicing kayıtları.
3. **Missense sonucu** (`docs/missense_sonuc.md`): 34.295 satırlık doygunluk tablosu üretildi, dokuz
   sağlama geçti. Ulaşılabilir 10.627 değişimden 1.499'u eşiği geçiyor, yalnız 176'sının güvenilir
   gnomAD taşıyıcısı var. 68 değişim KIF1A'da birebir patojenik, 46'sı hiç görülmemiş.
4. **Fonksiyon kaybı** (`docs/fonksiyon_kaybi_ekseni.md`): 204 güvenilir LoF'un 45'i NMD kaçışı ve allel
   yükünün %60'ını taşıyor; gerçek NMD-competent sayı 158. Motor domainde 29-31 tolere edilen kesik
   varyant var, bu katman A'nın gerekçesini baskın-negatif mekanizmaya indiriyor. KIF1A haployetmezlik
   geni DEĞİL (ClinGen HI 0). Hüküm: "olası doz duyarlılığı", kanıtlanmış haployetmezlik değil.
5. **Ölçekleme** (`docs/olcekleme_fizibilite.md`): 3.198 çift puanlandı, %18'i filtreleri geçiyor,
   lift medyanı 0,048, %35,5'inde kazanç yok. Paralog aktarımı genel bir özellik değil. KIF1A/KIF13A
   152/3.198 (üst çeyrek, sıfır zayıflık bayrağı). Sonraki çift önerisi NSD1/NSD3 ve PCDH19/PCDH10.
   Tuzak: ClinVar'da `pathogenic[clinsig]` diye alan yok, Entrez sessizce All Fields'e çeviriyor.

Ayrıca `HANDOFF.md` yazıldı (model bağımsız devir belgesi) ve `~/.claude/skills/kif13a-atlas` skill'i
oluşturuldu.

Kalan boşluklar: hücre tipi belgesi, splice ekseni (26 varyant), git commit yok, paketleme yok.

---

## 16 Eylül 2026, ikinci oturum: doğrulama turu

Kullanıcı "bu hastalarda sadece KIF13A bozuk değil mi, hataya yer yok" dedi. İki insan ipucu, onları
kırmakla görevlendirilmiş ajanlara verildi. **İkisi de en alt kanıt basamağına düştü ve düşen
argümanlardan biri bizimdi.**

### 253 kb delesyon (VCV000442594) -> `docs/dogrulama_delesyon.md`

- gnomAD v4 ekzom CNV `436229__DEL`: 448 kb, KIF13A'nın %99,2'si, 1/464.297. Ayrıca SV setinde
  `DEL_chr6_33713a79` (164 kb, %71,3, PASS, lof) ve 84 kb'lik bir üçüncüsü. ~527.343 kişide en az 3
  bağımsız gen yok edici delesyon, kabaca 1/175.000.
- Kalıtım `unknown`, son değerlendirme 2014, gönderen yorumu şablon, atıf verilen Kaminsky 2011 bu lokusu
  hiç incelememiş (yalnız 14 tekrarlayan bölgeye bakmış).
- Aralıkta 37 enhancer, 9 CTCF, 1 promotor. 52 kb ötedeki NUP153 KIF13A'dan daha kısıtlı (LOEUF 0,277).
- Düzeltmeler: ClinVar'da 20 değil **6** CNV kaybı; patojenik olanlar 23-237 değil **23-382** gen kapsıyor;
  313 kaydın 15 "Pathogenic"inin hepsi Mb ölçekli 6p CNV'si, **0 Likely pathogenic**, 0 gen içi patojenik.
- Güç uyarısı: gnomAD v4'te alt kohort **yok** (controls/non_neuro yalnız v2.1'de). v2.1'de 1/175.000 bir
  olay için beklenen sayı 0,03-0,06, yani oradaki sıfır güç yetersizliğidir, kanıt değil.
- dbVar kullanılamaz: nstd101 obsolete, yerine geçen nstd102 Case-Set ve 21/21 örnek alanı boş.
- DGV Gold: tam gen delesyonu yok (iddianın lehine tek bulgu).

### DECIPHER (aynı ajanın ikinci turu)

- **Hasta 451718**: 6:17.825.889-17.972.163 GRCh38, 146.275 bp, **tamamen gen içinde**, genin %63,6'sı,
  array mean ratio -0,999, Uncertain, kalıtım bilinmiyor, fenotip "Autistic behavior", 1 varyant 1 fenotip.
  Ensembl ile hem GRCh38 hem GRCh37'de doğrulandı: aralıkta yalnız 1 lncRNA ve 1 psödogen var.
  Tek genlik bu vakada geometrinin sonucu, veri tabanı iddiası değil.
- **Hasta 504413**: 182.672 bp, genin %73,0'ü, **maternal**, **likely benign**. Doğrudan karşı örnek.
- DDD araştırma verisi: KIF13A'da tam 2 varyant (stop_gained + frameshift), **ikisi de maternal, ikisi de
  Uncertain**.
- **DDG2P'de KIF13A yok** (2.869 kayıt, 12 kinesin var). KIF1A'nın iki kaydı: monoallelic = "altered gene
  product structure", biallelic = "absent gene product / LoF". Yani kardeş gende bile tek kopya kaybı
  hastalık yapmıyor.
- Haployetmezlik tahmincileri çelişiyor: pHaplo 0,99 lehte, sHet 0,049 ve %HI 42,3 aleyhte.

**Toplam desen: kalıtımı bilinen 4 LoF olayının 3'ü anneden. Tek de novo olan Lam 2021 ve o da mozaik.**

### p.Lys1122Glu (VCV002570672) -> `docs/dogrulama_hudsonalpha.md`

- `StudyName = CSER-HudsonAlpha`, 384 varyantlık partinin 124 VUS'undan biri. Partide 152 varyant de novo
  işaretli, bizimki `Origin: unknown`, `AffectedStatus: unknown`, gönderen yorumu boş.
- **Bizim hatamız:** "1.293 gönderim içinde KIF13A'da tam 1 kayıt, demek ki bilinçli aday çağrısı"
  çıkarımı yanlış. Tek gen kayıtlarındaki 607 genin 368'inde (%60,6) tam olarak 1 kayıt var; bu onların
  **en sık** deseni. Payda ölçülmeden oran kanıt değildir.
- Kalıntı 1122 gerçekten korunmuş (203 ortologun 186'sında Lys) ama coiled coil 1100-1138 içinde ve o
  bölgede gnomAD'de AC>0 olan 44 missense var, **p.Lys1104Glu AC 1 dahil (birebir aynı Lys->Glu)**.
- AlphaMissense 0,8044, kalibre eşiğin (0,968) 0,164 altında. MPC 0,573, REVEL 0,642, PrimateAI tolere.
- Chiari muhtemelen tesadüfi: herhangi bir endikasyonla görüntülenen 14.116 çocuğun %3,6'sı; gelişim
  geriliği kohortunda 6/6 etiyolojiyle ilgisiz.
- **Belge düzeltmesi:** `docs/fonksiyon_kaybi_ekseni.md` 4.4'teki "etkilenme durumu present" yanlış;
  XML'de `AffectedStatus = unknown`, "present" ayrı bir alandan (`ClinicalFeaturesAffectedStatus`) geliyor.

### Bu oturumda kendim ölçtüğüm: LoF toleransı karşılaştırma paneli

12 ClinGen HI-3 nörogelişimsel geniyle aynı boru hattı. KIF13A'nın gözlenen/beklenen LoF değeri **0,372**,
panelin hepsinden yüksek (0,037-0,287). ADNP'de NMD-yetkin HC LoF sayısı **sıfır**. Ham taşıyıcı frekansı
tek yüksek frekanslı allellerle kirleniyor (ARID1B c.611_612insGC AC 798, SCN1A c.1086C>A AC 773,
KIF1A c.4257-2A>C AC 824), o yüzden yalnız obs/exp okunmalı. Tablo skill'de `proje_durumu.md` 9d.
Lam 2021 canlı doğrulandı: ekzom **ve** genom yapılmış, de novo postzigotik mozaik, ama 5 yılda 3 atıf ve
hiçbiri bağımsız bir KIF13A vakası değil.

### Ortam

`mcp\atlas_server\.venv` arkasındaki uv CPython 3.12.14 Windows tarafından engellendi, venv hiç
başlamıyor. Sistem `python` 3.13.14 doğrulandı: httpx, numpy, pyyaml, biopython, pandas, scipy, mcp var,
`atlas_mcp.core` içe aktarılıyor, `.env` okunuyor. `CLAUDE.md` ve skill güncellendi; `HANDOFF.md`'deki
komutlar ve MCP kayıtları hâlâ ölü venv'i gösteriyor.

### Devam eden ve sıradaki

Komisyonlandı, sonucu bekleniyor: `docs/denovo_ekseni.md` (31.058 gelişimsel bozukluk üçlüsünde de novo
sayımı, ANKRD11/CHD8 pozitif kontrol, KIF13B null) ve `docs/bolgesel_kisitlilik.md` (motor domaine yerel
missense kısıtlılığı, KIF1A pozitif kontrol).

Yalnız bir insanın yapabileceği tek iş: DECIPHER 451718'in yatıran merkezine ve HudsonAlpha'ya
(örnek A00166) **kalıtım** sormak. İkisi de bir basamak yükseltebilir.

Skill baştan güncellendi (1.034 satır): `proje_durumu.md` 6, 6b ve 9d bölümleri, `tuzaklar.md` artık 45
tuzak (26-39 kanıt değerlendirme, 40-45 ortam), `yontem.md`'ye Adım 8 eklendi, `sonraki_adimlar.md`
yeniden yazıldı, `veri_katmani.md`'ye MCP'siz çalışan kaynaklar eklendi.

## 16 Eylül 2026, de novo ekseni tamamlandı

`docs/denovo_ekseni.md` + `scripts/denovo_axis.py` + `data/interim/denovo/` (11 TSV/JSON, hepsi
provenanslı). **Hüküm: KARŞI KANIT.**

Dört kaynak canlı bağlandı: Kaplanis 2020 (PMID 33057194, 31.058 DD triosu, ek tabloları Europe PMC
PMC7116826'dan), Gene4Denovo2 (PMID 41024706, arayüzünün kendi JSON API'si bulundu ve kullanıldı),
denovo-db v1.6.1 (PMID 27907889, SSC olmayan kısım), gnomAD v2.1.1 mutasyon hızları. SFARI Gene
okundu. DECIPHER CNV verisi anahtar istiyor, durduk; denovo-db SSC kısmı kullanım kısıtı beyan
ediyor, indirmedik.

**KIF13A: etkilenmiş hiçbir bireyde de novo LoF yok.** Kaplanis'te 0 (beklenen 0,367, denovoWEST
p=0,709); tüm katalogların toplamında 0 (beklenen 1,16, Gene4Denovo2'nin kendi p'si 0,975);
denovo-db'de 0. Bulunan tek de novo LoF `p.N1649Kfs*3`, SPARK'ın **sağlıklı kardeş** kolunda ve
ekzon 39'da (39 ekzonun son'u, Ensembl'den doğrulandı), yani NMD kaçış bölgesinde: iki ayrı nedenle
doz kanıtı değil. Missense de beklentinin altında (2'ye karşı 3,16).

Boru hattı doğrulandı: ANKRD11 111 LoF (beklenen 0,308), CHD8 19 (beklenen 0,497), toplamda 148x ve
47x zenginleşme. KIF1A tam beklenen deseni verdi: **53 de novo missense, 0 de novo LoF**, yani
"hastalık yapar ama LoF ile değil" durumunun nasıl göründüğü ölçüldü. Etkin trio sayısı sinonim de
novo'lardan kalibre edildi (34 sinonim / 2×mu_syn = 97.700, Gene4Denovo2'nin beyan ettiği 130.439
bireyle tutarlı).

Nicel sınır: sıfır gözlem %95 güvenle Kaplanis'te 8,2x, toplamda **2,6x** üzerini dışlıyor. Yani
ANKRD11/CHD8 tipi haployetmezlik elendi, hafif doz duyarlılığı elenmedi.

### Bu eksende bulunan yeni tuzaklar (HANDOFF 16-23)

Kaplanis ek tablosu **GRCh37**'dir ve bunu yazmaz (test edildi, 4/4; GRCh38 sanılsa dört KIF13A
satırının da sınıfı değişirdi). Gene4Denovo2'nin `AAChange_refGene` sütunu KIF13A için NM_001105567/8
kullanıp **13 kalıntı düşük** numara veriyor, doğrusu `AAChange_mrefGene` (NM_022113): 1770
tuzağından farklı, üçüncü bir izoform kayması. denovo-db her varyantı her transkript için ayrı satır
yazıyor (KIF13A'nın "101 satırı" 25 varyant). Gene4Denovo2'nin LoF etiketi güvenilmez (KIF1A'da 3
bazlık çerçeve içi insersiyon `stopgain` etiketli). `1 - CDF` Poisson'u pozitif kontrollerde taşıp
sahte `p = 0` basıyor, üst kuyruk log uzayında toplanmalı.

**Projeyi doğrudan etkileyen düzeltme: VCV000442594'ün ClinVar `Origin` alanı "unknown".** 253
kb'lik delesyon belgelerde sık sık doz kanıtı olarak anılıyor ama **de novo olduğu kayıtlı değil**;
de novo argümanına giremez.

Motor domainde üç gerçek insan de novo missense'i bulundu: **D65E** (Kaplanis, DD), **R173Q**
(Yuen 2017 + Cirnigliaro 2023, ASD, bir mi iki birey olduğu ayırt edilemedi), **Q237H** (Chesi 2013,
**ALS**). Üçü de doygunluk tablosuyla kesiştirilmeyi bekliyor; bu tablonun ilk gerçek sınavı.

Ortam: venv hâlâ ölü (Smart App Control), sistem Python 3.13 ile üretildi, `atlas_mcp.core` sorunsuz.

## 17 Eyl 2026: NSD1/NSD3 ile yöntemin örneklem dışı doğrulaması

Yöntem makalesinin merkezi iddiası ("paralog aktarımı genel değil ama hangi çiftte işe yarayacağı
ucuz ölçütlerle kestirilebilir") ilk kez örneklem dışı sınandı. Tarama NSD1/NSD3'ü 3.198 çiftin
1. sırasına koymuştu; boru hattı bu çift üzerinde, tahmin satırı kapalı tutularak baştan sona
çalıştırıldı. Betik `scripts/nsd_pair_calibration.py` (10 aşama), belge
`docs/nsd1_nsd3_dogrulama.md`, veri `data/interim/nsd_pair/`, aday listesi
`data/processed/NSD3_candidates_20260917.tsv`, kimlikler `config/targets.yaml`.

Hüküm: **öngörü tuttu.** Paylaşılan katalitik modül (NSD1 1539-2210, PHD1-CysHis, 672 aa, özdeşlik
0,646) içinde AlphaMissense AUROC 0,973 (0,942-0,995), %90 özgüllükte duyarlılık 0,95, gen-özgü eşik
0,9732; KIF1A'nın motor domainindeki 0,925 / 0,79 / 0,968 değerlerinin üzerinde ve aynı yön CADD
(0,781 -> 0,820), ESM1b (0,851 -> 0,902) ve PrimateAI (0,801 -> 0,857) ile de görülüyor. NSD3
aktarımı: 1.336 gnomAD missense, genel eşik 399 varyant çağırırken kalibre eşik 101; katmanlar
A 23 (19 güvenilir gözlemli), B 45 (43), C 33 (27), D 30 (22).

Kürasyon: ClinVar 2.694 tek gen kaydı, 158 P/LP missense, 125 kullanılabilir çekirdek pozitif
(87 kalıntı, 120'si modül içinde), 26 PP3'e dayalı pozitif ayrı tutuldu, negatifler modül içi gnomAD
AC >= 3 (72; AC = 2 için 47). ClinVar B/LB missense'in yalnız 6'sı modül içinde, yani konum
karışıklığı KIF1A'da olduğu gibi burada da geçerli ("modül içinde mi" tüm protein AUROC'u 0,964).

Tutmayan kısım: taramanın sıralamada kullandığı **kazanç** (lift) ölçütü. NSD1/NSD3'ün kazancı
KIF1A/KIF13A'nın iki katı (0,250'ye karşı 0,128) ama ölçülen paralog-özdeşlik referans çizgisi
ayırt edilemez (0,740'a karşı 0,744). Ayrıca kardeşin kısıtlılığı taramada gen düzeyinde ölçülüyor:
bölgesel hesap NSD3'ün paylaşılan modülünün proteinin geri kalanından daha kısıtlı **olmadığını**
gösteriyor (o/e missense 0,760'a karşı modül dışı 0,774), NSD1'de ise modül 0,433 ve modül dışı
1,032. Öğretmen/öğrenci derinlik oranı iki çiftte neredeyse aynı (1,75 ve 1,81 kat).

Sonraki adım (makale için zorunlu): taramanın **başarısız olacağını söylediği** bir çiftin aynı boru
hattından geçirilmesi. İki vaka da üst %5'te; yanlış-pozitif kontrolü olmadan kestirim iddiası
yalnız üst uçta desteklenmiş durumda.

Yeni tuzaklar: ClinVar'ın `pathogenic[clinsig]` şişmesi NSD1'de bağımsız olarak yeniden üretildi
(doğru Properties biçimi 158, yanlış biçim 371). NSD3'ün eski sembolü WHSC1L1'e Ensembl
`lookup/symbol` 400 döner. NSD3'ün çok çalışılan kısa izoformu Q9BZ95-3 (645 aa) SET domaini
taşımaz, numaralandırması MANE ile karşılaştırılamaz. Ensembl REST bu oturumda aralıklı 500/503
döndürdü, Ensembl çağrıları 9 denemeye çıkarıldı.

## 3-4 Ekim 2026: üç TRUBA ekseni, AlphaGenome değerlendirmesi, ilk commit

Üç ajan TRUBA'da paralel çalıştı (kullanım limiti yüzünden ikişer kez kesildiler; TRUBA işleri sürdü).

- **B, insan yaşam boyu** (`docs/insan_yasam_boyu.md`): faredeki yaşla artış insanda yok (335 donör,
  güç >0,99); oligodendrosit soyu birinciliği tekrarlandı ama insanda OPC birinci; OPC'de erken
  çocuklukta artış var ama KIF13A'ya özgü değil.
- **C, nakavt ve popEVE** (`docs/nakavt_ekseni.md`): bilgi vermiyor (KIF1A'yı bile ayırt edemiyor);
  RGC-ME'de 2 homozigot, büyük olasılıkla NMD kaçışı; popEVE KIF1A ile birebir aynı 68 patojenik
  değişimin hiçbirini KIF13A'da ağır saymıyor.
- **A, hücre tipi eQTL kolokalizasyonu** (`docs/eqtl_hucre_tipi.md`, ön sonuç): ak madde haplotipi
  oligodendrosit ve OPC'deki KIF13A eQTL'iyle kolokalize değil (PP.H3 0,78-1,00); mikroglia ipucu LD
  paneline duyarlı. Ajan kullanıcı tarafından durduruldu, belge tablolardan ana oturumda yazıldı.
- AlphaGenome Atlas / AVI incelendi (HANDOFF 3.15): splice ekseni için yararlı, motor missense için değil.
- Synapse token'ı kullanıcı TRUBA'ya kendisi koydu; ilk token yalnız `view` yetkiliydi (tuzak).
- Skill güncellendi (77 tuzak, skor tablosu) ve depoya `skill/kif13a-atlas/` olarak kopyalandı.
- Depo ilk kez commit edildi ve GitHub'a gönderildi; `data/processed` ve `data/interim` dahil,
  `data/raw` hariç. ClinVar başvuru sahiplerinin e-posta adresleri depodaki dosyalardan çıkarıldı.
